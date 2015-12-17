/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package picard.sam.markduplicates;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.SortingCollection;
import picard.util.BatchBuffer;
import picard.sam.markduplicates.util.OpticalDuplicateFinder;

import java.util.List;
import java.util.Map;

import static picard.sam.markduplicates.EstimateLibraryComplexity.*;

/**
 *@author Pavel Silin, Pavel_Silin@epam.com
 **/

public class PairReadSequenceAssembler implements Runnable {

    private final Log log = Log.getInstance(PairReadSequenceAssembler.class);
    private Map<String, PairedReadSequence> pendingByName;
    private List<SAMReadGroupRecord> readGroups;
    private SortingCollection<PairedReadSequence> sorter;
    private ProgressLogger progress;
    private final boolean useBarcodes;
    private OpticalDuplicateFinder opticalDuplicateFinder;
    private BatchBuffer<SAMRecord> batchBuffer;
    int minIdenticalBases;
    int minMeanQuality;
    String barcodeTag;
    String readOneBarcodeTag;
    String readTwoBarcodeTag;
    int readSize;

    /**
     * @param batchCapacity
     * @param pendingByName
     * @param readGroups
     * @param sorter
     * @param progress
     * @param useBarcode
     * @param opticalDuplicateFinder
     * @param minIdenticalBases
     * @param minMeanQuality
     * @param barcodeTag
     * @param readOneBarcodeTag
     * @param readTwoBarcodeTag
     * @param readSize
     */
    public PairReadSequenceAssembler(int batchCapacity, Map<String, PairedReadSequence> pendingByName,
                                     List<SAMReadGroupRecord> readGroups, SortingCollection<PairedReadSequence> sorter,
                                     ProgressLogger progress, boolean useBarcode,
                                     OpticalDuplicateFinder opticalDuplicateFinder, int minIdenticalBases,
                                     int minMeanQuality, String barcodeTag, String readOneBarcodeTag, String
                                             readTwoBarcodeTag, int readSize) {
        this.batchBuffer = new BatchBuffer<SAMRecord>(SAMRECORDS_PACK_SIZE, batchCapacity);
        this.pendingByName = pendingByName;
        this.readGroups = readGroups;
        this.sorter = sorter;
        this.progress = progress;
        this.useBarcodes = useBarcode;
        this.opticalDuplicateFinder = opticalDuplicateFinder;
        this.minIdenticalBases = minIdenticalBases;
        this.minMeanQuality = minMeanQuality;
        this.barcodeTag = barcodeTag;
        this.readOneBarcodeTag = readOneBarcodeTag;
        this.readTwoBarcodeTag = readTwoBarcodeTag;
        this.readSize = readSize;
    }

    @Override
    public void run() {
        while (true){
            try {
                for (final SAMRecord rec : batchBuffer.take()) {
                    if (rec == BatchBuffer.FINISH_RECORD) {
                        return;
                    }
                    PairedReadSequence prs = pendingByName.remove(rec.getReadName());

                    if (prs == null) {
                        prs = useBarcodes ? new PairedReadSequenceWithBarcodes(readSize) : new PairedReadSequence(readSize);
                        if (opticalDuplicateFinder.addLocationInformation(rec.getReadName(), prs)) {
                            final SAMReadGroupRecord rg = rec.getReadGroup();
                            if (rg != null) prs.setReadGroup((short) readGroups.indexOf(rg));
                        }

                        pendingByName.put(rec.getReadName(), prs);
                    }

                    // Read passes quality check if both ends meet the mean quality criteria
                    final boolean passesQualityCheck = passesQualityCheck(rec.getReadBases(), rec.getBaseQualities(),
                            minIdenticalBases, minMeanQuality);
                    prs.qualityOk = prs.qualityOk && passesQualityCheck;

                    // Get the bases and restore them to their original orientation if necessary
                    final byte[] bases = rec.getReadBases();
                    if (rec.getReadNegativeStrandFlag()) SequenceUtil.reverseComplement(bases);

                    final PairedReadSequenceWithBarcodes prsWithBarcodes = (useBarcodes)
                            ? (PairedReadSequenceWithBarcodes) prs : null;

                    if (rec.getFirstOfPairFlag()) {
                        prs.read1 = bases;
                        if (useBarcodes) {
                            prsWithBarcodes.barcode = getReadBarcodeValue(rec, barcodeTag);
                            prsWithBarcodes.readOneBarcode = getReadBarcodeValue(rec, readOneBarcodeTag);
                        }
                    } else {
                        prs.read2 = bases;
                        if (useBarcodes) {
                            prsWithBarcodes.readTwoBarcode = getReadBarcodeValue(rec, readTwoBarcodeTag);
                        }
                    }

                    if (prs.read1 != null && prs.read2 != null && prs.qualityOk) {
                        sorter.add(prs);
                    }
                    progress.record(rec);
                }
            } catch (InterruptedException e) {
                log.error("Failed to process read assembler", e);
                break;
            }
        }
    }

    /**
     * Checks that the average quality over the entire read is >= min, and that the first N bases do
     * not contain any no-calls.
     */
    private boolean passesQualityCheck(final byte[] bases, final byte[] quals, final int seedLength,
        final int minQuality) {
        if (bases.length < seedLength) return false;

        for (int i = 0; i < seedLength; ++i) {
            if (SequenceUtil.isNoCall(bases[i])) return false;
        }

        int total = 0;
        for (final byte b : quals) total += b;
        return total / quals.length >= minQuality;
    }


    public void add(SAMRecord rec) throws InterruptedException {
        batchBuffer.add(rec);
    }

    public void finish() throws InterruptedException {
        batchBuffer.add(BatchBuffer.FINISH_RECORD);
    }

}
