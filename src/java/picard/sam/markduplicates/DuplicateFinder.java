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
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.PeekableIterator;
import picard.sam.markduplicates.EstimateLibraryComplexity.PairedReadSequence;
import picard.sam.markduplicates.EstimateLibraryComplexity.PairedReadSequenceWithBarcodes;
import picard.sam.markduplicates.util.OpticalDuplicateFinder;
import picard.util.ArraysUtil;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.concurrent.BlockingQueue;

/**
 *@author Pavel Silin, Pavel_Silin@epam.com
 **/

public class DuplicateFinder implements Runnable {

    /**
     *This parameter determines the choice of the algorithm: if the group size > MIDDLE_SIZE_LIBRARY, the modified
     * algorithm applies
     */
    private static final int MIDDLE_SIZE_LIBRARY = 100;

    public static final String UNKNOWN_LIBRARY = "Unknown";
    private static final List FINISH_SIGNAL_FOR_DUPLICATION_FINDER = Collections.EMPTY_LIST;

    private final Log log = Log.getInstance(DuplicateFinder.class);

    private final int minIdenticalBases;
    private final double maxDiffRate;
    private final BlockingQueue<List<PairedReadSequence>> groupsQueue;
    private final Map<String, Histogram<Integer>> duplicationHistosByLibrary;
    private final Map<String, Histogram<Integer>> opticalHistosByLibrary;
    private int lengthOfShingle;
    private int numberOfBandsRead1;
    private int numberOfBandsRead2;
    private int groupsProcessed;
    private final List<SAMReadGroupRecord> readGroups;
    private OpticalDuplicateFinder opticalDuplicateFinder;
    private final boolean useBarcode;

    public DuplicateFinder(List<SAMReadGroupRecord> readGroups, BlockingQueue<List<PairedReadSequence>> groupsQueue,
                           OpticalDuplicateFinder opticalDuplicateFinder, final boolean useBarcode,
                           int groupsProcessed, int minIdenticalBases, double maxDiffRate) {
        this.readGroups = readGroups;
        this.groupsQueue = groupsQueue;
        this.groupsProcessed = groupsProcessed;
        this.useBarcode = useBarcode;
        this.opticalDuplicateFinder = opticalDuplicateFinder;
        this.minIdenticalBases = minIdenticalBases;
        this.maxDiffRate = maxDiffRate;
        this.duplicationHistosByLibrary = new HashMap<String, Histogram<Integer>>();
        this.opticalHistosByLibrary = new HashMap<String, Histogram<Integer>>();
    }

    public Map<String, Histogram<Integer>> getDuplicationHistosByLibrary(){
        return duplicationHistosByLibrary;
    }

    public Map<String, Histogram<Integer>> getOpticalHistosByLibrary(){
        return opticalHistosByLibrary;
    }

    @Override
    public void run() {
        long lastLogTime = System.currentTimeMillis();

        while (true) {
            try {
                List<PairedReadSequence> group = groupsQueue.take();

                if (isFinish(group)) {
                    return;
                }

                final Map<String, List<PairedReadSequence>> sequencesByLibrary = splitByLibrary(group, readGroups);

                // Now process the reads by library
                for (final Map.Entry<String, List<PairedReadSequence>> entry : sequencesByLibrary.entrySet()) {
                    final String library = entry.getKey();
                    final List<PairedReadSequence> seqs = entry.getValue();

                    Histogram<Integer> duplicationHisto = duplicationHistosByLibrary.get(library);
                    Histogram<Integer> opticalHisto = opticalHistosByLibrary.get(library);

                    if (duplicationHisto == null) {
                        duplicationHisto = new Histogram<Integer>("duplication_group_count", library);
                        opticalHisto = new Histogram<Integer>("duplication_group_count", "optical_duplicates");
                        duplicationHistosByLibrary.put(library, duplicationHisto);
                        opticalHistosByLibrary.put(library, opticalHisto);
                    }

                    if (seqs.size() < MIDDLE_SIZE_LIBRARY || useBarcode) {
                        searchSmallLibrary(seqs, duplicationHisto, opticalHisto, useBarcode);
                    } else {
                        fillHashValues(seqs);
                        searchBigLibrary(seqs, duplicationHisto, opticalHisto);

                    }
                }

                ++groupsProcessed;
                if (log.isEnabled(Log.LogLevel.INFO) && lastLogTime < (System.currentTimeMillis() - 60000)) {
                    log.info("Processed " + groupsProcessed + " groups.");
                    lastLogTime = System.currentTimeMillis();
                }
            } catch (InterruptedException e) {
                log.error("Failed to find duplicates", e);
                break;
            }
        }
    }

    public void initAlgorithmProperties(PeekableIterator<PairedReadSequence> iterator) {
        if (!iterator.hasNext()){
            return;
        }
        PairedReadSequence pairedRead = iterator.peek();
        int maxReadLength = Math.max(pairedRead.read1.length - minIdenticalBases,
                pairedRead.read2.length - minIdenticalBases);

        //if read.length % shingle.length != 0, we have tail of the read which is not included in the hashValues,
        // and we need extra hash value to get prs.hashValues.length = (number of errors) + 1
        int numberOfBands = (int) (maxReadLength * maxDiffRate) + 1;
        lengthOfShingle = maxReadLength / numberOfBands;
    }

    /**
     *Put in the queue EMPTY_LIST, indicating that it is time to shut down.
     */
    public void finish() throws InterruptedException {
        this.groupsQueue.put(FINISH_SIGNAL_FOR_DUPLICATION_FINDER);
    }

    private boolean isFinish(List<PairedReadSequence> group) {
        return group.equals(FINISH_SIGNAL_FOR_DUPLICATION_FINDER);
    }

    private void searchBigLibrary(List<PairedReadSequence> seqs, Histogram<Integer> duplicationHisto,
                                  Histogram<Integer> opticalHisto) {
        sortPairedReadsOnBands(seqs);

        HashSet<PairedReadSequence> somethingDuplicate = new HashSet<PairedReadSequence>();
        for (PairedReadSequence lhs : seqs) {
            if (lhs == null) continue;

            final List<PairedReadSequence> dupes = new ArrayList<PairedReadSequence>();

            for (PairedReadSequence rhs : lhs.possibleCopies) {
                if (somethingDuplicate.contains(rhs)) continue;

                if (!lhs.equals(rhs) && matchesByShingles(lhs, rhs, this.maxDiffRate)) {
                    dupes.add(rhs);
                    int index = rhs.indexInLibrary;
                    seqs.set(index, null);
                }

            }

            somethingDuplicate.addAll(dupes);
            fillHistogram(duplicationHisto, opticalHisto, lhs, dupes);
        }
    }

    private void sortPairedReadsOnBands(List<PairedReadSequence> seqs){
        Map<Integer, ArrayList<PairedReadSequence>> band = new HashMap<Integer, ArrayList<PairedReadSequence>>();

        for (PairedReadSequence prs : seqs) {
            putPairedReadOnBand(band, prs, prs.hashValuesRead1);
            putPairedReadOnBand(band, prs, prs.hashValuesRead2);
        }

        for (ArrayList<PairedReadSequence> copies : band.values()) {
            for (PairedReadSequence prs : copies) {
                prs.possibleCopies.addAll(copies);
            }
        }
    }

    private void putPairedReadOnBand(Map<Integer, ArrayList<PairedReadSequence>> band,
                                     PairedReadSequence prs, int[] readHashValue) {
        for (int key: readHashValue) {
            ArrayList<PairedReadSequence> possibleCopies = band.get(key);

            if (possibleCopies == null) {
                possibleCopies = new ArrayList<PairedReadSequence>();
                band.put(key, possibleCopies);
            }
            possibleCopies.add(prs);
        }
    }

    private boolean matchesByShingles(PairedReadSequence lhs, PairedReadSequence rhs, double maxDiffRate) {
        final int read1Length = Math.min(lhs.read1.length, rhs.read1.length);
        final int read2Length = Math.min(lhs.read2.length, rhs.read2.length);
        final int maxErrors = (int) Math.floor((read1Length + read2Length) * maxDiffRate);
        int errors = 0;

        errors += compareReadToRead(new HashedPairedReadEnd(lhs.read1, lhs.hashValuesRead1),
                new HashedPairedReadEnd(rhs.read1, rhs.hashValuesRead1), maxErrors);
        if (errors > maxErrors){
            return false;
        }

        errors += compareReadToRead(new HashedPairedReadEnd(lhs.read2, lhs.hashValuesRead2),
                new HashedPairedReadEnd(rhs.read2, rhs.hashValuesRead2), maxErrors);
        if (errors > maxErrors){
            return false;
        }

        return true;
    }

    private int compareReadToRead(HashedPairedReadEnd read1, HashedPairedReadEnd read2, int maxErrors){
        int errors = 0;
        final int read1Length = Math.min(read1.read.length, read2.read.length);
        int[] minHashValueRead = read1.shingleHashValue.length > read2.shingleHashValue.length ? read2.shingleHashValue:
                read1.shingleHashValue;

        for(int k = 0; k < minHashValueRead.length; k++){
            if (read1.shingleHashValue[k] != read2.shingleHashValue[k]){
                errors += compareShingles(read1.read, read2.read, getReadOffset(k), getReadOffset((k+1)));
                if(errors > maxErrors){
                    return errors;
                }
            }
        }

        if(read1Length > getReadOffset(minHashValueRead.length)){
            errors += compareShingles(read1.read, read2.read, getReadOffset(minHashValueRead.length), read1Length);
            if(errors > maxErrors){
                return errors;
            }
        }

        return errors;
    }

    private int getReadOffset(int k) {
        return k*lengthOfShingle + minIdenticalBases;
    }

    private int compareShingles(byte[] read1, byte[] read2, int start, int stop) {
        int errors = 0;
        for (int i = start; i < stop; ++i) {
            if (read1[i] != read2[i]) {
                errors++;
            }
        }
        return errors;
    }

    private void searchSmallLibrary(List<PairedReadSequence> seqs, Histogram<Integer> duplicationHisto,
                                    Histogram<Integer> opticalHisto, final boolean useBarcodes) {
        for (int i = 0; i < seqs.size(); ++i) {
            final PairedReadSequence lhs = seqs.get(i);
            if (lhs == null) continue;
            final List<PairedReadSequence> dupes = new ArrayList<PairedReadSequence>();

            for (int j = i + 1; j < seqs.size(); ++j) {
                final PairedReadSequence rhs = seqs.get(j);
                if (rhs == null) continue;

                if (matches(lhs, rhs, this.maxDiffRate, useBarcodes)) {
                    dupes.add(rhs);
                    seqs.set(j, null);
                }
            }

            fillHistogram(duplicationHisto, opticalHisto, lhs, dupes);
        }
    }

    private void fillHashValues(List<PairedReadSequence> sequences) {
        for (PairedReadSequence prs : sequences) {
            int amountShinglesInRead1 = (prs.read1.length - minIdenticalBases) / lengthOfShingle;
            int amountShinglesInRead2 = (prs.read2.length - minIdenticalBases) / lengthOfShingle;

            prs.possibleCopies = new HashSet<PairedReadSequence>();

            prs.hashValuesRead1 = new int[amountShinglesInRead1];
            for (int i = 0; i < amountShinglesInRead1; i++) {
                int st = getReadOffset(i);
                int end = getReadOffset((i + 1));
                prs.hashValuesRead1[i] = ArraysUtil.hash(prs.read1, st, end);
            }

            prs.hashValuesRead2 = new int[amountShinglesInRead2];
            for (int j = 0; j < amountShinglesInRead2; j++) {
                int st = getReadOffset(j);
                int end = getReadOffset((j + 1));
                prs.hashValuesRead2[j] = ArraysUtil.hash(prs.read2, st, end);
            }
        }
    }

    private void fillHistogram(Histogram<Integer> duplicationHisto, Histogram<Integer> opticalHisto,
                               PairedReadSequence lhs, List<PairedReadSequence> dupes) {
        if (dupes.size() > 0) {
            dupes.add(lhs);
            final int duplicateCount = dupes.size();
            duplicationHisto.increment(duplicateCount);

            final boolean[] flags = opticalDuplicateFinder.findOpticalDuplicates(dupes);
            for (final boolean b : flags) {
                if (b) opticalHisto.increment(duplicateCount);
            }
        } else {
            duplicationHisto.increment(1);
        }
    }

    /**
     * Checks to see if two reads pairs have sequence that are the same, give or take a few
     * errors/diffs as dictated by the maxDiffRate.
     */
    private boolean matches(final PairedReadSequence lhs, final PairedReadSequence rhs, final double maxDiffRate,
                            final boolean useBarcodes) {

        final int read1Length = Math.min(lhs.read1.length, rhs.read1.length);
        final int read2Length = Math.min(lhs.read2.length, rhs.read2.length);
        final int maxErrors = (int) Math.floor((read1Length + read2Length) * maxDiffRate);
        int errors = 0;

        if (useBarcodes) {
            final PairedReadSequenceWithBarcodes lhsWithBarcodes = (PairedReadSequenceWithBarcodes) lhs;
            final PairedReadSequenceWithBarcodes rhsWithBarcodes = (PairedReadSequenceWithBarcodes) rhs;
            if (lhsWithBarcodes.barcode != rhsWithBarcodes.barcode ||
                    lhsWithBarcodes.readOneBarcode != rhsWithBarcodes.readOneBarcode ||
                    lhsWithBarcodes.readTwoBarcode != rhsWithBarcodes.readTwoBarcode) {
                return false;
            }
        }

        // The loop can start from minIdenticalBases because we've already confirmed that
        // at least those first few bases are identical when sorting.
        for (int i = this.minIdenticalBases; i < read1Length; ++i) {
            if (lhs.read1[i] != rhs.read1[i]) {
                if (++errors > maxErrors) return false;
            }
        }

        for (int i = this.minIdenticalBases; i < read2Length; ++i) {
            if (lhs.read2[i] != rhs.read2[i]) {
                if (++errors > maxErrors) return false;
            }
        }

        return true;
    }

    /**
     * Takes a list of PairedReadSequence objects and splits them into lists by library.
     */
    Map<String, List<PairedReadSequence>> splitByLibrary(final List<PairedReadSequence> input,
                                                         final List<SAMReadGroupRecord> rgs) {
        final Map<String, List<PairedReadSequence>> out = new HashMap<String, List<PairedReadSequence>>();
        for (final PairedReadSequence seq : input) {
            String library;

            if (seq.getReadGroup() != -1) {
                library = rgs.get(seq.getReadGroup()).getLibrary();
                if (library == null) library = UNKNOWN_LIBRARY;
            } else {
                library = UNKNOWN_LIBRARY;
            }

            List<PairedReadSequence> librarySeqs = out.get(library);
            if (librarySeqs == null) {
                librarySeqs = new ArrayList<PairedReadSequence>();
                out.put(library, librarySeqs);
            }
            librarySeqs.add(seq);
            seq.indexInLibrary = librarySeqs.size() - 1;
        }
        return out;
    }

    private class HashedPairedReadEnd {
        byte[] read;
        int[] shingleHashValue;

        private HashedPairedReadEnd(byte[] read, int[] shingleHashValue){
            this.read = read;
            this.shingleHashValue = shingleHashValue;
        }
    }
}

