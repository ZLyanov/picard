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

import htsjdk.samtools.Defaults;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.samtools.util.StringUtil;
import picard.PicardException;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.Metrics;
import picard.sam.DuplicationMetrics;
import picard.sam.markduplicates.util.AbstractOpticalDuplicateFinderCommandLineProgram;
import picard.sam.markduplicates.util.OpticalDuplicateFinder;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.TimeUnit;

import static java.lang.Math.pow;

/**
 * <p>Attempts to estimate library complexity from sequence alone. Does so by sorting all reads
 * by the first N bases (5 by default) of each read and then comparing reads with the first
 * N bases identical to each other for duplicates.  Reads are considered to be duplicates if
 * they match each other with no gaps and an overall mismatch rate less than or equal to
 * MAX_DIFF_RATE (0.03 by default).</p>
 * <p/>
 * <p>Reads of poor quality are filtered out so as to provide a more accurate estimate. The filtering
 * removes reads with any no-calls in the first N bases or with a mean base quality lower than
 * MIN_MEAN_QUALITY across either the first or second read.</p>
 * <p/>
 * <p>The algorithm attempts to detect optical duplicates separately from PCR duplicates and excludes
 * these in the calculation of library size. Also, since there is no alignment to screen out technical
 * reads one further filter is applied on the data.  After examining all reads a Histogram is built of
 * [#reads in duplicate set -> #of duplicate sets]; all bins that contain exactly one duplicate set are
 * then removed from the Histogram as outliers before library size is estimated.</p>
 *
 * @author Tim Fennell
 */
@CommandLineProgramProperties(
        usage = "Attempts to estimate library complexity from sequence of read pairs alone. Does so by sorting all reads " +
                "by the first N bases (5 by default) of each read and then comparing reads with the first " +
                "N bases identical to each other for duplicates.  Reads are considered to be duplicates if " +
                "they match each other with no gaps and an overall mismatch rate less than or equal to " +
                "MAX_DIFF_RATE (0.03 by default).\n\n" +
                "Reads of poor quality are filtered out so as to provide a more accurate estimate. The filtering " +
                "removes reads with any no-calls in the first N bases or with a mean base quality lower than " +
                "MIN_MEAN_QUALITY across either the first or second read.\n\n" +
                "Unpaired reads are ignored in this computation.\n\n" +
                "The algorithm attempts to detect optical duplicates separately from PCR duplicates and excludes " +
                "these in the calculation of library size. Also, since there is no alignment to screen out technical " +
                "reads one further filter is applied on the data.  After examining all reads a Histogram is built of " +
                "[#reads in duplicate set -> #of duplicate sets] all bins that contain exactly one duplicate set are " +
                "then removed from the Histogram as outliers before library size is estimated.",
        usageShort = "Estimates library complexity from the sequence of read pairs",
        programGroup = Metrics.class
)
public class EstimateLibraryComplexity extends AbstractOpticalDuplicateFinderCommandLineProgram {

    @Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "One or more files to combine and " +
            "estimate library complexity from. Reads can be mapped or unmapped.")
    public List<File> INPUT;

    @Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output file to writes per-library metrics to.")
    public File OUTPUT;

    @Option(doc = "The minimum number of bases at the starts of reads that must be identical for reads to " +
            "be grouped together for duplicate detection.  In effect total_reads / 4^max_id_bases reads will " +
            "be compared at a time, so lower numbers will produce more accurate results but consume " +
            "exponentially more memory and CPU.")
    public int MIN_IDENTICAL_BASES = 5;

    @Option(doc = "The maximum rate of differences between two reads to call them identical.")
    public double MAX_DIFF_RATE = 0.03;

    @Option(doc = "The minimum mean quality of the bases in a read pair for the read to be analyzed. Reads with " +
            "lower average quality are filtered out and not considered in any calculations.")
    public int MIN_MEAN_QUALITY = 20;

    @Option(doc = "Do not process self-similar groups that are this many times over the mean expected group size. " +
            "I.e. if the input contains 10m read pairs and MIN_IDENTICAL_BASES is set to 5, then the mean expected " +
            "group size would be approximately 10 reads.")
    public int MAX_GROUP_RATIO = 500;

    @Option(doc = "Barcode SAM tag (ex. BC for 10X Genomics)", optional = true)
    public String BARCODE_TAG = null;

    @Option(doc = "Read one barcode SAM tag (ex. BX for 10X Genomics)", optional = true)
    public String READ_ONE_BARCODE_TAG = null;

    @Option(doc = "Read two barcode SAM tag (ex. BX for 10X Genomics)", optional = true)
    public String READ_TWO_BARCODE_TAG = null;

    @Option(doc = "Number of bases which contains in one read")
    public int READ_SIZE = 150;

    /**
     * Batch size SAMRecords
     */
    public static final int SAMRECORDS_PACK_SIZE = 1024;

    public int groupsQueueCapacity;

    private final Log log = Log.getInstance(EstimateLibraryComplexity.class);

    /**
     * Little class to hold the sequence of a pair of reads and tile location information.
     */
    static class PairedReadSequence implements OpticalDuplicateFinder.PhysicalLocation {
        int numberBasesInRead;
        short readGroup = -1;
        short tile = -1;
        short x = -1;
        short y = -1;
        boolean qualityOk = true;
        byte[] read1;
        byte[] read2;
        short libraryId;
        int[] hashValuesRead1;
        int[] hashValuesRead2;
        HashSet<PairedReadSequence> possibleCopies;
        int indexInLibrary;

        public PairedReadSequence(int numberBasesInRead){
            this.numberBasesInRead = numberBasesInRead;
        }

        public int getSizeInBytes() {
            // rough guess at memory footprint
            return 16 + 4 + (2 * 4) + 1 + 2 * (24 + 8 + numberBasesInRead) + 2 + (2 * (24 + 8)) + 8 + 4;
        }

        public short getReadGroup() { return this.readGroup; }

        public void setReadGroup(final short readGroup) { this.readGroup = readGroup; }

        public short getTile() { return this.tile; }

        public void setTile(final short tile) { this.tile = tile; }

        public short getX() { return this.x; }

        public void setX(final short x) { this.x = x; }

        public short getY() { return this.y; }

        public void setY(final short y) { this.y = y; }

        public short getLibraryId() { return this.libraryId; }

        public void setLibraryId(final short libraryId) { this.libraryId = libraryId; }

        public SortingCollection.Codec<PairedReadSequence> getCodec() {
            return new PairedReadCodec(numberBasesInRead);
        }
    }

    static class PairedReadSequenceWithBarcodes extends PairedReadSequence {
        int barcode; // primary barcode for this read (and pair)
        int readOneBarcode; // read one barcode, 0 if not present
        int readTwoBarcode; // read two barcode, 0 if not present or not paired

        public PairedReadSequenceWithBarcodes(int numberBasesInRead) {
            super(numberBasesInRead);
        }

        public PairedReadSequenceWithBarcodes(final PairedReadSequence val) {
            super(val.numberBasesInRead);
            if (null == val) throw new PicardException("val was null");
            this.readGroup = val.getReadGroup();
            this.tile = val.getTile();
            this.x = val.getX();
            this.y = val.getY();
            this.qualityOk = val.qualityOk;
            this.read1 = val.read1.clone();
            this.read2 = val.read2.clone();
            this.libraryId = val.getLibraryId();
        }

        public int getSizeInBytes() {
            return super.getSizeInBytes() + (3 * 4); // rough guess at memory footprint
        }
    }

    /**
     * Codec class for writing and read PairedReadSequence objects.
     */
    static class PairedReadCodec implements SortingCollection.Codec<PairedReadSequence> {
        protected int numberBasesInRead;
        protected DataOutputStream out;
        protected DataInputStream in;

        public PairedReadCodec(int numberBasesInRead){
            this.numberBasesInRead = numberBasesInRead;
        }

        public void setOutputStream(final OutputStream out) {
            this.out = new DataOutputStream(out);
        }

        public void setInputStream(final InputStream in) {
            this.in = new DataInputStream(in);
        }

        public void encode(final PairedReadSequence val) {
            try {
                this.out.writeShort(val.readGroup);
                this.out.writeShort(val.tile);
                this.out.writeShort(val.x);
                this.out.writeShort(val.y);
                this.out.writeInt(val.read1.length);
                this.out.write(val.read1);
                this.out.writeInt(val.read2.length);
                this.out.write(val.read2);
            } catch (final IOException ioe) {
                throw new PicardException("Error write out read pair.", ioe);
            }
        }

        public PairedReadSequence decode() {
            try {
                final PairedReadSequence val = new PairedReadSequence(numberBasesInRead);
                try {
                    val.readGroup = this.in.readShort();
                } catch (final EOFException eof) {
                    return null;
                }

                val.tile = this.in.readShort();
                val.x = this.in.readShort();
                val.y = this.in.readShort();

                int length = this.in.readInt();
                val.read1 = new byte[length];
                if (this.in.read(val.read1) != length) {
                    throw new PicardException("Could not read " + length + " bytes from temporary file.");
                }

                length = this.in.readInt();
                val.read2 = new byte[length];
                if (this.in.read(val.read2) != length) {
                    throw new PicardException("Could not read " + length + " bytes from temporary file.");
                }

                return val;
            } catch (final IOException ioe) {
                throw new PicardException("Exception reading read pair.", ioe);
            }
        }

        @Override
        public SortingCollection.Codec<PairedReadSequence> clone() { return new PairedReadCodec(numberBasesInRead); }
    }


    /**
     * Codec class for writing and read PairedReadSequence objects.
     */
    static class PairedReadWithBarcodesCodec extends PairedReadCodec {

        public PairedReadWithBarcodesCodec(int numberBasesInRead) {
            super(numberBasesInRead);
        }

        @Override
        public void encode(final PairedReadSequence val) {
            if (!(val instanceof PairedReadSequenceWithBarcodes)) {
                throw new PicardException("Val was not a PairedReadSequenceWithBarcodes");
            }
            final PairedReadSequenceWithBarcodes data = (PairedReadSequenceWithBarcodes) val;

            super.encode(val);

            try {
                this.out.writeInt(data.barcode);
                this.out.writeInt(data.readOneBarcode);
                this.out.writeInt(data.readTwoBarcode);
            } catch (final IOException ioe) {
                throw new PicardException("Error write out read pair.", ioe);
            }
        }

        @Override
        public PairedReadSequence decode() {
            try {
                final PairedReadSequence parentVal = super.decode();
                if (null == parentVal) return null; // EOF
                final PairedReadSequenceWithBarcodes val = new PairedReadSequenceWithBarcodes(parentVal);
                val.barcode = this.in.readInt();
                val.readOneBarcode = this.in.readInt();
                val.readTwoBarcode = this.in.readInt();

                return val;
            } catch (final IOException ioe) {
                throw new PicardException("Exception reading read pair.", ioe);
            }
        }

        @Override
        public SortingCollection.Codec<PairedReadSequence> clone() { return new PairedReadWithBarcodesCodec(numberBasesInRead); }
    }

    /**
     * Comparator that orders read pairs on the first N bases of both reads.
     */
    class PairedReadComparator implements Comparator<PairedReadSequence> {
        final int BASES = EstimateLibraryComplexity.this.MIN_IDENTICAL_BASES;

        public int compare(final PairedReadSequence lhs, final PairedReadSequence rhs) {
            // First compare the first N bases of the first read
            for (int i = 0; i < BASES; ++i) {
                final int retval = lhs.read1[i] - rhs.read1[i];
                if (retval != 0) return retval;
            }

            // Then compare the first N bases of the second read
            for (int i = 0; i < BASES; ++i) {
                final int retval = lhs.read2[i] - rhs.read2[i];
                if (retval != 0) return retval;
            }

            return 0;
        }
    }

    public int getBarcodeValue(final SAMRecord record) {
        return getReadBarcodeValue(record, BARCODE_TAG);
    }

    public static int getReadBarcodeValue(final SAMRecord record, final String tag) {
        if (null == tag) return 0;
        final String attr = record.getStringAttribute(tag);
        if (null == attr) return 0;
        else return attr.hashCode();
    }

    private int getReadOneBarcodeValue(final SAMRecord record) {
        return getReadBarcodeValue(record, READ_ONE_BARCODE_TAG);
    }

    private int getReadTwoBarcodeValue(final SAMRecord record) {
        return getReadBarcodeValue(record, READ_TWO_BARCODE_TAG);
    }

    /** Stock main method. */
    public static void main(final String[] args) {
        new EstimateLibraryComplexity().instanceMainWithExit(args);
    }

    /**
     * Method that does most of the work.  Reads through the input BAM file and extracts the
     * read sequences of each read pair and sorts them via a SortingCollection.  Then traverses
     * the sorted reads and looks at small groups at a time to find duplicates.
     */
    @Override
    protected int doWork() {
        initMetricProperties();

        for (final File f : INPUT) IOUtil.assertFileIsReadable(f);

        log.info("Will store " + MAX_RECORDS_IN_RAM + " read pairs in memory before sorting.");

        final List<SAMReadGroupRecord> readGroups = new ArrayList<SAMReadGroupRecord>();
        final int recordsRead = 0;
        final SortingCollection<PairedReadSequence> sorter;
        final boolean useBarcodes = (null != BARCODE_TAG || null != READ_ONE_BARCODE_TAG || null != READ_TWO_BARCODE_TAG);

        if (!useBarcodes) {
            sorter = SortingCollection.newInstance(PairedReadSequence.class,
                    new PairedReadCodec(READ_SIZE),
                    new PairedReadComparator(),
                    MAX_RECORDS_IN_RAM,
                    TMP_DIR);
        } else {
            sorter = SortingCollection.newInstance(PairedReadSequence.class,
                    new PairedReadWithBarcodesCodec(READ_SIZE),
                    new PairedReadComparator(),
                    MAX_RECORDS_IN_RAM,
                    TMP_DIR);
        }

        // Loop through the input files and pick out the read sequences etc.
        final ProgressLogger progress = new ProgressLogger(log, (int) 1e6, "Read");
        for (final File f : INPUT) {
            final Map<String, PairedReadSequence> pendingByName = new HashMap<String, PairedReadSequence>();
            final SamReader in = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).
                    async(Defaults.USE_ASYNC_IO).open(f);

            readGroups.addAll(in.getFileHeader().getReadGroups());

            final ExecutorService service = Executors.newSingleThreadExecutor();

            long possibleBatchCapacity = (Runtime.getRuntime().maxMemory()
                    / (100 * (new PairedReadSequence(READ_SIZE)).getSizeInBytes() * SAMRECORDS_PACK_SIZE));
            possibleBatchCapacity = possibleBatchCapacity > 64 ? 64: possibleBatchCapacity;

            int batchCapacity = (Defaults.MULTITHREAD_SORTING_COLLECTION == false)
                    || possibleBatchCapacity < 1 ? 1: (int) possibleBatchCapacity;

            PairReadSequenceAssembler pairReadSequenceAssembler = new PairReadSequenceAssembler(batchCapacity, pendingByName,
                    readGroups, sorter, progress, useBarcodes, opticalDuplicateFinder, MIN_IDENTICAL_BASES, MIN_MEAN_QUALITY,
                    BARCODE_TAG, READ_ONE_BARCODE_TAG, READ_TWO_BARCODE_TAG, READ_SIZE);
            service.execute(pairReadSequenceAssembler);

            try {
                collectReads(in, pairReadSequenceAssembler);

                service.shutdown();
                service.awaitTermination(1, TimeUnit.DAYS);
            } catch (InterruptedException e) {
                service.shutdownNow();
                // Restore the interrupted status.
                Thread.currentThread().interrupt();

                throw new PicardException("Failed to complete EstimateLibraryComplexity", e);
            }
        }


        log.info("Finished reading - moving on to scanning for duplicates.");

        // Now go through the sorted reads and attempt to find duplicates
        final PeekableIterator<PairedReadSequence> iterator = new PeekableIterator<PairedReadSequence>(sorter.iterator());
        int groupsProcessed = 0;
        final int meanGroupSize = Math.max(1, (recordsRead / 2) / (int) pow(4, MIN_IDENTICAL_BASES * 2));
        final int maxGroupSize = meanGroupSize * MAX_GROUP_RATIO;
        final BlockingQueue<List<PairedReadSequence>> groupsQueue =
                new LinkedBlockingQueue<List<PairedReadSequence>>(groupsQueueCapacity);

        final ExecutorService searchService = Executors.newSingleThreadExecutor();

        DuplicateFinder duplicateFinder = new DuplicateFinder(readGroups, groupsQueue,
                opticalDuplicateFinder, useBarcodes, groupsProcessed, MIN_IDENTICAL_BASES, MAX_DIFF_RATE);
        duplicateFinder.initAlgorithmProperties(iterator);

        searchService.execute(duplicateFinder);

        try {
            searchGroups(iterator, maxGroupSize, meanGroupSize, groupsQueue);
            completeSearch(duplicateFinder, sorter, iterator, searchService);
        } catch (InterruptedException e) {
            searchService.shutdownNow();
            // Restore the interrupted status.
            Thread.currentThread().interrupt();

            throw new PicardException("Failed to complete EstimateLibraryComplexity", e);
        }

        finalizeAndWriteMetrics(duplicateFinder);

        return 0;
    }

    private void collectReads(SamReader in, PairReadSequenceAssembler pairReadSequenceAssembler) throws InterruptedException {
        for (final SAMRecord rec : in) {
            if (!rec.getReadPairedFlag()
                    || !rec.getFirstOfPairFlag() && !rec.getSecondOfPairFlag()) {
                continue;
            }

            pairReadSequenceAssembler.add(rec);
        }

        pairReadSequenceAssembler.finish();
        CloserUtil.close(in);
    }

    /**
     * Pulls out of the iterator the next group of reads that can be compared to each other to
     * identify duplicates.
     */
    List<PairedReadSequence> getNextGroup(final PeekableIterator<PairedReadSequence> iterator,
                                          int maxGroupSize, int meanGroupSize) {
        final List<PairedReadSequence> group = new ArrayList<PairedReadSequence>();
        final PairedReadSequence first = iterator.next();
        group.add(first);

        int restSize = 0;

        while (iterator.hasNext()) {
            final PairedReadSequence next = iterator.peek();
            if(!isIdentical(first, next)) {
                break;
            }

            if (group.size() > maxGroupSize) {
                restSize++;
                iterator.next();
                continue;
            }
            group.add(iterator.next());
        }

        if (group.size() > meanGroupSize * this.MAX_GROUP_RATIO) {
            final PairedReadSequence prs = group.get(0);
            log.warn("Omitting group with over " + this.MAX_GROUP_RATIO + " times the expected mean number of read pairs. " +
                    "Mean=" + meanGroupSize + ", Actual=" + (group.size() + restSize) + ". Prefixes: " +
                    StringUtil.bytesToString(prs.read1, 0, this.MIN_IDENTICAL_BASES) +
                    " / " +
                    StringUtil.bytesToString(prs.read2, 0, this.MIN_IDENTICAL_BASES));
            return null;
        }

        return group;
    }

    private void searchGroups(PeekableIterator<PairedReadSequence> iterator, int maxGroupSize, int meanGroupSize,
                              BlockingQueue<List<PairedReadSequence>> groupsQueue) throws InterruptedException {
        while (iterator.hasNext()) {
            // Get the next group and split it apart by library
            final List<PairedReadSequence> group = getNextGroup(iterator, maxGroupSize, meanGroupSize);

            if (group != null) {
                groupsQueue.put(group);
            }
        }
    }

    private boolean isIdentical(PairedReadSequence first, PairedReadSequence second) {
        for (int i = 0; i < MIN_IDENTICAL_BASES; ++i) {
            if (first.read1[i] != second.read1[i] || first.read2[i] != second.read2[i])
                return false;
        }
        return true;
    }

    private void completeSearch(DuplicateFinder duplicateFinder, SortingCollection<PairedReadSequence> sorter,
                                PeekableIterator<PairedReadSequence> iterator,
                                ExecutorService service) throws InterruptedException {
        duplicateFinder.finish();

        service.shutdown();
        service.awaitTermination(1, TimeUnit.DAYS);

        iterator.close();
        sorter.cleanup();
    }

    private void finalizeAndWriteMetrics(DuplicateFinder duplicateFinder) {
        Map<String, Histogram<Integer>> duplicationHistosByLibrary = duplicateFinder.getDuplicationHistosByLibrary();
        Map<String, Histogram<Integer>> opticalHistosByLibrary = duplicateFinder.getOpticalHistosByLibrary();

        final MetricsFile<DuplicationMetrics, Integer> file = getMetricsFile();
        for (final Map.Entry<String, Histogram<Integer>> entry : duplicationHistosByLibrary.entrySet()) {
            final DuplicationMetrics metrics = new DuplicationMetrics();
            metrics.LIBRARY = entry.getKey();
            final Histogram<Integer> duplicationHisto = entry.getValue();
            final Histogram<Integer> opticalHisto = opticalHistosByLibrary.get(metrics.LIBRARY);

            // Filter out any bins that have only a single entry in them and calcu
            for (final Map.Entry<Integer, Histogram.Bin> entryBin : duplicationHisto.entrySet()) {
                final double duplicateGroups = entryBin.getValue().getValue();
                final Integer bin = entryBin.getKey();
                final Histogram.Bin od = opticalHisto.get(bin);
                final double opticalDuplicates = od == null ? 0 : od.getValue();

                if (duplicateGroups > 1) {
                    double value = bin * duplicateGroups;
                    metrics.READ_PAIRS_EXAMINED += value;
                    metrics.READ_PAIR_DUPLICATES += value - duplicateGroups;
                    metrics.READ_PAIR_OPTICAL_DUPLICATES += opticalDuplicates;
                }
            }

            metrics.calculateDerivedMetrics();
            file.addMetric(metrics);
            file.addHistogram(duplicationHisto);
        }

        file.write(OUTPUT);
    }

    /**
     * Initializes the settings of a metric such as the maximum number of reads in the memory.
     */
    private void initMetricProperties() {
        final int sizeInBytes;
        if (null != BARCODE_TAG || null != READ_ONE_BARCODE_TAG || null != READ_TWO_BARCODE_TAG) {
            sizeInBytes = new PairedReadSequenceWithBarcodes(READ_SIZE).getSizeInBytes();
        } else {
            sizeInBytes = new PairedReadSequence(READ_SIZE).getSizeInBytes();
        }
        long possibleGroupsCapacity = (Runtime.getRuntime().maxMemory() / (10 * sizeInBytes * MAX_GROUP_RATIO));
        possibleGroupsCapacity = possibleGroupsCapacity > 64 ? 64: possibleGroupsCapacity;

        groupsQueueCapacity = (Defaults.MULTITHREAD_SORTING_COLLECTION == false)
                || 1 > possibleGroupsCapacity ? 1: (int) possibleGroupsCapacity;

        MAX_RECORDS_IN_RAM = (int) (Runtime.getRuntime().maxMemory() / sizeInBytes) / 2;
    }
}
