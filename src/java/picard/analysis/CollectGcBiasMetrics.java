/*
 * The MIT License
 *
 * Copyright (c) 2015 The Broad Institute
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

package picard.analysis;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.IOUtil;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.programgroups.Metrics;
import picard.metrics.GcBiasMetrics;
import picard.util.RExecutor;

import java.io.File;
import java.text.NumberFormat;
import java.util.List;
import java.util.Set;

/**
 * Tool to collect information about GC bias in the reads in a given BAM file. Computes
 * the number of windows (of size specified by SCAN_WINDOW_SIZE) in the genome at each GC%
 * and counts the number of read starts in each GC bin.  What is output and plotted is
 * the "normalized coverage" in each bin - i.e. the number of reads per window normalized
 * to the average number of reads per window across the whole genome.
 *
 * @author Tim Fennell
 * edited by Kylee Bergin
 */
@CommandLineProgramProperties(
        usage = CollectGcBiasMetrics.USAGE_SUMMARY + CollectGcBiasMetrics.USAGE_DETAILS,
        usageShort = CollectGcBiasMetrics.USAGE_SUMMARY,
        programGroup = Metrics.class
)
public class CollectGcBiasMetrics extends SinglePassSamProgram {

    static final String USAGE_SUMMARY = "Collects information regarding GC bias from a SAM/BAM input file.  ";
    static final String USAGE_DETAILS = "Tool that collects information about the proportions of guanine (G) and cytosine (C)" +
            " nucleotides in a sample.  Regions of high and low G + C content have been shown to interfere with mapping/aligning," +
            " ultimately leading to low read depth and fragmented genome assemblies, a phenomenon known as \"GC bias\".  " +
            "Detailed information on the effects GC-bias on NGS data can be found at DOI: 10.1371/journal.pone.0062856/.<br /><br />." +
            "" +
            "For each run, the corresponding reference sequence is divided into bins or windows based on the percentage of G + C" +
            " content ranging from 0 - 100%.  The percentages of G + C are determined from a defined length of sequence, the default " +
            "value is set at 100 bases.   " +
            "Although the mean of the distribution will vary among organisms, human DNA has a mean GC-content of 40%, suggesting " +
            "preponderance of AT-rich regions.  <br /><br />" +

            "GC bias is calculated and output as both summary (optional) and detailed metrics (required).  The output tables for both the summary and the detailed metrics are \".txt\" files while the chart is a \".pdf\" file. <br /><br /> " +

            "The GcBiasSummaryMetrics provides high-level metrics that capture run-specific bias information including" +
            " WINDOW_SIZE, ALIGNED_READS, TOTAL_CLUSTERS, AT_DROPOUT, and GC_DROPOUT.  While WINDOW_SIZE refers to the" +
            "numbers of bases used for the distribution (see above), the ALIGNED_READS and" +
            " TOTAL_CLUSTERS are the total number of aligned reads and the total number of reads (after filtering) " +
            "produced in a run.   In addition, the tool produces both AT_DROPOUT and GC_DROPOUT metrics, which indicate the percentage of misaligned reads that correlate with low (%-GC is < 50%) or high (%-GC is > 50%) GC content respectively.  <br /><br />" +
            "" +
            "GcBiasDetailedMetrics produces both a chart (pdf) and a data table.  The table output includes GC percentages " +
            "for each bin (GC), the percentage of WINDOWS corresponding to each GC bin of the reference sequence, the numbers of reads that start within a particular %GC content bin (READ_STARTS), " +
            "and the mean base quality of the reads that correspond to a specific GC-content distribution window (MEAN_BASE_QUALITY).  NORMALIZED_COVERAGE is a relative measure of sequence coverage by the reads at a particular GC-content." +

            "  The percentage of \"coverage\" or depth in a GC bin is calculated by dividing the number of reads of a particular GC content, " +
            "by the mean number of reads of all GC bins.  A number of 1 represents mean coverage, a number less than " +
            "one represents lower than mean coverage (e.g. 0.5 means half as much coverage as average) while a " +
            "number greater than one represents higher than mean coverage (e.g. 3.1 means this GC bin has 3.1 times" +
            " more reads per window than average).  Tool also plots mean base-quality scores of the reads within each GC-content bin, enabling the user to determine how base-quality scores vary with GC-content.  <br /> <br />"+
            "The chart ouptut associated with this data table plots the NORMALIZED_COVERAGE, the distribution of WINDOWs corresponding to GC percentages, and base qualities corresponding to each %GC bin."+


            "<h4>Usage Example:</h4>"+
            "<pre>" +
            "java -jar picard.jar CollectGcBiasMetrics \\<br />"+
            "      I=input.bam \\<br />"+
            "      O=gc_bias_metrics.txt \\<br />"+
            "      CHART=gc_bias_metrics.pdf \\<br />"+
            "      R=reference_sequence.fasta"+
            "</pre>"+
            "For detailed explanations of the output metrics, please see: " +
            "https://broadinstitute.github.io/picard/picard-metric-definitions.html#GcBiasMetrics" +
            "<hr />";
    /** The location of the R script to do the plotting. */
    private static final String R_SCRIPT = "picard/analysis/gcBias.R";

    // Usage and parameters

    @Option(shortName = "CHART", doc = "The PDF file to render the chart to.")
    public File CHART_OUTPUT;

    @Option(shortName = "S", doc = "The text file to write summary metrics to.")
    public File SUMMARY_OUTPUT;

    @Option(shortName = "WINDOW_SIZE", doc = "The size of the scanning windows on the reference genome that are used to bin reads.")
    public int SCAN_WINDOW_SIZE = 100;

    @Option(shortName = "MGF", doc = "For summary metrics, exclude GC windows that include less than this fraction of the genome.")
    public double MINIMUM_GENOME_FRACTION = 0.00001;

    @Option(shortName = "BS", doc = "Whether the SAM or BAM file consists of bisulfite sequenced reads.")
    public boolean IS_BISULFITE_SEQUENCED = false;

    @Option(shortName = "LEVEL", doc = "The level(s) at which to accumulate metrics.")
    public Set<MetricAccumulationLevel> METRIC_ACCUMULATION_LEVEL = CollectionUtil.makeSet(MetricAccumulationLevel.ALL_READS);

    // Calculates GcBiasMetrics for all METRIC_ACCUMULATION_LEVELs provided
    private GcBiasMetricsCollector multiCollector;

    // Bins for the histograms to track the number of windows at each GC, and the number of read starts
    // at bins of each GC %. Need 101 to get from 0-100.
    private static final int BINS = 101;

    ////////////////////////////////////////////////////////////////////////////
    // Stock main method
    ////////////////////////////////////////////////////////////////////////////
    public static void main(final String[] args) {
        System.exit(new CollectGcBiasMetrics().instanceMain(args));
    }

    /////////////////////////////////////////////////////////////////////////////
    // Setup calculates windowsByGc for the entire reference. Must be done at
    // startup to avoid missing reference contigs in the case of small files
    // that may not have reads aligning to every reference contig.
    /////////////////////////////////////////////////////////////////////////////
    @Override
    protected void setup(final SAMFileHeader header, final File samFile) {
        IOUtil.assertFileIsWritable(CHART_OUTPUT);
        IOUtil.assertFileIsWritable(SUMMARY_OUTPUT);
        IOUtil.assertFileIsReadable(REFERENCE_SEQUENCE);

        //Calculate windowsByGc for the reference sequence
        final int[] windowsByGc = GcBiasUtils.calculateRefWindowsByGc(BINS, REFERENCE_SEQUENCE, SCAN_WINDOW_SIZE);

        //Delegate actual collection to GcBiasMetricCollector
        multiCollector = new GcBiasMetricsCollector(METRIC_ACCUMULATION_LEVEL, windowsByGc, header.getReadGroups(), SCAN_WINDOW_SIZE, IS_BISULFITE_SEQUENCED);
    }

    ////////////////////////////////////////////////////////////////////////////
    // MultiCollector acceptRead
    ////////////////////////////////////////////////////////////////////////////
    @Override
    protected void acceptRead(final SAMRecord rec, final ReferenceSequence ref) {
        multiCollector.acceptRecord(rec, ref);
    }

    /////////////////////////////////////////////////////////////////////////////
    // Write out all levels of normalized coverage metrics to a file
    /////////////////////////////////////////////////////////////////////////////
    @Override
    protected void finish() {
        multiCollector.finish();
        final MetricsFile<GcBiasMetrics, Integer> file = getMetricsFile();
        final MetricsFile<GcBiasDetailMetrics, ?> detailMetricsFile = getMetricsFile();
        final MetricsFile<GcBiasSummaryMetrics, ?> summaryMetricsFile = getMetricsFile();
        multiCollector.addAllLevelsToFile(file);
        final List<GcBiasMetrics> gcBiasMetricsList = file.getMetrics();
        for(final GcBiasMetrics gcbm : gcBiasMetricsList){
            final List<GcBiasDetailMetrics> gcDetailList = gcbm.DETAILS.getMetrics();
            for(final GcBiasDetailMetrics d : gcDetailList) {
                detailMetricsFile.addMetric(d);
            }
            summaryMetricsFile.addMetric(gcbm.SUMMARY);
        }
        detailMetricsFile.write(OUTPUT);
        summaryMetricsFile.write(SUMMARY_OUTPUT);

        final NumberFormat fmt = NumberFormat.getIntegerInstance();
        fmt.setGroupingUsed(true);
        RExecutor.executeFromClasspath(R_SCRIPT,
                OUTPUT.getAbsolutePath(),
                SUMMARY_OUTPUT.getAbsolutePath(),
                CHART_OUTPUT.getAbsolutePath(),
                String.valueOf(SCAN_WINDOW_SIZE));
    }
}


