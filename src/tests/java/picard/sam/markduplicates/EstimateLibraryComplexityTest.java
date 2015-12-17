/*
 * The MIT License
 *
 * Copyright (c) 2014 The Broad Institute
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

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.Histogram;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;
import picard.sam.DuplicationMetrics;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

public class EstimateLibraryComplexityTest extends CommandLineProgramTest {

    private static final File TEST_DATA_DIR = new File("testdata/picard/sam");

    @Override
    public String getCommandLineProgramName() {
        return EstimateLibraryComplexity.class.getSimpleName();
    }

    @Test
    public void testMetricsELC() throws IOException {
        final File input = new File(TEST_DATA_DIR, "forMetrics.sam");
        final File outfile   = File.createTempFile("metrics", ".txt");
        outfile.deleteOnExit();
        final String[] args = new String[] {
                "INPUT="  + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath()
        };

        HashMap<Integer, Double> expectedHisto = new HashMap<Integer, Double>();
        createExpectedHistogram(expectedHisto);

        Assert.assertEquals(runPicardCommandLine(args), 0);

        // Check the values written to metrics.txt against our input expectations
        final MetricsFile<DuplicationMetrics, Comparable<?>> metricsOutput = new MetricsFile<DuplicationMetrics, Comparable<?>>();
        metricsOutput.read(new FileReader(outfile));

        DuplicationMetrics metrics = metricsOutput.getMetrics().get(metricsOutput.getMetrics().size() - 1);
        Histogram<Comparable<?>> duplicationHisto = metricsOutput.getHistogram();

        Assert.assertEquals(metrics.UNPAIRED_READS_EXAMINED, 0);
        Assert.assertEquals(metrics.READ_PAIRS_EXAMINED, 4);
        Assert.assertEquals(metrics.UNMAPPED_READS, 0);
        Assert.assertEquals(metrics.UNPAIRED_READ_DUPLICATES, 0);
        Assert.assertEquals(metrics.READ_PAIR_DUPLICATES, 0);
        Assert.assertEquals(metrics.READ_PAIR_OPTICAL_DUPLICATES, 0);
        Assert.assertEquals(metrics.PERCENT_DUPLICATION, 0.0d);

        for(Map.Entry<Integer, Double> entry: expectedHisto.entrySet()){
            Assert.assertEquals(duplicationHisto.get(entry.getKey()).getValue(), entry.getValue());
        }
    }

    private void createExpectedHistogram(HashMap<Integer, Double> refHisto) {
        refHisto.put(1, 4.0d);
        refHisto.put(2, 1.0d);
    }
}
