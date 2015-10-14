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

import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import picard.PicardException;
import picard.util.MathUtil;

import java.io.IOException;
import java.util.Arrays;
import java.util.*;

/**
 * Created by David Benjamin on 5/13/15.
 */
public class TheoreticalSensitivity {
    /*
    @param the probability of depth n is depthDistribution[n] for n = 0, 1. . . N - 1
    @param the probability of quality q is qualityDistribution[q] for q = 0, 1. . . Q
    @param sample size is the number of random sums of quality scores for each m
    @param logOddsThreshold is the log_10 of the likelihood ratio required to call a SNP,
    for example 5 if the variant likelihood must be 10^5 times greater
     */
    private static final Log log = Log.getInstance(TheoreticalSensitivity.class);

    public static double hetSNPSensitivity(final long[] depthDistribution, final long[] qualityDistribution,
                                  final int sampleSize, final double logOddsThreshold) {
        final int N = depthDistribution.length;

        log.info("Creating Roulette Wheel");
        final RouletteWheel qualitySampler = new RouletteWheel(qualityDistribution);

        //qualitySums[m] is a random sample of sums of m quality scores, for m = 0, 1, N - 1
        log.info("Calculating quality sums from quality sampler");
        final List<ArrayList<Integer>> qualitySums = qualitySampler.sampleCumulativeSums(N, sampleSize);

        //if a quality sum of m qualities exceeds the quality sum threshold for n total reads, a SNP is called
        final ArrayList<Double> qualitySumThresholds = new ArrayList<Double>(N);
        for (int n = 0; n < N; n++) qualitySumThresholds.add(10*(n*Math.log10(2) + logOddsThreshold));


        //probabilityToExceedThreshold[m][n] is the probability that the sum of m quality score
        //exceeds the nth quality sum threshold
        log.info("Calculating theoretical het sensitivity");
        final List<ArrayList<Double>> probabilityToExceedThreshold = proportionsAboveThresholds(qualitySums, qualitySumThresholds);
        final List<ArrayList<Double>> altDepthDistribution = hetAltDepthDistribution(N);
        double result = 0.0;
        for (int n = 0; n < N; n++) {
            for (int m = 0; m <= n; m++) {
                result += depthDistribution[n] * altDepthDistribution.get(n).get(m)* probabilityToExceedThreshold.get(m).get(n);
            }
        }
        return result;
    }

    //given L lists of lists and N thresholds, count the proportion of each list above each threshold
    public static List<ArrayList<Double>> proportionsAboveThresholds(final List<ArrayList<Integer>> lists, final List<Double> thresholds) {
        final ArrayList<ArrayList<Double>> result = new ArrayList<ArrayList<Double>>();

        for (final ArrayList<Integer> list : lists) {
            final ArrayList<Double> newRow = new ArrayList<Double>(Collections.nCopies(thresholds.size(),0.0));
            Collections.sort(list);
            int n = 0;
            int j = 0;  //index within the ordered sample
            while (n < thresholds.size() && j < list.size()) {
                if (thresholds.get(n) > list.get(j)) j++;
                else newRow.set(n++, (double) (list.size() - j)/list.size());
            }
            result.add(newRow);
        }
        return result;

    }


    //Utility function for making table of binomial distribution probabilities nCm * (0.5)^n
    //for n = 0, 1 . . . N - 1 and m = 0, 1. . . n
    public static List<ArrayList<Double>> hetAltDepthDistribution(final int N) {
        final List<ArrayList<Double>> table = new ArrayList<ArrayList<Double>>();
        for (int n = 0; n < N; n++) {
            final ArrayList<Double> nthRow = new ArrayList<Double>();

            //add the 0th element, then elements 1 through n - 1, then the nth.
            //Note that nCm = (n-1)C(m-1) * (n/m)
            nthRow.add(Math.pow(0.5, n));
            for (int m = 1; m < n; m++) nthRow.add((n*0.5/m)*table.get(n - 1).get(m - 1));
            if (n > 0) nthRow.add(nthRow.get(0));

            table.add(nthRow);
        }

        return table;

    }


    /*
    Perform random draws from {0, 1. . . N - 1} according to a list of relative probabilities.

    We use an O(1) stochastic acceptance algorithm -- see Physica A, Volume 391, Page 2193 (2012) --
    which works well when the ratio of maximum weight to average weight is not large.
     */
    public static class RouletteWheel {
        final private List<Double> probabilities;
        final private int N;
        private int count=0;

        RouletteWheel(final long[] weights) {
            N = weights.length;

            probabilities = new ArrayList<Double>();
            final double wMax = (double)MathUtil.max(weights);
            if(wMax==0){throw new PicardException("Quality score distribution is empty.");}
            for (final long w : weights) {
                probabilities.add(w / wMax);
            }
        }

        public int draw() {
            while (true) {
                final int n = (int) (N * Math.random());
                count++;
                if (Math.random() < probabilities.get(n)) {
                    count = 0;
                    return n;
                }
                else if(count>=600){
                    return 0;
                }
            }
        }

        //get samples of sums of 0, 1, 2,. . .  N - 1 draws
        public List<ArrayList<Integer>> sampleCumulativeSums(final int maxNumberOfSummands, final int sampleSize) {
            final List<ArrayList<Integer>> result = new ArrayList<ArrayList<Integer>>();
            for (int m = 0; m < maxNumberOfSummands; m++) result.add(new ArrayList<Integer>());

            for (int iteration = 0; iteration < sampleSize; iteration++) {
                int cumulativeSum = 0;
                for (int m = 0; m < maxNumberOfSummands; m++) {
                    result.get(m).add(cumulativeSum);
                    cumulativeSum += draw();
                }
                if(iteration%1000==0){log.info(iteration + " sampling iterations completed");}
            }
            return result;
        }

    }

}
