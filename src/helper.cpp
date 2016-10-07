#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::export]]
DataFrame calculate_probabilities_helper(DataFrame sortedInteractions,
    NumericVector isOutlier, bool useBinning, int desiredPerBin,
    double distScaling, int observedIntraInRangeSum) {

    NumericVector c1 = sortedInteractions["interactionDistance"];
    IntegerVector c2 = sortedInteractions["hitCount"];

    int lcount = sortedInteractions.nrows();

    // the following five lists will be the print outputs
    NumericVector x(lcount); // avg genomic distances of bins
    NumericVector y(lcount); // avg interaction probabilities of bins
    NumericVector yerr(lcount); // stderrs of bins
    IntegerVector pairCounts(lcount); // number of pairs in bins
    // number of interactions (reads) in bins
    IntegerVector interactionTotals(lcount);

    // the following variables will be used to calculate the above five lists
    int noOfPairsForBin = 0;
    double meanCountPerPair = 0;
    double M2 = 0;
    int interactionTotalForBin = 0;
    int interactionTotalForBinTermination = 0;
    double distanceTotalForBin = 0;
    double lastDistanceForBin = -1;

    int index = 0;

    for (int i = 0; i < lcount; i++) {
        double interactionDistance = c1[i];
        int interactionCount = c2[i];

        // if one bin is full or it's the last bin
        if (noOfPairsForBin > 0 && ((!useBinning &&
            lastDistanceForBin != interactionDistance) || (useBinning
            && lastDistanceForBin != interactionDistance &&
            interactionTotalForBinTermination >= desiredPerBin))) {

            // calculate the things that need to be calculated
            double avgDistance = distanceTotalForBin / noOfPairsForBin *
                distScaling;
            double meanProbabilityObsv = meanCountPerPair /
                observedIntraInRangeSum;
            double se_p = meanProbabilityObsv;
            // update se_q if there are more than 1 pairs in the bin
            if (noOfPairsForBin > 1) {
                double var = M2 / (noOfPairsForBin - 1);
                double sd = sqrt(var);
                double se = sd / sqrt(noOfPairsForBin);
                se_p = se / observedIntraInRangeSum;
            }

            x[index] = avgDistance;
            y[index] = meanProbabilityObsv;
            yerr[index] = se_p;
            pairCounts[index] = noOfPairsForBin;
            interactionTotals[index] = interactionTotalForBin;

            index += 1;

            // now that we saved what we need
            // set the values back to defaults and go on to the next bin
            noOfPairsForBin = 0;
            meanCountPerPair = 0;
            M2 = 0;
            interactionTotalForBin = 0;
            interactionTotalForBinTermination = 0;
            distanceTotalForBin = 0;
            lastDistanceForBin = -1;
        }

        // Now go back to processing the read values of interactionDistance
        // and interactionCount
        // this check is necessary for the second pass of fit-hic
        // we want to only use the non-outlier interactions in our
        // probability calculation
        if (isOutlier[i] == 0) {
            distanceTotalForBin += interactionDistance / distScaling;
            interactionTotalForBin += interactionCount;
            noOfPairsForBin += 1;
            double delta = interactionCount - meanCountPerPair;
            meanCountPerPair += delta / noOfPairsForBin;
            M2 += delta * (interactionCount - meanCountPerPair);
        }

        interactionTotalForBinTermination += interactionCount;
        lastDistanceForBin = interactionDistance;
    }

    x = x[Range(0, index - 1)];
    y = y[Range(0, index - 1)];
    yerr = yerr[Range(0, index - 1)];
    pairCounts = pairCounts[Range(0, index - 1)];
    interactionTotals = interactionTotals[Range(0, index - 1)];

    return DataFrame::create(_["x"] = x,
        _["y"] = y,
        _["yerr"] = yerr,
        _["pairCounts"] = pairCounts,
        _["interactionTotals"] = interactionTotals);
}

// [[Rcpp::export]]
NumericVector benjamini_hochberg_correction_helper(NumericVector p_values,
    int num_total_tests, NumericVector sorted_pvals, IntegerVector order) {

    int size = p_values.size();
    NumericVector q_values(size);

    int i = 0;
    double prev_bh_value = 0;

    while(i < size) {
        double p_value = sorted_pvals[i];
        double bh_value;

        if (p_value == 1) {
            bh_value = 1;
        } else {
            bh_value = p_value * num_total_tests / (i + 1);
            bh_value = fmin(bh_value, 1);
        }

        bh_value = fmax(bh_value, prev_bh_value);
        prev_bh_value = bh_value;
        q_values[order[i]] = bh_value;

        i += 1;
    }

    return q_values;
}
