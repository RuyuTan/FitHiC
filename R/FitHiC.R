#' Fit-Hi-C
#'
#' Fit-Hi-C is a tool for assigning statistical confidence estimates
#' to intra-chromosomal contact maps produced by genome-wide genome
#' architecture assays such as Hi-C.
#'
#' @param fragsfile The path specifies where FRAGSFILE is located in the
#'                  file system. FRAGSFILE stores the information about
#'                  midpoints (or start indices) of the fragments. It should
#'                  consist of 5 columns: first column stands for chromosome
#'                  name; third column stands for the midPoint; fourth column
#'                  stands for the hitCount; second column and fifth column
#'                  can be arbitrary.
#' @param intersfile The path specifies where INTERSFILE is located in the
#'                   file sytem. INTERSFILE stores the information about
#'                   interactions between fragment pairs. It should consist
#'                   of 5 columns: first column and third column stand for
#'                   the chromosome names of the fragment pair; second column
#'                   and fourth column stand for midPoints of the fragment
#'                   pair; fifth column stands for hitCount.
#' @param outdir The path specifies where the output files will be stored
#'               in the file system. If the path does not exist, it will be
#'               automatically created.
#' @param biasfile The path specifies where BIASFILE is located in the file
#'                 system. BIASFILE stores the information about biases
#'                 calculated by ICE for each locus. It should consist of
#'                 3 columns: first column stands for chromosome name; second
#'                 column stands for the midPoint; third column stands for
#'                 the bias. This argument is OPTIONAL.
#' @param noOfPasses Number of passes after the initial (before) fit. DEFAULT
#'                   is 1 (after).
#' @param noOfBins Number of equal-occupancy (count) bins. Default is 100.
#' @param mappabilityThreshold Minimum number of hits per locus that has to
#'                             exist to call it mappable. DEFAULT is 1.
#' @param libname Name of the library that is analyzed to be used for plots.
#'                DEFAULT is empty.
#' @param distUpThres Upper bound on the intra-chromosomal distance range
#'                    (unit: base pairs). DEFAULT is no limit.
#' @param distLowThres Lower bound on the intra-chromosomal distance range
#'                     (unit: base pairs). DEFAULT is no limit.
#' @param visual Use this flag for generating plots. DEFAULT is False.
#'
#' @return None
#'
#' @author Ruyu Tan, \email{tanry33@gmail.com}
#'
#' @examples
#' FitHiC(system.file("extdata", "fragmentLists/Duan_yeast_EcoRI.gz",
#'     package = "FitHiC"), system.file("extdata",
#'     "contactCounts/Duan_yeast_EcoRI.gz", package = "FitHiC"), getwd(),
#'     libname="Duan_yeast_EcoRI", distUpThres=250000, distLowThres=10000,
#'     visual=TRUE)
#'
#' @import data.table
#' @import fdrtool
#' @importFrom grDevices dev.off png
#' @importFrom graphics arrows legend lines par plot
#' @importFrom Rcpp evalCpp
#' @importFrom stats pbinom predict smooth.spline
#' @importFrom utils read.table write.table
#' @useDynLib FitHiC, .registration=TRUE
#' @export
FitHiC <- function(fragsfile, intersfile, outdir, biasfile="none", noOfPasses=1,
noOfBins=100, mappabilityThreshold=1, libname="", distUpThres=-1,
distLowThres=-1, visual=FALSE) {

    distScaling <- 10000.0
    toKb <- 10**-3
    toMb <- 10**-6
    toProb <- 10**5

    useBinning <- TRUE # This is no more an option
    useInters <- FALSE # This is no more an option

    if (distUpThres == -1) {
        distUpThres <- Inf # -1 by default, means no upper bound
    }
    if (distLowThres == -1) {
        distLowThres <- (-Inf) # -1 by default, means no lower bound
    }

    print("Fit-Hi-C is processing ...")

    r1 <- generate_FragPairs(fragsfile, mappabilityThreshold, distUpThres,
        distLowThres)

    possibleInterAllCount <- r1[["possibleInterAllCount"]]
    possibleIntraAllCount <- r1[["possibleIntraAllCount"]]
    possibleIntraInRangeCount <- r1[["possibleIntraInRangeCount"]]

    listOfMappableFrags <- r1[["listOfMappableFrags"]]

    possiblePairsPerDistance <- r1[["possiblePairsPerDistance"]]

    baselineInterChrProb <- r1[["baselineInterChrProb"]]
    baselineIntraChrProb <- r1[["baselineIntraChrProb"]]

    biasDic <- NULL
    if (biasfile != "none") {
        biasDic <- read_ICE_biases(biasfile)
    }

    r2 <- read_All_Interactions(intersfile, biasDic, listOfMappableFrags,
        possiblePairsPerDistance, distUpThres, distLowThres)

    sortedInteractions <- r2[["sortedInteractions"]]

    observedInterAllSum <- r2[["observedInterAllSum"]]
    observedIntraAllSum <- r2[["observedIntraAllSum"]]
    observedIntraInRangeSum <- r2[["observedIntraInRangeSum"]]

    observedInterAllCount <- r2[["observedInterAllCount"]]
    observedIntraAllCount <- r2[["observedIntraAllCount"]]
    observedIntraInRangeCount <- r2[["observedIntraInRangeCount"]]

    minObservedGenomicDist <- r2[["minObservedGenomicDist"]]
    maxObservedGenomicDist <- r2[["maxObservedGenomicDist"]]

    possiblePairsPerDistance <- r2[["possiblePairsPerDistance"]]

    tempData <- calculate_Probabilities(sortedInteractions,
        rep(0, nrow(sortedInteractions)),
        paste(libname, ".fithic_pass1", sep=""), noOfBins, useBinning,
        distScaling, observedIntraInRangeSum, outdir, visual, libname, toKb,
        toProb)

    isOutlier <- fit_Spline(tempData$x, tempData$y, tempData$yerr, intersfile,
        sortedInteractions, biasDic, paste(libname, ".spline_pass1", sep=""),
        1, outdir, visual, distLowThres, distUpThres, toKb, toProb, useInters,
        baselineIntraChrProb, baselineInterChrProb, observedInterAllSum,
        observedIntraAllSum, observedIntraInRangeSum, possibleInterAllCount,
        possibleIntraAllCount, possibleIntraInRangeCount,
        maxObservedGenomicDist)

    for (i in seq(2, noOfPasses + 1)) {
        tempData <- calculate_Probabilities(sortedInteractions, isOutlier,
            paste(libname, ".fithic_pass", i, sep=""), noOfBins, useBinning,
            distScaling, observedIntraInRangeSum, outdir, visual, libname,
            toKb, toProb)

        isOutlier <- fit_Spline(tempData$x, tempData$y, tempData$yerr,
            intersfile, sortedInteractions, biasDic,
            paste(libname, ".spline_pass", i, sep=""), i, outdir, visual,
            distLowThres, distUpThres, toKb, toProb, useInters,
            baselineIntraChrProb, baselineInterChrProb, observedInterAllSum,
            observedIntraAllSum, observedIntraInRangeSum, possibleInterAllCount,
            possibleIntraAllCount, possibleIntraInRangeCount,
            maxObservedGenomicDist)
    }

    print("Execution of Fit-Hi-C completed successfully. [DONE]")

    return
}

generate_FragPairs <- function(infilename, mappabilityThreshold, distUpThres,
distLowThres) {

    print("Running generate_FragPairs method ...")

    chr <- hitCount <- mid1 <- mid2 <- NULL

    possibleInterAllCount <- 0
    possibleIntraAllCount <- 0
    possibleIntraInRangeCount <- 0

    baselineIntraChrProb <- 0 # 1.0/possibleIntraAllCount
    baselineInterChrProb <- 0 # 1.0/possibleInterAllCount

    # read the fragments file
    data <- data.table(read.table(gzfile(infilename), header=FALSE,
        col.names=c("chr", "C2", "mid", "hitCount", "C5")))
    data <- subset(data, hitCount >= mappabilityThreshold)

    # list of all chromosomes
    chrList <- unique(data$chr)
    # list of all mappable fragments
    listOfMappableFrags <- data.table(chr=data$chr, mid=data$mid)

    ### In Python, possiblePairsPerDistance is a dictionary(key-value pairs) ###
    ### In R, possiblePairsPerDistance is a data frame(for optimization)     ###
    ### Each column is chr, mid1, mid2, interactionDistance, hitCount, bias  ###
    possiblePairsPerDistance <- NULL

    totalNoOfFrags <- nrow(data)
    for (i in chrList) {
        fragsPerChr <- subset(data, chr == i)
        tempLen <- nrow(fragsPerChr)

        chr_mid1_data <- data.table(chr=fragsPerChr$chr, mid1=fragsPerChr$mid,
            dummy=rep(1, tempLen))
        mid2_data <- data.table(mid2=fragsPerChr$mid, dummy=rep(1, tempLen))

        intraPairs <- merge(chr_mid1_data, mid2_data, by="dummy", all=TRUE,
            allow.cartesian=TRUE)
        intraPairs <- subset(intraPairs, mid1 < mid2)
        intraPairs <- data.table(chr=intraPairs$chr, mid1=intraPairs$mid1,
            mid2=intraPairs$mid2,
            interactionDistance=abs(intraPairs$mid1 - intraPairs$mid2))

        intraInRangePairs <- subset(intraPairs,
            in_range_check(intraPairs$interactionDistance,
            distLowThres, distUpThres))
        possiblePairsPerDistance <- rbind(possiblePairsPerDistance,
            intraInRangePairs)

        possibleInterAllCount <- possibleInterAllCount +
            (totalNoOfFrags - tempLen) * tempLen
        possibleIntraInRangeCount <- possibleIntraInRangeCount +
            nrow(intraInRangePairs)
        possibleIntraAllCount <- possibleIntraAllCount +
            tempLen * (tempLen - 1) / 2
    }

    # divide the possibleInterAllCount by 2 so that
    # every inter-chr interaction is counted only once
    possibleInterAllCount <- possibleInterAllCount / 2
    # calculate inter-chr probabilities
    if (possibleInterAllCount > 0) {
        baselineInterChrProb <- 1.0 / possibleInterAllCount
    }
    baselineIntraChrProb <- 1.0 / possibleIntraAllCount

    print("Complete generate_FragPairs method [OK]")

    return(list(possibleInterAllCount=possibleInterAllCount,
        possibleIntraAllCount=possibleIntraAllCount,
        possibleIntraInRangeCount=possibleIntraInRangeCount,
        listOfMappableFrags=listOfMappableFrags,
        possiblePairsPerDistance=possiblePairsPerDistance,
        baselineInterChrProb=baselineInterChrProb,
        baselineIntraChrProb=baselineIntraChrProb))
}

read_ICE_biases <- function(infilename) {
    print("Running read_ICE_biases method ...")

    bias <- NULL

    data <- data.table(read.table(gzfile(infilename), header=FALSE,
        col.names=c("chr", "mid", "bias")))
    biasDic <- data[bias < 0.5 | bias > 2, bias := -1]

    print("Complete read_ICE_biases method [OK]")

    return(biasDic)
}

read_All_Interactions <- function(infilename, biasDic, listOfMappableFrags,
possiblePairsPerDistance, distUpThres, distLowThres) {

    print("Running read_All_Interactions method ...")

    chr1 <- chr2 <- NULL

    data <- data.table(read.table(gzfile(infilename), header=FALSE,
        col.names = c("chr1", "mid1", "chr2", "mid2", "hitCount")))

    # set bias1 and bias2
    if (!is.null(biasDic)) {
        names(biasDic) <- c("chr1", "mid1", "bias1")
        data <- merge(data, biasDic, by=c("chr1", "mid1"), all.x=TRUE)

        names(biasDic) <- c("chr2", "mid2", "bias2")
        data <- merge(data, biasDic, by=c("chr2", "mid2"), all.x=TRUE)

        data[is.na(data)] <- 1.0
    } else {
        data <- data.table(data, bias1=rep(1.0, nrow(data)),
            bias2=rep(1.0, nrow(data)))
    }

    # filter by listOfMappableFrags
    names(listOfMappableFrags) <- c("chr1", "mid1")
    data <- merge(data, listOfMappableFrags, by=c("chr1", "mid1"), all=FALSE)

    names(listOfMappableFrags) <- c("chr2", "mid2")
    data <- merge(data, listOfMappableFrags, by=c("chr2", "mid2"), all=FALSE)

    names(listOfMappableFrags) <- c("chr", "mid")

    # inter
    inter_data <- subset(data, chr1 != chr2)
    observedInterAllSum <- sum(inter_data$hitCount)
    observedInterAllCount <- nrow(inter_data)

    # intra
    intra_data <- subset(data, chr1 == chr2)
    observedIntraAllSum <- sum(intra_data$hitCount)
    observedIntraAllCount <- nrow(intra_data)

    # intraInRange
    intraInRange_data <- data.table(chr=intra_data$chr1,
        mid1=pmin(intra_data$mid1, intra_data$mid2),
        mid2=pmax(intra_data$mid1, intra_data$mid2),
        interactionDistance=abs(intra_data$mid1 - intra_data$mid2),
        hitCount=intra_data$hitCount,
        bias=intra_data$bias1 * intra_data$bias2)
    intraInRange_data <- subset(intraInRange_data,
        in_range_check(intraInRange_data$interactionDistance, distLowThres,
        distUpThres))
    minObservedGenomicDist <- min(intraInRange_data$interactionDistance)
    maxObservedGenomicDist <- max(intraInRange_data$interactionDistance)

    observedIntraInRangeSum <- sum(intraInRange_data$hitCount)
    observedIntraInRangeCount <- nrow(intraInRange_data)

    temp <- nrow(possiblePairsPerDistance)
    possiblePairsPerDistance <- merge(possiblePairsPerDistance,
        intraInRange_data, by=c("chr", "mid1", "mid2", "interactionDistance"),
        all=TRUE, allow.cartesian=TRUE)
    if (temp < nrow(possiblePairsPerDistance)) {
        print("Illegal fragment pair")
        q()
    }

    hitCount_data <- possiblePairsPerDistance$hitCount
    hitCount_data[is.na(hitCount_data)] <- 0

    bias_data <- possiblePairsPerDistance$bias
    bias_data[is.na(bias_data)] <- 1.0

    possiblePairsPerDistance <- data.table(chr=possiblePairsPerDistance$chr,
        mid1=possiblePairsPerDistance$mid1, mid2=possiblePairsPerDistance$mid2,
        interactionDistance=possiblePairsPerDistance$interactionDistance,
        hitCount=hitCount_data, bias=bias_data)

    sortedInteractions <- data.table(
        interactionDistance=possiblePairsPerDistance$interactionDistance,
        hitCount=possiblePairsPerDistance$hitCount,
        bias=possiblePairsPerDistance$bias)

    print("Complete read_All_Interactions method [OK]")

    return(list(sortedInteractions=
        sortedInteractions[order(sortedInteractions$interactionDistance), ],
        observedInterAllSum=observedInterAllSum,
        observedIntraAllSum=observedIntraAllSum,
        observedIntraInRangeSum=observedIntraInRangeSum,
        observedInterAllCount=observedInterAllCount,
        observedIntraAllCount=observedIntraAllCount,
        observedIntraInRangeCount=observedIntraInRangeCount,
        minObservedGenomicDist=minObservedGenomicDist,
        maxObservedGenomicDist=maxObservedGenomicDist,
        possiblePairsPerDistance=possiblePairsPerDistance))
}

calculate_Probabilities <- function(sortedInteractions, isOutlier, figname,
noOfBins, useBinning, distScaling, observedIntraInRangeSum, outdir, visual,
libname, toKb, toProb) {

    print("Running calculating_Probabilities method ...")

    desiredPerBin <- observedIntraInRangeSum / noOfBins

    data <- calculate_probabilities_helper(sortedInteractions, isOutlier,
        useBinning, desiredPerBin, distScaling, observedIntraInRangeSum)

    if (!file.exists(outdir)) {
        dir.create(outdir, recursive=TRUE)
    }

    if (visual) {
        print(paste("Plotting ", figname, ".png", sep=""))
        png(filename=paste(outdir, "/", figname, ".png", sep=""),
            width=800, height=600)
        titleStr <- paste(
            "Binning observed interactions using equal occupancy bins.\n",
            "No. of bins: ", noOfBins, ", Library: ", libname,
            ", No. of interactions: ", observedIntraInRangeSum, sep="")
        plot(data$x * toKb, data$y * toProb, pch=21, col="black", bg="red",
            cex=1.5, xlab="Genomic distance (kb)",
            ylab=expression(paste("Contact probability (x", 10^-5, ")",
            sep="")), main=titleStr)
        arrows(data$x * toKb, data$y * toProb - data$yerr * toProb,
            data$x * toKb, data$y * toProb + data$yerr * toProb,
            length=0.05, angle=90, code=3)
        legend("topright", legend=c("Mean", "Standard error"), pch=c(21, 91),
            pt.bg=c("red", "black"))
        dev.off()
    }

    print(paste("Writing ", figname, ".txt", sep=""))
    file.create(paste(outdir, "/", figname, ".txt", sep=""))
    outputData <- data.table(avgGenomicDist=trunc(data$x),
        contactProbability=data$y, standardError=data$yerr,
        noOfLocusPairs=data$pairCounts,
        totalOfContactCounts=data$interactionTotals)
    write.table(format(outputData, digits=3),
        file=paste(outdir, "/", figname, ".txt", sep=""), quote=FALSE, sep="\t",
        row.names=FALSE, col.names=TRUE)

    print("Complete calculating_Probabilities method [OK]")

    return(data.table(x=data$x, y=data$y, yerr=data$yerr))
}

fit_Spline <- function(x, y, yerr, infilename, sortedInteractions, biasDic,
figname, passNo, outdir, visual, distLowThres, distUpThres, toKb, toProb,
useInters, baselineIntraChrProb, baselineInterChrProb, observedInterAllSum,
observedIntraAllSum, observedIntraInRangeSum, possibleInterAllCount,
possibleIntraAllCount, possibleIntraInRangeCount, maxObservedGenomicDist) {

    print("Running fit_Spline method ...")

    bias <- bias1 <- bias2 <- chr1 <- chr2 <- hitCount <-
        interactionDistance <- isOutlier <- p_val <- NULL

    ius <- smooth.spline(x, y)

    tempMaxX <- max(x)
    tempMinX <- min(x)
    tempList <- unique(trunc(sortedInteractions$interactionDistance))

    splineX <- tempList[tempList >= tempMinX & tempList <= tempMaxX]
    splineY <- predict(ius, splineX)$y

    newSplineY <- monoreg(splineX, splineY, type="antitonic")$yf

    if (!file.exists(outdir)) {
        dir.create(outdir, recursive=TRUE)
    }

    if (visual) {
        print(paste("Plotting ", figname, ".png", sep=""))

        png(filename=paste(outdir, "/", figname, ".png", sep=""),
            width=800, height=600)

        par(mfrow=c(2, 1))

        if (distLowThres > -1 & distUpThres > -1) {
            plot(splineX * toKb, newSplineY * toProb, type="l", col="green",
                xlim=c(distLowThres * toKb, distUpThres * toKb),
                xlab="Genomic distance (kb)",
                ylab=expression(paste("Contact probability (x", 10^-5, ")",
                sep="")))
        } else {
            plot(splineX * toKb, newSplineY * toProb, type="l", col="green",
                xlab="Genomic distance (kb)",
                ylab=expression(paste("Contact probability (x", 10^-5, ")",
                sep="")))
        }
        arrows(x * toKb, y * toProb - yerr * toProb, x * toKb,
            y * toProb + yerr * toProb, col="red", length=0.05, angle=90,
            code=3)
        if (useInters) {
            lines(x * toKb, rep(baselineIntraChrProb, length(x)) * toProb,
                col="black")
            lines(x * toKb, rep(baselineInterChrProb, length(x)) * toProb,
                col="blue")
            legend("topright", legend=c(paste("spline-", passNo, sep=""),
                "Mean with std. error", "Baseline intra-chromosomal",
                "Baseline inter-chromosomal"), pch=c(NA, 91, NA, NA),
                lty=c(1, NA, 1, 1), col=c("green", "red", "black", "blue"))
        } else {
            legend("topright", legend=c(paste("spline-", passNo, sep=""),
                "Mean with std. error"), pch=c(NA, 91),
                lty=c(1, NA), col=c("green", "red"))
        }

        if (distLowThres > -1 & distUpThres > -1) {
            plot(splineX, newSplineY, type="l", log="xy", col="green",
                xlim=c(distLowThres, distUpThres),
                xlab="Genomic distance (log-scale)",
                ylab="Contact probability (log-scale)")
        } else {
            plot(splineX, newSplineY, type="l", log="xy", col="green",
                xlab="Genomic distance (log-scale)",
                ylab="Contact probability (log-scale)")
        }
        arrows(x, y - yerr, x, y + yerr, col="red", length=0.05, angle=90,
            code=3)
        if (useInters) {
            lines(x, rep(baselineIntraChrProb, length(x)), col="black")
            lines(x, rep(baselineInterChrProb, length(x)), col="blue")
        }

        dev.off()
    }

    data <- data.table(read.table(gzfile(infilename), header=FALSE,
        col.names = c("chr1", "mid1", "chr2", "mid2", "hitCount")))

    if (!is.null(biasDic)) {
        names(biasDic) <- c("chr1", "mid1", "bias1")
        data <- merge(data, biasDic, by=c("chr1", "mid1"), all.x=TRUE)

        names(biasDic) <- c("chr2", "mid2", "bias2")
        data <- merge(data, biasDic, by=c("chr2", "mid2"), all.x=TRUE)

        data[is.na(data)] <- 1.0
    } else {
        data <- data.table(data, bias1=rep(1.0, nrow(data)), bias2=rep(1.0,
            nrow(data)))
    }

    ### chr1 mid1 chr2 mid2 interactionDistance hitCount bias1 bias2 p_val ###
    data <- data.table(chr1=data$chr1, mid1=data$mid1, chr2=data$chr2,
        mid2=data$mid2, interactionDistance=abs(data$mid1 - data$mid2),
        hitCount=data$hitCount, bias1=data$bias1, bias2=data$bias2,
        p_val=rep(Inf, nrow(data)))
    data[chr1 == chr2 & useInters, p_val := 1 - pbinom(hitCount - 1,
        observedInterAllSum, baselineInterChrProb * bias1 * bias2)]
    data[chr1 == chr2 & interactionDistance > distUpThres,
        p_val := 1 - pbinom(hitCount - 1, observedIntraAllSum, 1)]
    data[chr1 == chr2 & interactionDistance <= distLowThres, p_val := 1]
    data[chr1 == chr2 & in_range_check(interactionDistance,
        distLowThres, distUpThres), p_val := 1 - pbinom(hitCount - 1,
        observedIntraInRangeSum, newSplineY[pmin(bisect_left(splineX,
        pmin(pmax(interactionDistance, tempMinX), tempMaxX)),
        length(splineX))] * bias1 * bias2)]
    data[chr1 == chr2 & (bias1 < 0 | bias2 < 0), p_val := 1]

    p_vals <- data$p_val
    p_vals <- p_vals[p_vals != Inf]

    if (useInters) {
        q_vals <- benjamini_hochberg_correction(p_vals,
            possibleInterAllCount + possibleIntraAllCount)
    } else {
        q_vals <- benjamini_hochberg_correction(p_vals,
            possibleIntraInRangeCount)
    }

    print(paste("Writing p-values to file ", figname, ".significances.txt.gz",
        sep=""))
    file.create(paste(outdir, "/", figname, ".significances.txt.gz", sep=""))
    gz <- gzfile(paste(outdir, "/", figname, ".significances.txt.gz", sep=""))
    outputData <- data.table(chr1=data$chr1, fragmentMid1=data$mid1,
        chr2=data$chr2, fragmentMid2=data$mid2, contactCount=data$hitCount,
        p_value=p_vals, q_value=q_vals)
    if (useInters) {
        outputData <- subset(outputData, chr1 != chr2)
    } else {
        outputData <- subset(outputData, chr1 == chr2 &
            in_range_check(abs(outputData$fragmentMid1 -
            outputData$fragmentMid2), distLowThres, distUpThres))
    }
    write.table(outputData, gz, quote=FALSE, sep="\t", row.names=FALSE,
        col.names=TRUE)

    outlierThres <- 1 / possibleIntraInRangeCount
    tempData <- data.table(
        interactionDistance=sortedInteractions$interactionDistance,
        hitCount=sortedInteractions$hitCount,
        bias=sortedInteractions$bias,
        p_val=rep(Inf, nrow(sortedInteractions)),
        isOutlier=rep(0, nrow(sortedInteractions)))
    tempData[, p_val := 1 - pbinom(hitCount - 1, observedIntraInRangeSum,
        newSplineY[pmin(bisect_left(splineX, pmin(pmax(interactionDistance,
        tempMinX), tempMaxX)), length(splineX))] * bias)]
    tempData[p_val < outlierThres, isOutlier := 1]

    belowData <- subset(tempData, isOutlier == 1)
    aboveData <- subset(tempData, isOutlier == 0)

    distsBelow <- belowData$interactionDistance
    distsAbove <- aboveData$interactionDistance

    intcountsBelow <- belowData$hitCount
    intcountsAbove <- aboveData$hitCount

    belowThresCount <- nrow(belowData)
    aboveThresCount <- nrow(aboveData)

    if (visual) {
        print(paste("Plotting results of extracting outliers to file ", figname,
            ".extractOutliers.png", sep=""))

        png(filename=paste(outdir, "/", figname, ".extractOutliers.png",
            sep=""), width=800, height=600)

        downsample <- 30
        randIndcsAbove <- sample(seq(1, length(intcountsAbove)),
            trunc(length(intcountsAbove) / downsample))
        randIndcsAbove <- randIndcsAbove[order(randIndcsAbove)]
        downsample <- 20
        randIndcsBelow <- sample(seq(1, length(intcountsBelow)),
            trunc(length(intcountsBelow) / downsample))
        randIndcsBelow <- randIndcsBelow[order(randIndcsBelow)]

        if (length(intcountsBelow) > 0 &
            (distLowThres > -1 & distUpThres > -1)) {
            plot(distsBelow[randIndcsBelow] * toKb,
                intcountsBelow[randIndcsBelow], pch=20, col="red",
                xlim=c(0, distUpThres * toKb),
                ylim=c(0, min(max(intcountsBelow), 1500)),
                xlab="Genomic distance (kb)", ylab="Contact counts")
        } else if (length(intcountsBelow) > 0) {
            plot(distsBelow[randIndcsBelow] * toKb,
                intcountsBelow[randIndcsBelow], pch=20, col="red",
                ylim=c(0, min(max(intcountsBelow), 1500)),
                xlab="Genomic distance (kb)", ylab="Contact counts")
        } else if (distLowThres > -1 & distUpThres > -1) {
            plot(distsBelow[randIndcsBelow] * toKb,
                intcountsBelow[randIndcsBelow], pch=20, col="red",
                xlim=c(0, distUpThres * toKb), xlab="Genomic distance (kb)",
                ylab="Contact counts")
        } else {
            plot(distsBelow[randIndcsBelow] * toKb,
                intcountsBelow[randIndcsBelow], pch=20, col="red",
                xlab="Genomic distance (kb)", ylab="Contact counts")
        }
        lines(append(splineX, maxObservedGenomicDist) * toKb, append(newSplineY,
            newSplineY[length(newSplineY)]) * observedIntraInRangeSum,
            col="green")
        legend("topright", legend=c("Outliers (p-value < 1/M)",
            paste("spline-", passNo, " (x N)", sep="")), pch=c(20, NA),
            lty=c(NA, 1), col=c("red", "green"))

        dev.off()
    }

    if (visual) {
        print(paste("Plotting q-values to file ", figname, ".qplot.png",
            sep=""))
    }
    minFDR <- 0
    maxFDR <- 0.05
    increment <- 0.001
    plot_qvalues(q_vals, minFDR, maxFDR, increment, paste(figname, ".qplot",
        sep=""), outdir, visual)

    print("Complete fit_Spline method [OK]")

    return(tempData$isOutlier)
}

plot_qvalues <- function(q_values, minFDR, maxFDR, increment, figname, outdir,
visual) {

    qvalTicks <- seq(minFDR, maxFDR + increment - increment / 2, increment)
    qvalBins <- floor(q_values / increment)

    data <- qvalBins[qvalBins < length(qvalTicks)] + 1
    significantTicks <- as.data.frame(table(append(seq(1, length(qvalTicks)),
        data)))$Freq - 1

    # make it cumulative
    significantTicks <- cumsum(significantTicks)
    # shift them by 1
    significantTicks <- append(0,
        significantTicks[1 : length(significantTicks) - 1])

    if (!file.exists(outdir)) {
        dir.create(outdir, recursive=TRUE)
    }

    if (visual) {
        png(filename=paste(outdir, "/", figname, ".png", sep=""),
            width=800, height=600)
        plot(qvalTicks, significantTicks, type="o", pch=17, col="blue", lty=1,
            xlab="FDR threshold", ylab="Significant contacts")
        dev.off()
    }
}

in_range_check <- function(interactionDistance, distLowThres, distUpThres) {
    return(interactionDistance > distLowThres &
        interactionDistance <= distUpThres)
}

bisect_left <- function(a, x) {
    return(length(a) + 1 - findInterval(-x, sort(-a)))
}

benjamini_hochberg_correction <- function(p_values, num_total_tests) {
    order <- order(p_values)
    sorted_pvals <- p_values[order]
    return(benjamini_hochberg_correction_helper(p_values, num_total_tests,
        sorted_pvals, order - 1))
}
