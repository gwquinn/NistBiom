# This software was developed at the National Institute of Standards and Technology (NIST) by
# employees of the Federal Government in the course of their official duties. Pursuant to title 17
# Section 105 of the United States Code, this software is not subject to copyright protection and
# is in the public domain. NIST assumes no responsibility whatsoever for its use by other parties,
# and makes no guarantees, expressed or implied, about its quality, reliability, or any other
# characteristic.

colors.15 <- c("blue4", "darkred", "orange3", "steelblue3", "purple3",
               "darkgreen", "red", "darkcyan", "orange4", "pink3",
               "grey40", "yellow3", "black", "olivedrab4", "seagreen1")

#' Constructs a DET object from a set of mated and nonmated scores.
#'
#' If num.points is a single numerical value, then it specifies the number of decision
#' thresholds to use when generating the DET curve. The selected thresholds will
#' correspond to evenly spaced FMRs in the log space. If num.points is a numeric vector,
#' then the elements of the vector will be used as the decision thresholds when
#' generating the DET curve.
#'
#' Failure-to-enroll (FTE) comparisons should be represented by NAs. When positive is
#' set to true, FTE comparisons contribute to the false non-match rate. When
#' positive is set to false, FTE comparisons contribute to the false match rate.
#'
#' @param gen.scores genuine (aka mated) comparison scores.
#' @param nonmated.scores impostor (aka nonmated) comparison scores.
#' @param FMR.range range over which to construct the DET.
#' @param num.points The number of decision thresholds to use when generating the DET
#'                   curve.
#' @param similarity Are the scores measures of similarity (TRUE) or distances (FALSE).
#' @param positive Is it a positive or negative system (determines how FTEs are treated).
#' @return A DET curve as a data.table with three columns: "t", "FMR", and "FNMR".
#' @export
#' @examples
#' gen.scores <- rnorm(1e5, mean=4)
#' nonmated.scores <- rnorm(1e5, mean=0)
#' DET <- make.DET(gen.scores, nonmated.scores)
#' get.FNMR(DET, FMR=1e-2)
make.DET <- function(mated.scores, nonmated.scores, FMR.range=c(1e-8, 1), num.points=4096,
                     similarity="auto", positive.system=TRUE, one.to.one=TRUE)
{
   # If caller wants us to figure out whether scores are similarities or distances.
   if (similarity == "auto")
   {
      mated    <- median(mated.scores,    na.rm=TRUE)
      nonmated <- median(nonmated.scores, na.rm=TRUE)
     
      # Assume scores are similarity measures if median is greater for mated scores.
      similarity <- mated > nonmated
   }

   # Convert dissimilarity scores to similarity scores.
   if (!similarity)
   {
      mated.scores    <- -mated.scores
      nonmated.scores <- -nonmated.scores
   }

   # Positive systems translate FTEs into low scores, negative systems into high scores.
   FTE.score <- ifelse(positive.system, -Inf, Inf)

   mated.scores[is.na(mated.scores)]       <- FTE.score
   nonmated.scores[is.na(nonmated.scores)] <- FTE.score

   # If thresholds specified explicitly.
   if (length(num.points) > 1)
   {
      thresholds <- num.points
      
      if (!similarity)
         thresholds <- -thresholds
   }
   else
   {
      # Minimum FMR cannot be lower than the FTE rate.
      if (sum(nonmated.scores == FTE.score))
         FMR.range[1] <- max(FMR.range[1], mean(nonmated.scores == FTE.score))
      
      # Find thresholds corresponding to equidistant FMRs on a log scale.
      FMR.range <- log10(FMR.range)
      quantiles <- 10^seq(FMR.range[1], FMR.range[2], length=num.points)
      
      thresholds <- quantile(nonmated.scores, quantiles, type=1)

      # Restrict to unique thresholds and cap off ends.
      thresholds <- unique(c(-Inf, sort(thresholds), Inf))
   }

   # Construct the DET (minus to make correct strict decisions).
   DET <- data.table(t    = as.double(thresholds),
                     FMR  = ecdf(-nonmated.scores)(-thresholds),
                     FNMR = 1 - ecdf(-mated.scores)(-thresholds))

   # Flip threshold sign if scores are distances.
   if (!similarity)
      DET[,t := -t]

   # Discard duplicate points.
   DET <- unique(DET, by=c("FMR", "FNMR"))
 
   # Change column names if requesting one-to-many DET. 
   if (!one.to.one)
      setnames(DET, c("FMR", "FNMR"), c("FPIR", "FNIR"))
   
   DET
}

# Similarity if FMR is positively correlated with threshold
similarity <- function(DET)
{
   col <- intersect(c("FMR", "FPIR"), names(DET))
   cor(DET[[col]], DET$t, method='kendall') < 0
}

#' Get FNMR at the specified FMR(s) or threshold(s).
#'
#' Must specify either FMR or t but not both.
#'
#' @param DET A data.table object or path to a DET file.
#' @param FMR A single false match rate (or vector of FMRs).
#' @param t A single decision threshold (or vector of thresholds).
#' @return A vector of false non-match rates.
#' @export
get.FNMR <- function(DET, ...)
{
   # If first argument is a path to a DET file, load it.
   if (is.character(DET))
      DET <- fread(DET)
   else
      DET <- copy(DET) # avoid altering the original passed object.
   
   args <- data.table(...)
   
   # Convert 1:N column names to 1:1 if necessary.
   setnames(args, c("FPIR", "FNIR"), c("FMR", "FNMR"), skip_absent=TRUE)
   setnames(DET,  c("FPIR", "FNIR"), c("FMR", "FNMR"), skip_absent=TRUE)

   # Set appropriate rolling join behavior.
   roll <- ifelse("t" %in% names(args) && similarity(DET), -Inf, Inf)

   DET[args, FNMR, on=names(args), roll=roll]
}

#' Get FMR at the specified threshold(s).
#'
#' @param DET A data.table object or path to a DET file.
#' @param t A single decision threshold (or vector of thresholds).
#' @return A vector of false match rates the same length as t.
#' @export
get.FMR <- function(DET, ...)
{
   # If first argument is a path to a DET file, load it.
   if ("character" %in% class(DET))
      DET <- fread(DET)
   else
      DET <- copy(DET) # avoid altering the original passed object.

   args <- list(...)

   # Convert 1:N column names to 1:1 if necessary.
   setnames(args, c("FPIR", "FNIR"), c("FMR", "FNMR"), skip_absent=TRUE)
   setnames(DET,  c("FPIR", "FNIR"), c("FMR", "FNMR"), skip_absent=TRUE)
   
   # Set appropriate rolling join behavior. 
   roll <- ifelse("t" %in% names(args) && !similarity(DET), -Inf, Inf)

   DET[args, FMR, on=names(args), roll=roll]
}

#' Get threshold at the specified FMR or FNMR.
#'
#' @param DET A data.table object or path to a DET file.
#' @param FMR A single false match rate (or vector of FMRs).
#' @param FNMR A singe false match rate (or vector of FNMRs).
#' @return A vector of thresholds.
#' @export
get.threshold <- function(DET, ...)
{
   # If first argument is a path, load the DET file.
   if (is.character(DET))
      DET <- fread(DET)

   args <- list(...)

   DET[args, t, on=names(args), roll=-Inf]
}

get.FNIR <- get.FNMR
get.FPIR <- get.FMR
