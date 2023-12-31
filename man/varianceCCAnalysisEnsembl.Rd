% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SNVAnalysis.R
\name{varianceCCAnalysisEnsembl}
\alias{varianceCCAnalysisEnsembl}
\title{Variance Case-Control Analysis where the Ensembl Database is the Control}
\usage{
varianceCCAnalysisEnsembl(variants, rsids, totalCaseSamples, useChi = FALSE)
}
\arguments{
\item{variants}{The data.frame containing all the variants, variant information, and rsIDs.}

\item{rsids}{The data.frame containing all the extra information regarding the known
variants.}

\item{totalCaseSamples}{The total number of case samples. This should be the number of rows
in the metadata file as each row represents one case, hence, one VCF.}

\item{useChi}{Instead of using Fisher's exact test, the chi-square test will be employed.
Default is FALSE.}
}
\value{
A data.frame containing all the p-values and the shape of the 2-way table
employed for each test.
}
\description{
This function makes it easy to test for potential associations between the
studied case, and the average population. To elaborate, this function takes
a data.frame of all the variants, and the rsIDs found to be linked with said
variants, and uses the pre-existing allele frequency information from the rsID
data.frame (fetched from Ensembl). Then, it moves over each mutation with population
information in public databases and runs a Fisher's exact test. The Fisher's exact
test is beneficial for low sample counts as it is an exact measure of correlation in
contrast to various other tests. It is more computational intensive however. Hence,
this function also allows the use of the chi-square test. The table formed is between
whether or not the position is mutated, and whether or not the mutation occurred in a
healthy individual (the control group), or the case study individual.
}
