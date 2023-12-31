% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SNVAnalysis.R
\name{varianceCCAnalysisPheno}
\alias{varianceCCAnalysisPheno}
\title{Variance Case-Control Analysis where the Case and Control are User Provided}
\usage{
varianceCCAnalysisPheno(
  variants,
  totalCaseSamples,
  phenotypeName,
  useChi = FALSE
)
}
\arguments{
\item{variants}{The data.frame containing the variant information.}

\item{totalCaseSamples}{The total number of case samples that goes into
the variants data.frame. In other words, this should equal the number of rows
found in the metadata file used to load the VCFs as each row in that file
should represent one case.}

\item{phenotypeName}{The name of the phenotype to categorize upon. E.g., if
each VCF had trinary column named "Type" and the values in the column are of
A, B, or C, then the categorical values this function would discover are A,
B, and C and therefore, would require the breaking down of the tables pairwise.}

\item{useChi}{Instead of using Fisher's exact test, the chi-square test will be employed.
Default is FALSE.}
}
\value{
A data.frame containing all the p-values and the shape of the 2-way table
employed for each test.
}
\description{
Similar to varianceCCAnalysisEnsembl, however, instead of using Ensembl,
the name of a column from the metadata used to load the VCFs is taken.
The function then finds all possible categorical values for the column
and built pairwise tables based off these.
}
