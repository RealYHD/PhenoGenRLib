
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PhenoGenRLib

<!-- badges: start -->
<!-- badges: end -->

The goal of PhenoGenRLib is to simplify nucleotide variant analysis.

As next generational sequencing (NGS) begins taking off, more and more
data is readily available to be used. Arguably, there is an
overabundance of data that has yet been used to it’s fullest potential.
PhenoGenRLib promises to provide simple ways of loading VCFs,
associating them with sample metadata, and lastly, running associative
studies by applying the metadata.

## Installation

You can install the development version of PhenoGenRLib like so:

``` r
require("devtools")
devtools::install github("RealYHD/PhenoGenRLib",
build vignettes = TRUE)
library("PhenoGenRLib")
```

## Getting Started

To get started, have a datasheet ready in the form of a `CSV`. This
datasheet should at the very least, contain one column, where each row
in that column contains the filename of the VCF including the `.vcf`.
For the following example, we will assume that such a file is called
`huntingtons_datasheet_shortened.csv` and is located at
`./inst/extdata/huntingtons_datasheet_shortened.csv` with the column
containing the VCF filenames being named `vcfs`. We will also need the
location of the VCFs. Let’s assume they can be found at the same place
as the metadata CSV `./inst/extdata/`. Then:

``` r
library(PhenoGenRLib)
variants <- PhenoGenRLib::linkVariantsWithMetadata(
  metadata = "inst/extdata/huntingtons_datasheet_shortened.csv",
  vcfDir = "inst/extdata/",
  vcfColName = "vcfs"
)
#> Rows: 8 Columns: 3
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: ","
#> chr (2): vcfs, dummy_pheno
#> dbl (1): chromosome
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
#> READING VCF
#>  * checking if file exists... PASS
#>  * Reading vcf header...
#>    Done
#>  * Reading vcf body...
#>    Done
#>  * Parse vcf header...
#>    Done
#>  * Split info...
#>  * Done
#>  * Split samples...
#>    Done
#> READING VCF
#>  * checking if file exists... PASS
#>  * Reading vcf header...
#>    Done
#>  * Reading vcf body...
#>    Done
#>  * Parse vcf header...
#>    Done
#>  * Split info...
#>  * Done
#>  * Split samples...
#>    Done
#> READING VCF
#>  * checking if file exists... PASS
#>  * Reading vcf header...
#>    Done
#>  * Reading vcf body...
#>    Done
#>  * Parse vcf header...
#>    Done
#>  * Split info...
#>  * Done
#>  * Split samples...
#>    Done
#> READING VCF
#>  * checking if file exists... PASS
#>  * Reading vcf header...
#>    Done
#>  * Reading vcf body...
#>    Done
#>  * Parse vcf header...
#>    Done
#>  * Split info...
#>  * Done
#>  * Split samples...
#>    Done
#> READING VCF
#>  * checking if file exists... PASS
#>  * Reading vcf header...
#>    Done
#>  * Reading vcf body...
#>    Done
#>  * Parse vcf header...
#>    Done
#>  * Split info...
#>  * Done
#>  * Split samples...
#>    Done
#> READING VCF
#>  * checking if file exists... PASS
#>  * Reading vcf header...
#>    Done
#>  * Reading vcf body...
#>    Done
#>  * Parse vcf header...
#>    Done
#>  * Split info...
#>  * Done
#>  * Split samples...
#>    Done
#> READING VCF
#>  * checking if file exists... PASS
#>  * Reading vcf header...
#>    Done
#>  * Reading vcf body...
#>    Done
#>  * Parse vcf header...
#>    Done
#>  * Split info...
#>  * Done
#>  * Split samples...
#>    Done
#> READING VCF
#>  * checking if file exists... PASS
#>  * Reading vcf header...
#>    Done
#>  * Reading vcf body...
#>    Done
#>  * Parse vcf header...
#>    Done
#>  * Split info...
#>  * Done
#>  * Split samples...
#>    Done
```

Checkout the documents and vignettes for where to go from here!

# Contributions

PhenoGenRLib stands on the shoulder of giants, and it would be a
disservice to not name them:

- Thank you Syed Haider et al. for providing bedr. It was greatly
  helpful in simplifying the data ingress features.
- Thanks to the entire Biomart Team for providing an awesome and easy to
  use interface to large public databases!
- ggplot2 was very helpful in generating figures. Thanks to Wickham et.
  al!
- This entire project wouldn’t have been possible without the help of
  the TidyVerse team. Despite not using every single package from their
  library, much work and diagnostic made use of their tools.
- Tibble helped simplify data storage and accession. Thanks Muller et.
  al!

No generative AI was used for this project directly, however, learning
about how R works and how some of the syntax differs from other
languages was aided by ChatGPT.

This was a BCB410H1 UofT Bioinformatics project by Harrison Deng.

# Citations

    Müller K, Wickham H (2023). _tibble: Simple Data Frames_. R package
      version 3.2.1, <https://CRAN.R-project.org/package=tibble>.

    Haider S, Waggott D, C. Boutros P (2019). _bedr: Genomic Region
      Processing using Tools Such as 'BEDTools', 'BEDOPS' and 'Tabix'_. R
      package version 1.0.7, <https://CRAN.R-project.org/package=bedr>.

    BioMart and Bioconductor: a powerful link between biological
      databases and microarray data analysis. Steffen Durinck, Yves Moreau,
      Arek Kasprzyk, Sean Davis, Bart De Moor, Alvis Brazma and Wolfgang
      Huber, Bioinformatics 21, 3439-3440 (2005).

    H. Wickham. ggplot2: Elegant Graphics for Data Analysis.
      Springer-Verlag New York, 2016.

    Wickham H, Hester J, Bryan J (2023). _readr: Read Rectangular Text
      Data_. R package version 2.1.4,
      <https://CRAN.R-project.org/package=readr>.

# Acknowledgements

This package was developed as part of an assessment for 2023 BCB410H:
Applied Bioinformatics course at the University of Toronto, Toronto,
CANADA. PhenoGenRLib welcomes issues, enhancement requests, and other
contributions. To submit an issue, use the GitHub issues.
