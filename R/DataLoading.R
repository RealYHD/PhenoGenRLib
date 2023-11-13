library(readr)
library(tibble)
library(VariantTools)
library(gmapR)

linkVariantsWithMetadata <- function(metadataFile, bamDir, refSeq) {
    metadata = readr::read_csv(metadataFile)
    for (bamFile in metadata["filename"]) {
        fastaFile <- rtracklayer::FastaFile(base::normalizePath(refSeq))
        gmapGenome <- gmapR::GmapGenome(fastaFile, create=TRUE)
        tally.param <- VariantTools::TallyVariantsParam(gmapGenome, high_base_quality = 56L)
        variants = VariantTools::callVariants(base::normalizePath(bamFile), tally.param=tally.param)
        # TODO add to tibble containing all known variants
    }
    # TODO return the tibble containing all known variants 
}
