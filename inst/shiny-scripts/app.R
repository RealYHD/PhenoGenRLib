library(shiny)
library(shinyWidgets)
library(PhenoGenRLib)
library(readr)
library(DT)

# Define UI for app that draws a histogram ----
ui <- fluidPage(

  # App title ----
  titlePanel("Variant Visualization from PhenoGenRLib"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      shiny::fileInput(
        inputId = "variantsMetadata",
        label = "A CSV of the variants containing the reference and alternate
        alleles for a given position with headers represented by VCF columns.
        Easily obtained by linking variants with a metadata file via function
        linkVariantsWithMetadata.",
        accept = c(".csv"),
        multiple = FALSE
      ),
      shiny::textInput(
        inputId = "variantsMetadataVCFCol",
        label = "VCF Filename Column",
        value = "vcfs"
      ),
      shiny::fileInput(
        inputId = "variants",
        label = "Variant Call Format (VCF) files to be used defined in the
        variants metadata file previously uploaded.",
        accept = c(".vcf"),
        multiple = TRUE
      ),
      shiny::checkboxGroupInput(
        inputId = "nucleoBases",
        label = "The nucleotide bases to sum the frequencies. Selecting just
        one would show the frequency of the single base.",
        choices = c("A", "T", "C", "G"),
        selected = c("A", "T", "C", "G"),
        inline = TRUE
      ),
      shiny::textInput(
        inputId = "chromCol",
        label = "Chromosome Column Name",
        value = "chromosome"
      ),
      shiny::numericInput(
        inputId = "posOffset",
        label = "The Offset of the Position",
        value = 0
      ),
      shiny::numericInput(
        inputId = "hgVersion",
        label = "Human Genome Version",
        value = 38,
        min = 37,
        max = 38
      ),
      shiny::checkboxInput(
        inputId = "ccWithEnsembl",
        label = "Case-control variants with Ensembl Data",
        value = FALSE
      )
    ),

    # Main panel for displaying outputs ----
    mainPanel(
      DT::dataTableOutput(outputId = "metadataTable"),
      DT::dataTableOutput(outputId = "variantsTable"),
      shiny::plotOutput(outputId = "variancePlot"),
      DT::dataTableOutput(outputId = "associationTable")
    )
  )
)

# Define server logic required to draw a histogram ----
server <- function(input, output) {
  metadata <- shiny::reactive({
    shiny::req(input$variantsMetadata)
    metadata <- readr::read_csv(input$variantsMetadata$datapath) # Read the metadata file
    return(metadata)
  })

  metadataWithPaths <- shiny::reactive({
    shiny::req(metadata)
    shiny::req(input$variantsMetadataVCFCol)
    shiny::req(input$variants)

    if (input$variantsMetadataVCFCol %in% colnames(metadata())) {
      metadata <- base::data.frame(metadata())
      metadata[[input$variantsMetadataVCFCol]] <- lapply( # Update the vcf columns to the correct paths
        metadata[[input$variantsMetadataVCFCol]],
        function(vcf) {
          path <- input$variants$datapath[input$variants$name == vcf]
          return(path)
        }
      )
    } else {
      stop(base::paste("No column named \"", input$variantsMetadataVCFCol, "\" found in: ", list(base::colnames(metadata()))), sep = "")
    }
    return(metadata)
  })

  variants <- shiny::reactive({
    shiny::req(metadataWithPaths)
    shiny::req(input$variantsMetadataVCFCol)

    shiny::withProgress(message = "Loading Variants and Linking Metadata", {
      variants <- PhenoGenRLib::linkVariantsWithMetadata( # Generate variants
        metadata = metadataWithPaths(),
        vcfColName = input$variantsMetadataVCFCol,
        progress = function(it, total) {
          shiny::setProgress(value = it/total)
        }
      )
      return(variants)
    })
  })

  mappedVariants <- shiny::reactive({
    shiny::req(variants)
    shiny::req(input$chromCol)
    shiny::req(input$posOffset)
    shiny::req(input$hgVersion)

    shiny::withProgress(message = "Mapping RSIDs", {
      mappedVariants <- PhenoGenRLib::mapRsidsForVariants(
        chromCol = input$chromCol,
        variants = variants(),
        offset = input$posOffset,
        hostGenVersion = input$hgVersion,
        progress = function(it, total) {
          shiny::setProgress(value = it/total)
        }
      )
      return(mappedVariants)
    })
  })

  heatmap <- shiny::reactive({
    shiny::req(variants)
    heatmap <- PhenoGenRLib::generatePositionHeatmap(
      variants = variants()
    )
    return(heatmap)
  })

  ensemblAssociationTable <- shiny::reactive({
    shiny::req(mappedVariants)
    shiny::req(input$ccWithEnsembl)

    shiny::withProgress(message = "Calculating Association P-Values", {
      ensemblAssocTable <- PhenoGenRLib::varianceCCAnalysisEnsembl(
        mappedVariants()$nvs,
        mappedVariants()$rsids,
        base::nrow(metadata()),
        progress = function(it, total) {
          shiny::setProgress(value = it/total)
        }
      )
      return(ensemblAssocTable)
    })
  })

  output$metadataTable <- DT::renderDataTable({
    shiny::req(metadata)
    return(metadata())
  })

  output$variantsTable <- DT::renderDataTable({
    shiny::req(variants)
    return(mappedVariants()$nvs[c(
      "CHROM",
      "POS",
      "REF",
      "ALT",
      "local_id",
      "refsnp_id"
    )])
  })

  output$variancePlot <- shiny::renderPlot({
    shiny::req(heatmap)

    plot <- PhenoGenRLib::visVariantHeatMap(
      heatmap = heatmap(),
      nucbases = input$nucleoBases
    )
    return(plot)
  })

  output$associationTable <- DT::renderDataTable({
    shiny::req(ensemblAssociationTable)
    return(ensemblAssociationTable())
  })

}

shinyApp(ui = ui, server = server)
