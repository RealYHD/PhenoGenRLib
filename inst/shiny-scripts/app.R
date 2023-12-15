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
        value = "vcf"
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
      )
    ),

    # Main panel for displaying outputs ----
    mainPanel(
      DT::dataTableOutput(outputId = "metadataTable"),
      DT::dataTableOutput(outputId = "variantsTable"),
      shiny::plotOutput(outputId = "variancePlot"),
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

  heatmap <- shiny::reactive({
    shiny::req(variants)
    heatmap <- PhenoGenRLib::generatePositionHeatmap(
      variants = variants()
    )
    return(heatmap)
  })

  output$metadataTable <- DT::renderDataTable({
    shiny::req(metadata)
    return(metadata())
  })

  output$variantsTable <- DT::renderDataTable({
    shiny::req(variants)
    return(variants()[c("CHROM", "POS", "REF", "ALT")])
  })

  output$variancePlot <- shiny::renderPlot({
    shiny::req(heatmap)

    plot <- PhenoGenRLib::visVariantHeatMap(
      heatmap = heatmap(),
      nucbases = input$nucleoBases
    )
    return(plot)
  })

}

shinyApp(ui = ui, server = server)
