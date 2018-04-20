# 
# Data science class
#
# University of Cincinnati/Cincinnati Children's
#
# Demonstrate GCT file upload, creating a signature, submit to iLincs for correlations
#
# TODO: Add error checking
#

library(shiny)
library(DT)
library(shinyWidgets)

shinyUI(fluidPage(

  # Application title
  titlePanel("GCT upload"),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      fileInput('file1', 'Choose GCT File', accept=c('text/gct','.gct')),   
      selectInput("variable", "Grouping Variable:", choices=c()),
      selectizeInput('group1', "Group1", choices = NULL, multiple = TRUE),
      selectizeInput('group2', "Group2", choices = NULL, multiple = TRUE),
      selectInput("id", "Identifier:", choices=c()),
      selectInput("idtype", "Identifier Type:", choices=c("GeneSymbol", "geneid")),
      selectInput("difffunction", "Differential function:", choices=c("t-test")),
      h4("Limit input genes to L1000 Set:"),
      switchInput("lk", onStatus = "success", offStatus = "danger"),      
      selectInput("top", label = h3("No. Diff Exp Genes"), 
                  choices = list("All" = 1, "50" = 2, 
                                 "100" = 3, "250" = 3,
                                 "500" = 4, "1000" = 5), 
                  selected = 1),
      
      actionButton("compute", "Compute and Submit Signature")
    ),

    # Show a plot of the generated distribution
    mainPanel(
      dataTableOutput("signature_data"),
      hr(),
      dataTableOutput("correlated_data"),
      hr(),
      dataTableOutput("gct_sample_data"),
      hr(),
      dataTableOutput("gct_probe_data")
    )
  )
))
