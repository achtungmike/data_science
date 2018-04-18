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
