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
library(httr)
library(dplyr)

options(shiny.maxRequestSize=30*1024^2)

choices = c(0,50,100,250,500,100)

#
# GCT file read function
#
read.gct13 <-
  function (file) {
    h <- readLines(file,n=2)
    
    if (h[1]!="#1.3") {
      stop(paste("Unknown file format version:", h[1], sep = " "))
    }
    
    h_arr <- strsplit(h[2], '\t')
    tryCatch ({
      num_samples <- as.numeric(h_arr[[1]][1])
      num_probes <- as.numeric(h_arr[[1]][2])
      num_probe_metadata <- as.numeric(h_arr[[1]][3])
      num_sample_metadata <- as.numeric(h_arr[[1]][4])
    }, error = function(err) {
      stop(paste("Error parsing GCT parameters:", err, sep = " "))
    })
    
    expr <- read.table(file, skip = 2, header = TRUE, sep = "\t", as.is=T, quote = "\"", check.names=FALSE, stringsAsFactors=FALSE)
    # , na.strings="", stringsAsFactors=FALSE)
    
    sample_metadata <- t(expr[1:num_sample_metadata,(num_probe_metadata+2):ncol(expr)])
    probe_metadata <- expr[(num_sample_metadata+1):nrow(expr),1:num_probe_metadata]
    
    # check uniqueness probes
    checkProbeName <- table(probe_metadata[,1])
    if(max(checkProbeName) > 1) {
      stop(paste("Probes/genes in gct file should be unique: ", names(which.max(checkProbeName)), sep = " "))
    }
    rownames(probe_metadata) <- probe_metadata[,1]
    
    # check uniqueness samples
    checkSampleName <- table(sample_metadata[,1])
    if(max(checkSampleName) > 1) {
      stop(paste("Samples in gct file should be unique (samples=",num_sample_metadata,", probes=", num_probe_metadata,"): ", names(which.max(checkSampleName)), sep = " "))
    }
    # rownames(sample_metadata) <- sample_metadata[,1]
    colnames(sample_metadata) <- expr[1:num_sample_metadata,1]
    
    # remove metadata from table
    values <- expr[(num_sample_metadata+1):nrow(expr),(num_probe_metadata+2):ncol(expr)]
    rownames(values) <- rownames(probe_metadata)
    
    # convert to numeric
    values[] <- lapply(as.data.frame(values), function(x) { as.numeric(as.character(x)) })
    
    return(list(values=values, 
                sample_metadata=as.data.frame(sample_metadata), 
                probe_metadata=as.data.frame(probe_metadata)))
  }

shinyServer(function(input, output, session) {
  
  values <- reactiveValues(gctdata=NULL)
  
  # handle file upload and show sample metadata
  output$gct_sample_data <- renderDataTable({
    
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, it will be a data frame with 'name',
    # 'size', 'type', and 'datapath' columns. The 'datapath'
    # column will contain the local filenames where the data can
    # be found.
    
    inFile <- input$file1
    
    if (is.null(inFile))
      return(NULL)
    
    values$gctdata <- read.gct13(inFile$datapath)
    updateSelectInput(session, "variable", choices=colnames(values$gctdata$sample_metadata))
    updateSelectInput(session, "id", choices=colnames(values$gctdata$probe_metadata))
    
    values$gctdata$sample_metadata
    
  }, caption = "GCT Sample Metadata")
  
  # show probe metadata
  output$gct_probe_data <- renderDataTable({
    values$gctdata$probe_metadata
  }, caption = "GCT Probe Metadata")
  
  # handle UI group values update based on selected variable
  observe({
    if (input$variable !="") {
      updateSelectizeInput(session, 'group1', choices=unique(values$gctdata$sample_metadata[,input$variable]), server = TRUE)
      updateSelectizeInput(session, 'group2', choices=unique(values$gctdata$sample_metadata[,input$variable]), server = TRUE)
    }
  })
  
  # Create signature, upload to API, display results
  observeEvent(input$compute, {
    withProgress(message = 'Creating signature', value = 0, {
      
      #
      # Filter into two groups
      #
      group1 <- values$gctdata$sample_metadata[(values$gctdata$sample_metadata[,input$variable] %in% input$group1) & (!is.na(values$gctdata$sample_metadata[,input$variable])),]
      group2 <- values$gctdata$sample_metadata[(values$gctdata$sample_metadata[,input$variable] %in% input$group2) & (!is.na(values$gctdata$sample_metadata[,input$variable])),]
      
      #
      # ensure the gct values are numeric
      #
      values$gctdata$values[] <- lapply(values$gctdata$values, function(x) { as.numeric(as.character(x)) })
      
      #
      # select differential function
      #
      incProgress(1/3, detail = paste0("Running ",input$difffunction))
      if (input$difffunction=="t-test") {
         diff_result <- as.data.frame(apply(values$gctdata$values, 1, 
                                            function(x) t.test(unlist(x[rownames(group1)], use.names = FALSE),
                                                               unlist(x[rownames(group2)], use.names = FALSE))$p.value))

         #Pulls t-statistic as measure of diff-exp
         diff_result$Value_LogDiffExp <- apply(values$gctdata$values, 1, 
                                               function(x) t.test(unlist(x[rownames(group1)], use.names = FALSE),
                                                                  unlist(x[rownames(group2)], use.names = FALSE))$statistic)
      }
    
      #
      # format signature output
      #
      output_id_column_name <- paste0(input$id,"_",input$idtype)
      diff_result <- data.frame(values$gctdata$probe_metadata[input$id], diff_result[,1:2])
      colnames(diff_result) <- c(output_id_column_name, "Significance_pvalue", "Value_LogDiffExp")

  
      # User is allowed to choose if they want to filter based upon L1000 or not (input$lk)
      # They want to filter, if they don't we don't have to do anythign special.
      if (input$lk == TRUE)
      {
        # get the list of L1000 genes
        t <- read.csv2('L1000.txt', sep = '\t')
        l_genes <- t %>% filter(pr_is_lm == "1") %>% select(pr_gene_symbol) %>% as.vector()
        l_genes <- l_genes$pr_gene_symbol
        
        # now filter diff_results
        diff_result <- diff_result[diff_result[[output_id_column_name]] %in% l_genes,]
      }
      
      # We now will select just the top 100 if that's what the user wanted.
      # They only want top 100 (we do this by abs(diff_exp)
      
      # We want to filter
      if (input$top != 1)
      {
        # Need to translate the filter they want
        tmp = as.integer(input$top)
        filt = choices[tmp]
        
        # Add column for absolute Diff Exp to sort on
        diff_result$abso = lapply(diff_result$Value_LogDiffExp, abs)
        diff_result <- as.data.frame(lapply(diff_result, unlist))
        
        # Grab top 100 based on absolute
        diff_result <- arrange(diff_result, desc(abso))
        
        diff_result <- head(diff_result, filt)
        
        # Now drop that column
        diff_result$abso <- NULL
      }      
      
      
      #
      # show signature in a table
      #
      output$signature_data <- DT::renderDataTable({
        diff_result
      }, caption = "Signature to submit to iLincs")
      
      incProgress(1/3, detail = paste("Submitting the signature to iLincs"))
      
      #
      # create temporary csv file to submit into API
      #
      ftemp <- tempfile(pattern = "file", fileext=".csv", tmpdir = tempdir())
      write.csv(diff_result, ftemp, row.names = F, quote=F)
      cat(ftemp)
      
      r <- POST("http://www.ilincs.org/api/SignatureMeta/uploadAndAnalyze?lib=LIB_5", body = list(file = upload_file(ftemp)))
      
      l <- lapply(content(r)$status$data, function(x) unlist(x))
      ilincs_result <- data.frame(t(sapply(l,c)))
      
      # Need to fix similarity which comes back as factor
      # but need to convert to numeric so it can be sorted
      # on in DataTable by user
      ilincs_result$similarity <- as.numeric(as.character(ilincs_result$similarity))
      
      
      #
      # show correlation results
      #
      output$correlated_data <- DT::renderDataTable({
        datatable( ilincs_result, rownames = TRUE, caption = "Correlated signatures")
      

    })
  })
})
  
})  
