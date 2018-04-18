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

