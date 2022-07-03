# delboy messages.
.db_message <- function(msg, color){
  msg <- paste(msg,"\n",sep="")
  switch(color,
         blue = cat(crayon::blue(msg)),
         red = cat(crayon::red(msg)),
         green = cat(crayon::green(msg)),
         magenta = cat(crayon::magenta(msg))
  )
}

# Check data inputs.
.check_data_inputs <- function(data, group_1, group_2, gene_column, grna_column = NULL,
                               bcorr_data_validation = NULL, batches = NULL){
  if(any(!group_1 %in% colnames(data)))
    stop(paste("unable to find the following group 1 columns:",setdiff(group_1,colnames(data))))
  if(any(!group_2 %in% colnames(data)))
    stop(paste("unable to find the following group 2 columns:",setdiff(group_2,colnames(data))))
  if(!gene_column %in% colnames(data))
    stop(paste("unable to find the gene column",gene_column))
  if(!is.null(grna_column)){
    if(!grna_column %in% colnames(data))
      stop(paste("unable to find the gRNA column",grna_column))
  }
  if(!is.null(bcorr_data_validation)){
    if(!is.data.frame(bcorr_data_validation))
      stop(paste("'bcorr_data_validation' should be a data frame, instead got",
                 class(bcorr_data_validation)))
  }
  if(!is.null(batches)){
    if(length(batches) != length(group_1) + length(group_2))
      stop(paste("the number of batched samples must equal the number of total samples across the two groups:\n",
                 "expected",length(group_1) + length(group_2),"but got",length(batches)))
    if(is.null(names(batches))) 
      stop("batches must be a named character vector with names referring to sample columns in the input data")
    if(length(setdiff(names(batches),colnames(data))))
      stop(paste("batch samples missing in data input:",setdiff(names(batches),colnames(data))))
  }
}
