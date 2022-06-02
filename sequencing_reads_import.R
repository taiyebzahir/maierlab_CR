sequencing_reads_import <- function(filename){
  abundance <- read.table(filename)
  abundance <- abundance %>% select(1,3) %>% spread(df.V6, miseq_percent)
  return(abundance)
}


