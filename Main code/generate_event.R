generate_event <- function (event_filename){
  library(R.matlab)
  event_func <- readMat(event_filename)
  peak_data <- event_func[[3]]
  ann_data <- event_func[[1]]
  peak_list <- list("peakval" = peak_data, "annotations" = ann_data)
  return(peak_list)
}