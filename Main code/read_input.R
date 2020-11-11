read_input <- function() {
  train_n <- readline(prompt = "Enter the training split:")
  test_n <- readline(prompt = "Enter the testing split:")
  val_n <- readline(prompt = "Enter the validation split:")
  Fs <- readline(prompt = "Enter the sampling frequency:")
  alpha <- readline(prompt = "What is the value of Pfa?")
  kern <- readline(prompt = "What kernel would you like to use?")
  input_list <- list( 
                     "train" = as.numeric(train_n), 
                     "test" = as.numeric(test_n), 
                     "validation" = as.numeric(val_n),
                     "freq" = as.integer(Fs),
                     "alpha" = as.numeric(alpha),
                     "kern" = kern)
  return(input_list)
}
