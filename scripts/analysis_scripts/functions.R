read_file <- function(x) {
  f <- paste(here("data", "input_data", "plot_coordinates"), x, sep = "/")
  g <- read.csv(f)
  return(g)
}

ddm_to_dd <- function(dms) {
  
  # split strings to separate degrees and minutes.m
  split1 <- str_split(as.vector(dms), pattern = " ")
  
  # extract degrees
  degrees <- substr(split1[[1]][2], 1, 2)
  degrees <- as.double(degrees)
  
  # extract minutes.m
  minutes.m <- str_replace(split1[[1]][3], pattern = "'", replace = "")
  minutes.m <- as.double(minutes.m)
  
  # convert to degrees
  dot_d <- minutes.m/60
  
  # calculate decimal degrees
  decimal_degrees = degrees + dot_d
  
  return(decimal_degrees)
}