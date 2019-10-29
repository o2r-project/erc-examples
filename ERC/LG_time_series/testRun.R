library("plumber")
path = paste("test.R", sep = "")
r <- plumb(path)
r$run(host = "0.0.0.0", port=8123)