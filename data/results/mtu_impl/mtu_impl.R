library(statsr)
library(dplyr)
library(ggplot2)

mtu <- read.csv("mtu_impl_desktop_uc5.csv", sep = ";")
#sapply(mtu, function(x) sum(is.na(x))) # used to check number of NAs
# Change the NAs on internal_time to 999 seconds (the timeout)
mtu$internal_time <- sapply(mtu$internal_time, function(x) if (is.na(x)) 999 else x)
# Remove empty column created by extra ';'
mtu$X <- NULL

mtu1 <- filter(mtu, algorithm == "cpp-mtu1_desktop_uc5" | algorithm == "fmtu1_desktop_uc5") # only mtu1
mtu2 <- filter(mtu, algorithm == "cpp-mtu2_desktop_uc5" | algorithm == "fmtu2_desktop_uc5") # only mtu2
cppmtu <- filter(mtu, algorithm == "cpp-mtu1_desktop_uc5" | algorithm == "cpp-mtu2_desktop_uc5") # only c++ code
fmtu <- filter(mtu, algorithm == "fmtu1_desktop_uc5" | algorithm == "fmtu2_desktop_uc5") # only fortran

algs <- c("cpp-mtu1_desktop_uc5", "cpp-mtu2_desktop_uc5", "fmtu1_desktop_uc5", "fmtu2_desktop_uc5")
mtus <- lapply(algs, function(alg) filter(mtu, algorithm == alg))
names(mtus) <- algs


ukp_time_comp_plot <- function (data) {
  data <- select(data, algorithm, filename, internal_time)
  data$inst_class <- sapply(data$filename, (function (f) { factor(strsplit(as.character(f), "_")[[1]][1])}))
  data <- data %>% group_by(filename) %>% mutate(fname_mean_time = mean(internal_time)) %>% arrange(fname_mean_time)
  data$filename <- factor(data$filename, levels = unique(data$filename))
  ggplot(data, aes(x = as.numeric(filename), y = internal_time, color = algorithm, shape = inst_class)) 
}

# Useful information about the data 
sapply(mtus, function(x) mean(x$internal_time))
summary(mtu)
ukp_time_comp_plot(mtu1) + geom_point() + scale_y_log10()
