library(statsr)
library(dplyr)
library(ggplot2)

mtu <- read.csv("mtu_impl_desktop_uc5.csv", sep = ";")
# Change the NAs on internal_time to 999 seconds (the timeout)
mtu$internal_time <- sapply(mtu$internal_time, function(x) if (is.na(x)) 999 else x)
# Remove empty column created by extra ';'
mtu$X <- NULL

mtu1 <- filter(mtu, algorithm == "cpp-mtu1_desktop_uc5" | algorithm == "fmtu1_desktop_uc5") # only mtu1
mtu2 <- filter(mtu, algorithm == "cpp-mtu2_desktop_uc5" | algorithm == "fmtu2_desktop_uc5") # only mtu2

ukp_time_comp_plot <- function (data) {
  data <- select(data, algorithm, filename, internal_time)
  data$abbr_alg <- sapply(data$algorithm, (function (f) { factor(strsplit(as.character(f), "_")[[1]][1])}))
  data <- data %>% group_by(filename) %>% mutate(fname_mean_time = mean(internal_time)) %>% arrange(fname_mean_time)
  data$filename <- factor(data$filename, levels = unique(data$filename))
  p <- ggplot(data, aes(x = as.numeric(filename), y = internal_time, color = abbr_alg)) 
  p <- p + geom_point() + scale_y_log10()
  p <- p + xlab("Instances ordered by the time taken to solve")
  p <- p + ylab("Time on seconds (log10 scale)")
}

p1 <- ukp_time_comp_plot(mtu1) + ggtitle("MTU1 (Fortran vs C++) over 454 instances")
ggsave(filename = "mtu1_comparison_y_log.pdf", plot = p1)
p2 <- ukp_time_comp_plot(mtu2) + ggtitle("MTU2 (Fortran vs C++) over 454 instances")
ggsave(filename = "mtu2_comparison_y_log.pdf", plot = p2)

