library(statsr)
library(dplyr)
library(ggplot2)
library(stringi)

env <- read.csv("env_influence.csv", sep = ";")
env$X <- NULL

# Checking the standard deviation difference between serial and parallel
# on the same computer, with the same algorithm.
alg_name = "ukp5"
com_name = "notebook"
serial_big_sd <- filter(env, algorithm == alg_name & computer == com_name & mode == "serial") %>% group_by(filename) %>% summarise(internal_time_sd = sd(internal_time)) %>% arrange(filename)
parallel_big_sd <- filter(env, algorithm == alg_name & computer == com_name & mode == "parallel") %>% group_by(filename) %>% summarise(internal_time_sd = sd(internal_time)) %>% arrange(filename)
com_sd <- data.frame(filename = serial_big_sd$filename, serial_sd = serial_big_sd$internal_time_sd, parallel_sd = parallel_big_sd$internal_time_sd)

com_sd2 <- com_sd %>% mutate(serial_parallel_sd_ratio = serial_sd / parallel_sd)# %>% mutate(inst_class = strsplit(as.character(filename), "_")[[1]][1])
com_sd2$inst_class <- sapply(com_sd2$filename, (function (f) { factor(strsplit(as.character(f), "_")[[1]][1])}))

p <- ggplot(data = com_sd2, aes(x = serial_sd, y = parallel_sd)) + ggtitle("SD Time (seconds, fixed scale)") + geom_point() + coord_fixed()
ggsave('sd_time_comp_parallel_serial.pdf', p)

# Table with the mean internal_time of all runs. This is made grouping the runs
# by algorithm, computer, mode and filename, removing the columns run_number
# (don't make sense anymore), opt, external_memory and external_time; adding
# inst_class and computer_mode; the column internal_time now has the mean of all
# similar runs (the ones with the same algorithm, computer, mode and file).
# TODO: Consider filtering by algorithm before generating chart and generating
# two charts (one for each algorithm).
env_file_mean_times <- env %>% select(algorithm, computer, mode, filename, internal_time, external_memory) %>% group_by(algorithm, computer, mode, filename) %>% summarise(internal_time = mean(internal_time), external_memory = mean(external_memory))
env_file_mean_times$inst_class <- sapply(env_file_mean_times$filename, (function (f) { factor(strsplit(as.character(f), "_")[[1]][1])}))
env_file_mean_times <- env_file_mean_times %>% mutate(computer_mode = factor(paste(sep = "_", as.character(algorithm), as.character(computer), as.character(mode)))) # TODO: reorder factor using the time maximums, so the colors will be in order in the figure
env_file_mean_times <- env_file_mean_times %>% group_by(filename) %>% mutate(mean_time = mean(internal_time)) %>% arrange(mean_time)
env_file_mean_times$filename <- factor(env_file_mean_times$filename, levels = unique(env_file_mean_times$filename))

# With the files (y axis) ordered in increasing order of ratio
# This should starts at ~1, stay at 1~1.5 for the most files, and end at 2~3
alg_name = "ukp5"
com_name = "notebook"
serial <- filter(env_file_mean_times, algorithm == alg_name & computer == com_name & mode == "serial") %>% arrange(filename) %>% select(filename,internal_time,external_memory,inst_class)
parallel <- filter(env_file_mean_times, algorithm == alg_name & computer == com_name & mode == "parallel") %>% arrange(filename) %>% select(filename,internal_time,external_memory)
ratio_parallel_serial <- data.frame(filename = serial$filename, inst_class = serial$inst_class, serial_time = serial$internal_time, parallel_time = parallel$internal_time, mem_use = serial$external_memory) %>% mutate(ratio = parallel_time/serial_time) %>% arrange(ratio)
ratio_parallel_serial$filename <- factor(ratio_parallel_serial$filename, levels = unique(ratio_parallel_serial$filename))
p <- ggplot(ratio_parallel_serial, aes(x = mem_use/1024, y = ratio, shape = inst_class)) + geom_point()
ggsave('ratio_parallel_serial_by_mem_use.pdf', p)
p <- ggplot(ratio_parallel_serial, aes(x = serial_time, y = ratio, shape = inst_class)) + geom_point()
ggsave('ratio_parallel_serial_by_serial_time.pdf', p)
p <- ggplot(ratio_parallel_serial, aes(x = serial_time, y = parallel_time, shape = inst_class)) + geom_point()
ggsave('parallel_by_serial_time.pdf', p)

