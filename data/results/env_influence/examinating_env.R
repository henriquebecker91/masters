env <- read.csv("env_influence.csv", sep = ";")

library(statsr)
library(dplyr)
library(ggplot2)

# Cheking the biggest standard deviations between runs by hand.
filter(env, algorithm == "ukp5" & computer == "notebook" & mode == "serial") %>% group_by(filename) %>% summarise(internal_time_sd = sd(internal_time)) %>% arrange(desc(internal_time_sd))
filter(env, algorithm == "ukp5" & computer == "notebook" & mode == "parallel") %>% group_by(filename) %>% summarise(internal_time_sd = sd(internal_time)) %>% arrange(desc(internal_time_sd))

# Biggest and smallest SD between same file on the notebook, serial, ukp5
filter(env, algorithm == "ukp5" & computer == "notebook" & mode == "serial" & filename == "saw_n100000wmin10000-18-s152335818c5364756.ukp")
filter(env, algorithm == "ukp5" & computer == "notebook" & mode == "serial" & filename == "ss2_wmin1000wmax500000n50000-0-s827183242c5970156.ukp")

# Biggest and smallest SD between same file on the notebook, parallel, ukp5
filter(env, algorithm == "ukp5" & computer == "notebook" & mode == "parallel" & filename == "nsds2_n50000wmin20000-19-s618866736c1469822.ukp")
filter(env, algorithm == "ukp5" & computer == "notebook" & mode == "parallel" & filename == "sc_a5n5000wmin10000-0-c591952.ukp")

# Checking the standard deviation difference between serial and parallel
# on the same computer, with the same algorithm.
alg_name = "pyasukpt"
com_name = "notebook"
serial_big_sd <- filter(env, algorithm == alg_name & computer == com_name & mode == "serial") %>% group_by(filename) %>% summarise(internal_time_sd = sd(internal_time)) %>% arrange(filename)
parallel_big_sd <- filter(env, algorithm == alg_name & computer == com_name & mode == "parallel") %>% group_by(filename) %>% summarise(internal_time_sd = sd(internal_time)) %>% arrange(filename)
com_sd <- data.frame(filename = serial_big_sd$filename, serial_sd = serial_big_sd$internal_time_sd, parallel_sd = parallel_big_sd$internal_time_sd)

com_sd2 <- com_sd %>% mutate(serial_parallel_sd_ratio = serial_sd / parallel_sd)# %>% mutate(inst_class = strsplit(as.character(filename), "_")[[1]][1])
com_sd2$inst_class <- sapply(com_sd2$filename, (function (f) { factor(strsplit(as.character(f), "_")[[1]][1])}))

ggplot(data = com_sd2, aes(x = serial_sd, y = parallel_sd, color = inst_class)) + ggtitle("SD Time (fixed scale)") + geom_point() + coord_fixed()
ggplot(data = com_sd2, aes(x = serial_sd, y = parallel_sd, color = inst_class)) + ggtitle("SD Time (non-fixed scale)") + geom_point()
ggplot(data = com_sd2, aes(x = serial_sd, y = parallel_sd, color = inst_class)) + ggtitle("SD Time (log10(y) scale)") + geom_point() + scale_y_log10()

# Table with the mean internal_time of all runs. This is made grouping the runs
# by algorithm, computer, mode and filename, removing the columns run_number
# (don't make sense anymore), opt, external_memory and external_time; adding
# inst_class and computer_mode; the column internal_time now has the mean of all
# similar runs (the ones with the same algorithm, computer, mode and file).
# TODO: Consider filtering by algorithm before generating chart and generating
# two charts (one for each algorithm).
env_file_mean_times <- env %>% select(algorithm, computer, mode, filename, internal_time) %>% group_by(algorithm, computer, mode, filename) %>% summarise(internal_time = mean(internal_time))
env_file_mean_times$inst_class <- sapply(env_file_mean_times$filename, (function (f) { factor(strsplit(as.character(f), "_")[[1]][1])}))
env_file_mean_times <- env_file_mean_times %>% mutate(computer_mode = factor(paste(sep = "_", as.character(algorithm), as.character(computer), as.character(mode)))) # TODO: reorder factor using the time maximums, so the colors will be in order in the figure
env_file_mean_times <- env_file_mean_times %>% group_by(filename) %>% mutate(mean_time = mean(internal_time)) %>% arrange(mean_time)
env_file_mean_times$filename <- factor(env_file_mean_times$filename, levels = unique(env_file_mean_times$filename))
ggplot(env_file_mean_times, aes(x = as.numeric(filename), y = internal_time, color = computer_mode, shape = inst_class)) + ggtitle("UKP5 Time Variation By Env") + geom_point()
ggplot(env_file_mean_times, aes(x = as.numeric(filename), y = internal_time, color = computer_mode, shape = inst_class)) + ggtitle("UKP5 Time Variation By Env") + geom_point() + scale_y_log10()
ggplot(env_file_mean_times, aes(x = as.numeric(filename), y = internal_time, color = computer_mode)) + ggtitle("UKP5 Time Variation By Env") + geom_line()
ggplot(env_file_mean_times, aes(x = as.numeric(filename), y = internal_time, color = computer_mode)) + ggtitle("UKP5 Time Variation By Env") + geom_line() + scale_y_log10()

# Ratio of the mean execution time on parallel by the mean execution time on serial
# With the files (y axis) ordered in increasing order of ratio
# This should starts at ~1, stay at 1~1.5 for the most files, and end at 2~3
alg_name = "pyasukpt"
com_name = "notebook"
serial <- filter(env_file_mean_times, algorithm == alg_name & computer == com_name & mode == "serial") %>% arrange(filename) %>% select(filename,internal_time,inst_class)
parallel <- filter(env_file_mean_times, algorithm == alg_name & computer == com_name & mode == "parallel") %>% arrange(filename) %>% select(filename,internal_time)
ratio_parallel_serial <- data.frame(filename = serial$filename, inst_class = serial$inst_class, serial_time = serial$internal_time, parallel_time = parallel$internal_time) %>% mutate(ratio = parallel_time/serial_time) %>% arrange(ratio)
ratio_parallel_serial$filename <- factor(ratio_parallel_serial$filename, levels = unique(ratio_parallel_serial$filename))
ggplot(ratio_parallel_serial, aes(x = as.numeric(filename), y = ratio, shape = inst_class)) + ggtitle("Some Algorithm Ratio Between Parallel and Serial in the same machine") + geom_point()
