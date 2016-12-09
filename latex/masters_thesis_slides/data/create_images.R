library(ggplot2)
library(statsr)
library(dplyr)
library(stringi)
require(graphics)

for (name in c('hi_n5000-0-s731232778c5052835',
              'nsds2_n50000wmin50000-0-s155213243c1070927',
              'saw_n10000wmin10000-0-s985850175c608451',
              'sc_a5n10000wmin10000-0-c6020768',
              'ss2_wmin1000wmax1000000n10000-0-s827183242c7523725',
              '128_16_std_breqd-n16384-s0')) {
  csv <- paste(name, '.csv', sep = '')
  image <- paste(name, '.pdf', sep = '')
  
  t <- read.csv(csv, sep = ';')
  ggplot(t, aes(x = w, y = p)) +
    geom_point() +
    xlab('weight') + 
    ylab('profit')
  ggsave(image, units = 'cm', width = 10, height = 7.5)
}

csv_dir <- '../../masters_thesis/data/'
add_base <- function (csv_name) { paste(csv_dir, csv_name, sep = '') }
csp_csv <- read.csv(add_base('cutstock_knap_solvers.csv'), sep = ";")
csp_csv$X <- NULL
# forgot to add sort_time to knapsack time in the chart generation codes, now
# changing it here to avoid changing many places
csp_csv <- mutate(csp_csv, hex_sum_knapsack_time = hex_sum_knapsack_time + hex_sum_sort_time)
breq_csv <- read.csv(add_base("128_16_std_breqd_all.csv"), sep = ";")
breq_csv $X <- NULL
fast <- read.csv(add_base("pya_ukp5.csv"), sep = ";")
fast$X <- NULL
slow <- read.csv(add_base("mtus.csv"), sep = ";")
slow$X <- NULL
mtu <- read.csv(add_base("mtu_impl_desktop_uc1.csv"), sep = ";")
mtu$X <- NULL
env_data <- read.csv(add_base("env_influence.csv"), sep = ";")
env_data$X <- NULL

compare_num_iter <- function(csv, first_m, second_m) {
  first_t <- csv %>% select(filename, algorithm, total_iter) %>%
    filter(algorithm == first_m) %>% arrange(filename)
  second_t <- csv %>% select(filename, algorithm, total_iter) %>%
    filter(algorithm == second_m) %>% arrange(filename)
  first_t$diff <- first_t$total_iter - second_t$total_iter
  first_t$norm_diff <- first_t$diff / second_t$total_iter
  first_t <- first_t %>% arrange(norm_diff)
  first_t$order <- 1:length(first_t$norm_diff)
  first_t
}

ukp_time_comp_plot <- function (data, legend) {
  data <- select(data, algorithm, filename, internal_time)
  data$algorithm <- sapply(data$algorithm, (function (f) { factor(strsplit(as.character(f), "_")[[1]][1])}))
  data$language <- data$algorithm
  data$language <- sapply(data$language, (function (f) { gsub("cpp-mtu[12]", "C++", f) }))
  data$language <- sapply(data$language, (function (f) { gsub("fmtu[12]", "Fortran", f) }))
  data$language <- factor(data$language, levels = c('C++', 'Fortran'))
  data$algorithm <- sapply(data$algorithm, (function (f) { gsub(".*1", "MTU1", f) }))
  data$algorithm <- sapply(data$algorithm, (function (f) { gsub(".*2", "MTU2", f) }))
  data$algorithm <- factor(data$algorithm, levels = c('MTU1', 'MTU2'))
  data <- data %>% group_by(filename) %>% mutate(fname_mean_time = mean(internal_time)) %>% arrange(fname_mean_time)
  data$filename <- factor(data$filename, levels = unique(data$filename))
  p <- ggplot(data, aes(x = as.numeric(filename), y = internal_time, color = language))
  p <- p + geom_point() + scale_y_continuous(
    trans = 'log10', limits = c(0.001, 1000),
    breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000),
    labels = c('0.001', '0.01', '0.1', '1', '10', '100', '1000')
  )
  p <- p + xlab("Instances ordered by mean time to solve")
  p <- p + ylab("Seconds to solve")
  p <- p + theme(legend.position = legend)
}

fix_names <- function(x) {
  x <- gsub('all', 'ALL', x)
  x <- gsub('ss2', 'S.S.', x)
  x <- gsub('saw', 'SAW', x)
  x <- gsub('nsds2', 'P.P.', x)
  x <- gsub('sc', 'S.C.', x)
  x <- gsub('hi', 'W.C.D.', x)
  x
}

fix_names2 <- function(x) {
  x <- gsub('_cutstock', '', x)
  x
}

csv_no_na <- fast[complete.cases(fast), ]
csv_no_na$type <- sapply(csv_no_na$filename, (function (f) {
  factor(strsplit(as.character(f), "_")[[1]][1],
         levels = c('all', 'hi', 'nsds2', 'saw', 'sc', 'ss2'))
}))
csv_no_na_mean <- csv_no_na %>%
  group_by(filename) %>%
  mutate(mean_methods_time = mean(internal_time))
#csv_no_na$n <- sapply(csv_no_na$filename, (function (f) { as.numeric(gsub("n", "", strsplit(as.character(f), "-")[[1]][2]))}))
dt_copy <- csv_no_na_mean
dt_copy$type <- factor('all', levels = c('all', 'hi', 'nsds2', 'saw', 'sc', 'ss2'))
csv_no_na_order <- rbind(csv_no_na_mean, dt_copy) %>%
  arrange(mean_methods_time)
csv_no_na_order$filename <- factor(csv_no_na_order$filename,
                                   levels = unique(csv_no_na_order$filename))

csv_no_na_order$type <- fix_names(csv_no_na_order$type)

red <- '#F8766D'
blue <- '#619CFF'
green <- '#00BA38'
myplot <- function(t, sel_type) {
  t2 <- filter(t, type == sel_type)
  ggplot(t2,            
         aes(x = as.numeric(filename),
             y = internal_time,
             color = algorithm)) +
    geom_point() +
    scale_y_continuous(trans = 'log10', limits = c(0.001, 1000),
                       breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000),
                       labels = c('0.001', '0.01', '0.1', '1', '10',
                                  '100', '1000')) +
    ylab('Seconds to solve') +
    xlab('Instance index when ordered by average time') +
    theme(legend.position = 'bottom', text = element_text(size=8)) +
    scale_color_manual(values = c('mgreendp' = red, 'pyasukp' = green, 'ukp5' = blue))
}

myplot(csv_no_na_order, 'ALL')
ggsave('ALL.pdf', units = 'cm', width = 10, height = 7.5)
myplot(csv_no_na_order, 'S.S.')
ggsave('SS.pdf', units = 'cm', width = 10, height = 7.5)
myplot(csv_no_na_order, 'SAW')
ggsave('SAW.pdf', units = 'cm', width = 10, height = 7.5)
myplot(csv_no_na_order, 'P.P.')
ggsave('PP.pdf', units = 'cm', width = 10, height = 7.5)
myplot(csv_no_na_order, 'S.C.')
ggsave('SC.pdf', units = 'cm', width = 10, height = 7.5)
myplot(csv_no_na_order, 'W.C.D.')
ggsave('WCD.pdf', units = 'cm', width = 10, height = 7.5)

csv_no_na <- breq_csv
#csv_no_na[is.na(csv_no_na$internal_time),]$internal_time <- 1000
csv_no_na$n <- sapply(csv_no_na$filename,
                      (function (f) {
                        as.numeric(gsub("n", "", strsplit(as.character(f), "-")[[1]][2]))
                      }))

csv_no_na <- breq_csv
#csv_no_na[is.na(csv_no_na$internal_time),]$internal_time <- 1000
csv_no_na$n <- sapply(csv_no_na$filename,
                      (function (f) {
                        as.numeric(gsub("n", "", strsplit(as.character(f), "-")[[1]][2]))
                      }))

mtu$internal_time <- sapply(mtu$internal_time, function(x) if (is.na(x)) 1000 else x)

mtu1 <- filter(mtu, algorithm == "cpp-mtu1_desktop_uc1" | algorithm == "fmtu1_desktop_uc1") # only mtu1
mtu2 <- filter(mtu, algorithm == "cpp-mtu2_desktop_uc1" | algorithm == "fmtu2_desktop_uc1") # only mtu2

ukp_time_comp_plot(mtu1, 'bottom')
ggsave('mtu1.pdf', units = 'cm', width = 10, height = 7.5)
ukp_time_comp_plot(mtu2, 'bottom')
ggsave('mtu2.pdf', units = 'cm', width = 10, height = 7.5)

ggplot(csv_no_na,
       aes(x = n,#* (1 + (as.numeric(algorithm) - 1)/10),
           y = internal_time,
           color = algorithm)) +
  geom_point() +
  scale_y_continuous(trans = "log10",
                     breaks = c(1000, 100, 10, 1, 0.1, 0.01, 0.001),
                     labels = c('1000', '100', '10', '1', '0.1', '0.01', '0.001'),
                     limits = c(0.001, 1000)) +
  scale_x_continuous(trans = "log2",
                     breaks = c(2^11, 2^12, 2^13, 2^14, 2^15, 2^16,
                                2^17, 2^18, 2^19, 2^20),
                     labels = c(bquote(2^11), bquote(2^12), bquote(2^13), bquote(2^14), bquote(2^15), bquote(2^16),
                                bquote(2^17), bquote(2^18), bquote(2^19), bquote(2^20))) +
  ylab('Seconds to solve') +
  xlab('Value of n (instance size)') +
  theme(legend.position = 'bottom', text = element_text(size=8))

ggsave('breq_exp.pdf', units = 'cm', width = 10, height = 7.5)

csv_no_na <- csp_csv
csv_no_na$algorithm <- sapply(csv_no_na$algorithm, fix_names2)
csv_no_na[is.na(csp_csv$hex_sum_knapsack_time), ]$hex_sum_knapsack_time <- 600
#csv_no_na <- csv[!is.na(csv$hex_sum_knapsack_time), ]
csv_no_na_order <- csv_no_na %>% group_by(filename) %>% mutate(mean_methods_time = mean(hex_sum_knapsack_time)) %>% arrange(mean_methods_time)
csv_no_na_order$filename <- factor(csv_no_na_order$filename, levels = unique(csv_no_na_order$filename))
ggplot(csv_no_na_order,
       aes(x = as.numeric(filename),
           y = hex_sum_knapsack_time,
           color = algorithm)) +
  xlab('Instance index when ordered by mean time') +
  ylab('Time spent in pricing (in seconds)') +
  geom_point() +
  scale_y_continuous(
    trans = 'log10', limits = c(0.0001, 1000),
    breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000),
    labels = c('0.001', '0.01', '0.1', '1', '10', '100', '1000')
  ) +
  theme(legend.position = 'right', text = element_text(size=8))
ggsave('pricing_time.pdf', units = 'cm', width = 10, height = 7.5)

corr_time <- csp_csv
corr_time$algorithm <- sapply(csv_no_na$algorithm, fix_names2)
corr_time$relative_time <- corr_time$hex_sum_knapsack_time / (corr_time$hex_sum_master_prob_time + corr_time$hex_sum_knapsack_time)
corr_time$tt <- corr_time$hex_sum_knapsack_time + corr_time$hex_sum_master_prob_time
corr_time[is.na(corr_time$tt), ]$tt <- 600
corr_time <- corr_time %>% group_by(filename) %>% mutate(mean_methods_time = mean(tt)) %>% arrange(mean_methods_time)
corr_time$filename <- factor(corr_time$filename, levels = unique(corr_time$filename))
ggplot(corr_time,
       aes(x = as.numeric(filename),
           y = relative_time * 100,
           color = algorithm)) +
  xlab('Instance index when ordered by the mean time') +
  ylab('Time spent in pricing (in % of total time)') +
  geom_point() +
  theme(legend.position = 'right', text = element_text(size=8))
ggsave('pricing_percentage.pdf', units = 'cm', width = 10, height = 7.5)

env_file_mean_times <- env_data %>%
  select(algorithm, computer, mode, filename, internal_time) %>%
  group_by(algorithm, computer, mode, filename) %>%
  summarise(mean_time = mean(internal_time))

create_facet <- function(env_file_mean_times, alg_name, com_name) {
  serial <- filter(env_file_mean_times, mode == "serial" &
                     algorithm == alg_name &
                     computer == com_name) %>% arrange(filename)
  parallel <- filter(env_file_mean_times, mode == "parallel" &
                       algorithm == alg_name &
                       computer == com_name) %>% arrange(filename)
  ratio_parallel_serial <- data.frame(filename = serial$filename,
                                      algorithm = serial$algorithm,
                                      computer = serial$computer,
                                      serial_time = serial$mean_time,
                                      parallel_time = parallel$mean_time) %>%
    mutate(ratio = parallel_time/serial_time)
  ratio_parallel_serial
}

t1 <- create_facet(env_file_mean_times, "ukp5", "notebook")
t2 <- create_facet(env_file_mean_times, "ukp5", "desktop")
t3 <- create_facet(env_file_mean_times, "pyasukpt", "notebook")
t4 <- create_facet(env_file_mean_times, "pyasukpt", "desktop")
t <- rbind(t1, t2, t3, t4)
t$inst_class <- sapply(t$filename, (function (f) { factor(strsplit(as.character(f), "_")[[1]][1])}))
t$env <- paste(t$algorithm, t$computer, sep = '_')
t$inst_class <- fix_names(t$inst_class)
ggplot(t, aes(x = serial_time, y = ratio)) + #, color = inst_class)) +
  xlab('Serial runs times (seconds, log10)') +
  ylab('parallel time / serial time') +
  geom_point() + facet_wrap(~ env) +
  theme(text = element_text(size=8)) + 
  #theme(legend.position = 'bottom', text = element_text(size=8)) + theme(legend.title = element_blank()) +
  scale_x_continuous(trans = 'log10', breaks = c(0.1, 10, 100))
ggsave('env_influence.pdf', units = 'cm', width = 10, height = 7.5)

