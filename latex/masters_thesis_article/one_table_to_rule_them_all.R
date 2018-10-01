library(dplyr)

pya_csv <- read.csv(
  '../data/pya_ukp5.csv', sep = ";", stringsAsFactors = FALSE
)
pya_csv$X <- NULL
pya_columns <- pya_csv %>% group_by(algorithm) %>% 
  filter(algorithm != 'cplex2nt' & algorithm != 'gurobi0nt') %>%
  summarise(
    finished = 4540 - sum(is.na(internal_time)),
    mean_time = mean(internal_time, na.rm = T)
  )

mtu_csv <- read.csv(
  '../data/mtus_pya.csv', sep = ";", stringsAsFactors = FALSE
) %>% filter(algorithm != 'fmtu1' & algorithm != 'fmtu2')
mtu_csv$X <- NULL
# get only the lines of pya_csv related to reduced dataset instances and
# add the lines of mtu_csv to them
reduced_pya <- rbind(semi_join(pya_csv, mtu_csv, by = 'filename'), mtu_csv)
# creates a column correct opt for the whole table
correct_opt <- filter(reduced_pya, algorithm == 'ordered_step_off') %>%
  select(filename, c_opt = opt)
reduced_pya <- inner_join(reduced_pya, correct_opt, by = 'filename')

reduced_pya <- reduced_pya %>% mutate(
  opt = ifelse(is.na(opt) | opt == 0, NA, opt)
) %>% mutate(
  # Some Gurobi Instances had rounding errors and finished with the
  # optimal solution but wrong opt (error of one unity in values greater
  # than 10^6), other had larger errors and wrong solutions, this is the
  # cuttoff to excluding only instances with wrong solutions.
  opt = ifelse(is.na(opt) | abs(opt - c_opt)/opt > 1e-7, NA, opt) 
) %>% mutate(
  internal_time = ifelse(is.na(opt) | internal_time > 1799, NA, internal_time)
) %>% filter(
  algorithm == 'cpp-mtu1' |
  algorithm == 'cpp-mtu2' |
  algorithm == 'cplex2nt' |
  algorithm == 'gurobi0nt'
)

reduced_columns <- reduced_pya %>% group_by(algorithm) %>% summarise(
  finished = 454 - sum(is.na(internal_time)),
  mean_time = mean(internal_time, na.rm = T)
)
both_pya <- rbind(pya_columns, reduced_columns)

rr_csv <- read.csv(
  '../data/realistic_random.csv', sep = ";", stringsAsFactors = FALSE
)
rr_csv$X <- NULL
rr_columns <- rr_csv %>% group_by(algorithm) %>% summarise(
  finished = 80 - sum(is.na(internal_time)),
  mean_time = mean(internal_time, na.rm = T)
)

breq_csv <- read.csv(
  '../data/128_16_std_breqd.csv', sep = ";", stringsAsFactors = FALSE
)
breq_csv$X <- NULL
breq_columns <- breq_csv %>% group_by(algorithm) %>% summarise(
  finished = 80 - sum(is.na(internal_time)),
  mean_time = mean(internal_time, na.rm = T)
)

csp_csv <- read.csv(
  '../data/csp_6195.csv', sep = ";", stringsAsFactors = FALSE
)
csp_csv$X <- NULL
csp_csv <- csp_csv %>% mutate(
  internal_time = hex_sum_knapsack_time + hex_sum_sort_time
)
csp_columns <- csp_csv %>% group_by(algorithm) %>% summarise(
  finished = 6195 - sum(is.na(internal_time)),
  mean_time = mean(internal_time, na.rm = T)
)

inner_join(both_pya, rr_columns, by = 'algorithm')
