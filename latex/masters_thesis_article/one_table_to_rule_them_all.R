library(dplyr)

pya_csv <- read.csv(
  '../data/pya_ukp5.csv', sep = ";", stringsAsFactors = FALSE
)
pya_csv$X <- NULL
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
# clean the cplex and gurobi timeouts, and the 
sum(is.na(reduced_pya$internal_time))
#t <- data.frame(diff = (reduced_pya$opt - reduced_pya$c_opt))
#t$r <- t$diff / reduced_pya$opt 
#t <- filter(t, diff != 0)
#filter(reduced_pya, opt == 0)
#reduced_pya <- reduced_pya %>% mutate(internal_time = case_when(
#  opt == 0 ~ as.numeric(NA),
#  opt != 0 ~ internal_time
#))

ifelse(reduced_pya$opt == 0, NA, reduced_pya$internal_time)

reduced_pya <- reduced_pya %>% mutate(internal_time = ifelse(
  (opt - c_opt)/opt < -1e-7, NA, internal_time
))
sum(is.na(reduced_pya$internal_time))
# GUROBI RETURNED THE WRONG ANSWER FOR NO 
reduced_pya <- filter(reduced_pya, internal_time < 1799)


breq_csv <- read.csv(
  '../data/128_16_std_breqd.csv', sep = ";", stringsAsFactors = FALSE
)
rr_csv <- read.csv(
  '../data/realistic_random.csv', sep = ";", stringsAsFactors = FALSE
)
csp_csv <- read.csv(
  '../data/csp_6195.csv', sep = ";", stringsAsFactors = FALSE
)

reduced_columns <- reduced_pya %>% group_by(algorithm) %>% summarise(
  mean = mean(internal_time)
)

