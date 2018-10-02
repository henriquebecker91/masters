library(dplyr)

pya_csv <- read.csv(
  '../data/pya_ukp5.csv', sep = ";", stringsAsFactors = FALSE
) %>% filter(algorithm != 'fmtu1' & algorithm != 'fmtu2')
pya_csv$X <- NULL

# Get only the algorithms executed over the reduced pyasukp dataset.
reduced_pya <- read.csv(
  '../data/mtus_pya.csv', sep = ";", stringsAsFactors = FALSE
) %>% filter(algorithm != 'fmtu1' & algorithm != 'fmtu2')
reduced_pya$X <- NULL
reduced_pya <- rbind(
  reduced_pya,
  filter(pya_csv, algorithm == 'cplex2nt' | algorithm == 'gurobi0nt')
)
# We need to get the runs of OSO for the reduced pyasukp dataset to
# have one algorithm with the correct instance values for all instances.
reduced_pya <- rbind(
  reduced_pya,
  semi_join(
    filter(pya_csv, algorithm == 'ordered_step_off'),
    reduced_pya,
    by = 'filename'
  )
)

clean_and_summarise <- function(dt, correct_alg = 'none') {
  qt_inst <- length(unique(dt$filename))
  has_opt <- correct_alg != 'none'
  
  # Creates a column with the correct opt.
  if (has_opt) {
    correct_opt <- dt %>% filter(algorithm == correct_alg) %>%
      select(filename, c_opt = opt)
    dt <- inner_join(dt, correct_opt, by = 'filename')
  }
  # Set internal_time to NA for instances in which Gurobi finished because
  # memory exhaustion (opt == 0), or instances for which Gurobi had wrong
  # results (abs(opt - c_opt)/opt > 1e-7), or instances for which
  # cplex or gurobi hit the timeout (internal_time > 1799).
  if (has_opt) {
    dt <- dt %>% mutate(
      internal_time = ifelse(
        (is.na(opt)                    |
           opt == 0                    |
           abs(opt - c_opt)/opt > 1e-7 |
           internal_time > 1799),
        NA,
        internal_time
      )
    )
  } else {
    dt <- dt %>% mutate(
      internal_time = ifelse(internal_time > 1799, NA, internal_time)
    )
  }
  
  dt <- dt %>% select(algorithm, internal_time)
  
  cols <- dt %>% group_by(algorithm) %>% summarise(
    finished = sprintf('%d', qt_inst - sum(is.na(internal_time))),
    mean_time = sprintf('%.2f', mean(internal_time, na.rm = T))
  )
  
  return(cols)
}

red_cols <- clean_and_summarise(reduced_pya, 'ordered_step_off') %>%
  filter(algorithm != 'ordered_step_off')
fast_cols <- clean_and_summarise(filter(
  pya_csv, algorithm != 'cplex2nt' & algorithm != 'gurobi0nt'
))
pya_cols <- rbind(red_cols, fast_cols)
pya_cols$algorithm <- recode(pya_cols$algorithm, pyasukp = 'eduk2')

rr_csv <- read.csv(
  '../data/realistic_random.csv', sep = ";", stringsAsFactors = FALSE
)
rr_cols <- clean_and_summarise(rr_csv, 'ordered_step_off')

breq_csv <- read.csv(
  '../data/128_16_std_breqd.csv', sep = ";", stringsAsFactors = FALSE
)
breq_cols <- clean_and_summarise(breq_csv, 'cpp-mtu2')

csp_csv <- read.csv(
  '../data/csp_6195.csv', sep = ";", stringsAsFactors = FALSE
) %>% filter(
  algorithm == 'cplex_cutstock' |
  algorithm == 'mtu1_cutstock'  |
  algorithm == 'ordso_int_ns'   |
  algorithm == 'terso_int'
) %>% mutate(
  internal_time = hex_sum_knapsack_time + hex_sum_sort_time
) %>% select(algorithm, filename, internal_time)
csp_cols <- clean_and_summarise(csp_csv)
csp_cols$algorithm <- recode(csp_cols$algorithm,
  # not really cplex2nt, but it is cplex, and this is the code
  # we refer to it in the other csvs
  cplex_cutstock = 'cplex2nt',
  mtu1_cutstock = 'cpp-mtu1',
  ordso_int_ns = 'ordered_step_off',
  terso_int = 'terminating_step_off'
)

create_empty_rows <- function (existing, universe) {
  forgotten <- vector(mode = 'character')
  for (u in universe) {
    if (!(u %in% existing)) { forgotten <- append(u, forgotten) }
  }
  return(data.frame(
    algorithm = forgotten,
    finished = '--',
    mean_time = '--'
  ))
}

add_missing <- function (dt) {
  all_algs = c(
    'cpp-mtu1', 'cpp-mtu2', 'cplex2nt', 'gurobi0nt', 'ordered_step_off',
    'terminating_step_off', 'mgreendp', 'eduk', 'eduk2'
  )
  mock <- create_empty_rows(dt$algorithm, all_algs)
  return(rbind(dt, mock))
}


pya_cols <- add_missing(pya_cols) %>% rename(
  pya_fin = finished,
  pya_mtm = mean_time
)
rr_cols <- add_missing(rr_cols) %>% rename(
  rr_fin = finished,
  rr_mtm = mean_time
)
breq_cols <- add_missing(breq_cols) %>% rename(
  breq_fin = finished,
  breq_mtm = mean_time
)
csp_cols <- add_missing(csp_cols) %>% rename(
  csp_fin = finished,
  csp_mtm = mean_time
)

one_table <- inner_join(pya_cols, rr_cols, by = 'algorithm')
one_table <- inner_join(one_table, breq_cols, by = 'algorithm')
one_table <- inner_join(one_table, csp_cols, by = 'algorithm')
one_table$algorithm <- recode(
  one_table$algorithm,
  eduk  = 'EDUK',
  eduk2 = 'EDUK2',
  `cpp-mtu1` = 'MTU1',
  `cpp-mtu2` = 'MTU2',
  gurobi0nt = 'Gurobi',
  cplex2nt = 'Cplex',
  ordered_step_off = 'OSO',
  terminating_step_off = 'TSO',
  mgreendp = 'GFDP'
)
one_table$algorithm <- factor(
  one_table$algorithm,
  levels = c(
    'TSO', 'OSO', 'GFDP', 'EDUK2', 'EDUK',
    'MTU2', 'MTU1', 'Cplex', 'Gurobi'
  )
)
one_table <- one_table %>% arrange(algorithm)

library(xtable)

print.xtable(
  xtable(one_table),
  only.contents = T,
  include.rownames = F,
  include.colnames = F,
  hline.after = NULL,
  comment = F
)
