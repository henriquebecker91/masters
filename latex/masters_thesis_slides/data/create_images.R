library(ggplot2)

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
