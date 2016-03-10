library(ggplot2)
library(grid)
library(gridExtra)

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

files <- c('postponed_per/nsds2_all.all', 'saw/saw_all.all', 'sc/sc_all.all', 'ss2/ss2_all.all', 'wcd/hi/wcd_hi_all.all', 'all.all')
title <- c('Postponed Periodicity', 'SAW', 'Strong Correlation', 'Subset Sum', 'No Collective Dominance', 'All Datasets')

plots <- list()
for (i in 1:6) { 
	t <- read.table(files[i], header=T, sep=";") 
	t$filename = NULL
	t$ukp5_ext_time = NULL
	t$ukp5_mem = NULL
	t$ukp5_opt = NULL
	t$pya_ext_time = NULL
	t$pya_mem = NULL
	t$pya_opt = NULL
	t <- t[order(t$pya_time), ]

	size <- dim(t)[1]
	ix <- c(1:size, 1:size)
	times <- c(t$pya_time, t$ukp5_time)
	col <- c(rep('PYAsUKP', size), rep('UKP5', size))
	t2 <- data.frame(ix, times, col)

	p = 	ggplot(t2, aes(x = ix, y = times, color = col)) +
		ggtitle(title[i]) +
		scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100, 1000), labels = c('0.01', '0.1', '1', '10', '100', '1000'), limits = c(0.01, 1000)) + 
		#scale_y_log10(0, 1000) + 
		#scale_y_continuous(labels = math_format()) +
		geom_point(aes(color = factor(col))) +
		scale_color_manual(values = c("black", "gray")) +
		theme(legend.position="none", axis.title.y = element_blank(), axis.title.x = element_blank(), plot.title = element_text(size = 10))
	plots[[i]] <- p

	ggsave(filename = paste(files[i], ".png", sep=""), plot = p, width = (210/3), height = (297/3), units = c('mm'), dpi = 300)
}

pn <- grid.arrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], ncol = 3)
ggsave(filename = './six_plots.png', plot = pn, width = 210*0.75, height = 297*0.4, units = c('mm'), dpi = 300)

