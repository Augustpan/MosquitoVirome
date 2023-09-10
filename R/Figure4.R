source("R/dev/Figure4_related_cst_profiling_dev.R")
source("R/dev/Figure4a_dev.R")
source("R/dev/Figure4b_dev.R")
source("R/dev/Figure4c_dev.R")
source("R/dev/Figure4de_dev.R")
source("R/dev/Figure4fg_dev.R")

library(cowplot)
library(ggplot2)

f4b_placeholder = grid::rectGrob()
f4c_placeholder = grid::rectGrob()

f4a_mod = f4a + 
  theme(legend.position = "none")

f4top = plot_grid(f4a, f4b_placeholder, nrow =1, labels=c("a", "b"))
fsub = plot_grid(f4d, f4e, f4f, f4g,
                 nrow = 2, ncol = 2, align="hv", labels=c("d", "e", "f", "g"))
f4left = plot_grid(f4top, fsub, nrow = 2, rel_heights = c(1,2))


f4 = plot_grid(f4left, f4c_placeholder, 
               rel_widths = c(2.25, 1), labels=c("", "c"))


f4ab_new = plot_grid(f4a_mod, f4b_placeholder, nrow = 2, labels=c("a", "b"))
f4abc_new = plot_grid(f4ab_new, f4b_placeholder, nrow = 1, labels=c("", "c"),  rel_widths = c(1.5, 1))
f4defg_new = plot_grid(f4d, f4e, f4f, f4g,
                       nrow = 1, align="hv", labels=c("d", "e", "f", "g"))
f4_new = plot_grid(f4abc_new, f4defg_new, rel_heights = c(3.8,1),
                   nrow=2)
ggsave(filename = "output/Figure4_merged.pdf", 
       plot = f4_new,
       width = 180,
       height = 225,
       units = "mm")