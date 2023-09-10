source("R/dev/Figure3ab_dev.R")
source("R/dev/Figure3cdhij_dev.R")
source("R/dev/Figure3c_dev.R")
source("R/dev/Figure3g_dev.R")
source("R/dev/Figure3ef_dev.R")

library(cowplot)
library(ggplot2)

#### merge plots ####
f3c_placeholder = grid::rectGrob()
f3b_placeholder = grid::rectGrob()
f3g_placeholder = grid::rectGrob()

f3c_mod = f3c + 
  theme(legend.position = "none", 
        legend.key.width = unit(10, units="mm"),
        legend.key.height = unit(2, units="mm"))

#f3c_mod = f3c_placeholder
f3cd = plot_grid(f3c_mod, f3d, nrow=1, rel_widths= c(1.2, 1), labels=c("c", "d"))
f3bcd = plot_grid(f3b_placeholder, f3cd, nrow=2, labels=c("b", ""))

f3a_mod = f3a + 
  theme(legend.position = "none")

f3abcd = plot_grid(f3a_mod, f3bcd, nrow=1, rel_widths = c(1.25, 1), labels=c("a", ""))

f3e_mod = f3e +
  theme(axis.text = element_text(size=6, color="black"),
        axis.text.x = element_text(face="italic"),
        axis.title.y = element_text(size=7),
        axis.title.x = element_blank())
f3f_mod = f3f +
  theme(axis.text = element_text(size=6, color="black"),
        axis.text.x = element_text(face="italic"),
        axis.title.y = element_text(size=7),
        axis.title.x = element_blank())

f3g_mod = f3g +
  theme(axis.text = element_text(size=6, color="black"),
        axis.title.y = element_text(size=7),
        axis.title.x = element_text(size=7))

f3efg = plot_grid(f3e_mod, f3f_mod, f3g_mod, 
                  nrow=1, rel_widths = c(1, 1, 0.75), labels=c("e", "f", "g"), align="hv")

f3h_mod = f3h +
  ylab("Viral richness per individual") +
  xlab("Annual Mean Temperature") +
  theme_bw()+
  theme(axis.text = element_text(size=6, color="black"),
        axis.title = element_text(size=7))
f3i_mod = f3i +
  ylab("Viral richness per individual") +
  xlab("Annual Precipitation") +
  theme_bw()+
  theme(axis.text = element_text(size=6, color="black"),
        axis.title = element_text(size=7))
f3j_mod = f3j +
  ylab("Viral richness per individual") +
  xlab("Mammal richness") +
  theme_bw()+
  theme(axis.text = element_text(size=6, color="black"),
        axis.title = element_text(size=7))
f3hij = plot_grid(f3h_mod, f3i_mod, f3j_mod, 
                  nrow=1, rel_widths = c(1, 1, 1), labels=c("h", "i", "j"))

f3 = plot_grid(f3abcd, f3efg, f3hij, nrow=3, rel_heights = c(1.5, 1, 0.7))

ggsave(filename = "output/Figure3_merged.pdf", 
       plot = f3,
       width = 180,
       height = 193,
       units = "mm")
