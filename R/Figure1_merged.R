library(cowplot)
library(ggplot2)

f1a =  tar_read("plot_figure_1a")
f1b =  tar_read("plot_figure_1b")
f1c =  tar_read("plot_figure_1c")
f1d =  tar_read("plot_figure_1d")
f1e =  tar_read("plot_figure_1e")
f1f =  tar_read("plot_figure_1f")

f1b_mod = f1b + 
  xlim(c(0, 0.55)) + 
  theme(plot.margin = margin(t = 0, r = 0, b = 5, l = 0, unit = "mm"))

f1c_mod = f1c + 
  theme(axis.text = element_text(size=5, color="black"),
        axis.title = element_text(size=6),
        plot.margin = margin(t = 0, r = 15, b = 0, l = 0, unit = "mm"))

f1ab = plot_grid(f1b_mod, f1c_mod, 
                 nrow = 2, 
                 rel_heights = c(1.8, 1), 
                 labels = c("b", "c"))

f1a_mod = f1a + 
  xlim(c(74, 134)) + 
  theme(legend.position = "none",
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "mm"))

f1abc = plot_grid(f1a_mod, f1ab, 
                  nrow = 1, 
                  rel_widths = c(2.25, 1), 
                  labels = c("a", ""))

f1d_mod = f1d + 
  theme(strip.text.x = element_text(size=6, face = "italic"), 
        axis.title.y = element_text(size=6, face = "plain"), 
        
        axis.text = element_text(size=5, color="black"),
        #legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size=6),
        legend.position = c(0, -0.075),
        legend.direction = "horizontal",
        legend.key.size = unit(3, "mm"),
        legend.justification = c("left", "top"),
        legend.box.just = "left",
        legend.margin = margin(0, 0, 0, 0),
        plot.margin = margin(t = 5, r = 7.5, b = 5, l = 5, unit = "mm"))

f1abcd = plot_grid(f1abc, f1d_mod, 
                   nrow = 2, 
                   rel_heights = c(2.5, 1), 
                   labels = c("", "d"))

f1e_placeholder = grid::rectGrob()
f1f_mod = f1f + 
  ylab("Prevalence") +
  guides(fill=guide_legend(ncol=2, by_row=F)) +
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size=5, color="black"),
        axis.title.y = element_text(size=6),
        legend.position = "right",
        legend.justification = c("left", "top"),
        legend.box.just = "left",
        legend.box.margin = margin(0, 0, 0, 0),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.background = element_blank(),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 5, unit = "mm"))

f1ef = plot_grid(f1e_placeholder, f1f_mod, 
                 nrow = 1, 
                 rel_widths = c(2, 1), 
                 labels = c("e", "f"))

f1 = plot_grid(f1abcd, f1ef, 
               nrow = 2, 
               rel_heights = c(4, 1))

ggsave(filename = "output/Figure1_merged.pdf", 
       plot = f1,
       width = 180,
       height = 160,
       units = "mm")
