library(tidyverse)

f3c_data = read_csv("output/Figure3c_data.csv") %>%
  mutate(model_rank = 1:nrow(.))

plot_data = f3c_data %>%
  pivot_longer(cols=2:ncol(.)-2) %>%
  select(-deltaAIC) %>%
  mutate(model_rank = factor(model_rank, levels=10:1)) %>%
  mutate(value = value * 100)

f3c = ggplot() + 
  geom_tile(aes(name, model_rank, fill=value), data=plot_data) +
  geom_text(aes(x=9.75, y=model_rank, label=round(deltaAIC, digits = 2)), 
            size=2, hjust = 0, data=f3c_data) + 
  annotate("text", x=9.75, y=11, label=expression(Delta*AIC), hjust=0, size=2) + 
  scale_y_discrete(name ="Rank",
                   labels = setNames(paste0("", 1:10), 1:10)) + 
  scale_fill_continuous(type = "viridis") +
  coord_cartesian(clip="off") +
  theme(text = element_text(size=6, color="black"), 
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 30, vjust = 1, hjust=1, 
                                   size=5, color="black"),
        axis.title = element_blank(),
        #legend.position = "none",
        plot.margin = margin(t = 2, r = 3, b = 0, l = 0, unit = "mm"))

ggsave(filename="output/Figure3c.pdf", plot=f3c, width=43, height=41, units="mm")