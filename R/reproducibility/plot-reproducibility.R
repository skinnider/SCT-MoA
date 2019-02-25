# Plot reproducibility test statistics. 
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(pheatmap)
source("R/theme.R")
source("R/Functions.R")

# read data
dat = read.delim("data/reproducibility/pancreas.txt")
# clean up coefficients
dat$coefficient %<>% clean_metric() %>% as.character()

# plot as heatmap instead 
breaks = c(-20, 0, 10, 20, 30, 40, 50, 100, 1000)
cut = dat %>%
  # group Z scores into categories
  mutate(z_cat = cut(z_score, breaks = breaks)) %>%
  # convert to matrix
  group_by(coefficient, z_cat) %>%
  summarise(n = n()) %>%
  ungroup() 
mat = cut %>%
  spread(z_cat, n) %>%
  column_to_rownames('coefficient') %>%
  as.matrix() %>%
  replace(., is.na(.), 0)

# add missing points
for (coef in unique(cut$coefficient)) {
  for (cat in unique(cut$z_cat)) {
    if (with(cut, sum(z_cat == cat & coefficient == coef)) == 0) {
      cut %<>% rbind(list(coefficient = coef, z_cat = cat, n = 0))
    }
  }
}

# plot heatmap
levels = with(dat, reorder(coefficient, z_score, median)) %>% levels()
p1 = ggplot(cut, aes(x = z_cat, y = factor(coefficient, levels = levels),
                    fill = n)) +
  geom_tile() +
  scale_x_discrete("Z score", expand = c(0, 0),
                   labels = c("< 0", levels(cut$z_cat)[-c(1,8)], ">100")) +
  scale_y_discrete(labels = c(
    levels[1:15], expression(rho[p]), expression(phi[s]))) +
  scale_fill_scico(palette = 'devon', limits = c(0, 20),
                   name = 'Network pairs', expand = c(0, 0)) +
  theme_sc + 
  theme(axis.title.y = element_blank(),
        axis.ticks = element_blank())
p1

# plot bar chart
summary = dat %>%
  group_by(coefficient) %>% 
  summarise(z_score = median(z_score))
p2 = ggplot(summary, aes(x = factor(coefficient, levels = levels), y = z_score,
                         fill = 0.8)) + 
  geom_col() +
  scale_fill_scico(palette = 'devon', limits = c(0, 1.35), guide = F) +
  coord_flip() +
  scale_y_continuous("median Z score", limits = c(0, 100),
                     expand = c(0, 0)) +
  theme_sc + 
  theme(axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.text.y = element_blank(),
        axis.line.x = element_line(color = 'grey50'))
p2

p = p1 + p2 + plot_layout(widths = c(1.6, 1))
ggsave("fig/Figure_3.pdf", p, width = 12, height = 9, units = 'cm')

# how many pairs were significant, per coefficient? 
threshold = abs(qnorm(0.05 / n_distinct(with(dat, paste(id1, id2, cell_type)))))
dat %>% 
  group_by(coefficient) %>%
  summarise(n_sig = sum(z_score > abs(threshold))) %>%
  ungroup() %>%
  arrange(n_sig)
