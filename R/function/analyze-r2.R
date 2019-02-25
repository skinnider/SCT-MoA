# Analyze the r2 of analytical and experimental factors in explaining variation
# in overall functional coherence.
setwd("~/git/SCT-MoA")
options(stringsAsFactors = F)
library(ontologyIndex)
library(tidyverse)
library(magrittr)
library(broom)
source("R/Functions.R")

# read AUROC data
auc = read.delim("data/function/auroc.txt.gz") %>%
  mutate(dataset = gsub("-expr.*$", "", dataset))

# calculate median AUC
med = auc %>%
  group_by(dataset, coefficient) %>%
  summarise(median_AUC = median(auroc)) %>%
  ungroup()

# tag protocol, and drop very rare ones
prop = read.csv("data/geo/datasets.csv") %>%
  filter(Protocol != "")

# filter to one dataset per publication
opp = read.delim("data/geo/one-dataset-per-publication.txt")
med0 = med %>%
  dplyr::rename(Dataset = dataset) %>%
  left_join(prop, by = 'Dataset') %>%
  filter(Dataset %in% opp$Dataset)

# drop protocols found in only a single dataset
protocol = med0 %>%
  group_by(Dataset) %>%
  dplyr::slice(1) %>%
  pull(Protocol)
one_off = names(which(table(protocol) == 1))
med0 %<>% filter(!Protocol %in% one_off)

# also tag transcript coverage and capture
## obtained from https://teichlab.github.io/scg_lib_structs/
# 3'/5' vs. full-length
med0 %<>% mutate(
  bias = ifelse(Protocol %in% c("Smart-seq2", "SMARTer C1", "SMARTer"), 
                'full-length', ifelse(Protocol == 'STRT-seq', '5\'', '3\'')),
  # capture method
  capture = ifelse(Protocol %in% c(
    "Smart-seq2", "SMARTer C1", "SMARTer", "STRT-seq", "CEL-Seq2"),
    'FACS', ifelse(Protocol %in% c("Drop-seq", "inDrop", "Chromium", "GemCode"),
                   "Droplet", "Nanowell")))

# fit univariate models
vars = c("coefficient", "Protocol", "Cells", "Dropouts", "bias", "capture")
models = map(vars, ~ lm(as.formula(paste("median_AUC ~", .)), dat = med0))
r2 = map(models, glance) %>% 
  bind_rows() %>%
  mutate(variable = vars) %>%
  mutate(x = fct_recode(
    variable, 
    "measure of association" = "coefficient",
    "sequencing protocol" = "Protocol",
    "# of cells" = "Cells",
    "% zeroes" = "Dropouts",
    "transcript coverage" = "bias",
    "cell isolation/capture" = "capture")) %>%
  mutate(fill = p.value < 0.05)

# plot
pal = colorRampPalette(rev(
  brewer.pal(n = 12, name = "Paired")))(17)
p = ggplot(r2, aes(x = reorder(x, -adj.r.squared), y = adj.r.squared,
                  fill = fill)) +
  geom_col() +
  scale_y_continuous(expression(adjusted~r^2), expand = c(0, 0), 
                     limits = c(0, 0.3001)) +
  scale_fill_manual(values = c('TRUE' = pal[17], 'FALSE' = 'grey80'),
                    guide = F) +
  clean_theme +
  theme(axis.title.x = element_blank(),
        axis.line.y = element_line(color = 'grey50'),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(5.5, 5.5, 5.5, 13), "pt"))
p
ggsave("fig/Figure_S5.pdf", p, width = 7, height = 8, units = 'cm',
       useDingbats = F)
