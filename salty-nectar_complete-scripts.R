# Plant-pollinator interactions
# Ethan VanValkeburg, Thiago Goncalves Souza, Nathan Sanders, Paul CaraDonna
# Updated 10/2023

# libraries
library(devtools) # tools and github package import
library(readxl) # read excel files
library(tidyverse) # plot, manipulate data
library(bipartite) # organize and analyze network data
library(lme4) # mixed effects modeling
library(car) # stats package
library(glmmTMB) #generalized mixed effects models with a different algorithm
library(vegan) # useful ecology functions — diversity indices
library(gridExtra) # plot side by side
library(lemon) # for plotting
library(patchwork) #organizing plots
library(magrittr) # handling data
library(cowplot) # plotting multiple plots together
library(benthos) # diversity indieces, PIE
library(ggtext) # editing ggplot text with markdown formatting
library(lubridate) # manipulate date data
library(plotrix) #std error and other functions
library(ggrepel) # repel ggplot objects
library(performance)

bee_obs <- readxl::read_excel("code/NaNectar_observations.xlsx") %>% # data import
  # remove extra flower species
  subset(flower.sp != "AQCO") %>%
  subset(flower.sp != "MECI") %>% 
  # update variable names
  mutate(treatment = replace(treatment, treatment == "C", "Control")) %>% 
  mutate(week = week(date) - 26) %>% 
  # remove notes
  dplyr::select(!notes) 
# final dataset
head(bee_obs)

# Variable names
# Labels of sp. name for the plant codes
plant_names <- as_labeller(c(
  'DEBA'="Delphinium barbeyi",
  'ERSP'="Erigeron speciosus",
  'HEQU'="Helianthella quinquenervis",
  'MECI'="Mertensia ciliata",
  'HEMU'='Heliomeris multiflora'))

plant_names_hash <- setNames(as.list(c("Delphinium barbeyi",
                                       "Erigeron speciosus",
                                       "Helianthella quinquenervis",
                                       "Mertensia ciliata",
                                       "Heliomeris multiflora")),
                             c('DEBA',
                               'ERSP',
                               'HEQU',
                               'MECI',
                               'HEMU'))

# treatment labels to reuse
treatment.labels = c("Control", "+ Na")

#saltycol <- c("#abafb3", "#0269c2")
saltycol <- c("#7e83b4", "#fbce24") # morton salt



##################
## Data Summary ##
##################

# total number of visits
sum(bee_obs$no.visits)
# total number of unique visitors
bee_obs %>% 
  filter(no.visits > 0) %>%
  nrow()

no_bombus <- bee_obs %>% 
  filter(!grepl("Bombus", pollinator.sp))
sum(no_bombus$no.visits)
length(no_bombus$no.visits)

# sample size per plant
sample_sizes <- bee_obs %>% 
  group_by(flower.sp, pair.num, date) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  group_by(flower.sp, pair.num) %>% 
  summarise(num = n())
sample_sizes  # sample sizes for each plant

mean(sample_sizes$num) # approximately 3.5 observations per plant. 
std.error(sample_sizes$num)

#total number of flowers spiked
bee_obs %>% 
  group_by(flower.sp, pair.num, treatment, date) %>% 
  summarise(count = mean(fl.count)) %>% 
  ungroup() %>% 
  pull(count) %>% 
  sum()

# calculate the visitation per plant per treatment
# NOT grouping across all observations of a single plant
# a 20 min observation period is considered a single point [FYI:each observation occurred on a different day]
# tot.visits are the total number of flowers visited during the observation periods
# n is the number of unique foraging bouts —> estimate of unique visitors
sum_visits <- bee_obs %>% 
  group_by(flower.sp, pair.num, treatment, date, week) %>% #only ever did a single plant on a day 
  summarise(tot.visits = sum(no.visits),
            n = n()) %>% # n() needs to be corrected to account for 0 visitors
  mutate(n = (tot.visits != 0) * n) %>%  # if there are 0 visits, then there are 0 visitors
  ungroup()

## NO. VISITS GLMM MODELS AND COMPARISONS
# GLMM models (test poisson and the negative binomial)
bee_glme.poisson <- glmmTMB(tot.visits ~ treatment * flower.sp + (1|pair.num), family = poisson, data = sum_visits)  
summary(bee_glme.poisson)
car::Anova(bee_glme.poisson, type = 3) #
check_overdispersion(bee_glme.poisson) # overdispersed

# nb is the favored model
bee_glme.nb.1 <- glmmTMB(tot.visits ~ treatment * flower.sp + treatment * week + (1|pair.num), family=nbinom2, data = sum_visits)  
summary(bee_glme.nb.1) # no significant interaction
car::Anova(bee_glme.nb.1, type = 2) #type = 3 for interaction term, can go to type = 2 if there is not signifigant interaction
# check for overdispersion
check_overdispersion(bee_glme.nb.1) #OK
AIC(bee_glme.nb.1) 


## Summary statistics for visits and visitation
# Determine the average and error of the visitation and visitors for the different treatments
sum_visits %>% 
  group_by(treatment) %>% 
  summarise(mean = mean(tot.visits),
            se = sd(tot.visits)/sqrt(n()),
            median = median(tot.visits),
            n.mean = mean(n),
            n.se = sd(n)/sqrt(n()),
            n.median = median(n))

# Average visitation and se by sp.
sum_visits %>% 
  group_by(flower.sp) %>% 
  summarise(mean = mean(tot.visits),
            se = sd(tot.visits/sqrt(n())))



## FIGURE 1a
plant_counts <- bee_obs %>% 
  group_by(treatment, flower.sp, pair.num, date, week) %>% 
  summarise(no.visits = sum(no.visits)) %>% 
  ungroup() %>% 
  group_by(treatment, flower.sp) %>% 
  summarise(no.visits.mean = mean(no.visits),
            no.visits.se = std.error(no.visits)) %>% 
  ungroup() %>% 
  mutate(flower.sp = str_replace(flower.sp, "DEBA", "Delphinium"),
         flower.sp = str_replace(flower.sp, "HEQU", "Helianthella"),
         flower.sp = str_replace(flower.sp, "HEMU", "Heliomeris"),
         flower.sp = str_replace(flower.sp, "ERSP", "Erigeron"))
plant_counts



# plotting for plant focus
plant_visits <- ggplot(data = plant_counts,
                       aes(x = treatment,
                           y = no.visits.mean, 
                           ymin=no.visits.mean-no.visits.se, 
                           ymax=no.visits.mean+no.visits.se,
                           group = flower.sp,
                           color = treatment,
                           fill = treatment,
                           label = flower.sp)) +
  # geom_jitter(data = bee_obs, mapping = aes(x = treatment, y = no.visits,
  #                                           ymin = NULL, ymax = NULL,
  #                                           color = flower.sp),
  #             alpha = 0.1, width = 0.05, height = 0.4) +
  geom_line(aes(group = flower.sp), color = "black", alpha = 0.5, linewidth = 0.75,
            position = position_dodge(width = 0.1)) +
  geom_errorbar(width = 0.1, position = position_dodge(width = 0.1)) +
  # geom_segment(aes(x = treatment, xend = trea+ 0.06, y = y, yend = y), size = 0.2) +
  geom_point(size = 6, shape = 21, color = "white", position = position_dodge(width = 0.1)) +
  geom_text_repel(aes(label=ifelse((treatment == 'Na'), as.character(flower.sp),'')),
                  hjust= 0,
                  vjust= 0.5,
                  size = 3,
                  color = "black",
                  segment.color = NA,
                  fontface = "italic",
                  direction = "y",
                  arrow = NULL,
                  seed = 3,
                  nudge_x = 0.12) +
  # NOTE: https://stackoverflow.com/questions/47492191/aligning-labels-with-ggrepel
  #ylim(0, 18)  +
  scale_x_discrete(expand = expansion(mult = c(0.3, .4)),
                   labels = treatment.labels) + 
  scale_color_manual(values = saltycol) +
  scale_fill_manual(values = saltycol) +
  xlab("")+
  ylab("Mean number of visits")+
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(linewidth = 0.25),
    legend.position = "none",
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18)
  )
plant_visits

## VISITATION BY POLLINATOR SPECIES
# alternative
bee_counts <- bee_obs %>% 
  group_by(treatment, date, flower.sp, pollinator.sp) %>% 
  summarise(no.visits = sum(no.visits)) %>% 
  ungroup() %>% 
  # make sure I have a record of zero interactions
  pivot_wider(names_from = pollinator.sp, values_from = no.visits, values_fill = 0) %>% 
  # pivot back
  pivot_longer(cols = !treatment:flower.sp, names_to = "pollinator.sp", values_to = "no.visits") %>% 
  # remove Bombus individuals not identified to species
  filter(!pollinator.sp %in% c("Bombus.sp", "0")) %>% 
  # group by pollinator species, calculate the total number of visits
  group_by(pollinator.sp) %>% 
  mutate(total = sum(no.visits)) %>% 
  ungroup() %>% 
  # group by pollinator species and treatment and calculate the total number of visits for each treatment 
  group_by(treatment, pollinator.sp, total) %>% 
  summarise(num.visitors = sum(no.visits > 0),
            no.visits.mean = mean(no.visits),
            no.visits.se = 0) %>% 
  ungroup() %>% 
  # filter out species with not enough observations
  filter(total > 20) %>% 
  # adjust zero values so they appear on the plot (NOTE: addition not in the analysis)
  mutate(no.visits.mean = no.visits.mean + 0.1) %>% 
  # modify names
  mutate(pollinator.sp = str_replace_all(pollinator.sp, "\\.", " "))

head(bee_counts)  

# prepare for the paired test (wilcoxon)
# calculate the difference between treatments
bee_counts_paired <- bee_counts %>%
  select(c(treatment, no.visits.mean, pollinator.sp)) %>% 
  pivot_wider(names_from = treatment, values_from = no.visits.mean) %>% 
  mutate(diff = Na - Control,
         prop = Na / Control) %>% 
  mutate(prop = ifelse(prop == Inf, 0, prop))

# lme to determine whether pollinator species have more visits on sodium-enriched plants
# repeat most of what is done above
bee_obs_filtered <- bee_obs %>% 
  group_by(treatment, date, flower.sp, pollinator.sp) %>% 
  # calculate number of visits by a species during an observation period
  summarise(no.visits = sum(no.visits)) %>% 
  ungroup() %>% 
  # make sure I have a record of NO interactions
  pivot_wider(names_from = pollinator.sp, values_from = no.visits, values_fill = 0) %>% 
  pivot_longer(cols = !treatment:flower.sp, names_to = "pollinator.sp", values_to = "no.visits") %>%
  # remove Bombus spp. not identified to species
  filter(!pollinator.sp %in% c("Bombus.sp", "0")) %>% 
  # calculate values (number of visits) for each pollinator species
  group_by(pollinator.sp) %>% 
  mutate(total = sum(no.visits)) %>% 
  # use species with sufficient visitation
  filter(total > 20) %>% 
  ungroup()

# glmer
m <- glmer(no.visits ~ treatment + (1 | pollinator.sp), family = poisson, data = bee_obs_filtered)
summary(m)
car::Anova(m, type = 2)

## FIGURE 1b
bombus_visits <- ggplot(data = bee_counts,
                        aes(x = treatment,
                            y = no.visits.mean,
                            ymin = no.visits.mean-no.visits.se,
                            ymax = no.visits.mean+no.visits.se,
                            group = pollinator.sp,
                            fill = treatment,
                            color = treatment)) +
  geom_line(aes(group = pollinator.sp), color = "black", alpha = 0.5, size = 0.75,
            position = position_dodge(width = 0.2)) +
  geom_point(size = 5, shape = 21, color = "white", position = position_dodge(width = 0.2)) +
  geom_errorbar(size = 0.5, width = 0.2, position = position_dodge(width = 0.2)) +
  geom_text_repel(aes(label=ifelse((treatment == 'Na'), as.character(pollinator.sp),'')),
                  hjust = 0,
                  vjust = 0.5,
                  nudge_x = 0.12,
                  size = 3,
                  color = "black",
                  segment.color = NA,
                  fontface = "italic",
                  direction = "y",
                  # position = position_dodge(width = 0.2),
                  arrow = NULL,
                  seed = 3) +
  scale_x_discrete(expand = expansion(mult = c(0.3, .6)),
                   labels = treatment.labels) + 
  scale_color_manual(values = saltycol) +
  scale_fill_manual(values = saltycol) +
  scale_y_continuous(trans='log2') +
  xlab("")+
  ylab("Mean number of visits")+
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(linewidth = 0.25),
    legend.position = "none",
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18)
  )

bombus_visits

## FIGURE 1
plot_grid(plant_visits, NULL, bombus_visits, ncol = 3,
          rel_widths = c(5,0.5,5),
          labels = c("a)", "", "b)"),
          label_size = 18)
          



### FIGURE 2

### NMDS among plant species
# create a community matrix with each cell representaing a single week of observations
comm_matrix <- bee_obs %>% 
  dplyr::select(flower.sp, pair.num, treatment, pollinator.sp, no.visits) %>% 
  pivot_wider(names_from = pollinator.sp,
              values_from = no.visits,
              values_fn = list(no.visits = sum),
              values_fill = 0) %>% 
  # remove the blank species
  dplyr::select(-'0') %>% 
  #mutate(sum = sum(.[5:ncol(.)])) %>% 
  filter_at(c(4:27), any_vars(. != 0)) # may need to be removed
comm_matrix

# simplify labels
comm_matrix1 <- comm_matrix %>% 
  unite("label" , flower.sp : treatment, remove = TRUE) %>% 
  tibble::column_to_rownames(var = "label") %>% 
  na.omit()
comm_matrix1
dim(comm_matrix1)

# Data check
sum(rowSums(comm_matrix1)==0)
glimpse(comm_matrix1)
nvisits <- rowSums(comm_matrix1); nvisits # total visits on each plant in each week
test <- data.frame(comm_matrix[,1:4], nvisits = nvisits);test # complete summary

# extract variables
comm_vars <- comm_matrix %>% 
  dplyr::select(c(1:4))

# create NMDS
NMDS <- metaMDS(comm_matrix1, k = 2, distance = "bray")
plot(NMDS) # preview

# create final NMDS plot
### PLOT NMDS with scores
scores <- as.data.frame(scores(NMDS)$sites)

#a dd columns to data frame 
# extract the scores for the plot
scores$flower.sp = comm_matrix$flower.sp
scores$pair.num = comm_matrix$pair.num
scores$treatment = comm_matrix$treatment
#scores$week = comm_matrix$week

# retrieve the variable names
scores$combined = comm_matrix %>% 
  unite("flower_treatment", c(flower.sp, treatment), remove = TRUE, sep = "_") %>% 
  .$flower_treatment


# plot the NMDS with ggplot
xx = ggplot(scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes( shape = flower.sp, colour = treatment, size = flower.sp)) + 
  theme(axis.text.y = element_text(colour = "black", size = 12), 
        axis.text.x = element_text(colour = "black",  size = 12), 
        legend.text = element_markdown(), 
        legend.position = c(0.81, 0.25), 
        axis.title.y = element_text(size = 16), 
        axis.title.x = element_text( size = 16, colour = "black"), 
        legend.title = element_text(size = 12, colour = "black"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "Treatment", y = "NMDS2", shape = "Plant Species")  + 
  scale_colour_manual(values = saltycol, labels = c("Control", "Na-enriched")) +
  scale_shape_manual(values = c(16, 17, 15, 18), 
                     labels = c("*D. babeyi*", "*E. speciosus*", "*H. multiflora*", "*H. quinquenervis*")) +
  scale_size_manual(values = c(4,4,4,5)) +
  guides(size = FALSE,
         shape = guide_legend(override.aes = list(size = 3)),
         color = guide_legend(override.aes = list(size = 3)))

xx

## BETADISPER 
# Multivariate homogeneity of groups dispersions (variances)
# multivaraite version of levene's test (checks whether several groups have the same variance in the population)
# summary: Betadisper first calculates the average distance of group members to the group centroid in multivariate 
# space (generated by a distance matrix). Then, an ANOVA is done to test if the dispersions (variances) of 
# groups are different.
# for adonis: when we know that groups have homogenous variances, we can test if compositions are different 

## PERMDISP
diss <- vegdist(comm_matrix1, method = "bray")

# groups are control/sodium enriched (from comm_matrix)
treatment.groups <- as.factor(comm_matrix$treatment)
species.groups <- as.factor(comm_matrix$flower.sp)
pair.groups <- as.factor(comm_matrix$pair.num)

beta_treatment <- betadisper(diss,
                             group = treatment.groups)

beta_sp <- betadisper(diss,
                      group = species.groups) 


# perform test
anova(beta_treatment) # p = 0.5 Homogenously dispersed
plot(beta_treatment)

anova(beta_sp) # p = 0.01 [[plant species]]
plot(beta_sp)



### ADONIS
# Permutational Multivariate Analysis of Variance Using Distance Matrices
# Summary: Adonis analyzes and partitions sums of squares using distance matrices. 
# It can be seen as an ANOVA using distance matrices (analogous to MANOVA - multivariate 
# analysis of variance). Therefore, it is used to test if two or more groups have similar compositions.
## using the same distance matrix as is listed above

# using strata; restrict permutations within the blocks of strata (examines the effect of treatment within flower.sp only)
mod1 <- adonis(comm_matrix1 ~ treatment * species.groups, method = "bray", data = comm_vars, 
               strata = pair.groups, permutations = 999)

mod1$aov.tab



## HURLBURT'S PROBABILITY OF INTERSPECIFIC ENCOUNTER
# for each pollinator species; compare between treatments
bee_hpie <- bee_obs %>%
  # remove null observations
  filter(pollinator.sp != 0) %>%
  # remove bombus not identified to species
  filter(pollinator.sp != "Bombus.sp") %>%
  # calculate HPIE for each pollinator species in each treatment
  select(pollinator.sp, flower.sp, treatment, no.visits) %>%
  group_by(pollinator.sp, treatment) %>%
  summarise(hpie = hpie(taxon = flower.sp, count = no.visits),
            num = sum(no.visits),
            div = n_distinct(flower.sp)) %>% 
  group_by(pollinator.sp) %>% 
  filter(max(div) != 1) %>% 
  ungroup() %>% 
  arrange(desc(hpie))

# dumbbell plot
bee_hpie_1 <- bee_hpie %>% 
  group_by(pollinator.sp) %>% 
  # calculate the total number of observations for each species
  mutate(num = sum(num)) %>% 
  select(-div) %>% 
  ungroup() %>% 
  # place them side by side
  pivot_wider(names_from = treatment,
              values_from = hpie,
              values_fill = 0) %>% 
  # Replace NA values with 0 evenness
  mutate(Control = replace(Control, Control == "NaN", 0)) %>% 
  # Fix names
  mutate(pollinator.sp = gsub(".", " ", pollinator.sp, fixed = TRUE))

#View(bee_hpie_1)
bee_hpie_1

mean(bee_hpie_1$Control)
std.error(bee_hpie_1$Control)
mean(bee_hpie_1$Na)
std.error(bee_hpie_1$Na)



# Wilcox test
hpie_wilc <- wilcox.test(x = bee_hpie_1$Control, y = bee_hpie_1$Na, 
                         alternative = "less", paired = TRUE)
hpie_wilc




#library(ggalt)
bee_hpie_bell <- ggplot(bee_hpie_1, aes(y = reorder(pollinator.sp, num))) +
  geom_segment(aes(x = -0.05, xend = 0.79,
                   y = reorder(pollinator.sp, num), yend = reorder(pollinator.sp, num)),
               size = 0.5, color = "lightgrey", alpha = 0.5) +
  geom_segment(aes(x = Control, xend = Na,
                   y = reorder(pollinator.sp, num), yend = reorder(pollinator.sp, num)),
               size = 8, lineend = "round", alpha = 0.5,
               color = "black") +
  geom_segment(aes(x = Control, xend = Na,
                   y = reorder(pollinator.sp, num), yend = reorder(pollinator.sp, num)),
               size = 7, lineend = "round",
               color = ifelse((bee_hpie_1$Na - bee_hpie_1$Control > 0), "darkgrey", "white")) +
  geom_point(aes(x = Control, y = reorder(pollinator.sp, num)), size = 7, color = saltycol[1]) +
  geom_point(aes(x = Na, y = reorder(pollinator.sp, num)), size = 7, color = saltycol[2]) +
  scale_x_continuous(limits = c(-0.05, 0.79), expand = c(0,0)) +
  ylab("") +
  xlab("Probability of Interspecific Encounter") +
  theme_bw() + 
  theme(
    panel.grid = element_blank(),
    #panel.border = element_blank(),
    axis.line = element_blank(),
    legend.position = "none",
    axis.text.y = element_text(face = "italic", size =12),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.ticks.y = element_blank()
  )

# FIGURE 2

plot_grid(xx, NULL, bee_hpie_bell, ncol = 3,
          rel_widths = c(5, 0.01, 4),
          labels = c("a)", "", "b)"),
          label_size = 18)



## SUPPLEMENTAL FIGURES
### Temporal visitation plot
## Pollinator visits by week and, treatment, and species
## Line plot

# Function to make a plot for each plot
# input is data for a single species
make_line_wk <- function(dat) {
  # summarize data for individual points and error bars (line plot)
  sum_dat <- dat %>% 
    group_by(treatment, week) %>% 
    summarise(mean = mean(tot.visits), 
              sd = sd(tot.visits), 
              se = (sd(tot.visits) / sqrt(n())))
  
  plant.sp <- plant_names_hash[[dat$flower.sp[[1]]]]
  
  # main plot by week for a single flower
  main <- ggplot(data = sum_dat, aes(x = week, y = mean, 
                                     group = treatment, color = treatment)) +
    geom_line() + # lines connecting means
    geom_point(aes()) + #points for means
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,
                  position=position_dodge(0.05)) + #standard error
    geom_jitter(data = dat, alpha = 0.25, 
                height = 0, width = 0.5,
                aes(x = week, y = tot.visits)) + #individual observations
    scale_color_manual(values = c(saltycol, "grey")) +
    xlab("Week") +
    xlim(1, 5) + # weeks since the start of the experiment
    ylim(min(dat$tot.visits),max(dat$tot.visits)) +
    ylab("Pollinator visits") +
    labs(title = plant.sp) +
    theme_bw() +
    theme(
      # axis.line = element_line(),
      panel.grid = element_blank(),
      strip.text = element_text(face = "italic", size = rel(0.7)),
      strip.background = element_blank(),
      legend.position = "None",
      axis.text.x = element_text(vjust = 0.5),
      axis.title.x = element_blank(),
      plot.title = element_text(face = "italic", size = 10, hjust = 0.5)
    )
  
  # summary plot for the entire sampling period
  summary <- ggplot(data = dat, aes(x = as.factor(treatment), tot.visits, fill = treatment)) +
    geom_boxplot() + # boxplot
    scale_fill_manual(values=saltycol, labels = treatment.labels) +
    scale_color_manual(values = c("black")) +
    ylim(min(dat$tot.visits),max(dat$tot.visits)) +
    theme_bw() +
    theme(#axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      # axis.line = element_line(),
      panel.grid = element_blank(),
      # axis.ticks.y = element_blank(),
      strip.text = element_text(face = "italic", size = rel(0.7)),
      strip.background = element_blank(),
      legend.position = "None",
      axis.text.x = element_text(angle = 90, vjust = 0.5),
      axis.title.x = element_blank())
  
  main + summary + plot_layout(ncol = 2, widths = c(3, 1)) # organize them side by side
} 

# make individual plots, store in a list  
wk_line <- sum_visits %>% 
  split(.$flower.sp) %>%  # seperate sum_visits  by flower.sp
  lapply(make_line_wk) # make a plot for each


# combine and layout plots
wk_line <- wk_line %>% 
  wrap_plots() + # collate the plots
  plot_layout(guides = "collect") # collect them into a single plot

wk_line # view


## NUMBER OF UNIQUE VISITORS
####
beevisits_glme.nb <- glmmTMB(n ~ treatment * flower.sp + treatment * week + (1|pair.num), family=nbinom2, data = sum_visits)  
summary(beevisits_glme.nb)
car::Anova(beevisits_glme.nb, type = 2) #type = 3 for interaction term, can go to type = 2 if there is not signifigant interaction
AIC(beevisits_glme.nb) #looks like negative binomial model is a better fit (lower AIC) -- marginal


visitors <- ggplot(sum_visits, aes(x = treatment, y = n, fill = treatment)) +
  geom_boxplot(width =0.3, outlier.color = 'darkgrey', outlier.size = 0.8) +
  scale_x_discrete(labels = treatment.labels) +
  stat_summary(fun.y="mean",color="black", shape=16, size = 0.75) +
  scale_fill_manual(values=saltycol, labels = treatment.labels) +
  scale_color_manual(values = c("black")) +
  facet_wrap(~ flower.sp, scales = "free_y", ncol = 4, labeller = plant_names) +
  xlab(NULL) +
  ylab("Distinct Foraging Bouts") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold.italic", size = rel(0.7)),
    strip.background = element_blank(),
    legend.position = "None",
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  )
visitors
####

#### DIVERSITY FIGURES
# create community matrix 
# each 20 minute observation period gets its own community matrix
bee_comm_mat <- bee_obs %>% 
  # group_by(flower.sp, date, treatment, pollinator.sp, week) %>% ## too uneven
  group_by(flower.sp, pair.num, treatment, date, pollinator.sp) %>%
  summarise(no.visits = sum(no.visits),
            here = 1) %>% 
  ungroup() %>% 
  spread(key = pollinator.sp, value = no.visits) %>% 
  select(-'0') %>% 
  replace(is.na(.), 0) 
bee_comm_mat # complete community matrix [ok]
# calculate the diversity indices for the bee communities
# using vegan for the div
bee_div <- bee_comm_mat %>% 
  mutate(shannon = diversity(bee_comm_mat[,6:ncol(bee_comm_mat)], index = "shannon"),
         simpson = diversity(bee_comm_mat[,6:ncol(bee_comm_mat)], index = "simpson"),
         richness = specnumber(bee_comm_mat[,6:ncol(bee_comm_mat)])) %>% 
  select(flower.sp, pair.num, date, treatment, shannon, simpson, richness) %>% 
  # add week back in from date
  mutate(week = week(date))
bee_div # with diversity metrices

# summarise the mean and se for diversity
bee_div %>% 
  group_by(treatment) %>% 
  summarise(mean = mean(richness),
            se = sd(richness)/sqrt(n()),
            simpmean = mean(simpson),
            simpse = sd(simpson)/sqrt(n()))



# Plot the differences in diversiyt
# make two plats, one for each, ad then combine.

# SPecies richness
ggrich <- ggplot(bee_div, aes(x = treatment, y = richness, fill = treatment)) +
  geom_boxplot(width =0.3, outlier.color = 'darkgrey', outlier.size = 0.8) + 
  scale_x_discrete(labels = treatment.labels) +
  scale_fill_manual(values=saltycol, labels = treatment.labels) +
  scale_color_manual(values = c("black")) +
  stat_summary(fun.y="mean",color="black", shape=16, size = 0.75) +
  #xlab("Treatment") +
  ylab("Species Richness") +
  theme_bw() +
  theme(strip.text = element_text(face = "bold.italic", size = rel(0.7)),
        strip.background = element_blank(),
        legend.position = "None",
        axis.text.x = element_text(angle = 90, vjust = 0.5),
        axis.title.x = element_blank(),
        panel.grid = element_blank())


# simpson's diversity index from vegan
ggsimp <- ggplot(bee_div, aes(x = treatment, y = simpson, fill = treatment)) +
  geom_boxplot(width =0.3, outlier.color = 'darkgrey', outlier.size = 0.8) + 
  scale_x_discrete(labels = treatment.labels) +
  scale_fill_manual(values=saltycol, labels = treatment.labels) +
  scale_color_manual(values = c("black")) +
  stat_summary(fun.y="mean",color="black", shape=16, size = 0.75) +
  #xlab("Treatment") +
  ylab("Simpson's Index") +
  theme_bw() +
  theme(strip.text = element_text(face = "bold.italic", size = rel(0.7)),
        strip.background = element_blank(),
        legend.position = "None",
        axis.text.x = element_text(angle = 90, vjust = 0.5),
        axis.title.x = element_blank(),
        panel.grid = element_blank())

# arrange them side-by-side
divplot <- grid.arrange(ggrich, ggsimp, ncol=2); divplot

ggrich


# Compare the diversity between the two treatments
# comparisons of richness

div_glme <- glmmTMB(richness ~ treatment + (1|pair.num), 
                    family=poisson, data = bee_div)  

summary(div_glme)
car::Anova(div_glme, type = 2) #type = 3 for interaction term, can go to type = 2 if there is not signifigant interaction

AIC(div_glme)


simpson_lme <- lmer(simpson ~ treatment * flower.sp + (1|date), data = bee_div)
summary(simpson_lme)
car::Anova(simpson_lme, type = 3)
AIC(simpson_lme)





## HANDLING TIME ##
# handling time is the length of foraging bout/number of flowers visited. ie. How long it spends on each flower.
#####
# Does handling time differ between treatments?
# total time on the flower unit
# time spent foraging on each flower
# time per flower
bee_obs_handling <- bee_obs %>% 
  filter(handling != "-") %>% 
  mutate(handling = as.numeric(handling)) %>% 
  mutate(time.on.flower = replace_na(handling / no.visits, 0)) # should I be replacing these with 0 or eliminating them?

# limited just to bumble bees on DEBA and HEQU
bee_obs_handling <- bee_obs_handling %>% 
  filter(grepl("Bombus", pollinator.sp)) %>% 
  filter(flower.sp != "ERSP", flower.sp != "MECI", flower.sp != "HEMU")


mod.handling <- lmer(time.on.flower ~ treatment * flower.sp + (1|pair.num), data = bee_obs_handling)
summary(mod.handling)
car::Anova(mod.handling, type = 2)


# PLOT HANDLING TIME
handling <- ggplot(bee_obs_handling, aes(x = treatment, y = time.on.flower, fill = treatment)) +
  geom_boxplot(aes(color = treatment), width = 0.35) +
  facet_wrap(.~ flower.sp, scales = "free_y", ncol = 2, labeller = plant_names) +
  scale_x_discrete(labels = treatment.labels) +
  scale_fill_manual(values=saltycol) +
  scale_color_manual(values=c("black", "black")) +
  #geom_point() +
  xlab("Treatment") +
  ylab("Handling time estimate (visits/time)") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold.italic", size = rel(0.7)),
    strip.background = element_blank(),
    legend.position = "None",
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  )
handling

ggsave("handling.png",
       handling,
       width = 10,
       height = 16,
       units = "cm",
       dpi = 300)

# can I do this and test if it only significant for DEBA? 
# It isn't for both species, but there does seem to be a distinct relationship
bee_obs_handling_DEBA <- bee_obs_handling %>% 
  filter(grepl("Bombus", pollinator.sp)) %>%
  filter(flower.sp != "ERSP", flower.sp != "MECI", flower.sp != "HEMU") %>% 
  filter(flower.sp == "DEBA")

mod.handling <- lmer(time.on.flower ~ treatment + (1|pair.num), data = bee_obs_handling_DEBA)
summary(mod.handling)
car::Anova(mod.handling, type = 2)


