### Read and evaluate the gabor data

### Libraries
library(tidyverse) # for data wrangling and plotting
library(here)
library(viridis)
library(patchwork)
source('rtfmri_theme.R')

# Path


path = 'E:/Gabor/'

bids_folder = 'C:/Users/Jeane/Documents/Behavior_Analysis/'


files = list.files(path = path, pattern = '.*staircase.csv', recursive = T, full.names = T)

data_df <- read.delim(file = paste(bids_folder,'participants.tsv', sep = ''), header = T, 
                      sep = '\t')

gabor_data <- files %>% 
  set_names() %>% 
  map_df(read_csv, .id = "source") %>% 
  mutate(subid = str_extract(source, ".P..."),
         session = str_extract(source, '(?<=_)[1-4]{1}(?=_)')) %>% #Look ahead and look behind! 
  select(-source,-cumulative_response_time,-iti_onset,-iti_dur, -stim_onset) %>% 
  mutate(session = case_when(session == '2' ~ '3',
                             TRUE ~ session))

data_df$subid <- data_df$participant_id

## Merge Gabor data with group, age, etc.
## Exclude VP001  because really bad. 

gabor_data <-  left_join(gabor_data, data_df, by = 'subid') %>% 
  select(-participant_id) %>% 
  mutate(group = case_when(group == 'lV1' ~ 'lV1',
                           group == 'rV1' ~ 'rV1',
                           group == 'rAIC' ~ 'rAIC',
                           group == 'noFeedback' ~ 'noFeedback'
                           ))

## Select only the final value of the staircases
gabor_final_value <- gabor_data %>% 
  filter(trial %in% c(96,97,98,99,100)) %>% 
  group_by(subid,visfield,session,group) %>% 
  summarise(opacity = mean(opacity)) %>% 
  ungroup()

gabor_final_value<-subset(gabor_final_value, opacity < 1.0)

ggplot(gabor_final_value, aes(x=group, y= opacity, fill=session))+geom_boxplot()


## summarise the final values Short Term effects
sum_gabor_final_value <- gabor_final_value %>% 
  filter(session <= 3) %>%
  group_by(visfield, subid, session, group) %>% 
  summarise(meanOpacity = mean(opacity)) %>% 
  ungroup()


## all sessions 1-4 
sum_gabor_final_value <- gabor_final_value %>% 
  group_by(visfield, subid, session, group) %>% 
  summarise(meanOpacity = mean(opacity)) %>% 
  ungroup()

# Left visual field
gabor_data %>% 
  filter(visfield == 'left') %>% 
  ggplot(aes(x=trial, y = opacity)) +
  geom_point(size = 0.3) +
  geom_line() + 
  facet_grid(subid ~ session) + 
  ggtitle('left visual field') + 
  theme_rtfmri() -> l_vis_field_staircase

# right visual field
gabor_data %>% 
  filter(visfield == 'right') %>% 
  ggplot(aes(x=trial, y = opacity)) +
  geom_point(size = 0.3) +
  geom_line() + 
  facet_grid(subid ~ session) + 
  ggtitle('right visual field') + 
  theme_rtfmri() -> r_vis_field_staircase


l_vis_field_staircase 

r_vis_field_staircase 

 plot_annotation('Gabor Staircase')

# Plot summarised staircases
## Left
sum_gabor_final_value %>% 
  drop_na() %>% 
  group_by(visfield, session, group) %>% 
  summarise(meanOpacity = mean(meanOpacity)) %>% 
  ungroup() %>% 
  ggplot() + 
  aes(x = session, y = meanOpacity, group = group, color = group) +
  geom_point(size = 3) +
  geom_line(size = 2) + 
  geom_point(data = sum_gabor_final_value,
             aes(x = session, y = meanOpacity, group = interaction(group, subid), color = group),
             alpha = 0.3) +
  geom_line(data = sum_gabor_final_value,
            aes(x = session, y = meanOpacity, group = interaction(group, subid), color = group),
            alpha = 0.3) + 
  facet_wrap(.~visfield, labeller = label_both) +
  theme_rtfmri() +
  ggtitle("Gabor Staircase values per group and sessions")


res.aov <- anova_test(
  data = gabor_final_value, dv = opacity, wid = subid,
  between = group, within = c(session,visfield)
)

get_anova_table(res.aov)




gaborRight <- gabor_data %>% 
  filter(visfield == 'right')%>%
  filter(trial %in% c(96,97,98,99,100)) %>% 
  group_by(subid,session,group) %>% 
  summarise(opacity = mean(opacity)) %>% 
  ungroup()

res.aov <- anova_test(
  data = gaborRight, dv = opacity, wid = subid,
  between = group, within = session
)

get_anova_table(res.aov)

gaborLeft <- gabor_data %>% 
  filter(visfield == 'left')%>%
  filter(trial %in% c(96,97,98,99,100)) %>% 
  group_by(subid,session,group) %>% 
  summarise(opacity = mean(opacity)) %>% 
  ungroup()

res.aov <- anova_test(
  data = gaborLeft, dv = opacity, wid = subid,
  between = group, within = session
)

get_anova_table(res.aov)

