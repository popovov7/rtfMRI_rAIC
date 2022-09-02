library(tidyverse)
library(ggpubr)
library(rstatix)
library(readxl)
source('rtfmri_theme2.R')

ConnrAICSession <- read_excel("D:/fMRI/AllData/final_Tabel/ConnrAICSession.xlsx")

ConnrAICSession$subid <- ConnrAICSession$id

ConnrAIC_long <- ConnrAICSession%>%
  pivot_longer(cols = r0:r13, names_to = "run", values_to ="nfb_value")

ConnrAIC_long_TB <- ConnrAIC_long%>%
  subset(run == 'r0' | run == 'r13')

ConnrAIC_long_TB<- ConnrAIC_long_TB%>%
  mutate(group = case_when(group == 'lV1' ~ 'lV1',
                             group == 'rV1' ~ 'rV1',
                             group == 'rAIC' ~ 'rAIC',
                             group == 'noFeedback' ~ 'mental-rehearsal'))


## Two-way ANOVA

## Assumptions

# Outliers

Outliers <- ConnrAIC_long_TB %>%
  group_by(run, group) %>%
  identify_outliers(nfb_value)

view(Outliers)

# Normality

ConnrAIC_long_TB %>%
  group_by(run, group) %>%
  shapiro_test(nfb_value)

ggqqplot(ConnrAIC_long_TB, "nfb_value", ggtheme = theme_bw()) +
  facet_grid(run ~ group)


## Homogenity

ConnrAIC_long_TB %>%
  group_by(run) %>%
  levene_test(nfb_value ~ group)

box_m(ConnrAIC_long_TB[, "nfb_value", drop = FALSE], ConnrAIC_long_TB$group)

#ANOVA

res.aov <- anova_test(
  data = ConnrAIC_long_TB, dv = nfb_value, wid = subid,
  between = group, within = run
)
get_anova_table(res.aov)


one.way2 <- ConnrAIC_long_TB %>%
  group_by(group) %>%
  anova_test(dv = nfb_value, wid = subid, within = run) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way2

pwc2 <- ConnrAIC_long_TB %>%
  group_by(group) %>%
  pairwise_t_test(
    nfb_value ~ run, paired = TRUE, 
    p.adjust.method = "bonferroni"
  ) %>%
  select(-df, -statistic, -p) # Remove details
pwc2

##Test normality of residuals saved as attributes

anova_attributes <- attributes(res.aov)
attrmodel <- anova_attributes$args$model
resid <- attrmodel$residuals
shapiro.test(resid)

## p value > 0.05 normal distributed residuals

##Main effect of run

ConnrAIC_long_TB$group<- as.factor(ConnrAIC_long_TB$group)

one.way2 <- ConnrAIC_long_TB %>%
  group_by(group) %>%
  anova_test(dv = nfb_value, wid = subid, within = group) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way2

## Pairwise comparison

pwc2 <- ConnrAIC_long_TB %>%
  group_by(group) %>%
  pairwise_t_test(
    nfb_value ~ run, paired = TRUE, 
    p.adjust.method = "bonferroni"
  )  
pwc2


## Boxplot

ggplot(ConnrAIC_long_TB, aes(x=group, y=nfb_value, color = run))+
  geom_boxplot(outlier.shape = NA)+
  labs( x="Group", y="Contrast value rAIC")+
  geom_point(position=position_jitterdodge(0.2), aes(color = run))+
  scale_y_continuous(breaks=seq(-2,2,by=0.5))+
  scale_color_discrete(name = "Run",labels=c('Baseline Session 02', 'Transfer Session 03'))+
  scale_x_discrete(limits = c("rAIC", "rV1", "lV1", "mental-rehearsal"))+
  theme_rtfmri2()

##Linear Regression just NFB runs

# rAIC group

ConnrAIC_long_NFBruns <- ConnrAIC_long%>%
  subset(run != 'r0' & run != 'r13' & run != 'r6' & run != 'r7')

ConnrAIC_long_NFBruns[c('r', 'Run')] <- str_split_fixed(ConnrAIC_long_NFBruns$run, '', 2)

ConnrAIC_long_NFBruns$Run <- as.numeric(ConnrAIC_long_NFBruns$Run)

ConnrAIC_long_NFBruns['Run'][ConnrAIC_long_NFBruns['Run'] == 8] <- 6
ConnrAIC_long_NFBruns['Run'][ConnrAIC_long_NFBruns['Run'] == 9] <- 7
ConnrAIC_long_NFBruns['Run'][ConnrAIC_long_NFBruns['Run'] == 10] <- 8
ConnrAIC_long_NFBruns['Run'][ConnrAIC_long_NFBruns['Run'] == 11] <- 9
ConnrAIC_long_NFBruns['Run'][ConnrAIC_long_NFBruns['Run'] == 12] <- 10

ConnrAIC_long_NFBruns$Run<-factor(ConnrAIC_long_NFBruns$Run, levels=c(1,2,3,4,5,6,7,8,9,10))

pd <- position_dodge(0.5)

ConnrAIC_long_NFBruns<- ConnrAIC_long_NFBruns%>%
  mutate(group = case_when(group == 'lV1' ~ 'lV1',
                             group == 'rV1' ~ 'rV1',
                             group == 'rAIC' ~ 'rAIC',
                             group == 'noFeedback' ~ 'mental-rehearsal'))

 ConnrAIC_long_NFBruns%>%
  group_by(group, Run)%>% 
  summarise(mean = mean(nfb_value), sd = sd(nfb_value), n = n(), se = sd/sqrt(n))%>%
  ungroup(Run)%>%
  ggplot( aes(x=Run, y = mean, group= group, color=group))+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width =.1, position = pd)+
  geom_line(position = pd, size = .5)+
  geom_point(position = pd, size = 3)+
  geom_smooth(method = 'lm', se=FALSE, size = .1)+
  theme_minimal()+ scale_y_continuous(breaks=seq(-2,2,by=0.2))+
   scale_fill_discrete(breaks=c("rAIC", "rV1", "lV1", "mental-rehearsal"))+
  labs( x="NFB run number", y="Contrast value rAIC")+ scale_color_discrete(name = "Group")+theme_rtfmri2() 
  

ConnrAIC_long_NFBruns_rAIC <- ConnrAIC_long_NFBruns%>%
  subset(group == 'rAIC')

## linear regression rAIC group
ConnrAIC_long_NFBruns_rAIC$Run <- as.numeric(ConnrAIC_long_NFBruns_rAIC$Run)
modelrAIC <- lm(nfb_value ~ Run, ConnrAIC_long_NFBruns_rAIC)
summary(modelrAIC)
cor.test(ConnrAIC_long_NFBruns_rAIC$Run, ConnrAIC_long_NFBruns_rAIC$nfb_value)

## linear regression V1 group
ConnrAIC_long_NFBruns_V1 <- ConnrAIC_long_NFBruns%>%
  subset(group == 'rV1')
ConnrAIC_long_NFBruns_V1$Run <- as.numeric(ConnrAIC_long_NFBruns_V1$Run)
modelV1 <- lm(nfb_value ~ Run, ConnrAIC_long_NFBruns_V1)
summary(modelV1)

## linear regression NF group
ConnrAIC_long_NFBruns_NF <- ConnrAIC_long_NFBruns%>%
  subset(group.y == 'noFeedback')
ConnrAIC_long_NFBruns_NF$Run <- as.numeric(ConnrAIC_long_NFBruns_NF$Run)
modelNF <- lm(nfb_value ~ Run, ConnrAIC_long_NFBruns_NF)
summary(modelNF)



### Whole V1

ConnWholeV1Session <- read_excel("E:/fMRI/AllData/final_Tabel/ConnWholeV1Session.xlsx")

ConnWholeV1Session$subid <- ConnWholeV1Session$id

ConnWholeV1_long <- ConnWholeV1Session%>%
  pivot_longer(cols = r0:r13, names_to = "run", values_to ="nfb_value")

ConnWholeV1_long_TB <- ConnWholeV1_long%>%
  subset(run == 'r0' | run == 'r13')

ConnWholeV1_long_TB<- ConnWholeV1_long_TB%>%
  mutate(group.y = case_when(group == 'lV1' ~ 'V1',
                             group == 'rV1' ~ 'V1',
                             group == 'rAIC' ~ 'rAIC',
                             group == 'noFeedback' ~ 'noFeedback'))


## Two-way ANOVA

## Assumptions

# Outliers

Outliers <-ConnWholeV1_long_TB %>%
  group_by(run, group.y) %>%
  identify_outliers(nfb_value)

view(Outliers)

# Normality

ConnWholeV1_long_TB %>%
  group_by(run, group.y) %>%
  shapiro_test(nfb_value)

ggqqplot(ConnWholeV1_long_TB, "nfb_value", ggtheme = theme_bw()) +
  facet_grid(run ~ group.y)

## Homogenity

ConnWholeV1_long_TB %>%
  group_by(run) %>%
  levene_test(nfb_value ~ group.y)

box_m(ConnWholeV1_long_TB[, "nfb_value", drop = FALSE], ConnWholeV1_long_TB$group.y)

#ANOVA

res.aov <- anova_test(
  data = ConnWholeV1_long_TB, dv = nfb_value, wid = subid,
  between = group, within = run
)
get_anova_table(res.aov)

##Main effect of run

one.way2 <- ConnWholeV1_long_TB %>%
  group_by(group) %>%
  anova_test(dv = nfb_value, wid = subid, within = run) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way2

## Pairwise comparison

pwc2 <- ConnrAIC_long_TB %>%
  group_by(group) %>%
  pairwise_t_test(
    nfb_value ~ run, paired = TRUE, 
    p.adjust.method = "bonferroni"
  )  %>%
  select(-df, -statistic, -p)
pwc2

## Boxplot
ggplot(ConnWholeV1_long_TB, aes(x=group.y, y=nfb_value, color = run))+geom_boxplot()+
  labs( x="Group", y="Contrast value")+geom_point(position=position_jitterdodge(0.2), aes(color = run))+
  scale_y_continuous(breaks=seq(-2,2,by=0.5))+
scale_color_discrete(name = "Run",labels=c('Baseline Session 02', 'Transfer Session 03'))+
  theme_rtfmri2()

##Linear Regression just NFB runs

# rAIC group

ConnWholeV1_long_NFBruns <- ConnWholeV1_long%>%
  subset(run != 'r0' & run != 'r13' & run != 'r6' & run != 'r7')

ConnWholeV1_long_NFBruns[c('r', 'Run')] <- str_split_fixed(ConnWholeV1_long_NFBruns$run, '', 2)
ConnWholeV1_long_NFBruns$Run <- as.numeric(ConnWholeV1_long_NFBruns$Run)
ConnWholeV1_long_NFBruns['Run'][ConnWholeV1_long_NFBruns['Run'] == 8] <- 6
ConnWholeV1_long_NFBruns['Run'][ConnWholeV1_long_NFBruns['Run'] == 9] <- 7
ConnWholeV1_long_NFBruns['Run'][ConnWholeV1_long_NFBruns['Run'] == 10] <- 8
ConnWholeV1_long_NFBruns['Run'][ConnWholeV1_long_NFBruns['Run'] == 11] <- 9
ConnWholeV1_long_NFBruns['Run'][ConnWholeV1_long_NFBruns['Run'] == 12] <- 10

ConnWholeV1_long_NFBruns$Run<-factor(ConnWholeV1_long_NFBruns$Run, levels=c(1,2,3,4,5,6,7,8,9,10))

pd <- position_dodge(0.1)

ConnWholeV1_long_NFBruns<- ConnWholeV1_long_NFBruns%>%
  mutate(group.y = case_when(group == 'lV1' ~ 'V1',
                             group == 'rV1' ~ 'V1',
                             group == 'rAIC' ~ 'rAIC',
                             group == 'noFeedback' ~ 'noFeedback'))

ConnWholeV1_long_NFBruns%>%
  group_by(group.y, Run)%>% 
  summarise(mean = mean(nfb_value), sd = sd(nfb_value), n = n(), se = sd/sqrt(n))%>%
  ungroup(Run)%>%
  ggplot( aes(x=Run, y = mean, group= group.y, color=group.y))+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width =.1, position = pd)+
  geom_line(position = pd, size = .5)+
  geom_point(position = pd, size = 3)+
  geom_smooth(method = 'lm', se=FALSE, size = .1)+
  theme_minimal()+ scale_y_continuous(breaks=seq(-2,2,by=0.1))+
  labs( x="NFB run number", y="Contrast value")+ scale_color_discrete(name = "Group") 


expand_limits(y=c(-0.7, 1))+ scale_y_continuous(breaks=c( -0.7, -0.6,-0.5, -0.4, -0.3, -0.2, -0.1,0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))

ConnWholeV1_long_NFBruns_rAIC <- ConnWholeV1_long_NFBruns%>%
  subset(group == 'rAIC')

## linear regression rAIC group
ConnWholeV1_long_NFBruns_rAIC$Run <- as.numeric(ConnWholeV1_long_NFBruns_rAIC$Run)
modelrAICforV1 <- lm(nfb_value ~ Run, ConnWholeV1_long_NFBruns_rAIC)
summary(modelrAICforV1)


## linear regression V1 group
ConnWholeV1_long_NFBruns_V1 <- ConnWholeV1_long_NFBruns%>%
  subset(group == 'rV1')
ConnWholeV1_long_NFBruns_V1$Run <- as.numeric(ConnWholeV1_long_NFBruns_V1$Run)
modelV1V1 <- lm(nfb_value ~ Run, ConnWholeV1_long_NFBruns_V1)
summary(modelV1V1)

## linear regression V1 group
ConnWholeV1_long_NFBruns_lV1 <- ConnWholeV1_long_NFBruns%>%
  subset(group == 'lV1')
ConnWholeV1_long_NFBruns_lV1$Run <- as.numeric(ConnWholeV1_long_NFBruns_lV1$Run)
modelV1lV1 <- lm(nfb_value ~ Run, ConnWholeV1_long_NFBruns_lV1)
summary(modelV1lV1)


## linear regression NF group
ConnWholeV1_long_NFBruns_NF <- ConnWholeV1_long_NFBruns%>%
  subset(group.y == 'noFeedback')
ConnWholeV1_long_NFBruns_NF$Run <- as.numeric(ConnWholeV1_long_NFBruns_NF$Run)
modelNF <- lm(nfb_value ~ Run, ConnWholeV1_long_NFBruns_NF)
summary(modelNF)



### right V1

ConnrV1Session <- read_excel("C:/Users/Jeane/Documents/PhD/Analysis_Paper1/NFB/ConnrV1Session.xlsx")

ConnrV1Session$subid <- ConnrV1Session$id

ConnrV1_long <- ConnrV1Session%>%
  pivot_longer(cols = r0:r13, names_to = "run", values_to ="nfb_value")

ConnrV1_long_TB <- ConnrV1_long%>%
  subset(run == 'r0' | run == 'r13')

ConnrV1_long_TB<- ConnrV1_long_TB%>%
  mutate(group = case_when(group == 'lV1' ~ 'lV1',
                             group == 'rV1' ~ 'rV1',
                             group == 'rAIC' ~ 'rAIC',
                             group == 'noFeedback' ~ 'mental-rehearsal'))


## Two-way ANOVA

## Assumptions

# Outliers

Outliers <-ConnrV1_long_TB %>%
  group_by(run, group.y) %>%
  identify_outliers(nfb_value)

view(Outliers)

# Normality

ConnrV1_long_TB %>%
  group_by(run, group.y) %>%
  shapiro_test(nfb_value)

ggqqplot(ConnrV1_long_TB, "nfb_value", ggtheme = theme_bw()) +
  facet_grid(run ~ group.y)

## Homogenity

ConnrV1_long_TB %>%
  group_by(run) %>%
  levene_test(nfb_value ~ group.y)

box_m(ConnrV1_long_TB[, "nfb_value", drop = FALSE], ConnrV1_long_TB$group.y)

#ANOVA

res.aov <- anova_test(
  data = ConnrV1_long_TB, dv = nfb_value, wid = subid,
  between = group, within = run
)
get_anova_table(res.aov)

##Main effect of run

one.way2 <- ConnWholeV1_long_TB %>%
  group_by(group) %>%
  anova_test(dv = nfb_value, wid = subid, within = run) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way2

## Pairwise comparison

pwc2 <- ConnrAIC_long_TB %>%
  group_by(group) %>%
  pairwise_t_test(
    nfb_value ~ run, paired = TRUE, 
    p.adjust.method = "bonferroni"
  )  %>%
  select(-df, -statistic, -p)
pwc2

## Boxplot
ggplot(ConnrV1_long_TB, aes(x=group, y=nfb_value, color = run))+geom_boxplot(outlier.shape = NA)+
  labs( x="Group", y="Contrast value rV1")+geom_point(position=position_jitterdodge(0.2), aes(color = run))+
  scale_y_continuous(breaks=seq(-2,2,by=0.5))+
  scale_color_discrete(name = "Run",labels=c('Baseline Session 02', 'Transfer Session 03'))+  scale_x_discrete(limits = c("rAIC", "rV1", "lV1", "mental-rehearsal"))+
  theme_rtfmri2()

##Linear Regression just NFB runs

# rAIC group

ConnrV1_long_NFBruns <- ConnrV1_long%>%
  subset(run != 'r0' & run != 'r13' & run != 'r6' & run != 'r7')

ConnrV1_long_NFBruns[c('r', 'Run')] <- str_split_fixed(ConnrV1_long_NFBruns$run, '', 2)

ConnrV1_long_NFBruns$Run <- as.numeric(ConnrV1_long_NFBruns$Run)

ConnrV1_long_NFBruns['Run'][ConnrV1_long_NFBruns['Run'] == 8] <- 6
ConnrV1_long_NFBruns['Run'][ConnrV1_long_NFBruns['Run'] == 9] <- 7
ConnrV1_long_NFBruns['Run'][ConnrV1_long_NFBruns['Run'] == 10] <- 8
ConnrV1_long_NFBruns['Run'][ConnrV1_long_NFBruns['Run'] == 11] <- 9
ConnrV1_long_NFBruns['Run'][ConnrV1_long_NFBruns['Run'] == 12] <- 10

ConnrV1_long_NFBruns$Run<-factor(ConnrV1_long_NFBruns$Run, levels=c(1,2,3,4,5,6,7,8,9,10))

pd <- position_dodge(0.1)

ConnrV1_long_NFBruns<- ConnrV1_long_NFBruns%>%
  mutate(group = case_when(group == 'lV1' ~ 'lV1',
                             group == 'rV1' ~ 'rV1',
                             group == 'rAIC' ~ 'rAIC',
                             group == 'noFeedback' ~ 'mental-rehearsal'))

ConnrV1_long_NFBruns%>%
  group_by(group, Run)%>% 
  summarise(mean = mean(nfb_value), sd = sd(nfb_value), n = n(), se = sd/sqrt(n))%>%
  ungroup(Run)%>%
  ggplot( aes(x=Run, y = mean, group= group, color=group))+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width =.1, position = pd)+
  geom_line(position = pd, size = .5)+
  geom_point(position = pd, size = 3)+
  geom_smooth(method = 'lm', se=FALSE, size = .1)+
  theme_minimal()+ scale_y_continuous(breaks=seq(-2,2,by=0.2))+
  labs( x="NFB run number", y="Contrast value rV1")+ scale_color_discrete(name = "Group")+theme_rtfmri2() 


expand_limits(y=c(-0.7, 1))+ scale_y_continuous(breaks=c( -0.7, -0.6,-0.5, -0.4, -0.3, -0.2, -0.1,0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))

ConnrV1_long_NFBruns_rAIC <- ConnrV1_long_NFBruns%>%
  subset(group == 'rAIC')

## linear regression rAIC group
ConnrV1_long_NFBruns_rAIC$Run <- as.numeric(ConnrV1_long_NFBruns_rAIC$Run)
modelrAICforrV1 <- lm(nfb_value ~ Run, ConnrV1_long_NFBruns_rAIC)
summary(modelrAICforrV1)


## linear regression rV1 group
ConnrV1_long_NFBruns_V1 <- ConnrV1_long_NFBruns%>%
  subset(group == 'rV1')
ConnrV1_long_NFBruns_V1$Run <- as.numeric(ConnrV1_long_NFBruns_V1$Run)
modelrV1V1 <- lm(nfb_value ~ Run, ConnrV1_long_NFBruns_V1)
summary(modelrV1V1)

## linear regression lV1 group
ConnrV1_long_NFBruns_lV1 <- ConnrV1_long_NFBruns%>%
  subset(group == 'lV1')
ConnrV1_long_NFBruns_lV1$Run <- as.numeric(ConnrV1_long_NFBruns_lV1$Run)
modelrV1lV1 <- lm(nfb_value ~ Run, ConnrV1_long_NFBruns_lV1)
summary(modelrV1lV1)

## linear regression NF group
ConnrV1_long_NFBruns_NF <- ConnrV1_long_NFBruns%>%
  subset(group.y == 'noFeedback')
ConnrV1_long_NFBruns_NF$Run <- as.numeric(ConnrV1_long_NFBruns_NF$Run)
modelNF <- lm(nfb_value ~ Run, ConnrV1_long_NFBruns_NF)
summary(modelNF)


### left V1


ConnlV1Session <- read_excel("C:/Users/Jeane/Documents/PhD/Analysis_Paper1/NFB/ConnlV1Session.xlsx")

ConnlV1Session$subid <- ConnlV1Session$id

ConnlV1_long <- ConnlV1Session%>%
  pivot_longer(cols = r0:r13, names_to = "run", values_to ="nfb_value")

ConnlV1_long_TB <- ConnlV1_long%>%
  subset(run == 'r0' | run == 'r13')

ConnlV1_long_TB<- ConnlV1_long_TB%>%
  mutate(group = case_when(group == 'lV1' ~ 'lV1',
                             group == 'rV1' ~ 'rV1',
                             group == 'rAIC' ~ 'rAIC',
                             group == 'noFeedback' ~ 'mental-rehearsal'))


## Two-way ANOVA

## Assumptions

# Outliers

Outliers <-ConnlV1_long_TB %>%
  group_by(run, group.y) %>%
  identify_outliers(nfb_value)

view(Outliers)

# Normality

ConnlV1_long_TB %>%
  group_by(run, group.y) %>%
  shapiro_test(nfb_value)

ggqqplot(ConnlV1_long_TB, "nfb_value", ggtheme = theme_bw()) +
  facet_grid(run ~ group.y)

## Homogenity

ConnlV1_long_TB %>%
  group_by(run) %>%
  levene_test(nfb_value ~ group.y)

box_m(ConnlV1_long_TB[, "nfb_value", drop = FALSE], ConnlV1_long_TB$group.y)

#ANOVA

res.aov <- anova_test(
  data = ConnlV1_long_TB, dv = nfb_value, wid = subid,
  between = group, within = run
)
get_anova_table(res.aov)

##Main effect of run

one.way2 <- ConnWholeV1_long_TB %>%
  group_by(group) %>%
  anova_test(dv = nfb_value, wid = subid, within = run) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way2

## Pairwise comparison

pwc2 <- ConnrAIC_long_TB %>%
  group_by(group) %>%
  pairwise_t_test(
    nfb_value ~ run, paired = TRUE, 
    p.adjust.method = "bonferroni"
  )  %>%
  select(-df, -statistic, -p)
pwc2

## Boxplot
ggplot(ConnlV1_long_TB, aes(x=group, y=nfb_value, color = run))+geom_boxplot(outlier.shape = NA)+
  labs( x="Group", y="Contrast value lV1")+geom_point(position=position_jitterdodge(0.2), aes(color = run))+
  scale_y_continuous(breaks=seq(-2,2,by=0.5))+
  scale_color_discrete(name = "Run",labels=c('Baseline Session 02', 'Transfer Session 03'))+  scale_x_discrete(limits = c("rAIC", "rV1", "lV1", "mental-rehearsal"))+
  theme_rtfmri2()

##Linear Regression just NFB runs

# rAIC group

ConnlV1_long_NFBruns <- ConnlV1_long%>%
  subset(run != 'r0' & run != 'r13' & run != 'r6' & run != 'r7')

ConnlV1_long_NFBruns[c('r', 'Run')] <- str_split_fixed(ConnlV1_long_NFBruns$run, '', 2)

ConnlV1_long_NFBruns$Run <- as.numeric(ConnlV1_long_NFBruns$Run)

ConnlV1_long_NFBruns['Run'][ConnlV1_long_NFBruns['Run'] == 8] <- 6
ConnlV1_long_NFBruns['Run'][ConnlV1_long_NFBruns['Run'] == 9] <- 7
ConnlV1_long_NFBruns['Run'][ConnlV1_long_NFBruns['Run'] == 10] <- 8
ConnlV1_long_NFBruns['Run'][ConnlV1_long_NFBruns['Run'] == 11] <- 9
ConnlV1_long_NFBruns['Run'][ConnlV1_long_NFBruns['Run'] == 12] <- 10

ConnlV1_long_NFBruns$Run<-factor(ConnlV1_long_NFBruns$Run, levels=c(1,2,3,4,5,6,7,8,9,10))

pd <- position_dodge(0.1)

ConnlV1_long_NFBruns<- ConnlV1_long_NFBruns%>%
  mutate(group = case_when(group == 'lV1' ~ 'lV1',
                             group == 'rV1' ~ 'rV1',
                             group == 'rAIC' ~ 'rAIC',
                             group == 'noFeedback' ~ 'mental-rehearsal'))

ConnlV1_long_NFBruns%>%
  group_by(group, Run)%>% 
  summarise(mean = mean(nfb_value), sd = sd(nfb_value), n = n(), se = sd/sqrt(n))%>%
  ungroup(Run)%>%
  ggplot( aes(x=Run, y = mean, group= group, color=group))+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width =.1, position = pd)+
  geom_line(position = pd, size = .5)+
  geom_point(position = pd, size = 3)+
  geom_smooth(method = 'lm', se=FALSE, size = .1)+
  theme_minimal()+ scale_y_continuous(breaks=seq(-2,2,by=0.2))+
  labs( x="NFB run number", y="Contrast value lV1")+ scale_color_discrete(name = "Group")+ theme_rtfmri2()


expand_limits(y=c(-0.7, 1))+ scale_y_continuous(breaks=c( -0.7, -0.6,-0.5, -0.4, -0.3, -0.2, -0.1,0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))

ConnlV1_long_NFBruns_rAIC <- ConnlV1_long_NFBruns%>%
  subset(group == 'rAIC')

## linear regression rAIC group
ConnlV1_long_NFBruns_rAIC$Run <- as.numeric(ConnlV1_long_NFBruns_rAIC$Run)
modelrAICforlV1 <- lm(nfb_value ~ Run, ConnlV1_long_NFBruns_rAIC)
summary(modelrAICforlV1)


## linear regression rV1 group
ConnlV1_long_NFBruns_rV1 <- ConnlV1_long_NFBruns%>%
  subset(group == 'rV1')
ConnlV1_long_NFBruns_rV1$Run <- as.numeric(ConnlV1_long_NFBruns_rV1$Run)
modellV1rV1 <- lm(nfb_value ~ Run, ConnlV1_long_NFBruns_rV1)
summary(modellV1rV1)

## linear regression lV1 group
ConnlV1_long_NFBruns_lV1 <- ConnlV1_long_NFBruns%>%
  subset(group == 'lV1')
ConnlV1_long_NFBruns_lV1$Run <- as.numeric(ConnlV1_long_NFBruns_lV1$Run)
modellV1lV1 <- lm(nfb_value ~ Run, ConnlV1_long_NFBruns_lV1)
summary(modellV1lV1)

## linear regression NF group
ConnlV1_long_NFBruns_NF <- ConnlV1_long_NFBruns%>%
  subset(group == 'no feedback')
ConnlV1_long_NFBruns_NF$Run <- as.numeric(ConnlV1_long_NFBruns_NF$Run)
modelNF <- lm(nfb_value ~ Run, ConnlV1_long_NFBruns_NF)
summary(modelNF)

