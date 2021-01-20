
################################################################
# SIMULATIONS RESULTS - MDA Research
# Marta Bofill and Guadalupe Gómez
################################################################

# DESCRIPTION
# This code creates the set of scenarios for the simulation. In this case, we assume  higher
# numbers of patients assigned to active treatment. In particular, we consider 2:1 ratio trials.

rm(list = ls())
load("C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survmixer/code_PAPER/Additional-Simulations/results_sim/RESULTS_sim.RData")

library(ggplot2)
library(gridExtra)
library(ggpubr)
library(hrbrthemes)
library(viridis)

data$cases = 4
for(i in 1:dim(data)[1]){
  if(data$Delta_r[i]==0 & data$Delta_nr[i]==0) data$cases[i]= 1
  if(data$Delta_r[i]==0 & data$Delta_nr[i]!=0) data$cases[i]= 2
  if(data$Delta_r[i]!=0 & data$Delta_nr[i]==0) data$cases[i]=3
}
data$cases = as.factor(data$cases)
summary(data$cases)

data$diff_power = data$Test_Reject - data$Test_Reject_LR
data$diff_alpha = data$Test_Reject_size - data$Test_Reject_LR_size


summary(data$os_samplesize)
data = subset(data, data$os_samplesize>100 & data$os_samplesize<5000)
summary(data$os_samplesize)
summary(data$os_effect)


############
# Boxplots alpha and power

windows(height = 14, width = 14)

p1 <- ggplot(data, aes(x=cases, y=Test_Reject,  color=cases)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="grey50", size=0.4, alpha=1)  + labs(y = "Power RMST test", x ="Settings", color ="Settings")+ coord_cartesian(ylim = c(0.65, 1))#+ ylim(0.65, 0.9)
p2 <- ggplot(data, aes(x=cases, y=Test_Reject_LR,  color=cases)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="grey50", size=0.4, alpha=1) + labs(y = "Power logrank test", x ="Settings", color ="Settings") + coord_cartesian(ylim = c(0.65, 1))#+ ylim(0.65, 0.9)
diff_p12 <- ggplot(data, aes(x=cases, y=diff_power,  color=cases)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="grey50", size=0.4, alpha=1) + labs(y = "Difference Power (RMST - logrank) ", x ="Settings", color ="Settings") #+ coord_cartesian(ylim = c(-0.1, 0.1))

p3 <- ggplot(data, aes(x=cases, y=Test_Reject_size, color=cases)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="grey50", size=0.4, alpha=1) + labs(y = "Significance level RMST test", x ="Settings", color ="Settings")  + coord_cartesian(ylim = c(0.035, 0.08)) # + ylim(0.035, 0.09)
p4 <- ggplot(data, aes(x=cases, y=Test_Reject_LR_size, color=cases)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="grey50", size=0.4, alpha=1)+ labs(y = "Significance level logrank test", x ="Settings", color ="Settings") + coord_cartesian(ylim = c(0.035, 0.08)) # + ylim(0.035, 0.09)
diff_p34 <- ggplot(data, aes(x=cases, y=diff_alpha,  color=cases)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="grey50", size=0.4, alpha=1) + labs(y = "Difference Significance level (RMST - logrank )", x ="Settings", color ="Settings")+ coord_cartesian(ylim = c(-0.03, 0.03))


figure <- ggarrange(p1,p2,diff_p12,p3,p4,diff_p34, ncol=3, nrow=2, common.legend = TRUE, legend="bottom")
annotate_figure(figure,
                top = text_grob("Power and Significance level",
                # top = text_grob(expression(paste("Power and Significance level (", bshape^{(0)}, "=", bshape^{(1)}, ")")),
                face = "bold", size = 14))

#

windows(height = 14, width = 12)
figure_mod <- ggarrange(p1,p2,p3,p4, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")
annotate_figure(figure_mod,
                top = text_grob("Power and Significance level",
                                face = "bold", size = 14))

############
# Sample size and effect size

windows(height = 14, width = 14)

p5 <- ggplot(data, aes(x=cases, y=os_samplesize,  color=cases)) +
  geom_boxplot(outlier.shape = NA) + labs(y = "Sample size", x ="Settings", color ="Settings") #+ scale_y_continuous(limits = quantile(data$os_samplesize, c(0.1, 0.9)))
p5_bis <-  ggplot(data, aes(x=NA., y=os_samplesize, shape=cases, color=cases)) + geom_point(size=2)+ labs(y = "Sample size",x = "Scenarios") +
  theme(legend.position = "none",  axis.text.x = element_blank())
p6 <- ggplot(data, aes(x=cases, y=os_effect,  color=cases)) +
  geom_boxplot() + labs(y = "Effect size (RMST difference)", x ="Settings", color ="Settings")
p6_bis <- ggplot(data, aes(x=NA., y=os_effect, shape=cases, color=cases)) + geom_point(size=2)+ labs(y = "Effect size",x = "Scenarios")+
  theme(legend.position = "none",  axis.text.x = element_blank())

figure2 <- ggarrange(p5,p5_bis,p6,p6_bis, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")
annotate_figure(figure2, top = text_grob("Sample size and Effect size",
                                         # expression(paste("Sample size and Effect size (", bshape^{(0)}, "=", bshape^{(1)}, ")")),
                                        face = "bold", size = 14))

windows(height = 7, width = 10)

figure2_bis <- ggarrange(p5,p6, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")
annotate_figure(figure2_bis, top = text_grob("Sample size and Effect size",
                                         # expression(paste("Sample size and Effect size (", bshape^{(0)}, "=", bshape^{(1)}, ")")),
                                         face = "bold", size = 14))


summary(subset(data,data$cases==1)$os_effect);summary(subset(data,data$cases==1)$os_samplesize)
summary(subset(data,data$cases==2)$os_effect);summary(subset(data,data$cases==2)$os_samplesize)
summary(subset(data,data$cases==3)$os_effect);summary(subset(data,data$cases==3)$os_samplesize)
summary(subset(data,data$cases==4)$os_effect);summary(subset(data,data$cases==4)$os_samplesize)

############
# Scatterplots alpha and power

windows(height = 14, width = 14)
par(mfrow = c(2, 2))

p1 <- ggplot(data, aes(x=NA., y=Test_Reject, shape=cases, color=cases)) +
  geom_point(size=2)+ ylim(0.65, 0.9) + labs(y = "Power RMST test")
p2 <- ggplot(data, aes(x=NA., y=Test_Reject_LR, shape=cases, color=cases)) +
  geom_point(size=2)+ ylim(0.65, 0.9) + labs(y = "Power logrank test")
p3 <- ggplot(data, aes(x=NA., y=Test_Reject_size  , shape=cases, color=cases)) +
  geom_point(size=2)+ ylim(0.035, 0.09) + labs(y = "Significance level RMST test")
p4 <- ggplot(data, aes(x=NA., y=Test_Reject_LR_size, shape=cases, color=cases)) +
  geom_point(size=2)+ ylim(0.035, 0.09) + labs(y = "Significance level logrank test")

grid.arrange(p1,p2,p3,p4, ncol=2)

################################################################
################################################################
