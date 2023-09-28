# Generate the visulization plot of the marignal mean for pcgs in CH study
# Sep 22nd;
library(readr)
library(tidyverse)
sample_data_PCGS_tutorial_paper_ICC01 <- read_csv("~/Downloads/sample_data_PCGS_tutorial_paper_ICC01.csv")
sample_data_PCGS_tutorial_paper_cp_mean = sample_data_PCGS_tutorial_paper_ICC01%>%
  group_by(clusters, periods)%>%summarise(means = mean(pcgs_integer), cp_size = n())%>%
  left_join(sample_data_PCGS_tutorial_paper_ICC01%>%select(clusters,periods,intervention_binary),by=c("clusters","periods"))
write.csv(sample_data_PCGS_tutorial_paper_cp_mean, "/Users/zying/Downloads/Tutorial paper CH pcgs simulated data codes and results v2/data/ClusterPeriodSize_CH_Simulated_PCGS.csv")
cohort.labs <- c("SNF 1", "SNF 2", "SNF 3", "SNF 4", "SNF 5", "SNF 6")
names(cohort.labs) <- c(1, 2, 3, 4, 5, 6)
options(scipen=10000)
ggplot(data = sample_data_PCGS_tutorial_paper_cp_mean, aes(periods, means, colour =factor( intervention_binary))) +
  geom_line() + geom_point(aes(size = cp_size)) + 
  scale_size_continuous(range = c(2, 5))+
  facet_wrap(~clusters,  ncol=3, labeller = labeller(clusters = cohort.labs)) +
  theme_bw() +
  scale_color_manual(values = c("grey69", "slategray2", "cornflowerblue"),
                     labels = c("Control", "Intervention", "Maintenance"),
                     breaks = c("0", "1", "2")) +
  xlab("Period") + ylab("Mean outcome") +
  scale_y_continuous(limits = c(16, 24), breaks = seq(16, 24, 0.5)) +
  scale_x_continuous(breaks = seq(1, 22, 1)) +
  theme(legend.position="bottom", legend.direction = "horizontal", legend.box = "vertical") +
  guides(color = guide_legend(title = "Phase"),
         size = guide_legend(title = "Cluster-period size"))

