# plot_model_selection.r
# Formerly: plotting_test_comparisons.r
# Created by John Trochta
# Date created:  11/06/2019
# Summary:
# This script runs the Bayesian ASA model for PWS using different covariates and 
# extracts the output
library(ggplot2)
library(dplyr)
library(R2admb)

setwd(here::here())
source(file=here::here("bin/data_reader.R"))
source(file=here::here("bin/data_header_reader.R"))

combine_effects <- function(folders,folder.names){
  modelPaths <- here::here(paste0("output/figs & tables/",folders,"_"))
  resultsPath <- here::here("output/figs & tables/")
  hyp_post <- read.csv(file=paste0(modelPaths[1],"effect posteriors.csv"),stringsAsFactors = FALSE)
  hyp_post <- bind_cols(run.set=rep(folder.names[1],nrow(hyp_post)),
                        hyp_post)
  if(length(folders)>1){
    for(i in 2:length(folders)){
      X <- read.csv(file=paste0(modelPaths[i],"effect posteriors.csv"),stringsAsFactors = FALSE)
      X <- bind_cols(run.set=rep(folder.names[i],nrow(X)),X)
      hyp_post <- rbind(hyp_post,X)
    }
  }
  return(hyp_post)
}

##############################################
# First plot estimated (fixed) effects from different time frames for SSB and R
folders <- c("2020_04_all fixed mortality effects",
             "2020_04_all fixed mortality effects_disease only",
             "2020_04_all fixed mortality effects & start from 1994",
             "2020_04_all fixed mortality effects & start from 2007",
             "2020_10_all fixed mortality effects & end 2009_2",
             "2020_04_all fixed mortality effects_disease only & start from 2007",
             "2020_10_all fixed mortality effects & start from 1994",
             "2020_10_all fixed mortality effects & start from 2007",
             "2020_10_all fixed mortality effects")
folder.names <- c("1980-2017","1994-2017","1994-2017","2007-2017","1980-2009","2007-2017","1994-2017","2007-2017","1980-2017")
# folder.names <- c("1980-2017","1994-2017","2007-2017","1980-2009")

hyp_post <- combine_effects(folders,folder.names)
hyp_post <- hyp_post %>% mutate(analysis.type="a) Mortality effects")
effects.df <- hyp_post

folders <- c("2020_04_all fixed age0 effects",
             "2020_04_all fixed age0 effects & start from 1994",
             "2020_04_all fixed age0 effects & start from 2007",
             "2020_10_all fixed age0 effects & end 2009_2")
folder.names <- c("1980-2017","1994-2017","2007-2017","1980-2009")

hyp_post <- combine_effects(folders,folder.names)
hyp_post <- hyp_post %>% mutate(analysis.type="b) Recruitment effects")
effects.df <- bind_rows(effects.df,hyp_post)

# REMOVE SPECIFIC COVARIATES:  Pacific Cod
effects.df.2 = filter(effects.df,variable!='GOA_Pacific_cod_SB_Ages_3.')
effects.df.2 = filter(effects.df.2,!(variable=='Summer_humpback_whale_estimates_Age_3.' & run.set!='1980-2009'))
effects.df.2 = filter(effects.df.2,!(variable=='PWS_Avg._Zooplankton_density_81_99'))
effects.df.2 = filter(effects.df.2,!(variable=='Disease_indices' & run.set=='2007-2017'))
effects.df.2 = filter(effects.df.2,!(variable=='Disease_indices_3' & run.set=='2007-2017'))
effects.df.2 = filter(effects.df.2,!(variable=='GOA_Arrowtooth_flounder_female_SB_Ages_3.'))


hypothesis_names <- c("I. hoferi pre-2007 (Ages 5+)",
                      "I. hoferi post-2007 (Ages 5+)",
                      "VHSV (Ages 3-4)",
                      "Summer upwelling (Ages 3+)",
                      "Summer PDO (Ages 3+)",
                      "Winter PDO (Ages 3+)",
                      'Summer NPGO (Ages 3+)',
                      'Winter NPGO (Ages 3+)',
                      'Total pink salmon returns (Ages 3+)',
                      'GOA walleye pollock SSB (Ages 3+)',
                      'GOA Pacific cod SB (Ages 3+)',
                      'GOA Arrowtooth flounder female SB (Ages 3+)',
                      'GOA Arrowtooth flounder total biomass (Ages 3+)',
                      "Summer humpback whale estimates (Age 3+)",
                      "Humpback whale counts (Age 3+)",
                      "Disease indices",
                      "Disease indices_2",
                      "Disease indices_3",
                      "Null",
                      "1989 regime shift",
                      "Hatchery-released juvenile pink salmon",
                      "GOA walleye pollock age 1",
                      'Summer NPGO',
                      'Summer PDO',
                      'Age 0 scale growth (lagged 1 year)',
                      'Freshwater discharge',
                      "PWS Avg. Zooplankton density '81-'99",
                      "Full")

hypothesis = matrix(NA,nrow=length(hypothesis_names),ncol=2)
hypothesis[,2] <- gsub(pattern=" |-|/",replacement="_",x=hypothesis_names)
hypothesis[,2] <- gsub(pattern="[()']",replacement="",x=hypothesis[,2])
hypothesis[,2] <- gsub(pattern="[+]",replacement=".",x=hypothesis[,2])
hypothesis[,1] = hypothesis_names

hypothesis[hypothesis[,2]=='1989_regime_shift',2] <- 'X1989_regime_shift'

# Re-assign names so they are appropriate for plotting
effects.df.2 <- effects.df.2 %>% 
  group_by(variable) %>%
  mutate(final.names = hypothesis[which(hypothesis[,2]%in%variable),1])

# Find "disease_indices_#" and rename to specific covariates
dis_var_names <- c("Disease_indices","Disease_indices_2","Disease_indices_3")
fin_dis_names <- c("I. hofer pre-2007 (Ages 5+)","I. hofer post-2007 (Ages 5+)","VHSV (Ages 3-4)")

effects.df.2 <- effects.df.2 %>% group_by(variable) %>%
  mutate(final.names = ifelse(variable%in%dis_var_names,
  fin_dis_names[which(dis_var_names%in%variable)],
  final.names))


# Assign numeric ID values to each covariate
effects.df.2 <- effects.df.2 %>% 
  group_by(analysis.type) %>% 
  mutate(id = as.numeric(as.factor(final.names))-1)

plot_labels <- ungroup(effects.df.2) %>% group_by(analysis.type,final.names) %>% 
  summarize(id=unique(id))


library(reshape)
library(gridExtra)
library(ggridges)
# hyp_post <- data.frame(hyp_post,
#                        id=as.integer(factor(hyp_post$covariates)))
# colors <- c( "firebrick1",  "#FFBA00FF","darkmagenta","#13F24AFF", "blue", "#00A4FFFF" )
font.size <- 14
plot.split <- c("a) Mortality effects","b) Recruitment effects")

sub.labels <- c("a)", "b)")

resultsPath <- here::here("output/figs & tables/")

pl.1 <- plot_labels
df.1 <- effects.df.2

ggplot(df.1,aes(x=value,y=final.names,height=..density..)) + 
  geom_density_ridges(stat="density",aes(fill=run.set),scale=1,alpha=0.4,trim=TRUE) +
  #geom_point(size=2.5,color="black",fill="white", shape=21,position=position_dodge(width=0.75)) +
  facet_wrap(.~analysis.type,scales="free",strip.position = "top",ncol=2)+
  labs(y=NULL,x=expression("Estimated effect ("~beta*")"))+
  scale_fill_manual(values=c("2007-2017"="#800000",
                              "1994-2017"="#FFBA00FF", 
                              "1980-2017"="#00BFFF",             #"#63ACBE",
                             "1980-2009" = "#191970"))+
  #geom_point(size=0,aes(x=0,y=0))+
  #theme_ridges()+
  theme_classic()+
  scale_y_discrete(expand = c(0.01,0)) +
  geom_vline(xintercept=0, 
             lwd=0.4, colour="grey30")+
  theme(strip.background = element_blank(),#element_rect(fill="white", colour=NULL),
        strip.text.x = element_text(size=font.size,face="bold"),
        #strip.placement = "none",#"outside",
        axis.text.y = element_text(size=font.size),
        axis.text.x = element_text(size=font.size),
        axis.title = element_text(size=font.size,face="bold"),
        legend.position = "none",#c(0.15, 0.9),
        legend.title=element_blank(),
        legend.text=element_blank(),#element_text(size=font.size),
        legend.key=element_blank(),
        legend.background = element_blank())+
        #panel.background = element_rect(color = "black", fill=NA, size=1,linetype="solid")) +
guides(colour = guide_legend(override.aes = list(size=2,linetype=0,fill=NA),reverse=TRUE)) +
    theme(legend.position = c(0.4, 0.9),
          legend.title=element_blank(),
          legend.text=element_text(size=font.size),
          legend.key=element_blank(),
          legend.background = element_blank())

#pdf(file=paste0(resultsPath,"/Figure_Comparing effects from different time frames.pdf"),width=14,height=6)
ggsave(filename=paste0(resultsPath,"/Figure_Comparing effect posteriors.png"),
       width=14, height=6, units="in",dpi=600)


# Alternative way of comparing - speading out time frames
df.mor <- filter(df.1,analysis.type=="a) Mortality effects")
df.mor$run.set <- factor(df.mor$run.set,levels= c("1980-2017","1994-2017","2007-2017","1980-2009"))
ggplot(df.mor,aes(x=value,y=final.names,height=..density..)) + 
  geom_density_ridges(stat="density",aes(fill=run.set),scale=1,alpha=0.4,trim=FALSE,size=0.1) +
  facet_wrap(.~run.set,strip.position = "top",nrow=1)+
  labs(y=NULL,x=expression("Estimated effect ("~beta*")"))+
  scale_fill_manual(values=c("2007-2017"="#800000",
                             "1994-2017"="#FFBA00FF", 
                             "1980-2017"="#00BFFF",     
                             "1980-2009" = "#191970"))+
  theme_classic()+
  scale_y_discrete(expand = c(0.01,0)) +
  geom_vline(xintercept=0, 
             lwd=0.4, colour="grey30",linetype="dashed")+
  theme(strip.background = element_blank(),#element_rect(fill="white", colour=NULL),
        strip.text.x = element_text(size=font.size),
        #strip.placement = "none",#"outside",
        axis.text.y = element_text(size=font.size),
        axis.text.x = element_text(size=font.size-4),
        axis.title = element_text(size=font.size,face="bold"),
        legend.position = "none",#c(0.15, 0.9),
        legend.title=element_blank(),
        legend.text=element_blank(),#element_text(size=font.size),
        legend.key=element_blank(),
        legend.background = element_blank())

ggsave(filename=paste0(resultsPath,"/Figure_Comparing mortality effect posteriors.png"),
       width=11, height=6, units="in",dpi=600)



df.rec <- filter(df.1,analysis.type=="b) Recruitment effects")
df.rec$run.set <- factor(df.rec$run.set,levels= c("1980-2017","1994-2017","2007-2017","1980-2009"))
ggplot(df.rec,aes(x=value,y=final.names,height=..density..)) + 
  geom_density_ridges(stat="density",aes(fill=run.set),scale=1,alpha=0.4,trim=FALSE,size=0.1) +
  #geom_point(size=2.5,color="black",fill="white", shape=21,position=position_dodge(width=0.75)) +
  facet_wrap(.~run.set,strip.position = "top",nrow=1)+
  labs(y=NULL,x=expression("Estimated effect ("~beta*")"))+
  scale_fill_manual(values=c("2007-2017"="#800000",
                             "1994-2017"="#FFBA00FF", 
                             "1980-2017"="#00BFFF",     
                             "1980-2009" = "#191970"))+
  #geom_point(size=0,aes(x=0,y=0))+
  #theme_ridges()+
  theme_classic()+
  scale_y_discrete(expand = c(0.01,0)) +
  geom_vline(xintercept=0, 
             lwd=0.4, colour="grey30",linetype="dashed")+
  theme(strip.background = element_blank(),#element_rect(fill="white", colour=NULL),
        strip.text.x = element_text(size=font.size),
        #strip.placement = "none",#"outside",
        axis.text.y = element_text(size=font.size),
        axis.text.x = element_text(size=font.size-4),
        axis.title = element_text(size=font.size,face="bold"),
        legend.position = "none",#c(0.15, 0.9),
        legend.title=element_blank(),
        legend.text=element_blank(),#element_text(size=font.size),
        legend.key=element_blank(),
        legend.background = element_blank())

#pdf(file=paste0(resultsPath,"/Figure_Comparing effects from different time frames.pdf"),width=14,height=6)
ggsave(filename=paste0(resultsPath,"/Figure_Comparing recruitment effect posteriors.png"),
       width=10, height=4, units="in",dpi=600)
