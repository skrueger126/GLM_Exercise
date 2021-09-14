#R code for Kuhlman et al. 2021. Relative bee abundance varies by collection method and flower richness:
# implications for understanding patterns in bee community data. Ecological Evidence and Solutions.
# Contact Phil Hahn (hahnp@ufl.edu) for questions about the code.

library(emmeans)
library(car)
library(lme4)
library(lmerTest)
library(patchwork)
library(RVAideMemoire)
library(DHARMa)
library(tidyverse)

bf <- read_csv('MK_ESE_BeeFlower_Phenology_20142017.csv')

### FIGURE 1 - Plot for flower richness, bee abundance, and bee richness across sampling dates
MK_ESE_Fig1 <- ggplot(data = bf, aes(x = WOY)) + # DOY on x-axis, number on y-axis
  geom_vline(xintercept=27.5, linetype=2, lwd=2, color="grey")+
  geom_smooth(data = bf, aes(x = WOY, y = Flower_Richness, color="Flower_Richness"), span=.5, size=1.5, fill="violet") +
  geom_smooth(data = bf, aes(x = WOY, y = bee_richness, color="bee_richness"), span=.5, size=1.5, fill="grey") +
  geom_smooth(data = bf, aes(x = WOY, y = bee_abund/100, color="bee_abund"), span=.5, size=1.5, fill="goldenrod2")+
  geom_point(data = bf, aes(y = Flower_Richness, shape=as_factor(Year)), color="violetred2", size=2) +
  geom_point(data = bf, aes(y = bee_richness, shape=as_factor(Year)), color="black", size=2) +
  geom_point(data = bf, aes(y = bee_abund/100, shape=as_factor(Year)), color="goldenrod4", size=2) +
  scale_x_continuous(name = "Sampling date", breaks = c(18,22,26,31,36), labels=c("May","June","July","August","Sept")) +
  scale_y_continuous(name = "Flower or bee richness", limits=c(-5,95), sec.axis = sec_axis(~.*100, name="Bee abundance"))+
  #scale_color_manual(values=c("violet"), labels=c("Flower richness")) +
  scale_color_manual(values=c("goldenrod2","grey","violet"), labels=c("Bee abund","Bee richness","Flower richness")) +
  scale_shape_manual(values=c(15,16,17,18)) +
  theme(axis.line = element_line(size = 1, colour = "black"))+
  theme(panel.background = element_blank())+
  theme(panel.border = element_rect(color="black", fill=NA, size=1.5)) +
  theme(text = element_text(size=18))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.ticks.length=unit(0.2, "cm"),  
        axis.text.x = element_text(margin=margin(5,5,5,5,"pt"),colour="black"),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt"),colour="black"))+
  theme(plot.title = element_text(hjust = 0.5,size = 10))+
  theme(strip.background = element_rect(colour=NA, fill=NA))+
  theme(strip.text = element_text(size=18))+
  theme(legend.position = "top", legend.title=element_blank(),
        legend.background = element_rect(fill="white", size=1, linetype="solid", color="black"))+
  guides(color=guide_legend(override.aes = list(fill=NA)))

tiff("MK_ESE_Fig1.tiff", width=10, height=6, units='in',
     res=600, compression='lzw')
MK_ESE_Fig1
dev.off()

####################################################################################################################
####### CORRELATIONS for Fig. 2 ###############################################################################################

### BEE ABUNDANCE - PANS CORRELATING WITH FLORAL RICHNESS ###
# all season
cP <- cor.test(bf$bee_abund, bf$Flower_Richness, method = "pearson")
cP   

# early flowering 
cE <- cor.test(bf$bee_abund[bf$WOY < 28], bf$Flower_Richness[bf$WOY < 28],
               method = "pearson")
cE

# late season 
cL <- cor.test(bf$bee_abund[bf$WOY >= 28], bf$Flower_Richness[bf$WOY >= 28],
               method = "pearson")
cL

### BEE RICHNESS - PANS CORRELATING WITH FLORAL RICHNESS ###
# all season
cR <- cor.test(bf$bee_richness, bf$Flower_Richness, method = "pearson")
cR  

# early season
cRE <- cor.test(bf$bee_richness[bf$WOY < 28], bf$Flower_Richness[bf$WOY < 28], method = "pearson")
cRE   

# late season
cRL <- cor.test(bf$bee_richness[bf$WOY >= 28], bf$Flower_Richness[bf$WOY >= 28], method = "pearson")
cRL


fig2atext <- data.frame(label=c("A","B"), Season=c('Early','Late'), x=c(17.5,37),y=7000)
fig2atextr <- data.frame(label=c("r = -0.80","r = 0.12"), Season=c('Early','Late'), x=c(25,30),y=1200)
fig2atextp <- data.frame(label=c("p < 0.001","p = 0.65"), Season=c('Early','Late'), x=c(25,30),y=700)

fig2a <- ggplot(data = bf, aes(x = Flower_Richness, y = bee_abund)) +
  geom_smooth(data = bf %>% filter(Season=="Early"), 
              aes(x = Flower_Richness, y = bee_abund), method='lm', se=T, color='black') +
  geom_point(aes(shape=as_factor(Year)), size=3) + facet_wrap(~Season, scales='free_x')+
  scale_y_continuous(name = "Bee abundance")+xlab("Floral richness")+
  scale_shape_manual(values=c(15,16,17,18)) +
  theme(axis.line = element_line(size = 1, colour = "black"))+
  theme(panel.background = element_blank())+
  theme(panel.border = element_rect(color="black", fill=NA, size=1.5)) +
  theme(text = element_text(size=18))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.ticks.length=unit(0.2, "cm"),  
        axis.text.x = element_text(margin=margin(5,5,5,5,"pt"),colour="black"),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt"),colour="black"))+
  theme(plot.title = element_text(hjust = 0.5,size = 10))+
  theme(strip.background = element_rect(colour=NA, fill=NA))+
  theme(strip.text = element_text(size=18))+
  theme(legend.position = c(.59,.85), legend.title=element_blank(),
        legend.background = element_rect(fill="white", size=1.5, linetype="solid", color="black"))+
  geom_text(data=fig2atext, aes(x=x, y=y, label=label), size=8)+
  geom_text(data=fig2atextr, aes(x=x, y=y, label=label), size=6) +
  geom_text(data=fig2atextp, aes(x=x, y=y, label=label), size=6)


fig2btext <- data.frame(label=c("C","D"), Season=c('Early','Late'), x=c(17.5,37),y=95)
fig2btextr <- data.frame(label=c("r = -0.20","r = 0.68"), Season=c('Early','Late'), x=c(25,35),y=45)
fig2btextp <- data.frame(label=c("p = 0.46","p < 0.001"), Season=c('Early','Late'), x=c(25,35),y=40)

fig2b <- ggplot(data = bf, aes(x = Flower_Richness, y = bee_richness)) +
  geom_smooth(data = bf %>% filter(Season=="Late"), 
              aes(x = Flower_Richness, y = bee_richness), method='lm', se=T, color='black') +
  geom_point(aes(shape=as_factor(Year)), size=3) + facet_wrap(~Season, scales='free_x')+
  scale_y_continuous(name = "Bee richness")+xlab("Floral richness")+
  scale_shape_manual(values=c(15,16,17,18)) +
  theme(axis.line = element_line(size = 1, colour = "black"))+
  theme(panel.background = element_blank())+
  theme(panel.border = element_rect(color="black", fill=NA, size=1.5)) +
  theme(text = element_text(size=18))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.ticks.length=unit(0.2, "cm"),  
        axis.text.x = element_text(margin=margin(5,5,5,5,"pt"),colour="black"),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt"),colour="black"))+
  theme(plot.title = element_text(hjust = 0.5,size = 10))+
  theme(strip.background = element_rect(colour=NA, fill=NA))+
  theme(strip.text = element_text(size=18))+
  theme(legend.position = c(.59,.85), legend.title=element_blank(),
        legend.background = element_rect(fill="white", size=1.5, linetype="solid", color="black"))+
  geom_text(data=fig2btext, aes(x=x, y=y, label=label), size=8) +
  geom_text(data=fig2btextr, aes(x=x, y=y, label=label), size=6) +
  geom_text(data=fig2btextp, aes(x=x, y=y, label=label), size=6)


tiff("MK_ESE_Fig2.tiff", width=10, height=10, units='in',
     res=600, compression='lzw')
(fig2a/fig2b)
dev.off()


#####################################################################################
### correlate Pan and Net abundances
# netted bees and flower richness data set
b3a <- read_csv('MK_ESE_BeePanNet_Abundances.csv')

cPnN_rich <- cor.test(b3a$bee_richness_pt,b3a$bee_richness_n, method = "pearson")
cPnN_rich

fig3r <- data.frame(label=c("r = -0.20","r = 0.68"), Season=c('Early','Late'), x=c(25,35),y=45)
fig3p <- data.frame(label=c("p = 0.46","p < 0.001"), Season=c('Early','Late'), x=c(25,35),y=40)

MK_ESE_Fig3 <- ggplot(data= b3a %>% filter(Year!=2014)) + 
  geom_smooth(data = b3a %>% filter(Year!=2014), 
              aes(x = bee_richness_pt, y = bee_richness_n), method='lm', se=T, color='black')+
  geom_point(aes(x=bee_richness_pt, y=bee_richness_n, shape=as.factor(Year)), size=3)+
  theme(axis.line = element_line(size = 1, colour = "black"))+
  xlab('Bee richness in pans') + ylab('Bee richness in nets')+
  scale_shape_manual(values=c(16,17,18)) +
  annotate(geom="text", x=60, y=8, label="r = 0.58")+
  annotate(geom="text", x=60, y=6, label="p = 0.004")+
  theme(panel.background = element_blank())+
  theme(panel.border = element_rect(color="black", fill=NA, size=1.5)) +
  theme(text = element_text(size=18))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.ticks.length=unit(0.2, "cm"),  
        axis.text.x = element_text(margin=margin(5,5,5,5,"pt"),colour="black"),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt"),colour="black"))+
  theme(plot.title = element_text(hjust = 0.5,size = 10))+
  theme(strip.background = element_rect(colour=NA, fill=NA))+
  theme(strip.text = element_text(size=18))+
  theme(legend.position = c(.9,.15), legend.title=element_blank(),
        legend.background = element_rect(fill="white", size=1.5, linetype="solid", color="black")) 

tiff("MK_ESE_Fig3.tiff", width=6, height=4, units='in',
     res=600, compression='lzw')
MK_ESE_Fig3
dev.off()

###############################################
################################################################
### PAN TRAP VS NETTING BY GENUS FIGURE #bfl
bfl <- read_csv('MK_ESE_BeeAbundFloral_Early20152017.csv')

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

## Fig 4 abundance figures
MK_ESE_Fig4a <- ggplot(data = bfl , 
       aes(x = Flower_Richness, y = log(bee_abund+1), color=CollMeth)) + # DOY on x-axis, number on y-axis
  geom_smooth(data = bfl , aes(x = Flower_Richness, y = log(bee_abund+1), group=CollMeth, linetype=CollMeth), 
              method='lm')+
  geom_point(data = bfl , aes(x = Flower_Richness, y = log(bee_abund+1), shape=CollMeth),size=.8,alpha=.1) +
  facet_wrap(~GenusName)+ aes(color=GenusName) +
  ggtitle("A")+
  scale_shape_manual(values=c(15,19), labels=c("Nets","Pans"))+scale_linetype_manual(values=c("dashed","solid"), labels=c("Nets","Pans"))+
  #scale_fill_manual(values=21)+
  scale_color_manual(values=cbPalette)+
  #scale_color_manual(values=c("black","blue"), labels=c("Nets","Pans"))+
  scale_y_continuous(name = "log(Bee abundance +1)")+scale_x_continuous(name = "Floral richness")+
  theme(axis.line = element_line(size = 1, colour = "black"))+
  theme(panel.background = element_blank())+
  theme(panel.border = element_rect(color="black", fill=NA, size=1.5)) +
  theme(text = element_text(size=18))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.ticks.length=unit(0.2, "cm"),  
        axis.text.x = element_text(margin=margin(5,5,5,5,"pt"),colour="black"),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt"),colour="black"),
        plot.title = element_text(hjust = 0,vjust=-8,size = 20))+
  theme(strip.background = element_rect(colour=NA, fill=NA))+
  theme(strip.text = element_text(size=18))+
  theme(legend.position = "top", legend.title=element_blank(),
        legend.background = element_rect(fill="white", size=1.5, linetype="solid", color="black"))+
  guides(color=F,linetype=F, shape=guide_legend(override.aes = list(size=3, alpha=1)))

MK_ESE_Fig4b <- ggplot(data = bfl , 
               aes(x = Flower_Richness, y = pres, color=CollMeth)) + # DOY on x-axis, number on y-axis
  geom_smooth(data = bfl , aes(x = Flower_Richness, y = pres, group=CollMeth, linetype=CollMeth), 
              method='glm', method.args=list(family="binomial"))+
  geom_point(data = bfl %>% filter(CollMeth=='pt'), 
              aes(x = Flower_Richness, y = pres+.05, shape=CollMeth),size=1,alpha=.1) +
  geom_point(data = bfl %>% filter(CollMeth=='n'), 
             aes(x = Flower_Richness, y = pres-.05, shape=CollMeth),size=1,alpha=.1) +
  facet_wrap(~GenusName)+ aes(color=GenusName) +
  ggtitle("B")+
  scale_shape_manual(values=c(15,19), labels=c("Nets","Pans"))+scale_linetype_manual(values=c("dashed","solid"), labels=c("Nets","Pans"))+
  #scale_fill_manual(values=21)+
  scale_color_manual(values=cbPalette)+
  #scale_color_manual(values=c("black","blue"), labels=c("Nets","Pans"))+
  scale_y_continuous(name = "Presence", limits=c(-.2,1.15),breaks=c(0,.5,1))+scale_x_continuous(name = "Floral richness")+
  theme(axis.line = element_line(size = 1, colour = "black"))+
  theme(panel.background = element_blank())+
  theme(panel.border = element_rect(color="black", fill=NA, size=1.5)) +
  theme(text = element_text(size=18))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.ticks.length=unit(0.2, "cm"),  
        axis.text.x = element_text(margin=margin(5,5,5,5,"pt"),colour="black"),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt"),colour="black"),
        plot.title = element_text(hjust = 0,vjust=-8,size = 20))+
  theme(strip.background = element_rect(colour=NA, fill=NA))+
  theme(strip.text = element_text(size=18))+
  theme(legend.position = "top", legend.title=element_blank(),
        legend.background = element_rect(fill="white", size=1.5, linetype="solid", color="black"))+
  theme(legend.key.width = unit(2,"cm"))+
  guides(color=F,shape=guide_legend(override.aes = list(size=1.5,alpha=1)), 
         linetype=guide_legend(override.aes=list(color="black")))

tiff("MK_ESE_Fig4.tiff", width=12, height=8, units='in',
     res=600, compression='lzw')
MK_ESE_Fig4a+MK_ESE_Fig4b
dev.off()

###################################################################################################################
### generalized linear model for bee abund ########################################################################
glm1 <- glmer(bee_abund~CollMeth*GenusName*scale(Flower_Richness)+(1|Location_Name/CollMeth/GenusName/Year)+(1|obs), 
              family="poisson", data=bfl , nAGQ = 0, 
              control = glmerControl(calc.derivs = FALSE,optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
summary(glm1)
Anova(glm1, type=2)

overdisp.glmer(glm1) # overdispersion ratio calculator from RVAideMemoire

sim_glm1 <- simulateResiduals(fittedModel = glm1, n = 250) ## simulate resids from DHARMa package
plot(sim_glm1) # residuals look ok, deviations sensative to sample size?

## extract emmeans for plotting
glm1mm <- as.data.frame(emtrends(glm1, pairwise~CollMeth|GenusName, var="Flower_Richness", adjust="tukey"))
emtrends(glm1, pairwise~CollMeth, var="Flower_Richness", adjust="tukey", transform='response')
tableS2 <- as.data.frame(emmeans(glm1, ~CollMeth|GenusName, transform='response'))

gendat <- as.data.frame(emmeans(glm1, ~CollMeth|GenusName, transform='response', at=list(Flower_Richness=1)))
gendat$div <- "Low"

gendat2 <- as.data.frame(emmeans(glm1, ~CollMeth|GenusName, transform='response', at=list(Flower_Richness=15)))
gendat2$div <- "High"

gendat_full <- rbind(gendat,gendat2)
gendat_full$CollMeth <- sub("n","Nets", gendat_full$CollMeth)
gendat_full$CollMeth <- sub("pt","Pans", gendat_full$CollMeth)
gendat_full$div <- factor(gendat_full$div, levels=c('Low','High'))

## construct figure 5a
MK_ESE_Fig5a <- ggplot(gendat_full, aes(fill=GenusName, x=div, y=rate)) + geom_bar(position="fill",stat="identity")+
  facet_wrap(~CollMeth) + xlab('Floral richness')+ ylab('Relative abundance')+
  scale_fill_manual()+
  ggtitle("A")+
  theme(panel.background = element_blank(), plot.margin=margin(.5,.5,.5,.5, "cm"))+
  theme(panel.border = element_rect(color="black", fill=NA, size=1.5)) +
  theme(text = element_text(size=18))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.ticks.length=unit(0.2, "cm"),  
        axis.text.x = element_text(margin=margin(5,5,5,5,"pt"),colour="black"),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt"),colour="black"),
        plot.title = element_text(hjust = 0,size = 20))+
  theme(strip.background = element_rect(colour=NA, fill=NA))+
  theme(strip.text = element_text(size=18)) +
  theme(legend.position = "none")
  #theme(legend.position = "top", legend.title=element_blank(),
   #     legend.background = element_rect(fill="white", size=1.5, linetype="solid", color="black"))


##########################################################################################################
## glmm binomial model ####################################################################################

glm2 <- glmer(pres~CollMeth*GenusName*scale(Flower_Richness)+(1|Location_Name/CollMeth/GenusName/Year), 
              family="binomial", data=bfl, nAGQ = 0,
              control = glmerControl(calc.derivs = FALSE,optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
summary(glm2)
Anova(glm2, type=2)

overdisp.glmer(glm2) # overdispersion ratio calculator from RVAideMemoire

sim_glm2 <- simulateResiduals(fittedModel = glm2, n = 250)
plot(sim_glm2) # residuals look good

## extract emmeans for plotting
gendatl <- as.data.frame(emmeans(glm2, ~CollMeth|GenusName, transform='response', at=list(Flower_Richness=1)))
gendatl$div <- "Low"

gendatl2 <- as.data.frame(emmeans(glm2, ~CollMeth|GenusName, transform='response', at=list(Flower_Richness=15)))
gendatl2$div <- "High"

gendatl_full <- rbind(gendatl,gendatl2)
gendatl_full$CollMeth <- sub("n","Nets", gendatl_full$CollMeth)
gendatl_full$CollMeth <- sub("pt","Pans", gendatl_full$CollMeth)
gendatl_full$div <- factor(gendatl_full$div, levels=c('Low','High'))

## construct figure 5b
MK_ESE_Fig5b <- ggplot(gendatl_full, aes(fill=GenusName, x=div, y=prob)) + geom_bar(position="fill",stat="identity")+
  facet_wrap(~CollMeth) + xlab('Floral richness')+ ylab('Relative presence')+
  ggtitle("B")+
  theme(panel.background = element_blank(), plot.margin=margin(.5,.5,.5,.5, "cm"))+
  theme(panel.border = element_rect(color="black", fill=NA, size=1.5)) +
  theme(text = element_text(size=18))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.ticks.length=unit(0.2, "cm"),  
        axis.text.x = element_text(margin=margin(5,5,5,5,"pt"),colour="black"),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt"),colour="black"),
        plot.title = element_text(hjust = 0,size = 20))+
  theme(strip.background = element_rect(colour=NA, fill=NA))+
  theme(strip.text = element_text(size=18))+
  theme(legend.position = "top", legend.title=element_blank(),
    legend.background = element_rect(fill="white", size=1.5, linetype="solid", color="black"))


tiff("MK_ESE_Fig5.tiff", width=10, height=15, units='in',
       res=600, compression='lzw')
MK_ESE_Fig5a/MK_ESE_Fig5b
dev.off()





