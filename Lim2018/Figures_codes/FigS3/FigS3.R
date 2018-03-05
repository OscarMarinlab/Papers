library (plotrix) 
library (ggplot2) 
library (plyr)

DATA <-read.csv("MZsvZdapi_ratio.csv", header = TRUE)
d1<-subset(DATA, cellmarker=="Tbr2")
d2<-subset(DATA, cellmarker == "Tdt")
p1<- ggplot(d1, aes(y= ratiotodapi, x = area)) +
  #facet_grid( . ~ area , scales="free_y") +
  geom_boxplot(color = "chartreuse4", fill = "chartreuse4", alpha=0.2) + 
  geom_jitter() +
  ylab("Fraction normalized to Dapi")+
  theme_classic() 

p2<- ggplot(d2, aes(y= ratiotodapi, x = area)) +
  #facet_grid( . ~ area , scales="free_y") +
  geom_boxplot(color = "red", fill = "red", alpha=0.2) + 
  geom_jitter() +
  ylab("Fraction normalized to Dapi")+
  theme_classic() 