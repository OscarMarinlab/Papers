
library(plyr)
library(ggplot2)
library (plotrix) #use for std.error
C1 <-read.csv("2016-03-28-trackcomplie.csv", header = TRUE, sep = ",")
t0 <- subset(C1, T.min ==0)
tsub<- t0[, c(1:3, 5, 6)]
ct0<-colnames(tsub)
ct0[2:3] <- c("sx", "sy")
colnames(tsub) <- ct0
C2 <-merge(C1, tsub)
Dist_C2 <-ddply(C2, .(cell_no, region, T.min, geno), summarize, xc = X.micron - sx, yc = Y.micron -sy)
Dist_C3 <-ddply(Dist_C2, .(cell_no, region, T.min, geno), summarize, 
      d2= sqrt((xc)^2 + (yc)^2))
library(ggplot2)
s1<- subset(Dist_C3, T.min>1000 & T.min<3000)

test1<- ddply(s1, .(region, cell_no, geno), summarize, 
              I_d2=(lm(d2~T.min))$coefficients [1],
             I_time=abs((lm(d2~T.min))$coefficients [2]))

test3<- ddply(test1, .(region, geno), summarize, mIt = mean(I_time), sem = std.error(I_time))
test1$geno <- factor(test1$geno, levels = c("nkx", "dlx"))
test1sub = subset(test1, region == "N-tip")
test1sub2 = subset(test1, region == "soma")

#statistical test
Ntip_pv <- wilcox.test(I_time ~ geno, data=test1sub) 
soma_pv <- wilcox.test(I_time ~ geno, data=test1sub2) 
p.adjust(0.01166, method = "bonferroni", n=2)
p.adjust(0.6334, method = "bonferroni", n=2)

#fig4b plot
gp2 <- ggplot(test1, aes(y= I_time, x = geno, colour = geno))+
  facet_grid(region ~. , scales="free_y") +
  geom_boxplot(outlier.shape = NA) + geom_jitter()+
  scale_fill_manual(values=c("#555555","#FF6600"))+
  scale_colour_manual(values=c("#555555","#FF6600"))+
  ylab("Average velocity (um per min)")+
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title=element_text(size=18,face="bold"),
        legend.key = element_rect(fill = "white"),
        legend.background = element_rect(fill = "white"),
        legend.text = element_text(size = 14),
        #panel.grid.major = element_line(colour = "grey40"),
        #panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold", size = 14, colour = "black")
  )
