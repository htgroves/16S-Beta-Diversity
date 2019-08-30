library(Biostrings) 
library(ggplot2)
library(plyr)
library(stringr)
library(gridExtra)
library(reshape2)
library(knitr)
library(vegan)
library(DESeq2)
library(phyloseq)
library(scales)
library(ggthemes)
library(xlsx)
library(rJava)

load("HG0115_rarefied")

# beta diversity following RSV infection

## splitting dates post infection

d0d4_a = subset_samples(RSV_a_rare, DatePostInfection%in%c("D0", "D4"))

d0d4.dis = distance(d0d4_a, "bray") # this creates a distance matrix

d0d4.ord <- ordinate(d0d4_a, "NMDS", "d0d4.dis") # can put just bray here but have replaced it with the distance matrix

d0d4.df = data.frame(sample_data(d0d4_a)) # have to make the sample data a dataframe

d0d4.adonis = adonis(d0d4.dis ~ DatePostInfection, d0d4.df)

print(d0d4.adonis) # 0.011

d0d7_a = subset_samples(RSV_a_rare, DatePostInfection%in%c("D0", "D7"))

d0d7.dis = distance(d0d7_a, "bray") # this creates a distance matrix

d0d7.ord <- ordinate(d0d7_a, "NMDS", "d0d7.dis") # can put just bray here but have replaced it with the distance matrix

d0d7.df = data.frame(sample_data(d0d7_a)) # have to make the sample data a dataframe

d0d7.adonis = adonis(d0d7.dis ~ DatePostInfection, d0d7.df)

print(d0d7.adonis) # 0.01

##

RSV.a.dis = distance(RSV_a_rare, "bray") 

RSV.ord.a <- ordinate(RSV_a_rare, "NMDS", "RSV.a.dis") # stress = 0.17

RSV.a.df = data.frame(sample_data(RSV_a_rare)) 

RSV.a.adonis = adonis(RSV.a.dis ~ DatePostInfection, RSV.a.df)

print(RSV.a.adonis) # p = 0.001

RSVa_NMDS_bray = plot_ordination(RSV_a_rare, RSV.ord.a, type="samples", color="DatePostInfection", title="RSV infection", justDF = TRUE)

plot.new()

dpiRa = factor(sample_data(RSV_a_rare)$DatePostInfection)# greats a factor which contains the grouping structure for the data

print(levels(dpiRa))

ord_RSV_a = ordiellipse(RSV.ord.a, dpiRa, display = "sites", kind = "se", conf = 0.95, label = T)


df_ell <- data.frame()
for(g in levels(RSVa_NMDS_bray$DatePostInfection)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(RSVa_NMDS_bray[RSVa_NMDS_bray$DatePostInfection==g,],
                                                   veganCovEllipse(ord_RSV_a[[g]]$cov,ord_RSV_a[[g]]$center,ord_RSV_a[[g]]$scale)))
                                ,DatePostInfection=g))
}

plot_ordination(RSV_a_rare, RSV.ord.a, type="samples", shape="DatePostInfection", title="RSV") +
  geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2, shape=DatePostInfection), size=1,linetype=1) +
  labs(shape="Day") +
  geom_point(size=3) +
  theme_base() +
  theme(strip.text.x = element_text(size=18)) +
  theme(plot.title = element_text(face = "bold", size = 25)) + 
  theme(axis.text = element_text(size = 18)) +
  theme(axis.title.x = element_text(size = 22)) +
  theme(axis.title.y = element_text(size = 22)) +
  theme(plot.title = element_text(margin=margin(0,0,20,0))) +
  theme(axis.title.y = element_text(margin=margin(0,20,0,0))) +
  theme(axis.title.x = element_text(margin=margin(20,0,0,0))) +
  theme(plot.title = element_text(hjust=0.5))  +
  theme(plot.margin = unit(c(0,4,0,0), "cm")) +
  theme(plot.background=element_blank()) 
