################script for random regression models for GWAs of longitundinal traits (greeness)
rm(list=ls())

pheno<- read.table("phen_rrm.txt") ##phenotype file for running the softwares 
colnames(pheno)[4]<-"AGE"
age_year <- seq(min(pheno$AGE),max(pheno$AGE))
age_year <- seq(19,122)
head(age_year)

max(pheno$AGE)
min(pheno$AGE)

#Packages
install.packages('orthopolynom')
#install.packages('matrixStats')
#install.packages("matrixStats")
library(matrixStats)
library(orthopolynom)

### creating a legendre
### The time values must be scaled to a range from -1 to +1. The scaling formula is
at <- as.matrix(-1+2*(age_year-min(age_year))/(max(age_year)-min(age_year)))
M <- cbind(as.matrix(rep(1,nrow(at))),at,at^2,at^3,at^4)
n <- 4 ##number of legendre polynomial 

# Coefficients of lengendre
b <- legendre.polynomials(n, normalized=T)
gama0 <- polynomial.coefficients(b)
head(gama0)

#gama0[[5]]
aux[upper.tri(aux)] <- 0
gamat <- aux
gama <- t(gamat)

# Creating the independent transformed variable
K <- M%*%gama
K<-K[,c(1:2)] ## because its for the L1 
head(K)
tail(K)
di
comp = read.table('solutions',skip=1)
head(comp)
coef0=comp[comp[,2]==3,][,3:4]
dim(coef0)
coef1=comp[comp[,2]==4,][,3:4]
dim(coef1)
#coef2=comp[comp[,2]==7,][,3:4]
#dim(coef2)
#coef3=comp[comp[,2]==8,][,3:4]
#dim(coef3)

coefs= cbind(coef0,coef1)#,coef2,coef3)
head(coefs)

coefs1=coefs[,c(2,4)]#6,8)]
head(coefs1)


VG=matrix(0,nrow(coef0),1)

for(i in 1:nrow(coef0))
{
  VG[i,] <- sum((coef0[i,2]*K[,1])+(coef1[i,2]*K[,2]))#+(coef2[i,2]*K[,3])+(coef3[i,2]*K[,4]))
}

head(VG)
dim(VG) 

###########this is to estimate a total GEVB for each line

## 1) i run blupf90 and then I run postGSf90. I used the GEVBs for each coefficient in the sollution file from blupf90 to run postGSf90  and obtain SNP effects for each coefficient 
##from postGSf90 we obtained the results across flight dates for each one of the coefficients there is even manhattan plots

##estiamating the SNP effects for each DAP
snp_sol<- read.table("snp_sol", header = TRUE)

###obtaining the matrix T 
x1<- 19:122  ####using DAP


leg5coef <- legendre.polynomials(n=1, normalized=TRUE) #choose degree and if it should be normalized
head(leg5coef)

leg5 <- as.matrix(as.data.frame(polynomial.values(polynomials=leg5coef, 
                                                  x=scaleX(x1, u=-1, v=1)))) #scaling x
colnames(leg5) <- c("leg0", "leg1") #, "leg2","leg3", "leg4")
head(leg5)

###the leg2 matrix have the coefficients of snp
a<- subset(snp_sol,snp_sol$effect==4)
b<- subset(snp_sol,snp_sol$effect==5)
a<- subset(snp_sol,snp_sol$effect==3)
b<- subset(snp_sol,snp_sol$effect==4)

c<- merge(a,b, by="snp")
c<- c[,c(1,6,14)]
rownames(c)<- c$snp
c<-c[,-1]
c<-t(c)
leg5<-as.matrix(leg5)
c<- as.matrix(c)



#########multiply both matrices######################
SNP_effectDAP<-leg5 %*% c

##now I merge with the SNP IDS that are in the map file
map<- read.table("map1.txt", header = TRUE)

colnames(SNP_effectDAP)<- map$MARKER
names<-as.data.frame(map$MARKER)
SNP_effectDAP[1:10,1:20]

########loop to pick the top 10 SNPs per day 
write.csv(SNP_effectDAP, file="SNP_effects.csv")
SNP_effectDAP<- t(SNP_effectDAP)

SNP_effectDAP<- cbind(names, SNP_effectDAP)
colnames(SNP_effectDAP)[1]<-"MARKER"



names<- c("MARKER", 19:122 )
colnames(SNP_effectDAP)<- names

##filter from 69 to 122
SNP_effectDAP<- SNP_effectDAP[,c(1,53:105)]


data_long <- gather(SNP_effectDAP, DAP, SNP_effect, '19':'122', factor_key=TRUE)
TopSNP<-data_long[order(data_long[,3], decreasing = TRUE),]
highSNP<- TopSNP[1:4000,]
write.csv(highSNP, file="highSNP.csv")


####filtering only the previouse SNPS)
list<-read.csv("GWAsSNPs.csv") ## Significant SNPs identified for single point gwas of related traits
list<- unique(list[,12])
list<- as.character(list)
data_long <- gather(SNP_effectDAP, DAP, SNP_effect, '19':'122', factor_key=TRUE)

data_long<- subset(data_long, data_long$MARKER %in% list)


data_long$MARKER <- factor(data_long$MARKER)
data_long$SNP_effect <- as.numeric(data_long$SNP_effect)
data_long$DAP <- as.numeric(as.character(data_long$DAP))
###removing the last letters 

marker<- data_long%>% separate(MARKER,c("x","chr","Position", "y","z")) 
marker$chr<- as.integer(marker$chr)
marker$MARKER<- paste(marker$chr,marker$Position,sep = "_")
data_long$MARKER<- marker$MARKER
data_long$DAP<- as.integer(data_long$DAP)


names<- unique(data_long$MARKER)

data_long$MARKER <- factor(data_long$MARKER, levels=names)
#display factor levels for region
levels(data_long$MARKER)

png("SNPsRRMGWAsSNPs.png", width = 800, height = 800)
ggplot(data_long, aes(DAP, y= MARKER, fill = SNP_effect)) + 
  geom_tile(colours="black") + 
  scale_fill_viridis_c(na.value="transparent") +
  scale_x_continuous(breaks = seq(19,122,by=20))+
  labs(x="Days After Planting (DAP)",
       y="SNPs",
       fill="SNP effect") +
  theme_classic(base_size = 18)
dev.off()















############################
TopSNP<- data.frame(matrix(ncol=1,nrow = 0))
colnames(TopSNP)<-"MARKER"
for (i in 2:ncol(SNP_effectDAP)){
  df<- SNP_effectDAP[order(SNP_effectDAP[,i], decreasing = TRUE),]
  df<- df[1:70,c(1,i)]
  TopSNP<- merge(TopSNP,df, by="MARKER", all=TRUE)
  }

write.csv(TopSNP, file="top_SNPs_1.csv")




#####################extracting h2 and variance  to to this we use the leg5 matrix 
data=read.table("postgibbs_samples")[,-(1:3)]
colnames(data) <- c("vgC0","covaC0C1","vaC1","ve")


a<-colMeans(data)              

##spliting the vector a in G and R
G<- a[1:3]
library(patr1ckm)

G<-vec2sym(G, diagonal = NULL, lower = FALSE, byrow = TRUE) ##variance covariance
G<-as.matrix(G)

G_cor<- cov2cor(G) ## additive genetic correlation 

residual <- a[4]

##leg 5 matrix of covariates for the function 
##covariance matrix for all DAP 
cov<- leg5 %*% G %*% t(leg5)
write.csv(cov,file="varcov.csv")

Gcor<- cov2cor(cov) #
colnames(Gcor)<- x1
rownames(Gcor)<-x1

write.csv(Gcor, file="Gcorr.csv")
###heritability
diag<- diag(cov)
diag<- as.data.frame(diag)
h2<- diag%>% mutate(h2=diag/(diag+residual))
h2<-h2[,-1]
h2<-cbind(x1,h2)
head(h2)

####genetic correaltions between days for both year and locations

library(reshape2)
Gcor<- read.csv("Gcorr70_122.csv")
names<- 70:122
Gcor<- Gcor[,-1]

colnames(Gcor)<- names
rownames(Gcor)<- names
Gcor<-as.matrix(Gcor)
melted_cormat <- melt(Gcor)
head(melted_cormat)

png
png("Gcor_DAP70_122.png", width = 1000, height = 1000)
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "black")+
  scale_fill_viridis_c(limit = c(-1,1), 
                       name="Genetic Correlation") +xlab("Days after Planting") + ylab("Days after Planting")+
  theme(legend.position = "bottom", text = element_text(size = 24,  color = "black"), axis.title = element_text(size = 26,color = "black"), axis.text = element_text(size = 20,color = "black"), legend.key.width= unit(4, 'cm'), plot.background=element_blank(),
        #remove plot border
        panel.border=element_blank())+
   scale_x_continuous(breaks = seq(70,122,by=5))+
  scale_y_continuous(breaks = seq(70,122,by=5))
 

dev.off()

##plot of heritability
library(ggplot2)

h2<- as.data.frame(h2)

png("H2_DAPbothyears.png", width = 800, height = 600)
 ggplot(h2, aes(x1, h2)) +
  geom_line(na.rm = TRUE) +  
   geom_point(na.rm= TRUE,color="purple")+
  ggtitle("") +
  xlab("Days After Planting (DAP)") + ylab("Narrow-sense heritability") +
               theme_classic(base_size = 20)+
     scale_x_continuous(breaks = seq(19,122,by=10))
 dev.off()

 
 
 
 
 
 ###graph across years 
 #load all heatmaps 
All<-read.csv("/home/descamil/Maturity/nlme/bothyears/L1/top_SNPs.csv")
All<-read.csv("/home/descamil/Maturity/nlme/bothyears/L1/SelectedSNPsRRM.csv")
All<- All[,-1]
names<- c("MARKER", 19:122 )
colnames(All)<- names
data_long <- gather(All, DAP, SNP_effect, '19':'122', factor_key=TRUE)
write.csv(data_long, file="data_long.csv")
data_long<-read.csv("data_long.csv")
data_long<-data_long[,-1]


library(ggplot2)
data_long<- read.csv("highSNP.csv") ###########using this one with the top 1% of SNP effect 
data_long<-data_long[,-1]
data_long$MARKER <- factor(data_long$MARKER)
data_long$SNP_effect <- as.numeric(data_long$SNP_effect)
data_long$DAP <- as.numeric(data_long$DAP)
###removing the last letters 
marker<- data_long%>% separate(MARKER,c("x","chr","Position", "y","z")) 
marker$chr<- as.integer(marker$chr)
marker$MARKER<- paste(marker$chr,marker$Position,sep = "_")
data_long$MARKER<- marker$MARKER

#both<- subset(data_long,data_long$group=='2019-20')
#write.csv(both, file="long_all.csv")  #=LEFT(B2,LEN(B2)-4) excel formula to organize marker name

All$MARKER<- paste(All$Chr, All$Position, sep="_")
data_long<- All[,c(5,10,6)]
write.csv(data_long, file = "selected_SNPs.csv")
a<- table(data_long$MARKER)
a<- as.data.frame(a)
a<-a[order(a[,2], decreasing = TRUE),]

list<- unique(a$Freq)
list1<- list[list<24]
list<- list[list>23]

#flight_list<- c(244, 245, 247, 249, 251, 253, 255, 258, 260, 262, 266, 268, 273, 280)
names<- subset(a, a$Freq %in% list1)
names<- names[,1]

names1<- subset(a, a$Freq %in% list)
names1<- names1[,1]
#both<-read.csv("/home/descamil/Maturity/nlme/bothyears/L1/long_all.csv")

##first plot 
dat<- subset(data_long, data_long$MARKER %in% names)
#dat$MARKER<- as.factor(dat$MARKER)
#levels(dat$MARKER)<- names
dat$MARKER <- factor(dat$MARKER, levels=names)
#display factor levels for region
levels(dat$MARKER)

png("SNPsRRMset1.png", width = 800, height = 800)
ggplot(dat, aes(DAP, y= MARKER, fill = SNP_effect)) + 
  geom_tile(colour = "white") + 
  scale_fill_viridis_c(na.value="transparent") +
   scale_x_continuous(breaks = seq(70,122,by=5))+
   labs(x="Days After Planting (DAP)",
       y="SNPs",
       fill="SNP effect") +
  theme_linedraw(base_size = 20)

dev.off()
 
 
##second plot 
dat1<- subset(data_long, data_long$MARKER %in% names1)
#dat$MARKER<- as.factor(dat$MARKER)
#levels(dat$MARKER)<- names
dat1$MARKER <- factor(dat1$MARKER, levels=names1)
#display factor levels for region
levels(dat1$MARKER)

png("SNPsRRMset2.png", width = 800, height = 800)
ggplot(dat1, aes(DAP, y= MARKER, fill = SNP_effect)) + 
  geom_tile(colour = "white") + 
  scale_fill_viridis_c(na.value="transparent") +
  scale_x_continuous(breaks = seq(70,122,by=5))+
  labs(x="Days After Planting (DAP)",
       y="SNPs",
       fill="SNP effect") +
  theme_linedraw(base_size = 20)

dev.off()


################
Effect<- rep(1,31,by=1)
Effect<- as.data.frame(Effect)
names<- as.data.frame(names)
names<- cbind(names,Effect)

Effect<- rep(2,33,by=1)
Effect<- as.data.frame(Effect)
names1<- as.data.frame(names1)
names1<- cbind(names1,Effect)
colnames(names1)[1]<- "names"
dat<- rbind(names,names1)


write.csv(dat, file = "SelectedSNP_list.csv")



###creating a table with the flights and the time period where R8 was measured 
phen<- read.table("phen_rrm.txt", header = FALSE)
a<- subset(phen,phen$V2=='ACRE2019')
a$V4<- a$V4-165

b<-subset(phen,phen$V2=='ACRE2020')
b$V4<- b$V4-144
phen<-rbind(a,b)
write.table(phen,file="phen_rrm.txt",sep=" ",row.names=F,col.names=F,quote=F)


phen<- phen[,c(2,4)]
phen<-distinct(phen)


library(matrixStats)
library(orthopolynom)
library(vistime)


data<- read.csv("/home/descamil/Maturity/Flightplotdata.csv")

data$event<- as.factor(data$event)
data$Group<- as.factor(data$Group)
data$start<- as.numeric(data$start)
dat<- data[-c(1:26),]
fly<- data[c(1:26),]
fly$start<-as.integer(fly$start)
#fly<-fly[,-1]

png("flights.png", width = 1000, height = 400)

ggplot(dat, aes(x=start, y= Group,  colour=as.factor(event))) +
  geom_line(size=15)+
  scale_color_brewer(palette="Dark2")+
  geom_point( data = fly, color = "black", aes( shape=event ),size=8,stroke=2)+
scale_shape_manual(values=c(18, 20))+
  labs(x="Days After Planting (DAP)",
       y="") +theme(legend.position = "bottom",panel.border = element_rect(colour = "black", fill=NA, size=1),  legend.title=element_blank(),text = element_text(size = 20,  color = "black"), axis.title = element_text(size = 20,face="bold",color = "black"), axis.text = element_text(size = 20,color = "black"),panel.background = element_rect(fill = NA),
                    panel.grid.major = element_line(colour = "grey"),
                    panel.ontop = TRUE, panel.grid.major.x = element_blank(),
                    panel.grid.minor.x = element_blank())+
guides(colour = guide_legend(override.aes = list(size=6, stroke=3)))+
  scale_x_continuous(breaks = seq(0,160,by=20))

    
dev.off()












