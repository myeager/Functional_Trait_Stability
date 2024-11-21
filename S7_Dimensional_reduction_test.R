library(ggplot2)
library(FD)
library(tidyr)
library(dplyr)

#spc.df <- read.csv("/Users/mallarie.yeager/Documents/R_projects/FT Stability/Ordination Analyses/Outputs/Species_comm_dat.csv")
#FT.df <- read.csv("/Users/mallarie.yeager/Documents/R_projects/FT Stability/Ordination Analyses/Outputs/Functional_trait_comm_dat.csv")
#comm.df <- read.csv("Ordination Analyses/Data/com.dat.fish.csv")
FT_comm_df <- read.csv("/Users/mallarie.yeager/Documents/R_projects/Functional_Trait_Stability/02_Ordination Analysis/Data/FT_community_df.csv")


#spc_fam <- unique(FT_comm_df[,1:2])
#comm.df <- merge(spc_fam, comm.df, by = "Common_name")

wes_palettes <- list(
  Darjeeling1 = c("#FF0000", "#00A08A", "#F2AD00", "#F98400", "#5BBCD6"),
  Cavalcanti1 = c("#D8B70A", "#02401B", "#A2A475", "#81A88D", "#972D15"),
  GrandBudapest1 = c("#F1BB7B", "#FD6467", "#5B1A18", "#D67236")
)

ponds <- c("NP","PP","PJ","WP", "GH","QP")
pal <- data.frame(Pond = ponds,col = c(wes_palettes$Darjeeling1[1], 
                                       wes_palettes$Darjeeling1[2], 
                                       wes_palettes$Darjeeling1[5],
                                       wes_palettes$Darjeeling1[4], 
                                       wes_palettes$Cavalcanti1[1],
                                       wes_palettes$GrandBudapest1[3]))



#new way to avoid demensional reduction: (1) sub sample 13 species, 
#(2) calculate FTs,(3) calculate instability for species and traits, (4) repeat
#8 species each

####13 species and FT each####
uni.fam <- unique(FT_comm_df$Family)

uni.spc <- unique(FT_comm_df$Common_name)
uni.fts <- colnames(FT_comm_df[,11:23])

MCMC_yr_centroid_dist_FT <- NULL
MCMC_yr_centroid_dist_COM <- NULL
MCMC_yr_yr_dist_FT <- NULL
MCMC_yr_yr_dist_COM <- NULL
spc_assemblges <- NULL
FT_assemblges <- NULL

for(j in 1:10000){
  smp.col <- sample(1:27, 13, replace = F)
  spc.smp <- uni.spc[smp.col]
  spc_assemblges <- rbind(spc_assemblges, cbind(spc.smp,j))
  
  #FT.col <- sample(11:23, 13, replace = F)
  #FT.smp <- uni.fts[FT.col-10]
  #FT_assemblges <- rbind(FT_assemblges, cbind(FT.smp,j))
  
  com.smp <- FT_comm_df %>%
    filter(Common_name %in% spc.smp)
  #com.smp <- cbind(com.smp[,1:10], com.smp[,FT.col])
  
  FT.mat1  <- com.smp[!duplicated(com.smp[,c(2,8)]),] #removes duplicates of spc*size combo
  FT.mat1 <- FT.mat1[order(FT.mat1$FL),] #reorders matrix by FL 
  FT.mat1 <- FT.mat1[order(FT.mat1$Common_name),] #then Common name
  
  FT.mat <- as.data.frame(FT.mat1[,c(11:ncol(FT.mat1))]) #Creates trait matrix with spc*size row names
  rownames(FT.mat) <- paste(FT.mat1$Common_name, FT.mat1$FL) #creates rownames
  #sapply(FT.mat[,1:12], FUN = min) #checks to make sure FTs are not negative
  
  fish.freq <- aggregate(FREQ ~ Common_name + FL + Pond + Year, data = com.smp, sum) #sums freq up for all spc * size combos
  fish.freq <- fish.freq[order(fish.freq$FL),] #reorders matrix by FL 
  fish.freq <- fish.freq[order(fish.freq$Common_name),] #then Common name
  fish.freq$spc_fl <- paste(fish.freq$Common_name, fish.freq$FL) #creates a species x length column
  fish.freq <- fish.freq[,-c(1:2)] #removes common name and length cols
  abund <- reshape(fish.freq, timevar = "spc_fl" , idvar = c("Pond", "Year"), direction = "wide") #reshapes data to wide format - spc*size on cols
  abund.mat <- abund[,-c(1:2)] #removes pond and year columns
  row.names(abund.mat) <- paste(abund$Pond, abund$Year, abund$month) #creates rownames
  names(abund.mat) <- substring(names(abund.mat), 6) #cleans up col names
  abund.mat[is.na(abund.mat)] <- 0 #places zeros were NAs were
  if (nrow(abund.mat) != 36){
    print("not all pond years represented")
  } else{
    
    #FT and Abundance Check#
    check <- cbind(rownames(FT.mat), colnames(abund.mat)) #checks to make sure FT rownames and abund colnames match
    
    
    FD.test <- FD::functcomp(x = as.matrix(FT.mat), a = as.matrix(abund.mat))
    
    FT.MC <- FD.test
    FT.MC <- cbind(rownames(FD.test), FT.MC)
    FT.MC <- cbind(substr(FT.MC[,1], 4, 7), FT.MC)
    FT.MC <- cbind(substr(FT.MC[,2], 1, 2), FT.MC)
    colnames(FT.MC)[1:2] <- c("Pond", "Year")
    FT.MC <- FT.MC[,-3]
    FT.MC <- FT.MC[order(FT.MC$Year),]
    FT.MC <- FT.MC[order(FT.MC$Pond),]
    
    FT.MC.mds <- metaMDS(FT.MC[,3:ncol(FT.MC)], k = 2, trymax = 10000, maxit = 10000)
    coords <- as.data.frame(FT.MC.mds$points)
    coords <- cbind(FT.MC[,1:2], coords)
    coords <- merge(coords, pal, by = "Pond")
    coords <- coords[order(coords$Year),]
    coords <- coords[order(coords$Pond),]
    
    plot(coords$MDS1, coords$MDS2,  type = "n", xlab = "nMDS1", xlim = c(-1,1), 
         ylim = c(-1,1), ylab = "nMDS2")
    abline(h = 0, lty = 2)
    abline(v = 0, lty = 2)
    #orditorp(com.mds, display = "species", air = 0.01, cex = 1)
    for (i in 1:nrow(coords)) {
      text(coords$MDS1[i], coords$MDS2[i], 
           labels = substr(as.character(coords$Year[i]), 
                           start = 3, stop = 4),
           col = c(as.character(coords$col[i])),
           pos = 1, cex = 0.75)
    }
    NP.e <- ordiellipse(FT.MC.mds, kind = "se", groups = FT.MC$Pond, show = "NP",
                        lwd = 2, col = as.character(coords$col[coords$Pond == "NP"]))
    PP.e <- ordiellipse(FT.MC.mds, kind = "se",groups = FT.MC$Pond, show = "PP", 
                        lwd = 2, col = as.character(coords$col[coords$Pond == "PP"]))
    PJ.e <- ordiellipse(FT.MC.mds, kind = "se",groups = FT.MC$Pond, show = "PJ", 
                        lwd = 2, col = as.character(coords$col[coords$Pond == "PJ"]))
    WP.e <- ordiellipse(FT.MC.mds, kind = "se",groups = FT.MC$Pond, show = "WP", 
                        lwd = 2, col= as.character(coords$col[coords$Pond == "WP"]))
    GH.e <- ordiellipse(FT.MC.mds, kind = "se",groups = FT.MC$Pond, show = "GH", 
                        lwd = 2, col = as.character(coords$col[coords$Pond == "GH"]))
    QP.e <- ordiellipse(FT.MC.mds, kind = "se", groups = FT.MC$Pond, show = "QP",
                        lwd = 2, col = as.character(coords$col[coords$Pond == "QP"]))
    
    
    tmp.MCMC_yr_centroid_dist <- data.frame(Pond = coords$Pond, Year = coords$Year,
                                            distance = NA, rel_distance = NA, 
                                            stress = FT.MC.mds$stress, converged = 
                                              FT.MC.mds$converged > 0, rep = j)
    
    tmp.MCMC_yr_yr_dist <- data.frame(Pond = c(rep("GH",5), rep("NP",5), rep("PJ",5),
                                               rep("PP",5), rep("QP",5), rep("WP",5)), 
                                      Year = rep(c("10_11", "11_12","12_13", "13_14", "14_15"),6),
                                      distance = NA, rel_distance = NA, 
                                      stress = FT.MC.mds$stress, converged = 
                                        FT.MC.mds$converged > 0, rep = j)
    
    
    #Calculate distance from centroid FT
    distmat.ft <- vegdist(FT.MC[,3:ncol(FT.MC)])
    com.beta <- betadisper(distmat.ft, group = FT.MC$Pond, type = "centroid")
    tmp.MCMC_yr_centroid_dist[,3] <-  com.beta$distances
    tmp.MCMC_yr_centroid_dist[,4] <- tmp.MCMC_yr_centroid_dist$distance/max(tmp.MCMC_yr_centroid_dist$distance)
    MCMC_yr_centroid_dist_FT <- rbind(MCMC_yr_centroid_dist_FT, tmp.MCMC_yr_centroid_dist)
    
    
    #calc yr-to-yr distance
    distmat.ft.mat <- as.matrix(distmat.ft)
    tmp.dist.ft <- NULL
    for(i in 1:nrow(distmat.ft.mat) - 1){ 
      tmp.dist.ft <- c(tmp.dist.ft, distmat.ft.mat[i,i+1])
    }
    tmp.dist.ft <- as.data.frame(tmp.dist.ft[-c(6,12,18,24,30)])
    tmp.MCMC_yr_yr_dist[,3] <-  tmp.dist.ft
    tmp.MCMC_yr_yr_dist[,4] <- tmp.dist.ft/max(tmp.dist.ft)
    MCMC_yr_yr_dist_FT <- rbind(MCMC_yr_yr_dist_FT, tmp.MCMC_yr_yr_dist)
    
    #Community comp
    avg.by.year <- aggregate(FREQ ~ Pond + Year + Common_name, FUN = mean, data = com.smp)
    sums.3 <- reshape(avg.by.year, timevar = c("Common_name"), ids = c(FREQ), idvar = c("Pond","Year"), direction = "wide")
    sums.2 <- sums.3[,3:ncol(sums.3)]
    #removes the "FREQ." before each name
    names(sums.2) <- substring(names(sums.2), 6)
    #replaces NA's with 0's
    sums.2[is.na(sums.2)] <- 0
    #sums.1 <- decostand(sums.2[,9:ncol(sums.2)], method = "hellinger")
    sums <- cbind(sums.3[,1:2], sums.2)
    sums <- sums[order(sums$Year),]
    sums <- sums[order(sums$Pond),]
   # sums$Pd_Yr <- paste(sums$Pond, sums$Year)
    #sums <- sums[,c(1:2,3:10)]
    
    COM.MC <- sums
    #write_csv(COM.MC, "CWM.df.csv")
    COM.MC <- COM.MC[order(COM.MC$Year),]
    COM.MC <- COM.MC[order(COM.MC$Pond),]
    
    COM.MC.mds <- metaMDS(COM.MC[,3:ncol(COM.MC)], k = 2, trymax = 10000, maxit = 10000)
    coords <- as.data.frame(COM.MC.mds$points)
    coords <- cbind(COM.MC[,1:2], coords)
    coords <- merge(coords, pal, by = "Pond")
    coords <- coords[order(coords$Year),]
    coords <- coords[order(coords$Pond),]
    
    plot(coords$MDS1, coords$MDS2,  type = "n", xlab = "nMDS1", xlim = c(-1,1), 
         ylim = c(-1,1), ylab = "nMDS2")
    abline(h = 0, lty = 2)
    abline(v = 0, lty = 2)
    #orditorp(com.mds, display = "species", air = 0.01, cex = 1)
    for (i in 1:nrow(coords)) {
      text(coords$MDS1[i], coords$MDS2[i], 
           labels = substr(as.character(coords$Year[i]), 
                           start = 3, stop = 4),
           col = c(as.character(coords$col[i])),
           pos = 1, cex = 0.75)
    }
    NP.e <- ordiellipse(COM.MC.mds, kind = "se", groups = COM.MC$Pond, show = "NP",
                        lwd = 2, col = as.character(coords$col[coords$Pond == "NP"]))
    PP.e <- ordiellipse(COM.MC.mds, kind = "se",groups = COM.MC$Pond, show = "PP", 
                        lwd = 2, col = as.character(coords$col[coords$Pond == "PP"]))
    PJ.e <- ordiellipse(COM.MC.mds, kind = "se",groups = COM.MC$Pond, show = "PJ", 
                        lwd = 2, col = as.character(coords$col[coords$Pond == "PJ"]))
    WP.e <- ordiellipse(COM.MC.mds, kind = "se",groups = COM.MC$Pond, show = "WP", 
                        lwd = 2, col= as.character(coords$col[coords$Pond == "WP"]))
    GH.e <- ordiellipse(COM.MC.mds, kind = "se",groups = COM.MC$Pond, show = "GH", 
                        lwd = 2, col = as.character(coords$col[coords$Pond == "GH"]))
    QP.e <- ordiellipse(COM.MC.mds, kind = "se", groups = COM.MC$Pond, show = "QP",
                        lwd = 2, col = as.character(coords$col[coords$Pond == "QP"]))
    
    tmp.MCMC_yr_centroid_dist <- data.frame(Pond = coords$Pond, Year = coords$Year,
                                            distance = NA, rel_distance = NA, 
                                            stress = COM.MC.mds$stress, converged = 
                                              COM.MC.mds$converged > 0, rep = j)
    
    tmp.MCMC_yr_yr_dist <- data.frame(Pond = c(rep("GH",5), rep("NP",5), rep("PJ",5),
                                               rep("PP",5), rep("QP",5), rep("WP",5)), 
                                      Year = rep(c("10_11", "11_12","12_13", "13_14", "14_15"),6),
                                      distance = NA, rel_distance = NA, 
                                      stress = COM.MC.mds$stress, converged = 
                                        COM.MC.mds$converged > 0, rep = j)
    
    #Calculate distance from yr to centroid com 
    distmat.com <- vegdist(COM.MC[,3:ncol(COM.MC)])
    com.beta <- betadisper(distmat.com, group = COM.MC$Pond, type = "centroid")
    tmp.MCMC_yr_centroid_dist[,3] <-  com.beta$distances
    tmp.MCMC_yr_centroid_dist[,4] <- tmp.MCMC_yr_centroid_dist$distance/max(tmp.MCMC_yr_centroid_dist$distance)
    MCMC_yr_centroid_dist_COM <- rbind(MCMC_yr_centroid_dist_COM, tmp.MCMC_yr_centroid_dist)
    
    #calc yr-to-yr distance
    dismat.com.mat <- as.matrix(distmat.com)
    tmp.dist.com <- NULL
    for(i in 1:nrow(dismat.com.mat) - 1){ 
      tmp.dist.com <- c(tmp.dist.com, dismat.com.mat[i,i+1])
    }
    tmp.dist.com <- as.data.frame(tmp.dist.com[-c(6,12,18,24,30)])
    tmp.MCMC_yr_yr_dist[,3] <-  tmp.dist.com
    tmp.MCMC_yr_yr_dist[,4] <- tmp.dist.com/max(tmp.dist.com)
    MCMC_yr_yr_dist_COM <- rbind(MCMC_yr_yr_dist_COM, tmp.MCMC_yr_yr_dist)
  }
}
MCMC_yr_yr_dist_FT_13 <-MCMC_yr_yr_dist_FT
MCMC_yr_centroid_dist_FT_13 <- MCMC_yr_centroid_dist_FT
MCMC_yr_yr_dist_COM_13 <- MCMC_yr_yr_dist_COM
MCMC_yr_centroid_dist_COM_13 <- MCMC_yr_centroid_dist_COM

avg.yr_centroid_dist_FT <- MCMC_yr_centroid_dist_FT_13 %>%
  group_by(Pond, Year) %>%
  summarise(mean_rel_distance = mean(rel_distance), mean_distance = mean(distance))%>%
  mutate(type = "traits", Year = as.numeric(Year))


avg.yr_centroid_dist_COM <- MCMC_yr_centroid_dist_COM_13 %>%
  group_by(Pond, Year) %>%
  summarise(mean_rel_distance = mean(rel_distance), mean_distance = mean(distance))%>%
  mutate(type = "species")


avg.yr_centroid_dist <- rbind(avg.yr_centroid_dist_COM, avg.yr_centroid_dist_FT)

beta.yr_to_cent_dist_13 <- merge(avg.yr_centroid_dist, pal, by = "Pond")
summary(aov(mean_rel_distance ~ Pond * type, beta.yr_to_cent_dist))
summary(aov(mean_rel_distance ~ Pond, beta.yr_to_cent_dist))
summary(aov(mean_rel_distance ~ type, beta.yr_to_cent_dist))
mod <- aov(mean_rel_distance ~ Pond, beta.yr_to_cent_dist)
TukeyHSD(mod)
pal <- pal[order(pal$Pond),]
summary(aov(mean_rel_distance ~ Pond * type, beta.yr_to_yr_dist))
par(mfrow = c(2,2))
##Main figs- Fig 1
#yr to centroid
##Pond
###Test
summary(aov.mod <- aov(mean_rel_distance ~ Pond, beta.yr_to_cent_dist_13))
TukeyHSD(aov.mod)
###Plot
boxplot(beta.yr_to_cent_dist_13$mean_rel_distance ~ Pond, beta.yr_to_cent_dist_13, col = pal$col, ylab = "Relative year to centroid distance")
points(aggregate(mean_rel_distance ~ Pond, data = beta.yr_to_cent_dist_13, mean)[2], pch = 16)
##Type
###Test
summary(aov.mod <- aov(mean_rel_distance ~ type, beta.yr_to_cent_dist_13))
TukeyHSD(aov.mod)
###Plot
boxplot(beta.yr_to_cent_dist_13$mean_rel_distance ~ type, beta.yr_to_cent_dist_13,
        col = c("#506E82", "#768655"), ylab = "Relative year to centroid distance",
        names = c("Species", "Functional Traits"), xlab = "")
points(aggregate(mean_rel_distance ~ type, data = beta.yr_to_cent_dist_13, mean)[2], pch = 16)


#yr to yr
avg.yr_yr_dist_FT <- MCMC_yr_yr_dist_FT_13 %>%
  group_by(Pond, Year) %>%
  summarise(mean_rel_distance = mean(rel_distance), mean_distance = mean(distance))%>%
  mutate(type = "traits")


avg.yr_yr_dist_COM <- MCMC_yr_yr_dist_COM_13 %>%
  group_by(Pond, Year) %>%
  summarise(mean_rel_distance = mean(rel_distance), mean_distance = mean(distance))%>%
  mutate(type = "species")


avg.yr_yr_dist <- rbind(avg.yr_yr_dist_COM, avg.yr_yr_dist_FT)

beta.yr_to_yr_dist_13 <- merge(avg.yr_yr_dist, pal, by = "Pond")
summary(aov(mean_rel_distance ~ Pond * type, beta.yr_to_yr_dist_13))
summary(aov(mean_rel_distance ~ Pond, beta.yr_to_yr_dist_13))
summary(aov(mean_rel_distance ~ type, beta.yr_to_yr_dist_13))
mod <- aov(mean_rel_distance ~ Pond, beta.yr_to_yr_dist_13)
TukeyHSD(mod)
##Pond
###Test
summary(aov(mean_rel_distance ~ Pond, beta.yr_to_yr_dist_13))
###Plot
boxplot(beta.yr_to_yr_dist_13$mean_rel_distance ~ Pond, beta.yr_to_yr_dist_13, col = pal$col, ylab = "Relative year to year distance")
points(aggregate(mean_rel_distance ~ Pond, data = beta.yr_to_yr_dist_13, mean)[2], pch = 16)

###Test
summary(aov.mod <- aov(mean_rel_distance ~ type, beta.yr_to_yr_dist_13))
TukeyHSD(aov.mod)
yr.yr.ty.pval <- round(summary(aov.mod)[[1]][["Pr(>F)"]][1],3)
boxplot(beta.yr_to_yr_dist$mean_rel_distance ~ type, beta.yr_to_yr_dist_13,
        col = c("#506E82", "#768655"), ylab = "Relative year to year distance",
        names = c("Species", "Functional Traits"), xlab = "")
points(aggregate(mean_rel_distance ~ type, data = beta.yr_to_yr_dist_13, mean)[2], pch = 16)


## species assemblages

spc_assemblges <- as.data.frame(spc_assemblges)
spc_counts_MCMC <- spc_assemblges %>%
  group_by(spc.smp)%>%
  summarise(prop_of_iter = length(j)/10000)


####8 species and FTs each####

uni.fam <- unique(FT_comm_df$Family)

uni.spc <- unique(FT_comm_df$Common_name)
uni.fts <- colnames(FT_comm_df[,11:23])

MCMC_yr_centroid_dist_FT <- NULL
MCMC_yr_centroid_dist_COM <- NULL
MCMC_yr_yr_dist_FT <- NULL
MCMC_yr_yr_dist_COM <- NULL
spc_assemblges <- NULL
FT_assemblges <- NULL

for(j in 1:10000){
smp.col <- sample(1:27, 8, replace = F)
spc.smp <- uni.spc[smp.col]
spc_assemblges <- rbind(spc_assemblges, cbind(spc.smp,j))

FT.col <- sample(11:23, 8, replace = F)
FT.smp <- uni.fts[FT.col]
FT_assemblges <- rbind(FT_assemblges, cbind(FT.smp,j))

com.smp <- FT_comm_df %>%
  filter(Common_name %in% spc.smp)
com.smp <- cbind(com.smp[,1:10], com.smp[,FT.col])

FT.mat1  <- com.smp[!duplicated(com.smp[,c(2,8)]),] #removes duplicates of spc*size combo
FT.mat1 <- FT.mat1[order(FT.mat1$FL),] #reorders matrix by FL 
FT.mat1 <- FT.mat1[order(FT.mat1$Common_name),] #then Common name

FT.mat <- as.data.frame(FT.mat1[,c(11:ncol(FT.mat1))]) #Creates trait matrix with spc*size row names
rownames(FT.mat) <- paste(FT.mat1$Common_name, FT.mat1$FL) #creates rownames
#sapply(FT.mat[,1:12], FUN = min) #checks to make sure FTs are not negative

fish.freq <- aggregate(FREQ ~ Common_name + FL + Pond + Year, data = com.smp, sum) #sums freq up for all spc * size combos
fish.freq <- fish.freq[order(fish.freq$FL),] #reorders matrix by FL 
fish.freq <- fish.freq[order(fish.freq$Common_name),] #then Common name
fish.freq$spc_fl <- paste(fish.freq$Common_name, fish.freq$FL) #creates a species x length column
fish.freq <- fish.freq[,-c(1:2)] #removes common name and length cols
abund <- reshape(fish.freq, timevar = "spc_fl" , idvar = c("Pond", "Year"), direction = "wide") #reshapes data to wide format - spc*size on cols
abund.mat <- abund[,-c(1:2)] #removes pond and year columns
row.names(abund.mat) <- paste(abund$Pond, abund$Year, abund$month) #creates rownames
names(abund.mat) <- substring(names(abund.mat), 6) #cleans up col names
abund.mat[is.na(abund.mat)] <- 0 #places zeros were NAs were
if (nrow(abund.mat) != 36){
  print("not all pond years represented")
} else{

#FT and Abundance Check#
check <- cbind(rownames(FT.mat), colnames(abund.mat)) #checks to make sure FT rownames and abund colnames match


FD.test <- FD::functcomp(x = as.matrix(FT.mat), a = as.matrix(abund.mat))

FT.MC <- FD.test
FT.MC <- cbind(rownames(FD.test), FT.MC)
FT.MC <- cbind(substr(FT.MC[,1], 4, 7), FT.MC)
FT.MC <- cbind(substr(FT.MC[,2], 1, 2), FT.MC)
colnames(FT.MC)[1:2] <- c("Pond", "Year")
FT.MC <- FT.MC[,-3]
FT.MC <- FT.MC[order(FT.MC$Year),]
FT.MC <- FT.MC[order(FT.MC$Pond),]

FT.MC.mds <- metaMDS(FT.MC[,3:ncol(FT.MC)], k = 2, trymax = 10000, maxit = 10000)
coords <- as.data.frame(FT.MC.mds$points)
coords <- cbind(FT.MC[,1:2], coords)
coords <- merge(coords, pal, by = "Pond")
coords <- coords[order(coords$Year),]
coords <- coords[order(coords$Pond),]

plot(coords$MDS1, coords$MDS2,  type = "n", xlab = "nMDS1", xlim = c(-1,1), 
     ylim = c(-1,1), ylab = "nMDS2")
abline(h = 0, lty = 2)
abline(v = 0, lty = 2)
#orditorp(com.mds, display = "species", air = 0.01, cex = 1)
for (i in 1:nrow(coords)) {
  text(coords$MDS1[i], coords$MDS2[i], 
       labels = substr(as.character(coords$Year[i]), 
                       start = 3, stop = 4),
       col = c(as.character(coords$col[i])),
       pos = 1, cex = 0.75)
}
NP.e <- ordiellipse(FT.MC.mds, kind = "se", groups = FT.MC$Pond, show = "NP",
                    lwd = 2, col = as.character(coords$col[coords$Pond == "NP"]))
PP.e <- ordiellipse(FT.MC.mds, kind = "se",groups = FT.MC$Pond, show = "PP", 
                    lwd = 2, col = as.character(coords$col[coords$Pond == "PP"]))
PJ.e <- ordiellipse(FT.MC.mds, kind = "se",groups = FT.MC$Pond, show = "PJ", 
                    lwd = 2, col = as.character(coords$col[coords$Pond == "PJ"]))
WP.e <- ordiellipse(FT.MC.mds, kind = "se",groups = FT.MC$Pond, show = "WP", 
                    lwd = 2, col= as.character(coords$col[coords$Pond == "WP"]))
GH.e <- ordiellipse(FT.MC.mds, kind = "se",groups = FT.MC$Pond, show = "GH", 
                    lwd = 2, col = as.character(coords$col[coords$Pond == "GH"]))
QP.e <- ordiellipse(FT.MC.mds, kind = "se", groups = FT.MC$Pond, show = "QP",
                    lwd = 2, col = as.character(coords$col[coords$Pond == "QP"]))


tmp.MCMC_yr_centroid_dist <- data.frame(Pond = coords$Pond, Year = coords$Year,
                                        distance = NA, rel_distance = NA, 
                                        stress = FT.MC.mds$stress, converged = 
                                          FT.MC.mds$converged > 0, rep = j)

tmp.MCMC_yr_yr_dist <- data.frame(Pond = c(rep("GH",5), rep("NP",5), rep("PJ",5),
                                           rep("PP",5), rep("QP",5), rep("WP",5)), 
                                  Year = rep(c("10_11", "11_12","12_13", "13_14", "14_15"),6),
                                  distance = NA, rel_distance = NA, 
                                  stress = FT.MC.mds$stress, converged = 
                                    FT.MC.mds$converged > 0, rep = j)


#Calculate distance from centroid FT
distmat.ft <- vegdist(FT.MC[,3:ncol(FT.MC)])
com.beta <- betadisper(distmat.ft, group = FT.MC$Pond, type = "centroid")
tmp.MCMC_yr_centroid_dist[,3] <-  com.beta$distances
tmp.MCMC_yr_centroid_dist[,4] <- tmp.MCMC_yr_centroid_dist$distance/max(tmp.MCMC_yr_centroid_dist$distance)
MCMC_yr_centroid_dist_FT <- rbind(MCMC_yr_centroid_dist_FT, tmp.MCMC_yr_centroid_dist)


#calc yr-to-yr distance
distmat.ft.mat <- as.matrix(distmat.ft)
tmp.dist.ft <- NULL
for(i in 1:nrow(distmat.ft.mat) - 1){ 
  tmp.dist.ft <- c(tmp.dist.ft, distmat.ft.mat[i,i+1])
}
tmp.dist.ft <- as.data.frame(tmp.dist.ft[-c(6,12,18,24,30)])
tmp.MCMC_yr_yr_dist[,3] <-  tmp.dist.ft
tmp.MCMC_yr_yr_dist[,4] <- tmp.dist.ft/max(tmp.dist.ft)
MCMC_yr_yr_dist_FT <- rbind(MCMC_yr_yr_dist_FT, tmp.MCMC_yr_yr_dist)

#Community comp
avg.by.year <- aggregate(FREQ ~ Pond + Year + Common_name, FUN = mean, data = com.smp)
sums.3 <- reshape(avg.by.year, timevar = c("Common_name"), ids = c(FREQ), idvar = c("Pond","Year"), direction = "wide")
sums.2 <- sums.3[,3:ncol(sums.3)]
#removes the "FREQ." before each name
names(sums.2) <- substring(names(sums.2), 6)
#replaces NA's with 0's
sums.2[is.na(sums.2)] <- 0
#sums.1 <- decostand(sums.2[,9:ncol(sums.2)], method = "hellinger")
sums <- cbind(sums.3[,1:2], sums.2)
sums <- sums[order(sums$Year),]
sums <- sums[order(sums$Pond),]
sums$Pd_Yr <- paste(sums$Pond, sums$Year)
sums <- sums[,c(1:2,3:10)]

COM.MC <- sums
#write_csv(COM.MC, "CWM.df.csv")
COM.MC <- COM.MC[order(COM.MC$Year),]
COM.MC <- COM.MC[order(COM.MC$Pond),]

COM.MC.mds <- metaMDS(COM.MC[,3:ncol(COM.MC)], k = 2, trymax = 10000, maxit = 10000)
coords <- as.data.frame(COM.MC.mds$points)
coords <- cbind(COM.MC[,1:2], coords)
coords <- merge(coords, pal, by = "Pond")
coords <- coords[order(coords$Year),]
coords <- coords[order(coords$Pond),]

plot(coords$MDS1, coords$MDS2,  type = "n", xlab = "nMDS1", xlim = c(-1,1), 
     ylim = c(-1,1), ylab = "nMDS2")
abline(h = 0, lty = 2)
abline(v = 0, lty = 2)
#orditorp(com.mds, display = "species", air = 0.01, cex = 1)
for (i in 1:nrow(coords)) {
  text(coords$MDS1[i], coords$MDS2[i], 
       labels = substr(as.character(coords$Year[i]), 
                       start = 3, stop = 4),
       col = c(as.character(coords$col[i])),
       pos = 1, cex = 0.75)
}
NP.e <- ordiellipse(COM.MC.mds, kind = "se", groups = COM.MC$Pond, show = "NP",
                    lwd = 2, col = as.character(coords$col[coords$Pond == "NP"]))
PP.e <- ordiellipse(COM.MC.mds, kind = "se",groups = COM.MC$Pond, show = "PP", 
                    lwd = 2, col = as.character(coords$col[coords$Pond == "PP"]))
PJ.e <- ordiellipse(COM.MC.mds, kind = "se",groups = COM.MC$Pond, show = "PJ", 
                    lwd = 2, col = as.character(coords$col[coords$Pond == "PJ"]))
WP.e <- ordiellipse(COM.MC.mds, kind = "se",groups = COM.MC$Pond, show = "WP", 
                    lwd = 2, col= as.character(coords$col[coords$Pond == "WP"]))
GH.e <- ordiellipse(COM.MC.mds, kind = "se",groups = COM.MC$Pond, show = "GH", 
                    lwd = 2, col = as.character(coords$col[coords$Pond == "GH"]))
QP.e <- ordiellipse(COM.MC.mds, kind = "se", groups = COM.MC$Pond, show = "QP",
                    lwd = 2, col = as.character(coords$col[coords$Pond == "QP"]))

tmp.MCMC_yr_centroid_dist <- data.frame(Pond = coords$Pond, Year = coords$Year,
                                        distance = NA, rel_distance = NA, 
                                        stress = COM.MC.mds$stress, converged = 
                                          COM.MC.mds$converged > 0, rep = j)

tmp.MCMC_yr_yr_dist <- data.frame(Pond = c(rep("GH",5), rep("NP",5), rep("PJ",5),
                                           rep("PP",5), rep("QP",5), rep("WP",5)), 
                                  Year = rep(c("10_11", "11_12","12_13", "13_14", "14_15"),6),
                                        distance = NA, rel_distance = NA, 
                                        stress = COM.MC.mds$stress, converged = 
                                          COM.MC.mds$converged > 0, rep = j)
                                                          
#Calculate distance from yr to centroid com 
distmat.com <- vegdist(COM.MC[,3:ncol(COM.MC)])
com.beta <- betadisper(distmat.com, group = COM.MC$Pond, type = "centroid")
tmp.MCMC_yr_centroid_dist[,3] <-  com.beta$distances
tmp.MCMC_yr_centroid_dist[,4] <- tmp.MCMC_yr_centroid_dist$distance/max(tmp.MCMC_yr_centroid_dist$distance)
MCMC_yr_centroid_dist_COM <- rbind(MCMC_yr_centroid_dist_COM, tmp.MCMC_yr_centroid_dist)

#calc yr-to-yr distance
dismat.com.mat <- as.matrix(distmat.com)
tmp.dist.com <- NULL
for(i in 1:nrow(dismat.com.mat) - 1){ 
  tmp.dist.com <- c(tmp.dist.com, dismat.com.mat[i,i+1])
}
tmp.dist.com <- as.data.frame(tmp.dist.com[-c(6,12,18,24,30)])
tmp.MCMC_yr_yr_dist[,3] <-  tmp.dist.com
tmp.MCMC_yr_yr_dist[,4] <- tmp.dist.com/max(tmp.dist.com)
MCMC_yr_yr_dist_COM <- rbind(MCMC_yr_yr_dist_COM, tmp.MCMC_yr_yr_dist)
}
}
MCMC_yr_centroid_dist_COM_8 <- MCMC_yr_centroid_dist_COM
MCMC_yr_yr_dist_COM_8 <- MCMC_yr_yr_dist_COM

MCMC_yr_yr_dist_FT_8 <- MCMC_yr_yr_dist_FT
MCMC_yr_centroid_dist_FT_8 <- MCMC_yr_centroid_dist_FT

avg.yr_centroid_dist_FT <- MCMC_yr_centroid_dist_FT_8 %>%
  group_by(Pond, Year) %>%
  summarise(mean_rel_distance = mean(rel_distance), mean_distance = mean(distance))%>%
  mutate(type = "traits", Year = as.numeric(Year))


avg.yr_centroid_dist_COM <- MCMC_yr_centroid_dist_COM_8 %>%
  group_by(Pond, Year) %>%
  summarise(mean_rel_distance = mean(rel_distance), mean_distance = mean(distance))%>%
  mutate(type = "species")


avg.yr_centroid_dist <- rbind(avg.yr_centroid_dist_COM, avg.yr_centroid_dist_FT)

beta.yr_to_cent_dist_8 <- merge(avg.yr_centroid_dist, pal, by = "Pond")
summary(aov(mean_rel_distance ~ Pond * type, beta.yr_to_cent_dist_8))
summary(aov(mean_rel_distance ~ Pond, beta.yr_to_cent_dist_8))
summary(aov(mean_rel_distance ~ type, beta.yr_to_cent_dist_8))
mod <- aov(mean_rel_distance ~ Pond, beta.yr_to_cent_dist_8)
TukeyHSD(mod)
pal <- pal[order(pal$Pond),]
summary(aov(mean_rel_distance ~ Pond * type, beta.yr_to_cent_dist_8))
par(mfrow = c(2,2))
##Main figs- Fig 1
#yr to centroid
##Pond
###Test
summary(aov.mod <- aov(mean_rel_distance ~ Pond, beta.yr_to_cent_dist_8))
TukeyHSD(aov.mod)
###Plot
boxplot(beta.yr_to_cent_dist_8$mean_rel_distance ~ Pond, beta.yr_to_cent_dist_8, col = pal$col, ylab = "Relative year to centroid distance")
points(aggregate(mean_rel_distance ~ Pond, data = beta.yr_to_cent_dist_8, mean)[2], pch = 16)
##Type
###Test
summary(aov.mod <- aov(mean_rel_distance ~ type, beta.yr_to_cent_dist_8))
TukeyHSD(aov.mod)
###Plot
boxplot(beta.yr_to_cent_dist$mean_rel_distance ~ type, beta.yr_to_cent_dist_8,
        col = c("#506E82", "#768655"), ylab = "Relative year to centroid distance",
        names = c("Species", "Functional Traits"), xlab = "")
points(aggregate(mean_rel_distance ~ type, data = beta.yr_to_cent_dist_8, mean)[2], pch = 16)


#yr to yr
avg.yr_yr_dist_FT <- MCMC_yr_yr_dist_FT_8 %>%
  group_by(Pond, Year) %>%
  summarise(mean_rel_distance = mean(rel_distance), mean_distance = mean(distance))%>%
  mutate(type = "traits")


avg.yr_yr_dist_COM <- MCMC_yr_yr_dist_COM_8 %>%
  group_by(Pond, Year) %>%
  summarise(mean_rel_distance = mean(rel_distance), mean_distance = mean(distance))%>%
  mutate(type = "species")


avg.yr_yr_dist <- rbind(avg.yr_yr_dist_COM, avg.yr_yr_dist_FT)

beta.yr_to_yr_dist_8 <- merge(avg.yr_yr_dist, pal, by = "Pond")
summary(aov(mean_rel_distance ~ Pond * type, beta.yr_to_yr_dist))
summary(aov(mean_rel_distance ~ Pond, beta.yr_to_yr_dist_8))
summary(aov(mean_rel_distance ~ type, beta.yr_to_yr_dist_8))
mod <- aov(mean_rel_distance ~ Pond, beta.yr_to_yr_dist_8)
TukeyHSD(mod)
##Pond
###Test
summary(aov(mean_rel_distance ~ Pond, beta.yr_to_yr_dist_8))
###Plot
boxplot(beta.yr_to_yr_dist_8$mean_rel_distance ~ Pond, beta.yr_to_yr_dist_8, col = pal$col, ylab = "Relative year to year distance")
points(aggregate(mean_rel_distance ~ Pond, data = beta.yr_to_yr_dist_8, mean)[2], pch = 16)

###Test
summary(aov.mod <- aov(mean_rel_distance ~ type, beta.yr_to_yr_dist_8))
TukeyHSD(aov.mod)
yr.yr.ty.pval <- round(summary(aov.mod)[[1]][["Pr(>F)"]][1],3)
boxplot(beta.yr_to_yr_dist_8$mean_rel_distance ~ type, beta.yr_to_yr_dist_8,
        col = c("#506E82", "#768655"), ylab = "Relative year to year distance",
        names = c("Species", "Functional Traits"), xlab = "")
points(aggregate(mean_rel_distance ~ type, data = beta.yr_to_yr_dist_8, mean)[2], pch = 16)


## species assemblages

spc_assemblges_8 <- as.data.frame(spc_assemblges)
spc_counts_MCMC <- spc_assemblges_8 %>%
  group_by(spc.smp)%>%
  summarise(prop_of_iter = length(j)/10000)









###Figure for supp

###yr-cent-13
#png("/Users/mallarie.yeager/Documents/R_projects/Functional_Trait_Stability/02_Ordination Analysis/Plots/deminsional_reduction_test_plot.png",
 #   width = 1200, height = 600, units = "px")
par(mfrow = c(2,4))
boxplot(beta.yr_to_cent_dist_13$mean_rel_distance ~ Pond, beta.yr_to_cent_dist_13, col = pal$col, ylab = "Relative year to centroid distance")
points(aggregate(mean_rel_distance ~ Pond, data = beta.yr_to_cent_dist_13, mean)[2], pch = 16)
###Plot
boxplot(beta.yr_to_cent_dist_13$mean_rel_distance ~ type, beta.yr_to_cent_dist_13,
        col = c("#506E82", "#768655"), ylab = "Relative year to centroid distance",
        names = c("Species", "Functional Traits"), xlab = "")
points(aggregate(mean_rel_distance ~ type, data = beta.yr_to_cent_dist_13, mean)[2], pch = 16)


#yr to yr- 13
boxplot(beta.yr_to_yr_dist_13$mean_rel_distance ~ Pond, beta.yr_to_yr_dist_13, col = pal$col, ylab = "Relative year to year distance")
points(aggregate(mean_rel_distance ~ Pond, data = beta.yr_to_yr_dist_13, mean)[2], pch = 16)

boxplot(beta.yr_to_yr_dist_13$mean_rel_distance ~ type, beta.yr_to_yr_dist_13,
        col = c("#506E82", "#768655"), ylab = "Relative year to year distance",
        names = c("Species", "Functional Traits"), xlab = "")
points(aggregate(mean_rel_distance ~ type, data = beta.yr_to_yr_dist_13, mean)[2], pch = 16)



###yr-cent-8
boxplot(beta.yr_to_cent_dist_8$mean_rel_distance ~ Pond, beta.yr_to_cent_dist_8, col = pal$col, ylab = "Relative year to centroid distance")
points(aggregate(mean_rel_distance ~ Pond, data = beta.yr_to_cent_dist_8, mean)[2], pch = 16)
###Plot
boxplot(beta.yr_to_cent_dist_8$mean_rel_distance ~ type, beta.yr_to_cent_dist_8,
        col = c("#506E82", "#768655"), ylab = "Relative year to centroid distance",
        names = c("Species", "Functional Traits"), xlab = "")
points(aggregate(mean_rel_distance ~ type, data = beta.yr_to_cent_dist_8, mean)[2], pch = 16)


#yr to yr- 8
boxplot(beta.yr_to_yr_dist_8$mean_rel_distance ~ Pond, beta.yr_to_yr_dist_8, col = pal$col, ylab = "Relative year to year distance")
points(aggregate(mean_rel_distance ~ Pond, data = beta.yr_to_yr_dist_8, mean)[2], pch = 16)

boxplot(beta.yr_to_yr_dist_8$mean_rel_distance ~ type, beta.yr_to_yr_dist_8,
        col = c("#506E82", "#768655"), ylab = "Relative year to year distance",
        names = c("Species", "Functional Traits"), xlab = "")
points(aggregate(mean_rel_distance ~ type, data = beta.yr_to_yr_dist_8, mean)[2], pch = 16)
#dev.off()

