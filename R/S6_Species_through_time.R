library(ggplot2)
library(dplyr)
library(patchwork)

sums <- read.csv("02_Ordination Analyses/Outputs/Species_comm_dat.csv")

col <- pnw_palette("Sunset",100)
col.detrit <- col[c(1,7,13)]

col.meio<- col[c(55,57,61,63,65,67,69,71,73)]
col.macro<- col[c(84,86,88,90,92,94,96,98,100)]

spc.long <- sums %>%
  pivot_longer(3:ncol(sums),
               values_to = "count",
               names_to = "spcs")

ggplot(spc.long, aes(x = as.numeric(Year), y = log10(count), color = spcs))+
  geom_line()+
  facet_wrap(.~Pond, scales = "free_y")
#subset ponds which are "most stable" in the system - NP and PJ
spc.long.NP_PJ <- spc.long %>%
  filter(Pond == "NP" | Pond == "PJ")

###### Planktivores#######
col.plankt <- col[c(12,18,28,38,48,58)]
plank <- subset(spc.long.NP_PJ, spcs == "Atlantic Menhaden" | spcs == "Atlantic Silverside"|
                  spcs == "Alewife" | spcs == "Bay Anchovy" |
                  spcs == "Northern Pipefish" | spcs == "White Perch")
sum.plank <- aggregate(count ~ Pond + Year, plank, sum)

sum.plank$spcs <- "Total sum"
sum.plank <- sum.plank[,c(1:2,4,3)]
plank <- rbind(plank, sum.plank)

plank$spcs <- factor(plank$spcs, levels = c("Total sum","Alewife","Atlantic Menhaden",
                                               "Atlantic Silverside", "Bay Anchovy",
                                               "Northern Pipefish","White Perch"))

(p1 <- ggplot(plank, aes(x = as.numeric(Year), y = log10(count), color = spcs))+
  geom_line(lwd = 1.2)+
  ylab("Counts (log10)")+
  xlab("Year")+
  scale_color_viridis_d(direction = 1, option = "magma", name = "Planktivores")+
  facet_wrap(.~Pond, scales = "free_y")+
  theme_bw())



##### Meiofauna #####
meiofan <- subset(spc.long.NP_PJ, spcs == "Mummichog" | spcs == "Striped Killifish"| spcs == "Atlantic Tomcod"
                  | spcs == "Fourspine Stickleback"| spcs == "Goby"| spcs == "Smallmouth Flounder"| 
                    spcs == "Spotfin Mojarra"| spcs == "Tautog"| spcs == "Winter Flounder")

sum.meiofan <- aggregate(count ~ Pond + Year, meiofan, sum)

sum.meiofan$spcs <- "Total sum"
sum.meiofan <- sum.meiofan[,c(1:2,4,3)]
meiofan <- rbind(meiofan, sum.meiofan)
meiofan$spcs <- factor(meiofan$spcs, levels = c("Total sum", "Atlantic Tomcod","Fourspine Stickleback",
                                                "Goby", "Mummichog", "Smallmouth Flounder",
                                                "Spotfin Mojarra","Striped Killifish","Tautog",
                                                "Winter Flounder"))

(p2 <- ggplot()+
  geom_line(meiofan, mapping = aes(x = as.numeric(Year), y = log10(count), color = spcs),lwd = 1.2)+
  ylab("Counts (log10)")+
  xlab("Year")+
  scale_color_viridis_d(direction = 1, option ="turbo", name = "Hunting Meiofauna")+
  facet_wrap(.~Pond, scales = "free_y")+
  theme_bw())



##### Macrofauna #####
macrofan <- subset(spc.long.NP_PJ, spcs == "Atlantic Needlefish" | spcs == "Black seabass"| spcs == "Crevalle Jack"
                  | spcs == "Grubby"| spcs == "Inshore Lizardfish"| spcs == "Northern Kingfish"| 
                    spcs == "Oyster Toadfish"| spcs == "Scup"| spcs == "Summer Flounder")

sum.macrofan <- aggregate(count ~ Pond + Year, macrofan, sum)

sum.macrofan$spcs <- "Total sum"
sum.macrofan <- sum.macrofan[,c(1:2,4,3)]
macrofan <- rbind(macrofan, sum.macrofan)
macrofan$spcs <- factor(macrofan$spcs, levels = c("Total sum","Atlantic Needlefish", "Black seabass","Crevalle Jack",
                                                  "Grubby", "Inshore Lizardfish", "Northern Kingfish",
                                                   "Oyster Toadfish", "Scup","Summer Flounder"))

(p3 <- ggplot()+
  geom_line(macrofan, mapping = aes(x = as.numeric(Year), y = log10(count), color = spcs),lwd = 1.2)+
  ylab("Counts (log10)")+
  xlab("Year")+
  scale_color_viridis_d(direction = 1, option ="mako", name = "Hunting Macrofauna")+
  facet_wrap(.~Pond, scales = "free_y")+
  theme_bw())


##### Detritivore #####
detrit <- subset(spc.long.NP_PJ, spcs == "Rainwater Killifish" | spcs == "Sheepshead Minnow"|spcs == "White Mullet" )

sum.detrit <- aggregate(count ~ Pond + Year, detrit, sum)

sum.detrit$spcs <- "Total sum"
sum.detrit <- sum.detrit[,c(1:2,4,3)]
detrit <- rbind(detrit, sum.detrit)
detrit$spcs <- factor(detrit$spcs, levels = c("Total sum","Rainwater Killifish", "Sheepshead Minnow","White Mullet"))

(p4 <- ggplot()+
    geom_line(detrit, mapping = aes(x = as.numeric(Year), y = log10(count), color = spcs),lwd = 1.2)+
    ylab("Counts (log10)")+
    xlab("Year")+
    scale_color_viridis_d(direction = 1, option ="cividis", name = "Detritivore")+
    facet_wrap(.~Pond, scales = "free_y")+
    theme_bw())



(plot <-(p4+p1)/(p2+p3))
#ggsave("Tests for reviewers/spc_tru_time.png",plot)
