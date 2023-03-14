

# Loading packages
library(tidyverse)
library(readxl)
library(grid)


# ------------------------------------------------------------------------------
##                          LOADING DATA


# Loading the master data sheet
Master <- read_excel("Murua_and_Gaxiola_MasterDataSheet.xlsx", 
                     sheet = "Master", 
                     range = "A1:S365")


# Setting treatments as factors
Master$Litter_spp <- factor(Master$Litter_spp, 
                     levels = c("Aristotelia", "Nothofagus", "Lomatia"))
Master$Soil_type <- factor(Master$Soil_type, 
                     levels = c("Aristotelia", "Nothofagus", "Lomatia"))
Master$ID <- as.character(Master$ID) # Setting ID as character


# Subset of master datasheet containing only replicates destined to mass loss
Master.m <- Master[Master$Sample_type == "ML", -12:-13]


# Importing litter elemental data measured over time
Sept.18.h <- read_excel("Chemical_analyses/Murua_and_Gaxiola_EA_LeafLitter_09-2018.xls",
                        sheet = 1, 
                        range = "B11:J47")
Jun.19 <- read_excel("Chemical_analyses/Murua_and_Gaxiola_EA_Soil_06-2019.xls",
                     sheet = 1, 
                     range = "B11:L54")
Jun.19 <- Jun.19[, c(-5, -8)]
Ag.19.h <- read_excel("Chemical_analyses/Murua_and_Gaxiola_EA_SoilandLeafLitter_08-2019.xls",
                      sheet = 1, 
                      range = "B27:J163")


# Combining litter elemental data
Elements <- rbind(Sept.18.h, Jun.19, Ag.19.h)
names(Elements) <- c("ID", "sample.mass", "N.mg", "N.ug", "C.mg", "C.ug", 
                     "N.perc", "C.perc", "CN")
Elements$`ID` <- 
  str_split(Elements$`ID`, "AG ", simplify = TRUE)[, 2] %>%
  str_replace_all(" ", "") %>%
  str_replace_all("V2", "") %>%
  str_remove_all("J")


# Combining master data set with litter elemental data 
nut <- inner_join(Master.m, Elements, by = "ID")


# Calculating the initial nitrogen for each litter type (harvest 1)
nut2 <- nut

Init.tab2 <-
nut[nut$Harvest == 1,] %>% group_by(Litter_spp) %>%
 summarise(Cp.i = mean(C.perc),
           Np.i = mean(N.perc),
           CNi = mean(CN),
           NCi = mean(N.perc/C.perc))


# Reading enzyme assays
Enzymes <- read_excel("Enzyme_assays/Murua_and_Gaxiola_MasterFluo.xlsx", 
                      sheet = "Master",
                      range = "A1:T1969")
names(Enzymes)[c(1, 3, 4, 10, 11, 15, 16, 19)] <- 
  c("Assay", "Enzyme", "Order", "mean.SC", "mean.SF", "Time", "Mass", "Harvest")


# Excluding null activities
Enzymes <- Enzymes[Enzymes$ACT != 0,]


# Averaging activities
Final <-
  Enzymes %>% group_by(Harvest, IDD, Enzyme) %>%
  summarise(Act = mean(ACT)) %>%
  pivot_wider(names_from = Enzyme, values_from = Act)


# Reading processed enzyme activities
pox<- read.csv("Pox.csv", header = T, sep = ";", dec = ",")


# Renaming variables
names(pox)[c(2, 4, 7, 31, 11)] <- c("Litter_spp", "Soil_type", "Harvest", "CN", "days")


# Merging the data set without 0s and the processed data
pox <- merge(pox[, c("ID", "Litter_spp", "Soil_type", "color", 
                     "Harvest", "CN", "days")],
             Final,                        
             by.x = c("ID", "color"), 
             by.y = c("IDD", "Harvest"))

# Setting the order of treatment levels
pox$Soil_type = factor(pox$Soil_type, levels = c("Aristotelia", "Nothofagus",  "Lomatia"))
pox$Litter_spp = factor(pox$Litter_spp, levels = c("Aristotelia",  "Nothofagus", "Lomatia"))
pox <- pox[order(pox$Harvest, pox$Litter_spp, pox$Soil_type),]


# Including CN measurements for the last harvest
New.cn <- nut[nut$Harvest == 6,] %>% group_by(Litter_spp, Soil_type, Harvest) %>%
 summarise(CN = mean(CN))                 

pox[pox$Harvest == 6,]$CN <- rep(New.cn$CN, each = 3)



# ------------------------------------------------------------------------------
#                        NITROGEN RELEASE CURVES


# Joining the datasets
nut2 <- inner_join(nut2,                             
                   Init.tab2,
                   by = "Litter_spp")


# Calculating absolute masses of C and N
nut2$C <- (nut2$C.perc/100)*nut2$Final_mass
nut2$N <- (nut2$N.perc/100)*nut2$Final_mass


# Absolute mass of initial carbon
nut2$Co <- (nut2$Cp.i/100)*nut2$Initial_mass


# Calculating mean absolute mass per treatment
mean.C <-
  nut2[, c("Litter_spp", "Soil_type", "Harvest", "C")] %>%
  group_by(Litter_spp, Soil_type, Harvest) %>%
  summarise(mean.C = mean(C))


# Calculating remaining carbon
nut2$c <- nut2$C/nut2$Co


# Absolute mass of initial nitrogen
nut2$No <- (nut2$Np.i/100)*nut2$Initial_mass


# Treatment means of absolute mass
mean.N <-
group_by(nut2[, c("N", "Litter_spp", "Soil_type", "Harvest")], 
         Litter_spp, Soil_type, Harvest) %>%
 summarise(mean.N = mean(N))


# Calculating remaining nitrogen
nut2$n <- nut2$N/nut2$No
nut3 <- split(nut2, interaction(nut2$Litter_spp, nut2$Soil_type))
Means2 <- 
  group_by(nut2, Litter_spp, Soil_type, Harvest) %>%
  summarise(c = mean(c), n = mean(n)) %>%
  filter(Harvest == 1)

Max <- sapply(lapply(nut3, "[[", "c"), 
               function(is){which(is == max (is))}) 

Max.df <- data.frame(Litter_spp = rep(levels(nut2$Litter_spp), 3),
                    Soil_type = rep(levels(nut2$Litter_spp), each = 3),
                    NCi = rep(NA, 9), 
                    c = rep(NA, 9), 
                    n = rep(NA, 9))    

for(i in 1:length(Max)){
  Max.df[i,] <- lapply(nut3, "[", c("Litter_spp", "Soil_type", "NCi", "c", "n"))[[i]][Max[i],]
}


# Modifying initial values
Max.df$c.dif <- 1 - Max.df$c
Max.df$n.dif <- 1 - Max.df$n

nut4 <- lapply(nut3, "[", c("Litter_spp", "Soil_type", "NCi", "c", "n"))

for(i in 1:9){
  nut4[[i]]$c <-  nut4[[i]]$c + Max.df[i, "c.dif"]
  nut4[[i]]$n <-  nut4[[i]]$n + Max.df[i,"n.dif"]
}

nut4 <- rbind(nut4[[1]], nut4[[2]], nut4[[3]], nut4[[4]], nut4[[5]], nut4[[6]],
              nut4[[7]], nut4[[8]], nut4[[9]])


# ------------------------------------------------------------------------------

# Fitting the models
Mr2 <- nls(n ~ (c*0.1/NCi) + (1-(0.1/NCi))*c^(1/(1-e)), 
           data = nut4[nut4$Litter_spp=="Aristotelia" & nut4$Soil_type=="Lomatia",], 
           start = list(e = 0.3), 
           trace = TRUE)

Mm2 <- nls(n ~ (c*0.1/NCi) + (1-(0.1/NCi))*c^(1/(1-e)), 
           data = nut4[nut4$Litter_spp=="Aristotelia" & nut4$Soil_type=="Aristotelia",], 
           start = list(e = 0.3), 
           trace = TRUE)

Mc2 <- nls(n ~ (c*0.1/NCi) + (1-(0.1/NCi))*c^(1/(1-e)), 
           data = nut4[nut4$Litter_spp=="Aristotelia" & nut4$Soil_type=="Nothofagus",], 
           start = list(e = 0.3), 
           trace = TRUE)

Cr2 <- nls(n ~ (c*0.1/NCi) + (1-(0.1/NCi))*c^(1/(1-e)), 
           data = nut4[nut4$Litter_spp=="Nothofagus" & nut4$Soil_type=="Lomatia",], 
           start = list(e = 0.3), 
           trace = TRUE)

Cm2 <- nls(n ~ (c*0.1/NCi) + (1-(0.1/NCi))*c^(1/(1-e)), 
           data = nut4[nut4$Litter_spp=="Nothofagus" & nut4$Soil_type=="Aristotelia",], 
           start = list(e = 0.3), 
           trace = TRUE)

Cc2 <- nls(n ~ (c*0.1/NCi) + (1-(0.1/NCi))*c^(1/(1-e)), 
           data = nut4[nut4$Litter_spp=="Nothofagus" & nut4$Soil_type=="Nothofagus",], 
           start = list(e = 0.3), 
           trace = TRUE)

Rr2 <- nls(n ~ (c*0.1/NCi) + (1-(0.1/NCi))*c^(1/(1-e)), 
           data = nut4[nut4$Litter_spp=="Lomatia" & nut4$Soil_type=="Lomatia",], 
           start = list(e = 0.3), 
           trace = TRUE)

Rm2 <- nls(n ~ (c*0.1/NCi) + (1-(0.1/NCi))*c^(1/(1-e)), 
           data = nut4[nut4$Litter_spp=="Lomatia" & nut4$Soil_type=="Aristotelia",], 
           start = list(e = 0.3), 
           trace = TRUE)

Rc2 <- nls(n ~ (c*0.1/NCi) + (1-(0.1/NCi))*c^(1/(1-e)), 
           data = nut4[nut4$Litter_spp=="Lomatia" & nut4$Soil_type=="Nothofagus",], 
           start = list(e = 0.3), 
           trace = TRUE)

xMr2 <- seq(0, 1, length.out = 2500)
xMm2 <- seq(0, 1, length.out = 2300)
xMc2 <- seq(0, 1, length.out = 2500)
xCr2 <- seq(0, 1, length.out = 2500)
xCm2 <- seq(0, 1, length.out = 2000)
xCc2 <- seq(0, 1, length.out = 2000)
xRr2 <- seq(0, 1, length.out = 2400)
xRm2 <- seq(0, 1, length.out = 2600)
xRc2 <- seq(0, 1, length.out = 2600)


yMr2 <-  predict(Mr2, list(c = xMr2))
yMm2 <-  predict(Mm2, list(c = xMm2))
yMc2 <-  predict(Mc2, list(c = xMc2))
yCr2 <-  predict(Cr2, list(c = xCr2))
yCm2 <-  predict(Cm2, list(c = xCm2))
yCc2 <-  predict(Cc2, list(c = xCc2))
yRr2 <-  predict(Rr2, list(c = xRr2))
yRm2 <-  predict(Rm2, list(c = xRm2))
yRc2 <-  predict(Rc2, list(c = xRc2))


DMr2 <- data.frame(carb = (1-xMr2), n = yMr2, Litter_spp = rep("Aristotelia", 2500), Soil_type = rep("Lomatia", 2500))
DMm2 <- data.frame(carb = (1-xMm2), n = yMm2, Litter_spp = rep("Aristotelia", 2300), Soil_type = rep("Aristotelia", 2300))
DMc2 <- data.frame(carb = (1-xMc2), n = yMc2, Litter_spp = rep("Aristotelia", 2500), Soil_type = rep("Nothofagus", 2500))
DCr2 <- data.frame(carb = (1-xCr2), n = yCr2, Litter_spp = rep("Nothofagus", 2500), Soil_type = rep("Lomatia", 2500))
DCm2 <- data.frame(carb = (1-xCm2), n = yCm2, Litter_spp = rep("Nothofagus", 2000), Soil_type = rep("Aristotelia", 2000))
DCc2 <- data.frame(carb = (1-xCc2), n = yCc2, Litter_spp = rep("Nothofagus", 2000), Soil_type = rep("Nothofagus", 2000))
DRr2 <- data.frame(carb = (1-xRr2), n = yRr2, Litter_spp = rep("Lomatia", 2400), Soil_type = rep("Lomatia", 2400))
DRm2 <- data.frame(carb = (1-xRm2), n = yRm2, Litter_spp = rep("Lomatia", 2600), Soil_type = rep("Aristotelia", 2600))
DRc2 <- data.frame(carb = (1-xRc2), n = yRc2, Litter_spp = rep("Lomatia", 2600), Soil_type = rep("Nothofagus", 2600))

Opt.Mr2 <- data.frame(carb = (1-xMr2))
Opt.Mm2 <- data.frame(carb = (1-xMm2))
Opt.Mc2 <- data.frame(carb = (1-xMc2))
Opt.Cr2 <- data.frame(carb = (1-xCr2))
Opt.Cm2 <- data.frame(carb = (1-xCm2))
Opt.Cc2 <- data.frame(carb = (1-xCc2))
Opt.Rr2 <- data.frame(carb = (1-xRr2))
Opt.Rm2 <- data.frame(carb = (1-xRm2))
Opt.Rc2 <- data.frame(carb = (1-xRc2))

Opt.Mr2$zoil = factor(DMr2$Soil_type, c("Aristotelia", "Nothofagus", "Lomatia"))
Opt.Mm2$zoil = factor(DMm2$Soil_type, c("Aristotelia", "Nothofagus", "Lomatia"))
Opt.Mc2$zoil = factor(DMc2$Soil_type, c("Aristotelia", "Nothofagus", "Lomatia"))
Opt.Cr2$zoil = factor(DCr2$Soil_type, c("Aristotelia", "Nothofagus", "Lomatia"))
Opt.Cm2$zoil = factor(DCm2$Soil_type, c("Aristotelia", "Nothofagus", "Lomatia"))
Opt.Cc2$zoil = factor(DCc2$Soil_type, c("Aristotelia", "Nothofagus", "Lomatia"))
Opt.Rr2$zoil = factor(DRr2$Soil_type, c("Aristotelia", "Nothofagus", "Lomatia"))
Opt.Rm2$zoil = factor(DRm2$Soil_type, c("Aristotelia", "Nothofagus", "Lomatia"))
Opt.Rc2$zoil = factor(DRc2$Soil_type, c("Aristotelia", "Nothofagus", "Lomatia"))

nut4$zoil = factor(nut4$Soil_type, c("Aristotelia", "Nothofagus", "Lomatia"))
DMr2$zoil = factor(DMr2$Soil_type, c("Aristotelia", "Nothofagus", "Lomatia"))
DMm2$zoil = factor(DMm2$Soil_type, c("Aristotelia", "Nothofagus", "Lomatia"))
DMc2$zoil = factor(DMc2$Soil_type, c("Aristotelia", "Nothofagus", "Lomatia"))
DCr2$zoil = factor(DCr2$Soil_type, c("Aristotelia", "Nothofagus", "Lomatia"))
DCm2$zoil = factor(DCm2$Soil_type, c("Aristotelia", "Nothofagus", "Lomatia"))
DCc2$zoil = factor(DCc2$Soil_type, c("Aristotelia", "Nothofagus", "Lomatia"))
DRr2$zoil = factor(DRr2$Soil_type, c("Aristotelia", "Nothofagus", "Lomatia"))
DRm2$zoil = factor(DRm2$Soil_type, c("Aristotelia", "Nothofagus", "Lomatia"))
DRc2$zoil = factor(DRc2$Soil_type, c("Aristotelia", "Nothofagus", "Lomatia"))

Opt.Mr2$zpp = factor(DMr2$Litter_spp, c("Aristotelia", "Nothofagus", "Lomatia"))
Opt.Mm2$zpp = factor(DMm2$Litter_spp, c("Aristotelia", "Nothofagus", "Lomatia"))
Opt.Mc2$zpp = factor(DMc2$Litter_spp, c("Aristotelia", "Nothofagus", "Lomatia"))
Opt.Cr2$zpp = factor(DCr2$Litter_spp, c("Aristotelia", "Nothofagus", "Lomatia"))
Opt.Cm2$zpp = factor(DCm2$Litter_spp, c("Aristotelia", "Nothofagus", "Lomatia"))
Opt.Cc2$zpp = factor(DCc2$Litter_spp, c("Aristotelia", "Nothofagus", "Lomatia"))
Opt.Rr2$zpp = factor(DRr2$Litter_spp, c("Aristotelia", "Nothofagus", "Lomatia"))
Opt.Rm2$zpp = factor(DRm2$Litter_spp, c("Aristotelia", "Nothofagus", "Lomatia"))
Opt.Rc2$zpp = factor(DRc2$Litter_spp, c("Aristotelia", "Nothofagus", "Lomatia"))

nut4$zpp = factor(nut4$Litter_spp, c("Aristotelia", "Nothofagus", "Lomatia"))
DMr2$zpp = factor(DMr2$Litter_spp, c("Aristotelia", "Nothofagus", "Lomatia"))
DMm2$zpp = factor(DMm2$Litter_spp, c("Aristotelia", "Nothofagus", "Lomatia"))
DMc2$zpp = factor(DMc2$Litter_spp, c("Aristotelia", "Nothofagus", "Lomatia"))
DCr2$zpp = factor(DCr2$Litter_spp, c("Aristotelia", "Nothofagus", "Lomatia"))
DCm2$zpp = factor(DCm2$Litter_spp, c("Aristotelia", "Nothofagus", "Lomatia"))
DCc2$zpp = factor(DCc2$Litter_spp, c("Aristotelia", "Nothofagus", "Lomatia"))
DRr2$zpp = factor(DRr2$Litter_spp, c("Aristotelia", "Nothofagus", "Lomatia"))
DRm2$zpp = factor(DRm2$Litter_spp, c("Aristotelia", "Nothofagus", "Lomatia"))
DRc2$zpp = factor(DRc2$Litter_spp, c("Aristotelia", "Nothofagus", "Lomatia"))


# Function to make theoretical curve
M.fun <- function(c, L, e){
  (c*0.1/L) + (1-(0.1/L))*c^(1/(1-e))
}

Opt.Mr2$n <- M.fun(c = Opt.Mr2$carb, L = 0.019, e = 0.197)
Opt.Mm2$n <- M.fun(c = Opt.Mm2$carb, L = 0.019, e = 0.197)
Opt.Mc2$n <- M.fun(c = Opt.Mc2$carb, L = 0.019, e = 0.197)
Opt.Cr2$n <- M.fun(c = Opt.Cr2$carb, L = 0.013, e = 0.130)
Opt.Cm2$n <- M.fun(c = Opt.Cm2$carb, L = 0.013, e = 0.130)
Opt.Cc2$n <- M.fun(c = Opt.Cc2$carb, L = 0.013, e = 0.130)
Opt.Rr2$n <- M.fun(c = Opt.Rr2$carb, L = 0.011, e = 0.111)
Opt.Rm2$n <- M.fun(c = Opt.Rm2$carb, L = 0.011, e = 0.111)
Opt.Rc2$n <- M.fun(c = Opt.Rc2$carb, L = 0.011, e = 0.111)



#Curves by litter species and soil type

Int <- interaction(nut4$Litter_spp, nut4$Soil_type)

Grids <- LETTERS[1:9] # In capitals
names(Grids) <- unique(Int)
nut4$Grid <- unname(Grids[Int])

grids <- letters[1:9] # Lower case
names(grids) <- unique(Int)
nut4$grid <- unname(grids[Int])

Lines <- c("Fitted" = "solid", "Expected" = "dashed")

pdf("Fig. 2.pdf", height = 6, width = 8)
ggplot() + geom_point(data = nut4, mapping = aes(x = (1-c)*100, y = n*100)) + 
  geom_line(data = DMr2, mapping = aes(x = carb*100, y = n*100, linetype = "Fitted")) +
  geom_line(data = DMm2, mapping = aes(x = carb*100, y = n*100, linetype = "Fitted")) + 
  geom_line(data = DMc2, mapping = aes(x = carb*100, y = n*100, linetype = "Fitted")) + 
  geom_line(data = DCr2, mapping = aes(x = carb*100, y = n*100, linetype = "Fitted")) + 
  geom_line(data = DCm2, mapping = aes(x = carb*100, y = n*100, linetype = "Fitted")) + 
  geom_line(data = DCc2, mapping = aes(x = carb*100, y = n*100, linetype = "Fitted")) + 
  geom_line(data = DRr2, mapping = aes(x = carb*100, y = n*100, linetype = "Fitted")) + 
  geom_line(data = DRm2, mapping = aes(x = carb*100, y = n*100, linetype = "Fitted")) + 
  geom_line(data = DRc2, mapping = aes(x = carb*100, y = n*100, linetype = "Fitted")) + 
  geom_line(data = Opt.Mr2, mapping = aes(x = (1 - carb)*100, y = n*100, linetype = "Expected"), alpha = 0.8, ) +
  geom_line(data = Opt.Mm2, mapping = aes(x = (1 - carb)*100, y = n*100, linetype = "Expected"), alpha = 0.8) +
  geom_line(data = Opt.Mc2, mapping = aes(x = (1 - carb)*100, y = n*100, linetype = "Expected"), alpha = 0.8) +
  geom_line(data = Opt.Cr2, mapping = aes(x = (1 - carb)*100, y = n*100, linetype = "Expected"), alpha = 0.8) +
  geom_line(data = Opt.Cm2, mapping = aes(x = (1 - carb)*100, y = n*100, linetype = "Expected"), alpha = 0.8) +
  geom_line(data = Opt.Cc2, mapping = aes(x = (1 - carb)*100, y = n*100, linetype = "Expected"), alpha = 0.8) +
  geom_line(data = Opt.Rr2, mapping = aes(x = (1 - carb)*100, y = n*100, linetype = "Expected"), alpha = 0.8) +
  geom_line(data = Opt.Rm2, mapping = aes(x = (1 - carb)*100, y = n*100, linetype = "Expected"), alpha = 0.8) +
  geom_line(data = Opt.Rc2, mapping = aes(x = (1 - carb)*100, y = n*100, linetype = "Expected"), alpha = 0.8) +
  geom_hline(alpha = 0.5, yintercept = 100, linetype = "dotted") +
  facet_grid(zoil~zpp) + expand_limits(y = 0) + 
  scale_linetype_manual(values = Lines) +
  geom_text(x = 0, y = 190, aes(label = grid), data = nut4) +
  labs(colour = "Soil type", x = "C mass loss (%)", y = "Remaining N (%)",
       subtitle = "Litter species") +
  guides(linetype = guide_legend(title = NULL)) +
  
  theme_classic() + 
  theme(strip.text = element_text(size = 10, face = "italic"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.background = element_rect(fill = "transparent"),
        legend.position = c(0.15, 0.67),
        plot.margin = unit(c(0.1,1,0.2,0.2), "cm")) 
grid.text(label = "Soil type", 
          x = unit(0.97, "npc"), y = 0.85, rot = 270, gp=gpar(fontsize=11))
dev.off()


# Table with Manzoni's parameters
sMr <- summary(Mr2)

models <- list(Mm2, Mc2, Mr2, Cm2, Cc2, Cr2, Rm2, Rc2, Rr2)
Coefs <- sapply(models, coef)
Sum.M <- lapply(models, summary)
Sum.M2 <- t(sapply(Sum.M, "[[", 10))

NCo2 <- aggregate(nut4$NCi, 
                  list(Litter_spp = nut4$Litter_spp, 
                       Soil_type = nut4$Soil_type), 
                  mean)


N_Imm_SE <- c((max(M.fun(c = Opt.Mm2$carb, 
                         L = unlist(Init.tab2[1, 5]), 
                         e = Sum.M2[1,1] + Sum.M2[1,2])) - 
                 max(M.fun(c = Opt.Mm2$carb, 
                           L = unlist(Init.tab2[1, 5]), 
                           e = Sum.M2[1,1] - Sum.M2[1,2])))*100,
              
              (max(M.fun(c = Opt.Mm2$carb, 
                         L = unlist(Init.tab2[1, 5]), 
                         e = Sum.M2[2,1] + Sum.M2[2,2])) -
                 max(M.fun(c = Opt.Mm2$carb, 
                           L = unlist(Init.tab2[1, 5]), 
                           e = Sum.M2[2,1] - Sum.M2[2,2])))*100,
              
              (max(M.fun(c = Opt.Mm2$carb, 
                         L = unlist(Init.tab2[1, 5]), 
                         e = Sum.M2[3,1] + Sum.M2[3,2])) -
                 max(M.fun(c = Opt.Mm2$carb, 
                           L = unlist(Init.tab2[1, 5]), 
                           e = Sum.M2[3,1] - Sum.M2[3,2])))*100,
              
              (max(M.fun(c = Opt.Mm2$carb, 
                         L = unlist(Init.tab2[3, 5]), 
                         e = Sum.M2[4,1] + Sum.M2[4,2])) -
                 max(M.fun(c = Opt.Mm2$carb, 
                           L = unlist(Init.tab2[3, 5]), 
                           e = Sum.M2[4,1] - Sum.M2[4,2])))*100,
              
              (max(M.fun(c = Opt.Mm2$carb, 
                         L = unlist(Init.tab2[3, 5]), 
                         e = Sum.M2[5,1] + Sum.M2[5,2])) -
                 max(M.fun(c = Opt.Mm2$carb, 
                           L = unlist(Init.tab2[3, 5]), 
                           e = Sum.M2[5,1] - Sum.M2[5,2])))*100,
              
              (max(M.fun(c = Opt.Mm2$carb, 
                         L = unlist(Init.tab2[3, 5]), 
                         e = Sum.M2[6,1] + Sum.M2[6,2])) -
                 max(M.fun(c = Opt.Mm2$carb, 
                           L = unlist(Init.tab2[3, 5]), 
                           e = Sum.M2[6,1] - Sum.M2[6,2])))*100,
              
              (max(M.fun(c = Opt.Mm2$carb, 
                         L = unlist(Init.tab2[2, 5]), 
                         e = Sum.M2[7,1] + Sum.M2[7,2])) -
                 max(M.fun(c = Opt.Mm2$carb, 
                           L = unlist(Init.tab2[2, 5]), 
                           e = Sum.M2[7,1] - Sum.M2[7,2])))*100,
              
              (max(M.fun(c = Opt.Mm2$carb, 
                         L = unlist(Init.tab2[2, 5]), 
                         e = Sum.M2[8,1] + Sum.M2[8,2])) -
                 max(M.fun(c = Opt.Mm2$carb, 
                           L = unlist(Init.tab2[2, 5]), 
                           e = Sum.M2[8,1] - Sum.M2[8,2])))*100,
              
              (max(M.fun(c = Opt.Mm2$carb, 
                         L = unlist(Init.tab2[2, 5]), 
                         e = Sum.M2[9,1] + Sum.M2[9,2])) -
                 max(M.fun(c = Opt.Mm2$carb, 
                           L = unlist(Init.tab2[2, 5]), 
                           e = Sum.M2[9,1] - Sum.M2[9,2])))*100)



DatM <- data.frame(Litter_spp = rep(c("Aristotelia", "Nothofagus", "Lomatia"), 
                             each = 3), 
                   Soil_type = rep(c("Aristotelia", "Nothofagus", "Lomatia"), 
                               3), 
                   N_immobilized = c(max(yMm2), max(yMc2), max(yMr2), 
                                     max(yCm2), max(yCc2), max(yCr2), 
                                     max(yRm2), max(yRc2), max(yRr2))*100 - 100, 
                   N_imm_SE = N_Imm_SE,
                   e = Coefs, 
                   e_SE = Sum.M2[,2], 
                   rcM = 1/(Coefs*0.1), 
                   CN.1 = 1/NCo2$x, 
                   e.esp = 8.6/(1/NCo2$x))


DatM$e.dif <- DatM$e.esp - DatM$e
DatM$Litter_spp <- factor(DatM$Litter_spp, levels = c("Aristotelia", "Nothofagus", "Lomatia"))
DatM$Soil_type <- factor(DatM$Soil_type, levels = c("Aristotelia", "Nothofagus", "Lomatia"))


 
 
 
 
 
# ------------------------------------------------------------------------------
#                         THRESHOLD ELEMENT RATIO
 
 
 # Calculating C:N threshold element ratio (TER)
 pox$TERcn <- (pox$BG/pox$NAG)*pox$CN  
 
 
 # Ranking TER values to identify outliers
 pox[order(pox$TERcn, decreasing = TRUE),    # Samples 71, 170, and 155 have extremely high values
     c("ID", "BG","NAG","TERcn")]            
 
 # Excluding samples 71, 170 y 155
 pox <- pox[pox$ID != 71 & pox$ID != 170 & pox$ID != 155,]
 
 
 #  Plotting threshold element ratio by litter species and soil type
 Int <- interaction(pox$Litter_spp, pox$Soil_type)
 grids <- letters[1:9]
 names(grids) <- unique(Int)
 pox$grid <- unname(grids[Int])
 
 Lines <- c("TER" = "solid", "Litter C:N" = "dashed")
 
 pdf("Fig. 3.pdf", height = 6, width = 8)
 ggplot(data = pox) + 
   stat_summary(aes(x = days, y = TERcn), fun.data = "mean_se", geom = "point") +
   stat_summary(aes(x = days, y = TERcn), fun.data = "mean_se", geom = "errorbar", width = 4) +
   stat_summary(aes(x = days, y = TERcn, linetype = "TER"), fun.data = "mean_se", geom = "line") + 
   stat_summary(aes(x = days, y = CN, group = interaction(Litter_spp, Soil_type), linetype = "Litter C:N"), fun.data = "mean_se", geom = "line", alpha = 0.5) + 
   facet_grid(Soil_type ~ Litter_spp) + 
   theme_classic() + 
   labs(x = "Time (days)", y = "TER", colour = "Soil type", subtitle = "Litter species") +
   theme(strip.text = element_text(size = 10, face = "italic"),
         legend.text = element_text(size = 10, face = "italic"),
         legend.position = c(0.15, 0.91),
         legend.background = element_rect(fill = "transparent"),
         plot.margin = unit(c(0.1,1,0.2,0.2), "cm")) +
   coord_cartesian(ylim = c(0, 300)) +
   geom_text(x = 30, y = 290, aes(label = grid)) +
   scale_linetype_manual(values = Lines) +
   guides(linetype = guide_legend(title = NULL)) 
 
 grid.text(label = "Soil type", 
           x = unit(0.97, "npc"), y = 0.85, rot = 270, gp = gpar(fontsize = 11)) 
 dev.off()
 
 
# -----------------------------------------------------------------------------
#                                 MASS LOSS
 
 
jos<- read.csv("Jos.csv", header = T, sep = ";", dec = ",")

names(jos)[c(2, 3, 4, 5, 6, 9, 10, 12, 14)] <- 
  c("Litter_spp", "spp", "Initial_mass", "Final_mass", "Soil_type", "days", "years", "Harvest", "Remaining_mass")

jos$Litter_spp = factor(jos$Litter_spp, c("Aristotelia", "Nothofagus", "Lomatia"))
jos$Soil_type = factor(jos$Soil_type, c("Aristotelia", "Nothofagus", "Lomatia"))

jos$spp = factor(jos$Litter_spp, c("Aristotelia", "Nothofagus", "Lomatia"))
jos$soil = factor(jos$Soil_type, c("Aristotelia", "Nothofagus", "Lomatia"))

joz <- jos

colnames(joz)[c(3, 7)] <- c("zpp", "zoil")   


cb_palette <- c("#000000", "#E69F00", "#56B4E9") 


#Remaining mass by litter species. Color is soil type
panels <- letters[1:3]
names(panels) <- unique(jos$Litter_spp)
jos$panel <- unname(panels[jos$Litter_spp])


# pdf("Fig. S1.pdf", height = 4, width = 10)
ggplot(data = jos) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 5, data = jos,
               linewidth = 0.4, mapping = aes(x = days, y = Remaining_mass*100, group = interaction(Soil_type, Litter_spp), color = soil)) + 
  stat_summary(fun.data = "mean_se", geom = "line", data = jos, 
               mapping = aes(x = days, y = Remaining_mass*100, group = interaction(Soil_type, Litter_spp), color = soil)) +
  stat_summary(fun.data = "mean_se", geom = "point", size = 2, data = joz, 
               mapping = aes(x = days, y = Remaining_mass*100, group = interaction(zoil, Litter_spp)), alpha = 0.1) +
  stat_summary(fun.data = "mean_se", geom = "line", data = joz, 
               mapping = aes(x = days, y = Remaining_mass*100, group = interaction(zoil, Litter_spp)), alpha = 0.1) + 
  stat_summary(fun.data = "mean_se",  geom = "point", size = 2, data = jos, 
               mapping = aes(x = days, y = Remaining_mass*100, colour = soil)) +
  facet_grid(~spp) + 
  theme_classic() + 
  geom_text(x = 30, y = 100, aes(label = panel)) +
  labs(y = "Mass remaining (%)", x = "Time (days)", subtitle =  "Litter species", colour = "Soil type")  +
  theme(legend.text = element_text(face = "italic", size = 8), legend.position = c(0.94, 0.81), 
        legend.title = element_text(size = 8), strip.text = element_text(size = 10, face = "italic"),) +
  scale_colour_manual(values = cb_palette)
# dev.off()


# Fitting decay models

# Single exponential model by litter species and soil type
SExp <- nls(Remaining_mass ~ exp(-k[interaction(Litter_spp, Soil_type)]*years), 
            data = jos, 
            start = list(k = rep(1, 9)), 
            trace = TRUE)

x <- seq(0, 0.745, length.out = 245)

dat <- expand.grid(years = x, 
                   Litter_spp = levels(jos$Litter_spp), 
                   Soil_type = levels(jos$Soil_type))

nModel <- predict(SExp, dat)
pred <- cbind(dat, nModel)


# Single exponential model excluding first two months
model2 <- nls(Remaining_mass ~ n[interaction(Litter_spp, Soil_type)]*exp(-k[interaction(Litter_spp, Soil_type)]*years), 
              data = jos[jos$Harvest > 2,], start = list(k = rep(1, 9), n = rep(1, 9)), 
              trace = TRUE)
nModel2 <- predict(model2, dat)
pred2 <- cbind(dat, nModel2)


# Asymptotic model (Wider and Lang, 1982)
# The result is identical as if using Berg and Ekbohm (1991)
Asym <- nls(Remaining_mass ~ C[interaction(Litter_spp, Soil_type)] + 
              (1 - C[interaction(Litter_spp, Soil_type)])*exp(-k[interaction(Litter_spp, Soil_type)]*years), 
            data = jos, 
            start = list(C = rep(0.4, 9), 
                         k = rep(3, 9)), 
            trace = TRUE)

nAsym <- predict(Asym, dat)
pred3 <- cbind(dat, nAsym)


# Data frame with parameters of all models
options(digits = 3)
models1 <- list(summary(SExp), summary(model2), summary(Asym))
est <- list(lapply(models1[1], coef), lapply(models1[2], coef), lapply(models1[3], coef))
fs <- expand.grid(Litter_spp = levels(jos$Litter_spp), Soil_type = levels(jos$Soil_type))
Dat <- data.frame(Litter_spp = fs[,1], Soil_type = fs[,2], k = est[[1]][[1]][,1],
                  k.init = est[[3]][[1]][10:18, 1],
                  k.2 = est[[2]][[1]][1:9,1], SE.k = est[[1]][[1]][,2],
                  SE.kinit = est[[3]][[1]][10:18, 2],
                  SE.k2 = est[[2]][[1]][1:9,2], Asym = est[[3]][[1]][1:9, 1], 
                  SE.Asym = est[[3]][[1]][1:9, 2])

Dat$Litter_spp <- factor(Dat$Litter_spp, levels = c("Aristotelia", "Nothofagus", "Lomatia"))
Dat$Soil_type <- factor(Dat$Soil_type, levels = c("Aristotelia", "Nothofagus", "Lomatia"))
Dat <- Dat[order(Dat$Litter_spp, Dat$Soil_type),]

Dat.t <- cbind(Dat, Tmed = (log(2)/Dat$k), Tmed.2 = (log(2)/Dat$k.2), 
               SE.tm = (log(2)/(Dat$k - Dat$SE.k)) - (log(2)/(Dat$k + Dat$SE.k)))
Dat.o <- arrange(Dat.t, Litter_spp, Soil_type)
Tmed <- data.frame(y = 0.5, x = Dat.t$Tmed)
Tmed.2 <- data.frame(y = 0.5, x = Dat.t$Tmed.2)

Dat$CN <- rep(c(13.3, 33.0, 63.1), each = 3)
Dat <- merge(Dat, DatM, by = c("Litter_spp", "Soil_type"), sort = FALSE)



#Dat$CN <- c(rep(18.6, 3), 
#            rep(32.8, 3), 
#            rep(60.6, 3))
#Dat$k <- c(2.501, 1.518, 2.28, 1.515, 1.142, 1.343, 0.827, 0.690, 0.856)


# Plotting decomposition rate (k) vs initial litter C:N
# pdf("Fig. 1.pdf", height = 6, width = 8)
ggplot(data = Dat, mapping = aes(x = CN, y = k)) + 
  geom_errorbar(aes(ymin = k - SE.k, ymax = k + SE.k, group = interaction(Soil_type, Litter_spp)), width = 1) +
  geom_point(aes(fill = Soil_type, shape = Litter_spp), size = 4) +
  theme_classic() +
  scale_y_continuous(breaks = c(1, 1.5, 2, 2.5)) +
  labs(x = "Litter C:N", y = expression(paste("k  ", (year^-1))), shape = "Litter species", fill = "Soil type") +
  scale_shape_manual(values=c(23, 24, 22)) +
  scale_fill_manual(values = c("white", "black", "gray"), guide = guide_legend(override.aes = list(shape = 21))) +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 18, face = "bold", hjust = 1),
        legend.text = element_text(size = 11, face = "italic"),
        legend.title = element_text(size = 12),
        legend.position = c(0.9, 0.75),
        panel.background = element_rect(colour = "black", linewidth = 0.7), 
        axis.line = element_line(NULL)) 
# dev.off()


# ------------------------------------------------------------------------------
#                        ENZYMATC ALLOCATION


# Joining activities in a single column
pox2 <- pivot_longer(data = pox, cols = c("BG", "NAG", "AP"), 
                     names_to = "Enzyme", values_to = "Activity") 

pox2$Litter_spp <- factor(pox2$Litter_spp, levels = c("Aristotelia", "Nothofagus", "Lomatia"))
pox2$Soil_type <- factor(pox2$Soil_type, levels = c("Aristotelia", "Nothofagus", "Lomatia"))

# Enzyme activity per litter x soil combination
Int <- interaction(pox2$Litter_spp, pox2$Soil_type)
grids <- letters[1:9]
names(grids) <- unique(Int)
pox2$grid <- unname(grids[Int])

pdf("Fig. S2.pdf", height = 6, width = 8)
ggplot(data = pox2, mapping = aes(x = days, y = Activity)) + 
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 4, aes(colour = Enzyme)) +
  stat_summary(fun.data = "mean_se", geom = "point", aes(colour = Enzyme)) + 
  stat_summary(fun.data = "mean_se", geom = "line", aes(colour = Enzyme)) + 
  labs(x = "Time (days)", 
       y = "Activity (nmol/h*g litter)", 
       colour = "Enzyme",
       subtitle = "Litter species") +
  scale_color_manual(values = cb_palette) +
  theme_classic() + 
  expand_limits( x = 0, y = c(0, 14)) +
  geom_text(x = 10, y = 5.6*10^6, aes(label = grid)) +
  theme(legend.position = c(0.9, 0.87),
        plot.margin = unit(c(0.1,1,0.2,0.2), "cm"),
        strip.text = element_text(size = 10, 
                                  face = "italic")) +
  facet_grid(Soil_type ~ Litter_spp)  

grid.text(label = "Soil type", 
          x = unit(0.97, "npc"), y = 0.85, rot = 270, gp=gpar(fontsize=11))
dev.off()


# log enzyme activity per litter x soil combination
pdf("Fig. S3.pdf", height = 6, width = 8)
ggplot(data = pox2, mapping = aes(x = days, y = log(Activity))) + 
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 4, aes(colour = Enzyme)) +
  stat_summary(fun.data = "mean_se", geom = "point", aes(colour = Enzyme)) + 
  stat_summary(fun.data = "mean_se", geom = "line", aes(colour = Enzyme)) + 
  labs(x = "Time (days)", 
       y = "ln Activity (nmol/h*g litter)", 
       colour = "Enzyme",
       subtitle = "Litter species") +
  scale_color_manual(values = cb_palette) +
  theme_classic() + 
  expand_limits( x = 0, y = c(0, 14)) +
  geom_text(x = 10, y = 15, aes(label = grid)) +
  theme(legend.position = c(0.9, 0.75),
        plot.margin = unit(c(0.1,1,0.2,0.2), "cm"),
        strip.text = element_text(size = 10, 
                                  face = "italic")) +
  facet_grid(Soil_type ~ Litter_spp)  
grid.text(label = "Soil type", 
          x = unit(0.97, "npc"), y = 0.85, rot = 270, gp=gpar(fontsize=11))
dev.off()



# BG/NAG in time by soil
grids <- letters[1:3]
names(grids) <- unique(pox$Soil_type)
pox$grid <- unname(grids[pox$Soil_type])

# pdf("Fig. S4.pdf", height = 4, width = 10)
ggplot(pox, aes(x = days, y = BG/NAG)) +
  geom_hline(yintercept = 1, alpha = 0.5, linetype = "dashed") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 10, aes(color = Litter_spp)) +
  stat_summary(fun.data = "mean_se", geom = "line", aes(color = Litter_spp)) +
  stat_summary(fun.data = "mean_se", geom = "point", aes(color = Litter_spp)) +
  facet_grid(~ Soil_type) + 
  geom_text(x = 30, y = 9.5, aes(label = grid)) +
  scale_color_manual(values = cb_palette) +
  labs(x = "Time (days)", y = "BG:NAG", color = "Litter species", subtitle = "Soil type") +
  theme_classic() +
  theme(legend.position = c(0.94, 0.8),
        strip.text = element_text(size = 10, face = "italic"),
        legend.text = element_text(size = 10, face = "italic")) 
# + ylim(c(0, 3.7))
# dev.off()


# Patterns of eznymatic allocation. Panels are soil type
ggplot(pox, aes(x = NAG/AP, y = BG/NAG)) +
  geom_point(aes(color = Litter_spp)) + 
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) + 
  geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.5) + 
  ylim(c(0, 6)) + xlim(0, 15) +
  facet_grid(~ Soil_type) +
  scale_color_manual(values = cb_palette) +
  theme_classic()


# Patterns of ezymatic allocation. Panels are litter spp
ggplot(pox, aes(x = NAG/AP, y = BG/NAG)) +
  geom_point(aes(color = Soil_type)) + 
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) + 
  geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.5) + 
  ylim(c(0, 6)) + xlim(0, 15) +
  facet_grid(~ Litter_spp) +
  scale_color_manual(values = cb_palette) +
  theme_classic()


# BG:NAG effect size relative to time 0
pox$BG_NAG <- pox$BG/pox$NAG

BG_NAG_t0 <- 
pox %>%
group_by(Litter_spp, Soil_type, Harvest) %>%
  summarise(BG_NAG_t0 = mean(BG_NAG)) %>%
  filter(Harvest == 1)

pox3 <- inner_join(pox, BG_NAG_t0[, -3], by = c("Litter_spp", "Soil_type"))

pox3$Enz_CN_change <- (pox3$BG_NAG_t0/pox3$BG_NAG)*100 # Initial over current OR
pox3$Enz_CN_change <- (pox3$BG_NAG/pox3$BG_NAG_t0)*100 # Current over initial


grids <- letters[1:3]
names(grids) <- unique(pox3$Soil_type)
pox3$grid <- unname(grids[pox3$Soil_type])

# pdf("Fig. S5.pdf", height = 4, width = 10)
ggplot(pox3, aes(x = days, y = Enz_CN_change)) +
  geom_hline(yintercept = 100, alpha = 0.5, linetype = "dashed") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 10, aes(color = Litter_spp)) +
  stat_summary(fun.data = "mean_se", geom = "line", aes(color = Litter_spp)) +
  stat_summary(fun.data = "mean_se", geom = "point", aes(color = Litter_spp)) +
  facet_grid(~ Soil_type) + 
  geom_text(x = 30, y = 860, aes(label = grid)) +
  scale_color_manual(values = cb_palette) +
  labs(x = "Time (days)", y = "BG:NAG (% of initial BG:NAG)", color = "Litter species", subtitle = "Soil type") +
  theme_classic() +
  theme(legend.position = c(0.94, 0.8),
        strip.text = element_text(size = 10, face = "italic"),
        legend.text = element_text(size = 10, face = "italic")) 
 #dev.off()


# BG effect size relative to time 0
BG_t0 <- 
  pox %>%
  group_by(Litter_spp, Soil_type, Harvest) %>%
  summarise(BG_t0 = mean(BG)) %>%
  filter(Harvest == 1)

pox3 <- inner_join(pox3, BG_t0[, -3], by = c("Litter_spp", "Soil_type"))

pox3$BG_change <- (pox3$BG/pox3$BG_t0)*100

# pdf("Fig. S5,pdf", height = 4, width = 10)
ggplot(pox3, aes(x = days, y = BG_change)) +
  geom_hline(yintercept = 100, alpha = 0.5, linetype = "dashed") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 10, aes(color = Litter_spp)) +
  stat_summary(fun.data = "mean_se", geom = "line", aes(color = Litter_spp)) +
  stat_summary(fun.data = "mean_se", geom = "point", aes(color = Litter_spp)) +
  facet_grid(~ Soil_type) + 
  geom_text(x = 30, y = 430, aes(label = grid)) +
  scale_color_manual(values = cb_palette) +
  labs(x = "Time (days)", y = "BG", color = "Litter species", subtitle = "Soil type") +
  theme_classic() +
  theme(legend.position = c(0.1, 0.8),
        strip.text = element_text(size = 10, face = "italic"),
        legend.text = element_text(size = 10, face = "italic")) 
# + ylim(c(0, 3.7))
# dev.off()






