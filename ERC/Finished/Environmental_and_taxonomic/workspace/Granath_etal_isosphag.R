###############################################################################################
# R-script to reproduce results presented in the paper:
# Granath et al. "Environmental and taxonomic controls of carbon and oxygen 
#   stable isotope composition in Sphagnum across broad climatic and geographic ranges"
# Journal: Biogeosciences, 15, 5189-5202, 2018
#          https://doi.org/10.5194/bg-15-5189-2018
#
# Contact: gustaf.granath@gmail.com
###############################################################################################

# tested under 
# R version 3.3.1 
# lme4_1.1-12 rworldmap_1.3-6 sp_1.2-3 piecewiseSEM_1.1.4 car_2.1-3           
# ggplot2_2.2.1  lattice_0.20-34

# load required packages
library(lme4)
library(car)
library(piecewiseSEM)
library(ggplot2)
library(rworldmap)
library(lattice)

# import data
dat <- read.csv("Granath_etal_gsp_isosphag.csv")

# Calculate ET/P variable
dat$ep <- dat$evap/dat$precip

# Correlations among variables
cor(dat[, c("prod", "HWT_end_season", "evap", "precip", "avt", "elevation")], use="complete.obs")

# Figure 1 - map and sampling summary ####

# Sampled sites and species
sisp <- aggregate(cbind(Patch_coord_lon, Patch_coord_lat) ~ Species + Site, data = dat, FUN=head, 1)
sisp <-reshape(sisp, timevar = "Species", v.names = "Species", idvar = "Site", direction = "wide")
sisp$samp <- factor(ifelse(!(is.na(sisp$Species.S.fuscum)) & is.na(sisp$Species.S.magellanicum), "Sf",
                           ifelse(is.na(sisp$Species.S.fuscum) & !(is.na(sisp$Species.S.magellanicum)),"Sm", "Sm+Sf")))
# make map
worldMap <- getMap()
world.points <- fortify(worldMap)
world.points$region <- world.points$id
world.df <- world.points[,c("long","lat","group", "region")]

worldmap <- ggplot() + 
  geom_polygon(data = world.df, aes(x = long, y = lat, group = group), fill = "darkgrey") +
  scale_y_continuous(breaks = (-2:2) * 30) +
  scale_x_continuous(breaks = (-4:4) * 45)+
  coord_map("ortho", orientation=c(90, 90, -90),xlim = c(-180, 175),  ylim = c(30, 80)) +
  theme(axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 15),
        legend.text = element_text(face = "italic", size = 14),
        legend.position = c(0.12,0.9),
        legend.key.size = unit(1, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_line(color="black", size=0.1),
        panel.ontop=TRUE,
        plot.background = element_rect(fill = "azure1", colour = "grey50")) +
  geom_point(data=sisp, aes(x=Patch_coord_lon, y=Patch_coord_lat, group = NULL, color = samp, shape = samp), size=2) + 
  scale_color_manual(name=c("Species"), labels= c("S. fuscum", "S.magellanicum", "S. fuscum & \nS. magellanicum"), 
                     values = c("black", "red", "blue")) +
  scale_shape_manual(name=c("Species"), labels= c("S. fuscum", "S.magellanicum", "S. fuscum & \nS. magellanicum"), 
                     values = c(1,2,3)) 

worldmap <- worldmap + # add lat-lon labels
  annotate("text", x = -45, y = 24, label =  as.character(expression(45*~degree*W)), size = 4, parse=T) +
  annotate("text", x = -23, y = 30, label = as.character(expression(30*~degree*N)), size = 4, parse=T) +
  annotate("text", x = -23, y = 60, label = as.character(expression(60*~degree*N)), size = 4, parse=T) +
  annotate("text", x = 45, y = 24, label = as.character(expression(45*~degree*E)), size = 4, parse=T) +
  annotate("text", x = 90, y = 37, label = as.character(expression(90*~degree*E)), size = 4, parse=T) +
  annotate("text", x = -90, y = 37, label = as.character(expression(90*~degree*W)), size = 4, parse=T)
worldmap

# save worldmap
ggsave(filename = "fig1_map.png", worldmap, width = 18, height = 18, units = "cm", dpi = 300)
#ggsave(filename = "fig1_map.eps", worldmap, width = 18, height = 18, units = "cm")

# Sampling summary
sum(data.frame(table(dat$Site, dat$Species))$Freq == 1) # count sites-species combinations with 1 sample

# sites with multiple samples at patch level
tt <- table(paste(dat$Site,dat$Patch_id))>1
sites <- do.call(rbind, strsplit(rownames(tt), ' (?=[^ ]+$)', perl=TRUE))[,1]
sum(table(data.frame(sites,tt))[,2]>0)
# Total number of patches
nlevels(factor(paste(dat$Site,dat$Patch_id)))

# Table 1 - Within-between site variation ####

# d13C
# S.fuscum
prop.var.fusc.mod <-lmer(iso_13C_all ~ 1 + 
                           (1|Site), dat[dat$Species == "S.fuscum",], na.action = "na.omit")
summary(prop.var.fusc.mod)
# between site
VarCorr(prop.var.fusc.mod)$Site[1] / (VarCorr(prop.var.fusc.mod)$Site[1] + attr(VarCorr(prop.var.fusc.mod), "sc")^2) 
# within-site
attr(VarCorr(prop.var.fusc.mod), "sc")^2 / (VarCorr(prop.var.fusc.mod)$Site[1] + attr(VarCorr(prop.var.fusc.mod), "sc")^2) 

# S.magellanicum
prop.var.mag.mod <-lmer(iso_13C_all ~ 1 + 
                          (1|Site), dat[dat$Species == "S.magellanicum",], na.action = "na.omit")
summary(prop.var.mag.mod)
# between site
VarCorr(prop.var.mag.mod)$Site[1] / (VarCorr(prop.var.mag.mod)$Site[1] + attr(VarCorr(prop.var.mag.mod), "sc")^2) 
# within-site
attr(VarCorr(prop.var.mag.mod), "sc")^2 / (VarCorr(prop.var.mag.mod)$Site[1] + attr(VarCorr(prop.var.mag.mod), "sc")^2) 

# d18O
# S.fuscum
prop.var.fusc.mod <-lmer(iso_18O_all ~ 1 + 
                           (1|Site), dat[dat$Species == "S.fuscum",], na.action = "na.omit")
summary(prop.var.fusc.mod)
# between site
VarCorr(prop.var.fusc.mod)$Site[1] / (VarCorr(prop.var.fusc.mod)$Site[1] + attr(VarCorr(prop.var.fusc.mod), "sc")^2) 
# within-site
attr(VarCorr(prop.var.fusc.mod), "sc")^2 / (VarCorr(prop.var.fusc.mod)$Site[1] + attr(VarCorr(prop.var.fusc.mod), "sc")^2) 

# S.magellanicum
prop.var.mag.mod <-lmer(iso_18O_all ~ 1 + 
                          (1|Site), dat[dat$Species == "S.magellanicum",], na.action = "na.omit")
summary(prop.var.mag.mod)
# between site
VarCorr(prop.var.mag.mod)$Site[1] / (VarCorr(prop.var.mag.mod)$Site[1] + attr(VarCorr(prop.var.mag.mod), "sc")^2) 
# within-site
attr(VarCorr(prop.var.mag.mod), "sc")^2 / (VarCorr(prop.var.mag.mod)$Site[1] + attr(VarCorr(prop.var.mag.mod), "sc")^2) 

# d13C analyses ####

# HWT is important for d13C but a few sites do not have HWT measurements and are removed
dat.sub <- dat[!(is.na(dat$HWT_end_season)),]

# Following part will build Table 2
H2.1 <- lmer(iso_13C_all ~ HWT_end_season * Species + 
               (HWT_end_season|Site), dat.sub, na.action = "na.omit", REML = TRUE)
H2.1.red <- update(H2.1, .~. -HWT_end_season:Species)
C1.res <- Anova(H2.1, type =2, test.statistic = "F")
C1.res 
# There is an interaction with HWT so HWT is centralised to make it easier to interpret the
# species effect (ie effect at mean HWT)

dat.sub$HWT_end_season.st <- scale(dat.sub$HWT_end_season, scale=FALSE) # centralise 
H2.1 <- lmer(iso_13C_all ~ HWT_end_season.st * Species + 
               (HWT_end_season.st|Site), dat.sub, na.action = "na.omit", REML = TRUE)
H2.1.red <- update(H2.1, .~. -HWT_end_season:Species)
C1.res <- Anova(H2.1, type =2, test.statistic = "F")
C1.res
C1.res$Df.res <- round(C1.res$Df.res, 2)
summary(H2.1)
# interaction effect so refit to get slopes for s.mag and s.fusc
betas <- do.call(paste, c(as.data.frame(round(summary(H2.1)$coefficients[,1:2], digits=5)), sep="+-"))
betas[4] <- betas[2] # S. fuscum slope
betas[2] <- NA # no main effect of HWT
H2.1.2 <- update(H2.1, .~ HWT_end_season.st*relevel(Species,ref="S.magellanicum")+(HWT_end_season.st|Site) )
betas[4] <- paste("S.fus: ", betas[4], " / ",
                  "S.mag: ", do.call(paste, c(as.data.frame(round(summary(H2.1.2)$coefficients[,1:2], digits=5)), sep="+-"))[2],
                  sep = "")

plot(H2.1)
summary(H2.1) 

H2.1.sp <- update(H2.1, .~ Species + (1|Site))
H2.1.hwt <- update(H2.1, .~ HWT_end_season + (HWT_end_season|Site))
H2.1.null <- update(H2.1, .~ 1 + (1|Site))

# within-site variation explained by HWT
(attr(VarCorr(H2.1.sp), "sc")^2-attr(VarCorr(H2.1), "sc")^2) / attr(VarCorr(H2.1.sp), "sc")^2 

# Variance explained by HWT * species model
Cexp.pred13C.int <- as.numeric((VarCorr(H2.1.null)$Site[1]-VarCorr(H2.1)$Site[1]) / VarCorr(H2.1.null)$Site[1])
Cexp.pred13C.int.tot <- rsquared(H2.1)[5]$Marginal

varExp <- round(cbind(c(Cexp.pred13C.int), 
                      c(Cexp.pred13C.int.tot)), digits = 4)

nm <- 19 # number of table rows
Ctable <- data.frame(model=rep(NA,nm), variable = rep(NA,nm), effect = rep(NA,nm), Fvalue = rep(NA,nm), 
                     df = rep(NA,nm), Pvalue = rep(NA,nm), N = rep(NA,nm),
                     varExpsite = rep(NA,nm), varExptotal = rep(NA,nm))

Ctable[1,"model"] <- paste(deparse(formula(H2.1))[1], deparse(formula(H2.1))[2],sep="")
Ctable$variable[1:(0+nrow(C1.res))] <- rownames(C1.res)
Ctable$effect[1:(0+nrow(C1.res))] <- betas[-1]
Ctable$Fvalue[1:(0+nrow(C1.res))] <- round(C1.res$F, digits = 2)
Ctable$df[1:(0+nrow(C1.res))] <- do.call(paste, c(C1.res[c("Df", "Df.res")], sep=", "))
Ctable$Pvalue[1:(0+nrow(C1.res))] <- round(C1.res$'Pr(>F)', digits = 4)
Ctable$N[1] <- nrow(H2.1@frame)
Ctable[3, c("varExpsite", "varExptotal")] <- varExp

# Check Table
Ctable

# S. magellanicum most important and a
# submodel is fitted for S.magellanicum
dat.sub.mag <- dat.sub[dat.sub$Species == "S.magellanicum",]
dat.sub.mag$iso_13C_all_Smag <- dat.sub.mag$iso_13C_all
H2.1.mag <- lmer(iso_13C_all_Smag ~ HWT_end_season.st + 
                   (1|Site), dat.sub.mag, 
                 na.action = "na.omit", REML = TRUE)
H2.1.mag.null <- lmer(iso_13C_all_Smag ~ 1 + 
                        (1|Site), dat.sub.mag, 
                      na.action = "na.omit", REML = TRUE)

C2.res <- Anova(H2.1.mag, type = 2, test.statistic = "F")
summary(H2.1.mag)
C2.res$Df.res <- round(C2.res$Df.res, 2)
betas <- do.call(paste, c(as.data.frame(round(summary(H2.1.mag)$coefficients[,1:2], digits=5)), sep="+-"))

# variance explained by HWT for S.mag
Cexp.pred13C.hwt.mag <- as.numeric((VarCorr(H2.1.mag.null)$Site[1]-VarCorr(H2.1.mag)$Site[1]) / VarCorr(H2.1.mag.null)$Site[1])
Cexp.pred13C.hwt.mag.tot <- rsquared(H2.1.mag)[5]$Marginal
varExp <- round(cbind(c(Cexp.pred13C.hwt.mag), 
                      c(Cexp.pred13C.hwt.mag.tot)), digits = 4)

# within-site variation explained by HWT for S.mag
(attr(VarCorr(H2.1.mag.null), "sc")^2-attr(VarCorr(H2.1.mag), "sc")^2) / attr(VarCorr(H2.1.mag.null), "sc")^2 

Ctable[4,"model"] <- deparse(formula(H2.1.mag))[1]
Ctable$variable[4:(3+nrow(C2.res))] <- rownames(C2.res)
Ctable$effect[4:(3+nrow(C2.res))] <- betas[-1]
Ctable$Fvalue[4:(3+nrow(C2.res))] <- round(C2.res$F, digits = 2)
Ctable$df[4:(3+nrow(C2.res))] <- do.call(paste, c(C2.res[c("Df", "Df.res")], sep=", "))
Ctable$Pvalue[4:(3+nrow(C2.res))] <- round(C2.res$'Pr(>F)', digits = 4)
Ctable$N[4] <- nrow(H2.1.mag@frame)
Ctable[4, c("varExpsite", "varExptotal")] <- varExp

# test NPP separately
dat.sub <-dat[!(is.na(dat$prod)),] # a few samples do not have NPP data
H2.2 <- lmer(iso_13C_all ~ prod * Species + 
               (1|Site), dat.sub, na.action = "na.omit", REML =TRUE)
H2.2.red <- update(H2.2, .~. -prod:Species) 
H2.2.red.sp <- update(H2.2.red, .~. -prod) 
H2.2.red.prod <- update(H2.2.red, .~. -Species) 
H2.2.null <- update(H2.2, .~.-prod*Species)
C3.res <- Anova(H2.2, type =2, test.statistic = "F")
C3.res
C3.res$Df.res <- round(C3.res$Df.res, 2)
summary(H2.2)
betas <- do.call(paste, c(as.data.frame(round(summary(H2.2.red)$coefficients[,1:2], digits=5)), sep="+-"))

# No interaction so calculate R2 for each main effect
# NPP
Cexp.pred13C.prod <- as.numeric((VarCorr(H2.2.red.sp)$Site[1]-VarCorr(H2.2.red)$Site[1]) / VarCorr(H2.2.null)$Site[1])
Cexp.pred13C.prod.tot <- rsquared(H2.2.red)[5]$Marginal-rsquared(H2.2.red.sp)[5]$Marginal

# Species
Cexp.pred13C.sp <- as.numeric((VarCorr(H2.2.red.prod)$Site[1]-VarCorr(H2.2.red)$Site[1]) / VarCorr(H2.2.null)$Site[1])
Cexp.pred13C.sp.tot <- rsquared(H2.2.red)[5]$Marginal-rsquared(H2.2.red.prod)[5]$Marginal

varExp <- round(cbind(c(Cexp.pred13C.prod, Cexp.pred13C.sp), 
                      c(Cexp.pred13C.prod.tot, Cexp.pred13C.sp.tot)), digits = 4)

# within-site variation explained by NPP
(attr(VarCorr(H2.2.red.sp), "sc")^2-attr(VarCorr(H2.2.red), "sc")^2) / attr(VarCorr(H2.2.red.sp), "sc")^2 

Ctable[5,"model"] <- deparse(formula(H2.2))[1]
Ctable$variable[5:(4+nrow(C3.res))] <- rownames(C3.res)
Ctable$effect[5:6] <- betas[-1]
Ctable$Fvalue[5:(4+nrow(C3.res))] <- round(C3.res$F, digits = 2)
Ctable$df[5:(4+nrow(C3.res))] <- do.call(paste, c(C3.res[c("Df", "Df.res")], sep=", "))
Ctable$Pvalue[5:(4+nrow(C3.res))] <- round(C3.res$'Pr(>F)', digits = 4)
Ctable$N[5] <- nrow(H2.2@frame)
Ctable[5:6, c("varExpsite", "varExptotal")] <- varExp

# Test HWT and NPP for S. mag
H2.2.mag <- lmer(iso_13C_all_Smag ~ HWT_end_season + prod + 
                   (1|Site), dat.sub.mag, 
                 na.action = "na.omit", REML = TRUE)
H2.2.mag.null <- lmer(iso_13C_all_Smag ~ 1 + 
                        (1|Site), dat.sub.mag, 
                      na.action = "na.omit", REML = TRUE)
# between site expl variation
as.numeric((VarCorr(H2.1.mag.null)$Site[1]-VarCorr(H2.2.mag)$Site[1]) / VarCorr(H2.1.mag.null)$Site[1])
# total
rsquared(H2.2.mag)[5]$Marginal
rsquared(H2.2.mag)[5]$Marginal - rsquared(H2.1.mag)[5]$Marginal #added by NPP

#  Test HWT and NPP for S. fuscum
dat.sub.fus <- dat.sub[dat.sub$Species == "S.fuscum",]
dat.sub.fus$iso_13C_all_Sfus <- dat.sub.fus$iso_13C_all
H2.1.fus <- lmer(iso_13C_all_Sfus ~ HWT_end_season  + 
                   (1|Site), dat.sub.fus, 
                 na.action = "na.omit", REML = TRUE)
H2.1.fus.null <- lmer(iso_13C_all_Sfus ~ 1 + 
                        (1|Site), dat.sub.fus, 
                      na.action = "na.omit", REML = TRUE)
H2.2.fus <- lmer(iso_13C_all_Sfus ~ HWT_end_season + prod + 
                   (1|Site), dat.sub.fus, 
                 na.action = "na.omit", REML = TRUE)

# between site expl variation
as.numeric((VarCorr(H2.1.fus.null)$Site[1]-VarCorr(H2.1.fus)$Site[1]) / VarCorr(H2.1.fus.null)$Site[1])
# total
rsquared(H2.2.fus)[5]$Marginal
rsquared(H2.2.fus)[5]$Marginal - rsquared(H2.1.fus)[5]$Marginal #added by NPP

# Test additional environmental variables 
# elevation + temperature + ET/P + precipitation

# centralise HWT
dat$HWT_end_season.st <- scale(dat$HWT_end_season, scale=FALSE) # centralise HWT
# make data set with S.mag as reference level for later coef estiamtes
dat2 <- dat 
dat2$Species <- relevel(dat2$Species, ref = "S.magellanicum")

# Additional variables were tested by adding them seperately to a "start model" that
# was considered the "best" model: HWT*Species + NPP.
# the HWT*NPP term was not included as it was not important
H2.3.full <- lmer(iso_13C_all ~ HWT_end_season.st*Species + HWT_end_season.st*prod +
       (1|Site), dat, na.action = "na.omit", REML = TRUE)
Anova(H2.3.full, type =2, test.statistic = "F") # not important
# continue with "best" model
H2.3.start <- lmer(iso_13C_all ~ HWT_end_season.st*Species + prod +
                     (1|Site), dat, na.action = "na.omit", REML = TRUE)
H2.3.start.no_prod <- update(H2.3.start, .~ . -prod)
H2.3.start.red <- update(H2.3.start, .~ . -HWT_end_season.st:Species)
H2.3.start.red.no_hwt <- update(H2.3.start.red, .~ . -HWT_end_season.st)
H2.3.start.red.no_sp <- update(H2.3.start.red, .~ . -Species)
H2.3.start.red.no_prod <- update(H2.3.start.red, .~ . -prod)
H2.3.start.only_prod <- update(H2.3.start, .~ . -HWT_end_season.st*Species)
H2.3.null <- update(H2.3.start, .~ . -HWT_end_season.st*Species -prod)

C4.res <- Anova(H2.3.start, type =2, test.statistic = "F")
C4.res$Df.res <- round(C4.res$Df.res, 2)
summary(H2.3.start)
betas <- do.call(paste, c(as.data.frame(round(summary(H2.3.start)$coefficients[,1:2], digits=5)), sep="+-"))

# interaction effect present so refit to get slope for s.mag and s.fusc
betas <- do.call(paste, c(as.data.frame(round(summary(H2.3.start)$coefficients[,1:2], digits=5)), sep="+-"))
betas[5] <- betas[2] # S. fuscum slope
betas[2] <- NA # no main effect of HWT
H2.3.start.mag <- update(H2.3.start, .~ HWT_end_season.st*relevel(Species,ref="S.magellanicum") + prod+(1|Site) )
betas[5] <- paste("S.fus: ", betas[5], " / ",
                  "S.mag: ", do.call(paste, c(as.data.frame(round(summary(H2.3.start.mag)$coefficients[,1:2], digits=5)), sep="+-"))[2],
                  sep = "")

# HWT*species interaction
Cexp.pred13C.int <- as.numeric((VarCorr(H2.3.start.only_prod)$Site[1]-VarCorr(H2.3.start)$Site[1]) / VarCorr(H2.3.null)$Site[1])
Cexp.pred13C.int.tot <- (rsquared(H2.3.start)[5] - rsquared(H2.3.start.only_prod)[5])$Marginal

varExp <- round(cbind(c(Cexp.pred13C.prod, Cexp.pred13C.int), 
                      c(Cexp.pred13C.prod.tot, Cexp.pred13C.int.tot)), digits = 4)

Ctable[8,"model"] <- deparse(formula(H2.3.start))[1]
Ctable$variable[8:(7+nrow(C4.res))] <- rownames(C4.res)
Ctable$effect[8:(7+nrow(C4.res))] <- betas[-1]
Ctable$Fvalue[8:(7+nrow(C4.res))] <- round(C4.res$F, digits = 2)
Ctable$df[8:(7+nrow(C4.res))] <- do.call(paste, c(C4.res[c("Df", "Df.res")], sep=", "))
Ctable$Pvalue[8:(7+nrow(C4.res))] <- round(C4.res$'Pr(>F)', digits = 4)
Ctable$N[8] <- nrow(H2.3.start@frame)
Ctable[10:11, c("varExpsite", "varExptotal")] <- varExp

# add elevation to start model
H2.3.start.ele <- update(H2.3.start, .~ . + elevation*Species)
C5.res <- Anova(H2.3.start.ele, type =2, test.statistic = "F")
C5.res$Df.res <- round(C5.res$Df.res, 2)
summary(H2.3.start.ele)

# no interaction so get main effect coef
H2.3.start.ele_main <- update(H2.3.start.ele, .~ . - elevation:Species)
H2.3.start.ele_main_no_ele <- update(H2.3.start.ele_main, .~ . - elevation)

betas <- do.call(paste, c(as.data.frame(round(summary(H2.3.start.ele_main)$coefficients[,1:2], digits=5)), sep="+-"))

# explained variance by elevation
Cexp.pred13C.ele <- as.numeric((VarCorr(H2.3.start.ele_main_no_ele)$Site[1]-VarCorr(H2.3.start.ele_main)$Site[1]) / 
                                 VarCorr(H2.3.null)$Site[1])
Cexp.pred13C.ele.tot <- (rsquared(H2.3.start.ele_main)[5] - rsquared(H2.3.start.ele_main_no_ele)[5])$Marginal

varExp <- round(cbind(c(Cexp.pred13C.ele), 
                      c(Cexp.pred13C.ele.tot)), digits = 4)

Ctable[12,"model"] <- "start+elevation*species"
Ctable$variable[12] <- "elevation"
Ctable$effect[12] <- betas[5]
Ctable$Fvalue[12:13] <- round(C5.res$F, digits = 2)[c(4,6)]
Ctable$df[12:13] <- do.call(paste, c(C5.res[c("Df", "Df.res")], sep=", "))[c(4,6)]
Ctable$Pvalue[12:13] <- round(C5.res$'Pr(>F)', digits = 4)[c(4,6)]
Ctable$N[12] <- nrow(H2.3.start@frame)
Ctable[12, c("varExpsite", "varExptotal")] <- varExp

# add ET/P to start model
H2.3.start.ep <- update(H2.3.start, .~ . + ep*Species)
C6.res <- Anova(H2.3.start.ep, type =2, test.statistic = "F")
C6.res$Df.res <- round(C6.res$Df.res, 2)
summary(H2.3.start.ep)

# interaction effect so refit to get slope for s.mag and s.fusc
betas <- do.call(paste, c(as.data.frame(round(summary(H2.3.start.ep)$coefficients[,1:2], digits=5)), sep="+-"))
betas[7] <- betas[5] # S. fuscum slope
betas[5] <- NA # no main effect of HWT
H2.3.start.ep2 <- update(H2.3.start.ep, .~. , data=dat2)
betas[7] <- paste("S.fus: ", betas[7], " / ",
                  "S.mag: ", do.call(paste, c(as.data.frame(round(summary(H2.3.start.ep2)$coefficients[,1:2], digits=5)), sep="+-"))[5],
                  sep = "")

# ET/P*speceis interaction
Cexp.pred13C.ep.int <- as.numeric((VarCorr(H2.3.start)$Site[1]-VarCorr(H2.3.start.ep)$Site[1]) / VarCorr(H2.3.null)$Site[1])
Cexp.pred13C.int.ep.tot <- rsquared(H2.3.start.ep)[5]$Marginal - rsquared(H2.3.start)[5]$Marginal

varExp <- round(cbind(c(Cexp.pred13C.ep.int), 
                      c(Cexp.pred13C.int.ep.tot)), digits = 4)

Ctable[14,"model"] <- "start+(evap/precip)*species"
Ctable$variable[14:15] <- c("evap/precip", "e/p * Species")
Ctable$effect[15] <- betas[7]
Ctable$Fvalue[14:15] <- round(C6.res$F, digits = 2)[c(4,6)]
Ctable$df[14:15] <- do.call(paste, c(C6.res[c("Df", "Df.res")], sep=", "))[c(4,6)]
Ctable$Pvalue[14:15] <- round(C6.res$'Pr(>F)', digits = 4)[c(4,6)]
Ctable$N[14] <- nrow(H2.3.start@frame)
Ctable[15, c("varExpsite", "varExptotal")] <- varExp

# add precipitation to start model
H2.3.start.pre <- update(H2.3.start, .~ . + precip*Species)
C7.res <- Anova(H2.3.start.pre, type =2, test.statistic = "F")
C7.res$Df.res <- round(C7.res$Df.res, 2)
summary(H2.3.start.pre)

# no interaction so get main effect coef
H2.3.start.pre_main <- update(H2.3.start.pre, .~ . - precip:Species)
H2.3.start.pre_main_no_pre <- update(H2.3.start.pre_main, .~ . - precip)

betas <- do.call(paste, c(as.data.frame(round(summary(H2.3.start.pre_main)$coefficients[,1:2], digits=5)), sep="+-"))

# Explained variance by precipitation
Cexp.pred13C.pre <- as.numeric((VarCorr(H2.3.start.pre_main_no_pre)$Site[1]-VarCorr(H2.3.start.pre_main)$Site[1]) / 
                                 VarCorr(H2.3.null)$Site[1])
Cexp.pred13C.pre.tot <- (rsquared(H2.3.start.pre_main)[5] - rsquared(H2.3.start.pre_main_no_pre)[5])$Marginal

varExp <- round(cbind(c(Cexp.pred13C.pre), 
                      c(Cexp.pred13C.pre.tot)), digits = 4)

Ctable[16,"model"] <- "start+precipitation*species"
Ctable$variable[16:17] <- c("precip", "precip*Species")
Ctable$effect[16] <- betas[5]
Ctable$Fvalue[16:17] <- round(C7.res$F, digits = 2)[c(4,6)]
Ctable$df[16:17] <- do.call(paste, c(C7.res[c("Df", "Df.res")], sep=", "))[c(4,6)]
Ctable$Pvalue[16:17] <- round(C7.res$'Pr(>F)', digits = 4)[c(4,6)]
Ctable$N[16] <- nrow(H2.3.start@frame)
Ctable[16, c("varExpsite", "varExptotal")] <- varExp

# add mean temperature to start model
H2.3.start.tmp <- update(H2.3.start, .~ . + avt*Species)
C8.res <- Anova(H2.3.start.tmp, type =2, test.statistic = "F")
C8.res$Df.res <- round(C8.res$Df.res, 2)
summary(H2.3.start.tmp)

# interaction effect so refit to get slope for s.mag and s.fusc
betas <- do.call(paste, c(as.data.frame(round(summary(H2.3.start.tmp)$coefficients[,1:2], digits=5)), sep="+-"))
betas[7] <- betas[5] # S. fuscum slope
betas[5] <- NA # no main effect of HWT
H2.3.start.tmp2 <- update(H2.3.start.tmp, .~. , data=dat2)
betas[7] <- paste("S.fus: ", betas[7], " / ",
                  "S.mag: ", do.call(paste, c(as.data.frame(round(summary(H2.3.start.tmp2)$coefficients[,1:2], digits=5)), sep="+-"))[5],
                  sep = "")

# Temp*species interaction
Cexp.pred13C.tmp.int <- as.numeric((VarCorr(H2.3.start)$Site[1]-VarCorr(H2.3.start.tmp)$Site[1]) / VarCorr(H2.3.null)$Site[1])
Cexp.pred13C.int.tmp.tot <- rsquared(H2.3.start.tmp)[5]$Marginal - rsquared(H2.3.start)[5]$Marginal

varExp <- round(cbind(c(Cexp.pred13C.tmp.int), 
                      c(Cexp.pred13C.int.tmp.tot)), digits = 4)

Ctable[18,"model"] <- "start+temp*species"
Ctable$variable[18:19] <- c("temp", "temp*Species")
Ctable$effect[19] <- betas[7]
Ctable$Fvalue[18:19] <- round(C8.res$F, digits = 2)[c(4,6)]
Ctable$df[18:19] <- do.call(paste, c(C8.res[c("Df", "Df.res")], sep=", "))[c(4,6)]
Ctable$Pvalue[18:19] <- round(C8.res$'Pr(>F)', digits = 4)[c(4,6)]
Ctable$N[18] <- nrow(H2.3.start@frame)
Ctable[19, c("varExpsite", "varExptotal")] <- varExp

write.csv(Ctable, file = "table2_C.csv")

# Figure 2 ####
# subset with 13C data and HWT
dat.sub.f2 <- dat[!(is.na(dat$iso_13C_all)) & !(is.na(dat$HWT_end_season)), ]
# get equations for Fig 2
fig2.mod <- lmer(iso_13C_all ~ HWT_end_season * Species + 
                   (HWT_end_season|Site), dat.sub.f2, na.action = "na.omit", REML = TRUE)
# S.fuscum
paste(fixef(fig2.mod)[1]," + ",fixef(fig2.mod)[2])
# S.magellanicum
paste((fixef(fig2.mod)[1] + fixef(fig2.mod)[3]), " + ", (fixef(fig2.mod)[2]+fixef(fig2.mod)[4]) )

# check site effect
rf <- ranef( fig2.mod, condVar = TRUE)
dotplot(rf)
dat.sub.f2$pred <- predict(fig2.mod)

# Plot
#newdat <- rbind(aggregate(HWT_end_season ~ Species, data=dat.sub.f2, FUN = min),
#                aggregate(HWT_end_season ~ Species, data=dat.sub.f2, FUN = max))

#
mima <- aggregate(HWT_end_season ~ Species, data=dat.sub.f2, function (x) c(min(x),max(x)) )$HWT_end_season
newdat <- rbind(data.frame(Species = "S.fuscum", HWT_end_season = seq(mima[1, 1],mima[1, 2], 0.25)),
      data.frame(Species = "S.magellanicum", HWT_end_season = seq(mima[2, 1],mima[2, 2], 0.25)))
newdat$predMean <- predict(fig2.mod, newdata=newdat, re.form=~0)

mm <- model.matrix(formula('~HWT_end_season * Species'),newdat)
newdat$distance <- mm %*% fixef(fig2.mod)
newdat$pvar1 <- diag(mm %*% tcrossprod(vcov(fig2.mod),mm))
newdat <- data.frame(
  newdat
  , plo = newdat$distance-2*sqrt(newdat$pvar1)
  , phi = newdat$distance+2*sqrt(newdat$pvar1) )
#


#newdat$predMean <- predict(fig2.mod, newdata=newdat, re.form=~0)
newdat$Site <- dat.sub.f2$Site[1] # hack to avoid ggplot error

fig2 <- ggplot(dat.sub.f2,aes(HWT_end_season, iso_13C_all, group=interaction(Site, Species), 
                           color=Species, shape = Species )) + 
  geom_ribbon(data=newdat, aes(y=plo,x=HWT_end_season, ymin = plo, ymax = phi, fill = Species), 
              alpha = .15, show.legend = FALSE, linetype=0) +
  geom_point(alpha = 1) +
  #geom_line(aes(y=pred, color=Species), size=0.8, alpha =0.3) + #if plot site-specific slopes
  geom_line(data=newdat, aes(y=predMean, x=HWT_end_season, color=Species), alpha =1, size = 1.2) +
  labs(y=expression(paste("Tissue ", delta^{13}, "C (\u2030)")), 
       x= "HWT (cm)") +
  scale_fill_manual(values=c("black", "red"), name="fill") +
  scale_color_manual(breaks = c("S.fuscum", "S.magellanicum"), values = c("S.fuscum" = "black", "S.magellanicum" = "red"),
                     labels = c("S. fuscum", "S. magellanicum")) +
  scale_shape_manual(breaks = c("S.fuscum", "S.magellanicum"), values = c(1,2),
                     labels = c("S. fuscum", "S. magellanicum")) +
  ylim(c(-35, -23)) + 
  theme_bw()+
  theme(legend.justification = c(0, 0), 
        legend.position = c(0.62, 0.02),
        legend.key.size = unit(2, "line"),
        legend.text = element_text(size=12, face = "italic"),
        legend.title = element_text(size=12)) +
  guides(shape = guide_legend(override.aes = list(size = 1)))

fig2
# legend line and point size are not equal and here is a hack
# to fix this
library(grid)
grid.ls(grid.force()) # check slots to edit 
grid.gedit("key-[-0-9]-1-1.4-2-4-2", size = unit(2, "mm"))    
grid.gedit("key-[-0-9]-1-1.5-2-5-2", size = unit(2, "mm"))    
fig2 <- grid.grab() # save edits
ggsave("fig2_13C_ci.png", fig2,  dpi = 300, height = 5, width=6, units="in")
# save as eps
cairo_pdf(filename = "fig2_13C_ci.eps", height = 5, width=6) 
grid.draw(fig2)
dev.off()

# d18O analyses ####
# First we remove the two samples without O isotopes. 
dat.sub <- dat[!(is.na(dat$iso_18O_all)),]

# The following code will produce Table 3.
# First, Mean annual precip d18O 
H1.1 <- lmer(iso_18O_all ~ bowen_18O_annual*Species + 
               (1|Site), dat.sub, na.action = "na.omit", REML = TRUE)
H1.1.red1 <- lmer(iso_18O_all ~ bowen_18O_annual + Species + 
                    (1|Site), dat.sub, na.action = "na.omit", REML = TRUE)
H1.1.red2 <- lmer(iso_18O_all ~  Species + 
                    (1|Site), dat.sub, na.action = "na.omit", REML = TRUE)
H1.1.red3 <- lmer(iso_18O_all ~  bowen_18O_annual + 
                    (1|Site), dat.sub, na.action = "na.omit", REML = TRUE)
H1.1.null <- lmer(iso_18O_all ~  1 + 
                    (1|Site), dat.sub, na.action = "na.omit", REML = TRUE)
O1.res <- Anova(H1.1, type = "2", test.statistic="F")
O1.res
O1.res$Df.res <- round(O1.res$Df.res, 2)

plot(H1.1)
summary(H1.1)
summary(H1.1.red1) # main effect of species
betas <- do.call(paste, c(as.data.frame(round(summary(H1.1.red1)$coefficients[,1:2], digits=5)), sep="+-"))

# variance explained
# modelled d18O explaining between site variation
Oexp.pred18O <- as.numeric((VarCorr(H1.1.red2)$Site-VarCorr(H1.1.red1)$Site) / VarCorr(H1.1.null)$Site)
Oexp.pred18O.tot <- rsquared(H1.1.red1)[5] - rsquared(H1.1.red2)[5]

# species explaining variation
Oexp.species <- rsquared(H1.1.red1)[5] - rsquared(H1.1.red3)[5]

nm <- 17 # number of table rows 
Otable <- data.frame(model=rep(NA,nm), variable = rep(NA,nm), effect = rep(NA,nm), Fvalue = rep(NA,nm), 
                     df = rep(NA,nm), Pvalue = rep(NA,nm), N = rep(NA,nm),
                     varExpsite = rep(NA,nm), varExptotal = rep(NA,nm))
Otable[1,"model"] <- deparse(formula(H1.1))
Otable$variable[1:nrow(O1.res)] <- rownames(O1.res)
Otable$effect[1:2] <- betas[-1] # for model w/o interaction (main effects)
Otable$Fvalue[1:nrow(O1.res)] <- round(O1.res$F, digits = 2)
Otable$df[1:nrow(O1.res)] <- do.call(paste, c(O1.res[c("Df", "Df.res")], sep=", "))
Otable$Pvalue[1:nrow(O1.res)] <- round(O1.res$'Pr(>F)', digits = 4)
Otable$N[1] <- nrow(H1.1@frame)
Otable$varExpsite[1] <- round(Oexp.pred18O, digits = 3)
Otable$varExptotal[1] <- round(Oexp.pred18O.tot, digits = 3)[[1]]
Otable$varExptotal[2] <- round(Oexp.species, digits = 3)[[1]]


# Second, growing season precip d18O 
H1.1 <- lmer(iso_18O_all ~ bowen_18O_AprSep*Species + 
               (1|Site), dat.sub, na.action = "na.omit", REML = TRUE)
H1.1.red1 <- lmer(iso_18O_all ~ bowen_18O_AprSep + Species + 
                    (1|Site), dat.sub, na.action = "na.omit", REML = TRUE)
H1.1.red2 <- lmer(iso_18O_all ~  Species + 
                    (1|Site), dat.sub, na.action = "na.omit", REML = TRUE)
H1.1.red3 <- lmer(iso_18O_all ~  bowen_18O_AprSep + 
                    (1|Site), dat.sub, na.action = "na.omit", REML = TRUE)
O1.res <- Anova(H1.1, type = "2", test.statistic="F")
O1.res
O1.res$Df.res <- round(O1.res$Df.res, 2)

plot(H1.1)
summary(H1.1)
summary(H1.1.red1) # main effect of species
betas <- do.call(paste, c(as.data.frame(round(summary(H1.1.red1)$coefficients[,1:2], digits=5)), sep="+-"))

# variance explained
# modelled d18O explaining between site variation
Oexp.pred18O <- as.numeric((VarCorr(H1.1.red2)$Site-VarCorr(H1.1.red1)$Site) / VarCorr(H1.1.null)$Site)
Oexp.pred18O.tot <- rsquared(H1.1.red1)[5] - rsquared(H1.1.red2)[5]

# species explaining variation
Oexp.species <- rsquared(H1.1.red1)[5] - rsquared(H1.1.red3)[5]

Otable[4,"model"] <- deparse(formula(H1.1))
Otable$variable[4:(3+nrow(O1.res))] <- rownames(O1.res)
Otable$effect[4:5] <- betas[-1] # for model w/o interaction (main effects)
Otable$Fvalue[4:(3+nrow(O1.res))] <- round(O1.res$F, digits = 2)
Otable$df[4:(3+nrow(O1.res))] <- do.call(paste, c(O1.res[c("Df", "Df.res")], sep=", "))
Otable$Pvalue[4:(3+nrow(O1.res))] <- round(O1.res$'Pr(>F)', digits = 4)
Otable$N[4] <- nrow(H1.1@frame)
Otable$varExpsite[4] <- round(Oexp.pred18O, digits = 3)
Otable$varExptotal[4] <- round(Oexp.pred18O.tot, digits = 3)[[1]]
Otable$varExptotal[5] <- round(Oexp.species, digits = 3)[[1]]


# Third, winter season precip d18O 
H1.1 <- lmer(iso_18O_all ~ bowen_18O_JanApr*Species + 
               (1|Site), dat.sub, na.action = "na.omit", REML = TRUE)
H1.1.red1 <- lmer(iso_18O_all ~ bowen_18O_JanApr + Species + 
                    (1|Site), dat.sub, na.action = "na.omit", REML = TRUE)
H1.1.red2 <- lmer(iso_18O_all ~  Species + 
                    (1|Site), dat.sub, na.action = "na.omit", REML = TRUE)
H1.1.red3 <- lmer(iso_18O_all ~  bowen_18O_JanApr + 
                    (1|Site), dat.sub, na.action = "na.omit", REML = TRUE)
O1.res <- Anova(H1.1, type = "2", test.statistic="F")
O1.res
O1.res$Df.res <- round(O1.res$Df.res, 2)

plot(H1.1)
summary(H1.1)
summary(H1.1.red1) # main effect of species
betas <- do.call(paste, c(as.data.frame(round(summary(H1.1.red1)$coefficients[,1:2], digits=5)), sep="+-"))

# variance explained
# modelled d18O explaining between site variation
Oexp.pred18O <- as.numeric((VarCorr(H1.1.red2)$Site-VarCorr(H1.1.red1)$Site) / VarCorr(H1.1.null)$Site)
Oexp.pred18O.tot <- rsquared(H1.1.red1)[5] - rsquared(H1.1.red2)[5]

# species explaining variation
Oexp.species <- rsquared(H1.1.red1)[5] - rsquared(H1.1.red3)[5]

Otable[7,"model"] <- deparse(formula(H1.1))
Otable$variable[7:(6+nrow(O1.res))] <- rownames(O1.res)
Otable$effect[7:8] <- betas[-1] # for model w/o interaction (main effects)
Otable$Fvalue[7:(6+nrow(O1.res))] <- round(O1.res$F, digits = 2)
Otable$df[7:(6+nrow(O1.res))] <- do.call(paste, c(O1.res[c("Df", "Df.res")], sep=", "))
Otable$Pvalue[7:(6+nrow(O1.res))] <- round(O1.res$'Pr(>F)', digits = 4)
Otable$N[7] <- nrow(H1.1@frame)
Otable$varExpsite[7] <- round(Oexp.pred18O, digits = 3)
Otable$varExptotal[7] <- round(Oexp.pred18O.tot, digits = 3)[[1]]
Otable$varExptotal[8] <- round(Oexp.species, digits = 3)[[1]]

# Check first part of Table 3
Otable

# Test with precipitation-weighted modelled d18O
winter.w.full <- lmer(iso_18O_all ~ bowen_18O_JanApr.w + Species + 
               (1|Site), dat.sub, na.action = "na.omit", REML = TRUE)
winter.w.red <- lmer(iso_18O_all ~ Species + 
                        (1|Site), dat.sub, na.action = "na.omit", REML = TRUE)
grow.w.full <- lmer(iso_18O_all ~ bowen_18O_AprSep + Species + 
                        (1|Site), dat.sub, na.action = "na.omit", REML = TRUE)
grow.w.red <- lmer(iso_18O_all ~ Species + 
                       (1|Site), dat.sub, na.action = "na.omit", REML = TRUE)
annual.w.full <- lmer(iso_18O_all ~ bowen_18O_annual.w + Species + 
                         (1|Site), dat.sub, na.action = "na.omit", REML = TRUE)
annual.w.red <- lmer(iso_18O_all ~ Species + 
                        (1|Site), dat.sub, na.action = "na.omit", REML = TRUE)

# variance explained by precip-weighted
# modelled d18O explaining between site variation
Oexp.pred18O.winter <- as.numeric((VarCorr(winter.w.red )$Site-VarCorr(winter.w.full)$Site) / VarCorr(H1.1.null)$Site)
Oexp.pred18O.grow <- as.numeric((VarCorr(grow.w.red )$Site-VarCorr(grow.w.full)$Site) / VarCorr(H1.1.null)$Site)
Oexp.pred18O.annual <- as.numeric((VarCorr(annual.w.red )$Site-VarCorr(annual.w.full)$Site) / VarCorr(H1.1.null)$Site)

Oexp.pred18O.tot <- rsquared(H1.1.red1)[5] - rsquared(H1.1.red2)[5]



# Annual precip d18O is the best predictor and will be used
# in subsequent models
# First, test effects of precip d18O + HWT

# a few sites without HWT measurements so they are removed
dat.sub <- dat[!(is.na(dat$HWT_end_season)),]

H1.3.full <- lmer(iso_18O_all ~ bowen_18O_annual + HWT_end_season + 
                    (1|Site), dat.sub, REML = TRUE)
H1.3.red <- lmer(iso_18O_all ~ bowen_18O_annual + 
                   (1|Site), dat.sub, REML = TRUE)
H1.3.red.sp <- lmer(iso_18O_all ~ HWT_end_season + 
                      (1|Site), dat.sub, REML = TRUE)
H1.3.null <- lmer(iso_18O_all ~ 1 + 
                    (1|Site), dat.sub, REML = TRUE)
O2.res <- Anova(H1.3.full, type = 2, test.statistic = "F")
O2.res$Df.res <- round(O2.res$Df.res, 2)

# Check model save coefs
plot(H1.3.full)
summary(H1.3.full)
betas <- do.call(paste, c(as.data.frame(round(summary(H1.3.full)$coefficients[,1:2], digits=5)), sep="+-"))

# variance explained
# modelled precip 18O explaining between site variation
Oexp.pred18O.hwt <- as.numeric((VarCorr(H1.3.red.sp)$Site-VarCorr(H1.3.full)$Site) / VarCorr(H1.3.null)$Site)
Oexp.pred18O.hwt.tot <- rsquared(H1.3.full)[5] - rsquared(H1.3.red.sp)[5]

# explained variance HWT
Oexp.hwt.site <-as.numeric((VarCorr(H1.3.red)$Site-VarCorr(H1.3.full)$Site) / VarCorr(H1.3.null)$Site)
Oexp.hwt.tot <- rsquared(H1.3.full)[5] - rsquared(H1.3.red)[5]
# explained variance within site by HWT
(attributes(VarCorr(H1.3.red))$sc-attributes(VarCorr(H1.3.full))$sc) / attributes(VarCorr(H1.3.red))$sc #almost zero, 3%

# Fill table
Otable[10,"model"] <- deparse(formula(H1.3.full))
Otable$variable[10:(9+nrow(O2.res))] <- rownames(O2.res)
Otable$effect[10:(9+nrow(O2.res))] <- betas[-1]
Otable$Fvalue[10:(9+nrow(O2.res))] <- round(O2.res$F, digits = 2)
Otable$df[10:(9+nrow(O2.res))] <- do.call(paste, c(O2.res[c("Df", "Df.res")], sep=", "))
Otable$Pvalue[10:(9+nrow(O2.res))] <- round(O2.res$'Pr(>F)', digits = 4)
Otable$N[10] <- nrow(H1.3.full@frame)
Otable$varExpsite[10] <- round(Oexp.pred18O.hwt, digits = 3)
Otable$varExptotal[10] <- round(Oexp.pred18O.hwt.tot, digits = 3)[[1]]
Otable$varExpsite[11] <- round(Oexp.hwt.site, digits = 3)
Otable$varExptotal[11] <- round(Oexp.hwt.tot, digits = 3)[[1]]

# Check HWT part
Otable

# Test additional environmental variables

# make data set with S.magellanicum as reference level
# for later use
dat2 <- dat 
dat2$Species <- relevel(dat2$Species, ref = "S.magellanicum")

# start model is the same as earlier(Species + bowen_18O_JanApr) and 
# each additional variable is added one by one to this model

# add ET/P to start model
H1.4.start.ep <- lmer(iso_18O_all ~ ep*Species + bowen_18O_JanApr +
                        (1|Site), dat, na.action = "na.omit", REML = TRUE)
H1.4.start.red.no_ep <- update(H1.4.start.ep, .~ . -ep:Species -ep)
H1.4.null <- update(H1.4.start.ep, .~ . -ep*Species -bowen_18O_JanApr)

C3.res <- Anova(H1.4.start.ep, type =2, test.statistic = "F")
C3.res$Df.res <- round(C3.res$Df.res, 2)
summary(H1.4.start.ep)

# Possibly an interaction effect so refit to get slope for s.mag and s.fusc
betas <- do.call(paste, c(as.data.frame(round(summary(H1.4.start.ep)$coefficients[,1:2], digits=5)), sep="+-"))
betas[5] <- betas[2] # S. fuscum slope
betas[2] <- NA # no main effect of HWT
H1.4.start.ep2 <- update(H1.4.start.ep, .~. , data=dat2)
betas[5] <- paste("S.fus: ", betas[5], " / ",
                  "S.mag: ", do.call(paste, c(as.data.frame(round(summary(H1.4.start.ep2)$coefficients[,1:2], digits=5)), sep="+-"))[2],
                  sep = "")

# ET/P * species interaction
Oexp.pred18O.ep.int <- as.numeric((VarCorr(H1.4.start.red.no_ep)$Site[1]-VarCorr(H1.4.start.ep)$Site[1]) / VarCorr(H1.4.null)$Site[1])
Oexp.pred18O.int.ep.tot <- rsquared(H1.4.start.ep)[5]$Marginal - rsquared(H1.4.start.red.no_ep)[5]$Marginal

varExp <- round(cbind(c(Oexp.pred18O.ep.int), 
                      c(Oexp.pred18O.int.ep.tot)), digits = 4)

Otable[12,"model"] <- "start+(evap/precip)*species"
Otable$variable[12:13] <- c("evap/precip", "e/p * Species")
Otable$effect[13] <- betas[5]
Otable$Fvalue[12:13] <- round(C3.res$F, digits = 2)[c(1,4)]
Otable$df[12:13] <- do.call(paste, c(C3.res[c("Df", "Df.res")], sep=", "))[c(1,4)]
Otable$Pvalue[12:13] <- round(C3.res$'Pr(>F)', digits = 4)[c(1,4)]
Otable$N[12] <- nrow(H1.4.start.ep@frame)
Otable[13, c("varExpsite", "varExptotal")] <- varExp

# add precip to start model
H1.4.start.pre <- lmer(iso_18O_all ~ precip*Species + bowen_18O_JanApr +
                         (1|Site), dat, na.action = "na.omit", REML = TRUE)
H1.4.start.red.no_pre <- update(H1.4.start.pre, .~ . -precip:Species -precip)
C4.res <- Anova(H1.4.start.pre, type =2, test.statistic = "F")
C4.res$Df.res <- round(C4.res$Df.res, 2)
summary(H1.4.start.pre)

# no interaction so get main effect coef of precip
H1.4.start.pre_main <- update(H1.4.start.pre, .~ . - precip:Species)
H1.4.start.pre_main_no_pre <- update(H1.4.start.pre_main, .~ . - precip)

betas <- do.call(paste, c(as.data.frame(round(summary(H1.4.start.pre_main)$coefficients[,1:2], digits=5)), sep="+-"))

# # explained variance by precipitation
Oexp.pred18O.pre <- as.numeric((VarCorr(H1.4.start.pre_main_no_pre)$Site[1]-VarCorr(H1.4.start.pre_main)$Site[1]) / 
                                 VarCorr(H1.4.null)$Site[1])
Oexp.pred18O.pre.tot <- (rsquared(H1.4.start.pre_main)[5] - rsquared(H1.4.start.pre_main_no_pre)[5])$Marginal

varExp <- round(cbind(c(Oexp.pred18O.pre), 
                      c(Oexp.pred18O.pre.tot)), digits = 4)

Otable[14,"model"] <- "start+precipitation*species"
Otable$variable[14:15] <- c("precip", "precip*Species")
Otable$effect[14] <- betas[2]
Otable$Fvalue[14:15] <- round(C4.res$F, digits = 2)[c(1,4)]
Otable$df[14:15] <- do.call(paste, c(C4.res[c("Df", "Df.res")], sep=", "))[c(1,4)]
Otable$Pvalue[14:15] <- round(C4.res$'Pr(>F)', digits = 4)[c(1,4)]
Otable$N[14] <- nrow(H1.4.start.ep@frame)
Otable[14, c("varExpsite", "varExptotal")] <- varExp

# add mean temperature to start model
H1.4.start.tmp <- lmer(iso_18O_all ~ avt*Species + bowen_18O_JanApr +
                         (1|Site), dat, na.action = "na.omit", REML = TRUE)
H1.4.start.red.no_tmp <- update(H1.4.start.tmp, .~ . -avt:Species -avt)
C5.res <- Anova(H1.4.start.tmp, type =2, test.statistic = "F")
C5.res$Df.res <- round(C5.res$Df.res, 2)
summary(H1.4.start.tmp)

# no interaction so get main effect coef for temp
H1.4.start.tmp_main <- update(H1.4.start.tmp, .~ . - avt:Species)
H1.4.start.tmp_main_no_tmp <- update(H1.4.start.tmp_main, .~ . - avt)

betas <- do.call(paste, c(as.data.frame(round(summary(H1.4.start.tmp_main)$coefficients[,1:2], digits=5)), sep="+-"))

# explained variance by temperature
Oexp.pred18O.tmp <- as.numeric((VarCorr(H1.4.start.tmp_main_no_tmp)$Site[1]-VarCorr(H1.4.start.tmp_main)$Site[1]) / 
                                 VarCorr(H1.4.null)$Site[1])
Oexp.pred18O.tmp.tot <- (rsquared(H1.4.start.tmp_main)[5] - rsquared(H1.4.start.tmp_main_no_tmp)[5])$Marginal

varExp <- round(cbind(c(Oexp.pred18O.tmp), 
                      c(Oexp.pred18O.tmp.tot)), digits = 4)

Otable[16,"model"] <- "start+temp*species"
Otable$variable[16:17] <- c("temp", "temp*Species")
Otable$effect[16] <- betas[2]
Otable$Fvalue[16:17] <- round(C5.res$F, digits = 2)[c(1,4)]
Otable$df[16:17] <- do.call(paste, c(C5.res[c("Df", "Df.res")], sep=", "))[c(1,4)]
Otable$Pvalue[16:17] <- round(C5.res$'Pr(>F)', digits = 4)[c(1,4)]
Otable$N[16] <- nrow(H1.4.start.tmp@frame)
Otable[16, c("varExpsite", "varExptotal")] <- varExp

# Save Table 3
write.csv(Otable, file = "table3_O.csv")


# Figure 3 ####

# Figure 3 shows model with main effect of Species
fig3.mod <- lmer(iso_18O_all ~ bowen_18O_annual + Species + 
                   (1|Site), dat, na.action = "na.omit")
plot(fig3.mod)
summary(fig3.mod)

# get equations for Fig 3
# S.fuscum
paste(fixef(fig3.mod)[1]," + ",fixef(fig3.mod)[2])
# S.magellanicum
paste((fixef(fig3.mod)[1] + fixef(fig3.mod)[3]), " + ", fixef(fig3.mod)[2])

################################
# Get BLUPs for plotting as an alternative. not shown in the paper.
#rf <- ranef(fig3.mod, condVar=TRUE)
#dotplot(rf)
#pred.df <- data.frame(aggregate(dat$bowen_18O_annual, by = list(dat$Site), mean))
#colnames(pred.df)  <- c("Site", "bowen_18O_annual")
#pred.df <- pred.df[pred.df$Site %in% rownames(rf$Site),]

# Make new data frame with BLUPs and their standard error.
#rff <- data.frame(pred.df, RF=rf$Site[,1],Var=attr(rf$Site,"postVar")[,,])
# compute the group level regression
# put NAs for species-site combination that werent sampled
#fusc.samp <- table(dat$Site, dat$Species)[,1]>0
#mag.samp <- table(dat$Site, dat$Species)[,2]>0

#rff$RFF.fusc <- with(rff, RF+fixef(fig3.mod)[1]+fixef(fig3.mod)[2]*bowen_18O_annual)
#rff$RFF.fusc <- ifelse(fusc.samp, rff$RFF.fusc, NA)
#rff$RFF.mag <- with(rff,RF+fixef(fig3.mod)[1]+fixef(fig3.mod)[2]*bowen_18O_annual+fixef(fig3.mod)[3])
#rff$RFF.mag <- ifelse(mag.samp, rff$RFF.mag, NA)
#rff <- reshape(rff, varying = c("RFF.fusc", "RFF.mag"), direction = "long")
################################

# plot raw data with SDs
agg.dat <- aggregate(cbind(iso_18O_all, bowen_18O_annual) ~ Site + Species, dat, FUN = function (x) c(mean(x), sd(x)))
agg.dat <- data.frame(agg.dat[,1:2], agg.dat[,3], agg.dat[,4][,1])
colnames(agg.dat)[3:5] <- c("iso_18O_all", "se", "bowen_18O_annual")

# To get confidence interval bands
# We have a lot of data so 2 x SE works well as 95% CIs
mima <- aggregate(bowen_18O_annual ~ Species, data=dat, function (x) c(min(x),max(x)) )$bowen_18O_annual
newdat <- rbind(data.frame(Species = "S.fuscum", bowen_18O_annual = seq(mima[1, 1],mima[1, 2], 0.025)),
                data.frame(Species = "S.magellanicum", bowen_18O_annual = seq(mima[2, 1],mima[2, 2], 0.025)))
newdat$predMean <- predict(fig3.mod, newdata=newdat, re.form=~0)

mm <- model.matrix(formula('~bowen_18O_annual + Species'),newdat)
newdat$distance <- mm %*% fixef(fig3.mod)
newdat$pvar1 <- diag(mm %*% tcrossprod(vcov(fig3.mod),mm))
newdat <- data.frame(
  newdat
  , plo = newdat$distance-2*sqrt(newdat$pvar1)
  , phi = newdat$distance+2*sqrt(newdat$pvar1) )
#

# Plot figure 3
fig3 <- ggplot(agg.dat, aes(x=bowen_18O_annual, y=iso_18O_all, ymin=iso_18O_all-se, ymax=iso_18O_all+se, color = Species, shape = Species)) +
  geom_point(size=2)+ geom_linerange(show.legend = FALSE) +
  geom_ribbon(data=newdat, aes(y=plo, x=bowen_18O_annual, ymin = plo, ymax = phi, fill = Species), 
              alpha = .15, show.legend = FALSE, linetype=0, inherit.aes = FALSE) +
  geom_line(data=newdat, aes(y=predMean, x=bowen_18O_annual, color=Species), 
            alpha =1, size = 1.2,inherit.aes = FALSE) +
  labs(y=expression(paste("Tissue ", delta^{18}, "O (\u2030)")), 
       x= expression(paste("Modelled precipitation ", delta^{18}, "O (\u2030)"))) +
  scale_fill_manual(values=c("black", "red"), name="fill") +
  scale_color_manual(name = "Species", breaks = c("S.fuscum", "S.magellanicum"), values = c("S.fuscum" = "black", "S.magellanicum" = "red"),
                     labels = c("S. fuscum", "S. magellanicum")) +
  scale_shape_manual(name=c("Species"), breaks = c("S.fuscum", "S.magellanicum"),values = c(1,2), labels = c("S. fuscum", "S. magellanicum")) +
  theme_bw()+
  theme(legend.justification = c(0, 0), 
        legend.position = c(0.65, 0.03),
        legend.key.size = unit(2, "line"),
        legend.text = element_text(size=12, face = "italic"),
        legend.title = element_text(size=12)) #+

fig3
# legend line and point size are not equal and here is a hack
# to fix this
library(grid)
grid.ls(grid.force()) # check slots to edit 
grid.gedit("key-[-0-9]-1-1.4-2-4-2", size = unit(3, "mm"))    
grid.gedit("key-[-0-9]-1-1.5-2-5-2", size = unit(3, "mm"))    
fig3 <- grid.grab() # save edits
ggsave("fig3_18o_ci_raw.png", fig3,  dpi = 300)
# save as eps
cairo_pdf(filename = "fig3_18o_ci_raw.eps", height = 5, width=6) 
grid.draw(fig3)
dev.off()

# Species differences in HWT ####
# use nlme for weights function (account for increasing variance)
detach(package:lme4)
library(nlme)
H1.2 <- lme(HWT_end_season ~ Species, 
            random = ~ 1|Site, dat, weights = varPower(),na.action = "na.omit", method = "REML")
anova(H1.2)
plot(H1.2)
summary(H1.2) # 11 cm lower for S.mag on average (33 cm vs 22 cm)
detach(package:nlme)

# Table S1 -  data overview ####
range(dat$elevation)
range(dat$Patch_coord_lat)

tables1 <- aggregate(cbind(Patch_coord_lat, Patch_coord_lon, elevation, samp_year, Species) ~ Region + Site, dat, mean)
colnames(tables1)[7] <- "species"
tables1$species <-ifelse(tables1$species == 1, "S.fuscum",
                         ifelse(tables1$species == 2, "S.magellanicum","S.fuscum+S.magellanicum")) 

# one site had samples from two years
tables1$samp_year <- ifelse(tables1$samp_year != 2013 & tables1$samp_year != 2014, "2013 & 2014", tables1$samp_year)
colnames(tables1) <- c("region", "site", "latitude_wgs84", "longitude_wgs84", "elevation_masl", "sample_year", "species") 
tables1$latitude_wgs84 <- round(tables1$latitude_wgs84, 5)
tables1 <- tables1[order(tables1$region),]
write.csv(tables1, file="tables1rr.csv", row.names = FALSE)

# Table S2 - Within-between site variation ####
cols <- c("iso_18O_all", "iso_13C_all", "bowen_18O_annual", "avt", 
          "evap", "precip", "ep", "prod", "HWT_end_season", "elevation")
sumstat <- aggregate(cbind(iso_18O_all, iso_13C_all, bowen_18O_annual, avt, 
                           evap, precip, ep, prod, HWT_end_season, elevation) ~ Species, data=dat, function (x) c(mean(x), sd(x), range(x)))
sumstat <- do.call(data.frame, sumstat)

tables2 <- reshape(sumstat, varying = c(colnames(sumstat)[-1]), direction="long", timevar="var", times=cols,
                   v.names=c("1mean", "2sd", "3min", "4max"), idvar= "Species")
rownames(tables2)<-NULL
tables2 <- tables2[,c(2,1,3:6)]
colnames(tables2)[3:6] <- c("mean", "sd", "min", "max")
tables2[,3:6] <- round(tables2[,3:6], digits = 2)
write.csv(tables2, file = "tables2.csv", row.names = FALSE)

# Figure S1 - d18O pred per month ####
# check variation in pred 18O over the year
predO <- read.csv("d18O_pred_months.csv")
predO <- predO[,-c(1:4)] # remove site info and coords
predO$site <- factor(1:nrow(predO))
# change to long data format
predO.long <- reshape(predO, idvar = "site", varying = list(1:12),
                    v.names = "mod18O", direction = "long")
predO.long$time <- as.Date(paste("2013",predO.long$time,"01", sep="-"))

# check SD for each months
apply(predO,2, sd, na.rm=TRUE)

# plot
preO <- ggplot(predO.long[order(predO.long$site),], aes(y=mod18O, x = time, group=site)) +
  geom_line() +
  geom_point() +
  labs(y = expression(paste("Modelled precipitation ", delta^{18}, "O (\u2030)")),
       x = "Month") +
  scale_x_date(date_labels = "%b",date_breaks = "months")
ggsave(preO, file = "pred18Oseason.png")