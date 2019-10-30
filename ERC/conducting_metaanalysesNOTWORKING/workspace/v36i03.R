### load the package and the data for the BCG vaccine meta-analysis 
#install.packages("metafor")
library("metafor")
data("dat.bcg", package = "metafor")

print(dat.bcg, row.names = FALSE)

### calculate the log relative risks and corresponding sampling variances

dat <- escalc(measure = "RR", ai = tpos, bi = tneg, ci = cpos, 
              di = cneg, data = dat.bcg, append = TRUE)
print(dat[,-c(4:7)], row.names = FALSE)

### demonstrate the formula interface

k <- length(dat.bcg$trial)
dat.fm <- data.frame(study = factor(rep(1:k, each = 4)))
dat.fm$grp <- factor(rep(c("T", "T", "C", "C"), k), levels = c("T", "C"))
dat.fm$out <- factor(rep(c("+", "-", "+", "-"), k), levels = c("+", "-"))
dat.fm$freq <- with(dat.bcg, c(rbind(tpos, tneg, cpos, cneg)))
dat.fm

escalc(out ~ grp | study, weights = freq, data = dat.fm, measure = "RR")

### random-effects model

res <- rma(yi, vi, data = dat)
res

res <- rma(ai = tpos, bi = tneg, ci = cpos, 
           di = cneg, data = dat, measure = "RR")
res

### confidence intervals for tau^2, I^2, and H^2

confint(res)

### forst plot

forest(res, slab = paste(dat$author, dat$year, sep = ", "), 
       xlim = c(-16, 6), at = log(c(.05, .25, 1, 4)), atransf = exp,
       ilab = cbind(dat$tpos, dat$tneg, dat$cpos, dat$cneg), 
       ilab.xpos = c(-9.5, -8, -6, -4.5), cex = .75)
op <- par(cex = .75, font = 2)
text(c(-9.5, -8, -6, -4.5), 15, c("TB+", "TB-", "TB+", "TB-"))
text(c(-8.75, -5.25),       16, c("Vaccinated", "Control"))
text(-16,                   15, "Author(s) and Year",     pos = 4)
text(6,                     15, "Relative Risk [95% CI]", pos = 2)
par(op)

### mixed-effects model with absolute latitude and publication year as moderators

res <- rma(yi, vi, mods = cbind(ablat, year), data = dat)
res

res <- rma(yi, vi, mods = ~ ablat + year, data = dat)
res

### predicted average relative risks for 10, 20, 30, 40, and 50 degrees 
### absolute latitude holding publication year constant at 1970

predict(res, newmods = cbind(seq(from = 10, to = 60, by = 10), 1970), 
        transf = exp, addx = TRUE)

### plot of the predicted average relative risk as a function of absolute latitude

preds <- predict(res, newmods = cbind(0:60, 1970), transf = exp)
wi    <- 1/sqrt(dat$vi)
size  <- 0.5 + 3 * (wi - min(wi))/(max(wi) - min(wi))
plot(dat$ablat, exp(dat$yi), pch = 19, cex = size, 
     xlab = "Absolute Latitude", ylab = "Relative Risk", 
     las = 1, bty = "l", log = "y")
lines(0:60, preds$pred)
lines(0:60, preds$ci.lb, lty = "dashed")
lines(0:60, preds$ci.ub, lty = "dashed")
abline(h = 1, lty = "dotted")

### separate random-effects model for the various levels of the allocation factor

rma(yi, vi, data = dat, subset = (alloc=="random"))
rma(yi, vi, data = dat, subset = (alloc=="alternate"))
rma(yi, vi, data = dat, subset = (alloc=="systematic"))

### dummy-coding of the allocation factor

dat$a.random     <- ifelse(dat$alloc == "random",     1, 0)
dat$a.alternate  <- ifelse(dat$alloc == "alternate",  1, 0)
dat$a.systematic <- ifelse(dat$alloc == "systematic", 1, 0)

### mixed-effects model with dummies for the allocation factor

rma(yi, vi, mods = cbind(a.random, a.alternate, a.systematic), 
    intercept = FALSE, data = dat)
rma(yi, vi, mods = ~ factor(alloc) - 1, data = dat)

### same model with different parameterization so that the allocation 
### factor can be tested with the omnibus test

rma(yi, vi, mods = cbind(a.alternate, a.systematic), data = dat)
rma(yi, vi, mods = ~ factor(alloc), data = dat)

### mixed-effects model with Knapp and Hartung adjustment

rma(yi, vi, mods = ~ factor(alloc) + ablat, data = dat, knha = TRUE)

### empirical Type I error rate for a random moderator without 
### and with the Knapp and Hartung adjustment

pval1 <- pval2 <- rep(NA, 10000)
for (i in 1:10000) {
   x <- rnorm(13)
   pval1[i] <- rma(yi, vi, mods = x, data = dat, method = "DL")$pval[2]
   pval2[i] <- rma(yi, vi, mods = x, data = dat, knha = TRUE, method = "DL")$pval[2]
}
mean(pval1 < .05)
mean(pval2 < .05)

### random-effects model using the log relative risks and then transformation
### of the estimated average log relative risk to the average relative risk

res <- rma(yi, vi, data = dat)
predict(res, transf = exp, digits = 2)

### externally standardized residuals

res <- rma(yi, vi, mods = cbind(ablat, year), data = dat)
rstudent(res)

### influential case diagnostics

res <- rma(yi, vi, mods = cbind(ablat, year), data = dat)
inf <- influence(res)
plot(inf, plotdfb = TRUE)

### leave-one-out analysis

res <- rma(yi, vi, data = dat)
leave1out(res, transf = exp, digits = 3)

### forest plot and addpoly function

forest(dat$yi, dat$vi, atransf = exp, ylim = c(-3.5, 16),
       at = log(c(.05, .25, 1, 4, 20)), xlim = c(-9, 7), 
       slab = paste(dat$author, dat$year, sep=", "))
res <- rma(yi, vi, mods = cbind(ablat), data = dat)
preds <- predict(res, newmods = c(10, 30, 50))
addpoly(preds$pred, sei = preds$se, atransf = exp, 
        mlab = c("10 Degrees", "30 Degrees", "50 Degrees"))
text(-9, 15, "Author(s) and Year",     pos = 4, font = 2)
text( 7, 15, "Relative Risk [95% CI]", pos = 2, font = 2)
abline(h = 0)

### funnel plots

res <- rma(yi, vi, data = dat)
funnel(res, main = "Random-Effects Model")
res <- rma(yi, vi, mods = cbind(ablat), data = dat)
funnel(res, main = "Mixed-Effects Model")

### radial plots

res <- rma(yi, vi, data = dat, method = "FE")
radial(res, main = "Fixed-Effects Model")
res <- rma(yi, vi, data = dat, method = "REML")
radial(res, main = "Random-Effects Model")

### qq-normal plots

res <- rma(yi, vi, data = dat)
qqnorm(res, main = "Random-Effects Model")
res <- rma(yi, vi, mods = cbind(ablat), data = dat)
qqnorm(res, main = "Mixed-Effects Model")

### regression tests for funnel plot asymmetry

res <- rma(yi, vi, data = dat)
regtest(res, model = "lm")
res <- rma(ai = tpos, bi = tneg, ci = cpos, di = cneg, data = dat, 
           measure = "RR", mods = cbind(ablat, year))
regtest(res, predictor = "ni")

### trim and fill method

res <- rma(yi, vi, data = dat, method = "FE")
rtf <- trimfill(res)
rtf
funnel(rtf)

### cumulative meta-analysis

res <- rma(yi, vi, data = dat, slab = paste(author, year, sep=", "))
rcu <- cumul(res, order = order(dat$year))
forest(rcu, xlim = c(-6, 3), atransf = exp, 
       at = log(c(.10, .25, .5, 1, 2)))
text(-6, 15, "Author(s) and Year",     pos = 4, font = 2)
text( 3, 15, "Relative Risk [95% CI]", pos = 2, font = 2)

### likelihood ratio test of the absolute latitude moderator

res1 <- rma(yi, vi, mods = cbind(ablat), data = dat, method = "ML")
res2 <- rma(yi, vi, data = dat, method = "ML")
anova(res1, res2)

### permutation tests

res <- rma(yi, vi, data = dat)
permutest(res, exact = TRUE)

res <- rma(yi, vi, mods = cbind(ablat, year), data = dat)
permres <- permutest(res, iter = 10000, retpermdist = TRUE)
permres

hist(permres$zval.perm[,2], breaks = 140, freq = FALSE, 
     xlim = c(-5, 5), ylim = c(0, .4), main = "", 
     xlab = "Value of Test Statistic")
abline(v = res$zval[2], lwd = 2, lty = "dashed")
abline(v = 0, lwd = 2)
curve(dnorm, from = -5, to = 5, add = TRUE, lwd = 2, 
      col = rgb(1, 0, 0, alpha = .7))
lines(density(permres$zval.perm[,2]), lwd = 2, 
      col = rgb(0, 0, 1, alpha = .7))

### best linear unbiased predictions in the random-effects model

res <- rma(yi, vi, data = dat)
blup(res)

### comparison of the general linear model approach with the 
### Mantel-Haenszel and Peto's method using odds ratios

res <- rma(ai = tpos, bi = tneg, ci = cpos, di = cneg, data = dat, measure = "OR", method = "FE")
res
predict(res, transf = exp)
rma.mh(ai = tpos, bi = tneg, ci = cpos, di = cneg, data = dat, measure = "OR")
rma.peto(ai = tpos, bi = tneg, ci = cpos, di = cneg, data = dat)

