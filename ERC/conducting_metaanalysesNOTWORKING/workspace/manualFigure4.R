options(repos = "https://cloud.r-project.org/")
library("metafor")
data("dat.bcg", package = "metafor")
dat <- escalc(measure = "RR", ai = tpos, bi = tneg, ci = cpos, 
              di = cneg, data = dat.bcg, append = TRUE)
k <- length(dat.bcg$trial)
dat.fm <- data.frame(study = factor(rep(1:k, each = 4)))
dat.fm$grp <- factor(rep(c("T", "T", "C", "C"), k), levels = c("T", "C"))
dat.fm$out <- factor(rep(c("+", "-", "+", "-"), k), levels = c("+", "-"))
dat.fm$freq <- with(dat.bcg, c(rbind(tpos, tneg, cpos, cneg)))
preds <- predict(res, newmods = cbind(0:60, 1970), transf = exp)
dat$a.random     <- ifelse(dat$alloc == "random",     1, 0)
dat$a.alternate  <- ifelse(dat$alloc == "alternate",  1, 0)
dat$a.systematic <- ifelse(dat$alloc == "systematic", 1, 0)
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