library("metafor")
data("dat.bcg", package = "metafor")
dat <- escalc(measure = "RR", ai = tpos, bi = tneg, ci = cpos, 
              di = cneg, data = dat.bcg, append = TRUE)
k <- length(dat.bcg$trial)
dat.fm <- data.frame(study = factor(rep(1:k, each = 4)))
dat.fm$grp <- factor(rep(c("T", "T", "C", "C"), k), levels = c("T", "C"))
dat.fm$out <- factor(rep(c("+", "-", "+", "-"), k), levels = c("+", "-"))
dat.fm$freq <- with(dat.bcg, c(rbind(tpos, tneg, cpos, cneg)))
dat$a.random     <- ifelse(dat$alloc == "random",     1, 0)
dat$a.alternate  <- ifelse(dat$alloc == "alternate",  1, 0)
dat$a.systematic <- ifelse(dat$alloc == "systematic", 1, 0)
res <- rma(yi, vi, mods = cbind(ablat, year), data = dat)
permres <- permutest(res, iter = 10000, retpermdist = TRUE)
hist(permres$zval.perm[,2], breaks = 140, freq = FALSE, 
     xlim = c(-5, 5), ylim = c(0, .4), main = "", 
     xlab = "Value of Test Statistic")
abline(v = res$zval[2], lwd = 2, lty = "dashed")
abline(v = 0, lwd = 2)
curve(dnorm, from = -5, to = 5, add = TRUE, lwd = 2, 
      col = rgb(1, 0, 0, alpha = .7))
lines(density(permres$zval.perm[,2]), lwd = 2, 
      col = rgb(0, 0, 1, alpha = .7))