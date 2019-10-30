library("metafor")
data("dat.bcg", package = "metafor")
dat <- escalc(measure = "RR", ai = tpos, bi = tneg, ci = cpos, 
              di = cneg, data = dat.bcg, append = TRUE)
k <- length(dat.bcg$trial)
dat.fm <- data.frame(study = factor(rep(1:k, each = 4)))
dat.fm$grp <- factor(rep(c("T", "T", "C", "C"), k), levels = c("T", "C"))
dat.fm$out <- factor(rep(c("+", "-", "+", "-"), k), levels = c("+", "-"))
dat.fm$freq <- with(dat.bcg, c(rbind(tpos, tneg, cpos, cneg)))
res <- rma(ai = tpos, bi = tneg, ci = cpos, 
           di = cneg, data = dat, measure = "RR")
forest(res, slab = paste(dat$author, dat$year, sep = ", "), 
       xlim = c(-16, 6), at = log(c(.05, .25, 1, 4)), atransf = exp,
       ilab = cbind(dat$tpos, dat$tneg, dat$cpos, dat$cneg), 
       ilab.xpos = c(-9.5, -8, -6, -4.5), cex = .75)
text(c(-9.5, -8, -6, -4.5), 15, c("TB+", "TB-", "TB+", "TB-"))
text(c(-8.75, -5.25),       16, c("Vaccinated", "Control"))
text(-16,                   15, "Author(s) and Year",     pos = 4)
text(6,                     15, "Relative Risk [95% CI]", pos = 2)