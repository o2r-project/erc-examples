#' @get /figure1
#' @png 
function(newValue){ 
  newValue = as.numeric(newValue)
  require(devtools)
  options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)
  library("mice")
  library("VIM")
  imp <- mice(nhanes, print=FALSE)
  pred <- imp$predictorMatrix
  pred[,"bmi"] <- 0
  imp <- mice(nhanes, pred=pred, pri=FALSE)
  ini <- mice(nhanes2, maxit=0, pri=FALSE)
  pred <- ini$pred
  pred[,"bmi"] <- 0
  meth <- ini$meth
  meth["bmi"] <- ""
  imp <- mice(nhanes2, meth=meth, pred=pred, pri=FALSE)
  pred <- ini$pred
  pred["bmi",] <- 0
  imp <- mice(nhanes2, pred=pred, pri=FALSE, seed=51162)
  ini <- mice(popmis, maxit=0)
  pred <- ini$pred
  pred["popular",] <- c(0, -2, 0, 2, 1, 2, 0)
  imp <- mice(popmis, meth=c("","","2l.norm","","","",""), pred=pred, maxit=1, seed=71152)
  imp <- mice(nhanes, pred=quickpred(nhanes, minpuc=0.25, include="age"))
  nhanes2.ext <- cbind(nhanes2, lchl=log(nhanes2$chl))
  ini <- mice(nhanes2.ext, max=0, print=FALSE)
  meth <- ini$meth
  meth["lchl"] <- "~log(chl)"
  pred <- ini$pred
  pred[c("hyp","chl"),"lchl"] <- 0
  pred["bmi","chl"] <- 0
  imp <- mice(nhanes2.ext, meth=meth, pred=pred, seed=38788, print=FALSE)
  nhanes2.ext <- cbind(nhanes2, lchl=NA)
  ini <- mice(boys,max=0,print=FALSE)
  meth <- ini$meth
  meth["bmi"] <- "~I(wgt/(hgt/100)^2)"
  pred <- ini$pred
  pred[c("wgt","hgt","hc","reg"),"bmi"] <- 0
  pred[c("gen","phb","tv"),c("hgt","wgt","hc")] <- 0
  imp.idx <- mice(boys, pred=pred, meth=meth, maxit=20, seed=09212, print=FALSE)
  meth["bmi"] <- "~round(wgt/(hgt/100)^2,dig=2)"
  ini <- mice(cbind(boys,mat=NA),max=0,print=FALSE)
  meth <- ini$meth
  meth["mat"] <- "~I(as.integer(gen) + as.integer(phb) +
  + as.integer(cut(tv,breaks=c(0,3,6,10,15,20,25))))"
  meth["bmi"] <- "~I(wgt/(hgt/100)^2)"
  pred <- ini$pred
  pred[c("bmi","gen","phb","tv"),"mat"] <- 0
  pred[c("hgt","wgt","hc","reg"),"mat"] <- 1
  pred[c("hgt","wgt","hc","reg"),c("gen","phb","tv")] <- 0
  pred[c("wgt","hgt","hc","reg"),"bmi"] <- 0
  pred[c("gen","phb","tv"),c("hgt","wgt","hc")] <- 0
  imp.sum <- mice(cbind(boys,mat=NA), pred=pred, meth=meth, maxit=20, seed=10948, print=FALSE)
  nhanes2.ext <- cbind(nhanes2, bmi.chl=NA)
  ini <- mice(nhanes2.ext, max=0, print=FALSE)
  meth <- ini$meth
  meth["bmi.chl"] <- "~I((bmi-25)*(chl-200))"
  pred <- ini$pred
  pred[c("bmi","chl"),"bmi.chl"] <- 0
  #imp <- mice(nhanes2.ext, meth=meth, pred=pred, seed=51600, print=FALSE)
  nhanes2.ext <- cbind(nhanes2,age.1.bmi=NA,age.2.bmi=NA)
  ini <- mice(nhanes2.ext, max=0, print=FALSE)
  meth <- ini$meth
  meth["age.1.bmi"] <- "~I(age.1*(bmi-25))"
  meth["age.2.bmi"] <- "~I(age.2*(bmi-25))"
  pred <- ini$pred
  pred[c("age","bmi"),c("age.1.bmi","age.2.bmi")] <- 0
  #imp <- mice(nhanes2.ext, meth=meth, pred=pred, maxit=10)
  nhanes2.ext <- cbind(nhanes2, lchl=NA)
  ini <- mice(nhanes2.ext,max=0,pri=FALSE)
  meth <- ini$meth
  meth[c("lchl","chl")] <- c("~log(chl)","norm")
  pred <- ini$pred
  pred[c("hyp","chl"),"lchl"] <- 0
  pred["bmi","chl"] <- 0
  #imp <- mice(nhanes2.ext, meth=meth, pred=pred, seed=1, maxit=100)
  meth["lchl"] <- "~log(squeeze(chl, bounds=c(100,300)))"
  #imp <- mice(nhanes2.ext, meth=meth, pred=pred, seed=1, maxit=100)
  nhanes2.ext <- cbind(nhanes2, lchl=NA)
  ini <- mice(nhanes2.ext,max=0,print=FALSE)
  meth <- ini$meth
  meth[c("lchl","chl")] <- c("~log(chl)","norm")
  pred <- ini$pred
  pred[c("hyp","chl"),"lchl"] <- 0
  pred["bmi","chl"] <- 0
  post <- ini$post
  post["chl"] <- "imp[[j]][,i] <- squeeze(imp[[j]][,i],c(100,300))"
  #imp <- mice(nhanes2.ext, meth=meth, pred=pred, post=post, 
  #            seed=30031, maxit=10, print=FALSE)
  ini <- mice(cbind(boys,mat=NA),max=0,print=FALSE)
  meth <- ini$meth
  meth["mat"] <- "~I(as.integer(gen) + as.integer(phb) +
  + as.integer(cut(tv,breaks=c(0,3,6,10,15,20,25))))"
  meth["bmi"] <- "~I(wgt/(hgt/100)^2)"
  pred <- ini$pred
  pred[c("bmi","gen","phb","tv"),"mat"] <- 0
  pred[c("hgt","wgt","hc","reg"),"mat"] <- 1
  pred[c("hgt","wgt","hc","reg"),c("gen","phb","tv")] <- 0
  pred[c("wgt","hgt","hc","reg"),"bmi"] <- 0
  pred[c("gen","phb","tv"),c("hgt","wgt","hc")] <- 0
  post <- ini$post
  post["gen"] <- "imp[[j]][p$data$age[!r[,j]]<5,i] <- levels(boys$gen)[1]"
  post["phb"] <- "imp[[j]][p$data$age[!r[,j]]<5,i] <- levels(boys$phb)[1]"
  post["tv"]  <- "imp[[j]][p$data$age[!r[,j]]<5,i] <- 1"
  #imp <- mice(cbind(boys,mat=NA), pred=pred, meth=meth, post=post, maxit=10, print=FALSE)
  nhanes2.ext <- cbind(nhanes2, bmi.chl=NA)
  ini <- mice(nhanes2.ext, max=0, print=FALSE)
  meth <- ini$meth
  meth["bmi.chl"] <- "~I((bmi-25)*(chl-200))"
  pred <- ini$pred
  pred[c("bmi","chl"),"bmi.chl"] <- 0
  imp <- mice(nhanes2.ext, meth=meth, pred=pred, seed=51600, print=FALSE)
  vis<-imp$vis
  vis<-append(vis,vis[4],1)
  imp <- mice(nhanes2.ext, meth=meth, pred=pred, vis=vis)
  imp <- mice(nhanes2.ext, meth=meth, pred=pred, vis=c(2,4,5,3))
  imp <- mice(nhanes2.ext, meth=meth, pred=pred, vis="monotone")
  ini <- mice(nhanes2, maxit=0)
  import <- matrix(c(30,30,30,29,25,21,25,25,22,33,27,22,27,35,27,20,27,30), byrow=TRUE,nr=9)
  imp <- mice(nhanes, print=FALSE, seed=77172)
  imp$imp$bmi[,1:2] <- import
  imp <- mice(nhanes, maxit=4, seed=44612, print=FALSE)
  ini <- mice(boys,max=0,print=FALSE)
  meth <- ini$meth
  meth["bmi"] <- "~I(wgt/(hgt/100)^2)"
  meth["wgt"] <- "~I(bmi*(hgt/100)^2)"
  meth["hgt"] <- "~I(100*sqrt(wgt/bmi))"
  ini <- mice(boys,max=0,print=FALSE)
  meth <- ini$meth
  meth["bmi"] <- "~I(wgt/(hgt/100)^2)"
  m <- 5
  T <- newValue
  imp.kendall <- mice(boys, m=m, meth=imp.idx$meth, pred=imp.idx$pred, maxit=0, print=FALSE)
  tau <- matrix(NA,nrow=T,ncol=m)
  for (i in 1:T) {
    if(i==1) set.seed(09212)
    imp.kendall <- mice.mids(imp.kendall, maxit=1, print=FALSE)
    x <- complete(imp.kendall,"repeated")[,paste("gen",1:m,sep=".")]
    y <- complete(imp.kendall,"repeated")[,paste("phb",1:m,sep=".")]
    xn <- as.data.frame(lapply(x,as.numeric))
    yn <- as.data.frame(lapply(y,as.numeric))
    tau[i,] <- diag(cor(xn,yn,method="kendall"))
  }
  matplot(x=1:T,y=tau,xlab="Iteration",type="l")
}