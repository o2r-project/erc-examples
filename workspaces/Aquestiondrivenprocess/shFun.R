#Function runs coupled version of reservoir operations model
#and allows users to test the impact of different degrees of
#hedging and different streamflow sequences

shFun <- function(Q, TP, SU, V0, P0, D0, M0, A1, sigmaT, 
                  Kp, etaH, deltaB, deltaD, deltaI, 
                  deltaE, muS, alpha, beta, Vmax, Dmin, Mp) {
 
  #Preallocate varaible memory  
  V <- numeric(length = TP)         #Reservoir storage vector
  P <- numeric(length = TP)         #Population vector
  D <- numeric(length = TP)         #Per capita demand vector
  A <- numeric(length = TP)         #Surface area vector
  W <- numeric(length = TP)         #Withdrawal vector
  M <- numeric(length = TP)         #Shortage Mem vector
  S <- numeric(length = TP)         #Shortage vector
  dVdt <- numeric(length = TP)      #Volume change vector
  dPdt <- numeric(length = TP)      #Pop change vector
  dDdt <- numeric(length = TP)      #Demand change vector
  dMdt <- numeric(length = TP)      #Shortage mem. change vector
  
  #Intial Conditions
  t <- 1                  #Intial time [T]
  dt <- 1                 #Model time step
  V[t:SU] <- V0           #Initial storage
  W[t:SU] <- 0            #Initial storag
  P[t:SU] <- P0           #Initial population
  D[t:SU] <- D0           #Initial per cap. demand
  A[t:SU] <- A1           #Reservoir base area [L^2]
  S[t:SU] <- 0            #Initial shortage
  
  
  #Set up Stage storage curve
  Elev = numeric(length = 101)
  Elev = seq(-0.05, 0.5, length.out = 100)
  Vol = numeric(length = 100)
  i = 1;
  while (i < 101) {
    Vol[i] = Elev[i]/3*(A[1] + pi/4*(sigmaT*Elev[i] + (A[1]*4/pi)^0.5)^2 +
             (A[1]*pi/4*(sigmaT*Elev[i] + (A[1]*4/pi)^0.5)^2)^0.5)
    i = i + 1
  }
  elevmax <- max(Elev)
  elevmin <- min(Elev)
  
  #Simulate system dynamics
  t <- SU
  while (t < TP + 1) {
    InterpRes <- approx(Vol, Elev, V[t], method = "linear", n = 50,
                        elevmin, elevmax, rule = 2, f = 0, ties = mean)
    H  <- InterpRes$y #Current water height  
    A[t]  <-  pi*((A[1]/pi)^0.5 + sigmaT*H)^2    #area given elevation
    
    #Determine withdrawal decision
    if (V[t] + 0.25*Q[t-1] > D[t]*P[t] + Vmax) {
      W[t] <- V[t] + 0.25*Q[t-1] - Vmax;
    } else {
      if (V[t] + 0.25*Q[t-1] > Kp*D[t]*P[t]) {
        W[t] <- D[t]*P[t]
      } else {
        W[t] <- (V[t]+0.25*Q[t-1])/Kp
      }
    }

    #Is there a shortage?
    if (D[t]*P[t] > W[t]) {
      S[t] <- D[t]*P[t] - W[t]
    } else {
      S[t] <- 0
    }
       
    dVdt[t] <- Q[t] - W[t] - etaH*A[t] - 0.25*Q[t] - 0.5*Q[t]
    if (M[t] < Mp) {
      dPdt[t] <- (deltaB - deltaD + deltaI - deltaE)*P[t] 
    } else {
      dPdt[t] <- (deltaB - deltaD)*P[t] + (deltaI*(1-M[t]) - deltaE*M[t])*P[t]
    }
    dMdt[t] <- (S[t]/(D[t]*P[t]))^2*(1-M[t]) - muS*M[t]
    dDdt[t] <- -D[t]*(M[t]*alpha*(1 - Dmin/D[t]) + beta)
    
    V[t+1] <- max((dVdt[t]*dt + V[t]),0)
    P[t+1] <- dPdt[t]*dt + P[t]
    M[t+1] <- dMdt[t]*dt + M[t]
    D[t+1] <- dDdt[t]*dt + D[t]
    t <- t + 1
    
  } #close the while loop
  #create result matrix
  Tot  <-  D*P
  #result = matrix(nrow = TP+1, ncol = 7)
  Year  <- c(-4:55)
  result <- data.frame(Year[1:TP],Q[1:TP],V[1:TP],M[1:TP],D[1:TP],P[1:TP],Tot[1:TP])
  colnames(result) <- c("Year","Stream","Vol","Mem","Dem","Pop","Tot")
  return(result) 
} #close the function loop
