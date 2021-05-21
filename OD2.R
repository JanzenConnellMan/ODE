filepath="~/Validation/Outputs"
setwd(filepath)

install.packages("deSolve")

# Comparison plots 

#Packages 
library(deSolve)

#ODE Models 

Model.Fixed.Area        <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    P = matrix(state[1:Sp], nrow = Sp, ncol=1)
    Y =  matrix(c(Fec.Vals), nrow = Sp, ncol=1) 
    m =  matrix(c(Mort.Rates), nrow = Sp, ncol=1) 
    d =  matrix(c(Disp.Rates), nrow = Sp, ncol=1)
    a =  matrix(c(Alpha.Vals), nrow = Sp, ncol=1)
    R = Neighborhood.Radius 
    
    C = D*pi*(R^2)
    
    
    P[P< exp(-12)] <- 0
    P <- P/sum(P)         
    
    
    Local.seeds         <- ((1-d)*Y + d*P*Y) *exp(-a)
    
    NLocal.seeds_NO     <- d*P*Y 
    
    NLocal.seeds_YES    <- d*P*Y*exp(-a)
    
    #Prob_NO             <- (1-P)^C
    Prob_NO             <- exp(-P*C)
    #Prob_YES            <- (1-(1-P)^C)
    Prob_YES            <- 1 - exp(-P*C)
    
    
    
    Sum.Nlocal          <- LotMat%*%(Prob_NO * NLocal.seeds_NO)  + LotMat%*%(Prob_YES * NLocal.seeds_YES)
    
    
    #Recolonization 
    Recol.Func          <- Local.seeds*P/ (Local.seeds  + Sum.Nlocal)
    
    #Colonization
    Col.Func_NO           <- LotMat%*%(P /(Local.seeds  + Sum.Nlocal)   ) * Prob_NO  * NLocal.seeds_NO
    Col.Func_YES          <- LotMat%*%(P /(Local.seeds  + Sum.Nlocal)   ) * Prob_YES * NLocal.seeds_YES
    
    
    
    
    dP <- m*(   Recol.Func  +  Col.Func_NO +  Col.Func_YES    - P  )
    
    
    return(list(c(dP)))
  })
}
Model.Decay.Add_NIr     <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # Set the important variables / and Species Proprtions 
    P =  matrix(state[1:Sp], nrow = Sp, ncol=1)
    Y =  matrix(c(Fec.Vals), nrow = Sp, ncol=1) 
    m =  matrix(c(Mort.Rates), nrow = Sp, ncol=1) 
    d =  matrix(c(Disp.Rates), nrow = Sp, ncol=1)
    a =  matrix(c(Alpha.Vals), nrow = Sp, ncol=1)
    
    # In case of ODE Errors, this normlizes Numbers       
    P[P< exp(-12)] <- 0
    P <- P/sum(P)         
    
    # Seed Morality Expressions 
    Local.Mort  <- (2*pi*D*P/(r^2)) +1
    NLocal.Mort <- (2*pi*D*P/(r^2)) 
    
    # Local seeds + mortality
    Local.seeds           <- Y*( (1-d) + P*d)*exp(-a*Local.Mort)
    
    # Nonlocal seeds + mortality 
    NLocal.seeds          <- d*P*Y*exp(-a*NLocal.Mort)
    
    # All Nonlocal Seeds 
    Sum.Nlocal            <- LotMat%*%(NLocal.seeds)
    
    # Recolonization
    Recol.Func            <- Local.seeds*P/ (Local.seeds  + Sum.Nlocal)
    
    # Colonization
    Col.Func              <- LotMat%*%(P /(Local.seeds  + Sum.Nlocal)   ) * NLocal.seeds
    
    # Differeces 
    dP <- m*(   Recol.Func  +  Col.Func - P  )
    
    
    return(list(c(dP)))
  })
}

####### Parameters #######
# (Parms for beta + R (Radius) Comparison )

# Number of Species 

Sp <- 300

# Matrix For Interactions
LotMat <- matrix(1,Sp,Sp)
diag(LotMat) <- 0

# Disturbance Rate
m <- 2.5
Mort.Rates  <- rep(m,Sp)  

# Dispersal Rate 
d <- 1
Disp.Rates  <- rep(d,Sp)


# Decay.Distance 
r <- .1

# Density of Species 
D <- .2

Neighborhood.Radius <- sqrt(24/(D*pi))

set.seed(150)
Alpha.List<- seq(.1,5,length=50)
SD.FEC <- seq(.1,1,length=50)

Divs_FA <- data.frame(matrix(NA,nrow=length(Alpha.List),ncol=length(SD.FEC)))
Divs_DD <- data.frame(matrix(NA,nrow=length(Alpha.List),ncol=length(SD.FEC)))

for(ii in 1:length(SD.FEC)){
  
  set.seed(100)
  Fec.Vals   <- rlnorm(Sp,mean=0,sd=SD.FEC[ii])  
  
  Divs_FAV <- rep(NA,length(Alpha.List))
  Divs_DDV <- rep(NA,length(Alpha.List))
  
  for(jj in 1:length(Alpha.List)){
    
    # Alphas 
    A.Max <- Alpha.List[jj]
    A.Min <- Alpha.List[jj]
    Alpha.Vals <- seq(A.Min,A.Max,length=Sp)
    
    byn        <- 5
    times      <- seq(0, 2500, by = byn)
    inits      <- rep(1/Sp,Sp)
    Ps <-  setNames(inits,  paste0(rep(LETTERS[16], each=Sp), rep(1:Sp, 1))   )
    init <- c(Ps)
    parameters <- c()
    
    #    out.MDA1 <- ode(y=init, times=times, func=Model.Fixed.Area, parms=parameters,method="ode45")
    out.MDA2 <- ode(y=init, times=times, func=Model.Decay.Add_NIr, parms=parameters,method="ode45")
    
    #    Finals1   <- out.MDA1[length(times),2:(Sp+1)]
    #    Finals1   <- Finals1[Finals1>exp(-11)]
    #    Finals1   <- Finals1[!is.na(Finals1)]
    #   Divs_FAV[jj] <- length(Finals1)
    
    
    Finals2   <- out.MDA2[length(times),2:(Sp+1)]
    Finals2   <- Finals2[Finals2> exp(-11)]
    #  Finals2   <- Finals2[!is.na(Finals2)]
    Divs_DDV[jj] <- length(Finals2)
  }
  
  Divs_FA[,ii]  <- Divs_FAV
  Divs_DD[,ii]  <- Divs_DDV
  
}


Divs_FA <- as.matrix(Divs_FA)
Divs_DD <- as.matrix(Divs_DD)



write.csv(Divs_FA,"ODE_Neigh_E24_d1.csv",quote=F,row.names=F)
write.csv(Divs_DD,"ODE_Add_v10_d1.csv",quote=F,row.names=F)



