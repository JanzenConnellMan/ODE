filepath="~/Validation/Outputs"
setwd(filepath)

install.packages("deSolve",dependencies = TRUE, INSTALL_opts = '--no-lock')

# Comparison plots 

#Packages 
library(deSolve)

#ODE Models 

ShDiv <- function(Dist){
  
  
  -sum( Dist*log(Dist) )
  
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
    Local.Mort  <- 1+  pi*P*(1 + exp((r/(sqrt(D))) ))/( (1 - exp( (r/(sqrt(D))) )) )^2
    NLocal.Mort <-  pi*P*(1 + exp((r/(sqrt(D))) ))/( (1 - exp((r/(sqrt(D))) )) )^2
    
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

# Density of Species 
D <- .2


Alpha.List  <- seq(.1,5,length=50)
vlist       <- seq(5,10,length=50)

vlist      <- vlist[31:40]

Divs_FA <- data.frame(matrix(NA,nrow=length(Alpha.List),ncol=length(vlist)))

Divs_FA_SH <- data.frame(matrix(NA,nrow=length(Alpha.List),ncol=length(vlist)))


set.seed(100)
Fec.Vals   <- rlnorm(Sp,mean=0,sd=1)  


for(ii in 1:length(vlist)){
  
  r                   <- 1/vlist[ii]
  Neighborhood.Radius <- vlist[ii]*sqrt(2)
  
  Divs_FAV <- rep(NA,length(Alpha.List))
  
  Divs_FAVSH  <- rep(NA,length(Alpha.List))
  
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
    
    out.MDA1 <- ode(y=init, times=times, func=Model.Decay.Add_NIr, parms=parameters,method="ode45")
    #out.MDA2 <- ode(y=init, times=times, func=Model.Decay.Add_NIr2, parms=parameters,method="ode45")
    #   out.MDA3 <- ode(y=init, times=times, func=Model.NearNeigh, parms=parameters,method="ode45")
    #    out.MDA4 <- ode(y=init, times=times, func=Model.Fixed.Add, parms=parameters,method="ode45")
    
    Finals1   <- out.MDA1[length(times),2:(Sp+1)]
    Finals1   <- Finals1[Finals1>exp(-11)]
    Divs_FAV[jj] <- length(Finals1)
    Divs_FAVSH[jj] <- ShDiv(Finals1)
    
    
    #  Finals2   <- out.MDA2[length(times),2:(Sp+1)]
    #  Finals2   <- Finals2[Finals2> exp(-11)]
    #  Divs_DDV[jj] <- length(Finals2)
    #  Divs_DDVSH[jj] <- ShDiv(Finals2)
    
    
    #    Finals3   <- out.MDA3[length(times),2:(Sp+1)]
    #    Finals3   <- Finals3[Finals3>exp(-11)]
    #    Divs_DDV2[jj] <- length(Finals3)
    #    Divs_DDV2SH[jj] <- ShDiv(Finals3)
    
    
    #    Finals4   <- out.MDA4[length(times),2:(Sp+1)]
    #    Finals4   <- Finals4[Finals4> exp(-11)]
    #    Divs_FAV2[jj] <- length(Finals4)
    #    Divs_FAV2SH[jj] <- ShDiv(Finals4)
    
    
  }
  
  Divs_FA[,ii]    <- Divs_FAV
  #  Divs_FA2[,ii]   <- Divs_FAV2
  # Divs_DD[,ii]    <- Divs_DDV
  #  Divs_DD2[,ii]   <- Divs_DDV2
  
  Divs_FA_SH[,ii]   <- Divs_FAVSH
  #  Divs_FA2_SH[,ii]   <- Divs_FAV2SH
  #  Divs_DD_SH[,ii]   <- Divs_DDVSH
  #  Divs_DD2_SH[,ii]  <- Divs_DDV2SH
  
  
}



Divs_FA <- as.matrix(Divs_FA)
Divs_FA_SH <- as.matrix(Divs_FA_SH)


write.csv(Divs_FA,"ODE_NewAdd_Y10_D.csv",quote=F,row.names=F)
write.csv(Divs_FA_SH,"ODE_NewAdd_Y10_SH_D.csv",quote=F,row.names=F)



