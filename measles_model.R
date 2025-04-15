# measles model

# stochastic SEIR  model 
SEIR.onestep <- function (x, params) { #function to calculate one step of stochastic SEIR
  S <- x[2] #local variable for susceptible
  E <- x[3] #local variable for exposed
  I <- x[4] #local variable for infected
  R <- x[5] #local variable for recovered
  ## does something need to happen in the initial Q to account for qi>0 --A: no. 
  Qs <- x[6] # quarantined people who won't end up in R (for now, completely) 
  Qr <- x[7] # quarantined people who WILL end up in R (infected) 
  #      N <- X+Y+Z+R #total population size (subject to demographic change)
  with( #use with to simplify code
    as.list(params),
    {
      total.rate <-  beta*(I+ c*E)*S + v*S+ qs*S + k*E+ qspep*E+ gamma*I+qi*I+l*Qs + l*Qr
      if (total.rate > 0) { # events : new infection (S to S-1, E to E+1), vax an S (S to S-1, R to R+1), quar an S (S-1, Q+1)
        # progress an E (E-1, I+1) , quar an E (E-1, Q+1), rec an I (I-1, R+1), quar an I (I-1, Q+1), 
        # release a Qs (Qs-1, S+1), release a Qr (Qr-1, S+1) 
        tau <- rexp(n=1,rate=total.rate) #inter-event time
        new.xyz <- c(S,E,I,R,Qs,Qr) #initialize a local variable at previous state variable values (not strictly necessary)
        U <- runif(1) # uniform random deviate
        #             new.xyz<-c(X,Y,Z-1) 
        ## the thing below is called the same as the thing above... which one is it?? -A: holdover from the code i startd with; redundant
        new.xyz <- c(S, E, I, R+1, Qs, Qr-1) # last event is release a Qr (to R) (I removed waning an R, immunity lasts too long) 
        # for each event, if U < (sum up to that one) we say we are going to do that one. If none of the other ifs are true
        # that's the one that happens. this results in each event having the correct probability. 
        if (U<=(beta*(I+ c*E)*S + v*S+ qs*S + k*E+ qspep*E+ gamma*I+qi*I +l*Qs)/total.rate) new.xyz <- c(S+1, E, I, R, Qs-1, Qr) # release a Qs
        if (U<=(beta*(I+ c*E)*S + v*S+ qs*S + k*E+ qspep*E+ gamma*I+qi*I )/total.rate) new.xyz <- c(S, E, I-1, R, Qs, Qr+1) # quar an I
        if (U<=(beta*(I+ c*E)*S + v*S+ qs*S + k*E+ qspep*E+ gamma*I)/total.rate) new.xyz <- c(S, E, I-1, R+1, Qs, Qr) # rec an I
        if (U<=(beta*(I+ c*E)*S + v*S+ qs*S + k*E+ qspep*E)/total.rate) new.xyz <- c(S, E-1, I, R, Qs, Qr+1) # quar or PEP an E 
        if (U<=(beta*(I+ c*E)*S + v*S+ qs*S + k*E)/total.rate) new.xyz <- c(S, E-1, I+1, R,  Qs, Qr) # progress an E
        if (U<=(beta*(I+ c*E)*S + v*S+ qs*S )/total.rate) new.xyz <- c(S-1, E, I, R, Qs+1, Qr) # quar an S
        if (U<=(beta*(I+ c*E)*S + v*S)/total.rate) new.xyz <- c(S-1, E, I, R+1,  Qs, Qr) # vax an S
        if (U<=(beta*(I+ c*E)*S )/total.rate) new.xyz <- c(S-1, E+1, I, R,  Qs, Qr) # new infection
        c(tau,new.xyz) #store result
        
      } else { 
        return(NA) } 
    }
  )
}

# iterate the onestep function to simulate an outbreak (once)  
SEIR.model <- function (x, params, nstep) { #function to simulate stochastic SIR
  output <- data.frame("time" = 0,"S"=x[2], "E"=x[3], "I"=x[4], "R"=x[5], "Qs"=x[6], "Qr"=x[7],"ctime"=0,"day"=0,"inci"=0)
  ctime = inci = day = vector(mode = "numeric","length"=nstep) #vectors to save our calculated things...maybe add more later like for the converted time
  for (k in 1:nstep) { #iterate for nstep steps -- k is inside SEIR.onestep.......
    if (x[3] == 0 & x[4] == 0){break} #if E and I are 0 then we done 
    x <- SEIR.onestep(x,params)
    ctime[k] <- x[1] + sum(output[1:k,1])  # compute cumulative time 
    day[k] <- floor(ctime[k]) #time in days
    inci[k] <- as.numeric(output$E[k] < x[3]) # x[3] = E; this flags if there was an incident we will add it later
    ## when Q is turned on the else gives an error on the first iter
    if (any(is.na(x)) ) {
      break 
    } else {output[k+1,] <- c(x,ctime[k],day[k],inci[k])} 
  }
  output = output[ which(rowSums(output)>0), ] # only keep rows where some state was nonzero
}

# runs the simulations 
measles.sim <- function(vax.rate,pop.size,I0,pars,nsims){
  #pop.size #total population size
  
  VacFraction = vax.rate # -- VAX DATA WILL GO HERE -
  S0 <- round((1-VacFraction)*pop.size) # initial number susceptible 
  xstart <- c(time=0, S=S0, E=0, I = I0, R = pop.size-S0-I0, Qs=0, Qr=0) #initial conditions
  # R0 should be 12-18 in the absence of any qs etc, let's use that to set beta 
  #  R0=15; and in my model, R0 = N beta (1/(gamma+qi) ( k/(k+qs)), or if q=0, simply R0=beta/gamma, 
  gamma = 1/4
  b= 15*gamma/pop.size # beta = R0*gamma/N i think
  nstep = 50000
  if (is.null(pars$k)) {pars$k = 1/10} # for backwards compatibility
  params <- list(beta = b,
                 c=pars$c,
                 v=pars$v, # if on: 0.05 ( 0.01-0.1) rate of vaccination of S . makes a diff! 
                 qs = pars$qs, # does not make much diff if on: 0.06 (qs = 0.014 - 0.125 ) rate we find and quarantine susceptible people 
                 qspep = pars$qspep, # if on: 2/3 qs quarantine and/or PEP for exposed people
                 qi=pars$qi, # if on: 0.45 ( 0.2-0.72)  quarantine for infectious people (send home/isolate)
                 l=1/15, # mean duration of quarantine is 21 days but people do it imperfectly but some are infectious, gah! 
                 k=1/10, # mean E duration of 6 days before infectiousness
                 gamma=gamma) # 8 day infectiousness wo the qi  ) # parameters
  
  data <- vector(mode='list',length=nsims) #initialize list to store the output
  set.seed(12345)
  for (k in 1:nsims) { 
    data[[k]] <- as.data.frame(SEIR.model(xstart,params,nstep))
  }
  
  data <- bind_rows(data, .id="simnum")
  data <- data.frame(data,vax = as.factor(rep(as.character(vax.rate))))
  
  return(data)
}


### dont use these anymore !!! #####
# here's a function that harmonizes the time to 0, 1, 2, .. n days
# this function works on one entry of 'data' , whose cols are "time"  "S"     "E"     "I"     "R"     "Qs" "Qr"  "ctime" 
# cols 2:6 of the input have to be S E I R Q anad there has to be a col called ctime
convtime <- function(outk) { 
  t=0:ceiling(max(outk$ctime)) 
  res=data.frame(time=t, S=0, E=0, I=0, R=0, Qs=0, Qr=0)
  for (k in 2:length(t)) { 
    ind = max(which(outk$ctime <= t[k]))
    #  ind = last(filter(outk,ctime <= t[k]))
    res[k, 2:(ncol(res))] = outk[which(outk$ctime <= t[k]), 2:7] # NOTE hard-coded variable cols , if 2:7 isn't SEIRQs Qr there is a problem
  }
  return(res)
}


# get the incidence from the outbreak sims. Incidence: new cases. 
# this works on the data[[k]] output, and on the analogous data frames 
# after we convert the time to be 0,1,2, 3, instead of event times 
addincidence <- function(outk) {
  outk$incid = 0
  newCases = diff(outk$E)
  whereIncid= 1+ which(newCases >= 1) # indices where there are incident cases
  outk$incid[whereIncid] = newCases[whereIncid-1] # put 'em in 
  return(outk)
}
