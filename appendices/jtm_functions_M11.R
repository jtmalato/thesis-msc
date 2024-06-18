###############################-
########## MODEL M11 ##########
###############################-

# 1) RCM model for restriction 1 -----
rcm.M11.equation <- function(t, par, tau) {
  
  lambda <- exp(par[1])
  rho1 <- exp(par[2])
  p <- exp(par[3])/(1+exp(par[3]))
  rho2 <- rho1*p
  
  theta1 <- lambda/(lambda+rho1)
  theta2 <- lambda/(lambda+rho2)
  
  if(t<=tau){
    
    x1 <- (1-exp(-(lambda+rho1)*t))
    prob <- theta1*x1
    
  } else {
    
    x1 <- (1-exp(-(lambda+rho1)*tau))
    x2 <- (1-exp(-(lambda+rho2)*(t-tau)))
    prob <- theta2*x2 + theta1*x1*(1-x2)
  }
  
  return(prob)
}


# 2) Likelihood function -----
loglikelihood.M11 <- function(pos, age, par, tau) {
  
  # seropos/seroneg by age -----
  tabela <- table(age,pos)
  
  # group m|n|t -----
  tabela.binom <- cbind(tabela[,2], rowSums(tabela), as.numeric(rownames(tabela)))
  
  likelihood <- apply(tabela.binom, 1, function(x,par,tau) dbinom(x=x[1], size=x[2], prob=rcm.M11.equation(t=x[3],par,tau),log=F), par=par,tau=tau)
  
  return(sum(log(likelihood)))
}


# 3) Profile MLE for the parameters -----
mle.M11.estimates <- function(pos, age, time.int, n.start, par) { 
  
  best.loglik <- (-1)*(10^6) # low loglikelihood initial value
  
  if(time.int[1] == time.int[2]) time.int <- time.int[1]
  
  output <- c()
  
  cat('\n1) Profile likelihood\n')
  
  for(tau in time.int) { # vary \tau -----
    
    # cat('tau=',tau,'\n',sep='')
    
    register <- c()
    
    for(i in 1:n.start) { # repeat n.start times each estimate of \tau -----
      
      # estimativas iniciais aleatorias -----
      par.ini <- par + runif(3,-0.5,0.5)
      
      # profile mle -----
      sol <- optim(par=par.ini, loglikelihood.M11, pos=pos, tau=tau, age=age, control=list(fnscale=(-1),maxit=1E+6))
      
      lambda <- exp(sol$par[1])
      rho1 <- exp(sol$par[2])
      p.rho2 <- exp(sol$par[3])
      rho2 <- rho1*p.rho2
      
      loglik.total <- sol$value
      
      register <- rbind(register, c(tau, loglik.total, lambda, rho1, rho2, sol$convergence))
      
      # select the best estimate for \tau -----
      aux <- which.max(register[,2]) # select the biggest loglikelihood value
      
      register <- register[aux,] # only that register
    }
    
    output <- rbind(output,register)
    
    if(register[2] > best.loglik) {
      
      output.vf <- register
      
      best.loglik <- register[2]	
    }
  }
  
  colnames(output) <- c('tau', 'loglik', 'lambda', 'rho1', 'rho2', 'convergence')
  
  return(output)
}


# 4) CI -----
rcm.M11.confidence.interval <- function(age, pos, par, tau, loglik.1){
  
  param.aux <- log(par)
  
  loglik1 <- loglikelihood.M11(age=age, pos=pos, par=param.aux, tau=tau)
  
  output <- c()
  
  for(j in 1:3) {
    if(j==1)cat('\n2) Confidence interval for lambda')
    if(j==2)cat('\n3) Confidence interval for rho1')
    if(j==3)cat('\n4) Confidence interval for rho2')
    
    
    param1 <- param.aux[j]
    param <- param.aux[-j]
    
    p <- qchisq(0.95, 1)
    
    loglik.0 <- loglik1 - p/2
    
    f2 <- function(param1, param, tau, age, pos, loglik.0, j){
        rcm.M11.estimates.one.par.fixed(age=age, pos=pos, param=param, param1=param1, tau=tau, j=j)-loglik.0
      }
     
    
    ##### CALCULATING LOWER BOUND #####
    
    if(param1 < (-10)) {
      
      lower.bound <- exp(param1)
      
    } else {
      
      sol <- rcm.M11.estimates.one.par.fixed(age=age, pos=pos, param=param, param1=-10, tau=tau, j=j)
      
      # cat('\nloglik(MLE) =',round(loglik1,2))
      # cat('\nCutoff(loglik) =',round(loglik.0,2))
      # cat('\nsol =',round(sol,2))
      
      if(sol > loglik.0){
        
        lower.bound <- exp(-100)
        
      } else {
        
        # print(param1)
        # print(param)
        # print(j)
        # print(rcm.M11.estimates.one.par.fixed(age=age, pos=pos, param=param, param1=param1, tau=tau, j=j))
        
        sol <- tryCatch(uniroot(f2,c(-10,param1), age=age, pos=pos, loglik.0=loglik.0, j=j, param=param, tau=tau),error=function(e){
          sol<-list(root=NA)
          return(sol)
        })
        
        lower.bound <- exp(sol$root)
      }
    }
    
    
    ##### CALCULATING UPPER BOUND #####
    
    sol <- rcm.M11.estimates.one.par.fixed(age=age, pos=pos, param=param, param1=20, tau=tau, j=j)
    
    if(sol > loglik.0){
      
      upper.bound <- exp(100)
      
    } else {
      
      sol <- tryCatch(uniroot(f2,c(param1,5), age=age, pos=pos, loglik.0=loglik.0, j=j, param=param, tau=tau), error=function(e){
        sol<-list(root=NA)
        return(sol)
      })
      
      upper.bound <- exp(sol$root)
    }
    
    cat('\nLower bound =',round(lower.bound,4))
    cat('\nUpper bound =',round(upper.bound,4),'\n')
    
    output <- rbind(output, c(lower.bound, upper.bound))
  }
  
  return(output)
}


# 5.1) Fix ONE parameter and estimate the remaining ones -----
rcm.M11.estimates.one.par.fixed <- function(age, pos, param, param1, tau, j) {
  
  if(j==1) fit <- optim(par=c(runif(1,-3,-1), runif(1,-3,-1)), fn = loglikelihood.M11.one.par.fixed, age=age, pos=pos, param1=param1, tau=tau, j=j, control=list(fnscale=-1,pgtol=1E-10))
  if(j==2) fit <- optim(par=c(runif(1,-3,-1), runif(1,-3,-1)), fn = loglikelihood.M11.one.par.fixed, age=age, pos=pos, param1=param1, tau=tau, j=j, control=list(fnscale=-1,pgtol=1E-10))
  if(j==3) fit <- optim(par=c(runif(1,-3,-1), log(1/runif(1,0,1))), fn = loglikelihood.M11.one.par.fixed, age=age, pos=pos, param1=param1, tau=tau, j=j, control=list(fnscale=-1,pgtol=1E-10))
  
  return(fit$value)
}

# 5.2) -----
loglikelihood.M11.one.par.fixed <- function(age, pos, param, param1, tau ,j) {
  
  if(j==1){
    all.param <- c(param1,param)
    sol <- loglikelihood.M11(pos=pos, age=age, par=all.param, tau=tau)
  }
  
  if(j==2){
    all.param <- c(param[1],param1,param[2])
    sol <- loglikelihood.M11(pos=pos, age=age, par=all.param, tau=tau)
  }
  
  if(j==3){
    all.param <- c(param,param1)
    sol <- loglikelihood.M11.Q(pos=pos, age=age, par=all.param, tau=tau)
  }
  
  return(sol)
}


# 5.3) -----
rcm.M11.equation.Q <- function(t, par, tau) {
  
  lambda <- exp(par[1])
  rho2 <- exp(par[3])
  q <- exp(par[2])
  rho1 <- rho2*q
  # p <- 1/q
  
  theta1 <- lambda/(lambda+rho1)
  theta2 <- lambda/(lambda+rho2)
  
  if(t<=tau){
    
    x1 <- (1-exp(-(lambda+rho1)*t))	
    prob <- theta1*x1
    
  } else {
    
    x1 <- (1-exp(-(lambda+rho1)*tau))
    x2 <- (1-exp(-(lambda+rho2)*(t-tau)))
    prob <- theta1*x1+theta2*x2*(1-x1)
  }
  
  return(prob)
}

# 5.4) -----
loglikelihood.M11.Q <- function(pos, age, par, tau) {
  
  # seropos/seroneg by age -----
  tabela <- table(age,pos)
  
  # group m|n|t -----
  tabela.binom <- cbind(tabela[,2], rowSums(tabela), as.numeric(rownames(tabela)))
  
  likelihood <- apply(tabela.binom, 1, function(x,par,tau) dbinom(x=x[1], size=x[2], prob=rcm.M11.equation.Q(t=x[3],par,tau),log=F), par=par,tau=tau)
  
  return(sum(log(likelihood)))
}


# 6) FINAL FUNCTION -----
M11.analysis <- function(age, pos, time.int=c(1:40), n.start=1, par=log(runif(3,0,1))){
  
  # Get Profile likelihood -----
  mle <- mle.M11.estimates(pos=pos, age=age, time.int=time.int, n.start=n.start, par=par)
  
  # Get parameter estimates -----
  aux <- which.max(mle[,2])
  
  estimates <- mle[aux,]
  
  print(estimates)
  
  tau <- as.integer(estimates[1])
  loglikelihood.total <- estimates[[2]]
  lambda <- estimates[[3]]
  rho1 <- estimates[[4]]
  rho2 <- estimates[[5]]
  convergence <- as.integer(estimates[6])
  p.rho2 <- rho2/rho1
  
  ##### confidence intervals #####
  
  conf.int <- rcm.M11.confidence.interval(age=age ,pos=pos, par=c(lambda,rho1,p.rho2), tau=tau, loglik.1=loglikelihood.total)
  
  parameters <- cbind(c(lambda,rho1,rho2),conf.int)
  colnames(parameters) <- c('estimates', 'lower', 'upper')
  rownames(parameters) <- c('lambda', 'rho1', 'rho2')
  
  fitted.values <- c(1:max(age))
  
  exp.seroprev <- sapply(1:max(age),rcm.M11.equation, par=c(log(lambda),log(rho1), log(rho2)), tau=tau)
  
  fitted.values<-cbind(fitted.values,exp.seroprev)
  
  # print(conf.int)				
  
  output <- list(loglikelihood.total=loglikelihood.total, estimates=parameters, tau=tau, df=4, expected.seroprevalence=fitted.values, proflik=mle, model='M1,1')
  
  return(output)
}
