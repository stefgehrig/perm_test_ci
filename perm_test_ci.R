perm_test_ci <- function(outcome,            # name of numeric vector that encodes the outcome which is compared 
                         group,              # name of numeric, factor, or character vector that encodes the two groups
                         data = NULL,        # data frame containing the variables, either tibble data.frame (mandatory)
                         conf.int = TRUE,    # whether confidence intervals for the difference in means should be computed
                         conf.level = 0.95,  # desired level of confidence
                         perms = 5000){      # number of permutations for testing (>1000 is mandatory)
                                             # search steps for confidence limits are perms * 1.5
  
  if(perms <= 1e3)  stop("Number of permutations must be larger than 1000")
  if(is.null(data)) stop("Please supply outcome and group vector in a data.frame or tibble")

  # remove rows with NA
  data <- as.data.frame(data)
  data <- data[!is.na(data[,outcome]) & !is.na(data[,group]),]
  
  # create shortcuts
  out  <- as.numeric(as.character(data[,outcome]))
  gro  <- as.factor(data[,group])
  gro1 <- levels(gro)[1]
  gro2 <- levels(gro)[2]
  
    # calculate observed empirical difference
    emp.diff  <- mean(out[gro==gro1])-mean(out[gro==gro2])
     
    # simulate differences via permutation
    sim.diff  <- rep(NA, perms)
    for(i in seq(1:perms)){                            
       shuffled_gro       <- sample(gro, nrow(data), replace=F)
       sim.diff[i]        <- mean(out[shuffled_gro==gro1])-mean(out[shuffled_gro==gro2])
    }
     
    if (conf.int == FALSE){
     
       # plot distribution and observed value
       hist(sim.diff, main = "Histogram of\n simulated differences", xlab = "Simulated mean differences",
            xlim=c(min(c(emp.diff,sim.diff)), max(c(emp.diff,sim.diff))))
       abline(v=emp.diff, lwd = 2, col = "red")
       text(emp.diff, 0, "observed difference", col = "red")
       
       # calculate two-sided p-value
       p <- sum(abs(sim.diff) >= abs(emp.diff))/perms
       
       return(cat(paste("Mean difference (", gro1, " - ", gro2, ") = ", round(emp.diff,4), "\n",
                    "P-value (", perms, " permutations) = ", round(p,6), sep="")))}
     
    else{
       
      # plot distribution and observed value
      par(mfrow=c(1,3))
      hist(sim.diff, main = "Histogram of\n simulated differences", xlab = "Simulated mean differences",
           xlim=c(min(c(emp.diff,sim.diff)), max(c(emp.diff,sim.diff))))
      abline(v=emp.diff, lwd = 2, col = "red")
      text(emp.diff, 0, "Observed difference", col = "red")
      
      # calculate two-sided p-value
      p <- sum(abs(sim.diff) >= abs(emp.diff))/perms
      
      # calcultate confidence intervals according to Garthwaite 1996
      # define parameters for search
      a <- (1-conf.level)/2
      z <- qnorm(a, lower.tail = FALSE)
      k <- 2/(z * (2*pi)^(-1/2) * exp((-z^2)/2))
      m <- round(0.3*(2-a) / (a),0) 
      n <- round(perms * 2.5,0) # number of search iterations should be larger than number of permutations for testing
      
      # create shifted group mean such that means are equal
      shifted_mean <- emp.diff + mean(out[gro==gro2])
      
      # define starting points for upper and lower search for CI endpoints
      sim.diff<-NULL
      shuffled_gro <- NULL
      for (i in seq(1,(2-a)/a,1)){
        shuffled_gro       <- sample(gro, nrow(data), replace=F)
        sim.diff[i]        <- mean(out[shuffled_gro==gro1])-shifted_mean}
        
      startsearch.lower <- emp.diff+(min(sim.diff[sim.diff!=min(sim.diff)]) - 
                             max(sim.diff[sim.diff!=max(sim.diff)]))/2
      startsearch.upper <- emp.diff-(min(sim.diff[sim.diff!=min(sim.diff)]) - 
                             max(sim.diff[sim.diff!=max(sim.diff)]))/2 
        
      # initiate vectors for search results
      U <- NULL
      L <- NULL
      
      # search lower CI
      L[1] <- startsearch.lower
      for(i in seq(1,n-1,1)){
        
        emp_delta    <- L[i] + out[gro==gro1]
        emp_selmean  <- mean(sample(c(emp_delta, out[gro==gro2]),
                                    length(out[gro==gro1]), replace=F))
        
        c <- k*(emp.diff - L[i])
        
        if (emp_selmean < mean(out[gro==gro1])){
          
          L[i+1] <- L[i] + c*a/(i+m-1)}
        
        else{
          
          L[i+1] <- L[i] - c*(1-a)/(i+m-1)}
      }
      
      # save mean of of last 100 iterations and plot convergence diagnostic
      result.lower <- round(mean(tail(L,1e2)),4)
      plot(L, type="l", ylab="Confidence limit", xlab = "Search step", 
           main = "Search diagnostics lower limit\n(line shows final value\n as mean of last 10 iterations)")
      abline(h = result.lower)
        
      # search upper CI
      U[1] <- startsearch.upper
      for(i in seq(1,n-1,1)){
        
        emp_delta    <- U[i] + out[gro==gro1]
        emp_selmean  <- mean(sample(c(emp_delta, out[gro==gro2]),
                                    length(out[gro==gro1]), replace=F))
        
        c <- k*(U[i] - emp.diff)
        
        if (emp_selmean > mean(out[gro==gro1])){
          
          U[i+1] <- U[i] - c*a/(i+m-1)}
        
        else{
          
          U[i+1] <- U[i] + c*(1-a)/(i+m-1)}
      }
      
      # save mean of of last 100 iterations and plot convergence diagnostic
      result.upper <- round(mean(tail(U,1e2)),4)
      plot(U, type="l", ylab="Confidence limit", xlab = "Search step", 
           main = "Search diagnostics upper limit\n(line shows final value\n as mean of last 100 iterations)")
      abline(h = result.upper)
      
      return(cat(paste("Observed mean difference (", gro1, " - ", gro2, ") = ", round(emp.diff,4), "\n",
                       "p-value (", perms, " permutations) = ", round(p,8), "\n",
                       "Confidence interval ", "(", conf.level*100, "%): ", 
                       result.lower, "; ", result.upper,
                       sep="")))
    }
}
    
#################
#### Example ####
#################
# set.seed(0)
# 
# df <- data.frame(
#   treatment_arm   = c( rep("control", 250), rep("intervention", 250) ),
#   outcome_measure = c( rnorm(250, 100, 10), rnorm (250, 102, 10) )
# )
# 
# perm_test_ci(
#   outcome    = "outcome_measure",
#   group      = "treatment_arm",
#   data       = df,
#   conf.int   = TRUE,
#   conf.level = 0.95,
#   perms      = 1e4
# )
