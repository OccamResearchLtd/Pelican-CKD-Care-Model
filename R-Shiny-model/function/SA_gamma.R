# The function for parameters following the gamma distribution
## to generate the lower and upper bounds for one-way sensitivity analysis
## and to generate the random value for probabilistic sensitivity analysis
SA_gamma <- function(vec.mean_se = NULL, 
                     vec.mean_l_u = NULL,
                     num.mean = NULL,
                     num.se = NULL,
                     num.lower = NULL,
                     num.upper = NULL,
                     n = 12,
                     tot_n = 12) {
  case <- NA
  
  if(is.vector(vec.mean_se) & length(vec.mean_se) == 2 & 
     is.null(vec.mean_l_u) & is.null(num.mean) & is.null(num.se) & 
     is.null(num.lower) & is.null(num.upper)) {
    case <- 1 
    # a vector with 2 elements for Mean and SE
    
    
  } else if(is.vector(vec.mean_l_u) & length(vec.mean_l_u) == 3 & 
            is.null(vec.mean_se) & is.null(num.mean) & is.null(num.se) & 
            is.null(num.lower) & is.null(num.upper)) {
    case <- 2 
    # a vector with 3 elements for Mean, Lower bound and Upper bound
    
  } else if(is.numeric(num.mean) & is.numeric(num.se) & 
            is.null(vec.mean_se) & is.null(vec.mean_l_u) & 
            is.null(num.lower) & is.null(num.upper)) {
    case <- 3
    # 2 separate numbers for Mean and SE
    
  } else if(is.numeric(num.mean) & 
            is.numeric(num.lower) & is.numeric(num.upper) &
            is.null(vec.mean_se) & is.null(vec.mean_l_u) & 
            is.null(num.se)) {
    case <- 4
    # 3 separate numbers for Mean, Lower bound and Upper bound
  }
  
  if(is.na(case)) {stop("Invalid input.")}
  
  mean <- switch(case, 
                 as.double(vec.mean_se[1]),
                 as.double(vec.mean_l_u[1]),
                 as.double(num.mean),
                 as.double(num.mean))
  
  se <- switch(case, 
               as.double(vec.mean_se[2]),
               NULL,
               as.double(num.se),
               NULL)
  
  lb <- switch(case, 
               NULL,
               as.double(vec.mean_l_u[2]),
               NULL,
               as.double(num.lower))
  
  ub <- switch(case, 
               NULL,
               as.double(vec.mean_l_u[3]),
               NULL,
               as.double(num.upper))
  
  
  output <- rep(NA,tot_n*2+2)
  names(output) <- c("base", "prob", rep(c("lower", "upper"),tot_n))
  output[1] <- mean
  
  # calculate the two parameters for gamma distribution
  se <- ifelse(is.null(se), (ub - lb) / (2 * qnorm(0.975)), se)
  alpha <- (mean ^ 2) / (se ^ 2)
  beta <- (se ^ 2) / mean
  
  output[2*n+1] <- ifelse(is.null(lb), qgamma(0.025, shape = alpha, scale = beta), lb)
  output[2*n+2] <- ifelse(is.null(ub), qgamma(0.975, shape = alpha, scale = beta), ub)
  
  # draw a random number from 0 to 1
  rand <- runif(1)
  # use "while" loop to avoid rand == 1, which is theoretically possible with 0 probability
  while(rand == 1) {
    rand <- runif(1)
  }
  output[2] <- qgamma(rand, shape = alpha, scale = beta)
  
  names(output)[is.na(output)] <- names(output)[1]
  output[is.na(output)] <- output[1]
  
  return(output)
  
}