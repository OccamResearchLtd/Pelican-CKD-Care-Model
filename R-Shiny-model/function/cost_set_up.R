## Function to set up the costs used in the core model for G3a to G5 ----
cost_set_up_g3a_g5 <- function(data,
                               cat1, 
                               cat2, 
                               cat3,
                               event) {
  ### No diabetes
  c_mat_temp1 <- array(dim = c(length(cat1), length(cat2), length(cat3)),
                       dimnames = list(eGFR = cat1,
                                       Age = cat2,
                                       Year = cat3))
  
  noevent_nodm_df <- subset(data, State == event & 
                              Diabetes == "no baseline t2d")
  #### Iterate over the array
  i_s <- i_a <- i_f <- 0
  for (s in levels(noevent_nodm_df$Stage)) {
    i_s <- i_s + 1
    i_a <- 0
    for (a in levels(noevent_nodm_df$Age)) {
      i_a <- i_a + 1
      i_f <- 0
      for (f in levels(noevent_nodm_df$FirstYear)) {
        i_f <- i_f + 1
        # Find the corresponding Live value in the dataframe
        live_value <- noevent_nodm_df[noevent_nodm_df$Stage == s & 
                                        noevent_nodm_df$Age == a & 
                                        noevent_nodm_df$FirstYear == f, 
                                      "Live"]
        # If a Live value is found, replace the NA with it
        if (length(live_value) > 0) {
          c_mat_temp1[i_s, i_a, i_f] <- live_value
        }
      }
    }
  }
  
  
  ### With diabetes
  c_mat_temp2 <- array(dim = c(length(cat1), length(cat2), length(cat3)),
                       dimnames = list(eGFR = cat1,
                                       Age = cat2,
                                       Year = cat3))
  
  noevent_dm_df <- subset(data, State == event & 
                            Diabetes == "baseline t2d")
  
  #### Iterate over the array
  i_s <- i_a <- i_f <- 0
  for (s in levels(noevent_dm_df$Stage)) {
    i_s <- i_s + 1
    i_a <- 0
    for (a in levels(noevent_dm_df$Age)) {
      i_a <- i_a + 1
      i_f <- 0
      for (f in levels(noevent_dm_df$FirstYear)) {
        i_f <- i_f + 1
        # Find the corresponding Live value in the dataframe
        live_value <- noevent_dm_df[noevent_dm_df$Stage == s & 
                                      noevent_dm_df$Age == a & 
                                      noevent_dm_df$FirstYear == f, 
                                    "Live"]
        # If a Live value is found, replace the NA with it
        if (length(live_value) > 0) {
          c_mat_temp2[i_s, i_a, i_f] <- live_value
        }
      }
    }
  }
  
  ### average over diabetes
  return(list(c_mat_temp1, c_mat_temp2))
}
