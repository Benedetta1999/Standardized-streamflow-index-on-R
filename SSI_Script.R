pacman::p_load(readr,SPEI,tidyverse,fitdistrplus,actuar,dplyr,future.apply, future,
               lubridate,gridExtra,FAdist,SCI,bbmle,lmom, PearsonDS, logspline)

#File input
{
  #Read .csv files of daily discharge in working directory
  file_list <- list.files(getwd(), pattern = "\\.csv$", full.names = TRUE)
}

#List of functions for fitting and K-S test calculation for Monte Carlo  method
{
  ks_test_mc <- list(
    #Pareto
    function(x){
      tryCatch({
        r <- rpareto(length(x$data), shape=x$estimate[1], scale=x$estimate[2])
        fit <- fitdist(r, "pareto")
        as.numeric(ks.test(r,"ppareto",
                           shape= fit$estimate["shape"],
                           scale = fit$estimate["scale"]
        )$statistic)
      }, error = function(e){
        NaN
      })},
    
    #Log-logistic
    function(x){
      tryCatch({
        r <- rllogis(length(x$data), shape=x$estimate[1], scale=x$estimate[2])
        fit <- fitdist(r, "llogis")
        as.numeric(ks.test(r,"pllogis",
                           shape= fit$estimate["shape"],
                           scale = fit$estimate["scale"]
        )$statistic)
      }, error = function(e){
        NaN
      })},
    
    #Log-normal
    function(x){
      tryCatch({
        r <- rlnorm(length(x$data), meanlog=x$estimate[1], sdlog=x$estimate[2])
        fit <- fitdist(r, "lnorm")
        as.numeric(ks.test(r,"plnorm",
                           meanlog = fit$estimate["meanlog"],
                           sdlog = fit$estimate["sdlog"]
        )$statistic)
      }, error = function(e){
        NaN
      })},
    
    #Pearson Type III
    function(x){
      tryCatch({
        r <- rgamma3(length(x$data), shape=x$estimate[1],
                     scale=x$estimate[2],
                     thres=x$estimate[3])
        fit <- fitdist(r, "gamma3",
                       start=list(
                         shape=mean(x$data)^2/var(x$data),
                         scale=var(x$data)/mean(x$data),
                         thres=0))
        
        as.numeric(ks.test(r,"pgamma3",
                           shape = fit$estimate["shape"],
                           scale = fit$estimate["scale"],
                           thres = fit$estimate["thres"]
        )$statistic)
      }, error = function(e){
        NaN
      })},
    
    #GEV
    function(x){
      tryCatch({
        r <- replicate(length(x$data),       #replicate is used due to an error in the rgev function
                       rgev(length(x$data), shape=x$estimate[1],
                            scale=x$estimate[2],
                            location=x$estimate[3]))
        fit <- fitdist(r, "gev",
                       start=list(shape = 1,
                                  scale = 1,
                                  location = 0))
        
        as.numeric(ks.test(r,"pgev",
                           shape = fit$estimate["shape"],
                           scale = fit$estimate["scale"],
                           location = fit$estimate["location"]
        )$statistic)
      }, error = function(e){
        NaN
      })},
    
    #Weibull
    function(x){
      tryCatch({
        r <- rweibull(length(x$data), shape=x$estimate[1], scale=x$estimate[2])
        fit <- fitdist(r, "weibull")
        as.numeric(ks.test(r,"pweibull",
                           shape= fit$estimate["shape"],
                           scale = fit$estimate["scale"]
        )$statistic)
      }, error = function(e){
        NaN
      })}
  ) 


#List of functions to obtain p-value from a vector of simulated K-S values
  p_val_ks_mc <- list(
    #Pareto
    function(x, spline_fit){
      1 - plogspline(ks.test(x$data,
                             "ppareto",
                             shape= x$estimate["shape"],
                             scale = x$estimate["scale"])$statistic,
                     spline_fit)
    },
    
    #Log-logistic
    function(x, spline_fit){
      1 - plogspline(ks.test(x$data,
                             "pllogis",
                             shape= x$estimate["shape"],
                             scale = x$estimate["scale"])$statistic,
                     spline_fit)
    },
    
    #Log-normal
    function(x,spline_fit){
      1 - plogspline(ks.test(x$data,
                             "plnorm",
                             meanlog= x$estimate["meanlog"],
                             sdlog = x$estimate["sdlog"])$statistic,
                     spline_fit)
    },
    
    #Pearson Type III
    function(x, spline_fit){
      1-plogspline(ks.test(x$data,
                           "pgamma3",
                           shape = x$estimate["shape"],
                           scale = x$estimate["scale"],
                           thres = x$estimate["thres"])$statistic,
                   spline_fit)
    },
    
    #GEV
    function(x, spline_fit){
      1-plogspline(ks.test(x$data,
                           "pgev",
                           shape = x$estimate["shape"],
                           scale = x$estimate["scale"],
                           location = x$estimate["location"])$statistic,
                   spline_fit)
    },
    
    #Weibull
    function(x,spline_fit){
      1 - plogspline(ks.test(x$data,
                             "pweibull",
                             shape= x$estimate["shape"],
                             scale = x$estimate["scale"])$statistic,
                     spline_fit)
    }
  )
}

#Repeat for all discharge files in working directory
{
  for(m in 1:length(file_list)){
    #Create dataframe from .csv file
    my_df <- read_csv(file_list[m])
    
    names(my_df) <- c("Date", "Discharge")
      
    #Name of location
    loc_name <- gsub(".csv","",gsub("^.*/","",file_list[[m]]))
    
    #Monthly mean
    {
      my_df$Year_month <- paste(year(my_df$Date), month(my_df$Date), sep = "-")
      
      months_unique <- unique(my_df$Year_month)
      
      #Check if enough values are present in each month - remove month if >20% values are missing
      for(i in length(months_unique)){
        if(sum(!is.na(my_df$Discharge[my_df$Year_month == months_unique[i]]))<
           0.8*days_in_month(ym(months_unique[i]))){
          my_df[my_df$Year_month == months_unique[i],] <- NA #Flag as NA
        }
      }
      rm(months_unique)
      
      #Monthly mean
      discharge <- aggregate(my_df$Discharge, by = list(ym(my_df$Year_month)),
                             FUN = function(x) mean(x, na.rm = TRUE))
      names(discharge) <- c("Date", "Discharge")
    }
    
    #Plot monthly discharge
    {
      ggplot(discharge, aes(x = Date, y = Discharge)) +
        geom_line(group = 1, color = "black") +
        labs(title = "Mean monthly discharge",
             x = "",
             y = "Mean discharge (m³/s)") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_x_date(date_breaks = "1 years", date_labels = "%Y")
    }
    
    #Fit distributions to each month of the year
    {
      #Create list for the fitting of the distributions
      distr <- list(fit_pareto = list(),
                    fit_ll = list(),
                    fit_ln = list(),
                    fit_pears3 = list(),
                    fit_gev = list(),
                    fit_weib = list())
    
      for(i in 1:12){
        #Extract values of discharge for a given month of the year
        month_disch <- discharge$Discharge[month(discharge$Date)==i]
        
        #Fitting
        distr$fit_pareto[[i]] <- tryCatch({fitdist(month_disch, "pareto")},
                                    error = function(e){return(list(aic=Inf))})
        
        distr$fit_ll[[i]] <- tryCatch({fitdist(month_disch, "llogis")},
                                error = function(e){return(list(aic=Inf))})
        
        distr$fit_ln[[i]] <- tryCatch({fitdist(month_disch, "lnorm")},
                                error = function(e){return(list(aic=Inf))})
        
        distr$fit_pears3[[i]] <- tryCatch({fitdist(month_disch, "gamma3",
                                                   start=list(
                                                     shape=mean(month_disch)^2/var(month_disch),
                                                     scale=var(month_disch)/mean(month_disch),
                                                     thres=0))},
                                          error = function(e){return(list(aic=Inf))})
        
        distr$fit_gev[[i]] <- tryCatch({fitdist(month_disch, "gev",
                                      start = list(
                                        shape=1,
                                        scale=1,
                                        location=0))},
                                      error = function(e){return(list(aic=Inf))})
        
        distr$fit_weib[[i]] <- tryCatch({fitdist(month_disch, "weibull")},
                                  error = function(e){return(list(aic=Inf))})
    
      }
      
      rm(month_disch)
    }
    
    #Test if distributions are to be rejected through K-S test with simulated p-value
    {
      n_sims <- 1e2
      
      #Multithread calculation
      plan(multisession, workers = 6)
      
    
      #Simulate p-value through MC method
      ks_p <- array(data = NA, dim = c(length(distr), length(distr[[1]])))
      for(i in 1:length(distr)){
        for(j in 1:12){
          if(is.finite(distr[[i]][[j]]$aic)){ #Checks if distribution could be evaluated or skips
            stat <- future_replicate(n_sims, ks_test_mc[[i]](distr[[i]][[j]]),
                                     future.packages = c("actuar",
                                                         "fitdistrplus",
                                                         "FAdist"))
            fit <- logspline(stat)
            
            ks_p[i,j] <- p_val_ks_mc[[i]](distr[[i]][[j]],fit)
            
            #The code can be slow, so indices are printed
            print(c(i,j))
          }
        }
      }
      
      #Remove distributions that do not pass the K-S test
      for(i in 1:6){
        for(j in 1:12){
          if(is.finite(ks_p[i,j]) && ks_p[i,j]<0.05){
            distr[[i]][[j]] <- list(aic=Inf)
            print(c(i,j))
          }
        }
      }
    }
    
    #Find best distribution for each month fo the year using AIC
    {
      best_fit <- list()
      for(i in 1:12){
        aic_best <- Inf
        for(j in 1:length(distr)){
          if(distr[[j]][[i]]$aic<aic_best){
            aic_best <- distr[[j]][[i]]$aic
            best_fit[[i]] <- distr[[j]][[i]]
          }
        }
        
      }
    }
    
    #Case where a month has no available distributions
    {
      for(i in 1:12){
        if(is.null(best_fit[[i]])){
          #A month without available distribution (excluded by K-S test) is given one
          #even if wouldn't pass the test
          
          #Distribution fitting
          {
            #Extract values of discharge for a given month of the year
            month_disch <- discharge$Discharge[month(discharge$Date)==i]
            
            #Fitting
            distr$fit_pareto[[i]] <- tryCatch({fitdist(month_disch, "pareto")},
                                              error = function(e){return(list(aic=Inf))})
            
            distr$fit_ll[[i]] <- tryCatch({fitdist(month_disch, "llogis")},
                                          error = function(e){return(list(aic=Inf))})
            
            distr$fit_ln[[i]] <- tryCatch({fitdist(month_disch, "lnorm")},
                                          error = function(e){return(list(aic=Inf))})
            
            distr$fit_pears3[[i]] <- tryCatch({fitdist(month_disch, "gamma3",
                                                       start=list(
                                                         shape=mean(month_disch)^2/var(month_disch),
                                                         scale=var(month_disch)/mean(month_disch),
                                                         thres=0))},
                                              error = function(e){return(list(aic=Inf))})
            
            distr$fit_gev[[i]] <- tryCatch({fitdist(month_disch, "gev",
                                                    start = list(
                                                      shape=1,
                                                      scale=1,
                                                      location=0))},
                                           error = function(e){return(list(aic=Inf))})
            
            distr$fit_weib[[i]] <- tryCatch({fitdist(month_disch, "weibull")},
                                            error = function(e){return(list(aic=Inf))})
          }
          
          #Choice of best distribution
          aic_best <- Inf
          for(j in 1:length(distr)){
            if(distr[[j]][[i]]$aic<aic_best){
              aic_best <- distr[[j]][[i]]$aic
              best_fit[[i]] <- distr[[j]][[i]]
            }
          }
          best_fit[[j]]
        }
      }
    }
    
    #Plot chosen distributions
    {
      png(filename = paste(loc_name,
                           "png",sep = "."),
          width = 3000, height = 2000)
      par(mfrow=c(3,8))
      for(i in 1:12){
        cdfcomp(best_fit[[i]])
        qqcomp(best_fit[[i]])
      }
      
      dev.off()
    }
    
    #Calculate the SSI using the best distributions
    {
      ssi <- data.frame(Date = discharge$Date,
                        SSI = rep(NA, length(discharge$Date)))
      for(i in 1:12){
        #Calculate quantiles
        switch(best_fit[[i]]$distname, 
               pareto={
                 p <- ppareto(discharge$Discharge[month(discharge$Date)==i],
                              shape=best_fit[[i]]$estimate[1],
                              scale=best_fit[[i]]$estimate[2])
               },
               llogis={
                 p <- pllogis(discharge$Discharge[month(discharge$Date)==i],
                              shape=best_fit[[i]]$estimate[1],
                              scale=best_fit[[i]]$estimate[2])    
               },
               lnorm={
                 p <- plnorm(discharge$Discharge[month(discharge$Date)==i],
                             meanlog=best_fit[[i]]$estimate[1],
                             sdlog=best_fit[[i]]$estimate[2])
               },
               gamma3={
                 p <- pgamma3(discharge$Discharge[month(discharge$Date)==i],
                              shape=best_fit[[i]]$estimate[1],
                              scale=best_fit[[i]]$estimate[2],
                              thres=best_fit[[i]]$estimate[3])
               },
               gev={
                 p <- pgev(discharge$Discharge[month(discharge$Date)==i],
                      shape=best_fit[[i]]$estimate[1],
                      scale=best_fit[[i]]$estimate[2],
                      location=best_fit[[i]]$estimate[3])
               },
               weibull={
                 p <- pweibull(discharge$Discharge[month(discharge$Date)==i],
                              shape=best_fit[[i]]$estimate[1],
                              scale=best_fit[[i]]$estimate[2])
               }
        )
        
        #Normalization
        ssi$SSI[month(ssi$Date)==i] <- qnorm(p, 0, 1)
      }  
    }
    
    #Plot SSI values
    {
      ggplot(ssi, aes(x = Date, y = SSI, fill = SSI)) +
        geom_col(group = 1, show.legend = FALSE) +
        scale_fill_stepsn(breaks = c(0),
                          colours = c("red", "blue"))+
        labs(title = paste("Standardized Streamflow Index -", loc_name),
             x = "",
             y = "SSI") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_x_date(date_breaks = "1 years", date_labels = "%Y")
    }
    
    #Save SSI data
    save(ssi,file = paste(loc_name,
                          "RData",sep = "."))
    
  } #Closes the loop through all files in the working directory
}
