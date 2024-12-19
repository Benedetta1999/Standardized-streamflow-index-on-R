pacman::p_load(readr,SPEI,openxlsx,tidyverse,fitdistrplus,actuar,evd,dplyr,
               lubridate,gridExtra,FAdist,SCI,bbmle,lmom)

#File input
{
  #Read .csv files of daily discharge in working directory
  file_list <- list.files(getwd(), pattern = "\\.csv$", full.names = TRUE)
  
  #Create dataframe from .csv file
  my_df <- read_csv(file_list[13])
  
  names(my_df) <- c("Date", "Discharge")
  
  rm(file_list)
}

#Monthly mean
{
  my_df$Year_month <- paste(year(my_df$Date), month(my_df$Date), sep = "-")
  
  months_unique <- unique(my_df$Year_month)
  
  #Check if enough values are present in each month - remove month if >20% values are missing
  for(i in length(months_unique)){
    if(sum(!is.na(my_df$Discharge[my_df$Year_month == months_unique[1]]))<
       0.8*days_in_month(ym(months_unique[1]))){
      my_df <- my_df[-(my_df$Year_month == months_unique[1]),] #Removes month
    }
  }
  rm(months_unique)
  
  #Monthly mean
  discharge <- aggregate(my_df$Discharge, by = list(ym(my_df$Year_month)), mean)
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
  
  #Formula to calculate Pearson Type III
  fitdist_pe3 <- function(x){tryCatch({
    lmom_disch <- as.vector(pelpe3(x))
    lmom_disch <- list(shape = lmom_disch[3],
                       scale = lmom_disch[2],
                       location = lmom_disch[1])
    
    fit <- suppressWarnings(fitdist(month_disch, "pe3",start = lmom_disch))
    return(fit)
  },
  error = function(e) {
    fit <- list(aic = Inf)
    return(fit)
  })}
  
  #Formula to calculate GEV
  fitdist_gev <- function(x){tryCatch({
    lmom_disch <- list(shape = 1,
                       scale = 1,
                       location = 0)
    fit <- fitdist(month_disch, "gev",
                       start = lmom_disch)
    return(fit)
  },
  error = function(e) {
    fit <- list(aic = Inf)
    return(fit)
  })}
  
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
    
    distr$fit_pears3[[i]] <- fitdist_pe3(month_disch)
    
    distr$fit_gev[[i]] <- fitdist_gev(month_disch)
    
    distr$fit_weib[[i]] <- tryCatch({fitdist(month_disch, "weibull")},
                              error = function(e){return(list(aic=Inf))})

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

#Plot chosen distributions
{
  par(mfrow=c(3,8))
  for(i in 1:12){
    cdfcomp(best_fit[[i]])
    qqcomp(best_fit[[i]])
  }
}