#' Interpolation for iAR, CiAR, and BiAR Classes
#'
#' Performs interpolation on time series with missing values. This method is implemented for:
#' 1. Irregular Autoregressive models (iAR)
#' 2. Complex Irregular Autoregressive models (CiAR)
#' 3. Bivariate Autoregressive models (BiAR)
#' 
#' @name interpolation
#' 
#' @param x An object of class \code{iAR}, \code{CiAR}, or \code{BiAR}, containing the model specification and parameters:
#'   \itemize{
#'     \item For \code{iAR}:
#'       \itemize{
#'         \item \code{family}: The distribution family of the iAR model (one of "norm", "t", or "gamma").
#'         \item \code{series}: A numeric vector representing the time series to interpolate.
#'         \item \code{coef}: The coefficients of the iAR model.
#'         \item \code{zero_mean}: Logical, whether the model assumes a zero-mean series.
#'         \item \code{standardized}: Logical, whether the model uses standardized data (only for "norm" family).
#'         \item \code{mean}: The mean parameter (only for "gamma" family).
#'       }
#'     \item For \code{CiAR}:
#'       \itemize{
#'         \item \code{coef}: The coefficients of the CiAR model.
#'         \item \code{series_esd}: The series of error standard deviations for the CiAR model.
#'         \item \code{zero_mean}: Logical, whether the model assumes a zero-mean series.
#'         \item \code{standardized}: Logical, whether the model uses standardized data.
#'         \item \code{seed}: Optional seed for random number generation.
#'       }
#'     \item For \code{BiAR}:
#'       \itemize{
#'         \item \code{coef}: The coefficients of the BiAR model.
#'         \item \code{series_esd}: The series of error standard deviations for the BiAR model.
#'         \item \code{zero_mean}: Logical, whether the model assumes a zero-mean series.
#'         \item \code{yini1}: Initial value for the first time series (for BiAR models).
#'         \item \code{yini2}: Initial value for the second time series (for BiAR models).
#'         \item \code{seed}: Optional seed for random number generation.
#'       }
#'   }
#'
#' @param ... Additional arguments (unused).
#'
#' @return An object of the same class as \code{x} with interpolated time series.
#'
#' @description
#' This method performs imputation of missing values in a time series using an autoregressive model.
#' The imputation is done iteratively for each missing value, utilizing available data and model coefficients.
#' Depending on the model family, the interpolation is performed differently:
#' - For \code{norm}: A standard autoregressive model for normally distributed data.
#' - For \code{t}: A model for time series with t-distributed errors.
#' - For \code{gamma}: A model for time series with gamma-distributed errors.
#' - For \code{CiAR}: A complex irregular autoregressive model.
#' - For \code{BiAR}: A bivariate autoregressive model.
#'
#' @details
#' The method handles missing values (NA) in the time series by imputing them iteratively. 
#' For each missing value, the available data is used to fit the autoregressive model, and the missing value is imputed based on the model's output.
#' For the \code{CiAR} and \code{BiAR} models, the error standard deviations and initial values are also considered during imputation.
#'
#' @seealso \code{\link{forecast}} for forecasting methods for these models.
#'
#' @references
#' \insertRef{Eyheramendy_2018}{iAR},\insertRef{Elorrieta_2019}{iAR},\insertRef{Elorrieta_2021}{iAR}
#'
#' @examples
#' # Interpolation for iAR model
#' library(iAR)
#' n=100
#' set.seed(6714)
#' o=iAR::utilities()
#' o<-gentime(o, n=n)
#' times=o@times
#' model_norm <- iAR(family = "norm", times = times, coef = 0.9)  
#' model_norm <- sim(model_norm)
#' y=model_norm@series
#' y1=y/sd(y)
#' model_norm@series=y1
#' model_norm@series_esd=rep(0,n)
#' model_norm <- kalman(model_norm) 
#' print(model_norm@coef)
#' napos=10
#' model_norm@series[napos]=NA
#' model_norm <- interpolation(model_norm)
#' interpolation=na.omit(model_norm@interpolated_values)
#' mse=as.numeric(y1[napos]-interpolation)^2
#' print(mse)
#' plot(times,y,type='l',xlim=c(times[napos-5],times[napos+5]))
#' points(times,y,pch=20)
#' points(times[napos],interpolation*sd(y),col="red",pch=20)
#' 
#' # Interpolation for CiAR model
#' model_CiAR <- CiAR(times = times,coef = c(0.9, 0))
#' model_CiAR <- sim(model_CiAR)
#' y=model_CiAR@series
#' y1=y/sd(y)
#' model_CiAR@series=y1
#' model_CiAR@series_esd=rep(0,n)
#' model_CiAR <- kalman(model_CiAR)
#' print(model_CiAR@coef)
#' napos=10
#' model_CiAR@series[napos]=NA
#' model_CiAR <- interpolation(model_CiAR)
#' interpolation=na.omit(model_CiAR@interpolated_values)
#' mse=as.numeric(y1[napos]-interpolation)^2
#' print(mse)
#' plot(times,y,type='l',xlim=c(times[napos-5],times[napos+5]))
#' points(times,y,pch=20)
#' points(times[napos],interpolation*sd(y),col="red",pch=20)

#' # Interpolation for BiAR model
#' model_BiAR <- BiAR(times = times,coef = c(0.9, 0.3), rho = 0.9)
#' model_BiAR <- sim(model_BiAR)
#' y=model_BiAR@series
#' y1=y/apply(y,2,sd)
#' model_BiAR@series=y1
#' model_BiAR@series_esd=matrix(0,n,2)
#' model_BiAR <- kalman(model_BiAR)
#' print(model_BiAR@coef) 
#' napos=10
#' model_BiAR@series[napos,1]=NA
#' model_BiAR@series_esd[napos,1]=NA
#' model_BiAR <- interpolation(model_BiAR)
#' interpolation=na.omit(model_BiAR@interpolated_values[,1])
#' mse=as.numeric(y1[napos,1]-interpolation)^2
#' print(mse)
#' par(mfrow=c(2,1))
#' plot(times,y[,1],type='l',xlim=c(times[napos-5],times[napos+5]))
#' points(times,y[,1],pch=20)
#' points(times[napos],interpolation*apply(y,1,sd)[1],col="red",pch=20)
#' plot(times,y[,2],type='l',xlim=c(times[napos-5],times[napos+5]))
#' points(times,y[,2],pch=20)
#' @export
interpolation <- S7::new_generic("interpolation", "x")
S7::method(generic = interpolation, signature = iAR) <- function(x, yini = 0) {
  if(length(x@series)==0) stop("The interpolation method needs a time series")
  if(x@family == "norm"){
    if(length(x@coef)==0) stop("The interpolation method needs the coefficient of the iAR model")
    # Original series
    vector <- x@series
    
    # Find positions of NA values in the vector
    nas <- which(is.na(vector))
    
    # If there are NA values, return the imputed vector
    if (length(nas) > 0) {
      for (i in seq_along(nas)) {
        # Current position of NA value
        na_pos <- nas[i]
        
        # Data available up to before the current NA value
        if((i+1) <= length(nas)){
          available_data <- vector[-nas[((i+1):length(nas))]]
          available_times <- x@times[-nas[((i+1):length(nas))]]
        }  else {
          available_data <- vector
          available_times <- x@times
        }
        if(is.integer(x@series_esd)){
          available_series_esd <- x@series_esd
        } else {
          if((i+1) <= length(nas)) available_series_esd <- x@series_esd[-nas[((i+1):length(nas))]] else available_series_esd <- x@series_esd
        }
        
        # Apply the imputation function to the available data
        imputed_value <- iARinterpolation(coef = x@coef,
                                          series = available_data,
                                          times = available_times,
                                          series_esd = available_series_esd,
                                          yini = yini,
                                          zero_mean = x@zero_mean,
                                          standardized = x@standardized)$fitted
        
        # Replace the NA value with the imputed value
        vector[na_pos] <- imputed_value
      }
      
      
      x@interpolated_values <- vector[nas]
      x@interpolated_times <- x@times[nas]
      x@interpolated_series <- vector
    }
    
    return(x)
  }
  if(x@family == "t"){
    if(length(x@coef)==0) stop("The interpolation method needs the coefficient of the iAR-T model")
    # Original series
    vector <- x@series
    
    # Find positions of NA values in the vector
    nas <- which(is.na(vector))
    
    # If there are NA values, return the imputed vector
    if (length(nas) > 0) {
      
      # Iterate over the positions of NA values
      for (i in seq_along(nas)) {
        # Current position of NA value
        na_pos <- nas[i]
        
        # Data available up to before the current NA value
        if((i+1) <= length(nas)){
          available_data <- vector[-nas[((i+1):length(nas))]]
          available_times <- x@times[-nas[((i+1):length(nas))]]
        }  else {
          available_data <- vector
          available_times <- x@times
        }
        if(is.integer(x@series_esd)){
          available_series_esd <- x@series_esd
        } else {
          if((i+1) <= length(nas)) available_series_esd <- x@series_esd[-nas[((i+1):length(nas))]] else available_series_esd <- x@series_esd
        }
        
        # Apply the imputation function to the available data
        imputed_value <- iARtinterpolation(coef = c(x@coef,x@sigma),
                                           series = available_data,
                                           times = available_times,
                                           df = x@df,
                                           # series_esd = available_series_esd, # es necesario aqui?
                                           # zero_mean = x@zero_mean, # no esta en la funcion original
                                           yini = yini)$fitted
        
        # Replace the NA value with the imputed value
        vector[na_pos] <- imputed_value
      }
      
      
      x@interpolated_values <- vector[nas]
      x@interpolated_times <- x@times[nas]
      x@interpolated_series <- vector
    }
    
    return(x)
  }
  
  if(x@family == "gamma"){
    if(length(x@coef)==0) stop("The interpolation method needs the coefficients of the iAR-Gamma model")
    # Original series
    vector <- x@series
    
    # Find positions of NA values in the vector
    nas <- which(is.na(vector))
    
    # If there are NA values, return the imputed vector
    if (length(nas) > 0) {
      
      # Iterate over the positions of NA values
      for (i in seq_along(nas)) {
        # Current position of NA value
        na_pos <- nas[i]
        
        # Data available up to before the current NA value
        if((i+1) <= length(nas)){
          available_data <- vector[-nas[((i+1):length(nas))]]
          available_times <- x@times[-nas[((i+1):length(nas))]]
        }  else {
          available_data <- vector
          available_times <- x@times
        }
        if(is.integer(x@series_esd)){
          available_series_esd <- x@series_esd
        } else {
          if((i+1) <= length(nas)) available_series_esd <- x@series_esd[-nas[((i+1):length(nas))]] else available_series_esd <- x@series_esd
        }
        
        # Apply the imputation function to the available data
        imputed_value <- iARginterpolation(coef = c(x@coef, x@mean, x@variance),
                                           series = available_data,
                                           times = available_times,
                                           # series_esd = available_series_esd, # preguntar?
                                           yini = yini)$fitted
        
        # Replace the NA value with the imputed value
        vector[na_pos] <- imputed_value
      }
      
      
      x@interpolated_values <- vector[nas]
      x@interpolated_times <- x@times[nas]
      x@interpolated_series <- vector
    }
    
    return(x)
  }
}
method(interpolation, CiAR) <- function(x, yini = 0, seed = NULL) {
  if(length(x@series)==0) stop("The interpolation method needs a time series")
  if(length(x@coef)==0) stop("The interpolation method needs the coefficients of the CiAR model")
  # Original series
  vector <- x@series
  
  # Find positions of NA values in the vector
  nas <- which(is.na(vector))
  
  # If there are NA values, return the imputed vector
  if (length(nas) > 0) {
    
    # Iterate over the positions of NA values
    for (i in seq_along(nas)) {
      # Current position of NA value
      na_pos <- nas[i]
      
      # Data available up to before the current NA value
      if((i+1) <= length(nas)){
        available_data <- vector[-nas[((i+1):length(nas))]]
        available_times <- x@times[-nas[((i+1):length(nas))]]
      }  else{
        available_data <- vector
        available_times <- x@times
      }
      if(is.integer(x@series_esd)){
        available_series_esd <- x@series_esd
      } else {
        if((i+1) <= length(nas)) available_series_esd <- x@series_esd[-nas[((i+1):length(nas))]] else available_series_esd <- x@series_esd
      }
      
      # Apply the imputation function to the available data
      imputed_value <- CiARinterpolation(coef = x@coef,
                                         series = available_data,
                                         times = available_times,
                                         series_esd = available_series_esd,
                                         yini = yini,
                                         zero_mean = x@zero_mean,
                                         standardized = x@standardized,
                                         seed = seed)$fitted
      
      # Replace the NA value with the imputed value
      vector[na_pos] <- imputed_value
    }
    
    
    x@interpolated_values <- vector[nas]
    x@interpolated_times <- x@times[nas]
    x@interpolated_series <- vector
  }
  
  return(x)
}
method(interpolation, BiAR) <- function(x, yini1 = 0, yini2 = 0, seed = NULL) {
  if(length(x@series)==0) stop("The interpolation method needs a bivariate time series")
  if(length(x@coef)==0) stop("The interpolation method needs the coefficients of the BiAR model")
  # Original series
  vector <- x@series
  # Find positions of NA values in the vector
  nas <- sort(unique(unlist(apply(is.na(vector),2,which))))
  auxtimes <- matrix(NA,dim(vector)[1],dim(vector)[2])
  auxseries <-  auxtimes
  # If there are NA values, return the imputed vector
  if (length(nas) > 0) {
    
    # Iterate over the positions of NA values
    for (i in seq_along(nas)) {
      # Current position of NA value
      na_pos <- nas[i]
      
      # Data available up to before the current NA value
      if((i+1) <= length(nas)){
        available_data <- vector[-nas[((i+1):length(nas))],]
        available_times <- x@times[-nas[((i+1):length(nas))]]
      }  else{
        available_data <- vector
        available_times <- x@times
      }
      na_serie=which(is.na(available_data[na_pos,]))
      nmiss=length(na_serie)
      if(is.integer(x@series_esd)) x@series_esd <- matrix(0, nrow=length(x@times),ncol = 2)
      if(is.integer(x@series_esd)){
        available_series_esd <- x@series_esd
      } 
      else {
        x@series_esd[na_pos,na_serie] <- 0
        if((i+1) <= length(nas)) available_series_esd <- x@series_esd[-nas[((i+1):length(nas))],] else available_series_esd <- x@series_esd
      }
      imputed_value <-   iAR:::BiARinterpolation(coef = x@coef,
                                           series1 = available_data[,1],
                                           series2 = available_data[,2],
                                           times = available_times,
                                           series_esd1 = available_series_esd[,1],
                                           series_esd2 = available_series_esd[,2],
                                           yini1 = yini1,
                                           yini2 = yini2,
                                           zero_mean = x@zero_mean,nmiss=nmiss)$fitted
      # Replace the NA value with the imputed value
      vector[na_pos,na_serie] <- imputed_value
      auxtimes[na_pos,na_serie] <- x@times[na_pos]
      auxseries[na_pos,na_serie] <- vector[na_pos,na_serie]
    }
    x@interpolated_values <- auxseries
    x@interpolated_times <- auxtimes
    x@interpolated_series <- vector
  }
  return(x)
}