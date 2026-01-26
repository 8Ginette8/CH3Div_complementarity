strip_glm.full = function(cm) {
  cm$y = c()
  cm$model = c()
  
  cm$residuals = c()
  cm$fitted.values = c()
  cm$effects = c()
  cm$qr$qr = c()  
  cm$linear.predictors = c()
  cm$weights = c()
  cm$prior.weights = c()
  cm$data = c()

  
  cm$family$variance = c()
  cm$family$dev.resids = c()
  cm$family$aic = c()
  cm$family$validmu = c()
  cm$family$simulate = c()
  attr(cm$terms,".Environment") = c()
  attr(cm$formula,".Environment") = c()
  
  cm
}


strip_glm.light = function(m1){
 m1$data <- NULL
 m1$y <- NULL
 m1$linear.predictors <- NULL
 m1$weights <- NULL
 m1$fitted.values <- NULL
 m1$model <- NULL
 m1$prior.weights <- NULL
 m1$residuals <- NULL
 m1$effects <- NULL
 m1$qr$qr <- NULL
 return(m1)
}