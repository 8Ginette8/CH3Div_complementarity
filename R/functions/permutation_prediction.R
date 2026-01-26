# ### =====================================================================
# ### Define speed.VarImp function
# ### =====================================================================

# ### This function is based on the principle of permutation importance. This metric is applicable to any type of model.
# ### The basic idea is to consider a variable important if it has a positive effect on the prediction accuracy (classification), or MSE (regression).
# ### The risk of using this metric is a potential bias towards collinear predictive variables.

speed.VarImp <- function (model, data, nperm, type, rescale) {
   if (type=="response") {
     if (class(model)[1] == "gbm") {
       ref <- predict(model, data, type="response", n.trees=model$n.trees)
     } else {
       ref <- predict(model, data, type="response")
     }
   }
   if (type=="prob") {
     ref <- predict(model, data, type="prob")[,2]
   }
   VarImp <- vector()
   for (i in 1:ncol(data)) {
     print(names(data)[i])
     refi <- vector()
     for (j in 1:nperm) {
       cali <- data
       cali[,i] <- cali[sample(1:nrow(cali), nrow(cali)), i]
       if (type=="response") {
         if (class(model)[1] == "gbm") {
           refi <- c(refi, 1 - cor(ref, predict(model, cali, type="response", n.trees=model$n.trees)))
         } else {
           refi <- c(refi, 1 - cor(ref, predict(model, cali, type="response")))
         }
       }
       if (type=="prob") {
         refi <- c(refi, 1 - cor(ref, predict(model, cali, type="prob")[,2]))
       }
     }
     VarImp <- c(VarImp, round(mean(refi), 3))
   }
   names(VarImp) <- names(data)
   if (rescale==TRUE) {
     VarImp <- VarImp/sum(VarImp)*100
   }
   return(VarImp)
}



# # Et voilà comment je l'appelle dans R:
# ### Variable Importance
# #VarImp.glm <- speed.VarImp(model=model.glm, data=Data_SDM[,Predictors], nperm=1, rescale=TRUE, type="response")
# #VarImp.gam <- speed.VarImp(model=model.gam, data=Data_SDM[,Predictors], nperm=1, rescale=TRUE, type="response")
# #VarImp.gbm <- speed.VarImp(model=model.gbm, data=Data_SDM[,Predictors], nperm=1, rescale=TRUE, type="response")
# #VarImp.rdf <- speed.VarImp(model=model.rdf, data=Data_SDM[,Predictors], nperm=1, rescale=TRUE, type="prob")
# #VarImp.max <- speed.VarImp(model=model.max, data=Data_SDM[,Predictors], nperm=1, rescale=TRUE, type="response")
# #VarImp <- list(glm=VarImp.glm, gam=VarImp.gam, gbm=VarImp.gbm, randomForest=VarImp.rdf, MaxEnt=VarImp.max)