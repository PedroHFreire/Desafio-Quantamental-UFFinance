#rollingForecast <- function(returns, chosen_specs) {
  cluster <- parallel::makeCluster(parallel::detectCores() - 1)
  parallel::clusterExport(cluster, varlist = c("chosen_specs",
                                               "returns"))
  
  
  init_smpl_size <- length(xts::first(returns, '2 years'))
  
  chosen_rolls <- lapply(chosen_specs,
                         FUN = rugarch::ugarchroll,
                         data = returns,
                         n.start = init_smpl_size,
                         refit.window = "moving",
                         refit.every = 1,
                         calculate.VaR = FALSE,
                         keep.coef = TRUE,
                         cluster = cluster)
  
  parallel::stopCluster(cluster)
  
  
  chosen_rolls <- chosen_rolls[sapply(chosen_rolls, rugarch::convergence) != 1]# Tirando modelos que nao convergiram EDIT colocar trycatch pq se tiver janela que nao converge dá ruim, pode usar o método resume para re-estimar também
  
  
  for (roll in chosen_rolls) {
    sigma <- as.data.frame(roll)$Sigma
    
    if(!exists("pred_interval")){
      pred_interval <- sigma
      
    } else {
      pred_interval <- cbind(pred_interval, sigma)
    }
    
  }
  
  colnames(pred_interval) <- names(chosen_rolls)
  
  
  pred_output <- cbind(apply(pred_interval, 1, min),
                       apply(pred_interval, 1, max))
  
  colnames(pred_output) <- c("Min pred.", "Max pred.")
  
  # Output
  # EDIT  numero do semestre
  out <- list("Intervalo de previsao"= xts::xts(pred_output, order.by = as.Date(rownames(as.data.frame(roll)))),
              "Modelos escolhidos" = names(chosen_rolls))
  
#  return(out)
#}
# Falta adicionar a re-selecao
# Vale comparar esses acima via MSE do garch roll aqui no forecast (real out of sample)
