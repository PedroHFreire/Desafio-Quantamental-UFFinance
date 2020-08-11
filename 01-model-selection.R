#modelSelection <- function(returns) {
  # Análise preliminar ------------------------------------------------------
  
  returns_presubset <- xts::first(returns, '2 years')
  
  # Teste de heteroscedasticidade stepwise
  
  wrapperArchTest <- function(lags, x) {
    
    test_out <- FinTS::ArchTest(x = x, lags = lags)
    
    return(test_out$p.value)
  }
  
  max_lag <- 5
  
  arch_tests_pvalues <- unlist(lapply(1:max_lag, wrapperArchTest, x = returns_presubset))
  
  if (sum(arch_tests_pvalues >= 0.05) != 0) {
    #stop( # Edit, voltar com o stop
    print("Hipótese nula de que não há efeito ARCH não foi rejeitada,
          não há heteroscedasticidade a ser modelada.") # EDIT Ta não rejeitando só pro primeiro e quinto lag, o que significa?
  }
  
  
  # Especificação do modelo para a variância --------------------------------
  
  arma_order <- c(0, 0)
  models <- c("sGARCH", "eGARCH", "gjrGARCH")
  distributions <- c("norm", "std", "ged", "snorm", "sstd", "sged")
  parameters <-  expand.grid(p = 0:2, q = 0:2)
  
  garch_specs <- vector(mode = "list",
                        length = length(models) * length(distributions) * nrow(parameters)) # + submodels
  counter <- 1
  
  for (model in models) {
    for (distribution in distributions) {
      for (i in 1:nrow(parameters)) {
      
        garch_specs[[counter]] <- rugarch::ugarchspec(
          
          mean.model         = list(armaOrder = arma_order),
          variance.model     = list(model = model, garchOrder = c(parameters[i, 1], parameters[i, 2])),
          distribution.model = distribution)
        
        counter <- counter + 1
      }
    }
  }
  
  # fGARCH
  submodels <- c("NGARCH", "TGARCH")
  for (submodel in submodels) {
    for (distribution in distributions) {
      for (i in 1:nrow(parameters)) {
        
        garch_specs[[counter]] <- rugarch::ugarchspec(
          
          mean.model         = list(armaOrder = arma_order),
          variance.model     = list(model    = "fGARCH",
                                    garchOrder = c(parameters[i, 1], parameters[i, 2]),
                                    submodel = submodel),
          distribution.model = distribution)
        
        counter <- counter + 1
      }
    }
  }
  
  names(garch_specs) <- as.character(1:(counter - 1))
  
  
  
  # Estimação dos modelos ---------------------------------------------------
  
  cl <- parallel::makeCluster(parallel::detectCores() - 1)
  parallel::clusterEvalQ(cl, { library("rugarch") })
  parallel::clusterExport(cl, varlist = c("garch_specs", "returns_presubset"))
  
  garch_fits <- parallel::parLapply(cl, 
                                    garch_specs, 
                                    fun  = rugarch::ugarchfit, 
                                    data = returns_presubset,
                                    USE.NAMES = TRUE)
  
  parallel::stopCluster(cl)
  
  
  
  
  # Escolhe modelo ----------------------------------------------------------
  
  ## Checa modelos 
  fits <- garch_fits # cria uma copia que sera filtrada
  
  
  # Modelos que não convergiram:
  
  convergence_indicator <- sapply(garch_fits, rugarch::convergence)
  
  fits <- garch_fits[convergence_indicator == 0]
  
  print(paste(length(convergence_indicator[convergence_indicator == 1]), "de", length(garch_fits), "modelos não convergiram."))
  
  
  
  # Teste de autocorrelação serial
  
  
  lm_discard <- c()
  for (i in 1:length(fits)){
    stdreturns <- residuals(fits[[i]], standardize = TRUE)
    if (Box.test(abs(stdreturns), 22, type = "Ljung-Box")$p.value < 0.05) {
      lm_discard <- c(lm_discard, i)
    }
  }
  
  print(paste(length(lm_discard), "de", length(fits), "modelos nao passaram no teste LM."))
  fits <- fits[-lm_discard]
  
  
  # Criterios de informação
  
  
  getBIC <- function(fit){ return(rugarch::infocriteria(fit)[2])} # Retorna o BIC
  bics <- sapply(fits, getBIC)
  
  bic_sup <- mean(bics) + 2 * sd(bics)
  
  bic_discard <- names(bics[bics > bic_sup])
  
  print(paste(length(bic_discard), "de", length(fits), "modelos tem BIC maior que intervalo aceitavel."))
  fits <- fits[!(names(fits) %in% bic_discard)]
  
  
  # Analise pseudo out of sample
  
  garch_specs <- garch_specs[names(fits)] # Seleciona somente os que passaram no teste anterior
  
  out_smpl_size <- length(xts::last(returns, '1 year'))
  window_size <- length(xts::first(returns, '2 years'))
  
  
  cluster <- parallel::makeCluster(parallel::detectCores() - 1)
  parallel::clusterExport(cluster, varlist = c("garch_specs", "returns", "out_smpl_size", "window_size"))
  
  garch_rolls <- lapply(garch_specs,
                        FUN = rugarch::ugarchroll,
                        data = returns, # Amostra completa
                        forecast.length = out_smpl_size, # Out of sample para testar
                        refit.every = 21 * 6, # Periodicidade da re-estimação
                        refit.window = 'moving',
                        window.size = window_size, # Tamanho da janela de dados a serem usadas na estimação
                        calculate.VaR = FALSE,
                        keep.coef = TRUE,
                        cluster = cluster)
  
  parallel::stopCluster(cluster)
  
  garch_rolls <- garch_rolls[sapply(garch_rolls, rugarch::convergence) != 1]
  
  # Critérios de informação e qualidade do ajuste ---------------------------
  
  ## Qualidade do ajuste para o modelo da variancia: erro medio de previsao ao quadrado
  
  garchVarianceMSE <- function(roll) {
    preds <- as.data.frame(roll)
    
    # Erro de previsão para a média
    e <- preds$Realized - preds$Mu
    
    # Erro de previsão para a variância
    d <- e ^ 2 - preds$Sigma ^ 2
    
    return(mean(d ^ 2))
  }
  
  rolls_mses <- vapply(garch_rolls, garchVarianceMSE, FUN.VALUE = numeric(1)) #EDIT resolver problema da convergencia, quando tem janela que nao converge dá erro e essa linha não funciona, tem que colocar um trycatch
  
  best_mses <- rolls_mses[rolls_mses <= summary(rolls_mses)["1st Qu."]]
  
  # Modelos escolhidos:
  chosen_models <- names(best_mses)
  
  if(length(chosen_models) > 10) {
    
    chosen_models <- names(sort(best_mses)[1:10])
    
  }
  
  # Especificacoes dos modelos escolhidos
  chosen_specs <- garch_specs[chosen_models]
  
  #return(chosen_specs)
#}
# EDIT Pesos dos modelos escolhidos: adicionar com otimizacao ou algoritmo de NN (https://www.rleripio.com.br/post/combinando-modelos-de-previsao/)


