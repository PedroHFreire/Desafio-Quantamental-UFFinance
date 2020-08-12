
# Preparando tudo ---------------------------------------------------------

# Carrega pacotes e seta diretorio
require("xts")
require("FinTS")
require("forecast")
require("rugarch")
require("parallel")
require("PerformanceAnalytics")
require("readxl")
require("xlsx")

setwd("D:/Google Drive/Desafio Quantamental/GARCH Vol. forecast/Scripts")

## Carrega funções 
wrapperArchTest <- function(lags, x) {
  
  test_out <- FinTS::ArchTest(x = x, lags = lags)
  
  return(test_out$p.value)
}

getBIC <- function(fit){
  bic <- rugarch::infocriteria(fit)[2]
  return(bic)
}

garchVarianceMSE <- function(roll) {
  preds <- as.data.frame(roll)
  
  e <- preds$Realized - preds$Mu
  
  d <- e ^ 2 - preds$Sigma ^ 2
  
  return(mean(d ^ 2))
}

## Especifica modelos GARCH

arma_order <- c(0, 0)
models <- c("sGARCH", "eGARCH", "gjrGARCH")
distributions <- c("norm", "std", "ged", "snorm", "sstd", "sged")
parameters <-  expand.grid(p = 0:1, q = 0:1)

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

## Lê e transforma dados

load("dados.RData")


# Código ------------------------------------------------------------------

prices <- dados$PETR3 # Pega um ativo, aqui que vai ser colocado o loop, variando a coluna


prices <- prices[paste(index(last(prices["2006"])), "2019", sep = "/")] # Amostra valida de dados
xts::endpoints(prices, on = "years")

all_returns <- PerformanceAnalytics::CalculateReturns(prices, method = "log")[-1]
sum(is.na(all_returns))
all_returns <- na.omit(all_returns)
sum(all_returns == 0)


all_returns <- all_returns["2013-07-23/2019-09-07"] # Simplificação para testar código
ep <- xts::endpoints(all_returns, on = "years")

## Loop 
y <- 1
for (z in 0:length(ep)) {
  
  first_obs <- ep[1 + z] + 1
  last_obs <- ep[4 + z]
  
  # Seleciona os modelos, dados:
  returns <- all_returns[first_obs:last_obs, ] # Primeiros 3 anos
  
  print(paste("Selecionando o modelo usando os dados:", 
              index(first(returns)),
              index(last(returns))))
  
  #####
  #####
  ### Model selection
  # Análise preliminar ------------------------------------------------------
  
  returns_presubset <- xts::first(returns, '2 years')
  
  # Teste de heteroscedasticidade stepwise
  
  max_lag <- 5
  
  arch_tests_pvalues <- unlist(lapply(1:max_lag, wrapperArchTest, x = returns_presubset))
  
  if (sum(arch_tests_pvalues >= 0.05) != 0) {
    #stop( # Edit, voltar com o stop
    print("Hipótese nula de que não há efeito ARCH não foi rejeitada,
          não há heteroscedasticidade a ser modelada.") # EDIT Ta não rejeitando só pro primeiro e quinto lag, o que significa?
  }
  
  
  # Estimação dos modelos ---------------------------------------------------
  
  cluster <- parallel::makeCluster(parallel::detectCores() - 1)
  parallel::clusterEvalQ(cluster, { library("rugarch") })
  parallel::clusterExport(cluster, varlist = c("garch_specs", "returns_presubset"))
  
  garch_fits <- parallel::parLapply(cluster, 
                                    garch_specs, 
                                    fun  = rugarch::ugarchfit, 
                                    data = returns_presubset,
                                    solver.control = list(tol = 1e-12),
                                    USE.NAMES = TRUE)
  
  parallel::stopCluster(cluster)
  
  
  
  
  # Escolhe modelo ----------------------------------------------------------
  fits <- garch_fits # cria uma copia que sera filtrada
  
  
  # Filtra modelos que não convergiram:
  convergence_indicator <- sapply(garch_fits, rugarch::convergence)
  
  fits <- garch_fits[convergence_indicator == 0]
  
  print(paste(length(convergence_indicator[convergence_indicator == 1]), "de", length(garch_fits), "modelos não convergiram."))
  
  
  # Filtra pelo teste de autocorrelação serial
  lm_discard <- c()
  for (i in 1:length(fits)){
    stdreturns <- residuals(fits[[i]], standardize = TRUE)
    if (Box.test(abs(stdreturns), 22, type = "Ljung-Box")$p.value < 0.05) {
      lm_discard <- c(lm_discard, i)
    }
  }
  
  print(paste(length(lm_discard), "de", length(fits), "modelos nao passaram no teste LM."))
  fits <- fits[-lm_discard]
  
  
  # Filtra pelos criterios de informação
  
  bics <- sapply(fits, getBIC)
  
  bic_sup <- mean(bics) + 2 * sd(bics)
  
  bic_discard <- names(bics[bics > bic_sup])
  
  print(paste(length(bic_discard), "de", length(fits), "modelos tem BIC maior que intervalo aceitavel."))
  fits <- fits[!(names(fits) %in% bic_discard)]
  
  
  # Filtra pelo MSE na analise pseudo out of sample
  
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
                        solver.control = list(tol = 1e-12),
                        cluster = cluster)
  
  parallel::stopCluster(cluster)
  
  garch_rolls <- garch_rolls[sapply(garch_rolls, rugarch::convergence) != 1]
  
  rolls_mses <- vapply(garch_rolls, garchVarianceMSE, FUN.VALUE = numeric(1)) 
  
  best_mses <- rolls_mses[rolls_mses <= summary(rolls_mses)["1st Qu."]]
  
  
  
  # Output selection: especificacoes dos modelos escolhidos
  chosen_models <- names(best_mses)
  
  if(length(chosen_models) > 5) {
    
    chosen_models <- names(sort(best_mses)[1:5])
    
  }
  
  chosen_specs <- garch_specs[chosen_models]
  
  
  #####
  #####
  ### Forecast
  
  # Realiza previsoes, dados:
  first_obs <- ep[2 + z] + 1
  last_obs <- ep[5 + z]
  returns <- all_returns[first_obs:last_obs, ] # Adiciona um ano onde durante vai ser feita a previsao com os modelos escolhidos
  
  print(paste("Criando previsoes usando os dados:", 
              index(first(returns)),
              index(last(returns))))
  
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
                         solver.control = list(tol = 1e-12),
                         cluster = cluster)
  
  parallel::stopCluster(cluster)
  
  
  chosen_rolls <- chosen_rolls[sapply(chosen_rolls, rugarch::convergence) != 1] 
  
  
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
  
  # Output forecast: maior e menor desvio padrao previsto pelos modelos
  out <- list("Intervalo de previsao"= xts::xts(pred_output, order.by = as.Date(rownames(as.data.frame(roll)))),
              "Modelos escolhidos" = names(chosen_rolls))
  
  
  ## Output do codigo
  if(!exists("output")){
    output <- list()
    output[[y]] <- out
    
  } else {
    output[[y]] <- out
    y <- y + 1
  }
  
}

