# Função que especifica modelo GARCH para uma série de retornos
# EDIT: colocar em função quando terminar


# Pacotes e inputs --------------------------------------------------------
require("xts")
require("FinTS")
require("forecast")
require("rugarch")
require("parallel")
require("PerformanceAnalytics")
library("readxl")
library("xlsx") # EDIT Isso soh precisa rodar uma vez, nao deve entrar no loop

ptm <- proc.time()
# Lê dados e transforma

load("dados.RData")

prices <- dados[, 1] # Pega um ativo, aqui que vai ser colocado o loop, variando a coluna

returns <- PerformanceAnalytics::CalculateReturns(prices, method = "log")[-1]
returns <- first(returns, '5 years') # Filtra a amostra de 5 anos que esta sendo trabalhada

returns <- na.omit(returns) #EDIT cuidado aqui, tem que ver com o backtest...


# Análise preliminar ------------------------------------------------------

returns_presubset <- first(returns, '4 years')

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


# Estima modelo para a média ----------------------------------------------

mean_fit <- forecast::auto.arima(returns_presubset,
                                 ic            = 'bic', # EDIT usar outros criterios?
                                 stepwise      = FALSE,
                                 approximation = FALSE,
                                 parallel      = TRUE,
                                 num.cores     = 7) #edit# tem o autoarfima tb e pode usar outros criterios ic


if (forecast::arimaorder(mean_fit)['d'] != 0) {
  stop('Serie tem componente sazonal, avaliar')
}

ar <- forecast::arimaorder(mean_fit)['p']
ma <- forecast::arimaorder(mean_fit)['q']
arma_order <- c(ar, ma)


# Especificação do modelo para a variância --------------------------------

# EDIT: colocar modelos simples: desvio padrao historico e janela de dados; EWMA; ARCH; GARCH(1, 1) média constante

## Modelos usados:
# Standard GARCH ("sGARCH");
# Exponential GARCH ("eGARCH");
# GJR-GARCH ("gjrGARCH");
# Asymmetric power ARCH ("apARCH");
# Nonlinear ARCH ("fGARCH", submodel: "NGARCH")
# Threshold GARCH ("fGARCH", submodel: "TGARCH")

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

proc.time() - ptm

# Estimação dos modelos ---------------------------------------------------

cl <- parallel::makeCluster(parallel::detectCores() - 1)
clusterEvalQ(cl, { library("rugarch") })
clusterExport(cl, varlist = c("garch_specs", "returns"))

garch_fits <- parallel::parLapply(cl, 
                                  garch_specs, 
                                  fun  = rugarch::ugarchfit, 
                                  data = returns_presubset,
                                  USE.NAMES = TRUE)

stopCluster(cl)




##########
### Choose
## Checa modelos 

# Checando convergência ---------------------------------------------------


# Modelos que não convergiram: # Checar na documentação se isso ta certo

convergence_discard <- c()
for (i in 1:length(garch_fits)) {
  if (garch_fits[[i]]@fit$convergence == 1) {
    convergence_discard <- c(convergence_discard, i)
  }
}

print(paste(length(convergence_discard), "de", length(garch_fits), "modelos não convergiram."))

# Removendo modelos que não convergiram:
fits <- garch_fits[-convergence_discard] 


#####
## Filtra os modelos pelos seguintes fatores:

# Significância dos parâmetros --------------------------------------------


## Analisando significancia dos coeficientes do modelo
significance_discard <- c()
for (i in 1:length(fits)) {
  zero_par <- sum(fits[[i]]@fit$matcoef[, 4] < 0.05) # coeficientes nao significantes
  
  if (zero_par != 0 || is.na(zero_par)) { # A matriz de coef resulta em NA's quando não consegue inverter a hessiana
    significance_discard <- c(significance_discard, i)
    
  }
}

#fits <- fits[-significance_discard] Ta descartando tudo kkk EDIT

# Teste de autocorrelação serial ------------------------------------------


lm_discard <- c()
for (i in 1:length(fits)){
  stdreturns <- residuals(fits[[i]], standardize = TRUE)
  if (Box.test(abs(stdreturns), 22, type = "Ljung-Box")$p.value < 0.05) {
    lm_discard <- c(lm_discard, i)
  }
}

fits <- fits[-lm_discard]


# Criterios de informação? EDIT


# Analise pseudo out of sample --------------------------------------------

garch_specs <- garch_specs[names(fits)] # Seleciona somente os que passaram no teste anterior

cluster <- parallel::makeCluster(parallel::detectCores() - 1)
clusterExport(cluster, varlist = c("garch_specs", "returns"))

garch_rolls <- vector(mode = "list", length = length(garch_specs))

i <- 1
for (garch_spec in garch_specs) {
  garch_rolls[[i]] <- ugarchroll(garch_spec,
                                 data = returns, # 5 anos de dados
                                 forecast.length = length(last(returns, '1 year')), # vai testar para um ano a frente
                                 refit.every = 21, # Re-estima todo mês
                                 refit.window = 'moving',
                                 window.size = 21 * 12 * 4, # Vai usar 4 anos de dados para estimar
                                 calculate.VaR = FALSE,
                                 keep.coef = TRUE,
                                 solver = "solnp",
                                 cluster = cluster)
  i <- i + 1
  
}

stopCluster(cluster)



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

rolls_mses <- vapply(garch_rolls, garchVarianceMSE, FUN.VALUE = numeric(1))
names(rolls_mses) <- names(garch_specs)


# Modelo escolhido:
chosen_one <- names(which.min(rolls_mses))

# Acha o spec que refere ao modelo escolhido
garch_spec <- garch_specs[[chosen_one]]

rm(list=ls()[! ls() %in% c("prices","garch_spec")])
