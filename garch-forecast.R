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
library("xlsx")

DadosIndustria363 <- readxl::read_excel("D:/Download/DadosIndustria363.xlsx", 
                                        col_types = c("date", "numeric"))
bova_asset <- DadosIndustria363

prices <- xts::xts(x = bova_asset$BOVA11, order.by = bova_asset$Data)

returns <- PerformanceAnalytics::CalculateReturns(prices, method = "log") 
returns <- returns["2010/2012-01-01"]
colnames(returns) <- "BOVA11"


# Análise preliminar ------------------------------------------------------

# Teste de heteroscedasticidade stepwise

wrapperArchTest <- function(lags, x) {
  
  test_out <- FinTS::ArchTest(x = x, lags = lags)
  
  return(test_out$p.value)
}

max_lag <- 5

arch_tests_pvalues <- unlist(lapply(1:max_lag, wrapperArchTest, x = returns))

if (sum(arch_tests_pvalues >= 0.05) != 0) {
  #stop( # Edit, voltar com o stop
  print("Hipótese nula de que não há efeito ARCH não foi rejeitada,
        não há heteroscedasticidade a ser modelada.") # EDIT Ta não rejeitando só pro primeiro e quinto lag, o que significa?
}


# Estima modelo para a média ----------------------------------------------


mean_fit <- forecast::auto.arima(returns,
                                 ic            = 'bic', # EDIT usar outros criterios?
                                 stepwise      = FALSE,
                                 approximation = FALSE,
                                 parallel      = TRUE,
                                 num.cores     = 7) #edit# tem o autoarfima tb e pode usar outros criterios ic

show(mean_fit)

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

# Estimação dos modelos ---------------------------------------------------

cl <- parallel::makeCluster(parallel::detectCores() - 1)
clusterEvalQ(cl, { library("rugarch") })
clusterExport(cl, varlist = c("garch_specs", "returns"))

garch_fits <- parallel::parLapply(cl, 
                                  garch_specs, 
                                  fun  = rugarch::ugarchfit, 
                                  data = returns,
                                  USE.NAMES = TRUE)

stopCluster(cl)




##########
## Choose
# Checa modelos -----------------------------------------------------------
## EDIT aqui coloca uma etapa para escolher os modelos
## Falta testar: stability (Nymblom); biasness e signs of estimated coefficients
# Baseado no paper do SeR # Rever isso
# Colocar model-averaging
# Falta considerar forecasting measures (backtest com ugarchroll?)
# Usar defined namespaces

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
    print(zero_par)
  }
}

fits <- fits[-significance_discard]

# Teste de autocorrelação serial ------------------------------------------


lm_discard <- c()
for (i in 1:length(fits)){
  stdreturns <- residuals(fits[[i]], standardize = TRUE)
  if (Box.test(abs(stdreturns), 22, type = "Ljung-Box")$p.value < 0.05) {
    lm_discard <- c(lm_discard, i)
  }
}

fits <- fits[-lm_discard]




# Critérios de informação e qualidade do ajuste ---------------------------

info_criterias <- vapply(fits,
                         rugarch::infocriteria,
                         FUN.VALUE = numeric(4),
                         USE.NAMES = FALSE)

## Qualidade do ajuste para o modelo da variancia: erro medio de previsao ao quadrado

garchVarianceMSE <- function(fit) {
  e <- residuals(fits[[1]])
  d <- e ^ 2 - sigma(fits[[1]]) ^ 2
  return(mean(d ^ 2))
}

fits_mses <- vapply(fits, garchVarianceMSE, FUN.VALUE = numeric(1))

comparison_mat <- rbind(fits_mses, info_criterias)

rownames(comparison_mat) <- c("Variance MSE", "Akaike", "Bayes", "Shibata", "Hannan-Quinn")

# Ordena cada um do menor para o maior e coloca o rank, a menor soma dos ranks é o melhor modelo
comparison_mat_ranks <- rbind(rank(comparison_mat[1, ]),
                              rank(comparison_mat[2, ]),
                              rank(comparison_mat[3, ]),
                              rank(comparison_mat[4, ]),
                              rank(comparison_mat[5, ]))

rownames(comparison_mat_ranks) <- rownames(comparison_mat)

# Modelo escolhido:
chosen_one <- names(which.min(colSums(comparison_mat_ranks))) # EDIT em caso de empate dá erro mas deveria escolher o mais simples


##########
## Forecast

# Previsão da volatilidade futura -----------------------------------------

# Os backtest começa em 2009-12-31, dai são usados log retornos dos últimos 4 anos (2006-2009) 
# para estimar os parâmetros e escolher o modelo, com isso esse modelo é usado pelos próximos 2 anos
# com os parâmetros re-estimados a cada trimestre, a previsão

# Acha o spec que refere ao modelo escolhido
garch_spec <- garch_specs[[chosen_one]]

# Filtra com os novos dados mantendo os parametros e preve dois dias a frente


# Re-estima os parametros a cada 6 meses
# Re-esolhe o modelo a cada dois anos

n_ahead <- 1 # Quer prever a volatilidade para quantos periodos a frente?


garch_boot <- rugarch::ugarchboot(garch_spec,
                                  data = returns,
                                  method = "Partial",
                                  n.ahead = n_ahead,
                                  n.bootpred = 10000)

vol_forecast <- rugarch::sigma(garch_boot@forc)

boot_pred_interval <- cbind(t(as.data.frame(garch_boot, which = "sigma", type = "summary")),  vol_forecast)


xlsx::write.xlsx(boot_pred_interval, "garch_vol.xlsx")


