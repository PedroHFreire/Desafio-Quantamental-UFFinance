# Previsão out-of-sample -------------------------------------------------- 
# Vai começar com os últimos 4 anos da primeira amostra de 5 anos
# Então vai prever um dia a frente, tirar o primeiro dia da amostra e adicionar o retorno realizado em d+1

setwd("D:/Google Drive/Desafio Quantamental/GARCH Vol. forecast/Scripts")
source("01-model-selection.R")

returns <- PerformanceAnalytics::CalculateReturns(prices, method = "log")[-1] 
colnames(returns) <- "BOVA11"

# Start (iteração 1)
returns_first <- first(returns, '5 years')
returns_subset <- last(returns_first, '4 years')


garch_fit <- ugarchfit(spec = garch_spec, returns_subset)

garch_boot <- rugarch::ugarchboot(garch_fit,
                                  data = returns_subset,
                                  method = "Partial",
                                  n.ahead = 1,
                                  n.bootpred = 10000)

vol_forecast <- rugarch::sigma(garch_boot@forc)

day_one <- which(index(returns) == index(first(returns_subset))) + 1 # Posição do segundo dia em "returns"
day_ahead <- which(index(returns) == index(last(returns_subset))) + 1 # posição do D + 1 em "returns"


# Outras iterações

# Rodar um for aqui para rodar a seleção a cada trimestre
forecast_subset <- returns[paste(index(returns[day_one]), "/")]


# Rodar um for aqui para cada dia de um trimestre

quarters <- endpoints(forecast_subset, on = "quarters")[-1]
quarters <- quarters[-length(quarters)]


j <- 1
for(quarter in quarters) {
  print(paste(j, quarter))
  quarter_subset <- forecast_subset[j:quarter]
  j <- quarter + 1
  
  for(i in 1:length(quarter_subset)) {
    returns_subset <- returns[day_one:day_ahead]
    
    garch_fit <- ugarchfit(spec = garch_spec, returns_subset)
    
    garch_boot <- rugarch::ugarchboot(garch_fit,
                                      data = returns_subset,
                                      method = "Partial",
                                      n.ahead = 1,
                                      n.bootpred = 10000)
    
    vol_forecast <- rugarch::sigma(garch_boot@forc)
    
    if(!exists("boot_pred_interval")){
      boot_pred_interval <- cbind(t(as.data.frame(garch_boot,
                                                  which = "sigma",
                                                  type = "summary")),
                                  vol_forecast)
      
    } else {
      boot_pred_interval <- rbind(boot_pred_interval,
                                  cbind(t(as.data.frame(garch_boot,
                                                        which = "sigma",
                                                        type = "summary")),
                                        vol_forecast))
    }
    
    day_one <- which(index(returns) == index(first(returns_subset))) + 1
    day_ahead <- which(index(returns) == index(last(returns_subset))) + 1
    
  }
  
  # re-select model
}

# se der erro lança um trycatch para usar o melhor modelo anterior a esse


