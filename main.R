# Pacotes e inputs --------------------------------------------------------
require("xts")
require("FinTS")
require("forecast")
require("rugarch")
require("parallel")
require("PerformanceAnalytics")
library("readxl")
library("xlsx") # EDIT Isso soh precisa rodar uma vez, nao deve entrar no loop

setwd("D:/Google Drive/Desafio Quantamental/GARCH Vol. forecast/Scripts")

# Carrega funções
#source("01-model-selection.R")
source("02-forecast.R")

# Lê dados e transforma

load("dados.RData")

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
  
  
  source("01-model-selection.R")
  
  # Realiza previsoes, dados:
  first_obs <- ep[2 + z] + 1
  last_obs <- ep[5 + z]
  returns <- all_returns[first_obs:last_obs, ] # Adiciona um ano onde durante vai ser feita a previsao com os modelos escolhidos
  
  print(paste("Criando previsoes usando os dados:", 
              index(first(returns)),
              index(last(returns))))
  
  source("02-forecast.R")
  
  if(!exists("output")){
    output <- list()
    output[[y]] <- out
    
  } else {
    output[[y]] <- out
    y <- y + 1
  }
  
}

