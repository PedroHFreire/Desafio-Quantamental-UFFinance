setwd("D:/Google Drive/Desafio Quantamental/GARCH Vol. forecast/Scripts")

dados <- read.csv("cotacoes_ativos_inicio_2000.csv",
                  sep = ";",
                  header = FALSE,
                  dec = ",",
                  stringsAsFactors = FALSE)


dados <- dados[, c(-3, -5)]
colnames(dados) <- c("Data", "Preco", "Ativo")

dados[1, "Data"] <- "2000-01-03 00:00:00.000"


dados$Data <- sapply(dados$Data, 
                     FUN = substr, 
                     start = 1, 
                     stop = 10,
                     USE.NAMES = FALSE)


dados <- reshape(dados,
                 idvar = "Data",
                 timevar = "Ativo",
                 direction = "wide")


colnames(dados)[-1] <- sapply(colnames(dados)[-1],
                              FUN = substr,
                              start = 7,
                              stop = 10000,
                              USE.NAMES = FALSE)


dados <- dados[order(dados$Data), ]


dados$Data <- as.Date(dados$Data)

library(xts)

dados <- xts::xts(x = dados[, -1], order.by = dados$Data)

dados <- dados["2004-12-31/"]

save(dados, file = "dados.RData")


# Fazendo double checks
