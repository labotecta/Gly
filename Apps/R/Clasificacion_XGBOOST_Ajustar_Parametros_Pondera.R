rm(list = ls())
library(xgboost)
library(caret)
library(data.table)

desamb <- "_Lucas"

disco   <- "C"
marcas  <- "444"
dmina   <- "600000"
dmaxa   <- "1000000"
anguloa <- "_8_16_15_70_"
dminb   <- "600000"
dmaxb   <- "1000000"
angulob <- "_8_16_15_70_"

casoa   <- paste(sep="",disco,":/Gly/G_1_10_",marcas,anguloa,dmina,"_",dmaxa)
casob   <- paste(sep="",disco,":/Gly/G_1_10_",marcas,angulob,dminb,"_",dmaxb)
porcentaje_entrenar <- ""

datosentrenar <- read.table(paste(sep="",casoa,"_muestra_",porcentaje_entrenar,"tot_meA.csv"), header=TRUE, sep=";", na.strings="NA", dec=",", strip.white=TRUE)
datosprobar   <- read.table(paste(sep="",casob,"_muestra_",porcentaje_entrenar,"tot_mpA.csv"), header=TRUE, sep=";", na.strings="NA", dec=",", strip.white=TRUE)
entrena  <- datosentrenar[,-c(5)]
prueba   <- datosprobar  [,-c(5)]
ne <- nrow(datosentrenar)
ponderae <- ifelse(datosentrenar$n == 0, 1, datosentrenar$n)
setDT(entrena)
setDT(prueba)
Eentrena <- entrena$E 
Eprueba <- prueba$E
coordenadasEntrena  <- model.matrix(~.+0, data = entrena [,-c("E"), with=F]) 
coordenadasPrueba   <- model.matrix(~.+0, data = prueba  [,-c("E"), with=F])
Eentrena  <- as.numeric(Eentrena)
Eprueba   <- as.numeric(Eprueba)
dentrena  <- xgb.DMatrix(data = coordenadasEntrena,  label = Eentrena) 
dprueba   <- xgb.DMatrix(data = coordenadasPrueba,   label = Eprueba)

logregobj <- function(preds, dtrain) {
  labels <- getinfo(dtrain, "label")
  preds  <- 1/(1 + exp(-preds))
  grad   <- (preds - labels) * ponderae
  hess   <- preds * (1 - preds) * ponderae
  return(list(grad = grad, hess = hess))
}
evalerror <- function(preds, dtrain) {
  labels <- getinfo(dtrain, "label")
  err <- as.numeric(sum(labels != (preds > 0)))/length(labels)
  return(list(metric = "error", value = err))
}

param <- list(
  objective    = logregobj
  ,eval_metric = evalerror
)

ss   <- 1
cb   <- 1
lnr  <- c(618,620,622)
leta <- c(0.083,0.084,0.085)
lcs  <- c(0.92,0.93,0.94,0.95,0.96)
lmd  <- c(12,13)

nl   <- 0
resultados <- data.frame(
	linea  = character(),
	stringsAsFactors = FALSE
)
fi <- paste(sep="",casoa,"_muestra_",porcentaje_entrenar,"Ajusta_Parametros_XGBoost",desamb,".txt")

for(m in 1:length(lnr)){
  nr <- lnr[m]
  for(i in 1:length(leta)){
    eta  <- leta[i]
    for(j in 1:length(lcs)){
      cs   <- lcs[j]
      for(k in 1:length(lmd)){
        md    <- lmd[k]
        set.seed(1235) 
        bst   <- xgb.train(subsample = ss, colsample_bylevel = cs, colsample_bytree = cb, eta = eta, max_depth = md, data = dentrena, params = param, nrounds = nr, print_every_n = 1000)
        pred  <- predict(bst, dprueba)
        pred  <- ifelse (pred > 0.5, 1, 0)
        mc    <- confusionMatrix(factor(pred,levels=0:1), factor(Eprueba,levels=0:1))
        linea <- paste(sep="",
          "Hyperleda",
          sprintf(" nr:%6d",nr),
          sprintf(" eta:%6.3f",eta),
          sprintf(" cs:%5.2f",cs),
          sprintf(" md:%5.1f",md),
          sprintf(" Kappa:%6.3f",unname(mc$overall[2])),
          sprintf(" McnemarPValue:%6.3f",unname(mc$overall[7])),
          sprintf(" Acierto:%7.3f",unname(mc$overall[1])*100),"%",
          sprintf(" Entrena:%6d No:%6d Si:%6d ( %7.3f%s )",nrow(entrena),nrow(entrena)-sum(entrena$E),sum(entrena$E),sum(entrena$E)*100/nrow(entrena),"%"),
          sprintf(" Prueba :%6d No:%6d Si:%6d ( %7.3f%s )",nrow(prueba),nrow(prueba)-sum(prueba$E),sum(prueba$E),sum(prueba$E)*100/nrow(prueba),"%"),
          sprintf(" %6.0f",mc$table[1,1]),
          sprintf(" %6.0f",mc$table[2,2])
        )
        nl <- nl + 1
        resultados[nl,"linea"] = linea
        cat(sprintf("%s", linea), file=fi, append=TRUE, sep="\n")
      }
    }
  }
}
cat(resultados[,"linea"],fill=1)

