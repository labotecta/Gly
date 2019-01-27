rm(list = ls())
library(xgboost)
library(caret)
library(data.table)

desamb  <- "_Lucas"

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

param <- list(
  objective   = "binary:logistic",
  eval_metric = "error")

lnr  <- c(625)
leta <- c(0.085)
lcs  <- c(0.93,0.94,0.95,0.96,0.97,0.98,0.99)
lmd  <- c(13)

nl   <- 0
resultados <- data.frame(
  linea  = character(),
  lineac = character(),
  stringsAsFactors = FALSE
)
fi <- paste(sep="",casoa,"_muestra_",porcentaje_entrenar,"Ajusta_Parametros_SP_XGBoost",desamb,".txt")

for(m in 1:length(lnr)){
  nr <- lnr[m]
  for(i in 1:length(leta)){
    eta  <- leta[i]
    for(j in 1:length(lcs)){
      cs   <- lcs[j]
      for(k in 1:length(lmd)){
        md    <- lmd[k]
        set.seed(1235) 
        bst   <- xgboost::xgboost(eta = eta, colsample_bylevel = cs, max_depth = md, data = dentrena, params = param, nrounds = nr, print_every_n = 1000)
        pred  <- predict(bst, dprueba)
        pred  <- ifelse (pred > 0.5, 1, 0)
        mc    <- confusionMatrix(factor(pred,levels=0:1), factor(Eprueba,levels=0:1))
        linea <- paste(sep="",
        "Hyperleda",
        sprintf(" eta:%4.2f",eta),
        sprintf(" cs:%5.2f",cs),
        sprintf(" md:%4.1f",md),
        sprintf(" Kappa:%6.3f",unname(mc$overall[2])),
        sprintf(" McnemarPValue:%5.3f",unname(mc$overall[7])),
        sprintf(" Acierto:%7.3f",unname(mc$overall[1])*100),"%",
        sprintf(" Entrena:%6d No:%6d Si:%6d (%7.3f%s)",nrow(entrena),nrow(entrena)-sum(entrena$E),sum(entrena$E),sum(entrena$E)*100/nrow(entrena),"%"),
        sprintf(" Prueba :%6d No:%6d Si:%6d (%7.3f%s)",nrow(prueba),nrow(prueba)-sum(prueba$E),sum(prueba$E),sum(prueba$E)*100/nrow(prueba),"%"),
        sprintf(" %6.0f",mc$table[1,1]),
        sprintf(" %6.0f",mc$table[2,2])
        )
        lineac  <- gsub('\\.',',',paste(sep="",
        "Hyperleda",
        sprintf(" %4.2f",eta),
        sprintf(" %5.2f",cs),
        sprintf(" %4.1f",md),
        sprintf(" %8.5f",unname(mc$overall[2])),
        sprintf(" %7.5f",unname(mc$overall[7])),
        sprintf(" %7.5f",unname(mc$overall[1])),
        sprintf(" %6d %6d %6d %7.5f",nrow(entrena),nrow(entrena)-sum(entrena$E),sum(entrena$E),sum(entrena$E)/nrow(entrena)),
        sprintf(" %6d %6d %6d %7.5f",nrow(prueba),nrow(prueba)-sum(prueba$E),sum(prueba$E),sum(prueba$E)/nrow(prueba)),
        sprintf(" %6d",mc$table[1,1]),
        sprintf(" %6d",mc$table[2,2])
        ))
        nl <- nl + 1
        resultados[nl,"linea"] = linea
        resultados[nl,"lineac"] = lineac
        cat(sprintf("%s", linea), file=fi, append=TRUE, sep="\n")
      }
    }
  }
}
cat(resultados[,"linea"],fill=1)
cat(resultados[,"lineac"],fill=1)

