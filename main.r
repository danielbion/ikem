setwd("C:/bitbucket/ikem")
source("functions_sda.r")
library('MASS')
library('Epi')
library('mclust')
source("Interval_EM.r")
source("IKEM_IS.r")
source("Interval_Kernel_EM.r")
source('Interval_Kernel_Fuzzy_C-Means.r')

#######################################################

x = as.matrix(read.table('Agaricus.txt', head = TRUE))
x = normalizarDadosSimbolicos(x)

alvo = rep(1, dim(x)[1])
alvo[c(7, 13, 18, 24)] = 0

#######################################################

x = as.matrix(read.table('temperatura.txt'))
x = x[-c(19,37),]
x = normalizarDadosSimbolicos(x)

alvo = rep(1, dim(x)[1])
alvo[c(1,2,8,10,11,14,15,17,20,21,24,25,26,27,28,30,32,33,34,35)] = 0

#######################################################

x = as.matrix(read.table('carros.txt', head = TRUE))
x = normalizarDadosSimbolicos(x)

alvo = rep(1, dim(x)[1])
alvo[c(4,11,15,22,23)] = 0

#######################################################


criterio1 = c()
while(length(criterio1) < 100){
	print(length(criterio1))
	tryCatch({
		retorno1 = IEM(x, tipoMatriz = 3, numRep = 3, plot = FALSE)
		roc = ROC(test = retorno1$Posteriori[,1], stat = alvo, PV = FALSE, MX = FALSE, MI = FALSE, )
		if(roc$AUC < 0.5){
			roc = ROC(test = retorno1$Posteriori[,2], stat = alvo, PV = FALSE, MX = FALSE, MI = FALSE, )
		}
		if(roc$AUC <= 1){
			criterio1 = c(criterio1, roc$AUC)	
		}
	}, error=function(e){})
}

criterio2 = c()
while(length(criterio2) < 100){
	print(length(criterio2))
	tryCatch({
		retorno1 = IKFCM(x, c = 2, m = 0.1, maxRep = 50, monteCarlo = 5)		
		roc = ROC(test = retorno1$U[,1], stat = alvo, PV = FALSE, MX = FALSE, MI = FALSE, )
		if(roc$AUC < 0.5){
			roc = ROC(test = retorno1$U[,2], stat = alvo, PV = FALSE, MX = FALSE, MI = FALSE, )
		}
		if(roc$AUC <= 1){
			criterio2 = c(criterio2, roc$AUC)	
		}
	}, error=function(e){})
}

criterio3 = c()
while(length(criterio3) < 100){
	print(length(criterio3))
	tryCatch({
		retorno1 = IKEM(x, maxRep = 50, kernelType = 0, D = 4, plot = FALSE, numRep = 10)
		roc = ROC(test = retorno1$Posteriori[,1], stat = alvo, PV = FALSE, MX = FALSE, MI = FALSE, )
		if(roc$AUC < 0.5){
			roc = ROC(test = retorno1$Posteriori[,2], stat = alvo, PV = FALSE, MX = FALSE, MI = FALSE, )
		}
		if(roc$AUC <= 1){
			criterio3 = c(criterio3, roc$AUC)	
		}
	}, error=function(e){})
}

criterio4 = c()
while(length(criterio4) < 100){
	print(length(criterio4))
	tryCatch({
		retorno1 = IKEMIS(x, repPerIteration = 50, plot = FALSE, iterations = 10)
		roc = ROC(test = retorno1$Posteriori[,1], stat = alvo, PV = FALSE, MX = FALSE, MI = FALSE, )
		if(roc$AUC < 0.5){
			roc = ROC(test = retorno1$Posteriori[,2], stat = alvo, PV = FALSE, MX = FALSE, MI = FALSE, )
		}
		if(roc$AUC <= 1){
			criterio4 = c(criterio4, roc$AUC)	
		}
	}, error=function(e){})
}