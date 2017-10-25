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


gerarAlvo = function(p){
	target = c()
	for(i in 1:dim(p)[1]){
		target = c(target, which(p[i,] == max(p[i,])))
	}
	return (target)
}

IEM_roc = c()
IEM_rand = c()
while(length(IEM_rand) < 100){
	print(length(IEM_rand))
	tryCatch({
		retorno1 = IEM(x, tipoMatriz = 3, numRep = 3, plot = FALSE)
		# roc = ROC(test = retorno1$Posteriori[,1], stat = alvo, PV = FALSE, MX = FALSE, MI = FALSE, )
		# if(roc$AUC < 0.5){
			# roc = ROC(test = retorno1$Posteriori[,2], stat = alvo, PV = FALSE, MX = FALSE, MI = FALSE, )
		# }
		# if(roc$AUC <= 1){
			# IEM_roc = c(IEM_roc, roc$AUC)	
		# }
		target = gerarAlvo(retorno1$Posteriori)
		IEM_rand = c(IEM_rand, adjustedRandIndex( target, alvo) )
	}, error=function(e){})
}

IKFCM_roc = c()
IKFCM_rand = c()
while(length(IKFCM_rand) < 100){
	print(length(IKFCM_rand))
	tryCatch({
		retorno1 = IKFCM(x, c = 2, m = 0.1, maxRep = 50, monteCarlo = 5)		
		# roc = ROC(test = retorno1$U[,1], stat = alvo, PV = FALSE, MX = FALSE, MI = FALSE, )
		# if(roc$AUC < 0.5){
			# roc = ROC(test = retorno1$U[,2], stat = alvo, PV = FALSE, MX = FALSE, MI = FALSE, )
		# }
		# if(roc$AUC <= 1){
			# IKFCM_roc = c(IKFCM_roc, roc$AUC)	
		# }
		target = gerarAlvo(retorno1$U)
		IKFCM_rand = c(IKFCM_rand, adjustedRandIndex( target, alvo) )		
	}, error=function(e){})
}

IKEM_roc = c()
IKEM_rand = c()
while(length(IKEM_rand) < 100){
	print(length(IKEM_rand))
	tryCatch({
		retorno1 = IKEM(x, maxRep = 50, kernelType = 0, D = 4, plot = FALSE, numRep = 10)
		# roc = ROC(test = retorno1$Posteriori[,1], stat = alvo, PV = FALSE, MX = FALSE, MI = FALSE, )
		# if(roc$AUC < 0.5){
			# roc = ROC(test = retorno1$Posteriori[,2], stat = alvo, PV = FALSE, MX = FALSE, MI = FALSE, )
		# }
		# if(roc$AUC <= 1){
			# IKEM_roc = c(IKEM_roc, roc$AUC)	
		# }
		target = gerarAlvo(retorno1$Posteriori)
		IKEM_rand = c(IKEM_rand, adjustedRandIndex( target, alvo) )
	}, error=function(e){})
}

IKEM_IS_roc = c()
IKEM_IS_rand = c()
while(length(IKEM_IS_rand) < 100){
	print(length(IKEM_IS_rand))
	tryCatch({
		retorno1 = IKEMIS(x, repPerIteration = 50, plot = FALSE, iterations = 10)
		# roc = ROC(test = retorno1$Posteriori[,1], stat = alvo, PV = FALSE, MX = FALSE, MI = FALSE, )
		# if(roc$AUC < 0.5){
			# roc = ROC(test = retorno1$Posteriori[,2], stat = alvo, PV = FALSE, MX = FALSE, MI = FALSE, )
		# }
		# if(roc$AUC <= 1){
			# IKEM_IS_roc = c(IKEM_IS_roc, roc$AUC)	
		# }
		target = gerarAlvo(retorno1$Posteriori)
		IKEM_IS_rand = c(IKEM_IS_rand, adjustedRandIndex( target, alvo) )		
	}, error=function(e){})
}
