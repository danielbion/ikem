setwd("C:/bitbucket/ikem")
source("functions_sda.r")
library('MASS')
library('Epi')
library('mclust')
source("Interval_EM.r")
source("IKEM_IS.r")
source("Interval_Kernel_EM.r")
source('Interval_Kernel_Fuzzy_C-Means.r')

getData = function(base, scenario){
		return(
		switch(base, 
			switch(scenario, 
				base1(n1=100, n2=100, var1=3, var2=3, range1=2, range2=2, offset2=c(10,5)),
				base1(100, 100, 3, 1, 2, 2, c(10,5)),
				base1(20, 20, 3, 3, 5, 1, c(0, 0))
			),
			switch(scenario, 
				base2(n1=300, n2=80, var2=0.1, range1=0.5, range2=0.5, offset2=c(-1.5, 1.5)),
				base2(n1=300, n2=80, var2=0.1, range1=0.5, range2=0.5, offset2=c(-2, 2)),
				base2(n1=300, n2=80, var2=0.1, range1=0.2, range2=0.5, offset2=c(-2, 2))
			),
			switch(scenario, 
				base3(n1 = 100, n2=80, range1=1, range2=1, var1=0.3),
				base3(n1 = 100, n2=80, range1=1, range2=0.5, var1=0.4),
				base3(n1 = 100, n2=80, range1=1, range2=0.5, var1=0.8)
			),
		))
}

MonteCarloIEM = function (mc, base, scenario){
	criterio = c()
	for(i in 1:mc){
		data1 = getData(base, scenario)
		x = data1[[1]]
		alvo = data1[[2]]
		
		retorno1 = IEM(x, tipoMatriz = 3, numRep = 1, plot = TRUE)
		
		roc = ROC(test = retorno1$Posteriori[,1], stat = alvo, PV = FALSE, MX = FALSE, MI = FALSE, )
		if(roc$AUC < 0.5){
			roc = ROC(test = retorno1$Posteriori[,2], stat = alvo, PV = FALSE, MX = FALSE, MI = FALSE, )
		}
		criterio = c(criterio, roc$AUC)	
	}
	return(criterio)
}

MonteCarloIKEM = function (mc, base, scenario){
	criterio = c()
	for(i in 1:mc){
		data1 = getData(base, scenario)
		x = data1[[1]]
		alvo = data1[[2]]
		
		retorno1 = IKEM(x, maxRep = 3, kernelType = 0, D = 13, plot = TRUE)
		
		roc = ROC(test = retorno1$Posteriori[,1], stat = alvo, PV = FALSE, MX = FALSE, MI = FALSE, )
		if(roc$AUC < 0.5){
			roc = ROC(test = retorno1$Posteriori[,2], stat = alvo, PV = FALSE, MX = FALSE, MI = FALSE, )
		}
		criterio = c(criterio, roc$AUC)	
	}
	return(criterio)
}

MonteCarloIKFCM = function (mc, base, scenario){
	criterio = c()
	for(i in 1:mc){
		data1 = getData(base, scenario)
		x = data1[[1]]
		alvo = data1[[2]]
		
		retorno1 = IKFCM(x, c = 2, m = 0.1, maxRep = 1, monteCarlo = 1)
		
		roc = ROC(test = retorno1$U[,1], stat = alvo, PV = FALSE, MX = FALSE, MI = FALSE, )
		if(roc$AUC < 0.5){
			roc = ROC(test = retorno1$U[,2], stat = alvo, PV = FALSE, MX = FALSE, MI = FALSE, )
		}
		criterio = c(criterio, roc$AUC)	
	}
	return(criterio)
}

criterio1_1 = MonteCarloIKFCM(100, 1, 1)
criterio1_2 = MonteCarloIKFCM(100, 1, 2)
criterio1_3 = MonteCarloIKFCM(100, 1, 3)
criterio2_1 = MonteCarloIKFCM(100, 2, 1)
criterio2_2 = MonteCarloIKFCM(100, 2, 2)
criterio2_3 = MonteCarloIKFCM(100, 2, 3)
criterio3_1 = MonteCarloIKFCM(100, 3, 1)
criterio3_2 = MonteCarloIKFCM(100, 3, 2)
criterio3_3 = MonteCarloIKFCM(100, 3, 3)


EMcriterio1_1 = MonteCarloIEM(100, 1, 1)
EMcriterio1_2 = MonteCarloIEM(100, 1, 2)
EMcriterio1_3 = MonteCarloIEM(100, 1, 3)
EMcriterio2_1 = MonteCarloIEM(100, 2, 1)
EMcriterio2_2 = MonteCarloIEM(100, 2, 2)
EMcriterio2_3 = MonteCarloIEM(100, 2, 3)
EMcriterio3_1 = MonteCarloIEM(100, 3, 1)
EMcriterio3_2 = MonteCarloIEM(100, 3, 2)
EMcriterio3_3 = MonteCarloIEM(100, 3, 3)

#######################################################

x = as.matrix(read.table('Agaricus.txt', head = TRUE))
x = normalizarDadosSimbolicos(x)

alvo = rep(1, dim(x)[1])
alvo[c(7, 13, 18, 24)] = 2

x_interval = interval(x[,1],x[,2])
y_interval = interval(x[,3],x[,4])	
z_interval = interval(x[,5],x[,6])	
intervalGraph3D(x_interval, y_interval, z_interval, alvo = alvo)

alvo[c(7, 13, 18, 24)] = 0

# retorno1 = IEM(x, tipoMatriz = 3, numRep = 1, plot = TRUE)
# ROC(test = retorno1$Posteriori[,1], stat = alvo, PV = FALSE, MX = FALSE, MI = FALSE, )

# retorno1 = IKEM(x, maxRep = 3, kernelType = 0, D = 13, plot = TRUE)
# ROC(test = retorno1$Posteriori[,1], stat = alvo, PV = FALSE, MX = FALSE, MI = FALSE, )

# retorno1 = IKFCM(x, c = 2, m = 0.1, maxRep = 1, monteCarlo = 1)
# ROC(test = retorno1$U[,1], stat = alvo, PV = FALSE, MX = FALSE, MI = FALSE, )


criterio = c()
for(i in 1:100){
	retorno1 = IEM(x, tipoMatriz = 3, numRep = 1, plot = FALSE)
	roc = ROC(test = retorno1$Posteriori[,1], stat = alvo, PV = FALSE, MX = FALSE, MI = FALSE, )
	if(roc$AUC < 0.5){
		roc = ROC(test = retorno1$Posteriori[,2], stat = alvo, PV = FALSE, MX = FALSE, MI = FALSE, )
	}
	criterio = c(criterio, roc$AUC)	
}

criterio = c()
for(i in 1:100){
	retorno1 = IKFCM(x, c = 2, m = 0.1, maxRep = 20, monteCarlo = 10)		
	roc = ROC(test = retorno1$U[,1], stat = alvo, PV = FALSE, MX = FALSE, MI = FALSE, )
	if(roc$AUC < 0.5){
		roc = ROC(test = retorno1$U[,2], stat = alvo, PV = FALSE, MX = FALSE, MI = FALSE, )
	}
	criterio = c(criterio, roc$AUC)	
}

criterio = c()
for(i in 1:100){
	retorno1 = IKEM(x, maxRep = 5, kernelType = 0, D = 10, plot = FALSE, numRep = 1)
	roc = ROC(test = retorno1$Posteriori[,1], stat = alvo, PV = FALSE, MX = FALSE, MI = FALSE, )
	if(roc$AUC < 0.5){
		roc = ROC(test = retorno1$Posteriori[,2], stat = alvo, PV = FALSE, MX = FALSE, MI = FALSE, )
	}
	criterio = c(criterio, roc$AUC)	
}

criterio = c()
for(i in 1:100){
	retorno1 = IKEMIS(x, repPerIteration = 50, plot = FALSE, iterations = 1)
	roc = ROC(test = retorno1$Posteriori[,1], stat = alvo, PV = FALSE, MX = FALSE, MI = FALSE, )
	if(roc$AUC < 0.5){
		roc = ROC(test = retorno1$Posteriori[,2], stat = alvo, PV = FALSE, MX = FALSE, MI = FALSE, )
	}
	criterio = c(criterio, roc$AUC)	
}


#######################################################

x = as.matrix(read.table('temperatura.txt'))
x = x[-c(19,37),]
x = normalizarDadosSimbolicos(x)

alvo = rep(1, dim(x)[1])
alvo[c(1,2,8,10,11,14,15,17,20,21,24,25,26,27,28,30,32,33,34,35)] = 2

x_interval = interval(x[,1],x[,2])
y_interval = interval(x[,3],x[,4])	
z_interval = interval(x[,5],x[,6])	
# intervalGraph3D(x_interval, y_interval, z_interval, alvo = alvo)

alvo[c(1,2,8,10,11,14,15,17,20,21,24,25,26,27,28,30,32,33,34,35)] = 0
# retorno1 = IEM(x, tipoMatriz = 3, numRep = 1, plot = FALSE)
# ROC(test = retorno1$Posteriori[,1], stat = alvo, PV = FALSE, MX = FALSE, MI = FALSE, )
# adjustedRandIndex( retorno1$Posteriori[,1], alvo)

# retorno1 = IKEM(x, maxRep = 3, kernelType = 0, D = 13, plot = FALSE)
# ROC(test = retorno1$Posteriori[,1], stat = alvo, PV = FALSE, MX = FALSE, MI = FALSE, )
# adjustedRandIndex( retorno1$Posteriori[,1], alvo)

# retorno1 = IKFCM(x, c = 2, m = 0.1, maxRep = 1, monteCarlo = 1)
# ROC(test = retorno1$U[,1], stat = alvo, PV = FALSE, MX = FALSE, MI = FALSE, )
# adjustedRandIndex( retorno1$U[,1], alvo)


criterio = c()
for(i in 1:100){
	retorno1 = IEM(x, tipoMatriz = 3, numRep = 1, plot = FALSE)
	roc = ROC(test = retorno1$Posteriori[,1], stat = alvo, PV = FALSE, MX = FALSE, MI = FALSE, )
	if(roc$AUC < 0.5){
		roc = ROC(test = retorno1$Posteriori[,2], stat = alvo, PV = FALSE, MX = FALSE, MI = FALSE, )
	}
	criterio = c(criterio, roc$AUC)	
}

criterio = c()
for(i in 1:100){
	retorno1 = IKFCM(x, c = 2, m = 0.1, maxRep = 20, monteCarlo = 10)		
	roc = ROC(test = retorno1$U[,1], stat = alvo, PV = FALSE, MX = FALSE, MI = FALSE, )
	if(roc$AUC < 0.5){
		roc = ROC(test = retorno1$U[,2], stat = alvo, PV = FALSE, MX = FALSE, MI = FALSE, )
	}
	criterio = c(criterio, roc$AUC)	
}

criterio = c()
for(i in 1:100){
	tryCatch({
		retorno1 = IKEM(x, maxRep = 1, kernelType = 0, D = 5, plot = FALSE, numRep = 10)
		roc = ROC(test = retorno1$Posteriori[,1], stat = alvo, PV = FALSE, MX = FALSE, MI = FALSE, )
		if(roc$AUC < 0.5){
			roc = ROC(test = retorno1$Posteriori[,2], stat = alvo, PV = FALSE, MX = FALSE, MI = FALSE, )
		}
		criterio = c(criterio, roc$AUC)	
	}, error=function(e){})
}

###########################################################################################################


x = as.matrix(read.table('carros.txt', head = TRUE))
x = normalizarDadosSimbolicos(x)

alvo = rep(1, dim(x)[1])
alvo[c(4,11,15,22,23)] = 2

x_interval = interval(x[,1],x[,2])
y_interval = interval(x[,3],x[,4])	
z_interval = interval(x[,5],x[,6])	
#intervalGraph3D(x_interval, y_interval, z_interval, alvo = alvo)

alvo[c(4,11,15,22,23)] = 0

criterio = c()
for(i in 1:100){
	retorno1 = IEM(x, tipoMatriz = 3, numRep = 10, plot = FALSE)
	roc = ROC(test = retorno1$Posteriori[,1], stat = alvo, PV = FALSE, MX = FALSE, MI = FALSE, )
	if(roc$AUC < 0.5){
		roc = ROC(test = retorno1$Posteriori[,2], stat = alvo, PV = FALSE, MX = FALSE, MI = FALSE, )
	}
	criterio = c(criterio, roc$AUC)	
}

criterio = c()
for(i in 1:100){
	retorno1 = IKFCM(x, c = 2, m = 0.1, maxRep = 20, monteCarlo = 10)		
	roc = ROC(test = retorno1$U[,1], stat = alvo, PV = FALSE, MX = FALSE, MI = FALSE, )
	if(roc$AUC < 0.5){
		roc = ROC(test = retorno1$U[,2], stat = alvo, PV = FALSE, MX = FALSE, MI = FALSE, )
	}
	criterio = c(criterio, roc$AUC)	
}

criterio = c()
for(i in 1:100){
	tryCatch({
		retorno1 = IKEM(x, maxRep = 50, kernelType = 0, D = 6, plot = FALSE, numRep = 10)
		roc = ROC(test = retorno1$Posteriori[,1], stat = alvo, PV = FALSE, MX = FALSE, MI = FALSE, )
		if(roc$AUC < 0.5){
			roc = ROC(test = retorno1$Posteriori[,2], stat = alvo, PV = FALSE, MX = FALSE, MI = FALSE, )
		}
		criterio = c(criterio, roc$AUC)	
	}, error=function(e){})
}