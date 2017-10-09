setwd("C:/Users/DanielBion/Documents/Mestrado/IKEM")
#options(digits = 20)
source("functions_sda.r")
library('MASS')
library('Epi')
source("Interval_EM.r")
source("Interval_Kernel_EM.r")
source('Interval_Kernel_Fuzzy_C-Means.r')

x = as.matrix(read.table('data.txt', head = TRUE, sep = ","))	
x = normalizarDadosSimbolicos(x)

############################################################################################################################################################
x = simulatedData1()
retorno = IEM(x, tipoMatriz = 3, numRep = 1, plot = TRUE)

retorno2 = IKEM(x, maxRep = 1, kernelType = 0, D = 13, plot = TRUE)

############################################################################################################################################################
	
#Simulação 1
dados = as.matrix(read.table('Target.txt', head = TRUE))	
x = matrix(0, dim(dados)[1], dim(dados)[2] * 2)
x[,c(1,2)] = dados[,1]
x[,c(3,4)] = dados[,2]
x[,2] = x[,2] + runif(length(x[,2]), 0.1, 0.3)
x[,4] = x[,4] + runif(length(x[,2]), 0.1, 0.3)
x = x[1:200,]
x = normalizarDadosSimbolicos(x)
x_interval = interval(x[,1],x[,2])
y_interval = interval(x[,3],x[,4])
intervalGraph2D(x_interval, y_interval)

alvo = rep(1, 200)
alvo[137:200] = 0

intervalGraph2D(x_interval, y_interval, alvo)
retorno2 = IKEM(x, maxRep = 1, kernelType = 0, D = 13, plot = TRUE)

ROC(test = retorno2$Posteriori[,1], stat=alvo)

############################################################################################################################################################

#Simulação 2
dados = as.matrix(read.table('banana.txt', head = TRUE))	
x = matrix(0, dim(dados)[1], dim(dados)[2] * 2)
x[,c(1,2)] = dados[,1]
x[,c(3,4)] = dados[,2]
x[,2] = x[,2] + runif(length(x[,2]), 0.1, 1)
x[,4] = x[,4] + runif(length(x[,2]), 0.1, 1)
x = normalizarDadosSimbolicos(x)
x_interval = interval(x[,1],x[,2])
y_interval = interval(x[,3],x[,4])
intervalGraph2D(x_interval, y_interval)

alvo = rep(1, 300)
alvo[200:300] = 2
intervalGraph2D(x_interval, y_interval, alvo)


alvo[200:300] = 0

idx = c();
for(i in 200:300){
	if(dist(rbind(c(x[i,1], x[i,3]), c(0.5, 0.5))) < 0.3){
		idx = c(idx, i)
	}
}
alvo[idx] = 1

retorno2 = IKEM(x, maxRep = 1, kernelType = 0, D = 13, plot = TRUE)
ROC(test = retorno2$Posteriori[,1], stat=alvo)

############################################################################################################################################################

retorno = IKFCM(x, c = 2, m = 0.1, maxRep = 5, monteCarlo = 1)

############################################################################################################################################################

#Simulação 3

dados = as.matrix(read.table('Target.txt', head = TRUE))	
x = matrix(0, dim(dados)[1], dim(dados)[2] * 2)
x[,c(1,2)] = dados[,1]
x[,c(3,4)] = dados[,2]
x[,2] = x[,2] + runif(length(x[,2]), 0.1, 0.3)
x[,4] = x[,4] + runif(length(x[,2]), 0.1, 0.3)
x = normalizarDadosSimbolicos(x)

x = x[1:200,]
y = x + 0.6;
x = rbind(x,y)

x = normalizarDadosSimbolicos(x)

x_interval = interval(x[,1],x[,2])
y_interval = interval(x[,3],x[,4])
intervalGraph2D(x_interval, y_interval)

retorno2 = IKEM(x, maxRep = 1, kernelType = 0, D = 15, plot = TRUE)

############################################################################################################################################################