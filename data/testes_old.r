
dados1 = function(){
	xoffset = 12
	yoffset = 4
	
	n = 100
	xa = runif(n, min = 5, max = 15)	
	ya = sin(xa*pi*2) * 3
	xa = cos(xa*pi*2) * 3
		
	yb = sin(xa*pi*2) * 3
	xb = cos(xa*pi*2) * 3 + xoffset
	
	yc = sin(xa*pi*2) * 2 + yoffset
	xc = cos(xa*pi*2) * 1.5 + xoffset/2
	
	yd = sin(xa*pi*2) * 2 - yoffset
	xd = cos(xa*pi*2) * 1.5 + xoffset/2
	
	xe = c(xa, xb, xc, xd)
	ye = c(ya, yb, yc ,yd)
	
	plot(xe, ye)
	
	xlimit = c(1,11)
	ylimit = c(-4.5,4.5)
		
	for(i in 1:length(xe)){
		if(xe[i] > xlimit[2] || xe[i] < xlimit[1]){
			xe[i] = -100
			ye[i] = -100			
		}
		if(ye[i] > ylimit[2] || ye[i] < ylimit[1]){
			xe[i] = -100
			ye[i] = -100			
		}
	}
	xe = xe[xe != -100]
	ye = ye[ye != -100]
	
	# xf = cbind(xe, ye)
		
	xr = runif(length(xe), min = 0.05, max = 0.5)
	yr = runif(length(xe), min = 0.04, max = 0.6)
		
	x = cbind(xe - xr/2, xe + xr/2, ye - yr/2, ye + yr/2)
	
	x = normalizarDadosSimbolicos(x)
	
	x_interval = interval(x[,1],x[,2])
	y_interval = interval(x[,3],x[,4])	
	intervalGraph2D(x_interval, y_interval)
	dim(x)
	return (x)
}



test2 = function(){
	xoffset = 8
	yoffset = 4
	
	n = 80	
	xa = runif(n, min = 5, max = 15)	
	ya = sin(xa*pi*2) * 2
	xa = cos(xa*pi*2) * 2
		
	yb = sin(xa*pi*2) * 2
	xb = cos(xa*pi*2) * 2 + xoffset
	
	yc = sin(xa*pi*2) * 2 + yoffset
	xc = cos(xa*pi*2) * 2 + xoffset/2
	
	yd = sin(xa*pi*2) * 2 - yoffset
	xd = cos(xa*pi*2) * 2 + xoffset/2
	
	xe = rbind(xa, xb, xc, xd)
	ye = rbind(ya, yb, yc ,yd)
	
	# plot(xe, ye)
	
	xlimit = c(0.5,7.5)
	ylimit = c(-4,4)
		
	for(i in 1:length(xe)){
		if(xe[i] > xlimit[2] || xe[i] < xlimit[1]){
			xe[i] = -100
			ye[i] = -100			
		}
		if(ye[i] > ylimit[2] || ye[i] < ylimit[1]){
			xe[i] = -100
			ye[i] = -100			
		}
	}
	xe = xe[xe != -100]
	ye = ye[ye != -100]
	
	x = cbind(xe, ye)
	plot(xe, ye)
	
	x = normalizarDados(x)
}

simulatedData1 = function(){
	xoffset = 10
	yoffset = 4.5
	
	n = 100
	xa = runif(n, min = 5, max = 15)	
	ya = sin(xa*pi*2) * 3
	xa = cos(xa*pi*2) * 3
		
	yb = sin(xa*pi*2) * 3
	xb = cos(xa*pi*2) * 3 + xoffset
	
	yc = sin(xa*pi*2) * 2 + yoffset
	xc = cos(xa*pi*2) * 1.5 + xoffset/2
	
	yd = sin(xa*pi*2) * 2 - yoffset
	xd = cos(xa*pi*2) * 1.5 + xoffset/2
	
	xe = c(xa, xb, xc, xd)
	ye = c(ya, yb, yc ,yd)
	
	plot(xe, ye)
	
	xlimit = c(1,9)
	ylimit = c(-4.5,4.5)
		
	for(i in 1:length(xe)){
		if(xe[i] > xlimit[2] || xe[i] < xlimit[1]){
			xe[i] = -100
			ye[i] = -100			
		}
		if(ye[i] > ylimit[2] || ye[i] < ylimit[1]){
			xe[i] = -100
			ye[i] = -100			
		}
	}
	xe = xe[xe != -100]
	ye = ye[ye != -100]
	
	# xf = cbind(xe, ye)
		
	xr = runif(length(xe), min = 0.05, max = 0.5)
	yr = runif(length(xe), min = 0.04, max = 0.6)
		
	x = cbind(xe - xr/2, xe + xr/2, ye - yr/2, ye + yr/2)
	
	x = normalizarDadosSimbolicos(x)
	
	x_interval = interval(x[,1],x[,2])
	y_interval = interval(x[,3],x[,4])	
	intervalGraph2D(x_interval, y_interval)
	dim(x)
	return (x)
}

simulatedData2 = function(){
	xoffset = 12
	yoffset = 4
	
	n = 100
	xa = runif(n, min = 5, max = 15)	
	ya = sin(xa*pi*2) * 3
	xa = cos(xa*pi*2) * 3
		
	yb = sin(xa*pi*2) * 3
	xb = cos(xa*pi*2) * 3 + xoffset
	
	yc = sin(xa*pi*2) * 2 + yoffset
	xc = cos(xa*pi*2) * 1.5 + xoffset/2
	
	yd = sin(xa*pi*2) * 2 - yoffset
	xd = cos(xa*pi*2) * 1.5 + xoffset/2
	
	xe = c(xa, xb, xc, xd)
	ye = c(ya, yb, yc ,yd)
	
	plot(xe, ye)
	
	xlimit = c(1,11)
	ylimit = c(-4.5,4.5)
		
	for(i in 1:length(xe)){
		if(xe[i] > xlimit[2] || xe[i] < xlimit[1]){
			xe[i] = -100
			ye[i] = -100			
		}
		if(ye[i] > ylimit[2] || ye[i] < ylimit[1]){
			xe[i] = -100
			ye[i] = -100			
		}
	}
	xe = xe[xe != -100]
	ye = ye[ye != -100]
	
	# xf = cbind(xe, ye)
		
	xr = runif(length(xe), min = 0.05, max = 0.5)
	yr = runif(length(xe), min = 0.04, max = 0.6)
		
	x = cbind(xe - xr/2, xe + xr/2, ye - yr/2, ye + yr/2)
	
	x = normalizarDadosSimbolicos(x)
	
	x_interval = interval(x[,1],x[,2])
	y_interval = interval(x[,3],x[,4])	
	intervalGraph2D(x_interval, y_interval)
	dim(x)
	return (x)
}



test = function(n){	
	n = 400
	noutliers = 100
	E <- rnorm(n, mean=0.01, sd=0.5)
	xc = runif(n, min = 5, max = 15)
		
	yc = sin(xc*pi*2) * 2 + E
	xc = cos(xc*pi*2) * 2
	
	noisexc = sin(xc*pi*2) * 4 + E*2
	noisexc = noisexc[1:noutliers]
	noiseyc = cos(xc*pi*2) * 4
	noiseyc = noiseyc[1:noutliers]
		
	for(i in 1:length(xc)){
		if(xc[i] > -0.5 && yc[i] < 0.5){
			xc[i] = -100
			yc[i] = -100			
		}
	}
	
	xc = xc[xc != -100]
	yc = yc[yc != -100]
		
	# noisexc = runif(noutliers, min = -5, max = 5)
	# noiseyc = runif(noutliers, min = -5, max = 5)
		
	xc = c(xc, noisexc)
	yc = c(yc, noiseyc)
		
	xr = runif(length(xc), min = 0.05, max = 0.5)
	yr = runif(length(xc), min = 0.04, max = 0.6)
		
	x = cbind(xc - xr/2, xc + xr/2, yc - yr/2, yc + yr/2)
	
	x_interval = interval(x[,1],x[,2])
	y_interval = interval(x[,3],x[,4])
	
	intervalGraph2D(x_interval, y_interval)
}

geraDadosSenoCenario1 = function(n){	
	n = 100
	E <- rnorm(n, mean=0.01, sd=0.03)
	xc = runif(n, min = 0, max = 0.75)
	xr = runif(n, min = 0.05, max = 0.5)
	yr = runif(n, min = 0.04, max = 0.6)
	yc = sin(xc*pi) + 10 + E
	
	o = order(xc, decreasing = FALSE)
	xc = xc[o]
	xr = xr[o]
	yc = yc[o]
	yr = yr[o]
	
	sdYC = sd(yc);
	idx = floor(runif(10, 0, n))
	for(i in 1:length(idx)){
		yc[idx[i]] = yc[idx[i]] - 6 * sdYC		
	}
		
	x = cbind(xc - xr/2, xc + xr/2, yc - yr/2, yc + yr/2)
	
	x_interval = interval(x[,1],x[,2])
	y_interval = interval(x[,3],x[,4])
	intervalGraph2D(x_interval, y_interval)
}

geraDadosSenoCenario2 = function(n){	
	n = 500
	E <- rnorm(n, mean=0.01, sd=0.03)
	xc = runif(n, min = 0, max = 0.75)
	xr = runif(n, min = 0.05, max = 0.5)
	yr = runif(n, min = 0.04, max = 0.6)
	yc = sin(xc*pi) + 10 + E
	
	o = order(xc, decreasing = FALSE)
	xc = xc[o]
	xr = xr[o]
	yc = yc[o]
	yr = yr[o]
	
	sdYC = sd(yc);
	idx = floor(runif(10, 0, n))
	for(i in 1:length(idx)){
		if(idx[i] < n/2)
		{
			yc[idx[i]] = yc[idx[i]] + 6 * sdYC
		}
		else
		{
			yc[idx[i]] = yc[idx[i]] - 6 * sdYC
		}
	}
		
	x = cbind(xc - xr/2, xc + xr/2, yc - yr/2, yc + yr/2)
	
	x_interval = interval(x[,1],x[,2])
	y_interval = interval(x[,3],x[,4])
	intervalGraph2D(x_interval, y_interval)
}

geraDadosSenoCenario3 = function(n){	
	n = 100
	E <- rnorm(n, mean=0.01, sd=0.03)
	xc = runif(n, min = 0, max = 0.75)
	xr = runif(n, min = 0.05, max = 0.5)
	yr = runif(n, min = 0.04, max = 0.6)
	yc = sin(xc*pi) + 10 + E
	
	o = order(xc, decreasing = FALSE)
	xc = xc[o]
	xr = xr[o]
	yc = yc[o]
	yr = yr[o]
	
	sdYR = sd(yr);
	idx = floor(runif(10, 0, n))
	for(i in 1:length(idx)){
		yr[idx[i]] = yr[idx[i]] + 12 * sdYR		
		
	}
		
	x = cbind(xc - xr/2, xc + xr/2, yc - yr/2, yc + yr/2)
	
	x_interval = interval(x[,1],x[,2])
	y_interval = interval(x[,3],x[,4])
	intervalGraph2D(x_interval, y_interval)
}

geraDadosSenoCenario4 = function(n){	
	n = 100
	E <- rnorm(n, mean=0.01, sd=0.03)
	xc = runif(n, min = 0, max = 0.75)
	xr = runif(n, min = 0.05, max = 0.5)
	yr = runif(n, min = 0.04, max = 0.6)
	yc = sin(xc*pi) + 10 + E
	
	o = order(xc, decreasing = FALSE)
	xc = xc[o]
	xr = xr[o]
	yc = yc[o]
	yr = yr[o]
	
	sdYC = sd(yc);
	sdYR = sd(yr);
	idxC = floor(runif(10, 0, n))
	idxC = c(10,20,30,40,50,60,70,80,90,100)
	idxR = sample(idxC, 5)
	
	for(i in 1:length(idxC)){
		random = runif(1);
		if(random > 0.5){
			yc[idxC[i]] = yc[idxC[i]] - 4 * sdYC		
		}else{
			yc[idxC[i]] = yc[idxC[i]] + 4 * sdYC
		}		
	}
	for(i in 1:length(idxR)){
		yr[idxR[i]] = yr[idxR[i]] + 12 * sdYR
	}

	x = cbind(xc - xr/2, xc + xr/2, yc - yr/2, yc + yr/2)
	
	x_interval = interval(x[,1],x[,2])
	y_interval = interval(x[,3],x[,4])
	intervalGraph2D(x_interval, y_interval)
}


geraDados = function(n, dimensao){	
	#Definir que o erro E é obtido de uma distribuiç˜ao normal com média 0 e desvio padr˜ao 1.
	E <- rnorm(n, mean=0, sd=1)
	
	beta0 = runif(1,min=-5,max=5)
	beta = runif(dimensao-1,min=-5,max=5)
	
	#Definir que o preditor de ponto xij (j = 1, 2) é gerado de uma distribuição uniforme com valores mínimo e máximo: 20 e 40 respectivamente
	X = matrix(0, nrow = n, ncol = dimensao-1)

	#Inicialização dos vetores de Centro e Range
	C=matrix(0, nrow = n, ncol = dimensao)
	R=matrix(0, nrow = n, ncol = dimensao)
	
	mu = runif(dimensao - 1,min=20,max=40)
	
	mcov = diag(dimensao-1)		
	diag(mcov) = 9
	
	p1 = 0.8
	for(i in 1:dimensao - 1){
		for(j in 1:dimensao - 1){
			if(i!=j){
				mcov[i,j] = p1*3*3
			}
		}
	}
	
	X = mvrnorm(n, mu, mcov)
	
	#Computar uma resposta de ponto yi = b0 + xi1*b1 + xi2*b2 + Ei.
	Y = beta0 + E
	
	for (p in 1:(dimensao-1)){
	#	X[,p] = runif(n,min=20,max=40)
		Y = Y + beta[p]*X[,p]
	}		
	
	#mu = cbind(X,Y)
	#mcov = diag(dimensao)		
	#diag(mcov) = 9
	
	#for (i in 1:n){
		#Obter uma amostra de 50 pontos no R3, de acordo com uma distribuiç˜ao normal
		#multivariada com vetor de média µ = (xi1, xi2, yi) e a matriz diagonal de covariância
		#com sjj = 9 (j = 1, 2, 3).
	#	INTERVAL = mvrnorm(n = 50, mu[i,], mcov)
		#min_X = apply(INTERVAL, 2, min)
		#max_X = apply(INTERVAL, 2, max)
	#	min_X = apply(INTERVAL, 2, quantile)[1,]
    #   max_X = apply(INTERVAL, 2, quantile)[4,] #3 Quartil
	#	C[i,] = (min_X + max_X)/2
	#	R[i,] = (max_X - min_X)#runif(1, 1, 10)
	#}
	
	C = cbind(X,Y)
	R = runif(n, 1, 10)
	for(i in 1:p-1){
		R = cbind(R, runif(n, 1, 10))
	}
	
	DataCenterRange = c()
	#Dados com centro e range
	for(p in 1:(dimensao)){
		DataCenterRange = cbind(DataCenterRange, C[,p],R[,p])
	}		
	return(DataCenterRange)
}


geraDadosOld2 = function(n, dimensao){	
	#Definir que o erro E é obtido de uma distribuiç˜ao normal com média 0 e desvio padr˜ao 1.
	E <- rnorm(n, mean=0, sd=1)
	
	beta0 = runif(1,min=-10,max=10)
	beta = runif(dimensao-1,min=-10,max=10)

	#B0 = beta[1]
	#B1 = beta[2]
	#B2 = beta[3]

	#Definir que o preditor de ponto xij (j = 1, 2) é gerado de uma distribuição uniforme com valores mínimo e máximo: 20 e 40 respectivamente
	X = matrix(0, nrow = n, ncol = dimensao-1)

	#Inicialização dos vetores de Centro e Range
	XC=matrix(0, nrow = n, ncol = dimensao-1)
	XR=matrix(0, nrow = n, ncol = dimensao-1)
	YC=c(rep(0,each=n))	
	YR=c(rep(0,each=n))

	#Computar uma resposta de ponto yi = b0 + xi1*b1 + xi2*b2 + Ei.
	Y = beta0 + E
	for (p in 1:(dimensao-1)){
		X[,p] = runif(n,min=20,max=40)
		Y = Y + beta[p]*X[,p]
	}
		
	for (i in 1:n){
		#Obter uma amostra de 50 pontos no R3, de acordo com uma distribuiç˜ao normal
		#multivariada com vetor de média µ = (xi1, xi2, yi) e a matriz diagonal de covariância
		#com sjj = 9 (j = 1, 2, 3).		
		
		for(p in 1:(dimensao-1)){			
			INTERVAL_X = rnorm (50,X[i,p],3)#sd(X)) 			
			min_X = quantile(INTERVAL_X)[1]
			max_X = quantile(INTERVAL_X)[5]
			XC[i,p] = (min_X + max_X)/2
			XR[i,p] = (max_X - min_X)
		}
				
		INTERVAL_Y = rnorm (50,Y[i],3)#sd(Y))
		min_Y = quantile(INTERVAL_Y)[1]
		max_Y = quantile(INTERVAL_Y)[5]		
		YC[i] = (min_Y + max_Y)/2
		YR[i] = (max_Y - min_Y)	
	}
	DataCenterRange = c()
	#Dados com centro e range
	for(p in 1:(dimensao-1)){
		DataCenterRange = cbind(DataCenterRange, XC[,p],XR[,p])
	}		
	DataCenterRange = cbind(DataCenterRange,YC,YR)	
}

geraOutlier = function (DataCenterRange, cenario, quantidade, deslocamentoC, deslocamentoR, retorno_min_max)
{
	YCIdx = dim(DataCenterRange)[2] - 1
	YRIdx = dim(DataCenterRange)[2]
	
	SDYC <- sd(DataCenterRange[,YCIdx])
	SDYR <- sd(DataCenterRange[,YRIdx])

	o = order(DataCenterRange[,YCIdx],decreasing = TRUE)
			
	for(i in 1:dim(DataCenterRange)[2]){
		DataCenterRange[,i] = DataCenterRange[,i][o]
	}
	
	Cenario = DataCenterRange
	
	if(quantidade > 0){
		outliers_idx = 1:quantidade
	
		#Cenario Center-Range      
		if(cenario == 1){
			Cenario[outliers_idx,YCIdx] = Cenario[outliers_idx,YCIdx] + (deslocamentoC*SDYC)
			
		} else if(cenario == 2){
			Cenario[outliers_idx[1:floor(quantidade/2)],YCIdx] = Cenario[outliers_idx[1:floor(quantidade/2)],YCIdx] + (deslocamentoC*SDYC)
			Cenario[outliers_idx[(ceiling(quantidade/2)+1):quantidade],YCIdx] = Cenario[(ceiling(quantidade/2)+1):quantidade,YCIdx] - (deslocamentoC*SDYC)
			
		} else if(cenario == 3){
			Cenario[outliers_idx,YRIdx] = Cenario[outliers_idx,YRIdx] * (deslocamentoR*SDYR)
			
		} else if(cenario == 4){
			Cenario[outliers_idx,YCIdx] = Cenario[outliers_idx,YCIdx] + (deslocamentoC*SDYC)
			#30% dos outliers de centro são escolhidos aleatoriamente para serem outliers de range
			idx = sample(outliers_idx, 30*length(outliers_idx)/100)
			Cenario[idx,YRIdx] = Cenario[idx,YRIdx] * (deslocamentoR*SDYR)		
			
		} else if(cenario == 5){
			#70% dos outliers são escolhidos aleatoriamente para serem outliers de centro
			idx70 = sample(outliers_idx, round(70*length(outliers_idx)/100))
			Cenario[idx70,YCIdx] = Cenario[idx70,YCIdx] + (deslocamentoC*SDYC)	
			#30% dos outliers são escolhidos aleatoriamente para serem outliers de range		
			idx30 = outliers_idx[-idx70]
			# idx30 = sample(outliers_idx, round(30*length(outliers_idx)/100))
			Cenario[idx30,YRIdx] = Cenario[idx30,YRIdx] * (deslocamentoR*SDYR)		
		}
	}
	#Cenário3 Min-Max
	Cenario_Min_Max = c()
	cont = 1
	
	for(i in seq(from = 1, to = (dim(DataCenterRange)[2] - 2), by = 2)){
		XMin = Cenario[,i] - ((Cenario[,i+1])/2)
		XMax = Cenario[,i] + ((Cenario[,i+1])/2)
		Cenario_Min_Max = cbind(Cenario_Min_Max, XMin, XMax)
		colnames(Cenario_Min_Max)[i] = paste0('XMin',cont)
		colnames(Cenario_Min_Max)[i+1] = paste0('XMax',cont)
		cont = cont + 1
	}
	
	YMin = Cenario[,YCIdx] - (Cenario[,YRIdx])/2	
	YMax = Cenario[,YCIdx] + (Cenario[,YRIdx])/2
		
	Cenario_Min_Max = cbind(Cenario_Min_Max,YMin,YMax)
	
	if(retorno_min_max){
		return(Cenario_Min_Max)
	}else{
		return(Cenario)
	}
}

#Só funcionava para 2 dimensões
geraDados_Old = function(quantidade){	
	n = quantidade

	#Definir que o erro E é obtido de uma distribuiç˜ao normal com média 0 e desvio padr˜ao 1.
	E <- rnorm(n, mean=0, sd=1)
	
	beta = runif(2,min=-10,max=10)

	B0 = beta[1]
	B1 = beta[2]

	#Definir que o preditor de ponto xij (j = 1, 2) é gerado de uma distribuição uniforme com valores mínimo e máximo: 20 e 40 respectivamente
	X = runif(n,min=20,max=40)

	#Inicialização dos vetores de Centro e Range
	XC=c(rep(0,each=n))
	YC=c(rep(0,each=n))
	XR=c(rep(0,each=n))
	YR=c(rep(0,each=n))

	#Computar uma resposta de ponto yi = b0 + xi1*b1 + xi2*b2 + Ei.
	Y = B0 + B1*X + E

	for (j in 1:n){
		#Obter uma amostra de 50 pontos no R3, de acordo com uma distribuiç˜ao normal
		#multivariada com vetor de média µ = (xi1, xi2, yi) e a matriz diagonal de covariância
		#com sjj = 9 (j = 1, 2, 3).		
		INTERVAL_X = rnorm (50,X[j],3) 
		INTERVAL_Y = rnorm (50,Y[j],3)
		
		#Computar os intervalos 
		min_X = quantile(INTERVAL_X)[1] #min(INTERVAL_X)
		max_X = quantile(INTERVAL_X)[5] #max(INTERVAL_X)
		min_Y = quantile(INTERVAL_Y)[1]
		max_Y = quantile(INTERVAL_Y)[5]

		XC[j] = (min_X + max_X)/2
		YC[j] = (min_Y + max_Y)/2
		XR[j] = (max_X - min_X)
		YR[j] = (max_Y - min_Y)	
	}
	#Dados com centro e range
	Data_Center_Range = cbind(XC,YC,XR,YR)	
}

#Só funcionava para 2 dimensões
geraOutlier_Old = function (Data_Center_Range, cenario, quantidade, deslocamento, retorno_min_max)
{
	SDYC <- sd(Data_Center_Range[,2])
	SDYR <- sd(Data_Center_Range[,4])

	o = order(Data_Center_Range[,2],decreasing = TRUE)
	Data_Center_Range[,1] = Data_Center_Range[,1][o]
	Data_Center_Range[,2] = Data_Center_Range[,2][o]
	Data_Center_Range[,3] = Data_Center_Range[,3][o]
	Data_Center_Range[,4] = Data_Center_Range[,4][o]	

	#temp = cbind(0,rep(0,each=quantidade),0,0)
	#Cenario = rbind(Data_Center_Range,temp) 
	Cenario = Data_Center_Range
	
	#size = dim(Data_Center_Range)[1]
	
	outliers_idx = 1:quantidade
	#outliers_idx = (size+1):(size+quantidade)
	#Cenario[outliers_idx,] = Data_Center_Range[1:quantidade,] 

	#Cenario Center-Range      
	if(cenario == 1){
		Cenario[outliers_idx,2] = Cenario[outliers_idx,2] + (deslocamento*SDYC)
		
	} else if(cenario == 2){
		Cenario[outliers_idx[1:floor(quantidade/2)],2] = Cenario[outliers_idx[1:floor(quantidade/2)],2] + (deslocamento*SDYC)
		Cenario[outliers_idx[ceiling(quantidade/2):quantidade],2] = Cenario[ceiling(quantidade/2):quantidade,2] - (deslocamento*SDYC)
		
	} else if(cenario == 3){
		Cenario[outliers_idx,4] = Cenario[outliers_idx,4] * (deslocamento*SDYR)
		
	} else if(cenario == 4){
		Cenario[outliers_idx,2] = Cenario[outliers_idx,2] + (deslocamento*SDYC)
		#30% dos outliers de centro são escolhidos aleatoriamente para serem outliers de range
		idx = sample(outliers_idx, 30*length(outliers_idx)/100)
		Cenario[idx,4] = Cenario[idx,4] * (deslocamento*SDYR)		
		
	} else if(cenario == 5){
		#70% dos outliers são escolhidos aleatoriamente para serem outliers de centro
		idx70 = sample(outliers_idx, 70*length(outliers_idx)/100)
		Cenario[idx70,2] = Cenario[idx70,2] + (deslocamento*SDYC)	
		#30% dos outliers são escolhidos aleatoriamente para serem outliers de range		
		idx30 = sample(outliers_idx, 30*length(outliers_idx)/100)
		Cenario[idx30,4] = Cenario[idx30,4] * (deslocamento*SDYR)		
	}

	#Cenário3 Min-Max
	XMin = Cenario[,1] - (Cenario[,3])/2
	YMin = Cenario[,2] - (Cenario[,4])/2
	XMax = Cenario[,1] + (Cenario[,3])/2
	YMax = Cenario[,2] + (Cenario[,4])/2
		
	Cenario_Min_Max = cbind(XMin,XMax,YMin,YMax)
	
	if(retorno_min_max){
		return(Cenario_Min_Max)
	}else{
		return(Cenario)
	}
}
