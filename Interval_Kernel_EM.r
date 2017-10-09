IKEM = function(x, maxRep, kernelType, D, plot, initialPosteriori,numRep){	
	IRBF = function(a_inf, a_sup, b_inf, b_sup, sigma){
		inf = exp(-sum((a_inf - b_inf)^2) / (2*sigma^2))
		sup = exp(-sum((a_sup - b_sup)^2) / (2*sigma^2))
		return (inf + sup)
	}
	
	IPF = function(a_inf, a_sup, b_inf, b_sup, d){
		inf = sum((a_inf * b_inf)^d)
		sup = sum((a_sup * b_sup)^d)
		return (inf + sup)
	}

	hausDist = function(xa, xb, H){
		infSeq =  seq(from = 1, to = pp, by = 2)
					
		soma = 0
		for(j in infSeq){
			difInf = abs(xa[j] - xb[j])
			difSup = abs(xa[j+1] - xb[j+1])
			maxDif = max(difInf, difSup)
			soma = soma + (maxDif/H[j])^2
		}		
		soma = soma^0.5
	}
	
	inicializarPosteriori = function(){		
		guess=sample(5:(n-5),2,replace=FALSE)
		
		if (guess[1] < guess[2]){
			guess[1]=guess[1] - 4
			guess[2]=guess[2] + 4
		}else{
			guess[2]=guess[2] - 4
			guess[1]=guess[1] + 4
		}
		
		#guess = c(37,47)
		#print(guess)
		#Escolher duas observações aleatoriamente
		g1 = x[guess[1],]       
		g2 = x[guess[2],]

		alvo = rep(0, n)

		#Definir o alvo inicial baseado na distância entre os dados e as duas observações escolhidas
		for(i in 1:n){
			d1 = 0
			d2 = 0
			for(j in 1:pp){
				d1 = d1 + (x[i,j] - g1[j])^2
				d2 = d2 + (x[i,j] - g2[j])^2
			}			
			if (d1 < d2){
				alvo[i]=1 
			}else{ 
				alvo[i]=2 
			}
		}
		
		infSeq =  seq(from = 1, to = pp, by = 2)
		
		H = c()		
		for(j in infSeq){
			soma = 0
			for(i in 1:n){
				for(h in 1:n){
					difInf = abs(x[i, j] - x[h, j])
					difSup = abs(x[i, j+1] - x[h, j+1])
					maxDif = max(difInf, difSup)
					soma = soma + maxDif^2					
				}
			}
			soma = soma/(2*(n^2))
			H[j] = soma^0.5
		}
		
		for(i in 1:n){
			d1 = hausDist(x[i,], g1, H)
			d2 = hausDist(x[i,], g2, H)
			
			if (d1 < d2){
				alvo[i]=1 
			}else{ 
				alvo[i]=2 
			}
		}
		
		theta = list(
			tau=c(0.5,0.5),
			media = matrix(			 
			  c(apply(x[alvo==1,1:pp],2,mean),
			    apply(x[alvo==2,1:pp],2,mean)),
			    nrow = c, ncol = pp, byrow=TRUE),	
			mcov = list(matrix(0,nrow=p,ncol=p,byrow=TRUE),
						matrix(0,nrow=p,ncol=p,byrow=TRUE))
		)		
				
		mcovA = matrix(0,nrow=p,ncol=p,byrow=TRUE)
		mcovB = matrix(0,nrow=p,ncol=p,byrow=TRUE)
		for(i in 1:p){
			infIdxi = (i*2) - 1
			supIdxi = (i*2)
			for(j in 1:p){
				infIdxj = (j*2) - 1
				supIdxj = (j*2)
				mcovA[i,j] = (cov(x[alvo==1, infIdxi], x[alvo==1, infIdxj]) + cov(x[alvo==1, supIdxi], x[alvo==1, supIdxj]))/2 
				mcovB[i,j] = (cov(x[alvo==2, infIdxi], x[alvo==2, infIdxj]) + cov(x[alvo==2, supIdxi], x[alvo==2, supIdxj]))/2 
			}			
		}		
				
		tipoMatriz = 3
		if(tipoMatriz == 1){
			# Matriz Identidade
			mcovA = diag(p)
			mcovB = diag(p)			
		}else if(tipoMatriz == 2){
			# Matriz diagonal
			for(i in 1:p){				
				for(j in 1:p){				
					if(i != j){
						mcovA[i,j] = 0
						mcovB[i,j] = 0
					}
				}
			}			
		}
		
		theta$mcov[[1]] = mcovA
		theta$mcov[[2]] = mcovB
					
		for (i in 1:n)
		{
			A = c()
			B = c()						
			
			for(k in 1:c){
				XMinf = c()
				XMsup = c()
				for(j in 1:p){
					infIdx = (2*j) - 1
					supIdx = 2*j
					XMinf = cbind(XMinf, x[i,infIdx] - theta$media[k,infIdx])	
					XMsup = cbind(XMsup, x[i,supIdx] - theta$media[k,supIdx])	
				}
				A[k] = XMinf %*% solve(theta$mcov[[k]]) %*% t(XMinf)
				B[k] = XMsup %*% solve(theta$mcov[[k]]) %*% t(XMsup)
				#A[k] = XMinf %*% ginv(theta$mcov[[k]]) %*% t(XMinf)
				#B[k] = XMsup %*% ginv(theta$mcov[[k]]) %*% t(XMsup)
			}			
					
			numA = theta$tau[1] * (exp(-1/2 * (A[1]+B[1])/2) / sqrt(((2*pi)^p) * det(theta$mcov[[1]])))	
			numB = theta$tau[2] * (exp(-1/2 * (A[2]+B[2])/2) / sqrt(((2*pi)^p) * det(theta$mcov[[2]])))
			den = numA + numB
			
			if(is.nan(den)){
				print('DEN == 0')
			}
			
			posteriori[i, 1] = numA/den
			posteriori[i, 2] = numB/den			
		}
		return(posteriori)
	}
	
	inicializarK = function(){
		kernelInicial = matrix(0, nrow=n,ncol=n)
		infSeq =  seq(from = 1, to = pp, by = 2)
		for(i in 1:n){
			for(j in 1:n){
				a_inf = x[i,infSeq]
				a_sup = x[i,infSeq+1]
				b_inf = x[j,infSeq]
				b_sup = x[j,infSeq+1]
				
				if(kernelType == 1){
					kernelInicial[i,j] = IPF(a_inf, a_sup, b_inf, b_sup, 2.0)
				}
				else{
					kernelInicial[i,j] = IRBF(a_inf, a_sup, b_inf, b_sup, 1.0)					
				}
			}
		}
		return (kernelInicial)
	}
	
	calcMatrixes = function(){
		kernelM = list(matrix(0,nrow=n,ncol=n),
				  matrix(0,nrow=n,ncol=n))
				  
		pkernelM = kernelM
		ckernelM = kernelM				  
		cpkernelM = kernelM
		W = kernelM
		Wlinha = kernelM
		wAux = as.matrix(rep(1, n))
		
		w = posteriori
		for(i in 1:n)
		{
			for(k in 1:c)
			{
				w[i,k] = sqrt(posteriori[i,k]/sum(posteriori[,k]))
			}
		}
				
		for(k in 1:c){				
			for(i in 1:n){
				for(j in 1:n){
					kernelM[[k]][i,j] = w[i,k] * w[j,k] * kernelInicial[i,j]
					pkernelM[[k]][i,j] = w[j,k] * kernelInicial[i,j]
				}
			}
		} 
		
		for(k in 1:c){				
			W[[k]] = w[,k] %*% t(w[,k])
			Wlinha[[k]] = wAux %*% t(w[,k])
		}
		for(k in 1:c){		
			ckernelM[[k]] = kernelM[[k]] - (W[[k]] %*% kernelM[[k]]) - (kernelM[[k]] %*% W[[k]]) + (W[[k]] %*% kernelM[[k]] %*% W[[k]])
			cpkernelM[[k]] = pkernelM[[k]] - (Wlinha[[k]] %*% kernelM[[k]]) - (pkernelM[[k]] %*% W[[k]]) + (Wlinha[[k]] %*% kernelM[[k]] %*% W[[k]])
		}
		return(list(ckernelM, cpkernelM))
	}
		
	E.step = function(){				
		fdp = matrix(0, n, c)
		for(k in 1:c){			
			lambda = eigen(ckernelM[[k]])[1][[1]]		
			#lambda2 = lambda[(D+1):n]			 ##
			lambda = lambda[1:D]
			# if(length(which(lambda < 0)) > 0){
				# print('lambda negativo, escolha um D menor')
			# }
			beta = eigen(ckernelM[[k]])[2][[1]][,1:D]			
			
			for(i in 1:n){
				gama = cpkernelM[[k]][,i]
				soma = 0
				
				#p = sum(lambda2)/(n-D)##
				# p = 0.01##
				
				#for(j in 1:D){					##
				#	soma = soma +(t(as.matrix(beta[,j])) %*% as.matrix(gama))^2##
				#}##
				#infSeq =  seq(from = 1, to = pp, by = 2)##
				#error = 1 - soma##
				##
				#ortogonal = 1/((2*pi*p)^((n-D)/2)) * exp(-error/2*p)##
				
				
				#soma = 0 ###
				for(j in 1:D){					
					soma = soma + ((t(as.matrix(beta[,j])) %*% as.matrix(gama))^2)/ (lambda[j])
				}
				num = exp(-1/2 * soma)				
				den = ((2*pi)^(D/2)) * prod(lambda^0.5)
				#fdp[i,k] = num*ortogonal/den								
				fdp[i,k] = num/den								
			}
		}
	
		numA = tau[1] * fdp[,1]
		numB = tau[2] * fdp[,2]
		den = numA + numB
		
		posteriori[,1] = numA / den
		posteriori[,2] = numB / den
		
		return(posteriori)
	}


	criterio = function(){
		ML = 0.0		
		for(i in 1:n){
			ML = ML + log(max(posteriori[i,]))
		}
		ML = -ML
		return(ML)
	}
	
	n = dim(x)[1]
	pp = dim(x)[2]
	p = pp/2
	c=2	
    epsilon = 0.000001
	
	melhorML = 10000000000;  
	
	for(ii in 1:numRep){
		tau = c(0.5, 0.5)			
		t = 0	
		MLatual = 1	
        MLdepois = -1
		kernelInicial = inicializarK()

		posteriori = matrix(0, nrow = n, ncol = c)		
		posteriori[,1] = runif(n, 0, 1)
		posteriori[,2] = 1 - posteriori[,1]
		posteriori = inicializarPosteriori()
		
		old_posteriori = posteriori
		
		while(t<maxRep && abs(MLatual - MLdepois)>epsilon){
			print(t)
			if(t == 0){
				MLdepois = 100
			}
			#	M.step()	
			retorno = calcMatrixes()
			
			ckernelM = retorno[[1]]
			cpkernelM = retorno [[2]]		
			tau = apply(posteriori, 2, sum)/n		
			#   E.step()			
			old_posteriori = posteriori
			
			posteriori = E.step()
			
			# print(criterio())
			t = t+1		
			
			alvo = rep(1, n)		
			if(length(which((posteriori[,1] > 0.5) == TRUE)) > length(which((posteriori[,1] <= 0.5) == TRUE)))
			{
				alvo[posteriori[,1]>0.5] = 1
				alvo[posteriori[,1]<=0.5] = 2
			}else{
				alvo[posteriori[,1]>0.5] = 2
				alvo[posteriori[,1]<=0.5] = 1
			}	
			
			if(plot == TRUE){
				if(p == 2){
					x_interval = interval(x[,1],x[,2])
					y_interval = interval(x[,3],x[,4])
					intervalGraph2D(x_interval,y_interval,alvo)			
				}
				if(p == 3){
					x_interval = interval(x[,1],x[,2])
					y_interval = interval(x[,3],x[,4])
					z_interval = interval(x[,5],x[,6])
					intervalGraph3D(x_interval, y_interval, z_interval, alvo = alvo)			
				}		
			}
			if(sum((posteriori - old_posteriori)^2) < 0.001){
				print("parou")
				break
			}
			
			MLatual = MLdepois
			MLdepois = criterio()
		}
		if(MLdepois < melhorML){
            melhorML = MLdepois
            melhorU = posteriori            
        }	
		if(melhorML < epsilon){
			break
		}
	}	

	retorno = list(melhorU,alvo, tau)
	names(retorno) = c('Posteriori','Alvo', 'Tau')
	return(retorno)	
}
