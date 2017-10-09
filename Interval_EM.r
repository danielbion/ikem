IEM = function(x, tipoMatriz, numRep, plot){	
	inicializarTheta = function(){		
		guess=sample(5:(n-5),2,replace=FALSE)
		
		# if (guess[1] < guess[2]){
			# guess[1]=guess[1] - 4
			# guess[2]=guess[2] + 4
		# }else{
			# guess[2]=guess[2] - 4
			# guess[1]=guess[1] + 4
		# }
		
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
		
		if(length(which(alvo==1)) == 0){
			alvo[1] = 1
		}else if(length(which(alvo==2)) == 0){
			alvo[1] = 2 
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
					
		return(theta)
	}
	
	E.step = function(){
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

	M.step = function(){
		for(i in 1:c){
			theta$tau[i] = sum(posteriori[,i])/n
			
			for(j in 1:pp){
				theta$media[i,j] = sum(posteriori[,i] * x[,j])/sum(posteriori[,i])
			}
		}
		
		mcovA = matrix(0, nrow = p, ncol = p)
		mcovB = matrix(0, nrow = p, ncol = p)
						
		for (i in 1:n){											
			for(p1 in 1:p){
				infIdxp1 = (p1*2) - 1
				supIdxp1 = (p1*2)
				for(p2 in 1:p){
					infIdxp2 = (p2*2) - 1
					supIdxp2 = (p2*2)
					
					W = (x[i,infIdxp1] - theta$media[1:2,infIdxp1]) * (x[i,infIdxp2] - theta$media[1:2,infIdxp2])
					V = (x[i,supIdxp1] - theta$media[1:2,supIdxp1]) * (x[i,supIdxp2] - theta$media[1:2,supIdxp2])
					
					mcovA[p1,p2] = mcovA[p1,p2] + (posteriori[i,1] * (W[1] + V[1]))
					mcovB[p1,p2] = mcovB[p1,p2] + (posteriori[i,2] * (W[2] + V[2]))									
				}			
			}			
		}	
		mcovA = mcovA/(2*sum(posteriori[,1]))
		mcovB = mcovB/(2*sum(posteriori[,2]))		
		
		if(tipoMatriz == 1){
			# Matriz Identidade
			mcovA = matrix(0,nrow=p,ncol=p,byrow=TRUE)
			mcovB = matrix(0,nrow=p,ncol=p,byrow=TRUE)
			for(i in 1:p){				
				mcovA[i,i] = 1
				mcovB[i,i] = 1
			}			
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
		
		return (theta)
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
	maxRep = 500
    epsilon = 0.000001
    
	#Comentar
	# tipoMatriz=3
	# numRep=20
	###########
		
	melhorML = 10000000000;  
		
	for(ii in 1:numRep){		
		cat(' - Modelo', ii, '\n')	
		t = 0
		MLatual = 1	
        MLdepois = -1
		
		posteriori = matrix(0, nrow = n, ncol = c)		
		theta = inicializarTheta()
			
		while(t<maxRep && abs(MLatual - MLdepois)>epsilon){
			if(t == 0){
				MLdepois = 100
			}
			cat(t ,'')
			if(is.nan(theta$mcov[[1]]) || is.nan(theta$mcov[[2]])){
				print('Matriz NaN')
				break
			}
			if(det(theta$mcov[[1]]) < (10^-10) || det(theta$mcov[[2]]) < (10^-10)){	
				print('Matriz não possui inversa')
				break
			}
			
			posteriori = E.step()		
			theta = M.step()			
				
			MLatual = MLdepois
			MLdepois = criterio()		
			
			#print(MLdepois)
			t = t+1		
		}
		cat('\nMLFinal',MLdepois, '\n\n')
		if(MLdepois < melhorML){
            melhorML = MLdepois
            melhorU = posteriori
            melhorTheta = theta
        }	
		if(melhorML < epsilon){
			break
		}
	}
	
	#Plot
	alvo = rep(1, n)		
	if(length(which((melhorU[,1] > 0.5) == TRUE)) > length(which((melhorU[,1] <= 0.5) == TRUE)))
	{
		alvo[melhorU[,1]>0.5] = 1
		alvo[melhorU[,1]<=0.5] = 2
	}else{
		alvo[melhorU[,1]>0.5] = 2
		alvo[melhorU[,1]<=0.5] = 1
	}	
	if(plot){
		if(p == 2){
			x_interval = interval(x[,1],x[,2])
			y_interval = interval(x[,3],x[,4])
			intervalGraph2D(x_interval,y_interval,alvo)
			points(mean(c(melhorTheta$media[1,1],melhorTheta$media[1,2])),mean(c(melhorTheta$media[1,3],melhorTheta$media[1,4])))
			points(mean(c(melhorTheta$media[2,1],melhorTheta$media[2,2])),mean(c(melhorTheta$media[2,3],melhorTheta$media[2,4])))
		}
		if(p == 3){
			x_interval = interval(x[,1],x[,2])
			y_interval = interval(x[,3],x[,4])
			z_interval = interval(x[,5],x[,6])
			intervalGraph3D(x_interval, y_interval, z_interval, alvo = alvo)
			points3d(mean(c(melhorTheta$media[1,1],melhorTheta$media[1,2])),mean(c(melhorTheta$media[1,3],melhorTheta$media[1,4])), mean(c(melhorTheta$media[1,5],melhorTheta$media[1,6])), color = 'black')
			points3d(mean(c(melhorTheta$media[2,1],melhorTheta$media[2,2])),mean(c(melhorTheta$media[2,3],melhorTheta$media[2,4])), mean(c(melhorTheta$media[2,5],melhorTheta$media[2,6])), color = 'black')
		}	
	}	
	retorno = list(melhorU, melhorTheta,alvo, melhorML)
	names(retorno) = c('Posteriori','Theta','Alvo','ML')
	return(retorno)	
}
