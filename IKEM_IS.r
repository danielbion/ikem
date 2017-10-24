IKEMIS = function(x, repPerIteration, D, plot, iterations){	
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
		
	E.step = function(){				
		for (i in 1:n)
		{
			A = c()
			B = c()						
			
			den = 0
			num = c()
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
				fdp = normal(A[k], B[k], theta$mcov[[k]])
				num = cbind(num, theta$tau[k] * fdp)				
			}	

			den = sum(num)
			
			if(is.nan(den)){
				cat('E-step: DEN == 0')
			}
			
			posteriori[i, ] = num/den			
		}
		return(posteriori)
	}
	
	normal = function(A, B, Sigma){
		num = exp(-0.5 * (A + B) * 0.5)
		den1 = ((2 * pi)^p)^0.5
		den2 = det(Sigma)
		return (num/(den1 * den2))
	}

	M.step = function(){	
		#tau = apply(posteriori, 2, sum)/n		
		for(k in 1:c){
			theta$tau[k] = sum(posteriori[,k])/n
		
			infSeq =  seq(from = 1, to = pp, by = 2)
			
			num = 0
			den = 0
			for(i in 1:n){
				a_inf = x[i,infSeq]
				a_sup = x[i,infSeq+1]
				b_inf = theta$media[k,infSeq]
				b_sup = theta$media[k,infSeq+1]
				
				kernelDist = IRBF(a_inf, a_sup, b_inf, b_sup, 1.0)				
				num = num + (posteriori[i,k] * x[i,] * kernelDist)	
				den = den + (posteriori[i,k] * kernelDist)	
			}
			theta$media[k,] = num / den		
		}		
				
		for(k in 1:c){
			mcov = matrix(0, nrow = p, ncol = p)
			for (i in 1:n){		
				m = kernelizedDistance(x, theta$media, i, k)
				mcov = mcov + (posteriori[i,k] * m)
			}
			theta$mcov[[k]] = mcov/(sum(posteriori[,k]))
		}
		
		return (theta)
	}
	
	kernelizedDistance = function(x, mu, i, k){	
		m = matrix(0, nrow = p, ncol = p)
		
		infSeq =  seq(from = 1, to = pp, by = 2)
		
		x_inf = x[i,infSeq]
		x_sup = x[i,infSeq+1]
		mu_inf = mu[k,infSeq]
		mu_sup = mu[k,infSeq+1]
				
		for(h in 1:p){	
			for(j in 1:p){				
				m[h, j] =  IRBF(x_inf[h], x_sup[h], x_inf[j], x_sup[j], 1.0) -
							  IRBF(x_inf[h], x_sup[h], mu_inf[j], mu_sup[j], 1.0) -
							  IRBF(mu_inf[h], mu_sup[h], x_inf[j], x_sup[j], 1.0) +
							  IRBF(mu_inf[h], mu_sup[h], mu_inf[j], mu_inf[j], 1.0)
			}
		}
		return (m)
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
		
		for(ii in 1:iterations){		
			t = 0	
			MLatual = 1	
			MLdepois = -1

			retorno = IEM(x, tipoMatriz = 3, numRep = 1, plot = FALSE)
			posteriori = retorno$Posteriori
			theta = retorno$Theta
			old_posteriori = posteriori
						
			while(t<repPerIteration && abs(MLatual - MLdepois)>epsilon){
				print(t)
				if(t == 0){
					MLdepois = 100
				}
				#	M.step()	
				#retorno = calcMatrixes()
				
				#ckernelM = retorno[[1]]
				#cpkernelM = retorno [[2]]		
				#tau = apply(posteriori, 2, sum)/n		
				#   E.step()			
				old_posteriori = posteriori
				
				#posteriori = E.step()
				theta = M.step()	
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

		retorno = list(melhorU,alvo, theta$tau)
		names(retorno) = c('Posteriori','Alvo', 'Tau')
		return(retorno)	
	
}
