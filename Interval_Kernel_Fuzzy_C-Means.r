IKFCM = function(x, c, m, maxRep, monteCarlo){
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
	
	inicializarTheta = function(){		
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
	
	ainicializarK = function(){
		K = matrix(0, nrow=n,ncol=n)
		for(i in 1:n){
			for(j in 1:n){
				a_inf = c(x[i,1],x[i,3])
				a_sup = c(x[i,2],x[i,4])
				b_inf = c(x[j,1],x[j,3])
				b_sup = c(x[j,2],x[j,4])
				K[i,j] = IRBF(a_inf, a_sup, b_inf, b_sup, 1.0)
				#K[i,j] = IPF(a_inf, a_sup, b_inf, b_sup, 2.0)
			}
		}
		return (K)
	}
	
	inicializarK = function(){
		K = matrix(0, nrow=n,ncol=n)
		infSeq =  seq(from = 1, to = pp, by = 2)
		for(i in 1:n){
			for(j in 1:n){
				a_inf = x[i,infSeq]
				a_sup = x[i,infSeq+1]
				b_inf = x[j,infSeq]
				b_sup = x[j,infSeq+1]
			
				K[i,j] = IRBF(a_inf, a_sup, b_inf, b_sup, 1.0)					
			}
		}
		return (K)
	}
	inicializarU = function(){
		U = matrix(runif(n*c), n, c)
		for(i in 1:n){
			U[i,] = U[i,]/sum(U[i,])
		}    
		return(U)
	}

	calcDist = function(){
		D = matrix(0,n,c)
		for(k in 1:c){		
			for(i in 1:n){
				# print(i)
								
				ckIdx = which(U[,k] > 0.5)
				# num2 = 0
				# for(l in ckIdx){	
					# num2 = num2 + U[l, k]^m * K[i, l]
				# }
				# num2 = num2 * 2
				
				num2 = 2 * sum(U[ckIdx, k]^m * K[i, ckIdx])
				
				# den2 = 0
				# for(l in ckIdx){	
					# den2 = den2 + U[l, k]^m
				# }
				den2 = sum(U[ckIdx, k]^m)
				
				num3 = 0
				for(i2 in ckIdx){
					for(l in ckIdx){
						num3 = num3 + U[i2, k]^m * U[l, k]^m * K[i2, l]
					}
				}
				
				# den3 = 0
				# for(l in ckIdx){
					# den3 = den3 + U[l,k]^m
				# }
				# den3 = den3^2
				
				den3 = sum(U[ckIdx, k]^m)^2
				
				D[i, k] = K[i, i] - (num2/den2) + (num3/den3)
			}			
		}
		return(D)
	}

	membership = function(){        
		U = matrix(0,n,c)
		for(k in 1:c){
			for(i in 1:n){		
				d1 = D[i,k] 
				sum = 0.0
				for(aa in 1:c){
					d2 = D[i,aa]
					sum = sum + (d1/d2)^(1.0/(m-1.0))
				}
				U[i,k] = sum^(-1.0) + 0.0000000001
			}
		}
		return(U)
	}
		
	criterio = function(){
		J = 0.0
		for(k in 1:n){
			for(i in 1:c){		
				d = D[k,i]
				J = J + U[k,i]^m * d
			}
		}
		return(J)
	}
	
	n = dim(x)[1]
	pp = dim(x)[2]
	p = pp/2
		
	#maxRep = 20
    epsilon = 0.000001
    numRep = monteCarlo
	
	K = inicializarK()
	
	for(ii in 1:numRep){
		Jc = c()
		t = 0;
		Jatual = 1;	
        Jdepois = -1;
        melhorJ = 10000000000;  
        melhorU = matrix(0, n, p)      
		
		U = inicializarU()
		
		while(t<maxRep && abs(Jatual - Jdepois)>epsilon){
			print(t)
			#P = prototipo()
			D = calcDist()
			U = membership()

			Jatual = Jdepois
			Jdepois = criterio()
			
			Jc = rbind(Jc, Jdepois)
			t = t+1
						
			alvo = rep(0, 499)
			alvo[U[,1]>0.5] = 1
			alvo[U[,2]>0.5] = 2
			# alvo[U[,3]>0.33] = 3

			x_interval = interval(x[,1],x[,2])
			y_interval = interval(x[,3],x[,4])
			intervalGraph2D(x_interval,y_interval,alvo)
        }
        
        if(Jdepois < melhorJ){
            melhorJ = Jdepois
            melhorU = U
			melhorJc = Jc
           # melhorP = P
		   	alvo = rep(0, 499)
			alvo[melhorU[,1]>0.5] = 1
			alvo[melhorU[,2]>0.5] = 2

			x_interval = interval(x[,1],x[,2])
			y_interval = interval(x[,3],x[,4])
			intervalGraph2D(x_interval,y_interval,alvo)
        }		
	}
	retorno = list(melhorU,melhorJ,melhorJc)
	names(retorno) = c('U','J','Jc')
	return(retorno)
}
