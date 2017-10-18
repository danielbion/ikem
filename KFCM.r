maxIteration = 50;
J = .Machine$integer.max;

KFCM <- function(data,centers,parM,sigma,tipo)
{
	count = 0;
	P = initializePrototypes(data,centers);
	Ubefore = initializeMembership(data,centers);
	Jbefore = J + 1.0;
	K = 0;
	parcel = 0;
	while(abs(Jbefore - J) > 0.0001 && count < maxIteration)
	{
		count = count+1;
		if(tipo == 1){
			K = updateKernels(data,P,sigma);			
		}
		if(tipo == 2){
			parcel = updateParcel(data,Ubefore,parM,sigma);
		}
		#print(K);
		DK = updateDistKernel(data,K,Ubefore,tipo,sigma,parM,parcel);
		#print(DK);
		U = updateMembership(DK,parM);
		#print(U);
		if(tipo == 1){
			P = updatePrototypes(data,U,K,parM);
		}
		#print(P);
		Jbefore = J;
		J = updateCriterion(U,DK,parM);
		#print(J);
		Ubefore = U;		
	}
	
	#print(Ubefore);
	
	L = getPartition(Ubefore);	
	resp = list(criterio = J, part = L, matriz = Ubefore);
	return(resp);
	#print(L);
}

initializePrototypes <- function(data,centers)
{
	nVar = ncol(data);
	nProt = length(centers);
	P = matrix(0,nProt,nVar);
	for(k in 1:nProt){
		P[k,] = data[centers[k],];		
	}
	return(P);
}

initializeMembership <- function(data,centers)
{
	nObj = nrow(data);
	nProt = length(centers);
	U = matrix(0,nObj,nProt);
	for(i in 1:nObj){
		U[i,] = runif(nProt, 0.0, 1.0);
		U[i,] = U[i,]/sum(U[i,]);
	}
	return(U);
}

updatePrototypes <- function(data,memberships,matrixKernel,parM)
{
	nObj = nrow(data);
	nVar = ncol(data);
	nProt = ncol(memberships);
	P = matrix(0,nProt,nVar);
	for(k in 1:nProt){
		s = 0;
		ss = 0
		for(j in 1:nObj){
			delta = (memberships[j,k])^(parM+0.0);
			ker = matrixKernel[j,k];
			obj = data[j,];
			s = s + ker*delta*obj;
			ss = ss + ker*delta;
		}
		s = s/ss;
		P[k,] = s;
	}
	return(P);
}

updateKernels <- function(data,prototypes,sigma)
{
	nObj = nrow(data);
	nProt = nrow(prototypes);
	K = matrix(0,nObj,nProt);
	for(j in 1:nObj){
		obj = data[j,];
		for(k in 1:nProt){
			prot = prototypes[k,];
			distance = sum((obj-prot)^2);
			kernelGaussian = exp(-distance/(2*sigma^2.0));
			K[j,k] = kernelGaussian;
		}
	}	
	return(K);
}

updateParcel <- function(data,U,parM,sigma)
{
	parcel = c();
	nObj = nrow(U);
	nProt = ncol(U);
	for(k in 1:nProt){
		soma = 0;
		p3 = 0;
		for(jj in 1:nObj){
			um1 = U[jj,k]^parM;
			soma = soma + um1;
			for(jjj in 1:nObj){
				um2 = U[jjj,k]^parM;
				p3 = p3 + um1*um2*calcKernel(data,jj,jjj,sigma);						
			}
		}
		p3 = p3/soma^2;
		parcel[k] = p3;
	}
	return(parcel);
}

updateDistKernel <- function(data,K,U,tipo,sigma,parM,parcel)
{
	nObj = nrow(U);
	nProt = ncol(U);
	DK = matrix(0,nObj,nProt);
	for(j in 1:nObj){
		for(k in 1:nProt){			
			distance = 0.0;
			if(tipo == 1){ #Input space (kernelizada)
				ker = K[j,k];
				distance = 2.0*(1.0 - ker);
			}
			else if(tipo == 2){ # Feature space
				p1 = calcKernel(data,j,j,sigma);				
				p2 = 0;
				soma = 0;
				for(jj in 1:nObj){
					um = U[jj,k]^parM;
					p2 = p2 + um*calcKernel(data,j,jj,sigma);
					soma = soma + um;
				}
				p2 = p2/soma;
				p3 = 0;
				soma = 0;
				#for(jj in 1:nObj){
				#	um1 = U[jj,k]^parM;
				#	soma = soma + um1;
				#	for(jjj in 1:nObj){
				#		um2 = U[jjj,k]^parM;
				#		p3 = p3 + um1*um2*calcKernel(data,jj,jjj,sigma);						
				#	}
				#}
				#p3 = p3/soma^2;
				p3 = parcel[k];
				distance = p1 - 2.0*p2 + p3;
			}
			DK[j,k] = distance + 0.00000001;
		}
	}	
	return(DK);
}

calcKernel <- function(data,x1,x2,sigma)
{
	obj1 = data[x1,];
	obj2 = data[x2,];
	distance = sum((obj1-obj2)^2);
	ker = exp(-distance/(2*sigma^2));
	
	return(ker);
}

updateMembership <-function(distKernel,parM)
{
	nObj = nrow(distKernel);
	nProt = ncol(distKernel);
	U = matrix(0,nObj,nProt);
	for(i in 1:nObj){
		for(k in 1:nProt){
			d1 = distKernel[i,k];
			soma = 0;
			for(kk in 1:nProt){
				d2 = distKernel[i,kk];				
				soma = soma + (d1/d2)^(1.0/(parM-1.0));
			}					
			soma = soma^(-1.0);			
			U[i,k] = soma;
		}					
	}
	return(U);	
}

updateCriterion <- function(memberships,distKernel,parM)
{
	J = 0;
	nObj = nrow(memberships);
	nProt = ncol(memberships);
	for(j in 1:nObj){
		for(k in 1:nProt){
			delta = (memberships[j,k])^(parM+0.0);
			distance = distKernel[j,k];
			J = J + delta*distance;
		}
	}
	return(J);
}


getPartition <- function(memberships)
{
	nObj = nrow(memberships);
	L = 0;
	for(i in 1:nObj){
		L[i] = which.max(memberships[i,]);
	}
	return(L);
}


