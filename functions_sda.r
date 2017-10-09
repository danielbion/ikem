library('rgl')

base1 = function(n1, n2, var1, var2, range1, range2, offset2){
	dimensao = 2
	mcov1 = diag(dimensao)		
	diag(mcov1) = var1*var1
	mcov2 = mcov1
	diag(mcov2) = var2*var2
	
	p1 = 0.7
	p2 = -0.7
	for(i in 1:dimensao){
		for(j in 1:dimensao){
			if(i!=j){
				mcov1[i,j] = p1*var1*var1
				mcov2[i,j] = p2*var2*var2
			}
		}
	}
	x1 = mvrnorm(n1, c(0,0), mcov1)
	x2 = mvrnorm(n2, offset2, mcov2)
	
	x = rbind(x1, x2)
	
	xr1 = runif(dim(x1)[1], min = 0.1, max = range1)
	yr1 = runif(dim(x1)[1], min = 0.1, max = range1)
	
	xr2 = runif(dim(x2)[1], min = 0.1, max = range2)
	yr2 = runif(dim(x2)[1], min = 0.1, max = range2)
	
	xr = c(xr1, xr2)
	yr = c(yr1, yr2)
	
	x = cbind(x[,1] - xr/2, x[,1] + xr/2, x[,2] - yr/2, x[,2] + yr/2)
	
	x = normalizarDadosSimbolicos(x)
	
	alvo = rep(1, n1 + n2)
	alvo[(n1+1):(n1+n2)] = 2
	
	x_interval = interval(x[,1],x[,2])
	y_interval = interval(x[,3],x[,4])	
	intervalGraph2D(x_interval, y_interval, alvo)
	
	alvo[alvo==2] = 0
	return(list(x, alvo))
}

base2=function(n1, n2, var2, range1, range2, offset2){
	a = runif(n1, min = 5, max = 15)	
	
	E <- rnorm(n1, mean=0, sd=0.1)
	
	ystart = sin(a*pi*2) * 3 + E
	xstart = cos(a*pi*2) * 3 + E
	
	x1 = cbind(xstart, ystart)
	for(i in 1:dim(x1)[1]){
		if(x1[i,1] > 0 || x1[i,2] < 0){
			x1[i,] = -100
		}		
	}
	x1 = x1[x1[,1] != -100,]
	
	
	dimensao = 2
	mcov2 = diag(dimensao)		
	diag(mcov2) = var2*var2	
	p2 = 0.01
	for(i in 1:dimensao){
		for(j in 1:dimensao){
			if(i!=j){				
				mcov2[i,j] = p2*var2*var2
			}
		}
	}	
	x2 = mvrnorm(n2, offset2, mcov2)	

	x = rbind(x1, x2)
	
	xr1 = runif(dim(x1)[1], min = 0.05, max = range1)
	yr1 = runif(dim(x1)[1], min = 0.05, max = range1)
	
	xr2 = runif(dim(x2)[1], min = 0.05, max = range2)
	yr2 = runif(dim(x2)[1], min = 0.05, max = range2)
	
	xr = c(xr1, xr2)
	yr = c(yr1, yr2)
	
	x = cbind(x[,1] - xr/2, x[,1] + xr/2, x[,2] - yr/2, x[,2] + yr/2)
	
	x = normalizarDadosSimbolicos(x)
	
	alvo = rep(1, dim(x1)[1] + n2)
	alvo[(dim(x1)[1]+1):(dim(x1)[1]+n2)] = 2
	
	x_interval = interval(x[,1],x[,2])
	y_interval = interval(x[,3],x[,4])	
	intervalGraph2D(x_interval, y_interval, alvo)
	
	alvo[alvo==2] = 0
	return(list(x, alvo))
}

base3=function(n1, n2, range1, range2, var1){
	dimensao = 2
	mcov1 = diag(dimensao)		
	diag(mcov1) = var1*var1
	
	p1 = 0.2
	for(i in 1:dimensao){
		for(j in 1:dimensao){
			if(i!=j){
				mcov1[i,j] = p1*var1*var1
			}
		}
	}
	x1 = mvrnorm(n1, c(0,0), mcov1)
	
	
	a = runif(n2, min = 5, max = 15)		
	E <- rnorm(n2, mean=0, sd=0.1)	
	ystart = sin(a*pi*2) * 3 + E
	xstart = cos(a*pi*2) * 3 + E	
	x2 = cbind(xstart, ystart)

	x = rbind(x1, x2)
	
	xr1 = runif(dim(x1)[1], min = 0.05, max = range1)
	yr1 = runif(dim(x1)[1], min = 0.05, max = range1)
	
	xr2 = runif(dim(x2)[1], min = 0.05, max = range2)
	yr2 = runif(dim(x2)[1], min = 0.05, max = range2)
	
	xr = c(xr1, xr2)
	yr = c(yr1, yr2)
	
	x = cbind(x[,1] - xr/2, x[,1] + xr/2, x[,2] - yr/2, x[,2] + yr/2)
	
	x = normalizarDadosSimbolicos(x)
	
	alvo = rep(1, dim(x1)[1] + n2)
	alvo[(dim(x1)[1]+1):(dim(x1)[1]+n2)] = 2
	
	x_interval = interval(x[,1],x[,2])
	y_interval = interval(x[,3],x[,4])	
	intervalGraph2D(x_interval, y_interval, alvo)
	
	alvo[alvo==2] = 0
	return(list(x, alvo))
}

normalizarDados=function(x){
	min1 = c()
	max1 = c()
	for(i in 1:dim(x)[2]){
		min1 = cbind(min1, min(x[,i]))		
		max1 = cbind(max1, max(x[,i]))
	}
	
	for(i in 1:dim(x)[1]){
		for(j in 1:dim(x)[2]){
			x[i,j] = (x[i,j] - min1[j]) / (max1[j] - min1[j])			
		}
	}
	return (x)
}

normalizarDadosSimbolicos=function(x){
	norm = c()
	for(i in seq(from=1, to=dim(x)[2], by=2)){
		norm = cbind(norm, min(x[,i]), max(x[,i+1]))		
	}
	
	for(i in 1:dim(x)[1]){
		for(j in seq(from=1, to=dim(x)[2], by=2)){
			x[i,j] = (x[i,j] - norm[j]) / (norm[j+1] - norm[j])
			x[i,j+1] = (x[i,j+1] - norm[j]) / (norm[j+1] - norm[j])
		}
	}
	return (x)
}

interval = function(min,max) {
	rval= list(
			minValue = min,
			maxValue = max
	)
	class(rval) = "interval"
	return(rval)
}

intervalGraph3D =
function(xIntervals,yIntervals,zIntervals, color=c('transparent','gray','red','blue','green'), alvo){	
	
	members <- match.arg(color)
	## get the vectors of the limits from the Intervals passed as parameter
	xmin = xIntervals$minValue
	xmax = xIntervals$maxValue
	ymin = yIntervals$minValue
	ymax = yIntervals$maxValue
	zmin = zIntervals$minValue
	zmax = zIntervals$maxValue
	## get the min and max limits for x, y and z variables
	minx = min(xmin)
	miny = min(ymin)
	minz = min(zmin)
	maxx = max(xmax)
	maxy = max(ymax)
	maxz = max(zmax)
	
	## drawing plot axes
	z <- seq(minz,maxz,abs(maxz - minz/4))
	x <- seq(minx,maxx,abs(maxx - minx/4))
	y <- seq(miny,maxy,abs(maxy - miny/4))

	open3d()
	plot3d(x, y, z, col = 'transparent',xlim = c(minx, maxx), ylim = c(miny, maxy), zlim = c(minz, maxz), box = FALSE)
		
	addSegment <- function(cubo, p1, p2){
		segments3d(c(cubo[p1,1],cubo[p2,1]),c(cubo[p1,2],cubo[p2,2]),c(cubo[p1,3],cubo[p2,3]))
	}
	
	addQuad <- function(cubo, p1, p2, p3, p4){
		quads3d(c(cubo[p1,1],cubo[p2,1],cubo[p3,1],cubo[p4,1]),c(cubo[p1,2],cubo[p2,2],cubo[p3,2],cubo[p4,2]),c(cubo[p1,3],cubo[p2,3],cubo[p3,3],cubo[p4,3]))
	}
		
	## drawing intervals
	createsRectangles <- function(n,...){
		cubo <- rbind(c(xmin[n],ymin[n],zmax[n]), c(xmin[n],ymin[n],zmin[n]), c(xmax[n],ymin[n],zmin[n]), c(xmax[n],ymax[n],zmin[n]), c(xmax[n],ymax[n],zmax[n]), c(xmin[n],ymax[n],zmax[n]), # < 6 outer
				c(xmax[n],ymin[n],zmax[n]), c(xmin[n],ymax[n],zmin[n]))		
		material3d(color = 'black', alpha = 1)
		addSegment(cubo, 1, 2)
		addSegment(cubo, 1, 7)
		addSegment(cubo, 1, 6)
		addSegment(cubo, 2, 3)
		addSegment(cubo, 2, 8)
		addSegment(cubo, 3, 4)
		addSegment(cubo, 3, 7)
		addSegment(cubo, 4, 5)
		addSegment(cubo, 4, 8)
		addSegment(cubo, 5, 6)
		addSegment(cubo, 5, 7)
		addSegment(cubo, 6, 8)
		
		if(color[1] != 'transparent'){
			material3d(color = color, alpha = 0.15)
			addQuad(cubo, 1,2,3,7)
			addQuad(cubo, 3,4,5,7)
			addQuad(cubo, 4,5,6,8)
			addQuad(cubo, 1,2,8,6)
			addQuad(cubo, 1,6,5,7)
			addQuad(cubo, 2,3,4,8)
		}
	}
	
	createsRectanglesWithTarget <- function(n,color,...){
		cubo <- rbind(c(xmin[n],ymin[n],zmax[n]), c(xmin[n],ymin[n],zmin[n]), c(xmax[n],ymin[n],zmin[n]), c(xmax[n],ymax[n],zmin[n]), c(xmax[n],ymax[n],zmax[n]), c(xmin[n],ymax[n],zmax[n]), # < 6 outer
				c(xmax[n],ymin[n],zmax[n]), c(xmin[n],ymax[n],zmin[n]))		
		material3d(color = color, alpha = 1)
		addSegment(cubo, 1, 2)
		addSegment(cubo, 1, 7)
		addSegment(cubo, 1, 6)
		addSegment(cubo, 2, 3)
		addSegment(cubo, 2, 8)
		addSegment(cubo, 3, 4)
		addSegment(cubo, 3, 7)
		addSegment(cubo, 4, 5)
		addSegment(cubo, 4, 8)
		addSegment(cubo, 5, 6)
		addSegment(cubo, 5, 7)
		addSegment(cubo, 6, 8)		
		
		
		material3d(color = color, alpha = c(0.1, 0.5)[alvo[i]])
		addQuad(cubo, 1,2,3,7)
		addQuad(cubo, 3,4,5,7)
		addQuad(cubo, 4,5,6,8)
		addQuad(cubo, 1,2,8,6)
		addQuad(cubo, 1,6,5,7)
		addQuad(cubo, 2,3,4,8)
		
	}
	
	if(missing(alvo)){
		for (i in 1:length(xmin)) {		
			createsRectangles(i)
		}
	}else{
		for (i in 1:length(xmin)) {				
			# createsRectanglesWithTarget(i, c("#92AF38","#F38211",'blue')[alvo[i]])
			createsRectanglesWithTarget(i, c("#666666","#333333",'blue')[alvo[i]])
			#createsRectanglesWithTarget(i, c("red","blue",'blue')[alvo[i]])
		}
	}
}

intervalGraph2D = function(xIntervals,yIntervals,alvo){	
	xmin = xIntervals$minValue
	xmax = xIntervals$maxValue
	ymin = yIntervals$minValue
	ymax = yIntervals$maxValue
			
	plot(c(xmin, xmax), c(ymin, ymax), type="n", xlab="X", ylab="Y",ylim=c(min(ymin), max(ymax)))
	if(missing(alvo)){
		if(xmin == xmax && ymin == ymax){
			points(xmin, ymin)
		}
		rect(xmin, ymin, xmax, ymax)			
	}else{
		if(xmin == xmax && ymin == ymax){
			# points(xmin, ymin, col = c("#666666","black")[alvo])
			# points(xmin, ymin, col = c("#92AF38","#F38211")[alvo])
			points(xmin, ymin, col = c("red","blue","green")[alvo])
		}
		rect(xmin, ymin, xmax, ymax, border = c("#666666","black", "green", "red")[alvo], lwd=c(1,2,1,1)[alvo])	
		#rect(xmin, ymin, xmax, ymax, border = c("#92AF38","#F38211")[alvo], lwd=c(1,2)[alvo])	
		#rect(xmin, ymin, xmax, ymax, border = c("red","blue","green")[alvo], lwd=c(1,2,2)[alvo])	
	}
}

CenterRangetoMinMax = function(x){
	y = cbind(x[,1] - x[,2], x[,1] + x[,2], x[,3] - x[,4], x[,3] + x[,4])
	return (y)
}