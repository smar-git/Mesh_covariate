
# This script belongs to Patric Brown: https://groups.google.com/g/r-inla-discussion-group/c/wjuadgCFnxE
simLGCP = function(plot.it=T, N=200){
library(RandomFields)
# try mu 3 and 4
param = list(mu=3, betax=1, betay=1.5, sigma=0.8, 
		range=0.3, roughness=2,model="matern", 
		gridsize=1/N)

mybbox = cbind(c(0,0), c(1,1))

dimnames(mybbox)=list(c("x","y"), c("min","max"))

x <- seq(mybbox["x","min"], mybbox["x","max"], by=param$gridsize) 
y <- seq(mybbox["y","min"], mybbox["y","max"], by=param$gridsize) 

f <- GaussRF(x=x, y=y, model="matern", grid=TRUE,
        param=c(mean=0, variance=1, nugget=0, 
                scale=param$range, alpha=param$roughness))

f = param$sigma*(f-mean(f))/sd(as.vector(f))

xmat = matrix(x, length(x), length(y))
ymat = matrix(y, length(x), length(y), byrow=T)

offset = matrix(seq(log(0.5), log(2), len=length(x)*length(y)), length(x), length(y))

lambda = exp(offset + f + param$mu + param$betax*xmat + param$betay*ymat)

pointsMat = matrix(rpois(length(lambda), param$gridsize^2*lambda), length(x), length(y))

coordsMat = outer(x, 1i*y, FUN="+")

Npoints = sum(pointsMat>0)

thepoints = coordsMat[pointsMat>0] +
		runif(Npoints, -param$gridsize, param$gridsize)+
		1i*runif(Npoints, -param$gridsize, param$gridsize)


Ndoubles = sum(pointsMat>1.5)
thepoints = c(thepoints, 
		coordsMat[pointsMat>1.5]+
				runif(Ndoubles, -param$gridsize, param$gridsize)+
				1i*runif(Ndoubles, -param$gridsize, param$gridsize)
)
if(plot.it) {
image(x, y, lambda)

points(thepoints)

}

return(list(points=thepoints, params=param, U=f,
				spPoints = 
						SpatialPointsDataFrame(
								cbind(Re(thepoints), Im(thepoints)),
								bbox = mybbox,
								data=data.frame(x=rep(1,length(thepoints)))) ,
			x = xmat, y=ymat, offset=offset)
)

}
