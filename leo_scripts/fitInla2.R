# This script belongs to Patric Brown: https://groups.google.com/g/r-inla-discussion-group/c/wjuadgCFnxE


library(spdep)
source("leo_scripts/sim.r")
thesim = simLGCP(F, N=100)


Napprox = 11
library(spdep)
library(raster)


thebb = bbox(thesim$sp)
rasterTemplate = raster(nrows=Napprox, ncols=Napprox, 
		xmn=thebb[1,1], xmx=thebb[1,2], 
		ymn=thebb[2,1], ymx=thebb[2,2])

library(raster)
x = rasterize(thesim$sp, rasterTemplate, field="x" , fun='sum')
x@data@values[is.na(x@data@values)] = 0

table(x@data@values)

if(F) {
	plot(x)
	points(thesim$sp)
}

reRast = function(x) {
 resample(
		raster(t(x), 	
				xmn=bbox(rasterTemplate)[1,1], xmx=bbox(rasterTemplate)[1,2], 
			ymn=bbox(rasterTemplate)[2,1], ymx=bbox(rasterTemplate)[2,2]),
rasterTemplate)@data@values
}

mydata = data.frame(count = x@data@values,
		offset = reRast(thesim$offset),
		x=reRast(thesim$x),
		y=reRast(thesim$y),
		cellID = seq(1, length(x@data@values)))

		
gridsize = diff(bbox(rasterTemplate)[1,])/(dim(x)[1]-1)		
		
		

	library(INLA)
	theformula = count ~ x + y + offset(offset) +
			f(cellID,model="matern2d", 
		    		ncol=dim(rasterTemplate)[1], nrow=dim(rasterTemplate)[2], 
					nu=2)

	inlaRes = inla(theformula, data=mydata, family="poisson", verbose=T)	

	thebb=bbox(rasterTemplate)
thePred = raster(
		matrix(inlaRes$summary.random$cellID$mean,
				dim(rasterTemplate)[1], dim(rasterTemplate)[2]),
		xmn=thebb[1,1], xmx=thebb[1,2], 
		ymn=thebb[2,1], ymx=thebb[2,2])

plot(thePred)

