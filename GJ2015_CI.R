	library(Rcpp)
	sourceCpp("GJ2015_simulation.cpp")
	
	sim = 101
  	set.seed(sim)
  	
  	# sample size
  	n=100
  	
  	#parameter truncated exponential distribution
  	B=2
  	
  	#generate truncated exponential distribution
	X=runif(n,0,1)
	X=-log(1-(1-exp(-B))*X)
	
	#compute confidence intervals
	intervals <- ComputeIntervals(X,B)
	
	# make a picture of the confidence intervals
	D <- intervals$CI
	
   	x <-D[,1]
   	y1<-D[,2]
   	y2<-D[,3]
   	y3<-D[,4]

   	plot(c(-10000,-10000),xlim=c(0,B), ylim=c(0.0,3), main= "", ylab="",xlab="",bty="n",las=1)
   	lines(x, y2,lty=1,lwd=2,col="red",type="s")
   	f <- function(x) {exp(-x)/(1-exp(-B))}
   	x0<-seq(0,max(x),by=0.01)
   	y0<-f(x0)
   	lines(x0,y0,lty=2,lwd=2,col="blue")
   	segments(x,y1,x,y3)

   


