# Generate n=101 equidistant points
m <- 1000
n <- 101
x <- 2*pi*seq(-1, 1, length=n)
# Initialize the matrix of fitted values 
fvlp <- fvnw <- fvss <- matrix(0, nrow= n, ncol= m)
#Generate data fit the data and store 
for (j in 1:m){
 
 # simulate y-values
 y <- (1-x^2) * exp(-0.5 * x^2) + rnorm(length(x), mean =0, sd = 0.2);
 # estimates finded
 fvlp[,j] <- predict(loess(y ~ x, span = 0.75), newdata = x);
 fvnw[,j] <- ksmooth(x, y, kernel="normal", bandwidth= 0.2, x.points=x)$y;
 fvss[,j] <- predict(smooth.spline(y ~ x), x=x)$y;
}
# mean estimates
meanlp = apply(fvlp,1,mean);
meannw = apply(fvnw,1,mean);
meanss = apply(fvss,1,mean);
dmin = min( meanlp, meannw, meanss);
dmax = max( meanlp, meannw, meanss);
matplot(x, meanlp, "l", ylim=c(dmin, dmax))
matlines(x, meannw, col="red")
matlines(x, meanss, col="blue")
legend("topright", legend=c("Loess", "NW Kernel", "Spline"), col=c("black", "
red", "blue"), lty=1, cex=0.6)
# Plot empirical Bias
ytrue = (1-x^2) * exp(-0.5 * x^2);
biaslp = apply(fvlp,1,mean) - ytrue;
biasnw = apply(fvnw,1, mean) - ytrue;
biasss = apply(fvss,1, mean) - ytrue;
d1min = min(biaslp, biasnw, biasss);
d1max = max(biaslp, biasnw, biasss);
matplot(x, biaslp, "l", ylim=c(d1min, d1max), ylab="Empirical Bias");
matlines(x, biasnw, col="red");
matlines(x, biasss, col="blue");
legend("topright", legend=c("Loess", "NW Kernel", "Spline"), col=c("black", "
red", "blue"), lty=1, cex=0.6)
# Plot empirical Variance
Varlp = apply( (fvlp - meanlp)^2, 1,sum) / m;
Varnw = apply( (fvnw - meannw)^2, 1,sum) / m;
Varss = apply( (fvss - meanss)^2, 1,sum) / m;
d2min = min(Varlp, Varnw, Varss);
d2max = max(Varlp, Varnw, Varss);
matplot(x, Varlp, "l", ylim=c(d2min, d2max), ylab="Empirical Variance");
matlines(x, Varnw, col="red");
matlines(x, Varss, col="blue");
legend("topright", legend=c("Loess", "NW Kernel", "Spline"), col=c("black", "
red", "blue"), lty=1, cex=0.6)
# Plot empirical MSE
MSElp = apply( (fvlp - ytrue)^2, 1,sum) / m;
MSEnw = apply( (fvnw - ytrue)^2, 1,sum) / m;
MSEss = apply( (fvss - ytrue)^2, 1,sum) / m;
d3min = min(MSElp, MSEnw, MSEss);
d3max = max(MSElp, MSEnw, MSEss);
matplot(x, MSElp, "l", ylim=c(d3min, d3max), ylab="Empirical MSE");
matlines(x, MSEnw, col="red");
matlines(x, MSEss, col="blue");
legend("topright", legend=c("Loess", "NW Kernel", "Spline"), col=c("black", "
red", "blue"), lty=1, cex=0.6)
## Non-equidistant points
set.seed(79)
x <- 2*pi*sort(c(0.5, -1 + rbeta(50,2,2), rbeta(50,2,2)))
fvlp <- fvnw <- fvss <- matrix(0, nrow= n, ncol= m)
#Generate and fit data
for (j in 1:m){
 # simulate y-values
 y <- (1-x^2) * exp(-0.5 * x^2) + rnorm(length(x), sd=0.2);
 fvlp[,j] <- predict(loess(y ~ x, span = 0.3365), newdata = x);
 fvnw[,j] <- ksmooth(x, y, kernel="normal", bandwidth= 0.2, x.points=x)$y;
 fvss[,j] <- predict(smooth.spline(y ~ x, spar= 0.7163), x=x)$y;
}
# mean estimator
meanlp = apply(fvlp,1,mean);
meannw = apply(fvnw,1,mean);
meanss = apply(fvss,1,mean);
dmin = min( meanlp, meannw, meanss);
dmax = max( meanlp, meannw, meanss);
matplot(x, meanlp, "l", ylim=c(dmin, dmax))
matlines(x, meannw, col="red")
matlines(x, meanss, col="blue")
legend("topright", legend=c("Loess", "NW Kernel", "Spline"), col=c("black", "
red", "blue"), lty=1, cex=0.6)
# Plot empirical Bias
ytrue = (1-x^2) * exp(-0.5 * x^2);
biaslp = apply(fvlp,1,mean) - ytrue;
biasnw = apply(fvnw,1, mean) - ytrue;
biasss = apply(fvss,1, mean) - ytrue;
d1min = min(biaslp, biasnw, biasss);
d1max = max(biaslp, biasnw, biasss);
matplot(x, biaslp, "l", ylim=c(d1min, d1max), ylab="Empirical Bias");
matlines(x, biasnw, col="red");
matlines(x, biasss, col="blue");
legend("topright", legend=c("Loess", "NW Kernel", "Spline"), col=c("black", "
red", "blue"), lty=1, cex=0.6)
# Plot empiricall Variance
Varlp = apply( (fvlp - meanlp)^2, 1,sum) / m;
Varnw = apply( (fvnw - meannw)^2, 1,sum) / m;
Varss = apply( (fvss - meanss)^2, 1,sum) / m;
d2min = min(Varlp, Varnw, Varss);
d2max = max(Varlp, Varnw, Varss);
matplot(x, Varlp, "l", ylim=c(d2min, d2max), ylab="Empirical Variance");
matlines(x, Varnw, col="red");
matlines(x, Varss, col="blue");
legend("topright", legend=c("Loess", "NW Kernel", "Spline"), col=c("black", "
red", "blue"), lty=1, cex=0.6)
# Plot empirical MSE
MSElp = apply( (fvlp - ytrue)^2, 1,sum) / m;
MSEnw = apply( (fvnw - ytrue)^2, 1,sum) / m;
MSEss = apply( (fvss - ytrue)^2, 1,sum) / m;
d3min = min(MSElp, MSEnw, MSEss);
d3max = max(MSElp, MSEnw, MSEss);
matplot(x, MSElp, "l", ylim=c(d3min, d3max), ylab="Empirical MSE");
matlines(x, MSEnw, col="red");
matlines(x, MSEss, col="blue");
legend("topright", legend=c("Loess", "NW Kernel", "Spline"), col=c("black", "
red", "blue"), lty=1, cex=0.6