qqpval <- function(p, cone=TRUE, log=TRUE, add=FALSE,...){
	p <- p[ 0< p & p < 1]
	p <- sort(p)
	if(log) p <- sort(-log10(p))
	N <- length(p)
	if(log){
		Ep <- sort(-log10(ppoints(N)))
		cone.lower <- -log10(qbeta(0.025, 1:N, N:1))
		cone.upper <- -log10(qbeta(0.975, 1:N, N:1))
		cone.x <- sort(-log10(ppoints(N)))
	} else {
		Ep <- sort(ppoints(N))
		cone.lower <- sort(qbeta(0.025, 1:N, N:1))
		cone.upper <- sort(qbeta(0.975, 1:N, N:1))
		cone.x <- rev(sort(ppoints(N)))
	}
	if(add){
		points(y=p, x=Ep,...)
	} else{
		if(log){
			plot(y=p,x=Ep,type="n", xlab = "Expected p-value, -log10", 
           ylab = "Observed p-value, -log10",...)
		} else {
		  plot(y=p,x=Ep,type="n", xlab = "Expected p-value", 
           ylab = "Observed p-value",...)
		}
		abline(0,1)
		polygon(y=c(rev(cone.lower),cone.upper),x=c(cone.x,rev(cone.x)), 
            col=1,lty=2,density=0)
		points(y=p, x=Ep,...)
	}
}


lasso <- function(y,X,lambda, beta=NULL, maxit=1000, tol = 2e-16, rescale=TRUE){
	if(is.null(beta)) beta <- rep(0,ncol(X))
	stopifnot(length(y) ==nrow(X) & ncol(X)==length(beta))
  if(rescale){
    X <- scale(X)*sqrt(nrow(X))/sqrt(nrow(X)-1)
    y <- scale(y,scale=FALSE)
  }
  res <- y-X%*%beta
	if(length(lambda)==1) lambda = rep(lambda,ncol(X))
  stopifnot(all(lambda >=0))
	stopifnot(all(!is.na(y)) & all(!is.na(X)) & all(!is.na(beta)))
	
  out <- .C("lasso",
		r=as.double(res), 
		X = as.double(X), 
		y=as.double(y), 
		beta=as.double(beta),
		p=as.integer(ncol(X)), 
		n=as.integer(nrow(X)), 
		lambda = as.double(lambda), 
		maxit= as.integer(maxit),
		eps = as.double(tol)
		)
	return(out)
}

lassoscore <- function(y,X, lambda=0, tol = .Machine$double.eps, maxit=1000, 
                       verbose=FALSE, subset = NULL, resvar=NULL){
	X <- apply(X,2,scale)*sqrt(nrow(X))/sqrt(nrow(X)-1)
  y <- scale(y,scale=FALSE)
  out0 <- lasso(y=y,X=X, lambda=lambda, tol = tol, maxit=maxit)
	wh <- out0$beta != 0
	if(is.null(resvar)) resvar <- sum(out0$r^2)/(out0$n-sum(wh)-1)
	if(is.null(subset)){
		subset <- 1:out0$p
	}
	if(verbose){
	  cat("\nProgress:\n")
	  pb <- txtProgressBar(min = 0, max = length(subset), style = 3)
	  pb.i <- 0
	}
	scores <- scorevar.cons <- scorevar.asm <- scorevar.sand <- numeric(out0$p)
  scores[-subset] <- scorevar.cons[-subset] <- scorevar.asm[-subset] <- scorevar.sand[-subset] <- NA
	if(length(intersect(which(!wh),subset)) > 0){
		scores[intersect(which(!wh),subset)] <- colSums(
      X[,intersect(which(!wh),subset),drop=FALSE]*out0$r)/sqrt(out0$n)
		Xs <- X[,wh,drop=FALSE]
		
		
		for(j in intersect(which(!wh),subset)){
			xx <- X[,j,drop=FALSE]
			scorevar.cons[j] <- resvar*crossprod(xx)/out0$n
			if(ncol(Xs) == 0){
				scorevar.asm[j] <- scorevar.sand[j] <- scorevar.cons[j]
			} else if(nrow(Xs) > ncol(Xs) & ncol(Xs) >0){
				Ui <- solve(crossprod(Xs)/out0$n)
				V <- with(out0,crossprod(as.vector(r)*Xs)/n - lambda[1]^2*outer(sign(beta[wh]),sign(beta[wh])))
				Va <- with(out0,crossprod(Xs,xx*r^2)/n -
				  lambda[1]*sign(beta[wh])*crossprod(xx*r)/n)
				Ua <-  crossprod(Xs,xx)/out0$n
        va <- crossprod(xx*out0$r)/out0$n
        
				scorevar.asm[j] <- resvar*(crossprod(xx)/out0$n - t(Ua)%*%Ui%*%Ua )
			  scorevar.sand[j] <- va + t(Ua)%*%(Ui)%*%(V%*%Ui%*%Ua - 2*Va) 
			}	
			if(verbose){
			  pb.i <- pb.i+1
			  setTxtProgressBar(pb, pb.i)
			}		
		}	
	}
	for(i in intersect(which(wh),subset)){
		out <- lasso(y=y,X=X[,-i], lambda=lambda, beta = out0$beta[-i], tol = tol, maxit=maxit)
		xx <- X[,i,drop=FALSE]
		Xs <- X[,-i,drop=FALSE][,out$beta !=0,drop=FALSE]

    scores[i] <- sum(out$r*X[,i])/sqrt(out$n)
		scorevar.cons[i] <- resvar*crossprod(xx)/out0$n
		if(ncol(Xs) == 0){
			scorevar.asm[i] <- scorevar.cons[i]
		} else if(nrow(Xs) > ncol(Xs) & ncol(Xs) >0){ 
		  Ui <- solve(crossprod(Xs)/out0$n)
		  V <- with(out,crossprod(as.vector(r)*Xs)/n - lambda[1]^2*outer(sign(beta[beta !=0]),sign(beta[beta !=0])))
		  Va <- with(out,crossprod(Xs,xx*r^2)/n -
		               lambda[1]*sign(beta[beta !=0])*crossprod(xx*r)/n)
		  Ua <-  crossprod(Xs,xx)/out0$n
		  va <- crossprod(xx*out0$r)/out0$n
		  
		  scorevar.asm[i] <- resvar*(crossprod(xx)/out0$n - t(Ua)%*%Ui%*%Ua )
		  scorevar.sand[i] <- va + t(Ua)%*%(Ui)%*%(V%*%Ui%*%Ua - 2*Va)
		}
		if(verbose){
		  pb.i <- pb.i+1
		  setTxtProgressBar(pb, pb.i)
		}
	}
	re <- list(
    "fit"=out0,
    "scores" = scores,
    "resvar" = resvar,
    "scorevar.cons" = scorevar.cons,
    "scorevar.asm" = scorevar.asm,
    "scorevar.sand" = scorevar.sand,	
	  "p.cons" = pchisq(scores^2/scorevar.cons,df=1,lower.tail=FALSE),
	  "p.asm" = pchisq(scores^2/scorevar.asm,df=1,lower.tail=FALSE),
	  "p.sand" = pchisq(scores^2/scorevar.sand,df=1,lower.tail=FALSE))
	
	class(re) <- "lassoscore"
	if(verbose) close(pb)
  return(re)
}

print.lassoscore <- function(x,...){
	cat("An object of class `lassoscore'\nbased on n =",x$fit$n, " observations on d =",x$fit$p, " features.",fill=TRUE)
}

summary.lassoscore <- function(object,...){
  cat("An object of class `lassoscore'\nbased on n =",object$fit$n, " observations on d =",object$fit$p, " features, with",sum(object$fit$beta !=0),"non-zero coefficients in regression of `y' on `X'",fill=TRUE)
  cat("\nAsymptotic p-values:\n")
      print(summary(object$p.asm))
  cat("\nConservative p-values:\n")
    print(summary(object$p.cons))
  cat("\nSandwich (model-agnostic) p-values:\n")
  summary(object$p.sand)
      }

plot.lassoscore <- function(x, type=c("asm","cons","sand","all"),cone=TRUE, log=TRUE, add=FALSE,...){
	switch(type[1],
		sand = qqpval(x$p.sand,cone,log,add, main = "QQ-plot of p-values",...),
		asm = qqpval(x$p.asm,cone,log,add,main = "QQ-plot of p-values",...),
		cons = qqpval(x$p.cons,cone,log,add,main = "QQ-plot of p-values",...),
		all = { qqpval(x$p.asm,cone,log,main = "QQ-plot of p-values",...);qqpval(x$p.cons,cone,log,add=TRUE,col=2,...);qqpval(x$p.sand,cone,log,add=TRUE,col=3,...);
			legend("topleft",col=1:3,legend=c("asymptotic", "conservative", "sandwich"),pch=1)}
		)		
}

