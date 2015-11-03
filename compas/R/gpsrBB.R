
gpsrBB <- function(y, A, dimx, tauVec, cont_steps = 4, initAty = TRUE, tolA = 0.0001, maxiter = 10000)
{
	stopCriterion <- 1			# 1 => check objective, 2 => check gradient
	found <- FALSE	
	
	# Barzilai-Borwein specific parameters	
	alphamin <- 1e-20
	alphamax <- 1e20
	
	Aty <- t(A) %*% y		# we will use this a lot, so we precompute

	matX <- matrix(0, dimx, 1)		# this will hold successive solutions
	
	x <- c(rep(0, dimx))		# default x	
	if(initAty)					# initialize x at A' * y
		x <- Aty
		
	#################################################################
	## NOTE: Since our vector x is known to be non-negative we 
	## will only use vector u, and omit vector v in our calculations
	#################################################################
	
	u <- x	
	
	ttau <- 0.5 * tauVec
	initialTau <- ttau
		
	cont_loop <- 0
	alpha <- 1.0
	keep_continuation <- TRUE
	while(keep_continuation)
	{
		cont_loop <- cont_loop + 1
		if(cont_loop > cont_steps)
			keep_continuation <- FALSE
		
		residual <- y - A%*%x
		
		gradq <- t(A) %*% residual				
		ttau <- pmax(initialTau, 0.5*max(abs(gradq)))
		
		diff_tau <- ttau - initialTau
		if(sqrt(sum(diff_tau * diff_tau)) < tolA)			# is the difference close to zero	
			keep_continuation <- FALSE		# this is adaptive continuation
				
		f <- 0.5*(t(residual) %*% residual) + sum(ttau * u)		

		resid_base <- y - residual
		
		iter <- 0
		keep_going <- TRUE		
		while(keep_going)
		{
			gradu <- t(A)%*%resid_base - Aty + ttau
						
			du <- pmax(u - alpha*gradu, 0.0) - u
			dx <- du
			old_u <- u
			
			auv <- A %*% dx
			dGd <- t(auv) %*% auv
			
			lambda <- 1			# default value								
			u <- old_u + lambda*du
			x <- u
			
			residual <- y - resid_base - lambda*auv
			prev_f <- f
			f <- 0.5*(t(residual) %*% residual) + sum(ttau * u)	

			# compute new alpha
			dd <- t(du) %*% du
			if(dGd <= 0)
				alpha <- alphamax			
			
			alpha <- min(alphamax, max(alphamin, dd/(1e-20 + dGd)))
				
			resid_base <- resid_base + lambda*auv
						
			iter <- iter + 1			
			
			if(stopCriterion == 1)	# stopping criterion based on the relative change in objective function
				keep_going <- (abs(f-prev_f)/(prev_f) > tolA)					
				
			if(stopCriterion == 2)	# stopping criterion based on relative norm of step taken						
				keep_going <- (sqrt(sum(dx^2))/sqrt(sum(x^2)) > tolA)			
			
			if(keep_going == FALSE)		# found a good solution
			{
				found <- TRUE
				if(cont_loop == 1)		# the first run
					matX[, 1] <- x					
				else
					matX <- cbind(matX, x)				
			}
				
			if(iter > maxiter)		# reached the maximum, so stop			
				keep_going <- FALSE				
		}		
	}	
	return(list(solFound=found, solMatrix=matX))
}

