
gpsrBasic <- function(y, A, dimx, tauVec, cont_steps = 3, initAty = TRUE, tolA = 0.0001, maxiter = 10000)
{	
	stopCriterion <- 1			# 1 => check objective, 2 => check gradient
	found <- FALSE
		
	mu <- 0.1					# sufficient decrease parameter for GP line search
	lambda_backtrack <- 0.5		# backtracking parameter for line search
	
	matX <- matrix(0, dimx, cont_steps)		# this matrix will hold successive solution vectors
	
	Aty <- t(A) %*% y
	
	x <- c(rep(0, dimx))		# default x	
	if(initAty)					# initialize x at Aty
		x <- Aty
	
	#################################################################
	## NOTE: Since our vector x is known to be non-negative we 
	## will only use vector u, and omit vector v in our calculations
	#################################################################
	
	u <- x			
	tauFactors <- 0.9/cont_steps * c(1:cont_steps)
	
	for(cont_loop in 1:cont_steps)
	{
		ttau <- tauFactors[cont_loop] * tauVec
		
		residual <- y - A%*%x
		f <- 0.5*(t(residual)%*%residual) + sum(ttau * u)		
		
		resid_base <- y - residual
		keep_going <- 1
		
		iter <- 0
		while(keep_going)
		{
			x_previous <- x
	
			# compute gradient						
			gradu <- (t(A) %*% resid_base) - Aty + ttau

			# set search direction			
			du <- 0 - gradu
			dx <- du
			old_u <- u

			# calculate useful matrix-vector product involving dx
			auv <- A %*% dx
			dGd <- t(auv) %*% auv
			
			lambda0 <- -(t(gradu) %*% du) / (1e-12 + dGd)
			lambda <- lambda0[1,1]		# for some reason, R thinks this is a matrix of 1-by-1
			
			while(1)	# backtrack search loop
			{
				# calculate step for this lambda and candidate point								
				du <- pmax(u - lambda*gradu, 0.0) - u
				u_new <- u + du
				
				dx <- du 
				x_new <- x + dx

				# evaluate function at the candidate point
				resid_base <- A %*% x_new
				residual <- y - resid_base
				
				f_new = 0.5*(t(residual)%*%residual) + sum(ttau * u_new)    
				
				# test sufficient decrease condition
				if(f_new <= f + mu * (t(gradu) %*% du))
					break
					
				lambda <- lambda * lambda_backtrack
			}
			
			u <- u_new
			prev_f <- f 
			f <- f_new
			x <- u

			iter <- iter + 1
			
			if(stopCriterion == 1)	# stopping criterion based on the relative change in objective function
				keep_going <- (prev_f>(tolA*tolA) && abs(f-prev_f)/(prev_f) > tolA)					
				
			if(stopCriterion == 2)	# stopping criterion based on relative norm of step taken						
				keep_going <- (x>(tolA*tolA) && sqrt(sum(dx^2))/sqrt(sum(x^2)) > tolA)			
			
			if(keep_going == FALSE)		# found a good solution
			{
				found <- TRUE
				matX[, cont_loop] <- x					
			}
				
			if(iter > maxiter)		# reached the maximum, so stop
				keep_going <- FALSE						
		}
	}
	return(list(solFound = found, solMatrix = matX))
}

