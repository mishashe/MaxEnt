fitLikelihood <- function(s)
{
  N <- ncol(s)
  nMFT <- fitIsing(s, "nMFT")
  # Inverse Ising Inference Using All the Data PRL 108, 090201 (2012)
  f <- function(par,ro)
  {
    hr <- par[1]
    Jr <- par[2:length(par)]
    Jr[ro] <- 0
    ex <- exp(-2*s[,ro]*(hr + s %*% Jr))
    LH <- sum(log(1+ex))
    return(LH)
  }
  gr <- function(par,ro)
  {
    hr <- par[1]
    Jr <- par[2:length(par)]
    Jr[ro] <- 0
    ex <- exp(-2*s[,ro]*(hr + s %*% Jr))
    gr <- par*0
    gr[1] <- sum(-2*s[,ro]*ex/(1+ex))
    gr[2:length(gr)] <- sapply(1:N,function(i){sum(-2*s[,ro]*s[,i]*ex/(ex+1))})
    gr[ro+1] <- 0
    return(gr)
  }
  parameters <- foreach(r = 1:N, .combine = rbind, .init = matrix(0, nrow=0, ncol=N+1), .inorder=TRUE) %dopar%
  {
    par_nMFT <- c(nMFT[[1]][r],nMFT[[2]][r,])
    par <- par_nMFT*0
    #par <- rnorm(N+1)
    # par <- optim(par, f, gr = gr, ro=r ,method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent")[2], lower = -Inf, upper = Inf,control = list(reltol=1e-15,abstol=1e-15,maxit=1e5,ndeps=1e-8), hessian = FALSE)$par
    # par <- maxNR(fn=f, grad = gr, hess = NULL, start=par, ro=r)$estimate
    # par <- nlm(f=f, p=par, ro=r, hessian = FALSE,
    #            fscale = 1, print.level = 0, ndigit = 12, gradtol = 1e-10,
    #            steptol = 1e-10, iterlim = 1e5, check.analyticals = TRUE)$estimate
    par <- nlminb(start=par, objective=f, gradient = gr, hessian = NULL, ro=r, control = list(iter.max=1e5, eval.max=1e5, x.tol=1e-14, rel.tol=1e-14), lower = -1-3*abs(par_nMFT), upper = 1+3*abs(par_nMFT))$par
    return(par)
  }
  h_est <- parameters[, 1]
  J_est <- parameters[, 2:(N+1)]
  diag(J_est) <- 0
  J_est <- (J_est+t(J_est))/2
  return(list(h_est,J_est))
}


fitIsing <- function(mut_matrix, method)
{
  N <- dim(mut_matrix)[2]
  if (method=="NoCor")
  {
    m <- colMeans(mut_matrix)
    h <- atanh(m)
    J <- matrix(0, ncol(mut_matrix), ncol(mut_matrix))
  }
  else if (method=="nMFT")
  {
    cov <- cov(mut_matrix)
    cov[is.na(cov)] <- 0
    cov <- (cov+t(cov))/2
    cov_inv <- ginv(cov)
    J <- -cov_inv
    J <- (J+t(J))/2
    diag(J) <- 0
    m <- colMeans(mut_matrix)
    m <- m + (m==-1)*1/dim(mut_matrix)[1] - (m==1)*1/dim(mut_matrix)[1]
    h <- atanh(m) - J %*% m
  }
  else if (method=="TAP")
  {
    browser()
    cov <- cov(mut_matrix)
    cov[is.na(cov)] <- 0
    cov <- (cov+t(cov))/2
    cov_inv <- ginv(cov)
    m <- colMeans(mut_matrix)
    J <- matrix(0,nrow=N,ncol=N)
    for (i in 1:N)    {
      for (j in (i+1):N)      {
        J[i,j] <- (-1+sqrt(1-8*cov_inv[i,j]*m[i]*m[j]))/4/m[i]/m[j]
        J[j,i] <- J[i,j]
      }
    }
    h <- atanh(m) - J %*% m - (1-m^2) %*% J %*% m
  }
  return(list(h,J))
}


sampleIsing <- function(h, J, nSamples, T)
{
  N <- length(h)
  if (all(J==0))
  {
    mut_matrix_combine <- foreach (z = 1:nSamples, .init = matrix(-1, nrow = 0, ncol = N), .combine=rbind, .inorder=FALSE) %dopar%
    {
      t(-1+2 * (runif(N, 0, 1) < exp(h)/(exp(h)+exp(-h))))
    }
    return(mut_matrix_combine)
  }
  
  
  mut_matrix_combine <- foreach (z = 1:nSamples, .init = matrix(-1, nrow = 0, ncol = N), .combine=rbind, .inorder=FALSE) %dopar%
  {
    # mut_matrix <- rep(-1, N)
    mut_matrix <- -1+2*(rand(1,N)<0.5)
    t <- 0
    while (t<T)
    {
      Ind <- sample(1:N,1)
      dE <- ( 2*h[Ind] + 2* mut_matrix %*% J[,Ind] ) * mut_matrix[Ind] 
      if (dE < 0 )
      {
        mut_matrix[Ind] <- -mut_matrix[Ind]
        t <- t+1
      }
      else if (runif(1, 0, 1) < exp(-dE))
      {
        mut_matrix[Ind] <- -mut_matrix[Ind]
        t <- t+1
      }
    }
    mut_matrix
  }
  return(mut_matrix_combine)
  # mut_matrix_dcis_sample <- foreach(t = 1:10, .combine = rbind, .init = matrix(0, nrow=0, ncol=N_dcis), .inorder=FALSE)%dopar%
  # {
  #   IsingSampler(100, J_dcis, h_dcis, beta = 1, nIter = 500, responses = c(-1L, 1L), method = c("MH", "CFTP", "direct")[1])
  # }
}
