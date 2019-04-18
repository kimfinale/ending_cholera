function (p, n, init, scale = rep(1, length(init)), adapt = !is.null(acc.rate), 
          acc.rate = NULL, gamma = 2/3, list = TRUE, showProgressBar = interactive(), 
          n.start = 0, ...) 
{
  if (adapt & !is.numeric(acc.rate)) 
    stop("Argument \"acc.rate\" is missing!")
  if (gamma <= 0.5 | gamma > 1) 
    stop("Argument \"gamma\" must be in (0.5, 1]!")
  if (is.numeric(adapt)) 
    n.adapt <- adapt
  if (adapt == TRUE) 
    n.adapt <- Inf
  if (adapt == FALSE) 
    n.adapt <- 0
  d <- length(init)
  X <- matrix(NA, ncol = d, nrow = n)
  colnames(X) <- names(init)
  X[1, ] <- init
  p.val <- rep(NA, n)
  val <- p(X[1, ], ...)
  if (is.list(val)) {
    returns.list <- TRUE
    extras <- list()
    if (!"log.density" %in% names(val)) {
      stop("The list returned by 'p' must contain an element named 'log.density!'")
    }
    if (length(val$log.density) > 1) 
      stop("The list element 'log.density' must be a scalar value!")
    p.val[1] <- val$log.density
    extras[[1]] <- val["log.density" != names(val)]
  }
  else {
    returns.list <- FALSE
    if (length(val) > 1) 
      stop("The function 'p' must return a scalar value or a named list! See ?MCMC.!")
    p.val[1] <- val
  }
  if (d > 1) {
    if (length(scale) == d) {
      M <- diag(scale)
    }
    else {
      M <- scale
    }
  }
  else {
    M <- matrix(scale)
  }
  if (ncol(M) != length(init)) 
    stop("Length or dimension of 'init' and 'scale' do not match!")
  S <- t(chol(M))
  cat("  generate", n, "samples \n")
  if (showProgressBar) {
    pb <- txtProgressBar(min = 0, max = n, style = 3)
  }
  update.step <- max(5, floor(n/100))
  k <- 0
  for (i in 2:n) {
    if (showProgressBar && i%%update.step == 0) {
      setTxtProgressBar(pb, i)
    }
    U <- rt(d, df = d)
    X.prop <- c(X[i - 1, ] + S %*% U)
    names(X.prop) <- names(init)
    val <- p(X.prop, ...)
    if (returns.list) {
      p.val.prop <- val$log.density
      extras.prop <- val["log.density" != names(val)]
    }
    else {
      p.val.prop <- val
    }
    alpha <- min(1, exp(p.val.prop - p.val[i - 1]))
    if (!is.finite(alpha)) 
      alpha <- 0
    if (runif(1) < alpha) {
      X[i, ] <- X.prop
      p.val[i] <- p.val.prop
      if (returns.list) {
        extras[[i]] <- extras.prop
      }
      k <- k + 1
    }
    else {
      X[i, ] <- X[i - 1, ]
      p.val[i] <- p.val[i - 1]
      if (returns.list) {
        extras[[i]] <- extras[[i - 1]]
      }
    }
    ii <- i + n.start
    if (ii < n.adapt) {
      adapt.rate <- min(1, d * ii^(-gamma))
      M <- S %*% (diag(d) + adapt.rate * (alpha - acc.rate) * 
                    U %*% t(U)/sum(U^2)) %*% t(S)
      eig <- eigen(M, only.values = TRUE)$values
      tol <- ncol(M) * max(abs(eig)) * .Machine$double.eps
      if (!isSymmetric(M) | is.complex(eig) | !all(Re(eig) > 
                                                   tol)) {
        M <- as.matrix(Matrix::nearPD(M)$mat)
      }
      S <- t(chol(M))
    }
  }
  if (showProgressBar) {
    close(pb)
  }
  acceptance.rate <- round(k/(n - 1), 3)
  if (list) {
    res <- list(samples = X, log.p = p.val, cov.jump = M, 
                n.sample = n, acceptance.rate = acceptance.rate, 
                adaption = adapt, sampling.parameters = list(sample.density = p, 
                                                             acc.rate = acc.rate, gamma = gamma))
    if (returns.list) {
      res$extra.values = extras
    }
    return(res)
  }
  else {
    cat("Acceptance rate:", acceptance.rate, "\n")
    return(X)
  }
}
<bytecode: 0x000000004271d828>
  <environment: namespace:adaptMCMC>