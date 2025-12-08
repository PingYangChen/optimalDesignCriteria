designV2M <- function(v_design, n, d) {
  tmp <- as.vector(t(v_design[,1:d]))		
  wt <- sqrt(as.vector(v_design[,ncol(v_design)]))
  ang <- numeric(n-1)
  for (i in 1:(n-1)) {
    ang[i] <- acos(wt[i]/sqrt(sum((wt[i:n]^2))))
    if (i < (n-1)) {
      ang[i] <- pi - ang[i]
    }
  }
  m_design <- c(tmp, ang)
  return(m_design)
}

designM2V <- function(m_design, n, d) {
  # weight of support points
  ang <- m_design[(n*d + 1):length(m_design)]
  wcumsin <- wcos <- numeric(n)
  wcumsin[1] <- 1; wcumsin[2:n] <- cumprod(sin(ang))
  wcos[1:(n-1)] <- cos(ang); wcos[n] <- 1
  wt <-	(wcumsin*wcos)^2
  tmp <- cbind(matrix(m_design[1:(n*d)], n, d, byrow = TRUE), wt)
  tmp <- as.matrix(tmp[do.call(order, as.data.frame(round(tmp, 4))),])
  dimnames(tmp) <- list(paste0("obs_", 1:n), c(paste0("dim_", 1:d), "weight"))
  v_design <- tmp
  return(v_design)
}
