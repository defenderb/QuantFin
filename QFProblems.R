mu <- c(0.08, 0.10, 0.13, 0.15, 0.20)
Sigma <- matrix(c(0.019600,-0.007560, 0.012880, 0.008750, -0.009800,
                  -0.007560, 0.032400, -0.004140, -0.009000, 0.009450,
                  0.012880, -0.004140, 0.052900, 0.020125, 0.020125,
                  0.008750, -0.009000, 0.020125, 0.062500, -0.013125,
                  -0.009800, 0.009450, 0.020125, -0.013125, 0.12250), 5,5)
sigmaP2 <- 0.0625
G <- function(x, mu, Sigma, sigmaP2)
{
  n <- length(mu)
  c(mu + rep(x[n+1], n) + 2*x[n+2]*(Sigma %*% x[1:n]),
    sum(x[1:n]) - 1,
    t(x[1:5]) %*% Sigma %*% x[1:5] - sigmaP2)
}

DG <- function(x, mu, Sigma, sigmaP2)
{
  n <- length(mu)
  grad <- matrix(0.0, n+2, n + 2)
  grad[1:n, 1:n] <- 2*x[n+2]*Sigma
  grad[1:n, n+1] <- 1
  grad[1:n, n+2] <- 2*(Sigma %*% x[1:5])
  grad[n+1, 1:n] <- 1
  grad[n+2, 1:n] <- 2*t(x[1:5]) %*% Sigma
  grad
}
x <- c(rep(0.5, 5), 1, 1)
for(i in 1:25)
  x <- x - solve(DG(x, mu, Sigma, 0.25^2),
                 G(x, mu, Sigma, 0.25^2))
x[1:length(mu)]
# risk minimization ===============================
G2 <- function(x, mu, Sigma, mu2)
  {
    n <- length(mu)
    c(2*(t(x[1:n]) %*% Sigma) +rep(x[n+1], n) + x[n+2]* mu,
      sum(x[1:n]) - 1,
      t(mu) %*% x[1:n] - mu2)
  }
DG2 <- function(x, mu, Sigma, mu2)
{
  n <- length(mu)
  grad <- matrix(0.0, n + 2, n + 2)
  grad[1:n, 1:n] <- 2*Sigma
  grad[1:n, n+1] <- 1
  grad[1:n, n+2] <- mu
  grad[n+1, 1:n] <- 1
  grad[n+2, 1:n] <- mu
  grad
}
mu2 <- 0.1991514
x <- c(rep(0.5, 5), 1, 1)
for(i in 1:25)
  x <- x - solve(DG2(x, mu, Sigma, mu2),
                 G2(x, mu, Sigma, mu2))
x[1:length(mu)]
