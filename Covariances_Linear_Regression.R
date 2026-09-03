# This script is provide some intuition on covariances and its relationship with linear regression.

# 1. Covariance between 2 variables is simply their linear regression coefficient

X <- rnorm(1000)
a <- 0.4
b <- 0.1

Y <- a * X + b + rnorm(1000, 0, 0.001)

print(a * var(X) - cov(X, Y))
print(a^2 * var(X) - var(Y))

