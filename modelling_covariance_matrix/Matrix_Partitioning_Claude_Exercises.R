## ============================================================
## Matrix partitioning: practice exercises
##
## Fill in each TODO, then run the check() line below it.
## Solutions are at the bottom of the file.
##
## Everything here is indexing + matrix multiplication.
## No calculus, no optimisation.
## ============================================================

check <- function(label, condition) {
  cat(sprintf("%-45s %s\n", label, if (isTRUE(condition)) "PASS" else "FAIL"))
}


## ------------------------------------------------------------
## 0. Setup: simulate 5 variables from a 2-factor structure
## ------------------------------------------------------------
## Vars 1-2 load on factor 1, vars 4-5 on factor 2,
## var 3 is a "bridge" loading on both.
## Nothing below depends on knowing how this was generated.

set.seed(42)
n <- 5000
F1 <- rnorm(n)
F2 <- rnorm(n)

Y <- cbind(
  V1 = 0.80 * F1              + rnorm(n, sd = sqrt(1 - 0.64)),
  V2 = 0.70 * F1              + rnorm(n, sd = sqrt(1 - 0.49)),
  V3 = 0.60 * F1 + 0.40 * F2  + rnorm(n, sd = sqrt(1 - 0.36 - 0.16)),
  V4 =              0.75 * F2 + rnorm(n, sd = sqrt(1 - 0.5625)),
  V5 =              0.65 * F2 + rnorm(n, sd = sqrt(1 - 0.4225))
)

S <- cov(Y)          # 5x5 observed covariance matrix
p <- ncol(S)

round(S, 3)


## ------------------------------------------------------------
## 1. Extracting blocks
## ------------------------------------------------------------
## Partition with group 1 = (V1, V2) and group 2 = (V3, V4, V5).
## So p1 = 2, p2 = 3.

i1 <- 1:2
i2 <- 3:5

S11 <- S[i1, i1]
S12 <- S[i1, i2]
# TODO: fill in S21 and S22
S21 <- S[i2, i1]
S22 <- S[i2,i2]

check("1a. S21 is the transpose of S12", isTRUE(all.equal(S21, t(S12))))
check("1b. S22 is 3x3",                  identical(dim(S22), c(3L, 3L)))


## ------------------------------------------------------------
## 2. Dimension arithmetic
## ------------------------------------------------------------
## Before running anything: write down on paper what dimension
## each of these should be. Then check.

check("2a. S11 is p1 x p1", identical(dim(S11), c(2L, 2L)))
check("2b. S12 is p1 x p2", identical(dim(S12), c(2L, 3L)))
check("2c. S21 is p2 x p1", identical(dim(S21), c(3L, 2L)))

## TODO: what dimension is S12 %*% solve(S22) ? S12 is p1 x p2, solve(S22) is p2 x p2 --> S12 %*% solve(S22) is p1 x p2 (2 x 3)
## Predict first, then uncomment:
dim(S12 %*% solve(S22))

## TODO: which of these two is a valid multiplication? Try both.
## One will error. Understand WHY before moving on. 

# Because of matrix multilication rule, number of columns in the first matrix must equal to the number of rows in the second matrix.

S12 %*% S22
# S22 %*% S12


## ------------------------------------------------------------
## 3. Regression coefficients live in the blocks of S
## ------------------------------------------------------------
## Claim: regressing group 1 on group 2 gives coefficients
##        B = S12 %*% solve(S22)
## Each ROW of B is one outcome variable's set of slopes.

# TODO: compute B
B <- S12 %*% solve(S22)

## Compare row 1 against an actual regression of V1 on V3,V4,V5.
fit1 <- lm(Y[, 1] ~ Y[, 3] + Y[, 4] + Y[, 5])
coef_lm <- unname(coef(fit1)[-1])   # drop the intercept

check("3a. B is 2x3", identical(dim(B), c(2L, 3L)))
check("3b. row 1 of B equals lm() slopes",
      isTRUE(all.equal(as.numeric(B[1, ]), coef_lm)))

## TODO: verify row 2 of B against a regression of V2 on V3,V4,V5.
## Write the lm() yourself.
lm(data = data.frame(Y), V2 ~ V3 + V4 + V5)
B

# Do the same for row 1
lm(data = data.frame(Y), V1 ~ V3 + V4 + V5)
B
## ------------------------------------------------------------
## 4. Schur complement = conditional covariance
## ------------------------------------------------------------
## Cov(group 1 | group 2) = S11 - S12 %*% solve(S22) %*% S21

# TODO: compute it
Schur <- S11 - S12 %*% solve(S22) %*% S21

check("4a. Schur is p1 x p1", identical(dim(Schur), c(2L, 2L)))

## The diagonal entries are residual variances. Check against lm:
resid_var1 <- var(residuals(fit1)) * (n - 1) / (n - 1)   # sample variance

check("4b. Schur[1,1] is close to residual variance of fit1",
      abs(Schur[1, 1] - var(residuals(fit1))) < 0.01)

## TODO: Schur[1,1] should be SMALLER than S11[1,1]. Why?
## Answer in one sentence in a comment here:
## Don't understand what the point of Schur component is.


## ------------------------------------------------------------
## 5. The trap: blocks of the inverse are not inverses of blocks
## ------------------------------------------------------------

K <- solve(S)        # precision matrix
K11 <- K[i1, i1]     # the "superscript" block, Sigma^{11}

check("5a. K11 equals inverse of the Schur complement",
      isTRUE(all.equal(K11, solve(Schur), check.attributes = FALSE)))

check("5b. K11 does NOT equal solve(S11)",
      !isTRUE(all.equal(K11, solve(S11), check.attributes = FALSE)))

## Look at both and see how far apart they are:
round(K11, 3)
round(solve(S11), 3)


## ------------------------------------------------------------
## 6. Non-contiguous partitions
## ------------------------------------------------------------
## Groups do not have to be adjacent columns. Put V1 and V4 in
## group 1, everything else in group 2.

j1 <- c(1, 4)
j2 <- c(2, 3, 5)

# TODO: extract all four blocks with these index vectors
T11 <- S[j1, j1]
T12 <- S[j1, j2]
T21 <- S[j2, j1]
T22 <- S[j2, j2]

check("6a. T12 is 2x3", identical(dim(T12), c(2L, 3L)))
check("6b. T11 diagonal is variances of V1 and V4",
      isTRUE(all.equal(diag(T11), c(S[1, 1], S[4, 4]), check.attributes = FALSE)))

## TODO: compute the regression of (V1, V4) on (V2, V3, V5)
## using the same formula as in section 3.
B2 <- T12 %*% solve(T22)

# Check
lm(data = data.frame(Y), V1 ~ V2 + V3 + V5)
lm(data = data.frame(Y), V4 ~ V2 + V3 + V5)

## ------------------------------------------------------------
## 7. From precision matrix to partial correlations
## ------------------------------------------------------------
## This is the step that turns K into a network.
## Partial correlation between i and j, given all others:
##      -K[i,j] / sqrt(K[i,i] * K[j,j])

# TODO: build the full matrix of partial correlations.
# Hint: D <- diag(1 / sqrt(diag(K))) then combine with K,
# and set the diagonal to 1.
D <- diag(1 / sqrt(diag(K)))
pcor <- -D %*% K %*% D

diag(pcor) <- 1

p_cor_boong <- -cov2cor(K)
check("7a. pcor is 5x5",        identical(dim(pcor), c(5L, 5L)))
check("7b. diagonal is all 1",  isTRUE(all.equal(diag(pcor), rep(1, 5),
                                                 check.attributes = FALSE)))
check("7c. pcor is symmetric",  isTRUE(all.equal(pcor, t(pcor))))

## Compare the marginal correlations against the partial ones:
round(cov2cor(S), 2)
round(pcor, 2)

## TODO: V1 and V5 have no direct link in the generating model
## (V1 loads only on F1, V5 only on F2). Look at cor(V1,V5)
## versus pcor(V1,V5). Which is near zero? Why is the other not?


## ------------------------------------------------------------
## 8. The GGM parameterisation
## ------------------------------------------------------------
## Delta = diag(K)^(-1/2), Omega = I - Delta %*% K %*% Delta
## Then Sigma should come back as Delta %*% solve(I - Omega) %*% Delta

# TODO: build Delta and Omega
Delta <- diag(diag(K) ^ (-1/2))
Omega <- diag(nrow = 5) - Delta %*% K %*% Delta

check("8a. Omega has zero diagonal",
      isTRUE(all.equal(diag(Omega), rep(0, 5), check.attributes = FALSE)))
check("8b. Omega equals pcor off the diagonal",
      isTRUE(all.equal(Omega + diag(5), pcor, check.attributes = FALSE)))

## TODO: reconstruct S from Delta and Omega and check it round-trips.
S_implied <- Delta %*% solve(diag(5) - Omega) %*% Delta
check("8c. reconstruction recovers S",
      isTRUE(all.equal(S_implied, S, check.attributes = FALSE)))



## ------------------------------------------------------------
## 9. Choose group 1 to be only 1 variable, and group 2 to be the rest
## ------------------------------------------------------------

l1 <- c(1)
l2 <- c(2,3,4,5)

V11 <- S[l1,l1]
V12 <- S[l1,l2]
V21 <- S[l2,l1]
V22 <- S[l2,l2]

partial_variance <- V12 %*% solve(V22)
partial_variance






