library(MASS)
library(gplots)

asd <- matrix(rnorm(1000),ncol=10)

eg <- eigen(cov(asd))

round(t(eg$vectors)%*%eg$vectors,5)

round(t(eg$vectors)%*%cov(asd)%*%eg$vectors,5)
eg$values

det(t(eg$vectors))
det(eg$vectors)

det(diag(eg$values))

sum(diag(ginv(cov(asd))))

prod(eg$values)
det(cov(asd))

dim(eg$vectors)
a <- rnorm(10)
B <- eg$vectors%*%diag(a)%*%t(eg$vectors)
a
1/a
sum(diag(ginv(B)))
sum(1/a)


X <- matrix(c(rep(1,24),rep(0,12),rep(1,12),rep(0:1,12)),ncol=3)
gamma <- matrix(rnorm(30),ncol=3)
X[3,]%*%gamma

A <- matrix(0,nrow=10,ncol=10)
for(j in 1:10){
    M <- matrix(0,nrow=10,ncol=10)
    M[j,j] <- 1
    v <- rep(0,10)
    v[j] <- 1
    A <- A + M%*%gamma%*%X[j,]%*%v
}
A


sdf <- eg$vectors[,1]
wer <- eg$vectors[,2]
sdft <- sdf%*%t(sdf)
wert <- wer%*%t(wer)
det(sdft)
det(wert + sdft)

########################
Q <- eg$vectors
answer <- ginv(Q%*%diag(a)%*%t(Q))
ginv(t(Q))%*%ginv(diag(a))%*%ginv(Q)

ginv(t(Q))%*%(ginv(Q)%*%ginv(diag(a)))


asd <- matrix(rnorm(100),ncol=10)
asd2 <- matrix(rnorm(100),ncol=10)
ginv(asd%*%asd2)
ginv(asd2)%*%ginv(asd
                  )



#####################################
gi <- matrix(rnorm(100),ncol=10)
sigma <- lapply(1:10, function(x){matrix(rnorm(100),ncol=10)})
trace <- function(x){sum(diag(x))}

sum(sapply(1:10, function(i){
    trace(t(gi[,i])%*%sigma[[i]]%*%(gi[,i]))
    t(gi[,i])%*%sigma[[i]]%*%gi[,i]
}))
t(gi)%*%sigma%*%gi
trace(t(gi)%*%sigma%*%gi)
trace(gi%*%t(gi)%*%sigma)

z <- rnorm(10)
trace(diag(z))
eg <- eigen(cov(sigma[[1]]))
trace(eg$vectors%*%diag(z)%*%t(eg$vectors))
      

##########################
numGenes <- 100
numSamples <-60
mu <- rnorm(numGenes,mean = 9)
batches <- c(rep(0,numSamples/2),rep(1,numSamples/2))
block1 <- c(rep(1,40),rep(0,60))
block2 <- c(rep(1,10),rep(0,90))
block3 <- c(rep(0,80),rep(.7,20))
block4 <- c(rep(0,40),rep(.9,30),rep(0,30))

batch1 <- cbind(block1, block2, block3, block4)
batch2 <- cbind(block1, block2, block3)      
Sigma1 <- batch1%*%t(batch1)
Sigma2 <- batch2%*%t(batch2)
diag(Sigma1) <- 4
diag(Sigma2) <- 4
heatmap.2(Sigma1, trace = "none")
heatmap.2(Sigma2, trace = "none")

data <- cbind(t(mvrnorm(sum(batches),mu=mu, Sigma = Sigma)),t(mvrnorm(sum(1-batches),mu=mu, Sigma = Sigma)))
heatmap.2(data, trace = "none", col = "bluered")

G_star <- data-rowMeans(data)
G_standard <- (G_star/sqrt(rowSums(G_star^2)))

eigenG <- eigen(G_standard%*%t(G_standard))
Q <- eigenG$vectors
D <- diag(eigenG$values)
plot(diag(D))
plot(Q[,1])
plot(Q[,2])
plot(Q[,3])
X <- cbind(rep(1,10))
X <- cbind(rep(1,10),c(rep(1,5),rep(0,5)))
X <- cbind(rep(1,10),c(rep(1,5),rep(0,5)), rnorm(10))











#########################################


G <- matrix(rnorm(1000),ncol=10)
G_star <- G-rowMeans(G)
G_standard <- (G_star/sqrt(rowSums(G_star^2)))

eigenG <- eigen(G_standard%*%t(G_standard))
Q <- eigenG$vectors
D <- diag(eigenG$values)
X <- cbind(rep(1,10))
X <- cbind(rep(1,10),c(rep(1,5),rep(0,5)))
X <- cbind(rep(1,10),c(rep(1,5),rep(0,5)), rnorm(10))

hat <- ginv(t(X)%*%X)%*%t(X)
hat%*%(2:11)
a <- Reduce('+', lapply(1:6, function(i){
    g <- G_star[,i,drop=F]
    ginv(Q)%*%(g)%*%t(g)%*%t(ginv(Q))
}))


est1 <- ginv(Q)%*%(G_star)%*%diag(hat[1,])%*%t(G_star)%*%t(ginv(Q))
est2 <- ginv(Q)%*%(G_star)%*%diag(hat[2,])%*%t(G_star)%*%t(ginv(Q))
est3 <- ginv(Q)%*%(G_star)%*%diag(hat[3,])%*%t(G_star)%*%t(ginv(Q))
gamma1 <- diag(est1)
gamma2 <- diag(est2)

fitted1 <- diag((X%*%t(cbind(gamma1,gamma2)))[1,])
fitted2 <- diag((X%*%t(cbind(gamma1,gamma2)))[2,])

# b <- ginv(Q)[2,,drop=F]%*%(G_star)%*%t(G_star)%*%t(ginv(Q)[2,,drop=F])
# ginv(Q)[1,,drop=F]%*%(G_star)%*%t(G_star)%*%t(ginv(Q)[1,,drop=F])
round(a[1:10,1:10],3)
round(b[1:10,1:10],3)
round(D[1:10,1:10],3)
round(est[1:10,1:10],3)
round(fitted1[1:10,1:10],3)
round(fitted2[1:10,1:10],3)

(Q%*%cov(t(G))%*%t(Q))[1:10,1:10]

Q%*%t(Q)

A <- cov(matrix(rnorm(1000),ncol=10))
B <- cov(matrix(rnorm(1000),ncol=10))
D <- diag(rnorm(10))

A%*%D%*%B
A%*%t(D)%*%t(B)
A%*%t(B%*%D)
