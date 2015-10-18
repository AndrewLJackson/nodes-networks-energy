library(sde)
library(viridis)

palette(viridis(8))

# -----------------------------------------------------------------
# Ornstein-Uhlenbeck process
# -----------------------------------------------------------------
set.seed(1)
d <- expression(0 - 20 * x)

s <- expression(0.1) 

time.max <- 100
N <- 10^3

y <- sde.sim(X0=0, drift=d, sigma=s, T=1, N=N, M=1)

x <- seq(0, time.max, length = N+1)


par(mfrow=c(1,2))
plot(x, y, main="Ornstein-Uhlenbeck", type="l")

print(var(y))


yy <- approxfun(x, y)

times <- seq(0, time.max, length = 100)

points(times, yy(times), pch = 19)

# -----------------------------------------------------------------
# same Ornstein-Uhlenbeck process with larger s
# -----------------------------------------------------------------
set.seed(1)
d <- expression(0 - 20 * x)

s <- expression(0.2) 

time.max <- 100
N <- 10^3

y2 <- sde.sim(X0=0, drift=d, sigma=s, T=1, N=N, M=1)

x <- seq(0, time.max, length = N+1)


plot(x, y2, main="Ornstein-Uhlenbeck", type="l")

print(var(y))


yy <- approxfun(x, y2)

times <- seq(0, time.max, length = 100)

points(times, yy(times), pch = 19)

# -----------------------------------------------------------------
# test correlation between the two processes
print(cov(cbind(y, y2)))
print(cor(y,y2))

# -----------------------------------------------------------------
# same Ornstein-Uhlenbeck process with larger s
# -----------------------------------------------------------------

d <- expression(0 - 100 * x)

time.max <- 500
N <- 10^4

s.list <- c(expression(0.01) , expression(0.1) , expression(1) , expression(10))

results <- matrix(NA, N+1, length(s.list))

for (i in 1:length(s.list)){
  
  set.seed(1)
  
  y2 <- sde.sim(X0=0, drift=d, sigma=s.list[i], T=1, N=N, M=1)
  
  results[, i] <- y2
   
  
  
}

x <- seq(0, time.max, length = N+1)

par(mfrow=c(1,1))
matplot(x, results, type="l", col = 1:ncol(results) , lty = 1)

