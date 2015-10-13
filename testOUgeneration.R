library(sde)


# Ornstein-Uhlenbeck process
d <- expression(0 - 20 * x)

s <- expression(1) 

time.max <- 100
N <- 10^3

y <- sde.sim(X0=0, drift=d, sigma=s, T=1, N=N, M=1)

x <- seq(0, time.max, length = N+1)


plot(x, y, main="Ornstein-Uhlenbeck", type="l")

print(var(y))


yy <- approxfun(x, y)

times <- seq(0, time.max, length = 100)

points(times, yy(times), pch = 19)


