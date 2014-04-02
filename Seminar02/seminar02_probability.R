# SE = Standard Error
rnorm(10)
rnorm(mean=100000, sd=2, n=10)

n <- 10
B <- 4

x <- matrix(data=rnorm(n*B), nrow=n)

rownames(x) <- (sprintf("obs%02d",1:n))
colnames(x) <- (sprintf("sam%02d",1:B))

colMeans(x)

mean(colMeans(x))

B <- 1000
x10 <- matrix(rnorm(10 * B), nrow = 10)
x100 <- matrix(rnorm(100 * B), nrow = 100)
x1000 <- matrix(rnorm(1000 * B), nrow = 1000)
x10000 <- matrix(rnorm(10000 * B), nrow = 10000)
xBar10 <- colMeans(x10)
xBar100 <- colMeans(x100)
xBar1000 <- colMeans(x1000)
xBar10000 <- colMeans(x10000)
xBarSd10 <- sd(colMeans(x10))
xBarSd100 <- sd(colMeans(x100))
xBarSd1000 <- sd(colMeans(x1000))
xBarSd10000 <- sd(colMeans(x10000))
IQR10 <- IQR(colMeans(x10))
IQR100 <- IQR(colMeans(x100))
IQR1000 <- IQR(colMeans(x1000))
IQR10000 <- IQR(colMeans(x10000))

cbind(sampSize = c(10, 100, 1000, 10000), trueSEM = 1 / sqrt(c(10, 100, 1000, 10000)), obsSEM = c(xBarSd10, xBarSd100, xBarSd1000, xBarSd10000), sampMeanIQR = c(IQR10, IQR100, IQR1000, IQR10000))
# Do the same as IQR() for mad()

x <- matrix(data=rnorm(40*100), nrow=40)
threshold  <- 0.1
mean(x <= threshold & x >= 0)
pnorm(threshold)

# x <- matrix(data=runif(4*10), nrow=4)
# threshold  <- 0.1
# mean(x <= threshold)
# punif(threshold)

min.x <- -5

max.x <- 5

num.samples <- 1000

x <- seq(from = min.x, to = max.x, length = num.samples)

# Open new blank plot with x limits from -5 to 5, and y limits from 0 to 1
plot(c(-5, 5), c(0, 1), xlab = "x", ylab = "f(x)", main = "Normal probability density function", 
     type = "n")

# Add each density plot one at a time
lines(x, dnorm(x, mean = 0, sd = 0.5), lwd = 2, col = "red")

lines(x, dnorm(x, mean = 0, sd = 1), lwd = 2, col = "green")

lines(x, dnorm(x, mean = 0, sd = 2), lwd = 2, col = "blue")

lines(x, dnorm(x, mean = -2, sd = 1), lwd = 2, col = "magenta")

# We can also add a legend to the plot
legend("topright", c("mean=0, sd=0.5", "mean=0, sd=1", "mean=0, sd=2", "mean=-2, sd=1"), 
       col = c("red", "green", "blue", "magenta"), lty = 1, lwd = 2)

normal.mean <- c(0, 0, 0, -2)

normal.sd <- c(0.5, 1, 2, 1)

colors <- c("red", "green", "blue", "magenta")

# Open new plot with x limits from -5 to 5, and y limits from 0 to 1
plot(c(-5, 5), c(0, 1), xlab = "x", ylab = "f(x)", main = "Normal probability density function", 
     type = "n")

# Add density plots with a for loop
for (i in 1:length(normal.mean)) {
  lines(x, dnorm(x, mean = normal.mean[i], sd = normal.sd[i]), lwd = 2, col = colors[i])
}

# Add a legend to the plot
legend("topright", paste0("mean=", normal.mean, ", sd=", normal.sd), col = colors, 
       lty = 1, lwd = 2)

# Open new plot with x limits from -5 to 5, and y limits from 0 to 1
plot(c(-5, 5), c(0, 1), xlab = "x", ylab = "f(x)", main = "Normal probability density function", 
     type = "n")

# Create our own user-defined function for plotting Normal probability
# density function
f <- function(col, ...) {
  lines(x, dnorm(x, ...), col = col, lwd = 2)
}

# apply this function with different parameters
plot.status <- mapply(f, mean = normal.mean, sd = normal.sd, col = colors)

# Add a legend to the plot
legend("topright", paste0("mean=", normal.mean, ", sd=", normal.sd), col = colors, 
       lty = 1, lwd = 2)

set.seed(1)

normal.mean <- 1
normal.sd <- 1

rnorm(n = 5, mean = normal.mean, sd = normal.sd)


# Draw a sample
set.seed(1)

y <- rnorm(n = 1000, mean = normal.mean, sd = normal.sd)

# Estimate sample density
estimated.density <- density(y)

# Plot the estimated density
plot(estimated.density, col = "blue", lwd = 2)

# Add the data points under the density plot and colors them organge
rug(y, col = "orange")

# Plot the true density lines
x <- seq(from = -5, to = 5, length = 1000)

true.density <- dnorm(x, mean = normal.mean, sd = normal.sd)

lines(x, true.density, col = "red", lwd = 2)

legend("topright", c("sample density", "true density"), col = c("blue", "red"), 
       lty = 1, lwd = 2)