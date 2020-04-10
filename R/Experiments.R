library(rbenchmark)

# Пример 1 ---------------------------------------------------------------
# Простой оператор, без истории

M1 <- cbind(c(3, 8, 0.2),
           c(-4 + 3i, 8, 0.1),
           c(9, -3, 1))
u1 <- c(1, 0.5, 1.7)
f1 <- c(1.8, 9.1, 1.3)

result11 <- GMSI(M1, f1, u1, 10^-3)
result12 <- solve(M1) %*% f1

print(result11)
cat("\n")
print(result12)
cat("\n")
print(all.equal(result11, result12))

result21 <- GMSI.history(M1, f1, u1, 10^-3)
result22 <- solve(M1) %*% f1

print(result21)
cat("\n")
print(result22)

plot(x = 1:(result21$iterations), y = abs(result21$var.hist[2,]), type = "l", lwd = 0.5, col = 3)
lines(x = 1:(result21$iterations), y = abs(result21$var.hist[1,]), type = "l", lwd = 0.5, col = 2)
lines(x = 1:(result21$iterations), y = abs(result21$var.hist[3,]), type = "l", lwd = 0.5, col = 4)
# Пример 2 ---------------------------------------------------------------
# Простой оператор без истории

M2 <- cbind(c(3, 4, 0.2),
            c(-4, 8, 0.1),
            c(9, -3, 1))
u2 <- c(1, 0.5, 1.7)
f2 <- c(1.8, 9.1, 1.3)

result21 <- GMSI.history(M2, f2, u2, 10^-3)
result22 <- solve(M2) %*% f2

print(result21)
cat("\n")
print(result22)

plot(x = 1:(result21$iterations), y = result21$var.hist[2,], type = "l", lwd = 0.5, col = 3, ylim = c(-2, 4))
lines(x = 1:(result21$iterations), y = result21$var.hist[1,], type = "l", lwd = 0.5, col = 2)
lines(x = 1:(result21$iterations), y = result21$var.hist[3,], type = "l", lwd = 0.5, col = 4)



# Пример 3 ----------------------------------------------------------------
# оператор 4х4

M3 <- cbind(c(.6, .91, -.013, .134),
            c(-.54, .01213, .341, -.231),
            c(.234, .198, .423, .013),
            c(.1231, 1.23, -.9412, .834))
u3 <- c(1, 1, 1, 1)
f3 <- c(1.8, 9.1, -5.4, 1.3)

result31 <- GMSI.history(M3, f3, u3, 10e-4)
result32 <- solve(M3) %*% f3 

print(result31$var)
print(result32)

rbenchmark::benchmark({
    result31 <- GMSI(M3, f3, u3, 10e-4)
}, replications = 1000)
rbenchmark::benchmark({
    result32 <- MatrixMethod(M3, f3) 
}, replications = 1000)
rbenchmark::benchmark({
    result33 <- Jacobi(M3, f3, u3, 10e-4) 
}, replications = 1000)


# Пример 5 ----------------------------------------------------------------
# Оператор 3x3 Jacobi
M4 <- cbind(c(1, 0.1, -0.1),
            c(0.1, 1, -0.1),
            c(-0.1, 0.1, 1))
f4 <- c(1.1, 1, 1)
u4 <- c(0, 0, 0)
lambs <- eigen(M4)$values

GMSI.history(M4,f4,u4)
Jacobi.history(M4,f4,u4)

rbenchmark::benchmark({
    result31 <- GMSI(M4, f4, u4, 10e-4)
}, replications = 1000)
rbenchmark::benchmark({
    result32 <- GDM(M4, f4, u4) 
}, replications = 1000)
rbenchmark::benchmark({
    result33 <- Jacobi(M4, f4, u4, 10e-4) 
}, replications = 1000)
rbenchmark::benchmark({
    result34 <- Chebishev(M4, f4, u4, lambs, 10L,10e-4) 
}, replications = 1000)


# Пример 6 ----------------------------------------------------------------
# Оператор 2х2
M5 <- cbind(c(0.75, 0), c(-1, 0.8))
f5 <- c(1, 1)
u5 <- c(4.4, 1.57)

rbenchmark::benchmark({
    result31 <- GMSI(M5, f5, u5, 10e-4)
}, replications = 1000)
rbenchmark::benchmark({
    result32 <- MatrixMethod(M5, f5) 
}, replications = 1000)
rbenchmark::benchmark({
    result33 <- Jacobi(M5, f5, u5, 10e-4) 
}, replications = 1000)
rbenchmark::benchmark({
    result34 <- SIM(M5, f5, u5, 10e-4) 
}, replications = 1000)
rbenchmark::benchmark({
    result35 <- IMRES(M5, f5, u5, 10e-4) 
}, replications = 1000)
rbenchmark::benchmark({
    result36 <- GDM(M5, f5, u5, 10e-8) 
}, replications = 1000)

print(result31)
print(result32)
print(result33)
print(result34)
print(result35)
print(result36)



result31 <- GMSI.history(M5, f5, u5, 10e-4)
result32 <- MatrixMethod(M5, f5) 
result33 <- Jacobi.history(M5, f5, u5, 10e-4) 

print(result31)
print(result32)
print(result33)


