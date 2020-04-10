
# MuFindScript ------------------------------------------------------------

muFind.jj <- function(lambs, draw = TRUE, path = "complexPlot") {
    # Проверка типа переменной - вектора
    stopifnot(is.numeric(lambs) || is.complex(lambs))
    stopifnot(is.logical(draw))
    n <- length(lambs) # Переменная количества считанных данных
    if (n <= 1) stop("dim of linear operator is emprty or not valid")
    # Если спектр оператора всё-таки действительная ось
    Rflash <- numeric(0)
    muflash <- complex(0)
    if (is.numeric(lambs) || all(Im(lambs) == 0)) {
        mu <- (max(Re(lambs)) + min(Re(lambs))) / 2 # считаем центр по двум крайним точкам
        R <- (max(Re(lambs)) - min(Re(lambs))) / 2 # считаем радиум круга по двум крайним точкам
        if (abs(mu) <= abs(R)) stop("Окружность оператора лежит на начале координат")
        muflash <- c(muflash, mu)
        Rflash <- c(Rflash, R)
    }
    # оказывается верным является выражение 5 < 5/0
    if (n == 2) {
        mu <- .mu2(labms[1], lambs[2])
        R <- .R2(lambs[1], lambs[2])
        if (abs(mu) <= abs(R)) stop("Окружность оператора лежит на начале координат")
        muflash <- c(muflash, mu)
        Rflash <- c(Rflash, R)
    }
    if (n >= 3) {
        
        for (i in (1:(n - 1))) {
            for (j in ((i + 1):n)) {
                mu <- .mu2(lambs[i], lambs[j])
                R <- .R2(lambs[i], lambs[j])
                if (abs(mu) <= abs(R)) next
                if (any(abs(mu - lambs) > R)) next
                muflash <- c(muflash, mu)
                Rflash <- c(Rflash, R)
            }
        }
        
        for (i in (1:(n - 2))) {
            for (j in ((i + 1):(n - 1))) {
                for (k in ((j + 1):n)) {
                    mu <- .mu3(lambs[i], lambs[j], lambs[k])
                    R <- .R3(mu, lambs[i])
                    if (abs(mu) <= abs(R)) next
                    if (any(abs(mu - lambs) > R)) next
                    muflash <- c(muflash, mu)
                    Rflash <- c(Rflash, R)
                }
            }
        }
    }
    mu <- muflash[which.min(Rflash)]
    R <- min(Rflash)
    if (draw) .drawComplex(lambs = lambs, mu = mu, R = R, fileName = path)
    return(mu)
}

lambs <- c(5 + 5i, 10 + 10i)
print(muFind(lambs))

muFind <- function(lambs, draw = TRUE, path = "complexPlot") {
    # Проверка типа переменной - вектора
    stopifnot((is.numeric(lambs) == TRUE) || (is.complex(lambs) == TRUE))
    stopifnot(is.logical(draw))
    n <- length(lambs) # Переменная количества считанных данных
    if (n <= 1) stop("dim of linear operator is emprty or not valid")
    # section of init z1, z2, z3 ----------------------------------------------
    flash <- FALSE # Логическая переменная для обозначения разрешения вывода результата
    # Конструкция if (...) {..} else if (...) {..} для проверки сложности алгоритма 
    if (n == 2) {
        # Если точки всего две, то определеяем центр окружности и радуис по первым формулам
        mu <- .mu2(lambs[1], lambs[2]) # Посчитали мю как параметр алгоритма
        R <- .R2(lambs[1], lambs[2]) # посчитали раидус окружности, описанной около двух точек
        if (is.infinite(R)) flash <- FALSE # Если радиус посчитается бесконечным, то это окружность с началом координат
        else if (abs(R) < abs(mu)) flash <- TRUE # Если радиус меньше размера вектора до центра, то посчитано верно, и 
        # проставляем флаг на разрешение вывода результата
        
    } else if (n > 2) {
        
        # Если точек три и больше, то считаем как по первым так и по вторым формулам
        U <- combn(1:n, m = 2) # combn() составляет из вектора все возможные комбинации по m членов
        mu <- .mu2(lambs[U[1,]], lambs[U[2,]])
        R <- .R2(lambs[U[1,]], lambs[U[2,]])
        Flags <- .FlagsCalc(mu, R, n, lambs)
        
        if (isTRUE(all.equal(rep(F, length(R)), as.logical(apply(Flags, 1, prod))))) {
            U <- combn(1:n, m = 3)
            mu <- .mu3(lambs[U[1,]], lambs[U[2,]], lambs[U[3,]])
            R <- .R3(mu, lambs[U[1,]])
            Flags <- .FlagsCalc(mu, R, n, lambs)
            if (!isTRUE(all.equal(rep(F, length(R)), as.logical(apply(Flags, 1, prod))))) flash <- TRUE
            
        } else flash <- TRUE
        
        if (any(is.infinite(R))) flash <- FALSE
    }
    
    # Условный оператор для вывода результата
    if (flash & n == 2) {
        if (draw) .drawComplex(lambs = lambs, mu = mu, R = R, fileName = path)
        cat("\n")
        res <- c(mu, R, sqrt(R^2 / abs(mu)^2))
        names(res) = c("mu", "R", "r(A)")
        return(res)
    } else if (flash & n > 2) {
        Number_of_Radius <- which(1 == apply(Flags, 1, prod))
        results <- matrix(nrow = length(Number_of_Radius), ncol = 3)
        for (i in 1:length(Number_of_Radius)) {
            j <- Number_of_Radius[i]
            results[i,] <- c(mu[j], R[j], sqrt(R[j]^2 / abs(mu[j])^2))
        }
        if (draw) .drawComplex(lambs = lambs, mu = results[which.min(results[,2]), 1], R = results[which.min(results[,2]), 2], fileName = path)
        cat("\n")
        res <- results[which.min(results[,2]), ]
        names(res) = c("mu", "R", "r(A)")
        return(res)
        
    } else {
        stop("Center of Decart plot belongs to area of linear operator's spectre")
        
    }
}

lambs <- rnorm(35, 35, 5) + 1i * rnorm(35, 35, 5)
muFind(lambs)

library()