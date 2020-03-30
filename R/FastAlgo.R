GMSI <- function(A = matrix(), f = numeric(), u = numeric(), eps = 10e-4) {
    # Определение размерности
    n <- ncol(A)
    m <- nrow(A)
    # Проверка на размерность матрицы, она должна быть больше или равна 2
    assertive.base::assert_all_are_true(n >= 2, severity = "stop")
    # Проверка на невырожденность
    assertive.base::assert_all_are_false(x = (m < n), severity = "stop")
    # Проверка невырожденности матрицы А: det(A) != 0
    # В общем случае, определитель надо расписать, но если программа позволяет, можно оставить
    #assertive.base::assert_all_are_true(x = (det(A) != 0), severity = "stop")
    # Поиск лямбд
    lambs <- eigen(A)$values
    # Определение итерационного члена
    n <- length(lambs)
    flash <- FALSE # Логическая переменная для обозначения разрешения записи в текстовый файл
    
    
    
    # Функция рассчитывает центр окружности минимального радиуса
    # на комплексной плоскости через две точки, которые принадлежат
    # множеству собственных чисел линейного оператора
    .mucalc2 <- function(z1, z2) {
        return(((z1 + z2) / 2) + 1i * ((Im(z1 * Conj(z2)) * (z2 - z1))/(2 * (abs(z1 * Conj(z2)) + Re(z1 * Conj(z2))))))
    }
    # Функция производит рассчет радиуса окружности на основе 
    # значений точек на комплексной плоскости, которые принадлежат
    # множеству собственных чисел линейного оператора
    .Rcalc2 <- function(z1, z2) {
        return(sqrt((abs(z1 - z2)^2 * abs(Conj(z1) * z2)) / (2 * (abs(Conj(z1) * z2) + Re(Conj(z1) * z2)))))
    }
    # Функция рассчитывает центр окружности минимального радиуса
    # на комплексной плоскости через три точки, которые принадлежат
    # множеству собственных чисел линейного оператора
    .mucalc3 <- function(z1, z2, z3) {
        return(1i * (abs(z1)^2 * (z2 - z3) + abs(z2)^2 * (z3 - z1) + abs(z3)^2 * (z1 - z2))/(2 * Im(z1 * Conj(z2) + z2 * Conj(z3) + z3 * Conj(z1))))
    }
    # Функция производит рассчет радиуса окружности на основе 
    # значений точек на комплексной плоскости, которые принадлежат
    # множеству собственных чисел линейного оператора
    .Rcalc3 <- function(MU, z) {
        return(sqrt(abs(MU - z)^2))
    }
    # Данная функция необходима для рассчета матрицы логических значений
    # Матрица на выходе отвечает за выявление такого мю и радиуса к ней,
    # который помещал бы все точки, принадлежащие множеству собственных значений матрицы оператора
    # Соответственно радиусам и центрам соотвествуют строки, а числам столбцы
    # делаем произведение матрицы по строкам для того, чтобы выявить истинный во всех реализациях вариант
    .FlagsCalc <- function(mu, R, n, lambs) {
        Flags <- matrix(FALSE, nrow = length(R), ncol = n + 1)
        for (j in 1:length(R)) {
            for (i in 1:n) {
                Flags[j, i] <- ((abs(R[j]) > abs(mu[j] - lambs[i])) || (abs(R[j]) - abs(mu[j] - lambs[i]) == 0))
            }
            Flags[j, n + 1] <- abs(R[j]) < abs(mu[j])
        }
        return(Flags)
    }
    
    
    
    # Конструкция if (...) {..} else if (...) {..} для проверки сложности алгоритма 
    if (n == 2) {
        
        # Если точки всего две, то определеяем центр окружности и радуис по первым формулам
        Z1 <- lambs[1] # Считали как первую точку в алгоритм
        Z2 <- lambs[2] # Считали как вторую точку в алгоритм
        mu <- .mucalc2(Z1, Z2) # Посчитали мю как параметр алгоритма
        R <- .Rcalc2(Z1, Z2) # посчитали раидус окружности, описанной около двух точек
        if (is.infinite(R)) flash <- FALSE # Если радиус посчитается бесконечным, то это окружность с началом координат
        else if (abs(R) < abs(mu)) flash <- TRUE # Если радиус меньше размера вектора до центра, то посчитано верно, и 
        # проставляем флаг на разрешение вывода результата
        
    } else if (n > 2) {
        
        # Если точек три и больше, то считаем как по первым так и по вторым формулам
        U <- combn(1:n, m = 2) # comdn() составляет из вектора все возможные комбинации по m членов
        Z1 <- lambs[U[1,]] # 
        Z2 <- lambs[U[2,]]
        mu <- .mucalc2(Z1, Z2)
        R <- .Rcalc2(Z1, Z2)
        Flags <- .FlagsCalc(mu, R, n, lambs)
        
        if (isTRUE(all.equal(rep(F, length(R)), as.logical(apply(Flags, 1, prod))))) {
            
            U <- combn(1:n, m = 3)
            Z1 <- lambs[U[1,]]
            Z2 <- lambs[U[2,]]
            Z3 <- lambs[U[3,]]
            mu <- .mucalc3(Z1, Z2, Z3)
            R <- .Rcalc3(mu, Z1)
            Flags <- .FlagsCalc(mu, R, n, lambs)
            if (!isTRUE(all.equal(rep(F, length(R)), as.logical(apply(Flags, 1, prod))))) flash <- TRUE
            
        } else flash <- TRUE
        
        if (any(is.infinite(R))) flash <- FALSE
    }
    # write results into txt --------------------------------------------------
    # Условный оператор для вывода результата в файл
    if (flash & n == 2) {
        results <- c(mu, R, sqrt(R^2 / abs(mu)^2))
        names(results) <- c("Center:", "Radius:", "Ro:")
        print(results)
        cat("\n")
        
    } else if (flash & n > 2) {
        Number_of_Radius <- which(1 == apply(Flags, 1, prod))
        results <- matrix(nrow = length(Number_of_Radius), ncol = 3)
        for (i in 1:length(Number_of_Radius)) {
            j <- Number_of_Radius[i]
            results[i,] <- c(mu[j], R[j], sqrt(R[j]^2 / abs(mu[j])^2))
        }
        Number_of_Radius <- which.min(results[,2])
        results <- results[Number_of_Radius, ]
        names(results) <- c("Center:", "Radius:", "Ro:")
        print(results)
        cat("\n")
        
    } else {
        results <- "Center of Decart plot belongs to area of linear operator's spectre"
        print(results)
        assertive.base::assert_all_are_false(is.character(results), severity = "stop")
    }
    
    repeat{
        ut = u
        u <- u - (1/results[1]) * (A %*% u - f)
        if (max(abs(u - ut)) < eps) break
    }
    return(u)
    
}



GMSI.history <- function(A = matrix(), f = numeric(), u = numeric(), eps = 10e-4) {
   
     # Определение размерности
    n <- ncol(A)
    m <- nrow(A)
    # Проверка на размерность матрицы, она должна быть больше или равна 2
    assertive.base::assert_all_are_true(n >= 2, severity = "stop")
    # Проверка на невырожденность
    assertive.base::assert_all_are_false(x = (m < n), severity = "stop")
    # Проверка невырожденности матрицы А: det(A) != 0
    # В общем случае, определитель надо расписать, но если программа позволяет, можно оставить
    #assertive.base::assert_all_are_true(x = (det(A) != 0), severity = "stop")
    # Поиск лямбд
    lambs <- eigen(A)$values
    # Определение итерационного члена
    n <- length(lambs)
    
    
    flash <- FALSE # Логическая переменная для обозначения разрешения записи в текстовый файл
    # Функция рассчитывает центр окружности минимального радиуса
    # на комплексной плоскости через две точки, которые принадлежат
    # множеству собственных чисел линейного оператора
    .mucalc2 <- function(z1, z2) {
        return(((z1 + z2) / 2) + 1i * ((Im(z1 * Conj(z2)) * (z2 - z1))/(2 * (abs(z1 * Conj(z2)) + Re(z1 * Conj(z2))))))
    }
    # Функция производит рассчет радиуса окружности на основе 
    # значений точек на комплексной плоскости, которые принадлежат
    # множеству собственных чисел линейного оператора
    .Rcalc2 <- function(z1, z2) {
        return(sqrt((abs(z1 - z2)^2 * abs(Conj(z1) * z2)) / (2 * (abs(Conj(z1) * z2) + Re(Conj(z1) * z2)))))
    }
    # Функция рассчитывает центр окружности минимального радиуса
    # на комплексной плоскости через три точки, которые принадлежат
    # множеству собственных чисел линейного оператора
    .mucalc3 <- function(z1, z2, z3) {
        return(1i * (abs(z1)^2 * (z2 - z3) + abs(z2)^2 * (z3 - z1) + abs(z3)^2 * (z1 - z2))/(2 * Im(z1 * Conj(z2) + z2 * Conj(z3) + z3 * Conj(z1))))
    }
    # Функция производит рассчет радиуса окружности на основе 
    # значений точек на комплексной плоскости, которые принадлежат
    # множеству собственных чисел линейного оператора
    .Rcalc3 <- function(MU, z) {
        return(sqrt(abs(MU - z)^2))
    }
    # Данная функция необходима для рассчета матрицы логических значений
    # Матрица на выходе отвечает за выявление такого мю и радиуса к ней,
    # который помещал бы все точки, принадлежащие множеству собственных значений матрицы оператора
    # Соответственно радиусам и центрам соотвествуют строки, а числам столбцы
    # делаем произведение матрицы по строкам для того, чтобы выявить истинный во всех реализациях вариант
    .FlagsCalc <- function(mu, R, n, lambs) {
        Flags <- matrix(FALSE, nrow = length(R), ncol = n + 1)
        for (j in 1:length(R)) {
            for (i in 1:n) {
                Flags[j, i] <- ((abs(R[j]) > abs(mu[j] - lambs[i])) || (abs(R[j]) - abs(mu[j] - lambs[i]) == 0))
            }
            Flags[j, n + 1] <- abs(R[j]) < abs(mu[j])
        }
        return(Flags)
    }
    
    
    
    # Конструкция if (...) {..} else if (...) {..} для проверки сложности алгоритма 
    if (n == 2) {
        
        # Если точки всего две, то определеяем центр окружности и радуис по первым формулам
        Z1 <- lambs[1] # Считали как первую точку в алгоритм
        Z2 <- lambs[2] # Считали как вторую точку в алгоритм
        mu <- .mucalc2(Z1, Z2) # Посчитали мю как параметр алгоритма
        R <- .Rcalc2(Z1, Z2) # посчитали раидус окружности, описанной около двух точек
        if (is.infinite(R)) flash <- FALSE # Если радиус посчитается бесконечным, то это окружность с началом координат
        else if (abs(R) < abs(mu)) flash <- TRUE # Если радиус меньше размера вектора до центра, то посчитано верно, и 
        # проставляем флаг на разрешение вывода результата
        
    } else if (n > 2) {
        
        # Если точек три и больше, то считаем как по первым так и по вторым формулам
        U <- combn(1:n, m = 2) # comdn() составляет из вектора все возможные комбинации по m членов
        Z1 <- lambs[U[1,]] # 
        Z2 <- lambs[U[2,]]
        mu <- .mucalc2(Z1, Z2)
        R <- .Rcalc2(Z1, Z2)
        Flags <- .FlagsCalc(mu, R, n, lambs)
        
        if (isTRUE(all.equal(rep(F, length(R)), as.logical(apply(Flags, 1, prod))))) {
            
            U <- combn(1:n, m = 3)
            Z1 <- lambs[U[1,]]
            Z2 <- lambs[U[2,]]
            Z3 <- lambs[U[3,]]
            mu <- .mucalc3(Z1, Z2, Z3)
            R <- .Rcalc3(mu, Z1)
            Flags <- .FlagsCalc(mu, R, n, lambs)
            if (!isTRUE(all.equal(rep(F, length(R)), as.logical(apply(Flags, 1, prod))))) flash <- TRUE
            
        } else flash <- TRUE
        
        if (any(is.infinite(R))) flash <- FALSE
    }
    
    # write results into txt --------------------------------------------------
    # Условный оператор для вывода результата в файл
   
     if (flash & n == 2) {
        results <- c(mu, R, sqrt(R^2 / abs(mu)^2))
        names(results) <- c("Center:", "Radius:", "Ro:")
        print(results)
        cat("\n")
        
    } else if (flash & n > 2) {
        Number_of_Radius <- which(1 == apply(Flags, 1, prod))
        results <- matrix(nrow = length(Number_of_Radius), ncol = 3)
        for (i in 1:length(Number_of_Radius)) {
            j <- Number_of_Radius[i]
            results[i,] <- c(mu[j], R[j], sqrt(R[j]^2 / abs(mu[j])^2))
        }
        Number_of_Radius <- which.min(results[,2])
        results <- results[Number_of_Radius, ]
        names(results) <- c("Center:", "Radius:", "Ro:")
        print(results)
        cat("\n")
        
    } else {
        results <- "Center of Decart plot belongs to area of linear operator's spectre"
        print(results)
        assertive.base::assert_all_are_false(is.character(results), severity = "stop")
    }
    

# graphics ----------------------------------------------------------------
    
    if ("plotrix" %in% rownames(installed.packages()) == FALSE) {install.packages("plotrix")}
    library(plotrix)
    jpeg("outputs/complexPlot.jpg", width = 800, height = 450)
    plot(x = c(Re(lambs), Re(results["Center:"])), asp = 1,
         y = c(Im(lambs), Im(results["Center:"])),
         xlim = c(Re(results["Center:"]) - Re(results["Radius:"]) - 1, Re(results["Center:"]) + Re(results["Radius:"]) + 1),
         ylim = c(Im(results["Center:"]) - Re(results["Radius:"]) - 1, Im(results["Center:"]) + Re(results["Radius:"]) + 1),
         main = "Спектр на комплексной плоскости линейного оператора")
    lines(x = Re(lambs), 
          y = Im(lambs), 
          type = "l")
    lines(x = c(Re(lambs[1]), Re(lambs[n])),
          y = c(Im(lambs[1]), Im(lambs[n])))
    abline(h = 0, lty = 1)
    abline(v = 0, lty = 1)
    # Функция из библиотеки 'plotrix' для рисования окружности с центром и радиусом
    draw.circle(x = Re(results["Center:"]), y = Im(results["Center:"]), radius = Re(results["Radius:"]), border = "red")
    grid()
    dev.off()
    

# Simple iterations -------------------------------------------------------

    
    
    i <- 1
    u.hist <- matrix(nrow = length(u))
    repeat{
        u.hist <- cbind(u.hist, u)
        ut <- u
        u <- u - (1/results[1]) * (A %*% u - f)
        i <- i + 1
        if (max(abs(u - ut)) < eps) break
    }
    result.list <- list(iterations = i, var = u, var.hist = u.hist, results = results)
    return(result.list)
}

