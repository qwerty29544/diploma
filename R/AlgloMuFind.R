library(plotrix)
# Functions init ----------------------------------------------------------


#' Title Данная функция возвращает результат вычисления центра окружности на комплексной плоскости для спектральной области линейного оператора
#' 
#' @param z1 is the first labmda of linear operator's Eigenvector 
#' @param z2 is the second labmda of linear operator's Eigenvector
#'
#' @return mu center 
#' @export
#'
#' @examples 
mucalc2 <- function(z1, z2) {
    return( ((z1 + z2) / 2) + 1i * ((Im(z1 * Conj(z2)) * (z2 - z1))/(2 * (abs(z1 * Conj(z2)) + Re(z1 * Conj(z2))))) )
}

#' Title считает радиус окружности по двум комплексным точкам спектра
#'
#' @param z1 
#' @param z2 
#'
#' @return R center
#' @export
#'
#' @examples
Rcalc2 <- function(z1, z2) {
    return(sqrt((abs(z1 - z2)^2 * abs(Conj(z1) * z2)) / (2 * (abs(Conj(z1) * z2) + Re(Conj(z1) * z2))) ))
}

#' Title считает центр окружности по трём комплексным точкам спектра
#'
#' @param z1 
#' @param z2 
#' @param z3 
#'
#' @return mu center
#' @export
#'
#' @examples
mucalc3 <- function(z1, z2, z3) {
    return(1i * (abs(z1)^2 * (z2 - z3) + abs(z2)^2 * (z3 - z1) + abs(z3)^2 * (z1 - z2))/(2 * Im(z1 * Conj(z2) + z2 * Conj(z3) + z3 * Conj(z1))))
}

#' Title считает радиус окружности по трём комплексным точкам спектра
#'
#' @param MU 
#' @param z 
#'
#' @return R
#' @export
#'
#' @examples
Rcalc3 <- function(MU, z) {
    return(sqrt(abs(MU - z)^2))
}

#' Title Данная функция реализует матрицу проверок для сочетаний центр многоугольника - точки многоугольника
#' Матрица проверок позволяет векторизовать процедуру проверки самого оптимального из предложенных центров окружности вокруг многоугольника. На основе генерации матрицы булевых значений, мы в строчках записываем мне сопоставления радиуса окружности с разностями межу центром окружности и крайними точками многоугольника.
#'
#' @param mu - центры окружностей, созданных треугольниками или отрезками внутри многоугольника
#' @param R - радиусы данных окружностей, которые создаются при помощи треугольков или отрезков многоугольника из mu
#' @param n - количество входных крайних точек спектральной области - создают многоугольник на комплексной плоскости
#' @param lambs - крайние точки спектра линейного оператора, создающие многоугольник на комплексной плоскости, для которого происходит поиск центра окружности и радиуса минимально возможного размера.
#'
#' @return Flags  - 
#' @export
#'
#' @examples
FlagsCalc <- function(mu, R, n, lambs) {
    Flags <- matrix(FALSE, nrow = length(R), ncol = n + 1)
    for (j in 1:length(R)) {
        for (i in 1:n) {
            Flags[j, i] <- ((abs(R[j]) > abs(mu[j] - lambs[i])) || (abs(R[j]) - abs(mu[j] - lambs[i]) == 0))
        }
        Flags[j, n + 1] <- abs(R[j]) < abs(mu[j])
    }
    return(Flags)
}


# section of initialization complex numbers -------------------------------

# Считывание данных о точках многоугольника из текстового документа
lambs <- read.table(file = "Docs/complexNumbeRS", header = T)[[1]]
cat(lambs) # Вывод считанных данных в консоль как есть
n <- length(lambs) # Переменная количества считанных данных

# section of init z1, z2, z3 ----------------------------------------------

flash <- FALSE # Логическая переменная для обозначения разрешения записи в текстовый файл

# Конструкция if (...) {..} else if (...) {..} для проверки сложности алгоритма 
if (n == 2) {
    
    # Если точки всего две, то определеяем центр окружности и радуис по первым формулам
    Z1 <- lambs[1] # Считали как первую точку в алгоритм
    Z2 <- lambs[2] # Считали как вторую точку в алгоритм
    mu <- mucalc2(Z1, Z2) # Посчитали мю как параметр алгоритма
    R <- Rcalc2(Z1, Z2) # посчитали раидус окружности, описанной около двух точек
    if (is.infinite(R)) flash <- FALSE # Если радиус посчитается бесконечным, то это окружность с началом координат
    else if (abs(R) < abs(mu)) flash <- TRUE # Если радиус меньше размера вектора до центра, то посчитано верно, и 
    # проставляем флаг на разрешение вывода результата
    
} else if (n > 2) {
    
    # Если точек три и больше, то считаем как по первым так и по вторым формулам
    U <- combn(1:n, m = 2) # comdn() составляет из вектора все возможные комбинации по m членов
    Z1 <- lambs[U[1,]] # 
    Z2 <- lambs[U[2,]]
    mu <- mucalc2(Z1, Z2)
    R <- Rcalc2(Z1, Z2)
    Flags <- FlagsCalc(mu, R, n, lambs)
    
    
    if (isTRUE(all.equal(rep(F, length(R)), as.logical(apply(Flags, 1, prod))))) {
        
        U <- combn(1:n, m = 3)
        Z1 <- lambs[U[1,]]
        Z2 <- lambs[U[2,]]
        Z3 <- lambs[U[3,]]
        mu <- mucalc3(Z1, Z2, Z3)
        R <- Rcalc3(mu, Z1)
        Flags <- FlagsCalc(mu, R, n, lambs)
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
    write.table(x = results, file = "Docs/results.txt", quote = FALSE)
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
    write.table(x = results, file = "Docs/results.txt", quote = FALSE)
} else {
    results <- "Center of Decart plot belongs to area of linear operator's spectre"
    write.table(x = results, file = "Docs/results.txt", quote = FALSE)
}


# Отрисовка результатов вычислений для спектра линейного оператора
jpeg("Docs/complexPlot.jpg")
plot(x = c(Re(lambs), Re(results["Center:"])), 
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

# Очистка памяти проекта
rm(list = ls())

