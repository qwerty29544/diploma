# Help func for circle ----------------------------------------------------

#' Title
#'
#' @return
#' @export
#'
#' @examples
getYmult <- function() {
    if (grDevices::dev.cur() == 1) {
        base::warning("No graphics device open.")
        ymult <- 1
    }
    else {
        xyasp <- graphics::par("pin")
        xycr <- base::diff(graphics::par("usr"))[c(1, 3)]
        ymult <- xyasp[1]/xyasp[2] * xycr[2]/xycr[1]
    }
    return(ymult)
}

# draw a circle func ------------------------------------------------------

#
#' Title
#'
#' @param x 
#' @param y 
#' @param radius 
#' @param nv 
#' @param border 
#' @param col 
#' @param lty 
#' @param density 
#' @param angle 
#' @param lwd 
#'
#' @return
#' @export
#'
#' @examples
circle <- function(x, y, radius, nv = 100, border = NULL, col = NA, lty = 1, density = NULL, angle = 45, lwd = 0.5) {
    # xylim <- graphics::par("usr")
    # plotdim <- graphics::par("pin")
    ymult <- getYmult()
    angle.inc <- 2 * pi/nv
    angles <- base::seq(0, 2 * pi - angle.inc, by = angle.inc)
    if (base::length(col) < base::length(radius)) 
        col <- base::rep(col, length.out = base::length(radius))
    for (circle in 1:base::length(radius)) {
        xv <- cos(angles) * radius[circle] + x
        yv <- sin(angles) * radius[circle] * ymult + y
        polygon(xv, yv, border = border, col = col[circle], lty = lty, 
                density = density, angle = angle, lwd = lwd)
    }
    invisible(list(x = xv, y = yv))
}

# mu with 2 cmplx ---------------------------------------------------------

# mu от 2 точек
#' Title
#'
#' @param z1 
#' @param z2 
#'
#' @return
#' @export
#'
#' @examples
mu2 <- function(z1, z2) {
    ((z1 + z2) / 2) + ((1i * Im(z1 * Conj(z2)) * (z2 - z1)) / (2 * (Mod(z1 * Conj(z2)) + Re(z1 * Conj(z2)))))
}

# R with 2 cmplx ----------------------------------------------------------

#
#' Title
#'
#' @param z1 
#' @param z2 
#'
#' @return
#' @export
#'
#' @examples
R2 <- function(z1, z2) {
    sqrt((Mod(z1 - z2) * Mod(z1 - z2) * Mod(Conj(z1) * z2)) / (2 * (Mod(Conj(z1) * z2) + Re(Conj(z1) * z2))))
}

# mu with 3 cmplx ---------------------------------------------------------

#
#' Title
#'
#' @param z1 
#' @param z2 
#' @param z3 
#'
#' @return
#' @export
#'
#' @examples
mu3 <- function(z1, z2, z3) {
    
}

# R with 3 cmplx  ---------------------------------------------------------

#
#' Title
#'
#' @param z1 
#' @param mu 
#'
#' @return
#' @export
#'
#' @examples
R3 <- function(z1, mu) {
    sqrt(Mod(z1 - mu)^2)
} 

# mu find -----------------------------------------------------------------

#' Title
#'
#' @param lambs 
#'
#' @return
#' @export
#'
#' @examples
muFind <- function(lambs = complex()) {
    stopifnot(!is.character(lambs))
    if (is.null(labms)) stop("lambs is NULL")
    if (length(lambs) == 1) stop("dim of operator is 1")
    dimA <- length(lambs)
    repeat{
        if (dimA == 2) {
            mu <- mu2(lambs[1], lambs[2])   # считаем значение центра для спектра
            R <- R2(lambs[1], lambs[2])     # считаем радиус для центра до значений точек спектра
            # Необходимо осуществить проверку условия на отрезок, принадлежащий началу координат
            if (is.infinite(R)) { # Если так, то выходим из алгоритма
                stop("Radius is infinite")
            }
            break        # выход из блока кода
            
        } else if (dimA >= 3) {
            
            U <- combn(1:dimA, m = 2)  # combn() составляет из вектора все возможные комбинации по m членов
            mu <- mu2(lambs[U[1,]], lambs[U[2,]])   # векторно считаем центры по двум точкам
            R <- R2(lambs[U[1,]], lambs[U[2,]])     # векторно считаем радиусы по двум точкам
            # Необходимо осуществить проверку условия на отрезок, принадлежащий началу координат
            if (is.infinite(R)) { # Если так, то выходим из алгоритма
                stop("Radius is infinite")
            }
            # В данном месте алгоритма необходимо перебрать все возможные варианты точек и центров,
            # в каждом варианте из которых мы получаем минимально покрывающие радиусы
            
            RFlash <- R[which.max(R)] # Я просто тут пытаюсь поэкспериментировать с перебором и не перебирать
            muFlash <- mu[which.max(R)]
            if ((prod(RFlash >= abs(muFlash - lambs)) == 1) && (RFlash < abs(muFlash))) {
                R <- RFlash
                mu <- muFlash
                break   # выход из блока кода
            } else {
                # переход на операции с 3 точками
                U <- combn(1:dimA, m = 3) # combn() составляет из вектора все возможные комбинации по m членов
                mu <- mu3(lambs[U[1,]], lambs[U[2,]], lambs[U[3,]])
                R <- R3(mu, lambs[U[1,]])
                if (is.infinite(R)) { # Если так, то выходим из алгоритма
                    stop("Radius is infinite")
                }
                
                RFlash <- R[which.max(R)] # Я просто тут пытаюсь поэкспериментировать с перебором и не перебирать
                muFlash <- mu[which.max(R)]
                if ((prod(RFlash >= abs(muFlash - lambs)) == 1) && (RFlash < abs(muFlash))) {
                    R <- RFlash
                    mu <- muFlash
                    break   # выход из блока кода
                } else {
                    stop("Спектр лежит на начале координат")
                }
            } 
        }
    }
    # Отрисовка ---------------------------------------------------------------
    
    grDevices::jpeg("complexPlot.jpg", width = 800, height = 800)
    plot(x = c(Re(lambs), Re(mu)), asp = 1,
         y = c(Im(lambs), Im(mu)),
         xlim = c(Re(mu) - 1.5*Re(R), Re(mu) + 1.5*Re(R)),
         ylim = c(Im(mu) - 1.5*Re(mu), Im(mu) + 1.5*Re(mu)),
         main = "Спектр на комплексной плоскости линейного оператора")
    lines(x = Re(lambs), 
          y = Im(lambs), 
          type = "l")
    lines(x = c(Re(lambs[1]), Re(lambs[dimA])),
          y = c(Im(lambs[1]), Im(lambs[dimA])))
    abline(h = 0, lty = 1)
    abline(v = 0, lty = 1)
    # Функция из библиотеки 'plotrix' для рисования окружности с центром и радиусом
    circle(x = Re(mu), y = Im(mu), radius = Re(R), border = "red")
    grid()
    grDevices::dev.off()

# возврат -----------------------------------------------------------------

    
    return(mu)
}


# muFind.File -------------------------------------------------------------

muFind.File <- function(path.input = "input.txt", path.output = "output.txt") {
    stopifnot(is.character(path.input), is.character(path.output))
    # Считывание данных о точках многоугольника из текстового документа
    lambs <- read.table(file = path.input, header = F)[[1]]
    if (is.null(labms)) stop("lambs is NULL")
    if (length(lambs) == 1) stop("dim of operator is 1")
    cat(lambs, "\n") # Вывод считанных данных в консоль как есть
    dimA <- length(lambs) # Переменная количества считанных данных

    repeat{
        if (dimA == 2) {
            mu <- mu2(lambs[1], lambs[2])   # считаем значение центра для спектра
            R <- R2(lambs[1], lambs[2])     # считаем радиус для центра до значений точек спектра
            # Необходимо осуществить проверку условия на отрезок, принадлежащий началу координат
            if (is.infinite(R)) { # Если так, то выходим из алгоритма
                write.table(x = "Radius is infinite", file = path.output, quote = FALSE)
                stop("Radius is infinite")
            }
            break        # выход из блока кода
            
        } else if (dimA >= 3) {
            
            U <- combn(1:dimA, m = 2)  # combn() составляет из вектора все возможные комбинации по m членов
            mu <- mu2(lambs[U[1,]], lambs[U[2,]])   # векторно считаем центры по двум точкам
            R <- R2(lambs[U[1,]], lambs[U[2,]])     # векторно считаем радиусы по двум точкам
            # Необходимо осуществить проверку условия на отрезок, принадлежащий началу координат
            if (is.infinite(R)) { # Если так, то выходим из алгоритма
                write.table(x = "Radius is infinite", file = path.output, quote = FALSE)
                stop("Radius is infinite")
            }
            # В данном месте алгоритма необходимо перебрать все возможные варианты точек и центров,
            # в каждом варианте из которых мы получаем минимально покрывающие радиусы
            
            RFlash <- R[which.max(R)] # Я просто тут пытаюсь поэкспериментировать с перебором и не перебирать
            muFlash <- mu[which.max(R)]
            if ((prod(RFlash >= abs(muFlash - lambs)) == 1) && (RFlash < abs(muFlash))) {
                R <- RFlash
                mu <- muFlash
                break   # выход из блока кода
            } else {
                # переход на операции с 3 точками
                U <- combn(1:dimA, m = 3) # combn() составляет из вектора все возможные комбинации по m членов
                mu <- mu3(lambs[U[1,]], lambs[U[2,]], lambs[U[3,]])
                R <- R3(mu, lambs[U[1,]])
                if (is.infinite(R)) { # Если так, то выходим из алгоритма
                    write.table(x = "Radius is infinite", file = path.output, quote = FALSE)
                    stop("Radius is infinite")
                }
                
                RFlash <- R[which.max(R)] # Я просто тут пытаюсь поэкспериментировать с перебором и не перебирать
                muFlash <- mu[which.max(R)]
                if ((prod(RFlash >= abs(muFlash - lambs)) == 1) && (RFlash < abs(muFlash))) {
                    R <- RFlash
                    mu <- muFlash
                    break   # выход из блока кода
                } else {
                    write.table(x = "Спектр лежит на начале координат", file = path.output, quote = FALSE)
                    stop("Спектр лежит на начале координат")
                }
            } 
        }
    }
    
    # Отрисовка ---------------------------------------------------------------
    
    grDevices::jpeg("complexPlot.jpg", width = 800, height = 800)
    plot(x = c(Re(lambs), Re(mu)), asp = 1,
         y = c(Im(lambs), Im(mu)),
         xlim = c(Re(mu) - 1.5*Re(R), Re(mu) + 1.5*Re(R)),
         ylim = c(Im(mu) - 1.5*Re(mu), Im(mu) + 1.5*Re(mu)),
         main = "Спектр на комплексной плоскости линейного оператора")
    lines(x = Re(lambs), 
          y = Im(lambs), 
          type = "l")
    lines(x = c(Re(lambs[1]), Re(lambs[dimA])),
          y = c(Im(lambs[1]), Im(lambs[dimA])))
    abline(h = 0, lty = 1)
    abline(v = 0, lty = 1)
    # Функция из библиотеки 'plotrix' для рисования окружности с центром и радиусом
    circle(x = Re(mu), y = Im(mu), radius = Re(R), border = "red")
    grid()
    grDevices::dev.off()

# возврат -----------------------------------------------------------------
    results <- c(mu, R, sqrt(R^2 / abs(mu)^2))
    names(results) <- c("Center:", "Radius:", "Ro:")
    print(results)
    cat("\n")
    write.table(x = results, file = path.output, quote = FALSE)
    return(mu)
}
# GMSI iterations ---------------------------------------------------------

#
#' Title
#'
#' @param A 
#' @param f 
#' @param u 
#' @param eps 
#'
#' @return
#' @export
#'
#' @examples
GMSI.muFind <- function(A = matrix(), f = numeric(), u = numeric(), eps = 10e-3) {
    stopifnot(is.matrix(A), is.numeric(A), is.numeric(f), is.numeric(u), is.numeric(eps))
# Блок приготовлений ------------------------------------------------------

    # Получение размерности матрицы А
    dimA <- dim(A)[1]
    
    # Проверка на n >= 2
    if (dimA[1] < 2) {
        warning("Operator must have dim >= 2")
        return(NULL)
    }
    
    # Проверка на размерность матрицы оператора
    if (dim(A)[1] != dim(A)[2]) {
        warning("Operator must be quadratic")
        return(NULL)
    }
    
    # Получаем собственные числа матрицы
    lambs <- eigen(A)$values

# Блок нахождения итерационного параметра ---------------------------------
    
    # Работаем далее по проверке условия для 2, 3, и более точек
    #'@TODO: МОЖНО ТУТ И ФУНКЦИЮ ОБЪЯВИТЬ ВООБЩЕТ
    repeat{
        if (dimA == 2) {
            mu <- mu2(lambs[1], lambs[2])   # считаем значение центра для спектра
            R <- R2(lambs[1], lambs[2])     # считаем радиус для центра до значений точек спектра
            # Необходимо осуществить проверку условия на отрезок, принадлежащий началу координат
            if (is.infinite(R)) { # Если так, то выходим из алгоритма
                warning("Radius is infinite")
                return(NULL)
            }
            break        # выход из блока кода
            
        } else if (dimA >= 3) {
            
            U <- combn(1:dimA, m = 2)  # combn() составляет из вектора все возможные комбинации по m членов
            mu <- mu2(lambs[U[1,]], lambs[U[2,]])   # векторно считаем центры по двум точкам
            R <- R2(lambs[U[1,]], lambs[U[2,]])     # векторно считаем радиусы по двум точкам
            # Необходимо осуществить проверку условия на отрезок, принадлежащий началу координат
            if (is.infinite(R)) { # Если так, то выходим из алгоритма
                warning("Radius is infinite")
                return(NULL)
            }
            # В данном месте алгоритма необходимо перебрать все возможные варианты точек и центров,
            # в каждом варианте из которых мы получаем минимально покрывающие радиусы
            
            RFlash <- R[which.max(R)] # Я просто тут пытаюсь поэкспериментировать с перебором и не перебирать
            muFlash <- mu[which.max(R)]
            if ((prod(RFlash >= abs(muFlash - lambs)) == 1) && (RFlash < abs(muFlash))) {
                R <- RFlash
                mu <- muFlash
                break   # выход из блока кода
            } else {
                # переход на операции с 3 точками
                U <- combn(1:dimA, m = 3) # combn() составляет из вектора все возможные комбинации по m членов
                mu <- mu3(lambs[U[1,]], lambs[U[2,]], lambs[U[3,]])
                R <- R3(mu, lambs[U[1,]])
                if (is.infinite(R)) { # Если так, то выходим из алгоритма
                    warning("Radius is infinite")
                    return(NULL)
                }
                
                RFlash <- R[which.max(R)] # Я просто тут пытаюсь поэкспериментировать с перебором и не перебирать
                muFlash <- mu[which.max(R)]
                if ((prod(RFlash >= abs(muFlash - lambs)) == 1) && (RFlash < abs(muFlash))) {
                    R <- RFlash
                    mu <- muFlash
                    break   # выход из блока кода
                } else {
                    warning("Спектр лежит на начале координат")
                    return(NULL)
                }
            } 
        }
    }
    # return(mu)
    
    repeat{
        ut = u
        u <- u - (1/mu) * (A %*% u - f)
        if (max(abs(u - ut)) < eps) break
    }
    return(u)
}

# GMSI iterations with history --------------------------------------------

#
#' Title
#'
#' @param A 
#' @param f 
#' @param u 
#' @param eps 
#'
#' @return
#' @export
#'
#' @examples
GMSI.muFind.history <- function(A = matrix(), f = numeric(), u = numeric(), eps = 10e-3) {
    stopifnot(is.matrix(A), is.numeric(A), is.numeric(f), is.numeric(u), is.numeric(eps))
    # Блок приготовлений ------------------------------------------------------
    
    # Получение размерности матрицы А
    dimA <- dim(A)[1]
    # Проверка на n >= 2
    if (dimA[1] < 2) {
        warning("Operator must have dim >= 2")
        return(NULL)
    }
    # Проверка на размерность матрицы оператора
    if (dim(A)[1] != dim(A)[2]) {
        warning("Operator must be quadratic")
        return(NULL)
    }
    # Получаем собственные числа матрицы
    lambs <- eigen(A)$values
    
    # Блок нахождения итерационного параметра ---------------------------------
    
    # Работаем далее по проверке условия для 2, 3, и более точек
    #'@TODO: МОЖНО ТУТ И ФУНКЦИЮ ОБЪЯВИТЬ ВООБЩЕТ
    repeat{
        if (dimA == 2) {
            mu <- mu2(lambs[1], lambs[2])   # считаем значение центра для спектра
            R <- R2(lambs[1], lambs[2])     # считаем радиус для центра до значений точек спектра
            # Необходимо осуществить проверку условия на отрезок, принадлежащий началу координат
            if (is.infinite(R)) { # Если так, то выходим из алгоритма
                warning("Radius is infinite")
                return(NULL)
            }
            break        # выход из блока кода
            
        } else if (dimA >= 3) {
            
            U <- combn(1:dimA, m = 2)  # combn() составляет из вектора все возможные комбинации по m членов
            mu <- mu2(lambs[U[1,]], lambs[U[2,]])   # векторно считаем центры по двум точкам
            R <- R2(lambs[U[1,]], lambs[U[2,]])     # векторно считаем радиусы по двум точкам
            # Необходимо осуществить проверку условия на отрезок, принадлежащий началу координат
            if (is.infinite(R)) { # Если так, то выходим из алгоритма
                warning("Radius is infinite")
                return(NULL)
            }
            # В данном месте алгоритма необходимо перебрать все возможные варианты точек и центров,
            # в каждом варианте из которых мы получаем минимально покрывающие радиусы
            
            RFlash <- R[which.max(R)] # Я просто тут пытаюсь поэкспериментировать с перебором и не перебирать
            muFlash <- mu[which.max(R)]
            if ((prod(RFlash >= abs(muFlash - lambs)) == 1) && (RFlash < abs(muFlash))) {
                R <- RFlash
                mu <- muFlash
                break   # выход из блока кода
            } else {
                # переход на операции с 3 точками
                U <- combn(1:dimA, m = 3) # combn() составляет из вектора все возможные комбинации по m членов
                mu <- mu3(lambs[U[1,]], lambs[U[2,]], lambs[U[3,]])
                R <- R3(mu, lambs[U[1,]])
                if (is.infinite(R)) { # Если так, то выходим из алгоритма
                    warning("Radius is infinite")
                    return(NULL)
                }
                
                RFlash <- R[which.max(R)] # Я просто тут пытаюсь поэкспериментировать с перебором и не перебирать
                muFlash <- mu[which.max(R)]
                if ((prod(RFlash >= abs(muFlash - lambs)) == 1) && (RFlash < abs(muFlash))) {
                    R <- RFlash
                    mu <- muFlash
                    break   # выход из блока кода
                } else {
                    warning("Спектр лежит на начале координат")
                    return(NULL)
                }
            } 
        }
    }
    
# Отрисовка ---------------------------------------------------------------

    grDevices::jpeg("complexPlot.jpg", width = 800, height = 800)
    plot(x = c(Re(lambs), Re(mu)), asp = 1,
         y = c(Im(lambs), Im(mu)),
         xlim = c(Re(mu) - Re(R) - 1, Re(mu) + Re(R) + 1),
         ylim = c(Im(mu) - Re(mu) - 1, Im(mu) + Re(mu) + 1),
         main = "Спектр на комплексной плоскости линейного оператора")
    lines(x = Re(lambs), 
          y = Im(lambs), 
          type = "l")
    lines(x = c(Re(lambs[1]), Re(lambs[dimA])),
          y = c(Im(lambs[1]), Im(lambs[dimA])))
    abline(h = 0, lty = 1)
    abline(v = 0, lty = 1)
    # Функция из библиотеки 'plotrix' для рисования окружности с центром и радиусом
    circle(x = Re(mu), y = Im(mu), radius = Re(R), border = "red")
    grid()
    grDevices::dev.off()

# итерации ----------------------------------------------------------------

    i <- 0
    u.hist <- matrix(u ,nrow = length(u))
    repeat{
        ut = u
        u <- u - (1/mu) * (A %*% u - f)
        u.hist <- cbind(u.hist, u)
        i <- i + 1
        if (max(abs(u - ut)) < eps) break
    }
    result.list <- list(iterations = i, var = u, var.hist = u.hist, results = c(Center = mu, Radius = R))
    return(result.list)
}


# GMSI Iterations ---------------------------------------------------------

#' Title
#'
#' @param A 
#' @param f 
#' @param u 
#' @param mu 
#' @param eps 
#'
#' @return
#' @export
#'
#' @examples
GMSI <- function(A = matrix(), f = numeric(), u = numeric(), mu = numeric(), eps = 10e-4) {
    stopifnot(is.matrix(A), is.numeric(A), is.numeric(f), is.numeric(u), is.numeric(eps))
    repeat{
        ut = u
        u <- u - (1/mu) * (A %*% u - f)
        if (max(abs(u - ut)) < eps) break
    }
    return(u)
}


# GMSI.history ------------------------------------------------------------

#' Title
#'
#' @param A 
#' @param f 
#' @param u 
#' @param mu 
#' @param eps 
#'
#' @return
#' @export
#'
#' @examples
GMSI.history <- function(A = matrix(), f = numeric(), u = numeric(), mu = numeric(), eps = 10e-4) {
    stopifnot(is.matrix(A), is.numeric(A), is.numeric(f), is.numeric(u), is.numeric(eps))
    i <- 0
    u.hist <- matrix(u, nrow = length(u))
    repeat{
        ut = u
        u <- u - (1/mu) * (A %*% u - f)
        i <- i + 1
        u.hist <- cbind(u.hist, u)
        if (max(abs(u - ut)) < eps) break
    }
    return.list <- list(iterations = i, var = u, var.hist = u.hist)
    return(return.list)
}

# GMSI.muFind.lambs -------------------------------------------------------

#' Title
#'
#' @param A 
#' @param f 
#' @param u 
#' @param mu 
#' @param eps 
#' @param lambs 
#'
#' @return
#' @export
#'
#' @examples
GMSI.muFind.lambs <- function(A = matrix(), f = numeric(), u = numeric(), eps = 10e-4, lambs = complex()) {
    stopifnot(is.matrix(A), is.numeric(A), is.numeric(f), is.numeric(u), is.numeric(eps), is.complex(lambs))
    # Получение размерности матрицы А
    dimA <- dim(A)[1]
    
    # Проверка на n >= 2
    if (dimA[1] < 2) {
        warning("Operator must have dim >= 2")
        return(NULL)
    }
    
    # Проверка на размерность матрицы оператора
    if (dim(A)[1] != dim(A)[2]) {
        warning("Operator must be quadratic")
        return(NULL)
    }
    # Блок нахождения итерационного параметра ---------------------------------
    
    # Работаем далее по проверке условия для 2, 3, и более точек
    #'@TODO: МОЖНО ТУТ И ФУНКЦИЮ ОБЪЯВИТЬ ВООБЩЕТ
    repeat{
        if (dimA == 2) {
            mu <- mu2(lambs[1], lambs[2])   # считаем значение центра для спектра
            R <- R2(lambs[1], lambs[2])     # считаем радиус для центра до значений точек спектра
            # Необходимо осуществить проверку условия на отрезок, принадлежащий началу координат
            if (is.infinite(R)) { # Если так, то выходим из алгоритма
                warning("Radius is infinite")
                return(NULL)
            }
            break        # выход из блока кода
            
        } else if (dimA >= 3) {
            
            U <- combn(1:dimA, m = 2)  # combn() составляет из вектора все возможные комбинации по m членов
            mu <- mu2(lambs[U[1,]], lambs[U[2,]])   # векторно считаем центры по двум точкам
            R <- R2(lambs[U[1,]], lambs[U[2,]])     # векторно считаем радиусы по двум точкам
            # Необходимо осуществить проверку условия на отрезок, принадлежащий началу координат
            if (is.infinite(R)) { # Если так, то выходим из алгоритма
                warning("Radius is infinite")
                return(NULL)
            }
            # В данном месте алгоритма необходимо перебрать все возможные варианты точек и центров,
            # в каждом варианте из которых мы получаем минимально покрывающие радиусы
            
            RFlash <- R[which.max(R)] # Я просто тут пытаюсь поэкспериментировать с перебором и не перебирать
            muFlash <- mu[which.max(R)]
            if ((prod(RFlash >= abs(muFlash - lambs)) == 1) && (RFlash < abs(muFlash))) {
                R <- RFlash
                mu <- muFlash
                break   # выход из блока кода
            } else {
                # переход на операции с 3 точками
                U <- combn(1:dimA, m = 3) # combn() составляет из вектора все возможные комбинации по m членов
                mu <- mu3(lambs[U[1,]], lambs[U[2,]], lambs[U[3,]])
                R <- R3(mu, lambs[U[1,]])
                if (is.infinite(R)) { # Если так, то выходим из алгоритма
                    warning("Radius is infinite")
                    return(NULL)
                }
                
                RFlash <- R[which.max(R)] # Я просто тут пытаюсь поэкспериментировать с перебором и не перебирать
                muFlash <- mu[which.max(R)]
                if ((prod(RFlash >= abs(muFlash - lambs)) == 1) && (RFlash < abs(muFlash))) {
                    R <- RFlash
                    mu <- muFlash
                    break   # выход из блока кода
                } else {
                    warning("Спектр лежит на начале координат")
                    return(NULL)
                }
            } 
        }
    }
    # return(mu)
    
    repeat{
        ut = u
        u <- u - (1/mu) * (A %*% u - f)
        if (max(abs(u - ut)) < eps) break
    }
    return(u)
}


# GMSI file ---------------------------------------------------------------

GMSI.File <- function(A = matrix(), f = numeric(), u = numeric(), 
                      eps = 10e-4, path.input = "input.txt", path.output = "output.txt") {
    stopifnot(is.character(path.input), is.character(path.output), 
              is.matrix(A), is.numeric(A), is.numeric(f), is.numeric(u), is.numeric(eps))
    # Считывание данных о точках многоугольника из текстового документа
    lambs <- read.table(file = path.input, header = F)[[1]]
    if (is.null(labms)) stop("lambs is NULL")
    if (length(lambs) == 1) stop("dim of operator is 1")
    cat(lambs, "\n") # Вывод считанных данных в консоль как есть
    dimA <- length(lambs) # Переменная количества считанных данных
    
    repeat{
        if (dimA == 2) {
            mu <- mu2(lambs[1], lambs[2])   # считаем значение центра для спектра
            R <- R2(lambs[1], lambs[2])     # считаем радиус для центра до значений точек спектра
            # Необходимо осуществить проверку условия на отрезок, принадлежащий началу координат
            if (is.infinite(R)) { # Если так, то выходим из алгоритма
                stop("Radius is infinite")
            }
            break        # выход из блока кода
            
        } else if (dimA >= 3) {
            
            U <- combn(1:dimA, m = 2)  # combn() составляет из вектора все возможные комбинации по m членов
            mu <- mu2(lambs[U[1,]], lambs[U[2,]])   # векторно считаем центры по двум точкам
            R <- R2(lambs[U[1,]], lambs[U[2,]])     # векторно считаем радиусы по двум точкам
            # Необходимо осуществить проверку условия на отрезок, принадлежащий началу координат
            if (is.infinite(R)) { # Если так, то выходим из алгоритма
                stop("Radius is infinite")
            }
            # В данном месте алгоритма необходимо перебрать все возможные варианты точек и центров,
            # в каждом варианте из которых мы получаем минимально покрывающие радиусы
            
            RFlash <- R[which.max(R)] # Я просто тут пытаюсь поэкспериментировать с перебором и не перебирать
            muFlash <- mu[which.max(R)]
            if ((prod(RFlash >= abs(muFlash - lambs)) == 1) && (RFlash < abs(muFlash))) {
                R <- RFlash
                mu <- muFlash
                break   # выход из блока кода
            } else {
                # переход на операции с 3 точками
                U <- combn(1:dimA, m = 3) # combn() составляет из вектора все возможные комбинации по m членов
                mu <- mu3(lambs[U[1,]], lambs[U[2,]], lambs[U[3,]])
                R <- R3(mu, lambs[U[1,]])
                if (is.infinite(R)) { # Если так, то выходим из алгоритма
                    stop("Radius is infinite")
                }
                
                RFlash <- R[which.max(R)] # Я просто тут пытаюсь поэкспериментировать с перебором и не перебирать
                muFlash <- mu[which.max(R)]
                if ((prod(RFlash >= abs(muFlash - lambs)) == 1) && (RFlash < abs(muFlash))) {
                    R <- RFlash
                    mu <- muFlash
                    break   # выход из блока кода
                } else {
                    stop("Спектр лежит на начале координат")
                }
            } 
        }
    }
    repeat{
        ut = u
        u <- u - (1/mu) * (A %*% u - f)
        if (max(abs(u - ut)) < eps) break
    }
    return(u)
    write.table(x = u, file = path.output, quote = FALSE)
}

# Chebishev's iterations --------------------------------------------------

chebishev <- function() {
    
}

# Gradient descent ---------------------------------------------------------------

#' Title
#'
#' @param A 
#' @param f 
#' @param u 
#' @param eps 
#'
#' @return
#' @export
#'
#' @examples
GDM <- function(A = matrix(), f = numeric(), u = numeric(), eps = 10e-4) {
    stopifnot(is.matrix(A), is.numeric(A), is.numeric(f), is.numeric(u), is.numeric(eps))
    dimA <- dim(A)[1]
    # Проверка на n >= 2
    if (dim(A)[1] < 2) {
        stop("Operator must have dim >= 2")
    }
    # Проверка на размерность матрицы оператора
    if (dim(A)[1] != dim(A)[2]) {
        stop("Operator must be quadratic")
    }
    r <- A %*% u - f
    chisl = (t(t(A) %*% r) %*% (t(A) %*% r))
    znam = t(A %*% t(A) %*% r) %*% (A %*% t(A) %*% r)
    u <- u - (chisl/znam)[1,1] * (t(A) %*% r)
    if ( abs((t(r) %*% (A %*% t(A) %*% r)))^2 / ((t(r) %*% r) * (t(A %*% t(A) %*% r)^2 %*% (A %*% t(A) %*% r)^2)) >= 1) {
        stop("Iterations do not sovling operator equation")
    }
    repeat{
        ut <- u
        chisl = (t(t(A) %*% r) %*% (t(A) %*% r))
        znam = t(A %*% t(A) %*% r) %*% (A %*% t(A) %*% r)
        r <- r - (chisl/znam)[1,1] * (A %*% t(A) %*% r)  
        u <- u - (chisl/znam)[1,1] * (t(A) %*% r)  
        if (max(abs(u - ut)) > eps) break
    }
    return(u)
}


# IMRES -------------------------------------------------------------------

#' Title
#'
#' @param A 
#' @param f 
#' @param u 
#' @param eps 
#'
#' @return
#' @export
#'
#' @examples
IMRES <- function(A = matrix(), f = numeric(), u = numeric(), eps = 10e-4) {
    stopifnot(is.matrix(A), is.numeric(A), is.numeric(f), is.numeric(u), is.numeric(eps))
    if (dim(A)[1] != dim(A)[2]) {
        stop("Operator must be a quadratic")
    }
    
    h <- A %*% u - f
    if ((1 - ((t(h) %*%  (A %*% h))^2 / ((t(h) %*% h) * (t(A %*% h) %*% (A %*% h))))) < 0) {
        stop("q >= 1, method is growing")
    }
    
    repeat {
        ut <- u
        h <- A %*% u - f
        tau <- (t(h) %*% (A %*% h)) / (t(A %*% h) %*% (A %*% h))
        u <- u - tau[1,1] * h
        if (max(abs(u - ut)) < eps) break
    }
    return(u)
}

# IMRES history -----------------------------------------------------------

#' Title
#'
#' @param A 
#' @param f 
#' @param u 
#' @param eps 
#'
#' @return
#' @export
#'
#' @examples
IMRES.history <- function(A = matrix(), f = numeric(), u = numeric(), eps = 10e-4) {
    stopifnot(is.matrix(A), is.numeric(A), is.numeric(f), is.numeric(u), is.numeric(eps))
    if (dim(A)[1] != dim(A)[2]) {
        stop("Operator must be a quadratic")
    }
    
    h <- A %*% u - f
    if ((1 - ((t(h) %*%  (A %*% h))^2 / ((t(h) %*% h) * (t(A %*% h) %*% (A %*% h))))) < 0) {
        stop("q >= 1, method is growing")
    }
    
    i = 0
    u.hist <- matrix(u, nrow = length(u))
    repeat {
        ut <- u
        h <- A %*% u - f
        tau <- (t(h) %*% (A %*% h)) / (t(A %*% h) %*% (A %*% h))
        u <- u - tau[1,1] * h
        i <- i + 1
        u.hist <- cbind(u.hist, u)
        if (max(abs(u - ut)) < eps) break
    }
    result.list <- list(iterations = i, var = u, var.hist = u.hist)
    return(result.list)
}

# MIMRES ------------------------------------------------------------------



# Matrix method -----------------------------------------------------------

#' Title
#'
#' @param A 
#' @param f 
#'
#' @return
#' @export
#'
#' @examples
MatrixMethod <- function(A = matirx(), f = numeric()) {
    stopifnot(is.matrix(A), is.numeric(A), is.numeric(f))
    return(solve(A) %*% f)
}


# Jacobi method -----------------------------------------------------------

#' Title
#'
#' @param A 
#' @param f 
#' @param u 
#' @param eps 
#'
#' @return
#' @export
#'
#' @examples
Jacobi <- function(A = matrix(), f = numeric(), u = numeric(), eps = 10e-4) {
    stopifnot(is.matrix(A), is.numeric(A), is.numeric(f), is.numeric(u), is.numeric(eps))
    dimA <- dim(A)[1]
    
    # Проверка на n >= 2
    if (dimA[1] < 2) {
        stop("Operator must have dim >= 2")

    }
    # Проверка на размерность матрицы оператора
    if (dim(A)[1] != dim(A)[2]) {
        stop("Operator must be quadratic")

    }
    
    # Проверка диагонального преобладания
    if (!all.equal(diag(A), apply(A, 1, max))) {
        stop("Operator must have a diagonal priority!")

    }

    
    repeat {
        TempX <- f
        for (i in 1:dimA) {
            for (g in 1:dimA) {
                if (i != g) {
                    TempX[i] <-  TempX[i] - A[i, g] * u[g]
                }
            }
        }
        TempX <- TempX / diag(A)
        norm <- max(abs(u - TempX))
        u <- TempX
        if (norm < eps) break
    }
    return(u)
}


# Jacobi history ----------------------------------------------------------

#' Title
#'
#' @param A 
#' @param f 
#' @param u 
#' @param eps 
#'
#' @return
#' @export
#'
#' @examples
Jacobi.history <- function(A = matrix(), f = numeric(), u = numeric(), eps = 10e-4) {
    stopifnot(is.matrix(A), is.numeric(A), is.numeric(f), is.numeric(u), is.numeric(eps))
    dimA <- dim(A)[1]
    
    # Проверка на n >= 2
    if (dimA[1] < 2) {
        stop("Operator must have dim >= 2")

    }
    # Проверка на размерность матрицы оператора
    if (dim(A)[1] != dim(A)[2]) {
        stop("Operator must be quadratic")

    }
    
    # Проверка диагонального преобладания
    if (!all.equal(diag(A), apply(A, 1, max))) {
        stop("Operator must have a diagonal priority!")

    }
    
    iter <- 0
    u.hist <- matrix(u ,nrow = length(u)) 
    repeat {
        TempX <- f
        for (i in 1:dimA) {
            for (g in 1:dimA) {
                if (i != g) {
                    TempX[i] <-  TempX[i] - A[i, g] * u[g]
                }
            }
        }
        TempX <- TempX / diag(A)
        norm <- max(abs(u - TempX))
        u <- TempX
        iter <- iter + 1
        u.hist <- cbind(u.hist, u)
        if (norm < eps) break

    }
    result.list <- list(iterations = iter, var = u, var.hist = u.hist)
    return(result.list)
}


# SIM ---------------------------------------------------------------------

#' Title
#'
#' @param A 
#' @param f 
#' @param u 
#' @param eps 
#'
#' @return
#' @export
#'
#' @examples
SIM <- function(A = matrix(), f = numeric(), u = numeric(), eps = 10e-4) {
    stopifnot(is.matrix(A), is.numeric(A), is.numeric(f), is.numeric(u), is.numeric(eps))
    dimA <- dim(A)[1]
    
    # Проверка на n >= 2
    if (dimA[1] < 2) {
        stop("Operator must have dim >= 2")
    }
    # Проверка на размерность матрицы оператора
    if (dim(A)[1] != dim(A)[2]) {
        stop("Operator must be quadratic")
    }

    B <- diag(1, nrow = dimA, ncol = dimA) - A
  
      repeat{
        ut <- u
        u <- B %*% u + f
        if (max(abs(u - ut)) < eps) break
    }
    return(u)
}


# SIM History -------------------------------------------------------------

SIM.history <- function(A = matrix(), f = numeric(), u = numeric(), eps = 10e-4) {
    stopifnot(is.matrix(A), is.numeric(A), is.numeric(f), is.numeric(u), is.numeric(eps))
    dimA <- dim(A)[1]
    
    # Проверка на n >= 2
    if (dimA[1] < 2) {
        stop("Operator must have dim >= 2")
    }
    # Проверка на размерность матрицы оператора
    if (dim(A)[1] != dim(A)[2]) {
        stop("Operator must be quadratic")
    }

    B <- diag(1, nrow = dimA, ncol = dimA) - A
    i <- 0
    u.hist <- matrix(u, nrow = length(u))
    repeat{
        ut <- u
        u <- B %*% u + f
        i <- i + 1
        u.hist <- cbind(u.hist, u)
        if (max(abs(u - ut)) < eps) break
    }
    result.list <- list(iterations = i, var = u, var.hist = u.hist)
    return(result.list)
}