
# libraries init ----------------------------------------------------------


# Проверка наличия пакета plotrix/ Если нет, то загрузка пакета 
if ("plotrix" %in% rownames(installed.packages()) == FALSE) {install.packages("plotrix")}
library(plotrix)


# read file with lambs and results ----------------------------------------

# Считывание данных о точках многоугольника из текстового документа
lambs <- read.table(file = "inputs/complexNumbeRS", header = T)[[1]]
cat(lambs) # Вывод считанных данных в консоль как есть
n <- length(lambs) # Переменная количества считанных данных

# Считывание данных из текстового файла вывода
results1 <- read.table(file = "outputs/results.txt", header = T)
results <- results1[[1]] # Запись в переменную объект вектора результатов
names(results) <- rownames(results1) # Название для вектора результатов
print(results) # Вывод результатов алгоритма


# Отрисовка результатов ---------------------------------------------------

# Отрисовка результатов вычислений для спектра линейного оператора
jpeg("outputs/complexPlot.jpg")
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


# очистка памяти ----------------------------------------------------------

rm(list = ls())
