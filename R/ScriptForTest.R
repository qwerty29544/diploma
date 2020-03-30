
# Section for text to complex ---------------------------------------------

# do not work
as.complex('5 + 1i*7')

# work
class(5 + 1i*7)

# work
class(5 + 7i)

# work
all.equal(5 + 1i*7, 5 + 7i)

# test
as.complex("5") + 1i*as.numeric("7")
'[:punct:]'


t <- read.table(file = "Docs/complexNumbeRS", header = T)
class(t[[1]])
t[[1]] <- as.complex(t[[1]])
t[[1]]
t

Conj(t[[1]])


library(Rcpp)
cppFunction('void Jacobi (int N, NumericMatrix A, NumericVector F, NumericVector X)
{
	const double eps = 0.001;
	double* TempX = new double[N];
	double norm; // норма, определяемая как наибольшая разность компонент столбца иксов соседних итераций.

	do {
		for (int i = 0; i < N; i++) {
			TempX[i] = F[i];
			for (int g = 0; g < N; g++) {
				if (i != g)
					TempX[i] -= A[i][g] * X[g];
			}
			TempX[i] /= A[i][i];
		}
        norm = fabs(X[0] - TempX[0]);
		for (int h = 0; h < N; h++) {
			if (fabs(X[h] - TempX[h]) > norm)
				norm = fabs(X[h] - TempX[h]);
			X[h] = TempX[h];
		}
	} while (norm > eps);
	delete[] TempX;
}')

