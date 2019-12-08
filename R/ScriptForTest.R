
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
