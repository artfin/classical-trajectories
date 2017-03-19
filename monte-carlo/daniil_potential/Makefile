default: test

test: test.c liblegendre.a
	clang $^ -o $@ -lm

liblegendre.a: legendre_polynomial.o
	ar -rcv $@ $^
