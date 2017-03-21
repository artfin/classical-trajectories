default: legendre

legendre: test.c liblegendre.a
	clang $^ -o $@ -lm

#legendre: legendre_polynomial_prb.c liblegendre.a
	#clang $^ -o $@ -lm

liblegendre.a: legendre_polynomial.o 
	ar -rcv $@ $^

 
