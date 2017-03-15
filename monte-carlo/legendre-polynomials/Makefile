default: legendre

legendre: test.c liblegendre.a
	clang $^ -o $@ 

#legendre: legendre_polynomial_prb.c liblegendre.a
#	clang $^ -o $@

liblegendre.a: legendre_polynomial.o
	ar -rcv $@ $^

 
