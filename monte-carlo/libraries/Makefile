default: wacky_dynamic 
# wacky_static

wacky_static: wacky.c libwacky.a 
	clang $^ -o $@

wacky_dynamic: wacky.c libwacky.so
	clang wacky.c -L. -lwacky -o $@

libwacky.a: wacky_math.o
	ar -rcv $@ $^  

libwacky.so: wacky_math.o
	clang -shared -o libwacky.so wacky_math.o
