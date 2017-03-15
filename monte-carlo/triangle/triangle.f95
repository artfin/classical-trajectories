program triangle_area

	implicit none
	real :: a, b, c, s, p, Area
	
	a = 3.1
	b = 4.2
	c = 5.5
	
	p = (a+b+c)/2

	Area = (p * (p - a) * (p - b) * (p - c))**0.5
	
	print *, 'Area: ', Area
	
end program triangle_area
