program first
	integer, parameter :: n = 25
	integer :: i
	real :: error, h, y
	x = 0.5
	h = 1.0
	print *, "i, h, y, error"
	do i = 1, n
		h = 0.25 * h
		y = (sin(x + h) - sin(x)) / h
		error = abs(cos(x) - y)
		print *, i, h, y, error
	end do
end program first
