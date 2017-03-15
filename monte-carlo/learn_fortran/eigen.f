Program Eigenvalue
c finding the eigenvalues of a complex matrix using LAPACK
	Implicit none
c declarations, notice double precision
	complex*16 A(3,3), b(3), DUMMY(1,1), WORK(6)
	integer i, ok
c define matrix A
	A(1,1)=(3.1, -1.8)
	A(1,2)=(1.3, 0.2) 
	A(1,3)=(-5.7, -4.3)
	A(2,1)=(1.0, 0)
	A(2,2)=(-6.9, 3.2)
	A(2,3)=(5.8, 2.2)
	A(3,1)=(3.4, -4.0)
	A(3,2)=(7.2, 2.9)
	A(3,3)=(-8.8, 3.2)
c
c find the solution using the LAPACK routine ZGEEV
	call ZGEEV('N', 'N', 3, A, 3, b, DUMMY, 1, DUMMY, 1, WORK, 6, WORK, ok)
c
c parameters in the order as they appear in the function call
c    no left eigenvectors, no right eigenvectors, order of input matrix A,
c    input matrix A, leading dimension of A, array for eigenvalues, 
c    array for left eigenvalue, leading dimension of DUMMY, 
c    array for right eigenvalues, leading dimension of DUMMY,
c    workspace array dim>=2*order of A, dimension of WORK
c    workspace array dim=2*order of A, return value 
c
c output of eigenvalues
	if (ok .eq. 0) then
	   do i=1, 3
	      write(*,*) b(i)
	   enddo
	else
	   write (*,*) "An error occured"
	endif
	end