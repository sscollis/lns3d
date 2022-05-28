program solve

      implicit none

      integer :: nx, ny, nz, n
      real, allocatable :: A(:,:), D1(:,:), D2(:,:), rhs(:,:,:)
      integer :: info, j
      integer, allocatable :: ipiv(:)
      real, external :: cpu, second

      integer :: task, ntasks, lot
      integer, allocatable :: taskarray(:,:), mark(:,:)
      external SGETRS

      write(*,"('Enter (nx,ny,nz) ==> ',$)")
      read(*,*) nx, ny, nz
      allocate( A(ny,ny), D1(ny,ny), D2(ny,ny), rhs(ny,nx,nz), ipiv(ny) )
 
      write(*,"('Enter ntasks ==> ',$)")
      read(*,*) ntasks
      allocate( taskarray(3,ntasks), mark(4,ntasks) )

      do task=1,ntasks
	taskarray(1,task)=3
	taskarray(3,task)=task
      end do

      lot = nx*nz / ntasks
      do task = 1, ntasks-1
	mark(1,task) = 1 + lot * (task-1)
	mark(2,task) = 1 + lot * ( task ) - 1
	mark(3,task) = mark(2,task) - mark(1,task) + 1
!	write(*,"(4(i5,1x))") task, mark(1,task), mark(2,task), mark(3,task)
      end do
      task = ntasks
      mark(1,ntasks) = 1 + lot * (ntasks-1)
      mark(2,ntasks) = nx*nz
      mark(3,ntasks) = mark(2,ntasks) - mark(1,ntasks) + 1
!     write(*,"(4(i5,1x))") task, mark(1,task), mark(2,task), mark(3,task)

      call chebyd(D1, ny-1)	! use ny-1 since this routine uses 0:ny
      D2 = matmul(D1, D1)	! compute the second derivative

      do n = 1, 50

      write(*,"(i5,1x,1pe13.6)") n, second()

      A = -2.0 * D2 + 2.0 * D1

      do j = 1, ny
	A(j,j) = A(j,j) + 1.0
      end do

      call random_number( rhs )

      call SGETRF( ny, ny, A, ny, ipiv, info )

      do task = 2, ntasks
	call tskstart( taskarray(1,task), SGETRS, 'N', ny, mark(3,task), A, ny, ipiv, rhs(1,mark(1,task),1), ny, mark(4,task) )
      end do
      call SGETRS( 'N', ny, mark(3,1), A, ny, ipiv, rhs(1,mark(1,1),1), ny, mark(4,1) )
      
      do task = 2, ntasks
	call tskwait(taskarray(1,task))
      end do

      end do

      stop
end program solve

real function cpu()
      real, save :: time = 0.0
      real, external :: second

      time = time + second()
      cpu = time
      return
end function cpu
!==============================================================================
	subroutine CHEBYD(D, N)
!==============================================================================
!
!     	Calculation the Chebyshev collocation derivative matrix of order N
!	at the Gauss-Labatto points
!
!	This is a brute force calculation and is actually valid for arbitrary
!	quadrature points.  Just change the mesh definition.
!
!	Since this routine is brute force, the chance for round-off errors
!	is significantly increased.  So watch out for large N
!
!	Revised: 2-9-96
!==============================================================================
	implicit none
	
	integer :: N, I, J, K
	real :: D(0:N,0:N), PI, X(0:N), a(0:N)
	real, parameter :: zero=0.0, pt5=0.5, one=1.0, two=2.0, six=6.0
!==============================================================================
!
!.... compute PI
!
	pi = acos(-one)
!
!.... form the mesh for the Chebyshev-Gauss-Lobatto points
!
	do j = 0, N
	  x(j) = COS( PI * REAL(j) / REAL(N) )
	end do
!
!.... compute the a(j)
!
	do j = 0, N
	  a(j) = one
	  do k = 0, N
	    if (k.ne.j) a(j) = a(j) * ( x(j) - x(k) )
	  end do
	end do
!
!.... form the diagonal of the matrix
!
	do j = 0, N
	  D(j,j) = zero
	  do k = 0, N
	    if (k.ne.j) D(j,j) = D(j,j) + one / ( x(j) - x(k) )
	  end do
	end do
!
!.... Form the off-diagonal terms
!            
	do j = 0, N
	  do k = 0, N
	    if (k.ne.j) D(j,k) = a(j) / ( a(k) * (x(j) - x(k)) )
	  end do
	end do
	
	return
	end

!==============================================================================
	subroutine CHEBYD_NEW(D, N)
!==============================================================================
!
!     	Calculation the Chebyshev collocation derivative matrix of order N
!
!	This routine uses the matrix definition given by:
!
!	W.S. Don and A. Solomonoff, (1995) Accuracy and Speed in
!	Computing the Chebyshev Collocation Derivative, SIAM J. Sci. Comput.,
!	Vol 16, No. 6. pp. 1253-1268.
!
!	This method uses trigonometric functions and flipping to reduce
!	roundoff error in the matrix.  They show that with this method
!	the derivative is as accurate as that using FFT's.
!
!	In a practical problem (computing the derivatives of a boundary
!	layer profile) I have found no perceivable accuracy difference 
!	between this method and the original method.
!
!	Revised: 2-9-96
!==============================================================================
	implicit none
	
	integer :: N, I, J, K
	real :: D(0:N,0:N), PI, C(0:N), X(0:N)
	real, parameter :: zero=0.0, pt5=0.5, one=1.0, two=2.0, six=6.0
!==============================================================================
!
!.... compute PI
!
	PI = acos(-one)
!
!.... compute normalization
!      
	do I = 1, N-1
	  C(I) = one
	end do
	C(0) = two
	C(N) = two
!
!.... form the mesh
!
	DO J = 0, N
	  X(J) = COS( PI * REAL(J) / REAL(N) )
	END DO
!
!.... form the derivative
!            
	do J = 0, N
	  do K = 0, N
	    if ( J.EQ.0 .AND. K.EQ.0) then
  
	      D(J,K) = (two * real(N)**2 + one) / six
  
	    else if ( J.EQ.N .AND. K.EQ.N) then
  
	      D(J,K) = -(two * real(N)**2 + one) / six
  
	    else if (J.EQ.K) then
  
	      D(J,K) = -pt5 * x(j) / (sin(pi*real(j)/real(N)))**2
  
	    else
  
	      D(J,K) = -pt5 * C(J) / C(k) * (-one)**(J+K) / &
			( sin(pt5*pi*real(j+k)/real(n)) * &
			  sin(pt5*pi*real(j-k)/real(n)) )
  
	    end if
	  end do
	end do
!
!.... correct the bottom half of the matrix by flipping
!
	do i = N/2+1, N
	  do j = 1, N
	    D(i,j) = -D(N-i,N-j)
	  end do
	end do
	
	return
	end

!==============================================================================
	subroutine CHEBYINT(GLtoG, GtoGL, N)
!==============================================================================
!
!	Calculation the Chebyshev interpolation matrices from Gauss-Lobatto 
!	to Gauss points, and from Gauss to Gauss-Lobatto points.
!
!==============================================================================
	implicit none
  
	integer :: N
	real :: GLtoG(0:N-1,0:N), GtoGL(0:N,0:N-1)
	real :: CGL(0:N,0:N),     CGLinv(0:N,0:N-1)
	real :: CG(0:N-1,0:N-1),  CGinv(0:N-1,0:N)
  
	real :: PI, C(0:N), Cb(0:N)
	real :: fact1, fact2
  
	real, parameter :: zero=0.0, pt5=0.5, one=1.0, two=2.0
  
	integer :: i, j
!==============================================================================
!
!.... compute PI
!
	PI = acos(-one)
!
!.... compute normalization coefficients
!      
	do I = 1, N-1
	  C(I)  = one
	  Cb(I) = one
	end do
	C(0)  = two
	Cb(0) = two
	Cb(N) = two
!
!.... Form the Chebyschev coefficients at the Gauss-Lobatto points
!
	fact1 = two / real(N)
	do i = 0, N
	  fact2 = fact1 / Cb(i)
	  do j = 0, N
	    CGL(i,j) = fact2 / Cb(j) * cos(pi*real(i)*real(j)/real(N))
	  end do
	end do
!
!.... Form the Inverse Chebyschev operator at the Gauss-Lobatto points
!.... Given the Chebyschev Gauss coefficients
!      
	do i = 0, N
	  do j = 0, N-1
	    CGLinv(i,j) = cos(pi*real(i)*real(j)/real(N))
	  end do
	end do
!
!.... Form the Chebyschev coefficients at the Gauss points
!
	fact1 = two / real(N)
	do i = 0, N-1
	  fact2 = fact1 / C(i)
	  do j = 0, N-1
	    CG(i,j) = fact2 * cos(pi*real(i)*(real(j)+pt5)/real(N))
	  end do
	end do
!
!.... Form the Inverse Chebyschev operator at the Gauss points
!.... given the Chebyschev Gauss-Labatto coefficients
!      
	do i = 0, N-1
	  do j = 0, N
	    CGinv(i,j) = cos(pi*(real(i)+pt5)*real(j)/real(N))
	  end do
	end do
!
!.... Interpolation matrix from Gauss-Lobatto to Gauss points
!
	GLtoG = matmul( CGinv, CGL )

!	write(*,*) 'GLtoG'
!	do i = 0, N-1
!	  write(*,10) (GLtoG(i,j),j=0,N)
!	end do
!
!.... Interpolation matrix from Gauss to Gauss-Lobatto points
!
	GtoGL = matmul( CGLinv, CG )

!	write(*,*) 'GtoGL'
!	do i = 0, N
!	  write(*,10) (GtoGL(i,j),j=0,N-1)
!	end do

	return
 10	format( 10(1pe13.6,1x) )
	end

!==============================================================================
	subroutine CHEBYDG(D, N)
!==============================================================================
!
!     	Calculation the Chebyshev collocation derivative matrix of order N
!	at the Gauss collocation points
!
!	This is a brute force calculation and is actually valid for arbitrary
!	quadrature points.  Just change the mesh definition.
!
!	Since this routine is brute force, the chance for round-off errors
!	is significantly increased.  So watch out for large N
!
!	Revised: 2-9-96
!==============================================================================
	implicit none
	
	integer :: N, I, J, K
	real :: D(0:N,0:N), PI, X(0:N), a(0:N)
	real, parameter :: zero=0.0, pt5=0.5, one=1.0, two=2.0, six=6.0
!==============================================================================
!
!.... compute PI
!
	pi = acos(-one)
!
!.... form the mesh for the Chebyshev-Gauss points
!
	do j = 0, N
	  x(j) = COS( pi * (two * real(j) + one) / (two * real(N) + 2) )
	end do
!
!.... form the mesh for the Chebyshev-Gauss-Lobatto points (test)
!
!	do j = 0, N
!	  x(j) = COS( PI * REAL(j) / REAL(N) )
!	end do
!
!.... compute the a(j)
!
	do j = 0, N
	  a(j) = one
	  do k = 0, N
	    if (k.ne.j) a(j) = a(j) * ( x(j) - x(k) )
	  end do
	end do
!
!.... form the diagonal of the matrix
!
	do j = 0, N
	  D(j,j) = zero
	  do k = 0, N
	    if (k.ne.j) D(j,j) = D(j,j) + one / ( x(j) - x(k) )
	  end do
	end do
!
!.... Form the off-diagonal terms
!            
	do j = 0, N
	  do k = 0, N
	    if (k.ne.j) D(j,k) = a(j) / ( a(k) * (x(j) - x(k)) )
	  end do
	end do
!
!.... correct the bottom half of the matrix by flipping.  (Is this valid?)
!
!	do i = N/2+1, N
!	  do j = 1, N
!	    D(i,j) = -D(N-i,N-j)
!	  end do
!	end do
	
	return
	end

