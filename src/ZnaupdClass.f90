MODULE ZnaupdClass 
TYPE :: ZnaupdIO
INTEGER :: ido           ! Reverse Communication Flag
CHARACTER :: bmat*1      ! Specifies type of eigenvalue problem
INTEGER :: n             ! Dimension of eigenvalue problem
CHARACTER :: which*2     ! Type of eigenvalues sought
INTEGER :: nev           ! Number of eigenvalues to be computed
REAL :: tol              ! Stopping criteria
COMPLEX, POINTER, &
     DIMENSION(:) :: resid ! final residual vector
INTEGER :: ncv           ! Number of columns of the matrix V
COMPLEX, POINTER, &
     DIMENSION(:,:) :: V ! Contains final set of Arnoldi Basis Vectors
INTEGER :: ldv           ! leading dimension of v
INTEGER :: iparam(11)    ! settings
INTEGER :: ipntr(14)
COMPLEX, POINTER, DIMENSION(:) :: workd ! work array of length 3*N
COMPLEX, POINTER, DIMENSION(:) :: workl ! work array of length lworkl
INTEGER :: lworkl
REAL, POINTER, DIMENSION(:) :: rwork    ! work array
INTEGER :: info          ! error flag on output
END TYPE ZnaupdIO
CONTAINS

!... Setup for Standard Eigenvalue Problem

SUBROUTINE StandardZnaupdSetup(myZnaupd, n)
IMPLICIT NONE
TYPE(ZnaupdIO), INTENT(INOUT) :: myZnaupd
INTEGER :: n, maxn, maxnev, maxncv
INTEGER :: mode, maxitr, ishfts, nev, ncv
CHARACTER :: bmat
CHARACTER(len=2) :: which
REAL :: tol

NAMELIST /inputparam/ mode, maxitr, bmat, which, nev, ncv, tol

write(*,*) 'Beginning ZnaupdClass Setup...'
write(*,*)
write(*,*)
write(*,*) '**************************************************'
write(*,*) 'Welcome and Enjoy Using and Re-using ZnaupdClass.'
write(*,*) '**************************************************'
write(*,*)

mode = 1  
maxitr = 100
ishfts = 1

myZnaupd%ido = 0
myZnaupd%bmat = 'I' 
myZnaupd%n = n
myZnaupd%which = 'LI'
myZnaupd%nev = 8 
myZnaupd%tol = 0.000005
myZnaupd%ncv = 40 
myZnaupd%ldv = n
myZnaupd%iparam(:) = 0
myZnaupd%iparam(1) = ishfts
myZnaupd%iparam(3) = maxitr
myZnaupd%iparam(7) = mode
myZnaupd%lworkl =  3*(myZnaupd%ncv)**2+5*(myZnaupd%ncv)
myZnaupd%info = 0 ! use random initial vector

OPEN(7, FILE='eig.inp', DELIM='APOSTROPHE')
READ(7, NML=inputparam)
WRITE(*,*)
WRITE(*,*) '***************'
WRITE(*,*) 'Namelist Input:'
WRITE(*,*) '***************'
WRITE(*,*)
WRITE(*,*)
WRITE(UNIT=*, NML=inputparam)
WRITE(*,*)
WRITE(*,*)

myZnaupd%bmat = bmat
myZnaupd%which = which
myZnaupd%nev = nev
myZnaupd%tol = tol
myZnaupd%ncv = ncv
myZnaupd%iparam(3) = maxitr
myZnaupd%iparam(7) = mode
myZnaupd%lworkl =  3*(myZnaupd%ncv)**2+5*(myZnaupd%ncv)

ALLOCATE( myZnaupd%resid(myZnaupd%n),                       &
          myZnaupd%V(myZnaupd%n,myZnaupd%ncv),              &
          myZnaupd%workd(3*myZnaupd%n),                     &
          myZnaupd%workl(3*myZnaupd%ncv**2+5*myZnaupd%ncv), &
          myZnaupd%rwork(myZnaupd%ncv) )

END SUBROUTINE StandardZnaupdSetup

SUBROUTINE CallZnaupd(myZnaupd)
IMPLICIT NONE
TYPE(ZnaupdIO), INTENT(INOUT) :: myZnaupd 
CALL znaupd (myZnaupd%ido, myZnaupd%bmat, myZnaupd%n, myZnaupd%which, &
     myZnaupd%nev, myZnaupd%tol, myZnaupd%resid, myZnaupd%ncv,        &
     myZnaupd%v, myZnaupd%ldv, myZnaupd%iparam, myZnaupd%ipntr,       &
     myZnaupd%workd, myZnaupd%workl, myZnaupd%lworkl,                 &
     myZnaupd%rwork, myZnaupd%info )
END SUBROUTINE CallZnaupd

END MODULE  ZnaupdClass
