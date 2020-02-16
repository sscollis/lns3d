SUBROUTINE eigdrv3d()

  USE ZnaupdClass 
  USE ZneupdClass
  USE global

  implicit none

!=============================================================================!

    INTEGER :: j, itr 
    TYPE(ZnaupdIO) :: myZnaupd
    TYPE(ZneupdIO) :: myZneupd
    COMPLEX, ALLOCATABLE :: ax(:)
    REAL, ALLOCATABLE :: rwork(:), rd(:,:)
    REAL :: dznrm2, dlapy2
    EXTERNAL dznrm2, zaxpy, dlapy2

!    include 'debug.h'
!     ndigit = -3
!     logfil = 6
!     mcaitr = 0 
!     mcapps = 0
!     mcaupd = 0
!     mcaup2 = 2
!     mceigh = 0
!     mceupd = 0

!... Setup myZnaupd & myZneupd

       CALL StandardZnaupdSetup(myZnaupd, ndof*nx*ny)
       CALL StandardZneupdSetup(myZneupd, myZnaupd)

       ALLOCATE( rwork(myZnaupd%ncv), rd(myZnaupd%ncv, 3), ax(myZnaupd%n) )

!... Reverse Communication (i.e. houston problem a have we)

       itr = 0

 10    CONTINUE 
         itr = itr+1

         CALL callZnaupd(myZnaupd)

         write(*,*) itr ! , myZnaupd%iparam

         IF (myZnaupd%ido .eq. -1 .or. myZnaupd%ido .eq. 1) THEN 

         !   myZnaupd%workd(myZnaupd%ipntr(2)) = v (output)
         !   myZnaupd%workd(myZnaupd%ipntr(1)) = r (input)

!         CALL itrbc3d( myZnaupd%workd(myZnaupd%ipntr(1)), vm )

         CALL lrhs3d(myZnaupd%workd(myZnaupd%ipntr(1)), &
                     myZnaupd%workd(myZnaupd%ipntr(2)), vm, x, y)

         CALL rhsbc3d( myZnaupd%workd(myZnaupd%ipntr(2)), &
                       myZnaupd%workd(myZnaupd%ipntr(1)), vm )

!         CALL eigbc3d( myZnaupd%workd(myZnaupd%ipntr(2)), &
!                       myZnaupd%workd(myZnaupd%ipntr(1)), vm )

!        CALL av (nx, myZnaupd%workd(myZnaupd%ipntr(1)), &
!                     myZnaupd%workd(myZnaupd%ipntr(2)))
         GO TO 10
         END IF

!... Post-Processing 
         IF (myZnaupd%info .lt. 0) THEN 
           PRINT *, ' '
           PRINT *, ' Error with _naupd, info = ', myZnaupd%info 
           PRINT *, ' Check the documentation of _naupd'
           PRINT *, ' '
         ELSE
           CALL callZneupd(myZneupd, myZnaupd) 
           IF (myZnaupd%info .ne. 0) THEN 
             PRINT *, ' '
             PRINT *, ' Error with _neupd, info = ', myZnaupd%info 
             PRINT *, ' Check the documentation of _neupd. '
             PRINT *, ' '
           ELSE
           DO 20 j = 1, myZnaupd%iparam(5)

!            CALL av(nx, myZnaupd%v(1,j), ax)

             CALL lrhs3d(myZnaupd%v(1,j), ax, vm, x, y)
             CALL eigBC3d(ax, myZnaupd%v(1,j), vm )

             CALL zaxpy(myZnaupd%n, -myZneupd%d(j),   & 
                        myZnaupd%v(1,j), 1, ax, 1)
             rd(j,1) = dble(myZneupd%d(j))
             rd(j,2) = dimag(myZneupd%d(j))
             rd(j,3) = dznrm2(myZnaupd%n, ax, 1)
             rd(j,3) = rd(j,3) / dlapy2(rd(j,1),rd(j,2))
 20        CONTINUE 

             CALL dmout(10, myZnaupd%iparam(5), 3, rd, myZnaupd%ncv, -6, &
                      'Ritz values (Real, Imag) and relative residuals')
         END IF 

         IF ( myZnaupd%info .eq. 1) THEN 
             PRINT *, ' '
             PRINT *, ' Maximum number of iterations reached.'
             PRINT *, ' '
         ELSE IF ( myZnaupd%info .eq. 3) THEN
             PRINT *, ' '
             PRINT *, ' No shifts could be applied during implicit', &
                     ' Arnoldi update, try increasing NCV.'
             PRINT *, ' '
         END IF 

         PRINT *, ' '
         PRINT *, '_NDRV1'
         PRINT *, '====== '
         PRINT *, ' '
         PRINT *, ' Size of the matrix is ', myZnaupd%n
         PRINT *, ' The number of Ritz values requested is ', myZnaupd%nev
         PRINT *, ' The number of Arnoldi vectors generated', &
                ' (NCV) is ', myZnaupd%ncv
         PRINT *, ' What portion of the spectrum: ', myZnaupd%which
         PRINT *, ' The number of converged Ritz values is ', &
                   myZnaupd%iparam(5)
         PRINT *, ' The number of Implicit Arnoldi update', &
                  ' iterations taken is ', myZnaupd%iparam(3)
         PRINT *, ' The number of OP*x is ', myZnaupd%iparam(9)
         PRINT *, ' The convergence criterion is ', myZnaupd%tol
         PRINT *, ' '
      END IF

!.... Output the eigenvalues and eigenvectors

      open(11,file='eig.dat', form='unformatted', status='unknown')
      write(11) myZnaupd%nev, 0.0, nx, ny, nz, ndof, &
                Re, Ma, Pr, gamma, cv
      write(11) myZneupd%d(1:myZnaupd%nev) 
      write(11) myZnaupd%v(:,:)
      close(11)

END SUBROUTINE 

!==========================================================================
!
!     matrix vector subroutine
!
!     The matrix used is the convection-diffusion operator
!     discretized using centered difference.
!
      subroutine av (nx, v, w)
      integer           nx, j, lo
      Complex*16        v(nx*nx), w(nx*nx), one, h2
      parameter         (one = (1.0D+0, 0.0D+0))
      external          zaxpy, tv
!
!     Computes w <--- OP*v, where OP is the nx*nx by nx*nx block
!     tridiagonal matrix
!
!                  | T -I          |
!                  |-I  T -I       |
!             OP = |   -I  T       |
!                  |        ...  -I|
!                  |           -I T|
!
!     derived from the standard central difference  discretization
!     of the convection-diffusion operator (Laplacian u) + rho*(du/dx)
!     with zero boundary condition.
!
!     The subroutine TV is called to computed y<---T*x.
!
!
      h2 = one / dcmplx((nx+1)*(nx+1))
!
      call tv(nx,v(1),w(1))
      call zaxpy(nx, -one/h2, v(nx+1), 1, w(1), 1)
!
      do 10 j = 2, nx-1
         lo = (j-1)*nx
         call tv(nx, v(lo+1), w(lo+1))
         call zaxpy(nx, -one/h2, v(lo-nx+1), 1, w(lo+1), 1)
         call zaxpy(nx, -one/h2, v(lo+nx+1), 1, w(lo+1), 1)
  10  continue
!
      lo = (nx-1)*nx
      call tv(nx, v(lo+1), w(lo+1))
      call zaxpy(nx, -one/h2, v(lo-nx+1), 1, w(lo+1), 1)
!
      return
      end subroutine
!=========================================================================

!=========================================================================
      subroutine tv (nx, x, y)
!
      integer           nx, j
      Complex*16  x(nx), y(nx), h, h2, dd, dl, du
!
      Complex*16  one, rho
      parameter         (one = (1.0D+0, 0.0D+0), rho = (1.0D+2, 0.0D+0))
!
!     Compute the matrix vector multiplication y<---T*x
!     where T is a nx by nx tridiagonal matrix with DD on the
!     diagonal, DL on the subdiagonal, and DU on the superdiagonal
!
      h   = one / dcmplx(nx+1)
      h2  = h*h
      dd  = (4.0D+0, 0.0D+0) / h2
      dl  = -one/h2 - (5.0D-1, 0.0D+0)*rho/h
      du  = -one/h2 + (5.0D-1, 0.0D+0)*rho/h
!
      y(1) =  dd*x(1) + du*x(2)
      do 10 j = 2,nx-1
         y(j) = dl*x(j-1) + dd*x(j) + du*x(j+1)
 10   continue
      y(nx) =  dl*x(nx-1) + dd*x(nx)
      return
      end subroutine


