SUBROUTINE eigdrv3d()

  USE ZnaupdClass 
  USE ZneupdClass
  use local2d
  use local3d
  USE global

  implicit none

!=============================================================================!
    LOGICAL :: doit=.false.
    INTEGER :: debug=1
    INTEGER :: i, j, idof, itr 
    TYPE(ZnaupdIO) :: myZnaupd
    TYPE(ZneupdIO) :: myZneupd
    REAL, ALLOCATABLE :: rwork(:), rd(:,:)
    REAL :: dznrm2, dlapy2, sor
    EXTERNAL dznrm2, zaxpy, dlapy2
    character*80 :: code='EigDrv3D$'

    complex, allocatable :: v(:,:,:), vold(:,:,:), v2old(:,:,:), r(:,:,:)
    real, allocatable :: r_r(:,:,:), r_im(:,:,:)
    real, allocatable :: Uc_r(:,:,:), Uc_im(:,:,:)
    real, allocatable :: Uh_r(:,:,:), Uh_im(:,:,:)
    real, allocatable :: rmatx(:,:,:,:,:), rmaty(:,:,:,:,:)
    real, allocatable :: MiUc_r(:,:,:), MiUh_im(:,:,:)
!    complex, allocatable :: matx(:,:,:,:,:)
    real, allocatable :: perx(:,:,:,:,:), per2x(:,:,:,:,:)
!    complex, allocatable :: maty(:,:,:,:,:)
    real, allocatable :: pery(:,:,:,:,:), per2y(:,:,:,:,:)
    real, allocatable    :: dtl(:,:)
    complex, allocatable :: bmat(:,:,:,:)
    integer, allocatable :: ipiv(:,:,:), info(:,:)

    integer :: ier, i1, i2, mem=0
    complex :: dummy
    real :: resnrm, resnrm1, resold=1.0e10, fact, ifact, resmax, resnod
    integer :: ixrmax, iyrmax

!.... allocate local time step

    allocate (dtl(nx,ny), STAT=ier)
    if (ier .ne. 0) call error(code,'Insufficient Memory for dtl$')
    dtl  = one
    alfa = one

!.... allocate the storage area for the disturbance field

    allocate (v(ndof,nx,ny),    vold(ndof,nx,ny), r(ndof,nx,ny), &
            r_r(ndof,nx,ny),    r_im(ndof,nx,ny),   &
           Uc_r(ndof,nx,ny),   Uc_im(ndof,nx,ny), &
           Uh_r(ndof,nx,ny),   Uh_im(ndof,nx,ny), &
         MiUc_r(ndof,nx,ny), MiUh_im(ndof,nx,ny), STAT=ier)
    if (ier .ne. 0) call error(code,'Insufficient Memory for v$')
    v = zero

!.... allocate memory for the LHS

    allocate (rmatx(5,ndof,ndof,nx,ny), rmaty(5,ndof,ndof,nx,ny), STAT=ier)
    if (ier .ne. 0) call error(code,'Insufficient Memory for mat$')

!.... allocate memory for the per

    if (xper .or. yper) then
      allocate (perx(2,ndof,ndof,nx,ny), per2x(2,ndof,ndof,nx,ny), &
                pery(2,ndof,ndof,nx,ny), per2y(2,ndof,ndof,nx,ny), &
                STAT=ier)
      if (ier .ne. 0) call error(code, 'Insufficient Memory for per$')
    end if
    
!... Setup myZnaupd & myZneupd

       CALL StandardZnaupdSetup(myZnaupd, ndof*nx*ny)
       CALL StandardZneupdSetup(myZneupd, myZnaupd)

       ALLOCATE( rwork(myZnaupd%ncv), rd(myZnaupd%ncv, 3) )

!... Change the initial residual so that it is smooth and satisfies the BC's

       if (.false.) then
         myZnaupd%info = 1
         !$omp parallel do private(i,idof,i1)
         do j = 1, ny
           do i = 1, nx
             do idof = 1, ndof
               i1 = idof + ndof*(i-1) + ndof*nx*(j-1)
               myZnaupd%resid(i1) = exp(-((x(i,j)-0.0)**2+(y(i,j)-10.0)**2)/3**2 )
!                                   exp( im * x(i,j) )
             end do
           end do
         end do
       end if

!... Reverse Communication (i.e. houston problem a have we)

       itr = 0

 10    CONTINUE 
         itr = itr+1

         CALL callZnaupd(myZnaupd)

         write(*,"(/,'==> ARPack Iteration ',i6)") itr

         IF (myZnaupd%ido .eq. -1 .or. myZnaupd%ido .eq. 1) THEN 

           if (myZnaupd%iparam(7) == 1) then  ! regular mode

!....  Regular mode:  A x = y

             CALL lrhs3d(myZnaupd%workd(myZnaupd%ipntr(1)), &
                         myZnaupd%workd(myZnaupd%ipntr(2)), vm, x, y)
             CALL rhsbc3d( myZnaupd%workd(myZnaupd%ipntr(2)), &
                           myZnaupd%workd(myZnaupd%ipntr(1)), vm )

           else if (myZnaupd%iparam(7) == 3) then ! shift and invert mode
             
!.... Iteratively solve (A - sigma*I) y = x
!....                   (A - i omega I) y = x
!....                   y = (A - i omega I)^-1 x

!.... compute the CFL

             resold = 1.0e10
             doit = .false.
             istep = 0
             call set_time_step(vm, dtl)
             alfa = one

!.... Set the initial condition (zero seems best)

!!$             do i = 1, myZnaupd%n
!!$               i1 = myZnaupd%ipntr(1) + (i-1)
!!$               i2 = myZnaupd%ipntr(2) + (i-1)
!!$               v(i,1,1) = myZnaupd%workd(i1)
!!$             end do

             !$omp parallel do private(i,idof)
             do j = 1, ny
               do i = 1, nx
                 do idof = 1, ndof
                   v(idof,i,j) = zero ! cos( x(i,j) )
                 end do
               end do
             end do

             call itrbc3D(v,vm)
             time = zero

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  Approximate Factorization
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

             call adi(v, MyZnaupd%workd(MyZnaupd%ipntr(1)), dtl) 
           end if
           go to 10
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

             CALL lrhs3d(myZnaupd%v(1,j), r, vm, x, y)
             CALL eigBC3d(r, myZnaupd%v(1,j), vm )

             CALL zaxpy(myZnaupd%n, -myZneupd%d(j),   & 
                        myZnaupd%v(1,j), 1, r, 1)
             rd(j,1) = dble(myZneupd%d(j))
             rd(j,2) = dimag(myZneupd%d(j))
             rd(j,3) = dznrm2(myZnaupd%n, r, 1)
             rd(j,3) = rd(j,3) / dlapy2(rd(j,1),rd(j,2))
 20        CONTINUE 
           CALL dmout(10, myZnaupd%iparam(5), 3, rd, myZnaupd%ncv, -6, &
                      'Ritz values (Real, Imag) and relative residuals')
           CALL dmout( 6, myZnaupd%iparam(5), 3, rd, myZnaupd%ncv, -6, &
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

!.... If shift and invert then lambda = sigma + 1/amu

      if (myZnaupd%iparam(7) == 3) then
        myZneupd%d(1:myZnaupd%nev) = im*omega + myZneupd%d(1:myZnaupd%nev)
      end if

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


!=============================================================================
        subroutine ebe3D( rl, vl, voldl, dtl, fl )
!
!  Backward Euler RHS:    r = vold - v - delt * (r - f)
!
!=============================================================================
        use global
        implicit none
        
        complex :: rl(ndof,nx,ny), vl(ndof,nx,ny), voldl(ndof,nx,ny)
        complex :: fl(ndof,nx,ny)
        real    :: dtl(nx,ny)
        integer :: i,j,idof
!=============================================================================

        !$omp parallel do private(i,idof)
        do j = 1, ny
          do i = 1, nx
            do idof = 1, ndof
             !rl(idof,i,j) = - dtl(i,j) * (rl(idof,i,j) - fl(idof,i,j)) 
             rl(idof,i,j) = voldl(idof,i,j) - vl(idof,i,j) - &
                            dtl(i,j) * (rl(idof,i,j) - fl(idof,i,j))
            end do
          end do
        end do
        
        return
        end

!=============================================================================!
        subroutine blkdiag3D( mat, Dh, Dhi, Vh11, Vh22 )
!  
!  For the block diagonal matrix
!
!=============================================================================!
        use global
        use stencil
        implicit none

        complex :: mat(ndof,ndof,nx,ny)
        real    :: Dh(ndof,ndof,nx,ny), Dhi(ndof,ndof,nx,ny)
        real    :: Vh11(6,nx,ny), Vh22(6,nx,ny)

        real    :: a3, b3

        integer :: i, j, idof, jdof
!=============================================================================!

!.... fourth-order stencil

        a3 = da3 / dxi**2
        b3 = da3 / deta**2

        !$omp parallel do private (i,idof,jdof)
        loop_j: do j = 1, ny
        loop_i: do i = 1, nx

!.... \hat{D} term

        do idof = 1, ndof
          do jdof = 1, ndof
            mat(idof,jdof,i,j) = Dh(idof,jdof,i,j) + im * Dhi(idof,jdof,i,j)
          end do
        end do

!.... \hat{V}_{\xi\xi} term

        mat(2,2,i,j) = mat(2,2,i,j) - a3 * Vh11(1,i,j)
        mat(2,3,i,j) = mat(2,3,i,j) - a3 * Vh11(5,i,j)
        mat(3,2,i,j) = mat(3,2,i,j) - a3 * Vh11(6,i,j)
        mat(3,3,i,j) = mat(3,3,i,j) - a3 * Vh11(2,i,j)
        mat(4,4,i,j) = mat(4,4,i,j) - a3 * Vh11(3,i,j)
        mat(5,5,i,j) = mat(5,5,i,j) - a3 * Vh11(4,i,j)

!.... \hat{v}_{\xi\eta} term

!.... \hat{V}_{\eta\eta} term

        mat(2,2,i,j) = mat(2,2,i,j) - b3 * Vh22(1,i,j)
        mat(2,3,i,j) = mat(2,3,i,j) - b3 * Vh22(5,i,j)
        mat(3,2,i,j) = mat(3,2,i,j) - b3 * Vh22(6,i,j)
        mat(3,3,i,j) = mat(3,3,i,j) - b3 * Vh22(2,i,j)
        mat(4,4,i,j) = mat(4,4,i,j) - b3 * Vh22(3,i,j)
        mat(5,5,i,j) = mat(5,5,i,j) - b3 * Vh22(4,i,j)

        end do loop_i
        end do loop_j

        end subroutine blkdiag3D

!=============================================================================
        subroutine MiMultBy(vl,rl, dtl)
!
!
!=============================================================================
        use global
        use local3d
        use local2d
        implicit none

        real :: rl(ndof,nx,ny), vl(ndof,nx,ny)
        real :: dtl(nx,ny)
        integer :: i,j
!=============================================================================

        rl(:,:,:) = zero

!.... compute first derivatives

        call grad( ndof, nx, ny, vl, g1v, g2v, dxi, deta, optx, opty, &
                    xper, yper, lsym, rsym, bsym, tsym, carp )

        loop_j: do j = 1, ny
        loop_i: do i = 1, nx

!.... \hat{A}_i, \hat{B}_i terms

        rl(2,i,j) = rl(2,i,j) - ( ABhi(1,i,j) * g1v(4,i,j) + &
                                  ABhi(3,i,j) * g2v(4,i,j) )
        rl(3,i,j) = rl(3,i,j) - ( ABhi(2,i,j) * g1v(4,i,j) + &
                                  ABhi(4,i,j) * g2v(4,i,j) )
        rl(4,i,j) = rl(4,i,j) - ( ABhi(1,i,j) * g1v(2,i,j) + &
                                  ABhi(2,i,j) * g1v(3,i,j) + &
                    ABhi(3,i,j) * g2v(2,i,j) + ABhi(4,i,j) * g2v(3,i,j) )


!.... \hat{D}_i term

        rl(1,i,j) = rl(1,i,j) + ( Dhi(1,1,i,j) * vl(1,i,j)      &
                                 + Dhi(1,2,i,j) * vl(2,i,j)     &
                                 + Dhi(1,3,i,j) * vl(3,i,j)     &
                                 + Dhi(1,4,i,j) * vl(4,i,j)     &
                                 + Dhi(1,5,i,j) * vl(5,i,j) )

        rl(2,i,j) = rl(2,i,j) + ( Dhi(2,1,i,j) * vl(1,i,j)      &
                                 + Dhi(2,2,i,j) * vl(2,i,j)     &
                                 + Dhi(2,3,i,j) * vl(3,i,j)     &
                                 + Dhi(2,4,i,j) * vl(4,i,j)     &
                                 + Dhi(2,5,i,j) * vl(5,i,j) )

        rl(3,i,j) = rl(3,i,j) + ( Dhi(3,1,i,j) * vl(1,i,j)      &
                                 + Dhi(3,2,i,j) * vl(2,i,j)     &
                                 + Dhi(3,3,i,j) * vl(3,i,j)     &
                                 + Dhi(3,4,i,j) * vl(4,i,j)     &
                                 + Dhi(3,5,i,j) * vl(5,i,j) )

        rl(4,i,j) = rl(4,i,j) + ( Dhi(4,1,i,j) * vl(1,i,j)      &
                                 + Dhi(4,2,i,j) * vl(2,i,j)     &
                                 + Dhi(4,3,i,j) * vl(3,i,j)     &
                                 + Dhi(4,4,i,j) * vl(4,i,j)     &
                                 + Dhi(4,5,i,j) * vl(5,i,j) )

        rl(5,i,j) = rl(5,i,j) + ( Dhi(5,1,i,j) * vl(1,i,j)      &
                                 + Dhi(5,2,i,j) * vl(2,i,j)     &
                                 + Dhi(5,3,i,j) * vl(3,i,j)     &
                                 + Dhi(5,4,i,j) * vl(4,i,j)     &
                                 + Dhi(5,5,i,j) * vl(5,i,j) )

!... add the time term

        rl(1,i,j) = dtl(i,j) * rl(1,i,j)
        rl(2,i,j) = dtl(i,j) * rl(2,i,j)
        rl(3,i,j) = dtl(i,j) * rl(3,i,j)
        rl(4,i,j) = dtl(i,j) * rl(4,i,j)
        rl(5,i,j) = dtl(i,j) * rl(5,i,j)

        end do loop_i
        end do loop_j

        return
        end
