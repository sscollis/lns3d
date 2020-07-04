SUBROUTINE eigdrv3d()

  USE ZnaupdClass 
  USE ZneupdClass
  use local2d
  use local3d
  USE global

  implicit none

!=============================================================================!
    LOGICAL :: gdoit= .true., doit=.true.
    INTEGER :: debug=1
    INTEGER :: i, j, idof, itr 
    TYPE(ZnaupdIO) :: myZnaupd
    TYPE(ZneupdIO) :: myZneupd
    REAL, ALLOCATABLE :: rwork(:), rd(:,:)
    REAL :: dznrm2, dlapy2, sor
    EXTERNAL dznrm2, zaxpy, dlapy2
    character(80) :: code='EigDrv3D$'

    complex, allocatable :: v(:,:,:), vold(:,:,:), v2old(:,:,:), r(:,:,:)
    complex, allocatable :: matx(:,:,:,:,:)
    complex, allocatable :: perx(:,:,:,:,:), per2x(:,:,:,:,:)
    complex, allocatable :: maty(:,:,:,:,:)
    complex, allocatable :: pery(:,:,:,:,:), per2y(:,:,:,:,:)
    real, allocatable    :: dtl(:,:)
    complex, allocatable :: bmat(:,:,:,:)
    integer, allocatable :: ipiv(:,:,:), info(:,:)

    complex, pointer :: vt(:,:,:), rt(:,:,:)
    integer :: ns, ne

    integer :: ier, i1, i2
    complex, pointer :: dummy(:,:,:)
    real :: resnrm, resnrm1, resold=1.0e10, fact, ifact, resmax, resnod
    integer :: ixrmax, iyrmax

!.... allocate local time step

    allocate (dtl(nx,ny), STAT=ier)
    if (ier .ne. 0) call error(code,'Insufficient Memory for dtl$')
    dtl  = one
    alfa = one

!.... allocate the storage area for the disturbance field

    allocate (v(ndof,nx,ny), vold(ndof,nx,ny), r(ndof,nx,ny), STAT=ier)
    if (ier .ne. 0) call error(code,'Insufficient Memory for v$')
    v = zero

!.... allocate memory for the LHS

    allocate (matx(5,ndof,ndof,nx,ny), maty(5,ndof,ndof,nx,ny), STAT=ier)
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
               if (idof.eq.ndof) then
                 myZnaupd%resid(i1) = exp(-((x(i,j)-0.0)**2+ &
                                            (y(i,j)-10.0)**2)/3**2 )
               else
                 myZnaupd%resid(i1) = zero
               end if
             end do
           end do
         end do
       end if

!.... write out restart file

       open(unit=10, file='resid.dat', form='unformatted', &
            status='unknown')
       write(10) itr, time, nx, ny, nz, ndof, Re, Ma, Pr, gamma, cv
       write(10) (myZnaupd%resid(i), i = 1,ndof*nx*ny)
       close(10)

!... Reverse Communication (i.e. houston problem a have we)

       itr = 0

 10    CONTINUE 
         itr = itr+1

         CALL callZnaupd(myZnaupd)

         write(*,"(/,'==> ARPack Iteration ',i6)") itr

         IF (myZnaupd%ido .eq. -1 .or. myZnaupd%ido .eq. 1) THEN 

           if (myZnaupd%iparam(7) == 1) then  ! regular mode

!....  Regular mode:  A x = y

! SSC 200704:  added pointers to remove mismatch
#if 1
             ns = myZnaupd%ipntr(1)
             ne = myZnaupd%ipntr(1)+myZnaupd%n
             vt(1:ndof,1:nx,1:ny) => myZnaupd%workd(ns:ne)
             ns = myZnaupd%ipntr(2)
             ne = myZnaupd%ipntr(2)+myZnaupd%n
             rt(1:ndof,1:nx,1:ny) => myZnaupd%workd(ns:ne)

             CALL lrhs3d(vt, rt, vm, x, y)
             CALL rhsbc3d(rt, vt, vm)
#else
             CALL lrhs3d(myZnaupd%workd(myZnaupd%ipntr(1)), &
                         myZnaupd%workd(myZnaupd%ipntr(2)), vm, x, y)
             CALL rhsbc3d( myZnaupd%workd(myZnaupd%ipntr(2)), &
                           myZnaupd%workd(myZnaupd%ipntr(1)), vm )
#endif

           else if (myZnaupd%iparam(7) == 3) then ! shift and invert mode
             
!.... Iteratively solve (A - sigma*I) y = x
!....                   (A - i omega I) y = x
!....                   y = (A - i omega I)^-1 x

!.... compute the CFL

             resold = 1.0e10
             istep = 0
             doit = gdoit
             call set_time_step(vm, dtl)
             alfa = one

!.... Build and factor the matrix

             call lhs1f3D( matx, Ah, Dh, Dhi, Vh11, ABhi, dtl, dummy )
             if (xper) then
               call cpenta2p( nx, ny, ndof, matx, dummy, perx, per2x, 0)
             else
               call cpenta2bc( nx, ny, ndof, matx, dummy, 0 )
             end if
             
             call lhs2f3D( maty, Bh, Dh, Dhi, Vh22, ABhi, dtl, dummy )
             if (yper) then
               call cpenta1p( nx, ny, ndof, maty, dummy, pery, per2y, 0 )
             else
               call cpenta1bc( nx, ny, ndof, maty, dummy, 0 )
             end if

!.... Set the initial condition (zero seems best)

!!$          do i = 1, myZnaupd%n
!!$            i1 = myZnaupd%ipntr(1) + (i-1)
!!$            i2 = myZnaupd%ipntr(2) + (i-1)
!!$            v(i,1,1) = myZnaupd%workd(i1)
!!$          end do

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
!  Point Jacobi Iteration
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

             if (.false.) then

!.... Form the block diagonal and factor

               allocate( bmat(ndof,ndof,nx,ny), ipiv(ndof,nx,ny), info(nx,ny) )

               call blkdiag3d( bmat, Dh, Dhi, Vh11, Vh22 )
               
               !$omp parallel do private(i), shared(ndof)
               do j = 1, ny
                 do i = 1, nx
                   call ZGETRF( ndof, ndof, bmat(1,1,i,j), ndof, ipiv(1,i,j), &
                                info(i,j) )
                   if (info(i,j) .ne. 0) then
                     write(*,*) i, j, info(i,j)
                   end if
                 end do
               end do

!.... Iteration loop

               do istep = 1, nstep

!.... Evaluate residual

               call lrhs3D( v, r, vm, x, y )
               !$omp parallel do private(i,idof,i1)
               do j = 1, ny
                 do i = 1, nx
                   do idof = 1, ndof
                     i1 = myZnaupd%ipntr(1) + &
                          (idof + ndof*(i-1) + ndof*nx*(j-1)) - 1
                     r(idof,i,j) = r(idof,i,j) - myZnaupd%workd(i1)
                   end do
                 end do
               end do
               call rhsBC3D( r, v, vm )

!.... Compute the norm of the residual (method 2)

               resnrm = zero
               !$omp parallel do private(i,idof) reduction(+: resnrm)
               do j = 1, ny
                 do i = 1, nx
                   do idof = 1, ndof
                     resnrm = resnrm + r(idof,i,j)*conjg(r(idof,i,j))
                   end do
                 end do
               end do
               resnrm = sqrt(resnrm)
               
               write(*,110) istep, resnrm
               write(80,110) istep, resnrm
               
               if ( resnrm .le. myZnaupd%tol * 0.1 ) exit

!.... Multiply through by the inverse block diagonal

               !$omp parallel do private(i), shared(ndof)
               do j = 1, ny
                 do i = 1, nx
                   call ZGETRS( 'N', ndof, 1, bmat(1,1,i,j), ndof, &
                                ipiv(1,i,j), r(1,i,j), ndof, info(i,j) )
                   if (info(i,j) .ne. 0) then
                     write(*,*) i, j, info(i,j)
                   end if
                 end do
               end do

!.... update solution using SOR
               
               sor = 1.0
               !$omp parallel do private(i,idof)
               do j = 1, ny
                 do i = 1, nx
                   do idof = 1, ndof
                     v(idof,i,j) = v(idof,i,j) - sor * r(idof,i,j)
                   end do
                 end do
               end do
               call itrbc3D(v,vm)

               if (debug.eq.2) then
                 j = 1
                 do i = 1, nx
                   write(20,50) x(i,j), &
                        real(v(1,i,j)), imag(v(1,i,j)), &
                        real(v(2,i,j)), imag(v(2,i,j)), &
                        real(v(3,i,j)), imag(v(3,i,j)), &
                        real(v(4,i,j)), imag(v(4,i,j)), &
                        real(v(5,i,j)), imag(v(5,i,j))
                 end do
                 close(20)

                 j = 1
                 do i = 1, nx
                   write(21,50) x(i,j), &
                        real(r(1,i,j)), imag(r(1,i,j)), &
                        real(r(2,i,j)), imag(r(2,i,j)), &
                        real(r(3,i,j)), imag(r(3,i,j)), &
                        real(r(4,i,j)), imag(r(4,i,j)), &
                        real(r(5,i,j)), imag(r(5,i,j))
                 end do
                 close(21)
                 if (istep.eq.nstep) call exit(1)
               end if

             end do

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  Approximate Factorization
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

             else     

!.... March in pseudo-time

             do istep = 1, nstep

!.... -------------------------> predictor phase <-------------------------

               !$omp parallel do private(i,idof)
               do j = 1, ny
                 do i = 1, nx
                   do idof = 1, ndof
                     vold(idof,i,j) = v(idof,i,j)
                   end do
                 end do
               end do

               time = time + Delt
               call itrbc3D(v,vm)
               
!.... ----------------------> multi-corrector phase <----------------------

               do iter = 1, niter

                 call lrhs3D( v, r, vm, x, y )
                 call ebe3D( r, v, vold, dtl, &
                             myZnaupd%workd(myZnaupd%ipntr(:)))
                 call rhsBC3D( r, v, vm )
                 
                 if (xper) then
                   call cpenta2p( nx, ny, ndof, matx, r, perx, per2x, 1 )
                 else
                   call cpenta2bc( nx, ny, ndof, matx, r, 1 )
                 end if
                 
                 if (yper) then
                   call cpenta1p( nx, ny, ndof, maty, r, pery, per2y, 1 )
                 else
                   call cpenta1bc( nx, ny, ndof, maty, r, 1 )
                 end if

!.... update the solution
                 
                 !$omp parallel do private(i,idof)
                 do j = 1, ny
                   do i = 1, nx
                     do idof = 1, ndof
                       v(idof,i,j) = v(idof,i,j) + r(idof,i,j)
                     end do
                   end do
                 end do
                 call itrbc3D(v,vm)

                 call lrhs3D( v, r, vm, x, y )
                 call ebe3D( r, v, vold, dtl, &
                             myZnaupd%workd(myZnaupd%ipntr(:)))
                 call rhsBC3D( r, v, vm )

!.... compute the norm of r

                 resnrm1 = zero
                 !$omp parallel do private(i,idof) reduction(+: resnrm1)
                 do j = 1, ny
                   do i = 1, nx
                     do idof = 1, ndof
                       resnrm1 = resnrm1 + &
                            r(idof,i,j)*conjg(r(idof,i,j))/dtl(i,j)**2
                     end do
                   end do
                 end do
!                resnrm1 = sqrt(resnrm1/real(nx*ny))
                 resnrm1 = sqrt(resnrm1)

                 if (debug.gt.0) then
!                  write(*,100) itr, iter, resnrm1
                   write(81,100) itr, iter, resnrm1
                 end if
100              format(i4,1x,i4,8(1x,1pe13.6))

!... Convergence criterion for inner iteration

                 if ( resnrm1 .le. 1.0e-8 ) exit

               end do  ! iter                 

!.... Compute the norm of the residual (method 1)

!!$               !$omp parallel do private(i,idof)
!!$               do j = 1, ny
!!$                 do i = 1, nx
!!$                   do idof = 1, ndof
!!$                     r(idof,i,j) = vold(idof,i,j) - v(idof,i,j)
!!$                   end do
!!$                 end do
!!$               end do
!!$
!!$               resnrm = zero
!!$               !$omp parallel do private(i,idof) reduction(+: resnrm)
!!$               do j = 1, ny
!!$                 do i = 1, nx
!!$                   do idof = 1, ndof
!!$                     resnrm = resnrm + &
!!$                          r(idof,i,j)*conjg(r(idof,i,j))/dtl(i,j)**2
!!$                   end do
!!$                 end do
!!$               end do
!!$               resnrm = sqrt(resnrm/real(nx*ny))

!.... Compute the norm of the residual (method 2)

               !$omp parallel do private(i,idof)
               do j = 1, ny
                 do i = 1, nx
                   do idof = 1, ndof
                     vold(idof,i,j) = v(idof,i,j)
                   end do
                 end do
               end do
               call lrhs3D( v, r, vm, x, y )
               call ebe3D( r, v, vold, dtl, &
                           myZnaupd%workd(myZnaupd%ipntr(:)))
               call rhsBC3D( r, v, vm )

               resnrm = zero; resmax = zero
!              !$omp parallel do private(i,idof) reduction(+: resnrm)
               do j = 1, ny
                 do i = 1, nx
                   do idof = 1, ndof
                     resnod = r(idof,i,j)*conjg(r(idof,i,j))/dtl(i,j)**2
                   end do
                   resnrm = resnrm + resnod
                   if (resnod .gt. resmax) then
                     resmax = resnod
                     ixrmax = i; iyrmax = j
                   end if
                  end do
               end do
               resnrm = sqrt(resnrm/real(nx*ny))
!              resnrm = sqrt(resnrm)
               resmax = sqrt(resmax)

               write(*,111) istep, delt, time, resnrm1, resnrm, &
                    abs((resnrm-resold)/resold), resmax, ixrmax, iyrmax
               write(80,110) istep, delt, time, resnrm1, resnrm, &
                    abs((resnrm-resold)/resold)
110            format(i4,8(1pe13.6,1x))
111            format(i4,6(1pe10.3,1x),'(',i3,1x,i3,')')

               if ( resnrm .le. myZnaupd%tol*0.1 ) exit

               fact = 1000.0
               ifact = one/fact

               if (gdoit) then
               if ( abs((resnrm-resold)/resold) .le. 5.0e-3) then
                 if (doit) then
                   delt = fact * delt
                   dtl =  fact * dtl
                   doit = .false.
                 else
                   doit = .true.
                   delt = ifact * delt
                   dtl =  ifact * dtl
                 end if
                 call lhs1f3D( matx, Ah, Dh, Dhi, Vh11, ABhi, dtl, dummy )
                 if (xper) then
                   call cpenta2p( nx, ny, ndof, matx, dummy, perx, per2x, 0)
                 else
                   call cpenta2bc( nx, ny, ndof, matx, dummy, 0 )
                 end if
                 
                 call lhs2f3D( maty, Bh, Dh, Dhi, Vh22, ABhi, dtl, dummy )
                 if (yper) then
                   call cpenta1p( nx, ny, ndof, maty, dummy, pery, per2y, 0 )
                 else
                   call cpenta1bc( nx, ny, ndof, maty, dummy, 0 )
                 end if
                 resold = 1.0e10
               else
                 resold = resnrm
               end if
             end if

             end do    ! istep

             if (debug.eq.3) then
             j = 1
             do i = 1, nx
               write(20,50) x(i,j), &
                    real(v(1,i,j)), imag(v(1,i,j)), &
                    real(v(2,i,j)), imag(v(2,i,j)), &
                    real(v(3,i,j)), imag(v(3,i,j)), &
                    real(v(4,i,j)), imag(v(4,i,j)), &
                    real(v(5,i,j)), imag(v(5,i,j))
             end do
             close(20)
50           format(12(e12.5,1x))

             do i = 1, myZnaupd%n
               i1 = myZnaupd%ipntr(1) + (i-1)
               i2 = myZnaupd%ipntr(2) + (i-1)
               r(i,1,1) = myZnaupd%workd(i1)
             end do
       
             j = 1
             do i = 1, nx
               write(21,50) x(i,j), &
                    real(r(1,i,j)), imag(r(1,i,j)), &
                    real(r(2,i,j)), imag(r(2,i,j)), &
                    real(r(3,i,j)), imag(r(3,i,j)), &
                    real(r(4,i,j)), imag(r(4,i,j)), &
                    real(r(5,i,j)), imag(r(5,i,j))
             end do
             close(21)

!.... write out restart file

             open(unit=10, file='soln.dat', form='unformatted', &
                  status='unknown')
             write(10) itr, time, nx, ny, nz, ndof, Re, Ma, Pr, gamma, cv
             write(10) v
             close(10)

!.... write out restart file

             open(unit=10, file='force.dat', form='unformatted', &
                  status='unknown')
             write(10) itr, time, nx, ny, nz, ndof, Re, Ma, Pr, gamma, cv
             write(10) r
             close(10)

             if (itr.eq.1) call exit(1)
             end if

             end if    ! iteration type

             do i = 1, myZnaupd%n
               i1 = myZnaupd%ipntr(1) + (i-1)
               i2 = myZnaupd%ipntr(2) + (i-1)
               myZnaupd%workd(i2) = v(i,1,1)
             end do

           else
             call error('eigdrv3d$','Illegal ARpack mode$')
           end if
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
! SSC 200704:  added pointer to fixe mismatched ranks
             vt(1:ndof,1:nx,1:ny) => myZnaupd%v(:,j)
             CALL lrhs3d(vt, r, vm, x, y)
             CALL eigBC3d(r, myZnaupd%v(1,j), vm )

             CALL zaxpy(myZnaupd%n, -myZneupd%d(j),   & 
                        myZnaupd%v(1,j), 1, r, 1)
             rd(j,1) = dble(myZneupd%d(j))
! SSC 200203:  Changed to avoid error
!            rd(j,2) = dimag(myZneupd%d(j))
             rd(j,2) = imagpart(myZneupd%d(j))
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
!!$              rl(idof,i,j) = -dtl(i,j) * (rl(idof,i,j) - fl(idof,i,j))
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
