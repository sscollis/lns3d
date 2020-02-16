

subroutine adi(v,workd,dtl)
 use global
 use local2d
 implicit none

    integer :: ier, i1, i2, mem=0
    complex :: dummy
    real :: resnrm, resnrm1, resold=1.0e10, fact, ifact, resmax, resnod
    integer :: ixrmax, iyrmax

    LOGICAL :: doit=.false.
    INTEGER :: debug=1
    INTEGER :: i, j, idof, itr

    complex :: v(ndof,nx,ny), vold(ndof,nx,ny), r(ndof,nx,ny)
    real :: r_r(ndof,nx,ny), r_im(ndof,nx,ny), workd(ndof,nx,ny)
    real :: Uc_r(ndof,nx,ny), Uc_im(ndof,nx,ny)
    real :: Uh_r(ndof,nx,ny), Uh_im(ndof,nx,ny)
    real :: rmatx(5,ndof,ndof,nx,ny), rmaty(5,ndof,ndof,nx,ny)
    real :: MiUc_r(ndof,nx,ny), MiUh_im(ndof,nx,ny)
    real :: perx(2,ndof,ndof,nx,ny), per2x(2,ndof,ndof,nx,ny) 
    real :: pery(2,ndof,ndof,nx,ny), per2y(2,ndof,ndof,nx,ny)
    real :: dtl(nx,ny)

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

!... First compute the residual

                 call lrhs3D( v, r, vm, x, y )
                 call ebe3D( r, v, vold, dtl, workd)
                 call rhsBC3D( r, v, vm )

!*****************************************************************!
!... New ADI technique derived by Dr. Collis and Mr. Bou-Rabee ...!
!*****************************************************************!
!
!    The derivation of this technique is available at:
!
!       http://www.owlnet.rice.edu/~q8gulf
!
!    Click on "Shift-And-Invert Applied."  In the next few lines
!    we will compute:
!
!     Uc_r  = rmaty^{-1} rmatx^{-1} r_r
!     Uc_im = r_im + omega * delt * Uc_r
!     Uh_im = rmaty^{-1} rmatx^{-1} Uc_im
!     Uh_r  = Uc_r - omega * delt * rmaty^{-1} rmatx^{-1} Uh_im
!
!    where r_im and r_r are the real and imaginary parts of the
!    residual respectively, and Uh_r and Uh_im are the real and
!    and imaginary parts of the update respectively.
!
!*****************************************************************!

!... Separate the real and imaginary parts of r (r_r and r_im)

            !$omp parallel do private(i,idof)
            do j = 1, ny
              do i = 1, nx
                do idof = 1, ndof
                    r_r(idof,i,j) = real(r(idof,i,j))
                  r_im(idof,i,j) = imag(r(idof,i,j))
                end do
              end do
            end do

!... Compute real part of Intermediate U (denoted Uc_r)

     !... Recall: Uc_r = rmaty^{-1} rmatx^{-1} r_r

     !... Call lhs1f to form the real part of matx (rmatx)

                 call lhs1f(rmatx, Ah, Dh, Vh11, dtl, v)

     !... Solve for rmatx^{-1} r_r

                 if (xper) then
                   call rpenta2p( nx, ny, ndof, rmatx, r_r, perx, per2x, 2 )
                 else
                   call rpenta2bc( nx, ny, ndof, rmatx, r_r, 2 )
                 end if

     !... Call lhs2f to form the real part of maty (rmaty)

                call lhs2f(rmaty, Bh, Dh, Vh22, dtl, v)

     !... Solve for rmaty^{-1} rmatx^{-1} r_r

                 if (yper) then
                   call rpenta1p( nx, ny, ndof, rmaty, r_r, pery, per2y, 2 )
                 else
                   call rpenta1bc( nx, ny, ndof, rmaty, r_r, 2 )
                 end if

!... Compute the imaginary part of Intermediate U (Uc_im)

     !... Recall: Uc_im = r_im + omega * delt * Uc_r
     !... Recall: Uc_im = r_im + Mi * Uc_r

                 !$omp parallel do private(i,idof)
                 do j = 1, ny
                   do i = 1, nx
                     do idof = 1, ndof
!               Uc_im(idof,i,j) = r_im(idof,i,j) + omega*dtl(i,j)*r_r(idof,i,j)
                        Uc_r(idof,i,j) = r_r(idof,i,j)
                     end do
                   end do
                 end do

                call MiMultBy(r_r, MiUc_r, dtl)

              !$omp parallel do private(i,idof)
               do j = 1, ny
                 do i = 1, nx
                   do idof = 1, ndof
                    Uc_im(idof,i,j) = r_im(idof,i,j) - MiUc_r(idof,i,j)
                   end do
                 end do
               end do

!... Compute the imaginary part of Uh, Uh_im

     !... Recall: Uh_im = rmaty^{-1} rmatx^{-1} Uc_im

     !... Compute rmaty^{-1} rmatx^{-1} Uc_im

                if (xper) then
                   call rpenta2p( nx, ny, ndof, rmatx, Uc_im, perx, per2x, 1 )
                 else
                   call rpenta2bc( nx, ny, ndof, rmatx, Uc_im, 1 )
                 end if

                 if (yper) then
                   call rpenta1p( nx, ny, ndof, rmaty, Uc_im, pery, per2y, 1 )
                 else
                   call rpenta1bc( nx, ny, ndof, rmaty, Uc_im, 1 )
                 end if

                !$omp parallel do private(i,idof)
                do j = 1, ny
                   do i = 1, nx
                     do idof = 1, ndof
                        Uh_im(idof,i,j) = Uc_im(idof,i,j)
                         r_im(idof,i,j) = Uh_im(idof,i,j)
                     end do
                   end do
                 end do

         call MiMultBy(Uh_im, MiUh_im, dtl)

     !... Compute rmaty^{-1} rmatx^{-1} Uh_im

                if (xper) then
                   call rpenta2p( nx, ny, ndof, rmatx, MiUh_im, perx, per2x, 1 )
                 else
                   call rpenta2bc( nx, ny, ndof, rmatx, MiUh_im, 1 )
                 end if

                 if (yper) then
                   call rpenta1p( nx, ny, ndof, rmaty, MiUh_im, pery, per2y, 1 )
                 else
                   call rpenta1bc( nx, ny, ndof, rmaty, MiUh_im, 1 )
                 end if

!... Compute the real part of Uh, Uh_r

     !... Recall: Uh_r = Uc_r - omega * delt * rmaty^{-1} rmatx^{-1} Uh_im

     !... Restore r by returning Uh_im and Uh_r into r

                 !$omp parallel do private(i,idof)
                 do j = 1, ny
                   do i = 1, nx
                     do idof = 1, ndof
                        r_r(idof,i,j) = Uc_r(idof,i,j) + MiUh_im(idof,i,j)
                          r(idof,i,j) = r_r(idof,i,j) + im * r_im(idof,i,j)
                          v(idof,i,j) = v(idof,i,j) + r(idof,i,j)
                     end do
                   end do
                 end do
                 call itrbc3D(v,vm)

!... compute the residual (again)

                 call lrhs3D( v, r, vm, x, y )
                 call ebe3D( r, v, vold, dtl, workd)
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
                if(debug.eq.1) then
                  write(*,100) itr, iter, resnrm1
                  write(81,100) itr, iter, resnrm1
                end if
100              format(i4,1x,i4,8(1x,1pe13.6))

!... Convergence criterion for inner iteration

!                 if ( resnrm1 .le. 1.0e-5 ) exit

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

               !$omp parallel do private(i,idof)
               do j = 1, ny
                 do i = 1, nx
                   do idof = 1, ndof
                     vold(idof,i,j) = v(idof,i,j)
                   end do
                 end do
               end do

!... compute the residual (for the third time)

               call lrhs3D( v, r, vm, x, y )
               call ebe3D( r, v, vold, dtl, workd)
               call rhsBC3D( r, v, vm )
!.... As you can probably tell, the residual is stored in r

!... Now compute the norm of the residual

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


!.... write some residual statistics to screen

               write(*,111) istep, delt, time, resnrm1, resnrm, &
                    abs((resnrm-resold)/resold), resmax, ixrmax, iyrmax

!.... also write the residual statistics to fort.80

               write(80,110) istep, delt, time, resnrm1, resnrm, &
                    abs((resnrm-resold)/resold)

110            format(i4,8(1pe13.6,1x))
111            format(i4,6(1pe10.3,1x),'(',i3,1x,i3,')')

               if ( resnrm .le. 1e-5 ) exit

               fact = 100.0
               ifact = one/fact

               if (.true.) then
               if ( abs((resnrm-resold)/resold) .le. 5.0e-6) then
                 if (.false.) then
                   delt = fact * delt
                   dtl =  fact * dtl
                   doit = .false.
!                 else
                   doit = .true.
                   delt = ifact * delt
                   dtl =  ifact * dtl
                 end if

                 resold = 1.0e10
               else
                 resold = resnrm
               end if
             end if

             end do    ! istep

             workd = v
return
end subroutine 
