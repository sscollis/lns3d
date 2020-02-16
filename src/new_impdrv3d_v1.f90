!=============================================================================!
        subroutine ImpDrv3D()
!  
!  Driver for implicit time advancement of the 3-d compressible 
!  Navier-Stokes solver. 
!
!  Modified to use LU factorization of LHS for a 5x speedup
!
!  Modified to filter the field in y if desired
!
!   7-24-98:  Removed SSD
!  12-18-00:  Changed IJ indices [SSC] 
!=============================================================================!
        use global
        use local2d
        use local3d
        implicit none
!=============================================================================!
        complex, allocatable :: v(:,:,:), vold(:,:,:), v2old(:,:,:), r(:,:,:)

        real, allocatable :: r_r(:,:,:), r_im(:,:,:)
        real, allocatable :: Uc_r(:,:,:), Uc_im(:,:,:)
        real, allocatable :: Uh_r(:,:,:), Uh_im(:,:,:)
        real, allocatable :: rmatx(:,:,:,:,:), bcx(:,:,:,:)
        real, allocatable :: perx(:,:,:,:,:), per2x(:,:,:,:,:)
        real, allocatable :: rmaty(:,:,:,:,:), bcy(:,:,:,:)
        real, allocatable :: pery(:,:,:,:,:), per2y(:,:,:,:,:)
        real, allocatable    :: dtl(:,:)
        
        real    :: rtime, cpu
        integer :: i, j, idof, jdof
        integer :: ier, mem = 0
        integer, external :: igetver
        real, external :: second
        character*80 :: name, code='ImpDrv3D$'
        
        logical :: lfilter = .false.
!=============================================================================!

!.... allocate local time step

        mem = mem + ny*nx
        allocate (dtl(nx,ny), STAT=ier)
        if (ier .ne. 0) call error(code,'Insufficient Memory for dtl$')

!.... compute the CFL

        call set_time_step(vm, dtl)

!.... allocate the storage area for the disturbance field

        mem = mem + 2*ny*nx*ndof
        allocate (v(ndof,nx,ny), vold(ndof,nx,ny), STAT=ier)
        if (ier .ne. 0) call error(code,'Insufficient Memory for v$')

!.... generate the initial condition

        call genini3d(v, vold)

!.... allocate memory for the RHS

        mem = mem + 2*ny*nx*ndof
        allocate (r(ndof,nx,ny), STAT=ier)
        if (ier .ne. 0) call error(code,'Insufficient Memory for r$')

!.... allocate memory for the previous, previous disturbance field

        if (impl.eq.2) then
          mem = mem + 2*ny*nx*ndof
          allocate (v2old(ndof,nx,ny), STAT=ier)
          if (ier .ne. 0) call error(code, 'Insufficient Memory for v2old$')
        end if

!.... allocate memory for the LHS

        mem = mem + ny*nx*ndof*ndof*5
        allocate (rmatx(5,ndof,ndof,nx,ny), rmaty(5,ndof,ndof,nx,ny), STAT=ier)
        if (ier .ne. 0) call error(code,'Insufficient Memory for mat$')

!.... allocate memory for the per

        if (xper .or. yper) then
          mem = mem + 2*2*ny*nx*ndof*ndof*2
          allocate (perx(2,ndof,ndof,nx,ny), per2x(2,ndof,ndof,nx,ny), &
                    pery(2,ndof,ndof,nx,ny), per2y(2,ndof,ndof,nx,ny), &
                    STAT=ier)
          if (ier .ne. 0) call error(code, 'Insufficient Memory for per$')
        end if
        
        allocate (r_r(ndof,nx,ny), r_im(ndof,nx,ny),   &
                  Uc_r(ndof,nx,ny), Uc_im(ndof,nx,ny), &
                  Uh_r(ndof,nx,ny), Uh_im(ndof,nx,ny), STAT=ier)
        if (ier .ne. 0) call error(code,'Insufficient Memory for matrices$')

!       write(*,"(' ImpDrv3D allocated ===> ',1pe13.6,' words')") float(mem)
!=============================================================================!

!.... issue a warning if Ny is even

#ifdef CRAY
        if (mod(ny,2).eq.0) then
          write(*,"('W A R N I N G:  Ny is even ==> bank conflicts on Cray')")
        end if
#endif

!.... issue a warning if filter is true

        if (lfilter) then
          write(*,"('W A R N I N G:  Filter is ON')")
        end if

!.... open the history file

        iver = igetver( base, 'h'//char(0) ) + 1
        call addver(base, 'h'//char(0), iver, filen, lfile)
        open(unit=ihist, file=filen, form='formatted', status='unknown')

!.... satisfy boundary conditions on v

        istep = 0
        call itrbc3D(v,vm)

!.... loop through the time steps

        do istep = 1, nstep
          
!.... -------------------------> predictor phase <-------------------------

          if (impl.eq.2 .and. lstep .gt. 0) then
            !$doacross local(i,idof)
            !$omp parallel do private(i,idof)
            do j = 1, ny
              do i = 1, nx
                do idof = 1, ndof
                  v2old(idof,i,j) = vold(idof,i,j)
                end do
              end do
            end do
          end if

          !$doacross local(i,idof)
          !$omp parallel do private(i,idof)
          do j = 1, ny
            do i = 1, nx
              do idof = 1, ndof
                vold(idof,i,j) = v(idof,i,j)
              end do
            end do
          end do

!.... try starting from zero every time

!!$          !$doacross local(i,idof)
!!$          !$omp parallel do private(i,idof)
!!$          do j = 1, ny
!!$            do i = 1, nx
!!$              do idof = 1, ndof
!!$                v(idof,i,j) = zero
!!$              end do
!!$            end do
!!$          end do

          time  = time + Delt
          call itrbc3D(v,vm)
          
!.... ----------------------> multi-corrector phase <----------------------

          do iter = 1, niter+1
          
            if (iter.eq.niter+1 .and. ires.eq.0) goto 200

!.... Backward Euler

            if (impl.eq.1 .or. (lstep .eq. 0 .and. impl.eq.2 )) then
              call lrhs3D(v, r, vm, x, y)
              call rbe3D( r, v, vold, dtl )
              alfa = one
            end if

!.... Two-step backward

            if ( lstep .gt. 0 .and. impl.eq.2 ) then
              call lrhs3D(v, r, vm, x, y)
              call rtb3D( r, v, vold, v2old, dtl)
              alfa = pt66
            end if

!.... Midpoint rule

            if ( impl.eq.3 ) then
              call lrhs3D( pt5*(v+vold), r, vm, x, y)
              call rbe( r, v, vold, dtl )
              alfa = pt5
            end if
            
!.... satisfy BC's on the RHS

            if (iter.le.niter+1) call rhsBC3D(r, v, vm)

!.... output RHS statistics

            call resstat3D(r,dtl)
            
            if (iter.eq.niter+1) then        ! when you want to know |res|
              if (ires.eq.2 .and. istep.eq.nstep) then
                name = 'resid.dat'
                call wdata(name,abs(v-vold),nx,ny,nz,ndof,Ma,zero,Re,time)
              end if
              if (ires.eq.3 .and. istep.eq.nstep) then
                name = 'resid.dat'
                call wdata(name,abs(r),nx,ny,nz,ndof,Ma,zero,Re,time)
              end if
              goto 200
            end if

!==============================================================================
!       S O L V E   F O U R T H   O R D E R   L H S 
!==============================================================================
            if (iLHS.eq.0) call error(code,'iLHS = 0 not supported$')

!.... form the LHS for the xi direction
            
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

              !$omp parallel do private(i,idof)
              do j = 1, ny
                do i = 1, nx
                  do idof = 1, ndof
                    Uc_im(idof,i,j) = r_im(idof,i,j) + &
                                      omega*dtl(i,j)*r_r(idof,i,j)
                    Uc_r(idof,i,j)  = r_r(idof,i,j)
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

     !... Compute rmaty^{-1} rmatx^{-1} Uh_im

              if (xper) then
                call rpenta2p( nx, ny, ndof, rmatx, Uh_im, perx, per2x, 1 )
              else
                call rpenta2bc( nx, ny, ndof, rmatx, Uh_im, 1 )
              end if
              
              if (yper) then
                call rpenta1p( nx, ny, ndof, rmaty, Uh_im, pery, per2y, 1 )
              else
                call rpenta1bc( nx, ny, ndof, rmaty, Uh_im, 1 )
              end if

!... Compute the real part of Uh, Uh_r 

     !... Recall: Uh_r = Uc_r - omega * delt * rmaty^{-1} rmatx^{-1} Uh_im

     !... Restore r by returning Uh_im and Uh_r into r

              !$omp parallel do private(i,idof)
              do j = 1, ny
                do i = 1, nx
                  do idof = 1, ndof
                    r_r(idof,i,j) = Uc_r(idof,i,j) - &
                                    omega * dtl(i,j) * Uh_im(idof,i,j)
                    r(idof,i,j) = r_r(idof,i,j) + im * r_im(idof,i,j)
                  end do
                end do
              end do

!==============================================================================

!.... Inviscid: rotate back to global coordinates on the body

            if ((.not. Navier) .and. (.not. yper)) then
              call bcfix( nx, ny, ndof, bnb, r )
            end if
              
!.... update v

            !$doacross local(i,idof)
            !$omp parallel do private(i,idof)
            do j = 1, ny
              do i = 1, nx
                do idof = 1, ndof
                  v(idof,i,j) = v(idof,i,j) + r(idof,i,j)
                end do
              end do
            end do

!.... filter if desired

            if (lfilter) call cfilter( ny, nx, ndof, v )

!.... satisfy boundary conditions on v

            call itrbc3D(v,vm)
            
          end do        ! loop on iter

 200      continue
          lstep = lstep + 1

!.... quit if less than 10 seconds remain

#ifdef CRAY
          call tremain(rtime)
          if (rtime .le. 10) goto 5000
#endif

!.... output time traces

!         if (tflag) call traces3D(v)
          
!.... compute temporal statistics

          if(fflag .and. xper .and. istep .ne. nstep .and. &
            mod(istep-1,itout) .eq. 0 ) then
            call fstat3D(.false., v, vm)        ! don't print profiles
          end if

!.... compute spatial statistics

!         if (sflag .and. linear.eq.1) then
!           call sstat3D(v, vm, x, y)
!         end if

!.... write out restart file

          if ( mod(lstep,ntout) .eq. 0 .and. istep .ne. nstep) then
            
            iver = igetver( base, 'R'//char(0) ) + 1
            call addver(base, 'R'//char(0), iver, filen, lfile)
            open(unit=ifile, file=filen, form='unformatted', &
                  status='unknown')
            write(ifile) lstep, time, nx, ny, nz, ndof, &
                         Re, Ma, Pr, gamma, cv
            write(ifile) v
            if (impl.eq.2) write(ifile) vold
            close(ifile)

            write(*,10) filen(1:index(filen,char(0))-1), lstep, time
            write(ihist,10) filen(1:index(filen,char(0))-1), lstep, time

          end if
          
        end do          ! loop on istep

!.... ---------------------->  Post Processing  <----------------------

5000    continue

!.... compute final statistics

        if (fflag .and. xper) then
          call fstat3D(.true., v, vm)           ! print profiles
        end if

!.... write the restart file

        iver = igetver( base, 'R'//char(0) ) + 1
        call addver(base, 'R'//char(0), iver, filen, lfile)
        open(unit=ifile, file=filen, form='unformatted', &
             status='unknown')
        write(ifile) lstep, time, nx, ny, nz, ndof, &
                     Re, Ma, Pr, gamma, cv
        write(ifile) v
        if (impl.eq.2) write(ifile) vold
        close(ifile)

        write(*,10) filen(1:index(filen,char(0))-1), lstep, time
        write(ihist,10) filen(1:index(filen,char(0))-1), lstep, time

!.... close the history file

        close(ihist)

        return
        
  10    format(' Wrote ',a,' at step = ',i7,' time = ',1pe13.6)

        end

!=============================================================================
        subroutine rbe3D( rl, vl, voldl, dtl )
!
!  Backward Euler RHS:    r = vold - v - delt * r
!
!  This also works for midpoint rule.
!=============================================================================
        use global
        implicit none
        
        complex :: rl(ndof,nx,ny), vl(ndof,nx,ny), voldl(ndof,nx,ny)
        real    :: dtl(nx,ny)
        integer :: i,j,idof
!=============================================================================

        !$doacross local(i,idof)
        !$omp parallel do private(i,idof)
        do j = 1, ny
          do i = 1, nx
            do idof = 1, ndof
              rl(idof,i,j) = voldl(idof,i,j) - vl(idof,i,j) - &
                             dtl(i,j) * rl(idof,i,j)
            end do
          end do
        end do
        
        return
        end

!=============================================================================
        subroutine rtb3D( rl, vl, voldl, v2oldl, dtl )
!
!  Two-step backward RHS:  r = onept33 * vold - pt33 *
!                              v2old - v - pt66 * delt * r
!
!=============================================================================
        use global
        implicit none
        
        complex :: rl(ndof,nx,ny), vl(ndof,nx,ny), voldl(ndof,nx,ny)
        complex :: v2oldl(ndof,nx,ny)
        real    :: dtl(nx,ny)
        integer :: i,j,idof
!=============================================================================

        !$doacross local(i,idof)
        !$omp parallel do private(i,idof)
        do j = 1, ny
          do i = 1, nx
            do idof = 1, ndof
              rl(idof,i,j) = onept33 * voldl(idof,i,j) - &
                             pt33 * v2oldl(idof,i,j) - vl(idof,i,j) - &
                             pt66 * dtl(i,j) * rl(idof,i,j)
            end do
          end do
        end do
        
        return
        end
