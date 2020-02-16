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
        complex, allocatable :: matx(:,:,:,:,:), bcx(:,:,:,:)
        complex, allocatable :: perx(:,:,:,:,:), per2x(:,:,:,:,:)
        complex, allocatable :: maty(:,:,:,:,:), bcy(:,:,:,:)
        complex, allocatable :: pery(:,:,:,:,:), per2y(:,:,:,:,:)
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
        allocate (matx(5,ndof,ndof,nx,ny), maty(5,ndof,ndof,nx,ny), STAT=ier)
        if (ier .ne. 0) call error(code,'Insufficient Memory for mat$')

!.... allocate memory for the per

        if (xper .or. yper) then
          mem = mem + 2*2*ny*nx*ndof*ndof*2
          allocate (perx(2,ndof,ndof,nx,ny), per2x(2,ndof,ndof,nx,ny), &
                    pery(2,ndof,ndof,nx,ny), per2y(2,ndof,ndof,nx,ny), &
                    STAT=ier)
          if (ier .ne. 0) call error(code, 'Insufficient Memory for per$')
        end if
        
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
            
            if ( (iter.eq.1 .and. istep.eq.1) .or. &
                 (iter.eq.1 .and. lstep.eq.1 .and. impl.eq.2) ) then

              call lhs1f3D( matx, Ah, Dh, Dhi, Vh11, ABhi, dtl, v )
  
              if (xper) then
                call cpenta2p( nx, ny, ndof, matx, r, perx, per2x, 0)
              else
                call cpenta2bc( nx, ny, ndof, matx, r, 0 )
              end if
  
            end if
                    
            if (xper) then
              call cpenta2p( nx, ny, ndof, matx, r, perx, per2x, 1)
            else
              call cpenta2bc( nx, ny, ndof, matx, r, 1 )
            end if

!.... form the LHS for the eta direction

            if ((iter.eq.1 .and. istep.eq.1) .or. &
                 (iter.eq.1 .and. lstep.eq.1 .and. impl.eq.2)) then

              call lhs2f3D( maty, Bh, Dh, Dhi, Vh22, ABhi, dtl, v )

              if (yper) then
                call cpenta1p( nx, ny, ndof, maty, r, pery, per2y, 0 )
              else
                call cpenta1bc( nx, ny, ndof, maty, r, 0 )
              end if

            end if

            if (yper) then
              call cpenta1p( nx, ny, ndof, maty, r, pery, per2y, 1)
            else
              call cpenta1bc( nx, ny, ndof, maty, r, 1 )
            end if

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
