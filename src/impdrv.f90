!=============================================================================!
        subroutine impdrv()
!  
!  Driver for implicit time advancement of the compressible Navier-Stokes
!  solver. 
!
!  This is a new version which uses either the 2nd or 4th order LHS
!
!  Revised: 7-15-96
!           12-15-00  Changed to IJ indices for Cached based computers [SSC]
!=============================================================================!
        use global
        use local2d
        implicit none
!=============================================================================!
        real, allocatable :: v(:,:,:), vold(:,:,:), v2old(:,:,:), r(:,:,:)
        real, allocatable :: mat(:,:,:,:,:), bc(:,:,:,:)
        real, allocatable :: per(:,:,:,:,:), per2(:,:,:,:)
        real, allocatable :: dtl(:,:)
        
        real :: rtime, cpu
        integer :: i, j, idof, jdof
        integer :: ier, mem = 0
        integer, external :: igetver
        real*4, external :: second
        character(80) :: name, code='ImpDrv$'
!=============================================================================!

!.... allocate the storage area for the disturbance field

        mem = mem + ny*nx*ndof
        allocate (v(ndof,nx,ny), vold(ndof,nx,ny), STAT=ier)
        if (ier .ne. 0) call error(code,'Insufficient Memory for v$')
        !$omp parallel do private(i,idof)
        do j = 1, ny
          do i = 1, nx
            do idof = 1, ndof
              v(idof,i,j) = zero; vold(idof,i,j) = zero
            end do
          end do
        end do

!.... generate the initial condition

        call genini(v, vold)

!.... allocate memory for the RHS

        mem = mem + ny*nx*ndof
        allocate (r(ndof,nx,ny), STAT=ier)
        if (ier .ne. 0) call error(code,'Insufficient Memory for r$')
        !$omp parallel do private(i,idof)
        do j = 1, ny
          do i = 1, nx
            do idof = 1, ndof
              r(idof,i,j) = zero
            end do
          end do
        end do

!.... allocate memory for the previous, previous disturbance field

        if (impl.eq.2) then
          mem = mem + ny*nx*ndof
          allocate (v2old(ndof,nx,ny), STAT=ier)
          if (ier .ne. 0) then
            call error(code, 'Insufficient Memory for v2old$')
          endif
          !$omp parallel do private(i,idof)
          do j = 1, ny
            do i = 1, nx
              do idof = 1, ndof
                v2old(idof,i,j) = zero
              end do
            end do
          end do
        end if

!.... allocate memory for the LHS

        mem = mem + ny*nx*ndof*5
        allocate (mat(5,ndof,ndof,nx,ny), STAT=ier)
        if (ier .ne. 0) then
          call error(code,'Insufficient Memory for mat$')
        endif
        !$omp parallel do private(i,idof,jdof)
        do j = 1, ny
          do i = 1, nx
            do jdof = 1, ndof
              do idof = 1, ndof
                mat(:,idof,jdof,i,j) = zero
              end do
            end do
          end do
        end do

!.... allocate memory for the per

        if (xper .or. yper) then
          mem = mem + ny*nx*ndof*ndof*2 + max(nx,ny)*ndof*ndof*6
          allocate (per(2,ndof,ndof,nx,ny), &
                    per2(6,ndof,ndof,max(nx,ny)), STAT=ier)
          if (ier .ne. 0) call error(code,'Insufficient Memory for per$')
          !$omp parallel do private(i,idof,jdof)
          do j = 1, ny
            do i = 1, nx
              do jdof = 1, ndof
                do idof = 1, ndof
                  per(:,idof,jdof,i,j) = zero
                end do
              end do
            end do
            do jdof = 1, ndof
              do idof = 1, ndof
                per2(:,idof,jdof,j) = zero
              end do
            end do
          end do
        end if
        
!.... allocate memory for the BC (for second order LHS)

        mem = mem + max(nx,ny)*ndof*ndof*14
        allocate (bc(14,ndof,ndof,max(nx,ny)), STAT=ier)
        if (ier .ne. 0) call error(code,'Insufficient Memory for bc$')

!.... allocate local time step

        mem = mem + ny*nx
        allocate (dtl(nx,ny), STAT=ier)
        if (ier .ne. 0) call error(code,'Insufficient Memory for dtl$')
        !$omp parallel do private(i)
        do j = 1, ny
          do i = 1, nx
            dtl(i,j) = zero
          end do
        end do

!       write(*,"(' ImpDrv allocated ===> ',1pe13.6,' words')") float(mem)
!=============================================================================!

!.... issue a warning if Ny is even

#ifdef CRAY
        if (mod(ny,2).eq.0) then
          write(*,"('WARNING: Ny is even ==> bank conflicts on Cray')")
        end if
#endif

!.... open the history file

        iver = igetver( base, 'h'//char(0) ) + 1
        call addver(base, 'h'//char(0), iver, filen, lfile)
        open(unit=ihist, file=filen, form='formatted', status='unknown')

!.... satisfy boundary conditions on v

        istep = 0
        call itrbc(v,vm)

!.... compute the CFL

        if (linear.eq.1) then
          call set_time_step(vm, dtl)
        else
          call set_time_step(v, dtl)
        end if

!.... loop through the time steps

        do istep = 1, nstep
          
!          cpu = second()

          if (linear.eq.0) call set_time_step(v, dtl)
          
!.... -------------------------> predictor phase <-------------------------

!.... quadratic prediction (makes things go unstable?)

!         if (impl.eq.2 .and. istep .gt. 2) then
!           r = 3.0 * v - 3.0 * vold + v2old
!         end if

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

!.... quick test

!!$          !$doacross local(i,idof)
!!$          !$omp parallel do private(i,idof)
!!$          do j = 1, ny
!!$            do i = 1, nx
!!$              do idof = 1, ndof
!!$                v(idof,i,j) = zero
!!$              end do
!!$            end do
!!$          end do
          
!!$          if (impl.eq.2 .and. istep .gt. 2) then
!!$            !$doacross local(i,idof)
!!$            !$omp parallel do private(i,idof)
!!$            do j = 1, ny
!!$              do i = 1, nx
!!$                do idof = 1, ndof
!!$                  v(idof,i,j) = r(idof,i,j)
!!$                end do
!!$              end do
!!$            end do
!!$          end if

          time  = time + Delt
          call itrbc(v,vm)
          
!         write(*,*) 'Predictor ', second()-cpu

!.... ----------------------> multi-corrector phase <----------------------

          do iter = 1, niter+1
          
            if (iter.eq.niter+1 .and. ires.eq.0) goto 200

!           cpu = second()
            
!.... Backward Euler

            if (impl.eq.1 .or. (lstep .eq. 0 .and. impl.eq.2 )) then
              if (linear.eq.0) then
                call genmtrx( v, zero, r)
              else
                call lrhs( v, r, vm, x, y)
              end if
              call rbe( r, v, vold, dtl )
              alfa = one
            end if

!.... Two-step backward

            if ( lstep .gt. 0 .and. impl.eq.2 ) then
              if (linear.eq.0) then
                call genmtrx( v, zero, r )
              else
                call lrhs(v, r, vm, x, y)
              end if
              call rtb( r, v, vold, v2old, dtl)
              alfa = pt66
            end if

!.... Midpoint rule

            if ( impl.eq.3 ) then
              if (linear.eq.0) then
                call genmtrx( pt5*(v+vold), zero, r )
              else
                call lrhs( pt5*(v+vold), r, vm, x, y)
              end if
              call rbe( r, v, vold, dtl )
              alfa = pt5
            end if
            
!           write(*,*) 'Build RHS ', second()-cpu

!.... satisfy BC's on the RHS

!           cpu = second()
            if (iter.le.niter+1) call rhsBC(r, v, vm)
!           write(*,*) 'Build rhsBC ', second()-cpu

!.... output RHS statistics

!           cpu = second()
            call resstat(r,dtl)
!           write(*,*) 'Resstat ', second()-cpu

            if (iter.eq.niter+1) then        ! when you want to know |res|
              if (ires.eq.2 .and. istep.eq.nstep) then
                name = 'resid.dat'
                if(impl.eq.1) then
                  call wdata(name,v-vold,nx,ny,nz,ndof,Ma,zero,Re,time)
                else
                  call wdata(name,v-vold,nx,ny,nz,ndof,Ma,zero,Re,time)
                end if
              end if
              goto 200
            end if

            if (iLHS.eq.0) then

!==============================================================================
!       S O L V E   S E C O N D   O R D E R   L H S 
!==============================================================================

!.... form the LHS for the xi direction
            
!           cpu = second()
            call lhs1( mat, Ah, Dh, Vh11, bc, dtl, v )
!           write(*,*) 'Build LHS1 ', second()-cpu

!           cpu = second()
            if (xper) then
              call solve2p( nx, ny, ndof, mat, r, per, per2 )
            else
              call solve2bc( nx, ny, ndof, mat, r, bc )
            end if
!           write(*,*) 'Solve LHS1 ', second()-cpu

!.... form the LHS for the eta direction

!           cpu = second()
            call lhs2( mat, Bh, Dh, Vh22, bc, dtl, v )
!           write(*,*) 'Build LHS2 ', second()-cpu
          
!           cpu = second()
            if (yper) then
              call solve1p( nx, ny, ndof, mat, r, per, per2 )
            else
              call solve1bc( nx, ny, ndof, mat, r , bc )
            end if
!           write(*,*) 'Solve LHS2 ', second()-cpu

            else
!==============================================================================
!       S O L V E   F O U R T H   O R D E R   L H S 
!==============================================================================

!.... form the LHS for the xi direction
            
!           cpu = second()
            call lhs1f( mat, Ah, Dh, Vh11, dtl, v )
!           write(*,*) 'Build LHS1 ', second()-cpu

!           cpu = second()
            if (xper) then
              call penta2p( nx, ny, ndof, mat, r, per, per2 )
            else
              call penta2bc( nx, ny, ndof, mat, r )
            end if
!           write(*,*) 'Solve LHS1 ', second()-cpu

!.... form the LHS for the eta direction

!           cpu = second()
            call lhs2f( mat, Bh, Dh, Vh22, dtl, v )
!           write(*,*) 'Build LHS2 ', second()-cpu
          
!           cpu = second()
            if (yper) then
              call penta1p( nx, ny, ndof, mat, r, per, per2 )
            else
              call penta1bc( nx, ny, ndof, mat, r )
            end if
!           write(*,*) 'Solve LHS2 ', second()-cpu
            
            end if       ! iLHS
!==============================================================================

!.... Inviscid: rotate back to global coordinates on the body

            if ((.not. Navier) .and. (.not. yper)) then
              call rbcfix( nx, ny, ndof, bnb, r )
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
            
!.... satisfy boundary conditions on v

            call itrbc(v,vm)
            
          end do        ! loop on iter

 200      continue
          lstep = lstep + 1

!.... quit if less than 10 seconds remain

#ifdef CRAY
          call tremain(rtime)
          if (rtime .le. 10) goto 5000
#endif

!.... output time traces

          if (tflag) call traces(v)
          
!.... compute temporal statistics

          if(fflag .and. xper .and. istep .ne. nstep .and. &
            mod(istep-1,itout) .eq. 0 ) then
            call fstat2d(.false., v, vm)        ! don't print profiles
          end if

!.... compute spatial statistics

          if (sflag .and. linear.eq.1) then
            call sstat(v, vm, x, y)
          end if

!.... write out restart and q files 

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
          call fstat2d(.true., v, vm)   ! print profiles
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

        end subroutine impdrv

!=============================================================================
        subroutine rbe( rl, vl, voldl, dtl )
!
!  Backward Euler RHS:    r = vold - v - delt * r
!
!  This also works for midpoint rule.
!=============================================================================
        use global
        implicit none
        
        real    :: rl(ndof,nx,ny), vl(ndof,nx,ny), voldl(ndof,nx,ny)
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
        end subroutine rbe

!=============================================================================
        subroutine rtb( rl, vl, voldl, v2oldl, dtl )
!
!  Two-step backward RHS:  r = onept33 * vold - pt33 *
!                              v2old - v - pt66 * delt * r
!
!=============================================================================
        use global
        implicit none
        
        real    :: rl(ndof,nx,ny), vl(ndof,nx,ny), voldl(ndof,nx,ny)
        real    :: v2oldl(ndof,nx,ny)
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
        end subroutine rtb
