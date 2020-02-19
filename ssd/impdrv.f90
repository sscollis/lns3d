!=============================================================================!
        subroutine impdrv()
!  
!  Driver for implicit time advancement of the compressible Navier-Stokes
!  solver. 
!
!  This is a new version which uses either the 2nd or 4th order LHS
!
!  Revised: 7-15-96
!
!=============================================================================!
        use stuff
        use local
        implicit none
!=============================================================================!
        real, allocatable :: v(:), vold(:), v2old(:), r(:)
        real, allocatable :: mat(:), bc(:), per(:), per2(:)
        real, allocatable :: mult(:), fact(:)
        real, allocatable :: dtl(:)
#ifdef __GFORTRAN__
  external :: fstat
#endif  
        real :: rtime, cpu
        integer :: i, j, idof
        integer :: ier, mem = 0
        integer, external :: igetver
        real, external :: second
        character*80 :: name, code='ImpDrv$'
!=============================================================================!

!.... allocate the storage area for the disturbance field

        mem = mem + ny*nx*ndof
        allocate (v(ny*nx*ndof), vold(ny*nx*ndof), STAT=ier)
        if (ier .ne. 0) call error(code,'Insufficient Memory for v$')

!.... generate the initial condition

        call genini(v, vold)

!.... allocate memory for the RHS

        mem = mem + ny*nx*ndof
        allocate (r(ny*nx*ndof), STAT=ier)
        if (ier .ne. 0) call error(code,'Insufficient Memory for r$')

!.... allocate memory for the previous, previous disturbance field

        if (impl.eq.2) then
          mem = mem + ny*nx*ndof
          allocate (v2old(ny*nx*ndof), STAT=ier)
          if (ier .ne. 0) call error(code, 'Insufficient Memory for v2old$')
        end if

!.... allocate memory for the LHS

        mem = mem + ny*nx*ndof*5
        allocate (mat(ny*nx*ndof*ndof*5), STAT=ier)
        if (ier .ne. 0) call error(code,'Insufficient Memory for mat$')

!.... allocate memory for the per

        if (xper .or. yper) then
          mem = mem + ny*nx*ndof*ndof*2 + max(nx,ny)*ndof*ndof*6
          allocate (per(ny*nx*ndof*ndof*2), &
                    per2(max(nx,ny)*ndof*ndof*6), STAT=ier)
          if (ier .ne. 0) call error(code,'Insufficient Memory for per$')
        end if
        
!.... allocate memory for the BC (for second order LHS)

        mem = mem + max(nx,ny)*ndof*ndof*14
        allocate (bc(max(nx,ny)*ndof*ndof*14), STAT=ier)
        if (ier .ne. 0) call error(code,'Insufficient Memory for bc$')

!.... allocate work space for implicit solver

        mem = mem + max(ny,nx)*2
        allocate (mult(max(ny,nx)), fact(max(ny,nx)), STAT=ier)
        if (ier .ne. 0) call error(code,'Insufficient Memory for work$')
        
!.... allocate local time step

        mem = mem + ny*nx
        allocate (dtl(ny*nx), STAT=ier)
        if (ier .ne. 0) call error(code,'Insufficient Memory for dtl$')

        write(*,"(' ImpDrv allocated ===> ',1pe13.6,' words')") float(mem)
!=============================================================================!

!.... issue a warning if Ny is even

        if (mod(ny,2).eq.0) then
          write(*,"('WARNING: Ny is even ==> bank conflicts')")
        end if

!.... open the history file

        iver = igetver( base, 'h'//char(0) ) + 1
        call addver(base, 'h'//char(0), iver, filen, lfile)
        open(unit=ihist, file=filen, form='formatted', status='unknown')

!.... satisfy boundary conditions on v

        istep = 0
        call itrbc(v,vm)

!.... compute the CFL

        if (linear.eq.1) then
          call dtcfl(vm, dtl)
          call mlocal('deallocate')
        end if

!.... loop through the time steps

        do istep = 1, nstep
          
          if (linear.eq.0) call dtcfl(v, dtl)
          
!.... -------------------------> predictor phase <-------------------------

!.... quadratic prediction (makes things go unstable?)

!         if (impl.eq.2 .and. istep .gt. 2) then
!           r = 3.0 * v - 3.0 * vold + v2old
!         end if

          if (impl.eq.2 .and. lstep .gt. 0) v2old = vold

          vold  = v

!         if (impl.eq.2 .and. istep .gt. 2) v = r 

          time  = time + Delt
          call itrbc(v,vm)
          
!.... ----------------------> multi-corrector phase <----------------------

          do iter = 1, niter+1
          
            if (iter.eq.niter+1 .and. ires.eq.0) goto 200

!           cpu = second()
            
!.... Backward Euler

            if (impl.eq.1 .or. (lstep .eq. 0 .and. impl.eq.2 )) then
              if (linear.eq.0) then
                call genmtrx( v, zero, .false. )
                call rhs( r, v )
              else
                call lrhs( v, r, vm, x, y)
              end if
              call rbe( r, v, vold, dtl )
              alfa = one
            end if

!.... Two-step backward

            if ( lstep .gt. 0 .and. impl.eq.2 ) then
              if (linear.eq.0) then
                call genmtrx( v, zero, .false. )
                call rhs( r, v )
              else
                call lrhs(v, r, vm, x, y)
              end if
              call rtb( r, v, vold, v2old, dtl)
              alfa = pt66
            end if

!.... Midpoint rule

            if ( impl.eq.3 ) then
              if (linear.eq.0) then
                call genmtrx( pt5*(v+vold), zero, .false. )
                call rhs( r, pt5*(v+vold) )
              else
                call lrhs( pt5*(v+vold), r, vm, x, y)
              end if
              call rbe( r, v, vold, dtl )
              alfa = pt5
            end if
            
!.... satisfy BC's on the RHS

            if (iter.le.niter+1) call rhsBC(r, v, vm)

!.... output RHS statistics

            call resstat(r,dtl)
            
!           write(*,*) 'Build RHS ', second()-cpu

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
            call lhs1( mat, Q1, Q2, Vh1, bc, dtl, v )
!           write(*,*) 'Build LHS1 ', second()-cpu

!           cpu = second()
            if (xper) then
              call solve1p( ny, nx, ndof, mat, r, per, per2, mult, fact )
            else
              call solve1bc( ny, nx, ndof, mat, r, bc, mult, fact )
!             call solve1( ny, nx, ndof, mat, r, mult, fact )
            end if
!           write(*,*) 'Solve LHS1 ', second()-cpu

!.... form the LHS for the eta direction

!           cpu = second()
            call lhs2( mat, Q1, Q2, Vh1, bc, dtl, v )
!           write(*,*) 'Build LHS2 ', second()-cpu
          
!           cpu = second()
            call solve2bc( ny, nx, ndof, mat, r , bc, mult, fact)
!           write(*,*) 'Solve LHS2 ', second()-cpu

            else
!==============================================================================
!       S O L V E   F O U R T H   O R D E R   L H S 
!==============================================================================

!.... form the LHS for the xi direction
            
!           cpu = second()
            call lhs1f( mat, Q1, Q2, Vh1, dtl, v )
!           write(*,*) 'Build LHS1 ', second()-cpu

!           cpu = second()
            if (xper) then
              call penta1p( ny, nx, ndof, mat, r, per, per2, mult, fact )
            else
              call penta1bc( ny, nx, ndof, mat, r, mult, fact )
            end if
!           write(*,*) 'Solve LHS1 ', second()-cpu

!.... form the LHS for the eta direction

!           cpu = second()
            call lhs2f( mat, Q1, Q2, Vh1, dtl, v )
!           write(*,*) 'Build LHS2 ', second()-cpu
          
!           cpu = second()
            if (yper) then
              call penta2p( ny, nx, ndof, mat, r, per, per2, mult, fact )
            else
              call penta2bc( ny, nx, ndof, mat, r, mult, fact)
            end if
!           write(*,*) 'Solve LHS2 ', second()-cpu
            
            end if       ! iLHS
!==============================================================================

!.... Inviscid: rotate back to global coordinates on the body

            if ((.not. Navier) .and. (.not. yper)) then
              call rbcfix( nx, ny, ndof, bnb, r )
            end if

!.... update v

            v = v + r
            
!.... satisfy boundary conditions on v

            call itrbc(v,vm)
            
          end do        ! loop on iter

 200      continue
          lstep = lstep + 1

!.... quit if less than 10 seconds remain

          call tremain(rtime)
          if (rtime .le. 10) goto 5000

!.... output time traces

          if (tflag) call traces(v)
          
!.... compute temporal statistics

          if(fflag .and. xper .and. istep .ne. nstep .and. &
            mod(istep-1,itout) .eq. 0 ) then
            call fstat(.false., v, vm)          ! don't print profiles
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
          call fstat(.true., v, vm)     ! print profiles
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
        subroutine rbe( rl, vl, voldl, dtl )
!
!  Backward Euler RHS:    r = vold - v - delt * r
!
!  This also works for midpoint rule.
!=============================================================================
        use stuff
        implicit none
        
        real    :: rl(ny*nx,ndof), vl(ny*nx,ndof), voldl(ny*nx,ndof)
        real    :: dtl(ny*nx)
        integer :: idof
!=============================================================================

        do idof = 1, ndof
          rl(:,idof) = voldl(:,idof) - vl(:,idof) - dtl(:) * rl(:,idof)
        end do
        
        return
        end

!=============================================================================
        subroutine rtb( rl, vl, voldl, v2oldl, dtl )
!
!  Two-step backward RHS:  r = onept33 * vold - pt33 *
!                              v2old - v - pt66 * delt * r
!
!=============================================================================
        use stuff
        implicit none
        
        real    :: rl(ny*nx,ndof), vl(ny*nx,ndof), voldl(ny*nx,ndof)
        real    :: v2oldl(ny*nx,ndof)
        real    :: dtl(ny*nx)
        integer :: idof
!=============================================================================

        do idof = 1, ndof
          rl(:,idof) = onept33 * voldl(:,idof) - &
                       pt33 * v2oldl(:,idof) - vl(:,idof) - &
                       pt66 * dtl(:) * rl(:,idof)
        end do
        
        return
        end
