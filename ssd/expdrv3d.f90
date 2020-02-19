!=============================================================================!
        subroutine expdrv3D()
!  
!  Driver for explicit time advancement of the three-dimensinoal 
!  compressible Navier-Stokes eqn.
!  
!  Revised:  4-19-96
!=============================================================================!
        use stuff
        use local
        implicit none
        
        integer, external :: igetver
        
        integer :: ier
        character*80 :: name, code='ExpDrv3D$'
        real :: rtime, told
        
        complex, allocatable :: v(:), vold(:), vint(:)
        complex, allocatable :: r1(:), r2(:), r3(:), r4(:)
        real, allocatable    :: dtl(:)
!=============================================================================!
!.... allocate storage for the local time step

        istep = 0
        allocate (dtl(ny*nx), STAT=ier)
        if (ier .ne. 0) call error(code,'Insufficient Memory for dtl$')

!.... compute the CFL

        call dtcfl(vm, dtl)
        call mlocal('deallocate')

!.... allocate the storage area for the disturbance field

        allocate (v(ny*nx*ndof), vold(ny*nx*ndof), vint(ny*nx*ndof), &
                  r1(ny*nx*ndof), r2(ny*nx*ndof), &
                  r3(ny*nx*ndof), r4(ny*nx*ndof), STAT=ier)
        if (ier .ne. 0) call error(code,'Insufficient Memory for v$')

!.... generate the initial condition

        call genini3d(v, vold)

!.... open the history file

        iver = igetver( base, 'h'//char(0) ) + 1
        call addver(base, 'h'//char(0), iver, filen, lfile)
        open(unit=ihist, file=filen, form='formatted', status='unknown')

!.... satisfy boundary conditions on v

        call itrbc3d(v,vm)

!.... loop through the time steps

        do istep = 1, nstep
          
!.... -----------------------> predictor phase <-----------------------

          vold = v
          
!.... --------------------> multi-corrector phase <--------------------

!.... Fourth-Order Runge-Kutta time advancement

          iter = 1
          vint = v
          told = time
          call rk3d(vint, r1, vm, x, y, dtl)

          iter = 2
          vint = v + pt5 * r1
          time = told + pt5 * Delt
          call itrbc3d(vint,vm)
          call rk3d(vint, r2, vm, x, y, dtl)
          
          iter = 3
          vint = v + pt5 * r2
          call itrbc3d(vint,vm)
          call rk3d(vint, r3, vm, x, y, dtl)
          
          iter = 4
          vint = v + r3
          time = told + Delt
          call itrbc3d(vint,vm)
          call rk3d(vint, r4, vm, x, y, dtl)
            
          v = vold + &
              0.1666666666667 * r1 + &
              0.3333333333333 * (r2 + r3) + &
              0.1666666666667 * r4

!.... satisfy boundary conditions on v

          call itrbc3d(v,vm)

          lstep = lstep + 1

!.... quit if less than 10 seconds remain

          call tremain(rtime)
          if (rtime .le. 10) goto 5000
          
!.... compute temporal statistics

!         if(fflag .and. xper .and. mod(lstep-1,itout) .eq. 0 .and. &
!             istep .ne. nstep) then
!           call fstat(.false., v, vm)          ! don't print profiles
!          end if

!.... compute spatial statistics

!         if (sflag) then
!           call sstat(v, vm, x, y)
!         end if
          
!.... write out restart and q files 

          if ( mod(lstep,ntout) .eq. 0 .and. istep .ne. nstep) then
            
            iver = igetver( base, 'R'//char(0) ) + 1
            call addver(base, 'R'//char(0), iver, filen, lfile)
            open(unit=ifile, file=filen, form='unformatted', &
                  status='unknown')
            write(ifile) lstep, time, nx, ny, nz, ndof, &
                         Re, Ma, Pr, gamma, cv
            write(ifile) v
            close(ifile)

            write(*,10) filen(1:index(filen,char(0))-1), lstep, time
            write(ihist,10) filen(1:index(filen,char(0))-1), lstep, time

          end if
          
        end do

!.... ---------------------->  Post Processing  <----------------------

5000    continue

!.... compute final statistics

!       if (fflag .and. xper) then
!         call fstat(.true., v, vm)     ! print profiles
!       end if

!.... write the restart file

        iver = igetver( base, 'R'//char(0) ) + 1
        call addver(base, 'R'//char(0), iver, filen, lfile)
        open(unit=ifile, file=filen, form='unformatted', &
             status='unknown')
        write(ifile) lstep, time, nx, ny, nz, ndof, &
                     Re, Ma, Pr, gamma, cv
        write(ifile) v
        close(ifile)

        write(*,10) filen(1:index(filen,char(0))-1), lstep, time
        write(ihist,10) filen(1:index(filen,char(0))-1), lstep, time

!.... close the history file

        close(ihist)

        return
        
  10    format(' Wrote ',a,' at step = ',i7,' time = ',1pe13.6)

        end
