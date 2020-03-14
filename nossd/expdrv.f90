!=============================================================================!
        subroutine expdrv()
!  
!  Driver for explicit time advancement of the compressible Navier-Stokes eqn.
!  
!  This is the routine is capable of 505 MFLOPS 
!  with ~ 7.0 micro-sec per/node/step on a 384 x 127 mesh
!
!  Revised:  6-11-96
!=============================================================================!
        use global
        use local
        implicit none
        
        integer, external :: igetver
        
        integer :: ier
        character*80 :: name
        real :: rtime, told, cpu
        real, external :: second
        
        real, allocatable :: v(:), vold(:)
        real :: vint(ny*nx*ndof)
        real :: r1(ny*nx*ndof), r2(ny*nx*ndof)
        real :: r3(ny*nx*ndof), r4(ny*nx*ndof)
        real :: dtl(ny*nx)
!=============================================================================!

!.... allocate the storage area for the disturbance field

        allocate (v(ny*nx*ndof), vold(ny*nx*ndof), STAT=ier)
        if (ier .ne. 0) call error('expdrv$','Insufficient Memory for v$')

!.... generate the initial condition

        call genini(v, vold)

!.... open the history file

        iver = igetver( base, 'h'//char(0) ) + 1
        call addver(base, 'h'//char(0), iver, filen, lfile)
        open(unit=ihist, file=filen, form='formatted', status='unknown')

!.... satisfy boundary conditions on v

        call itrbc(v,vm)

!.... compute the CFL

        if (linear.eq.1) then
          call dtcfl(vm, dtl)
          call mlocal('deallocate')
        end if

!.... loop through the time steps

        do istep = 1, nstep
          
          if (linear.eq.0) call dtcfl(v, dtl)

!.... -----------------------> predictor phase <-----------------------

          vold = v
          
!.... --------------------> multi-corrector phase <--------------------

!.... Fourth-Order Runge-Kutta time advancement

          iter = 1
          vint = v
          told = time
          if (linear.eq.1) then
            call rk(vint, r1, vm, x, y, dtl)
          else
            call nrk(vint, r1, dtl)
          end if

          iter = 2
          vint = v + pt5 * r1
          time = told + pt5 * Delt
          call itrbc(vint,vm)
          if (linear.eq.1) then
            call rk(vint, r2, vm, x, y, dtl)
          else
            call nrk(vint, r2, dtl)
          end if
          
          iter = 3
          vint = v + pt5 * r2
          call itrbc(vint,vm)
          if (linear.eq.1) then
            call rk(vint, r3, vm, x, y, dtl)
          else
            call nrk(vint, r3, dtl)
          end if
          
          iter = 4
          vint = v + r3
          time = told + Delt
          call itrbc(vint,vm)
          if (linear.eq.1) then
            call rk(vint, r4, vm, x, y, dtl)
          else
            call nrk(vint, r4, dtl)
          end if
            
          v = vold + &
              0.1666666666667 * r1 + &
              0.3333333333333 * (r2 + r3) + &
              0.1666666666667 * r4

!.... satisfy boundary conditions on v

          call itrbc(v,vm)

          lstep = lstep + 1

!.... quit if less than 10 seconds remain

#ifdef CRAY
          call tremain(rtime)
          if (rtime .le. 10) goto 5000
#endif

!.... compute temporal statistics

          if(fflag .and. xper .and. mod(lstep-1,itout) .eq. 0 .and. &
             istep .ne. nstep) then
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
            close(ifile)

            write(*,10) filen(1:index(filen,char(0))-1), lstep, time
            write(ihist,10) filen(1:index(filen,char(0))-1), lstep, time

          end if
          
        end do

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
        close(ifile)

        write(*,10) filen(1:index(filen,char(0))-1), lstep, time
        write(ihist,10) filen(1:index(filen,char(0))-1), lstep, time

!.... close the history file

        close(ihist)

        return
        
  10    format(' Wrote ',a,' at step = ',i7,' time = ',1pe13.6)

        end
