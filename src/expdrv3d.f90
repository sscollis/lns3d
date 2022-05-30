!=============================================================================!
        subroutine expdrv3D()
!  
!  Driver for explicit time advancement of the three-dimensinoal 
!  compressible Navier-Stokes eqn.
!  
!  Revised:  4-19-96
!            12-15-00  Changed index order for execution on SGI O2k
!=============================================================================!
        use global
        use local2d
        use local3d
        implicit none
        
        integer, external :: igetver
        
        integer :: i, j, idof, ier
        character(80) :: name, code='ExpDrv3D$'
        real :: rtime, told
        
        complex, allocatable :: v(:,:,:), vold(:,:,:), vint(:,:,:)
        complex, allocatable :: r1(:,:,:), r2(:,:,:), r3(:,:,:), r4(:,:,:)
        real, allocatable    :: dtl(:,:)
!=============================================================================!
!.... allocate storage for the local time step

        istep = 0
        allocate (dtl(nx,ny), STAT=ier)
        if (ier .ne. 0) call error(code,'Insufficient Memory for dtl$')

!.... compute the CFL

        call set_time_step(vm, dtl)

!.... allocate the storage area for the disturbance field

        allocate (v(ndof,nx,ny), vold(ndof,nx,ny), vint(ndof,nx,ny), &
                  r1(ndof,nx,ny), r2(ndof,nx,ny), &
                  r3(ndof,nx,ny), r4(ndof,nx,ny), STAT=ier)
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

          !$doacross local(i,idof)
          !$omp parallel do private(i,idof)
          do j = 1, ny
            do i = 1, nx
              do idof = 1, ndof
                vold(idof,i,j) = v(idof,i,j)
              end do
            end do
          end do
          
!.... --------------------> multi-corrector phase <--------------------

!.... Fourth-Order Runge-Kutta time advancement

          iter = 1
          !$doacross local(i,idof)
          !$omp parallel do private(i,idof)
          do j = 1, ny
            do i = 1, nx
              do idof = 1, ndof
                vint(idof,i,j) = v(idof,i,j)
              end do
            end do
          end do
          told = time
          call rk3d(vint, r1, vm, x, y, dtl)

          iter = 2
          !$doacross local(i,idof)
          !$omp parallel do private(i,idof)
          do j = 1, ny
            do i = 1, nx
              do idof = 1, ndof
                vint(idof,i,j) = v(idof,i,j) + pt5 * r1(idof,i,j)
              end do
            end do
          end do
          time = told + pt5 * Delt
          call itrbc3d(vint,vm)
          call rk3d(vint, r2, vm, x, y, dtl)
          
          iter = 3
          !$doacross local(i,idof)
          !$omp parallel do private(i,idof)
          do j = 1, ny
            do i = 1, nx
              do idof = 1, ndof
                vint(idof,i,j) = v(idof,i,j) + pt5 * r2(idof,i,j)
              end do
            end do
          end do
          call itrbc3d(vint,vm)
          call rk3d(vint, r3, vm, x, y, dtl)
          
          iter = 4
          !$doacross local(i,idof)
          !$omp parallel do private(i,idof)
          do j = 1, ny
            do i = 1, nx
              do idof = 1, ndof
                vint(idof,i,j) = v(idof,i,j) + r3(idof,i,j)
              end do
            end do
          end do
          time = told + Delt
          call itrbc3d(vint,vm)
          call rk3d(vint, r4, vm, x, y, dtl)
            
          !$doacross local(i,idof)
          !$omp parallel do private(i,idof)
          do j = 1, ny
            do i = 1, nx
              do idof = 1, ndof
                v(idof,i,j) = vold(idof,i,j) + &
                     0.16666666666666667 * r1(idof,i,j) + &
                     0.33333333333333333 * (r2(idof,i,j) + r3(idof,i,j)) + &
                     0.16666666666666667 * r4(idof,i,j)
              end do
            end do
          end do

!.... satisfy boundary conditions on v

          call itrbc3d(v,vm)

          lstep = lstep + 1

!.... output RHS statistics

        if (mod(lstep,itout).eq.0 .and. iter.eq.4) then
!         call resstat3D(r4, dtl)
          call vstat3d( v, dtl )
        end if
        
!.... quit if less than 10 seconds remain

#ifdef CRAY
          call tremain(rtime)
          if (rtime .le. 10) goto 5000
#endif

!.... compute temporal statistics

!         if(fflag .and. xper .and. mod(lstep-1,itout) .eq. 0 .and. &
!             istep .ne. nstep) then
!           call fstat2d(.false., v, vm)        ! don't print profiles
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
!         call fstat2d(.true., v, vm)   ! print profiles
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
