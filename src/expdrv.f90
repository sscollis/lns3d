!=============================================================================!
        subroutine expdrv()
!  
!  Driver for explicit time advancement of the compressible Navier-Stokes eqn.
!  
!  This is the routine is capable of 505 MFLOPS 
!  with ~ 7.0 micro-sec per/node/step on a 384 x 127 mesh
!
!  Revised:  6-11-96
!            11-10-99  Changed index order for execution on SGI O2k
!=============================================================================!
        use global
        use local2d
        implicit none
        
        integer, external :: igetver
        
        integer :: ier, i, j, idof, mem = 0
        character(80) :: name, code='ImpDrv$'
        real :: rtime, told

        real*4 :: cpul
        real*4, external :: second

#ifdef USE_STACK        
        real :: v(ndof,nx,ny), vold(ndof,nx,ny)
        !$sgi distribute v(*,*,block), vold(*,*,block)
        real :: vint(ndof,nx,ny)
        !$sgi distribute vint(*,*,block)
        real :: r1(ndof,nx,ny), r2(ndof,nx,ny)
        !$sgi distribute r1(*,*,block), r2(*,*,block)
        real :: r3(ndof,nx,ny), r4(ndof,nx,ny)
        !$sgi distribute r3(*,*,block), r4(*,*,block)
        real :: dtl(nx,ny)
        !$sgi distribute dtl(*,block)
#else
        real, allocatable :: v(:,:,:), vold(:,:,:)
        real, allocatable :: vint(:,:,:)
        real, allocatable :: r1(:,:,:), r2(:,:,:)
        real, allocatable :: r3(:,:,:), r4(:,:,:)
        real, allocatable :: dtl(:,:)
#endif
!=============================================================================!
        istep = 0

        cpul = second()

!.... allocate the storage area for the disturbance field

        mem = mem + 7*ndof*nx*ny + nx*ny
        allocate (v(ndof,nx,ny), vold(ndof,nx,ny), STAT=ier)
        if (ier .ne. 0) call error(code,'Insufficient Memory for v$')
        allocate (vint(ndof,nx,ny), STAT=ier)
        if (ier .ne. 0) call error(code,'Insufficient Memory for vint$')
        allocate (r1(ndof,nx,ny), r2(ndof,nx,ny), STAT=ier)
        if (ier .ne. 0) call error(code,'Insufficient Memory for r1$')
        allocate (r3(ndof,nx,ny), r4(ndof,nx,ny), STAT=ier)
        if (ier .ne. 0) call error(code,'Insufficient Memory for r3$')
        allocate (dtl(nx,ny), STAT=ier)
        if (ier .ne. 0) call error(code,'Insufficient Memory for dtl$')
        write(*,"(' ExpDrv allocated: ',1pe13.6,' words')") float(mem)

!.... initialize for 1st touch allocation on SGI

        !$doacross local(i,idof)
        !$omp parallel do private(i,idof)
        do j = 1, ny
          do i = 1, nx
            do idof = 1, ndof
              v(idof,i,j) = zero;    vold(idof,i,j) = zero
              vint(idof,i,j) = zero; r1(idof,i,j) = zero; 
              r2(idof,i,j) = zero;   r3(idof,i,j) = zero;   
              r4(idof,i,j) = zero;   dtl(i,j) = zero
            end do
          end do
        end do

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
          call set_time_step(vm, dtl)
        else
          call set_time_step(v, dtl)
        end if

!.... output statistics

!       call resstat(r1, dtl)
!       call vstat(v, dtl)

!       write(*,*) 'Initialize ',second()-cpul
        cpul = second()
        
!.... loop through the time steps

        do istep = 1, nstep

          if (linear.eq.0) call set_time_step(v, dtl)

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
          !write(*,*) 'Predictor ',second() - cpul
          cpul = second()
          
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
          if (linear.eq.1) then
            call rk(vint, r1, vm, x, y, dtl)
          else
            call nrk(vint, r1, dtl)
          end if
          !write(*,*) 'Stage 1 ',second() - cpul
          cpul = second()

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
          call itrbc(vint,vm)
          if (linear.eq.1) then
            call rk(vint, r2, vm, x, y, dtl)
          else
            call nrk(vint, r2, dtl)
          end if
          !write(*,*) 'Stage 2 ',second() - cpul
          cpul = second()
          
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
          call itrbc(vint,vm)
          if (linear.eq.1) then
            call rk(vint, r3, vm, x, y, dtl)
          else
            call nrk(vint, r3, dtl)
          end if
          !write(*,*) 'Stage 3 ',second() - cpul
          cpul = second()
          
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
          call itrbc(vint,vm)
          if (linear.eq.1) then
            call rk(vint, r4, vm, x, y, dtl)
          else
            call nrk(vint, r4, dtl)
            cpul = second()
          end if
          !write(*,*) 'Stage 4 ',second() - cpul
!         cpul = second()
            
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
!         write(*,*) 'Update ',second() - cpul
!         cpul = second()

!.... satisfy boundary conditions on v

          call itrbc(v,vm)
!         write(*,*) 'Stage 4 ',second() - cpul
          cpul = second()

          lstep = lstep + 1

!.... output statistics

          if (mod(lstep,itout).eq.0 .and. iter.eq.4) then
!           call resstat(r4, dtl)
            call vstat(v, dtl)
          end if
        
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
