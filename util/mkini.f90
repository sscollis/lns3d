!=============================================================================!
        program mkini
!
!  Generate a simple initial conditions
!
!  11-14-99:  Switched i,j indices [SSC]
!=============================================================================!
        use const 
        implicit none
        
!.... flow variables

        integer :: nx, ny, nz=1, ndof=5, ier
        real    :: Re, Ma, Pr, gamma=1.4, cv=716.5
        real    :: time=0.0
        integer :: lstep=0

        real, allocatable :: v(:,:,:)
        
!.... mesh

        real, allocatable :: x(:,:), y(:,:)
        real :: tmp, fact

        integer i, j, k, node, ians
        character(80) line
        character(1) ans
        integer conserve
        
!.... variables for the acoustic pulse

        real :: c1, c2, c3, c4
        real :: rho, u1, u2, u3, t, p
        real :: rhom, tm, cm, eps
        real :: amp, x0, sigma
        
!=============================================================================!
        write(*,"(/,'MakeIni ',/)") 
        write(*,*) 'WARNING:  Indices are I,J'
        write(*,"('Enter Re, Ma, Pr ==> ',$)")
        read (*,*) Re, Ma, Pr
        write(*,"('Write Primative(0) or Conservative(1) ==> ',$)")
        read (*,*) conserve

!.... read in the grid file

        open(unit=10,file='grid.dat',form='unformatted',status='old',err=100)
        read(10) nx, ny, nz
        allocate( x(nx,ny), y(nx,ny) )
        read(10) (((x(i,j), i = 1, nx), j = 1, ny), k = 1, nz), &
                 (((y(i,j), i = 1, nx), j = 1, ny), k = 1, nz), &
                 (((   tmp, i = 1, nx), j = 1, ny), k = 1, nz)
        close(10)

!.... allocate the storage area for the field

        allocate (v(ndof,nx,ny), STAT=ier)
        if (ier .ne. 0) &
          call error('makeini$','Insufficient Memory for v$')

        write(*,"('[E]NSA file, [N]o flow, [U]niform flow, ',&
       &          '[G]aussian pulse ==> ',$)")
        read(*,*) ans
        
        if (ans.eq.'E' .or. ans.eq.'e') then

!.... read and ENSA file

          open(10,file='pert.dat',status='old')
          read(10,"(a80)") line
          read(10,*) lstep, time, gamma, cv
          do i = 1, nx
            do j = 1, ny
              read(10,*) node, v(1,i,j), v(2,i,j), v(3,i,j), v(5,i,j)
              v(4,i,j) = zero
            end do
          end do
        
        else if (ans.eq.'N' .or. ans.eq.'n') then

!.... no flow

          v(1,:,:) = one
          v(2,:,:) = zero
          v(3,:,:) = zero
          v(4,:,:) = zero
          v(5,:,:) = one

        else if (ans.eq.'U' .or. ans.eq.'u') then

!.... uniform flow in the +x direction

          v(1,:,:) = one
          v(2,:,:) = one
          v(3,:,:) = zero
          v(4,:,:) = zero
          v(5,:,:) = one

        else if (ans.eq.'G' .or. ans.eq.'g') then

!.... Gaussian acoustic pulse:  set the initial field to no mean flow

          v(1,:,:) = one
          v(2,:,:) = zero
          v(3,:,:) = zero
          v(4,:,:) = zero
          v(5,:,:) = one

!.... add on an acoustic pulse

          write(*,"('amp, x0, and sigma ==> ',$)")
          read(*,*) amp, x0, sigma
           do i = 1, nx
             do j = 1, ny
               rhom = v(1,i,j)
               tm   = v(5,i,j)
               cm   = one / Ma * sqrt( tm )

               c1  = zero
               c2  = zero
               c3  = exp( -pt5 * ( (x(i,j) - x0)/sigma )**2 )
               c4  = zero

               rho = ( -c1 + pt5 * ( c3 + c4 ) ) / cm**2
               u1  = ( c3 - c4 ) * pt5 / ( rhom * cm )
               u2  = c2 / ( rhom * cm )
               p   = ( c3 + c4 ) * pt5
               t   = ( gamma * Ma**2 * p - tm * rho ) / rhom

               v(1,i,j) = v(1,i,j) + amp * rho
               v(2,i,j) = v(2,i,j) + amp * u1
               v(3,i,j) = v(3,i,j) + amp * u2
               v(4,i,j) = v(4,i,j)
               v(5,i,j) = v(5,i,j) + amp * t
             end do
           end do

        end if

!.... convert to conservative form

        if (conserve.eq.1) then
          fact = one / (gamma * Ma**2)
          do j = 1, ny
            do i = 1, nx
              v(2,i,j) = v(2,i,j) * v(1,i,j)
              v(3,i,j) = v(3,i,j) * v(1,i,j)
              v(4,i,j) = v(4,i,j) * v(1,i,j) 
              v(5,i,j) = v(1,i,j) * v(5,i,j) * fact / (gamma-one) +  &
                         pt5 * ( v(2,i,j)**2 + v(3,i,j)**2 + &
                                 v(4,i,j)**2) / &
                         v(1,i,j)
            end do
          end do
        end if

!.... write the restart file

        open(unit=10, file='output.R.0', form='unformatted', &
             status='unknown')
        write(10,err=1000) lstep, time, nx, ny, nz, ndof, &
                           Re, Ma, Pr, gamma, cv
        write(10,err=1000) v
        close(10)

        stop

100     call error('makeini$','Reading grid.dat file$')
1000    call error('makeini$','writing restart file$')
        
        end
        
!=============================================================================!
        subroutine error(name,msg)
!  
!       Generic error handler 
!  
!=============================================================================!
        implicit none
        
        integer loc
        character(*) name, msg

        loc = index(name,'$')-1
        write(*,"(/,'===============================================')")
        write(*,"('Error in --> ',a)") name(1:loc)
        loc = index(msg,'$')-1
        write(*,"('-----------> ',a)") msg(1:loc)
        write(*,"('===============================================',/)")
        
        call exit(1)
        
        stop
        end
