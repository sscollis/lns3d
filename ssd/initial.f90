!=============================================================================!
        module constants

          real, parameter :: zero    = 0.0000000000000000000d+0
          real, parameter :: pt25    = 2.5000000000000000000d-1
          real, parameter :: pt33    = 3.3333333333333333333d-1
          real, parameter :: pt5     = 5.0000000000000000000d-1
          real, parameter :: pt66    = 6.6666666666666666666d-1
          real, parameter :: one     = 1.0000000000000000000d+0
          real, parameter :: onept25 = 1.2500000000000000000d+0   
          real, parameter :: onept33 = 1.3333333333333333333d+0
          real, parameter :: onept5  = 1.5000000000000000000d+0
          real, parameter :: two     = 2.0000000000000000000d+0
          real, parameter :: three   = 3.0000000000000000000d+0
          real, parameter :: four    = 4.0000000000000000000d+0
          real, parameter :: pi      = 3.1415926535897932385d+0
          complex, parameter :: im   = (0.0,1.0)
          
        end module constants
!=============================================================================!
        program initial
!
!  Make an initial condition for the 3-D linearized Navier-Stokes solver
!
!=============================================================================!
        use constants
        implicit none
        
        integer ::  nx, ny, nz, ndof = 5, ier
        real    ::  Re, Ma, Pr, gamma = 1.4, cv = 716.5
        real    ::  time = 0.0, Lx, Ly, Lz, kx, ky, kz, kk, omega
        integer ::  lstep = 0

!.... data
        
        complex, allocatable :: v(:,:,:)
        real, allocatable :: vm(:,:,:)

!.... mesh

        real, allocatable :: x(:,:), y(:,:)
        real :: dxi, deta

        real :: tmp

!.... variables for the acoustic pulse

        complex :: c1, c2, c3, c4, rho, u1, u2, u3, t, p
        real    :: rhom, tm, um, cm
        integer :: i, j, k

        character*80 :: line, code='Initial$'
        character*1 ans
!=============================================================================!

!.... read in the grid file

        open(unit=10,file='grid.dat',form='unformatted',status='old',err=100)
        read(10) nx, ny, nz
        allocate( x(ny,nx), y(ny,nx) )
        read(10) (((x(j,i), i = 1, nx), j = 1, ny), k = 1, nz), &
                 (((y(j,i), i = 1, nx), j = 1, ny), k = 1, nz), &
                 (((   tmp, i = 1, nx), j = 1, ny), k = 1, nz)
        close(10)

!.... read in the mean field

        open(unit=10,file='mean.dat',form='unformatted',status='old',err=200)
        read(10) lstep, time, nx, ny, nz, ndof, &
                 Re, Ma, Pr, gamma, cv
        write(*,*) 'Re = ',Re,' Ma = ',Ma,' Pr = ',Pr,' Gamma = ',gamma
        allocate( vm(ny,nx,ndof) )
        read(10) vm
        close(10)

!.... allocate the storage area for the field

        allocate (v(ny,nx,ndof), STAT=ier)
        if (ier .ne. 0) call error(code,'Insufficient Memory for v$')

        v = zero

        write(*,"('[Z]ero disturbance, [A]coustic wave ==> ',$)")
        read(*,"(a1)") ans

        if (ans.eq.'A' .or. ans.eq.'a') then
          write(*,"('Enter omega ==> ',$)")
          read(*,*) omega
            
          rhom = one
          um   = one
          tm   = one
          cm   = sqrt( tm ) / Ma
          kk   = omega / (cm + um)
              
          do i = 1, nx
            do j = 1, ny
              
              c3   = cos( kk * x(j,i) ) + im * sin( kk * x(j,i) )
              
              v(j,i,1) = pt5 * c3 / cm**2
              v(j,i,2) = c3 * pt5 / ( rhom * cm )
              v(j,i,3) = zero
              v(j,i,4) = zero
              v(j,i,5) = (gamma*Ma**2 * c3 * pt5 - &
                          tm * pt5 * c3 / cm**2) / rhom
            end do
!           write(11,10) x(1,i), v(1,i,1), v(1,i,2), v(1,i,3), v(1,i,5)
10          format(8(1pe13.6,1x))
          end do
        end if
        
!.... write the restart file

        open(unit=10,file='output.R.0',form='unformatted',status='unknown')
        write(10,err=300) 0, zero, nx, ny, nz, ndof, &
                          Re, Ma, Pr, gamma, cv
        write(10,err=300) v
        close(10)

        stop

100     call error(code,'Reading grid.dat file$')
200     call error(code,'Reading mean.dat file$')
300     call error(code,'Writing restart file$')
        
        end
        
!=============================================================================!
        subroutine error(name,msg)
!  
!       Generic error handler 
!  
!=============================================================================!
        implicit none
        
        integer loc
        character*80 name, msg

        loc = index(name,'$')-1
        write(*,"(/,'*****************************************************')")
        write(*,"('Error in --> ',a)") name(1:loc)
        loc = index(msg,'$')-1
        write(*,"('-----------> ',a)") msg(1:loc)
        write(*,"('*****************************************************',/)")
        
        call exit(1)
        
        stop
        end
