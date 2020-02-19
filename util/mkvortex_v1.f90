!=============================================================================!
        program mkvortex
!
!  Generate an initial condition with one or more Oseen type vortices
!
!  Written:  11-15-99 
!=============================================================================!
        use const 
        implicit none
        
!.... flow variables

        integer :: nx, ny, nz=1, ndof=5, ier
        real    :: Re, Ma, Pr, gamma=1.4, cv=716.5
        real    :: time=0.0
        integer :: lstep=0

        real, allocatable :: v(:,:,:)
!$sgi distribute v(*,*,block)
        
!.... mesh

        real, allocatable :: x(:,:), y(:,:)
!$sgi distribute x(*,block), y(*,block)
        real :: tmp

        integer i, j, k
        character*80 line
        character*1 ans
        
!.... variables for the Oseen Vortex

        real :: L, circ, x0, y0, r, alpha = 1.2564312086261696770
        real :: rho, t, p, a_inf, rho_inf, v_theta, v_r, theta, fact

        real, external :: f
!=============================================================================!
        write(*,"(/,'MkVortex',/)") 
        write(*,"('Enter Re, Ma, Pr ==> ',$)")
        read (*,*) Re, Ma, Pr
        write(*,"('Enter x0, y0, L, Circ ==> ',$)")
        read (*,*) x0, y0, l, circ

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
        if (ier .ne. 0) call error('MkVortex$','Insufficient Memory for v$')

!$omp parallel do private(i)
        do j = 1, ny
          do i = 1, ny
            v(:,i,j) = zero
          end do
        end do

!.... Oseen vortex

        fact = one / (gamma * Ma**2)

!$omp parallel do private( i, r, theta, a_inf, rho_inf, v_theta, v_r, p, &
!$omp&                     rho, t )
        do j = 1, ny
          do i = 1, ny
            r = sqrt( (x(i,j)-x0)**2 + (y(i,j)-y0)**2 )
            theta = atan2( y(i,j)-y0 , x(i,j)-x0 )
            
            a_inf   = sqrt(one)/Ma
            rho_inf = one

            if (r.eq.0) then
              v_theta = zero
              v_r = zero
              p = rho_inf * a_inf**2 / gamma
            else
              v_theta = pt5 * circ / (pi * r) * &
                        ( one - exp( -alpha*(r/L)**2 ) )
              v_r = zero
              p = rho_inf * a_inf**2 / gamma * ( one - (gamma-one)*circ**2/ &
                  (two*a_inf*pi*r)**2 * f(alpha*r**2/L**2) )** &
                  (gamma/(gamma-one))
            end if

            rho = ( p / fact ) ** (-gamma)

            t = p / ( fact * rho )

            v(1,i,j) = rho
            v(2,i,j) = cos(theta) * v_r - sin(theta) * v_theta
            v(3,i,j) = sin(theta) * v_r + cos(theta) * v_theta
            v(4,i,j) = zero
            v(5,i,j) = t

          end do
        end do

        circ = -circ
        x0 = -x0

!$omp parallel do private( i, r, theta, a_inf, rho_inf, v_theta, v_r, p, &
!$omp&                     rho, t )
        do j = 1, ny
          do i = 1, ny
            r = sqrt( (x(i,j)-x0)**2 + (y(i,j)-y0)**2 )
            theta = atan2( y(i,j)-y0 , x(i,j)-x0 )
            
            a_inf   = sqrt(one)/Ma
            rho_inf = one

            if (r.eq.0) then
              v_theta = zero
              v_r = zero
              p = rho_inf * a_inf**2 / gamma
            else
              v_theta = pt5 * circ / (pi * r) * &
                        ( one - exp( -alpha*(r/L)**2 ) )
              v_r = zero
              p = rho_inf * a_inf**2 / gamma * ( one - (gamma-one)*circ**2/ &
                  (two*a_inf*pi*r)**2 * f(alpha*r**2/L**2) )** &
                  (gamma/(gamma-one))
            end if

            rho = ( p / fact ) ** (-gamma)

            t = p / ( fact * rho )

            v(1,i,j) = v(1,i,j) + rho
            v(2,i,j) = v(2,i,j) + cos(theta) * v_r - sin(theta) * v_theta
            v(3,i,j) = v(3,i,j) + sin(theta) * v_r + cos(theta) * v_theta
            v(4,i,j) = zero
            v(5,i,j) = v(5,i,j) + t

          end do
        end do

!.... write the restart file

        open(unit=10, file='output.R.0', form='unformatted', status='unknown')
        write(10,err=1000) lstep, time, nx, ny, nz, ndof, &
                           Re, Ma, Pr, gamma, cv
        write(10,err=1000) v
        close(10)

        stop

100     call error('MkVortex$','Reading grid.dat file$')
1000    call error('MkVortex$','writing restart file$')
        
        end

!=============================================================================!
        real function f(x)
!
!=============================================================================!
        implicit none

        real*8, external :: dEi
        real :: x

        if (x.gt.sqrt(200.0)) then
          f = 0.5
        else
          f = 0.5 - exp(-x) + 0.5*exp(-2.0*x) + x * dEi(-2.0*x) - x * dEi(-x)
        end if

        return
        end function f

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
        write(*,"(/,'=====================================================')")
        write(*,"('Error in --> ',a)") name(1:loc)
        loc = index(msg,'$')-1
        write(*,"('-----------> ',a)") msg(1:loc)
        write(*,"('=====================================================',/)")
        
        call exit(1)
        
        stop
        end
