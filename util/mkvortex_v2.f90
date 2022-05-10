!=============================================================================!
        program mkvortex
!
!  Generate an initial condition with one or more Oseen type vortices
!
!  Written:  11-15-99 
!  
!  Author :  Dr. Collis
! 
!  Modified by Guoquan Chen  (03-15-2001)
!
!  1. Add the radius of the vortex
!  2. Modify the expression of pressure term p
!     Taking use of the isentropy relation ,modify the expression of
!     the expression of density term Rho
!         
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
        logical :: conservative = .false.
        
!.... variables for the Oseen Vortex

        real :: xi, x0, y0, r, alpha = 1.2564312086261696770
        real :: rho, t, p, v_theta, v_r, theta, fact
        real :: x0_1, x0_2, y0_1, y0_2, amp_1, amp_2, amp, r1, r2

        real, external :: f
!=============================================================================!
        write(*,"(/,'MkVortex v2',/)") 
        write(*,*) 'WARNING:  Indices are I,J'
        write(*,"('(P)rimitive or (C)onservative ==> ',$)")
        read (*,*) ans
        if (ans.eq.'c' .or. ans.eq.'C') conservative = .true.
        write(*,"('Enter Re, Ma, Pr ==> ',$)")
        read (*,*) Re, Ma, Pr
        write(*,"('Enter amp, x0, y0, r for first vortex ==> ',$)")
        read (*,*) amp_1, x0_1, y0_1,r1
!   Judge the inputing radius 
        if (r1.le.zero) call error('MkVortex$','Incorrect r$')
        write(*,"('Enter amp, x0, y0, r  for second vortex ==> ',$)")
        read (*,*) amp_2, x0_2, y0_2,r2
        if (r2.le.zero) call error('MkVortex$','Incorrect r$')

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

!.... First Oseen vortex

        x0 = x0_1
        y0 = y0_1
        r  = r1
        amp = amp_1

        fact = one / (gamma * Ma**2)
        xi = one - exp( -alpha )

        if (amp.ne.zero) then

        !$omp parallel do private(i,r,theta,v_theta,v_r, p, &
        !$omp&                    rho, t )
        do j = 1, ny
          do i = 1, nx
            r = sqrt( (x(i,j)-x0)**2 + (y(i,j)-y0)**2 )/r1
            theta = atan2( y(i,j)-y0 , x(i,j)-x0 )
            
            if (r.eq.0) then
              v_theta = zero
              v_r = zero
              p = fact
            else
              v_theta = amp * ( one - exp( -alpha * r**2 ) ) / ( r * xi )
              v_r = zero
!              p = fact * ( one - (gamma-one) * Ma**2 / &
!                  (xi*r)**2 * f(alpha * r**2 ) )**(gamma/(gamma-one))
               p = fact * (one-amp**2 / ( fact * (r*xi)**2 )*f(alpha * r**2 ))
            end if

!            rho = ( p / fact ) ** (-gamma)
            rho = ( p / fact ) ** (one/gamma)

            t = p / ( fact * rho )

            v(1,i,j) = rho
            v(2,i,j) = cos(theta) * v_r - sin(theta) * v_theta
            v(3,i,j) = sin(theta) * v_r + cos(theta) * v_theta
            v(4,i,j) = zero
            v(5,i,j) = t

          end do
        end do

        end if

        x0 = x0_2
        y0 = y0_2
        r  = r2
        amp =  amp_2

        if (amp.ne.zero) then

        !$omp parallel do private(i,r,theta,v_theta,v_r,p,rho,t)
        do j = 1, ny
          do i = 1, nx
            r = sqrt( (x(i,j)-x0)**2 + (y(i,j)-y0)**2 )/r2
            theta = atan2( y(i,j)-y0 , x(i,j)-x0 )
            
            if (r.eq.0) then
              v_theta = zero
              v_r = zero
              p = fact
            else
              v_theta = amp * ( one - exp( -alpha* r**2 ) ) / ( r * xi )
              v_r = zero

!              p = fact * ( one - (gamma-one) * Ma**2 / &
!                  (xi*r)**2 * f(alpha * r**2 ) )**(gamma/(gamma-one))
               p = fact * (one - amp**2 / (fact * (r*xi)**2 )*f(alpha * r**2 ))
            end if

!            rho = ( p / fact ) ** (-gamma)
            rho = ( p / fact ) ** (one/gamma)

            t = p / ( fact * rho )

            v(1,i,j) = v(1,i,j) + rho
            v(2,i,j) = v(2,i,j) + cos(theta) * v_r - sin(theta) * v_theta
            v(3,i,j) = v(3,i,j) + sin(theta) * v_r + cos(theta) * v_theta
            v(4,i,j) = zero
            v(5,i,j) = v(5,i,j) + t

          end do
        end do

        end if

!.... convert to conservative variables if needed

        if (conservative) then
          fact = one / (gamma * Ma**2)
          do j = 1, ny
            do i = 1, nx
              v(5,i,j) = fact * v(1,i,j) * v(5,i,j) / (gamma-one) + &
                         pt5 * v(1,i,j) * ( v(2,i,j)**2 + v(3,i,j)**2 + &
                         v(4,i,j)**2 )
              v(2,i,j) = v(1,i,j) * v(2,i,j) 
              v(3,i,j) = v(1,i,j) * v(3,i,j)
              v(4,i,j) = v(1,i,j) * v(4,i,j)
            end do
          end do
        end if

!.... write the restart file

        open(unit=10, file='output.R.0', form='unformatted', status='unknown')
        write(10,err=1000) lstep, time, nx, ny, nz, ndof, &
                           Re, Ma, Pr, gamma, cv
        write(10,err=1000) v
        close(10)

!.... write out the potential boundary values

        open(unit=10, file='top.pot')
        do i = 1, nx
          write(10,*,err=110) v(:,i,ny)
        end do
        close(10)

        open(unit=10, file='bottom.pot')
        do i = 1, nx
          write(10,*,err=110) v(:,i,1)
        end do
        close(10)

        open(unit=10, file='left.pot')
        do j = 1, ny
          write(10,*,err=110) v(:,1,j)
        end do
        close(10)
        
        open(unit=10, file='right.pot')
        do j = 1, ny
          write(10,*,err=120) v(:,nx,j)
        end do
        close(10)

        stop

100     call error('MkVortex$','Reading grid.dat file$')
110     call error('MkVortex$','Writing left.pot file$')
120     call error('MkVortex$','Writing right.pot file$')
1000    call error('MkVortex$','Writing restart file$')
        
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
        character(*) name, msg

        loc = index(name,'$')-1
        write(*,"(/,'=====================================================')")
        write(*,"('Error in --> ',a)") name(1:loc)
        loc = index(msg,'$')-1
        write(*,"('-----------> ',a)") msg(1:loc)
        write(*,"('=====================================================',/)")
        
        call exit(1)
        
        stop
        end
