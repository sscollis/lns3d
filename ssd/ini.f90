!=============================================================================!
        program makedist
!
!  Make a disturbance file for the linearized Navier-Stokes solver
!
!=============================================================================!
        use constants
        implicit none
        
        integer ::  nx, ny, nz, ndof = 5, ier
        real    ::  Re, Ma, Pr, gamma = 1.4, cv = 716.5
        real    ::  time = 0.0
        integer ::  lstep = 0

!.... data
        
        real, allocatable :: v(:,:,:), vm(:,:,:)

!.... mesh

        real, allocatable :: x(:,:), y(:,:)
        real :: dxi, deta

        real :: tmp

!.... variables for the acoustic pulse

        real    :: c1, c2, c3, c4, rho, u1, u2, u3, t, p, rhom, tm, cm, eps
        integer :: i, j, k

        character*80 line
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
        allocate( vm(ny,nx,ndof) )
        read(10) vm
        close(10)

!.... allocate the storage area for the field

        allocate (v(ny,nx,ndof), STAT=ier)
        if (ier .ne. 0) call error('makeini$','Insufficient Memory for v$')

        v = zero

        write(*,"('[Z]ero disturbance, [G]aussian pulse, [S]pike ==> ',$)")
        read(*,"(a1)") ans

        if (ans.eq.'G' .or. ans.eq.'g') then
          write(*,"('[X]-pulse or [Y]-pulse ==> ',$)")
          read(*,"(a1)") ans
          if (ans.eq.'X' .or. ans.eq.'x') then
            do i = 1, nx
              do j = 1, ny
    
                rhom = vm(j,i,1)
                tm   = vm(j,i,5)
                cm   = one / Ma * sqrt( tm )
      
                c1  =  zero
                c2  =  zero
                c3  =  exp( -pt5 * ( (x(j,i) - 0.5)/0.05 )**2 )
                c4  =  zero ! exp( -pt5 * ( (x(j,i) - 0.5)/0.05 )**2 )
      
                rho = ( -c1 + pt5 * ( c3 + c4 ) ) / cm**2
                u1  = ( c3 - c4 ) * pt5 / ( rhom * cm )
                u2  = c2 / ( rhom * cm )
                p   = ( c3 + c4 ) * pt5
                t   = ( gamma * Ma**2 * p - tm * rho ) / rhom
                
                v(j,i,1) = rho
                v(j,i,2) = u1
                v(j,i,3) = u2
                v(j,i,4) = zero
                v(j,i,5) = t
              end do
              write(11,10) x(1,i), v(1,i,1), v(1,i,2), v(1,i,3), v(1,i,5)
 10           format(8(1pe13.6,1x))
            end do
          else 
!           eps = 0.0001
!           write(*,*) "Adding to zero mean flow with amplitude ", eps
            do i = 1, nx
              do j = 1, ny
                rhom = vm(j,i,1)
                tm   = vm(j,i,5)
                cm   = one / Ma * sqrt( tm )
      
                c1 =  zero
                c2 =  zero
                c3 =  zero
                c4 =  exp( -pt5 * ( (y(j,i) - 0.5) / 0.1 )**2 )

!               if ( mod(j-1,2).eq.0 ) then
!                 c3 =  exp( -pt5 * ( (y(j,i) - 20.0) / 2.0 )**2 )
!               else
!                 c3 = -exp( -pt5 * ( (y(j,i) - 20.0) / 2.0 )**2 )
!               end if
    
                rho = ( -c1 + pt5 * ( c3 + c4 ) ) / cm**2
                u1  = c2 / ( rhom * cm ) 
                u2  = ( c3 - c4 ) * pt5 / ( rhom * cm )
                p   = ( c3 + c4 ) * pt5
                t   = ( gamma * Ma**2 * p - tm * rho ) / rhom
                  
                v(j,i,1) = rho
                v(j,i,2) = u1
                v(j,i,3) = u2
                v(j,i,4) = zero
                v(j,i,5) = t

!               v(j,i,1) = one + eps * rho
!               v(j,i,2) = eps * u1
!               v(j,i,3) = eps * u2
!               v(j,i,4) = zero
!               v(j,i,5) = one + eps * t
              end do
              write(11,10) x(1,i), v(1,i,1), v(1,i,2), v(1,i,3), v(1,i,5)
            end do
          end if
        else if (ans.eq.'S' .or. ans.eq.'s') then
          write(*,"('[X]-spike, [Y]-spike, or [P]oint spike ==> ',$)")
          read(*,"(a1)") ans
          if (ans.eq.'X' .or. ans.eq.'x') then
            write(*,"('Enter (i) for spike ==> ',$)")
            read(*,*) i
            do j = 1, ny
              v(j,i,:) = one
            end do
          else if (ans.eq.'Y' .or. ans.eq.'y') then
            write(*,"('Enter (j) for spike ==> ',$)")
            read(*,*) j
            do i = 1, nx
              v(j,i,:) = one
            end do
          else
            write(*,"('Enter (i,j) for spike ==> ',$)")
            read(*,*) i,j
            v(j,i,:) = one
          end if
        end if
        
!.... write the restart file

        open(unit=10,file='output.R.0',form='unformatted',status='unknown')
        write(10,err=300) 0, zero, nx, ny, nz, ndof, &
                          Re, Ma, Pr, gamma, cv
        write(10,err=300) v
        close(10)

        stop

100     call error('makedist$','Reading grid.dat file$')
200     call error('makedist$','Reading mean.dat file$')
300     call error('makedist$','Writing restart file$')
        
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
