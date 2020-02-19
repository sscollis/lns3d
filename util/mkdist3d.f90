!=============================================================================!
        program mkdist3d
!
!  Make an initial condition for the 3-D linearized Navier-Stokes solver
!
!  Revised: 3-15-2001  Changed to IJ format
!
!=============================================================================!
        use const
        implicit none
        
        integer ::  nx, ny, nz, ndof = 5, ier
        real    ::  Re, Ma, Pr, gamma = 1.4, cv = 716.5
        real    ::  time = 0.0, Lx, Ly, Lz, kx, ky, kz, kk
        integer ::  lstep = 0

!.... data
        
        complex, allocatable :: v(:,:,:)
        real, allocatable :: vm(:,:,:)

!.... mesh

        real, allocatable :: x(:,:), y(:,:)
        real :: dxi, deta

        real :: tmp

!.... variables for the acoustic pulse

        complex :: c1, c2, c3, c4, rho, u1, u2, u3, t, p, rhom, tm, cm
        integer :: i, j, k

        character*80 :: line, code='Initial$'
        character*1 ans
!=============================================================================!

!.... read in the grid file

        open(unit=10,file='grid.dat',form='unformatted',status='old',err=100)
        read(10) nx, ny, nz
        allocate( x(nx,ny), y(nx,ny) )
        read(10) (((x(i,j), i = 1, nx), j = 1, ny), k = 1, nz), &
                 (((y(i,j), i = 1, nx), j = 1, ny), k = 1, nz), &
                 (((   tmp, i = 1, nx), j = 1, ny), k = 1, nz)
        close(10)

!.... read in the mean field

        open(unit=10,file='mean.dat',form='unformatted',status='old',err=200)
        read(10) lstep, time, nx, ny, nz, ndof, &
                 Re, Ma, Pr, gamma, cv
        write(*,10) Re, Ma, Pr, gamma
10      format('Re = ',1e11.3e3,' Ma = ',1e11.3e3,' Pr = ',1e11.3e3, &
               ' Gamma = ',1e11.3e3)
        allocate( vm(ndof,nx,ny) )
        read(10) vm
        close(10)

!.... allocate the storage area for the field

        allocate (v(ndof,nx,ny), STAT=ier)
        if (ier .ne. 0) call error(code,'Insufficient Memory for v$')

        v = zero

        write(*,"('[Z]ero disturbance, [W]ave, [S]pike ==> ',$)")
        read(*,"(a1)") ans

        if (ans.eq.'W' .or. ans.eq.'w') then
          write(*,"('Enter Lx, Ly, Lz ==> ',$)")
          read(*,*) Lx, Ly, Lz

          if (Lx.eq.zero) then
            kx = zero
          else
            kx = two * pi / Lx
          end if

          if (Ly.eq.zero) then
            ky = zero
          else
            ky = two * pi / Ly
          end if

          if (Lz.eq.zero) then
            kz = zero
          else  
            kz = two * pi / Lz
          end if

!         kk = sqrt( kx**2 + kz**2 )
          kk = sqrt( ky**2 + kz**2 )
          do i = 1, nx
            do j = 1, ny
              rhom = vm(1,i,j)
              tm   = vm(5,i,j)
              cm   = sqrt( tm ) / Ma 
              
              c1 = zero
              c2 = zero
!!$           c3 = cos( kx * x(i,j) ) * cos( ky * y(i,j) ) - &
!!$                sin( kx * x(i,j) ) * sin( ky * y(i,j) ) + &
!!$                im * ( cos( kx * x(i,j) ) * sin( ky * y(i,j) ) + &
!!$                          sin( kx * x(i,j) ) * cos( ky * y(i,j) ) )
              c3 = cos( kx * x(i,j) )
              c4 =  zero

              rho = ( -c1 + pt5 * ( c3 + c4 ) ) / cm**2

              u1 = ( c3 - c4 ) * pt5 / ( rhom * cm )
              u2 = c2 / ( rhom * cm )
              u3 = zero

!!$              u1  = c2 / ( rhom * cm )
!!$              if (ky.eq.zero) then
!!$                u2 = zero
!!$              else
!!$                u2  = ky / kk * ( c3 - c4 ) * pt5 / ( rhom * cm )
!!$              end if
!!$              if (kz.eq.zero) then
!!$                u3 = zero
!!$              else
!!$                u3  = kz / kk * ( c3 - c4 ) * pt5 / ( rhom * cm )
!!$              end if

              p   = ( c3 + c4 ) * pt5
              t   = ( gamma * Ma**2 * p - tm * rho ) / rhom
                  
              v(1,i,j) = rho
              v(2,i,j) = u1
              v(3,i,j) = u2
              v(4,i,j) = u3
              v(5,i,j) = t
            end do
!           write(11,10) x(1,i), v(1,i,1), v(1,i,2), v(1,i,3), v(1,i,5)
          end do
        else if (ans.eq.'S' .or. ans.eq.'s') then
          write(*,"('[X]-spike, [Y]-spike, or [P]oint spike ==> ',$)")
          read(*,"(a1)") ans
          if (ans.eq.'X' .or. ans.eq.'x') then
            write(*,"('Enter (i) for spike ==> ',$)")
            read(*,*) i
            do j = 1, ny
              v(:,i,j) = one
            end do
          else if (ans.eq.'Y' .or. ans.eq.'y') then
            write(*,"('Enter (j) for spike ==> ',$)")
            read(*,*) j
            do i = 1, nx
              v(:,i,j) = one
            end do
          else
            write(*,"('Enter (i,j) for spike ==> ',$)")
            read(*,*) i,j
            v(:,i,j) = one
          end if
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
        character(*) name, msg

        loc = index(name,'$')-1
        write(*,"(/,'*****************************************************')")
        write(*,"('Error in --> ',a)") name(1:loc)
        loc = index(msg,'$')-1
        write(*,"('-----------> ',a)") msg(1:loc)
        write(*,"('*****************************************************',/)")
        
        call exit(1)
        
        stop
        end
