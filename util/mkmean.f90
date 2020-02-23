!=============================================================================!
        program main
!
!  MKMEAN:  Make a mean file for the linearized Navier-Stokes solver
!
!  Revised:  3-15-01  Changed to IJ on input and output format, still JI
!                     internally
!
!=============================================================================!
        use const
        implicit none
        
        real, allocatable :: vm(:,:,:), x(:,:), y(:,:)
        real    :: tmp
        integer :: nx, ny, nz, ndof=5
        
        real    :: Ma, Re, Pr, gamma=1.4, gamma1=0.4, cv=716.5
        integer :: i, j, k, idof
!=============================================================================!
        write(*,"(/,'MakeMean ',/)") 
        write(*,"('Enter Re, Ma, Pr ==> ',$)")
        read (*,*) Re, Ma, Pr

!.... read in the grid file

        open(unit=10,file='grid.dat',form='unformatted',status='old')
        read(10) nx, ny, nz
        write(*,"('Read grid with dimensions (',i3,',',i3,',',i3,')')") &
                nx,ny,nz
        allocate( x(ny,nx), y(ny,nx) )
        read(10) (((x(j,i), i = 1, nx), j = 1, ny), k = 1, nz), &
                 (((y(j,i), i = 1, nx), j = 1, ny), k = 1, nz), &
                 (((   tmp, i = 1, nx), j = 1, ny), k = 1, nz)
        close(10)
        
        allocate( vm(ny,nx,ndof) )

        call mean( vm, x, y, nx, ny, ndof, Ma, Re, Pr )
                
        open(unit=10, file='mean.dat', form='unformatted', status='unknown')
        write(10) 0, zero, nx, ny, nz, ndof, Re, Ma, Pr, gamma, cv
        write(10) (((vm(j,i,idof), idof=1,ndof), i=1,nx), j=1,ny)
        close(10)
        
        stop
        end

!=============================================================================!
        module mkmean
!=============================================================================!
        use const
        real, allocatable :: xb(:), yb(:), el(:), r(:)
        real, allocatable :: bn1(:), bn2(:)
        real :: xi, yi
        integer :: jel

        contains

!=============================================================================!
        function uint(s)
!=============================================================================!
        implicit none
        real :: uint, xj, yj, s
!=============================================================================!
        xj = xb(jel) + bn2(jel) * s
        yj = yb(jel) - bn1(jel) * s
        
        uint = ( (xi-xj) ) /  &
     &         ( (xi-xj)**2 + (yi-yj)**2 )
        
        return
        end function uint

!=============================================================================!
        function vint(s)
!=============================================================================!
        implicit none
        real :: vint, xj, yj, s
!=============================================================================!
        xj = xb(jel) + bn2(jel) * s
        yj = yb(jel) - bn1(jel) * s
        
        vint = ( (yi-yj) ) /  &
     &         ( (xi-xj)**2 + (yi-yj)**2 )
        
        return
        end function vint

        end module mkmean

!=============================================================================!
        subroutine mean(vml, x, y, nx, ny, ndof, Ma, Re, Pr)
!=============================================================================!
        use const
        use mkmean
        implicit none
        
        integer :: nx, ny, ndof
        real    :: vml(ny,nx,ndof), x(ny,nx), y(ny,nx)
        real    :: Ma, Re, Pr, gamma=1.4, gamma1=0.4, cv=716.5
        
        real, allocatable :: ym(:,:), vt(:,:,:), vs(:,:,:)
        real :: ymaxm
        integer :: i, j, k, nxm, nym, ndofm, ier, idof
        
        real :: b, rbeta, u, v, t, rho, p, c, alpha

        character*1 ans
        character*80 name

        real :: result, errest, dummy
        integer :: np, nel

        real :: xii, eta, u1, u2, x0, y0, th, tmp
!=============================================================================!

        write(*,"('Use [M]ean flow profile,         [N]o flow, ', / , &
                & '    [U]nit flow,                 [P]otential flow,' , / , &
                & '    [C]urved mesn flow profile   ==> ',$)")
        read(*,"(a1)") ans
        
!.... input from file

        if (ans.eq.'M' .or. ans.eq.'m') then

!.... read the mean profile and spline to the computational grid

        write(*,"('Enter profile name ==> ',$)")
        read(*,*) name
        open (unit=10, file=name, form='formatted', status='old', err=1000)
        
        nxm = 1
        nym = 0
 20     continue
        read(10,*,end=30) tmp
        nym = nym + 1
        goto 20
 30     continue
        rewind(10)

!       read (10,*,err=1000,end=1010) nxm, nym, ndofm, ymaxm
!       if (ndofm .ne. ndof) then
!         write (*,*) 'NDOF in mean field is incorrect ',ndofm
!         call exit(1)
!       end if
        
        allocate( ym(nym,nxm), vt(nym,nxm,ndof), vs(nym,nxm,ndof), STAT=ier )
        if (ier .ne. 0) then
          write(*,*) 'Error allocating mean field'
          call exit(1)
        end if
        
        do i = 1, nxm
          do j = 1, nym
            read (10,*,err=1000,end=1010) ym(j,i), (vt(j,i,k),k=1,ndof)
          end do
          call SPLINE(nym, ym(1,i), vt(1,i,1), vs(1,i,1))
          call SPLINE(nym, ym(1,i), vt(1,i,2), vs(1,i,2))
          call SPLINE(nym, ym(1,i), vt(1,i,3), vs(1,i,3))
          call SPLINE(nym, ym(1,i), vt(1,i,4), vs(1,i,4))
          call SPLINE(nym, ym(1,i), vt(1,i,5), vs(1,i,5))
        end do
        
        close (10)
        
        ymaxm = ym(nym,1)

!.... Evaluate the mean field on the disturbance grid

        do i = 1, nx
          do j = 1, ny
            if (y(j,i) .lt. ymaxm) then
              call SPEVAL(nym, ym(1,1), vt(1,1,1), vs(1,1,1), &
                          y(j,i), vml(j,i,1))
              call SPEVAL(nym, ym(1,1), vt(1,1,2), vs(1,1,2), &
                          y(j,i), vml(j,i,2))
              call SPEVAL(nym, ym(1,1), vt(1,1,3), vs(1,1,3), &
                          y(j,i), vml(j,i,3))
              call SPEVAL(nym, ym(1,1), vt(1,1,4), vs(1,1,4), &
                          y(j,i), vml(j,i,4))
              call SPEVAL(nym, ym(1,1), vt(1,1,5), vs(1,1,5), &
                          y(j,i), vml(j,i,5))
            else
              vml(j,i,1) = vt(nym,1,1)
              vml(j,i,2) = vt(nym,1,2)
              vml(j,i,3) = vt(nym,1,3)
              vml(j,i,4) = vt(nym,1,4)
              vml(j,i,5) = vt(nym,1,5)
            end if
!           write (*,"(7(e13.6,1x))") y(j,i), (vml(j,i,k), k = 1, ndof)
          end do
        end do

!.... make the flow parallel

        write(*,"('WARNING: using parallel flow assumption...')")
        vml(:,:,3) = zero

        deallocate( ym, vt, vs )

        else if (ans.eq.'P' .or. ans.eq.'p') then       ! potential soln
        
          write(*,"('Enter the sweep angle ==> ',$)")
          read(*,*) alpha

          b = one / pi
          if (Ma.lt.one) then
            rbeta = sqrt(one - Ma**2)
          else
            rbeta = sqrt(Ma**2 - one)
          end if
          
          do i = 1, nx
            do j = 1, ny

              xi = x(j,i)
              yi = y(j,i)

!.... potential flow for a Rankine Oval

!             u   = one + (x(j,i)-b)/pi/rbeta**2/((x(j,i)-b)**2 + &
!                   rbeta**2*y(j,i)**2)
!             v   = rbeta*y(j,i)/pi/((x(j,i)-b)**2+rbeta**2*y(j,i)**2)
!             t   = one - pt5*gamma1*Ma**2*( u**2 + v**2 - one )
!             rho = t**(one/gamma1)
!             p   = rho * t / ( gamma * Ma**2 )
!             c   = sqrt(t)/Ma

!.... potential flow for parabolic cylinder (currently incompressible)

              eta = sqrt(pt5 - xi + pt5*sqrt((two*xi-one)**2 + four*yi**2))
              xii  = yi / eta
              u  = (eta**2 - eta + xii**2) / (xii**2 + eta**2)
              v  = xii / (xii**2 + eta**2)
              t   = one - pt5*gamma1*Ma**2*( u**2 + v**2 - one )
              rho = t**(one/gamma1)
              p   = rho * t / ( gamma * Ma**2 )
              c   = sqrt(t)/Ma

              vml(j,i,1) = rho
              vml(j,i,2) = u
              vml(j,i,3) = v
              vml(j,i,4) = tan(alpha)
              vml(j,i,5) = t
            end do
          end do
        
        else if (ans.eq.'U' .or. ans.eq.'u') then
        
          write(*,"('Enter the sweep angle ==> ',$)")
          read(*,*) alpha

          vml(:,:,1) = one
          vml(:,:,2) = one
          vml(:,:,3) = zero
          vml(:,:,4) = tan(alpha)
          vml(:,:,5) = one

        else if (ans.eq.'C' .or. ans.eq.'c') then

        write(*,"('Enter profile name ==> ',$)")
        read(*,*) name
        open (unit=10, file=name, form='formatted', status='old', err=1000)
        
        nxm = 1
        nym = 0
 40     continue
        read(10,*,end=50) tmp
        nym = nym + 1
        goto 40
 50     continue
        rewind(10)

!       read (10,*,err=1000,end=1010) nxm, nym, ndofm, ymaxm
!       if (ndofm .ne. ndof) then
!         write (*,*) 'NDOF in mean field is incorrect ',ndofm
!         call exit(1)
!       end if
        
        allocate( ym(nym,nxm), vt(nym,nxm,ndof), vs(nym,nxm,ndof), STAT=ier )
        if (ier .ne. 0) then
          write(*,*) 'Error allocating mean field'
          call exit(1)
        end if
        
        write(*,*) nxm, nym
        do i = 1, nxm
          do j = 1, nym
            read (10,*,err=1000,end=1010) ym(j,i), (vt(j,i,k),k=1,ndof)
          end do
          call SPLINE(nym, ym(1,i), vt(1,i,1), vs(1,i,1))
          call SPLINE(nym, ym(1,i), vt(1,i,2), vs(1,i,2))
          call SPLINE(nym, ym(1,i), vt(1,i,3), vs(1,i,3))
          call SPLINE(nym, ym(1,i), vt(1,i,4), vs(1,i,4))
          call SPLINE(nym, ym(1,i), vt(1,i,5), vs(1,i,5))
        end do
        close (10)

        ymaxm = ym(nym,1)
        
!.... Evaluate the mean field on the disturbance grid

        do i = 1, 1
          do j = 1, ny
            if (y(j,i) .lt. ymaxm) then
              call SPEVAL(nym, ym(1,1), vt(1,1,1), vs(1,1,1), &
                          y(j,1), vml(j,i,1))
              call SPEVAL(nym, ym(1,1), vt(1,1,2), vs(1,1,2), &
                          y(j,1), vml(j,i,2))
              call SPEVAL(nym, ym(1,1), vt(1,1,3), vs(1,1,3), &
                          y(j,1), vml(j,i,3))
              call SPEVAL(nym, ym(1,1), vt(1,1,4), vs(1,1,4), &
                          y(j,1), vml(j,i,4))
              call SPEVAL(nym, ym(1,1), vt(1,1,5), vs(1,1,5), &
                          y(j,1), vml(j,i,5))
            else
              vml(j,i,1) = vt(nym,1,1)
              vml(j,i,2) = vt(nym,1,2)
              vml(j,i,3) = vt(nym,1,3)
              vml(j,i,4) = vt(nym,1,4)
              vml(j,i,5) = vt(nym,1,5)
            end if
          end do
        end do

        write(*,"('WARNING: using parallel flow assumption...')")
        vml(:,:,3) = zero

!.... rotate to the new coordinate system

        write(*,"('Enter (x0, y0) ==> ',$)")
        read(*,*) x0, y0
        do i = 1, nx
          th = atan2( y(1,i)+y0, x(1,i)+x0 )
          do j = 1, ny
            vml(j,i,1) = vml(j,1,1)
            vml(j,i,2) =  sin(th) * vml(j,1,2) + cos(th) * vml(j,1,3)
            vml(j,i,3) = -cos(th) * vml(j,1,2) + sin(th) * vml(j,1,3)
            vml(j,i,4) = vml(j,1,4)
            vml(j,i,5) = vml(j,1,5)
          end do
        end do
        
        deallocate( ym, vt, vs )

        else            ! no flow

          vml(:,:,1) = one
          vml(:,:,2) = zero
          vml(:,:,3) = zero
          vml(:,:,4) = zero
          vml(:,:,5) = one
        
        end if
        
        return

1000    write(*,"('Error reading file mean.dat\')")
        call exit(1)
1010    write(*,"('End of file encountered in profile.dat\')")
        call exit(1)

        end

!-------------------------------------------------------------------------
!         b = one / pi
!         if (Ma.lt.one) then
!           rbeta = sqrt(one - Ma**2)
!         else
!           rbeta = sqrt(Ma**2 - one)
!         end if
!         
!         do i = 1, nx
!           do j = 1, ny
!             u   = one + (x(j,i)-b)/pi/rbeta**2/((x(j,i)-b)**2 + &
!                   rbeta**2*y(j,i)**2)
!             v   = rbeta*y(j,i)/pi/((x(j,i)-b)**2+rbeta**2*y(j,i)**2)
!             t   = one - pt5*gamma1*Ma**2*( u**2 + v**2 - one )
!             rho = t**(one/gamma1)
!             p   = rho * t / ( gamma * Ma**2 )
!             c   = sqrt(t)/Ma
!
!             vml(j,i,1) = rho
!             vml(j,i,2) = u
!             vml(j,i,3) = v
!             vml(j,i,4) = zero
!             vml(j,i,5) = t
!           end do
!         end do
!-------------------------------------------------------------------------
!         open(unit=10,file='bndry.dat')
!         read(10,*) np
!         nel = np - 1
!         allocate( xb(np), yb(np), el(nel), r(nel), bn1(nel), bn2(nel) )
!         do i = 1, np
!           read(10,*) xb(i), yb(i)
!         end do
!         close(10)
!
!         open(10,file='pot.dat')
!         do jel = 1, nel
!           read(10,*) dummy, dummy, el(jel), bn1(jel), bn2(jel), r(jel)
!         end do
!         close(10)
!
!         do i = 250, 250
!           do j = 10,ny
!             write(*,*) i, j
!
!             xi = x(j,i)
!             yi = y(j,i)
!
!             u = one
!             v = zero
!             do jel = 1, nel
!               write(*,*) jel
!               call QDNG( uint, 0, el(jel), 1.0e-10, 0.0, &
!     &                     result, errest)
!               u = u + r(jel)*pt5/pi * result
!               call QDNG( vint, 0, el(jel), 1.0e-10, 0.0, &
!     &                     result, errest)
!               v = v +  r(jel)*pt5/pi * result
!             end do
!
!             t   = one - pt5*gamma1*Ma**2*( u**2 + v**2 - one )
!             rho = t**(one/gamma1)
!             p   = rho * t / ( gamma * Ma**2 )
!             c   = sqrt(t)/Ma
!
!             vml(j,i,1) = rho
!             vml(j,i,2) = u
!             vml(j,i,3) = v
!             vml(j,i,4) = zero
!             vml(j,i,5) = t
!             write(11,"(8(1pe13.6,1x))") xi, rho, u, v, t
!           end do
!         end do
!-------------------------------------------------------------------------
!         open(unit=10,file='bndry.dat')
!         read(10,*) np
!         nel = np - 1
!         allocate( xb(np), yb(np), el(nel), r(nel), bn1(nel), bn2(nel) )
!         do i = 1, np
!           read(10,*) xb(i), yb(i)
!         end do
!         close(10)
!
!         open(10,file='pot.dat')
!         do jel = 1, nel
!           read(10,*) dummy, dummy, el(jel), bn1(jel), bn2(jel), r(jel)
!         end do
!         close(10)
!
!         do i = 1, nx
!           do j = 1,ny
!
!             xi = x(j,i)
!             yi = y(j,i)
!
!             u = one
!             v = zero
!             do jel = 1, nel
!               call QDNG( uint, 0, el(jel), 1.0e-10, 0.0, &
!     &                     result, errest)
!               u = u + r(jel)*pt5/pi * result
!               call QDNG( vint, 0, el(jel), 1.0e-10, 0.0, &
!     &                     result, errest)
!               v = v +  r(jel)*pt5/pi * result
!             end do
!
!             t   = one - pt5*gamma1*Ma**2*( u**2 + v**2 - one )
!             rho = t**(one/gamma1)
!             p   = rho * t / ( gamma * Ma**2 )
!             c   = sqrt(t)/Ma
!
!             vml(j,i,1) = rho
!             vml(j,i,2) = u
!             vml(j,i,3) = v
!             vml(j,i,4) = zero
!             vml(j,i,5) = t
!
!             write(11,"(8(1pe13.6,1x))") xi, rho, u, v, t
!           end do
!         end do

