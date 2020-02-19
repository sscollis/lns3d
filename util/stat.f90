!=============================================================================!
        module bspline

        integer           :: ny, kyord
        real, allocatable :: yknot(:)
        real, allocatable :: bs(:,:)

        real, external :: bsder

        end module bspline
!=============================================================================!
        module circle

        integer           :: idx
        real, allocatable :: fr(:), fi(:)

        end module circle
!=============================================================================!
        program statistics
!
!  This program computes the spatial statistics for a receptivity calculation
!
!=============================================================================!
        use const
        use bspline
        use circle
        implicit none

        real, allocatable :: v(:,:,:,:), x(:,:), y(:,:), xi(:), eta(:)
        complex, allocatable :: q(:,:,:,:), s(:,:,:), evec(:,:)
        complex :: scale
        real    :: dxi, deta, dummy
        integer :: i, j, k, idof, ier
        integer :: lstep, nx, nz, ndof
        real, allocatable :: time(:)
        real    :: Re, Ma, Pr, gamma, cv
        real    :: u(3)
        
        integer :: nver, num, n, n0, n1
        integer :: iout=10, iin=11
        
        character*80 base, name
        character*80, allocatable :: fname(:)
        character*1 ans

        real :: vmax
        integer :: jmax, ilen

!.... Cray FFT stuff
        
        real, allocatable :: table(:), work(:)

        real :: yloc

        real, external :: rtsafe, func, grad

        external funcd
!=============================================================================!
        write(*,"(/,'Compute Receptivity Stats',/)")
        write(*,"('Enter the base name ==> ',$)") 
        read(*,"(a80)") base
        write(*,"('Enter the starting version number ==> ',$)")
        read(*,*) nver
        write(*,"('Enter the number of files ==> ',$)")
        read(*,*) num

        allocate ( fname(num) )
        
!.... read in the grid file

        write(*,"('Reading grid file')")
        open(unit=10,file='grid.dat',form='unformatted',status='old',err=1000)
        read(10) nx, ny, nz
        allocate( x(ny,nx), y(ny,nx), xi(nx), eta(ny), stat=ier )
        if (ier .ne. 0) goto 2000
        read(10) (((x(j,i), i = 1, nx), j = 1, ny), k = 1, nz), &
                 (((y(j,i), i = 1, nx), j = 1, ny), k = 1, nz), &
                 ((( dummy, i = 1, nx), j = 1, ny), k = 1, nz)
        close(10)

        dxi  = one / real(nx-1)
        do i = 1, nx
          xi(i) = real(i-1) * dxi
        end do

        deta = one / real(ny-1)
        do j = 1, ny
          eta(j) = real(j-1) * deta
        end do

!.... allocate FFT stuff

        allocate ( table(100 + 4*num), work(4 + 4*num), time(num), stat=ier )
        if (ier .ne. 0) goto 2000
        
!.... read all data

        do n = nver, nver+num-1
          n0 = n - nver
          n1 = n - nver + 1
          call makename(base,n,fname(n1))
          open(unit=iin+n0,file=fname(n1),status='old',form='unformatted', &
               err=1000)
          read(iin+n0,err=1000) lstep, time(n1), nx, ny, nz, ndof, &
                                Re, Ma, Pr, gamma, cv
          ilen=index(fname(n1),' ')
          write(*,"('Reading file: ',a,' at t = ',1pe13.6,' n = ',i5)") &
                fname(n1)(1:ilen), time(n1), lstep
          if (n.eq.nver) then
            allocate( v(num,ny,nx,ndof), q(num/2+1,ny,nx,ndof), &
                      s(ny,nx,ndof), stat=ier )
            if (ier .ne. 0) goto 2000
          end if
          read(iin+n0) v(n1,:,:,:)
          close(iin+n0)
        end do

!.... write out a time-sequence

        j = 30
        i = 140
        open(iout,file='time.dat',form='formatted',status='unknown')
        do n = 1, num
          write(iout,10) time(n), ( v(n,j,i,idof), idof=1,ndof )
        end do
        close(iout)
         
!.... initialized the FFT tables

#ifdef CRAY
        call SCFFT(0, num, one/real(num), dummy, dummy, table, work, 0)
#else
        write(*,*) "Need to update FFT libraries"
        call exit(1)
#endif

!.... compute the FFT
        
        do idof = 1, ndof
          write(*,*) 'Computing FFT of dof = ',idof
          do i = 1, nx
            do j = 1, ny
#ifdef CRAY
              call SCFFT(1, num, one/real(num), v(:,j,i,idof), q(:,j,i,idof), &
                         table, work, 0)
#else
              write(*,*) "Need upt update FFT libraries"
              call exit(1)
#endif
            end do
          end do
        end do

!.... read/write raw FFT data

        write(*,"('[R]ead, [W]rite, or [I]gnore FFT data ==> ',$)")
        read(*,"(a1)") ans
        if (ans.eq.'R' .or. ans.eq.'r') then
          open(iout,file='raw.dat',form='unformatted',status='old')
          read(iout) s(:,:,:)
          close(iout)     
          q(2,:,:,:) = q(2,:,:,:) - s(:,:,:)
        else if (ans.eq.'W' .or. ans.eq.'w') then
          open(iout,file='raw.dat',form='unformatted',status='unknown')
          write(iout) q(2,:,:,:)
          close(iout)
        end if

!.... write out a profile

        if (.false.) then
        i = 250

        scale = zero
        do idof = 1, 5
          do j = 1, ny
            if ( abs(q(2,j,i,idof)) .gt. abs(scale) ) then
              scale = q(2,j,i,idof)
            end if
          end do
        end do

        scale = q(2,66,i,3) / (-2.2395809561954E-001,-2.7801657602073E-001)

        allocate( evec(ny,ndof) )
        do idof = 1, ndof
          do j = 1, ny
            evec(j,idof) = q(2,j,i,idof) / scale
          end do
        end do
        
        open(iout,file='efun.250')
        do j = 1, ny
          write(iout,10) y(j,i)-one, &
                         real(evec(j,1)), aimag(evec(j,1)), &
                         real(evec(j,2)), aimag(evec(j,2)), &
                         real(evec(j,3)), aimag(evec(j,3)), &
                         real(evec(j,5)), aimag(evec(j,5))
        end do
        close(iout)
        deallocate(evec)
        end if

!.... Wlezien's analysis

        allocate( fr(nx), fi(nx) )

        q(3,:,:,:) = 0.0

        do idof = 1, ndof
          do j = 1, ny
            u(1) = 0.0
            u(2) = 0.0
            u(3) = 0.0005
            fr = real( q(2,j,:,idof) )
            fi = aimag( q(2,j,:,idof) )
            do idx = 1, nx-2
              if (x(j,idx) .ge. 6.0 .and. x(j,idx) .le. 20.0) then
                call mnewt(10, u, 3, 1.0d-20, 1.0d-20)
                q(3,j,idx,idof) = cmplx(fr(idx)-u(1),fi(idx)-u(2))
              end if
            end do
          end do
        end do

!.... write out a profile

!       i = 222
!       i = 230
!       i = 240
        i = 250
!       i = 260
!       i = 300

        scale = zero
        do idof = 1, ndof
          do j = 1, ny
            if ( abs(q(3,j,i,idof)) .gt. abs(scale) ) then
              scale = q(3,j,i,idof)
            end if
          end do
        end do

!... M = 0.1
!
!       scale = q(3,59,i,3) / (-2.8263063565126E-002,-2.3881728178823E-001)
!       scale = q(3,59,i,3) / (-2.6039505121869E-002,-2.5712185713743E-001)
!       scale = q(3,61,i,3) / (-2.9500789369712E-002,-2.6853268873163E-001)
!       scale = q(3,61,i,3) / (-2.7472796946116E-002,-2.8111314786326E-001)
!       scale = q(3,61,i,3) / (-2.7307145524452E-002,-2.9063831841306E-001)
!       scale = q(3,63,i,3) / (-3.1945866777527E-002,-3.2989977048419E-001)

!.... M = 0.8
!
        scale = q(3,65,i,3) / (-2.2395809561954E-001,-2.7801657602073E-001)

        allocate( evec(ny,ndof) )
        do idof = 1, ndof
          do j = 1, ny
            evec(j,idof) = q(3,j,i,idof) / scale
          end do
        end do
        
        open(iout,file='efun.250')
        do j = 1, ny
          write(iout,10) y(j,i)-one, &
                         real(evec(j,1)), aimag(evec(j,1)), &
                         real(evec(j,2)), aimag(evec(j,2)), &
                         real(evec(j,3)), aimag(evec(j,3)), &
                         real(evec(j,5)), aimag(evec(j,5))
        end do
        close(iout)

        call exit(0)

!.... write out the polar plots

!       do j = 1, ny
          j = 60
          base = 'rpolar'
          call makename(base,j,name)
          open(iout,file=name,form='formatted',status='unknown')
          do i = 1, nx
            write(iout,10) x(j,i), real(q(2,j,i,1)), aimag(q(2,j,i,1))
          end do
          close(iout)

          j = 48
          base = 'upolar'
          call makename(base,j,name)
          open(iout,file=name,form='formatted',status='unknown')
          do i = 1, nx
            write(iout,10) x(j,i), real(q(2,j,i,2)), aimag(q(2,j,i,2))
          end do
          close(iout)

          j = 60
          base = 'vpolar'
          call makename(base,j,name)
          open(iout,file=name,form='formatted',status='unknown')
          do i = 1, nx
            write(iout,10) x(j,i), real(q(2,j,i,3)), aimag(q(2,j,i,3))
          end do
          close(iout)

          j = 60
          base = 'tpolar'
          call makename(base,j,name)
          open(iout,file=name,form='formatted',status='unknown')
          do i = 1, nx
            write(iout,10) x(j,i), real(q(2,j,i,5)), aimag(q(2,j,i,5))
          end do
          close(iout)
!       end do

!.... now write out the results

        open(iout,file='stat1.dat',form='unformatted',status='unknown')
        write(iout) nx, ny, nz
        write(iout) Ma, Pr, Re, time(1)
        write(iout) (((( abs(q(1,j,i,idof)), i = 1, nx), j = 1, ny),  &
                    k = 1, nz), idof = 1, ndof)
        close(iout)

        open(iout,file='stat2.dat',form='unformatted',status='unknown')
        write(iout) nx, ny, nz
        write(iout) Ma, Pr, Re, time(1)
        write(iout) (((( abs(q(2,j,i,idof)), i = 1, nx), j = 1, ny),  &
                    k = 1, nz), idof = 1, ndof)
        close(iout)
                
!.... Setup the B-spline interpolation
            
        write(*,"('Enter korder ==> ',$)")
        read(*,*) kyord
        if (kyord.eq.0) call exit(0)

        allocate ( yknot(ny+kyord), bs(ny,2) )
        call BSNAK( ny, eta, kyord, yknot)

!.... Find the maximum v at each x-station

        open(iout,file='vmax.dat',form='formatted',status='unknown')
        do i = 215, 300
          call BSINT( ny, eta,  real(q(2,:,i,3)), kyord, yknot, bs(:,1) )
          call BSINT( ny, eta, aimag(q(2,:,i,3)), kyord, yknot, bs(:,2) )

          yloc = rtsafe( funcd, eta(50), eta(70), 1.0e-14 )

!         write(*,"(3(1pe13.6,1x))") yloc, func(yloc), grad(yloc)
          write(iout,10) x(j,i), &
                         BSDER( 0, yloc, kyord, yknot, ny, bs(:,1) ), &
                         BSDER( 0, yloc, kyord, yknot, ny, bs(:,2) )
        end do
        close(iout)

!.... Find the maximum u at each x-station

        if (.true.) then
          open(iout,file='umax.dat',form='formatted',status='unknown')
          do i = 225, 300 
            call BSINT( ny, eta,  real(q(2,:,i,2)), kyord, yknot, bs(:,1) )
            call BSINT( ny, eta, aimag(q(2,:,i,2)), kyord, yknot, bs(:,2) )
  
!           do j = 2, ny
!             write (*,"(3(1pe13.6,1x))") eta(j), func(eta(j)), grad(eta(j))
!           end do
            
            yloc = rtsafe( funcd, eta(40), eta(52), 1.0e-14 )
  
!           write(*,"(3(1pe13.6,1x))") yloc, func(yloc), grad(yloc)
            write(iout,10) x(j,i), &
                            BSDER( 0, yloc, kyord, yknot, ny, bs(:,1) ), &
                            BSDER( 0, yloc, kyord, yknot, ny, bs(:,2) )
          end do
          close(iout)
        end if

        call exit(0)

 10     format(10(1pe20.13,1x))

 1000   write(*,*) 'File I/O error'
        call exit(1)
 
 2000   write(*,*) 'Insufficient memory'
        call exit(2)
        
        end
        
!=============================================================================!
        function func( x )

          use bspline

          real :: x, f, func, f1, f2
!=============================================================================!

          f1 = BSDER( 0, x, kyord, yknot, ny, bs(:,1) )
          f2 = BSDER( 0, x, kyord, yknot, ny, bs(:,2) )

          func = sqrt( f1**2 + f2**2 )

          return
        end 

!=============================================================================!
        function grad( x )

          use bspline

          real :: x, f, grad, f1, f2, g1, g2
!=============================================================================!

          f1 = BSDER( 0, x, kyord, yknot, ny, bs(:,1) )
          f2 = BSDER( 0, x, kyord, yknot, ny, bs(:,2) )

          f = sqrt( f1**2 + f2**2 )

          g1 = BSDER( 1, x, kyord, yknot, ny, bs(:,1) )
          g2 = BSDER( 1, x, kyord, yknot, ny, bs(:,2) )

          grad = ( f1*g1 + f2*g2 ) / f

          return
        end 

!=============================================================================!
        subroutine funcd( x, g, d )

          use bspline

          real :: x, f, g, d, f1, f2, g1, g2, d1, d2
!=============================================================================!

          f1 = BSDER( 0, x, kyord, yknot, ny, bs(:,1) )
          f2 = BSDER( 0, x, kyord, yknot, ny, bs(:,2) )

          f = sqrt( f1**2 + f2**2 )

          g1 = BSDER( 1, x, kyord, yknot, ny, bs(:,1) )
          g2 = BSDER( 1, x, kyord, yknot, ny, bs(:,2) )

          g = ( f1*g1 + f2*g2 ) / f

          d1 = BSDER( 2, x, kyord, yknot, ny, bs(:,1) )
          d2 = BSDER( 2, x, kyord, yknot, ny, bs(:,2) )

          d = (g1**2 + f1 * d1 + g2**2 + f2 * d2) / f - &
              ( f1*g1 + f2*g2 )**2 / f**3

          return
        end 

!=============================================================================!
      subroutine makename(base,iver,fname)
!
!.... put a version number on the filename
!
!=============================================================================!
      character*80 base, fname

      length = index(base,' ')
      fname = base
      if (iver .lt. 10) then
        write(fname(length:80),"('.R.',i1)") iver
      else if (iver .lt. 100) then
        write(fname(length:80),"('.R.',i2)") iver
      else
        write(fname(length:80),"('.R.',i3)") iver
      end if

      return
      end

!=============================================================================!
        subroutine usrfun(u,alpha,beta)
!=============================================================================!
        use const
        use circle
        implicit none

        integer, parameter :: np=3
        real u(np), alpha(np,np), beta(np)
!=============================================================================!

        alpha(1,1) =  two * ( u(1) - fr(idx) )
        alpha(1,2) =  two * ( u(2) - fi(idx) )
        alpha(1,3) = -two * u(3)
        
        alpha(2,1) =  two * ( u(1) - fr(idx+1) )
        alpha(2,2) =  two * ( u(2) - fi(idx+1) )
        alpha(2,3) = -two * u(3)

        alpha(3,1) =  two * ( u(1) - fr(idx+2) )
        alpha(3,2) =  two * ( u(2) - fi(idx+2) )
        alpha(3,3) = -two * u(3)
        
        beta(1) = u(3)**2 - (fr(idx)  -u(1))**2 - (fi(idx)  -u(2))**2
        beta(2) = u(3)**2 - (fr(idx+1)-u(1))**2 - (fi(idx+1)-u(2))**2
        beta(3) = u(3)**2 - (fr(idx+2)-u(1))**2 - (fi(idx+2)-u(2))**2

        return
        end
