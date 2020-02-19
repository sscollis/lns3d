!=============================================================================!
      program csubwave
!  
!  This program removes a forcing wave of a given frequency from a restart
!  file and outputs a plot3d file.
!  
!=============================================================================!
      use const
      implicit none

      real, allocatable :: vm(:,:,:)
      complex, allocatable :: v(:,:,:), q(:,:,:)

      real :: Ma, Re, Pr, gamma, cv, Rgas, time
      integer :: lstep, nx, ny, nz, ndof
      integer :: nxm, nym, nzm, ndofm
      integer :: i, j, k, idof, ix

      character*80 file, temp
      integer :: iloc, iend
#ifndef __GFORTRAN__
      integer, external :: iargc
#endif
      complex, allocatable :: rho(:,:), u1(:,:), u2(:,:), u3(:,:)
      complex, allocatable :: t(:,:), p(:,:)
      complex, allocatable :: c1(:,:), c2(:,:), c3(:,:), c4(:,:)

      real, allocatable :: rhom(:,:), um(:,:), tm(:,:), cm(:,:)
      real, allocatable :: x(:,:), y(:,:)

      real, allocatable :: wamp(:)
      logical :: lamp

      real :: tmp
      integer :: itmp

      real :: kk, alpha, delta, x0, r, rho0, c0, omega, th, tl
!=============================================================================!

!.... read the disturbance file

      if ( iargc() .gt. 0 ) then
        call getarg(1,temp)
        file = temp
      else
        file = 'output.R.1'
        write (*,"('Enter file name [',a,']? ',$)") file(1:index(file,' '))
        read (*,"(a80)") temp
        if (temp(1:1) .ne. ' ') file = temp
      end if
10    open(unit=10, file=file, form='unformatted', status='old', err=20)
      goto 30
20    write (*,"('>> Error opening [',a,'] ',/)") file(1:index(file,' '))
      write (*,"('Enter file name [',a,']? ',$)") file(1:index(file,' '))
      read (*,"(a80)") temp
      if (temp(1:1) .ne. ' ') file = temp
      goto 10
30    continue

      read(10) lstep, time, nx, ny, nz, ndof, Re, Ma, Pr, gamma, cv
      if (nz .ne. 1) then
        write(*,"(' Error:  nz <> 1 ')")
        call exit(1)
      end if
      allocate( v(ny,nx,ndof) )
      read(10) v
      close(10)
      write(*,"('Read disturbance field for ',a)") file(1:index(file,' '))
      write(*,"('  Time = ',1pe10.3,'  step = ',i6)") time, lstep

!.... read in the mean flow file

      open(unit=10, file='mean.dat', form='unformatted', status='old', &
        err=1000)
      read(10,err=1000) itmp, tmp, nxm, nym, nzm, ndofm, &
                        tmp, tmp, tmp, tmp, tmp
      if( nx.ne.nxm .or. ny.ne.nym .or. ndof.ne.ndofm) then
        write(*,"('>> Error:  Mean and disturbance fields do not agree')")
        call exit(1)
      end if
      allocate( vm(ny,nx,ndof) )
      read(10,err=1000) vm
      close(10)

!.... read in the grid file

      open (unit=10, file='grid.dat', form='unformatted', status='unknown')
      read(10) nx, ny, nz
      allocate (x(ny,nx), y(ny,nx))
      read(10) (((x(j,i), i = 1, nx), j = 1, ny), k = 1, nz), &
               (((y(j,i), i = 1, nx), j = 1, ny), k = 1, nz), &
               (((   tmp, i = 1, nx), j = 1, ny), k = 1, nz)
      close(10)

!.... read in the forcing wave amplitude

!      inquire( file='amp.dat', exist=lamp )
!      if (lamp) then
!        allocate( wamp(nx) )
!        write(*,*) 'Reading amp.dat'
!       open(66,file='amp.dat',form='formatted',status='unknown')
!       do i = 1, nx
!         read(66,*) tmp, wamp(i)
!       end do
!       close(66)
!      end if
      
!.... compute the acoustic wave

      allocate( rho(ny,nx), u1(ny,nx), u2(ny,nx), u3(ny,nx), &
                t(ny,nx), p(ny,nx) )
      allocate( rhom(ny,nx), um(ny,nx), tm(ny,nx), cm(ny,nx) )
      allocate( c1(ny,nx), c2(ny,nx), c3(ny,nx), c4(ny,nx) )

      rhom = vm(:,:,1)
      um   = vm(:,:,2)
      tm   = vm(:,:,5)
      cm   = sqrt( tm ) / Ma

      rho = v(:,:,1)
      u1  = v(:,:,2)
      u2  = v(:,:,3)
      u3  = v(:,:,4)
      t   = v(:,:,5)
      p   = ( rhom * t + tm * rho ) / (gamma * Ma**2)

!      if (.not. lamp) then
        write(*,"('Correcting for Re = ',1pe13.6,', Pr ',1pe13.6)") Re, Pr
        rho0  = one
        delta = (one/Re) * pt5 / rho0 * ( onept33 + (gamma - one) / Pr )
!      end if

      x0 = x(ny,1)
      write(*,"('Origin, X0 = ',1pe13.6)") x0

      write(*,"('Enter the angular frequency ==> ',$)")
      read(*,*) omega

      do i = 1, nx
        do j = 1, ny
          kk      = omega / ( cm(j,i) + um(j,i) )
          alpha   = kk**2 * delta / ( cm(j,i) + um(j,i) )
          c1(j,i) = zero
          c2(j,i) = zero
!         if (lamp) then
!           c3(j,i) = wamp(i)*( cos( kk * x(j,i) ) + im * sin( kk * x(j,i) ) )
!         else
            c3(j,i) = exp( -alpha * (x(j,i) - x0) ) * &
                      ( cos( kk * x(j,i) ) + im * sin( kk * x(j,i) ) )
!         end if
          c4(j,i) = zero
        end do
      end do

!     j = ny
!     do i = 1, nx
!       kk      = omega / ( cm(j,i) + um(j,i) )
!       alpha   = kk**2 * delta / ( cm(j,i) + um(j,i) )
!       write(90,"(2(1pe20.13,1x))") x(j,i), exp( -alpha * (x(j,i) - x0) )
!     end do

      rho = rho - ( -c1 + pt5 * ( c3 + c4 ) ) / cm**2
      u1  = u1 - ( c3 - c4 ) * pt5 / ( rhom * cm )
      u2  = u2 -  c2 / ( rhom * cm )
      p   = p - ( c3 + c4 ) * pt5
      t   = ( gamma * Ma**2 * p - tm * rho ) / rhom

!      rho = ( -c1 + pt5 * ( c3 + c4 ) ) / cm**2
!      u1  = ( c3 - c4 ) * pt5 / ( rhom * cm )
!      u2  = c2 / ( rhom * cm )
!      p   = ( c3 + c4 ) * pt5
!      t   = ( gamma * Ma**2 * p - tm * rho ) / rhom

      deallocate( rhom, tm, cm )
      deallocate( c1, c2, c3, c4 )

!.... form the output variables

      allocate( q(ny,nx,ndof) )

      q(:,:,1) = rho
      q(:,:,2) = u1
      q(:,:,3) = u2
      q(:,:,4) = p
      q(:,:,5) = t

      deallocate( rho, u1, u2, u3, t, p )

!.... write out a plot3d file

      iloc = index(file,'.R.')
      iend = index(file,' ')-1
      if (iloc .ne. 0) then
        temp = file(1:iloc)//'d'//file(iloc+2:iend)
      else
        temp = file(1:iend)//'.d'
      end if

      if ( iargc() .eq. 0 ) then
        file = temp
        write (*,"(/,'Enter PLOT3D file name [',a,']? ',$)")    &
          file(1:index(file,' '))
        read (*,"(a80)") temp
        if (temp(1:1) .ne. ' ') file = temp
      else if ( iargc() .eq. 1 ) then
        file = temp
      else if ( iargc() .eq. 2 ) then
        call getarg(2,temp)
        file = temp
      end if
40    open(unit=10, file=file, form='unformatted', status='unknown', err=50)
      goto 60
50    write (*,"('>> Error opening [',a,'] ',/)") file(1:index(file,' '))
      write (*,"(/,'Enter PLOT3D file name [',a,']? ',$)") &
        file(1:index(file,' '))
      read (*,"(a80)") temp
      if (temp(1:1) .ne. ' ') file = temp
      goto 40
60    continue

      write(*,"('Enter the nondimensional time ==> ',$)")
      read(*,*) time

      write(10,err=1010) nx, ny, nz
      write(10,err=1010) Ma, 0.0, Re, time
      write(10,err=1010) (((( real(q(j,i,idof) * exp(-im*two*pi*time)), &
                         i = 1, nx), j = 1, ny), k = 1, nz), idof = 1, ndof)
      close(10,err=1010)

!.... write out a slice

      i = 1
      do j = 1, ny
        write(10,"(6(1e13.6,1x))") x(j,i), abs(q(j,i,1)), abs(q(j,i,2)), &
                                   abs(q(j,i,3)), abs(q(j,i,4)), abs(q(j,i,5))
        write(11,"(6(1e13.6,1x))") x(j,i), real(q(j,i,1)), real(q(j,i,2)), &
                                   real(q(j,i,3)), real(q(j,i,4)), &
                                   real(q(j,i,5))
        write(12,"(6(1e13.6,1x))") x(j,i), aimag(q(j,i,1)), aimag(q(j,i,2)), &
                                   aimag(q(j,i,3)), aimag(q(j,i,4)), &
                                   aimag(q(j,i,5))
      end do
      close(10)
      close(11)
      close(12)

!.... write out a slice

      write(*,"('Enter j ==> ',$)")
      read(*,*) j
      
      do i = 1, nx
        write(16,"(6(1e13.6,1x))") x(j,i), abs(q(j,i,1)), abs(q(j,i,2)), &
                                   abs(q(j,i,3)), abs(q(j,i,4)), abs(q(j,i,5))
        write(17,"(6(1e13.6,1x))") x(j,i), real(q(j,i,1)), real(q(j,i,2)), &
                                   real(q(j,i,3)), real(q(j,i,4)), &
                                   real(q(j,i,5))
        write(18,"(6(1e13.6,1x))") x(j,i), aimag(q(j,i,1)), aimag(q(j,i,2)), &
                                   aimag(q(j,i,3)), aimag(q(j,i,4)), &
                                   aimag(q(j,i,5))
      end do
      close(16)
      close(17)
      close(18)
      do i = 1, nx
        r = sqrt( x(j,i)**2 + y(j,i)**2 )
        th = atan2( y(j,i), x(j,i) )

!       write(13,"(6(1e13.6,1x))") th, abs(q(j,i,1)), abs(q(j,i,2)), &
!                                  abs(q(j,i,3)), abs(q(j,i,4)), abs(q(j,i,5))

        write(13,"(6(1e13.6,1x))") r*cos(th)*abs(q(j,i,4))**2, &
                                   r*sin(th)*abs(q(j,i,4))**2
      end do
      close(13)

!.... write out a slice at a particular time

      write(*,"('Enter time ==> ',$)")
      read(*,*) tl
      
      q(:,:,:) = q(:,:,:) * exp(-im*omega*tl)
      
      i = 1
      do j = 1, ny
        r = sqrt( x(j,i)**2 + y(j,i)**2 )
        write(14,"(6(1e13.6,1x))") r, real(q(j,i,1)), &
                                   real(q(j,i,2)), real(q(j,i,3)), &
                                   real(q(j,i,4)), real(q(j,i,5))
        write(15,"(6(1e13.6,1x))") r, aimag(q(j,i,1)), &
                                   aimag(q(j,i,2)), aimag(q(j,i,3)), &
                                   aimag(q(j,i,4)), aimag(q(j,i,5))
      end do
      close(14)
      close(15)

      call exit(0)      

1000  write(*,"('>> Error reading mean field:  mean.dat')")
      call exit(1)

1010  write(*,"('>> Error writing PLOT3D file')")
      call exit(1)

      end program csubwave
