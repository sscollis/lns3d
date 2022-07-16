!=============================================================================!
module lpost3d

      integer :: nbs, kord=5
      real, allocatable :: knot(:)
      real, allocatable :: bs(:)

      real, external :: BSDER

contains

      function fun(x)
        real :: x, fun
        fun = BSDER( 0, x, kord, knot, nbs, bs )
      end function fun

end module lpost3d

!=============================================================================!
        program main
!
!  lpost3d:  Post processor for three-dimensional linear calculations
!
!=============================================================================!
        use const
        use lpost3d
        use fmax
        implicit none

!.... flow data

        complex, allocatable :: v(:,:,:), q(:,:,:), base(:,:,:), evec(:,:)
        real, allocatable :: u(:,:,:,:), xyz(:,:,:,:), ynn(:)

!.... mean flow data

        real, allocatable :: vm(:,:,:), umb(:), vmb(:), wmb(:)
        integer :: nxm, nym, nzm, ndofm, lstepm
        real :: timem, Rem, Mam, Prm, gammam, cvm

!.... mesh

        real, allocatable :: x(:,:), y(:,:), xi(:), eta(:), z(:), s(:)
        real :: dxi, deta, dz
        real :: wn1, wn2

!.... metrics

        real, allocatable :: m1(:,:),  m2(:,:),  n1(:,:),  n2(:,:),     &
                             m11(:,:), m12(:,:), m22(:,:),              &
                             n11(:,:), n12(:,:), n22(:,:)
        real :: m1l, m2l

        real, allocatable :: bn1(:), bn2(:)
        complex, allocatable :: u1(:), u2(:)

!.... parameters

        real    :: Ma, Re, Pr, gamma, cv, Rgas, time
        integer :: lstep, nx, ny, nz, ndof, Nkz

        integer :: i, j, k, idof

        character(80) file, temp
        integer :: iloc, iend
#ifndef __GFORTRAN__
        integer, external :: iargc
#endif
        integer :: nx2, ny2, nz2

        integer :: itmp
        real :: tmp, arc, amp
        complex :: p

        real :: b, beta, Lz, kz

        logical :: xper=.false., yper=.false.

        integer :: npi
        real, allocatable :: rlog(:), ilog(:)

        complex :: scale = (0.0,0.0)
        real :: smin, smax, savg

        real :: dke, dke1, dke2, dke3
        real, allocatable :: dkeint(:), dkei1(:), dkei2(:), dkei3(:), fint(:)
        complex :: ke, te, mf1
        complex,allocatable :: teint(:)

        real, allocatable :: umax(:), vmax(:), wmax(:), pro(:)
        real :: urmax, uimax
        complex, allocatable :: ucmax(:)
        real :: yumax, yvmax, ywmax, errest

!.... flags

        logical :: threed=.false., debug=.false., rotate=.false.
        logical :: sub=.false., growth=.false., lfilter=.false.
        logical :: add=.false., ltime=.false., switch_ij=.true.
        logical :: calcp=.false.
        integer, parameter :: mfile = 5
        integer :: body=0, emax=0
        integer :: narg, iarg, nfile=0, ifile(mfile)
        character(80) :: arg
        logical :: oreal=.true.
!=============================================================================!
!.... parse the argument list

        narg = iargc()
        do iarg = 1, narg
          call getarg(iarg,arg)
          if (arg(1:1) .ne. '-') then
            nfile = nfile + 1
            if (nfile .gt. mfile) then
              write(*,*) '>> Error in argument list, too many file names'
              call exit(1)
            end if
            ifile(nfile) = iarg
          else
            select case (arg(1:3))
            case ('-3 ')                ! three-dimensional output
              threed = .true.
            case ('-d ')                ! debug
              debug = .true.
            case ('-g ')                ! compute growth rate
              growth = .true.
            case ('-m ')                ! compute efun maxima
              emax = 1
            case ('-r ')                ! rotate to mesh coordinates
              rotate = .true.
            case ('-i ')                ! ouput imaginary component
              oreal = .false.
            case ('-s ')                ! Subtract a base solution
              sub = .true.
            case ('-t ')                ! Write solution at a certain time
              ltime = .true.
            case ('-a ')                ! add a base solution
              add = .true.
            case ('-f ')                ! Filter the field
              lfilter = .true.
            case ('-p ')                ! assume it a parabolic cylinder
              body = 1
            case ('-c ')                ! assume it's a circular cylinder
              body = 2
            case ('-ji')
              switch_ij = .false.
            case('-cp')
              calcp = .true.
            case ('-h ')
              write(*,"('-----------------------------------------------')")
              write(*,"('Usage:  lpost3d [options] [file1]')")
              write(*,"('-----------------------------------------------')")
              write(*,"('    -h:  this help')")
              write(*,"('    -d:  debug')")
              write(*,"('    -g:  Compute growth rate')")
              write(*,"('    -3:  3-D output')")
              write(*,"('    -i:  output imaginary component')")
              write(*,"('    -m:  efun maxima')")
              write(*,"('    -r:  rotate to mesh coordinates')")
              write(*,"('    -s:  subtract base.dat')")
              write(*,"('    -t:  write solution at a certain time')")
              write(*,"('    -a:  add mean.dat')")
              write(*,"('    -f:  filter the field')")
              write(*,"('    -p:  assume a parabolic cylinder')")
              write(*,"('    -c:  assume a circular cylinder')")
              write(*,"('   -cp:  calculate pressure in 4th slot')")
              write(*,"('   -ji:  read metric and restart in ij format')")
              write(*,"('-----------------------------------------------')")
              call exit(0)
            case default
              write(*,"('Argument ',i2,' ignored.')") iarg
            end select
          end if
        end do

!.... read the restart file

        if ( nfile .gt. 0 ) then
          call getarg(ifile(1),temp)
          file = temp
        else
          file = 'output.R.1'
          write (*,"('Enter file name [',a,']? ',$)") file(1:index(file,' '))
          read (*,"(a80)") temp
          if (temp(1:1) .ne. ' ') file = temp
        end if
 10     open(unit=10, file=file, form='unformatted', status='old', err=20)
        goto 30
 20     write (*,"('>> Error opening [',a,'] ',/)") file(1:index(file,' '))
        write (*,"('Enter file name [',a,']? ',$)") file(1:index(file,' '))
        read (*,"(a80)") temp
        if (temp(1:1) .ne. ' ') file = temp
        goto 10
 30     continue

        read(10) lstep, time, nx, ny, nz, ndof, Re, Ma, Pr, gamma, cv
        if (nz .ne. 1) then
          write(*,"(' Error:  nz <> 1 ')")
          call exit(1)
        end if
        allocate( v(ny,nx,ndof) )
        if (switch_ij) then
          read(10) (((v(j,i,idof), idof=1,ndof), i=1,nx), j=1,ny)
        else
          read(10) v
        end if
        close(10)
        write(*,"('Read disturbance field for ',a)") file(1:index(file,' '))
        write(*,"('  Time = ',1pe10.3,'  step = ',i6)") time, lstep

        if (sub) then
 11       open(unit=10,file='base.dat',form='unformatted',status='old',err=21)
          goto 31
 21       write (*,"('>> Error opening [',a,'] ',/)") 'base.dat'
          call exit(1)
 31       continue
          read(10) itmp, tmp, itmp, itmp, itmp, itmp, tmp, tmp, tmp, tmp, tmp
          allocate( base(ny,nx,ndof) )
          if (switch_ij) then
            read(10) (((base(j,i,idof), idof=1,ndof), i=1,nx), j=1,ny)
          else
            read(10) base
          end if
          close(10)
          v = v - base
          deallocate(base)
        end if

!.... read in the grid file

        allocate( x(ny,nx), y(ny,nx), xi(nx), eta(ny), s(nx) )
        open(unit=10,file='grid.dat',form='unformatted',status='old')
        read(10) nx2, ny2, nz2
        if (nx2.ne.nx .or. ny2.ne.ny .or. nz2.ne.nz) then
          write(*,*) 'Grid and data file dimensions do not match'
          call exit(1)
        end if
        read(10) (((x(j,i), i = 1, nx), j = 1, ny), k = 1, nz), &
                 (((y(j,i), i = 1, nx), j = 1, ny), k = 1, nz), &
                 (((   tmp, i = 1, nx), j = 1, ny), k = 1, nz)
        close(10)

!.... allocate storage for metrics

        allocate (m1(ny,nx),  m2(ny,nx),  n1(ny,nx),  n2(ny,nx), &
                  m11(ny,nx), m12(ny,nx), m22(ny,nx),            &
                  n11(ny,nx), n12(ny,nx), n22(ny,nx) )

!.... read in the metric file

        open (unit=10,file='metric.dat',form='unformatted', status='old')
        if (switch_ij) then
          read(10) (( m1(j,i), i=1,nx),j=1,ny), &
                   (( m2(j,i), i=1,nx),j=1,ny), &
                   (( n1(j,i), i=1,nx),j=1,ny), &
                   (( n2(j,i), i=1,nx),j=1,ny), &
                   ((m11(j,i), i=1,nx),j=1,ny), &
                   ((m12(j,i), i=1,nx),j=1,ny), &
                   ((m22(j,i), i=1,nx),j=1,ny), &
                   ((n11(j,i), i=1,nx),j=1,ny), &
                   ((n12(j,i), i=1,nx),j=1,ny), &
                   ((n22(j,i), i=1,nx),j=1,ny)
        else
          read(10) m1, m2, n1, n2, m11, m12, m22, n11, n12, n22
        end if
        close(10)

!.... read in the mean flow file

        open(unit=10, file='mean.dat', form='unformatted', status='old')
        read(10) lstepm, timem, nxm, nym, nzm, ndofm, &
                 Rem, Mam, Prm, gammam, cvm
        if (nxm.ne.nx .or. nym.ne.ny .or. nzm.ne.nz .or. ndofm.ne.ndof) then
          write(*,*) 'Mean and disturbance file dimensions do not match'
          call exit(1)
        end if
        allocate( vm(ny,nx,ndof) )
        if (switch_ij) then
          read(10) (((vm(j,i,idof), idof=1,ndof), i=1,nx), j=1,ny)
        else
          read(10) vm
        end if
        close(10)

!.... read in the body discription file

        if (body .eq. 0) then
          open(unit=10,file='body.dat',form='formatted',status='old',err=100)
          do i = 1, nx
            read(10,*,err=100) tmp, s(i)
          end do
          close(10)
          goto 101
        end if
100     continue

        if (body.eq.1) then        ! parabolic cylinder
          if (x(1,1).lt.zero) x(1,1) = zero
          s(:) = Sqrt(x(1,:) + two*x(1,:)**2)/Sqrt(two) + &
                 Log(one + four*x(1,:) + two**(onept5)*Sqrt(x(1,:) + &
                 two*x(1,:)**2)) / four
        else if (body.eq.2) then   ! circular cylinder
          s(:) = pi - atan2( y(1,:), x(1,:) )
        else
          s(:) = x(1,:)
          write(*,*) 'WARNING:  Setting s = x'
        end if
 101    continue

!.... make the xi grid

        dxi = one / float(nx-1)
        do i = 1, nx
          xi(i) = zero + (i-1) * dxi
        end do

!.... make the eta grid

        deta = one / float(ny-1)
        do j = 1, ny
          eta(j) = float(j-1) * deta
        end do

!.... filter the field?

        if (lfilter) then
          write(*,*) 'Filtering in Y'
          call cfilter( ny, nx, v(:,:,1) )
          call cfilter( ny, nx, v(:,:,2) )
          call cfilter( ny, nx, v(:,:,3) )
          call cfilter( ny, nx, v(:,:,4) )
          call cfilter( ny, nx, v(:,:,5) )
        end if

!.... form the output variables

        allocate( q(ny,nx,ndof) )

        q(:,:,1) = v(:,:,1)
        q(:,:,2) = v(:,:,2)
        q(:,:,3) = v(:,:,3)
        if (threed) then
          q(:,:,4) = v(:,:,4)
        else if (calcp) then                        ! pressure
          q(:,:,4) = ( v(:,:,1) * vm(:,:,5) + &
                      vm(:,:,1) *  v(:,:,5) ) / (gamma * Ma**2)
        else
          q(:,:,4) = v(:,:,4)
        end if
        q(:,:,5) = v(:,:,5)

!.... amp.dat file

        j = 1
        open(unit=98,file='amp.dat')
        do i = 1, nx
          p = ( v(j,i,1) * vm(j,i,5) + &
               vm(j,i,1) *  v(j,i,5) ) / (gamma * Ma**2)
          write(98,*) x(j,i), abs(2.0*p)
        end do
        close(98)

!.... rotate to the conformal coordinate system

        if (rotate) then
          allocate( bn1(ny), bn2(ny), u1(ny), u2(ny) )
          do i = 1, nx
            bn1      = n1(:,i) / sqrt( n1(:,i)**2 + n2(:,i)**2 )
            bn2      = n2(:,i) / sqrt( n1(:,i)**2 + n2(:,i)**2 )
            u1       = bn2 * q(:,i,2) - bn1 * q(:,i,3)
            u2       = bn1 * q(:,i,2) + bn2 * q(:,i,3)
            q(:,i,2) = u1
            q(:,i,3) = u2
          end do
          deallocate( bn1, bn2, u1, u2 )
        end if

!.... form the 2D field

        if (ltime) then
          write(*,"('Enter the nondimensional time ==> ',$)")
          read(*,*) time
        else
          time = zero
        end if

        allocate( u(nx,ny,1,ndof) )

        if (oreal) then
          do idof = 1, ndof
            do k = 1, 1
              do i = 1, nx
                do j = 1, ny
                  u(i,j,k,idof) = real( q(j,i,idof) * exp(-im*two*pi*time) )
                end do
              end do
            end do
          end do
        else
          do idof = 1, ndof
            do k = 1, 1
              do i = 1, nx
                do j = 1, ny
                  u(i,j,k,idof) = imag( q(j,i,idof) * exp(-im*two*pi*time) )
                end do
              end do
            end do
          end do
        end if

        iloc = index(file,'.R.')
        iend = index(file,' ')-1
        if (iloc .ne. 0) then
          temp = file(1:iloc)//'q'//file(iloc+2:iend)
        else
          temp = file(1:iend)//'.q'
        end if
        if ( nfile .eq. 0 ) then
          file = temp
          write (*,"(/,'Enter PLOT3D file name [',a,']? ',$)")  &
                  file(1:index(file,' '))
          read (*,"(a80)") temp
          if (temp(1:1) .ne. ' ') file = temp
        else if ( nfile .eq. 1 ) then
          file = temp
        else if ( nfile .eq. 2 ) then
          call getarg(ifile(2),temp)
          file = temp
        end if
        open(unit=10, file=file, form='unformatted', status='unknown')
        write(10) nx, ny, 1
        write(10) Ma, Pr, Re, time
        write(10) u
        close(10)
        deallocate( u )

!.... write out the magnitude into a separate (real) LNS file
!.... for interpolation to another mesh

        open (10,file='mag.dat',form='unformatted',status='unknown')
        write(10) lstep, time, nx, ny, nz, ndof, Re, Ma, Pr, gamma, cv
        write(10) abs( q(:,:,:) )
        close(10)

!.... form the 3D field

        if (.true. .and. threed) then

!.... make the z grid

        write(*,"('Enter Nz ==> ',$)")
        read(*,*) nz
        write(*,"('Enter Lz, Nkz (periodic extension) ==> ',$)")
        read(*,*) Lz, Nkz

        if (Lz.eq.zero) then
          kz = zero
        else
          kz = two * pi / Lz
        end if
        dz = Nkz * Lz / float(nz-1)
        allocate( z(nz) )
        do k = 1, nz
          z(k) = float(k-1) * dz
        end do

        allocate( u(nx,ny,nz,ndof) )

        if (add) then
          write(*,"('Enter amp ==> ',$)")
          read(*,*) amp
          do idof = 1, ndof
            do k = 1, nz
              do i = 1, nx
                do j = 1, ny
                  u(i,j,k,idof) = vm(j,i,idof) + &
                                  amp * real( q(j,i,idof) * exp(im*kz*z(k)) )
                end do
              end do
            end do
          end do
        else
          do idof = 1, ndof
            do k = 1, nz
              do i = 1, nx
                do j = 1, ny
                  u(i,j,k,idof) = real( q(j,i,idof) * exp(im*kz*z(k)) )
                end do
              end do
            end do
          end do
        end if

        write(*,*) 'Writing 3D output file'
        open(unit=10, file='out3d.q', form='unformatted', status='unknown')
        write(10) nx, ny, nz
        write(10) Ma, Pr, Re, time
        write(10) u
        close(10)
        deallocate( u)

        allocate( xyz(nx,ny,nz,3) )
        do k = 1, nz
          do i = 1, nx
            do j = 1, ny
              xyz(i,j,k,1) = x(j,i)
              xyz(i,j,k,2) = y(j,i)
              xyz(i,j,k,3) = z(k)
            end do
          end do
        end do

        write(*,*) 'Writing 3D grid file'
        open(unit=10,file='grid3d.xyz',form='unformatted',status='unknown')
        write(10) nx, ny, nz
        write(10) xyz
        close(10)
        deallocate( xyz )

        end if

        if (debug) then

!.... Magnitude at a fixed j off the wall

        write(*,"('Enter j ==> ',$)")
        read(*,*) j
        do i = 1, nx
          write(10,"(6(1pe21.12E4,1x))") s(i),   &
                                        abs(q(j,i,1)), &
                                        abs(q(j,i,2)), &
                                        abs(q(j,i,3)), &
                                        abs(q(j,i,4)), &
                                        abs(q(j,i,5))
        end do
        close(10)
        do i = 1, nx
          write(11,"(6(1pe21.12E4,1x))") s(i),   &
                                        real(q(j,i,1)), &
                                        real(q(j,i,2)), &
                                        real(q(j,i,3)), &
                                        real(q(j,i,4)), &
                                        real(q(j,i,5))
        end do
        close(11)
        do i = 1, nx
          write(12,"(6(1pe21.12E4,1x))") s(i),   &
                                        aimag(q(j,i,1)), &
                                        aimag(q(j,i,2)), &
                                        aimag(q(j,i,3)), &
                                        aimag(q(j,i,4)), &
                                        aimag(q(j,i,5))
        end do
        close(12)

!.... On left boundary

        i = 1
        do j = 1, ny
          write(13,"(6(1pe21.12E4,1x))") x(j,i),   &
                                        real(q(j,i,1)), &
                                        real(q(j,i,2)), &
                                        real(q(j,i,3)), &
                                        real(q(j,i,4)), &
                                        real(q(j,i,5))
        end do
        close(13)
        do j = 1, ny
          write(14,"(6(1pe21.12E4,1x))") x(j,i),   &
                                        aimag(q(j,i,1)), &
                                        aimag(q(j,i,2)), &
                                        aimag(q(j,i,3)), &
                                        aimag(q(j,i,4)), &
                                        aimag(q(j,i,5))
        end do
        close(14)

!.... Output mode shapes at index i normal to the wall

        write(*,"('Enter i ==> ',$)")
        read(*,*) i

!.... Rotate to Body Fixed Coordinates

        allocate( evec(ny,ndof), ynn(ny), umb(ny), vmb(ny), wmb(ny) )

        wn1 = n1(1,i) / sqrt( n1(1,i)**2 + n2(1,i)**2 )
        wn2 = n2(1,i) / sqrt( n1(1,i)**2 + n2(1,i)**2 )
        arc = zero
        do j = 1, ny
          ynn(j) = arc
          if (j.ne.ny) arc = arc + sqrt( (x(j+1,i)-x(j,i))**2 + &
                                         (y(j+1,i)-y(j,i))**2 )
          evec(j,1) = q(j,i,1)
          evec(j,2) = wn2 * q(j,i,2) - wn1 * q(j,i,3)
          evec(j,3) = wn1 * q(j,i,2) + wn2 * q(j,i,3)
          evec(j,4) = q(j,i,4)
          evec(j,5) = q(j,i,5)

          umb(j) = wn2 * vm(j,i,2) - wn1 * vm(j,i,3)
          vmb(j) = wn1 * vm(j,i,2) + wn2 * vm(j,i,3)
          wmb(j) = vm(j,i,4)
        end do

        do j = 1, ny
          write(15,"(6(1pe21.12E4,1x))") ynn(j),   &
                                        abs(q(j,i,1)), &
                                        abs(q(j,i,2)), &
                                        abs(q(j,i,3)), &
                                        abs(q(j,i,4)), &
                                        abs(q(j,i,5))
        end do
        close(15)

        do j = 1, ny
          write(16,"(6(1pe21.12E4,1x))") ynn(j),   &
                                        real(q(j,i,1)), &
                                        real(q(j,i,2)), &
                                        real(q(j,i,3)), &
                                        real(q(j,i,4)), &
                                        real(q(j,i,5))
        end do
        close(16)

        do j = 1, ny
          write(32,"(6(1pe21.12E4,1x))") ynn(j),   &
                                        aimag(q(j,i,1)), &
                                        aimag(q(j,i,2)), &
                                        aimag(q(j,i,3)), &
                                        aimag(q(j,i,4)), &
                                        aimag(q(j,i,5))
        end do
        close(32)

!==============================================================================
!.... output some other important stuff like:
!==============================================================================
!.... Disturbance Kinetic energy
!.... Total energy
!.... Kinetic energy
!.... streamwise mass flux
!==============================================================================

        do j = 1, ny
          dke = vm(j,i,1) * ( abs(evec(j,2))**2 + abs(evec(j,3))**2 + &
                              abs(evec(j,4))**2 )
          ke  = ( umb(j) * evec(j,2) + vmb(j) * evec(j,3) + &
                  wmb(j) * evec(j,4) )
          te  = evec(j,5) / (gamma * (gamma-one) * Ma**2) + ke
          mf1 = vm(j,i,1) * evec(j,2) + umb(j) * evec(j,1)

!         dke = vm(j,i,1) * ( abs(q(j,i,2))**2 + abs(q(j,i,3))**2 + &
!                             abs(q(j,i,4))**2 )
!         ke  = ( vm(j,i,2) * q(j,i,2) + vm(j,i,3) * q(j,i,3) + &
!                 vm(j,i,4) * q(j,i,4) )
!         te  = q(j,i,5) / (gamma * (gamma-one) * Ma**2) + ke
!         mf1 = vm(j,i,1) * q(j,i,2) + vm(j,i,2) * q(j,i,1)

          write(18,"(11(1pe21.12E4,1x))") ynn(j), dke, &
                                          abs(ke), abs(te), abs(mf1)
        end do
        close(18)

!.... normalize the mode shape

        do idof = 1, ndof
          do j = 1, ny
            if ( abs(evec(j,idof)) .gt. abs(scale) ) scale = evec(j,idof)
          end do
        end do
        write(*,*) 'Scale factor = ', scale, abs(scale)
        if (scale .ne. zero) then
          evec = evec / scale
        end if

!.... output a normalized mode-shape

        do j = 1, ny
          write(17,"(11(1pe21.12E4,1x))") ynn(j),   &
                     real( evec(j,1) ), aimag( evec(j,1) ), &
                     real( evec(j,2) ), aimag( evec(j,2) ), &
                     real( evec(j,3) ), aimag( evec(j,3) ), &
                     real( evec(j,4) ), aimag( evec(j,4) ), &
                     real( evec(j,5) ), aimag( evec(j,5) )
        end do
        close(17)

        deallocate( evec, ynn, umb, vmb, wmb )

        end if                  ! debug

!==============================================================================
!.... Compute the growth-rates
!==============================================================================

        smin = -one
        smax = -one
        if (growth) then
          write(*,"('Enter smin, smax ==> ',$)")
          read(*,*) smin, smax
        end if

        allocate( evec(ny,ndof), ynn(ny), umax(nx), &
                  vmax(nx), wmax(nx), pro(ny), ucmax(nx) )

!.... compute the local max of (u,v,w) in body-normal coordinates

        if (emax.eq.1) then

        open(unit=50,file='emax.dat',status='unknown',form='formatted')
        umax  = one
        vmax  = one
        wmax  = one
        ucmax = one

        do i = 1, nx
          if ( s(i) .ge. smin .and. s(i) .le. smax ) then
            wn1 = n1(1,i) / sqrt( n1(1,i)**2 + n2(1,i)**2 )
            wn2 = n2(1,i) / sqrt( n1(1,i)**2 + n2(1,i)**2 )
            arc = zero
            do j = 1, ny
              ynn(j) = arc
              if (j.ne.ny) arc = arc + sqrt( (x(j+1,i)-x(j,i))**2 + &
                                             (y(j+1,i)-y(j,i))**2 )
              evec(j,1) = q(j,i,1)
              evec(j,2) = wn2 * q(j,i,2) - wn1 * q(j,i,3)
              evec(j,3) = wn1 * q(j,i,2) + wn2 * q(j,i,3)
              evec(j,4) = q(j,i,4)
              evec(j,5) = q(j,i,5)
            end do

!           pro = abs( vm(:,i,1) * evec(:,2) )
            pro = abs( evec(:,2) )
            umax(i) = findmax( ny, ynn, pro, yumax )

            pro = real( evec(:,2) )
            urmax = getval( ny, ynn, pro, yumax )
            pro = aimag( evec(:,2) )
            uimax = getval( ny, ynn, pro, yumax )
            ucmax(i) = cmplx(urmax,uimax)

!           pro = abs( vm(:,i,1) * evec(:,3) )
            pro = abs( evec(:,3) )
!           vmax(i) = findmax( ny, ynn, pro, yvmax )
            yvmax = yumax
            vmax(i) = getval( ny, ynn, pro, yvmax )

!           pro = abs( vm(:,i,1) * evec(:,4) )
            pro = abs( evec(:,4) )
!           wmax(i) = findmax( ny, ynn, pro, ywmax )
            ywmax = yumax
            wmax(i) = getval( ny, ynn, pro, ywmax )

            write(50,"(8(1pe21.12E4,1x))") s(i), umax(i), vmax(i), wmax(i), &
                                           yumax, yvmax, ywmax
          end if
        end do
        close(50)

!.... now compute the growth rates

        open(unit=51,file='emax-grw.dat',status='unknown',form='formatted')
        do i = 2, nx
          if ( s(i) .ge. smin .and. s(i) .le. smax ) then
            savg = (s(i)+s(i-1))*pt5
            write(51,"(8(1pe21.12E4,1x))") (s(i)+s(i-1))*pt5,       &
                -(log(umax(i))-log(umax(i-1)))/(s(i)-s(i-1)),       &
                -(log(vmax(i))-log(vmax(i-1)))/(s(i)-s(i-1)),       &
                -(log(wmax(i))-log(wmax(i-1)))/(s(i)-s(i-1))
          end if
        end do
        close(51)

        open(unit=52,file='ucmax-grw.dat',status='unknown',form='formatted')
        npi = 0
        allocate( rlog(nx), ilog(nx) )
        do i = 2, nx
          savg = (s(i)+s(i-1))*pt5
          if ( savg .gt. smin .and. savg .lt. smax ) then
            if (aimag(log(ucmax(i))) .gt. aimag(log(ucmax(i-1))) ) then
              npi = npi + 2
            end if
            rlog(i) = real(log(ucmax(i)))
            ilog(i) = aimag(log(ucmax(i))) - npi * pi
            write(52,"(6(1pe21.12E4,1x))") (s(i)+s(i-1))*pt5,   &
                    (ilog(i)-ilog(i-1))/(s(i)-s(i-1)), &
                   -(rlog(i)-rlog(i-1))/(s(i)-s(i-1)), &
                    (rlog(i)+rlog(i-1))*pt5, &
                    (ilog(i)+ilog(i-1))*pt5
          end if
        end do
        close(52)
        deallocate( rlog, ilog )

        deallocate( umax, vmax, wmax, pro, ucmax )

        end if

!.... compute the disturbance kinetic energy integral

        allocate( dkeint(nx), dkei1(nx), dkei2(nx), dkei3(nx) )
        dkeint = zero
        dkei1  = zero
        dkei2  = zero
        dkei3  = zero

        nbs = ny
        allocate( knot(nbs+kord), bs(nbs), fint(nbs) )

        do i = 1, nx

!.... transform to body-normal coordinate system

          wn1 = n1(1,i) / sqrt( n1(1,i)**2 + n2(1,i)**2 )
          wn2 = n2(1,i) / sqrt( n1(1,i)**2 + n2(1,i)**2 )
          arc = zero
          do j = 1, ny
            ynn(j) = arc
            if (j.ne.ny) arc = arc + sqrt( (x(j+1,i)-x(j,i))**2 + &
                                           (y(j+1,i)-y(j,i))**2 )
            evec(j,1) = q(j,i,1)
            evec(j,2) = wn2 * q(j,i,2) - wn1 * q(j,i,3)
            evec(j,3) = wn1 * q(j,i,2) + wn2 * q(j,i,3)
            evec(j,4) = q(j,i,4)
            evec(j,5) = q(j,i,5)
          end do

!.... Integrate in physical space with B-splines

          if (.false.) then
          call BSNAK( nbs, ynn, kord, knot)

          fint(:) = (abs(evec(:,2))**2+abs(evec(:,3))**2+abs(evec(:,4))**2)
          call BSINT( nbs, ynn, fint, kord, knot, bs )
          call QDAG( fun, zero, ynn(ny), 1.0e-8, 1.0e-8, 2, dkeint(i), errest )

          fint(:) = abs(evec(:,2))**2
          call BSINT( nbs, ynn, fint, kord, knot, bs )
          call QDAG( fun, zero, ynn(ny), 1.0e-8, 1.0e-8, 2, dkei1(i), errest )

          fint(:) = abs(evec(:,3))**2
          call BSINT( nbs, ynn, fint, kord, knot, bs )
          call QDAG( fun, zero, ynn(ny), 1.0e-8, 1.0e-8, 2, dkei2(i), errest )

          fint(:) = abs(evec(:,4))**2
          call BSINT( nbs, ynn, fint, kord, knot, bs )
          call QDAG( fun, zero, ynn(ny), 1.0e-8, 1.0e-8, 2, dkei3(i), errest )
          end if

!.... integral in physical space using trapezoid

          open(unit=20,file='dke.dat',status='unknown',form='formatted')
          if (.true.) then
            do j = 1, ny-1
              dke1 = abs(evec(j,2))**2
              dke2 = abs(evec(j,3))**2
              dke3 = abs(evec(j,4))**2
              dke  = ( abs(evec(j,2))**2 + abs(evec(j,3))**2 + &
                       abs(evec(j,4))**2 )
              dkeint(i) = dkeint(i) + pt5 * (ynn(j+1)-ynn(j)) * dke
              dkei1(i)  = dkei1(i)  + pt5 * (ynn(j+1)-ynn(j)) * dke1
              dkei2(i)  = dkei2(i)  + pt5 * (ynn(j+1)-ynn(j)) * dke2
              dkei3(i)  = dkei3(i)  + pt5 * (ynn(j+1)-ynn(j)) * dke3
              dke1 = abs(evec(j+1,2))**2
              dke2 = abs(evec(j+1,3))**2
              dke3 = abs(evec(j+1,4))**2
              dke  = ( abs(evec(j+1,2))**2 + abs(evec(j+1,3))**2 + &
                       abs(evec(j+1,4))**2 )
              dkeint(i) = dkeint(i) + pt5 * (ynn(j+1)-ynn(j)) * dke
              dkei1(i)  = dkei1(i)  + pt5 * (ynn(j+1)-ynn(j)) * dke1
              dkei2(i)  = dkei2(i)  + pt5 * (ynn(j+1)-ynn(j)) * dke2
              dkei3(i)  = dkei3(i)  + pt5 * (ynn(j+1)-ynn(j)) * dke3
            end do
          end if
          write(20,"(6(1pe21.12E4,1x))") s(i), sqrt(dkeint(i)), &
                          sqrt(dkei1(i)), sqrt(dkei2(i)), sqrt(dkei3(i))
        end do
        close(20)
        deallocate( evec, ynn )
        deallocate( knot, bs, fint )

        if (growth) then
          open(unit=21,file='dke-grw.dat',status='unknown',form='formatted')
          do i = 2, nx
            savg = (s(i)+s(i-1))*pt5
            if ( savg .gt. smin .and. savg .lt. smax ) then
              write(21,"(6(1pe21.12E4,1x))") (s(i)+s(i-1))*pt5,       &
                -pt5*(log(dkeint(i))-log(dkeint(i-1)))/(s(i)-s(i-1)), &
                -pt5*(log(dkei1(i))-log(dkei1(i-1)))/(s(i)-s(i-1)),   &
                -pt5*(log(dkei2(i))-log(dkei2(i-1)))/(s(i)-s(i-1)),   &
                -pt5*(log(dkei3(i))-log(dkei3(i-1)))/(s(i)-s(i-1))
            end if
          end do
          close(21)
        end if

!.... compute the total energy integral

        allocate( teint(nx) )
        open(unit=22,file='te.dat',status='unknown',form='formatted')
        do i = 1, nx
          j = 1
          ke  = ( vm(j,i,2) * q(j,i,2) + vm(j,i,3) * q(j,i,3) + &
                  vm(j,i,4) * q(j,i,4) )
          te  = q(j,i,5) / (gamma * (gamma-one) * Ma**2) + ke
          teint(i) = pt5 * deta * te
          do j = 2, ny-1
            ke  = ( vm(j,i,2) * q(j,i,2) + vm(j,i,3) * q(j,i,3) + &
                    vm(j,i,4) * q(j,i,4) )
            te  = q(j,i,5) / (gamma * (gamma-one) * Ma**2) + ke
            teint(i) = teint(i) + deta * te
          end do
          j = ny
          ke  = ( vm(j,i,2) * q(j,i,2) + vm(j,i,3) * q(j,i,3) + &
                  vm(j,i,4) * q(j,i,4) )
          te  = q(j,i,5) / (gamma * (gamma-one) * Ma**2) + ke
          teint(i) = teint(i) + pt5 * deta * te

          write(22,"(6(1pe21.12E4,1x))") s(i), abs(teint(i)), &
                                         real(teint(i)), aimag(teint(i))
        end do
        close(22)

        if (growth) then
          open(unit=23,file='te-grw.dat',status='unknown',form='formatted')
          npi = 0
          allocate( rlog(nx), ilog(nx) )
          do i = 2, nx
            savg = (s(i)+s(i-1))*pt5
            if ( savg .gt. smin .and. savg .lt. smax ) then
              if (aimag(log(teint(i))) .gt. aimag(log(teint(i-1))) ) then
                npi = npi + 2
              end if
              rlog(i) = real(log(teint(i)))
              ilog(i) = aimag(log(teint(i))) - npi * pi
              write(23,"(6(1pe21.12E4,1x))") (s(i)+s(i-1))*pt5,   &
                      (ilog(i)-ilog(i-1))/(s(i)-s(i-1)), &
                     -(rlog(i)-rlog(i-1))/(s(i)-s(i-1)), &
                      (rlog(i)+rlog(i-1))*pt5, &
                      (ilog(i)+ilog(i-1))*pt5
            end if
          end do
          close(23)
          deallocate( rlog, ilog )
        end if

!.... compute the kinetic energy integral

        open(unit=24,file='ke.dat',status='unknown',form='formatted')
        do i = 1, nx
          j = 1
          ke  = ( vm(j,i,2) * q(j,i,2) + vm(j,i,3) * q(j,i,3) + &
                  vm(j,i,4) * q(j,i,4) )
          te  = ke
          teint(i) = pt5 * deta * te
          do j = 2, ny-1
            ke  = ( vm(j,i,2) * q(j,i,2) + vm(j,i,3) * q(j,i,3) + &
                    vm(j,i,4) * q(j,i,4) )
            te  = ke
            teint(i) = teint(i) + deta * te
          end do
          j = ny
          ke  = ( vm(j,i,2) * q(j,i,2) + vm(j,i,3) * q(j,i,3) + &
                  vm(j,i,4) * q(j,i,4) )
          te  = ke
          teint(i) = teint(i) + pt5 * deta * te

          write(24,"(6(1pe21.12E4,1x))") s(i), abs(teint(i)), &
                                         real(teint(i)), aimag(teint(i))
        end do
        close(24)

        if (growth) then
          open(unit=25,file='ke-grw.dat',status='unknown',form='formatted')
          npi = 0
          allocate( rlog(nx), ilog(nx) )
          do i = 2, nx
            savg = (s(i)+s(i-1))*pt5
            if ( savg .gt. smin .and. savg .lt. smax ) then
              if (aimag(log(teint(i))) .gt. aimag(log(teint(i-1))) ) then
                npi = npi + 2
              end if
              rlog(i) = real(log(teint(i)))
              ilog(i) = aimag(log(teint(i))) - npi * pi
              write(25,"(6(1pe21.12E4,1x))") (s(i)+s(i-1))*pt5,   &
                      (ilog(i)-ilog(i-1))/(s(i)-s(i-1)), &
                     -(rlog(i)-rlog(i-1))/(s(i)-s(i-1)), &
                      (rlog(i)+rlog(i-1))*pt5, &
                      (ilog(i)+ilog(i-1))*pt5
            end if
          end do
          close(25)
          deallocate( rlog, ilog )
        end if

!.... compute the u energy integral

        open(unit=26,file='ue.dat',status='unknown',form='formatted')
        do i = 1, nx
          j = 1
          ke  = ( vm(j,i,2) * q(j,i,2) )
          te  = ke
          teint(i) = pt5 * deta * te
          do j = 2, ny-1
            ke  = ( vm(j,i,2) * q(j,i,2) )
            te  = ke
            teint(i) = teint(i) + deta * te
          end do
          j = ny
          ke  = ( vm(j,i,2) * q(j,i,2) )
          te  = ke
          teint(i) = teint(i) + pt5 * deta * te
          write(26,"(6(1pe21.12E4,1x))") s(i), abs(teint(i)), &
                                         real(teint(i)), aimag(teint(i))
        end do
        close(26)

        if (growth) then
          open(unit=27,file='ue-grw.dat',status='unknown',form='formatted')
          npi = 0
          allocate( rlog(nx), ilog(nx) )
          do i = 2, nx
            savg = (s(i)+s(i-1))*pt5
            if ( savg .gt. smin .and. savg .lt. smax ) then
              if (aimag(log(teint(i))) .gt. aimag(log(teint(i-1))) ) then
                npi = npi + 2
              end if
              rlog(i) = real(log(teint(i)))
              ilog(i) = aimag(log(teint(i))) - npi * pi
              write(27,"(6(1pe21.12E4,1x))") (s(i)+s(i-1))*pt5,   &
                      (ilog(i)-ilog(i-1))/(s(i)-s(i-1)), &
                     -(rlog(i)-rlog(i-1))/(s(i)-s(i-1)), &
                      (rlog(i)+rlog(i-1))*pt5, &
                      (ilog(i)+ilog(i-1))*pt5
            end if
          end do
          close(27)
          deallocate( rlog, ilog )
        end if

!.... compute the v energy integral

#ifdef ADDITIONAL_VELOCITY_INTEGRALS

        open(unit=28,file='ve.dat',status='unknown',form='formatted')
        do i = 1, nx
          j = 1
          ke  = ( vm(j,i,3) * q(j,i,3) )
          te  = ke
          teint(i) = pt5 * deta * te
          do j = 2, ny-1
            ke  = ( vm(j,i,3) * q(j,i,3) )
            te  = ke
            teint(i) = teint(i) + deta * te
          end do
          j = ny
          ke  = ( vm(j,i,3) * q(j,i,3) )
          te  = ke
          teint(i) = teint(i) + pt5 * deta * te

          write(28,"(6(1pe21.12E4,1x))") s(i), abs(teint(i)), &
                                         real(teint(i)), aimag(teint(i))
        end do
        close(28)

        if (growth) then
          open(unit=29,file='ve-grw.dat',status='unknown',form='formatted')
          npi = 0
          allocate( rlog(nx), ilog(nx) )
          do i = 2, nx
            savg = (s(i)+s(i-1))*pt5
            if ( savg .gt. smin .and. savg .lt. smax ) then
              if (aimag(log(teint(i))) .gt. aimag(log(teint(i-1))) ) then
                npi = npi + 2
              end if
              rlog(i) = real(log(teint(i)))
              ilog(i) = aimag(log(teint(i))) - npi * pi
              write(29,"(6(1pe21.12E4,1x))") (s(i)+s(i-1))*pt5,   &
                      (ilog(i)-ilog(i-1))/(s(i)-s(i-1)), &
                     -(rlog(i)-rlog(i-1))/(s(i)-s(i-1)), &
                      (rlog(i)+rlog(i-1))*pt5, &
                      (ilog(i)+ilog(i-1))*pt5
            end if
          end do
          close(29)
          deallocate( rlog, ilog )
        end if

!.... compute the w energy integral

        open(unit=30,file='we.dat',status='unknown',form='formatted')
        do i = 1, nx
          j = 1
          ke  = ( vm(j,i,4) * q(j,i,4) )
          te  = ke
          teint(i) = pt5 * deta * te
          do j = 2, ny-1
            ke  = ( vm(j,i,4) * q(j,i,4) )
            te  = ke
            teint(i) = teint(i) + deta * te
          end do
          j = ny
          ke  = ( vm(j,i,4) * q(j,i,4) )
          te  = ke
          teint(i) = teint(i) + pt5 * deta * te

          write(30,"(6(1pe21.12E4,1x))") s(i), abs(teint(i)), &
                                         real(teint(i)), aimag(teint(i))
        end do
        close(30)

        if (growth) then
          open(unit=31,file='we-grw.dat',status='unknown',form='formatted')
          npi = 0
          allocate( rlog(nx), ilog(nx) )
          do i = 2, nx
            savg = (s(i)+s(i-1))*pt5
            if ( savg .gt. smin .and. savg .lt. smax ) then
              if (aimag(log(teint(i))) .gt. aimag(log(teint(i-1))) ) then
                npi = npi + 2
              end if
              rlog(i) = real(log(teint(i)))
              ilog(i) = aimag(log(teint(i))) - npi * pi
              write(31,"(6(1pe21.12E4,1x))") (s(i)+s(i-1))*pt5,   &
                      (ilog(i)-ilog(i-1))/(s(i)-s(i-1)), &
                     -(rlog(i)-rlog(i-1))/(s(i)-s(i-1)), &
                      (rlog(i)+rlog(i-1))*pt5, &
                      (ilog(i)+ilog(i-1))*pt5
            end if
          end do
          close(31)
          deallocate( rlog, ilog )
        end if

#endif 

        call exit(0)
        end
