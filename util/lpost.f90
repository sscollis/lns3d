!=============================================================================!
        program lpost 
!  
!  Post processor for linear calculations
!  
!=============================================================================!
        use const
        implicit none

!.... flow data

        real, allocatable :: v(:,:,:), q(:,:,:)
        
        real, allocatable :: g1v(:,:,:), g2v(:,:,:)
        real, allocatable :: g11v(:,:,:), g12v(:,:,:), g22v(:,:,:)
        
        real, allocatable :: g1vl(:), g2vl(:), g11vl(:), g12vl(:), g22vl(:)

!.... mean flow data

        real, allocatable :: vm(:,:,:)
        integer :: nxm, nym, nzm, ndofm, lstepm
        real :: timem, Rem, Mam, Prm, gammam, cvm
        
!.... mesh

        real, allocatable :: x(:,:), y(:,:), xi(:), eta(:)
        real :: dxi, deta

!.... metrics

        real, allocatable :: m1(:,:),  m2(:,:),  n1(:,:),  n2(:,:),     &
                             m11(:,:), m12(:,:), m22(:,:),              &
                             n11(:,:), n12(:,:), n22(:,:)
        real, allocatable :: bn1(:), bn2(:)

        real :: m1l, m2l

!.... parameters

        real    :: Ma, Re, Pr, gamma, cv, Rgas, time
        integer :: lstep, nx, ny, nz, ndof

        integer :: i, j, k, idof

        character(80) file, temp
        integer :: iloc, iend
#ifndef __GFORTRAN__    
        integer, external :: iargc
#endif

        integer :: nx2, ny2, nz2

        real :: tmp
        
        real :: b, beta

        logical :: xper=.false., yper=.false.
        
!.... temporary variables

        real :: rhom, tm, cm, c1, c2, c3, c4, rho, u1, u2, p, t

!.... flags

        logical :: sub=.false., deriv=.false.
        integer, parameter :: mfile = 5
        integer :: narg, iarg, nfile=0, ifile(mfile)
        character(80) :: arg

        logical :: switch_ij = .true.
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
            select case (arg(1:2))
            case ('-s')         ! Subtract a gaussian pulse
              sub = .true.
            case ('-d')         ! Output derivative file
              deriv = .true.
            case ('-h')
              write(*,"('-----------------------------------------------')")
              write(*,"('Usage:  lpost [options] [file1]')")
              write(*,"('-----------------------------------------------')")
              write(*,"('   -h:  this help')")
              write(*,"('   -s:  subtract Gaussian pulse')")
              write(*,"('   -d:  derivative Plot3d file')")
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
        
        read(10) lstep, time, nx, ny, nz, ndof, &
                 Re, Ma, Pr, gamma, cv
        if (nz .ne. 1) then
          write(*,"(' Error:  nz <> 1 ')")
          call exit(1)
        end if
        allocate( v(ny,nx,ndof) )
        if (switch_ij) then
          read(10) (((v(j,i,k), k=1,ndof), i=1,nx), j=1,ny)
        else
          read(10) v
        endif
        close(10)
        write(*,"('Read disturbance field for ',a)") file(1:index(file,' '))
        write(*,"('  Time = ',1pe10.3,'  step = ',i6)") time, lstep

!.... read in the grid file

        allocate( x(ny,nx), y(ny,nx), xi(nx), eta(ny) )
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
          read(10) (((vm(j,i,k), k=1,ndof), i=1,nx), j=1,ny)
        else
          read(10) vm
        endif
        close(10)

!.... allocate storage for metrics

        allocate (m1(ny,nx),  m2(ny,nx),  n1(ny,nx),  n2(ny,nx), &
                  m11(ny,nx), m12(ny,nx), m22(ny,nx),            &
                  n11(ny,nx), n12(ny,nx), n22(ny,nx) )

!.... read in the metric file

        open (unit=10,file='metric.dat',form='unformatted', status='old')
        if (switch_ij) then
          read(10) &
            (( m1(j,i), i=1,nx), j=1,ny), &
            (( m2(j,i), i=1,nx), j=1,ny), &
            (( n1(j,i), i=1,nx), j=1,ny), &
            (( n2(j,i), i=1,nx), j=1,ny), &
            ((m11(j,i), i=1,nx), j=1,ny), &
            ((m12(j,i), i=1,nx), j=1,ny), &
            ((m22(j,i), i=1,nx), j=1,ny), &
            ((n11(j,i), i=1,nx), j=1,ny), &
            ((n12(j,i), i=1,nx), j=1,ny), &
            ((n22(j,i), i=1,nx), j=1,ny)
        else
          read(10) m1, m2, n1, n2, m11, m12, m22, n11, n12, n22
        endif
        close(10)

!.... compute the wall normals

        allocate( bn1(nx), bn2(nx) )
        j = 1
        do i = 1, nx
          bn1(i) = n1(j,i) / sqrt( n1(j,i)**2 + n2(j,i)**2 )
          bn2(i) = n2(j,i) / sqrt( n1(j,i)**2 + n2(j,i)**2 )
        end do

!.... Compute first derivatives of field in the mapped space (4th order)

        allocate( g1v(ny,nx,ndof), g2v(ny,nx,ndof) )
        call grad(ndof, nx, ny, v, g1v, g2v, dxi, deta, -1, -1, &
                  xper, yper, .false., .false., .false., .false.)
        
!.... Compute second derivatives of field (4th order)
        
        allocate( g11v(ny,nx,ndof), g12v(ny,nx,ndof), g22v(ny,nx,ndof) )
        call grad2(ndof, nx, ny, v, g1v, g11v, g12v, g22v, dxi, deta, &
                   -1, -1, xper, yper, .false., .false., &
                   .false., .false.)

!.... transform the gradients to physical space

        allocate(g1vl(ny), g2vl(ny), g11vl(ny), g12vl(ny), g22vl(ny))
        do idof = 1, ndof
          do i = 1, nx
            g1vl  = g1v(:,i,idof)*m1(:,i) + g2v(:,i,idof)*n1(:,i)
            g2vl  = g1v(:,i,idof)*m2(:,i) + g2v(:,i,idof)*n2(:,i)
              
            g11vl = g11v(:,i,idof)       * m1(:,i)*m1(:,i)      + &
                    two * g12v(:,i,idof) * m1(:,i)*n1(:,i)      + &
                    g22v(:,i,idof)       * n1(:,i)*n1(:,i)      + &
                    g1v(:,i,idof)        * m11(:,i)             + &
                    g2v(:,i,idof)        * n11(:,i)
  
            g12vl = g11v(:,i,idof)       * m1(:,i)*m2(:,i)      + &
                    g12v(:,i,idof)       * m1(:,i)*n2(:,i)      + &
                    g12v(:,i,idof)       * m2(:,i)*n1(:,i)      + &
                    g22v(:,i,idof)       * n1(:,i)*n2(:,i)      + &
                    g1v(:,i,idof)        * m12(:,i)             + &
                    g2v(:,i,idof)        * n12(:,i)
  
            g22vl = g11v(:,i,idof)       * m2(:,i)*m2(:,i)      + &
                    two * g12v(:,i,idof) * m2(:,i)*n2(:,i)      + &
                    g22v(:,i,idof)       * n2(:,i)*n2(:,i)      + &
                    g1v(:,i,idof)        * m22(:,i)             + &
                    g2v(:,i,idof)        * n22(:,i)
  
            g1v(:,i,idof)  = g1vl
            g2v(:,i,idof)  = g2vl
            g11v(:,i,idof) = g11vl
            g12v(:,i,idof) = g12vl
            g22v(:,i,idof) = g22vl
          end do
        end do
        deallocate( g1vl, g2vl, g11vl, g12vl, g22vl )
        
!.... form the output variables

        allocate( q(ny,nx,ndof) )

        q(:,:,1) = v(:,:,1)                                     ! rho
        q(:,:,2) = v(:,:,2)                                     ! u-velocity
        q(:,:,3) = v(:,:,3)                                     ! v-velocity
        q(:,:,4) = ( v(:,:,1) * vm(:,:,5) + &
                    vm(:,:,1) *  v(:,:,5) ) / (gamma * Ma**2)   ! pressure
        q(:,:,5) = v(:,:,5)                                     ! temperature

!.... subtract off a Gaussian Pulse to form the error

        if (sub) then
          do i = 1, nx
            do j = 1, ny
              rhom = vm(j,i,1)
              tm   = vm(j,i,5)
              cm   = sqrt( tm ) / Ma
              
              c1 =  zero
              c2 =  zero
              c3 =  exp( -pt5 * ( ((y(j,i)-time) + 2.0) / 0.05 )**2 )
              c4 =  exp( -pt5 * ( ((y(j,i)+time) - 2.0) / 0.05 )**2 )
              
              rho = ( -c1 + pt5 * ( c3 + c4 ) ) / cm**2
              u1  = c2 / ( rhom * cm ) 
              u2  = ( c3 - c4 ) * pt5 / ( rhom * cm )
              p   = ( c3 + c4 ) * pt5
              t   = ( gamma * Ma**2 * p - tm * rho ) / rhom
              
              q(j,i,1) = q(j,i,1) - rho
              q(j,i,2) = q(j,i,2) - u1
              q(j,i,3) = q(j,i,3) - u2
              q(j,i,4) = q(j,i,4) - p
              q(j,i,5) = q(j,i,5) - t
            end do
          end do
        end if

!.... write out data on a ray from the LE

        i = 1
        do j = 1, ny
          write(11,"(6(1pe13.6,1x))") x(j,i), q(j,i,1), q(j,i,2), q(j,i,3), &
                                      q(j,i,4), q(j,i,5)
        end do
        close(11)

!.... write out data along a line in x

        j = 1
!       write(*,"('Enter j ==> ',$)")
!       read(*,*) j
        do i = 1, nx
          write(12,"(6(1pe14.6E3,1x))") x(j,i), q(j,i,1), q(j,i,2), q(j,i,3), &
                                        q(j,i,4), q(j,i,5)
        end do
        close(12)

!.... write out data on the top boundary

        j = ny-5
        do i = 1, nx
          write(13,"(6(1pe14.6E3,1x))") x(j,i), q(j,i,1), q(j,i,2), q(j,i,3), &
                                        q(j,i,4), q(j,i,5)
        end do
        close(13)

!.... write out data on the inflow boundary

        i = 1
!       write(*,"('Enter i ==> ',$)")
!       read(*,*) i
        do j = 1, ny
          write(14,"(6(1pe14.6E3,1x))") y(j,i), q(j,i,1), q(j,i,2), q(j,i,3), &
                                        q(j,i,4), q(j,i,5)
        end do
        close(14)

!.... write out the wall quantities:
!....           x station
!....           wall vorticity
!....           rho
!....           t
!....           ut
!....           un

        if (.true.) then
        open(10,file='wall.dat',form='formatted',status='unknown')
        do i = 1, nx
          m1l = n1(1,i) / sqrt( n1(1,i)**2 + n2(1,i)**2 )
          m2l = n2(1,i) / sqrt( n1(1,i)**2 + n2(1,i)**2 )
          write(10,11) x(1,i),                                          &
                       g1v(1,i,3) - g2v(1,i,2),                         &
                       v(1,i,1), v(1,i,5),                              &
                       m1l * v(1,i,2) + m2l * v(1,i,3),                 &
                      -m2l * v(1,i,2) + m1l * v(1,i,3)
        end do
        close(10)
        end if

!.... write out a plot3d file

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
 40     open(unit=10, file=file, form='unformatted', status='unknown', err=50)
        goto 60
 50     write (*,"('>> Error opening [',a,'] ',/)") file(1:index(file,' '))
        write (*,"(/,'Enter PLOT3D file name [',a,']? ',$)") &
              file(1:index(file,' '))
        read (*,"(a80)") temp
        if (temp(1:1) .ne. ' ') file = temp
        goto 40
 60     continue
              
        write(10,err=1010) nx, ny, nz
        write(10,err=1010) Ma, Pr, Re, time
        write(10,err=1010) ((((q(j,i,idof), i = 1, nx), j = 1, ny),  &
                            k = 1, nz), idof = 1, ndof)
        close(10,err=1010)
        write(*,"('Wrote first Plot3d file')")

        if (deriv) then

        q(:,:,1) = g1v(:,:,3) - g2v(:,:,2)      ! vorticity
        q(:,:,2) = g1v(:,:,2) + g2v(:,:,3)      ! dilatation
        q(:,:,3) = g1v(:,:,5) * n1(:,:) / sqrt( n1(:,:)**2 + n2(:,:)**2 ) + &
                   g2v(:,:,5) * n2(:,:) / sqrt( n1(:,:)**2 + n2(:,:)**2 )
        q(:,:,4) = ( v(:,:,1) * vm(:,:,5) + &
                    vm(:,:,1) *  v(:,:,5) ) / (gamma * Ma**2)
        q(:,:,5) = (g1v(:,:,2) * n1(:,:) / sqrt( n1(:,:)**2 + n2(:,:)**2 )  + &
                    g2v(:,:,2) * n2(:,:) / sqrt( n1(:,:)**2 + n2(:,:)**2 )) * &
                    n2(:,:) / sqrt( n1(:,:)**2 + n2(:,:)**2 ) - &
                   (g1v(:,:,3) * n1(:,:) / sqrt( n1(:,:)**2 + n2(:,:)**2 )  + &
                    g2v(:,:,3) * n2(:,:) / sqrt( n1(:,:)**2 + n2(:,:)**2 )) * &
                    n1(:,:) / sqrt( n1(:,:)**2 + n2(:,:)**2 )

!.... write out a plot3d file

        iloc = index(file,'.q.')
        iend = index(file,' ')-1
        if (iloc .ne. 0) then
          temp = file(1:iloc)//'d'//file(iloc+2:iend)
        else
          temp = file(1:iend)//'.d'
        end if
        if ( nfile .eq. 0 ) then
          file = temp
          write (*,"(/,'Enter PLOT3D file name [',a,']? ',$)")  &
                  file(1:index(file,' '))
          read (*,"(a80)") temp
          if (temp(1:1) .ne. ' ') file = temp
        else if ( nfile .eq. 1 .or. nfile .eq. 2) then
          file = temp
        else if ( nfile .eq. 3 ) then
          call getarg(ifile(3),temp)
          file = temp
        end if
 400    open(unit=10, file=file, form='unformatted', status='unknown', err=500)
        goto 600
 500    write (*,"('>> Error opening [',a,'] ',/)") file(1:index(file,' '))
        write (*,"(/,'Enter PLOT3D file name [',a,']? ',$)") &
              file(1:index(file,' '))
        read (*,"(a80)") temp
        if (temp(1:1) .ne. ' ') file = temp
        goto 400
 600    continue
              
        write(10,err=1010) nx, ny, nz
        write(10,err=1010) Ma, Pr, Re, time
        write(10,err=1010) ((((q(j,i,idof), i = 1, nx), j = 1, ny),  &
                            k = 1, nz), idof = 1, ndof)
        close(10,err=1010)
        write(*,"('Wrote second Plot3d file')")
        
        end if

        call exit(0)    

 11     format(9(1pe20.11e4,1x))

 1010   write(*,"('>> Error writing PLOT3D file')")
        call exit(1)

        end
