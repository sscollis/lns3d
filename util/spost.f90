!=============================================================================!
        module bspline

        integer           :: ny, kyord
        real, allocatable :: yknot(:)
        real, allocatable :: bs(:,:)

        real, external :: bsder

        real :: rhoe, ue, ve, we, te, pe

        end module bspline
!============================================================================!
        program spost 
!  
!  Post processor for swept, nonlinear calculations
!  
!============================================================================!
        use bspline
        use const
        implicit none

!.... flow data

        real, allocatable :: v(:,:,:), q(:,:,:), p(:,:)
        
        real, allocatable :: g1v(:,:,:), g2v(:,:,:)
        real, allocatable :: g11v(:,:,:), g12v(:,:,:), g22v(:,:,:)
        
        real, allocatable :: g1vl(:), g2vl(:), g11vl(:), g12vl(:), &
                             g22vl(:)

!.... mesh

        real, allocatable :: x(:,:), y(:,:), xi(:), eta(:), s(:)
        real :: dxi, deta, dz, z(10)

!.... metrics

        real, allocatable :: m1(:,:),  m2(:,:),  n1(:,:),  n2(:,:),   &
                             m11(:,:), m12(:,:), m22(:,:),            &
                             n11(:,:), n12(:,:), n22(:,:)
        real :: m1l, m2l, bn1, bn2

!.... parameters

        real    :: Ma, Re, Pr, gamma, cv, Rgas, time, gamma1, vinf, &
                   alpha
        integer :: lstep, nx, nz, ndof

        integer :: i, j, k, idof, itmp

        character(80) file, temp, base, fname
        integer :: iloc, iend, imin, imax, inc, iver
#ifndef __GFORTRAN__    
        integer, external :: iargc
#endif

        integer :: nx2, ny2, nz2

        real :: tmp, arc
        
        real :: b, beta, eps

!.... variables for thickness calculation

        real :: rho, u1, u2, t, delta, theta, us, un, edge, errest
        integer :: ierr
        real, allocatable :: ynn(:), uss(:), unn(:)

        real, external :: rtsafe
        external funcd, dfun, thfun
!============================================================================!
!.... read the restart file

        if ( iargc() .gt. 0 ) then
          call getarg(1,temp)
          file = temp
        else
          file = 'output.R.1'
          write (*,"('Enter file name [',a,']? ',$)") &
            file(1:index(file,' '))
          read (*,"(a80)") temp
          if (temp(1:1) .ne. ' ') file = temp
        end if
 10     open(unit=10, file=file, form='unformatted', status='old', err=20)
        goto 30
 20     write (*,"('>> Error opening [',a,'] ',/)") &
          file(1:index(file,' '))
        write (*,"('Enter file name [',a,']? ',$)") &
          file(1:index(file,' '))
        read (*,"(a80)") temp
        if (temp(1:1) .ne. ' ') file = temp
        goto 10
 30     continue
        
        read(10) lstep, time, nx, ny, nz, ndof, &
                 Re, Ma, Pr, gamma, cv
!       write(*,*) lstep, time, nx, ny, nz, ndof, &
!                  Re, Ma, Pr, gamma, cv
        gamma1 = gamma - one
        Rgas = gamma1 * cv
        if (nz .ne. 1) then
          write(*,"(' Error:  nz <> 1 ')")
          call exit(1)
        end if
        allocate( v(ny,nx,ndof) )
        read(10) v
        close(10)
        write(*,"('Read flow field for ',a)") file(1:index(file,' '))
        write(*,"('  Time = ',1pe10.3,'  step = ',i6)") time, lstep

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
 
!.... read in the body discription file (if available)

        s = x(1,:)
        open(unit=10,file='body.dat',form='formatted',status='old',err=100)
        do i = 1, nx
           read(10,*) tmp, s(i)
        end do
        close(10)
 100    continue

!.... make the xi grid
        
        dxi = one / float(nx-1)
        
        do i = 1, nx
          xi(i) = real(i-1) * dxi
        end do

!.... make the eta grid

        deta = one / float(ny-1)
        
        do j = 1, ny
          eta(j) = float(j-1) * deta
        end do

!.... filter the field?

        if (.false.) then
          write(*,*) 'Filtering in Y'
          call filter( ny, nx, v(:,:,1) )
          call filter( ny, nx, v(:,:,2) )
          call filter( ny, nx, v(:,:,3) )
          call filter( ny, nx, v(:,:,4) )
          call filter( ny, nx, v(:,:,5) )
          open(unit=10,file='filter.dat',form='unformatted', &
               status='unknown')
          write(10) lstep, time, nx, ny, nz, ndof, &
                    Re, Ma, Pr, gamma, cv
          write(10) v
          close(10)
        end if
        
!.... allocate storage for metrics

        allocate (m1(ny,nx),  m2(ny,nx),  n1(ny,nx),  n2(ny,nx), &
                  m11(ny,nx), m12(ny,nx), m22(ny,nx),            &
                  n11(ny,nx), n12(ny,nx), n22(ny,nx) )

!.... read in the metric file

        open (unit=10,file='metric.dat',form='unformatted', &
              status='old')
        read(10) m1, m2, n1, n2, m11, m12, m22, n11, n12, n22
        close(10)
        
!.... Compute first derivatives of field in the mapped space

        allocate( g1v(ny,nx,ndof), g2v(ny,nx,ndof) )
        call grad(ndof, nx, ny, v, g1v, g2v, dxi, deta, -1, -1, &
                  .false., .false.)
        
!.... Compute second derivatives of field
                
        allocate( g11v(ny,nx,ndof), g12v(ny,nx,ndof), g22v(ny,nx,ndof) )
        call grad2(ndof, nx, ny, v, g1v, g11v, g12v, g22v, dxi, deta, &
                   -1, -1, .false., .false.)

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

        allocate( q(ny,nx,ndof), p(ny,nx) )

        P    = v(:,:,1) * v(:,:,5) / (gamma * Ma**2)

        write(*,"('Enter sweep angle (deg) ==> ',$)")
        read(*,*) alpha
        alpha = alpha * pi / 180.0
        vinf = sqrt( one + tan(alpha)**2 )

!.... write out primitive variables or conservative variables in Plot3d
!.... normalization

        if (.false.) then
          q(:,:,1) = v(:,:,1)                                   ! rho
          q(:,:,2) = v(:,:,2)                                   ! u-velocity
          q(:,:,3) = v(:,:,3)                                   ! v-velocity
          q(:,:,4) = two*(p(:,:)-one/(gamma*Ma**2))             ! Cp
          q(:,:,5) = v(:,:,5)                                   ! temperature
        else
          q(:,:,1) = v(:,:,1)
          q(:,:,2) = v(:,:,1) * v(:,:,2) * Ma
          q(:,:,3) = v(:,:,1) * v(:,:,3) * Ma
          q(:,:,4) = v(:,:,1) * v(:,:,4) * Ma
          q(:,:,5) = v(:,:,1) * ( one/gamma1/gamma * v(:,:,5) + pt5 * &
                     Ma**2 * (v(:,:,2)**2+v(:,:,3)**2+v(:,:,4)**2) )  
        end if

!.... compute potential flow quantities

!       b = 1.0 / pi
!       beta = sqrt(one - Ma**2)
!       q(:,:,1) = one + (x(:,:)-b)/pi/beta**2/((x(:,:)-b)**2+&
!                         beta**2*y(:,:)**2)
!       q(:,:,2) = v(:,:,2)                                     ! u-velocity
!       q(:,:,3) = v(:,:,3)                                     ! v-velocity
!       q(:,:,4) = beta*y(:,:)/pi/((x(:,:)-b)**2+beta**2*y(:,:)**2)
!       q(:,:,5) = v(:,:,5)                                     ! temperature

!.... write out quantities along the ray from the nose

        i = 1
        do j = 1, ny
          write(11,"(7(1pe16.9,1x))") x(j,i), v(j,i,1), v(j,i,2), &
                                      v(j,i,3), &
                                      v(j,i,4), v(j,i,5), p(j,i)
        end do

!.... write out quantities along the far-field boundary

        j = ny
!       write(*,"('Enter j ==> ',$)")
!       read(*,*) j
        do i = 1, nx
          m1l = n1(j,i) / sqrt( n1(j,i)**2 + n2(j,i)**2 )
          m2l = n2(j,i) / sqrt( n1(j,i)**2 + n2(j,i)**2 )
          un  = v(j,i,2) * m1l + v(j,i,3) * m2l 
          write(12,"(8(1pe16.9,1x))") &
            x(j,i), y(j,i), v(j,i,1), v(j,i,2), &
            v(j,i,3), v(j,i,4), v(j,i,5), p(j,i)
        end do

!.... write out quantities along the outflow boundnary

        i = nx
        do j = 1, ny
          write(13,"(7(1pe16.9,1x))") &
            y(j,i), v(j,i,1), v(j,i,2), v(j,i,3), &
            v(j,i,4), v(j,i,5), p(j,i)
        end do

!.... write out quantities along the 1 node in from the outflow boundnary

        i = nx-1
        do j = 1, ny
          write(14,"(7(1pe16.9,1x))") y(j,i), v(j,i,1), v(j,i,2), v(j,i,3), &
                                      v(j,i,4), v(j,i,5), p(j,i)
        end do

!.... write out quantities along the 10 nodes in from the outflow boundnary

        i = nx-10
        do j = 1, ny
          write(15,"(7(1pe16.9,1x))") y(j,i), v(j,i,1), v(j,i,2), v(j,i,3), &
                                      v(j,i,4), v(j,i,5), p(j,i)
        end do

!.... write out quantities along the 20 nodes in from the outflow boundnary

        i = nx-20
        do j = 1, ny
          write(16,"(7(1pe16.9,1x))") y(j,i), v(j,i,1), v(j,i,2), v(j,i,3), &
                                      v(j,i,4), v(j,i,5), p(j,i)
        end do

!.... write out data on a ray from the LE

        eps = 0.0001
        i = 1
        do j = 1, ny
          write(17,"(7(1pe14.6E3,1x))") y(j,i), &
                                        (v(j,i,1)-one)/eps, &
                                        v(j,i,2)/eps, &
                                        v(j,i,3)/eps, &
                                        v(j,i,4)/eps, &
                                        (v(j,i,5)-one)/eps, &
                                        (p(j,i) - one/(gamma*Ma**2))/eps
        end do
        close(17)

!.... write out profiles at user supplied locations

        if (.false.) then
        write(*,"('Enter the x-stations for profiles [min,max,inc] ==> ',$)")
        read(*,*) imin, imax, inc
        base = 'profile'
        if (imin.ne.0) then
        open(11,file='loc.dat')
        do i = imin, imax, inc
          call makename(base,i,fname)
          iloc = index(fname,' ')
          write(*,110) s(i), x(1,i), fname(1:iloc)
 110      format('Saving profile at s, x_b = ',2(1pe20.13,1x),'in file: ',a)
          write(11,"(i4,1x,2(1pe13.6,1x))") i, s(i), x(1,i)
          arc = zero
          bn1 = n1(1,i) / sqrt( n1(1,i)**2 + n2(1,i)**2 )
          bn2 = n2(1,i) / sqrt( n1(1,i)**2 + n2(1,i)**2 )
          open(10,file=fname)
          do j = 1, ny
            us = bn2 * v(j,i,2) - bn1 * v(j,i,3)
            un = bn1 * v(j,i,2) + bn2 * v(j,i,3)
!           write(10,"(7(1pe16.9,1x))") arc, v(j,i,1), us, un, v(j,i,4), &
!                                       v(j,i,5), p(j,i)
            write(10,"(7(1pe16.9,1x))") arc, v(j,i,1), v(j,i,2), v(j,i,3), &
                                        v(j,i,4), v(j,i,5), p(j,i)
            if (j.ne.ny) arc = arc + sqrt( (x(j+1,i)-x(j,i))**2 + &
                                           (y(j+1,i)-y(j,i))**2 )
          end do
          close(10)
        end do
        close(11)
        end if
        end if

!.... write out a plot3d file

        iloc = index(file,'.R.')
        iend = index(file,' ')-1
        if (iloc .ne. 0) then
          temp = file(1:iloc)//'e'//file(iloc+2:iend)
        else
          temp = file(1:iend)//'.e'
        end if

        if ( iargc() .eq. 0 ) then
          file = temp
          write (*,"(/,'Enter PLOT3D file name [',a,']? ',$)")  &
                  file(1:index(file,' '))
          read (*,"(a80)") temp
          if (temp(1:1) .ne. ' ') file = temp
        else if ( iargc() .eq. 1 ) then
          file = temp
        else if ( iargc() .eq. 2 ) then
          call getarg(2,temp)
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

        write(*,*) 'Writing data file...'
        nz = 10
        write(10,err=1010) nx, ny, nz
        write(10,err=1010) Ma * vinf, alpha*180.0/pi, Re * vinf, time
        write(10,err=1010) ((((q(j,i,idof), i = 1, nx), j = 1, ny),  &
                            k = 1, nz), idof = 1, ndof)
        close(10,err=1010)

!.... write out a 3d grid file

        dz = 500.0/float(nz-1)
        do k = 1, nz
          z(k) = float(k-1)*dz
        end do
        write(*,*) 'Writing 3D grid file...'
        open(unit=10,file='grid3d.dat',form='unformatted',status='unknown')
        write(10) nx, ny, nz
        write(10) (((x(j,i), i = 1, nx), j = 1, ny), k = 1, nz), &
                  (((y(j,i), i = 1, nx), j = 1, ny), k = 1, nz), &
                  (((  z(k), i = 1, nx), j = 1, ny), k = 1, nz)
        close(10)

!.... write out a second plot3d file with other stuff

        if (.false.) then    ! don't mess with the gradients

        if (.false.) then

!.... look at the second derivatives
        
!         q(:,:,1) = g22v(:,:,1)
!         q(:,:,2) = g22v(:,:,2)
!         q(:,:,3) = g22v(:,:,3)
!         q(:,:,4) = g22v(:,:,4)
!         q(:,:,5) = g22v(:,:,5)

          q(:,:,1) = sqrt( v(:,:,5) ) / Ma                      ! c
          q(:,:,2) = v(:,:,2) / q(:,:,1)                        ! u-Mach
          q(:,:,3) = v(:,:,3) / q(:,:,1)                        ! v-Mach
          q(:,:,4) = sqrt(v(:,:,2)**2 + v(:,:,3)**2) / q(:,:,1) ! Mach
          q(:,:,5) = log( v(:,:,1) * v(:,:,5) / &
                          (gamma * Ma**2) / v(:,:,1) ** gamma ) ! entropy

        else

          q(:,:,1) = g1v(:,:,3) - g2v(:,:,2)
          q(:,:,2) = g1v(:,:,2) + g2v(:,:,3)
          q(:,:,3) = g1v(:,:,1) * v(:,:,5) + v(:,:,1) * g1v(:,:,5)
          q(:,:,4) = g2v(:,:,1) * v(:,:,5) + v(:,:,1) * g2v(:,:,5)
          q(:,:,5) = v(:,:,1) * v(:,:,5) / (gamma * Ma**2) + &
                     pt5 * v(:,:,1) * ( v(:,:,2)**2 + &
                     v(:,:,3)**2 + v(:,:,4)**2 )

        end if
        
!.... write out a plot3d file

        iloc = index(file,'.q.')
        iend = index(file,' ')-1
        if (iloc .ne. 0) then
          temp = file(1:iloc)//'d'//file(iloc+2:iend)
        else
          temp = file(1:iend)//'.d'
        end if

        if ( iargc() .eq. 0 ) then
          file = temp
          write (*,"(/,'Enter PLOT3D file name [',a,']? ',$)")  &
                  file(1:index(file,' '))
          read (*,"(a80)") temp
          if (temp(1:1) .ne. ' ') file = temp
        else if ( iargc() .eq. 1 ) then
          file = temp
        else if ( iargc() .eq. 2 ) then
          call getarg(2,temp)
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

        end if

!.... write out the wall quantities:
!....   1       x station
!....   2       arc length
!....   3       wall vorticity
!....   4       -Cp
!....   5       dp/ds
!....   6       du/dn
!....   7       dt/dn
!....   8       us
!....   9       un
!....   0       t

        open(unit=10, file='wall.dat', form='formatted', status='unknown')
        do i = 1, nx
          m1l = n1(1,i) / sqrt( n1(1,i)**2 + n2(1,i)**2 )
          m2l = n2(1,i) / sqrt( n1(1,i)**2 + n2(1,i)**2 )
          write(10,11) x(1,i), s(i),                                    &
                       g1v(1,i,3) - g2v(1,i,2),                         &
                       -two * (v(1,i,1) * v(1,i,5) - one) /             &
                          (gamma * Ma**2),                              &
                       (m2l * (g1v(1,i,1) *   v(1,i,5)  +               &
                                 v(1,i,1) * g1v(1,i,5)) -               &
                        m1l * (g2v(1,i,1) *   v(1,i,5)  +               &
                                 v(1,i,1) * g2v(1,i,5))) /              &
                          (gamma * Ma**2),                              &
                       (g1v(1,i,2) * m1l + g2v(1,i,2) * m2l) * m2l -    &
                       (g1v(1,i,3) * m1l + g2v(1,i,3) * m2l) * m1l,     &
                       g1v(1,i,5) * m1l + g2v(1,i,5) * m2l,             &
                       v(1,i,2)*m2l-v(1,i,3)*m1l,                       &
                       v(1,i,2)*m1l+v(1,i,3)*m2l,                       &
                       v(1,i,5)
        end do
        close(10)
        
!.... Setup the B-spline interpolation
            
        if (.false.) then
        
        write(*,"('Enter korder (0 = terminate) ==> ',$)")
        read(*,*) kyord
        if (kyord.eq.0) call exit(0)

        allocate ( yknot(ny+kyord), bs(ny,ndof), ynn(ny), uss(ny), unn(ny) )

!.... compute the boundary layer parameters

        open(10,file='delta.dat',status='unknown')
        do i = 2, nx

          bn1 = n1(1,i) / sqrt( n1(1,i)**2 + n2(1,i)**2 )
          bn2 = n2(1,i) / sqrt( n1(1,i)**2 + n2(1,i)**2 )
          
          arc = zero
          do j = 1, ny
            ynn(j) = arc
            if (j.ne.ny) arc = arc + sqrt( (x(j+1,i)-x(j,i))**2 + &
                                           (y(j+1,i)-y(j,i))**2 )
            uss(j) = bn2 * v(j,i,2) - bn1 * v(j,i,3)
            unn(j) = bn1 * v(j,i,2) + bn2 * v(j,i,3)
          end do
          
          call BSNAK( ny, ynn, kyord, yknot)
          call BSINT( ny, ynn, v(:,i,1), kyord, yknot, bs(:,1) )
          call BSINT( ny, ynn,      uss, kyord, yknot, bs(:,2) )
          call BSINT( ny, ynn,      unn, kyord, yknot, bs(:,3) )
          call BSINT( ny, ynn, v(:,i,4), kyord, yknot, bs(:,4) )
          call BSINT( ny, ynn, v(:,i,5), kyord, yknot, bs(:,5) )

          do j = 1, ny-1
            if ( uss(j+1) .lt. uss(j) ) goto 70         ! rough estimate
          end do
 70       edge = rtsafe( funcd, ynn(j-1), ynn(j+2), 1.0e-14 )
 
          rhoe = BSDER( 0, edge, kyord, yknot, ny, bs(:,1) )
          ue   = BSDER( 0, edge, kyord, yknot, ny, bs(:,2) )
          ve   = BSDER( 0, edge, kyord, yknot, ny, bs(:,3) )
          we   = BSDER( 0, edge, kyord, yknot, ny, bs(:,4) )
          te   = BSDER( 0, edge, kyord, yknot, ny, bs(:,5) )
          pe   = rhoe * te / (gamma * Ma**2)
          write(*,"(8(1pe13.6,1x))") x(1,i), edge, rhoe, ue, ve, pe, te
          
          call qdag(  dfun, zero, edge, 1.0e-12, 1.0e-12, 2, delta, errest)
          call qdag( thfun, zero, edge, 1.0e-12, 1.0e-12, 2, theta, errest)

          write(10,"(8(1pe13.6,1x))") x(1,i), delta, theta

          if(.false.) then
            do j = 1, ny
              rho = BSDER( 0, ynn(j), kyord, yknot, ny, bs(:,1) )
              us  = BSDER( 0, ynn(j), kyord, yknot, ny, bs(:,2) )
              un  = BSDER( 0, ynn(j), kyord, yknot, ny, bs(:,3) )
              t   = BSDER( 0, ynn(j), kyord, yknot, ny, bs(:,5) )
              write(*,"(8(1pe13.6,1x))") ynn(j), rho, us, un, t
            end do
          end if
          
        end do
        close(10)

        end if

        call exit(0)    

 11     format(12(1pe16.9,1x))

 1010   write(*,"('>> Error writing PLOT3D file')")
        call exit(1)

!.... grid diagnostic (not used)

        if (.false.) then
        j = 1
        do i = 1, nx
          write(30,"(8(1pe13.6,1x))")                                    &
              (-m1(j,i)*m2(j,i)*m11(j,i) + m1(j,i)*m1(j,i)*m12(j,i) -    &
              m2(j,i)*m2(j,i)*m12(j,i) + m2(j,i)*m1(j,i)*m22(j,i) ) *    &
              sqrt(n1(j,i)**2 + n2(j,i)**2) /                            &
              sqrt(m1(j,i)**2 + m2(j,i)**2) /                            &
              (m1(j,i)*n2(j,i)-m2(j,i)*n1(j,i)),                         &
            ( n1(j,i)*n2(j,i)*n11(j,i) - n1(j,i)*n1(j,i)*n12(j,i) +      &
              n2(j,i)*n2(j,i)*n12(j,i) - n2(j,i)*n1(j,i)*n22(j,i) ) *    &
              sqrt(m1(j,i)**2 + m2(j,i)**2) /                            &
              sqrt(n1(j,i)**2 + n2(j,i)**2) /                            &
              (m1(j,i)*n2(j,i)-m2(j,i)*n1(j,i)),                         &
            ( n1(j,i)*n2(j,i)*m11(j,i) - n1(j,i)*n1(j,i)*m12(j,i) +      &
              m1(j,i)*n2(j,i)*n11(j,i) - m1(j,i)*n1(j,i)*n12(j,i) +      &
              n2(j,i)*n2(j,i)*m12(j,i) - n2(j,i)*n1(j,i)*m22(j,i) +      &
              m2(j,i)*n2(j,i)*n12(j,i) - m2(j,i)*n1(j,i)*n22(j,i) ) *    &
              sqrt(m1(j,i)**2 + m2(j,i)**2) /                            &
              sqrt(n1(j,i)**2 + n2(j,i)**2) /                            &
              (m1(j,i)*n2(j,i)-m2(j,i)*n1(j,i)),                         &
               m1(j,i)*n1(j,i) + m2(j,i)*n2(j,i)
        end do
        end if

        end

!=============================================================================!
        function dfun(x)

        use bspline

        real :: x, dfun, rho, u
!=============================================================================!

        rho = BSDER( 0, x, kyord, yknot, ny, bs(:,1) )
        u   = BSDER( 0, x, kyord, yknot, ny, bs(:,2) )

        dfun = 1.0 - rho * u / (rhoe * ue)
        
        return
        end 
!=============================================================================!
        function thfun(x)

        use bspline

        real :: x, thfun, rho, u
!=============================================================================!

        rho = BSDER( 0, x, kyord, yknot, ny, bs(:,1) )
        u   = BSDER( 0, x, kyord, yknot, ny, bs(:,2) )

        thfun = rho * u / (rhoe * ue) * ( 1.0 - u / ue )
        
        return
        end 
!=============================================================================!
        subroutine funcd( x, g, d )

          use bspline

          real :: x, f, g, d
!=============================================================================!

          f = BSDER( 0, x, kyord, yknot, ny, bs(:,2) )
          g = BSDER( 1, x, kyord, yknot, ny, bs(:,2) )
          d = BSDER( 2, x, kyord, yknot, ny, bs(:,2) )

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
        write(fname(length:80),"('.',i1)") iver
      else if (iver .lt. 100) then
        write(fname(length:80),"('.',i2)") iver
      else
        write(fname(length:80),"('.',i3)") iver
      end if

      return
      end
