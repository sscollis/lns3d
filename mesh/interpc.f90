!=============================================================================!
        module constants

          real, parameter :: zero    = 0.0000000000000000000d+0
          real, parameter :: pt25    = 2.5000000000000000000d-1
          real, parameter :: pt33    = 3.3333333333333333333d-1
          real, parameter :: pt5     = 5.0000000000000000000d-1
          real, parameter :: pt66    = 6.6666666666666666666d-1
          real, parameter :: one     = 1.0000000000000000000d+0
          real, parameter :: onept25 = 1.2500000000000000000d+0   
          real, parameter :: onept33 = 1.3333333333333333333d+0
          real, parameter :: onept5  = 1.5000000000000000000d+0
          real, parameter :: two     = 2.0000000000000000000d+0
          real, parameter :: three   = 3.0000000000000000000d+0
          real, parameter :: four    = 4.0000000000000000000d+0
          real, parameter :: pi      = 3.1415926535897932385d+0

        end module constants
!=============================================================================!
        program inter
!
!  This program interpolates from one computational mesh to another
!  using B-Spline basis function with the order set by the user.
!
!  This requires that you specify the analytical mapping functions in xi,eta
!  for the original mesh.
!  
!  Revised:  6/14/2022
!
!=============================================================================!
        use constants
        implicit none

!.... first grid and data

        integer :: nx1, ny1, nz1
        real, allocatable :: xy1(:,:,:), v1(:,:,:), xi1(:), eta1(:)
        real :: dxi1, deta1
        
!.... second grid and data

        integer :: nx2, ny2, nz2
        real, allocatable :: xy2(:,:,:), v2(:,:,:), xi2(:), eta2(:)
        real :: dxi2, deta2
        real, allocatable :: p(:,:), g1p(:,:), g2p(:,:)

!.... stuff for B-spline

        integer           :: kxord, kyord
        real, allocatable :: xknot(:), yknot(:)
        real, allocatable :: bsv(:,:,:)

!.... stuff for file I/O

        character(80) :: name, filen

!.... local vars

        real :: Re, Ma, Pr, gamma, cv
        integer :: i, j, k, idof, ndof=5, lstep
        real :: tmp, time

        real, external :: arc, calcxi
        
!.... stuff for mapping spline

        integer, parameter :: npts = 500
        real :: nn(npts), rr(npts), n2(npts), dnn
        real :: ss(npts), xx(npts), s2(npts), dss
        real :: rmax, rd1, rd2, dd, ximin, ximax
        real :: sx, sy, dxmin, dymin, c1, c2
        real, external :: calcdd, calcs
#if defined(USE_IMSL) || defined(USE_BSLIB)
        real, external :: BS2DR
        external :: BS2IN
#endif        
        real :: x, y, xi, eta, s, r, n

        real :: drmin, dsmin, b, rc, sc, cm, smin, smax

!.... argument parameters

        integer :: iarg, narg
        character(80) :: arg
        integer :: yflag=1, xflag=0
        logical :: plot3d = .false., pot=.false., useIJ=.true.
#ifndef __GFORTRAN__
        integer, external :: iargc
#endif
!=============================================================================!

!.... parse the argument list

        narg = iargc()
        do iarg = 1, narg
          call getarg(iarg,arg)
          select case (arg(1:2))
          case ('-y')
            select case (arg(1:3))
            case ('-y1')                ! exponential
              yflag = 1
            case ('-y2')                ! hyperbolic tangent
              yflag = 2
            case ('-y3')                ! Mahesh's mapping
              yflag = 3
            case default                ! uniform
              yflag = 0
            end select
          case ('-x')
            select case (arg(1:3))
            case ('-x1')                ! exponential
              xflag = 1
            case ('-x3')                ! Mahesh' mapping in xi
              xflag = 3
            case ('-x4')                ! Mahesh' mapping in s
              xflag = 4
            case default                ! uniform
              xflag = 0
            end select
          case ('-i')
            select case (arg(1:3))
            case ('-ij')
              useIJ = .true.            ! use IJ format (default)
            case ('-ji')
              useIJ = .false.           ! use JI format
            case default
              write(*,"('Argument ',i2,' ignored.')") iarg
             end select
          case ('-o')                   ! Output a Plot3d file
            plot3d = .true.
          case ('-p')                   ! generate potential BC data
            pot = .true.
          case ('-h')
            write(*,"('Usage:  interpc [options]')")
            write(*,"('------------------------------------------')")
            write(*,"('   -h:  this help')")
            write(*,"('------------------------------------------')")
            write(*,"('  -x1:  exponential stretching in x')")
            write(*,"('  -x3:  Mahesh''s mapping in xi')")
            write(*,"('  -x4:  Mahesh''s mapping in s')")
            write(*,"('------------------------------------------')")
            write(*,"('  -y1:  exponential stretching in y')")
            write(*,"('  -y2:  hyperbolic tangent stretching in y')")
            write(*,"('  -y3:  Mahesh''s mapping in y')")
            write(*,"('------------------------------------------')")
            write(*,"('   -o:  Output a Plot3d file')")
            write(*,"('   -p:  Generate potential info')")
            write(*,"('  -ij:  Solution in IJ format (default)   ')")
            write(*,"('  -ji:  Solution in JI format')")
            write(*,"('------------------------------------------')")
            call exit(0)
          case default
            write(*,"('Argument ',i2,' ignored.')") iarg
          end select
        end do

!.... read the first grid file

        filen = 'grid.dat'
 200    write (*,"('Enter grid to interpolate from [',a,']? ',$)") &
                    filen(1:index(filen,' '))
        read (*,"(a)") name
        if (name(1:1) .ne. ' ') filen = name

!.... allocate memory for the physical grid

        open (unit=10,file=filen,form='unformatted',status='old',err=200)
        read(10,err=200) nx1, ny1, nz1
        write(*,"('Nx = ',i4,', Ny = ',i4)") nx1, ny1
        allocate (xy1(ny1,nx1,2))
        read(10,err=200) &
                 (((xy1(j,i,1), i = 1, nx1), j = 1, ny1), k = 1, nz1), &
                 (((xy1(j,i,2), i = 1, nx1), j = 1, ny1), k = 1, nz1), &
                 (((       tmp, i = 1, nx1), j = 1, ny1), k = 1, nz1)
        close(10)

!.... allocate memory for the computational grid

        allocate (xi1(nx1), eta1(ny1))

!.... make the computational grid

        dxi1  = one / float(nx1-1)
        do i = 1, nx1
          xi1(i) = (i-1) * dxi1
        end do
        
        deta1 = one / float(ny1-1)
        do j = 1, ny1
          eta1(j) = (j-1) * deta1
        end do

!.... read in the first data file

        allocate(v1(ny1,nx1,ndof))
        filen = 'output.dat'
 201    write (*,"('Enter data filename to interpolate from [',a,']? ',$)") &
                    filen(1:index(filen,' '))
        read (*,"(a)") name
        if (name(1:1) .ne. ' ') filen = name

        open(unit=10, file=filen, form='unformatted', status='old', err=201)
        read(10,err=201) lstep, time, nx1, ny1, nz1, ndof, &
                         Re, Ma, Pr, gamma, cv
        if (useIJ) then
          read(10,err=201) (((v1(j,i,idof), idof=1,ndof), i=1,nx1), j=1,ny1)
        else
          read(10,err=201) v1
        endif
        close(10)

!.... read the second grid file

        filen = 'grid.dat'
 202    write (*,"('Enter grid to interpolate too [',a,']? ',$)") &
                    filen(1:index(filen,' '))
        read (*,"(a)") name
        if (name(1:1) .ne. ' ') filen = name

!.... allocate memory for the physical grid

        open (unit=10, file=filen, form='unformatted', status='old', err=202)
        read(10,err=202) nx2, ny2, nz2
        write(*,"('Nx = ',i4,', Ny = ',i4)") nx2, ny2
        allocate (xy2(ny2,nx2,2))
        read(10,err=202) &
                 (((xy2(j,i,1), i = 1, nx2), j = 1, ny2), k = 1, nz2), &
                 (((xy2(j,i,2), i = 1, nx2), j = 1, ny2), k = 1, nz2), &
                 (((       tmp, i = 1, nx2), j = 1, ny2), k = 1, nz2)
        close(10)

!.... allocate memory for the computational grid

        allocate (xi2(nx2), eta2(ny2))

!.... make the computational grid

        dxi2  = one / float(nx2-1)
        do i = 1, nx2
          xi2(i) = (i-1) * dxi2
        end do
        
        deta2 = one / float(ny2-1)
        do j = 1, ny2
          eta2(j) = (j-1) * deta2
        end do

!.... spline the tangent mapping function to get the inverse

        if(xy1(1,1,1) .le. zero) then
          ximin = zero
        else
          ximin = sqrt( xy1(1,1,1) )
        end if
        ximax = sqrt( xy1(1,nx1,1) )

        if (xflag.eq.1) then
          write(*,*) 'WARNING:  not corrected for nonzero origin'
          ximax = sqrt( xy1(1,nx1,1) )
          write(*,"('Enter dxmin ==> ',$)") 
          read(*,*) dxmin
          dss = one / dble(npts-1)
          sx = calcs( ximax, dxmin, dxi1 )
          write(*,*) 'Sx = ', sx
          c2 = log( dxmin / dxi1 )
          c1 = (log( sx*dxmin / dxi1) - c2) / dxi1
          do i = 1, npts
            ss(i) = dble(i-1) * dss
            xx(i) = one / c1 * ( exp(c1 * ss(i) + c2) - exp(c2) ) 
          end do
          call NR_SPLINE( xx, ss, npts, 1.0e31, 1.0e31, s2)
        else if (xflag.eq.3) then
          write(*,*) 'WARNING:  not corrected for nonzero origin'
          write(*,"('Enter dsmin, b, sc ==> ',$)") 
          read(*,*) dsmin, b, sc
          dsmin = dsmin / sqrt(two)
          smax = sqrt( xy1(1,nx1,1) )
          cm = ( two * b * tanh(b*sc) + (nx1-1)*dsmin/smax * &
                 log( cosh(b*(one-sc)) / cosh(b*(one+sc)) ) ) / &
               ( one - (nx1-1)*dsmin/smax )
          dss = one / dble(npts-1)
          do i = 1, npts
            ss(i) = dble(i-1) * dss
            xx(i) = smax*( cm * ss(i) + log( cosh(b*(ss(i)-sc)) / &
                    cosh(b*(ss(i)+sc)) ) ) / &
                    (cm + log( cosh(b*(one-sc)) / cosh(b*(one+sc)) ) )
          end do
          call NR_SPLINE( xx, ss, npts, 1.0e31, 1.0e31, s2)
        else if (xflag.eq.4) then
          write(*,"('Enter dsmin, b, sc ==> ',$)") 
          read(*,*) dsmin, b, sc
          smin = arc(ximin)
          smax = arc(ximax)
          write(*,"('Smin = ',1pe13.6,'  Smax = ',1pe13.6)") smin, smax
          cm = ( two * b * tanh(b*sc) + (nx1-1)*dsmin/(smax-smin) * &
                 log( cosh(b*(one-sc)) / cosh(b*(one+sc)) ) ) / &
               ( one - (nx1-1)*dsmin/(smax-smin) )
          dss = one / dble(npts-1)
          do i = 1, npts
            ss(i) = dble(i-1) * dss
            xx(i) = smin + (smax-smin)*( cm * ss(i) + &
                    log( cosh(b*(ss(i)-sc)) / &
                    cosh(b*(ss(i)+sc)) ) ) / &
                    (cm + log( cosh(b*(one-sc)) / cosh(b*(one+sc)) ) )
            xx(i) = calcxi( xx(i) )
          end do
          call NR_SPLINE( xx, ss, npts, 1.0e31, 1.0e31, s2)
        end if

!.... spline the normal mapping function to get the inverse

        eta = sqrt( pt5 - xy1(ny1,1,1) + &
                    pt5 * sqrt( (two*xy1(ny1,1,1)-one)**2 + &
                    four * xy1(ny1,1,2)**2) )
        rmax = abs( (one - eta**2) * pt5 )
        write(*,*) 'Rmax = ',rmax
        if (yflag.eq.1) then
          write(*,"('Enter dymin ==> ',$)") 
          read(*,*) dymin
          dnn = one / dble(npts-1)
          sy = calcs( rmax, dymin, deta1 )
          write(*,*) 'Sy = ', sy
          c2 = log( dymin / deta1 )
          c1 = (log( sy*dymin / deta1) - c2) / deta1
          do j = 1, npts
            nn(j) = dble(j-1) * dnn
            rr(j) = one / c1 * ( exp(c1 * nn(j) + c2) - exp(c2) ) 
            write(*,*) j, c1, nn(j), rr(j)
          end do
          write(*,*) 'Calling NR_SPLINE'
          call NR_SPLINE( rr, nn, npts, 1.0e31, 1.0e31, n2)
        else if (yflag.eq.2) then
          rd1  = 0.00028                ! set for R=1000 PCYL, r=500
          rd2  = 5.5
          dd   = calcdd(rd1,rd2)
          dnn = one / dble(npts-1)
          do j = 1, npts
            nn(j) = dble(j-1) * dnn
            rr(j) = rmax*(pt5 + one/Tanh(dd*pt5)*Tanh(dd*(-pt5 + &
                    nn(j)))*pt5)/(Sqrt(rd2/rd1) + &
                    (one - Sqrt(rd2/rd1))* &
                    (pt5 + one/Tanh(dd*pt5)*Tanh(dd*(-pt5 + &
                    nn(j)))*pt5))
          end do
          call NR_SPLINE( rr, nn, npts, 1.0e31, 1.0e31, n2)
        else if (yflag.eq.3) then
          write(*,"('Enter drmin, b, rc ==> ',$)") 
          read(*,*) drmin, b, rc
          cm = ( two * b * tanh(b*rc) + (ny1-1)*drmin/rmax * &
                 log( cosh(b*(one-rc)) / cosh(b*(one+rc)) ) ) / &
               ( one - (ny1-1)*drmin/rmax )
          dnn = one / dble(npts-1)
          do j = 1, npts
            nn(j) = dble(j-1) * dnn
            rr(j) = rmax*( cm * nn(j) + log( cosh(b*(nn(j)-rc)) / &
                    cosh(b*(nn(j)+rc)) ) ) / &
                    (cm + log( cosh(b*(one-rc)) / cosh(b*(one+rc)) ) )
          end do
          call NR_SPLINE( rr, nn, npts, 1.0e31, 1.0e31, n2)
        end if

        allocate(v2(ny2,nx2,ndof))

!.... form the B-Spline Interpolant (assume order 5 for 4th order accuracy)

        kxord = 5
!       write(*,"('Enter korder ==> ',$)")
!       read(*,*) kxord
        kyord = kxord

        allocate (xknot(nx1+kxord), yknot(ny1+kyord), bsv(ny1,nx1,ndof))

        write(*,*) "Starting B-splines"
        call BSNAK( nx1, xi1,  kxord, xknot)
        call BSNAK( ny1, eta1, kyord, yknot)
#if defined(USE_IMSL) || defined(USE_BSLIB)
        do idof = 1, ndof
          write(*,*) idof
          call BS2IN( ny1, eta1, nx1, xi1, v1(:,:,idof), ny1, kyord, &
                      kxord, yknot, xknot, bsv(:,:,idof))
        end do
        write(*,*) 'B-spline complete...begin interpolating'
#else
        write(*,*) 'Must use IMSL'
        stop 1
#endif
        if(xy1(1,1,1) .le. zero) then
          ximin = zero
        else
          ximin = sqrt( xy1(1,1,1) )
        end if
        ximax = sqrt( xy1(1,nx1,1) )

!.... now interpolate to the new mesh
        
        do i = 1, nx2
          write(*,"('i = ',i4)") i
          do j = 1, ny2
            x   = xy2(j,i,1)
            y   = xy2(j,i,2)
            
            eta = sqrt( pt5 - x + pt5 * sqrt( (two*x-one)**2 + &
                    four * y**2) )
            xi  = y / eta / sqrt(two)

            s   = ( xi - ximin ) / ( ximax - ximin )
            r   = pt5 * ( eta**2 - one )
            
            if (xflag.ne.0) &
              call NR_SPLINT( xx, ss, s2, npts, xi, s, tmp, tmp )
            if (yflag.ne.0) &
              call NR_SPLINT( rr, nn, n2, npts,  r, n, tmp, tmp )

!.... make sure that the computational coordinates are in range

            if (s.lt.zero) s = zero
            if (n.lt.zero) n = zero

!           write (*,"(7(1pe13.6,1x))") x, y, xi, eta, s, r, n
            
            if ( s .le. (one + 1.0e-8) ) then
              if (s .gt. one) s = one
              if (n .gt. one) n = one
#if defined(USE_IMSL) || defined(USE_BSLIB)
              v2(j,i,1) = BS2DR( 0, 0, n, s, kyord, kxord, yknot, &
                                 xknot, ny1, nx1, bsv(:,:,1) )
              v2(j,i,2) = BS2DR( 0, 0, n, s, kyord, kxord, yknot, &
                                 xknot, ny1, nx1, bsv(:,:,2) )
              v2(j,i,3) = BS2DR( 0, 0, n, s, kyord, kxord, yknot, &
                                 xknot, ny1, nx1, bsv(:,:,3) )
              v2(j,i,4) = BS2DR( 0, 0, n, s, kyord, kxord, yknot, &
                                 xknot, ny1, nx1, bsv(:,:,4) )
              v2(j,i,5) = BS2DR( 0, 0, n, s, kyord, kxord, yknot, &
                                 xknot, ny1, nx1, bsv(:,:,5) )
#else
              write(*,*) 'Must use IMSL'
              stop 1
#endif
            else
              write(*,"(2(i4,1x),4(1pe13.6,1x))") i,j,x,y,s,r
              v2(j,i,1) = one
              v2(j,i,2) = zero
              v2(j,i,3) = zero
              v2(j,i,4) = zero
              v2(j,i,5) = one
            end if
          end do
        end do              

!.... write out the new data file

        filen = 'output.R.0'
        write (*,"('Enter interpolated data filename [',a,']? ',$)") &
                    filen(1:index(filen,' '))
        read (*,"(a)") name
        if (name(1:1) .ne. ' ') filen = name
        open(unit=10, file=filen, form='unformatted', status='unknown')
        write(10) lstep, time, nx2, ny2, nz2, ndof, &
                  Re, Ma, Pr, gamma, cv
        if (useIJ) then
          write(10) (((v2(j,i,idof), idof=1,ndof), i=1,nx2), j=1,ny2)
        else
          write(10) v2
        endif
        close(10)
        write(*,*) "Finished writing data file"

        if (pot) then

!.... write out boundary values

        open(20,file='top.pot',form='formatted')
        do i = 1, nx2
          write(20,*) (v2(ny2,i,idof), idof = 1, ndof)
        end do
        close(20)

        open(20,file='right.pot',form='formatted')
        do j = 1, ny2
          write(20,*) (v2(j,nx2,idof), idof = 1, ndof)
        end do
        close(20)

!.... write out the pressure gradient field
!.... this is the gradient in the first generalized coordinate
!.... direction for the second mesh.
        
        allocate( p(ny2,nx2), g1p(ny2,nx2), g2p(ny2,nx2) )

        p = v2(:,:,1) * v2(:,:,5) / (gamma * Ma**2)
        
        call grad(1, nx2, ny2, p, g1p, g2p, dxi2, deta2, -1, -1, &
                  .false., .false., .false., .false., &
                  .false., .false., .false.)

        open(20,file='pg.pot',form='unformatted')
        write(20) g1p
        close(20)

        end if

!.... write out a plot3d file
        
        if (plot3d) then
          open(10, file='int.dat', form='unformatted', status='unknown')
          write(10) nx2, ny2, nz2
          write(10) Ma, Pr, Re, time
          if (pot) then
            write(10) ((((v2(j,i,idof), i = 1, nx2), j = 1, ny2),  &
                                        k = 1, nz2), idof = 1, 3), &
                      (((       p(j,i), i = 1, nx2), j = 1, ny2),  &
                                        k = 1, nz2), &
                      (((     g1p(j,i), i = 1, nx2), j = 1, ny2),  &
                                        k = 1, nz2)
          else
            write(10) ((((v2(j,i,idof),i=1,nx2),j=1,ny2),k=1,nz2),idof=1,5)
          endif
          close(10)
        end if
        
        stop
        end 

!=============================================================================!
        function arc(xi)
!
!  Compute the arclength along the parabolic cylinder
!
!=============================================================================!
        implicit none

        real :: arc, xi
!=============================================================================!

        arc = sqrt(2.0)*0.5*xi*sqrt(1.0+2.0*xi**2) + &
              0.5*log(sqrt(2.0)*xi + sqrt(1.0+2.0*xi**2))

        return
        end

!=============================================================================!
        function calcxi(s)
!=============================================================================!
        implicit none

        real :: calcxi, s
        real, external :: funcxi, zbrent

        real :: sloc
        common /xistuff/ sloc
!=============================================================================!
        sloc = s

        calcxi = zbrent(funcxi,0.0,2000.0,1.0e-12)

        return
        end

!=============================================================================!
        function funcxi(xi)
!=============================================================================!
        implicit none
        
        real funcxi, xi
        real, external :: arc

        real :: sloc
        common /xistuff/ sloc
!=============================================================================!
        funcxi = sloc - arc(xi)

        return
        end 
