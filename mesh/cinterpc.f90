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
        program cinter
!
!  This program interpolates from one computational mesh to another
!  using B-Spline basis function with the order set by the user.
!  Since this routine works in computational space, you cannot
!  change the mapping function between meshes.
!
!=============================================================================!
        use constants
        implicit none

!.... first grid and data

        integer :: nx1, ny1, nz1
        real, allocatable :: xy1(:,:,:), xi1(:), eta1(:)
        complex, allocatable :: v1(:,:,:)
        real :: dxi1, deta1
        
!.... second grid and data

        integer :: nx2, ny2, nz2
        real, allocatable :: xy2(:,:,:), xi2(:), eta2(:)
        complex, allocatable :: v2(:,:,:)
        real :: dxi2, deta2

!.... stuff for B-spline

        integer           :: kxord, kyord
        real, allocatable :: xknot(:), yknot(:)
        real, allocatable :: bsv(:,:,:)
        real, allocatable :: bsi(:,:,:)

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
#ifdef USE_IMSL
        real, external :: BS2DR
#endif
        real :: x, y, xi, eta, s, r, n

        real :: drmin, dsmin, b, rc, sc, cm, smin, smax

!.... argument parameters

        integer :: iarg, narg
        character(80) :: arg
        integer :: yflag=1, xflag=0
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
          case ('-h')
            write(*,"('Usage:  interpc [options]')")
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

        open (unit=10, file=filen, form='unformatted', status='old',err=200)
        read(10) nx1, ny1, nz1
        write(*,"('Nx = ',i4,', Ny = ',i4)") nx1, ny1
        allocate (xy1(ny1,nx1,2))
        read(10) (((xy1(j,i,1), i = 1, nx1), j = 1, ny1), k = 1, nz1), &
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
        open(unit=10, file=filen, form='unformatted', status='old',err=201)
        read(10) lstep, time, nx1, ny1, nz1, ndof, &
                 Re, Ma, Pr, gamma, cv
        read(10) v1
        close(10)

!.... read the second grid file

        filen = 'grid.dat'
 202    write (*,"('Enter grid to interpolate too [',a,']? ',$)") &
                    filen(1:index(filen,' '))
        read (*,"(a)") name
        if (name(1:1) .ne. ' ') filen = name

!.... allocate memory for the physical grid

        open (unit=10, file=filen, form='unformatted', status='old',err=202)
        read(10) nx2, ny2, nz2
        write(*,"('Nx = ',i4,', Ny = ',i4)") nx2, ny2
        allocate (xy2(ny2,nx2,2))
        read(10) (((xy2(j,i,1), i = 1, nx2), j = 1, ny2), k = 1, nz2), &
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

        allocate(v2(ny2,nx2,ndof))

        if(xy1(1,1,1) .le. zero) then
          ximin = zero
        else
          ximin = sqrt( xy1(1,1,1) )
        end if
        ximax = sqrt( xy1(1,nx1,1) )

!.... spline the tangent mapping function to get the inverse

        if (xflag.eq.1) then
          write(*,*) 'WARNING:  not corrected for nonzero origin'
          write(*,"('Enter dxmin ==> ',$)") 
          read(*,*) dxmin
          dss = one / real(npts-1)
          sx = calcs( ximax, dxmin, dxi1 )
          write(*,*) 'Sx = ', sx
          c2 = log( dxmin / dxi1 )
          c1 = (log( sx*dxmin / dxi1) - c2) / dxi1
          do i = 1, npts
            ss(i) = real(i-1) * dss
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
          dss = one / real(npts-1)
          do i = 1, npts
            ss(i) = real(i-1) * dss
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
          dss = one / real(npts-1)
          do i = 1, npts
            ss(i) = real(i-1) * dss
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
          dnn = one / real(npts-1)
          sy = calcs( rmax, dymin, deta1 )
          write(*,*) 'Sy = ', sy
          c2 = log( dymin / deta1 )
          c1 = (log( sy*dymin / deta1) - c2) / deta1
          do j = 1, npts
            nn(j) = real(j-1) * dnn
            rr(j) = one / c1 * ( exp(c1 * nn(j) + c2) - exp(c2) ) 
          end do
          call NR_SPLINE( rr, nn, npts, 1.0e31, 1.0e31, n2)
        else if (yflag.eq.2) then
          rd1  = 0.00028                ! set for R=1000 PCYL, r=500
          rd2  = 5.5
          dd   = calcdd(rd1,rd2)
          
          dnn = one / real(npts-1)
          do j = 1, npts
            nn(j) = real(j-1) * dnn
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
          
          dnn = one / real(npts-1)
          do j = 1, npts
            nn(j) = real(j-1) * dnn
            rr(j) = rmax*( cm * nn(j) + log( cosh(b*(nn(j)-rc)) / &
                    cosh(b*(nn(j)+rc)) ) ) / &
                    (cm + log( cosh(b*(one-rc)) / cosh(b*(one+rc)) ) )
          end do
          call NR_SPLINE( rr, nn, npts, 1.0e31, 1.0e31, n2)
        end if

!.... form the B-Spline Interpolant

        kxord = 5
!       write(*,"('Enter korder ==> ',$)")
!       read(*,*) kxord
        kyord = kxord

        allocate (xknot(nx1+kxord), yknot(ny1+kyord), &
                  bsv(ny1,nx1,ndof), bsi(ny1,nx1,ndof) )

        call BSNAK( nx1, xi1,  kxord, xknot)
        call BSNAK( ny1, eta1, kyord, yknot)
#ifdef USE_IMSL
        do idof = 1, ndof
          call BS2IN( ny1, eta1, nx1, xi1,  real(v1(:,:,idof)), ny1, &
                      kyord, kxord, yknot, xknot, bsv(:,:,idof))
          call BS2IN( ny1, eta1, nx1, xi1, aimag(v1(:,:,idof)), ny1, &
                      kyord, kxord, yknot, xknot, bsi(:,:,idof))
        end do
        write(*,*) 'B-spline complete...begin interpolating'
#else
        write(*,*) 'Must use IMSL currently'
        stop 1
#endif        

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
            
            if (xflag.ne.0) call NR_SPLINT( xx, ss, s2, npts, xi, s, tmp, tmp )
            if (yflag.ne.0) call NR_SPLINT( rr, nn, n2, npts,  r, n, tmp, tmp )

!.... make sure that the computational coordinates are in range

            if (s.lt.zero) s = zero
            if (n.lt.zero) n = zero

!           if (i.eq.nx2) write (*,"(7(1pe13.6,1x))") x, y, xi, eta, s, r, n
            
            if ( s .le. (one + 1.0e-8) ) then
              if (s .gt. one) s = one
              if (n .gt. one) n = one
#ifdef USE_IMSL
              v2(j,i,1) = cmplx(BS2DR( 0, 0, n, s, kyord, kxord,  &
                                yknot, xknot, ny1, nx1, bsv(:,:,1) ), &
                          BS2DR(0, 0, n, s, kyord, kxord, yknot, xknot,&
                                ny1, nx1, bsi(:,:,1) ) )
              v2(j,i,2) = cmplx(BS2DR( 0, 0, n, s, kyord, kxord, &
                                yknot, xknot, ny1, nx1, bsv(:,:,2) ), &
                          BS2DR(0, 0, n, s, kyord, kxord, yknot, xknot,&
                                ny1, nx1, bsi(:,:,2) ) )
              v2(j,i,3) = cmplx(BS2DR( 0, 0, n, s, kyord, kxord, &
                                yknot, xknot, ny1, nx1, bsv(:,:,3) ), &
                          BS2DR(0, 0, n, s, kyord, kxord, yknot, xknot,&
                                ny1, nx1, bsi(:,:,3) ) )
              v2(j,i,4) = cmplx(BS2DR( 0, 0, n, s, kyord, kxord, &
                                yknot, xknot, ny1, nx1, bsv(:,:,4) ), &
                          BS2DR(0, 0, n, s, kyord, kxord, yknot, xknot,&
                                ny1, nx1, bsi(:,:,4) ) )
              v2(j,i,5) = cmplx(BS2DR( 0, 0, n, s, kyord, kxord, &
                                yknot, xknot, ny1, nx1, bsv(:,:,5) ), &
                          BS2DR(0, 0, n, s, kyord, kxord, yknot, xknot,&
                                ny1, nx1, bsi(:,:,5) ) )
#else
              write(*,*) 'Must use IMSL currently'
              stop 1
#endif        
            else
              v2(j,i,:) = zero
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
        write(10) 0, 0, nx2, ny2, nz2, ndof, &
                  Re, Ma, Pr, gamma, cv
        write(10) v2
        close(10)

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
