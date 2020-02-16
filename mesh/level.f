c=============================================================================c
        program level
c=============================================================================c
c
c  Purpose:  Given a closed body, this program generates level sets
c            In has the option of performing a least squares B-spline to
c            smooth the input data.
c
c  Author:   S. Scott Collis
c
c  Revised:  6-17-97
c
c=============================================================================c
        implicit real (a-h,o-z)

        logical IRIS
        parameter (IRIS = .false.)

        parameter (zero=0.0d0,   pt5=0.5d0,    pt25=0.25d0,     one=1.0d0, 
     &             onept5=1.5d0, two=2.0d0,    twopt5=2.5d0, 
     &             three=3.0d0,  BIG = 1.0e10, pi = 3.1415926535897932385d+0)

        parameter (iout=10, iin=11, ibs=12)

        parameter (mpts = 1024)
        parameter (nsd = 2)

        dimension  b(mpts,nsd), r(nsd), dr(nsd)
        dimension  bigb(mpts,nsd), bigt(mpts,nsd), bigs(mpts)
        dimension  arcl(mpts)
        dimension  x1(mpts), x2(mpts)
        dimension  rr(mpts), ss(mpts), zz(mpts), cur(mpts)
        dimension  xi(mpts), eta(mpts)
c
c.... needed for plot3d files
c
        real xx(mpts,mpts),  yy(mpts,mpts)
c
c.... B-splines data structure
c
        parameter (kmax = 8)
        dimension xknot(kmax+mpts), yknot(kmax+mpts)
        dimension s(mpts), t(mpts,nsd)
        common /bspline/ npts, t, darc, s1, s2, s, nknot, xknot, yknot, 
     &                   korder, ncoef
c
c.... external functions
c
        external crvarc
c
c.... argument list parameters
c
        integer iarg, narg
        character*80  arg
        integer yflag, xflag
        logical plot3d, debug, period, sym, ellip, field
#ifndef __GFORTRAN__
        integer, external iargc
#endif
c=============================================================================c
        yflag = 0
        xflag = 0
        plot3d = .false.
        debug = .false.
        period = .false.
        sym = .false.
        ellip = .false.
        field = .false.
c
c.... parse the argument list
c
        narg = iargc()
        do iarg = 1, narg
          call getarg(iarg,arg)
          if (arg(1:2) .eq. '-y') then
            if (arg(1:3) .eq. '-y1') then               ! exponential
              yflag = 1
            else if (arg(1:3) .eq. '-y2') then          ! hyperbolic tangent
              yflag = 2
            else if (arg(1:3) .eq. '-y3') then          ! Mahesh
              yflag = 3
            else if (arg(1:3) .eq. '-y4') then          ! Algebraic
              yflag = 4
            else                                        ! uniform
              yflag = 0
            end if
          else if (arg(1:2) .eq. '-x') then
            if (arg(1:3) .eq. '-x1') then               ! exponential
              xflag = 1
            else if (arg(1:3) .eq. '-x2') then          ! hyperbolic tangent
              xflag = 2
            else if (arg(1:3) .eq. '-x3') then          ! Mahesh in xi
              xflag = 3
            else if (arg(1:3) .eq. '-x4') then          ! Mahesh in s
              xflag = 4
            else                                        ! uniform in xi
              xflag = 0
            end if
          else if (arg(1:2) .eq. '-c') then
            period = .true.
          else if (arg(1:2) .eq. '-d') then
            debug = .true.
          else if (arg(1:2) .eq. '-e') then
            ellip = .true.
          else if (arg(1:2) .eq. '-f') then
            field = .true.
          else if (arg(1:2) .eq. '-p') then
            plot3d = .true.
          else if (arg(1:2) .eq. '-s') then
            sym = .true.
          else if (arg(1:2) .eq. '-h') then
            write(*,"('Usage:  level [options]')")
            write(*,"('   -h:  this help')")
            write(*,"('   -d:  debug output')")
            write(*,"('   -e:  use the elliptic smoother')")
            write(*,"('   -f:  write field file format')")
            write(*,"('   -c:  Harwire periodicity')")
            write(*,"('   -s:  Hardwire LE symmetry')")
            write(*,"('  -x1:  exponential stretching in x')")
            write(*,"('  -x2:  hyperbolic tangent stretching in x')")
            write(*,"('  -x3:  Mahesh''s mapping in x')")
            write(*,"('  -x4:  Mahesh''s mapping in s')")
            write(*,"('  -y1:  exponential stretching in y')")
            write(*,"('  -y2:  hyperbolic tangent stretching in y')")
            write(*,"('  -y3:  Mahesh''s mapping in y')")
            write(*,"('  -y4:  Algebraic mapping in y')")
            call exit(0)
          else
            write(*,"('Argument ',i2,' ignored.')") iarg
          end if
        end do
c
c.... get grid dimensions
c
        write(*,"(/,'Enter Nx, Ny ==> ',$)")
        read(*,*) nx, ny
        if (nx .gt. mpts .or. ny .gt. mpts) then
          write(*,*) 'Increase mpts... ', mpts
          call exit(1)
        endif
c
c.... B-spline should be at least 4th order
c
        korder = 5
c       write(*,"('Enter korder ==> ',$)")
c       read(*,*) korder
        if (korder.gt.kmax) then
          write(*,"('Error: Increase kmax')")
          call exit(1)
        end if

c=============================================================================
c       A p p r o x i m a t e   t h e   B o d y 
c=============================================================================
c
c.... open the input file
c       
        open(iin, file='body.dat', status='unknown')
c
c.... read the boundary file
c
        read (iin, *) npts
        if (npts.gt.mpts) then
          write(*,*) 'Increase mpts in level', mpts
          call exit(1)
        end if
        do i = 1, npts
          read(iin,*) b(i,1), b(i,2)
        end do
        close(iin)
c
c.... initialize the arclength using a linear approximation
c
        arcmin = BIG
        s(1) = zero
        do i = 2, npts
           arc = sqrt((b(i,1)-b(i-1,1))**2 + (b(i,2)-b(i-1,2))**2)
           arcmin = min(arc,arcmin)
           s(i) = s(i-1) + arc
c          write (14,100) s(i), b(i,1), b(i,2)
        end do
        arctot = s(npts)
        write (*,"(/,'Linear Estimation')")
        write (*,"('Minimum arclength on body = ', 1pe20.13)") arcmin
        write (*,"('Total arclength on body   = ', 1pe20.13)") arctot
c
c.... Compute a smooth B-spine
c
        ncoef = npts
        write(*,"(/,'Least squares spline')")
        write(*,"('Npts = ',i4,'  Enter Ncoef ( <= Npts ) ==> ',$)") npts
        read(*,*) ncoef
        nknot = korder + ncoef
        call splcrv(npts, s, b(1,1), b(1,2), t(1,1), t(1,2), 
     &              nknot, xknot, yknot, korder, ncoef)
c
c.... compute the total arclength of the B-splined body
c
        arctot = crvarc(npts, s, t(1,1), t(1,2), nknot, xknot, yknot, 
     &                  korder, ncoef, arcl)
        write (*,"(/,'After smoothing Spline')")
        write (*,"('Total arclength on body = ', 1pe20.13)") arctot

        if (.false.) then
c
c.... Oversample the smooth spline
c.... It turns out that you don't have to oversample
c
        itime = 3
        if (itime*(npts-1)+1 .gt. mpts) then
           write(*,"('ERROR: increase mpts in level')")
           call exit(1)
        end if
        do i = 1, npts-1
           i0 = (i-1)*itime+1
           ds = (s(i+1)-s(i))/float(itime)
           do j = 0, itime-1
              sint = s(i) + float(j) * ds
              bigs(i0+j) = sint
              b(i0+j,1)  = BSDER( 0, sint, korder, xknot, ncoef, t(1,1) )
              b(i0+j,2)  = BSDER( 0, sint, korder, yknot, ncoef, t(1,2) )
           end do
        end do
        i0 = itime * (npts-1) + 1
        sint = s(npts)
        bigs(i0) = sint
        b(i0,1)  = BSDER( 0, sint, korder, xknot, ncoef, t(1,1) )
        b(i0,2)  = BSDER( 0, sint, korder, yknot, ncoef, t(1,2) )
        npts = itime * (npts-1) + 1
        do i = 1, npts
           s(i) = bigs(i)
        end do
c
c.... With the oversampled smooth data iterate to find the spline based on
c.... arclength.
c
        icnt = 0
        ncoef = npts
        nknot = korder + ncoef
 200    icnt = icnt + 1
        call splcrv(npts, s, b(1,1), b(1,2), t(1,1), t(1,2), 
     &              nknot, xknot, yknot, korder, ncoef)
        arcold = arctot
        arctot = crvarc(npts, s, t(1,1), t(1,2), nknot, xknot, yknot, 
     &                  korder, ncoef, arcl)
c
c.... find the max error
c
        errmax = 0.0
        do i = 2, npts
          err = abs( (arcl(i)-arcl(i-1)) - (s(i)-s(i-1)) )
          errmax = max(err,errmax)
        end do
        imax = 10
        if ( errmax .gt. 1.0e-10 .and. icnt .le. imax ) then
           write(*,*) icnt, errmax
           do i = 1, npts
              s(i) = arcl(i)
           end do
           goto 200
        end if

        write (*,"(/,'After re-Spline')")
        write (*,"('Total arclength on body = ', 1pe20.13)") arctot
        write (*,"('Converged to ', 1pe20.13)") abs(arctot-arcold)
c
c.... end of oversample
c
        end if
c
c.... write out the body as a check
c
        if (debug) then
          do i = 1, npts
            write (15,100) s(i), b(i,1), b(i,2)
          end do
        end if
c
c.... save the B-spline representation
c
        open(unit=ibs, file='bs.dat', status='unknown')
        write(ibs,40) korder, ncoef
        do i = 1, ncoef + korder
           write(ibs,50) xknot(i), yknot(i)
        end do
        do i = 1, ncoef
           write(ibs,50) t(i,1), t(i,2)
        end do
 40     format(2(i4,1x))
 50     format(2(1pe22.15,1x))
        close(ibs)

c=============================================================================
c       B o d y   T a n g e n t   D i r e c t i o n
c=============================================================================

        smin = zero
        smax = arctot
        savg = (smin + smax) * pt5
        write(*,"(/,'Smin = ',1pe13.6,' Savg = ',1pe13.6,' Smax = ',1pe13.6)")
     &        smin, savg, smax
c
c.... arbitrary mapping function or uniform mesh
c
        dxi = one / real(nx-1)
        ds = (smax-smin) / real(nx-1)
        write(*,"(/,'Uniform dS = ',1pe20.13)") ds
c
c.... Uniform arc-length
c
        if (xflag.eq.0) then
          write(*,"(/,'Uniform grid in S')")
          do i = 1, nx
            xi(i) = real(i-1) * dxi
            ss(i) = smin + real(i-1) * ds
          end do
        end if

        if (xflag.eq.1) then
           write(*,"('ERROR:  xflag = 1 is not supported')")
           call exit(1)
        end if

        if (xflag.eq.2) then
           write(*,"('ERROR:  xflag = 2 is not supported')")
           call exit(1)
        end if

        if (xflag.eq.3) then
           write(*,"('ERROR:  xflag = 3 is not supported')")
           call exit(1)
        end if
c
c.... Uniform arc-length + mapping
c
        if (xflag.eq.4) then
          write(*,"(/,'Mahesh map in S')")
          write(*,"('Enter dsmin, b, sc ==> ',$)") 
          read(*,*) dsmin, rb, sc
          cm = ( two * rb * tanh(rb*sc) + real(nx-1)*dsmin/smax * 
     &         log( cosh(rb*(one-sc)) / cosh(rb*(one+sc)) ) ) /
     &         ( one - real(nx-1)*dsmin/smax )
          do i = 1, nx
            xi(i)   = dble(i-1) * dxi
            ss(i)   = smax*( cm * xi(i) + log( cosh(rb*(xi(i)-sc)) /
     &                cosh(rb*(xi(i)+sc)) ) ) / 
     &                (cm + log( cosh(rb*(one-sc)) / cosh(rb*(one+sc)) ) )
          end do
        end if
c
c.... spline the mapping function for the body tangent dirction
c
c       ds = two / (500-1)
c       c  = 0.3333333333333333d0
c       c1 = 4.0d0
c       c2 = 20.0d0
c       c3 = 0.0
c       do i = 1, 500
c         s1(i) = -one + (i-1) * ds
c         f = pt5 * c1 * tanh(c2*(s1(i)-c3))
c         zz(i) = c * ( s1(i) + f + 1.0 - pt5 * c1 * tanh(c2*(-1-c3))) - 1.0
c       end do
c       call spline(500, zz, s1, s2)
c       do i = 1, nx
c         xi(i) = real(i-1) * dxi
c         call speval(500, zz, s1, s2, xi(i), ss(i))
c         ss(i) = ss(i) * savg + savg
c         write(19,110) i, xi(i), ss(i)
c       end do

c=============================================================================
c       B o d y   N o r m a l   D i r e c t i o n
c=============================================================================

        deta = one / real(ny-1)
        rmin = zero
        write(*,"(/,'Enter Rmax ==> ',$)")
        read(*,*) rmax
        drr = (rmax-rmin) / real(ny-1)
        write(*,"(/,'Uniform dR = ',1pe20.13)") drr
c
c.... uniform mesh in the wall normal direction
c
        if (yflag.eq.0) then
          write(*,"(/,'Uniform grid in R')")
          do j = 1, ny
            eta(j) = real(j-1) * deta
            rr(j)  = rmin + dble(j-1) * drr
c           write(22,110) j, rr(j), eta(j)
          end do
        end if
c
c.... exponential map in the wall normal direction
c
        if (yflag.eq.1) then
          write(*,"(/,'Exponential map in R')")
          deta2 = one / real(128-1)
          dymin = 5.0d-5
          sy = 1.088d0
          c2 = log( dymin / deta2 )
          c1 = (log( sy*dymin / deta2) - c2) / deta2
          do j = 1, ny
            eta(j)  = dble(j-1) * deta
            rr(j)   = one / c1 * ( exp(c1 * eta(j) + c2) - exp(c2) ) 
c           write(22,110) j, rr(j), eta(j)
          end do
        end if
c
c.... Hyperbolic tangent stretching
c
        if (yflag .eq. 2) then
          write(*,"(/,'Hyperbolic stretching in R')")
          ds1 = 0.1                     ! set for R=1000 cylinder
          ds2 = 2.0
          dd  = 2.36921620400439        ! from Mathematica (grid2.ma)
c         ds1 = 0.01                    ! set for R=2400 MSE, r=26
c         ds2 = 2.5
c         dd  = 3.898620676884022               
c         ds1 = 0.005                   ! set for R=2400 MSE, r=50
c         ds2 = 5.0
c         dd  = 3.898620676884022               
c         ds1 = 0.0005                  ! set for R=2400 MSE, r=50
c         ds2 = 5.0
c         dd  = 5.36966703089523                
          do j = 1, ny
            eta(j) = float(j-1) * deta
            rr(j) =  rmax*(0.5 + 1.0/Tanh(dd/2.0)*Tanh(dd*(-0.5 +
     &               eta(j)))/2.0)/(Sqrt(ds2/ds1) +
     &               (1.0 - Sqrt(ds2/ds1))*
     &               (0.5 + 1.0/Tanh(dd/2.0)*Tanh(dd*(-0.5 +
     &               eta(j)))/2.0))
c           write(22,110) j, rr(j), eta(j)
          end do
        end if
c
c.... Mahesh's mapping, personal communication (6-5-96)
c
        if (yflag.eq.3) then
          write(*,"(/,'Mahesh map in R')")
          write(*,"('Enter drmin, b, rc ==> ',$)") 
          read(*,*) drmin, rb, rc
          cm = ( two * rb * tanh(rb*rc) + (ny-1)*drmin/rmax * 
     &           log( cosh(rb*(one-rc)) / cosh(rb*(one+rc)) ) ) /
     &         ( one - (ny-1)*drmin/rmax )
          do j = 1, ny
            eta(j)  = real(j-1) * deta
            rr(j)   = rmax*( cm * eta(j) + log( cosh(rb*(eta(j)-rc)) /
     &                cosh(rb*(eta(j)+rc)) ) ) / 
     &                (cm + log( cosh(rb*(one-rc)) / cosh(rb*(one+rc)) ) )
c           write(22,110) j, rr(j), eta(j)
          end do
        end if
c
c.... Algebraic stretched mesh
c
        if (yflag.eq.4) then
          write(*,"(/,'Algebraic grid in R')")
          write(*,"('Enter ri ==> ',$)")
          read(*,*) ri
          aa = rmax * ri / ( rmax - two * ri )
          bb = one + aa / rmax
          do j = 1, ny
            eta(j) = (j-1) * deta
            rr(j)   = aa * eta(j) / (bb - eta(j))
c           write(22,110) j, rr(j), eta(j)
 110        format(i5,6(1x,1pe13.6))
          end do
        end if

c=============================================================================
c       F o r m   L e v e l   S e t s
c=============================================================================
        write(*,"(/,'Discretizing boundary...')")
        do i = 1, nx
          if (i.eq.1) then
            s1 = s(1)
            sint = 0.0
            nold = 1
          else
            s1 = sint
            darc = ss(i) - ss(i-1)
            do n = nold, npts
               if ( arcl(n) .ge. ss(i) ) then
                 s2 = s(n)
                 nold = n
                 goto 11
              end if
            end do
 11         sint = gets()
          end if
          r(1)  = BSDER( 0, sint, korder, xknot, ncoef, t(1,1) )
          r(2)  = BSDER( 0, sint, korder, yknot, ncoef, t(1,2) )
          dr(1) = BSDER( 1, sint, korder, xknot, ncoef, t(1,1) )
          dr(2) = BSDER( 1, sint, korder, yknot, ncoef, t(1,2) )
          do k = 1, nsd
            bigb(i,k) = r(k)
            bigt(i,k) = dr(k)
          end do
          bigs(i) = ss(i)
        end do
c
c.... periodicity
c
        if (period) then
           write(*,"(/,'Hardwire periodicity (1,0) ==> ',$)")
           do k = 1, nsd
              bigb(nx,k) = bigb(1,k)
              bigt(nx,k) = bigt(1,k)
           end do
        end if
c       
c.... symmetry
c
        if (sym) then
           write(*,"(/,'Hardwire LE symmetry (1,0) ==> ')")
           bigb(1,1) = zero
           bigb(1,2) = zero
           bigt(1,1) = zero
           bigt(1,2) = one
        end if
c
        nptsbig = nx
c
c.... exchange the curves
c
        do i = 1, nptsbig
          b(i,1) = bigb(i,1)
          b(i,2) = bigb(i,2)
          t(i,1) = bigt(i,1)
          t(i,2) = bigt(i,2)
          s(i)   = bigs(i)
        end do
        npts = nptsbig
c
c.... find the minimum arc length
c.... NB: s(:) is really the local arclength on the body now
c
        arcmin = BIG
        arcl(1) = zero
        do j = 2, npts
          arc    = s(j) - s(j-1)
          arcmin = min(arc,arcmin)
          arcl(j) = arcl(j-1) + arc
        end do
        write (*,"(/,'Discretized Body')")
        write (*,"('Minimum arclength on body = ', 1pe20.13)") arcmin
        write (*,"('Total   arclength on body = ', 1pe20.13)") arcl(npts)
c
c.... reset nx just in case
c
        nx = npts
c
c.... compute the outward normal using the spline data
c
        do i = 1, npts
          x1(i) = -t(i,2)
          x2(i) =  t(i,1) 
          rx    =  sqrt(x1(i)**2 + x2(i)**2)
          x1(i) =  x1(i) / rx
          x2(i) =  x2(i) / rx
          if (i .gt. 2) then
            if ( (x2(i) .gt. zero .and. x2(i-2) .lt. zero) .or.
     &           (x1(i) .eq. -one) ) then
              write(*,"(/,'Estimate LE at',i4,3(1x,1pe13.6))")
     &           i-1, b(i-1,1), b(i-1,2), s(i-1)
            end if
          end if
        end do
c
c.... Compute a B-spline through the boundary normal vectors
c
        call BSNAK(npts, s, korder, xknot)
        call BSNAK(npts, s, korder, yknot)
        call BSINT(npts, s, x1, korder, xknot, t(1,1) )
        call BSINT(npts, s, x2, korder, yknot, t(1,2) )
c
c.... write out the body, the surface normal vector, and the curvature
c
        open(20,file='body.out')
        open(21,file='curv.out')
        do i = 1, npts
           write (20,100) b(i,1), b(i,2)
           dr(1) = BSDER( 1, s(i), korder, xknot, npts, t(1,1) )
           dr(2) = BSDER( 1, s(i), korder, yknot, npts, t(1,2) )
           cur(i) = sqrt( dr(1)**2 + dr(2)**2 )
           write (21,100) s(i), b(i,1), b(i,2), x1(i), x2(i), cur(i)
        end do
        close(20)
        close(21)
c
c.... conformal mapping
c
        call conformal( npts, s, b(1,1), b(1,2), cur )
c
c.... make the level sets
c
        do i = 1, npts
          xx(i,1) = b(i,1)
          yy(i,1) = b(i,2)
        end do

        do j = 2, ny
          dy = rr(j) - rr(j-1)
          do i = 1, npts
            b(i,1) = b(i,1) + x1(i) * dy
            b(i,2) = b(i,2) + x2(i) * dy
            xx(i,j) = b(i,1)
            yy(i,j) = b(i,2)
          end do
c
c.... output each level to a separate file
c
          if (.false.) then
            do i = 1, npts
              write (29+iy-1,*) b(i,1), b(i,2)
            end do
          end if
        end do
c
c.... try a simple elliptic smoother
c
        if (ellip) then
           write(*,"(/,'Using elliptic smoother...')")
           write(*,"('Enter eps, nmax ==> ',$)")
           read(*,*) eps, nmax
           call elliptic( nx, ny, xx, mpts, yy, mpts, eps, nmax)
        end if

c=============================================================================
c       O u t p u t   R e s u l t s
c=============================================================================
c
c.... output in Streett's format
c
        if (field) then
           open(10,file='field.mean',form='unformatted')
           write(10) nx, ny
           write(10) (((0.0, j=1, ny), i=1,nx), k=1,4)
           write(10) (s(i), i=1,nx)
           write(10) (rr(j), j=1,ny)
           write(10) (cur(i), i=1,nx)
           close(10)
        end if
c
c.... output grid
c
        nz = 1
        write(*,"(/,'Writing plot3d files for (',i4,',',i4,',',i4,')')") 
     &    nx, ny, nz
c
c.... write out binary plot3d files if on the iris
c
        if (IRIS) then

        call wgrid(xx, yy, nx, ny, nz, mpts, mpts, 'grid.bin'//char(0))
        call wdata(xx, yy, nx, ny, nz, mpts, mpts, 'q.bin'//char(0))

        if (plot3d) then
c
c.... body fitted coordinates
c
        do i = 1, nx
          do j = 1, ny
            xx(i,j) = ss(i)
            yy(i,j) = rr(j)
          end do
        end do
        call wgrid(xx, yy, nx, ny, nz, mpts, mpts, 'grid1.bin'//char(0))
c
c.... computational coordinats
c
        do i = 1, nx
          do j = 1, ny
            xx(i,j) = xi(i)
            yy(i,j) = eta(j)
          end do
        end do
        call wgrid(xx, yy, nx, ny, nz, mpts, mpts, 'grid2.bin'//char(0))

        end if

        else
c
c.... write out unformatted plot3d files
c
        open(unit=iout,file='grid.dat',form='unformatted')
        write(iout) nx, ny, 1
        write(iout) (((xx(i,j), i=1,nx), j=1,ny), k = 1, nz),
     &              (((yy(i,j), i=1,nx), j=1,ny), k = 1, nz),
     &              (((    0.0, i=1,nx), j=1,ny), k = 1, nz)
        close(iout)

        if (plot3d) then

        open(unit=iout,file='q.dat',form='unformatted')
        write(iout) nx, ny, 1
        write(iout) 0.0, 0.0, 0.0, 0.0
        write(iout) (((    1.0, i=1,nx), j=1,ny), k = 1, nz),
     &              (((xx(i,j), i=1,nx), j=1,ny), k = 1, nz),
     &              (((yy(i,j), i=1,nx), j=1,ny), k = 1, nz),
     &              (((    0.0, i=1,nx), j=1,ny), k = 1, nz),
     &              (((    1.0, i=1,nx), j=1,ny), k = 1, nz)
        close(iout)
c
c.... body fitted coordinates
c
        do i = 1, nx
          do j = 1, ny
            xx(i,j) = ss(i)
            yy(i,j) = rr(j)
          end do
        end do
        open(unit=iout,file='grid1.dat',form='unformatted')
        write(iout) nx, ny, 1
        write(iout) (((xx(i,j), i=1,nx), j=1,ny), k = 1, nz),
     &              (((yy(i,j), i=1,nx), j=1,ny), k = 1, nz),
     &              (((    0.0, i=1,nx), j=1,ny), k = 1, nz)
        close(iout)
c
c.... computational coordinates
c
        do i = 1, nx
          do j = 1, ny
            xx(i,j) = xi(i)
            yy(i,j) = eta(j)
          end do
        end do
        open(unit=iout,file='grid2.dat',form='unformatted')
        write(iout) nx, ny, 1
        write(iout) (((xx(i,j), i=1,nx), j=1,ny), k = 1, nz),
     &              (((yy(i,j), i=1,nx), j=1,ny), k = 1, nz),
     &              (((    0.0, i=1,nx), j=1,ny), k = 1, nz)
        close(iout)

        endif
c
        endif
c
        call exit(0)
 100    format(8(1pe13.6,1x))
        end
