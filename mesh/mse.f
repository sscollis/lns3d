c=============================================================================c
        program MSE
c=============================================================================c
c
c  Purpose:  make a body-fitted Super-ellipse mesh
c
c  Author:   S. Scott Collis
c
c  Revised:  9-6-96
c  Revised:  3-11-2020
c
c=============================================================================c
c       implicit double precision (a-h,o-z)

        logical IRIS, DIAG
        parameter (IRIS = .false.)
        parameter (DIAG = .false.)

        parameter (zero=0.0d0,    pt5=0.5d0,  one=1.0d0, 
     &             onept5=1.5d0,  two=2.0d0,  twopt5=2.5d0, 
     &             three=3.0d0)
c
c.... must use a large mpts to make error in spline very small
c
        parameter (mpts=50000)
        dimension zz(mpts), s1(mpts), s2(mpts), ds1(mpts), ds2(mpts),
     &            dss1(mpts), dss2(mpts)

        parameter (mx=512,my=512)
        
        dimension xb(mx), yb(mx)
        dimension bn1(mx), bn2(mx)
        
        dimension xi(mx),  eta(my)

        dimension ss(mx), rr(my), dss(mx), drr(my), d2ss(mx), d2rr(my)

        dimension x(mx,my),    y(mx,my)
        dimension rm1(mx,my),  rm2(mx,my)
        dimension rn1(mx,my),  rn2(mx,my)
        dimension rjac(mx,my)
        dimension rm11(mx,my), rm12(mx,my), rm21(mx,my), rm22(mx,my)
        dimension rn11(mx,my), rn12(mx,my), rn21(mx,my), rn22(mx,my)
        dimension g11x(mx,my), g11y(mx,my), g12x(mx,my), g12y(mx,my)
        dimension g22x(mx,my), g22y(mx,my)
        
        dimension d(10), e(10), f1(10), f2(10), g0(10)
        
        common /stuff/ AR, rm, rn, xmin, xmax

        external xloc, yloc, calcs

c.... fourth order central difference ( 1 2 x 4 5 )

        parameter ( ga1 =  8.333333333333333333333E-02 )
        parameter ( ga2 = -6.666666666666666666667E-01 )
        parameter ( ga3 =  6.666666666666666666667E-01 )
        parameter ( ga4 = -8.333333333333333333333E-02 )

c.... fourth order one-sided ( x 2 3 4 5 )

        parameter ( gc1 = -2.083333333333333333333E+00 )
        parameter ( gc2 =  4.000000000000000000000E+00 )
        parameter ( gc3 = -3.000000000000000000000E+00 )
        parameter ( gc4 =  1.333333333333333333333E+00 )
        parameter ( gc5 = -2.500000000000000000000E-01 )

c.... fourth order biased difference ( 1 x 2 3 4 5 )

        parameter ( gb1 = -2.500000000000000000000E-01 )
        parameter ( gb2 = -8.333333333333333333333E-01 )
        parameter ( gb3 =  1.500000000000000000000E+00 )
        parameter ( gb4 = -5.000000000000000000000E-01 )
        parameter ( gb5 =  8.333333333333333333333E-02 )

        integer narg, iarg
        character(80) arg

        logical metric_ji
        integer xflag, yflag
c=============================================================================c
        metric_ji = .false.
        xflag = 1
        yflag = 1
c
c.... parse the argument list
c
        narg = iargc()
        do iarg = 1, narg
          call getarg(iarg,arg)
          select case (arg(1:3))
          case ('-ms')
            metric_ji = .true.
          case ('-ji')
            metric_ji = .true.
          case ('-x1')           ! exponential (default)
            xflag = 1
c         case ('-x2')           ! hyperbolic tangent
c           xflag = 2
          case ('-x3')           ! Mahesh in xi
            xflag = 3
c         case  ('-x4')          ! Mahesh in s
c           xflag = 4
          case ('-y1')           ! exponential (default)
            yflag = 1
          case ('-y2')           ! hyperbolic tangent
            yflag = 2
c         case ('-y3')           ! Mahesh
c           yflag = 3
          case ('-h')
            write(*,"('-----------------------------------------------')")
            write(*,"('Usage:  mse [options] ')")
            write(*,"('-----------------------------------------------')")
            write(*,"('   -h:  this help')")
            write(*,"('-----------------------------------------------')")
            write(*,"('  -ms:  write metrics in ji format')")
            write(*,"('  -ji:  write metrics in ji format')")
            write(*,"('  -x1:  exponential stretching in x (default)')")
c           write(*,"('  -x2:  hyperbolic tangent map in x')")
            write(*,"('  -x3:  Mahesh''s mapping in x')")
c           write(*,"('  -x4:  Mahesh''s mapping in s')")
            write(*,"('  -y1:  exponential stretching in y (default)')")
            write(*,"('  -y2:  hyperbolic tangent map in y')")
c           write(*,"('  -y3:  Mahesh''s mapping in y')")
            write(*,"('-----------------------------------------------')")
            call exit(0)
          case default
            write(*,"('Argument ',i2,' ignored.')") iarg
          end select
        end do
c
c.... keyboard input
c
        write(*,"('Enter AR (often 6) ==> ',$)")
        read(*,*) AR
        xmin = zero
        rmin = zero
        write(*,"('Enter xmax, rmax ==> ',$)")
        read(*,*) xmax, rmax
        write(*,"('Enter m (usually 4) ==> ',$)")
        read(*,*) rm
        rn = two
        write(*,"('Enter Nx, Ny ==> ',$)")
        read(*,*) nx, ny
c
c.... discretize in the body tangent direction
c
        smin = zero
        smax = arc(xmin,xmax)
        savg = (smin + smax) * pt5
        write(*,10) xmin, xmax, smax
        if (nx.gt.mx) stop 'nx > mx'
        dxi = one / real(nx-1)
        ds  = (smax-smin) / real(nx-1)

        write(*,*) 'Smax = ',smax,', xloc(smax) = ',xloc(xmin,smax)
c
c.... exponential map in the streamwise direction
c
        if (xflag.eq.1) then
          write(*,"('Enter dsmin ==> ',$)") 
          read(*,*) dsmin
          sx = calcs( smin, smax, dsmin, dxi )
          write(*,*) 'Sx = ', sx
          c2 = log( dsmin / dxi )
          c1 = (log( sx * dsmin / dxi) - c2) / dxi
          do i = 1, nx
            xi(i)   = real(i-1) * dxi
            ss(i)   = smin + one / c1 * ( exp(c1 * xi(i) + c2) - exp(c2) )
            dss(i)  = exp(c1 * xi(i) + c2)
            d2ss(i) = c1 * exp(c1 * xi(i) + c2)
            if (i.ne.1) then
              write(10,30) xi(i), ss(i), dss(i), d2ss(i),
     &                     (ss(i)-ss(i-1))/dxi
            else
              write(10,30) xi(i), ss(i), dss(i), d2ss(i), zero
            end if
          end do
        end if
c
c.... Mahesh mapping, personal communication (6-5-96)
c
        if (xflag.eq.3) then
          write(*,"('Enter dsmin, b, sc ==> ',$)") 
          read(*,*) dsmin, b, sc
          cm = ( two * b * tanh(b*sc) + (nx-1)*dsmin/smax * 
     &           log( cosh(b*(one-sc)) / cosh(b*(one+sc)) ) ) /
     &         ( one - (nx-1)*dsmin/smax )
          do i = 1, nx
            xi(i)   = dble(i-1) * dxi
            ss(i)   = smin + smax*( cm * xi(i) + log( cosh(b*(xi(i)-sc)) /
     &                cosh(b*(xi(i)+sc)) ) ) / 
     &                (cm + log( cosh(b*(one-sc)) / cosh(b*(one+sc)) ) )
            dss(i)  = smax*(cm + b*tanh(b*(xi(i)-sc)) - b*tanh(b*(xi(i)+sc)))/
     &                (cm + log( cosh(b*(one-sc)) / cosh(b*(one+sc)) ) )
            d2ss(i) = smax*(-b**2*(tanh(b*(xi(i)-sc)))**2 + 
     &                b**2*(tanh(b*(xi(i)+sc)))**2)/
     &                (cm + log( cosh(b*(one-sc)) / cosh(b*(one+sc)) ) )
          end do
        end if
c
c.... spline the tangent mapping function
c
        if (.false.) then
          npts = mpts
          if (npts.gt.mpts) stop 'npts > mpsts'
          ds = two / real(npts-1)
          c  = 0.333333333333333333333d0
          c1 = 4.0d0
          c2 = 20.0d0
          c3 = zero
          do i = 1, npts
            s1(i) = -one + real(i-1) * ds
            f = pt5 * c1 * tanh(c2*(s1(i)-c3))
            zz(i) = c * ( s1(i) + f + one - pt5 * c1 * 
     &            tanh( c2*(-one-c3) ) ) - one
          end do
        end if
c
c.... spline the tangent mapping function
c
        if (.false.) then
          write(*,*) 'Discretize in the wall tangent direction'
          npts = mpts
          if (npts.gt.mpts) stop 'npts > mpsts'
          ds = one / real(npts-1)
  
          c  = 0.50251256281407d0       ! WARNING: Must be updated
  
          c1 = 1.0d0
          c2 = 150.0d0
          c3 = zero
          f0 = pt5 * c1 * tanh(c2*(zero-c3))
  
          d(1)  = 50.0d0
          e(1)  = 30.0d0
          f1(1) = -0.01d0
          f2(1) = 0.01d0
  
          d(2)  = 1.0d0
          e(2)  = 20.0d0
          f1(2) = 0.12d0
          f2(2) = 0.17d0
  
          d(3)  = 1.0d0
          e(3)  = 20.0d0
          f1(3) = 0.03d0
          f2(3) = 0.10d0
  
          kmax  = 3
          do k = 1, kmax
            g0(k) = (d(k)-one)/(two*e(k)) * log( Cosh(e(k)*(zero-f1(k)))/
     &            Cosh(e(k)*(zero-f2(k))) )
          end do
          
          do i = 1, npts
            s1(i) = zero + real(i-1) * ds
            f = pt5 * c1 * tanh(c2*(s1(i)-c3))
            zz(i) = s1(i) + f - f0
            do k = 1, kmax
              g = (d(k)-one)/(two*e(k)) * log( Cosh(e(k)*(s1(i)-f1(k)))/
     &           Cosh(e(k)*(s1(i)-f2(k))) )
              zz(i) = zz(i) + g - g0(k)
            end do
            zz(i) = c * ( zz(i) )
          end do
c
c.... compute the derivative of the mapping function and spline it
c
          ds1(1) = ( gc1 * zz(1) + gc2 * zz(2) + gc3 * zz(3) +
     &               gc4 * zz(4) + gc5 * zz(5) ) / ds
          ds1(2) = ( gb1 * zz(1) + gb2 * zz(2) + gb3 * zz(3) +
     &               gb4 * zz(4) + gb5 * zz(5) ) / ds
          do i = 3, npts-2
             ds1(i) = ( ga1 * zz(i-2) + ga2 * zz(i-1) + 
     &                  ga3 * zz(i+1) + ga4 * zz(i+2) ) / ds
          end do
          ds1(npts-1) = -( gb1 * zz(npts) + gb2 * zz(npts-1) + 
     &                     gb3 * zz(npts-2) +
     &                     gb4 * zz(npts-3) + gb5 * zz(npts-4) ) / ds
          ds1(npts) = -( gc1 * zz(npts) + gc2 * zz(npts-1) + 
     &                   gc3 * zz(npts-2) +
     &                   gc4 * zz(npts-3) + gc5 * zz(npts-4) ) / ds
          do i = 1, npts
             ds1(i) = one / ds1(i)
          end do
c
c.... compute the second derivative of the mapping function
c
          dss1(1) = ( gc1 * ds1(1) + gc2 * ds1(2) + gc3 * ds1(3) +
     &                gc4 * ds1(4) + gc5 * ds1(5) ) / ds
          dss1(2) = ( gb1 * ds1(1) + gb2 * ds1(2) + gb3 * ds1(3) +
     &                gb4 * ds1(4) + gb5 * ds1(5) ) / ds
          do i = 3, npts-2
             dss1(i) = ( ga1 * ds1(i-2) + ga2 * ds1(i-1) + 
     &                   ga3 * ds1(i+1) + ga4 * ds1(i+2) ) / ds
          end do
          dss1(npts-1) = -( gb1 * ds1(npts) + gb2 * ds1(npts-1) + 
     &                      gb3 * ds1(npts-2) + gb4 * ds1(npts-3) + 
     &                      gb5 * ds1(npts-4) ) / ds
          dss1(npts) = -( gc1 * ds1(npts) + gc2 * ds1(npts-1) + 
     &                    gc3 * ds1(npts-2) +
     &                    gc4 * ds1(npts-3) + gc5 * ds1(npts-4) ) / ds
          do i = 1, npts
             dss1(i) = dss1(i) * ds1(i)
          end do

        end if
        
        call NR_SPLINE(zz,  s1,npts,1.0d11,1.0d11,  s2)
        call NR_SPLINE(zz, ds1,npts,1.0d11,1.0d11, ds2)
        call NR_SPLINE(zz,dss1,npts,1.0d11,1.0d11,dss2)
c
c.... arbitrary mapping function or uniform mesh
c
        if (.false.) then
          do i = 1, nx
            xi(i) = zero + real(i-1) * dxi
c           call NR_SPLINT(zz,s1,s2,npts, xi(i), ss(i), dss(i), d2ss(i))
            if (.true.) then
              call NR_SPLINT(zz,  s1,  s2,npts, xi(i),   ss(i), tmp, tmp)
              call NR_SPLINT(zz, ds1, ds2,npts, xi(i),  dss(i), tmp, tmp)
              call NR_SPLINT(zz,dss1,dss2,npts, xi(i), d2ss(i), tmp, tmp)
              ss(i)  = smin + ss(i) * smax
              dss(i) = dss(i) * smax
              d2ss(i) = d2ss(i) * smax
            else
              ss(i)   = smin + float(i-1) * ds
              dss(i)  = ds/dxi
              d2ss(i) = zero
            end if
            if (i.ne.1) then
              write(10,30) xi(i), ss(i), dss(i), d2ss(i),
     &                     (ss(i)-ss(i-1))/dxi
            else
              write(10,30) xi(i), ss(i), dss(i), d2ss(i), zero
            end if
          end do
        end if
c
c=============================================================================
c
c.... discretize in the body normal direction
c
        if (ny.gt.my) stop 'ny > my'
        deta = one / real(ny-1)
        dr   = (rmax-rmin) / real(ny-1)
        write(*,15) rmin, rmax
c
c.... exponential map in the wall normal direction
c
        if (yflag.eq.1) then
          write(*,"('Enter drmin ==> ',$)") 
          read(*,*) drmin
          sr = calcs( rmin, rmax, drmin, deta )
          write(*,*) 'Sr = ', sr
          c2 = log( drmin / deta )
          c1 = (log( sr*drmin / deta) - c2) / deta
          do j = 1, ny
            eta(j)  = real(j-1) * deta
            rr(j)   = rmin + one / c1 * ( exp(c1 * eta(j) + c2) - exp(c2) ) 
            drr(j)  = exp(c1 * eta(j) + c2)
            d2rr(j) = c1 * exp(c1 * eta(j) + c2)
            if (j.ne.1) then
              write(20,30) eta(j), rr(j), drr(j), d2rr(j),
     &                     (rr(j)-rr(j-1))/deta
            else
              write(20,30) eta(j), rr(j), drr(j), d2rr(j), zero
            end if
          end do
        end if
c
c.... Hyperbolic tangent stretching
c
        if (yflag.eq.2) then
        rd1 = 0.0005d0                  ! set for R=2400 MSE, r=50
        rd2 = 5.0d0
        dd  = 5.36966703089523d0        

c       rd1 = 0.0050d0                  ! set for R=2400 MSE, r=50
c       rd2 = 5.0d0
c       dd  = 3.898620676884022 

c       rd1 = 0.0100d0                  ! set for R=2400 MSE, r=50
c       rd2 = 5.0d0
c       dd  = 3.422429472016738d0       

        do j = 1, ny
          eta(j) = real(j-1) * deta
          rr(j) =  rmax*(pt5 + one/Tanh(dd*pt5)*Tanh(dd*(-pt5 +
     &             eta(j)))*pt5)/(Sqrt(rd2/rd1) +
     &             (one - Sqrt(rd2/rd1))*
     &             (pt5 + one/Tanh(dd*pt5)*Tanh(dd*(-pt5 +
     &             eta(j)))*pt5))
          drr(j) = dd*Sqrt(rd2/rd1)*rmax*Cosh(dd*pt5)*
     &             one/Cosh(dd*(-pt5 + eta(j)))*
     &             (Sinh(dd*(one - eta(j))) + Sinh(dd*eta(j)))/
     &             (Sqrt(rd2/rd1)*Sinh(dd*(one - eta(j))) + 
     &             Sinh(dd*eta(j)))**2
          d2rr(j) = dd**2*Sqrt(rd2/rd1)*rmax*
     &              (-Cosh(dd*(one - three*eta(j))) + 
     &              Sqrt(rd2/rd1)*Cosh(dd*(two - three*eta(j))) - 
     &              Cosh(dd*(one - eta(j))) + 
     &              two*Sqrt(rd2/rd1)*Cosh(dd*(one - eta(j))) - 
     &              two*Cosh(dd*eta(j)) + 
     &              Sqrt(rd2/rd1)*Cosh(dd*eta(j)))*
     &              (one/Cosh(dd*(-pt5 + eta(j))))**2*(Sinh(dd))/
     &              (two*(Sqrt(rd2/rd1)*(Sinh(dd*(one - eta(j)))) + 
     &              (Sinh(dd*eta(j))))**3)
          write(20,30) eta(j), rr(j), drr(j), d2rr(j)
        end do
        end if
c
c=============================================================================c
c.... Build the mesh
c=============================================================================c
c
c.... Compute the boundary and outward normal
c
        j = 1
        do i = 1, nx
          if (i.eq.1) then
             xl = xloc( xmin, ss(i) )
          else
             xl = xloc( xold, ss(i) - ss(i-1) )
          end if
          if (xl.lt.AR) then
             yl = yloc( xl )
             th = atan2( one/AR * rm/rn * (one - xl/AR)**(rm-one),
     &            (one - (one-xl/AR)**rm)**((rn-one)/rn) )
          else
             yl = one
             th = zero
          end if
          bn1(i) = -sin(th)
          bn2(i) =  cos(th)
          
          xb(i) = xl
          yb(i) = yl
          
          write(30,30) xl, yl, th, bn1(i), bn2(i)
          xold = xl
        end do
        write(*,*) 'Making the mesh'
c
c... Quick diagnostic
c
#ifdef MSE_DEBUG
        do i = 2, nx
          write(*,*) xb(i-1), xb(i), arc(xb(i-1),xb(i)), ss(i), arc(0.0,xb(i))
        end do
#endif
c
c.... Make the mesh
c
        do j = 1, ny
          do i = 1, nx
            x(i,j) = xb(i) + rr(j) * bn1(i)
            y(i,j) = yb(i) + rr(j) * bn2(i)
          end do
        end do
        write(*,*) 'Making the metrics'
c
c.... write out the body.dat file
c
        open(9,file='body.dat',status='unknown')
        j = 1
        do i = 1, nx
           write(9,"(1(i5,1x),5(1pe20.13,1x))") i, ss(i), xb(i), yb(i),
     &                                          bn1(i), bn2(i)
        end do
        close(9)
c
c.... write out a coord.dat file for the output plane
c
        open(9,file='coord.dat',status='unknown')
        i = nx
        do j = 1, ny
           write(9,"(2(i5,1x),2(1pe13.6,1x))") 0, 0, x(i,j), y(i,j)
        end do
        close(9)
c
c.... write out a top file
c
        open(9,file='top.dat',status='unknown')
        j = ny
        do i = 1, nx
           write(9,"(1(i5,1x),2(1pe13.6,1x))") i, x(i,j), y(i,j)
        end do
        close(9)
c
c.... Make the metrics
c
        do j = 1, ny
c         write(*,*) 'j = ',j
          do i = 1, nx
c==============================================================================
c.... first derivative metrics
c==============================================================================
            xl = xb(i) 
            yl = yb(i)

            if (xl .eq. zero) then
              dydx = 1.0d30             ! one / zero
            else if (xl .lt. AR) then
              dydx = one / AR * rm/rn * ( one - xl/AR )**(rm-one) *
     &               ( one - (one - xl/AR)**rm )**((one-rn)/rn)
            else
              dydx = zero
            end if
            dxbds  = one / sqrt(one + dydx**2)
  
            if (xl .eq. zero) then
              dxdy = zero
            else if (xl .lt. AR) then
              dxdy = AR * rn/rm * ( one - xl/AR )**(one-rm) *
     &               (one - (one - xl/AR)**rm)**((rn-one)/rn)
            else
              dxdy = 1.0d30             ! one / zero
            end if
            dybds = one / sqrt( dxdy**2 + one )
  
            if (xl .eq. zero) then
              dx = zero
              dy = one/AR * rm/rn * (one)**(rm-one)
              
              ddxdx = 1.0d30
              ddydx = one/AR * rm*(rm-one)/rn * (one)**(rm-two)*
     &                (-one/AR)

              if (rn.eq.two) then
                ddxdy = one
              else
                ddxdy = (rn-one)*yl**(rn-two)
              end if
              ddydy = (one-rm)/AR*(one-yl**rn)**(-one/rm)*yl**(rn-one)

              d2dxdx2 = 1.0d30
              d2dydx2 = rm*(two - three*rm + rm**2)*(one)**rm/
     &                  (rn*(AR)**3)

              if (rn.eq.two) then
                d2dxdy2 = zero
                d2dydy2 = (-one+rm)*(one-yl**rn)**(-one-one/rm)*
     &                    (rm-rm*rn-rm*yl**rn-rn*yl**rn+rm*rn*yl**rn)/(AR*rm)
              else
                d2dxdy2 = (-two + rn)*(-one + rn)*yl**(-three + rn)
                d2dydy2 = (-one+rm)*yl**(-two+rn)*(one-yl**rn)**(-one-one/rm)*
     &                    (rm-rm*rn-rm*yl**rn-rn*yl**rn+rm*rn*yl**rn)/(AR*rm)
              endif
            else if (xl .lt. AR) then
              dx = (one - (one - xl/AR)**rm)**((rn-one)/rn)
              dy = one/AR * rm/rn * (one - xl/AR)**(rm-one)

              ddxdx = (rn-one)/rn * (one - (one - xl/AR)**rm)**(-one/rn) *
     &                (rm/AR * (one - xl/AR)**(rm-one))
              ddydx = one/AR * rm*(rm-one)/rn * (one - xl/AR)**(rm-two)*
     &                (-one/AR)

              ddxdy = (rn-one)*yl**(rn-two)
              ddydy = (one-rm)/AR*(one-yl**rn)**(-one/rm)*yl**(rn-one)

              d2dxdx2 = rm*(-one + rn)*(one - xl/AR)**rm*
     &                  (one - (one - xl/AR)**rm)**(-one - one/rn)*
     &                  (rn - rm*rn - rm*(one - xl/AR)**rm - 
     &                  rn*(one - xl/AR)**rm + rm*rn*(one - xl/AR)**rm)/
     &                  (rn**2*(-AR + xl)**2)
              d2dydx2 = rm*(two - three*rm + rm**2)*(one - xl/AR)**rm/
     &                  (rn*(AR - xl)**3)

              if (rn.eq.two) then
                d2dxdy2 = zero
                d2dydy2 = (-one+rm)*(one-yl**rn)**(-one-one/rm)*
     &                    (rm-rm*rn-rm*yl**rn-rn*yl**rn+rm*rn*yl**rn)/(AR*rm)
              else
                d2dxdy2 = (-two + rn)*(-one + rn)*yl**(-three + rn)
                d2dydy2 = (-one+rm)*yl**(-two+rn)*(one-yl**rn)**(-one-one/rm)*
     &                    (rm-rm*rn-rm*yl**rn-rn*yl**rn+rm*rn*yl**rn)/(AR*rm)
              endif
            end if
            
            if ( xl .ge. AR ) then
              dbn1 = zero
              dbn2 = zero
            else if ( abs(bn1(i)) .gt. abs(bn2(i)) ) then
              dbn1 = ( -ddydy/(dx**2 + dy**2)**pt5 + pt5*dy*(two*dx*ddxdy + 
     &                  two*dy*ddydy)/(dx**2 + dy**2)**onept5 ) * dybds
              dbn2 = ( ddxdy/(dx**2 + dy**2)**pt5 - pt5*dx*(two*dx*ddxdy + 
     &                 two*dy*ddydy)/(dx**2 + dy**2)**onept5 ) * dybds
            else
              dbn1  = ( -ddydx/(dx**2 + dy**2)**pt5 + pt5*dy*(two*dx*ddxdx + 
     &                   two*dy*ddydx)/(dx**2 + dy**2)**onept5 ) * dxbds
              dbn2  = ( ddxdx/(dx**2 + dy**2)**pt5 - pt5*dx*(two*dx*ddxdx + 
     &                  two*dy*ddydx)/(dx**2 + dy**2)**onept5 ) * dxbds
            end if

            dsdxi  = dss(i)
            drdxi  = zero

            dsdeta = zero
            drdeta = drr(j)

            dxds  = dxbds + rr(j) * dbn1
            dxdr  = bn1(i)
            
            dyds  = dybds + rr(j) * dbn2
            dydr  = bn2(i)
            
            dxdxi  = dxds * dsdxi  + dxdr * drdxi
            dxdeta = dxds * dsdeta + dxdr * drdeta

            dydxi  = dyds * dsdxi  + dydr * drdxi
            dydeta = dyds * dsdeta + dydr * drdeta
           
            rjac(i,j) = one / ( dxdxi * dydeta - dxdeta * dydxi )
            
            rm1(i,j) =  dydeta * rjac(i,j)
            rm2(i,j) = -dxdeta * rjac(i,j)
            
            rn1(i,j) = -dydxi * rjac(i,j)
            rn2(i,j) =  dxdxi * rjac(i,j)
            
            if (j .eq. 1)  write(60,30) ss(i), bn1(i), dbn1, bn2(i), dbn2,
     &                     sqrt(dbn1**2 + dbn2**2)

            if (DIAG) then
            if (j .eq. 1)  write(40,30) ss(i), x(i,j), dxds, 
     &                     (x(i+1,j)-x(i-1,j))/(ss(i+1)-ss(i-1))
            if (j .eq. 1)  write(50,30) ss(i), y(i,j), dyds,
     &                     (y(i+1,j)-y(i-1,j))/(ss(i+1)-ss(i-1))
            if (j .eq. 1)  write(60,30) ss(i), bn1(i), dbn1, bn2(i), dbn2,
     &                     sqrt(dbn1**2 + dbn2**2)
            if (j .eq. ny) write(65,30) ss(i), rm1(i,j), rn1(i,j), rm2(i,j), 
     &                     rn2(i,j), rjac(i,j)
            if (j .eq. ny) write(66,30) ss(i), dxds, dsdxi, dxbds, dbn1
            if (j .eq. ny) write(67,30) ss(i), dydeta, dxdeta, dydxi, dxdxi
            end if

            if (j.eq.1) write(35,30) ss(i), sqrt(dbn1**2+dbn2**2)
c==============================================================================
c.... second derivative metrics
c==============================================================================
            if (xl .ge. AR) then
              d2xdy2  = 1.0d99 ! one / zero
              d2ydx2  = zero
              d2xbds2 = zero
              d2ybds2 = zero
            else if ( abs(bn1(i)) .gt. abs(bn2(i)) ) then
              if (rn.eq.two) then
                d2xdy2 = AR*rn*(one - yl**rn)**(one/rm)*
     &                   (-rm + rm*rn + rm*yl**rn - rn*yl**rn)/
     &                   (rm**2*(-one + yl**rn)**2)
              else
                d2xdy2 = AR*rn*yl**(-two + rn)*(one - yl**rn)**(one/rm)*
     &                   (-rm + rm*rn + rm*yl**rn - rn*yl**rn)/
     &                   (rm**2*(-one + yl**rn)**2)
              end if
              d2ybds2 = -(one + dxdy**2)**(-onept5) * dxdy *
     &                  d2xdy2 * dybds
              d2xbds2 = d2xdy2*(dybds)**2 + dxdy*d2ybds2
            else
              d2ydx2 = rm*(one - xl/AR)**rm*(one-(one-xl/AR)**rm)**(one/rn)*
     &                  (rn - rm*rn + rm*(one - xl/AR)**rm - 
     &                  rn*(one - xl/AR)**rm)/
     &                  (rn**2*(-AR + xl)**2*(-one + (one - xl/AR)**rm)**2)
              d2xbds2 = -(one + dydx**2)**(-onept5) * dydx *
     &                  d2ydx2 * dxbds
              d2ybds2 = d2ydx2*(dxbds)**2 + dydx*d2xbds2
            end if    

            if ( xl .ge. AR ) then
              d2bn1 = zero
              d2bn2 = zero
            else ! if ( abs(bn1(i)) .gt. abs(bn2(i)) ) then
              d2bn1 = ((ddxdy*(dy*ddxdy-dx*ddydy)+dx*(dy*d2dxdy2-dx*d2dydy2))/
     &                (dx**2+dy**2)**(onept5) -
     &                (three*dx*(dy*ddxdy-dx*ddydy)*(dx*ddxdy+dy*ddydy))/
     &                (dx**2+dy**2)**(twopt5))*(dybds)**2 +
     &                (-ddydy/(dx**2 + dy**2)**pt5 + pt5*dy*(two*dx*ddxdy + 
     &                  two*dy*ddydy)/(dx**2 + dy**2)**onept5) * d2ybds2
              d2bn2 = ((ddydy*(dy*ddxdy-dx*ddydy)+dy*(dy*d2dxdy2-dx*d2dydy2))/
     &                (dx**2+dy**2)**(onept5) -
     &                (three*dy*(dy*ddxdy-dx*ddydy)*(dx*ddxdy+dy*ddydy))/
     &                (dx**2+dy**2)**(twopt5))*(dybds)**2 +
     &                (ddxdy/(dx**2 + dy**2)**pt5 - pt5*dx*(two*dx*ddxdy + 
     &                 two*dy*ddydy)/(dx**2 + dy**2)**onept5) * d2ybds2
c           else
c             d2bn1 = (three*dx*(ddydx*dx - ddxdx*dy)*(ddxdx*dx + ddydx*dy)/
c     &         (dx**2 + dy**2)**(twopt5) + 
c     &         ddxdx*(-(ddydx*dx) + ddxdx*dy)/
c     &         (dx**2 + dy**2)**(onept5) + 
c     &         dx*(dy*d2dxdx2 - dx*d2dydx2)/
c     &         (dx**2 + dy**2)**(onept5))*(dxbds)**2 +
c     &         (dx*(-(ddydx*dx) + ddxdx*dy)/(dx**2 + dy**2)**(onept5)) *
c     &          d2xbds2
c
c                     ((ddxdx*(dy*ddxdx-dx*ddydx)+dx*(dy*d2dxdx2-dx*d2dydx2))/
c     &                (dx**2+dy**2)**(onept5) -
c     &                (three*dx*(dy*ddxdx-dx*ddydx)*(dx*ddxdx+dy*ddydx))/
c     &                (dx**2+dy**2)**(twopt5))*(dxbds)**2 +
c     &               (-ddydx/(dx**2 + dy**2)**pt5 + pt5*dy*(two*dx*ddxdx + 
c     &                  two*dy*ddydx)/(dx**2 + dy**2)**onept5) * d2xbds2
c
c             d2bn2 = (three*dy*(ddydx*dx - ddxdx*dy)*(ddxdx*dx + ddydx*dy)/
c     &         (dx**2 + dy**2)**(twopt5) + 
c     &         ddydx*(-(ddydx*dx) + ddxdx*dy)/
c     &         (dx**2 + dy**2)**(onept5) + 
c     &         dy*(dy*d2dxdx2 - dx*d2dydx2)/
c     &         (dx**2 + dy**2)**(onept5))*(dxbds)**2 +
c     &         (dy*(-(ddydx*dx) + ddxdx*dy)/(dx**2 + dy**2)**(onept5)) *
c     &         d2xbds2
c             
c             ((ddydx*(dy*ddxdx-dx*ddydx)+dy*(dy*d2dxdx2-dx*d2dydx2))/
c     &                (dx**2+dy**2)**(onept5) -
c     &                (three*dy*(dy*ddxdx-dx*ddydx)*(dx*ddxdx+dy*ddydx))/
c     &                (dx**2+dy**2)**(twopt5))*(dxbds)**2 +
c     &               (ddxdx/(dx**2 + dy**2)**pt5 - pt5*dx*(two*dx*ddxdx + 
c     &                 two*dy*ddydx)/(dx**2 + dy**2)**onept5) * d2xbds2
            end if
            
            d2sdxi2    = d2ss(i)
            d2sdxideta = zero
            d2sdeta2   = zero

            d2rdxi2    = zero
            d2rdxideta = zero
            d2rdeta2   = d2rr(j)

            d2xds2  = d2xbds2 + rr(j) * d2bn1
            d2xdsdr = dbn1
            d2xdr2  = zero
            
            d2yds2  = d2ybds2 + rr(j) * d2bn2
            d2ydsdr = dbn2
            d2ydr2  = zero
            
            d2xdxi2 = d2xds2 * dsdxi**2 + dxds * d2sdxi2 + 
     &                d2xdsdr * drdxi * dsdxi + d2xdsdr * drdxi * dsdxi +
     &                d2xdr2 * drdxi**2 + dxdr * d2rdxi2
     
            d2xdxideta = d2xds2 * dsdeta * dsdxi + dxds * d2sdxideta + 
     &                   d2xdsdr * drdeta * dsdxi + d2xdsdr * drdxi * dsdeta +
     &                   d2xdr2 * drdeta * drdxi + dxdr * d2rdxideta

            d2xdeta2 = d2xds2 * dsdeta**2 + dxds * d2sdeta2 + 
     &                 d2xdsdr * drdeta * dsdeta + d2xdsdr * drdeta * dsdeta +
     &                 d2xdr2 * drdeta**2 + dxdr * d2rdeta2

            d2ydxi2 = d2yds2 * dsdxi**2 + dyds * d2sdxi2 + 
     &                d2ydsdr * drdxi * dsdxi + d2ydsdr * drdxi * dsdxi +
     &                d2ydr2 * drdxi**2 + dydr * d2rdxi2
            
            d2ydxideta = d2yds2 * dsdeta * dsdxi + dyds * d2sdxideta + 
     &                   d2ydsdr * drdeta * dsdxi + d2ydsdr * drdxi * dsdeta +
     &                   d2ydr2 * drdeta * drdxi + dydr * d2rdxideta

            d2ydeta2 = d2yds2 * dsdeta**2 + dyds * d2sdeta2 + 
     &                 d2ydsdr * drdeta * dsdeta + d2ydsdr * drdeta * dsdeta +
     &                 d2ydr2 * drdeta**2 + dydr * d2rdeta2

            drjacdxi = -(d2xdxi2*dydeta + dxdxi*d2ydxideta - d2xdxideta*dydxi - 
     &                   dxdeta*d2ydxi2) * rjac(i,j)**2
     
            drjacdeta = -(d2xdxideta*dydeta + dxdxi*d2ydeta2 - d2xdeta2*dydxi -
     &                    dxdeta*d2ydxideta) * rjac(i,j)**2
     
            rm11(i,j) = drjacdxi * dydeta * rm1(i,j) +
     &                  d2ydxideta * rjac(i,j) * rm1(i,j) +
     &                  drjacdeta * dydeta * rn1(i,j) +
     &                  d2ydeta2 * rjac(i,j) * rn1(i,j)
     
            rm12(i,j) = drjacdxi * dydeta * rm2(i,j) +
     &                  d2ydxideta * rjac(i,j) * rm2(i,j) +
     &                  drjacdeta * dydeta * rn2(i,j) +
     &                  d2ydeta2 * rjac(i,j) * rn2(i,j)

            rm21(i,j) = -( drjacdxi * dxdeta * rm1(i,j) +
     &                     d2xdxideta * rjac(i,j) * rm1(i,j) +
     &                     drjacdeta * dxdeta * rn1(i,j) +
     &                     d2xdeta2 * rjac(i,j) * rn1(i,j) )

            rm22(i,j) = -( drjacdxi * dxdeta * rm2(i,j) +
     &                     d2xdxideta * rjac(i,j) * rm2(i,j) +
     &                     drjacdeta * dxdeta * rn2(i,j) +
     &                     d2xdeta2 * rjac(i,j) * rn2(i,j) )

            rn11(i,j) = -( drjacdxi * dydxi * rm1(i,j) +
     &                     d2ydxi2 * rjac(i,j) * rm1(i,j) +
     &                     drjacdeta * dydxi * rn1(i,j) +
     &                     d2ydxideta * rjac(i,j) * rn1(i,j) )
     
            rn12(i,j) = -( drjacdxi * dydxi * rm2(i,j) +
     &                     d2ydxi2 * rjac(i,j) * rm2(i,j) +
     &                     drjacdeta * dydxi * rn2(i,j) +
     &                     d2ydxideta * rjac(i,j) * rn2(i,j) )
     
            rn21(i,j) = drjacdxi * dxdxi * rm1(i,j) +
     &                  d2xdxi2 * rjac(i,j) * rm1(i,j) +
     &                  drjacdeta * dxdxi * rn1(i,j) +
     &                  d2xdxideta * rjac(i,j) * rn1(i,j)
     
            rn22(i,j) = drjacdxi * dxdxi * rm2(i,j) +
     &                  d2xdxi2 * rjac(i,j) * rm2(i,j) +
     &                  drjacdeta * dxdxi * rn2(i,j) +
     &                  d2xdxideta * rjac(i,j) * rn2(i,j)

            if (DIAG) then
            if (j .eq. 1) write(70,30) ss(i), dxbds, d2xbds2, dybds, d2ybds2
            if (j .eq. 1) write(71,30) ss(i), bn1(i), dbn1, d2bn1
            if (j .eq. 1) write(72,30) ss(i), bn2(i), dbn2, d2bn2
            end if

            if (j .eq. ny) write(74,30) x(i,j), rm1(i,j), rn1(i,j), rm2(i,j), 
     &                                  rn2(i,j), rjac(i,j)
            if (j .eq. ny) write(75,30) x(i,j), rm11(i,j), rm12(i,j), 
     &                                  rm21(i,j), rm22(i,j)
            if (j .eq. ny) write(76,30) x(i,j), rn11(i,j), rn12(i,j), 
     &                                  rn21(i,j), rn22(i,j)
            if (j .eq. ny) write(77,30) ss(i), d2xds2, d2xbds2, d2bn1
            if (j .eq. ny) write(78,30) ss(i), d2yds2, d2ybds2, d2bn2

            g11x(i,j) = d2xdxi2
            g11y(i,j) = d2ydxi2
            g12x(i,j) = d2xdxideta
            g12y(i,j) = d2ydxideta
            g22x(i,j) = d2xdeta2
            g22y(i,j) = d2ydeta2
          end do
        end do
c
c.... compute the error in the cross derivatives
c
        errn = zero
        errm = zero
        do i = 1, nx
          do j = 1, ny
            errm = errm + ( rm12(i,j) - rm21(i,j) )**2
            errn = errn + ( rn12(i,j) - rn21(i,j) )**2
          end do
        end do
        write(*,*) 'Error in m ',sqrt(errm)/real(nx*ny)
        write(*,*) 'Error in n ',sqrt(errn)/real(nx*ny)
c
c..... diagnostic
c
        j = 1
        do i = 1, nx
          write(80,30) ss(i),
     &               (-rm1(i,j)*rm2(i,j)*rm11(i,j) + rm1(i,j)**2*rm12(i,j) -
     &               rm2(i,j)**2*rm12(i,j) + rm1(i,j)*rm2(i,j)*rm22(i,j) ) *
     &               sqrt(rn1(i,j)**2 + rn2(i,j)**2) /
     &               sqrt(rm1(i,j)**2 + rm2(i,j)**2) /
     &               (rm1(i,j)*rn2(i,j)-rm2(i,j)*rn1(i,j)),
     &               (rn1(i,j)*rn2(i,j)*rn11(i,j) - rn1(i,j)**2*rn12(i,j) +
     &               rn2(i,j)**2*rn12(i,j) - rn1(i,j)*rn2(i,j)*rn22(i,j) ) *
     &               sqrt(rm1(i,j)**2 + rm2(i,j)**2) /
     &               sqrt(rn1(i,j)**2 + rn2(i,j)**2) /
     &               (rm1(i,j)*rn2(i,j)-rm2(i,j)*rn1(i,j))
        end do
c
c.... write out some Plot3d files
c    
        nz = 1
        write(*,*) 'Writing plot3d files for ',nx, ny, nz
c
c.... write out binary if on an IRIS
c
        if (IRIS) then
          call wgrid(x, y, nx, ny, nz, mx, my, 'grid.bin'//char(0))
          call wdata(rm1, rn1, rm2, rn2, rjac, nx, ny, nz, mx, my, 
     &               'm1.bin'//char(0))
          call wdata(rm11, rm12, rm21, rm22, rjac, nx, ny, nz, mx, my, 
     &               'm11.bin'//char(0))
          call wdata(rn11, rn12, rn21, rn22, rjac, nx, ny, nz, mx, my, 
     &               'n11.bin'//char(0))
          call wdata(g11x, g11y, g12x, g12y, g22x, nx, ny, nz, mx, my, 
     &               'grad.bin'//char(0))
        else
c
c.... write out unformatted if on a CRAY
c
          open(unit=10,file='grid.dat',form='unformatted')
          write(10) nx, ny, nz
          write(10) (((real(x(i,j)), i=1,nx), j=1,ny), k = 1, nz),
     &              (((real(y(i,j)), i=1,nx), j=1,ny), k = 1, nz),
     &              (((real(zero),   i=1,nx), j=1,ny), k = 1, nz)
          close(10)
        
          open(unit=10,file='m1.dat',form='unformatted')
          write(10) nx, ny, nz
          write(10) 0.0, 0.0, 0.0, 0.0
          write(10) (((real( rm1(i,j)), i=1,nx), j=1,ny), k = 1, nz),
     &              (((real( rn1(i,j)), i=1,nx), j=1,ny), k = 1, nz),
     &              (((real( rm2(i,j)), i=1,nx), j=1,ny), k = 1, nz),
     &              (((real( rn2(i,j)), i=1,nx), j=1,ny), k = 1, nz),
     &              (((real(rjac(i,j)), i=1,nx), j=1,ny), k = 1, nz)
          close(10)

          open(unit=10,file='m11.dat',form='unformatted')
          write(10) nx, ny, nz
          write(10) 0.0, 0.0, 0.0, 0.0
          write(10) (((real(rm11(i,j)), i=1,nx), j=1,ny), k = 1, nz),
     &              (((real(rm12(i,j)), i=1,nx), j=1,ny), k = 1, nz),
     &              (((real(rm21(i,j)), i=1,nx), j=1,ny), k = 1, nz),
     &              (((real(rm22(i,j)), i=1,nx), j=1,ny), k = 1, nz),
     &              (((real(rm1(i,j)*rn1(i,j)+rm2(i,j)*rn2(i,j)), 
     &              i=1,nx), j=1,ny), k = 1, nz)
          close(10)

          open(unit=10,file='n11.dat',form='unformatted')
          write(10) nx, ny, nz
          write(10) 0.0, 0.0, 0.0, 0.0
          write(10) (((real(rn11(i,j)), i=1,nx), j=1,ny), k = 1, nz),
     &              (((real(rn12(i,j)), i=1,nx), j=1,ny), k = 1, nz),
     &              (((real(rn21(i,j)), i=1,nx), j=1,ny), k = 1, nz),
     &              (((real(rn22(i,j)), i=1,nx), j=1,ny), k = 1, nz),
     &              (((real(rjac(i,j)), i=1,nx), j=1,ny), k = 1, nz)
          close(10)

          if (metric_ji) then
c
c.... write out the metric file (Note the order of i and j are reversed)
c
            open (unit=10, file='metric.dat', form='unformatted',
     &            status='unknown')
            write(10) (( real( rm1(i,j)), j=1,ny), i=1,nx),
     &                (( real( rm2(i,j)), j=1,ny), i=1,nx),
     &                (( real( rn1(i,j)), j=1,ny), i=1,nx),
     &                (( real( rn2(i,j)), j=1,ny), i=1,nx),
     &                (( real(rm11(i,j)), j=1,ny), i=1,nx),
     &                (( real(rm12(i,j)), j=1,ny), i=1,nx),
     &                (( real(rm22(i,j)), j=1,ny), i=1,nx),
     &                (( real(rn11(i,j)), j=1,ny), i=1,nx),
     &                (( real(rn12(i,j)), j=1,ny), i=1,nx),
     &                (( real(rn22(i,j)), j=1,ny), i=1,nx)
            close(10)
          else
c
c.... write out the metric file in IJ ordering
c
            open (unit=10, file='metric.dat', form='unformatted',
     &            status='unknown')
            write(10) (( real( rm1(i,j)), i=1,nx), j=1,ny),
     &                (( real( rm2(i,j)), i=1,nx), j=1,ny),
     &                (( real( rn1(i,j)), i=1,nx), j=1,ny),
     &                (( real( rn2(i,j)), i=1,nx), j=1,ny),
     &                (( real(rm11(i,j)), i=1,nx), j=1,ny),
     &                (( real(rm12(i,j)), i=1,nx), j=1,ny),
     &                (( real(rm22(i,j)), i=1,nx), j=1,ny),
     &                (( real(rn11(i,j)), i=1,nx), j=1,ny),
     &                (( real(rn12(i,j)), i=1,nx), j=1,ny),
     &                (( real(rn22(i,j)), i=1,nx), j=1,ny)
            close(10)
          end if
        end if

        stop      
 10     format('xmin = ',1pe13.6,' xmax = ',1pe13.6,' S = ',1pe13.6)
 15     format('rmin = ',1pe13.6,' rmax = ',1pe13.6)
 20     format(i5,6(1x,1pe13.6))
 30     format(7(1x,1pe20.13))
        end

c=============================================================================c
        function calcs( ximin1, ximax1, dxmin1, ds1)
c
c  Calculate a stretching ratio to fill a domain using exponential
c  stretching.
c
c=============================================================================c
c       implicit double precision (a-h,o-z)
        parameter (zero=0.0d0, pt5=0.5d0, one=1.0d0, two=2.0d0)
        
        common /stuffs/ ximin, ximax, dxmin, ds

        external funcs, rtflsp
c=============================================================================c
        ximin = ximin1
        ximax = ximax1
        dxmin = dxmin1
        ds = ds1
        
c       calcs = rtflsp(funcs, 1.01e0, 1.02e0, 1.0e-11)
        calcs = zbrent(funcs, 1.0001e0, 1.5e0, 1.0e-11)

        return
        end
c=============================================================================c
        function funcs(sx)
c=============================================================================c
c       implicit double precision (a-h,o-z)
        parameter (zero=0.0d0, pt5=0.5d0, one=1.0d0, two=2.0d0)
        
        common /stuffs/ ximin, ximax, dxmin, ds
c=============================================================================c
        c2 = log( dxmin / ds )
        c1 = (log( sx*dxmin / ds) - c2) / ds
             
        funcs = ximax - ximin - one / c1 * ( exp(c1 + c2) - exp(c2) )
        
        return
        end 
c=============================================================================c
        function yloc(x)
c=============================================================================c
c       implicit double precision (a-h,o-z)
        parameter (zero=0.0d0, pt5=0.5d0, one=1.0d0, two=2.0d0)

        common /stuff/ AR, rm, rn, xmin, xmax
c=============================================================================c
        yloc = ( one - ( one - x/AR )**rm )**(one/rn)

        return
        end
c=============================================================================c
        function xloc(x,ds)
c=============================================================================c
c       implicit double precision (a-h,o-z)
        parameter (zero=0.0d0, pt5=0.5d0, one=1.0d0, two=2.0d0)

        external rtsec, rtflsp, func, zbrent

        common /stuff/ AR, rm, rn, xmin, xmax
        common /distance/ darc, x1
c=============================================================================c
        darc = ds
        x1   = x

c       xloc = rtsec(func,x1,xmax,1.0e-11)
        xloc = rtflsp(func,x1,x1+two*ds,1.0e-11)
c       xloc = zbrent(func,x1,x1+two*ds,1.0e-11)

        return
        end 
c=============================================================================c
        function func(x)
c=============================================================================c
c       implicit double precision (a-h,o-z)
        parameter (zero=0.0d0, pt5=0.5d0, one=1.0d0, two=2.0d0)

        external arc

        common /distance/ darc, x1
c=============================================================================c
        func = darc - arc(x1,x)

        return
        end 
c=============================================================================c
        function arc(x1,x2)
c=============================================================================c
c       implicit double precision (a-h,o-z)
        parameter (zero=0.0d0, pt5=0.5d0, one=1.0d0, two=2.0d0)

        parameter EPS = 1.0e-10
        parameter HMIN = 1.0e-8

        external derivs1, derivs2, RKQCR

        common /stuff/ AR, rm, rn, xmin, xmax
c=============================================================================c
        s = zero
        xtmp = AR * pt5

        if (x1 .eq. x2 ) then
           arc = zero
           return
        end if

!       if (x1 .gt. x2) then
!          write(*,*) x1, x2
!          stop 'Error in arc:  x1 > x2'
!       end if
!       if (x1 .lt. xmin) stop 'Error in arc:  x1 < xmin'
!       if (x2 .gt. xmax) stop 'Error in arc:  x2 > xmax'

        if (x1 .le. xtmp) then
           y1 = ( one - ( one - x1/AR )**rm )**(one/rn)
           if (x2 .le. xtmp) then
             y2 = ( one - ( one - x2/AR )**rm )**(one/rn)
           else
             y2 = ( one - ( one - xtmp/AR )**rm )**(one/rn)
           end if
           call ODEINTR(s,1,y1,y2,EPS,(y2-y1)*pt5,
     &                  HMIN,nok,nbad,derivs2,RKQCR)
c          write(*,*) "1: NOK = ", nok, ', NBAD = ', nbad
           if (x2 .gt. xtmp) then
             call ODEINTR(s,1,xtmp,x2,EPS,(x2-xtmp)*pt5,
     &                    HMIN,nok,nbad,derivs1,RKQCR)
c            write(*,*) "2: NOK = ", nok, ', NBAD = ', nbad
           end if
        else
           call ODEINTR(s,1,x1,x2,EPS,(x2-x1)*pt5,
     &                  HMIN,nok,nbad,derivs1,RKQCR)
c          write(*,*) "3: NOK = ", nok, ', NBAD = ', nbad
        end if

        arc = s

        return
        end
c=============================================================================c
        subroutine derivs1(n,s,x,ds)
c=============================================================================c
c       implicit double precision (a-h,o-z)
        parameter (zero=0.0d0, pt5=0.5d0, one=1.0d0, two=2.0d0)

        common /stuff/ AR, rm, rn, xmin, xmax
c=============================================================================c
        if (x.lt.AR) then
          dydx = one / AR * rm/rn * ( one - x/AR )**(rm-one) *
     &           ( one - (one - x/AR)**rm )**((one-rn)/rn)
        else
          dydx = zero
        end if

        ds = sqrt( one + dydx**2 )

        return
        end
c=============================================================================c
        subroutine derivs2(n,s,y,ds)
c=============================================================================c
c       implicit double precision (a-h,o-z)
        parameter (zero=0.0d0, pt5=0.5d0, one=1.0d0, two=2.0d0)

        common /stuff/ AR, rm, rn, xmin, xmax
c=============================================================================c
        x = AR * ( one - (one - y**rn)**(one/rm) )

c       if (x .gt. pt5*AR) stop 'Error in derivs2'

        dxdy = AR * rn/rm * ( one - x/AR )**(one-rm) *
     &         (one - (one - x/AR)**rm)**((rn-one)/rn)

        ds = sqrt( dxdy**2 + one )

        return
        end
