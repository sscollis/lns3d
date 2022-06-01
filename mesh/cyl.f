c=============================================================================c
        program cyl
c=============================================================================c
c
c  Purpose:  make a body-fitted circular-cylinder mesh
c
c  Inputs:   This uses Maheshs mapping in x and y
c
c  Author:   S. Scott Collis
c
c  Revised:  8-21-96
c            3-36-01  Switched Metric file to IJ order
c
c=============================================================================c
c       implicit double precision (a-h,o-z)

        logical IRIS, DIAG
        parameter (IRIS = .false.)
        parameter (DIAG = .false.)

        parameter (zero=0.0d0,    pt5=0.5d0,  one=1.0d0, 
     &             onept5=1.5d0,  two=2.0d0,  twopt5=2.5d0, 
     &             three=3.0d0)

        parameter (mx=1024,my=255)
        
        dimension xb(mx), yb(mx)
        dimension bn1(mx), bn2(mx)
        
        dimension xi(mx), eta(my)

        dimension ss(mx), rr(my), dss(mx), drr(my), d2ss(mx), d2rr(my)

        dimension x(mx,my),    y(mx,my)
        dimension rm1(mx,my),  rm2(mx,my)
        dimension rn1(mx,my),  rn2(mx,my),  rjac(mx,my)
        dimension rm11(mx,my), rm12(mx,my), rm21(mx,my), rm22(mx,my)
        dimension rn11(mx,my), rn12(mx,my), rn21(mx,my), rn22(mx,my)
        dimension g11x(mx,my), g11y(mx,my), g12x(mx,my), g12y(mx,my)
        dimension g22x(mx,my), g22y(mx,my)
        
        dimension d(10), e(10), f1(10), f2(10), g0(10)
        
        common /stuff/ radius

        external xloc, yloc

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
c=============================================================================c
        xmin = -one
        xmax =  one

        write(*,*) 'Make a body fitted Mesh on a Circular Cylinder'
        write(*,*) '*** metric.dat file is in IJ format ***'
                                              
c....   Radius = radius of cylinder
c....   xmin, xmax = are the start and end stations, e.g. (-Radius,Radius)
c....   x0, y0 = a shift in the origin
c
        write(*,"('Enter Radius, Xmin, Xmax, x0, y0 ==> ',$)")
        read(*,*) radius, xmin, xmax, x0, y0
c
c.... discretize in the body tangent direction
c
        smin = zero
        smax = arc(xmin,xmax)
        savg = (smin + smax) * pt5
        write(*,10) xmin, xmax, smax

        write(*,"('Enter Nx ==> ',$)")
        read(*,*) nx
        if (nx.gt.mx) stop 'nx > mx'

        dxi  = one / real(nx-1)
        ds   = (smax-smin) / real(nx-1)

        write(*,"('ds uniform = ',1pe20.13)") (smax-smin)/real(nx-1)
c
c.... Uniform arc-length
c
        if (.true.) then
          do i = 1, nx
            xi(i)   = real(i-1) * dxi
            ss(i)   = smin + float(i-1) * ds
            dss(i)  = ds / dxi
            d2ss(i) = zero
          end do
        end if
c
c.... Uniform arc-length + mapping
c
        if (.false.) then
        write(*,"('Enter dsmin, b, sc ==> ',$)") 
        read(*,*) dsmin, b, sc

        cm = ( two * b * tanh(b*sc) + real(nx-1)*dsmin/smax * 
     &         log( cosh(b*(one-sc)) / cosh(b*(one+sc)) ) ) /
     &       ( one - real(nx-1)*dsmin/smax )

        do i = 1, nx
           xi(i)   = real(i-1) * dxi
           ss(i)   = smax*( cm * xi(i) + log( cosh(b*(xi(i)-sc)) /
     &               cosh(b*(xi(i)+sc)) ) ) / 
     &               (cm + log( cosh(b*(one-sc)) / cosh(b*(one+sc)) ) )
           dss(i)  = smax*(cm - b*(1.0/Cosh(b*(xi(i) - sc)))*
     &               (1.0/Cosh(b*(xi(i) + sc)))*Sinh(2.0*b*sc)) /
     &               (cm + log( cosh(b*(one-sc)) / cosh(b*(one+sc)) ) )
           d2ss(i) = smax*(-b**2*(tanh(b*(xi(i)-sc)))**2 + 
     &               b**2*(tanh(b*(xi(i)+sc)))**2)/
     &               (cm + log( cosh(b*(one-sc)) / cosh(b*(one+sc)) ) )
        end do
        do i = 1, nx
           if (i.eq.1) then
             write(10,30) xi(i), ss(i), dss(i), d2ss(i), 
     &                    (dss(i+1)-dss(i))/dxi
           else if (i.eq.nx) then
             write(10,30) xi(i), ss(i), dss(i), d2ss(i), 
     &                    (dss(i)-dss(i-1))/dxi
           else
             write(10,30) xi(i), ss(i), dss(i), d2ss(i), 
     &                    (dss(i+1)-dss(i-1))/dxi*pt5
           end if
        end do
        end if
c
c.... discretize in the body normal direction
c
        write(*,*) 'Discretize in the wall normal direction'
        write(*,"('Enter Ny ==> ',$)")
        read(*,*) ny
        if (ny.gt.my) stop 'ny > my'
        rmin = zero
        write(*,"('Enter Rmax ==> ',$)")
        read(*,*) rmax

        dr = (rmax-rmin)/real(ny-1)
        deta = one / real(ny-1)
        write(*,"('dr uniform = ',1pe20.13)") (rmax-rmin)/real(ny-1)
c
c.... Uniform spacing
c
        if (.true.) then
          do j = 1, ny
            eta(j)  = real(j-1) * deta
            rr(j)   = rmin + real(j-1) * dr
            drr(j)  = dr / deta
            d2rr(j) = zero
            write(20,30) eta(j), rr(j), drr(j), d2rr(j)
          end do
        end if
c
c.... exponential map in the wall normal direction
c
        if (.false.) then
        deta2 = one / real(128-1)
        dymin = 5.0d-5
        sy = 1.088d0
        c2 = log( dymin / deta2 )
        c1 = (log( sy*dymin / deta2) - c2) / deta2
        do j = 1, ny
          eta(j)  = real(j-1) * deta
          rr(j)   = one / c1 * ( exp(c1 * eta(j) + c2) - exp(c2) ) 
          drr(j)  = exp(c1 * eta(j) + c2)
          d2rr(j) = c1 * exp(c1 * eta(j) + c2)
          write(20,30) eta(j), rr(j), drr(j), d2rr(j)
        end do
        end if
c
c.... Hyperbolic tangent stretching
c
        if (.false.) then       
        rd1 = 0.00028d0                 ! set for R=1000 PCYL, r=500
        rd2 = 5.5d0
        dd  = 5.665549493763092d0
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
          write(20,30) eta(j), rr(j), drr(j), d2rr(j), d2rr(j)/drr(j)
        end do
        end if
c
c.... Maheshs mapping, personal communication (6-5-96)
c
        if (.false.) then
          write(*,"('Enter drmin, b, rc ==> ',$)") 
          read(*,*) drmin, b, rc
          cm = ( two * b * tanh(b*rc) + (ny-1)*drmin/rmax * 
     &           log( cosh(b*(one-rc)) / cosh(b*(one+rc)) ) ) /
     &         ( one - (ny-1)*drmin/rmax )
          do j = 1, ny
            eta(j)  = real(j-1) * deta
            rr(j)   = rmax*( cm * eta(j) + log( cosh(b*(eta(j)-rc)) /
     &                cosh(b*(eta(j)+rc)) ) ) / 
     &                (cm + log( cosh(b*(one-rc)) / cosh(b*(one+rc)) ) )
            drr(j)  = rmax*(cm + b*tanh(b*(eta(j)-rc)) - 
     &                b*tanh(b*(eta(j)+rc))) /
     &                (cm + log( cosh(b*(one-rc)) / cosh(b*(one+rc)) ) )
            d2rr(j) = rmax*(-b**2*(tanh(b*(eta(j)-rc)))**2 + 
     &                b**2*(tanh(b*(eta(j)+rc)))**2)/
     &                (cm + log( cosh(b*(one-rc)) / cosh(b*(one+rc)) ) )
            write(20,30) eta(j), rr(j), drr(j), d2rr(j)
          end do
        end if
c
c.... Compute the boundary and outward normal
c
        j = 1
        do i = 1, nx
          if (i.eq.1) then
             xl = xmin
          else
             xl = xloc( xold, ss(i) - ss(i-1) )
          end if
          yl = yloc( xl )
          th = atan2( -xl, sqrt( radius**2 - xl**2 ) ) 
          bn1(i) = -sin(th)
          bn2(i) =  cos(th)
          
          xb(i) = xl
          yb(i) = yl
          
          write(30,30) xl, yl, th, bn1(i), bn2(i)
          xold = xl
        end do
        write(*,*) 'Making the mesh'
c
c.... Make the mesh
c
        do j = 1, ny
          do i = 1, nx
            x(i,j) = xb(i) - x0 + rr(j) * bn1(i) 
            y(i,j) = yb(i) - y0 + rr(j) * bn2(i)
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
     &                                         bn1(i), bn2(i)
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
c.... write out the nose element spacing
c
c       dr = x(1,1) - x(1,2)
c       ds = atan2( sqrt(1.0-x(2,1)**2), x(2,1) ) -
c     &      atan2( sqrt(1.0-x(1,1)**2), x(1,1) )
c
c       write(*,11) ds, dr, (x(1,2)-x(1,3))/dr
c 11    format('ds_min = ',1pe20.13,' dr_min = ',1pe20.13,
c     &         ' dr(2)/dr(1)     = ',1pe20.13)
c
c.... write out the maximum spacing
c
c       dr = x(1,ny-1) - x(1,ny)
c       yl = y(nx,1)
c       ds = yl * 0.5 * sqrt( 1.0 + yl**2) + 
c     &       0.5 * log( abs( yl + sqrt(1.0 + yl**2)))
c        yl = y(nx-1,1)
c       ds = ds - (yl * 0.5 * sqrt( 1.0 + yl**2) + 
c     &       0.5 * log( abs( yl + sqrt(1.0 + yl**2))))
c       write(*,21) ds, dr, (x(1,ny-2)-x(1,ny-1))/dr
c 21    format('ds_max = ',1pe20.13,' dr_max = ',1pe20.13,
c     &         ' dr(ny)/dr(ny-1) = ',1pe20.13)
c
c.... Make the metrics
c
        do j = 1, ny
          do i = 1, nx
c==============================================================================
c.... first derivative metrics
c==============================================================================
            xl = xb(i) 
            yl = yb(i)
            xtmp = one

            if (yl .eq. zero) then
              dydx = -1.0d30
            else
              dydx = (-xl)/yl
            end if
            dxbds  = one / sqrt(one + dydx**2)

            if (xl .eq. zero) then
              dxdy = -1.0d30
            else
              dxdy = yl/(-xl)
            end if
            if (xl .lt. zero) then
              dybds = one / sqrt( dxdy**2 + one )
            else if (xl .eq. zero) then
              dybds = zero
            else
              dybds = -one / sqrt( dxdy**2 + one )
            end if

            if (yl .eq. zero) then
              dx = zero
              dy = -xl
              
              ddxdx = -1.0d30
              ddydx = -one

              ddxdy = one
              ddydy = zero

              d2dxdx2 = -1.0d30
              d2dydx2 = zero

              d2dxdy2 = zero
              d2dydy2 = (xl**2 + yl**2)/xl**3
            else if (xl .eq. zero) then
              dx = yl
              dy = zero

              ddxdx = zero
              ddydx = -one

              ddxdy = one
              ddydy = 1.0d30

              d2dxdx2 = -(yl**2 + xl**2)/yl**3
              d2dydx2 = zero

              d2dxdy2 = zero
              d2dydy2 = 1.0d30  
            else
              dx = yl
              dy = -xl

              ddxdx = -xl / yl
              ddydx = -one

              ddxdy = one
              ddydy = yl / xl

              d2dxdx2 = -(yl**2 + xl**2)/yl**3
              d2dydx2 = zero

              d2dxdy2 = zero
              d2dydy2 = (xl**2 + yl**2)/xl**3
            end if
            
            if ( abs(bn1(i)) .gt. abs(bn2(i)) ) then
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
     &                       (x(i+1,j)-x(i-1,j))/(ss(i+1)-ss(i-1))
              if (j .eq. 1)  write(50,30) ss(i), y(i,j), dyds,
     &                       (y(i+1,j)-y(i-1,j))/(ss(i+1)-ss(i-1))
              if (j .eq. 1)  write(60,30) ss(i), bn1(i), dbn1, bn2(i), dbn2,
     &                       sqrt(dbn1**2 + dbn2**2)
              if (j .eq. ny) write(65,30) ss(i), rm1(i,j), rn1(i,j), rm2(i,j), 
     &                       rn2(i,j), rjac(i,j)
              if (j .eq. ny) write(66,30) ss(i), dxds, dsdxi, dxbds, dbn1
              if (j .eq. ny) write(67,30) ss(i), dydeta, dxdeta, dydxi, dxdxi
            end if

            if (j.eq.1) write(35,30) ss(i), sqrt(dbn1**2+dbn2**2)
c==============================================================================
c.... second derivative metrics
c==============================================================================
            if ( abs(bn1(i)) .gt. abs(bn2(i)) ) then
              d2xdy2  = -(xl**2+yl**2)/xl**3
              if (xl.lt.zero) then
                d2ybds2 = -(one + dxdy**2)**(-onept5) * dxdy * d2xdy2 * dybds
              else
                d2ybds2 = (one + dxdy**2)**(-onept5) * dxdy * d2xdy2 * dybds
              end if
              d2xbds2 = d2xdy2*(dybds)**2 + dxdy*d2ybds2
            else
              d2ydx2  = -(xl**2+yl**2)/yl**3
              d2xbds2 = -(one + dydx**2)**(-onept5) * dydx * d2ydx2 * dxbds
              d2ybds2 = d2ydx2*(dxbds)**2 + dydx*d2xbds2
            end if    

            if (j .eq. 1)  write(61,30) ss(i), dxbds, d2xbds2, dybds, d2ybds2

            if ( abs(bn1(i)) .gt. abs(bn2(i)) ) then
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
            else
              d2bn1 = ((ddxdx*(dy*ddxdx-dx*ddydx)+dx*(dy*d2dxdx2-dx*d2dydx2))/
     &                (dx**2+dy**2)**(onept5) -
     &                (three*dx*(dy*ddxdx-dx*ddydx)*(dx*ddxdx+dy*ddydx))/
     &                (dx**2+dy**2)**(twopt5))*(dxbds)**2 +
     &                (-ddydx/(dx**2 + dy**2)**pt5 + pt5*dy*(two*dx*ddxdx + 
     &                  two*dy*ddydx)/(dx**2 + dy**2)**onept5) * d2xbds2
              d2bn2 = ((ddydx*(dy*ddxdx-dx*ddydx)+dy*(dy*d2dxdx2-dx*d2dydx2))/
     &                (dx**2+dy**2)**(onept5) -
     &                (three*dy*(dy*ddxdx-dx*ddydx)*(dx*ddxdx+dy*ddydx))/
     &                (dx**2+dy**2)**(twopt5))*(dxbds)**2 +
     &                (ddxdx/(dx**2 + dy**2)**pt5 - pt5*dx*(two*dx*ddxdx + 
     &                 two*dy*ddydx)/(dx**2 + dy**2)**onept5) * d2xbds2
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

            if (j .eq. 1) then
              a = dxbds + 0.5 * dbn1
              b = dybds + 0.5 * dbn2
              h = sqrt( a**2 + b**2 )
              dads = d2xbds2 + 0.5 * d2bn1
              dbds = d2ybds2 + 0.5 * d2bn2
              dhds = ( a * dads + b * dbds ) / h
              write(44,30) ss(i), h, a, b, dads, dbds, dhds
              write(45,30) ss(i), h, d2xbds2, d2bn1, d2ybds2, d2bn2
            end if

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
c         call wgrid(x, y, nx, ny, nz, mx, my, 'grid.bin'//char(0))
c         call wdata(rm1, rn1, rm2, rn2, rjac, nx, ny, nz, mx, my, 
c     &               'm1.bin'//char(0))
c         call wdata(rm11, rm12, rm21, rm22, rjac, nx, ny, nz, mx, my, 
c     &               'm11.bin'//char(0))
c         call wdata(rn11, rn12, rn21, rn22, rjac, nx, ny, nz, mx, my, 
c     &               'n11.bin'//char(0))
c         call wdata(g11x, g11y, g12x, g12y, g22x, nx, ny, nz, mx, my, 
c     &               'grad.bin'//char(0))
        else
c
c.... write out unformatted if on a CRAY
c
          open(unit=10,file='grid.dat',form='unformatted')
          write(10) nx, ny, nz
          write(10) ((((x(i,j)), i=1,nx), j=1,ny), k = 1, nz),
     &              ((((y(i,j)), i=1,nx), j=1,ny), k = 1, nz),
     &              ((((zero),   i=1,nx), j=1,ny), k = 1, nz)
          close(10)
        
          open(unit=10,file='m1.dat',form='unformatted')
          write(10) nx, ny, nz
          write(10) 0.0, 0.0, 0.0, 0.0
          write(10) (((( rm1(i,j)), i=1,nx), j=1,ny), k = 1, nz),
     &              (((( rn1(i,j)), i=1,nx), j=1,ny), k = 1, nz),
     &              (((( rm2(i,j)), i=1,nx), j=1,ny), k = 1, nz),
     &              (((( rn2(i,j)), i=1,nx), j=1,ny), k = 1, nz),
     &              ((((rjac(i,j)), i=1,nx), j=1,ny), k = 1, nz)
          close(10)

          open(unit=10,file='m11.dat',form='unformatted')
          write(10) nx, ny, nz
          write(10) 0.0, 0.0, 0.0, 0.0
          write(10) ((((rm11(i,j)), i=1,nx), j=1,ny), k = 1, nz),
     &              ((((rm12(i,j)), i=1,nx), j=1,ny), k = 1, nz),
     &              ((((rm21(i,j)), i=1,nx), j=1,ny), k = 1, nz),
     &              ((((rm22(i,j)), i=1,nx), j=1,ny), k = 1, nz),
     &              ((((rjac(i,j)), i=1,nx), j=1,ny), k = 1, nz)
          close(10)

          open(unit=10,file='n11.dat',form='unformatted')
          write(10) nx, ny, nz
          write(10) 0.0, 0.0, 0.0, 0.0
          write(10) ((((rn11(i,j)), i=1,nx), j=1,ny), k = 1, nz),
     &              ((((rn12(i,j)), i=1,nx), j=1,ny), k = 1, nz),
     &              ((((rn21(i,j)), i=1,nx), j=1,ny), k = 1, nz),
     &              ((((rn22(i,j)), i=1,nx), j=1,ny), k = 1, nz),
     &              ((((rjac(i,j)), i=1,nx), j=1,ny), k = 1, nz)
          close(10)
c
c.... write out the metric file (Note the order of i and j are reversed)
c
          open (unit=10, file='metric.ji', form='unformatted', 
     &          status='unknown')
          write(10) (( ( rm1(i,j)), j=1,ny), i=1,nx),
     &              (( ( rm2(i,j)), j=1,ny), i=1,nx),
     &              (( ( rn1(i,j)), j=1,ny), i=1,nx),
     &              (( ( rn2(i,j)), j=1,ny), i=1,nx),
     &              (( (rm11(i,j)), j=1,ny), i=1,nx),
     &              (( (rm12(i,j)), j=1,ny), i=1,nx),
     &              (( (rm22(i,j)), j=1,ny), i=1,nx),
     &              (( (rn11(i,j)), j=1,ny), i=1,nx),
     &              (( (rn12(i,j)), j=1,ny), i=1,nx),
     &              (( (rn22(i,j)), j=1,ny), i=1,nx)
          close(10)
c
c.... Now write out in IJ order
c
          open (unit=10, file='metric.dat', form='unformatted', 
     &          status='unknown')
          write(10) (( ( rm1(i,j)), i=1,nx), j=1,ny),
     &              (( ( rm2(i,j)), i=1,nx), j=1,ny),
     &              (( ( rn1(i,j)), i=1,nx), j=1,ny),
     &              (( ( rn2(i,j)), i=1,nx), j=1,ny),
     &              (( (rm11(i,j)), i=1,nx), j=1,ny),
     &              (( (rm12(i,j)), i=1,nx), j=1,ny),
     &              (( (rm22(i,j)), i=1,nx), j=1,ny),
     &              (( (rn11(i,j)), i=1,nx), j=1,ny),
     &              (( (rn12(i,j)), i=1,nx), j=1,ny),
     &              (( (rn22(i,j)), i=1,nx), j=1,ny)
          close(10)

        end if

c.... write out the derivatives on the body

        write(*,"('Enter j ==> ',$)") 
        read(*,*) j
        do i = 1, nx
          write(50,25) xi(i), rm1(i,j), rm2(i,j)
          write(51,25) xi(i), rn1(i,j), rn2(i,j)
          write(52,25) xi(i), rm11(i,j), rn11(i,j), rm22(i,j), 
     .                 rm12(i,j), rn12(i,j)
          write(53,25) xi(i), rjac(i,j), one/rjac(i,j)
        end do

        stop      
 10     format('xmin = ',1pe13.6,' xmax = ',1pe13.6,' S = ',1pe13.6)
 20     format(i5,6(1x,1pe13.6))
 30     format(7(1x,1pe13.6))
 25     format(8(1pe13.6,1x))
        end

c=============================================================================c
        function yloc(x)
c=============================================================================c
c       implicit double precision (a-h,o-z)
        parameter (zero=0.0d0, pt5=0.5d0, one=1.0d0, two=2.0d0)

        common /stuff/ radius
c=============================================================================c
        yloc = sqrt( radius**2 - x**2 )

        return
        end

c=============================================================================c
        function xloc(x,ds)
c=============================================================================c
c       implicit double precision (a-h,o-z)
        parameter (zero=0.0d0, pt5=0.5d0, one=1.0d0, two=2.0d0)

        common /stuff/ radius
        common /distance/ darc, x1
c=============================================================================c
        darc = ds
        x1   = x

        th1 = atan2( sqrt( radius**2 - x1**2 ), x1 )
        th  = th1 - ds/radius
        xloc = radius * cos(th)

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
        common /stuff/ radius
c=============================================================================c
        th1 = atan2( sqrt( radius**2 - x1**2 ), x1 )
        th2 = atan2( sqrt( radius**2 - x2**2 ), x2 )

        arc = -radius*(th2 - th1)

        return
        end
