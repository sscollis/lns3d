c=============================================================================c
        program pc2
c=============================================================================c
c
c  Purpose:  make a body-fitted parabolic-cylinder mesh
c
c  Inputs:   This uses Mahesh mapping in x and y
c
c  Author:   S. Scott Collis
c
c  Revised:  10-1-96
c
c=============================================================================c
        implicit double precision (a-h,o-z)

        logical IRIS, DIAG
        parameter (IRIS = .false.)
        parameter (DIAG = .false.)

        parameter (zero=0.0d0,    pt5=0.5d0,  one=1.0d0, 
     &             onept5=1.5d0,  two=2.0d0,  twopt5=2.5d0, 
     &             three=3.0d0)
c
c.... must use a large mpts to make error in spline very small
c
        parameter (mpts=5000)
        dimension zz(mpts), s1(mpts), s2(mpts), ds1(mpts), ds2(mpts),
     &            dss1(mpts), dss2(mpts)

        parameter (mx=1024,my=127)
        
        dimension xb(mx), yb(mx), cur(mx)
        dimension bn1(mx), bn2(mx)
        
        dimension xi(mx),  eta(my)

        dimension ss(mx), rr(my), dss(mx), drr(my), d2ss(mx), d2rr(my)

        dimension x(mx,my),    y(mx,my)
        dimension rm1(mx,my),  rm2(mx,my)
        dimension rn1(mx,my),  rn2(mx,my),  rjac(mx,my)
        dimension rm11(mx,my), rm12(mx,my), rm21(mx,my), rm22(mx,my)
        dimension rn11(mx,my), rn12(mx,my), rn21(mx,my), rn22(mx,my)
        dimension g11x(mx,my), g11y(mx,my), g12x(mx,my), g12y(mx,my)
        dimension g22x(mx,my), g22y(mx,my)
        
        dimension d(10), e(10), f1(10), f2(10), g0(10)
        
        common /stuff/ AR, rm, rn, xmin, xmax

        external xloc, yloc

c.... argument parameters

        integer iarg, narg
        character(80)  arg
        integer yflag, xflag
        logical plot3d, debug
#ifndef __GFORTRAN__
        integer, external iargc
#endif
c=============================================================================c
        yflag = 1
        xflag = 0
        plot3d = .false.
        debug = .false.
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
          else if (arg(1:2) .eq. '-p') then
            plot3d = .true.
          else if (arg(1:2) .eq. '-d') then
            debug = .true.
          else if (arg(1:2) .eq. '-h') then
            write(*,"('Usage:  pc [options]')")
            write(*,"('   -h:  this help')")
            write(*,"('   -d:  debug output')")
            write(*,"('  -x1:  exponential stretching in x')")
            write(*,"('  -x2:  hyperbolic tangent stretching in x')")
            write(*,"('  -x3:  Mahesh''s mapping in x')")
            write(*,"('  -x4:  Mahesh''s mapping in s')")
            write(*,"('  -y1:  exponential stretching in y')")
            write(*,"('  -y2:  hyperbolic tangent stretching in y')")
            write(*,"('  -y3:  Mahesh''s mapping in y')")
            call exit(0)
          else
            write(*,"('Argument ',i2,' ignored.')") iarg
          end if
        end do
c
c.... get user input
c
        write(*,"('Enter (Nx, Ny) ==> ',$)")
        read(*,*) nx, ny
        if (nx.gt.mx) stop 'nx > mx'
        if (ny.gt.my) stop 'ny > my'

        write(*,*) 'Discretize in the wall tangent direction'

        write(*,"('Enter xmin, xmax ==> ',$)")
        read(*,*) xmin, xmax
c
c.... discretize in the body tangent direction
c
        smin = arc(zero,xmin)
        smax = arc(zero,xmax)
        savg = (smin + smax) * pt5
        write(*,10) xmin, xmax, smin, smax, smax-smin

        dxi  = one / dble(nx-1)
        ds   = (smax-smin) / dble(nx-1)

        write(*,"('ds uniform = ',1pe20.13)") ds
c
c.... Uniform arc-length
c
        if (xflag.eq.0) then
          do i = 1, nx
            xi(i)   = dble(i-1) * dxi
            ss(i)   = smin + float(i-1) * ds
            dss(i)  = ds / dxi
            d2ss(i) = zero
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
        write(*,"('Enter dsmin, b, sc ==> ',$)") 
        read(*,*) dsmin, b, sc

        cm = ( two * b * tanh(b*sc) + real(nx-1)*dsmin/smax * 
     &         log( cosh(b*(one-sc)) / cosh(b*(one+sc)) ) ) /
     &       ( one - real(nx-1)*dsmin/smax )

        do i = 1, nx
           xi(i)   = dble(i-1) * dxi
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

        deta = one / dble(ny-1)
        rmin = zero
        write(*,"('Enter Rmax ==> ',$)")
        read(*,*) rmax
        dr = (rmax-rmin) / real(ny-1)
c
c.... uniform mesh in the wall normal direction
c
        if (yflag.eq.0) then
          do j = 1, ny
             eta(j)  = dble(j-1) * deta
             rr(j)   = rmin + dble(j-1) * dr
             drr(i)  = dr / deta
             d2rr(i) = zero
          end do
        end if
c
c.... exponential map in the wall normal direction
c
        if (yflag.eq.1) then
        deta2 = one / real(128-1)
        dymin = 5.0d-5
        sy = 1.088d0
        c2 = log( dymin / deta2 )
        c1 = (log( sy*dymin / deta2) - c2) / deta2
        do j = 1, ny
          eta(j)  = dble(j-1) * deta
          rr(j)   = one / c1 * ( exp(c1 * eta(j) + c2) - exp(c2) ) 
          drr(j)  = exp(c1 * eta(j) + c2)
          d2rr(j) = c1 * exp(c1 * eta(j) + c2)
          write(20,30) eta(j), rr(j), drr(j), d2rr(j)
        end do
        end if
c
c.... Hyperbolic tangent stretching
c
        if (yflag.eq.2) then
        
        rd1 = 0.00028d0                 ! set for R=1000 PCYL, r=500
        rd2 = 5.5d0
        dd  = 5.665549493763092d0       

        do j = 1, ny
          eta(j) = dble(j-1) * deta
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
c.... Mahesh mapping, personal communication (6-5-96)
c
        if (yflag.eq.3) then
          write(*,"('Enter drmin, b, rc ==> ',$)") 
          read(*,*) drmin, b, rc
          cm = ( two * b * tanh(b*rc) + (ny-1)*drmin/rmax * 
     &           log( cosh(b*(one-rc)) / cosh(b*(one+rc)) ) ) /
     &         ( one - (ny-1)*drmin/rmax )
          do j = 1, ny
            eta(j)  = dble(j-1) * deta
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
!            xl = xloc( xmin, ss(i) )
             xl = xmin
          else
             xl = xloc( xold, ss(i) - ss(i-1) )
          end if
          yl = yloc( xl )
          th = atan2( one, yl )
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
            x(i,j) = xb(i) + rr(j) * bn1(i)
            y(i,j) = yb(i) + rr(j) * bn2(i)
          end do
        end do
        write(*,*) 'Making the metrics'
c
c.... write out the body.dat file
c.... rk is the curvature for a parabolic cylinder
c
        open(9,file='body.dat',status='unknown')
        j = 1
        do i = 1, nx
           rk = (two * xb(i) + one)**(-onept5)
           write(9,"(1(i5,1x),6(1pe20.13,1x))") i, ss(i), xb(i), yb(i),
     &                                          bn1(i), bn2(i), rk
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
        dr = sqrt( (x(1,1) - x(1,2))**2 + (y(1,1) - y(1,2))**2 )
        yl = y(2,1)
        ds = yl * 0.5 * sqrt( 1.0 + yl**2) + 
     &       0.5 * log( abs( yl + sqrt(1.0 + yl**2)))

        write(*,11) ds, dr, sqrt( (x(1,2)-x(1,3))**2 + (y(1,2)-y(1,3))**2 )/dr
 11     format('ds_min = ',1pe20.13,' dr_min = ',1pe20.13,
     &         ' dr(2)/dr(1)     = ',1pe20.13)
c
c.... write out the maximum spacing
c
        dr = x(1,ny-1) - x(1,ny)
        yl = y(nx,1)
        ds = yl * 0.5 * sqrt( 1.0 + yl**2) + 
     &       0.5 * log( abs( yl + sqrt(1.0 + yl**2)))
        yl = y(nx-1,1)
        ds = ds - (yl * 0.5 * sqrt( 1.0 + yl**2) + 
     &       0.5 * log( abs( yl + sqrt(1.0 + yl**2))))
        write(*,21) ds, dr, (x(1,ny-2)-x(1,ny-1))/dr
 21     format('ds_max = ',1pe20.13,' dr_max = ',1pe20.13,
     &         ' dr(ny)/dr(ny-1) = ',1pe20.13)
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
            xtmp = one

            if (xl .eq. zero) then
              dydx = 1.0d30
            else
              dydx = one / sqrt( two * xl )
            end if
            dxbds  = one / sqrt(one + dydx**2)
  
            if (xl .eq. zero) then
              dxdy = zero
            else
              dxdy = sqrt( two * xl )
            end if
            dybds = one / sqrt( dxdy**2 + one )
  
            if (xl .eq. zero) then
              dx = zero
              dy = one
              
              ddxdx = 1.0d30
              ddydx = zero

              ddxdy = one
              ddydy = zero

              d2dxdx2 = 1.0d30
              d2dydx2 = zero

              d2dxdy2 = zero
              d2dydy2 = zero
            else
              dx = yl
              dy = one

              ddxdx = one / yl
              ddydx = zero

              ddxdy = one
              ddydy = zero

              d2dxdx2 = -sqrt(two)/(two**2 * xl**onept5)
              d2dydx2 = zero

              d2dxdy2 = zero
              d2dydy2 = zero
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

            if (j.eq.1) then
               cur(i) = sqrt(dbn1**2+dbn2**2)
               write(35,30) ss(i), sqrt(dbn1**2+dbn2**2)
            end if
c==============================================================================
c.... second derivative metrics
c==============================================================================
            if ( abs(bn1(i)) .gt. abs(bn2(i)) ) then
              d2xdy2  = one
              d2ybds2 = -(one + dxdy**2)**(-onept5) * dxdy * d2xdy2 * dybds
              d2xbds2 = d2xdy2*(dybds)**2 + dxdy*d2ybds2
            else
              d2ydx2  = -one / ( two * sqrt(two) * xl**onept5 )
              d2xbds2 = -(one + dydx**2)**(-onept5) * dydx * d2ydx2 * dxbds
              d2ybds2 = d2ydx2*(dxbds)**2 + dydx*d2xbds2
            end if    

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
              a = dxbds + one * dbn1
              b = dybds + one * dbn2
              h = sqrt( a**2 + b**2 )
              dads = d2xbds2 + one * d2bn1
              dbds = d2ybds2 + one * d2bn2
              dhds = ( a * dads + b * dbds ) / h
              write(44,30) ss(i), h, a, b, dads, dbds, dhds
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
        write(*,*) 'Error in m ',sqrt(errm)/dble(nx*ny)
        write(*,*) 'Error in n ',sqrt(errn)/dble(nx*ny)
c
c..... diagnostic
c
        if (debug) then
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
        end if

c=============================================================================
c       O u t p u t   R e s u l t s
c=============================================================================
c
c.... output in Streett format
c
        open(10,file='field.mean',form='unformatted')
        write(10) nx, ny
        write(10) (((0.0, j=1, ny), i=1,nx), k=1,4)
        write(10) (sngl(ss(i)), i=1,nx)
        write(10) (sngl(rr(j)), j=1,ny)
        write(10) (sngl(cur(i)), i=1,nx)
        close(10)
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
          write(10) (((sngl(x(i,j)), i=1,nx), j=1,ny), k = 1, nz),
     &              (((sngl(y(i,j)), i=1,nx), j=1,ny), k = 1, nz),
     &              (((sngl(zero),   i=1,nx), j=1,ny), k = 1, nz)
          close(10)
        
          if (plot3d) then

          open(unit=10,file='m1.dat',form='unformatted')
          write(10) nx, ny, nz
          write(10) 0.0, 0.0, 0.0, 0.0
          write(10) (((sngl( rm1(i,j)), i=1,nx), j=1,ny), k = 1, nz),
     &              (((sngl( rn1(i,j)), i=1,nx), j=1,ny), k = 1, nz),
     &              (((sngl( rm2(i,j)), i=1,nx), j=1,ny), k = 1, nz),
     &              (((sngl( rn2(i,j)), i=1,nx), j=1,ny), k = 1, nz),
     &              (((sngl(rjac(i,j)), i=1,nx), j=1,ny), k = 1, nz)
          close(10)

          open(unit=10,file='m11.dat',form='unformatted')
          write(10) nx, ny, nz
          write(10) 0.0, 0.0, 0.0, 0.0
          write(10) (((sngl(rm11(i,j)), i=1,nx), j=1,ny), k = 1, nz),
     &              (((sngl(rm12(i,j)), i=1,nx), j=1,ny), k = 1, nz),
     &              (((sngl(rm21(i,j)), i=1,nx), j=1,ny), k = 1, nz),
     &              (((sngl(rm22(i,j)), i=1,nx), j=1,ny), k = 1, nz),
     &              (((sngl(rm1(i,j)*rn1(i,j)+rm2(i,j)*rn2(i,j)), 
     &              i=1,nx), j=1,ny), k = 1, nz)
          close(10)

          open(unit=10,file='n11.dat',form='unformatted')
          write(10) nx, ny, nz
          write(10) 0.0, 0.0, 0.0, 0.0
          write(10) (((sngl(rn11(i,j)), i=1,nx), j=1,ny), k = 1, nz),
     &              (((sngl(rn12(i,j)), i=1,nx), j=1,ny), k = 1, nz),
     &              (((sngl(rn21(i,j)), i=1,nx), j=1,ny), k = 1, nz),
     &              (((sngl(rn22(i,j)), i=1,nx), j=1,ny), k = 1, nz),
     &              (((sngl(rjac(i,j)), i=1,nx), j=1,ny), k = 1, nz)
          close(10)
          
          end if
c
c.... write out the metric file (Note the order of i and j are reversed)
c
          open (unit=10, file='metric.dat', form='unformatted', 
     &          status='unknown')
          write(10) (( sngl( rm1(i,j)), j=1,ny), i=1,nx),
     &              (( sngl( rm2(i,j)), j=1,ny), i=1,nx),
     &              (( sngl( rn1(i,j)), j=1,ny), i=1,nx),
     &              (( sngl( rn2(i,j)), j=1,ny), i=1,nx),
     &              (( sngl(rm11(i,j)), j=1,ny), i=1,nx),
     &              (( sngl(rm12(i,j)), j=1,ny), i=1,nx),
     &              (( sngl(rm22(i,j)), j=1,ny), i=1,nx),
     &              (( sngl(rn11(i,j)), j=1,ny), i=1,nx),
     &              (( sngl(rn12(i,j)), j=1,ny), i=1,nx),
     &              (( sngl(rn22(i,j)), j=1,ny), i=1,nx)
          close(10)
        end if

        stop      
 10     format('Xmin = ',1pe13.6,' Xmax = ',1pe13.6,
     &         ' Smin = ',1pe13.6,' Smax = ',1pe13.6,' S = ',1pe13.6)
 20     format(i5,6(1x,1pe13.6))
 30     format(7(1x,1pe20.13))
        end

c=============================================================================c
        function yloc(x)
c=============================================================================c
        implicit double precision (a-h,o-z)
        parameter (zero=0.0d0, pt5=0.5d0, one=1.0d0, two=2.0d0)

        common /stuff/ AR, rm, rn, xmin, xmax
c=============================================================================c
        yloc = sqrt( two * x )

        return
        end

c=============================================================================c
        function xloc(x,ds)
c=============================================================================c
        implicit double precision (a-h,o-z)
        parameter (zero=0.0d0, pt5=0.5d0, one=1.0d0, two=2.0d0)

        external rtsec, rtflsp, func

        common /stuff/ AR, rm, rn, xmin, xmax
        common /distance/ darc, x1
c=============================================================================c
        darc = ds
        x1   = x

c       xloc = rtsec(func,x1,xmax,1.0d-8)

        xloc = rtflsp(func,x1,x1+two*ds,1.0d-8)

c       xloc = zbrent(func,x1,x1+two*ds,1.0d-8)

        return
        end 

c=============================================================================c
        function func(x)
c=============================================================================c
        implicit double precision (a-h,o-z)
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
        implicit double precision (a-h,o-z)
c=============================================================================c
        xi1 = sqrt(x1)
        xi2 = sqrt(x2)


        arc = sqrt(2.0)*0.5*xi2*sqrt(1.0+2.0*xi2**2) +
     &        0.5*log(sqrt(2.0)*xi2 + sqrt(1.0+2.0*xi2**2))

        arc = arc - ( sqrt(2.0)*0.5*xi1*sqrt(1.0+2.0*xi1**2) +
     &        0.5*log(sqrt(2.0)*xi1 + sqrt(1.0+2.0*xi1**2)) )

        return
        end

c=============================================================================c
c       function arc(x1,x2)
c=============================================================================c
c       implicit double precision (a-h,o-z)
c       parameter (zero=0.0d0, pt5=0.5d0, one=1.0d0, two=2.0d0)
c
c       external derivs1, derivs2, RKQCR
c
c       common /stuff/ AR, rm, rn, xmin, xmax
c=============================================================================c
c       s = zero
c       xtmp = one
c
c       if (x1 .eq. x2 ) then
c          arc = zero
c          return
c       end if
c
c       if (x1 .gt. x2) then
c          write(*,*) x1, x2
c          stop 'Error in arc:  x1 > x2'
c       end if
c       if (x1 .lt. xmin) stop 'Error in arc:  x1 < xmin'
c       if (x2 .gt. xmax) stop 'Error in arc:  x2 > xmax'
c
c       if (x1 .le. xtmp) then
c           y1 = sqrt( two * x1 )
c          if (x2 .le. xtmp) then
c             y2 = sqrt( two * x2 )
c          else
c             y2 = sqrt( two * xtmp )
c          end if
c           call ODEINTR(s,1,y1,y2,1.0d-15,(y2-y1)*pt5,
c     &                  1.d-15,nok,nbad,derivs2,RKQCR)
c          if (x2 .gt. xtmp) then
c             call ODEINTR(s,1,xtmp,x2,1.0d-15,(x2-x1)*pt5,
c     &                    1.d-15,nok,nbad,derivs1,RKQCR)
c          end if
c       else
c          call ODEINTR(s,1,x1,x2,1.0d-15,(x2-x1)*pt5,
c     &                 1.d-15,nok,nbad,derivs1,RKQCR)
c       end if
c
c       arc = s
c
c       return
c       end
c
c=============================================================================c
        subroutine derivs1(n,s,x,ds)
c=============================================================================c
        implicit double precision (a-h,o-z)
        parameter (zero=0.0d0, pt5=0.5d0, one=1.0d0, two=2.0d0)

        common /stuff/ AR, rm, rn, xmin, xmax
c=============================================================================c
        dydx = one / sqrt( two * x )

        ds = sqrt( one + dydx**2 )

        return
        end

c=============================================================================c
        subroutine derivs2(n,s,y,ds)
c=============================================================================c
        implicit double precision (a-h,o-z)
        parameter (zero=0.0d0, pt5=0.5d0, one=1.0d0, two=2.0d0)

        common /stuff/ AR, rm, rn, xmin, xmax
c=============================================================================c
        x = pt5 * y**2
        
        dxdy = y

        ds = sqrt( dxdy**2 + one )

        return
        end
