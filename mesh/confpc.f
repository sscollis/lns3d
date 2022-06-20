c=============================================================================c
        program confpc
c=============================================================================c
c
c  Purpose:  make parabolic-cylinder mesh using conformal mapping
c
c  Author:   S. Scott Collis
c
c  Revised:  8-7-96
c
c  (1) Fixed a bug in the second derivative metric for Mahesh mapping
c
c  (2) Added the capability of doing mapping in s as opposed to \xi
c
c  Revised:  8-12-96
c
c  (1) Added a virtual orgin option so that I can do a downstream section
c      of the plate.
c
c  (2) Cleaned up some extaneous stuff
c
c   Make sure to compile with -N80 on the Cray
c=============================================================================c
        implicit double precision (a-h,o-z)

        logical IRIS, DIAG
        parameter (IRIS = .false.)
        
        parameter (mx=1024,my=1024)


        parameter (zero   = 0.0d0,  pt5 = 0.5d0,  one    = 1.0d0, 
     &             onept5 = 1.5d0,  two = 2.0d0,  twopt5 = 2.5d0, 
     &             three  = 3.0d0)
        
        dimension x(mx,my), y(mx,my), xi(mx), eta(my), s(mx)
        dimension ss(mx), dss(mx), d2ss(mx)
        dimension rr(my), drr(my), d2rr(my)
        dimension rm1(mx,my),  rm2(mx,my)
        dimension rn1(mx,my),  rn2(mx,my)
        dimension rjac(mx,my)
        dimension rm11(mx,my), rm12(mx,my), rm21(mx,my), rm22(mx,my)
        dimension rn11(mx,my), rn12(mx,my), rn21(mx,my), rn22(mx,my)
        dimension g11x(mx,my), g11y(mx,my), g12x(mx,my), g12y(mx,my)
        dimension g22x(mx,my), g22y(mx,my)

        external calcdd, calcs

!.... argument parameters

        integer narg
        character(80) arg
        integer yflag, xflag
        logical plot3d, debug, useIJ
#ifndef __GFORTRAN__
        integer, external iargc
#endif
        external arc, calcxi
c=============================================================================c
        yflag = 1
        xflag = 0
        plot3d = .false.
        debug = .false.
        useIJ = .true.
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
          else if (arg(1:3) .eq. '-ji') then
            useIJ = .false.
          else if (arg(1:2) .eq. '-h') then
            write(*,"('Usage:  confpc [options]')")
            write(*,"('   -h:  this help')")
            write(*,"('   -d:  debug output')")
            write(*,"('  -x1:  exponential stretching in x')")
            write(*,"('  -x2:  hyperbolic tangent stretching in x')")
            write(*,"('  -x3:  Mahesh''s mapping in x')")
            write(*,"('  -x4:  Mahesh''s mapping in s')")
            write(*,"('  -y1:  exponential stretching in y')")
            write(*,"('  -y2:  hyperbolic tangent stretching in y')")
            write(*,"('  -y3:  Mahesh''s mapping in y')")
            write(*,"('  -ji:  Output metrics in JI ordering')")
            call exit(0)
          else
            write(*,"('Argument ',i2,' ignored.')") iarg
          end if
        end do
c
c.... get user input
c
        write(*,"('Enter (nx,ny) ==> ',$)") 
        read(*,*) nx, ny
        if (nx.gt.mx) then
           write(*,*) 'Nx is greater than Mx ',mx
           call exit(1)
        end if
        if (ny.gt.my) then
           write(*,*) 'Ny is greater than My ',my
           call exit(1)
        end if
        write(*,"('Enter (Xmin, Xmax) ==> ',$)")
        read(*,*) qLx, rLx
        write(*,"('Enter (Ymax) ==> ',$)")
        read(*,*) rLy
        qLy = zero
c
c.... uniform grid in computational space
c       
        ximax = sqrt(rLx)
        ximin = sqrt(qLx)
        dxi = (ximax - ximin) / dble(nx-1)
        write(*,*) 'Uniform dxi = ',dxi
        ds  = one / dble(nx-1)
        do i = 1, nx
          ss(i)   = dble(i-1) * ds
          xi(i)   = ximin + dble(i-1) * dxi
          dss(i)  = (ximax - ximin)
          d2ss(i) = zero
        end do

        etamax = sqrt(one + two * rLy)
        etamin = sqrt(one + two * qLy)
        deta = (etamax - etamin) / dble(ny-1)
        write(*,*) 'Uniform deta = ', deta
        write(*,*) 'Uniform dy = ', pt5*((one+deta)**2-one)
        dr   = one / dble(ny-1)
        do j = 1, ny
          rr(j)   = dble(j-1) * dr
          eta(j)  = etamin + dble(j-1) * deta
          drr(j)  = (etamax - etamin)
          d2rr(j) = zero
        end do
c
c.... exponential map in the streamwise direction
c
        if (xflag.eq.1) then
          write(*,"('Enter dxmin ==> ',$)") 
          read(*,*) dxmin
          if (dxmin.ne.zero) then
            sx = calcs( ximin, ximax, dxmin, ds )
            write(*,*) 'Sx = ', sx
            c2 = log( dxmin / ds )
            c1 = (log( sx * dxmin / ds) - c2) / ds
            do i = 1, nx
              ss(i)   = dble(i-1) * ds
              xi(i)   = ximin + one / c1 * ( exp(c1 * ss(i) + c2) - 
     &                  exp(c2) ) 
              dss(i)  = exp(c1 * ss(i) + c2)
              d2ss(i) = c1 * exp(c1 * ss(i) + c2)
            end do
            write(*,*) 'Xmax = ', xi(nx)**2
          else
            write(*,*) 'Uniform mesh in xi'
          end if
        end if
c
c.... hyperbolic tangent mapping in the streamwise direction
c.... this increases resolution in the sponge
c
        if (xflag.eq.2) then
        smax = (ximax - ximin)
        sd1 = 0.4d0             ! set for R=1000 PCYL, Lx=500
        sd2 = 0.6d0
        dd = calcdd(sd1,sd2)
        write(*,*) 'dd = ',dd
        dxi = one / dble(nx-1)
        do i = 1, nx
          ss(i) =  dble(i-1) * ds
          xi(i) =  ximin + smax*(pt5 + one/Tanh(dd*pt5)*Tanh(dd*(-pt5 +
     &             ss(i)))*pt5)/(Sqrt(sd2/sd1) +
     &             (one - Sqrt(sd2/sd1))*
     &             (pt5 + one/Tanh(dd*pt5)*Tanh(dd*(-pt5 +
     &             ss(i)))*pt5))
          dss(i) = dd*Sqrt(sd2/sd1)*smax*Cosh(dd*pt5)*
     &             one/Cosh(dd*(-pt5 + ss(i)))*
     &             (Sinh(dd*(one - ss(i))) + Sinh(dd*ss(i)))/
     &             (Sqrt(sd2/sd1)*Sinh(dd*(one - ss(i))) + 
     &             Sinh(dd*ss(i)))**2
          d2ss(i) = dd**2*Sqrt(sd2/sd1)*smax*
     &              (-Cosh(dd*(one - three*ss(i))) + 
     &              Sqrt(sd2/sd1)*Cosh(dd*(two - three*ss(i))) - 
     &              Cosh(dd*(one - ss(i))) + 
     &              two*Sqrt(sd2/sd1)*Cosh(dd*(one - ss(i))) - 
     &              two*Cosh(dd*ss(i)) + 
     &              Sqrt(sd2/sd1)*Cosh(dd*ss(i)))*
     &              (one/Cosh(dd*(-pt5 + ss(i))))**2*(Sinh(dd))/
     &              (two*(Sqrt(sd2/sd1)*(Sinh(dd*(one - ss(i)))) + 
     &              (Sinh(dd*ss(i))))**3)
        end do
        end if
c
c.... Mahesh's mapping, personal communication (6-5-96)
c
        if (xflag.eq.3) then
          write(*,"('Enter dsmin, b, sc ==> ',$)") 
          read(*,*) dsmin, b, sc
          dsmin = dsmin / sqrt(two)
          smax = (ximax - ximin)
          cm = ( two * b * tanh(b*sc) + (nx-1)*dsmin/smax * 
     &           log( cosh(b*(one-sc)) / cosh(b*(one+sc)) ) ) /
     &         ( one - (nx-1)*dsmin/smax )
          do i = 1, nx
            ss(i)   = dble(i-1) * ds
            xi(i)   = ximin + smax*( cm * ss(i) + 
     &                log( cosh(b*(ss(i)-sc)) /
     &                cosh(b*(ss(i)+sc)) ) ) / 
     &                (cm + log( cosh(b*(one-sc)) / cosh(b*(one+sc)) ) )
            dss(i)  = smax*(cm + b*tanh(b*(ss(i)-sc)) - 
     &                b*tanh(b*(ss(i)+sc)))/
     &                (cm + log( cosh(b*(one-sc)) / cosh(b*(one+sc)) ) )
            d2ss(i) = smax*(-b**2*(tanh(b*(ss(i)-sc)))**2 + 
     &                b**2*(tanh(b*(ss(i)+sc)))**2)/
     &                (cm + log( cosh(b*(one-sc)) / cosh(b*(one+sc)) ) )
          end do
        end if
c
c.... Uniform arc-length mapping
c
        if (xflag.eq.4) then
          smin  = arc(ximin)
          smax  = arc(ximax)
          write(*,"('smin, smax, ds uniform = ',3(1pe20.13,1x))") 
     &      smin, smax, (smax-smin)/real(nx-1)

          write(*,"('Enter dsmin, b, sc ==> ',$)") 
          read(*,*) dsmin, b, sc

          cm = ( two * b * tanh(b*sc) + (nx-1)*dsmin/(smax-smin) * 
     &           log( cosh(b*(one-sc)) / cosh(b*(one+sc)) ) ) /
     &         ( one - (nx-1)*dsmin/(smax-smin) )
          do i = 1, nx
            ss(i)   = dble(i-1) * ds
            s(i)    = smin + (smax-smin)*( cm * ss(i) + 
     &                log( cosh(b*(ss(i)-sc)) /
     &                cosh(b*(ss(i)+sc)) ) ) / 
     &                (cm + log( cosh(b*(one-sc)) / cosh(b*(one+sc)) ) )
            dss(i)  = (smax-smin)*(cm - b*(1.0/Cosh(b*(ss(i) - sc)))*
     &                (1.0/Cosh(b*(ss(i) + sc)))*Sinh(2.0*b*sc)) /
     &                (cm + log( cosh(b*(one-sc)) / cosh(b*(one+sc)) ) )
            d2ss(i) = (smax-smin)*(-b**2*(tanh(b*(ss(i)-sc)))**2 + 
     &                b**2*(tanh(b*(ss(i)+sc)))**2)/
     &                (cm + log( cosh(b*(one-sc)) / cosh(b*(one+sc)) ) )
          end do
c
c.... transform from s to \xi
c
          do i = 1, nx
            xi(i)   = calcxi( s(i) )
            dxids   = 1.0 / (
     &                Sqrt(2.)*xi(i)**2/Sqrt(1. + 2.*xi(i)**2) + 
     &                Sqrt(1. + 2.*xi(i)**2)/Sqrt(2.) +
     &                (Sqrt(2.) + 2.*xi(i)/Sqrt(1. + 2.*xi(i)**2))/
     &                (2.*(Sqrt(2.)*xi(i) + Sqrt(1. + 2.*xi(i)**2))) )
            d2xids2 = (
     &                -2.*xi(i)*(Sqrt(2.) + 2.**(5./2.)*xi(i)**2 + 
     &                4.*xi(i)*Sqrt(1.0 + 2.0*xi(i)**2))/
     &                ((1.0 + 2.0*xi(i)**2)**(3.0/2.0)*(2.0*xi(i) + 
     &                Sqrt(2.0)*Sqrt(1.0 + 2.0*xi(i)**2))**2) ) * dxids
            d2ss(i) = dss(i)**2 * d2xids2 + d2ss(i) * dxids
            dss(i)  = dss(i) * dxids
          end do
        end if
c
c.... output the tangent mapping
c
        if (debug) then
          open(10,file='tangent.dat')
          do i = 1, nx
            write(10,30) ss(i), xi(i), dss(i), d2ss(i), d2ss(i)/dss(i)
          end do
          close(10)
        end if
c
c.... exponential map in the wall normal direction
c
        if (yflag.eq.1) then
          rmin = qLy
          rmax = rLy
          dymin = 1.0d-2
          sy = 1.01d0
          write(*,"('Enter dymin ==> ',$)") 
          read(*,*) dymin
          if (dymin.ne.zero) then
            sy = calcs( rmin, rmax, dymin, dr )
            write(*,*) 'Sy = ', sy
            c2 = log( dymin / dr )
            c1 = (log( sy*dymin / dr) - c2) / dr
            do j = 1, ny
              rr(j)   = dble(j-1) * dr
              eta(j)  = one / c1 * ( exp(c1 * rr(j) + c2) - exp(c2) ) 
              drr(j)  = exp(c1 * rr(j) + c2)
              d2rr(j) = c1 * exp(c1 * rr(j) + c2)
            end do
            write(*,*) 'Rmax = ', eta(ny)
            do j = 1, ny
              d2rr(j) = -(drr(j)**2)/( 1.0 + 2.0 * eta(j) )**(1.5) + 
     &                  d2rr(j) / sqrt( 1.0 + 2.0 * eta(j) )
              drr(j)  = drr(j) / sqrt( 1.0 + 2.0 * eta(j) )
              eta(j)  = sqrt( 1.0 + 2.0 * eta(j) )
            end do
          else
            write(*,*) 'Uniform mesh in eta'
          end if
        end if
c
c.... hyperbolic tangent stretching in wall normal direction
c
        if (yflag.eq.2) then
        rmin = qLy
        rmax = rLy
        rd1 = 0.00028d0         ! set for R=1000 PCYL, r=500
        rd2 = 5.5d0
        dd = calcdd(rd1,rd2)
        write(*,*) 'dd = ',dd
        deta = one / dble(ny-1)
        do j = 1, ny
          rr(j) = dble(j-1) * dr
          eta(j) = rmax*(pt5 + one/Tanh(dd*pt5)*Tanh(dd*(-pt5 +
     &             rr(j)))*pt5)/(Sqrt(rd2/rd1) +
     &             (one - Sqrt(rd2/rd1))*
     &             (pt5 + one/Tanh(dd*pt5)*Tanh(dd*(-pt5 +
     &             rr(j)))*pt5))
          drr(j) = dd*Sqrt(rd2/rd1)*rmax*Cosh(dd*pt5)*
     &             one/Cosh(dd*(-pt5 + rr(j)))*
     &             (Sinh(dd*(one - rr(j))) + Sinh(dd*rr(j)))/
     &             (Sqrt(rd2/rd1)*Sinh(dd*(one - rr(j))) + 
     &             Sinh(dd*rr(j)))**2
          d2rr(j) = dd**2*Sqrt(rd2/rd1)*rmax*
     &              (-Cosh(dd*(one - three*rr(j))) + 
     &              Sqrt(rd2/rd1)*Cosh(dd*(two - three*rr(j))) - 
     &              Cosh(dd*(one - rr(j))) + 
     &              two*Sqrt(rd2/rd1)*Cosh(dd*(one - rr(j))) - 
     &              two*Cosh(dd*rr(j)) + 
     &              Sqrt(rd2/rd1)*Cosh(dd*rr(j)))*
     &              (one/Cosh(dd*(-pt5 + rr(j))))**2*(Sinh(dd))/
     &              (two*(Sqrt(rd2/rd1)*(Sinh(dd*(one - rr(j)))) + 
     &              (Sinh(dd*rr(j))))**3)
        end do
        write(*,*) 'Rmax = ', eta(ny)
        do j = 1, ny
          d2rr(j) = -(drr(j)**2)/( 1.0 + 2.0 * eta(j) )**(1.5) + 
     &              d2rr(j) / sqrt( 1.0 + 2.0 * eta(j) )
          drr(j)  = drr(j) / sqrt( 1.0 + 2.0 * eta(j) )
          eta(j)  = sqrt( 1.0 + 2.0 * eta(j) )
        end do
        end if
c
c.... Mahesh's mapping, personal communication (6-5-96)
c
        if (yflag.eq.3) then
          write(*,"('Enter drmin, b, rc ==> ',$)") 
          read(*,*) drmin, b, rc
          rmax = rLy
          cm = ( two * b * tanh(b*rc) + (ny-1)*drmin/rmax * 
     &           log( cosh(b*(one-rc)) / cosh(b*(one+rc)) ) ) /
     &         ( one - (ny-1)*drmin/rmax )
          do j = 1, ny
            rr(j)   = dble(j-1) * dr
            eta(j)  = rmax*( cm * rr(j) + log( cosh(b*(rr(j)-rc)) /
     &                cosh(b*(rr(j)+rc)) ) ) / 
     &                (cm + log( cosh(b*(one-rc)) / cosh(b*(one+rc)) ) )
            drr(j)  = rmax*(cm + b*tanh(b*(rr(j)-rc)) - 
     &                b*tanh(b*(rr(j)+rc))) /
     &                (cm + log( cosh(b*(one-rc)) / cosh(b*(one+rc)) ) )
            d2rr(j) = rmax*(-b**2*(tanh(b*(rr(j)-rc)))**2 + 
     &                b**2*(tanh(b*(rr(j)+rc)))**2)/
     &                (cm + log( cosh(b*(one-rc)) / cosh(b*(one+rc)) ) )
          end do
c
c.... transform from r to \eta
c
          do j = 1, ny
            d2rr(j) = -(drr(j)**2)/( one + two * eta(j) )**onept5 + 
     &                d2rr(j) / sqrt( one + two * eta(j) )
            drr(j)  = drr(j) / sqrt( one + two * eta(j) )
            eta(j)  = sqrt( one + two * eta(j) )
          end do
        end if
c
c.... output the normal mapping
c
        if (debug) then
          open(20,file='normal.dat')
          do j = 1, ny
             write(20,30) rr(j), eta(j), drr(j), d2rr(j), d2rr(j)/drr(j)
          end do
          close(20)
        end if
c       
c.... build the mesh
c       
        do i = 1, nx
          do j = 1, ny
            x(i,j) = xi(i)**2 + (one - eta(j)**2) * pt5
            y(i,j) = sqrt(two) * xi(i) * eta(j)

            if (i.eq.1 .and. j.eq.1 .and. qLx.eq.zero) then
              x(i,j) = zero
              y(i,j) = zero
            end if

            dsdxi  = dss(i)
            drdxi  = zero

            dsdeta = zero
            drdeta = drr(j)

            dxds  = two * xi(i)
            dxdr  = -eta(j)
            
            dyds  = sqrt(two) * eta(j)
            dydr  = sqrt(two) * xi(i)
            
            dxdxi  = dxds * dsdxi  + dxdr * drdxi
            dxdeta = dxds * dsdeta + dxdr * drdeta

            dydxi  = dyds * dsdxi  + dydr * drdxi
            dydeta = dyds * dsdeta + dydr * drdeta

            rjac(i,j) = one / ( dxdxi * dydeta - dxdeta * dydxi )
            
            rm1(i,j) =  dydeta * rjac(i,j)
            rm2(i,j) = -dxdeta * rjac(i,j)
            
            rn1(i,j) = -dydxi * rjac(i,j)
            rn2(i,j) =  dxdxi * rjac(i,j)

            d2sdxi2    = d2ss(i)
            d2sdxideta = zero
            d2sdeta2   = zero

            d2rdxi2    = zero
            d2rdxideta = zero
            d2rdeta2   = d2rr(j)

            d2xds2  = two
            d2xdsdr = zero
            d2xdr2  = -one
            
            d2yds2  = zero
            d2ydsdr = sqrt(two)
            d2ydr2  = zero
            
            d2xdxi2 = d2xds2 * dsdxi**2 + dxds * d2sdxi2 + 
     &                d2xdsdr * drdxi * dsdxi + 
     &                d2xdsdr * drdxi * dsdxi +
     &                d2xdr2 * drdxi**2 + dxdr * d2rdxi2
     
            d2xdxideta = d2xds2 * dsdeta * dsdxi + dxds * d2sdxideta + 
     &                   d2xdsdr * drdeta * dsdxi + 
     &                   d2xdsdr * drdxi * dsdeta +
     &                   d2xdr2 * drdeta * drdxi + dxdr * d2rdxideta

            d2xdeta2 = d2xds2 * dsdeta**2 + dxds * d2sdeta2 + 
     &                 d2xdsdr * drdeta * dsdeta + 
     &                 d2xdsdr * drdeta * dsdeta +
     &                 d2xdr2 * drdeta**2 + dxdr * d2rdeta2

            d2ydxi2 = d2yds2 * dsdxi**2 + dyds * d2sdxi2 + 
     &                d2ydsdr * drdxi * dsdxi + 
     &                d2ydsdr * drdxi * dsdxi +
     &                d2ydr2 * drdxi**2 + dydr * d2rdxi2
            
            d2ydxideta = d2yds2 * dsdeta * dsdxi + dyds * d2sdxideta + 
     &                   d2ydsdr * drdeta * dsdxi + 
     &                   d2ydsdr * drdxi * dsdeta +
     &                   d2ydr2 * drdeta * drdxi + dydr * d2rdxideta

            d2ydeta2 = d2yds2 * dsdeta**2 + dyds * d2sdeta2 + 
     &                 d2ydsdr * drdeta * dsdeta + 
     &                 d2ydsdr * drdeta * dsdeta +
     &                 d2ydr2 * drdeta**2 + dydr * d2rdeta2

            drjacdxi = -(d2xdxi2*dydeta + dxdxi*d2ydxideta - 
     &                   d2xdxideta*dydxi -
     &                   dxdeta*d2ydxi2) * rjac(i,j)**2
     
            drjacdeta = -(d2xdxideta*dydeta + dxdxi*d2ydeta2 - 
     &                    d2xdeta2*dydxi -
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
c.... write out the nose element spacing
c
        dr = x(1,1) - x(1,2)
        yl = y(2,1)
        ds = yl * 0.5 * sqrt( 1.0 + yl**2) + 
     &       0.5 * log( abs( yl + sqrt(1.0 + yl**2)))

        write(*,10) ds, dr, (x(1,2)-x(1,3))/dr
 10     format('ds_min = ',1pe20.13,' dr_min = ',1pe20.13,
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
        write(*,20) ds, dr, (x(1,ny-2)-x(1,ny-1))/dr
 20     format('ds_max = ',1pe20.13,' dr_max = ',1pe20.13,
     &         ' dr(ny)/dr(ny-1) = ',1pe20.13)
c
c.... write out the wall spacing
c
        if (debug) then
        do i = 2, nx
          yl = y(i,1)
          sl = yl * 0.5 * sqrt( 1.0 + yl**2) + 
     &         0.5 * log( abs( yl + sqrt(1.0 + yl**2)))
          yl = y(i-1,1)
          ds = sl - (yl * 0.5 * sqrt( 1.0 + yl**2) + 
     &         0.5 * log( abs( yl + sqrt(1.0 + yl**2))))
          write(70,"(3(1pe20.13,1x))") sl, ds
        end do
        close(70)
        end if
c
c.... write out the body.dat file
c
        open(9,file='body.dat',status='unknown')
        j = 1
        do i = 1, nx
          yl = y(i,j)
          sl = yl * 0.5 * sqrt( 1.0 + yl**2) + 
     &         0.5 * log( abs( yl + sqrt(1.0 + yl**2)))
          bn1 = rn1(i,j) / sqrt( rn1(i,j)**2 + rn2(i,j)**2 )
          bn2 = rn2(i,j) / sqrt( rn1(i,j)**2 + rn2(i,j)**2 )
          write(9,"(1(i5,1x),5(1pe20.13,1x))") i, sl, x(i,j), yl, 
     &                                         bn1, bn2
        end do
        close(9)
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
        else
c
c.... write out unformatted if on a CRAY
c
#ifdef USE_SNGL
          open(unit=10,file='grid.dat',form='unformatted')
          write(10) nx, ny, nz
          write(10) (((sngl(x(i,j)), i=1,nx), j=1,ny), k = 1, nz),
     &              (((sngl(y(i,j)), i=1,nx), j=1,ny), k = 1, nz),
     &              (((sngl(zero),   i=1,nx), j=1,ny), k = 1, nz)
          close(10)
#else
          open(unit=10,file='grid.dat',form='unformatted')
          write(10) nx, ny, nz
          write(10) (((x(i,j), i=1,nx), j=1,ny), k = 1, nz),
     &              (((y(i,j), i=1,nx), j=1,ny), k = 1, nz),
     &              (((zero,   i=1,nx), j=1,ny), k = 1, nz)
          close(10)
#endif
        
          if (debug) then

          open(unit=10,file='m1.dat',form='unformatted')
          write(10) nx, ny, nz
          write(10) 0.0, 0.0, 0.0, 0.0
          write(10) (((rm1(i,j), i=1,nx), j=1,ny), k = 1, nz),
     &              (((rn1(i,j), i=1,nx), j=1,ny), k = 1, nz),
     &              (((rm2(i,j), i=1,nx), j=1,ny), k = 1, nz),
     &              (((rn2(i,j), i=1,nx), j=1,ny), k = 1, nz),
     &              (((rjac(i,j), i=1,nx), j=1,ny), k = 1, nz)
          close(10)

          open(unit=10,file='m11.dat',form='unformatted')
          write(10) nx, ny, nz
          write(10) 0.0, 0.0, 0.0, 0.0
          write(10) (((rm11(i,j), i=1,nx), j=1,ny), k = 1, nz),
     &              (((rm12(i,j), i=1,nx), j=1,ny), k = 1, nz),
     &              (((rm21(i,j), i=1,nx), j=1,ny), k = 1, nz),
     &              (((rm22(i,j), i=1,nx), j=1,ny), k = 1, nz),
     &              (((rm1(i,j)*rn1(i,j)+rm2(i,j)*rn2(i,j), 
     &              i=1,nx), j=1,ny), k = 1, nz)
          close(10)

          open(unit=10,file='n11.dat',form='unformatted')
          write(10) nx, ny, nz
          write(10) 0.0, 0.0, 0.0, 0.0
          write(10) (((rn11(i,j), i=1,nx), j=1,ny), k = 1, nz),
     &              (((rn12(i,j), i=1,nx), j=1,ny), k = 1, nz),
     &              (((rn21(i,j), i=1,nx), j=1,ny), k = 1, nz),
     &              (((rn22(i,j), i=1,nx), j=1,ny), k = 1, nz),
     &              (((rjac(i,j), i=1,nx), j=1,ny), k = 1, nz)
          close(10)
          
          end if  ! debug

          if (useIJ) then

          open (unit=10, file='metric.dat', form='unformatted', 
     &          status='unknown')
          write(10) (( rm1(i,j), i=1,nx), j=1,ny),
     &              (( rm2(i,j), i=1,nx), j=1,ny),
     &              (( rn1(i,j), i=1,nx), j=1,ny),
     &              (( rn2(i,j), i=1,nx), j=1,ny),
     &              ((rm11(i,j), i=1,nx), j=1,ny),
     &              ((rm12(i,j), i=1,nx), j=1,ny),
     &              ((rm22(i,j), i=1,nx), j=1,ny),
     &              ((rn11(i,j), i=1,nx), j=1,ny),
     &              ((rn12(i,j), i=1,nx), j=1,ny),
     &              ((rn22(i,j), i=1,nx), j=1,ny)
          close(10)

          else
c
c.... write out the metric file (Note the order of i and j are reversed)
c
          open (unit=10, file='metric.dat', form='unformatted', 
     &          status='unknown')
          write(10) (( rm1(i,j), j=1,ny), i=1,nx),
     &              (( rm2(i,j), j=1,ny), i=1,nx),
     &              (( rn1(i,j), j=1,ny), i=1,nx),
     &              (( rn2(i,j), j=1,ny), i=1,nx),
     &              ((rm11(i,j), j=1,ny), i=1,nx),
     &              ((rm12(i,j), j=1,ny), i=1,nx),
     &              ((rm22(i,j), j=1,ny), i=1,nx),
     &              ((rn11(i,j), j=1,ny), i=1,nx),
     &              ((rn12(i,j), j=1,ny), i=1,nx),
     &              ((rn22(i,j), j=1,ny), i=1,nx)
          close(10)
   
          endif

        end if
        
        stop
 30     format(7(1x,1pe20.13))
        end

c=============================================================================c
        function arc(xi)
c=============================================================================c
        implicit double precision (a-h,o-z)
c=============================================================================c
        arc = sqrt(2.0)*0.5*xi*sqrt(1.0+2.0*xi**2) +
     &        0.5*log(sqrt(2.0)*xi + sqrt(1.0+2.0*xi**2))
        return
        end

c=============================================================================c
        function calcxi(s)
c=============================================================================c
        implicit double precision (a-h,o-z)

        common /xistuff/ sloc

        external funcxi, zbrent
c=============================================================================c
        sloc = s

        calcxi = zbrent(funcxi,0.0d0,2000.0d0,1.0d-14)

        return
        end

c=============================================================================c
        function funcxi(xi)
c=============================================================================c
        implicit double precision (a-h,o-z)
        
        common /xistuff/ sloc
c=============================================================================c
        funcxi = sloc - arc(xi)

        return
        end 

c=============================================================================c
        function calcdd(d1, d2)
c=============================================================================c
        implicit double precision (a-h,o-z)
        parameter (zero=0.0d0, pt5=0.5d0, one=1.0d0, two=2.0d0)
        
        common /stuff/ rd1, rd2
        
        external func, rtflsp
c=============================================================================c
        rd1 = d1
        rd2 = d2
        
        dd1 = 1.0
        dd2 = 7.0
        
c       calcdd = rtflsp(func,dd1,dd2,1.0d-14)
        calcdd = zbrent(func,dd1,dd2,1.0d-14)
        
        return
        end

c=============================================================================c
        function func(x)
c=============================================================================c
        implicit double precision (a-h,o-z)
        parameter (zero=0.0d0, pt5=0.5d0, one=1.0d0, two=2.0d0)
        
        common /stuff/ rd1, rd2
c=============================================================================c
        func = sinh(x)/x - one / sqrt(rd1 * rd2)

        return
        end 

c=============================================================================c
        function calcs( ximin1, ximax1, dxmin1, ds1)
c
c  Calculate a stretching ratio to fill a domain using exponential
c  stretching.
c
c=============================================================================c
        implicit double precision (a-h,o-z)
        parameter (zero=0.0d0, pt5=0.5d0, one=1.0d0, two=2.0d0)
        
        common /stuffs/ ximin, ximax, dxmin, ds

        external funcs, rtflsp
c=============================================================================c
        ximin = ximin1
        ximax = ximax1
        dxmin = dxmin1
        ds = ds1
        
c       calcs = rtflsp(funcs, 1.01d0, 1.02d0, 1.0d-14)
        calcs = zbrent(funcs, 1.0001d0, 1.2d0, 1.0d-14)

        return
        end

c=============================================================================c
        function funcs(sx)
c=============================================================================c
        implicit double precision (a-h,o-z)
        parameter (zero=0.0d0, pt5=0.5d0, one=1.0d0, two=2.0d0)
        
        common /stuffs/ ximin, ximax, dxmin, ds
c=============================================================================c
        c2 = log( dxmin / ds )
        c1 = (log( sx*dxmin / ds) - c2) / ds
             
        funcs = ximax - ximin - one / c1 * ( exp(c1 + c2) - exp(c2) )
        
        return
        end 
