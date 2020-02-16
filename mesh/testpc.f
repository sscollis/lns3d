c=============================================================================c
        program testpc
c=============================================================================c
        implicit double precision (a-h,o-z)

        logical IRIS, DIAG
        parameter (IRIS = .false.)
        
        parameter (mx=512,my=512)

        parameter (zero=0.0d0,    pt5=0.5d0,  one=1.0d0, 
     &             onept5=1.5d0,  two=2.0d0,  twopt5=2.5d0, 
     &             three=3.0d0)
        
        dimension x(mx,my), y(mx,my), xi(mx), eta(my)
        dimension ss(mx), dss(mx), d2ss(mx)
        dimension rr(my), drr(my), d2rr(my)
        dimension rm1(mx,my),  rm2(mx,my)
        dimension rn1(mx,my),  rn2(mx,my)
        dimension rjac(mx,my)
        dimension rm11(mx,my), rm12(mx,my), rm21(mx,my), rm22(mx,my)
        dimension rn11(mx,my), rn12(mx,my), rn21(mx,my), rn22(mx,my)
        dimension g11x(mx,my), g11y(mx,my), g12x(mx,my), g12y(mx,my)
        dimension g22x(mx,my), g22y(mx,my)

        common /stuff/ rd1, rd2
        
        external calcdd
c=============================================================================c
        write(*,"('Enter (nx,ny) ==> ',$)") 
        read(*,*) nx, ny
        write(*,"('Enter (Xmax, Ymax) ==> ',$)")
        read(*,*) rLx, rLy
c
c.... uniform grid in computational space
c       
        dxi = sqrt(rLx) / dble(nx-1)
        do i = 1, nx
          xi(i) = dble(i-1)*dxi
          dss(i) = sqrt(rLx)
          d2ss(i) = zero
        end do

        etamax = sqrt(one + two * rLy)
        etamin = one
        deta = (etamax - etamin) / dble(ny-1)
        do j = 1, ny
          eta(j) = one + dble(j-1)*deta
          drr(j) = (etamax - etamin)
          d2rr(j) = zero
        end do
c
c.... Hyperbolic tangent stretching in wall normal direction
c
        if (.false.) then
        rmax = rLy
        
        rd1 = 0.0045d0
        rd2 = 5.0d0
        dd = calcdd()

        write(*,*) 'dd = ',dd
        
        deta = 1.0 / dble(ny-1)
        
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
        write(*,*) 'Rmax = ',rr(ny)
        do j = 1, ny
          eta(j) = sqrt( 1.0 + 2.0 * rr(j) )
        end do
        end if
c
c.... exponential map in the wall normal direction
c
        if (.false.) then
        deta2 = 1.0 / dble(ny-1)
        dymin = 1.0d-2
        sy = 1.02d0
        c2 = log( dymin / deta2 )
        c1 = (log( sy*dymin / deta2) - c2) / deta2
        do j = 1, ny
          eta(j) = dble(j-1) * deta2
          rr(j)  = 1.0 / c1 * ( exp(c1 * eta(j) + c2) - exp(c2) ) 
          drr(j)  = exp(c1 * eta(j) + c2)
          d2rr(j) = c1 * exp(c1 * eta(j) + c2)
          write(20,30) eta(j), rr(j), drr(j), d2rr(j), d2rr(j)/drr(j)
        end do
        write(*,*) 'Rmax = ',rr(ny)
        do j = 1, ny
          eta(j) = sqrt( 1.0 + 2.0 * rr(j) )
        end do
        end if
                
        do i = 1, nx
          do j = 1, ny
            x(i,j) = xi(i)**2 + (1.0 - eta(j)**2)/2.0
            y(i,j) = sqrt(2.0) * xi(i) * eta(j)

            dsdxi  = dss(i)
            drdxi  = zero

            dsdeta = zero
            drdeta = drr(j)

            dxds  = 2.0 * xi(i)
            dxdr  = -eta(j)
            
            dyds  = sqrt(2.0) * eta(j)
            dydr  = sqrt(2.0) * xi(i)
            
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

        write(*,*) ds, dr, (x(1,2)-x(1,3))/dr
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
          open(unit=10,file='grid.dat',form='unformatted')
          write(10) nx, ny, nz
          write(10) (((sngl(x(i,j)), i=1,nx), j=1,ny), k = 1, nz),
     &              (((sngl(y(i,j)), i=1,nx), j=1,ny), k = 1, nz),
     &              (((sngl(zero),   i=1,nx), j=1,ny), k = 1, nz)
          close(10)
        
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
 30     format(7(1x,1pe20.13))
        end

c=============================================================================c
        function calcdd()
c=============================================================================c
        implicit double precision (a-h,o-z)
        parameter (zero=0.0d0, pt5=0.5d0, one=1.0d0, two=2.0d0)
        
        common /stuff/ rd1, rd2
        
        external func
c=============================================================================c
        dd1 = 3.0
        dd2 = 7.0
        
        calcdd = rtflsp(func,dd1,dd2,1.0d-14)
        
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
