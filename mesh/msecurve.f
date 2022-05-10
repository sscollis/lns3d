c=============================================================================c
        program msecurve
c=============================================================================c
c
c  WARNING:  Only partially coded .... not useful yet!
c
c  Purpose:  Compute the "optimal" streamwise point distribution based on the
c            curvature of the MSE 
c
c  Author:   S. Scott Collis
c
c  Revised:  7-15-96
c
c=============================================================================c
        implicit double precision (a-h,o-z)

        parameter (mx=5000)
        dimension eta(mx), s(mx), ds1(mx), dss1(mx)
        common /NR_RUNGE_PATH/ eta, s

        common /xloc/ s1
        
        common /stuff/ AR, rm, rn, xmin, c1, c2, dsmin, xmax, smax, nx

        external derivs, func, rtflsp, func2, curv

        parameter (zero=0.0d0,    pt5=0.5d0,  one=1.0d0, 
     &             onept5=1.5d0,  two=2.0d0,  twopt5=2.5d0, 
     &             three=3.0d0)

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
        AR = 6.0
        rm = 4.0
        rn = 2.0
        xmin = 0.0
        
        s(1) = 0.0
        eta1 = 0.0
        eta2 = 1.0

        dsmin = 0.05
        dsmax = 2.0
        
        write(*,"('Enter nx1, dsmin, Lx ==> ',$)")
        read(*,*) nx, dsmin, xmax
        if (nx .gt. mx) then
          write(*,*) 'ERROR: nx > mx in pcurve'
          call exit(1)
        endif
        smax  = sloc(xmax)
        write(*,*) 'smax = ',smax
        
!       dx = xmax / real(nx-1)
!       do i = 1, nx
!         x = real(i-1)*dx
!         write(*,*) x, curv(x)
!       end do
!       
!       stop
        
        
        dsmax = rtflsp( func2, 0.1d0, 5.0d0, 1.0d-14 )
        write(*,*) 'dsmax = ',dsmax

        c2 = (1.0 / dsmax) / real(nx-1)
        c1 = (1.0 / dsmin) / real(nx-1) - c2
        
        write(*,"('Enter nx2 ==> ',$)")
        read(*,*) nx
        if (nx .gt. mx) then
          write(*,*) 'ERROR: nx > mx in pcurve'
          call exit(1)
        endif
        
        call NR_RUNGE( s(1), 1, eta1, eta2, nx-1, derivs )
c
c.... compute the derivative of the mapping function and spline it
c
c       do i = 1, nx
c         if (i.eq.1) then
c           dsl = (s(i+1)-s(i))/(eta(i+1)-eta(i))
c           dssl = (ds1(i+1)-ds1(i))/(eta(i+1)-eta(i))
c         else if (i.eq.nx) then
c           dsl = (s(i)-s(i-1))/(eta(i)-eta(i-1))
c           dssl = (ds1(i)-ds1(i-1))/(eta(i)-eta(i-1))
c         else
c           dsl = (s(i+1)-s(i-1))/(eta(i+1)-eta(i-1))
c           dssl = (ds1(i+1)-ds1(i-1))/(eta(i+1)-eta(i-1))
c         end if
c       end do
        
        ds = one / real(nx-1)
          
        ds1(1) = ( gc1 * s(1) + gc2 * s(2) + gc3 * s(3) +
     &             gc4 * s(4) + gc5 * s(5) ) / ds
        ds1(2) = ( gb1 * s(1) + gb2 * s(2) + gb3 * s(3) +
     &             gb4 * s(4) + gb5 * s(5) ) / ds
        do i = 3, nx-2
           ds1(i) = ( ga1 * s(i-2) + ga2 * s(i-1) + 
     &                ga3 * s(i+1) + ga4 * s(i+2) ) / ds
        end do
        ds1(nx-1) = -( gb1 * s(nx) + gb2 * s(nx-1) + 
     &                 gb3 * s(nx-2) +
     &                 gb4 * s(nx-3) + gb5 * s(nx-4) ) / ds
        ds1(nx) = -( gc1 * s(nx) + gc2 * s(nx-1) + 
     &               gc3 * s(nx-2) +
     &               gc4 * s(nx-3) + gc5 * s(nx-4) ) / ds
c       do i = 1, nx
c          ds1(i) = 1.0 / ds1(i)
c       end do
c
c.... compute the second derivative of the mapping function
c
        dss1(1) = ( gc1 * ds1(1) + gc2 * ds1(2) + gc3 * ds1(3) +
     &              gc4 * ds1(4) + gc5 * ds1(5) ) / ds
        dss1(2) = ( gb1 * ds1(1) + gb2 * ds1(2) + gb3 * ds1(3) +
     &              gb4 * ds1(4) + gb5 * ds1(5) ) / ds
        do i = 3, nx-2
           dss1(i) = ( ga1 * ds1(i-2) + ga2 * ds1(i-1) + 
     &                 ga3 * ds1(i+1) + ga4 * ds1(i+2) ) / ds
        end do
        dss1(nx-1) = -( gb1 * ds1(nx) + gb2 * ds1(nx-1) + 
     &                  gb3 * ds1(nx-2) + gb4 * ds1(nx-3) + 
     &                  gb5 * ds1(nx-4) ) / ds
        dss1(nx) = -( gc1 * ds1(nx) + gc2 * ds1(nx-1) + 
     &                gc3 * ds1(nx-2) +
     &                gc4 * ds1(nx-3) + gc5 * ds1(nx-4) ) / ds
c       do i = 1, nx
c          dss1(i) = dss1(i) * ds1(i)
c       end do

        open(10,file='tangent.map',status='unknown')
        do i = 1, nx
          s1 = s(i)
          x = rtflsp( func, 0.0d0, 1000.0d0, 1.0d-14)
          write(10,10) eta(i), s(i), ds1(i), dss1(i)
 10       format(8(1pe13.6,1x))
        end do
        close(10)
        
        stop
        end
c=============================================================================c
        function func2(dsmax)
c=============================================================================c
        implicit double precision (a-h,o-z)

        parameter (mx=5000)
        dimension eta(mx), s(mx)
        common /NR_RUNGE_PATH/ eta, s

        common /stuff/ AR, rm, rn, xmin, c1, c2, dsmin, xmax, smax, nx
        
        external derivs
c=============================================================================c
        s(1) = 0.0
        eta1 = 0.0
        eta2 = 1.0

        c2 = (1.0 / dsmax) / real(nx-1)
        c1 = (1.0 / dsmin) / real(nx-1) - c2
        
        call NR_RUNGE( s(1), 1, eta1, eta2, nx-1, derivs )
        
        func2 = smax - s(nx)
         
        write(*,*) dsmax, func2
        
        return
        end

        function curv( x )
        implicit double precision (a-h,o-z)
        common /stuff/ AR, rm, rn, xmin, c1, c2, dsmin, xmax, smax, nx
        
        dydx = rm/(rn*AR) * (1.0 - x/AR)**(rm-1.0)*
     &         (1.0-(1.0-x/AR)**rm)**((1.0-rn)/(rn))
        d2ydx2 = 1.0/AR**2 * rm**2*(1.0-rn)/rn**2 * (1.0-x/AR)**(2.0*(rm-1)) *
     &           (1.0-(1.0-x/AR)**rm)**((1.0-2.0*rn)/rn) -
     &           (rm*(rm-1.0)/(AR**2*rn)*(1.0-x/AR)**(rm-2)*
     &           (1.0-(1.0-x/AR)**rm)**((1.0-rn)/rn))
        curv = abs(d2ydx2) / (1.0 + dydx**2)**(3.0/2.0)

        if (x .eq. 0.0) then
          curv = 3.0
        else if ( x .ge. AR ) then
          curv = 0.0
        end if
        
        return
        end

c=============================================================================c
        subroutine derivs(eta, s, ds)
c=============================================================================c
        implicit double precision (a-h,o-z)
        external rtflsp, func, curv
        
        common /xloc/ s1
        common /stuff/ AR, rm, rn, xmin, c1, c2, dsmin, xmax, smax, nx
c=============================================================================c
        s1 = s
        x  = rtflsp( func, 0.0d0, s1, 1.0d-14)
        write(*,*) 'derivs ', x
        rk = curv(x)

        ds = 1.0 / ( c1 * abs(rk)/3.0 + c2 )
        
        return
        end
c=============================================================================c
        function func(x)
c=============================================================================c
        implicit double precision (a-h,o-z)
        common /xloc/ s
c=============================================================================c
        func = s - sloc(x)
        
        return
        end
c=============================================================================c
        function sloc(x)
c=============================================================================c
        implicit double precision (a-h,o-z)
c------------------------------------------------------------------------------
        sloc = arc(0.0d0,x)
        
        return
        end

c=============================================================================c
        function arc(x1,x2)
c=============================================================================c
        implicit double precision (a-h,o-z)
        parameter (zero=0.0d0, pt5=0.5d0, one=1.0d0, two=2.0d0)

        external derivs1, derivs2, RKQCR

        common /stuff/ AR, rm, rn, xmin, c1, c2, dsmin, xmax, smax, nx
c=============================================================================c
        s = zero
        xtmp = AR * pt5

c       write(*,"(8(1pe13.6,1x))") x1, x2, AR, rm, rn
        
        if (x1 .eq. x2 ) then
           arc = zero
           return
        end if

        if (x1 .le. xtmp) then
           y1 = ( one - ( one - x1/AR )**rm )**(one/rn)
           if (x2 .le. xtmp) then
              y2 = ( one - ( one - x2/AR )**rm )**(one/rn)
           else
              y2 = ( one - ( one - xtmp/AR )**rm )**(one/rn)
           end if

           call ODEINTR(s,1,y1,y2,1.0d-13,(y2-y1)*pt5,
     &                  1.d-14,nok,nbad,derivs2,RKQCR)

           if (x2 .gt. xtmp) then
              call ODEINTR(s,1,xtmp,x2,1.0d-13,(x2-x1)*pt5,
     &                     1.d-14,nok,nbad,derivs1,RKQCR)
           end if

        else

           call ODEINTR(s,1,x1,x2,1.0d-13,(x2-x1)*pt5,
     &                  1.d-14,nok,nbad,derivs1,RKQCR)

        end if

        arc = s

        return
        end
c=============================================================================c
        subroutine derivs1(n,s,x,ds)
c=============================================================================c
        implicit double precision (a-h,o-z)
        parameter (zero=0.0d0, pt5=0.5d0, one=1.0d0, two=2.0d0)

        common /stuff/ AR, rm, rn, xmin, c1, c2, dsmin, xmax, smax, nx
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
        implicit double precision (a-h,o-z)
        parameter (zero=0.0d0, pt5=0.5d0, one=1.0d0, two=2.0d0)

        common /stuff/ AR, rm, rn, xmin, c1, c2, dsmin, xmax, smax, nx
c=============================================================================c
        x = AR * ( one - (one - y**rn)**(one/rm) )

        dxdy = AR * rn/rm * ( one - x/AR )**(one-rm) *
     &         (one - (one - x/AR)**rm)**((rn-one)/rn)

        ds = sqrt( dxdy**2 + one )

        return
        end
