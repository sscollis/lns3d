c=============================================================================c
        program pcurve
c=============================================================================c
c
c  Purpose:  Compute the "optimal" streamwise point distribution based on the
c            curvature of the parabolic cylinder.
c
c  Author:   S. Scott Collis
c
c  Revised:  7-15-96
c
c=============================================================================c
        implicit double precision (a-h,o-z)

        parameter (mx=5000)
        dimension eta(mx), s(mx), ds1(mx), dss1(mx)
        common /path/ eta, s

        common /xloc/ s1
        
        common /stuff/ c1, c2, dsmin, xmax, smax, nx

        external derivs, func, rtflsp, func2

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
        dsmax = rtflsp( func2, 0.1d0, 5.0d0, 1.0d-14 )
        write(*,*) 'smax = ',smax,' dsmax = ',dsmax

        c2 = (1.0 / dsmax) / real(nx-1)
        c1 = (1.0 / dsmin) / real(nx-1) - c2
        
        write(*,"('Enter nx2 ==> ',$)")
        read(*,*) nx
        if (nx .gt. mx) then
          write(*,*) 'ERROR: nx > mx in pcurve'
          call exit(1)
        endif
        
        call RUNGE( s(1), 1, eta1, eta2, nx-1, derivs )
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
        common /path/ eta, s

        common /stuff/ c1, c2, dsmin, xmax, smax, nx
        
        external derivs
c=============================================================================c
        s(1) = 0.0
        eta1 = 0.0
        eta2 = 1.0

        c2 = (1.0 / dsmax) / real(nx-1)
        c1 = (1.0 / dsmin) / real(nx-1) - c2
        
        call RUNGE( s(1), 1, eta1, eta2, nx-1, derivs )
        
        func2 = smax - s(nx)
        
        return
        end
c=============================================================================c
        subroutine derivs(eta, s, ds)
c=============================================================================c
        implicit double precision (a-h,o-z)
        external rtflsp, func
        
        common /xloc/ s1
        common /stuff/ c1, c2, dsmin, xmax, smax, nx
c=============================================================================c
        s1 = s
        x  = rtflsp( func, 0.0d0, s1, 1.0d-14)
        
        rk = -(2.0 * x + 1.0)**(-1.5)

        ds = 1.0 / ( c1 * abs(rk) + c2 )
        
        return
        end
c=============================================================================c
        function func(x)
c=============================================================================c
        implicit double precision (a-h,o-z)
        common /xloc/ s
c=============================================================================c
        y = sqrt(2.0 * x)
        func = s - 0.5 * (y  * sqrt(1.0 + y**2) + 
     &                    log( abs(y + sqrt(1.0 + y**2)) ) )
        
        return
        end
c=============================================================================c
        function sloc(x)
c=============================================================================c
        implicit double precision (a-h,o-z)
c------------------------------------------------------------------------------
        y = sqrt(2.0 * x)
        sloc =  0.5 * (y  * sqrt(1.0 + y**2) + 
     &                 log( abs(y + sqrt(1.0 + y**2)) ) )
        
        return
        end
