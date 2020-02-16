        subroutine allocwrk
#ifdef USE_IMSL
        COMMON /WORKSP/  RWKSP
        REAL RWKSP(58388)
        CALL IWKIN(58388)
#endif
        return
        end
c=============================================================================c
        function calcs( ximax1, dxmin1, ds1)
c=============================================================================c
c       implicit double precision (a-h,o-z)
        parameter (zero=0.0d0, pt5=0.5d0, one=1.0d0, two=2.0d0)
        
        common /stuffs/ ximax, dxmin, ds

        external funcs, rtflsp
c=============================================================================c
        ximax = ximax1
        dxmin = dxmin1
        ds = ds1
        
c       calcs = rtflsp(funcs, 1.01d0, 1.02d0, 1.0d-14)
        calcs = zbrent(funcs, 1.0001e0, 1.2e0, 1.0e-14)

        return
        end
c=============================================================================c
        function funcs(sx)
c=============================================================================c
c       implicit double precision (a-h,o-z)
        parameter (zero=0.0d0, pt5=0.5d0, one=1.0d0, two=2.0d0)
        
        common /stuffs/ ximax, dxmin, ds
c=============================================================================c
        c2 = log( dxmin / ds )
        c1 = (log( sx*dxmin / ds) - c2) / ds
             
        funcs = ximax - one / c1 * ( exp(c1 + c2) - exp(c2) )
        
        return
        end 
c=============================================================================c
        function calcdd(d1,d2)
c=============================================================================c
c       implicit double precision (a-h,o-z)
        parameter (zero=0.0d0, pt5=0.5d0, one=1.0d0, two=2.0d0)
        
        common /stuff/ rd1, rd2
        
        external func
c=============================================================================c
        rd1 = d1
        rd2 = d2
        
        dd1 = 3.0
        dd2 = 7.0
        
c       calcdd = rtflsp(func,dd1,dd2,1.0d-14)
        calcdd = zbrent(func,dd1,dd2,1.0e-14)
        
        return
        end
c=============================================================================c
        function func(x)
c=============================================================================c
c       implicit double precision (a-h,o-z)
        parameter (zero=0.0d0, pt5=0.5d0, one=1.0d0, two=2.0d0)
        
        common /stuff/ rd1, rd2
c=============================================================================c
        func = sinh(x)/x - one / sqrt(rd1 * rd2)

        return
        end 
