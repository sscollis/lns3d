c-----------------------------------------------------------------------
        function crvdist(npts, bsx, bsy, s1, s2, 
     &                   nknot, xknot, yknot, korder, ncoef)
c-----------------------------------------------------------------------
        implicit real (a-h,o-z)
        dimension bsx(ncoef), bsy(ncoef)
        dimension xknot(nknot), yknot(nknot)
        dimension ystart(1)

        external derivs, RKQCR
c-----------------------------------------------------------------------
        crvdist = 0.0
        eps = 1.0e-10

        if (s1 .eq. s2) then
           crvdist = 0.0
           return
        end if

        ystart(1) = crvdist
        call ODEINTR(ystart, 1, s1, s2, eps, (s2-s1)*0.5, eps, nok,
     &               nbad, derivs, RKQCR)
        crvdist = ystart(1)

c       numquad = korder * 10
c       ds = (s2 - s1) / real(numquad)
c       
c       x1 = BSDER( 0, s1, korder, xknot, ncoef, bsx )
c       y1 = BSDER( 0, s1, korder, yknot, ncoef, bsy )
c
c       do i = 1, numquad
c         sint = s1 + (i * ds)
c         x2 = BSDER( 0, sint, korder, xknot, ncoef, bsx )
c         y2 = BSDER( 0, sint, korder, yknot, ncoef, bsy )
c         crvdist = crvdist + sqrt( (x2 - x1)**2 + (y2 - y1)**2 )
c         x1 = x2
c         y1 = y2
c       end do
        
        return
        end

c-----------------------------------------------------------------------
        subroutine derivs(nl, arcl, sl, darcl)
c-----------------------------------------------------------------------
        parameter (mpts=1024)
        dimension s(mpts), t(mpts,2)
        parameter (kmax = 8)
        dimension xknot(kmax+mpts), yknot(kmax+mpts)
        common /bspline/  npts, t, darc, s1, s2, s, nknot, xknot, yknot, 
     &                    korder, ncoef

        external bsder
c-----------------------------------------------------------------------
        dxds = BSDER( 1, sl, korder, xknot, ncoef, t(1,1) )
        dyds = BSDER( 1, sl, korder, yknot, ncoef, t(1,2) )

        darcl = sqrt( dxds**2 + dyds**2 )

        return
        end

c-----------------------------------------------------------------------
        function crvarc(npts, s, bsx, bsy, nknot, xknot, yknot, 
     &                  korder, ncoef, arcl)
c-----------------------------------------------------------------------
        implicit real (a-h,o-z)
        dimension s(npts), bsx(ncoef), bsy(ncoef), arcl(npts)
        dimension xknot(nknot), yknot(nknot)
c-----------------------------------------------------------------------
        numquad = korder * 10
        crvarc = 0.0
        arcl(1) = 0.0

        do i = 2, npts
          arcl(i) = arcl(i-1)
          s1 = s(i-1)
          s2 = s(i)
          ds = (s2 - s1) / real(numquad)

          arcl(i) = arcl(i-1) + crvdist(npts, bsx, bsy, s1, s2, 
     &              nknot, xknot, yknot, korder, ncoef)

c         x1 = BSDER( 0, s1, korder, xknot, ncoef, bsx )
c         y1 = BSDER( 0, s1, korder, yknot, ncoef, bsy )
c         do iquad = 1, numquad
c           sint = s1 + real(iquad)*ds
c           x2 = BSDER( 0, sint, korder, xknot, ncoef, bsx )
c           y2 = BSDER( 0, sint, korder, yknot, ncoef, bsy )
c           arcl(i) = arcl(i) + sqrt( (x2 - x1)**2 + (y2 - y1)**2 )
c           x1 = x2
c           y1 = y2
c         end do

        end do

        crvarc = arcl(npts)

        return
        end

c---------------------------------------------------------------------
        function gets()
c---------------------------------------------------------------------
        implicit real (a-h,o-z)
        
        parameter (mpts=1024)
        dimension s(mpts), t(mpts,2)
        parameter (kmax = 8)
        dimension xknot(kmax+mpts), yknot(kmax+mpts)
        common /bspline/ npts, t, darc, s1, s2, s, nknot, xknot, yknot, 
     &                   korder, ncoef
        
        external rtsec, func
c---------------------------------------------------------------------

        gets = rtsec(func,s1,s2,1.0e-10)

        return
        end

c---------------------------------------------------------------------
        function func(sint)
c---------------------------------------------------------------------
        implicit real (a-h,o-z)
        
        parameter (mpts=1024)
        dimension s(mpts), t(mpts,2)
        parameter (kmax = 8)
        dimension xknot(kmax+mpts), yknot(kmax+mpts)
        common /bspline/ npts, t, darc, s1, s2, s, nknot, xknot, yknot, 
     &                   korder, ncoef

        external crvdist
c---------------------------------------------------------------------
        
        func = darc - crvdist(npts, t(1,1), t(1,2), s1, sint, 
     &                        nknot, xknot, yknot, korder, ncoef)

        return
        end
