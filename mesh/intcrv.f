        subroutine intcrv( npts, nsd, s, b, t, sint, r, dr, mpts )
        
        implicit real*8 (a-h,o-z)
        dimension s(mpts), b(mpts,nsd), t(mpts,nsd), r(nsd), dr(nsd)
     
        klo = 1
        khi = npts
 10     if ( khi-klo .gt. 1 ) then
          k = (khi+klo)/2
          if (s(k) .gt. sint) then
            khi=k
          else
            klo=k
          end if
          goto 10
        end if
        
        u  = (sint - s(klo))

        ub = (s(khi) - s(klo))

        do j = 1, nsd
          ba = b(klo,j)
          bb = b(khi,j)
          ta = t(klo,j)
          tb = t(khi,j)
          
          c1 = -(2.0*ub*(bb-ta*ub-ba) - ub**2*(tb-ta))/ub**4
          
          c2 = -(ub**3*(tb-ta) - 3.0*ub**2*(bb-ta*ub-ba))/ub**4

          r(j)  = c1 * u**3 + c2 * u**2 + u*ta + ba
          dr(j) = 3.0 * c1 * u**2 + 2.0 * c2 * u + ta
        end do
        
        return
        end
