!=============================================================================!
        subroutine g2(v, g2v, ny, dy, opty, yper)
!
!  Take the gradient in the wall normal direction
!  updated to fourth order accurate differencing
!
!=============================================================================!
        use stencil
        implicit none

        integer ny, opty
        logical :: yper
        real v(ny), g2v(ny), dy

        real,parameter :: one = 1.0
        real dyinv
        real a, b, c, w
        real gy1, gy2, gy3, gy4, gy5, gy6
!=============================================================================!

        dyinv  = one / dy

!.... seven point stencil in y

        if (opty.eq.0) then
          c = 1.0 / 60.0
        else
          w = 2.0 * 3.1415926535897932385e+0 / 12.0
          c = (w/2.0 - 2.0*Sin(w)/3.0 + Sin(2.0*w)/12.0) / &
              (5.0*Sin(w) - 4.0*Sin(2.0*w) + Sin(3.0*w))
        end if
        
        a = (2.0 + 15.0 * c) / 3.0
        b = -(1.0 + 48.0 * c) / 12.0

        gy1 =  -c * dyinv
        gy2 =  -b * dyinv
        gy3 =  -a * dyinv
        gy4 =   a * dyinv
        gy5 =   b * dyinv
        gy6 =   c * dyinv

!=============================================================================!
!.... compute the gradient in y
!=============================================================================!

        if (yper) then
        
          g2v(1)       = ( gy1 * v(ny-3)        + &
                                gy2 * v(ny-2)   + &
                                gy3 * v(ny-1)   + &
                                gy4 * v(2)      + &
                                gy5 * v(3)      + &
                                gy6 * v(4)      ) 
  
          g2v(2)       = ( gy1 * v(ny-2)        + &
                                gy2 * v(ny-1)   + &
                                gy3 * v(1)      + &
                                gy4 * v(3)      + &
                                gy5 * v(4)      + &
                                gy6 * v(5)      ) 
  
          g2v(3)       = ( gy1 * v(ny-1)        + &
                                gy2 * v(1)      + &
                                gy3 * v(2)      + &
                                gy4 * v(4)      + &
                                gy5 * v(5)      + &
                                gy6 * v(6)      ) 
  
          g2v(ny-2)    = ( gy1 * v(ny-5)        + &
                                gy2 * v(ny-4)   + &
                                gy3 * v(ny-3)   + &
                                gy4 * v(ny-1)   + &
                                gy5 * v(1)      + &
                                gy6 * v(2)      ) 
  
          g2v(ny-1)    = ( gy1 * v(ny-4)        + &
                                gy2 * v(ny-3)   + &
                                gy3 * v(ny-2)   + &
                                gy4 * v(1)      + &
                                gy5 * v(2)      + &
                                gy6 * v(3)      ) 
  
          g2v(ny) = g2v(1)

        else

          g2v(1)     = ( gc1 * v(1)  + &
                              gc2 * v(2)  + &
                              gc3 * v(3)  + &
                              gc4 * v(4)  + &
                              gc5 * v(5)  ) * dyinv
  
          g2v(2)     = ( gb1 * v(1)  + &
                              gb2 * v(2)  + &
                              gb3 * v(3)  + &
                              gb4 * v(4)  + &
                              gb5 * v(5)  ) * dyinv
  
          g2v(3)     = ( ga1 * v(1)  + &
                              ga2 * v(2)  + &
                              ga3 * v(4)  + &
                              ga4 * v(5)  ) * dyinv
  
          g2v(ny-2)  = ( ga1 * v(ny-4)  + &
                              ga2 * v(ny-3)  + &
                              ga3 * v(ny-1)  + &
                              ga4 * v(ny  )  ) * dyinv
  
          g2v(ny-1) = -( gb1 * v(ny  )  + &
                              gb2 * v(ny-1)  + &
                              gb3 * v(ny-2)  + &
                              gb4 * v(ny-3)  + &
                              gb5 * v(ny-4)  ) * dyinv
  
          g2v(ny)   = -( gc1 * v(ny  )  + &
                              gc2 * v(ny-1)  + &
                              gc3 * v(ny-2)  + &
                              gc4 * v(ny-3)  + &
                              gc5 * v(ny-4)  ) * dyinv
                            
        end if
        
!.... interior

        g2v(4:ny-3) = ( gy1 * v(1:ny-6) + &
                             gy2 * v(2:ny-5)    + &
                             gy3 * v(3:ny-4)    + &
                             gy4 * v(5:ny-2)    + &
                             gy5 * v(6:ny-1)    + &
                             gy6 * v(7:ny  )    ) 

        return
        end
