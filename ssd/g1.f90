!=============================================================================!
        subroutine g1(v, g1v, nx, dx, optx, xper)
!
!  Take the gradient in the streamwise direction
!  updated to fourth order accurate differencing
!
!=============================================================================!
        use diff
        implicit none
        
        integer nx, optx
        logical :: xper
        real v(nx), g1v(nx), dx

        real,parameter :: one = 1.0
        real dxinv
        real a, b, c, w
        real gx1, gx2, gx3, gx4, gx5, gx6
!=============================================================================!

        dxinv  = one / dx

!.... seven point stencil in x

        if (optx.eq.0) then
          c = 1.0 / 60.0
        else
          w = 2.0 * 3.1415926535897932385e+0 / 12.0
          c = (w/2.0 - 2.0*Sin(w)/3.0 + Sin(2.0*w)/12.0) / &
              (5.0*Sin(w) - 4.0*Sin(2.0*w) + Sin(3.0*w))
        end if
        
        a = (2.0 + 15.0 * c) / 3.0
        b = -(1.0 + 48.0 * c) / 12.0

        gx1 =  -c * dxinv
        gx2 =  -b * dxinv
        gx3 =  -a * dxinv
        gx4 =   a * dxinv
        gx5 =   b * dxinv
        gx6 =   c * dxinv

!=============================================================================!
!.... compute the gradient in x
!=============================================================================!

        if (xper) then

          g1v(1)        = ( gx1 * v(nx-3)       + &
                                gx2 * v(nx-2)   + &
                                gx3 * v(nx-1)   + &
                                gx4 * v(2)      + &
                                gx5 * v(3)      + &
                                gx6 * v(4)      ) 
  
          g1v(2)        = ( gx1 * v(nx-2)       + &
                                gx2 * v(nx-1)   + &
                                gx3 * v(1)      + &
                                gx4 * v(3)      + &
                                gx5 * v(4)      + &
                                gx6 * v(5)      ) 
  
          g1v(3)        = ( gx1 * v(nx-1)       + &
                                gx2 * v(1)      + &
                                gx3 * v(2)      + &
                                gx4 * v(4)      + &
                                gx5 * v(5)      + &
                                gx6 * v(6)      ) 
  
          g1v(nx-2)     = ( gx1 * v(nx-5)       + &
                                gx2 * v(nx-4)   + &
                                gx3 * v(nx-3)   + &
                                gx4 * v(nx-1)   + &
                                gx5 * v(1)      + &
                                gx6 * v(2)      ) 
  
          g1v(nx-1)     = ( gx1 * v(nx-4)       + &
                                gx2 * v(nx-3)   + &
                                gx3 * v(nx-2)   + &
                                gx4 * v(1)      + &
                                gx5 * v(2)      + &
                                gx6 * v(3)      ) 
  
          g1v(nx) = g1v(1)
          
        else
        
          g1v(1)      = ( gc1 * v(1)  + &
                              gc2 * v(2)  + &
                              gc3 * v(3)  + &
                              gc4 * v(4)  + &
                              gc5 * v(5)  ) * dxinv
  
          g1v(2)      = ( gb1 * v(1)  + &
                              gb2 * v(2)  + &
                              gb3 * v(3)  + &
                              gb4 * v(4)  + &
                              gb5 * v(5)  ) * dxinv
  
          g1v(3)      = ( ga1 * v(1)  + &
                              ga2 * v(2)  + &
                              ga3 * v(4)  + &
                              ga4 * v(5)  ) * dxinv
          
          g1v(nx-2)   = ( ga1 * v(nx-4)  + &
                              ga2 * v(nx-3)  + &
                              ga3 * v(nx-1)  + &
                              ga4 * v(nx  )  ) * dxinv
  
          g1v(nx-1)  = -( gb1 * v(nx  )  + &
                              gb2 * v(nx-1)  + &
                              gb3 * v(nx-2)  + &
                              gb4 * v(nx-3)  + &
                              gb5 * v(nx-4)  ) * dxinv
  
          g1v(nx)    = -( gc1 * v(nx  )  + &
                              gc2 * v(nx-1)  + &
                              gc3 * v(nx-2)  + &
                              gc4 * v(nx-3)  + &
                              gc5 * v(nx-4)  ) * dxinv

        end if

!.... interior

        g1v(4:nx-3)  = ( gx1 * v(1:nx-6)        + &
                             gx2 * v(2:nx-5)    + &
                             gx3 * v(3:nx-4)    + &
                             gx4 * v(5:nx-2)    + &
                             gx5 * v(6:nx-1)    + &
                             gx6 * v(7:nx  )    ) 

        return
        end
