!=============================================================================!
        subroutine cgrad( ndof, nx, ny, v, g1v, g2v, dx, dy, optx, opty, &
                          xper, yper, lsym, rsym, bsym, tsym, carp)
!
!  Take the gradient of a complex 2-D field.
!  Updated to sixth-order accurate differencing on the interior with
!  the option for optimized fourth-order differencing.
!
!  Revised: 4-17-96
!=============================================================================!
        use stencil
        implicit none
        
        integer :: ndof, nx, ny, optx, opty
        logical :: xper, yper
        logical :: lsym, rsym, bsym, tsym, carp
        complex :: v(ny,nx,ndof), g1v(ny,nx,ndof), g2v(ny,nx,ndof)
        real    :: dx, dy
        
        real, parameter :: zero = 0.0, one = 1.0, pt5 = 0.5
        real dxinv, dyinv
        real a, b, c, w
        real gx1, gx2, gx3, gx4, gx5, gx6
        real gy1, gy2, gy3, gy4, gy5, gy6
        
        integer :: i, j, idof
        real :: eps = 1.0e-12, isign
!=============================================================================!
        dxinv  = one / dx
        dyinv  = one / dy

!.... seven point stencil in x

        if (optx.eq.0) then
          c = 1.0 / 60.0
        else if (optx.eq.-1) then
          c = 0.0
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

!.... seven point stencil in y

        if (opty.eq.0) then
          c = 1.0 / 60.0
        else if (opty.eq.-1) then
          c = 0.0
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
!.... compute the gradient in x
!=============================================================================!

        if (xper) then

          g1v(:,1,:)        = ( gx1 * v(:,nx-3,:)       + &
                                gx2 * v(:,nx-2,:)       + &
                                gx3 * v(:,nx-1,:)       + &
                                gx4 * v(:,2,:)          + &
                                gx5 * v(:,3,:)          + &
                                gx6 * v(:,4,:)  ) 
  
          g1v(:,2,:)        = ( gx1 * v(:,nx-2,:)       + &
                                gx2 * v(:,nx-1,:)       + &
                                gx3 * v(:,1,:)          + &
                                gx4 * v(:,3,:)          + &
                                gx5 * v(:,4,:)          + &
                                gx6 * v(:,5,:)  ) 
  
          g1v(:,3,:)        = ( gx1 * v(:,nx-1,:)       + &
                                gx2 * v(:,1,:)          + &
                                gx3 * v(:,2,:)          + &
                                gx4 * v(:,4,:)          + &
                                gx5 * v(:,5,:)          + &
                                gx6 * v(:,6,:)  ) 
  
          g1v(:,nx-2,:)     = ( gx1 * v(:,nx-5,:)       + &
                                gx2 * v(:,nx-4,:)       + &
                                gx3 * v(:,nx-3,:)       + &
                                gx4 * v(:,nx-1,:)       + &
                                gx5 * v(:,1,:)          + &
                                gx6 * v(:,2,:)  ) 
  
          g1v(:,nx-1,:)     = ( gx1 * v(:,nx-4,:)       + &
                                gx2 * v(:,nx-3,:)       + &
                                gx3 * v(:,nx-2,:)       + &
                                gx4 * v(:,1,:)          + &
                                gx5 * v(:,2,:)          + &
                                gx6 * v(:,3,:)  ) 
  
          g1v(:,nx,:) = g1v(:,1,:)
          
        else
        
          g1v(:,1,:)      = ( gc1 * v(:,1,:)  + &
                              gc2 * v(:,2,:)  + &
                              gc3 * v(:,3,:)  + &
                              gc4 * v(:,4,:)  + &
                              gc5 * v(:,5,:)  ) * dxinv

          g1v(:,2,:)      = ( gb1 * v(:,1,:)  + &
                              gb2 * v(:,2,:)  + &
                              gb3 * v(:,3,:)  + &
                              gb4 * v(:,4,:)  + &
                              gb5 * v(:,5,:)  ) * dxinv

          g1v(:,3,:)      = ( ga1 * v(:,1,:)  + &
                              ga2 * v(:,2,:)  + &
                              ga3 * v(:,4,:)  + &
                              ga4 * v(:,5,:)  ) * dxinv
          
          g1v(:,nx-2,:)   = ( ga1 * v(:,nx-4,:)  + &
                              ga2 * v(:,nx-3,:)  + &
                              ga3 * v(:,nx-1,:)  + &
                              ga4 * v(:,nx  ,:)  ) * dxinv
  
          g1v(:,nx-1,:)  = -( gb1 * v(:,nx  ,:)  + &
                              gb2 * v(:,nx-1,:)  + &
                              gb3 * v(:,nx-2,:)  + &
                              gb4 * v(:,nx-3,:)  + &
                              gb5 * v(:,nx-4,:)  ) * dxinv
  
          g1v(:,nx,:)    = -( gc1 * v(:,nx  ,:)  + &
                              gc2 * v(:,nx-1,:)  + &
                              gc3 * v(:,nx-2,:)  + &
                              gc4 * v(:,nx-3,:)  + &
                              gc5 * v(:,nx-4,:)  ) * dxinv

        end if

!.... implement the symmetric conditions

        if (lsym) then
          do idof = 1, ndof
            if (idof .eq. 3) then
              isign = -one
            else
              isign = one
            end if
  
            g1v(:,1,idof)     = ( isign * gx1 * v(:,4,idof)     + &
                                  isign * gx2 * v(:,3,idof)     + &
                                  isign * gx3 * v(:,2,idof)     + &
                                  gx4 * v(:,2,idof)             + &
                                  gx5 * v(:,3,idof)             + &
                                  gx6 * v(:,4,idof)             ) 
    
            g1v(:,2,idof)     = ( isign * gx1 * v(:,3,idof)     + &
                                  isign * gx2 * v(:,2,idof)     + &
                                  gx3 * v(:,1,idof)             + &
                                  gx4 * v(:,3,idof)             + &
                                  gx5 * v(:,4,idof)             + &
                                  gx6 * v(:,5,idof)             ) 
    
            g1v(:,3,idof)     = ( isign * gx1 * v(:,2,idof)     + &
                                  gx2 * v(:,1,idof)             + &
                                  gx3 * v(:,2,idof)             + &
                                  gx4 * v(:,4,idof)             + &
                                  gx5 * v(:,5,idof)             + &
                                  gx6 * v(:,6,idof)             ) 
          end do
        end if
        
        if (rsym) then
          do idof = 1, ndof
            if (idof .eq. 3) then
              isign = -one
            else
              isign = one
            end if
  
            g1v(:,nx-2,idof)  = ( gx1 * v(:,nx-5,idof)          + &
                                  gx2 * v(:,nx-4,idof)          + &
                                  gx3 * v(:,nx-3,idof)          + &
                                  gx4 * v(:,nx-1,idof)          + &
                                  gx5 * v(:,nx,idof)            + &
                                  isign * gx6 * v(:,nx-1,idof)          ) 
    
            g1v(:,nx-1,idof)  = ( gx1 * v(:,nx-4,idof)          + &
                                  gx2 * v(:,nx-3,idof)          + &
                                  gx3 * v(:,nx-2,idof)          + &
                                  gx4 * v(:,nx,idof)            + &
                                  isign * gx5 * v(:,nx-1,idof)  + &
                                  isign * gx6 * v(:,nx-2,idof)  ) 
    
            g1v(:,nx,idof)    = ( gx1 * v(:,nx-3,idof)          + &
                                  gx2 * v(:,nx-2,idof)          + &
                                  gx3 * v(:,nx-1,idof)          + &
                                  isign * gx4 * v(:,nx-1,idof)  + &
                                  isign * gx5 * v(:,nx-2,idof)  + &
                                  isign * gx6 * v(:,nx-3,idof)  ) 
          end do
        end if

!.... interior

        g1v(:,4:nx-3,:)  = ( gx1 * v(:,1:nx-6,:)        + &
                             gx2 * v(:,2:nx-5,:)        + &
                             gx3 * v(:,3:nx-4,:)        + &
                             gx4 * v(:,5:nx-2,:)        + &
                             gx5 * v(:,6:nx-1,:)        + &
                             gx6 * v(:,7:nx  ,:)        ) 

!=============================================================================!
!.... compute the gradient in y
!=============================================================================!

        if (yper) then
        
          g2v(1,:,:)        = ( gy1 * v(ny-3,:,:)       + &
                                gy2 * v(ny-2,:,:)       + &
                                gy3 * v(ny-1,:,:)       + &
                                gy4 * v(2,:,:)          + &
                                gy5 * v(3,:,:)          + &
                                gy6 * v(4,:,:)          ) 
  
          g2v(2,:,:)        = ( gy1 * v(ny-2,:,:)       + &
                                gy2 * v(ny-1,:,:)       + &
                                gy3 * v(1,:,:)          + &
                                gy4 * v(3,:,:)          + &
                                gy5 * v(4,:,:)          + &
                                gy6 * v(5,:,:)  ) 
  
          g2v(3,:,:)        = ( gy1 * v(ny-1,:,:)       + &
                                gy2 * v(1,:,:)          + &
                                gy3 * v(2,:,:)          + &
                                gy4 * v(4,:,:)          + &
                                gy5 * v(5,:,:)          + &
                                gy6 * v(6,:,:)  ) 
  
          g2v(ny-2,:,:)     = ( gy1 * v(ny-5,:,:)       + &
                                gy2 * v(ny-4,:,:)       + &
                                gy3 * v(ny-3,:,:)       + &
                                gy4 * v(ny-1,:,:)       + &
                                gy5 * v(1,:,:)          + &
                                gy6 * v(2,:,:)  ) 
  
          g2v(ny-1,:,:)     = ( gy1 * v(ny-4,:,:)       + &
                                gy2 * v(ny-3,:,:)       + &
                                gy3 * v(ny-2,:,:)       + &
                                gy4 * v(1,:,:)          + &
                                gy5 * v(2,:,:)          + &
                                gy6 * v(3,:,:)  ) 
  
          g2v(ny,:,:) = g2v(1,:,:)

        else

          g2v(1,:,:)     = ( gc1 * v(1,:,:)     + &
                             gc2 * v(2,:,:)     + &
                             gc3 * v(3,:,:)     + &
                             gc4 * v(4,:,:)     + &
                             gc5 * v(5,:,:)  ) * dyinv

          g2v(2,:,:)     = ( gb1 * v(1,:,:)     + &
                             gb2 * v(2,:,:)     + &
                             gb3 * v(3,:,:)     + &
                             gb4 * v(4,:,:)     + &
                             gb5 * v(5,:,:)  ) * dyinv
  
          g2v(3,:,:)     = ( ga1 * v(1,:,:)     + &
                             ga2 * v(2,:,:)     + &
                             ga3 * v(4,:,:)     + &
                             ga4 * v(5,:,:)  ) * dyinv
  
          g2v(ny-2,:,:)  = ( ga1 * v(ny-4,:,:)  + &
                             ga2 * v(ny-3,:,:)  + &
                             ga3 * v(ny-1,:,:)  + &
                             ga4 * v(ny  ,:,:)  ) * dyinv
  
          g2v(ny-1,:,:) = -( gb1 * v(ny  ,:,:)  + &
                             gb2 * v(ny-1,:,:)  + &
                             gb3 * v(ny-2,:,:)   + &
                             gb4 * v(ny-3,:,:)  + &
                             gb5 * v(ny-4,:,:)  ) * dyinv
  
          g2v(ny,:,:)   = -( gc1 * v(ny  ,:,:)  + &
                             gc2 * v(ny-1,:,:)  + &
                             gc3 * v(ny-2,:,:)  + &
                             gc4 * v(ny-3,:,:)  + &
                             gc5 * v(ny-4,:,:)  ) * dyinv
                            
        end if
        
!.... interior

        g2v(4:ny-3,:,:) = ( gy1 * v(1:ny-6,:,:) + &
                            gy2 * v(2:ny-5,:,:) + &
                            gy3 * v(3:ny-4,:,:) + &
                            gy4 * v(5:ny-2,:,:) + &
                            gy5 * v(6:ny-1,:,:) + &
                            gy6 * v(7:ny  ,:,:) ) 

!.... Implement Carpenter's boundary stencil

        if (carp .and. (.not. yper) ) then
          g2v(1,:,:)     = ( gg1 * v(1,:,:)     + &
                             gg2 * v(2,:,:)     + &
                             gg3 * v(3,:,:)     + &
                             gg4 * v(4,:,:)     + &
                             gg5 * v(5,:,:)     + &
                             gg6 * v(6,:,:)  ) * dyinv
  
          g2v(2,:,:)     = ( gh1 * v(1,:,:)     + &
                             gh2 * v(2,:,:)     + &
                             gh3 * v(3,:,:)     + &
                             gh4 * v(4,:,:)     + &
                             gh5 * v(5,:,:)     + &
                             gh6 * v(6,:,:)  ) * dyinv

          g2v(3,:,:)     = ( gi1 * v(1,:,:)     + &
                             gi2 * v(2,:,:)     + &
                             gi3 * v(3,:,:)     + &
                             gi4 * v(4,:,:)     + &
                             gi5 * v(5,:,:)     + &
                             gi6 * v(6,:,:)  ) * dyinv

          g2v(4,:,:)     = ( gj1 * v(1,:,:)     + &
                             gj2 * v(2,:,:)     + &
                             gj3 * v(3,:,:)     + &
                             gj4 * v(4,:,:)     + &
                             gj5 * v(5,:,:)     + &
                             gj6 * v(6,:,:)  ) * dyinv

          g2v(ny-3,:,:) = -( gj1 * v(ny,:,:)    + &
                             gj2 * v(ny-1,:,:)  + &
                             gj3 * v(ny-2,:,:)  + &
                             gj4 * v(ny-3,:,:)  + &
                             gj5 * v(ny-4,:,:)  + &
                             gj6 * v(ny-5,:,:)  ) * dyinv

          g2v(ny-2,:,:) = -( gi1 * v(ny,:,:)    + &
                             gi2 * v(ny-1,:,:)  + &
                             gi3 * v(ny-2,:,:)  + &
                             gi4 * v(ny-3,:,:)  + &
                             gi5 * v(ny-4,:,:)  + &
                             gi6 * v(ny-5,:,:)  ) * dyinv

          g2v(ny-1,:,:) = -( gh1 * v(ny,:,:)    + &
                             gh2 * v(ny-1,:,:)  + &
                             gh3 * v(ny-2,:,:)  + &
                             gh4 * v(ny-3,:,:)  + &
                             gh5 * v(ny-4,:,:)  + &
                             gh6 * v(ny-5,:,:)  ) * dyinv

          g2v(ny,:,:)   = -( gg1 * v(ny,:,:)    + &
                             gg2 * v(ny-1,:,:)  + &
                             gg3 * v(ny-2,:,:)  + &
                             gg4 * v(ny-3,:,:)  + &
                             gg5 * v(ny-4,:,:)  + &
                             gg6 * v(ny-5,:,:)  ) * dyinv
        end if
        
!.... implement a filter of roundoff noise

!       do idof = 1, ndof
!         do i = 1, nx
!           do j = 1, ny
!             if ( abs(g1v(j,i,idof)) .lt. eps ) g1v(j,i,idof) = zero
!             if ( abs(g2v(j,i,idof)) .lt. eps ) g2v(j,i,idof) = zero
!           end do
!         end do
!       end do

        return
        end
