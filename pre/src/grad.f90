!=============================================================================!
        module diff
!
!  Finite Difference coefficients
!
!=============================================================================!
!.... First derivatives
!=============================================================================!

!.... second order central difference

        real, parameter :: gd1 = -5.000000000000000E-01
        real, parameter :: gd2 =  5.000000000000000E-01
        
!.... second order one-sided difference

        real, parameter :: ge1 = -1.500000000000000E+00
        real, parameter :: ge2 =  2.000000000000000E+00
        real, parameter :: ge3 = -5.000000000000000E-01

!.... first order one-sided difference

        real, parameter :: gf1 = -1.000000000000000E+00
        real, parameter :: gf2 =  1.000000000000000E+00

        end module diff

!=============================================================================!
        subroutine grad( ndof, nx, ny, v, g1v, g2v, dx, dy, optx, opty, &
                         xper, yper)
!
!  Take the gradient of a 2-D field.
!  updated to fourth order accurate differencing
!
!=============================================================================!
        use diff
        implicit none
        
        integer ndof, nx, ny, optx, opty
        logical :: xper, yper
        real v(ny,ndof,nx), g1v(ny,ndof,nx), g2v(ny,ndof,nx)
        real dx, dy
        
        real, parameter :: one = 1.0
        real dxinv, dyinv
!=============================================================================!

        dxinv  = one / dx
        dyinv  = one / dy

!=============================================================================!
!.... compute the gradient in x
!=============================================================================!

        if (xper) then

          g1v(:,:,1)        = ( gd1 * v(:,:,nx-1)       + &
                                gd2 * v(:,:,2   )       ) * dxinv

          g1v(:,:,nx) = g1v(:,:,1)
          
        else
        
          g1v(:,:,1)      = ( gf1 * v(:,:,1)  + &
                              gf2 * v(:,:,2)  ) * dxinv
  
          g1v(:,:,nx)    = -( gf1 * v(:,:,nx  )  + &
                              gf2 * v(:,:,nx-1)  ) * dxinv

        end if

!.... interior

        g1v(:,:,2:nx-1)  = ( gd1 * v(:,:,1:nx-2)        + &
                             gd2 * v(:,:,3:nx  )        ) * dxinv

!=============================================================================!
!.... compute the gradient in y
!=============================================================================!

        if (yper) then
        
          g2v(1,:,:)       = ( gd1 * v(ny-1,:,:)        + &
                               gd2 * v(2,:,:)           ) * dyinv
  
          g2v(ny,:,:) = g2v(1,:,:)

        else

          g2v(1,:,:)     = ( gf1 * v(1,:,:)  + &
                             gf2 * v(2,:,:)  ) * dyinv
    
          g2v(ny,:,:)   = -( gf1 * v(ny  ,:,:)  + &
                             gf2 * v(ny-1,:,:)  ) * dyinv
                            
        end if
        
!.... interior

        g2v(2:ny-1,:,:) = ( gd1 * v(1:ny-2,:,:) + &
                            gd2 * v(3:ny  ,:,:) ) * dyinv 

        return
        end
