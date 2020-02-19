!=============================================================================!
        subroutine grad( ndof, nx, ny, v, g1v, g2v, dx, dy, optx, opty, &
                         xper, yper, lsym, rsym, bsym, tsym, carp)
!
!  Take the gradient of a 2-D field.
!  Updated to sixth-order accurate differencing on the interior with
!  the option for optimized fourth-order differencing.
!
!  Revised: 6-28-95
!
!=============================================================================!
        use stencil
        implicit none
        
        integer :: ndof, nx, ny, optx, opty
        logical :: xper, yper
        logical :: lsym, rsym, bsym, tsym, carp
        real    :: v(ndof,nx,ny), g1v(ndof,nx,ny), g2v(ndof,nx,ny)
        !$sgi distribute v(*,*,block), g1v(*,*,block), g2v(*,*,block)

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

!$doacross
!$omp parallel do
          do j = 1, ny

          g1v(:,1,j)        = ( gx1 * v(:,nx-3,j)       + &
                                gx2 * v(:,nx-2,j)       + &
                                gx3 * v(:,nx-1,j)       + &
                                gx4 * v(:,2,j)          + &
                                gx5 * v(:,3,j)          + &
                                gx6 * v(:,4,j)  ) 
  
          g1v(:,2,j)        = ( gx1 * v(:,nx-2,j)       + &
                                gx2 * v(:,nx-1,j)       + &
                                gx3 * v(:,1,j)          + &
                                gx4 * v(:,3,j)          + &
                                gx5 * v(:,4,j)          + &
                                gx6 * v(:,5,j)  ) 
  
          g1v(:,3,j)        = ( gx1 * v(:,nx-1,j)       + &
                                gx2 * v(:,1,j)          + &
                                gx3 * v(:,2,j)          + &
                                gx4 * v(:,4,j)          + &
                                gx5 * v(:,5,j)          + &
                                gx6 * v(:,6,j)  ) 
  
          g1v(:,nx-2,j)     = ( gx1 * v(:,nx-5,j)       + &
                                gx2 * v(:,nx-4,j)       + &
                                gx3 * v(:,nx-3,j)       + &
                                gx4 * v(:,nx-1,j)       + &
                                gx5 * v(:,1,j)          + &
                                gx6 * v(:,2,j)  ) 
  
          g1v(:,nx-1,j)     = ( gx1 * v(:,nx-4,j)       + &
                                gx2 * v(:,nx-3,j)       + &
                                gx3 * v(:,nx-2,j)       + &
                                gx4 * v(:,1,j)          + &
                                gx5 * v(:,2,j)          + &
                                gx6 * v(:,3,j)  ) 
  
          g1v(:,nx,j) = g1v(:,1,j)

          end do
!$omp end parallel do 
          
        else  if (carp ) then       !  Implement Carpenter's boundary stencil

!$doacross
!$omp parallel do
          do j = 1, ny

          g1v(:,1,j)     = ( gg1 * v(:,1,j)     + &
                             gg2 * v(:,2,j)     + &
                             gg3 * v(:,3,j)     + &
                             gg4 * v(:,4,j)     + &
                             gg5 * v(:,5,j)     + &
                             gg6 * v(:,6,j)  ) * dxinv
  
          g1v(:,2,j)     = ( gh1 * v(:,1,j)     + &
                             gh2 * v(:,2,j)     + &
                             gh3 * v(:,3,j)     + &
                             gh4 * v(:,4,j)     + &
                             gh5 * v(:,5,j)     + &
                             gh6 * v(:,6,j)  ) * dxinv

          g1v(:,3,j)     = ( gi1 * v(:,1,j)     + &
                             gi2 * v(:,2,j)     + &
                             gi3 * v(:,3,j)     + &
                             gi4 * v(:,4,j)     + &
                             gi5 * v(:,5,j)     + &
                             gi6 * v(:,6,j)  ) * dxinv

          g1v(:,4,j)     = ( gj1 * v(:,1,j)     + &
                             gj2 * v(:,2,j)     + &
                             gj3 * v(:,3,j)     + &
                             gj4 * v(:,4,j)     + &
                             gj5 * v(:,5,j)     + &
                             gj6 * v(:,6,j)  ) * dxinv

          g1v(:,nx-3,j) = -( gj1 * v(:,nx,j)    + &
                             gj2 * v(:,nx-1,j)  + &
                             gj3 * v(:,nx-2,j)  + &
                             gj4 * v(:,nx-3,j)  + &
                             gj5 * v(:,nx-4,j)  + &
                             gj6 * v(:,nx-5,j)  ) * dxinv

          g1v(:,nx-2,j) = -( gi1 * v(:,nx,j)    + &
                             gi2 * v(:,nx-1,j)  + &
                             gi3 * v(:,nx-2,j)  + &
                             gi4 * v(:,nx-3,j)  + &
                             gi5 * v(:,nx-4,j)  + &
                             gi6 * v(:,nx-5,j)  ) * dxinv

          g1v(:,nx-1,j) = -( gh1 * v(:,nx,j)    + &
                             gh2 * v(:,nx-1,j)  + &
                             gh3 * v(:,nx-2,j)  + &
                             gh4 * v(:,nx-3,j)  + &
                             gh5 * v(:,nx-4,j)  + &
                             gh6 * v(:,nx-5,j)  ) * dxinv

          g1v(:,nx,j)   = -( gg1 * v(:,nx,j)    + &
                             gg2 * v(:,nx-1,j)  + &
                             gg3 * v(:,nx-2,j)  + &
                             gg4 * v(:,nx-3,j)  + &
                             gg5 * v(:,nx-4,j)  + &
                             gg6 * v(:,nx-5,j)  ) * dxinv
          end do
!$omp end parallel do 

        else                  ! boundary differences
        
!$doacross
!$omp parallel do
          do j = 1, ny

!         g1v(:,1,j)      = ( v(:,2,j) - v(:,1,j) ) * dxinv
          g1v(:,1,j)      = ( gc1 * v(:,1,j)  + &
                              gc2 * v(:,2,j)  + &
                              gc3 * v(:,3,j)  + &
                              gc4 * v(:,4,j)  + &
                              gc5 * v(:,5,j)  ) * dxinv
!         g1v(:,1,j)      = ( ge1 * v(:,1,j)  + &
!                             ge2 * v(:,2,j)  + &
!                             ge3 * v(:,3,j)  + &
!                             ge4 * v(:,4,j)  + &
!                             ge5 * v(:,5,j)  + &
!                             ge6 * v(:,6,j)  + &
!                             ge7 * v(:,7,j)  ) * dxinv
  
!         g1v(:,2,j)      = ( v(:,3,j) - v(:,1,j) ) * pt5 * dxinv
          g1v(:,2,j)      = ( gb1 * v(:,1,j)  + &
                              gb2 * v(:,2,j)  + &
                              gb3 * v(:,3,j)  + &
                              gb4 * v(:,4,j)  + &
                              gb5 * v(:,5,j)  ) * dxinv
!         g1v(:,2,j)      = ( gf1 * v(:,1,j)  + &
!                             gf2 * v(:,2,j)  + &
!                             gf3 * v(:,3,j)  + &
!                             gf4 * v(:,4,j)  + &
!                             gf5 * v(:,5,j)  + &
!                             gf6 * v(:,6,j)  + &
!                             gf7 * v(:,7,j)  ) * dxinv

          g1v(:,3,j)      = ( ga1 * v(:,1,j)  + &
                              ga2 * v(:,2,j)  + &
                              ga3 * v(:,4,j)  + &
                              ga4 * v(:,5,j)  ) * dxinv
          
          g1v(:,nx-2,j)   = ( ga1 * v(:,nx-4,j)  + &
                              ga2 * v(:,nx-3,j)  + &
                              ga3 * v(:,nx-1,j)  + &
                              ga4 * v(:,nx  ,j)  ) * dxinv
  
          g1v(:,nx-1,j)  = -( gb1 * v(:,nx  ,j)  + &
                              gb2 * v(:,nx-1,j)  + &
                              gb3 * v(:,nx-2,j)  + &
                              gb4 * v(:,nx-3,j)  + &
                              gb5 * v(:,nx-4,j)  ) * dxinv
  
          g1v(:,nx,j)    = -( gc1 * v(:,nx  ,j)  + &
                              gc2 * v(:,nx-1,j)  + &
                              gc3 * v(:,nx-2,j)  + &
                              gc4 * v(:,nx-3,j)  + &
                              gc5 * v(:,nx-4,j)  ) * dxinv

          end do
!$omp end parallel do 

        end if

!.... implement the symmetric conditions

        if (lsym) then
!$doacross local( idof, isign )
!$omp parallel do private(idof, isign)
          do j = 1, ny
          do idof = 1, ndof
            if (idof .eq. 3) then
              isign = -one
            else
              isign = one
            end if
  
            g1v(idof,1,j)     = ( isign * gx1 * v(idof,4,j)     + &
                                  isign * gx2 * v(idof,3,j)     + &
                                  isign * gx3 * v(idof,2,j)     + &
                                  gx4 * v(idof,2,j)             + &
                                  gx5 * v(idof,3,j)             + &
                                  gx6 * v(idof,4,j)             ) 
    
            g1v(idof,2,j)     = ( isign * gx1 * v(idof,3,j)     + &
                                  isign * gx2 * v(idof,2,j)     + &
                                  gx3 * v(idof,1,j)             + &
                                  gx4 * v(idof,3,j)             + &
                                  gx5 * v(idof,4,j)             + &
                                  gx6 * v(idof,5,j)             ) 
    
            g1v(idof,3,j)     = ( isign * gx1 * v(idof,2,j)     + &
                                  gx2 * v(idof,1,j)             + &
                                  gx3 * v(idof,2,j)             + &
                                  gx4 * v(idof,4,j)             + &
                                  gx5 * v(idof,5,j)             + &
                                  gx6 * v(idof,6,j)             ) 
          end do
          end do
!$omp end parallel do 

        end if
        
        if (rsym) then
!$doacross local( idof, isign )
!$omp parallel do private(idof, isign)
          do j = 1, ny
          do idof = 1, ndof
            if (idof .eq. 3) then
              isign = -one
            else
              isign = one
            end if
  
            g1v(idof,nx-2,j)  = ( gx1 * v(idof,nx-5,j)          + &
                                  gx2 * v(idof,nx-4,j)          + &
                                  gx3 * v(idof,nx-3,j)          + &
                                  gx4 * v(idof,nx-1,j)          + &
                                  gx5 * v(idof,nx,  j)          + &
                                  isign * gx6 * v(idof,nx-1,j)  ) 
    
            g1v(idof,nx-1,j)  = ( gx1 * v(idof,nx-4,j)          + &
                                  gx2 * v(idof,nx-3,j)          + &
                                  gx3 * v(idof,nx-2,j)          + &
                                  gx4 * v(idof,nx,  j)          + &
                                  isign * gx5 * v(idof,nx-1,j)  + &
                                  isign * gx6 * v(idof,nx-2,j)  ) 
    
            g1v(idof,nx,j)    = ( gx1 * v(idof,nx-3,j)          + &
                                  gx2 * v(idof,nx-2,j)          + &
                                  gx3 * v(idof,nx-1,j)          + &
                                  isign * gx4 * v(idof,nx-1,j)  + &
                                  isign * gx5 * v(idof,nx-2,j)  + &
                                  isign * gx6 * v(idof,nx-3,j)  ) 
          end do
          end do
!$omp end parallel do 

        end if

!.... interior

!$doacross local(i,idof)
!$omp parallel do private(i,idof)
        do j = 1, ny
          do i = 4, nx-3
            do idof = 1, ndof
              g1v(idof,i,j)  = ( gx1 * v(idof,i-3,j)    + &
                                 gx2 * v(idof,i-2,j)    + &
                                 gx3 * v(idof,i-1,j)    + &
                                 gx4 * v(idof,i+1,j)    + &
                                 gx5 * v(idof,i+2,j)    + &
                                 gx6 * v(idof,i+3,j)    ) 
            end do
          end do
        end do
!$omp end parallel do 

!=============================================================================!
!.... compute the gradient in y
!=============================================================================!

        if (yper) then
        
!$doacross
!$omp parallel do
          do i = 1, nx

          g2v(:,i,1)        = ( gy1 * v(:,i,ny-3)       + &
                                gy2 * v(:,i,ny-2)       + &
                                gy3 * v(:,i,ny-1)       + &
                                gy4 * v(:,i,2)          + &
                                gy5 * v(:,i,3)          + &
                                gy6 * v(:,i,4)          ) 
  
          g2v(:,i,2)        = ( gy1 * v(:,i,ny-2)       + &
                                gy2 * v(:,i,ny-1)       + &
                                gy3 * v(:,i,1)          + &
                                gy4 * v(:,i,3)          + &
                                gy5 * v(:,i,4)          + &
                                gy6 * v(:,i,5)  ) 
  
          g2v(:,i,3)        = ( gy1 * v(:,i,ny-1)       + &
                                gy2 * v(:,i,1)          + &
                                gy3 * v(:,i,2)          + &
                                gy4 * v(:,i,4)          + &
                                gy5 * v(:,i,5)          + &
                                gy6 * v(:,i,6)  ) 
  
          g2v(:,i,ny-2)     = ( gy1 * v(:,i,ny-5)       + &
                                gy2 * v(:,i,ny-4)       + &
                                gy3 * v(:,i,ny-3)       + &
                                gy4 * v(:,i,ny-1)       + &
                                gy5 * v(:,i,1)          + &
                                gy6 * v(:,i,2)  ) 
  
          g2v(:,i,ny-1)     = ( gy1 * v(:,i,ny-4)       + &
                                gy2 * v(:,i,ny-3)       + &
                                gy3 * v(:,i,ny-2)       + &
                                gy4 * v(:,i,1)          + &
                                gy5 * v(:,i,2)          + &
                                gy6 * v(:,i,3)  ) 
  
          g2v(:,i,ny) = g2v(:,i,1)

          end do
!$omp end parallel do 

        else  if (carp ) then       !  Implement Carpenter's boundary stencil

!$doacross
!$omp parallel do
          do i = 1, nx

          g2v(:,i,1)     = ( gg1 * v(:,i,1)     + &
                             gg2 * v(:,i,2)     + &
                             gg3 * v(:,i,3)     + &
                             gg4 * v(:,i,4)     + &
                             gg5 * v(:,i,5)     + &
                             gg6 * v(:,i,6)  ) * dyinv
  
          g2v(:,i,2)     = ( gh1 * v(:,i,1)     + &
                             gh2 * v(:,i,2)     + &
                             gh3 * v(:,i,3)     + &
                             gh4 * v(:,i,4)     + &
                             gh5 * v(:,i,5)     + &
                             gh6 * v(:,i,6)  ) * dyinv

          g2v(:,i,3)     = ( gi1 * v(:,i,1)     + &
                             gi2 * v(:,i,2)     + &
                             gi3 * v(:,i,3)     + &
                             gi4 * v(:,i,4)     + &
                             gi5 * v(:,i,5)     + &
                             gi6 * v(:,i,6)  ) * dyinv

          g2v(:,i,4)     = ( gj1 * v(:,i,1)     + &
                             gj2 * v(:,i,2)     + &
                             gj3 * v(:,i,3)     + &
                             gj4 * v(:,i,4)     + &
                             gj5 * v(:,i,5)     + &
                             gj6 * v(:,i,6)  ) * dyinv

          g2v(:,i,ny-3) = -( gj1 * v(:,i,ny)    + &
                             gj2 * v(:,i,ny-1)  + &
                             gj3 * v(:,i,ny-2)  + &
                             gj4 * v(:,i,ny-3)  + &
                             gj5 * v(:,i,ny-4)  + &
                             gj6 * v(:,i,ny-5)  ) * dyinv

          g2v(:,i,ny-2) = -( gi1 * v(:,i,ny)    + &
                             gi2 * v(:,i,ny-1)  + &
                             gi3 * v(:,i,ny-2)  + &
                             gi4 * v(:,i,ny-3)  + &
                             gi5 * v(:,i,ny-4)  + &
                             gi6 * v(:,i,ny-5)  ) * dyinv

          g2v(:,i,ny-1) = -( gh1 * v(:,i,ny)    + &
                             gh2 * v(:,i,ny-1)  + &
                             gh3 * v(:,i,ny-2)  + &
                             gh4 * v(:,i,ny-3)  + &
                             gh5 * v(:,i,ny-4)  + &
                             gh6 * v(:,i,ny-5)  ) * dyinv

          g2v(:,i,ny)   = -( gg1 * v(:,i,ny)    + &
                             gg2 * v(:,i,ny-1)  + &
                             gg3 * v(:,i,ny-2)  + &
                             gg4 * v(:,i,ny-3)  + &
                             gg5 * v(:,i,ny-4)  + &
                             gg6 * v(:,i,ny-5)  ) * dyinv
          end do
!$omp end parallel do 

        else     ! normal boundary differences (must use for implicit)

!$doacross
!$omp parallel do
          do i = 1, nx

          g2v(:,i,1)     = ( gc1 * v(:,i,1)     + &
                             gc2 * v(:,i,2)     + &
                             gc3 * v(:,i,3)     + &
                             gc4 * v(:,i,4)     + &
                             gc5 * v(:,i,5)     ) * dyinv

          g2v(:,i,2)     = ( gb1 * v(:,i,1)     + &
                             gb2 * v(:,i,2)     + &
                             gb3 * v(:,i,3)     + &
                             gb4 * v(:,i,4)     + &
                             gb5 * v(:,i,5)     ) * dyinv
  
          g2v(:,i,3)     = ( ga1 * v(:,i,1)     + &
                             ga2 * v(:,i,2)     + &
                             ga3 * v(:,i,4)     + &
                             ga4 * v(:,i,5)     ) * dyinv
  
          g2v(:,i,ny-2)  = ( ga1 * v(:,i,ny-4)  + &
                             ga2 * v(:,i,ny-3)  + &
                             ga3 * v(:,i,ny-1)  + &
                             ga4 * v(:,i,ny  )  ) * dyinv
  
          g2v(:,i,ny-1) = -( gb1 * v(:,i,ny  )  + &
                             gb2 * v(:,i,ny-1)  + &
                             gb3 * v(:,i,ny-2)  + &
                             gb4 * v(:,i,ny-3)  + &
                             gb5 * v(:,i,ny-4)  ) * dyinv
  
          g2v(:,i,ny)   = -( gc1 * v(:,i,ny  )  + &
                             gc2 * v(:,i,ny-1)  + &
                             gc3 * v(:,i,ny-2)  + &
                             gc4 * v(:,i,ny-3)  + &
                             gc5 * v(:,i,ny-4)  ) * dyinv
                            
          end do
!$omp end parallel do 

        end if
        
!.... interior

!$doacross local(j,idof)
!$omp parallel do private(j,idof)
        do i = 1, nx
          do j = 4, ny-3
            do idof = 1, ndof
              g2v(idof,i,j) = ( gy1 * v(idof,i,j-3)     + &
                                gy2 * v(idof,i,j-2)     + &
                                gy3 * v(idof,i,j-1)     + &
                                gy4 * v(idof,i,j+1)     + &
                                gy5 * v(idof,i,j+2)     + &
                                gy6 * v(idof,i,j+3)     ) 
            end do
          end do
        end do
!$omp end parallel do 
        
!.... implement a filter of roundoff noise

        if (.false.) then

!$omp parallel do private(i,idof)
          do j = 1, ny
            do i = 1, nx
              do idof = 1, ndof
                if ( abs(g1v(idof,i,j)) .lt. eps ) g1v(idof,i,j) = zero
                if ( abs(g2v(idof,i,j)) .lt. eps ) g2v(idof,i,j) = zero
              end do
            end do
          end do

        end if

        return
        end
