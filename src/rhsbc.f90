!=============================================================================!
        subroutine rhsBC(rl, vl, vml) 
!  
!  Satisfy the boundary conditions 
!  
!=============================================================================!
        use global
        use stencil
        use pot
        implicit none

        real :: rl(ndof,nx,ny), vl(ndof,nx,ny), vml(ndof,nx,ny)
        
        real, allocatable :: p(:), pnorm(:)

        integer :: i, j
!=============================================================================!
!       L I N E A R   B O U N D A R Y   C O N D I T I O N S
!=============================================================================!
        if (linear.eq.1) then

        if (yper) then
        
          rl(:,:,ny) = zero
        
        else                    ! yper

!=============================================================================!
!       W a l l
!=============================================================================!
        if (Navier) then

!.... density boundary condition

          if (wall.eq.1) then
            rl(1,:,1) = gc1 * vl(1,:,1) + gc2 * vl(1,:,2) + &
                        gc3 * vl(1,:,3) + gc4 * vl(1,:,4) + &
                        gc5 * vl(1,:,5)
          end if
          
!.... linear extrapolation of rho at the wall

          if (wall.eq.2) then
            rl(1,:,1) = vl(1,:,1) - two * vl(1,:,2) + vl(1,:,3)
          end if

!.... no-slip boundary condition

          rl(2:ndof-1,:,1)  = zero

!.... isothermal wall

          if (wallt.eq.0) rl(ndof,:,1) = zero

!.... adiabatic boundary condition

          if (wallt.eq.1) then
            rl(ndof,:,1) = gc1 * vl(ndof,:,1) + gc2 * vl(ndof,:,2) + &
                           gc3 * vl(ndof,:,3) + gc4 * vl(ndof,:,4) + &
                           gc5 * vl(ndof,:,5)
          end if
        
        else            ! inviscid wall
        
          allocate( p(nx), pnorm(nx) )
          call lwallbc(vl, vml, p, pnorm)
          do i = 1, nx
            rl(1,i,1) = ( gc1 * vl(1,i,1) + gc2 * vl(1,i,2) + &
                          gc3 * vl(1,i,3) + gc4 * vl(1,i,4) + &
                          gc5 * vl(1,i,5) ) / deta - Pnorm(i)
          end do
          deallocate( p, pnorm )

          do i = 1, nx
            if ( abs(bnb(i,2)) .ge. abs(bnb(i,1)) ) then
              rl(2,i,1) = bnb(i,2) * rl(2,i,1) - bnb(i,1) * rl(3,i,1)
              rl(3,i,1) = -( bnb(i,1) * vl(2,i,1) + bnb(i,2) * vl(3,i,1) )
            else
              rl(3,i,1) = bnb(i,2) * rl(2,i,1) - bnb(i,1) * rl(3,i,1)
              rl(2,i,1) = -( bnb(i,1) * vl(2,i,1) + bnb(i,2) * vl(3,i,1) )
            end if
          end do

        end if          ! Navier
!=============================================================================!
!       T o p
!=============================================================================!

!.... freestream zero disturbance boundary conditions

          if (top.eq.0) then
            rl(:,:,ny) = zero
          end if
        
        end if                  ! yper

        if (xper) then
          rl(:,nx,:) = zero
        else                    ! xper
        
!=============================================================================!
!       L e f t   B o u n d a r y
!=============================================================================!
        if (left.eq.0) then

          rl(:,1,:)  = zero

        else if (left.eq.1) then

          rl(1,1,nbl+1:ny) = vl(1,1,nbl+1:ny) - &
                  (gamma * Ma * vml(1,1,nbl+1:ny) * &
                  sqrt(vml(5,1,nbl+1:ny)) * vl(2,1,nbl+1:ny) - &
                  vml(1,1,nbl+1:ny) * vl(5,1,nbl+1:ny)) / vml(5,1,nbl+1:ny)

!         rl(:,1,1:nbl) = vl(:,1,1:nbl) - vl(:,2,1:nbl)

        else if (left.eq.2) then

          do j = 1, ny
!           rl(1,1,j) = gc1 * vl(1,1,j) + gc2 * vl(1,2,j)  + &
!                       gc3 * vl(1,3,j) + gc4 * vl(1,4,j)  + &
!                       gc5 * vl(1,5,j)

            rl(2,1,j) = gc1 * vl(2,1,j) + gc2 * vl(2,2,j)  + &
                        gc3 * vl(2,3,j) + gc4 * vl(2,4,j)  + &
                        gc5 * vl(2,5,j)
                            
            rl(3,1,j) = zero

            rl(4,1,j) = gc1 * vl(4,1,j) + gc2 * vl(4,2,j)  + &
                        gc3 * vl(4,3,j) + gc4 * vl(4,4,j)  + &
                        gc5 * vl(4,5,j)

            rl(5,1,j) = gc1 * vl(5,1,j) + gc2 * vl(5,2,j)  + &
                        gc3 * vl(5,3,j) + gc4 * vl(5,4,j)  + &
                        gc5 * vl(5,5,j)
          end do

        else if (left.eq.7) then           ! symmetry boundary
          rl(3,1,:) = zero
        end if

!=============================================================================!
!       R i g h t   B o u n d a r y
!=============================================================================!
        if (right.eq.0) then

          rl(:,nx,:) = zero

        else if (right.eq.1) then

          rl(1,nx,nbl+1:ny) = vl(1,nx,nbl+1:ny) - &
            (gamma * Ma * vml(1,nx,nbl+1:ny) * &
            sqrt(vml(5,nx,nbl+1:ny)) * vl(2,nx,nbl+1:ny) - &
            vml(1,nx,nbl+1:ny) * vl(5,nx,nbl+1:ny)) /vml(5,nx,nbl+1:ny)

!         rl(:,nx,1:nbl) = vl(:,nx,1:nbl) - vl(:,nx-1,1:nbl)

        end if

        end if                  ! xper
                        
!=============================================================================!
!       N O N L I N E A R   B O U N D A R Y   C O N D I T I O N S
!=============================================================================!
        else                    ! linear = 0    

        if (yper) then
        
          rl(:,:,ny) = zero
        
        else                    ! yper

!=============================================================================!
!       W a l l
!=============================================================================!
        if (Navier) then

!.... no-slip boundary condition

          rl(2:ndof-1,:,1) = zero

!.... density boundary condition

          if (wall.eq.1) then
            rl(1,:,1) = gc1 * vl(1,:,1) + gc2 * vl(1,:,2) + &
                        gc3 * vl(1,:,3) + gc4 * vl(1,:,4) + &
                        gc5 * vl(1,:,5) 
          end if

!.... Wall normal momentum equation

          if (wall.eq.2) then
            allocate( p(nx), pnorm(nx) )
            call wallbc(vl, p, pnorm)
            do i = 1, nx
              rl(1,i,1) = ( gc1 * vl(1,i,1) + gc2 * vl(1,i,2) + &
                            gc3 * vl(1,i,3) + gc4 * vl(1,i,4) + &
                            gc5 * vl(1,i,5) ) / deta -          &
                            gamma * Ma**2 * Pnorm(i) / vl(5,i,1)
            end do
            deallocate( p, pnorm )
          end if
          
!.... Third order extrapolation to the wall for density
!.... (in computational space)

          if (wall.eq.3) then
            rl(1,:,1) = vl(1,:,1) - (4.0*vl(1,:,2) - 6.0*vl(1,:,3) + &
                                     4.0*vl(1,:,4) -     vl(1,:,5) )
          end if

!.... isothermal wall

          if (wallt.eq.0) then
            rl(ndof,:,1) = vl(ndof,:,1) - ( one + pt5 * gamma1 * Ma**2 * &
                           (one + tan(theta)**2) * sqrt(Pr) )
          end if

!.... adiabatic boundary condition

          if (wallt.eq.1) then
            rl(ndof,:,1) = gc1 * vl(ndof,:,1) + gc2 * vl(ndof,:,2) + &
                           gc3 * vl(ndof,:,3) + gc4 * vl(ndof,:,4) + &
                           gc5 * vl(ndof,:,5)
          end if

        else    ! inviscid
        
!.... rl(2,:,1) = wall tangent momentum
!.... rl(3,:,1) = wall normal velocity
          
          if (.false.) then
            rl(2,:,1) =  ( bnb(:,2) * rl(2,:,1) - bnb(:,1) * rl(3,:,1) )
            rl(3,:,1) = -( bnb(:,1) * rl(2,:,1) + bnb(:,2) * rl(3,:,1) )
          end if

!.... Add density extrapolation

          if (wall.eq.3) then

          end if

        end if  ! Navier
!=============================================================================!
!       T o p
!=============================================================================!
        if (Ma.gt.one) then

          rl(1,:,ny) = zero             ! set all variables
          rl(2,:,ny) = zero
          rl(3,:,ny) = zero
          rl(4,:,ny) = zero
          rl(5,:,ny) = zero

        else    ! Ma

          if (top.eq.0) then
            call ReimannRHS( nx, vl(:,:,ny), vl(:,:,ny-1), vl(:,:,ny-2), &
                             bnt, x(:,ny), y(:,ny), rl(:,:,ny), &
                             rhobt, ubt, vbt, wbt, tbt, pbt, cbt )
          end if
        
        end if                  ! Ma
        
        end if                  ! yper
        
!=============================================================================!

        if (xper) then
        
          rl(:,nx,:) = zero
        
        else                    ! xper
        
!=============================================================================!
!       L e f t  B o u n d a r y
!=============================================================================!

        if (Ma.lt.one) then

!.... first-order Riemann

          if (left.eq.0) then
            call ReimannRHS( ny, vl(:,1,:), vl(:,2,:), vl(:,3,:), &
                             bnl, x(1,:), y(1,:), rl(:,1,:), &
                             rhobl, ubl, vbl, wbl, tbl, pbl, cbl )
          end if

!.... Symmetry plane

          if (left.eq.2) then
            do j = 1, ny
              rl(1,1,j) = gc1 * vl(1,1,j) + gc2 * vl(1,2,j)  + &
                          gc3 * vl(1,3,j) + gc4 * vl(1,4,j)  + &
                          gc5 * vl(1,5,j)
  
              rl(2,1,j) = gc1 * vl(2,1,j) + gc2 * vl(2,2,j)  + &
                          gc3 * vl(2,3,j) + gc4 * vl(2,4,j)  + &
                          gc5 * vl(2,5,j)
                              
              rl(3,1,j) = zero
              rl(4,1,j) = gc1 * vl(4,1,j) + gc2 * vl(4,2,j)  + &
                          gc3 * vl(4,3,j) + gc4 * vl(4,4,j)  + &
                          gc5 * vl(4,5,j)
              rl(5,1,j) = gc1 * vl(5,1,j) + gc2 * vl(5,2,j)  + &
                          gc3 * vl(5,3,j) + gc4 * vl(5,4,j)  + &
                          gc5 * vl(5,5,j)
            end do
          end if

!.... Symmetry boundary

          if (left.eq.7) then
            rl(3,1,:) = zero
          end if
          
        end if  ! Ma

        if (.false.) then
          
!.... Zero'th order extrapolation in the viscous layers

          if (extrap.eq.0) then
            rl(1,1,1:nbl) = vl(1,1,1:nbl) - vl(1,2,1:nbl)
            rl(2,1,1:nbl) = vl(2,1,1:nbl) - vl(2,2,1:nbl)
            rl(3,1,1:nbl) = vl(3,1,1:nbl) - vl(3,2,1:nbl)
            rl(4,1,1:nbl) = vl(4,1,1:nbl) - vl(4,2,1:nbl)
            rl(5,1,1:nbl) = vl(5,1,1:nbl) - vl(5,2,1:nbl)
          
!.... First-order extrapolation in the viscous layers

          else if (extrap.eq.1) then
            rl(:,1,1:nbl) = vl(:,1,1:nbl) - &
                            ( two * vl(:,2,1:nbl) - vl(:,3,1:nbl) )
          
!.... Second-order extrapolation in the viscous layers (no good)

          else if (extrap.eq.2) then
            rl(:,1,1:nbl) = vl(:,1,1:nbl) - &
                     ( three * vl(:,2,1:nbl) - three * vl(:,3,1:nbl) + &
                       vl(:,4,1:nbl) )
          end if
          
        end if
!=============================================================================!
!       R i g h t  B o u n d a r y
!=============================================================================!

        if (Ma.lt.one) then

!.... First-order Riemann

          if (right.eq.0) then
            call ReimannRHS( ny, vl(:,nx,:), vl(:,nx-1,:), vl(:,nx-2,:), &
                             bnr, x(nx,:), y(nx,:), rl(:,nx,:), &
                             rhobr, ubr, vbr, wbr, tbr, pbr, cbr )
          end if

!.... Symmetry boundary

          if (right.eq.7) then
            rl(3,nx,:) = zero
          end if

        end if          ! Ma

!.... Zero'th order extrapolation in the viscous layers

        if (extrap.eq.0) then
          rl(:,nx,1:nbl) = vl(:,nx,1:nbl) - vl(:,nx-1,1:nbl)
          
!.... First-order extrapolation in the viscous layers

        else if (extrap.eq.1) then
          rl(:,nx,1:nbl) = vl(:,nx,1:nbl) - &
                         ( two * vl(:,nx-1,1:nbl) - vl(:,nx-2,1:nbl) )
        end if
        
        if (right.eq.8) then            ! hold initial condition
          rl(:,nx,:) = zero
        end if

        end if                  ! xper
        
        end if                  ! linear

        return
        end
