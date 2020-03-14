!=============================================================================!
        subroutine lhsbc1( mat, bc, vl, vml)
!  
!  Correct the LHS for boundary conditions in the \xi direction
!
!  This version supports the 2nd order LHS.
!
!  Revised: 4-16-96
!
!=============================================================================!
        use global
        use stencil
        implicit none
        
        real :: mat(ny,nx,ndof,ndof,3), bc(ny,ndof,ndof,14)
        real :: vl(ny,nx,ndof), vml(ny,nx,ndof)

!.... local variables for Riemann invariants

        real :: m1l(ny), m2l(ny), R1(ny), R2(ny), fact(ny), vint(ny,ndof)
        real :: rho(ny), p(ny), fact1(ny), fact2(ny), vn(ny), vt(ny)
        real :: ub(ny), vb(ny), tb(ny), rhob(ny), pb(ny), cb(ny)
        real :: xil(ny), etal(ny)
        real :: b, rbeta
        
        integer :: i, j, ij 
!=============================================================================!
!       L I N E A R   B O U N D A R Y   C O N D I T I O N S
!=============================================================================!
        if (linear.eq.1) then

        if (.not. yper) then

        if (Navier) then

!.... wall boundary condition

        if (wall.eq.1) then
          mat(1,:,1,:,:) = zero
          bc(1,1,:,:)    = zero
          mat(1,:,1,1,2) = one
        end if

!.... linear extrapolation of rho at the wall

        if (wall.eq.2) then
          mat(1,:,1,:,:) = zero
          bc(1,1,:,:)    = zero
          mat(1,:,1,1,3) = one
        end if

!.... no-slip

        mat(1,:,2:4,:,:) = zero
        bc(1,2:4,:,:)    = zero
        mat(1,:,2,2,2)   = one
        mat(1,:,3,3,2)   = one
        mat(1,:,4,4,2)   = one

!.... wall temperature boundary condition

        if (wallt.eq.0 .or. wallt.eq.1) then
          mat(1,:,5,:,:) = zero
          bc(1,5,:,:)    = zero
          mat(1,:,5,5,2) = one
        end if

        else
          do i = 1, nx
            mat(1,i,2,:,:) = bnb(i,2) * mat(1,i,2,:,:) - &
                             bnb(i,1) * mat(1,i,3,:,:)
          end do
          mat(1,:,3,:,:) = zero
          mat(1,:,3,2,2) = bnb(:,1)
          mat(1,:,3,3,2) = bnb(:,2)
        end if

!.... freestream hard boundary condition

        if (top.eq.0) then
          mat(ny,:,:,:,:) = zero
          bc(ny,:,:,:)    = zero
          mat(ny,:,1,1,2) = one
          mat(ny,:,2,2,2) = one
          mat(ny,:,3,3,2) = one
          mat(ny,:,4,4,2) = one
          mat(ny,:,5,5,2) = one
        end if
        
        end if          ! yper

!=============================================================================!
        
        if (.not. xper) then

!.... left boundary

          if (left.eq.0) then           ! zero disturbance
            mat(:,1,:,:,:) = zero
            bc(:,:,:,1:3)  = zero
            mat(:,1,1,1,2) = one
            mat(:,1,2,2,2) = one
            mat(:,1,3,3,2) = one
            mat(:,1,4,4,2) = one
            mat(:,1,5,5,2) = one
          else if (left.eq.1) then      ! nonreflecting BC
            mat(:,1,1,:,:) = zero
            bc(:,1,:,1:3)  = zero
            mat(:,1,1,1,2) = -one
            mat(:,1,1,2,2) = gamma * Ma * vml(:,1,1) * &
                             sqrt(vml(:,1,5)) / vml(:,1,5)
            mat(:,1,1,5,2) = -vml(:,1,1) / vml(:,1,5)
          endif        ! left

!.... Right boundary

          if (right.eq.0) then          ! zero disturbance
            mat(:,nx,:,:,:) = zero
            bc(:,:,:,12:14) = zero
            mat(:,nx,1,1,2) = one
            mat(:,nx,2,2,2) = one
            mat(:,nx,3,3,2) = one
            mat(:,nx,4,4,2) = one
            mat(:,nx,5,5,2) = one
          else if (right.eq.1) then     ! nonreflecting
            mat(:,nx,1,:,:) = zero
            bc(:,1,:,12:14) = zero
            mat(:,nx,1,1,2) = -one
            mat(:,nx,1,2,2) = gamma * Ma * vml(:,nx,1) * &
                              sqrt(vml(:,nx,5)) / vml(:,nx,5)
            mat(:,nx,1,5,2) = -vml(:,nx,1) / vml(:,nx,5)
          end if        ! right
          
        end if          ! xper

!=============================================================================!
!       N O N L I N E A R   B O U N D A R Y   C O N D I T I O N S
!=============================================================================!
        else            ! nonlinear

        if (.not. yper) then

        if (Navier) then

!.... wall boundary condition

        if (wall.eq.1 .or. wall.eq.2) then
          mat(1,:,1,:,:) = zero
          bc(1,1,:,:)    = zero
          mat(1,:,1,1,2) = one
        end if

!.... no-slip

        mat(1,:,2:4,:,:) = zero
        bc(1,2:4,:,:)    = zero
        mat(1,:,2,2,2)   = one
        mat(1,:,3,3,2)   = one
        mat(1,:,4,4,2)   = one

!.... adiabatic & isothermal

        if (wallt.eq.0 .or. wallt.eq.1) then
          mat(1,:,5,:,:) = zero
          bc(1,5,:,:)    = zero
          mat(1,:,5,5,2) = one
        end if

        else
          do i = 1, nx
            mat(1,i,2,:,:) = bnb(i,2) * mat(1,i,2,:,:) - &
                             bnb(i,1) * mat(1,i,3,:,:)
          end do
          mat(1,:,3,:,:) = zero
          mat(1,:,3,2,2) = bnb(:,1)
          mat(1,:,3,3,2) = bnb(:,2)
        end if
        
!.... Riemann Invariants on top boundary

        if (top.eq.0) then
          mat(ny,:,:,:,:) = zero
          bc(ny,:,:,:)    = zero
          mat(ny,:,1,1,2) = one
          mat(ny,:,2,2,2) = one
          mat(ny,:,3,3,2) = one
          mat(ny,:,4,4,2) = one
          mat(ny,:,5,5,2) = one
        end if
        
        end if          ! yper

!=============================================================================!
!       L e f t   B o u n d a r y
!=============================================================================!
        if (.not. xper) then

!.... symmetry plane

        if (left.eq.2) then
          mat(:,1,:,:,:) = zero
          bc(:,:,:,1:3)  = zero

          mat(:,1,1,1,2) = -gc1
          mat(:,1,1,1,3) = -gc2
          bc(:,1,1,1)    = -gc3
          bc(:,1,1,2)    = -gc4
          bc(:,1,1,3)    = -gc5

          mat(:,1,2,2,2) = -gc1
          mat(:,1,2,2,3) = -gc2
          bc(:,2,2,1)    = -gc3
          bc(:,2,2,2)    = -gc4
          bc(:,2,2,3)    = -gc5

          mat(:,1,3,3,2) = one

          mat(:,1,4,4,2) = -gc1
          mat(:,1,4,4,3) = -gc2
          bc(:,4,4,1)    = -gc3
          bc(:,4,4,2)    = -gc4
          bc(:,4,4,3)    = -gc5

          mat(:,1,5,5,2) = -gc1
          mat(:,1,5,5,3) = -gc2
          bc(:,5,5,1)    = -gc3
          bc(:,5,5,2)    = -gc4
          bc(:,5,5,3)    = -gc5
        end if

        if (left.eq.7) then             ! symmetry plane
          mat(:,1,3,:,:) = zero
          bc(:,3,:,1:3)  = zero
          mat(:,1,3,3,2) = one
        endif

!.... Left Riemann invariant boundary condition (first-order)
        
        if (left.eq.0 .and. Ma.lt.one) then

          mat(:,1,:,:,:) = zero
          bc(:,:,:,1:3)  = zero
        
          call ReimannLHS( ny, vl(:,1,:), vl(:,2,:), vl(:,3,:), &
                           bnl, x(:,1), y(:,1), mat(:,1,:,:,2), &
                           mat(:,1,:,:,3), bc(:,:,:,1) )

        end if
        
        if (.false.) then

          if (extrap.eq.0) then

!.... Zero-order extrapolation in the viscous layer

          mat(1:nbl,1,:,:,:) = zero
          bc(1:nbl,:,:,1:3)  = zero

          mat(1:nbl,1,1,1,2) = -one
          mat(1:nbl,1,2,2,2) = -one
          mat(1:nbl,1,3,3,2) = -one
          mat(1:nbl,1,4,4,2) = -one
          mat(1:nbl,1,5,5,2) = -one

          mat(1:nbl,1,1,1,3) = one
          mat(1:nbl,1,2,2,3) = one
          mat(1:nbl,1,3,3,3) = one
          mat(1:nbl,1,4,4,3) = one
          mat(1:nbl,1,5,5,3) = one
          
          end if
          
          if (extrap.eq.1) then

!.... First-order extrapolation in the viscous layer

          mat(1:nbl,1,:,:,:) = zero
          bc(1:nbl,:,:,1:3)  = zero

          mat(1:nbl,1,1,1,2) = -one
          mat(1:nbl,1,2,2,2) = -one
          mat(1:nbl,1,3,3,2) = -one
          mat(1:nbl,1,4,4,2) = -one
          mat(1:nbl,1,5,5,2) = -one

          mat(1:nbl,1,1,1,3) = two
          mat(1:nbl,1,2,2,3) = two
          mat(1:nbl,1,3,3,3) = two
          mat(1:nbl,1,4,4,3) = two
          mat(1:nbl,1,5,5,3) = two
          
          bc(1:nbl,1,1,1) = -one
          bc(1:nbl,2,2,1) = -one
          bc(1:nbl,3,3,1) = -one
          bc(1:nbl,4,4,1) = -one
          bc(1:nbl,5,5,1) = -one

          end if
        
          if (extrap.eq.2) then
          
!.... second-order extrapolation in the viscous layer

          mat(1:nbl,1,:,:,:) = zero
          bc(1:nbl,:,:,1:3)  = zero

          mat(1:nbl,1,1,1,2) = -one
          mat(1:nbl,1,2,2,2) = -one
          mat(1:nbl,1,3,3,2) = -one
          mat(1:nbl,1,4,4,2) = -one
          mat(1:nbl,1,5,5,2) = -one

          mat(1:nbl,1,1,1,3) = three
          mat(1:nbl,1,2,2,3) = three
          mat(1:nbl,1,3,3,3) = three
          mat(1:nbl,1,4,4,3) = three
          mat(1:nbl,1,5,5,3) = three
          
          bc(1:nbl,1,1,1) = -three
          bc(1:nbl,2,2,1) = -three
          bc(1:nbl,3,3,1) = -three
          bc(1:nbl,4,4,1) = -three
          bc(1:nbl,5,5,1) = -three

          bc(1:nbl,1,1,2) = one
          bc(1:nbl,2,2,2) = one
          bc(1:nbl,3,3,2) = one
          bc(1:nbl,4,4,2) = one
          bc(1:nbl,5,5,2) = one

          end if                ! extrapolation type

        end if
!=============================================================================!
!       R i g h t  B o u n d a r y
!=============================================================================!

!.... Right Riemann invariant boundary condition (first-order)
        
        if (right.eq.0 .and. Ma.lt.one) then

          mat(:,nx,:,:,:) = zero
          bc(:,:,:,12:14) = zero
                
          call ReimannLHS( ny, vl(:,nx,:), vl(:,nx-1,:), vl(:,nx-2,:), &
                           bnr, x(:,nx), y(:,nx), mat(:,nx,:,:,2), &
                           mat(:,nx,:,:,1), bc(:,:,:,14) )
        end if
          
!.... Zero'th order extrapolation in the viscous layer

          if (extrap.eq.0) then
            mat(1:nbl,nx,:,:,:) = zero
            bc(1:nbl,:,:,12:14) = zero
          
            mat(1:nbl,nx,1,1,2) = -one
            mat(1:nbl,nx,2,2,2) = -one
            mat(1:nbl,nx,3,3,2) = -one
            mat(1:nbl,nx,4,4,2) = -one
            mat(1:nbl,nx,5,5,2) = -one

            mat(1:nbl,nx,1,1,1) = one
            mat(1:nbl,nx,2,2,1) = one
            mat(1:nbl,nx,3,3,1) = one
            mat(1:nbl,nx,4,4,1) = one
            mat(1:nbl,nx,5,5,1) = one
          end if
          
!.... First-order extrapolation in the viscous layer

          if (extrap.eq.1) then
            mat(1:nbl,nx,:,:,:) = zero
            bc(1:nbl,:,:,12:14) = zero

            mat(1:nbl,nx,1,1,2) = -one
            mat(1:nbl,nx,2,2,2) = -one
            mat(1:nbl,nx,3,3,2) = -one
            mat(1:nbl,nx,4,4,2) = -one
            mat(1:nbl,nx,5,5,2) = -one

            mat(1:nbl,nx,1,1,1) = two
            mat(1:nbl,nx,2,2,1) = two
            mat(1:nbl,nx,3,3,1) = two
            mat(1:nbl,nx,4,4,1) = two
            mat(1:nbl,nx,5,5,1) = two
          
            bc(1:nbl,1,1,14) = -one
            bc(1:nbl,2,2,14) = -one
            bc(1:nbl,3,3,14) = -one
            bc(1:nbl,4,4,14) = -one
            bc(1:nbl,5,5,14) = -one
          end if                ! extrapolation type

!.... hold the IC

          if (right.eq.8) then
            mat(:,nx,:,:,:) = zero
            bc(:,:,:,12:14) = zero
            mat(:,nx,1,1,2) = -one
            mat(:,nx,2,2,2) = -one
            mat(:,nx,3,3,2) = -one
            mat(:,nx,4,4,2) = -one
            mat(:,nx,5,5,2) = -one
          end if

        end if                  ! xper
        
        end if                  ! linear
        
        return
        end
