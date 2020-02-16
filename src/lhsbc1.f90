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
        use pot
        implicit none
        
        real :: mat(3,ndof,ndof,nx,ny), bc(14,ndof,ndof,ny)
        real :: vl(ndof,nx,ny), vml(ndof,nx,ny)

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
          mat(:,1,:,:,1) = zero
          bc(:,1,:,1)    = zero
          mat(2,1,1,:,1) = one
        end if

!.... linear extrapolation of rho at the wall

        if (wall.eq.2) then
          mat(:,1,:,:,1) = zero
          bc(:,1,:,1)    = zero
          mat(3,1,1,:,1) = one
        end if

!.... no-slip

        mat(:,2:4,:,:,1) = zero
        bc(:,2:4,:,1)    = zero
        mat(2,2,2,:,1)   = one
        mat(2,3,3,:,1)   = one
        mat(2,4,4,:,1)   = one

!.... wall temperature boundary condition

        if (wallt.eq.0 .or. wallt.eq.1) then
          mat(:,5,:,:,1) = zero
          bc(:,5,:,1)    = zero
          mat(2,5,5,:,1) = one
        end if

        else
          do i = 1, nx
            mat(:,2,:,i,1) = bnb(i,2) * mat(:,2,:,i,1) - &
                             bnb(i,1) * mat(:,3,:,i,1)
          end do
          mat(:,3,:,:,1) = zero
          mat(2,3,2,:,1) = bnb(:,1)
          mat(2,3,3,:,1) = bnb(:,2)
        end if

!.... freestream hard boundary condition

        if (top.eq.0) then
          mat(:,:,:,:,ny) = zero
          bc(:,:,:,ny)    = zero
          mat(2,1,1,:,ny) = one
          mat(2,2,2,:,ny) = one
          mat(2,3,3,:,ny) = one
          mat(2,4,4,:,ny) = one
          mat(2,5,5,:,ny) = one
        end if
        
        end if          ! yper

!=============================================================================!
        
        if (.not. xper) then

!.... left boundary

          if (left.eq.0) then           ! zero disturbance
            mat(:,:,:,1,:) = zero
            bc(1:3,:,:,:)  = zero
            mat(2,1,1,1,:) = one
            mat(2,2,2,1,:) = one
            mat(2,3,3,1,:) = one
            mat(2,4,4,1,:) = one
            mat(2,5,5,1,:) = one
          else if (left.eq.1) then      ! nonreflecting BC
            mat(:,1,:,1,:) = zero
            bc(1:3,1,:,:)  = zero
            mat(2,1,1,1,:) = -one
            mat(2,1,2,1,:) = gamma * Ma * vml(1,1,:) * &
                             sqrt(vml(5,1,:)) / vml(5,1,:)
            mat(2,1,5,1,:) = -vml(1,1,:) / vml(5,1,:)
          endif        ! left

!.... Right boundary

          if (right.eq.0) then          ! zero disturbance
            mat(:,:,:,nx,:) = zero
            bc(12:14,:,:,:) = zero
            mat(2,1,1,nx,:) = one
            mat(2,2,2,nx,:) = one
            mat(2,3,3,nx,:) = one
            mat(2,4,4,nx,:) = one
            mat(2,5,5,nx,:) = one
          else if (right.eq.1) then     ! nonreflecting
            mat(:,1,:,nx,:) = zero
            bc(12:14,1,:,:) = zero
            mat(2,1,1,nx,:) = -one
            mat(2,1,2,nx,:) = gamma * Ma * vml(1,nx,:) * &
                              sqrt(vml(5,nx,:)) / vml(5,nx,:)
            mat(2,1,5,nx,:) = -vml(1,nx,:) / vml(5,nx,:)
          end if        ! right
          
        end if          ! xper

!=============================================================================!
!       N O N L I N E A R   B O U N D A R Y   C O N D I T I O N S
!=============================================================================!
        else            ! nonlinear

        if (.not. yper) then

        if (Navier) then

!.... wall boundary condition

        if (wall.eq.1 .or. wall.eq.2 .or. wall.eq.3) then
          mat(:,1,:,:,1) = zero
          bc(:,1,:,1)    = zero
          mat(2,1,1,:,1) = one
        end if

!.... no-slip

        mat(:,2:4,:,:,1) = zero
        bc(:,2:4,:,1)    = zero
        mat(2,2,2,:,1)   = one
        mat(2,3,3,:,1)   = one
        mat(2,4,4,:,1)   = one

!.... adiabatic & isothermal

        if (wallt.eq.0 .or. wallt.eq.1) then
          mat(:,5,:,:,1) = zero
          bc(:,5,:,1)    = zero
          mat(2,5,5,:,1) = one
        end if

        else     ! inviscid wall

          do i = 1, nx
            mat(:,2,:,i,1) = bnb(i,2) * mat(:,2,:,i,1) - &
                             bnb(i,1) * mat(:,3,:,i,1)
          end do
          mat(:,3,:,:,1) = zero
          mat(2,3,2,:,1) = bnb(:,1)
          mat(2,3,3,:,1) = bnb(:,2)

        end if
        
!.... Riemann Invariants on top boundary

        if (top.eq.0) then
          mat(:,:,:,:,ny) = zero
          bc(:,:,:,ny)    = zero
          mat(2,1,1,:,ny) = one
          mat(2,2,2,:,ny) = one
          mat(2,3,3,:,ny) = one
          mat(2,4,4,:,ny) = one
          mat(2,5,5,:,ny) = one
        end if
        
        end if          ! yper

!=============================================================================!
!       L e f t   B o u n d a r y
!=============================================================================!
        if (.not. xper) then

!.... symmetry plane

        if (left.eq.2) then
          mat(:,:,:,1,:) = zero
          bc(1:3,:,:,:)  = zero

          mat(2,1,1,1,:) = -gc1
          mat(3,1,1,1,:) = -gc2
          bc(1,1,1,:)    = -gc3
          bc(2,1,1,:)    = -gc4
          bc(3,1,1,:)    = -gc5

          mat(2,2,2,1,:) = -gc1
          mat(3,2,2,1,:) = -gc2
          bc(1,2,2,:)    = -gc3
          bc(2,2,2,:)    = -gc4
          bc(3,2,2,:)    = -gc5

          mat(2,3,3,1,:) = one

          mat(2,4,4,1,:) = -gc1
          mat(3,4,4,1,:) = -gc2
          bc(1,4,4,:)    = -gc3
          bc(2,4,4,:)    = -gc4
          bc(3,4,4,:)    = -gc5

          mat(2,5,5,1,:) = -gc1
          mat(3,5,5,1,:) = -gc2
          bc(1,5,5,:)    = -gc3
          bc(2,5,5,:)    = -gc4
          bc(3,5,5,:)    = -gc5
        end if

        if (left.eq.7) then             ! symmetry plane
          mat(:,3,:,1,:) = zero
          bc(1:3,3,:,:)  = zero
          mat(2,3,3,1,:) = one
        endif

!.... Left Riemann invariant boundary condition (first-order)
        
        if (left.eq.0 .and. Ma.lt.one) then

          mat(:,:,:,1,:) = zero
          bc(1:3,:,:,:)  = zero
        
          call ReimannLHS( ny, vl(:,1,:), vl(:,2,:), vl(:,3,:), &
                           bnl, x(1,:), y(1,:), mat(2,:,:,1,:), &
                           mat(3,:,:,1,:), bc(1,:,:,:), &
                           rhobl, ubl, vbl, wbl, tbl, pbl, cbl )

        end if
        
        if (.false.) then

          if (extrap.eq.0) then

!.... Zero-order extrapolation in the viscous layer

          mat(:,:,:,1,1:nbl) = zero
          bc(1:3,:,:,1:nbl)  = zero

          mat(2,1,1,1,1:nbl) = -one
          mat(2,2,2,1,1:nbl) = -one
          mat(2,3,3,1,1:nbl) = -one
          mat(2,4,4,1,1:nbl) = -one
          mat(2,5,5,1,1:nbl) = -one

          mat(3,1,1,1,1:nbl) = one
          mat(3,2,2,1,1:nbl) = one
          mat(3,3,3,1,1:nbl) = one
          mat(3,4,4,1,1:nbl) = one
          mat(3,5,5,1,1:nbl) = one
          
          end if
          
          if (extrap.eq.1) then

!.... First-order extrapolation in the viscous layer

          mat(:,:,:,1,1:nbl) = zero
          bc(1:3,:,:,1:nbl)  = zero

          mat(2,1,1,1,1:nbl) = -one
          mat(2,2,2,1,1:nbl) = -one
          mat(2,3,3,1,1:nbl) = -one
          mat(2,4,4,1,1:nbl) = -one
          mat(2,5,5,1,1:nbl) = -one

          mat(3,1,1,1,1:nbl) = two
          mat(3,2,2,1,1:nbl) = two
          mat(3,3,3,1,1:nbl) = two
          mat(3,4,4,1,1:nbl) = two
          mat(3,5,5,1,1:nbl) = two
          
          bc(1,1,1,1:nbl) = -one
          bc(1,2,2,1:nbl) = -one
          bc(1,3,3,1:nbl) = -one
          bc(1,4,4,1:nbl) = -one
          bc(1,5,5,1:nbl) = -one

          end if
        
          if (extrap.eq.2) then
          
!.... second-order extrapolation in the viscous layer

          mat(:,:,:,1,1:nbl) = zero
          bc(1:3,:,:,1:nbl)  = zero

          mat(2,1,1,1,1:nbl) = -one
          mat(2,2,2,1,1:nbl) = -one
          mat(2,3,3,1,1:nbl) = -one
          mat(2,4,4,1,1:nbl) = -one
          mat(2,5,5,1,1:nbl) = -one

          mat(3,1,1,1,1:nbl) = three
          mat(3,2,2,1,1:nbl) = three
          mat(3,3,3,1,1:nbl) = three
          mat(3,4,4,1,1:nbl) = three
          mat(3,5,5,1,1:nbl) = three
          
          bc(1,1,1,1:nbl) = -three
          bc(1,2,2,1:nbl) = -three
          bc(1,3,3,1:nbl) = -three
          bc(1,4,4,1:nbl) = -three
          bc(1,5,5,1:nbl) = -three

          bc(2,1,1,1:nbl) = one
          bc(2,2,2,1:nbl) = one
          bc(2,3,3,1:nbl) = one
          bc(2,4,4,1:nbl) = one
          bc(2,5,5,1:nbl) = one

          end if                ! extrapolation type

        end if
!=============================================================================!
!       R i g h t  B o u n d a r y
!=============================================================================!

!.... Right Riemann invariant boundary condition (first-order)
        
        if (right.eq.0 .and. Ma.lt.one) then

          mat(:,:,:,nx,:) = zero
          bc(12:14,:,:,:) = zero
                
          call ReimannLHS( ny, vl(:,nx,:), vl(:,nx-1,:), vl(:,nx-2,:), &
                           bnr, x(nx,:), y(nx,:), mat(2,:,:,nx,:), &
                           mat(1,:,:,nx,:), bc(14,:,:,:), &
                           rhobr, ubr, vbr, wbr, tbr, pbr, cbr )
        end if
          
!.... Zero'th order extrapolation in the viscous layer

          if (extrap.eq.0) then
            mat(:,:,:,nx,1:nbl) = zero
            bc(12:14,:,:,1:nbl) = zero
          
            mat(2,1,1,nx,1:nbl) = -one
            mat(2,2,2,nx,1:nbl) = -one
            mat(2,3,3,nx,1:nbl) = -one
            mat(2,4,4,nx,1:nbl) = -one
            mat(2,5,5,nx,1:nbl) = -one

            mat(1,1,1,nx,1:nbl) = one
            mat(1,2,2,nx,1:nbl) = one
            mat(1,3,3,nx,1:nbl) = one
            mat(1,4,4,nx,1:nbl) = one
            mat(1,5,5,nx,1:nbl) = one
          end if
          
!.... First-order extrapolation in the viscous layer

          if (extrap.eq.1) then
            mat(:,:,:,nx,1:nbl) = zero
            bc(12:14,:,:,1:nbl) = zero

            mat(2,1,1,nx,1:nbl) = -one
            mat(2,2,2,nx,1:nbl) = -one
            mat(2,3,3,nx,1:nbl) = -one
            mat(2,4,4,nx,1:nbl) = -one
            mat(2,5,5,nx,1:nbl) = -one

            mat(1,1,1,nx,1:nbl) = two
            mat(1,2,2,nx,1:nbl) = two
            mat(1,3,3,nx,1:nbl) = two
            mat(1,4,4,nx,1:nbl) = two
            mat(1,5,5,nx,1:nbl) = two
          
            bc(14,1,1,1:nbl) = -one
            bc(14,2,2,1:nbl) = -one
            bc(14,3,3,1:nbl) = -one
            bc(14,4,4,1:nbl) = -one
            bc(14,5,5,1:nbl) = -one
          end if                ! extrapolation type

!.... hold the IC

          if (right.eq.8) then
            mat(:,:,:,nx,:) = zero
            bc(12:14,:,:,:) = zero
            mat(2,1,1,nx,:) = -one
            mat(2,2,2,nx,:) = -one
            mat(2,3,3,nx,:) = -one
            mat(2,4,4,nx,:) = -one
            mat(2,5,5,nx,:) = -one
          end if

        end if                  ! xper
        
        end if                  ! linear
        
        return
        end
