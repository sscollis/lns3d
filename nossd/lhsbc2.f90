!=============================================================================!
        subroutine lhsbc2( mat, bc, vl )
!  
!  Correct the LHS for boundary conditions in the \eta direction
!
!  This version supports the 2nd order LHS
!
!  Revised: 4-16-96
!=============================================================================!
        use global
        use stencil
        use pot
        implicit none
        
        real :: mat(ny,nx,ndof,ndof,3), bc(nx,ndof,ndof,14)
        real :: vl(ny,nx,ndof)

        real :: detainv, detasinv
        
        integer :: iv, idof, jdof
        
!.... local variables for Riemann invariants

        real :: m1l(nx), m2l(nx), R1(nx), R2(nx), fact(nx), vint(nx,ndof)
        real :: rho(nx), p(nx), fact1(nx), fact2(nx), vn(nx), vt(nx), pnorm(nx)
        real :: ub(nx), vb(nx), tb(nx), rhob(nx), pb(nx), cb(nx)
        real :: xil(nx), etal(nx)
        real :: b, rbeta
        
        integer :: i, j, ij
!=============================================================================!
        detainv  = one / deta
        detasinv = one / deta**2
!=============================================================================!
!       L I N E A R   B O U N D A R Y   C O N D I T I O N S
!=============================================================================!
        if (linear.eq.1) then

        if (.not. yper) then

        if (Navier) then

!.... wall boundary condition

        if (wall.eq.1) then
          mat(1,:,1,:,:) = zero
          bc(:,1,:,1:3)  = zero
          mat(1,:,1,1,2) = -gc1
          mat(1,:,1,1,3) = -gc2
          bc(:,1,1,1)    = -gc3
          bc(:,1,1,2)    = -gc4
          bc(:,1,1,3)    = -gc5
        end if

!.... no-slip

        mat(1,:,2:4,:,:) = zero
        bc(:,2:4,:,1:3)  = zero
        mat(1,:,2,2,2)   = one
        mat(1,:,3,3,2)   = one
        mat(1,:,4,4,2)   = one

!.... isothermal

        if (wallt.eq.0) then
          mat(1,:,5,:,:) = zero
          bc(:,5,:,1:3) = zero
          mat(1,:,5,5,2) = one
        end if

!.... adiabatic

        if (wallt.eq.1) then
          mat(1,:,5,:,:) = zero
          bc(:,5,:,1:3)  = zero
          mat(1,:,5,5,2) = -gc1
          mat(1,:,5,5,3) = -gc2
          bc(:,5,5,1)    = -gc3
          bc(:,5,5,2)    = -gc4
          bc(:,5,5,3)    = -gc5
        end if

        else                    ! inviscid wall
          do i = 1, nx
            mat(1,i,2,:,:) = bnb(i,2) * mat(1,i,2,:,:) - &
                             bnb(i,1) * mat(1,i,3,:,:)
          end do
          mat(1,:,3,:,:) = zero
          mat(1,:,3,2,2) = bnb(:,1)
          mat(1,:,3,3,2) = bnb(:,2)
        end if

!.... freestream zero disturbance boundary condition

        if (top.eq.0) then
          mat(ny,:,:,:,:) = zero
          bc(:,:,:,12:14) = zero
          mat(ny,:,1,1,2) = one
          mat(ny,:,2,2,2) = one
          mat(ny,:,3,3,2) = one
          mat(ny,:,4,4,2) = one
          mat(ny,:,5,5,2) = one
        end if
        
        end if          ! yper

!=============================================================================!

        if (.not. xper) then

!.... Left boundary

          if (left.eq.0) then           ! zero disturbance BC
            mat(:,1,:,:,:) = zero
            bc(1,:,:,:)    = zero
            mat(:,1,1,1,2) = one
            mat(:,1,2,2,2) = one
            mat(:,1,3,3,2) = one
            mat(:,1,4,4,2) = one
            mat(:,1,5,5,2) = one
          else if (left.eq.1) then      ! nonreflecting BC
            mat(:,1,1,:,:) = zero
            bc(1,1,:,:)    = zero
            mat(:,1,1,1,2) = one
          end if

!.... Right boundary

          if (right.eq.0) then          ! zero disturbance
            mat(:,nx,:,:,:) = zero
            mat(:,nx,1,1,2) = one
            mat(:,nx,2,2,2) = one
            mat(:,nx,3,3,2) = one
            mat(:,nx,4,4,2) = one
            mat(:,nx,5,5,2) = one
          else if (right.eq.1) then     ! nonreflecting BC
            mat(:,nx,1,:,:) = zero
            bc(nx,1,:,:)    = zero
            mat(:,nx,1,1,2) = one
          end if

        end if          ! xper

!=============================================================================!
!       N O N L I N E A R   B O U N D A R Y   C O N D I T I O N S
!=============================================================================!
        else            ! nonlinear

        if (.not. yper) then

        if (Navier) then

!.... wall boundary condition

        if (wall.eq.1) then
          mat(1,:,1,:,:) = zero
          bc(:,1,:,1:3)  = zero
          mat(1,:,1,1,2) = -gc1
          mat(1,:,1,1,3) = -gc2
          bc(:,1,1,1)    = -gc3
          bc(:,1,1,2)    = -gc4
          bc(:,1,1,3)    = -gc5
        end if

!.... linear extrapolation of rho at the wall

        if (wall.eq.2) then
          mat(1,:,1,:,:) =  zero
          bc(:,1,:,1:3) =  zero
          mat(1,:,1,1,2) = -one
          mat(1,:,1,1,3) =  two
          bc(:,1,1,1)    = -one
          bc(:,1,1,2)    =  zero
          bc(:,1,1,3)    =  zero
        end if

!.... no-slip

        mat(1,:,2:4,:,:) = zero
        bc(:,2:4,:,1:3)  = zero
        mat(1,:,2,2,2)   = one
        mat(1,:,3,3,2)   = one
        mat(1,:,4,4,2)   = one

!.... isothermal

        if (wallt.eq.0) then
          mat(1,:,5,:,:) = zero
          bc(:,5,:,1:3)  = zero
          mat(1,:,5,5,2) = -one
        end if

!.... adiabatic

        if (wallt.eq.1) then
          mat(1,:,5,:,:) = zero
          bc(:,5,:,1:3)  = zero
          mat(1,:,5,5,2) = -gc1
          mat(1,:,5,5,3) = -gc2
          bc(:,5,5,1)    = -gc3
          bc(:,5,5,2)    = -gc4
          bc(:,5,5,3)    = -gc5
        end if

!.... partial tangent for normal momentum equation

        if (wall.eq.2) then
          call wallbc(vl, p, pnorm)
          mat(1,:,1,:,:) = zero
          bc(:,1,:,1:3)  = zero
          mat(1,:,1,1,2) = -gc1 * detainv
          mat(1,:,1,5,2) = -gamma * Ma**2 * Pnorm / vl(1,:,5)**2
          mat(1,:,1,1,3) = -gc2 * detainv
          bc(:,1,1,1)    = -gc3 * detainv
          bc(:,1,1,2)    = -gc4 * detainv
          bc(:,1,1,3)    = -gc5 * detainv
        end if
        
        else                    ! inviscid wall
          do i = 1, nx
            mat(1,i,2,:,:) = bnb(i,2) * mat(1,i,2,:,:) - &
                             bnb(i,1) * mat(1,i,3,:,:)
          end do
          mat(1,:,3,:,:) = zero
          mat(1,:,3,2,2) = bnb(:,1)
          mat(1,:,3,3,2) = bnb(:,2)
        end if
!=============================================================================!
!       T o p
!=============================================================================!
        if (Ma.lt.one) then

        if (top.eq.0) then
                
          mat(ny,:,:,:,:) = zero
          bc(:,:,:,12:14) = zero

          call ReimannLHS( nx, vl(ny,:,:), vl(ny-1,:,:), vl(ny-2,:,:), &
                           bnt, x(ny,:), y(ny,:), mat(ny,:,:,:,2), &
                           mat(ny,:,:,:,1), bc(:,:,:,14) )
                
        end if

!.... freestream hard boundary condition

        else
        
          mat(ny,:,:,:,:) = zero
          bc(:,:,:,12:14) = zero
          mat(ny,:,1,1,2) = one
          mat(ny,:,2,2,2) = one
          mat(ny,:,3,3,2) = one
          mat(ny,:,4,4,2) = one
          mat(ny,:,5,5,2) = one
        
        end if          ! Mach

        end if          ! yper
        
!=============================================================================!

        if (.not. xper) then

        if (Ma.lt.one) then

!.... Riemann invariants or extrapolation on the x-boundaries

          if (left.eq.0) then
            mat(:,1,:,:,:) = zero
            bc(1,:,:,:)    = zero
            mat(:,1,1,1,2) = one
            mat(:,1,2,2,2) = one
            mat(:,1,3,3,2) = one
            mat(:,1,4,4,2) = one
            mat(:,1,5,5,2) = one
          end if

!.... symmetry plan

          if (left.eq.2) then
            mat(:,1,:,:,:) = zero
            bc(1,:,:,:)    = zero
            mat(:,1,1,1,2) = one
            mat(:,1,2,2,2) = one
            mat(:,1,3,3,2) = one
            mat(:,1,4,4,2) = one
            mat(:,1,5,5,2) = one
          end if

          if (left.eq.7) then           ! symmetry plane
            mat(:,1,3,:,:) = zero
            bc(1,3,:,:)    = zero
            mat(:,1,3,3,2) = one
          end if

          if (right.eq.0 .or. right.eq.8) then
            mat(:,nx,1:5,:,:) = zero
            bc(nx,:,:,:)      = zero
            mat(:,nx,1,1,2)   = one
            mat(:,nx,2,2,2)   = one
            mat(:,nx,3,3,2)   = one
            mat(:,nx,4,4,2)   = one
            mat(:,nx,5,5,2)   = one
          end if

        else
        
          mat(nbl:ny,1,1:5,:,:) = zero
          bc(1,:,:,12:14)       = zero
          mat(nbl:ny,1,1,1,2)   = one
          mat(nbl:ny,1,2,2,2)   = one
          mat(nbl:ny,1,3,3,2)   = one
          mat(nbl:ny,1,4,4,2)   = one
          mat(nbl:ny,1,5,5,2)   = one

          mat(nbl:ny,nx,1:5,:,:) = zero
          bc(nx,:,:,12:14)       = zero
          mat(nbl:ny,nx,1,1,2)   = one
          mat(nbl:ny,nx,2,2,2)   = one
          mat(nbl:ny,nx,3,3,2)   = one
          mat(nbl:ny,nx,4,4,2)   = one
          mat(nbl:ny,nx,5,5,2)   = one

        end if

        end if          ! xper
        
        end if          ! linear

        return
        end
