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
        
        real :: mat(3,ndof,ndof,nx,ny), bc(14,ndof,ndof,nx)
        real :: vl(ndof,nx,ny)

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
          mat(:,1,:,:,1) = zero
          bc(1:3,1,:,:)  = zero
          mat(2,1,1,:,1) = -gc1
          mat(3,1,1,:,1) = -gc2
          bc(1,1,1,:)    = -gc3
          bc(2,1,1,:)    = -gc4
          bc(3,1,1,:)    = -gc5
        end if

!.... no-slip

        mat(:,2:4,:,:,1) = zero
        bc(1:3,2:4,:,:)  = zero
        mat(2,2,2,:,1)   = one
        mat(2,3,3,:,1)   = one
        mat(2,4,4,:,1)   = one

!.... isothermal

        if (wallt.eq.0) then
          mat(:,5,:,:,1) = zero
          bc(1:3,5,:,:) = zero
          mat(2,5,5,:,1) = one
        end if

!.... adiabatic

        if (wallt.eq.1) then
          mat(:,5,:,:,1) = zero
          bc(1:3,5,:,:)  = zero
          mat(2,5,5,:,1) = -gc1
          mat(3,5,5,:,1) = -gc2
          bc(1,5,5,:)    = -gc3
          bc(2,5,5,:)    = -gc4
          bc(3,5,5,:)    = -gc5
        end if

        else                    ! inviscid wall
          do i = 1, nx
            mat(:,2,:,i,1) = bnb(i,2) * mat(:,2,:,i,1) - &
                             bnb(i,1) * mat(:,3,:,i,1)
          end do
          mat(:,3,:,:,1) = zero
          mat(2,3,2,:,1) = bnb(:,1)
          mat(2,3,3,:,1) = bnb(:,2)
        end if

!.... freestream zero disturbance boundary condition

        if (top.eq.0) then
          mat(:,:,:,:,ny) = zero
          bc(12:14,:,:,:) = zero
          mat(2,1,1,:,ny) = one
          mat(2,2,2,:,ny) = one
          mat(2,3,3,:,ny) = one
          mat(2,4,4,:,ny) = one
          mat(2,5,5,:,ny) = one
        end if
        
        end if          ! yper

!=============================================================================!

        if (.not. xper) then

!.... Left boundary

          if (left.eq.0) then           ! zero disturbance BC
            mat(:,:,:,1,:) = zero
            bc(:,:,:,1)    = zero
            mat(2,1,1,1,:) = one
            mat(2,2,2,1,:) = one
            mat(2,3,3,1,:) = one
            mat(2,4,4,1,:) = one
            mat(2,5,5,1,:) = one
          else if (left.eq.1) then      ! nonreflecting BC
            mat(:,1,:,1,:) = zero
            bc(:,1,:,1)    = zero
            mat(2,1,1,1,:) = one
          end if

!.... Right boundary

          if (right.eq.0) then          ! zero disturbance
            mat(:,:,:,nx,:) = zero
            mat(2,1,1,nx,:) = one
            mat(2,2,2,nx,:) = one
            mat(2,3,3,nx,:) = one
            mat(2,4,4,nx,:) = one
            mat(2,5,5,nx,:) = one
          else if (right.eq.1) then     ! nonreflecting BC
            mat(:,1,:,nx,:) = zero
            bc(:,1,:,nx)    = zero
            mat(2,1,1,nx,:) = one
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
          mat(:,1,:,:,1) = zero
          bc(1:3,1,:,:)  = zero
          mat(2,1,1,:,1) = -gc1
          mat(3,1,1,:,1) = -gc2
          bc(1,1,1,:)    = -gc3
          bc(2,1,1,:)    = -gc4
          bc(3,1,1,:)    = -gc5
        end if

!.... linear extrapolation of rho at the wall

        if (wall.eq.2) then
          mat(:,1,:,:,1) =  zero
          bc(1:3,1,:,:) =  zero
          mat(2,1,1,:,1) = -one
          mat(3,1,1,:,1) =  two
          bc(1,1,1,:)    = -one
          bc(2,1,1,:)    =  zero
          bc(3,1,1,:)    =  zero
        end if

!.... Third-order extrapolation of rho at the wall

        if (wall.eq.3) then
          mat(:,1,:,:,1) =  zero
          mat(2,1,1,:,1) = -one
          mat(3,1,1,:,1) =  4.0
          bc(1,1,1,:)    = -6.0
          bc(2,1,1,:)    =  4.0
          bc(3,1,1,:)    = -one
        end if

!.... no-slip

        mat(:,2:4,:,:,1) = zero
        bc(1:3,2:4,:,:)  = zero
        mat(2,2,2,:,1)   = one
        mat(2,3,3,:,1)   = one
        mat(2,4,4,:,1)   = one

!.... isothermal

        if (wallt.eq.0) then
          mat(:,5,:,:,1) = zero
          bc(1:3,5,:,:)  = zero
          mat(2,5,5,:,1) = -one
        end if

!.... adiabatic

        if (wallt.eq.1) then
          mat(:,5,:,:,1) = zero
          bc(1:3,5,:,:)  = zero
          mat(2,5,5,:,1) = -gc1
          mat(3,5,5,:,1) = -gc2
          bc(1,5,5,:)    = -gc3
          bc(2,5,5,:)    = -gc4
          bc(3,5,5,:)    = -gc5
        end if

!.... partial tangent for normal momentum equation

        if (wall.eq.2) then
          call wallbc(vl, p, pnorm)
          mat(:,1,:,:,1) = zero
          bc(1:3,1,:,:)  = zero
          mat(2,1,1,:,1) = -gc1 * detainv
          mat(2,1,5,:,1) = -gamma * Ma**2 * Pnorm / vl(5,:,1)**2
          mat(3,1,1,:,1) = -gc2 * detainv
          bc(1,1,1,:)    = -gc3 * detainv
          bc(2,1,1,:)    = -gc4 * detainv
          bc(3,1,1,:)    = -gc5 * detainv
        end if
        
        else                    ! inviscid wall

          do i = 1, nx
            mat(:,2,:,i,1) = bnb(i,2) * mat(:,2,:,i,1) - &
                             bnb(i,1) * mat(:,3,:,i,1)
          end do
          mat(:,3,:,:,1) = zero
          mat(2,3,2,:,1) = bnb(:,1)
          mat(2,3,3,:,1) = bnb(:,2)

        end if
!=============================================================================!
!       T o p
!=============================================================================!
        if (Ma.lt.one) then

        if (top.eq.0) then
                
          mat(:,:,:,:,ny) = zero
          bc(12:14,:,:,:) = zero

          call ReimannLHS( nx, vl(:,:,ny), vl(:,:,ny-1), vl(:,:,ny-2), &
                           bnt, x(:,ny), y(:,ny), mat(2,:,:,:,ny), &
                           mat(1,:,:,:,ny), bc(14,:,:,:), &
                           rhobt, ubt, vbt, wbt, tbt, pbt, cbt )
                
        end if

!.... freestream hard boundary condition

        else
        
          mat(:,:,:,:,ny) = zero
          bc(12:14,:,:,:) = zero
          mat(2,1,1,:,ny) = one
          mat(2,2,2,:,ny) = one
          mat(2,3,3,:,ny) = one
          mat(2,4,4,:,ny) = one
          mat(2,5,5,:,ny) = one
        
        end if          ! Mach

        end if          ! yper
        
!=============================================================================!

        if (.not. xper) then

        if (Ma.lt.one) then

!.... Riemann invariants or extrapolation on the x-boundaries

          if (left.eq.0) then
            mat(:,:,:,1,:) = zero
            bc(:,:,:,1)    = zero
            mat(2,1,1,1,:) = one
            mat(2,2,2,1,:) = one
            mat(2,3,3,1,:) = one
            mat(2,4,4,1,:) = one
            mat(2,5,5,1,:) = one
          end if

!.... symmetry plan

          if (left.eq.2) then
            mat(:,:,:,1,:) = zero
            bc(:,:,:,1)    = zero
            mat(2,1,1,1,:) = one
            mat(2,2,2,1,:) = one
            mat(2,3,3,1,:) = one
            mat(2,4,4,1,:) = one
            mat(2,5,5,1,:) = one
          end if

          if (left.eq.7) then           ! symmetry plane
            mat(:,3,:,1,:) = zero
            bc(:,3,:,1)    = zero
            mat(2,3,3,1,:) = one
          end if

          if (right.eq.0 .or. right.eq.8) then
            mat(:,1:5,:,nx,:) = zero
            bc(:,:,:,nx)      = zero
            mat(2,1,1,nx,:)   = one
            mat(2,2,2,nx,:)   = one
            mat(2,3,3,nx,:)   = one
            mat(2,4,4,nx,:)   = one
            mat(2,5,5,nx,:)   = one
          end if

        else
        
          mat(:,1:5,:,1,nbl:ny) = zero
          bc(12:14,:,:,1)       = zero
          mat(2,1,1,1,nbl:ny)   = one
          mat(2,2,2,1,nbl:ny)   = one
          mat(2,3,3,1,nbl:ny)   = one
          mat(2,4,4,1,nbl:ny)   = one
          mat(2,5,5,1,nbl:ny)   = one

          mat(:,1:5,:,nx,nbl:ny) = zero
          bc(12:14,:,:,nx)       = zero
          mat(2,1,1,nx,nbl:ny)   = one
          mat(2,2,2,nx,nbl:ny)   = one
          mat(2,3,3,nx,nbl:ny)   = one
          mat(2,4,4,nx,nbl:ny)   = one
          mat(2,5,5,nx,nbl:ny)   = one

        end if

        end if          ! xper
        
        end if          ! linear

        return
        end
