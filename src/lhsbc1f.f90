!=============================================================================!
        subroutine lhsbc1f( mat, vl, vml)
!  
!  Correct the LHS for boundary conditions in the \xi direction.
!
!  This version supports the 4th order LHS. 
!
!  Revised: 10-16-95
!
!=============================================================================!
        use global
        use stencil
        use pot
        implicit none
        
        real :: mat(5,ndof,ndof,nx,ny)
        real :: vl(ndof,nx,ny), vml(ndof,nx,ny)

        real :: matl(5,2)

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
          mat(3,1,1,:,1) = one
        end if

!.... extrapolation of rho at the wall

        if (wall.eq.2) then
          mat(:,1,:,:,1) = zero
          mat(3,1,1,:,1) = one
        end if

!.... no-slip

        mat(:,2:4,:,:,1) = zero
        mat(3,2,2,:,1)   = one
        mat(3,3,3,:,1)   = one
        mat(3,4,4,:,1)   = one

!.... wall temperature boundary condition

        if (wallt.eq.0 .or. wallt.eq.1) then
          mat(:,5,:,:,1) = zero
          mat(3,5,5,:,1) = one
        end if
        
        else                    ! inviscid wall

          mat(:,1,:,:,1) = zero
          mat(3,1,1,:,1) = one
          
          do i = 1, nx
            if ( abs(bnb(i,2)) .ge. abs(bnb(i,1)) .and. .false.) then
              mat(:,2,:,i,1) = bnb(i,2) * mat(:,2,:,i,1) - &
                               bnb(i,1) * mat(:,3,:,i,1)
              mat(:,3,:,i,1) = zero
              mat(3,3,2,i,1) = bnb(i,1)
              mat(3,3,3,i,1) = bnb(i,2)
            else
              mat(:,3,:,i,1) = bnb(i,2) * mat(:,2,:,i,1) - &
                               bnb(i,1) * mat(:,3,:,i,1)
              mat(:,2,:,i,1) = zero
              mat(3,2,2,i,1) = bnb(i,1)
              mat(3,2,3,i,1) = bnb(i,2)
            end if
          end do
        end if                  ! Navier

!.... freestream zero disturbance boundary condition

        if (top.eq.0.or.top.eq.1) then
          mat(:,:,:,:,ny) = zero
          mat(3,1,1,:,ny) = one
          mat(3,2,2,:,ny) = one
          mat(3,3,3,:,ny) = one
          mat(3,4,4,:,ny) = one
          mat(3,5,5,:,ny) = one
        end if
        
        end if          ! yper

!=============================================================================!
        
        if (.not. xper) then

!.... left boundary

          if (left.eq.0.or.left.eq.4) then ! zero disturbance

            mat(:,:,:,1,:) = zero
            mat(3,1,1,1,:) = one
            mat(3,2,2,1,:) = one
            mat(3,3,3,1,:) = one
            mat(3,4,4,1,:) = one
            mat(3,5,5,1,:) = one

          else if (left.eq.1) then      ! nonreflecting BC

            mat(:,1,:,1,nbl+1:ny) = zero
            mat(3,1,1,1,nbl+1:ny) = -one
            mat(3,1,2,1,nbl+1:ny) = gamma * Ma * vml(1,1,nbl+1:ny) * &
                                    sqrt(vml(5,1,nbl+1:ny)) / vml(5,1,nbl+1:ny)
            mat(3,1,5,1,nbl+1:ny) = -vml(1,1,nbl+1:ny) / vml(5,1,nbl+1:ny)

            if (.false.) then
              mat(:,:,:,1,1:nbl) = zero

              mat(3,1,1,1,1:nbl) = -one
              mat(3,2,2,1,1:nbl) = -one
              mat(3,3,3,1,1:nbl) = -one
              mat(3,4,4,1,1:nbl) = -one
              mat(3,5,5,1,1:nbl) = -one

              mat(4,1,1,1,1:nbl) = one
              mat(4,2,2,1,1:nbl) = one
              mat(4,3,3,1,1:nbl) = one
              mat(4,4,4,1,1:nbl) = one
              mat(4,5,5,1,1:nbl) = one
            end if

          else if (left.eq.2) then      ! symmetry plane

!           mat(:,1,:,1,:) = zero
!           mat(3,1,1,1,:) = -gc1
!           mat(4,1,1,1,:) = -gc2
!           mat(5,1,1,1,:) = -gc3
!           mat(1,1,1,1,:) = -gc4
!           mat(2,1,1,1,:) = -gc5
  
            mat(:,2,:,1,:) = zero
            mat(3,2,2,1,:) = -gc1
            mat(4,2,2,1,:) = -gc2
            mat(5,2,2,1,:) = -gc3
            mat(1,2,2,1,:) = -gc4
            mat(2,2,2,1,:) = -gc5

            mat(:,3,:,1,:) = zero
            mat(3,3,3,1,:) = one

            mat(:,4,:,1,:) = zero
            mat(3,4,4,1,:) = -gc1
            mat(4,4,4,1,:) = -gc2
            mat(5,4,4,1,:) = -gc3
            mat(1,4,4,1,:) = -gc4
            mat(2,4,4,1,:) = -gc5

            mat(:,5,:,1,:) = zero
            mat(3,5,5,1,:) = -gc1
            mat(4,5,5,1,:) = -gc2
            mat(5,5,5,1,:) = -gc3
            mat(1,5,5,1,:) = -gc4
            mat(2,5,5,1,:) = -gc5

          else if (left.eq.7) then      ! symmetry plane

!           mat(:,3,:,1,:) = zero
!           mat(3,3,3,1,:) = one

          endif        ! left

!.... Right boundary

          if (right.eq.0) then          ! zero disturbance

            mat(:,:,:,nx,:) = zero
            mat(3,1,1,nx,:) = one
            mat(3,2,2,nx,:) = one
            mat(3,3,3,nx,:) = one
            mat(3,4,4,nx,:) = one
            mat(3,5,5,nx,:) = one

          else if (right.eq.1) then     ! nonreflecting

            mat(:,1,:,1,nbl+1:ny) = zero
            mat(3,1,1,1,nbl+1:ny) = -one
            mat(3,1,2,1,nbl+1:ny) = gamma * Ma * vml(1,1,nbl+1:ny) * &
                                    sqrt(vml(5,1,nbl+1:ny)) / vml(5,1,nbl+1:ny)
            mat(3,1,5,1,nbl+1:ny) = -vml(1,1,nbl+1:ny) / vml(5,1,nbl+1:ny)

            mat(:,1,:,nx,nbl+1:ny) = zero
            mat(3,1,1,nx,nbl+1:ny) = -one
            mat(3,1,2,nx,nbl+1:ny) = gamma * Ma * vml(1,nx,nbl+1:ny) * &
                                     sqrt(vml(5,nx,nbl+1:ny)) / &
                                     vml(5,nx,nbl+1:ny)
            mat(3,1,5,nx,nbl+1:ny) = -vml(1,nx,nbl+1:ny) / vml(5,nx,nbl+1:ny)

            if (.false.) then
              mat(:,:,:,nx,1:nbl) = zero

              mat(3,1,1,nx,1:nbl) = -one
              mat(3,2,2,nx,1:nbl) = -one
              mat(3,3,3,nx,1:nbl) = -one
              mat(3,4,4,nx,1:nbl) = -one
              mat(3,5,5,nx,1:nbl) = -one

              mat(2,1,1,nx,1:nbl) = one
              mat(2,2,2,nx,1:nbl) = one
              mat(2,3,3,nx,1:nbl) = one
              mat(2,4,4,nx,1:nbl) = one
              mat(2,5,5,nx,1:nbl) = one
            end if

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
          mat(3,1,1,:,1) = one
        end if

!.... Lele-Poinsot

        if (wall.eq.4) then

          ! LHS must be modified

        end if

!.... no-slip

        mat(:,2:4,:,:,1) = zero
        mat(3,2,2,:,1)   = one
        mat(3,3,3,:,1)   = one
        mat(3,4,4,:,1)   = one

!.... adiabatic & isothermal

        if (wallt.eq.0 .or. wallt.eq.1) then
          mat(:,5,:,:,1) = zero
          mat(3,5,5,:,1) = one
        end if

        else      ! inviscid:  rotate to body normal coordinates

          i = 1

          matl(1,1) = mat(3,1,2,i,1) * bnb(i,2) - mat(3,1,3,i,1) * bnb(i,1)
          matl(1,2) = mat(3,1,2,i,1) * bnb(i,1) + mat(3,1,3,i,1) * bnb(i,2)
          matl(2,1) = mat(3,2,2,i,1) * bnb(i,2) - mat(3,2,3,i,1) * bnb(i,1)
          matl(2,2) = mat(3,2,2,i,1) * bnb(i,1) + mat(3,2,3,i,1) * bnb(i,2)
          matl(3,1) = mat(3,3,2,i,1) * bnb(i,2) - mat(3,3,3,i,1) * bnb(i,1)
          matl(3,2) = mat(3,3,2,i,1) * bnb(i,1) + mat(3,3,3,i,1) * bnb(i,2)
          matl(4,1) = mat(3,4,2,i,1) * bnb(i,2) - mat(3,4,3,i,1) * bnb(i,1)
          matl(4,2) = mat(3,4,2,i,1) * bnb(i,1) + mat(3,4,3,i,1) * bnb(i,2)
          matl(5,1) = mat(3,5,2,i,1) * bnb(i,2) - mat(3,5,3,i,1) * bnb(i,1)
          matl(5,2) = mat(3,5,2,i,1) * bnb(i,1) + mat(3,5,3,i,1) * bnb(i,2)
          mat(3,:,2:3,i,1) = matl
          
          matl(1,1) = mat(4,1,2,i,1) * bnb(i+1,2) - mat(4,1,3,i,1) * bnb(i+1,1)
          matl(1,2) = mat(4,1,2,i,1) * bnb(i+1,1) + mat(4,1,3,i,1) * bnb(i+1,2)
          matl(2,1) = mat(4,2,2,i,1) * bnb(i+1,2) - mat(4,2,3,i,1) * bnb(i+1,1)
          matl(2,2) = mat(4,2,2,i,1) * bnb(i+1,1) + mat(4,2,3,i,1) * bnb(i+1,2)
          matl(3,1) = mat(4,3,2,i,1) * bnb(i+1,2) - mat(4,3,3,i,1) * bnb(i+1,1)
          matl(3,2) = mat(4,3,2,i,1) * bnb(i+1,1) + mat(4,3,3,i,1) * bnb(i+1,2)
          matl(4,1) = mat(4,4,2,i,1) * bnb(i+1,2) - mat(4,4,3,i,1) * bnb(i+1,1)
          matl(4,2) = mat(4,4,2,i,1) * bnb(i+1,1) + mat(4,4,3,i,1) * bnb(i+1,2)
          matl(5,1) = mat(4,5,2,i,1) * bnb(i+1,2) - mat(4,5,3,i,1) * bnb(i+1,1)
          matl(5,2) = mat(4,5,2,i,1) * bnb(i+1,1) + mat(4,5,3,i,1) * bnb(i+1,2)
          mat(4,:,2:3,i,1) = matl
          
          matl(1,1) = mat(5,1,2,i,1) * bnb(i+2,2) - mat(5,1,3,i,1) * bnb(i+2,1)
          matl(1,2) = mat(5,1,2,i,1) * bnb(i+2,1) + mat(5,1,3,i,1) * bnb(i+2,2)
          matl(2,1) = mat(5,2,2,i,1) * bnb(i+2,2) - mat(5,2,3,i,1) * bnb(i+2,1)
          matl(2,2) = mat(5,2,2,i,1) * bnb(i+2,1) + mat(5,2,3,i,1) * bnb(i+2,2)
          matl(3,1) = mat(5,3,2,i,1) * bnb(i+2,2) - mat(5,3,3,i,1) * bnb(i+2,1)
          matl(3,2) = mat(5,3,2,i,1) * bnb(i+2,1) + mat(5,3,3,i,1) * bnb(i+2,2)
          matl(4,1) = mat(5,4,2,i,1) * bnb(i+2,2) - mat(5,4,3,i,1) * bnb(i+2,1)
          matl(4,2) = mat(5,4,2,i,1) * bnb(i+2,1) + mat(5,4,3,i,1) * bnb(i+2,2)
          matl(5,1) = mat(5,5,2,i,1) * bnb(i+2,2) - mat(5,5,3,i,1) * bnb(i+2,1)
          matl(5,2) = mat(5,5,2,i,1) * bnb(i+2,1) + mat(5,5,3,i,1) * bnb(i+2,2)
          mat(5,:,2:3,i,1) = matl
          
          matl(1,1) = mat(1,1,2,i,1) * bnb(i+3,2) - mat(1,1,3,i,1) * bnb(i+3,1)
          matl(1,2) = mat(1,1,2,i,1) * bnb(i+3,1) + mat(1,1,3,i,1) * bnb(i+3,2)
          matl(2,1) = mat(1,2,2,i,1) * bnb(i+3,2) - mat(1,2,3,i,1) * bnb(i+3,1)
          matl(2,2) = mat(1,2,2,i,1) * bnb(i+3,1) + mat(1,2,3,i,1) * bnb(i+3,2)
          matl(3,1) = mat(1,3,2,i,1) * bnb(i+3,2) - mat(1,3,3,i,1) * bnb(i+3,1)
          matl(3,2) = mat(1,3,2,i,1) * bnb(i+3,1) + mat(1,3,3,i,1) * bnb(i+3,2)
          matl(4,1) = mat(1,4,2,i,1) * bnb(i+3,2) - mat(1,4,3,i,1) * bnb(i+3,1)
          matl(4,2) = mat(1,4,2,i,1) * bnb(i+3,1) + mat(1,4,3,i,1) * bnb(i+3,2)
          matl(5,1) = mat(1,5,2,i,1) * bnb(i+3,2) - mat(1,5,3,i,1) * bnb(i+3,1)
          matl(5,2) = mat(1,5,2,i,1) * bnb(i+3,1) + mat(1,5,3,i,1) * bnb(i+3,2)
          mat(1,:,2:3,i,1) = matl
          
          matl(1,1) = mat(2,1,2,i,1) * bnb(i+4,2) - mat(2,1,3,i,1) * bnb(i+4,1)
          matl(1,2) = mat(2,1,2,i,1) * bnb(i+4,1) + mat(2,1,3,i,1) * bnb(i+4,2)
          matl(2,1) = mat(2,2,2,i,1) * bnb(i+4,2) - mat(2,2,3,i,1) * bnb(i+4,1)
          matl(2,2) = mat(2,2,2,i,1) * bnb(i+4,1) + mat(2,2,3,i,1) * bnb(i+4,2)
          matl(3,1) = mat(2,3,2,i,1) * bnb(i+4,2) - mat(2,3,3,i,1) * bnb(i+4,1)
          matl(3,2) = mat(2,3,2,i,1) * bnb(i+4,1) + mat(2,3,3,i,1) * bnb(i+4,2)
          matl(4,1) = mat(2,4,2,i,1) * bnb(i+4,2) - mat(2,4,3,i,1) * bnb(i+4,1)
          matl(4,2) = mat(2,4,2,i,1) * bnb(i+4,1) + mat(2,4,3,i,1) * bnb(i+4,2)
          matl(5,1) = mat(2,5,2,i,1) * bnb(i+4,2) - mat(2,5,3,i,1) * bnb(i+4,1)
          matl(5,2) = mat(2,5,2,i,1) * bnb(i+4,1) + mat(2,5,3,i,1) * bnb(i+4,2)
          mat(2,:,2:3,i,1) = matl
          
          i = 2
          
          matl(1,1) = mat(2,1,2,i,1) * bnb(i-1,2) - mat(2,1,3,i,1) * bnb(i-1,1)
          matl(1,2) = mat(2,1,2,i,1) * bnb(i-1,1) + mat(2,1,3,i,1) * bnb(i-1,2)
          matl(2,1) = mat(2,2,2,i,1) * bnb(i-1,2) - mat(2,2,3,i,1) * bnb(i-1,1)
          matl(2,2) = mat(2,2,2,i,1) * bnb(i-1,1) + mat(2,2,3,i,1) * bnb(i-1,2)
          matl(3,1) = mat(2,3,2,i,1) * bnb(i-1,2) - mat(2,3,3,i,1) * bnb(i-1,1)
          matl(3,2) = mat(2,3,2,i,1) * bnb(i-1,1) + mat(2,3,3,i,1) * bnb(i-1,2)
          matl(4,1) = mat(2,4,2,i,1) * bnb(i-1,2) - mat(2,4,3,i,1) * bnb(i-1,1)
          matl(4,2) = mat(2,4,2,i,1) * bnb(i-1,1) + mat(2,4,3,i,1) * bnb(i-1,2)
          matl(5,1) = mat(2,5,2,i,1) * bnb(i-1,2) - mat(2,5,3,i,1) * bnb(i-1,1)
          matl(5,2) = mat(2,5,2,i,1) * bnb(i-1,1) + mat(2,5,3,i,1) * bnb(i-1,2)
          mat(2,:,2:3,i,1) = matl
          
          matl(1,1) = mat(3,1,2,i,1) * bnb(i,2) - mat(3,1,3,i,1) * bnb(i,1)
          matl(1,2) = mat(3,1,2,i,1) * bnb(i,1) + mat(3,1,3,i,1) * bnb(i,2)
          matl(2,1) = mat(3,2,2,i,1) * bnb(i,2) - mat(3,2,3,i,1) * bnb(i,1)
          matl(2,2) = mat(3,2,2,i,1) * bnb(i,1) + mat(3,2,3,i,1) * bnb(i,2)
          matl(3,1) = mat(3,3,2,i,1) * bnb(i,2) - mat(3,3,3,i,1) * bnb(i,1)
          matl(3,2) = mat(3,3,2,i,1) * bnb(i,1) + mat(3,3,3,i,1) * bnb(i,2)
          matl(4,1) = mat(3,4,2,i,1) * bnb(i,2) - mat(3,4,3,i,1) * bnb(i,1)
          matl(4,2) = mat(3,4,2,i,1) * bnb(i,1) + mat(3,4,3,i,1) * bnb(i,2)
          matl(5,1) = mat(3,5,2,i,1) * bnb(i,2) - mat(3,5,3,i,1) * bnb(i,1)
          matl(5,2) = mat(3,5,2,i,1) * bnb(i,1) + mat(3,5,3,i,1) * bnb(i,2)
          mat(3,:,2:3,i,1) = matl
          
          matl(1,1) = mat(4,1,2,i,1) * bnb(i+1,2) - mat(4,1,3,i,1) * bnb(i+1,1)
          matl(1,2) = mat(4,1,2,i,1) * bnb(i+1,1) + mat(4,1,3,i,1) * bnb(i+1,2)
          matl(2,1) = mat(4,2,2,i,1) * bnb(i+1,2) - mat(4,2,3,i,1) * bnb(i+1,1)
          matl(2,2) = mat(4,2,2,i,1) * bnb(i+1,1) + mat(4,2,3,i,1) * bnb(i+1,2)
          matl(3,1) = mat(4,3,2,i,1) * bnb(i+1,2) - mat(4,3,3,i,1) * bnb(i+1,1)
          matl(3,2) = mat(4,3,2,i,1) * bnb(i+1,1) + mat(4,3,3,i,1) * bnb(i+1,2)
          matl(4,1) = mat(4,4,2,i,1) * bnb(i+1,2) - mat(4,4,3,i,1) * bnb(i+1,1)
          matl(4,2) = mat(4,4,2,i,1) * bnb(i+1,1) + mat(4,4,3,i,1) * bnb(i+1,2)
          matl(5,1) = mat(4,5,2,i,1) * bnb(i+1,2) - mat(4,5,3,i,1) * bnb(i+1,1)
          matl(5,2) = mat(4,5,2,i,1) * bnb(i+1,1) + mat(4,5,3,i,1) * bnb(i+1,2)
          mat(4,:,2:3,i,1) = matl
          
          matl(1,1) = mat(5,1,2,i,1) * bnb(i+2,2) - mat(5,1,3,i,1) * bnb(i+2,1)
          matl(1,2) = mat(5,1,2,i,1) * bnb(i+2,1) + mat(5,1,3,i,1) * bnb(i+2,2)
          matl(2,1) = mat(5,2,2,i,1) * bnb(i+2,2) - mat(5,2,3,i,1) * bnb(i+2,1)
          matl(2,2) = mat(5,2,2,i,1) * bnb(i+2,1) + mat(5,2,3,i,1) * bnb(i+2,2)
          matl(3,1) = mat(5,3,2,i,1) * bnb(i+2,2) - mat(5,3,3,i,1) * bnb(i+2,1)
          matl(3,2) = mat(5,3,2,i,1) * bnb(i+2,1) + mat(5,3,3,i,1) * bnb(i+2,2)
          matl(4,1) = mat(5,4,2,i,1) * bnb(i+2,2) - mat(5,4,3,i,1) * bnb(i+2,1)
          matl(4,2) = mat(5,4,2,i,1) * bnb(i+2,1) + mat(5,4,3,i,1) * bnb(i+2,2)
          matl(5,1) = mat(5,5,2,i,1) * bnb(i+2,2) - mat(5,5,3,i,1) * bnb(i+2,1)
          matl(5,2) = mat(5,5,2,i,1) * bnb(i+2,1) + mat(5,5,3,i,1) * bnb(i+2,2)
          mat(5,:,2:3,i,1) = matl

          matl(1,1) = mat(1,1,2,i,1) * bnb(i+3,2) - mat(1,1,3,i,1) * bnb(i+3,1)
          matl(1,2) = mat(1,1,2,i,1) * bnb(i+3,1) + mat(1,1,3,i,1) * bnb(i+3,2)
          matl(2,1) = mat(1,2,2,i,1) * bnb(i+3,2) - mat(1,2,3,i,1) * bnb(i+3,1)
          matl(2,2) = mat(1,2,2,i,1) * bnb(i+3,1) + mat(1,2,3,i,1) * bnb(i+3,2)
          matl(3,1) = mat(1,3,2,i,1) * bnb(i+3,2) - mat(1,3,3,i,1) * bnb(i+3,1)
          matl(3,2) = mat(1,3,2,i,1) * bnb(i+3,1) + mat(1,3,3,i,1) * bnb(i+3,2)
          matl(4,1) = mat(1,4,2,i,1) * bnb(i+3,2) - mat(1,4,3,i,1) * bnb(i+3,1)
          matl(4,2) = mat(1,4,2,i,1) * bnb(i+3,1) + mat(1,4,3,i,1) * bnb(i+3,2)
          matl(5,1) = mat(1,5,2,i,1) * bnb(i+3,2) - mat(1,5,3,i,1) * bnb(i+3,1)
          matl(5,2) = mat(1,5,2,i,1) * bnb(i+3,1) + mat(1,5,3,i,1) * bnb(i+3,2)
          mat(1,:,2:3,i,1) = matl
          
          do i = 3, nx-2
            matl(1,1) = mat(1,1,2,i,1) * bnb(i-2,2) - mat(1,1,3,i,1)*bnb(i-2,1)
            matl(1,2) = mat(1,1,2,i,1) * bnb(i-2,1) + mat(1,1,3,i,1)*bnb(i-2,2)
            matl(2,1) = mat(1,2,2,i,1) * bnb(i-2,2) - mat(1,2,3,i,1)*bnb(i-2,1)
            matl(2,2) = mat(1,2,2,i,1) * bnb(i-2,1) + mat(1,2,3,i,1)*bnb(i-2,2)
            matl(3,1) = mat(1,3,2,i,1) * bnb(i-2,2) - mat(1,3,3,i,1)*bnb(i-2,1)
            matl(3,2) = mat(1,3,2,i,1) * bnb(i-2,1) + mat(1,3,3,i,1)*bnb(i-2,2)
            matl(4,1) = mat(1,4,2,i,1) * bnb(i-2,2) - mat(1,4,3,i,1)*bnb(i-2,1)
            matl(4,2) = mat(1,4,2,i,1) * bnb(i-2,1) + mat(1,4,3,i,1)*bnb(i-2,2)
            matl(5,1) = mat(1,5,2,i,1) * bnb(i-2,2) - mat(1,5,3,i,1)*bnb(i-2,1)
            matl(5,2) = mat(1,5,2,i,1) * bnb(i-2,1) + mat(1,5,3,i,1)*bnb(i-2,2)
            mat(1,:,2:3,i,1) = matl

            matl(1,1) = mat(2,1,2,i,1) * bnb(i-1,2) - mat(2,1,3,i,1)*bnb(i-1,1)
            matl(1,2) = mat(2,1,2,i,1) * bnb(i-1,1) + mat(2,1,3,i,1)*bnb(i-1,2)
            matl(2,1) = mat(2,2,2,i,1) * bnb(i-1,2) - mat(2,2,3,i,1)*bnb(i-1,1)
            matl(2,2) = mat(2,2,2,i,1) * bnb(i-1,1) + mat(2,2,3,i,1)*bnb(i-1,2)
            matl(3,1) = mat(2,3,2,i,1) * bnb(i-1,2) - mat(2,3,3,i,1)*bnb(i-1,1)
            matl(3,2) = mat(2,3,2,i,1) * bnb(i-1,1) + mat(2,3,3,i,1)*bnb(i-1,2)
            matl(4,1) = mat(2,4,2,i,1) * bnb(i-1,2) - mat(2,4,3,i,1)*bnb(i-1,1)
            matl(4,2) = mat(2,4,2,i,1) * bnb(i-1,1) + mat(2,4,3,i,1)*bnb(i-1,2)
            matl(5,1) = mat(2,5,2,i,1) * bnb(i-1,2) - mat(2,5,3,i,1)*bnb(i-1,1)
            matl(5,2) = mat(2,5,2,i,1) * bnb(i-1,1) + mat(2,5,3,i,1)*bnb(i-1,2)
            mat(2,:,2:3,i,1) = matl

            matl(1,1) = mat(3,1,2,i,1) * bnb(i,2) - mat(3,1,3,i,1) * bnb(i,1)
            matl(1,2) = mat(3,1,2,i,1) * bnb(i,1) + mat(3,1,3,i,1) * bnb(i,2)
            matl(2,1) = mat(3,2,2,i,1) * bnb(i,2) - mat(3,2,3,i,1) * bnb(i,1)
            matl(2,2) = mat(3,2,2,i,1) * bnb(i,1) + mat(3,2,3,i,1) * bnb(i,2)
            matl(3,1) = mat(3,3,2,i,1) * bnb(i,2) - mat(3,3,3,i,1) * bnb(i,1)
            matl(3,2) = mat(3,3,2,i,1) * bnb(i,1) + mat(3,3,3,i,1) * bnb(i,2)
            matl(4,1) = mat(3,4,2,i,1) * bnb(i,2) - mat(3,4,3,i,1) * bnb(i,1)
            matl(4,2) = mat(3,4,2,i,1) * bnb(i,1) + mat(3,4,3,i,1) * bnb(i,2)
            matl(5,1) = mat(3,5,2,i,1) * bnb(i,2) - mat(3,5,3,i,1) * bnb(i,1)
            matl(5,2) = mat(3,5,2,i,1) * bnb(i,1) + mat(3,5,3,i,1) * bnb(i,2)
            mat(3,:,2:3,i,1) = matl
            
            matl(1,1) = mat(4,1,2,i,1) * bnb(i+1,2) - mat(4,1,3,i,1)*bnb(i+1,1)
            matl(1,2) = mat(4,1,2,i,1) * bnb(i+1,1) + mat(4,1,3,i,1)*bnb(i+1,2)
            matl(2,1) = mat(4,2,2,i,1) * bnb(i+1,2) - mat(4,2,3,i,1)*bnb(i+1,1)
            matl(2,2) = mat(4,2,2,i,1) * bnb(i+1,1) + mat(4,2,3,i,1)*bnb(i+1,2)
            matl(3,1) = mat(4,3,2,i,1) * bnb(i+1,2) - mat(4,3,3,i,1)*bnb(i+1,1)
            matl(3,2) = mat(4,3,2,i,1) * bnb(i+1,1) + mat(4,3,3,i,1)*bnb(i+1,2)
            matl(4,1) = mat(4,4,2,i,1) * bnb(i+1,2) - mat(4,4,3,i,1)*bnb(i+1,1)
            matl(4,2) = mat(4,4,2,i,1) * bnb(i+1,1) + mat(4,4,3,i,1)*bnb(i+1,2)
            matl(5,1) = mat(4,5,2,i,1) * bnb(i+1,2) - mat(4,5,3,i,1)*bnb(i+1,1)
            matl(5,2) = mat(4,5,2,i,1) * bnb(i+1,1) + mat(4,5,3,i,1)*bnb(i+1,2)
            mat(4,:,2:3,i,1) = matl

            matl(1,1) = mat(5,1,2,i,1) * bnb(i+2,2) - mat(5,1,3,i,1)*bnb(i+2,1)
            matl(1,2) = mat(5,1,2,i,1) * bnb(i+2,1) + mat(5,1,3,i,1)*bnb(i+2,2)
            matl(2,1) = mat(5,2,2,i,1) * bnb(i+2,2) - mat(5,2,3,i,1)*bnb(i+2,1)
            matl(2,2) = mat(5,2,2,i,1) * bnb(i+2,1) + mat(5,2,3,i,1)*bnb(i+2,2)
            matl(3,1) = mat(5,3,2,i,1) * bnb(i+2,2) - mat(5,3,3,i,1)*bnb(i+2,1)
            matl(3,2) = mat(5,3,2,i,1) * bnb(i+2,1) + mat(5,3,3,i,1)*bnb(i+2,2)
            matl(4,1) = mat(5,4,2,i,1) * bnb(i+2,2) - mat(5,4,3,i,1)*bnb(i+2,1)
            matl(4,2) = mat(5,4,2,i,1) * bnb(i+2,1) + mat(5,4,3,i,1)*bnb(i+2,2)
            matl(5,1) = mat(5,5,2,i,1) * bnb(i+2,2) - mat(5,5,3,i,1)*bnb(i+2,1)
            matl(5,2) = mat(5,5,2,i,1) * bnb(i+2,1) + mat(5,5,3,i,1)*bnb(i+2,2)
            mat(5,:,2:3,i,1) = matl

          end do
            
          i = nx - 1
          
          matl(1,1) = mat(5,1,2,i,1) * bnb(i-3,2) - mat(5,1,3,i,1) * bnb(i-3,1)
          matl(1,2) = mat(5,1,2,i,1) * bnb(i-3,1) + mat(5,1,3,i,1) * bnb(i-3,2)
          matl(2,1) = mat(5,2,2,i,1) * bnb(i-3,2) - mat(5,2,3,i,1) * bnb(i-3,1)
          matl(2,2) = mat(5,2,2,i,1) * bnb(i-3,1) + mat(5,2,3,i,1) * bnb(i-3,2)
          matl(3,1) = mat(5,3,2,i,1) * bnb(i-3,2) - mat(5,3,3,i,1) * bnb(i-3,1)
          matl(3,2) = mat(5,3,2,i,1) * bnb(i-3,1) + mat(5,3,3,i,1) * bnb(i-3,2)
          matl(4,1) = mat(5,4,2,i,1) * bnb(i-3,2) - mat(5,4,3,i,1) * bnb(i-3,1)
          matl(4,2) = mat(5,4,2,i,1) * bnb(i-3,1) + mat(5,4,3,i,1) * bnb(i-3,2)
          matl(5,1) = mat(5,5,2,i,1) * bnb(i-3,2) - mat(5,5,3,i,1) * bnb(i-3,1)
          matl(5,2) = mat(5,5,2,i,1) * bnb(i-3,1) + mat(5,5,3,i,1) * bnb(i-3,2)
          mat(5,:,2:3,i,1) = matl
          
          matl(1,1) = mat(1,1,2,i,1) * bnb(i-2,2) - mat(1,1,3,i,1) * bnb(i-2,1)
          matl(1,2) = mat(1,1,2,i,1) * bnb(i-2,1) + mat(1,1,3,i,1) * bnb(i-2,2)
          matl(2,1) = mat(1,2,2,i,1) * bnb(i-2,2) - mat(1,2,3,i,1) * bnb(i-2,1)
          matl(2,2) = mat(1,2,2,i,1) * bnb(i-2,1) + mat(1,2,3,i,1) * bnb(i-2,2)
          matl(3,1) = mat(1,3,2,i,1) * bnb(i-2,2) - mat(1,3,3,i,1) * bnb(i-2,1)
          matl(3,2) = mat(1,3,2,i,1) * bnb(i-2,1) + mat(1,3,3,i,1) * bnb(i-2,2)
          matl(4,1) = mat(1,4,2,i,1) * bnb(i-2,2) - mat(1,4,3,i,1) * bnb(i-2,1)
          matl(4,2) = mat(1,4,2,i,1) * bnb(i-2,1) + mat(1,4,3,i,1) * bnb(i-2,2)
          matl(5,1) = mat(1,5,2,i,1) * bnb(i-2,2) - mat(1,5,3,i,1) * bnb(i-2,1)
          matl(5,2) = mat(1,5,2,i,1) * bnb(i-2,1) + mat(1,5,3,i,1) * bnb(i-2,2)
          mat(1,:,2:3,i,1) = matl
          
          matl(1,1) = mat(2,1,2,i,1) * bnb(i-1,2) - mat(2,1,3,i,1) * bnb(i-1,1)
          matl(1,2) = mat(2,1,2,i,1) * bnb(i-1,1) + mat(2,1,3,i,1) * bnb(i-1,2)
          matl(2,1) = mat(2,2,2,i,1) * bnb(i-1,2) - mat(2,2,3,i,1) * bnb(i-1,1)
          matl(2,2) = mat(2,2,2,i,1) * bnb(i-1,1) + mat(2,2,3,i,1) * bnb(i-1,2)
          matl(3,1) = mat(2,3,2,i,1) * bnb(i-1,2) - mat(2,3,3,i,1) * bnb(i-1,1)
          matl(3,2) = mat(2,3,2,i,1) * bnb(i-1,1) + mat(2,3,3,i,1) * bnb(i-1,2)
          matl(4,1) = mat(2,4,2,i,1) * bnb(i-1,2) - mat(2,4,3,i,1) * bnb(i-1,1)
          matl(4,2) = mat(2,4,2,i,1) * bnb(i-1,1) + mat(2,4,3,i,1) * bnb(i-1,2)
          matl(5,1) = mat(2,5,2,i,1) * bnb(i-1,2) - mat(2,5,3,i,1) * bnb(i-1,1)
          matl(5,2) = mat(2,5,2,i,1) * bnb(i-1,1) + mat(2,5,3,i,1) * bnb(i-1,2)
          mat(2,:,2:3,i,1) = matl
          
          matl(1,1) = mat(3,1,2,i,1) * bnb(i,2) - mat(3,1,3,i,1) * bnb(i,1)
          matl(1,2) = mat(3,1,2,i,1) * bnb(i,1) + mat(3,1,3,i,1) * bnb(i,2)
          matl(2,1) = mat(3,2,2,i,1) * bnb(i,2) - mat(3,2,3,i,1) * bnb(i,1)
          matl(2,2) = mat(3,2,2,i,1) * bnb(i,1) + mat(3,2,3,i,1) * bnb(i,2)
          matl(3,1) = mat(3,3,2,i,1) * bnb(i,2) - mat(3,3,3,i,1) * bnb(i,1)
          matl(3,2) = mat(3,3,2,i,1) * bnb(i,1) + mat(3,3,3,i,1) * bnb(i,2)
          matl(4,1) = mat(3,4,2,i,1) * bnb(i,2) - mat(3,4,3,i,1) * bnb(i,1)
          matl(4,2) = mat(3,4,2,i,1) * bnb(i,1) + mat(3,4,3,i,1) * bnb(i,2)
          matl(5,1) = mat(3,5,2,i,1) * bnb(i,2) - mat(3,5,3,i,1) * bnb(i,1)
          matl(5,2) = mat(3,5,2,i,1) * bnb(i,1) + mat(3,5,3,i,1) * bnb(i,2)
          mat(3,:,2:3,i,1) = matl
          
          matl(1,1) = mat(4,1,2,i,1) * bnb(i+1,2) - mat(4,1,3,i,1) * bnb(i+1,1)
          matl(1,2) = mat(4,1,2,i,1) * bnb(i+1,1) + mat(4,1,3,i,1) * bnb(i+1,2)
          matl(2,1) = mat(4,2,2,i,1) * bnb(i+1,2) - mat(4,2,3,i,1) * bnb(i+1,1)
          matl(2,2) = mat(4,2,2,i,1) * bnb(i+1,1) + mat(4,2,3,i,1) * bnb(i+1,2)
          matl(3,1) = mat(4,3,2,i,1) * bnb(i+1,2) - mat(4,3,3,i,1) * bnb(i+1,1)
          matl(3,2) = mat(4,3,2,i,1) * bnb(i+1,1) + mat(4,3,3,i,1) * bnb(i+1,2)
          matl(4,1) = mat(4,4,2,i,1) * bnb(i+1,2) - mat(4,4,3,i,1) * bnb(i+1,1)
          matl(4,2) = mat(4,4,2,i,1) * bnb(i+1,1) + mat(4,4,3,i,1) * bnb(i+1,2)
          matl(5,1) = mat(4,5,2,i,1) * bnb(i+1,2) - mat(4,5,3,i,1) * bnb(i+1,1)
          matl(5,2) = mat(4,5,2,i,1) * bnb(i+1,1) + mat(4,5,3,i,1) * bnb(i+1,2)
          mat(4,:,2:3,i,1) = matl
          
          i = nx
          
          matl(1,1) = mat(4,1,2,i,1) * bnb(i-4,2) - mat(4,1,3,i,1) * bnb(i-4,1)
          matl(1,2) = mat(4,1,2,i,1) * bnb(i-4,1) + mat(4,1,3,i,1) * bnb(i-4,2)
          matl(2,1) = mat(4,2,2,i,1) * bnb(i-4,2) - mat(4,2,3,i,1) * bnb(i-4,1)
          matl(2,2) = mat(4,2,2,i,1) * bnb(i-4,1) + mat(4,2,3,i,1) * bnb(i-4,2)
          matl(3,1) = mat(4,3,2,i,1) * bnb(i-4,2) - mat(4,3,3,i,1) * bnb(i-4,1)
          matl(3,2) = mat(4,3,2,i,1) * bnb(i-4,1) + mat(4,3,3,i,1) * bnb(i-4,2)
          matl(4,1) = mat(4,4,2,i,1) * bnb(i-4,2) - mat(4,4,3,i,1) * bnb(i-4,1)
          matl(4,2) = mat(4,4,2,i,1) * bnb(i-4,1) + mat(4,4,3,i,1) * bnb(i-4,2)
          matl(5,1) = mat(4,5,2,i,1) * bnb(i-4,2) - mat(4,5,3,i,1) * bnb(i-4,1)
          matl(5,2) = mat(4,5,2,i,1) * bnb(i-4,1) + mat(4,5,3,i,1) * bnb(i-4,2)
          mat(4,:,2:3,i,1) = matl
          
          matl(1,1) = mat(5,1,2,i,1) * bnb(i-3,2) - mat(5,1,3,i,1) * bnb(i-3,1)
          matl(1,2) = mat(5,1,2,i,1) * bnb(i-3,1) + mat(5,1,3,i,1) * bnb(i-3,2)
          matl(2,1) = mat(5,2,2,i,1) * bnb(i-3,2) - mat(5,2,3,i,1) * bnb(i-3,1)
          matl(2,2) = mat(5,2,2,i,1) * bnb(i-3,1) + mat(5,2,3,i,1) * bnb(i-3,2)
          matl(3,1) = mat(5,3,2,i,1) * bnb(i-3,2) - mat(5,3,3,i,1) * bnb(i-3,1)
          matl(3,2) = mat(5,3,2,i,1) * bnb(i-3,1) + mat(5,3,3,i,1) * bnb(i-3,2)
          matl(4,1) = mat(5,4,2,i,1) * bnb(i-3,2) - mat(5,4,3,i,1) * bnb(i-3,1)
          matl(4,2) = mat(5,4,2,i,1) * bnb(i-3,1) + mat(5,4,3,i,1) * bnb(i-3,2)
          matl(5,1) = mat(5,5,2,i,1) * bnb(i-3,2) - mat(5,5,3,i,1) * bnb(i-3,1)
          matl(5,2) = mat(5,5,2,i,1) * bnb(i-3,1) + mat(5,5,3,i,1) * bnb(i-3,2)
          mat(5,:,2:3,i,1) = matl
          
          matl(1,1) = mat(1,1,2,i,1) * bnb(i-2,2) - mat(1,1,3,i,1) * bnb(i-2,1)
          matl(1,2) = mat(1,1,2,i,1) * bnb(i-2,1) + mat(1,1,3,i,1) * bnb(i-2,2)
          matl(2,1) = mat(1,2,2,i,1) * bnb(i-2,2) - mat(1,2,3,i,1) * bnb(i-2,1)
          matl(2,2) = mat(1,2,2,i,1) * bnb(i-2,1) + mat(1,2,3,i,1) * bnb(i-2,2)
          matl(3,1) = mat(1,3,2,i,1) * bnb(i-2,2) - mat(1,3,3,i,1) * bnb(i-2,1)
          matl(3,2) = mat(1,3,2,i,1) * bnb(i-2,1) + mat(1,3,3,i,1) * bnb(i-2,2)
          matl(4,1) = mat(1,4,2,i,1) * bnb(i-2,2) - mat(1,4,3,i,1) * bnb(i-2,1)
          matl(4,2) = mat(1,4,2,i,1) * bnb(i-2,1) + mat(1,4,3,i,1) * bnb(i-2,2)
          matl(5,1) = mat(1,5,2,i,1) * bnb(i-2,2) - mat(1,5,3,i,1) * bnb(i-2,1)
          matl(5,2) = mat(1,5,2,i,1) * bnb(i-2,1) + mat(1,5,3,i,1) * bnb(i-2,2)
          mat(1,:,2:3,i,1) = matl
          
          matl(1,1) = mat(2,1,2,i,1) * bnb(i-1,2) - mat(2,1,3,i,1) * bnb(i-1,1)
          matl(1,2) = mat(2,1,2,i,1) * bnb(i-1,1) + mat(2,1,3,i,1) * bnb(i-1,2)
          matl(2,1) = mat(2,2,2,i,1) * bnb(i-1,2) - mat(2,2,3,i,1) * bnb(i-1,1)
          matl(2,2) = mat(2,2,2,i,1) * bnb(i-1,1) + mat(2,2,3,i,1) * bnb(i-1,2)
          matl(3,1) = mat(2,3,2,i,1) * bnb(i-1,2) - mat(2,3,3,i,1) * bnb(i-1,1)
          matl(3,2) = mat(2,3,2,i,1) * bnb(i-1,1) + mat(2,3,3,i,1) * bnb(i-1,2)
          matl(4,1) = mat(2,4,2,i,1) * bnb(i-1,2) - mat(2,4,3,i,1) * bnb(i-1,1)
          matl(4,2) = mat(2,4,2,i,1) * bnb(i-1,1) + mat(2,4,3,i,1) * bnb(i-1,2)
          matl(5,1) = mat(2,5,2,i,1) * bnb(i-1,2) - mat(2,5,3,i,1) * bnb(i-1,1)
          matl(5,2) = mat(2,5,2,i,1) * bnb(i-1,1) + mat(2,5,3,i,1) * bnb(i-1,2)
          mat(2,:,2:3,i,1) = matl
          
          matl(1,1) = mat(3,1,2,i,1) * bnb(i,2) - mat(3,1,3,i,1) * bnb(i,1)
          matl(1,2) = mat(3,1,2,i,1) * bnb(i,1) + mat(3,1,3,i,1) * bnb(i,2)
          matl(2,1) = mat(3,2,2,i,1) * bnb(i,2) - mat(3,2,3,i,1) * bnb(i,1)
          matl(2,2) = mat(3,2,2,i,1) * bnb(i,1) + mat(3,2,3,i,1) * bnb(i,2)
          matl(3,1) = mat(3,3,2,i,1) * bnb(i,2) - mat(3,3,3,i,1) * bnb(i,1)
          matl(3,2) = mat(3,3,2,i,1) * bnb(i,1) + mat(3,3,3,i,1) * bnb(i,2)
          matl(4,1) = mat(3,4,2,i,1) * bnb(i,2) - mat(3,4,3,i,1) * bnb(i,1)
          matl(4,2) = mat(3,4,2,i,1) * bnb(i,1) + mat(3,4,3,i,1) * bnb(i,2)
          matl(5,1) = mat(3,5,2,i,1) * bnb(i,2) - mat(3,5,3,i,1) * bnb(i,1)
          matl(5,2) = mat(3,5,2,i,1) * bnb(i,1) + mat(3,5,3,i,1) * bnb(i,2)
          mat(3,:,2:3,i,1) = matl
          
!.... now rotate the equations to the normal direction

          do i = 1, nx
            mat(:,2,:,i,1) = bnb(i,2) * mat(:,2,:,i,1) - &
                             bnb(i,1) * mat(:,3,:,i,1)
            mat(:,3,:,i,1) = bnb(i,1) * mat(:,2,:,i,1) + &
                             bnb(i,2) * mat(:,3,:,i,1)
          end do

!.... constrain the normal velocity to be zero

          mat(:,3,:,:,1) = zero
          mat(3,3,3,:,1) = one

        end if                  ! Navier
          
!.... Riemann Invariants on top boundary

        if (top.eq.0) then
          mat(:,:,:,:,ny) = zero
          mat(3,1,1,:,ny) = one
          mat(3,2,2,:,ny) = one
          mat(3,3,3,:,ny) = one
          mat(3,4,4,:,ny) = one
          mat(3,5,5,:,ny) = one
        end if
        
        end if          ! yper

!=============================================================================!
!       L e f t   B o u n d a r y
!=============================================================================!
        if (.not. xper) then

!.... constant pressure

!       mat(:,1,:,1,:) = zero
!       mat(3,1,1,1,:) = one

!.... symmetry plane

        if (left.eq.2) then
          mat(:,1,:,1,:) = zero
          mat(3,1,1,1,:) = -gc1
          mat(4,1,1,1,:) = -gc2
          mat(5,1,1,1,:) = -gc3
          mat(1,1,1,1,:) = -gc4
          mat(2,1,1,1,:) = -gc5

          mat(:,2,:,1,:) = zero
          mat(3,2,2,1,:) = -gc1
          mat(4,2,2,1,:) = -gc2
          mat(5,2,2,1,:) = -gc3
          mat(1,2,2,1,:) = -gc4
          mat(2,2,2,1,:) = -gc5

          mat(:,3,:,1,:) = zero
          mat(3,3,3,1,:) = one

          mat(:,4,:,1,:) = zero
          mat(3,4,4,1,:) = -gc1
          mat(4,4,4,1,:) = -gc2
          mat(5,4,4,1,:) = -gc3
          mat(1,4,4,1,:) = -gc4
          mat(2,4,4,1,:) = -gc5

          mat(:,5,:,1,:) = zero
          mat(3,5,5,1,:) = -gc1
          mat(4,5,5,1,:) = -gc2
          mat(5,5,5,1,:) = -gc3
          mat(1,5,5,1,:) = -gc4
          mat(2,5,5,1,:) = -gc5
        end if

        if (left.eq.7) then             ! symmetry plane
          mat(:,3,:,1,:) = zero
          mat(3,3,3,1,:) = one
        endif

!.... Left Riemann invariant boundary condition (first-order)
        
        if (left.eq.0 .and. Ma.lt.one) then
          mat(:,:,:,1,:) = zero
          call ReimannLHS( ny, vl(:,1,:), vl(:,2,:), vl(:,3,:), &
                           bnl, x(1,:), y(1,:), mat(3,:,:,1,:), &
                           mat(4,:,:,1,:), mat(5,:,:,1,:), &
                           rhobl, ubl, vbl, wbl, tbl, pbl, cbl )
        end if

!.... Left Side:  Riemann invariants (zero-order)

          if (.false.) then
          
          mat(:,:,:,1,:) = zero

          R1(:) = one - two / Ma / gamma1
          R2(:) = vl(2,2,:) + &
                  two * sqrt( vl(5,2,:) ) / Ma / gamma1
                
          rho(:) = vl(1,2,:)
          p(:)   = vl(1,2,:) * vl(5,2,:) / (gamma * Ma**2)

          fact1(:) = (rho(:)**gamma/(gamma*p(:)))**(one/gamma1) * &
                     (gamma1/four)**(two/gamma1) * (two/gamma1) * &
                     (R2(:) - R1(:))**((two-gamma1)/gamma1)
          fact2(:) = (gamma1/four * (R2(:) - R1(:)))**(two/gamma1) * &
                     one/(gamma*gamma1) * &
                     (rho(:)**gamma/(gamma*p(:)))**((two-gamma)/gamma1)
          mat(3,1,1,1,:) = -one
          mat(4,1,2,1,:) = fact1(:) + &
                           gamma1 * rho(:)**gamma1 / p(:)
          mat(4,1,5,1,:) = fact1(:) / ( gamma1 * Ma * sqrt(vl(5,2,:)) ) - &
                           fact2(:) * rho(:)**(gamma+1)/(gamma*Ma**2*p(:)**2)
          
          mat(3,2,2,1,:) = -one
          mat(4,2,2,1,:) = pt5
          mat(4,2,5,1,:) = pt5 / ( gamma1 * Ma * sqrt(vl(5,2,:)) )
  
          mat(3,3,3,1,:) = -one
          mat(4,3,3,1,:) = one
  
          mat(3,4,4,1,:) = -one
          mat(4,4,4,1,:) = one
  
          fact(:) = (Ma * gamma1 / four)**2 * two * ( R2(:) - R1(:) )
          mat(3,5,5,1,:) = -one
          mat(4,5,2,1,:) = fact(:)
          mat(4,5,3,1,:) = zero
          mat(4,5,5,1,:) = fact(:) / ( gamma1 * Ma * sqrt(vl(5,2,:)) )
          
          end if                ! Riemann boundary
        
          if (.false.) then

!.... Zero-order extrapolation in the viscous layer

          if (extrap.eq.0) then
            mat(:,:,:,1,1:nbl) = zero
  
            mat(3,1,1,1,1:nbl) = -one
            mat(3,2,2,1,1:nbl) = -one
            mat(3,3,3,1,1:nbl) = -one
            mat(3,4,4,1,1:nbl) = -one
            mat(3,5,5,1,1:nbl) = -one
  
            mat(4,1,1,1,1:nbl) = one
            mat(4,2,2,1,1:nbl) = one
            mat(4,3,3,1,1:nbl) = one
            mat(4,4,4,1,1:nbl) = one
            mat(4,5,5,1,1:nbl) = one
          end if
          
!.... First-order extrapolation in the viscous layer

          if (extrap.eq.1) then
            mat(:,:,:,1,1:nbl) = zero
  
            mat(3,1,1,1,1:nbl) = -one
            mat(3,2,2,1,1:nbl) = -one
            mat(3,3,3,1,1:nbl) = -one
            mat(3,4,4,1,1:nbl) = -one
            mat(3,5,5,1,1:nbl) = -one
  
            mat(4,1,1,1,1:nbl) = two
            mat(4,2,2,1,1:nbl) = two
            mat(4,3,3,1,1:nbl) = two
            mat(4,4,4,1,1:nbl) = two
            mat(4,5,5,1,1:nbl) = two
            
            mat(5,1,1,1,1:nbl) = -one
            mat(5,2,2,1,1:nbl) = -one
            mat(5,3,3,1,1:nbl) = -one
            mat(5,4,4,1,1:nbl) = -one
            mat(5,5,5,1,1:nbl) = -one
          end if                ! extrapolation type

          end if
!=============================================================================!
!       R i g h t  B o u n d a r y
!=============================================================================!

!.... constant pressure

!       mat(:,1,:,nx,:) = zero
!       mat(3,1,1,nx,:) = one

!.... Right Riemann invariant boundary condition (first-order)
        
        if (right.eq.0 .and. Ma.lt.one) then

          mat(:,:,:,nx,:) = zero
          call ReimannLHS( ny, vl(:,nx,:), vl(:,nx-1,:), vl(:,nx-2,:), &
                           bnr, x(nx,:), y(nx,:), mat(3,:,:,nx,:), &
                           mat(2,:,:,nx,:), mat(1,:,:,nx,:), &
                           rhobr, ubr, vbr, wbr, tbr, pbr, cbr )
                
        end if

!.... Right Side:  Riemann invariants (zero-order)

          if (.false.) then
          
          mat(:,:,:,nx,:) = zero

          R1(:) = one - two / Ma / gamma1
          R2(:) = vl(2,nx-1,:) + &
                  two * sqrt( vl(5,nx-1,:) ) / Ma / gamma1
                
          rho(:) = vl(1,nx-1,:)
          p(:)   = vl(1,nx-1,:) * vl(5,nx-1,:) / (gamma * Ma**2)

          fact1(:) = (rho(:)**gamma/(gamma*p(:)))**(one/gamma1) * &
                     (gamma1/four)**(two/gamma1) * (two/gamma1) * &
                     (R2(:) - R1(:))**((two-gamma1)/gamma1)
          fact2(:) = (gamma1/four * (R2(:) - R1(:)))**(two/gamma1) * &
                     one/(gamma*gamma1) * &
                     (rho(:)**gamma/(gamma*p(:)))**((two-gamma)/gamma1)
          mat(3,1,1,nx,:) = -one
          mat(2,1,2,nx,:) = fact1(:) + &
                            gamma1 * rho(:)**gamma1 / p(:)
          mat(2,1,5,nx,:) = fact1(:) / ( gamma1 * Ma * sqrt(vl(5,2,:)) ) - &
                            fact2(:) * rho(:)**(gamma+1)/(gamma*Ma**2*p(:)**2)
          
          mat(3,2,2,nx,:) = -one
          mat(2,2,2,nx,:) = pt5
          mat(2,2,5,nx,:) = pt5 / ( gamma1 * Ma * sqrt(vl(5,2,:)) )
  
          mat(3,3,3,nx,:) = -one
          mat(2,3,3,nx,:) = one
  
          mat(3,4,4,nx,:) = -one
          mat(2,4,4,nx,:) = one
  
          fact(:) = (Ma * gamma1 / four)**2 * two * ( R2(:) - R1(:) )
          mat(3,5,5,nx,:) = -one
          mat(2,5,2,nx,:) = fact(:)
          mat(2,5,3,nx,:) = zero
          mat(2,5,5,nx,:) = fact(:) / ( gamma1 * Ma * sqrt(vl(5,2,:)) )

          end if                ! Riemann boundary
          
!.... Zero'th order extrapolation in the viscous layer

          if (extrap.eq.0) then
            mat(:,:,:,nx,1:nbl) = zero
            
            mat(3,1,1,nx,1:nbl) = -one
            mat(3,2,2,nx,1:nbl) = -one
            mat(3,3,3,nx,1:nbl) = -one
            mat(3,4,4,nx,1:nbl) = -one
            mat(3,5,5,nx,1:nbl) = -one
  
            mat(2,1,1,nx,1:nbl) = one
            mat(2,2,2,nx,1:nbl) = one
            mat(2,3,3,nx,1:nbl) = one
            mat(2,4,4,nx,1:nbl) = one
            mat(2,5,5,nx,1:nbl) = one
          end if
          
!.... First-order extrapolation in the viscous layer

          if (extrap.eq.1) then
            mat(:,:,:,nx,1:nbl) = zero
  
            mat(3,1,1,nx,1:nbl) = -one
            mat(3,2,2,nx,1:nbl) = -one
            mat(3,3,3,nx,1:nbl) = -one
            mat(3,4,4,nx,1:nbl) = -one
            mat(3,5,5,nx,1:nbl) = -one
  
            mat(2,1,1,nx,1:nbl) = two
            mat(2,2,2,nx,1:nbl) = two
            mat(2,3,3,nx,1:nbl) = two
            mat(2,4,4,nx,1:nbl) = two
            mat(2,5,5,nx,1:nbl) = two
            
            mat(1,1,1,nx,1:nbl) = -one
            mat(1,2,2,nx,1:nbl) = -one
            mat(1,3,3,nx,1:nbl) = -one
            mat(1,4,4,nx,1:nbl) = -one
            mat(1,5,5,nx,1:nbl) = -one
          end if                ! extrapolation type

!.... Second-order extrapolation in the viscous layer

          if (extrap.eq.2) then

            call error("itrbc1f$","extrap.eq.2 not supported$")

            mat(:,:,:,nx,1:nbl) = zero

            mat(3,1,1,nx,1:nbl) = -one
            mat(3,2,2,nx,1:nbl) = -one
            mat(3,3,3,nx,1:nbl) = -one
            mat(3,4,4,nx,1:nbl) = -one
            mat(3,5,5,nx,1:nbl) = -one

            mat(2,1,1,nx,1:nbl) = three 
            mat(2,2,2,nx,1:nbl) = three 
            mat(2,3,3,nx,1:nbl) = three 
            mat(2,4,4,nx,1:nbl) = three 
            mat(2,5,5,nx,1:nbl) = three 

            mat(1,1,1,nx,1:nbl) = -three
            mat(1,2,2,nx,1:nbl) = -three
            mat(1,3,3,nx,1:nbl) = -three
            mat(1,4,4,nx,1:nbl) = -three
            mat(1,5,5,nx,1:nbl) = -three
          end if                ! extrapolation type

!.... hold the IC

          if (right.eq.8) then
            mat(:,:,:,nx,:) = zero
            mat(3,1,1,nx,:) = -one
            mat(3,2,2,nx,:) = -one
            mat(3,3,3,nx,:) = -one
            mat(3,4,4,nx,:) = -one
            mat(3,5,5,nx,:) = -one
          end if

        end if                  ! xper
        
        end if                  ! linear
        
        return
        end
