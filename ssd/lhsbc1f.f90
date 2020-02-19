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
        use stuff
        use diff
        use pot
        implicit none
        
        real :: mat(ny,nx,ndof,ndof,5)
        real :: vl(ny,nx,ndof), vml(ny,nx,ndof)

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
          mat(1,:,1,:,:) = zero
          mat(1,:,1,1,3) = one
        end if

!.... linear extrapolation of rho at the wall

        if (wall.eq.2) then
          mat(1,:,1,:,:) = zero
          mat(1,:,1,1,3) = one
        end if

!.... no-slip

        mat(1,:,2:4,:,:) = zero
        mat(1,:,2,2,3)   = one
        mat(1,:,3,3,3)   = one
        mat(1,:,4,4,3)   = one

!.... wall temperature boundary condition

        if (wallt.eq.0 .or. wallt.eq.1) then
          mat(1,:,5,:,:) = zero
          mat(1,:,5,5,3) = one
        end if
        
        else                    ! inviscid wall

          mat(1,:,1,:,:) = zero
          mat(1,:,1,1,3) = one
          
          do i = 1, nx
            if ( abs(bnb(i,2)) .ge. abs(bnb(i,1)) .and. .false.) then
              mat(1,i,2,:,:) = bnb(i,2) * mat(1,i,2,:,:) - &
                               bnb(i,1) * mat(1,i,3,:,:)
              mat(1,i,3,:,:) = zero
              mat(1,i,3,2,3) = bnb(i,1)
              mat(1,i,3,3,3) = bnb(i,2)
            else
              mat(1,i,3,:,:) = bnb(i,2) * mat(1,i,2,:,:) - &
                               bnb(i,1) * mat(1,i,3,:,:)
              mat(1,i,2,:,:) = zero
              mat(1,i,2,2,3) = bnb(i,1)
              mat(1,i,2,3,3) = bnb(i,2)
            end if
          end do
        end if                  ! Navier

!.... freestream zero disturbance boundary condition

        if (top.eq.0) then
          mat(ny,:,:,:,:) = zero
          mat(ny,:,1,1,3) = one
          mat(ny,:,2,2,3) = one
          mat(ny,:,3,3,3) = one
          mat(ny,:,4,4,3) = one
          mat(ny,:,5,5,3) = one
        end if
        
        end if          ! yper

!=============================================================================!
        
        if (.not. xper) then

!.... left boundary

          if (left.eq.0) then           ! zero disturbance

            mat(:,1,:,:,:) = zero
            mat(:,1,1,1,3) = one
            mat(:,1,2,2,3) = one
            mat(:,1,3,3,3) = one
            mat(:,1,4,4,3) = one
            mat(:,1,5,5,3) = one

          else if (left.eq.1) then      ! nonreflecting BC

            mat(nbl+1:ny,1,1,:,:) = zero
            mat(nbl+1:ny,1,1,1,3) = -one
            mat(nbl+1:ny,1,1,2,3) = gamma * Ma * vml(nbl+1:ny,1,1) * &
                                    sqrt(vml(nbl+1:ny,1,5)) / vml(nbl+1:ny,1,5)
            mat(nbl+1:ny,1,1,5,3) = -vml(nbl+1:ny,1,1) / vml(nbl+1:ny,1,5)

            if (.false.) then
              mat(1:nbl,1,:,:,:) = zero

              mat(1:nbl,1,1,1,3) = -one
              mat(1:nbl,1,2,2,3) = -one
              mat(1:nbl,1,3,3,3) = -one
              mat(1:nbl,1,4,4,3) = -one
              mat(1:nbl,1,5,5,3) = -one

              mat(1:nbl,1,1,1,4) = one
              mat(1:nbl,1,2,2,4) = one
              mat(1:nbl,1,3,3,4) = one
              mat(1:nbl,1,4,4,4) = one
              mat(1:nbl,1,5,5,4) = one
            end if

          else if (left.eq.2) then      ! symmetry plane

!           mat(:,1,1,:,:) = zero
!           mat(:,1,1,1,3) = -gc1
!           mat(:,1,1,1,4) = -gc2
!           mat(:,1,1,1,5) = -gc3
!           mat(:,1,1,1,1) = -gc4
!           mat(:,1,1,1,2) = -gc5
  
            mat(:,1,2,:,:) = zero
            mat(:,1,2,2,3) = -gc1
            mat(:,1,2,2,4) = -gc2
            mat(:,1,2,2,5) = -gc3
            mat(:,1,2,2,1) = -gc4
            mat(:,1,2,2,2) = -gc5

            mat(:,1,3,:,:) = zero
            mat(:,1,3,3,3) = one

            mat(:,1,4,:,:) = zero
            mat(:,1,4,4,3) = -gc1
            mat(:,1,4,4,4) = -gc2
            mat(:,1,4,4,5) = -gc3
            mat(:,1,4,4,1) = -gc4
            mat(:,1,4,4,2) = -gc5

            mat(:,1,5,:,:) = zero
            mat(:,1,5,5,3) = -gc1
            mat(:,1,5,5,4) = -gc2
            mat(:,1,5,5,5) = -gc3
            mat(:,1,5,5,1) = -gc4
            mat(:,1,5,5,2) = -gc5

          else if (left.eq.7) then      ! symmetry plane

!           mat(:,1,3,:,:) = zero
!           mat(:,1,3,3,3) = one

          endif        ! left

!.... Right boundary

          if (right.eq.0) then          ! zero disturbance

            mat(:,nx,:,:,:) = zero
            mat(:,nx,1,1,3) = one
            mat(:,nx,2,2,3) = one
            mat(:,nx,3,3,3) = one
            mat(:,nx,4,4,3) = one
            mat(:,nx,5,5,3) = one

          else if (right.eq.1) then     ! nonreflecting

            mat(nbl+1:ny,1,1,:,:) = zero
            mat(nbl+1:ny,1,1,1,3) = -one
            mat(nbl+1:ny,1,1,2,3) = gamma * Ma * vml(nbl+1:ny,1,1) * &
                                    sqrt(vml(nbl+1:ny,1,5)) / vml(nbl+1:ny,1,5)
            mat(nbl+1:ny,1,1,5,3) = -vml(nbl+1:ny,1,1) / vml(nbl+1:ny,1,5)

            mat(nbl+1:ny,nx,1,:,:) = zero
            mat(nbl+1:ny,nx,1,1,3) = -one
            mat(nbl+1:ny,nx,1,2,3) = gamma * Ma * vml(nbl+1:ny,nx,1) * &
                                     sqrt(vml(nbl+1:ny,nx,5)) / &
                                     vml(nbl+1:ny,nx,5)
            mat(nbl+1:ny,nx,1,5,3) = -vml(nbl+1:ny,nx,1) / vml(nbl+1:ny,nx,5)

            if (.false.) then
              mat(1:nbl,nx,:,:,:) = zero

              mat(1:nbl,nx,1,1,3) = -one
              mat(1:nbl,nx,2,2,3) = -one
              mat(1:nbl,nx,3,3,3) = -one
              mat(1:nbl,nx,4,4,3) = -one
              mat(1:nbl,nx,5,5,3) = -one

              mat(1:nbl,nx,1,1,2) = one
              mat(1:nbl,nx,2,2,2) = one
              mat(1:nbl,nx,3,3,2) = one
              mat(1:nbl,nx,4,4,2) = one
              mat(1:nbl,nx,5,5,2) = one
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

        if (wall.eq.1 .or. wall.eq.2) then
          mat(1,:,1,:,:) = zero
          mat(1,:,1,1,3) = one
        end if

!.... no-slip

        mat(1,:,2:4,:,:) = zero
        mat(1,:,2,2,3)   = one
        mat(1,:,3,3,3)   = one
        mat(1,:,4,4,3)   = one

!.... adiabatic & isothermal

        if (wallt.eq.0 .or. wallt.eq.1) then
          mat(1,:,5,:,:) = zero
          mat(1,:,5,5,3) = one
        end if

        else      ! inviscid:  rotate to body normal coordinates

            i = 1

            matl(1,1) = mat(1,i,1,2,3) * bnb(i,2) - mat(1,i,1,3,3) * bnb(i,1)
            matl(1,2) = mat(1,i,1,2,3) * bnb(i,1) + mat(1,i,1,3,3) * bnb(i,2)
            matl(2,1) = mat(1,i,2,2,3) * bnb(i,2) - mat(1,i,2,3,3) * bnb(i,1)
            matl(2,2) = mat(1,i,2,2,3) * bnb(i,1) + mat(1,i,2,3,3) * bnb(i,2)
            matl(3,1) = mat(1,i,3,2,3) * bnb(i,2) - mat(1,i,3,3,3) * bnb(i,1)
            matl(3,2) = mat(1,i,3,2,3) * bnb(i,1) + mat(1,i,3,3,3) * bnb(i,2)
            matl(4,1) = mat(1,i,4,2,3) * bnb(i,2) - mat(1,i,4,3,3) * bnb(i,1)
            matl(4,2) = mat(1,i,4,2,3) * bnb(i,1) + mat(1,i,4,3,3) * bnb(i,2)
            matl(5,1) = mat(1,i,5,2,3) * bnb(i,2) - mat(1,i,5,3,3) * bnb(i,1)
            matl(5,2) = mat(1,i,5,2,3) * bnb(i,1) + mat(1,i,5,3,3) * bnb(i,2)
            mat(1,i,:,2:3,3) = matl
            
            matl(1,1) = mat(1,i,1,2,4) * bnb(i+1,2) - mat(1,i,1,3,4) * bnb(i+1,1)
            matl(1,2) = mat(1,i,1,2,4) * bnb(i+1,1) + mat(1,i,1,3,4) * bnb(i+1,2)
            matl(2,1) = mat(1,i,2,2,4) * bnb(i+1,2) - mat(1,i,2,3,4) * bnb(i+1,1)
            matl(2,2) = mat(1,i,2,2,4) * bnb(i+1,1) + mat(1,i,2,3,4) * bnb(i+1,2)
            matl(3,1) = mat(1,i,3,2,4) * bnb(i+1,2) - mat(1,i,3,3,4) * bnb(i+1,1)
            matl(3,2) = mat(1,i,3,2,4) * bnb(i+1,1) + mat(1,i,3,3,4) * bnb(i+1,2)
            matl(4,1) = mat(1,i,4,2,4) * bnb(i+1,2) - mat(1,i,4,3,4) * bnb(i+1,1)
            matl(4,2) = mat(1,i,4,2,4) * bnb(i+1,1) + mat(1,i,4,3,4) * bnb(i+1,2)
            matl(5,1) = mat(1,i,5,2,4) * bnb(i+1,2) - mat(1,i,5,3,4) * bnb(i+1,1)
            matl(5,2) = mat(1,i,5,2,4) * bnb(i+1,1) + mat(1,i,5,3,4) * bnb(i+1,2)
            mat(1,i,:,2:3,4) = matl

            matl(1,1) = mat(1,i,1,2,5) * bnb(i+2,2) - mat(1,i,1,3,5) * bnb(i+2,1)
            matl(1,2) = mat(1,i,1,2,5) * bnb(i+2,1) + mat(1,i,1,3,5) * bnb(i+2,2)
            matl(2,1) = mat(1,i,2,2,5) * bnb(i+2,2) - mat(1,i,2,3,5) * bnb(i+2,1)
            matl(2,2) = mat(1,i,2,2,5) * bnb(i+2,1) + mat(1,i,2,3,5) * bnb(i+2,2)
            matl(3,1) = mat(1,i,3,2,5) * bnb(i+2,2) - mat(1,i,3,3,5) * bnb(i+2,1)
            matl(3,2) = mat(1,i,3,2,5) * bnb(i+2,1) + mat(1,i,3,3,5) * bnb(i+2,2)
            matl(4,1) = mat(1,i,4,2,5) * bnb(i+2,2) - mat(1,i,4,3,5) * bnb(i+2,1)
            matl(4,2) = mat(1,i,4,2,5) * bnb(i+2,1) + mat(1,i,4,3,5) * bnb(i+2,2)
            matl(5,1) = mat(1,i,5,2,5) * bnb(i+2,2) - mat(1,i,5,3,5) * bnb(i+2,1)
            matl(5,2) = mat(1,i,5,2,5) * bnb(i+2,1) + mat(1,i,5,3,5) * bnb(i+2,2)
            mat(1,i,:,2:3,5) = matl

            matl(1,1) = mat(1,i,1,2,1) * bnb(i+3,2) - mat(1,i,1,3,1) * bnb(i+3,1)
            matl(1,2) = mat(1,i,1,2,1) * bnb(i+3,1) + mat(1,i,1,3,1) * bnb(i+3,2)
            matl(2,1) = mat(1,i,2,2,1) * bnb(i+3,2) - mat(1,i,2,3,1) * bnb(i+3,1)
            matl(2,2) = mat(1,i,2,2,1) * bnb(i+3,1) + mat(1,i,2,3,1) * bnb(i+3,2)
            matl(3,1) = mat(1,i,3,2,1) * bnb(i+3,2) - mat(1,i,3,3,1) * bnb(i+3,1)
            matl(3,2) = mat(1,i,3,2,1) * bnb(i+3,1) + mat(1,i,3,3,1) * bnb(i+3,2)
            matl(4,1) = mat(1,i,4,2,1) * bnb(i+3,2) - mat(1,i,4,3,1) * bnb(i+3,1)
            matl(4,2) = mat(1,i,4,2,1) * bnb(i+3,1) + mat(1,i,4,3,1) * bnb(i+3,2)
            matl(5,1) = mat(1,i,5,2,1) * bnb(i+3,2) - mat(1,i,5,3,1) * bnb(i+3,1)
            matl(5,2) = mat(1,i,5,2,1) * bnb(i+3,1) + mat(1,i,5,3,1) * bnb(i+3,2)
            mat(1,i,:,2:3,1) = matl

            matl(1,1) = mat(1,i,1,2,2) * bnb(i+4,2) - mat(1,i,1,3,2) * bnb(i+4,1)
            matl(1,2) = mat(1,i,1,2,2) * bnb(i+4,1) + mat(1,i,1,3,2) * bnb(i+4,2)
            matl(2,1) = mat(1,i,2,2,2) * bnb(i+4,2) - mat(1,i,2,3,2) * bnb(i+4,1)
            matl(2,2) = mat(1,i,2,2,2) * bnb(i+4,1) + mat(1,i,2,3,2) * bnb(i+4,2)
            matl(3,1) = mat(1,i,3,2,2) * bnb(i+4,2) - mat(1,i,3,3,2) * bnb(i+4,1)
            matl(3,2) = mat(1,i,3,2,2) * bnb(i+4,1) + mat(1,i,3,3,2) * bnb(i+4,2)
            matl(4,1) = mat(1,i,4,2,2) * bnb(i+4,2) - mat(1,i,4,3,2) * bnb(i+4,1)
            matl(4,2) = mat(1,i,4,2,2) * bnb(i+4,1) + mat(1,i,4,3,2) * bnb(i+4,2)
            matl(5,1) = mat(1,i,5,2,2) * bnb(i+4,2) - mat(1,i,5,3,2) * bnb(i+4,1)
            matl(5,2) = mat(1,i,5,2,2) * bnb(i+4,1) + mat(1,i,5,3,2) * bnb(i+4,2)
            mat(1,i,:,2:3,2) = matl

            i = 2

            matl(1,1) = mat(1,i,1,2,2) * bnb(i-1,2) - mat(1,i,1,3,2) * bnb(i-1,1)
            matl(1,2) = mat(1,i,1,2,2) * bnb(i-1,1) + mat(1,i,1,3,2) * bnb(i-1,2)
            matl(2,1) = mat(1,i,2,2,2) * bnb(i-1,2) - mat(1,i,2,3,2) * bnb(i-1,1)
            matl(2,2) = mat(1,i,2,2,2) * bnb(i-1,1) + mat(1,i,2,3,2) * bnb(i-1,2)
            matl(3,1) = mat(1,i,3,2,2) * bnb(i-1,2) - mat(1,i,3,3,2) * bnb(i-1,1)
            matl(3,2) = mat(1,i,3,2,2) * bnb(i-1,1) + mat(1,i,3,3,2) * bnb(i-1,2)
            matl(4,1) = mat(1,i,4,2,2) * bnb(i-1,2) - mat(1,i,4,3,2) * bnb(i-1,1)
            matl(4,2) = mat(1,i,4,2,2) * bnb(i-1,1) + mat(1,i,4,3,2) * bnb(i-1,2)
            matl(5,1) = mat(1,i,5,2,2) * bnb(i-1,2) - mat(1,i,5,3,2) * bnb(i-1,1)
            matl(5,2) = mat(1,i,5,2,2) * bnb(i-1,1) + mat(1,i,5,3,2) * bnb(i-1,2)
            mat(1,i,:,2:3,2) = matl

            matl(1,1) = mat(1,i,1,2,3) * bnb(i,2) - mat(1,i,1,3,3) * bnb(i,1)
            matl(1,2) = mat(1,i,1,2,3) * bnb(i,1) + mat(1,i,1,3,3) * bnb(i,2)
            matl(2,1) = mat(1,i,2,2,3) * bnb(i,2) - mat(1,i,2,3,3) * bnb(i,1)
            matl(2,2) = mat(1,i,2,2,3) * bnb(i,1) + mat(1,i,2,3,3) * bnb(i,2)
            matl(3,1) = mat(1,i,3,2,3) * bnb(i,2) - mat(1,i,3,3,3) * bnb(i,1)
            matl(3,2) = mat(1,i,3,2,3) * bnb(i,1) + mat(1,i,3,3,3) * bnb(i,2)
            matl(4,1) = mat(1,i,4,2,3) * bnb(i,2) - mat(1,i,4,3,3) * bnb(i,1)
            matl(4,2) = mat(1,i,4,2,3) * bnb(i,1) + mat(1,i,4,3,3) * bnb(i,2)
            matl(5,1) = mat(1,i,5,2,3) * bnb(i,2) - mat(1,i,5,3,3) * bnb(i,1)
            matl(5,2) = mat(1,i,5,2,3) * bnb(i,1) + mat(1,i,5,3,3) * bnb(i,2)
            mat(1,i,:,2:3,3) = matl
            
            matl(1,1) = mat(1,i,1,2,4) * bnb(i+1,2) - mat(1,i,1,3,4) * bnb(i+1,1)
            matl(1,2) = mat(1,i,1,2,4) * bnb(i+1,1) + mat(1,i,1,3,4) * bnb(i+1,2)
            matl(2,1) = mat(1,i,2,2,4) * bnb(i+1,2) - mat(1,i,2,3,4) * bnb(i+1,1)
            matl(2,2) = mat(1,i,2,2,4) * bnb(i+1,1) + mat(1,i,2,3,4) * bnb(i+1,2)
            matl(3,1) = mat(1,i,3,2,4) * bnb(i+1,2) - mat(1,i,3,3,4) * bnb(i+1,1)
            matl(3,2) = mat(1,i,3,2,4) * bnb(i+1,1) + mat(1,i,3,3,4) * bnb(i+1,2)
            matl(4,1) = mat(1,i,4,2,4) * bnb(i+1,2) - mat(1,i,4,3,4) * bnb(i+1,1)
            matl(4,2) = mat(1,i,4,2,4) * bnb(i+1,1) + mat(1,i,4,3,4) * bnb(i+1,2)
            matl(5,1) = mat(1,i,5,2,4) * bnb(i+1,2) - mat(1,i,5,3,4) * bnb(i+1,1)
            matl(5,2) = mat(1,i,5,2,4) * bnb(i+1,1) + mat(1,i,5,3,4) * bnb(i+1,2)
            mat(1,i,:,2:3,4) = matl

            matl(1,1) = mat(1,i,1,2,5) * bnb(i+2,2) - mat(1,i,1,3,5) * bnb(i+2,1)
            matl(1,2) = mat(1,i,1,2,5) * bnb(i+2,1) + mat(1,i,1,3,5) * bnb(i+2,2)
            matl(2,1) = mat(1,i,2,2,5) * bnb(i+2,2) - mat(1,i,2,3,5) * bnb(i+2,1)
            matl(2,2) = mat(1,i,2,2,5) * bnb(i+2,1) + mat(1,i,2,3,5) * bnb(i+2,2)
            matl(3,1) = mat(1,i,3,2,5) * bnb(i+2,2) - mat(1,i,3,3,5) * bnb(i+2,1)
            matl(3,2) = mat(1,i,3,2,5) * bnb(i+2,1) + mat(1,i,3,3,5) * bnb(i+2,2)
            matl(4,1) = mat(1,i,4,2,5) * bnb(i+2,2) - mat(1,i,4,3,5) * bnb(i+2,1)
            matl(4,2) = mat(1,i,4,2,5) * bnb(i+2,1) + mat(1,i,4,3,5) * bnb(i+2,2)
            matl(5,1) = mat(1,i,5,2,5) * bnb(i+2,2) - mat(1,i,5,3,5) * bnb(i+2,1)
            matl(5,2) = mat(1,i,5,2,5) * bnb(i+2,1) + mat(1,i,5,3,5) * bnb(i+2,2)
            mat(1,i,:,2:3,5) = matl

            matl(1,1) = mat(1,i,1,2,1) * bnb(i+3,2) - mat(1,i,1,3,1) * bnb(i+3,1)
            matl(1,2) = mat(1,i,1,2,1) * bnb(i+3,1) + mat(1,i,1,3,1) * bnb(i+3,2)
            matl(2,1) = mat(1,i,2,2,1) * bnb(i+3,2) - mat(1,i,2,3,1) * bnb(i+3,1)
            matl(2,2) = mat(1,i,2,2,1) * bnb(i+3,1) + mat(1,i,2,3,1) * bnb(i+3,2)
            matl(3,1) = mat(1,i,3,2,1) * bnb(i+3,2) - mat(1,i,3,3,1) * bnb(i+3,1)
            matl(3,2) = mat(1,i,3,2,1) * bnb(i+3,1) + mat(1,i,3,3,1) * bnb(i+3,2)
            matl(4,1) = mat(1,i,4,2,1) * bnb(i+3,2) - mat(1,i,4,3,1) * bnb(i+3,1)
            matl(4,2) = mat(1,i,4,2,1) * bnb(i+3,1) + mat(1,i,4,3,1) * bnb(i+3,2)
            matl(5,1) = mat(1,i,5,2,1) * bnb(i+3,2) - mat(1,i,5,3,1) * bnb(i+3,1)
            matl(5,2) = mat(1,i,5,2,1) * bnb(i+3,1) + mat(1,i,5,3,1) * bnb(i+3,2)
            mat(1,i,:,2:3,1) = matl

          do i = 3, nx-2
            matl(1,1) = mat(1,i,1,2,1) * bnb(i-2,2) - mat(1,i,1,3,1) * bnb(i-2,1)
            matl(1,2) = mat(1,i,1,2,1) * bnb(i-2,1) + mat(1,i,1,3,1) * bnb(i-2,2)
            matl(2,1) = mat(1,i,2,2,1) * bnb(i-2,2) - mat(1,i,2,3,1) * bnb(i-2,1)
            matl(2,2) = mat(1,i,2,2,1) * bnb(i-2,1) + mat(1,i,2,3,1) * bnb(i-2,2)
            matl(3,1) = mat(1,i,3,2,1) * bnb(i-2,2) - mat(1,i,3,3,1) * bnb(i-2,1)
            matl(3,2) = mat(1,i,3,2,1) * bnb(i-2,1) + mat(1,i,3,3,1) * bnb(i-2,2)
            matl(4,1) = mat(1,i,4,2,1) * bnb(i-2,2) - mat(1,i,4,3,1) * bnb(i-2,1)
            matl(4,2) = mat(1,i,4,2,1) * bnb(i-2,1) + mat(1,i,4,3,1) * bnb(i-2,2)
            matl(5,1) = mat(1,i,5,2,1) * bnb(i-2,2) - mat(1,i,5,3,1) * bnb(i-2,1)
            matl(5,2) = mat(1,i,5,2,1) * bnb(i-2,1) + mat(1,i,5,3,1) * bnb(i-2,2)
            mat(1,i,:,2:3,1) = matl

            matl(1,1) = mat(1,i,1,2,2) * bnb(i-1,2) - mat(1,i,1,3,2) * bnb(i-1,1)
            matl(1,2) = mat(1,i,1,2,2) * bnb(i-1,1) + mat(1,i,1,3,2) * bnb(i-1,2)
            matl(2,1) = mat(1,i,2,2,2) * bnb(i-1,2) - mat(1,i,2,3,2) * bnb(i-1,1)
            matl(2,2) = mat(1,i,2,2,2) * bnb(i-1,1) + mat(1,i,2,3,2) * bnb(i-1,2)
            matl(3,1) = mat(1,i,3,2,2) * bnb(i-1,2) - mat(1,i,3,3,2) * bnb(i-1,1)
            matl(3,2) = mat(1,i,3,2,2) * bnb(i-1,1) + mat(1,i,3,3,2) * bnb(i-1,2)
            matl(4,1) = mat(1,i,4,2,2) * bnb(i-1,2) - mat(1,i,4,3,2) * bnb(i-1,1)
            matl(4,2) = mat(1,i,4,2,2) * bnb(i-1,1) + mat(1,i,4,3,2) * bnb(i-1,2)
            matl(5,1) = mat(1,i,5,2,2) * bnb(i-1,2) - mat(1,i,5,3,2) * bnb(i-1,1)
            matl(5,2) = mat(1,i,5,2,2) * bnb(i-1,1) + mat(1,i,5,3,2) * bnb(i-1,2)
            mat(1,i,:,2:3,2) = matl

            matl(1,1) = mat(1,i,1,2,3) * bnb(i,2) - mat(1,i,1,3,3) * bnb(i,1)
            matl(1,2) = mat(1,i,1,2,3) * bnb(i,1) + mat(1,i,1,3,3) * bnb(i,2)
            matl(2,1) = mat(1,i,2,2,3) * bnb(i,2) - mat(1,i,2,3,3) * bnb(i,1)
            matl(2,2) = mat(1,i,2,2,3) * bnb(i,1) + mat(1,i,2,3,3) * bnb(i,2)
            matl(3,1) = mat(1,i,3,2,3) * bnb(i,2) - mat(1,i,3,3,3) * bnb(i,1)
            matl(3,2) = mat(1,i,3,2,3) * bnb(i,1) + mat(1,i,3,3,3) * bnb(i,2)
            matl(4,1) = mat(1,i,4,2,3) * bnb(i,2) - mat(1,i,4,3,3) * bnb(i,1)
            matl(4,2) = mat(1,i,4,2,3) * bnb(i,1) + mat(1,i,4,3,3) * bnb(i,2)
            matl(5,1) = mat(1,i,5,2,3) * bnb(i,2) - mat(1,i,5,3,3) * bnb(i,1)
            matl(5,2) = mat(1,i,5,2,3) * bnb(i,1) + mat(1,i,5,3,3) * bnb(i,2)
            mat(1,i,:,2:3,3) = matl
            
            matl(1,1) = mat(1,i,1,2,4) * bnb(i+1,2) - mat(1,i,1,3,4) * bnb(i+1,1)
            matl(1,2) = mat(1,i,1,2,4) * bnb(i+1,1) + mat(1,i,1,3,4) * bnb(i+1,2)
            matl(2,1) = mat(1,i,2,2,4) * bnb(i+1,2) - mat(1,i,2,3,4) * bnb(i+1,1)
            matl(2,2) = mat(1,i,2,2,4) * bnb(i+1,1) + mat(1,i,2,3,4) * bnb(i+1,2)
            matl(3,1) = mat(1,i,3,2,4) * bnb(i+1,2) - mat(1,i,3,3,4) * bnb(i+1,1)
            matl(3,2) = mat(1,i,3,2,4) * bnb(i+1,1) + mat(1,i,3,3,4) * bnb(i+1,2)
            matl(4,1) = mat(1,i,4,2,4) * bnb(i+1,2) - mat(1,i,4,3,4) * bnb(i+1,1)
            matl(4,2) = mat(1,i,4,2,4) * bnb(i+1,1) + mat(1,i,4,3,4) * bnb(i+1,2)
            matl(5,1) = mat(1,i,5,2,4) * bnb(i+1,2) - mat(1,i,5,3,4) * bnb(i+1,1)
            matl(5,2) = mat(1,i,5,2,4) * bnb(i+1,1) + mat(1,i,5,3,4) * bnb(i+1,2)
            mat(1,i,:,2:3,4) = matl

            matl(1,1) = mat(1,i,1,2,5) * bnb(i+2,2) - mat(1,i,1,3,5) * bnb(i+2,1)
            matl(1,2) = mat(1,i,1,2,5) * bnb(i+2,1) + mat(1,i,1,3,5) * bnb(i+2,2)
            matl(2,1) = mat(1,i,2,2,5) * bnb(i+2,2) - mat(1,i,2,3,5) * bnb(i+2,1)
            matl(2,2) = mat(1,i,2,2,5) * bnb(i+2,1) + mat(1,i,2,3,5) * bnb(i+2,2)
            matl(3,1) = mat(1,i,3,2,5) * bnb(i+2,2) - mat(1,i,3,3,5) * bnb(i+2,1)
            matl(3,2) = mat(1,i,3,2,5) * bnb(i+2,1) + mat(1,i,3,3,5) * bnb(i+2,2)
            matl(4,1) = mat(1,i,4,2,5) * bnb(i+2,2) - mat(1,i,4,3,5) * bnb(i+2,1)
            matl(4,2) = mat(1,i,4,2,5) * bnb(i+2,1) + mat(1,i,4,3,5) * bnb(i+2,2)
            matl(5,1) = mat(1,i,5,2,5) * bnb(i+2,2) - mat(1,i,5,3,5) * bnb(i+2,1)
            matl(5,2) = mat(1,i,5,2,5) * bnb(i+2,1) + mat(1,i,5,3,5) * bnb(i+2,2)
            mat(1,i,:,2:3,5) = matl

            end do

            i = nx - 1

            matl(1,1) = mat(1,i,1,2,5) * bnb(i-3,2) - mat(1,i,1,3,5) * bnb(i-3,1)
            matl(1,2) = mat(1,i,1,2,5) * bnb(i-3,1) + mat(1,i,1,3,5) * bnb(i-3,2)
            matl(2,1) = mat(1,i,2,2,5) * bnb(i-3,2) - mat(1,i,2,3,5) * bnb(i-3,1)
            matl(2,2) = mat(1,i,2,2,5) * bnb(i-3,1) + mat(1,i,2,3,5) * bnb(i-3,2)
            matl(3,1) = mat(1,i,3,2,5) * bnb(i-3,2) - mat(1,i,3,3,5) * bnb(i-3,1)
            matl(3,2) = mat(1,i,3,2,5) * bnb(i-3,1) + mat(1,i,3,3,5) * bnb(i-3,2)
            matl(4,1) = mat(1,i,4,2,5) * bnb(i-3,2) - mat(1,i,4,3,5) * bnb(i-3,1)
            matl(4,2) = mat(1,i,4,2,5) * bnb(i-3,1) + mat(1,i,4,3,5) * bnb(i-3,2)
            matl(5,1) = mat(1,i,5,2,5) * bnb(i-3,2) - mat(1,i,5,3,5) * bnb(i-3,1)
            matl(5,2) = mat(1,i,5,2,5) * bnb(i-3,1) + mat(1,i,5,3,5) * bnb(i-3,2)
            mat(1,i,:,2:3,5) = matl

            matl(1,1) = mat(1,i,1,2,1) * bnb(i-2,2) - mat(1,i,1,3,1) * bnb(i-2,1)
            matl(1,2) = mat(1,i,1,2,1) * bnb(i-2,1) + mat(1,i,1,3,1) * bnb(i-2,2)
            matl(2,1) = mat(1,i,2,2,1) * bnb(i-2,2) - mat(1,i,2,3,1) * bnb(i-2,1)
            matl(2,2) = mat(1,i,2,2,1) * bnb(i-2,1) + mat(1,i,2,3,1) * bnb(i-2,2)
            matl(3,1) = mat(1,i,3,2,1) * bnb(i-2,2) - mat(1,i,3,3,1) * bnb(i-2,1)
            matl(3,2) = mat(1,i,3,2,1) * bnb(i-2,1) + mat(1,i,3,3,1) * bnb(i-2,2)
            matl(4,1) = mat(1,i,4,2,1) * bnb(i-2,2) - mat(1,i,4,3,1) * bnb(i-2,1)
            matl(4,2) = mat(1,i,4,2,1) * bnb(i-2,1) + mat(1,i,4,3,1) * bnb(i-2,2)
            matl(5,1) = mat(1,i,5,2,1) * bnb(i-2,2) - mat(1,i,5,3,1) * bnb(i-2,1)
            matl(5,2) = mat(1,i,5,2,1) * bnb(i-2,1) + mat(1,i,5,3,1) * bnb(i-2,2)
            mat(1,i,:,2:3,1) = matl

            matl(1,1) = mat(1,i,1,2,2) * bnb(i-1,2) - mat(1,i,1,3,2) * bnb(i-1,1)
            matl(1,2) = mat(1,i,1,2,2) * bnb(i-1,1) + mat(1,i,1,3,2) * bnb(i-1,2)
            matl(2,1) = mat(1,i,2,2,2) * bnb(i-1,2) - mat(1,i,2,3,2) * bnb(i-1,1)
            matl(2,2) = mat(1,i,2,2,2) * bnb(i-1,1) + mat(1,i,2,3,2) * bnb(i-1,2)
            matl(3,1) = mat(1,i,3,2,2) * bnb(i-1,2) - mat(1,i,3,3,2) * bnb(i-1,1)
            matl(3,2) = mat(1,i,3,2,2) * bnb(i-1,1) + mat(1,i,3,3,2) * bnb(i-1,2)
            matl(4,1) = mat(1,i,4,2,2) * bnb(i-1,2) - mat(1,i,4,3,2) * bnb(i-1,1)
            matl(4,2) = mat(1,i,4,2,2) * bnb(i-1,1) + mat(1,i,4,3,2) * bnb(i-1,2)
            matl(5,1) = mat(1,i,5,2,2) * bnb(i-1,2) - mat(1,i,5,3,2) * bnb(i-1,1)
            matl(5,2) = mat(1,i,5,2,2) * bnb(i-1,1) + mat(1,i,5,3,2) * bnb(i-1,2)
            mat(1,i,:,2:3,2) = matl

            matl(1,1) = mat(1,i,1,2,3) * bnb(i,2) - mat(1,i,1,3,3) * bnb(i,1)
            matl(1,2) = mat(1,i,1,2,3) * bnb(i,1) + mat(1,i,1,3,3) * bnb(i,2)
            matl(2,1) = mat(1,i,2,2,3) * bnb(i,2) - mat(1,i,2,3,3) * bnb(i,1)
            matl(2,2) = mat(1,i,2,2,3) * bnb(i,1) + mat(1,i,2,3,3) * bnb(i,2)
            matl(3,1) = mat(1,i,3,2,3) * bnb(i,2) - mat(1,i,3,3,3) * bnb(i,1)
            matl(3,2) = mat(1,i,3,2,3) * bnb(i,1) + mat(1,i,3,3,3) * bnb(i,2)
            matl(4,1) = mat(1,i,4,2,3) * bnb(i,2) - mat(1,i,4,3,3) * bnb(i,1)
            matl(4,2) = mat(1,i,4,2,3) * bnb(i,1) + mat(1,i,4,3,3) * bnb(i,2)
            matl(5,1) = mat(1,i,5,2,3) * bnb(i,2) - mat(1,i,5,3,3) * bnb(i,1)
            matl(5,2) = mat(1,i,5,2,3) * bnb(i,1) + mat(1,i,5,3,3) * bnb(i,2)
            mat(1,i,:,2:3,3) = matl
            
            matl(1,1) = mat(1,i,1,2,4) * bnb(i+1,2) - mat(1,i,1,3,4) * bnb(i+1,1)
            matl(1,2) = mat(1,i,1,2,4) * bnb(i+1,1) + mat(1,i,1,3,4) * bnb(i+1,2)
            matl(2,1) = mat(1,i,2,2,4) * bnb(i+1,2) - mat(1,i,2,3,4) * bnb(i+1,1)
            matl(2,2) = mat(1,i,2,2,4) * bnb(i+1,1) + mat(1,i,2,3,4) * bnb(i+1,2)
            matl(3,1) = mat(1,i,3,2,4) * bnb(i+1,2) - mat(1,i,3,3,4) * bnb(i+1,1)
            matl(3,2) = mat(1,i,3,2,4) * bnb(i+1,1) + mat(1,i,3,3,4) * bnb(i+1,2)
            matl(4,1) = mat(1,i,4,2,4) * bnb(i+1,2) - mat(1,i,4,3,4) * bnb(i+1,1)
            matl(4,2) = mat(1,i,4,2,4) * bnb(i+1,1) + mat(1,i,4,3,4) * bnb(i+1,2)
            matl(5,1) = mat(1,i,5,2,4) * bnb(i+1,2) - mat(1,i,5,3,4) * bnb(i+1,1)
            matl(5,2) = mat(1,i,5,2,4) * bnb(i+1,1) + mat(1,i,5,3,4) * bnb(i+1,2)
            mat(1,i,:,2:3,4) = matl

            i = nx

            matl(1,1) = mat(1,i,1,2,4) * bnb(i-4,2) - mat(1,i,1,3,4) * bnb(i-4,1)
            matl(1,2) = mat(1,i,1,2,4) * bnb(i-4,1) + mat(1,i,1,3,4) * bnb(i-4,2)
            matl(2,1) = mat(1,i,2,2,4) * bnb(i-4,2) - mat(1,i,2,3,4) * bnb(i-4,1)
            matl(2,2) = mat(1,i,2,2,4) * bnb(i-4,1) + mat(1,i,2,3,4) * bnb(i-4,2)
            matl(3,1) = mat(1,i,3,2,4) * bnb(i-4,2) - mat(1,i,3,3,4) * bnb(i-4,1)
            matl(3,2) = mat(1,i,3,2,4) * bnb(i-4,1) + mat(1,i,3,3,4) * bnb(i-4,2)
            matl(4,1) = mat(1,i,4,2,4) * bnb(i-4,2) - mat(1,i,4,3,4) * bnb(i-4,1)
            matl(4,2) = mat(1,i,4,2,4) * bnb(i-4,1) + mat(1,i,4,3,4) * bnb(i-4,2)
            matl(5,1) = mat(1,i,5,2,4) * bnb(i-4,2) - mat(1,i,5,3,4) * bnb(i-4,1)
            matl(5,2) = mat(1,i,5,2,4) * bnb(i-4,1) + mat(1,i,5,3,4) * bnb(i-4,2)
            mat(1,i,:,2:3,4) = matl

            matl(1,1) = mat(1,i,1,2,5) * bnb(i-3,2) - mat(1,i,1,3,5) * bnb(i-3,1)
            matl(1,2) = mat(1,i,1,2,5) * bnb(i-3,1) + mat(1,i,1,3,5) * bnb(i-3,2)
            matl(2,1) = mat(1,i,2,2,5) * bnb(i-3,2) - mat(1,i,2,3,5) * bnb(i-3,1)
            matl(2,2) = mat(1,i,2,2,5) * bnb(i-3,1) + mat(1,i,2,3,5) * bnb(i-3,2)
            matl(3,1) = mat(1,i,3,2,5) * bnb(i-3,2) - mat(1,i,3,3,5) * bnb(i-3,1)
            matl(3,2) = mat(1,i,3,2,5) * bnb(i-3,1) + mat(1,i,3,3,5) * bnb(i-3,2)
            matl(4,1) = mat(1,i,4,2,5) * bnb(i-3,2) - mat(1,i,4,3,5) * bnb(i-3,1)
            matl(4,2) = mat(1,i,4,2,5) * bnb(i-3,1) + mat(1,i,4,3,5) * bnb(i-3,2)
            matl(5,1) = mat(1,i,5,2,5) * bnb(i-3,2) - mat(1,i,5,3,5) * bnb(i-3,1)
            matl(5,2) = mat(1,i,5,2,5) * bnb(i-3,1) + mat(1,i,5,3,5) * bnb(i-3,2)
            mat(1,i,:,2:3,5) = matl

            matl(1,1) = mat(1,i,1,2,1) * bnb(i-2,2) - mat(1,i,1,3,1) * bnb(i-2,1)
            matl(1,2) = mat(1,i,1,2,1) * bnb(i-2,1) + mat(1,i,1,3,1) * bnb(i-2,2)
            matl(2,1) = mat(1,i,2,2,1) * bnb(i-2,2) - mat(1,i,2,3,1) * bnb(i-2,1)
            matl(2,2) = mat(1,i,2,2,1) * bnb(i-2,1) + mat(1,i,2,3,1) * bnb(i-2,2)
            matl(3,1) = mat(1,i,3,2,1) * bnb(i-2,2) - mat(1,i,3,3,1) * bnb(i-2,1)
            matl(3,2) = mat(1,i,3,2,1) * bnb(i-2,1) + mat(1,i,3,3,1) * bnb(i-2,2)
            matl(4,1) = mat(1,i,4,2,1) * bnb(i-2,2) - mat(1,i,4,3,1) * bnb(i-2,1)
            matl(4,2) = mat(1,i,4,2,1) * bnb(i-2,1) + mat(1,i,4,3,1) * bnb(i-2,2)
            matl(5,1) = mat(1,i,5,2,1) * bnb(i-2,2) - mat(1,i,5,3,1) * bnb(i-2,1)
            matl(5,2) = mat(1,i,5,2,1) * bnb(i-2,1) + mat(1,i,5,3,1) * bnb(i-2,2)
            mat(1,i,:,2:3,1) = matl

            matl(1,1) = mat(1,i,1,2,2) * bnb(i-1,2) - mat(1,i,1,3,2) * bnb(i-1,1)
            matl(1,2) = mat(1,i,1,2,2) * bnb(i-1,1) + mat(1,i,1,3,2) * bnb(i-1,2)
            matl(2,1) = mat(1,i,2,2,2) * bnb(i-1,2) - mat(1,i,2,3,2) * bnb(i-1,1)
            matl(2,2) = mat(1,i,2,2,2) * bnb(i-1,1) + mat(1,i,2,3,2) * bnb(i-1,2)
            matl(3,1) = mat(1,i,3,2,2) * bnb(i-1,2) - mat(1,i,3,3,2) * bnb(i-1,1)
            matl(3,2) = mat(1,i,3,2,2) * bnb(i-1,1) + mat(1,i,3,3,2) * bnb(i-1,2)
            matl(4,1) = mat(1,i,4,2,2) * bnb(i-1,2) - mat(1,i,4,3,2) * bnb(i-1,1)
            matl(4,2) = mat(1,i,4,2,2) * bnb(i-1,1) + mat(1,i,4,3,2) * bnb(i-1,2)
            matl(5,1) = mat(1,i,5,2,2) * bnb(i-1,2) - mat(1,i,5,3,2) * bnb(i-1,1)
            matl(5,2) = mat(1,i,5,2,2) * bnb(i-1,1) + mat(1,i,5,3,2) * bnb(i-1,2)
            mat(1,i,:,2:3,2) = matl

            matl(1,1) = mat(1,i,1,2,3) * bnb(i,2) - mat(1,i,1,3,3) * bnb(i,1)
            matl(1,2) = mat(1,i,1,2,3) * bnb(i,1) + mat(1,i,1,3,3) * bnb(i,2)
            matl(2,1) = mat(1,i,2,2,3) * bnb(i,2) - mat(1,i,2,3,3) * bnb(i,1)
            matl(2,2) = mat(1,i,2,2,3) * bnb(i,1) + mat(1,i,2,3,3) * bnb(i,2)
            matl(3,1) = mat(1,i,3,2,3) * bnb(i,2) - mat(1,i,3,3,3) * bnb(i,1)
            matl(3,2) = mat(1,i,3,2,3) * bnb(i,1) + mat(1,i,3,3,3) * bnb(i,2)
            matl(4,1) = mat(1,i,4,2,3) * bnb(i,2) - mat(1,i,4,3,3) * bnb(i,1)
            matl(4,2) = mat(1,i,4,2,3) * bnb(i,1) + mat(1,i,4,3,3) * bnb(i,2)
            matl(5,1) = mat(1,i,5,2,3) * bnb(i,2) - mat(1,i,5,3,3) * bnb(i,1)
            matl(5,2) = mat(1,i,5,2,3) * bnb(i,1) + mat(1,i,5,3,3) * bnb(i,2)
            mat(1,i,:,2:3,3) = matl
            
!.... now rotate the equations to the normal direction

            do i = 1, nx
              mat(1,i,2,:,:) = bnb(i,2) * mat(1,i,2,:,:) - &
                               bnb(i,1) * mat(1,i,3,:,:)
              mat(1,i,3,:,:) = bnb(i,1) * mat(1,i,2,:,:) + &
                               bnb(i,2) * mat(1,i,3,:,:)
            end do

!.... constrain the normal velocity to be zero

            mat(1,:,3,:,:) = zero
            mat(1,:,3,3,3) = one

        end if                  ! Navier
          
!.... Riemann Invariants on top boundary

        if (top.eq.0) then
          mat(ny,:,:,:,:) = zero
          mat(ny,:,1,1,3) = one
          mat(ny,:,2,2,3) = one
          mat(ny,:,3,3,3) = one
          mat(ny,:,4,4,3) = one
          mat(ny,:,5,5,3) = one
        end if
        
        end if          ! yper

!=============================================================================!
!       L e f t   B o u n d a r y
!=============================================================================!
        if (.not. xper) then

!.... constant pressure

!       mat(:,1,1,:,:) = zero
!       mat(:,1,1,1,3) = one

!.... symmetry plane

        if (left.eq.2) then
          mat(:,1,1,:,:) = zero
          mat(:,1,1,1,3) = -gc1
          mat(:,1,1,1,4) = -gc2
          mat(:,1,1,1,5) = -gc3
          mat(:,1,1,1,1) = -gc4
          mat(:,1,1,1,2) = -gc5

          mat(:,1,2,:,:) = zero
          mat(:,1,2,2,3) = -gc1
          mat(:,1,2,2,4) = -gc2
          mat(:,1,2,2,5) = -gc3
          mat(:,1,2,2,1) = -gc4
          mat(:,1,2,2,2) = -gc5

          mat(:,1,3,:,:) = zero
          mat(:,1,3,3,3) = one

          mat(:,1,4,:,:) = zero
          mat(:,1,4,4,3) = -gc1
          mat(:,1,4,4,4) = -gc2
          mat(:,1,4,4,5) = -gc3
          mat(:,1,4,4,1) = -gc4
          mat(:,1,4,4,2) = -gc5

          mat(:,1,5,:,:) = zero
          mat(:,1,5,5,3) = -gc1
          mat(:,1,5,5,4) = -gc2
          mat(:,1,5,5,5) = -gc3
          mat(:,1,5,5,1) = -gc4
          mat(:,1,5,5,2) = -gc5
        end if

        if (left.eq.7) then             ! symmetry plane
          mat(:,1,3,:,:) = zero
          mat(:,1,3,3,3) = one
        endif

!.... Left Riemann invariant boundary condition (first-order)
        
        if (left.eq.0 .and. Ma.lt.one) then
          mat(:,1,:,:,:) = zero
          call ReimannLHS( ny, vl(:,1,:), vl(:,2,:), vl(:,3,:), &
                           bnl, x(:,1), y(:,1), mat(:,1,:,:,3), &
                           mat(:,1,:,:,4), mat(:,1,:,:,5), &
                           rhobl, ubl, vbl, wbl, tbl, pbl, cbl )
        end if

!.... Left Side:  Riemann invariants (zero-order)

          if (.false.) then
          
          mat(:,1,:,:,:) = zero

          R1(:) = one - two / Ma / gamma1
          R2(:) = vl(:,2,2) + &
                  two * sqrt( vl(:,2,5) ) / Ma / gamma1
                
          rho(:) = vl(:,2,1)
          p(:)   = vl(:,2,1) * vl(:,2,5) / (gamma * Ma**2)

          fact1(:) = (rho(:)**gamma/(gamma*p(:)))**(one/gamma1) * &
                     (gamma1/four)**(two/gamma1) * (two/gamma1) * &
                     (R2(:) - R1(:))**((two-gamma1)/gamma1)
          fact2(:) = (gamma1/four * (R2(:) - R1(:)))**(two/gamma1) * &
                     one/(gamma*gamma1) * &
                     (rho(:)**gamma/(gamma*p(:)))**((two-gamma)/gamma1)
          mat(:,1,1,1,3) = -one
          mat(:,1,1,2,4) = fact1(:) + &
                           gamma1 * rho(:)**gamma1 / p(:)
          mat(:,1,1,5,4) = fact1(:) / ( gamma1 * Ma * sqrt(vl(:,2,5)) ) - &
                           fact2(:) * rho(:)**(gamma+1)/(gamma*Ma**2*p(:)**2)
          
          mat(:,1,2,2,3) = -one
          mat(:,1,2,2,4) = pt5
          mat(:,1,2,5,4) = pt5 / ( gamma1 * Ma * sqrt(vl(:,2,5)) )
  
          mat(:,1,3,3,3) = -one
          mat(:,1,3,3,4) = one
  
          mat(:,1,4,4,3) = -one
          mat(:,1,4,4,4) = one
  
          fact(:) = (Ma * gamma1 / four)**2 * two * ( R2(:) - R1(:) )
          mat(:,1,5,5,3) = -one
          mat(:,1,5,2,4) = fact(:)
          mat(:,1,5,3,4) = zero
          mat(:,1,5,5,4) = fact(:) / ( gamma1 * Ma * sqrt(vl(:,2,5)) )
          
          end if                ! Riemann boundary
        
          if (.false.) then

!.... Zero-order extrapolation in the viscous layer

          if (extrap.eq.0) then
            mat(1:nbl,1,:,:,:) = zero
  
            mat(1:nbl,1,1,1,3) = -one
            mat(1:nbl,1,2,2,3) = -one
            mat(1:nbl,1,3,3,3) = -one
            mat(1:nbl,1,4,4,3) = -one
            mat(1:nbl,1,5,5,3) = -one
  
            mat(1:nbl,1,1,1,4) = one
            mat(1:nbl,1,2,2,4) = one
            mat(1:nbl,1,3,3,4) = one
            mat(1:nbl,1,4,4,4) = one
            mat(1:nbl,1,5,5,4) = one
          end if
          
!.... First-order extrapolation in the viscous layer

          if (extrap.eq.1) then
            mat(1:nbl,1,:,:,:) = zero
  
            mat(1:nbl,1,1,1,3) = -one
            mat(1:nbl,1,2,2,3) = -one
            mat(1:nbl,1,3,3,3) = -one
            mat(1:nbl,1,4,4,3) = -one
            mat(1:nbl,1,5,5,3) = -one
  
            mat(1:nbl,1,1,1,4) = two
            mat(1:nbl,1,2,2,4) = two
            mat(1:nbl,1,3,3,4) = two
            mat(1:nbl,1,4,4,4) = two
            mat(1:nbl,1,5,5,4) = two
            
            mat(1:nbl,1,1,1,5) = -one
            mat(1:nbl,1,2,2,5) = -one
            mat(1:nbl,1,3,3,5) = -one
            mat(1:nbl,1,4,4,5) = -one
            mat(1:nbl,1,5,5,5) = -one
          end if                ! extrapolation type

          end if
!=============================================================================!
!       R i g h t  B o u n d a r y
!=============================================================================!

!.... constant pressure

!       mat(:,nx,1,:,:) = zero
!       mat(:,nx,1,1,3) = one

!.... Right Riemann invariant boundary condition (first-order)
        
        if (right.eq.0 .and. Ma.lt.one) then

          mat(:,nx,:,:,:) = zero
          call ReimannLHS( ny, vl(:,nx,:), vl(:,nx-1,:), vl(:,nx-2,:), &
                           bnr, x(:,nx), y(:,nx), mat(:,nx,:,:,3), &
                           mat(:,nx,:,:,2), mat(:,nx,:,:,1), &
                           rhobr, ubr, vbr, wbr, tbr, pbr, cbr )
                
        end if

!.... Right Side:  Riemann invariants (zero-order)

          if (.false.) then
          
          mat(:,nx,:,:,:) = zero

          R1(:) = one - two / Ma / gamma1
          R2(:) = vl(:,nx-1,2) + &
                  two * sqrt( vl(:,nx-1,5) ) / Ma / gamma1
                
          rho(:) = vl(:,nx-1,1)
          p(:)   = vl(:,nx-1,1) * vl(:,nx-1,5) / (gamma * Ma**2)

          fact1(:) = (rho(:)**gamma/(gamma*p(:)))**(one/gamma1) * &
                     (gamma1/four)**(two/gamma1) * (two/gamma1) * &
                     (R2(:) - R1(:))**((two-gamma1)/gamma1)
          fact2(:) = (gamma1/four * (R2(:) - R1(:)))**(two/gamma1) * &
                     one/(gamma*gamma1) * &
                     (rho(:)**gamma/(gamma*p(:)))**((two-gamma)/gamma1)
          mat(:,nx,1,1,3) = -one
          mat(:,nx,1,2,2) = fact1(:) + &
                            gamma1 * rho(:)**gamma1 / p(:)
          mat(:,nx,1,5,2) = fact1(:) / ( gamma1 * Ma * sqrt(vl(:,2,5)) ) - &
                            fact2(:) * rho(:)**(gamma+1)/(gamma*Ma**2*p(:)**2)
          
          mat(:,nx,2,2,3) = -one
          mat(:,nx,2,2,2) = pt5
          mat(:,nx,2,5,2) = pt5 / ( gamma1 * Ma * sqrt(vl(:,2,5)) )
  
          mat(:,nx,3,3,3) = -one
          mat(:,nx,3,3,2) = one
  
          mat(:,nx,4,4,3) = -one
          mat(:,nx,4,4,2) = one
  
          fact(:) = (Ma * gamma1 / four)**2 * two * ( R2(:) - R1(:) )
          mat(:,nx,5,5,3) = -one
          mat(:,nx,5,2,2) = fact(:)
          mat(:,nx,5,3,2) = zero
          mat(:,nx,5,5,2) = fact(:) / ( gamma1 * Ma * sqrt(vl(:,2,5)) )

          end if                ! Riemann boundary
          
!.... Zero'th order extrapolation in the viscous layer

          if (extrap.eq.0) then
            mat(1:nbl,nx,:,:,:) = zero
            
            mat(1:nbl,nx,1,1,3) = -one
            mat(1:nbl,nx,2,2,3) = -one
            mat(1:nbl,nx,3,3,3) = -one
            mat(1:nbl,nx,4,4,3) = -one
            mat(1:nbl,nx,5,5,3) = -one
  
            mat(1:nbl,nx,1,1,2) = one
            mat(1:nbl,nx,2,2,2) = one
            mat(1:nbl,nx,3,3,2) = one
            mat(1:nbl,nx,4,4,2) = one
            mat(1:nbl,nx,5,5,2) = one
          end if
          
!.... First-order extrapolation in the viscous layer

          if (extrap.eq.1) then
            mat(1:nbl,nx,:,:,:) = zero
  
            mat(1:nbl,nx,1,1,3) = -one
            mat(1:nbl,nx,2,2,3) = -one
            mat(1:nbl,nx,3,3,3) = -one
            mat(1:nbl,nx,4,4,3) = -one
            mat(1:nbl,nx,5,5,3) = -one
  
            mat(1:nbl,nx,1,1,2) = two
            mat(1:nbl,nx,2,2,2) = two
            mat(1:nbl,nx,3,3,2) = two
            mat(1:nbl,nx,4,4,2) = two
            mat(1:nbl,nx,5,5,2) = two
            
            mat(1:nbl,nx,1,1,1) = -one
            mat(1:nbl,nx,2,2,1) = -one
            mat(1:nbl,nx,3,3,1) = -one
            mat(1:nbl,nx,4,4,1) = -one
            mat(1:nbl,nx,5,5,1) = -one
          end if                ! extrapolation type

!.... hold the IC

          if (right.eq.8) then
            mat(:,nx,:,:,:) = zero
            mat(:,nx,1,1,3) = -one
            mat(:,nx,2,2,3) = -one
            mat(:,nx,3,3,3) = -one
            mat(:,nx,4,4,3) = -one
            mat(:,nx,5,5,3) = -one
          end if

        end if                  ! xper
        
        end if                  ! linear
        
        return
        end
