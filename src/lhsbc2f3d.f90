!=============================================================================!
        subroutine lhsbc2f3d( mat, vl, vml)
!  
!  Correct the LHS for boundary conditions in the \eta direction
!
!  This version supports the 4th order LHS
!
!  Revised: 10-17-95
!=============================================================================!
        use global
        use stencil
        implicit none
        
        complex :: mat(5,ndof,ndof,nx,ny), vl(ndof,nx,ny)
        real    :: vml(ndof,nx,ny)

        complex :: matl(5,2)
        
        integer :: i, j

!=============================================================================!
!       L I N E A R   B O U N D A R Y   C O N D I T I O N S
!=============================================================================!
        if (.not. yper) then

        if (Navier) then
        
!.... wall boundary condition

        if (wall.eq.1) then
          mat(:,1,:,is:ie,1) = zero
          mat(3,1,1,is:ie,1) = -gc1 
          mat(4,1,1,is:ie,1) = -gc2
          mat(5,1,1,is:ie,1) = -gc3
          mat(1,1,1,is:ie,1) = -gc4
          mat(2,1,1,is:ie,1) = -gc5
        end if

!.... Third-order extrapolation of rho at the wall

        if (wall.eq.3) then
          mat(:,1,:,is:ie,1) =  zero
          mat(3,1,1,is:ie,1) = -one
          mat(4,1,1,is:ie,1) =  4.0
          mat(5,1,1,is:ie,1) = -6.0
          mat(1,1,1,is:ie,1) =  4.0
          mat(2,1,1,is:ie,1) = -one
        end if

!.... no-slip

        mat(:,2:4,:,is:ie,1) = zero
        mat(3,2,2,is:ie,1)   = one
        mat(3,3,3,is:ie,1)   = one
        mat(3,4,4,is:ie,1)   = one

!.... isothermal

        if (wallt.eq.0) then
          mat(:,5,:,is:ie,1) = zero
          mat(3,5,5,is:ie,1) = one
        end if

!.... adiabatic

        if (wallt.eq.1) then
          mat(:,5,:,is:ie,1) = zero
          mat(3,5,5,is:ie,1) = -gc1
          mat(4,5,5,is:ie,1) = -gc2
          mat(5,5,5,is:ie,1) = -gc3
          mat(1,5,5,is:ie,1) = -gc4
          mat(2,5,5,is:ie,1) = -gc5
        end if

        else            ! inviscid:  rotate to body normal coordinates
        
          !call warning('lhsbc2f3d$','Inviscid BCs need to have is:ie$')

          !do i = 1, nx
          do i = is, ie 
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
            
            matl(1,1) = mat(2,1,2,i,2) * bnb(i,2) - mat(2,1,3,i,2) * bnb(i,1)
            matl(1,2) = mat(2,1,2,i,2) * bnb(i,1) + mat(2,1,3,i,2) * bnb(i,2)
            matl(2,1) = mat(2,2,2,i,2) * bnb(i,2) - mat(2,2,3,i,2) * bnb(i,1)
            matl(2,2) = mat(2,2,2,i,2) * bnb(i,1) + mat(2,2,3,i,2) * bnb(i,2)
            matl(3,1) = mat(2,3,2,i,2) * bnb(i,2) - mat(2,3,3,i,2) * bnb(i,1)
            matl(3,2) = mat(2,3,2,i,2) * bnb(i,1) + mat(2,3,3,i,2) * bnb(i,2)
            matl(4,1) = mat(2,4,2,i,2) * bnb(i,2) - mat(2,4,3,i,2) * bnb(i,1)
            matl(4,2) = mat(2,4,2,i,2) * bnb(i,1) + mat(2,4,3,i,2) * bnb(i,2)
            matl(5,1) = mat(2,5,2,i,2) * bnb(i,2) - mat(2,5,3,i,2) * bnb(i,1)
            matl(5,2) = mat(2,5,2,i,2) * bnb(i,1) + mat(2,5,3,i,2) * bnb(i,2)
            mat(2,:,2:3,i,2) = matl
            
            matl(1,1) = mat(1,1,2,i,3) * bnb(i,2) - mat(1,1,3,i,3) * bnb(i,1)
            matl(1,2) = mat(1,1,2,i,3) * bnb(i,1) + mat(1,1,3,i,3) * bnb(i,2)
            matl(2,1) = mat(1,2,2,i,3) * bnb(i,2) - mat(1,2,3,i,3) * bnb(i,1)
            matl(2,2) = mat(1,2,2,i,3) * bnb(i,1) + mat(1,2,3,i,3) * bnb(i,2)
            matl(3,1) = mat(1,3,2,i,3) * bnb(i,2) - mat(1,3,3,i,3) * bnb(i,1)
            matl(3,2) = mat(1,3,2,i,3) * bnb(i,1) + mat(1,3,3,i,3) * bnb(i,2)
            matl(4,1) = mat(1,4,2,i,3) * bnb(i,2) - mat(1,4,3,i,3) * bnb(i,1)
            matl(4,2) = mat(1,4,2,i,3) * bnb(i,1) + mat(1,4,3,i,3) * bnb(i,2)
            matl(5,1) = mat(1,5,2,i,3) * bnb(i,2) - mat(1,5,3,i,3) * bnb(i,1)
            matl(5,2) = mat(1,5,2,i,3) * bnb(i,1) + mat(1,5,3,i,3) * bnb(i,2)
            mat(1,:,2:3,i,3) = matl

            mat(:,2,:,i,1) = bnb(i,2) * mat(:,2,:,i,1) - &
                             bnb(i,1) * mat(:,3,:,i,1)
            mat(:,3,:,i,1) = bnb(i,1) * mat(:,2,:,i,1) + &
                             bnb(i,2) * mat(:,3,:,i,1)

!.... constrain the normal velocity to be zero

            mat(:,3,:,i,1) = zero
            mat(3,3,3,i,1) = one

          end do

        end if          ! Navier

!.... freestream boundary conditions

        if (top.eq.0 .or. top.eq.1) then
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

!.... Left boundary

          if (left.eq.0 .or. left.eq.4 .or. left.eq.5) then

            mat(:,:,:,1,:) = zero
            mat(3,1,1,1,:) = one
            mat(3,2,2,1,:) = one
            mat(3,3,3,1,:) = one
            mat(3,4,4,1,:) = one
            mat(3,5,5,1,:) = one

          end if

!.... Right boundary

          if (right.eq.0 .or. right.eq.4 .or. right.eq.5) then

            mat(:,:,:,nx,:) = zero
            mat(3,1,1,nx,:) = one
            mat(3,2,2,nx,:) = one
            mat(3,3,3,nx,:) = one
            mat(3,4,4,nx,:) = one
            mat(3,5,5,nx,:) = one
        
          end if

        end if          ! xper

        return
        end
