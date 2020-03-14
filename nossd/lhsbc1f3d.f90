!=============================================================================!
        subroutine lhsbc1f3d( mat, vl, vml)
!  
!  Correct the LHS for boundary conditions in the \xi direction.
!
!  This version supports the 4th order LHS. 
!
!  Revised: 10-16-95
!=============================================================================!
        use global
        use stencil
        implicit none
        
        complex :: mat(ny,nx,ndof,ndof,5), vl(ny,nx,ndof)
        real    :: vml(ny,nx,ndof)

        complex :: matl(5,2)

        integer :: i, j, ij
!=============================================================================!
!       L I N E A R   B O U N D A R Y   C O N D I T I O N S
!=============================================================================!
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
        
        else            ! inviscid:  rotate to body normal coordinates
        
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

!.... watch out for symetric boundaries and periodic boundaries

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

        end if          ! Navier
        
!.... freestream hard boundary condition

        if (top.eq.0 .or. top.eq.1) then
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

          if (left.eq.0 .or. left.eq.4 .or. left.eq.5) then
            mat(:,1,:,:,:) = zero
            mat(:,1,1,1,3) = one
            mat(:,1,2,2,3) = one
            mat(:,1,3,3,3) = one
            mat(:,1,4,4,3) = one
            mat(:,1,5,5,3) = one
         endif

!.... Right boundary

          if (right.eq.0 .or. right.eq.4 .or. right.eq.5) then
            mat(:,nx,:,:,:) = zero
            mat(:,nx,1,1,3) = one
            mat(:,nx,2,2,3) = one
            mat(:,nx,3,3,3) = one
            mat(:,nx,4,4,3) = one
            mat(:,nx,5,5,3) = one
          end if
          
!.... hold the IC

          if (right.eq.8) then
            mat(:,nx,:,:,:) = zero
            mat(:,nx,1,1,3) = -one
            mat(:,nx,2,2,3) = -one
            mat(:,nx,3,3,3) = -one
            mat(:,nx,4,4,3) = -one
            mat(:,nx,5,5,3) = -one
          end if

        end if          ! xper

        return
        end
