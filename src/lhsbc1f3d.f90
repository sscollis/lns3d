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

        if (wall.eq.1 .or. wall.eq.3) then
          mat(:,1,:,is:ie,1) = zero
          mat(3,1,1,is:ie,1) = one
        end if

!.... linear extrapolation of rho at the wall

        if (wall.eq.2) then
          mat(:,1,:,is:ie,1) = zero
          mat(3,1,1,is:ie,1) = one
        end if

!.... no-slip

        mat(:,2:4,:,is:ie,1) = zero
        mat(3,2,2,is:ie,1)   = one
        mat(3,3,3,is:ie,1)   = one
        mat(3,4,4,is:ie,1)   = one

!.... wall temperature boundary condition

        if (wallt.eq.0 .or. wallt.eq.1) then
          mat(:,5,:,is:ie,1) = zero
          mat(3,5,5,is:ie,1) = one
        end if
        
        else            ! inviscid:  rotate to body normal coordinates
        
            call warning('lhsbc1f3d$','Inviscid BCs need to have is:ie$')

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

!.... watch out for symetric boundaries and periodic boundaries

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

        end if          ! Navier
        
!.... freestream hard boundary condition

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

!.... left boundary

          if (left.eq.0 .or. left.eq.4 .or. left.eq.5) then
            mat(:,:,:,1,:) = zero
            mat(3,1,1,1,:) = one
            mat(3,2,2,1,:) = one
            mat(3,3,3,1,:) = one
            mat(3,4,4,1,:) = one
            mat(3,5,5,1,:) = one
         endif

!.... Right boundary

          if (right.eq.0 .or. right.eq.4 .or. right.eq.5) then
            mat(:,:,:,nx,:) = zero
            mat(3,1,1,nx,:) = one
            mat(3,2,2,nx,:) = one
            mat(3,3,3,nx,:) = one
            mat(3,4,4,nx,:) = one
            mat(3,5,5,nx,:) = one
          end if
          
!.... hold the IC

          if (right.eq.8) then
            mat(:,:,:,nx,:) = zero
            mat(3,1,1,nx,:) = -one
            mat(3,2,2,nx,:) = -one
            mat(3,3,3,nx,:) = -one
            mat(3,4,4,nx,:) = -one
            mat(3,5,5,nx,:) = -one
          end if

        end if          ! xper

        return
        end
