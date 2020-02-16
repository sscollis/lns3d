!=============================================================================!
        subroutine bcfix( nx, ny, ndof, bnb, r )
        
        integer :: ny, nx, ndof
        real    :: bnb(nx,2)
        complex :: r(ndof,nx,ny)
        complex :: tmp(nx)

        integer :: i
!=============================================================================!

        do i = 1, nx
          tmp(i)   =  bnb(i,2) * r(2,i,1) + bnb(i,1) * r(3,i,1)
          r(3,i,1) = -bnb(i,1) * r(2,i,1) + bnb(i,2) * r(3,i,1)
          r(2,i,1) = tmp(i)
        end do
        
        return
        end

!.... Gives wrong solution on San Diego?

!       tmp(:)   =  bnb(:,2) * r(2,:,1) + bnb(:,1) * r(3,:,1)
!       r(3,:,1) = -bnb(:,1) * r(2,:,1) + bnb(:,2) * r(3,:,1)
!       r(2,:,1) = tmp(:)


!=============================================================================!
        subroutine rbcfix( nx, ny, ndof, bnb, r )
        
        integer :: ny, nx, ndof
        real    :: bnb(nx,2)
        real    :: r(ndof,nx,ny)
        real    :: tmp(nx)

        integer :: i
!=============================================================================!

        do i = 1, nx
          tmp(i)   =  bnb(i,2) * r(2,i,1) + bnb(i,1) * r(3,i,1)
          r(3,i,1) = -bnb(i,1) * r(2,i,1) + bnb(i,2) * r(3,i,1)
          r(2,i,1) = tmp(i)
        end do
        
        return
        end
