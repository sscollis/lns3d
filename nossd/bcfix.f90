!=============================================================================!
	subroutine bcfix( nx, ny, ndof, bnb, r )
	
	integer :: ny, nx, ndof
	real    :: bnb(nx,2)
	complex :: r(ny,nx,ndof)
	complex :: tmp(nx)

	integer :: i
!=============================================================================!

	do i = 1, nx
	  tmp(i)   =  bnb(i,2) * r(1,i,2) + bnb(i,1) * r(1,i,3)
	  r(1,i,3) = -bnb(i,1) * r(1,i,2) + bnb(i,2) * r(1,i,3)
	  r(1,i,2) = tmp(i)
	end do
	
	return
	end

!.... Gives wrong solution on San Diego?

!	tmp(:)   =  bnb(:,2) * r(1,:,2) + bnb(:,1) * r(1,:,3)
!	r(1,:,3) = -bnb(:,1) * r(1,:,2) + bnb(:,2) * r(1,:,3)
!	r(1,:,2) = tmp(:)


!=============================================================================!
	subroutine rbcfix( nx, ny, ndof, bnb, r )
	
	integer :: ny, nx, ndof
	real    :: bnb(nx,2)
	real    :: r(ny,nx,ndof)
	real    :: tmp(nx)

	integer :: i
!=============================================================================!

	do i = 1, nx
	  tmp(i)   =  bnb(i,2) * r(1,i,2) + bnb(i,1) * r(1,i,3)
	  r(1,i,3) = -bnb(i,1) * r(1,i,2) + bnb(i,2) * r(1,i,3)
	  r(1,i,2) = tmp(i)
	end do
	
	return
	end
