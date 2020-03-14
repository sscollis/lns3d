!=============================================================================!
	subroutine csolve1bc( nsys, n, nblk, mat, rhs, bc, mult, fact )
!=============================================================================!
!  
!  Solves a block tridiagonal system of equations without pivoting. 
!  Vectorized over the first index.
!
!  This routine has been modified to allow for high-order boundary
!  conditions.
!
!  WARNING: Since pivoting is not done, this algorithm could fail
!
!  WARNING: This routine works for nblk = 5 only!
!
!  input:
!    nsys  	: number of systems to vectorize over
!    n 	   	: number of rows in the systems
!    nblk  	: dimension of the blocks in each row
!
!  inout:
!    mat	: the block tridiagonal matrix to be solved
!    rhs 	: the right hand side
!    bc  	: the bc array
!=============================================================================!
	implicit none

	integer, intent(in) :: n, nblk, nsys
	complex, intent(inout) :: mat(nsys,n,nblk,nblk*3), rhs(nsys,n,nblk)
	complex, intent(inout) :: bc(nsys,nblk,nblk,14)
	complex, intent(inout) :: mult(nsys), fact(nsys)
	
!.... useful constants

	real, parameter :: zero = 0.0, one = 1.0

!.... local variables
	
	integer :: nc, nr, i, j, l, m, p, q, s, r, t, iv
!=============================================================================!
!                     F O R W A R D   E L I M I N A T I O N 
!=============================================================================!

!.... first 3 rows are special due to boundary conditions

	i = 1
	
	do nc = 1, nblk - 1

	  j = nblk + nc

	  do m = nc, nblk - 1
	    l = m + 1
	    mult(:) = -mat(:,i,l,j) / mat(:,i,nc,j)
	    do q = nc + nblk + 1, nblk * 2
	      mat(:,i,l,q) = mat(:,i,l,q) + mult(:) * mat(:,i,nc,q)
	    end do
	    
	    do iv = 1, nsys
	      mat(iv,i,l,11) = mat(iv,i,l,11) + mult(iv) * mat(iv,i,nc,11)
	      mat(iv,i,l,12) = mat(iv,i,l,12) + mult(iv) * mat(iv,i,nc,12)
	      mat(iv,i,l,13) = mat(iv,i,l,13) + mult(iv) * mat(iv,i,nc,13)
	      mat(iv,i,l,14) = mat(iv,i,l,14) + mult(iv) * mat(iv,i,nc,14)
	      mat(iv,i,l,15) = mat(iv,i,l,15) + mult(iv) * mat(iv,i,nc,15)
  
	      bc(iv,l,1,1) = bc(iv,l,1,1) + mult(iv) * bc(iv,nc,1,1)
	      bc(iv,l,2,1) = bc(iv,l,2,1) + mult(iv) * bc(iv,nc,2,1)
	      bc(iv,l,3,1) = bc(iv,l,3,1) + mult(iv) * bc(iv,nc,3,1)
	      bc(iv,l,4,1) = bc(iv,l,4,1) + mult(iv) * bc(iv,nc,4,1)
	      bc(iv,l,5,1) = bc(iv,l,5,1) + mult(iv) * bc(iv,nc,5,1)
  
	      bc(iv,l,1,2) = bc(iv,l,1,2) + mult(iv) * bc(iv,nc,1,2)
	      bc(iv,l,2,2) = bc(iv,l,2,2) + mult(iv) * bc(iv,nc,2,2)
	      bc(iv,l,3,2) = bc(iv,l,3,2) + mult(iv) * bc(iv,nc,3,2)
	      bc(iv,l,4,2) = bc(iv,l,4,2) + mult(iv) * bc(iv,nc,4,2)
	      bc(iv,l,5,2) = bc(iv,l,5,2) + mult(iv) * bc(iv,nc,5,2)
  
	      bc(iv,l,1,3) = bc(iv,l,1,3) + mult(iv) * bc(iv,nc,1,3)
	      bc(iv,l,2,3) = bc(iv,l,2,3) + mult(iv) * bc(iv,nc,2,3)
	      bc(iv,l,3,3) = bc(iv,l,3,3) + mult(iv) * bc(iv,nc,3,3)
	      bc(iv,l,4,3) = bc(iv,l,4,3) + mult(iv) * bc(iv,nc,4,3)
	      bc(iv,l,5,3) = bc(iv,l,5,3) + mult(iv) * bc(iv,nc,5,3)
	      
	      rhs(iv,i,l) = rhs(iv,i,l) + mult(iv) * rhs(iv,i,nc)
	    end do

	  end do	      

	end do	! end of loop on nc

!.... for all columns
	  
	p = i + 1
	r = i + 2
	  
	do nc = 1, nblk

	  j = nc + nblk
	  fact(:) = one / mat(:,i,nc,j)
	  
	  do m = 1, nblk

!.... second row

	    mult(:) = -mat(:,p,m,nc) * fact(:)

	    do q = nc + 1, nblk
	      l = q + nblk
	      mat(:,p,m,q) = mat(:,p,m,q) + mult(:) * mat(:,i,nc,l)
	    end do

	    do iv = 1, nsys
	      mat(iv,p,m, 6) = mat(iv,p,m, 6) + mult(iv) * mat(iv,i,nc,11)
	      mat(iv,p,m, 7) = mat(iv,p,m, 7) + mult(iv) * mat(iv,i,nc,12)
	      mat(iv,p,m, 8) = mat(iv,p,m, 8) + mult(iv) * mat(iv,i,nc,13)
	      mat(iv,p,m, 9) = mat(iv,p,m, 9) + mult(iv) * mat(iv,i,nc,14)
	      mat(iv,p,m,10) = mat(iv,p,m,10) + mult(iv) * mat(iv,i,nc,15)
  
	      mat(iv,p,m,11) = mat(iv,p,m,11) + mult(iv) * bc(iv,nc,1,1)
	      mat(iv,p,m,12) = mat(iv,p,m,12) + mult(iv) * bc(iv,nc,2,1)
	      mat(iv,p,m,13) = mat(iv,p,m,13) + mult(iv) * bc(iv,nc,3,1)
	      mat(iv,p,m,14) = mat(iv,p,m,14) + mult(iv) * bc(iv,nc,4,1)
	      mat(iv,p,m,15) = mat(iv,p,m,15) + mult(iv) * bc(iv,nc,5,1)
  
	      bc(iv,m,1,4)  = bc(iv,m,1,4)  + mult(iv) * bc(iv,nc,1,2)
	      bc(iv,m,2,4)  = bc(iv,m,2,4)  + mult(iv) * bc(iv,nc,2,2)
	      bc(iv,m,3,4)  = bc(iv,m,3,4)  + mult(iv) * bc(iv,nc,3,2)
	      bc(iv,m,4,4)  = bc(iv,m,4,4)  + mult(iv) * bc(iv,nc,4,2)
	      bc(iv,m,5,4)  = bc(iv,m,5,4)  + mult(iv) * bc(iv,nc,5,2)
  
	      bc(iv,m,1,5)  = bc(iv,m,1,5)  + mult(iv) * bc(iv,nc,1,3)
	      bc(iv,m,2,5)  = bc(iv,m,2,5)  + mult(iv) * bc(iv,nc,2,3)
	      bc(iv,m,3,5)  = bc(iv,m,3,5)  + mult(iv) * bc(iv,nc,3,3)
	      bc(iv,m,4,5)  = bc(iv,m,4,5)  + mult(iv) * bc(iv,nc,4,3)
	      bc(iv,m,5,5)  = bc(iv,m,5,5)  + mult(iv) * bc(iv,nc,5,3)
  
	      rhs(iv,p,m) = rhs(iv,p,m) + mult(iv) * rhs(iv,i,nc)
	    end do

!.... third row

	    mult(:) = -bc(:,m,nc,7) * fact(:)

	    do q = nc + 1, nblk
	      l = q + nblk
	      bc(:,m,q,7) = bc(:,m,q,7) + mult(:) * mat(:,i,nc,l)
	    end do

	    do iv = 1, nsys
	      mat(iv,r,m, 1) = mat(iv,r,m, 1) + mult(iv) * mat(iv,i,nc,11)
	      mat(iv,r,m, 2) = mat(iv,r,m, 2) + mult(iv) * mat(iv,i,nc,12)
	      mat(iv,r,m, 3) = mat(iv,r,m, 3) + mult(iv) * mat(iv,i,nc,13)
	      mat(iv,r,m, 4) = mat(iv,r,m, 4) + mult(iv) * mat(iv,i,nc,14)
	      mat(iv,r,m, 5) = mat(iv,r,m, 5) + mult(iv) * mat(iv,i,nc,15)
  
	      mat(iv,r,m, 6) = mat(iv,r,m, 6) + mult(iv) * bc(iv,nc,1,1)
	      mat(iv,r,m, 7) = mat(iv,r,m, 7) + mult(iv) * bc(iv,nc,2,1)
	      mat(iv,r,m, 8) = mat(iv,r,m, 8) + mult(iv) * bc(iv,nc,3,1)
	      mat(iv,r,m, 9) = mat(iv,r,m, 9) + mult(iv) * bc(iv,nc,4,1)
	      mat(iv,r,m,10) = mat(iv,r,m,10) + mult(iv) * bc(iv,nc,5,1)
  
	      mat(iv,r,m,11) = mat(iv,r,m,11) + mult(iv) * bc(iv,nc,1,2)
	      mat(iv,r,m,12) = mat(iv,r,m,12) + mult(iv) * bc(iv,nc,2,2)
	      mat(iv,r,m,13) = mat(iv,r,m,13) + mult(iv) * bc(iv,nc,3,2)
	      mat(iv,r,m,14) = mat(iv,r,m,14) + mult(iv) * bc(iv,nc,4,2)
	      mat(iv,r,m,15) = mat(iv,r,m,15) + mult(iv) * bc(iv,nc,5,2)
  
	      bc(iv,m,1,6)  = bc(iv,m,1,6)  + mult(iv) * bc(iv,nc,1,3)
	      bc(iv,m,2,6)  = bc(iv,m,2,6)  + mult(iv) * bc(iv,nc,2,3)
	      bc(iv,m,3,6)  = bc(iv,m,3,6)  + mult(iv) * bc(iv,nc,3,3)
	      bc(iv,m,4,6)  = bc(iv,m,4,6)  + mult(iv) * bc(iv,nc,4,3)
	      bc(iv,m,5,6)  = bc(iv,m,5,6)  + mult(iv) * bc(iv,nc,5,3)
  
	      rhs(iv,r,m) = rhs(iv,r,m) + mult(iv) * rhs(iv,i,nc)
	    end do
	  end do

	end do	! end of loop on nc

!.... second special row

	i = 2
	
	do nc = 1, nblk - 1

	  j = nblk + nc

	  do m = nc, nblk - 1
	    l = m + 1
	    mult(:) = -mat(:,i,l,j) / mat(:,i,nc,j)
	    do q = nc + nblk + 1, nblk * 2
	      mat(:,i,l,q) = mat(:,i,l,q) + mult(:) * mat(:,i,nc,q)
	    end do
	    
	    do iv = 1, nsys
	      mat(iv,i,l,11) = mat(iv,i,l,11) + mult(iv) * mat(iv,i,nc,11)
	      mat(iv,i,l,12) = mat(iv,i,l,12) + mult(iv) * mat(iv,i,nc,12)
	      mat(iv,i,l,13) = mat(iv,i,l,13) + mult(iv) * mat(iv,i,nc,13)
	      mat(iv,i,l,14) = mat(iv,i,l,14) + mult(iv) * mat(iv,i,nc,14)
	      mat(iv,i,l,15) = mat(iv,i,l,15) + mult(iv) * mat(iv,i,nc,15)
  
	      bc(iv,l,1,4) = bc(iv,l,1,4) + mult(iv) * bc(iv,nc,1,4)
	      bc(iv,l,2,4) = bc(iv,l,2,4) + mult(iv) * bc(iv,nc,2,4)
	      bc(iv,l,3,4) = bc(iv,l,3,4) + mult(iv) * bc(iv,nc,3,4)
	      bc(iv,l,4,4) = bc(iv,l,4,4) + mult(iv) * bc(iv,nc,4,4)
	      bc(iv,l,5,4) = bc(iv,l,5,4) + mult(iv) * bc(iv,nc,5,4)
	      
	      bc(iv,l,1,5) = bc(iv,l,1,5) + mult(iv) * bc(iv,nc,1,5)
	      bc(iv,l,2,5) = bc(iv,l,2,5) + mult(iv) * bc(iv,nc,2,5)
	      bc(iv,l,3,5) = bc(iv,l,3,5) + mult(iv) * bc(iv,nc,3,5)
	      bc(iv,l,4,5) = bc(iv,l,4,5) + mult(iv) * bc(iv,nc,4,5)
	      bc(iv,l,5,5) = bc(iv,l,5,5) + mult(iv) * bc(iv,nc,5,5)
  
	      rhs(iv,i,l) = rhs(iv,i,l) + mult(iv) * rhs(iv,i,nc)
	    end do
	  end do	      

	end do	! end of loop on nc

!.... for all columns
	  
	p  = i + 1
	  
	do nc = 1, nblk

	  j = nc + nblk
	  fact(:) = one / mat(:,i,nc,j)
	  
	  do m = 1, nblk
	    mult(:) = -mat(:,p,m,nc) * fact(:)

	    do q = nc + 1, nblk
	      l = q + nblk
	      mat(:,p,m,q) = mat(:,p,m,q) + mult(:) * mat(:,i,nc,l)
	    end do

	    do iv = 1, nsys
	      mat(iv,p,m, 6) = mat(iv,p,m, 6) + mult(iv) * mat(iv,i,nc,11)
	      mat(iv,p,m, 7) = mat(iv,p,m, 7) + mult(iv) * mat(iv,i,nc,12)
	      mat(iv,p,m, 8) = mat(iv,p,m, 8) + mult(iv) * mat(iv,i,nc,13)
	      mat(iv,p,m, 9) = mat(iv,p,m, 9) + mult(iv) * mat(iv,i,nc,14)
	      mat(iv,p,m,10) = mat(iv,p,m,10) + mult(iv) * mat(iv,i,nc,15)
  
	      mat(iv,p,m,11) = mat(iv,p,m,11) + mult(iv) * bc(iv,nc,1,4)
	      mat(iv,p,m,12) = mat(iv,p,m,12) + mult(iv) * bc(iv,nc,2,4)
	      mat(iv,p,m,13) = mat(iv,p,m,13) + mult(iv) * bc(iv,nc,3,4)
	      mat(iv,p,m,14) = mat(iv,p,m,14) + mult(iv) * bc(iv,nc,4,4)
	      mat(iv,p,m,15) = mat(iv,p,m,15) + mult(iv) * bc(iv,nc,5,4)
  
	      bc(iv,m,1,6)  = bc(iv,m,1,6)  + mult(iv) * bc(iv,nc,1,5)
	      bc(iv,m,2,6)  = bc(iv,m,2,6)  + mult(iv) * bc(iv,nc,2,5)
	      bc(iv,m,3,6)  = bc(iv,m,3,6)  + mult(iv) * bc(iv,nc,3,5)
	      bc(iv,m,4,6)  = bc(iv,m,4,6)  + mult(iv) * bc(iv,nc,4,5)
	      bc(iv,m,5,6)  = bc(iv,m,5,6)  + mult(iv) * bc(iv,nc,5,5)
  
	      rhs(iv,p,m) = rhs(iv,p,m) + mult(iv) * rhs(iv,i,nc)
	    end do
	  end do

	end do	! end of loop on nc

!.... the third and final special row

	i = 3
	
	do nc = 1, nblk - 1

	  j = nblk + nc

	  do m = nc, nblk - 1
	    l = m + 1
	    mult(:) = -mat(:,i,l,j) / mat(:,i,nc,j)
	    do q = nc + nblk + 1, nblk * 2
	      mat(:,i,l,q) = mat(:,i,l,q) + mult(:) * mat(:,i,nc,q)
	    end do
	    
	    do iv = 1, nsys
	      mat(iv,i,l,11) = mat(iv,i,l,11) + mult(iv) * mat(iv,i,nc,11)
	      mat(iv,i,l,12) = mat(iv,i,l,12) + mult(iv) * mat(iv,i,nc,12)
	      mat(iv,i,l,13) = mat(iv,i,l,13) + mult(iv) * mat(iv,i,nc,13)
	      mat(iv,i,l,14) = mat(iv,i,l,14) + mult(iv) * mat(iv,i,nc,14)
	      mat(iv,i,l,15) = mat(iv,i,l,15) + mult(iv) * mat(iv,i,nc,15)
  
	      bc(iv,l,1,6) = bc(iv,l,1,6) + mult(iv) * bc(iv,nc,1,6)
	      bc(iv,l,2,6) = bc(iv,l,2,6) + mult(iv) * bc(iv,nc,2,6)
	      bc(iv,l,3,6) = bc(iv,l,3,6) + mult(iv) * bc(iv,nc,3,6)
	      bc(iv,l,4,6) = bc(iv,l,4,6) + mult(iv) * bc(iv,nc,4,6)
	      bc(iv,l,5,6) = bc(iv,l,5,6) + mult(iv) * bc(iv,nc,5,6)
	      
	      rhs(iv,i,l) = rhs(iv,i,l) + mult(iv) * rhs(iv,i,nc)
	    end do
	  end do	      

	end do	! end of loop on nc

!.... for all columns
	  
	p  = i + 1
	  
	do nc = 1, nblk

	  j = nc + nblk
	  fact(:) = one / mat(:,i,nc,j)
	  
	  do m = 1, nblk
	    mult(:) = -mat(:,p,m,nc) * fact(:)

	    do q = nc + 1, nblk
	      l = q + nblk
	      mat(:,p,m,q) = mat(:,p,m,q) + mult(:) * mat(:,i,nc,l)
	    end do

	    do iv = 1, nsys
	      mat(iv,p,m, 6) = mat(iv,p,m, 6) + mult(iv) * mat(iv,i,nc,11)
	      mat(iv,p,m, 7) = mat(iv,p,m, 7) + mult(iv) * mat(iv,i,nc,12)
	      mat(iv,p,m, 8) = mat(iv,p,m, 8) + mult(iv) * mat(iv,i,nc,13)
	      mat(iv,p,m, 9) = mat(iv,p,m, 9) + mult(iv) * mat(iv,i,nc,14)
	      mat(iv,p,m,10) = mat(iv,p,m,10) + mult(iv) * mat(iv,i,nc,15)
  
	      mat(iv,p,m,11) = mat(iv,p,m,11) + mult(iv) * bc(iv,nc,1,6)
	      mat(iv,p,m,12) = mat(iv,p,m,12) + mult(iv) * bc(iv,nc,2,6)
	      mat(iv,p,m,13) = mat(iv,p,m,13) + mult(iv) * bc(iv,nc,3,6)
	      mat(iv,p,m,14) = mat(iv,p,m,14) + mult(iv) * bc(iv,nc,4,6)
	      mat(iv,p,m,15) = mat(iv,p,m,15) + mult(iv) * bc(iv,nc,5,6)
  
	      rhs(iv,p,m) = rhs(iv,p,m) + mult(iv) * rhs(iv,i,nc)
	    end do
	  end do

	end do	! end of loop on nc

!.... the rest of the matrix is a standard tridiagonal up to the other boundary

 	do i = 4, n - 5

!.... do the first four columns

	  do nc = 1, nblk - 1
	  
	    j = nblk + nc

	    do m = nc, nblk - 1
	      l = m + 1
	      mult(:) = -mat(:,i,l,j) / mat(:,i,nc,j)
	      do q = nc + nblk + 1, nblk * 2
		mat(:,i,l,q) = mat(:,i,l,q) + mult(:) * mat(:,i,nc,q)
	      end do

	      do iv = 1, nsys
		mat(iv,i,l,11) = mat(iv,i,l,11) + mult(iv) * mat(iv,i,nc,11)
		mat(iv,i,l,12) = mat(iv,i,l,12) + mult(iv) * mat(iv,i,nc,12)
		mat(iv,i,l,13) = mat(iv,i,l,13) + mult(iv) * mat(iv,i,nc,13)
		mat(iv,i,l,14) = mat(iv,i,l,14) + mult(iv) * mat(iv,i,nc,14)
		mat(iv,i,l,15) = mat(iv,i,l,15) + mult(iv) * mat(iv,i,nc,15)
  
		rhs(iv,i,l) = rhs(iv,i,l) + mult(iv) * rhs(iv,i,nc)
	      end do
	    end do	      

	  end do	! end of loop on nc

!.... for all columns
	  
	  p  = i + 1
	    
	  do nc = 1, nblk

	    j = nc + nblk
	    fact(:) = one / mat(:,i,nc,j)
	    
	    do m = 1, nblk
	      mult(:) = -mat(:,p,m,nc) * fact(:)
	      do q = nc + 1, nblk
		l = q + nblk
		mat(:,p,m,q) = mat(:,p,m,q) + mult(:) * mat(:,i,nc,l)
	      end do

	      do iv = 1, nsys
		mat(iv,p,m, 6) = mat(iv,p,m, 6) + mult(iv) * mat(iv,i,nc,11)
		mat(iv,p,m, 7) = mat(iv,p,m, 7) + mult(iv) * mat(iv,i,nc,12)
		mat(iv,p,m, 8) = mat(iv,p,m, 8) + mult(iv) * mat(iv,i,nc,13)
		mat(iv,p,m, 9) = mat(iv,p,m, 9) + mult(iv) * mat(iv,i,nc,14)
		mat(iv,p,m,10) = mat(iv,p,m,10) + mult(iv) * mat(iv,i,nc,15)
  
		rhs(iv,p,m) = rhs(iv,p,m) + mult(iv) * rhs(iv,i,nc)
	      end do
	    end do

	  end do	! end of loop on nc
	  
	end do		! end of loop on i

!.... do the FIFTH to the last row

	i = n - 4
	
!.... first four columns

	  do nc = 1, nblk - 1
	  
	    j = nblk + nc

	    do m = nc, nblk - 1
	      l = m + 1
	      mult(:) = -mat(:,i,l,j) / mat(:,i,nc,j)
	      do q = nc + nblk + 1, nblk * 2
		mat(:,i,l,q) = mat(:,i,l,q) + mult(:) * mat(:,i,nc,q)
	      end do

	      do iv = 1, nsys
		mat(iv,i,l,11) = mat(iv,i,l,11) + mult(iv) * mat(iv,i,nc,11)
		mat(iv,i,l,12) = mat(iv,i,l,12) + mult(iv) * mat(iv,i,nc,12)
		mat(iv,i,l,13) = mat(iv,i,l,13) + mult(iv) * mat(iv,i,nc,13)
		mat(iv,i,l,14) = mat(iv,i,l,14) + mult(iv) * mat(iv,i,nc,14)
		mat(iv,i,l,15) = mat(iv,i,l,15) + mult(iv) * mat(iv,i,nc,15)
  
		rhs(iv,i,l) = rhs(iv,i,l) + mult(iv) * rhs(iv,i,nc)
	      end do
	    end do	      

	  end do	! end of loop on nc

!.... for all columns
	  
	  p = i + 1
	  r = i + 2
	  s = i + 3
	  t = i + 4
	  	    
	  do nc = 1, nblk

	    j = nc + nblk
	    fact(:) = one / mat(:,i,nc,j)
	    
	    do m = 1, nblk
	      mult(:) = -mat(:,p,m,nc) * fact(:)
	      do q = nc + 1, nblk
		l = q + nblk
		mat(:,p,m,q) = mat(:,p,m,q) + mult(:) * mat(:,i,nc,l)
	      end do

	      do iv = 1, nsys
		mat(iv,p,m, 6) = mat(iv,p,m, 6) + mult(iv) * mat(iv,i,nc,11)
		mat(iv,p,m, 7) = mat(iv,p,m, 7) + mult(iv) * mat(iv,i,nc,12)
		mat(iv,p,m, 8) = mat(iv,p,m, 8) + mult(iv) * mat(iv,i,nc,13)
		mat(iv,p,m, 9) = mat(iv,p,m, 9) + mult(iv) * mat(iv,i,nc,14)
		mat(iv,p,m,10) = mat(iv,p,m,10) + mult(iv) * mat(iv,i,nc,15)
  
		rhs(iv,p,m) = rhs(iv,p,m) + mult(iv) * rhs(iv,i,nc)
	      end do

!.... eliminate the boundary conditions: bc(8)

	      mult(:) = -bc(:,m,nc,8) * fact(:)
	      do q = nc + 1, nblk
		l = q + nblk
		bc(:,m,q,8) = bc(:,m,q,8) + mult(:) * mat(:,i,nc,l)
	      end do
  
	      do iv = 1, nsys
		mat(iv,r,m, 1) = mat(iv,r,m, 1) + mult(iv) * mat(iv,i,nc,11)
		mat(iv,r,m, 2) = mat(iv,r,m, 2) + mult(iv) * mat(iv,i,nc,12)
		mat(iv,r,m, 3) = mat(iv,r,m, 3) + mult(iv) * mat(iv,i,nc,13)
		mat(iv,r,m, 4) = mat(iv,r,m, 4) + mult(iv) * mat(iv,i,nc,14)
		mat(iv,r,m, 5) = mat(iv,r,m, 5) + mult(iv) * mat(iv,i,nc,15)
      
		rhs(iv,r,m) = rhs(iv,r,m) + mult(iv) * rhs(iv,i,nc)
	      end do

!.... eliminate the boundary conditions: bc(10)

	      mult(:) = -bc(:,m,nc,10) * fact(:)
	      do q = nc + 1, nblk
		l = q + nblk
		bc(:,m,q,10) = bc(:,m,q,10) + mult(:) * mat(:,i,nc,l)
	      end do
 
 	      do iv = 1, nsys
		bc(iv,m, 1,11) = bc(iv,m, 1,11) + mult(iv) * mat(iv,i,nc,11)
		bc(iv,m, 2,11) = bc(iv,m, 2,11) + mult(iv) * mat(iv,i,nc,12)
		bc(iv,m, 3,11) = bc(iv,m, 3,11) + mult(iv) * mat(iv,i,nc,13)
		bc(iv,m, 4,11) = bc(iv,m, 4,11) + mult(iv) * mat(iv,i,nc,14)
		bc(iv,m, 5,11) = bc(iv,m, 5,11) + mult(iv) * mat(iv,i,nc,15)
      
		rhs(iv,s,m) = rhs(iv,s,m) + mult(iv) * rhs(iv,i,nc)
	      end do
		
!.... eliminate the boundary conditions: bc(12)

	      mult(:) = -bc(:,m,nc,12) * fact(:)
	      do q = nc + 1, nblk
		l = q + nblk
		bc(:,m,q,12) = bc(:,m,q,12) + mult(:) * mat(:,i,nc,l)
	      end do
  
	      do iv = 1, nsys
		bc(iv,m, 1,13) = bc(iv,m, 1,13) + mult(iv) * mat(iv,i,nc,11)
		bc(iv,m, 2,13) = bc(iv,m, 2,13) + mult(iv) * mat(iv,i,nc,12)
		bc(iv,m, 3,13) = bc(iv,m, 3,13) + mult(iv) * mat(iv,i,nc,13)
		bc(iv,m, 4,13) = bc(iv,m, 4,13) + mult(iv) * mat(iv,i,nc,14)
		bc(iv,m, 5,13) = bc(iv,m, 5,13) + mult(iv) * mat(iv,i,nc,15)
      
		rhs(iv,t,m) = rhs(iv,t,m) + mult(iv) * rhs(iv,i,nc)
	      end do

	    end do

	  end do	! end of loop on nc

!.... do the FOURTH to the last row

	i = n - 3
	
!.... first four columns

	  do nc = 1, nblk - 1
	  
	    j = nblk + nc

	    do m = nc, nblk - 1
	      l = m + 1
	      mult(:) = -mat(:,i,l,j) / mat(:,i,nc,j)
	      do q = nc + nblk + 1, nblk * 2
		mat(:,i,l,q) = mat(:,i,l,q) + mult(:) * mat(:,i,nc,q)
	      end do

	      do iv = 1, nsys
		mat(iv,i,l,11) = mat(iv,i,l,11) + mult(iv) * mat(iv,i,nc,11)
		mat(iv,i,l,12) = mat(iv,i,l,12) + mult(iv) * mat(iv,i,nc,12)
		mat(iv,i,l,13) = mat(iv,i,l,13) + mult(iv) * mat(iv,i,nc,13)
		mat(iv,i,l,14) = mat(iv,i,l,14) + mult(iv) * mat(iv,i,nc,14)
		mat(iv,i,l,15) = mat(iv,i,l,15) + mult(iv) * mat(iv,i,nc,15)
  
		rhs(iv,i,l) = rhs(iv,i,l) + mult(iv) * rhs(iv,i,nc)
	      end do
	    end do	      

	  end do	! end of loop on nc

!.... for all columns
	  
	  p = i + 1
	  r = i + 2
	  s = i + 3
	  	    
	  do nc = 1, nblk

	    j = nc + nblk
	    fact(:) = one / mat(:,i,nc,j)
	    
	    do m = 1, nblk
	      mult(:) = -mat(:,p,m,nc) * fact(:)
	      do q = nc + 1, nblk
		l = q + nblk
		mat(:,p,m,q) = mat(:,p,m,q) + mult(:) * mat(:,i,nc,l)
	      end do

	      do iv = 1, nsys
		mat(iv,p,m, 6) = mat(iv,p,m, 6) + mult(iv) * mat(iv,i,nc,11)
		mat(iv,p,m, 7) = mat(iv,p,m, 7) + mult(iv) * mat(iv,i,nc,12)
		mat(iv,p,m, 8) = mat(iv,p,m, 8) + mult(iv) * mat(iv,i,nc,13)
		mat(iv,p,m, 9) = mat(iv,p,m, 9) + mult(iv) * mat(iv,i,nc,14)
		mat(iv,p,m,10) = mat(iv,p,m,10) + mult(iv) * mat(iv,i,nc,15)
  
		rhs(iv,p,m) = rhs(iv,p,m) + mult(iv) * rhs(iv,i,nc)
	      end do
	      
!.... eliminate the boundary conditions: bc(11)

	      mult(:) = -bc(:,m,nc,11) * fact(:)
	      do q = nc + 1, nblk
		l = q + nblk
		bc(:,m,q,11) = bc(:,m,q,11) + mult(:) * mat(:,i,nc,l)
	      end do
  
	      do iv = 1, nsys
		mat(iv,r,m, 1) = mat(iv,r,m, 1) + mult(iv) * mat(iv,i,nc,11)
		mat(iv,r,m, 2) = mat(iv,r,m, 2) + mult(iv) * mat(iv,i,nc,12)
		mat(iv,r,m, 3) = mat(iv,r,m, 3) + mult(iv) * mat(iv,i,nc,13)
		mat(iv,r,m, 4) = mat(iv,r,m, 4) + mult(iv) * mat(iv,i,nc,14)
		mat(iv,r,m, 5) = mat(iv,r,m, 5) + mult(iv) * mat(iv,i,nc,15)
      
		rhs(iv,r,m) = rhs(iv,r,m) + mult(iv) * rhs(iv,i,nc)
	      end do
	      
!.... eliminate the boundary conditions: bc(13)

	      mult(:) = -bc(:,m,nc,13) * fact(:)
	      do q = nc + 1, nblk
		l = q + nblk
		bc(:,m,q,13) = bc(:,m,q,13) + mult(:) * mat(:,i,nc,l)
	      end do
  
	      do iv = 1, nsys
		bc(iv,m, 1,14) = bc(iv,m, 1,14) + mult(iv) * mat(iv,i,nc,11)
		bc(iv,m, 2,14) = bc(iv,m, 2,14) + mult(iv) * mat(iv,i,nc,12)
		bc(iv,m, 3,14) = bc(iv,m, 3,14) + mult(iv) * mat(iv,i,nc,13)
		bc(iv,m, 4,14) = bc(iv,m, 4,14) + mult(iv) * mat(iv,i,nc,14)
		bc(iv,m, 5,14) = bc(iv,m, 5,14) + mult(iv) * mat(iv,i,nc,15)
      
		rhs(iv,s,m) = rhs(iv,s,m) + mult(iv) * rhs(iv,i,nc)
	      end do
	    end do

	  end do	! end of loop on nc

!.... do the THIRD to the last row

	i = n - 2
	
!.... first four columns

	  do nc = 1, nblk - 1
	  
	    j = nblk + nc

	    do m = nc, nblk - 1
	      l = m + 1
	      mult(:) = -mat(:,i,l,j) / mat(:,i,nc,j)
	      do q = nc + nblk + 1, nblk * 2
		mat(:,i,l,q) = mat(:,i,l,q) + mult(:) * mat(:,i,nc,q)
	      end do

	      do iv = 1, nsys
		mat(iv,i,l,11) = mat(iv,i,l,11) + mult(iv) * mat(iv,i,nc,11)
		mat(iv,i,l,12) = mat(iv,i,l,12) + mult(iv) * mat(iv,i,nc,12)
		mat(iv,i,l,13) = mat(iv,i,l,13) + mult(iv) * mat(iv,i,nc,13)
		mat(iv,i,l,14) = mat(iv,i,l,14) + mult(iv) * mat(iv,i,nc,14)
		mat(iv,i,l,15) = mat(iv,i,l,15) + mult(iv) * mat(iv,i,nc,15)
  
		bc(iv,l, 1,9) = bc(iv,l, 1,9) + mult(iv) * bc(iv,nc,1,9)
		bc(iv,l, 2,9) = bc(iv,l, 2,9) + mult(iv) * bc(iv,nc,2,9)
		bc(iv,l, 3,9) = bc(iv,l, 3,9) + mult(iv) * bc(iv,nc,3,9)
		bc(iv,l, 4,9) = bc(iv,l, 4,9) + mult(iv) * bc(iv,nc,4,9)
		bc(iv,l, 5,9) = bc(iv,l, 5,9) + mult(iv) * bc(iv,nc,5,9)
  
		rhs(iv,i,l) = rhs(iv,i,l) + mult(iv) * rhs(iv,i,nc)
	      end do
	    end do	      

	  end do	! end of loop on nc

!.... for all columns
	  
	  p = i + 1
	  r = i + 2
	    
	  do nc = 1, nblk

	    j = nc + nblk
	    fact(:) = one / mat(:,i,nc,j)
	    
	    do m = 1, nblk
	      mult(:) = -mat(:,p,m,nc) * fact(:)
	      do q = nc + 1, nblk
		l = q + nblk
		mat(:,p,m,q) = mat(:,p,m,q) + mult(:) * mat(:,i,nc,l)
	      end do

	      do iv = 1, nsys
		mat(iv,p,m, 6) = mat(iv,p,m, 6) + mult(iv) * mat(iv,i,nc,11)
		mat(iv,p,m, 7) = mat(iv,p,m, 7) + mult(iv) * mat(iv,i,nc,12)
		mat(iv,p,m, 8) = mat(iv,p,m, 8) + mult(iv) * mat(iv,i,nc,13)
		mat(iv,p,m, 9) = mat(iv,p,m, 9) + mult(iv) * mat(iv,i,nc,14)
		mat(iv,p,m,10) = mat(iv,p,m,10) + mult(iv) * mat(iv,i,nc,15)
  
		mat(iv,p,m,11) = mat(iv,p,m,11) + mult(iv) * bc(iv,nc,1,9)
		mat(iv,p,m,12) = mat(iv,p,m,12) + mult(iv) * bc(iv,nc,2,9)
		mat(iv,p,m,13) = mat(iv,p,m,13) + mult(iv) * bc(iv,nc,3,9)
		mat(iv,p,m,14) = mat(iv,p,m,14) + mult(iv) * bc(iv,nc,4,9)
		mat(iv,p,m,15) = mat(iv,p,m,15) + mult(iv) * bc(iv,nc,5,9)
  
		rhs(iv,p,m) = rhs(iv,p,m) + mult(iv) * rhs(iv,i,nc)
	      end do
	      
!.... eliminate the boundary condition:  bc(14)

	      mult(:) = -bc(:,m,nc,14) * fact(:)
	      do q = nc + 1, nblk
		l = q + nblk
		bc(:,m,q,14) = bc(:,m,q,14) + mult(:) * mat(:,i,nc,l)
	      end do
  
	      do iv = 1, nsys
		mat(iv,r,m, 1) = mat(iv,r,m, 1) + mult(iv) * mat(iv,i,nc,11)
		mat(iv,r,m, 2) = mat(iv,r,m, 2) + mult(iv) * mat(iv,i,nc,12)
		mat(iv,r,m, 3) = mat(iv,r,m, 3) + mult(iv) * mat(iv,i,nc,13)
		mat(iv,r,m, 4) = mat(iv,r,m, 4) + mult(iv) * mat(iv,i,nc,14)
		mat(iv,r,m, 5) = mat(iv,r,m, 5) + mult(iv) * mat(iv,i,nc,15)
      
		mat(iv,r,m, 6) = mat(iv,r,m, 6) + mult(iv) * bc(iv,nc,1,9)
		mat(iv,r,m, 7) = mat(iv,r,m, 7) + mult(iv) * bc(iv,nc,2,9)
		mat(iv,r,m, 8) = mat(iv,r,m, 8) + mult(iv) * bc(iv,nc,3,9)
		mat(iv,r,m, 9) = mat(iv,r,m, 9) + mult(iv) * bc(iv,nc,4,9)
		mat(iv,r,m,10) = mat(iv,r,m,10) + mult(iv) * bc(iv,nc,5,9)
  
		rhs(iv,r,m) = rhs(iv,r,m) + mult(iv) * rhs(iv,i,nc)
	      end do
	    end do

	  end do	! end of loop on nc

!.... do the SECOND to the last row (which is normal tridiagonal)

	i = n - 1
	
!.... first four columns

	  do nc = 1, nblk - 1
	  
	    j = nblk + nc

	    do m = nc, nblk - 1
	      l = m + 1
	      mult(:) = -mat(:,i,l,j) / mat(:,i,nc,j)
	      do q = nc + nblk + 1, nblk * 2
		mat(:,i,l,q) = mat(:,i,l,q) + mult(:) * mat(:,i,nc,q)
	      end do

	      do iv = 1, nsys
		mat(iv,i,l,11) = mat(iv,i,l,11) + mult(iv) * mat(iv,i,nc,11)
		mat(iv,i,l,12) = mat(iv,i,l,12) + mult(iv) * mat(iv,i,nc,12)
		mat(iv,i,l,13) = mat(iv,i,l,13) + mult(iv) * mat(iv,i,nc,13)
		mat(iv,i,l,14) = mat(iv,i,l,14) + mult(iv) * mat(iv,i,nc,14)
		mat(iv,i,l,15) = mat(iv,i,l,15) + mult(iv) * mat(iv,i,nc,15)
  
		rhs(iv,i,l) = rhs(iv,i,l) + mult(iv) * rhs(iv,i,nc)
	      end do
	    end do	      

	  end do	! end of loop on nc

!.... for all columns
	  
	  p = i + 1
	    
	  do nc = 1, nblk

	    j = nc + nblk
	    fact(:) = one / mat(:,i,nc,j)
	    
	    do m = 1, nblk
	      mult(:) = -mat(:,p,m,nc) * fact(:)
	      do q = nc + 1, nblk
		l = q + nblk
		mat(:,p,m,q) = mat(:,p,m,q) + mult(:) * mat(:,i,nc,l)
	      end do

	      do iv = 1, nsys
		mat(iv,p,m, 6) = mat(iv,p,m, 6) + mult(iv) * mat(iv,i,nc,11)
		mat(iv,p,m, 7) = mat(iv,p,m, 7) + mult(iv) * mat(iv,i,nc,12)
		mat(iv,p,m, 8) = mat(iv,p,m, 8) + mult(iv) * mat(iv,i,nc,13)
		mat(iv,p,m, 9) = mat(iv,p,m, 9) + mult(iv) * mat(iv,i,nc,14)
		mat(iv,p,m,10) = mat(iv,p,m,10) + mult(iv) * mat(iv,i,nc,15)
  
		rhs(iv,p,m) = rhs(iv,p,m) + mult(iv) * rhs(iv,i,nc)
	      end do
	    end do

	  end do	! end of loop on nc

!.... do the last row
	  
	i = n
	
	do nc = 1, nblk - 1
	
	  j = nblk + nc

	  do m = nc, nblk - 1
	    l = m + 1
	    mult(:) = -mat(:,i,l,j) / mat(:,i,nc,j)
	    do q = nc + nblk + 1, nblk * 2
	      mat(:,i,l,q) = mat(:,i,l,q) + mult(:) * mat(:,i,nc,q)
	    end do
	    rhs(:,i,l) = rhs(:,i,l) + mult(:) * rhs(:,i,nc)
	  end do	      

	end do	! end of loop on nc

!=============================================================================!
!                   B A C K W A R D   S U B S T I T U T I O N 
!=============================================================================!

!.... do the LAST row first

	i = n
	
	do m = nblk, 1, -1
	  j = m + nblk
	  do q = m+1, nblk
	    l = q + nblk
	    rhs(:,i,m) = rhs(:,i,m) - rhs(:,i,q) * mat(:,i,m,l)
	  end do
	  rhs(:,i,m) = rhs(:,i,m) / mat(:,i,m,j)
	end do
	
!.... do the SECOND to the last row (normal row)

	i = n - 1
	p = i + 1
	
	do m = nblk, 1, -1
	  j = m + nblk
	  do q = m+1, nblk
	    l = q + nblk
	    rhs(:,i,m) = rhs(:,i,m) - rhs(:,i,q) * mat(:,i,m,l)
	  end do

	  do iv = 1, nsys
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,p,1) * mat(iv,i,m,11)
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,p,2) * mat(iv,i,m,12)
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,p,3) * mat(iv,i,m,13)
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,p,4) * mat(iv,i,m,14)
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,p,5) * mat(iv,i,m,15)
  
	    rhs(iv,i,m) = rhs(iv,i,m) / mat(iv,i,m,j)
	  end do
	end do

!.... do the THIRD to the last row

	i = n - 2
	p = i + 1
	r = i + 2
	
	do m = nblk, 1, -1
	  j = m + nblk
	  do q = m+1, nblk
	    l = q + nblk
	    rhs(:,i,m) = rhs(:,i,m) - rhs(:,i,q) * mat(:,i,m,l)
	  end do

	  do iv = 1, nsys
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,p,1) * mat(iv,i,m,11)
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,p,2) * mat(iv,i,m,12)
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,p,3) * mat(iv,i,m,13)
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,p,4) * mat(iv,i,m,14)
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,p,5) * mat(iv,i,m,15)
  
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,r,1) * bc(iv,m,1,9)
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,r,2) * bc(iv,m,2,9)
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,r,3) * bc(iv,m,3,9)
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,r,4) * bc(iv,m,4,9)
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,r,5) * bc(iv,m,5,9)
  
	    rhs(iv,i,m) = rhs(iv,i,m) / mat(iv,i,m,j)
	  end do
	  
	end do

!.... now do the rest of the rows in reverse order upto the boundary

	do i = n - 3, 4, -1
	
	  p = i + 1
	  
	  do m = nblk, 1, -1
	    j = m + nblk
	    do q = m+1, nblk
	      l = q + nblk
	      rhs(:,i,m) = rhs(:,i,m) - rhs(:,i,q) * mat(:,i,m,l)
	    end do

	    do iv = 1, nsys
	      rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,p,1) * mat(iv,i,m,11)
	      rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,p,2) * mat(iv,i,m,12)
	      rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,p,3) * mat(iv,i,m,13)
	      rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,p,4) * mat(iv,i,m,14)
	      rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,p,5) * mat(iv,i,m,15)
  
	      rhs(iv,i,m) = rhs(iv,i,m) / mat(iv,i,m,j)
	    end do
	  end do
				  
	end do
	
!.... now do the boundary rows

	i = 3
	
	p = i + 1
	r = i + 2
	
	do m = nblk, 1, -1
	  j = m + nblk
	  do q = m+1, nblk
	    l = q + nblk
	    rhs(:,i,m) = rhs(:,i,m) - rhs(:,i,q) * mat(:,i,m,l)
	  end do

	  do iv = 1, nsys
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,p,1) * mat(iv,i,m,11)
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,p,2) * mat(iv,i,m,12)
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,p,3) * mat(iv,i,m,13)
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,p,4) * mat(iv,i,m,14)
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,p,5) * mat(iv,i,m,15)
  
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,r,1) * bc(iv,m,1,6)
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,r,2) * bc(iv,m,2,6)
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,r,3) * bc(iv,m,3,6)
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,r,4) * bc(iv,m,4,6)
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,r,5) * bc(iv,m,5,6)
  
	    rhs(iv,i,m) = rhs(iv,i,m) / mat(iv,i,m,j)
	  end do
	end do

	i = 2
	
	p = i + 1
	r = i + 2
	s = i + 3
	
	do m = nblk, 1, -1
	  j = m + nblk
	  do q = m+1, nblk
	    l = q + nblk
	    rhs(:,i,m) = rhs(:,i,m) - rhs(:,i,q) * mat(:,i,m,l)
	  end do

	  do iv = 1, nsys
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,p,1) * mat(iv,i,m,11)
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,p,2) * mat(iv,i,m,12)
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,p,3) * mat(iv,i,m,13)
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,p,4) * mat(iv,i,m,14)
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,p,5) * mat(iv,i,m,15)
  
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,r,1) * bc(iv,m,1,4)
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,r,2) * bc(iv,m,2,4)
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,r,3) * bc(iv,m,3,4)
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,r,4) * bc(iv,m,4,4)
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,r,5) * bc(iv,m,5,4)
  
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,s,1) * bc(iv,m,1,5)
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,s,2) * bc(iv,m,2,5)
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,s,3) * bc(iv,m,3,5)
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,s,4) * bc(iv,m,4,5)
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,s,5) * bc(iv,m,5,5)
  
	    rhs(iv,i,m) = rhs(iv,i,m) / mat(iv,i,m,j)
	  end do
	end do

	i = 1
	
	p = i + 1
	r = i + 2
	s = i + 3
	t = i + 4
	
	do m = nblk, 1, -1
	  j = m + nblk
	  do q = m+1, nblk
	    l = q + nblk
	    rhs(:,i,m) = rhs(:,i,m) - rhs(:,i,q) * mat(:,i,m,l)
	  end do

	  do iv = 1, nsys
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,p,1) * mat(iv,i,m,11)
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,p,2) * mat(iv,i,m,12)
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,p,3) * mat(iv,i,m,13)
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,p,4) * mat(iv,i,m,14)
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,p,5) * mat(iv,i,m,15)
  
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,r,1) * bc(iv,m,1,1)
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,r,2) * bc(iv,m,2,1)
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,r,3) * bc(iv,m,3,1)
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,r,4) * bc(iv,m,4,1)
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,r,5) * bc(iv,m,5,1)
  
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,s,1) * bc(iv,m,1,2)
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,s,2) * bc(iv,m,2,2)
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,s,3) * bc(iv,m,3,2)
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,s,4) * bc(iv,m,4,2)
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,s,5) * bc(iv,m,5,2)
  
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,t,1) * bc(iv,m,1,3)
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,t,2) * bc(iv,m,2,3)
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,t,3) * bc(iv,m,3,3)
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,t,4) * bc(iv,m,4,3)
	    rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,t,5) * bc(iv,m,5,3)
  
	    rhs(iv,i,m) = rhs(iv,i,m) / mat(iv,i,m,j)
	  end do
	end do
!=============================================================================!
	return
	end
