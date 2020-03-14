!=============================================================================!
	subroutine cpenta2p( np, nsys, nblk, mat, rhs, per, per2, 
     &                       mult, fact, code )
!=============================================================================!
!  
!  Solves a block pentadiagonal system of equations without pivoting.
!  Vectorized over the first index.
!
!  This routine solves a PERIODIC block pentadiagonal system.
!
!  WARNING: Since pivoting is not done, this algorithm could fail
!
!  input:
!    nsys  	: number of systems to vectorize over
!    n 	   	: number of rows in the systems
!    nblk  	: dimension of the blocks in each row
!    code	: 0 = LU factorization only
!               : 1 = Forward and back solves given LU
!		: 2 = LU and solve
!
!  inout:
!    mat	: the block pentadiagonal matrix to be solved
!    rhs 	: the right hand side
!=============================================================================!
	implicit none

	integer np, nblk, nsys, code
	complex mat(np,nsys,nblk,nblk,5), rhs(np,nsys,nblk)
	complex per(np,nsys,nblk,nblk,2), per2(np-2:np-1,nsys,nblk,nblk,np)
	complex mult(nsys), fact(nsys)
	
!.... useful constants

	real zero, one
	parameter (zero = 0.0, one = 1.0)

!.... local variables
	
	integer i, n, p, q, r, iv

	integer nc, l, m
!=============================================================================!

!.... don't solve for the redundant nodes

	n = np - 1

!=============================================================================!
!                       L U   F A C T O R I Z A T I O N 
!=============================================================================!
	if (code .eq. 0 .or. code .eq. 2) then

!	per  = zero
!	per2 = zero

      	call rzero( nsys*np*nblk*nblk*2, per )
      	call rzero( nsys*2*nblk*nblk*np, per2 )

!.... setup the periodic fillin vector

	do q = 1, nblk
	  do p = 1, nblk
	    do iv = 1, nsys
	      per(1,iv,p,q,1) = mat(1,iv,p,q,1)
	      per(1,iv,p,q,2) = mat(1,iv,p,q,2)
	      per(2,iv,p,q,1) = zero
	      per(2,iv,p,q,2) = mat(2,iv,p,q,1)
	      
	      per2(n-1,iv,p,q,1) = mat(n-1,iv,p,q,5)
	      per2(n-1,iv,p,q,2) = zero
	      per2(n,iv,p,q,1)   = mat(n,iv,p,q,4)
	      per2(n,iv,p,q,2)   = mat(n,iv,p,q,5)
	    end do
          end do
	end do

 	do i = 1, n - 6

!.... do the first four columns

	  do nc = 1, nblk - 1
	  
	    do m = nc, nblk - 1
	      l = m + 1
	      do iv = 1, nsys
		mult(iv) = -mat(i,iv,l,nc,3) / mat(i,iv,nc,nc,3)
		mat(i,iv,l,nc,3) = mult(iv)
	      end do
	      do q = nc + 1, nblk
	        do iv = 1, nsys
		  mat(i,iv,l,q,3) = mat(i,iv,l,q,3) + mult(iv) * mat(i,iv,nc,q,3)
		end do
	      end do
	      do q = 1, nblk
	        do iv = 1, nsys
		  mat(i,iv,l,q,4) = mat(i,iv,l,q,4) + mult(iv) * mat(i,iv,nc,q,4)
		  mat(i,iv,l,q,5) = mat(i,iv,l,q,5) + mult(iv) * mat(i,iv,nc,q,5)
	          per(i,iv,l,q,1) = per(i,iv,l,q,1) + mult(iv) * per(i,iv,nc,q,1)
	          per(i,iv,l,q,2) = per(i,iv,l,q,2) + mult(iv) * per(i,iv,nc,q,2)
		end do
	      end do
	    end do	      

	  end do	! end of loop on nc

!.... for all columns
	  
	  nc = 1
	    
	  p = i + 1
	  r = i + 2
	  
	  do iv = 1, nsys
	    fact(iv) = one / mat(i,iv,nc,nc,3)
	  end do
	  do m = 1, nblk
	    do iv = 1, nsys
	      mult(iv) = -mat(p,iv,m,nc,2) * fact(iv)
	      mat(p,iv,m,nc,2) = mult(iv)
	    end do
	    do q = nc + 1, nblk
	      do iv = 1, nsys
		mat(p,iv,m,q,2) = mat(p,iv,m,q,2) + mult(iv) * mat(i,iv,nc,q,3)
	      end do
	    end do
	    do q = 1, nblk
	      do iv = 1, nsys
		mat(p,iv,m,q,3) = mat(p,iv,m,q,3) + mult(iv) * mat(i,iv,nc,q,4)
		mat(p,iv,m,q,4) = mat(p,iv,m,q,4) + mult(iv) * mat(i,iv,nc,q,5)
		per(p,iv,m,q,1) = per(p,iv,m,q,1) + mult(iv) * per(i,iv,nc,q,1)
		per(p,iv,m,q,2) = per(p,iv,m,q,2) + mult(iv) * per(i,iv,nc,q,2)
	      end do
	    end do
	    
	    do iv = 1, nsys
	      mult(iv) = -mat(r,iv,m,nc,1) * fact(iv)
	      mat(r,iv,m,nc,1) = mult(iv)
	    end do
	    do q = nc + 1, nblk
	      do iv = 1, nsys
		mat(r,iv,m,q,1) = mat(r,iv,m,q,1) + mult(iv) * mat(i,iv,nc,q,3)
	      end do
	    end do
	    do q = 1, nblk
	      do iv = 1, nsys
		mat(r,iv,m,q,2) = mat(r,iv,m,q,2) + mult(iv) * mat(i,iv,nc,q,4)
		mat(r,iv,m,q,3) = mat(r,iv,m,q,3) + mult(iv) * mat(i,iv,nc,q,5)
		per(r,iv,m,q,1) = mult(iv) * per(i,iv,nc,q,1)
		per(r,iv,m,q,2) = mult(iv) * per(i,iv,nc,q,2)
	      end do
	    end do

!.... the periodic stuff on the next-to-the-last row

	    do iv = 1, nsys
	      mult(iv) = -per2(n-1,iv,m,nc,i) * fact(iv)
	      per2(n-1,iv,m,nc,i) = mult(iv)
	    end do
	    do q = nc + 1, nblk
	      do iv = 1, nsys
		per2(n-1,iv,m,q,i) = per2(n-1,iv,m,q,i) + mult(iv) * mat(i,iv,nc,q,3)
	      end do
	    end do
	    do q = 1, nblk
	      do iv = 1, nsys
		per2(n-1,iv,m,q,p) = per2(n-1,iv,m,q,p) + mult(iv) * mat(i,iv,nc,q,4)
		per2(n-1,iv,m,q,r) = mult(iv) * mat(i,iv,nc,q,5)
		mat(n-1,iv,m,q,3) = mat(n-1,iv,m,q,3) + mult(iv) * per(i,iv,nc,q,1)
		mat(n-1,iv,m,q,4) = mat(n-1,iv,m,q,4) + mult(iv) * per(i,iv,nc,q,2)
	      end do
	    end do
      
!.... the periodic stuff on the last row

	    do iv = 1, nsys
	      mult(iv) = -per2(n,iv,m,nc,i) * fact(iv)
	      per2(n,iv,m,nc,i) = mult(iv)
	    end do
	    do q = nc + 1, nblk
	      do iv = 1, nsys
		per2(n,iv,m,q,i) = per2(n,iv,m,q,i) + mult(iv) * mat(i,iv,nc,q,3)
	      end do
	    end do
	    do q = 1, nblk
	      do iv = 1, nsys
		per2(n,iv,m,q,p) = per2(n,iv,m,q,p) + mult(iv) * mat(i,iv,nc,q,4)
		per2(n,iv,m,q,r) = mult(iv) * mat(i,iv,nc,q,5)
		mat(n,iv,m,q,2) = mat(n,iv,m,q,2) + mult(iv) * per(i,iv,nc,q,1)
		mat(n,iv,m,q,3) = mat(n,iv,m,q,3) + mult(iv) * per(i,iv,nc,q,2)
	      end do
	    end do

	  end do	! end of loop on m

	  do nc = 2, nblk

	    p = i + 1
	    r = i + 2
	    
	    do iv = 1, nsys
	      fact(iv) = one / mat(i,iv,nc,nc,3)
	    end do
	    do m = 1, nblk
	      do iv = 1, nsys
		mult(iv) = -mat(p,iv,m,nc,2) * fact(iv)
		mat(p,iv,m,nc,2) = mult(iv)
	      end do
	      do q = nc + 1, nblk
                do iv = 1, nsys
		  mat(p,iv,m,q,2) = mat(p,iv,m,q,2) + mult(iv) * mat(i,iv,nc,q,3)
		end do
	      end do
	      do q = 1, nblk
	        do iv = 1, nsys
		  mat(p,iv,m,q,3) = mat(p,iv,m,q,3) + mult(iv) * mat(i,iv,nc,q,4)
		  mat(p,iv,m,q,4) = mat(p,iv,m,q,4) + mult(iv) * mat(i,iv,nc,q,5)
		  per(p,iv,m,q,1) = per(p,iv,m,q,1) + mult(iv) * per(i,iv,nc,q,1)
		  per(p,iv,m,q,2) = per(p,iv,m,q,2) + mult(iv) * per(i,iv,nc,q,2)
		end do
	      end do
	      
	      do iv = 1, nsys
		mult(iv) = -mat(r,iv,m,nc,1) * fact(iv)
		mat(r,iv,m,nc,1) = mult(iv)
	      end do
	      do q = nc + 1, nblk
	        do iv = 1, nsys
		  mat(r,iv,m,q,1) = mat(r,iv,m,q,1) + mult(iv) * mat(i,iv,nc,q,3)
		end do
	      end do
	      do q = 1, nblk
	        do iv = 1, nsys
		  mat(r,iv,m,q,2) = mat(r,iv,m,q,2) + mult(iv) * mat(i,iv,nc,q,4)
		  mat(r,iv,m,q,3) = mat(r,iv,m,q,3) + mult(iv) * mat(i,iv,nc,q,5)
		  per(r,iv,m,q,1) = per(r,iv,m,q,1) + mult(iv) * per(i,iv,nc,q,1)
		  per(r,iv,m,q,2) = per(r,iv,m,q,2) + mult(iv) * per(i,iv,nc,q,2)
		end do
	      end do

!.... the periodic stuff on the next-to-the-last row

	      do iv = 1, nsys
		mult(iv) = -per2(n-1,iv,m,nc,i) * fact(iv)
		per2(n-1,iv,m,nc,i) = mult(iv)
	      end do
	      do q = nc + 1, nblk
	        do iv = 1, nsys
		  per2(n-1,iv,m,q,i) = per2(n-1,iv,m,q,i) + mult(iv) * mat(i,iv,nc,q,3)
		end do
	      end do
	      do q = 1, nblk
	        do iv = 1, nsys
		  per2(n-1,iv,m,q,p) = per2(n-1,iv,m,q,p) + mult(iv) * mat(i,iv,nc,q,4)
		  per2(n-1,iv,m,q,r) = per2(n-1,iv,m,q,r) + mult(iv) * mat(i,iv,nc,q,5)
		  mat(n-1,iv,m,q,3) = mat(n-1,iv,m,q,3) + mult(iv) * per(i,iv,nc,q,1)
		  mat(n-1,iv,m,q,4) = mat(n-1,iv,m,q,4) + mult(iv) * per(i,iv,nc,q,2)
		end do
	      end do
	
!.... the periodic stuff on the last row

	      do iv = 1, nsys
		mult(iv) = -per2(n,iv,m,nc,i) * fact(iv)
		per2(n,iv,m,nc,i) = mult(iv)
	      end do
	      do q = nc + 1, nblk
	        do iv = 1, nsys
		  per2(n,iv,m,q,i) = per2(n,iv,m,q,i) + mult(iv) * mat(i,iv,nc,q,3)
		end do
	      end do
	      do q = 1, nblk
	        do iv = 1, nsys
		  per2(n,iv,m,q,p) = per2(n,iv,m,q,p) + mult(iv) * mat(i,iv,nc,q,4)
		  per2(n,iv,m,q,r) = per2(n,iv,m,q,r) + mult(iv) * mat(i,iv,nc,q,5)
		  mat(n,iv,m,q,2) = mat(n,iv,m,q,2) + mult(iv) * per(i,iv,nc,q,1)
		  mat(n,iv,m,q,3) = mat(n,iv,m,q,3) + mult(iv) * per(i,iv,nc,q,2)
		end do
	      end do

	    end do	! end of loop on m

	  end do	! end of loop on nc
	  
	end do		! end of loop on i

!====================================================================================
 	i = n - 5

!.... do the first four columns

	  do nc = 1, nblk - 1
	  
	    do m = nc, nblk - 1
	      l = m + 1
	      do iv = 1, nsys
		mult(iv) = -mat(i,iv,l,nc,3) / mat(i,iv,nc,nc,3)
		mat(i,iv,l,nc,3) = mult(iv)
	      end do
	      do q = nc + 1, nblk
	        do iv = 1, nsys
		  mat(i,iv,l,q,3) = mat(i,iv,l,q,3) + mult(iv) * mat(i,iv,nc,q,3)
		end do
	      end do
	      do q = 1, nblk
	        do iv = 1, nsys
		  mat(i,iv,l,q,4) = mat(i,iv,l,q,4) + mult(iv) * mat(i,iv,nc,q,4)
		  mat(i,iv,l,q,5) = mat(i,iv,l,q,5) + mult(iv) * mat(i,iv,nc,q,5)
	          per(i,iv,l,q,1) = per(i,iv,l,q,1) + mult(iv) * per(i,iv,nc,q,1)
	          per(i,iv,l,q,2) = per(i,iv,l,q,2) + mult(iv) * per(i,iv,nc,q,2)
		end do
	      end do
	    end do	      

	  end do	! end of loop on nc

!.... for all columns
	  
	  nc = 1
	  
	  p = i + 1
	  r = i + 2
	  
	  do iv = 1, nsys
	    fact(iv) = one / mat(i,iv,nc,nc,3)
	  end do
	  do m = 1, nblk
	    do iv = 1, nsys
	      mult(iv) = -mat(p,iv,m,nc,2) * fact(iv)
	      mat(p,iv,m,nc,2) = mult(iv)
	    end do
	    do q = nc + 1, nblk
	      do iv = 1, nsys
		mat(p,iv,m,q,2) = mat(p,iv,m,q,2) + mult(iv) * mat(i,iv,nc,q,3)
	      end do
	    end do
	    do q = 1, nblk
	      do iv = 1, nsys
		mat(p,iv,m,q,3) = mat(p,iv,m,q,3) + mult(iv) * mat(i,iv,nc,q,4)
		mat(p,iv,m,q,4) = mat(p,iv,m,q,4) + mult(iv) * mat(i,iv,nc,q,5)
		per(p,iv,m,q,1) = per(p,iv,m,q,1) + mult(iv) * per(i,iv,nc,q,1)
		per(p,iv,m,q,2) = per(p,iv,m,q,2) + mult(iv) * per(i,iv,nc,q,2)
	      end do
	    end do
	    
	    do iv = 1, nsys
	      mult(iv) = -mat(r,iv,m,nc,1) * fact(iv)
	      mat(r,iv,m,nc,1) = mult(iv)
	    end do
	    do q = nc + 1, nblk
	      do iv = 1, nsys
		mat(r,iv,m,q,1) = mat(r,iv,m,q,1) + mult(iv) * mat(i,iv,nc,q,3)
	      end do
	    end do
	    do q = 1, nblk
	      do iv = 1, nsys
		mat(r,iv,m,q,2) = mat(r,iv,m,q,2) + mult(iv) * mat(i,iv,nc,q,4)
		mat(r,iv,m,q,3) = mat(r,iv,m,q,3) + mult(iv) * mat(i,iv,nc,q,5)
		mat(r,iv,m,q,4) = mat(r,iv,m,q,4) + mult(iv) * per(i,iv,nc,q,1)
		per(r,iv,m,q,2) = mult(iv) * per(i,iv,nc,q,2)
	      end do
	    end do

!.... the periodic stuff on the next-to-the-last row

	    do iv = 1, nsys
	      mult(iv) = -per2(n-1,iv,m,nc,i) * fact(iv)
	      per2(n-1,iv,m,nc,i) = mult(iv)
	    end do
	    do q = nc + 1, nblk
	      do iv = 1, nsys
		per2(n-1,iv,m,q,i) = per2(n-1,iv,m,q,i) + mult(iv) * mat(i,iv,nc,q,3)
	      end do
	    end do
	    do q = 1, nblk
	      do iv = 1, nsys
		per2(n-1,iv,m,q,p) = per2(n-1,iv,m,q,p) + mult(iv) * mat(i,iv,nc,q,4)
		mat(n-1,iv,m,q,1) = mat(n-1,iv,m,q,1) + mult(iv) * mat(i,iv,nc,q,5)
		mat(n-1,iv,m,q,3) = mat(n-1,iv,m,q,3) + mult(iv) * per(i,iv,nc,q,1)
		mat(n-1,iv,m,q,4) = mat(n-1,iv,m,q,4) + mult(iv) * per(i,iv,nc,q,2)
	      end do
	    end do
      
!.... the periodic stuff on the last row

	    do iv = 1, nsys
	      mult(iv) = -per2(n,iv,m,nc,i) * fact(iv)
	      per2(n,iv,m,nc,i) = mult(iv)
	    end do
	    do q = nc + 1, nblk
	      do iv = 1, nsys
		per2(n,iv,m,q,i) = per2(n,iv,m,q,i) + mult(iv) * mat(i,iv,nc,q,3)
	      end do
	    end do
	    do q = 1, nblk
	      do iv = 1, nsys
		per2(n,iv,m,q,p) = per2(n,iv,m,q,p) + mult(iv) * mat(i,iv,nc,q,4)
		per2(n,iv,m,q,r) = mult(iv) * mat(i,iv,nc,q,5)
		mat(n,iv,m,q,2) = mat(n,iv,m,q,2) + mult(iv) * per(i,iv,nc,q,1)
		mat(n,iv,m,q,3) = mat(n,iv,m,q,3) + mult(iv) * per(i,iv,nc,q,2)
	      end do
	    end do

	  end do	! end of loop on m

	  do nc = 2, nblk

	    p = i + 1
	    r = i + 2
	    
	    do iv = 1, nsys
	      fact(iv) = one / mat(i,iv,nc,nc,3)
	    end do
	    do m = 1, nblk
	      do iv = 1, nsys
		mult(iv) = -mat(p,iv,m,nc,2) * fact(iv)
		mat(p,iv,m,nc,2) = mult(iv)
	      end do
	      do q = nc + 1, nblk
                do iv = 1, nsys
		  mat(p,iv,m,q,2) = mat(p,iv,m,q,2) + mult(iv) * mat(i,iv,nc,q,3)
		end do
	      end do
	      do q = 1, nblk
	        do iv = 1, nsys
		  mat(p,iv,m,q,3) = mat(p,iv,m,q,3) + mult(iv) * mat(i,iv,nc,q,4)
		  mat(p,iv,m,q,4) = mat(p,iv,m,q,4) + mult(iv) * mat(i,iv,nc,q,5)
		  per(p,iv,m,q,1) = per(p,iv,m,q,1) + mult(iv) * per(i,iv,nc,q,1)
		  per(p,iv,m,q,2) = per(p,iv,m,q,2) + mult(iv) * per(i,iv,nc,q,2)
		end do
	      end do
	      
	      do iv = 1, nsys
		mult(iv) = -mat(r,iv,m,nc,1) * fact(iv)
		mat(r,iv,m,nc,1) = mult(iv)
	      end do
	      do q = nc + 1, nblk
	        do iv = 1, nsys
		  mat(r,iv,m,q,1) = mat(r,iv,m,q,1) + mult(iv) * mat(i,iv,nc,q,3)
		end do
	      end do
	      do q = 1, nblk
	        do iv = 1, nsys
		  mat(r,iv,m,q,2) = mat(r,iv,m,q,2) + mult(iv) * mat(i,iv,nc,q,4)
		  mat(r,iv,m,q,3) = mat(r,iv,m,q,3) + mult(iv) * mat(i,iv,nc,q,5)
		  mat(r,iv,m,q,4) = mat(r,iv,m,q,4) + mult(iv) * per(i,iv,nc,q,1)
		  per(r,iv,m,q,2) = per(r,iv,m,q,2) + mult(iv) * per(i,iv,nc,q,2)
		end do
	      end do

!.... the periodic stuff on the next-to-the-last row

	      do iv = 1, nsys
		mult(iv) = -per2(n-1,iv,m,nc,i) * fact(iv)
		per2(n-1,iv,m,nc,i) = mult(iv)
	      end do
	      do q = nc + 1, nblk
	        do iv = 1, nsys
		  per2(n-1,iv,m,q,i) = per2(n-1,iv,m,q,i) + mult(iv) * mat(i,iv,nc,q,3)
		end do
	      end do
	      do q = 1, nblk
	        do iv = 1, nsys
		  per2(n-1,iv,m,q,p) = per2(n-1,iv,m,q,p) + mult(iv) * mat(i,iv,nc,q,4)
		  mat(n-1,iv,m,q,1) = mat(n-1,iv,m,q,1) + mult(iv) * mat(i,iv,nc,q,5)
		  mat(n-1,iv,m,q,3) = mat(n-1,iv,m,q,3) + mult(iv) * per(i,iv,nc,q,1)
		  mat(n-1,iv,m,q,4) = mat(n-1,iv,m,q,4) + mult(iv) * per(i,iv,nc,q,2)
		end do
	      end do
	
!.... the periodic stuff on the last row

	      do iv = 1, nsys
		mult(iv) = -per2(n,iv,m,nc,i) * fact(iv)
		per2(n,iv,m,nc,i) = mult(iv)
	      end do
	      do q = nc + 1, nblk
	        do iv = 1, nsys
		  per2(n,iv,m,q,i) = per2(n,iv,m,q,i) + mult(iv) * mat(i,iv,nc,q,3)
		end do
	      end do
	      do q = 1, nblk
	        do iv = 1, nsys
		  per2(n,iv,m,q,p) = per2(n,iv,m,q,p) + mult(iv) * mat(i,iv,nc,q,4)
		  per2(n,iv,m,q,r) = per2(n,iv,m,q,r) + mult(iv) * mat(i,iv,nc,q,5)
		  mat(n,iv,m,q,2) = mat(n,iv,m,q,2) + mult(iv) * per(i,iv,nc,q,1)
		  mat(n,iv,m,q,3) = mat(n,iv,m,q,3) + mult(iv) * per(i,iv,nc,q,2)
		end do
	      end do

	    end do	! end of loop on m

	  end do	! end of loop on nc
	  
!====================================================================================
 	i = n - 4

!.... do the first four columns

	  do nc = 1, nblk - 1
	  
	    do m = nc, nblk - 1
	      l = m + 1
	      do iv = 1, nsys
		mult(iv) = -mat(i,iv,l,nc,3) / mat(i,iv,nc,nc,3)
		mat(i,iv,l,nc,3) = mult(iv)
	      end do
	      do q = nc + 1, nblk
	        do iv = 1, nsys
		  mat(i,iv,l,q,3) = mat(i,iv,l,q,3) + mult(iv) * mat(i,iv,nc,q,3)
		end do
	      end do
	      do q = 1, nblk
	        do iv = 1, nsys
		  mat(i,iv,l,q,4) = mat(i,iv,l,q,4) + mult(iv) * mat(i,iv,nc,q,4)
		  mat(i,iv,l,q,5) = mat(i,iv,l,q,5) + mult(iv) * mat(i,iv,nc,q,5)
	          per(i,iv,l,q,1) = per(i,iv,l,q,1) + mult(iv) * per(i,iv,nc,q,1)
	          per(i,iv,l,q,2) = per(i,iv,l,q,2) + mult(iv) * per(i,iv,nc,q,2)
		end do
	      end do
	    end do	      

	  end do	! end of loop on nc

!.... for all columns
	  
	  do nc = 1, nblk

	    p = i + 1
	    r = i + 2
	    
	    do iv = 1, nsys
	      fact(iv) = one / mat(i,iv,nc,nc,3)
	    end do
	    do m = 1, nblk
	      do iv = 1, nsys
		mult(iv) = -mat(p,iv,m,nc,2) * fact(iv)
		mat(p,iv,m,nc,2) = mult(iv)
	      end do
	      do q = nc + 1, nblk
                do iv = 1, nsys
		  mat(p,iv,m,q,2) = mat(p,iv,m,q,2) + mult(iv) * mat(i,iv,nc,q,3)
		end do
	      end do
	      do q = 1, nblk
	        do iv = 1, nsys
		  mat(p,iv,m,q,3) = mat(p,iv,m,q,3) + mult(iv) * mat(i,iv,nc,q,4)
		  mat(p,iv,m,q,4) = mat(p,iv,m,q,4) + mult(iv) * mat(i,iv,nc,q,5)
		  mat(p,iv,m,q,5) = mat(p,iv,m,q,5) + mult(iv) * per(i,iv,nc,q,1)
		  per(p,iv,m,q,2) = per(p,iv,m,q,2) + mult(iv) * per(i,iv,nc,q,2)
		end do
	      end do
	      
	      do iv = 1, nsys
		mult(iv) = -mat(r,iv,m,nc,1) * fact(iv)
		mat(r,iv,m,nc,1) = mult(iv)
	      end do
	      do q = nc + 1, nblk
	        do iv = 1, nsys
		  mat(r,iv,m,q,1) = mat(r,iv,m,q,1) + mult(iv) * mat(i,iv,nc,q,3)
		end do
	      end do
	      do q = 1, nblk
	        do iv = 1, nsys
		  mat(r,iv,m,q,2) = mat(r,iv,m,q,2) + mult(iv) * mat(i,iv,nc,q,4)
		  mat(r,iv,m,q,3) = mat(r,iv,m,q,3) + mult(iv) * mat(i,iv,nc,q,5)
		  mat(r,iv,m,q,4) = mat(r,iv,m,q,4) + mult(iv) * per(i,iv,nc,q,1)
		  mat(r,iv,m,q,5) = mat(r,iv,m,q,5) + mult(iv) * per(i,iv,nc,q,2)
		end do
	      end do

!.... the periodic stuff on the next-to-the-last row

	      do iv = 1, nsys
		mult(iv) = -per2(n-1,iv,m,nc,i) * fact(iv)
		per2(n-1,iv,m,nc,i) = mult(iv)
	      end do
	      do q = nc + 1, nblk
	        do iv = 1, nsys
		  per2(n-1,iv,m,q,i) = per2(n-1,iv,m,q,i) + mult(iv) * mat(i,iv,nc,q,3)
		end do
	      end do
	      do q = 1, nblk
	        do iv = 1, nsys
		  mat(n-1,iv,m,q,1) = mat(n-1,iv,m,q,1) + mult(iv) * mat(i,iv,nc,q,4)
		  mat(n-1,iv,m,q,2) = mat(n-1,iv,m,q,2) + mult(iv) * mat(i,iv,nc,q,5)
		  mat(n-1,iv,m,q,3) = mat(n-1,iv,m,q,3) + mult(iv) * per(i,iv,nc,q,1)
		  mat(n-1,iv,m,q,4) = mat(n-1,iv,m,q,4) + mult(iv) * per(i,iv,nc,q,2)
		end do
	      end do
	
!.... the periodic stuff on the last row

	      do iv = 1, nsys
		mult(iv) = -per2(n,iv,m,nc,i) * fact(iv)
		per2(n,iv,m,nc,i) = mult(iv)
	      end do
	      do q = nc + 1, nblk
	        do iv = 1, nsys
		  per2(n,iv,m,q,i) = per2(n,iv,m,q,i) + mult(iv) * mat(i,iv,nc,q,3)
		end do
	      end do
	      do q = 1, nblk
	        do iv = 1, nsys
		  per2(n,iv,m,q,p) = per2(n,iv,m,q,p) + mult(iv) * mat(i,iv,nc,q,4)
		  mat(n,iv,m,q,1) = mat(n,iv,m,q,1) + mult(iv) * mat(i,iv,nc,q,5)
		  mat(n,iv,m,q,2) = mat(n,iv,m,q,2) + mult(iv) * per(i,iv,nc,q,1)
		  mat(n,iv,m,q,3) = mat(n,iv,m,q,3) + mult(iv) * per(i,iv,nc,q,2)
		end do
	      end do

	    end do	! end of loop on m

	  end do	! end of loop on nc
	  
!====================================================================================
 	i = n - 3

!.... do the first four columns

	  do nc = 1, nblk - 1
	  
	    do m = nc, nblk - 1
	      l = m + 1
	      do iv = 1, nsys
		mult(iv) = -mat(i,iv,l,nc,3) / mat(i,iv,nc,nc,3)
		mat(i,iv,l,nc,3) = mult(iv)
	      end do
	      do q = nc + 1, nblk
	        do iv = 1, nsys
		  mat(i,iv,l,q,3) = mat(i,iv,l,q,3) + mult(iv) * mat(i,iv,nc,q,3)
		end do
	      end do
	      do q = 1, nblk
	        do iv = 1, nsys
		  mat(i,iv,l,q,4) = mat(i,iv,l,q,4) + mult(iv) * mat(i,iv,nc,q,4)
		  mat(i,iv,l,q,5) = mat(i,iv,l,q,5) + mult(iv) * mat(i,iv,nc,q,5)
	          per(i,iv,l,q,2) = per(i,iv,l,q,2) + mult(iv) * per(i,iv,nc,q,2)
		end do
	      end do
	    end do	      

	  end do	! end of loop on nc

!.... for all columns
	  
	  do nc = 1, nblk

	    p = i + 1
	    r = i + 2
	    
	    do iv = 1, nsys
	      fact(iv) = one / mat(i,iv,nc,nc,3)
	    end do
	    do m = 1, nblk
	      do iv = 1, nsys
		mult(iv) = -mat(p,iv,m,nc,2) * fact(iv)
		mat(p,iv,m,nc,2) = mult(iv)
	      end do
	      do q = nc + 1, nblk
                do iv = 1, nsys
		  mat(p,iv,m,q,2) = mat(p,iv,m,q,2) + mult(iv) * mat(i,iv,nc,q,3)
		end do
	      end do
	      do q = 1, nblk
	        do iv = 1, nsys
		  mat(p,iv,m,q,3) = mat(p,iv,m,q,3) + mult(iv) * mat(i,iv,nc,q,4)
		  mat(p,iv,m,q,4) = mat(p,iv,m,q,4) + mult(iv) * mat(i,iv,nc,q,5)
		  mat(p,iv,m,q,5) = mat(p,iv,m,q,5) + mult(iv) * per(i,iv,nc,q,2)
		end do
	      end do
	      
	      do iv = 1, nsys
		mult(iv) = -mat(r,iv,m,nc,1) * fact(iv)
		mat(r,iv,m,nc,1) = mult(iv)
	      end do
	      do q = nc + 1, nblk
	        do iv = 1, nsys
		  mat(r,iv,m,q,1) = mat(r,iv,m,q,1) + mult(iv) * mat(i,iv,nc,q,3)
		end do
	      end do
	      do q = 1, nblk
	        do iv = 1, nsys
		  mat(r,iv,m,q,2) = mat(r,iv,m,q,2) + mult(iv) * mat(i,iv,nc,q,4)
		  mat(r,iv,m,q,3) = mat(r,iv,m,q,3) + mult(iv) * mat(i,iv,nc,q,5)
		  mat(r,iv,m,q,4) = mat(r,iv,m,q,4) + mult(iv) * per(i,iv,nc,q,2)
		end do
	      end do

!.... the periodic stuff on the last row

	      do iv = 1, nsys
		mult(iv) = -per2(n,iv,m,nc,i) * fact(iv)
		per2(n,iv,m,nc,i) = mult(iv)
	      end do
	      do q = nc + 1, nblk
	        do iv = 1, nsys
		  per2(n,iv,m,q,i) = per2(n,iv,m,q,i) + mult(iv) * mat(i,iv,nc,q,3)
		end do
	      end do
	      do q = 1, nblk
	        do iv = 1, nsys
		  mat(n,iv,m,q,1) = mat(n,iv,m,q,1) + mult(iv) * mat(i,iv,nc,q,4)
		  mat(n,iv,m,q,2) = mat(n,iv,m,q,2) + mult(iv) * mat(i,iv,nc,q,5)
		  mat(n,iv,m,q,3) = mat(n,iv,m,q,3) + mult(iv) * per(i,iv,nc,q,2)
		end do
	      end do

	    end do	! end of loop on m

	  end do	! end of loop on nc
	  
!====================================================================================
 	i = n - 2

!.... do the first four columns

	  do nc = 1, nblk - 1
	  
	    do m = nc, nblk - 1
	      l = m + 1
	      do iv = 1, nsys
		mult(iv) = -mat(i,iv,l,nc,3) / mat(i,iv,nc,nc,3)
		mat(i,iv,l,nc,3) = mult(iv)
	      end do
	      do q = nc + 1, nblk
	        do iv = 1, nsys
		  mat(i,iv,l,q,3) = mat(i,iv,l,q,3) + mult(iv) * mat(i,iv,nc,q,3)
		end do
	      end do
	      do q = 1, nblk
	        do iv = 1, nsys
		  mat(i,iv,l,q,4) = mat(i,iv,l,q,4) + mult(iv) * mat(i,iv,nc,q,4)
		  mat(i,iv,l,q,5) = mat(i,iv,l,q,5) + mult(iv) * mat(i,iv,nc,q,5)
		end do
	      end do
	    end do	      

	  end do	! end of loop on nc

!.... for all columns
	  
	  do nc = 1, nblk

	    p = i + 1
	    r = i + 2
	    
	    do iv = 1, nsys
	      fact(iv) = one / mat(i,iv,nc,nc,3)
	    end do
	    do m = 1, nblk
	      do iv = 1, nsys
		mult(iv) = -mat(p,iv,m,nc,2) * fact(iv)
		mat(p,iv,m,nc,2) = mult(iv)
	      end do
	      do q = nc + 1, nblk
                do iv = 1, nsys
		  mat(p,iv,m,q,2) = mat(p,iv,m,q,2) + mult(iv) * mat(i,iv,nc,q,3)
		end do
	      end do
	      do q = 1, nblk
	        do iv = 1, nsys
		  mat(p,iv,m,q,3) = mat(p,iv,m,q,3) + mult(iv) * mat(i,iv,nc,q,4)
		  mat(p,iv,m,q,4) = mat(p,iv,m,q,4) + mult(iv) * mat(i,iv,nc,q,5)
		end do
	      end do
	      
	      do iv = 1, nsys
		mult(iv) = -mat(r,iv,m,nc,1) * fact(iv)
		mat(r,iv,m,nc,1) = mult(iv)
	      end do
	      do q = nc + 1, nblk
	        do iv = 1, nsys
		  mat(r,iv,m,q,1) = mat(r,iv,m,q,1) + mult(iv) * mat(i,iv,nc,q,3)
		end do
	      end do
	      do q = 1, nblk
	        do iv = 1, nsys
		  mat(r,iv,m,q,2) = mat(r,iv,m,q,2) + mult(iv) * mat(i,iv,nc,q,4)
		  mat(r,iv,m,q,3) = mat(r,iv,m,q,3) + mult(iv) * mat(i,iv,nc,q,5)
		end do
	      end do

	    end do	! end of loop on m

	  end do	! end of loop on nc

!====================================================================================
	i = n - 1

!.... do the first four columns

	  do nc = 1, nblk - 1
	  
	    do m = nc, nblk - 1
	      l = m + 1
	      do iv = 1, nsys
		mult(iv) = -mat(i,iv,l,nc,3) / mat(i,iv,nc,nc,3)
		mat(i,iv,l,nc,3) = mult(iv)
	      end do
	      do q = nc + 1, nblk
                do iv = 1, nsys
		  mat(i,iv,l,q,3) = mat(i,iv,l,q,3) + mult(iv) * mat(i,iv,nc,q,3)
		end do
	      end do
	      do q = 1, nblk
	        do iv = 1, nsys
		  mat(i,iv,l,q,4) = mat(i,iv,l,q,4) + mult(iv) * mat(i,iv,nc,q,4)
		end do
	      end do
	    end do	      

	  end do	! end of loop on nc

!.... for all columns
	  
	  do nc = 1, nblk

	    p = i + 1
	    
	    do iv = 1, nsys
	      fact(iv) = one / mat(i,iv,nc,nc,3)
	    end do
	    
	    do m = 1, nblk
	      do iv = 1, nsys
		mult(iv) = -mat(p,iv,m,nc,2) * fact(iv)
		mat(p,iv,m,nc,2) = mult(iv)
	      end do
	      do q = nc + 1, nblk
	        do iv = 1, nsys
		  mat(p,iv,m,q,2) = mat(p,iv,m,q,2) + mult(iv) * mat(i,iv,nc,q,3)
		end do
	      end do
	      do q = 1, nblk
	        do iv = 1, nsys
		  mat(p,iv,m,q,3) = mat(p,iv,m,q,3) + mult(iv) * mat(i,iv,nc,q,4)
		end do
	      end do
	    end do

	  end do	! end of loop on nc

!.... do the last row

	i = n
	
	do nc = 1, nblk - 1
	
	  do m = nc, nblk - 1
	    l = m + 1
	    do iv = 1, nsys
	      mult(iv) = -mat(i,iv,l,nc,3) / mat(i,iv,nc,nc,3)
	      mat(i,iv,l,nc,3) = mult(iv)
	    end do
	    do q = nc + 1, nblk
	      do iv = 1, nsys
		mat(i,iv,l,q,3) = mat(i,iv,l,q,3) + mult(iv) * mat(i,iv,nc,q,3)
	      end do
	    end do
	  end do	      

	end do	! end of loop on nc

	end if
!=============================================================================!
!                     F O R W A R D   E L I M I N A T I O N 
!=============================================================================!
	if (code.eq.1 .or. code.eq.2) then
	
 	do i = 1, n - 6

!.... do the first four columns

	  do nc = 1, nblk - 1
	  
	    do m = nc, nblk - 1
	      l = m + 1
	      do iv = 1, nsys
		mult(iv) = mat(i,iv,l,nc,3)
		rhs(i,iv,l) = rhs(i,iv,l) + mult(iv) * rhs(i,iv,nc)
	      end do
	    end do	      

	  end do	! end of loop on nc

!.... for all columns
	  
	  nc = 1
	    
	  p = i + 1
	  r = i + 2
	  
	  do m = 1, nblk
	    do iv = 1, nsys
	      mult(iv) = mat(p,iv,m,nc,2)
	      rhs(p,iv,m) = rhs(p,iv,m) + mult(iv) * rhs(i,iv,nc)
	    end do
	    
	    do iv = 1, nsys
	      mult(iv) = mat(r,iv,m,nc,1)
	      rhs(r,iv,m) = rhs(r,iv,m) + mult(iv) * rhs(i,iv,nc)
	    end do

!.... the periodic stuff on the next-to-the-last row

	    do iv = 1, nsys
	      mult(iv) = per2(n-1,iv,m,nc,i)
	      rhs(n-1,iv,m) = rhs(n-1,iv,m) + mult(iv) * rhs(i,iv,nc)
	    end do
      
!.... the periodic stuff on the last row

	    do iv = 1, nsys
	      mult(iv) = per2(n,iv,m,nc,i)
	      rhs(n,iv,m) = rhs(n,iv,m) + mult(iv) * rhs(i,iv,nc)
	    end do

	  end do	! end of loop on m

	  do nc = 2, nblk

	    p = i + 1
	    r = i + 2
	    
	    do m = 1, nblk
	      do iv = 1, nsys
		mult(iv) = mat(p,iv,m,nc,2)
		rhs(p,iv,m) = rhs(p,iv,m) + mult(iv) * rhs(i,iv,nc)
	      end do
	      
	      do iv = 1, nsys
		mult(iv) = mat(r,iv,m,nc,1)
		rhs(r,iv,m) = rhs(r,iv,m) + mult(iv) * rhs(i,iv,nc)
	      end do

!.... the periodic stuff on the next-to-the-last row

	      do iv = 1, nsys
		mult(iv) = per2(n-1,iv,m,nc,i)
		rhs(n-1,iv,m) = rhs(n-1,iv,m) + mult(iv) * rhs(i,iv,nc)
	      end do
	
!.... the periodic stuff on the last row

	      do iv = 1, nsys
		mult(iv) = per2(n,iv,m,nc,i)
		rhs(n,iv,m) = rhs(n,iv,m) + mult(iv) * rhs(i,iv,nc)
	      end do

	    end do	! end of loop on m

	  end do	! end of loop on nc
	  
	end do		! end of loop on i

!====================================================================================
 	i = n - 5

!.... do the first four columns

	  do nc = 1, nblk - 1
	  
	    do m = nc, nblk - 1
	      l = m + 1
	      do iv = 1, nsys
		mult(iv) = mat(i,iv,l,nc,3)
		rhs(i,iv,l) = rhs(i,iv,l) + mult(iv) * rhs(i,iv,nc)
	      end do
	    end do	      

	  end do	! end of loop on nc

!.... for all columns
	  
	  nc = 1
	  
	  p = i + 1
	  r = i + 2
	  
	  do m = 1, nblk
	    do iv = 1, nsys
	      mult(iv) = mat(p,iv,m,nc,2)
	      rhs(p,iv,m) = rhs(p,iv,m) + mult(iv) * rhs(i,iv,nc)
	    end do
	    
	    do iv = 1, nsys
	      mult(iv) = mat(r,iv,m,nc,1)
	      rhs(r,iv,m) = rhs(r,iv,m) + mult(iv) * rhs(i,iv,nc)
	    end do

!.... the periodic stuff on the next-to-the-last row

	    do iv = 1, nsys
	      mult(iv) = per2(n-1,iv,m,nc,i)
	      rhs(n-1,iv,m) = rhs(n-1,iv,m) + mult(iv) * rhs(i,iv,nc)
	    end do
      
!.... the periodic stuff on the last row

	    do iv = 1, nsys
	      mult(iv) = per2(n,iv,m,nc,i)
	      rhs(n,iv,m) = rhs(n,iv,m) + mult(iv) * rhs(i,iv,nc)
	    end do

	  end do	! end of loop on m

	  do nc = 2, nblk

	    p = i + 1
	    r = i + 2
	    
	    do m = 1, nblk
	      do iv = 1, nsys
		mult(iv) = mat(p,iv,m,nc,2)
		rhs(p,iv,m) = rhs(p,iv,m) + mult(iv) * rhs(i,iv,nc)
	      end do
	      
	      do iv = 1, nsys
		mult(iv) = mat(r,iv,m,nc,1)
		rhs(r,iv,m) = rhs(r,iv,m) + mult(iv) * rhs(i,iv,nc)
	      end do

!.... the periodic stuff on the next-to-the-last row

	      do iv = 1, nsys
		mult(iv) = per2(n-1,iv,m,nc,i)
		rhs(n-1,iv,m) = rhs(n-1,iv,m) + mult(iv) * rhs(i,iv,nc)
	      end do
	
!.... the periodic stuff on the last row

	      do iv = 1, nsys
		mult(iv) = per2(n,iv,m,nc,i)
		rhs(n,iv,m) = rhs(n,iv,m) + mult(iv) * rhs(i,iv,nc)
	      end do

	    end do	! end of loop on m

	  end do	! end of loop on nc
	  
!====================================================================================
 	i = n - 4

!.... do the first four columns

	  do nc = 1, nblk - 1
	  
	    do m = nc, nblk - 1
	      l = m + 1
	      do iv = 1, nsys
		mult(iv) = mat(i,iv,l,nc,3)
		rhs(i,iv,l) = rhs(i,iv,l) + mult(iv) * rhs(i,iv,nc)
	      end do
	    end do	      

	  end do	! end of loop on nc

!.... for all columns
	  
	  do nc = 1, nblk

	    p = i + 1
	    r = i + 2
	    
	    do m = 1, nblk
	      do iv = 1, nsys
		mult(iv) = mat(p,iv,m,nc,2)
		rhs(p,iv,m) = rhs(p,iv,m) + mult(iv) * rhs(i,iv,nc)
	      end do
	      
	      do iv = 1, nsys
		mult(iv) = mat(r,iv,m,nc,1)
		rhs(r,iv,m) = rhs(r,iv,m) + mult(iv) * rhs(i,iv,nc)
	      end do

!.... the periodic stuff on the next-to-the-last row

	      do iv = 1, nsys
		mult(iv) = per2(n-1,iv,m,nc,i)
		rhs(n-1,iv,m) = rhs(n-1,iv,m) + mult(iv) * rhs(i,iv,nc)
	      end do
	
!.... the periodic stuff on the last row

	      do iv = 1, nsys
		mult(iv) = per2(n,iv,m,nc,i)
		rhs(n,iv,m) = rhs(n,iv,m) + mult(iv) * rhs(i,iv,nc)
	      end do

	    end do	! end of loop on m

	  end do	! end of loop on nc
	  
!====================================================================================
 	i = n - 3

!.... do the first four columns

	  do nc = 1, nblk - 1
	  
	    do m = nc, nblk - 1
	      l = m + 1
	      do iv = 1, nsys
		mult(iv) = mat(i,iv,l,nc,3)
		rhs(i,iv,l) = rhs(i,iv,l) + mult(iv) * rhs(i,iv,nc)
	      end do
	    end do	      

	  end do	! end of loop on nc

!.... for all columns
	  
	  do nc = 1, nblk

	    p = i + 1
	    r = i + 2
	    
	    do m = 1, nblk
	      do iv = 1, nsys
		mult(iv) = mat(p,iv,m,nc,2)
		rhs(p,iv,m) = rhs(p,iv,m) + mult(iv) * rhs(i,iv,nc)
	      end do
	      
	      do iv = 1, nsys
		mult(iv) = mat(r,iv,m,nc,1)
		rhs(r,iv,m) = rhs(r,iv,m) + mult(iv) * rhs(i,iv,nc)
	      end do

!.... the periodic stuff on the last row

	      do iv = 1, nsys
		mult(iv) = per2(n,iv,m,nc,i)
		rhs(n,iv,m) = rhs(n,iv,m) + mult(iv) * rhs(i,iv,nc)
	      end do

	    end do	! end of loop on m

	  end do	! end of loop on nc
	  
!====================================================================================
 	i = n - 2

!.... do the first four columns

	  do nc = 1, nblk - 1
	  
	    do m = nc, nblk - 1
	      l = m + 1
	      do iv = 1, nsys
		mult(iv) = mat(i,iv,l,nc,3)
		rhs(i,iv,l) = rhs(i,iv,l) + mult(iv) * rhs(i,iv,nc)
	      end do
	    end do	      

	  end do	! end of loop on nc

!.... for all columns
	  
	  do nc = 1, nblk

	    p = i + 1
	    r = i + 2
	    
	    do m = 1, nblk
	      do iv = 1, nsys
		mult(iv) = mat(p,iv,m,nc,2)
		rhs(p,iv,m) = rhs(p,iv,m) + mult(iv) * rhs(i,iv,nc)
	      end do
	      
	      do iv = 1, nsys
		mult(iv) = mat(r,iv,m,nc,1)
		rhs(r,iv,m) = rhs(r,iv,m) + mult(iv) * rhs(i,iv,nc)
	      end do

	    end do	! end of loop on m

	  end do	! end of loop on nc

!====================================================================================
	i = n - 1

!.... do the first four columns

	  do nc = 1, nblk - 1
	  
	    do m = nc, nblk - 1
	      l = m + 1
	      do iv = 1, nsys
		mult(iv) = mat(i,iv,l,nc,3)
		rhs(i,iv,l) = rhs(i,iv,l) + mult(iv) * rhs(i,iv,nc)
	      end do
	    end do	      

	  end do	! end of loop on nc

!.... for all columns
	  
	  do nc = 1, nblk

	    p = i + 1
	    
	    do m = 1, nblk
	      do iv = 1, nsys
		mult(iv) = mat(p,iv,m,nc,2)
		rhs(p,iv,m) = rhs(p,iv,m) + mult(iv) * rhs(i,iv,nc)
	      end do
	    end do

	  end do	! end of loop on nc

!.... do the last row

	i = n
	
	do nc = 1, nblk - 1
	
	  do m = nc, nblk - 1
	    l = m + 1
	    do iv = 1, nsys
	      mult(iv) = mat(i,iv,l,nc,3)
	      rhs(i,iv,l) = rhs(i,iv,l) + mult(iv) * rhs(i,iv,nc)
	    end do
	  end do	      

	end do	! end of loop on nc

!=============================================================================!
!                   B A C K W A R D   S U B S T I T U T I O N 
!=============================================================================!

!.... do the last row first

	i = n
	
	do m = nblk, 1, -1
	  do q = m+1, nblk
	    do iv = 1, nsys
	      rhs(i,iv,m) = rhs(i,iv,m) - rhs(i,iv,q) * mat(i,iv,m,q,3)
	    end do
	  end do
	  do iv = 1, nsys
	    rhs(i,iv,m) = rhs(i,iv,m) / mat(i,iv,m,m,3)
	  end do
	end do

!.... do the next-to-the-last row

	i = n - 1
	p = i + 1
	
	do m = nblk, 1, -1
	  do q = m+1, nblk
	    do iv = 1, nsys
	      rhs(i,iv,m) = rhs(i,iv,m) - rhs(i,iv,q) * mat(i,iv,m,q,3)
	    end do
	  end do
	  do q = 1, nblk
	    do iv = 1, nsys
	      rhs(i,iv,m) = rhs(i,iv,m) - rhs(p,iv,q) * mat(i,iv,m,q,4)
	    end do
	  end do
	  do iv = 1, nsys
	    rhs(i,iv,m) = rhs(i,iv,m) / mat(i,iv,m,m,3)
	  end do
	end do

	i = n - 2
	p = i + 1
	r = i + 2
	
	do m = nblk, 1, -1
	  do q = m+1, nblk
	    do iv = 1, nsys
	      rhs(i,iv,m) = rhs(i,iv,m) - rhs(i,iv,q) * mat(i,iv,m,q,3)
	    end do
	  end do
	  do q = 1, nblk
	    do iv = 1, nsys
	      rhs(i,iv,m) = rhs(i,iv,m) - rhs(p,iv,q) * mat(i,iv,m,q,4)
	      rhs(i,iv,m) = rhs(i,iv,m) - rhs(r,iv,q) * mat(i,iv,m,q,5)
	    end do
	  end do
	  do iv = 1, nsys
	    rhs(i,iv,m) = rhs(i,iv,m) / mat(i,iv,m,m,3)
	  end do
	end do

	i = n - 3
	p = i + 1
	r = i + 2
	
	do m = nblk, 1, -1
	  do q = m+1, nblk
	    do iv = 1, nsys
	      rhs(i,iv,m) = rhs(i,iv,m) - rhs(i,iv,q) * mat(i,iv,m,q,3)
	    end do
	  end do
	  do q = 1, nblk
	    do iv = 1, nsys
	      rhs(i,iv,m) = rhs(i,iv,m) - rhs(p,iv,q) * mat(i,iv,m,q,4)
	      rhs(i,iv,m) = rhs(i,iv,m) - rhs(r,iv,q) * mat(i,iv,m,q,5)
	      rhs(i,iv,m) = rhs(i,iv,m) - rhs(n,iv,q) * per(i,iv,m,q,2)
	    end do
	  end do
	  do iv = 1, nsys
	    rhs(i,iv,m) = rhs(i,iv,m) / mat(i,iv,m,m,3)
	  end do
	end do

!.... now do the rest of the rows in reverse order

	do i = n - 4, 1, -1
	
	  p = i + 1
	  r = i + 2
	  
	  do m = nblk, 1, -1
	    do q = m+1, nblk
	      do iv = 1, nsys
		rhs(i,iv,m) = rhs(i,iv,m) - rhs(i,iv,q) * mat(i,iv,m,q,3)
	      end do
	    end do
	    do q = 1, nblk
	      do iv = 1, nsys
		rhs(i,iv,m) = rhs(i,iv,m) - rhs(p,iv,q) * mat(i,iv,m,q,4)
		rhs(i,iv,m) = rhs(i,iv,m) - rhs(r,iv,q) * mat(i,iv,m,q,5)
		rhs(i,iv,m) = rhs(i,iv,m) - rhs(n-1,iv,q) * per(i,iv,m,q,1)
		rhs(i,iv,m) = rhs(i,iv,m) - rhs(n,iv  ,q) * per(i,iv,m,q,2)
	      end do
	    end do
	    do iv = 1, nsys
	      rhs(i,iv,m) = rhs(i,iv,m) / mat(i,iv,m,m,3)
	    end do
	  end do
				  
	end do

!.... account for periodicity

	do q = 1, nblk
	  do iv = 1, nsys
	    rhs(np,iv,q) = rhs(1,iv,q)
	  end do
	end do
	
	end if
!=============================================================================!
	return
	end
