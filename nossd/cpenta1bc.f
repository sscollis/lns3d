!=============================================================================!
	subroutine cpenta1bc( nsys, n, nblk, mat, rhs, mult, fact, code )
!=============================================================================!
!  
!  Solves a block pentadiagonal system of equations without pivoting.
!  Vectorized over the first index.
!
!  This version includes the Boundary Treatment, and complex math
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
!
!  Revised:  8-6-96
!=============================================================================!
	implicit none

	integer n, nblk, nsys, code
	complex mat(nsys,n,nblk,nblk,5), rhs(nsys,n,nblk)
	complex mult(nsys), fact(nsys)
	
!.... useful constants

	real zero, one
	parameter (zero = 0.0, one = 1.0)

!.... local variables
	
	integer nc, i, j, l, m, p, q, r, s, t, iv

!=============================================================================!
!                       L U   F A C T O R I Z A T I O N 
!=============================================================================!

	if (code .eq. 0 .or. code .eq. 2) then

!.... do the first row

	i = 1
	
!.... do the first four columns

	  do nc = 1, nblk - 1
	  
	    do m = nc, nblk - 1
	      l = m + 1
	      do iv = 1, nsys
		mult(iv) = -mat(iv,i,l,nc,3) / mat(iv,i,nc,nc,3)
		mat(iv,i,l,nc,3) = mult(iv)
	      end do
	      do q = nc + 1, nblk
		do iv = 1, nsys
		  mat(iv,i,l,q,3) = mat(iv,i,l,q,3) + mult(iv) * mat(iv,i,nc,q,3)
		end do
	      end do
	      do q = 1, nblk
		do iv = 1, nsys
		  mat(iv,i,l,q,4) = mat(iv,i,l,q,4) + mult(iv) * mat(iv,i,nc,q,4)
		  mat(iv,i,l,q,5) = mat(iv,i,l,q,5) + mult(iv) * mat(iv,i,nc,q,5)
		  mat(iv,i,l,q,1) = mat(iv,i,l,q,1) + mult(iv) * mat(iv,i,nc,q,1)
		  mat(iv,i,l,q,2) = mat(iv,i,l,q,2) + mult(iv) * mat(iv,i,nc,q,2)
		end do
	      end do
	    end do	      

	  end do	! end of loop on nc

!.... for all columns
	  
	  do nc = 1, nblk

	    p = i + 1
	    r = i + 2
	    
	    do iv = 1, nsys
	      fact(iv) = one / mat(iv,i,nc,nc,3)
	    end do
	    
	    do m = 1, nblk
	      do iv = 1, nsys
		mult(iv) = -mat(iv,p,m,nc,2) * fact(iv)
		mat(iv,p,m,nc,2) = mult(iv)
	      end do
	      do q = nc + 1, nblk
		do iv = 1, nsys
		  mat(iv,p,m,q,2) = mat(iv,p,m,q,2) + mult(iv) * mat(iv,i,nc,q,3)
		end do
	      end do
	      do q = 1, nblk
		do iv = 1, nsys
		  mat(iv,p,m,q,3) = mat(iv,p,m,q,3) + mult(iv) * mat(iv,i,nc,q,4)
		  mat(iv,p,m,q,4) = mat(iv,p,m,q,4) + mult(iv) * mat(iv,i,nc,q,5)
		  mat(iv,p,m,q,5) = mat(iv,p,m,q,5) + mult(iv) * mat(iv,i,nc,q,1)
		  mat(iv,p,m,q,1) = mat(iv,p,m,q,1) + mult(iv) * mat(iv,i,nc,q,2)
		end do
	      end do
	      
	      do iv = 1, nsys
		mult(iv) = -mat(iv,r,m,nc,1) * fact(iv)
		mat(iv,r,m,nc,1) = mult(iv)
	      end do
	      do q = nc + 1, nblk
		do iv = 1, nsys
		  mat(iv,r,m,q,1) = mat(iv,r,m,q,1) + mult(iv) * mat(iv,i,nc,q,3)
		end do
	      end do
	      do q = 1, nblk
	        do iv = 1, nsys
		  mat(iv,r,m,q,2) = mat(iv,r,m,q,2) + mult(iv) * mat(iv,i,nc,q,4)
		  mat(iv,r,m,q,3) = mat(iv,r,m,q,3) + mult(iv) * mat(iv,i,nc,q,5)
		  mat(iv,r,m,q,4) = mat(iv,r,m,q,4) + mult(iv) * mat(iv,i,nc,q,1)
		  mat(iv,r,m,q,5) = mat(iv,r,m,q,5) + mult(iv) * mat(iv,i,nc,q,2)
		end do
	      end do
	      
	    end do

	  end do	! end of loop on nc

!.... do the second row

	i = 2
	
!.... do the first four columns

	  do nc = 1, nblk - 1
	  
	    do m = nc, nblk - 1
	      l = m + 1
	      do iv = 1, nsys
	        mult(iv) = -mat(iv,i,l,nc,3) / mat(iv,i,nc,nc,3)
		mat(iv,i,l,nc,3) = mult(iv)
	      end do
	      do q = nc + 1, nblk
		do iv = 1, nsys
		  mat(iv,i,l,q,3) = mat(iv,i,l,q,3) + mult(iv) * mat(iv,i,nc,q,3)
		end do
	      end do
	      do q = 1, nblk
		do iv = 1, nsys
		  mat(iv,i,l,q,4) = mat(iv,i,l,q,4) + mult(iv) * mat(iv,i,nc,q,4)
		  mat(iv,i,l,q,5) = mat(iv,i,l,q,5) + mult(iv) * mat(iv,i,nc,q,5)
		  mat(iv,i,l,q,1) = mat(iv,i,l,q,1) + mult(iv) * mat(iv,i,nc,q,1)
		end do
	      end do
	    end do	      

	  end do	! end of loop on nc

!.... for all columns
	  
	  do nc = 1, nblk

	    p = i + 1
	    r = i + 2
	    
	    do iv = 1, nsys
	      fact(iv) = one / mat(iv,i,nc,nc,3)
	    end do
	    
	    do m = 1, nblk
	      do iv = 1, nsys
		mult(iv) = -mat(iv,p,m,nc,2) * fact(iv)
		mat(iv,p,m,nc,2) = mult(iv)
	      end do
	      do q = nc + 1, nblk
		do iv = 1, nsys
		  mat(iv,p,m,q,2) = mat(iv,p,m,q,2) + mult(iv) * mat(iv,i,nc,q,3)
		end do
	      end do
	      do q = 1, nblk
		do iv = 1, nsys
		  mat(iv,p,m,q,3) = mat(iv,p,m,q,3) + mult(iv) * mat(iv,i,nc,q,4)
		  mat(iv,p,m,q,4) = mat(iv,p,m,q,4) + mult(iv) * mat(iv,i,nc,q,5)
		  mat(iv,p,m,q,5) = mat(iv,p,m,q,5) + mult(iv) * mat(iv,i,nc,q,1)
		end do
	      end do
	      
	      do iv = 1, nsys
		mult(iv) = -mat(iv,r,m,nc,1) * fact(iv)
		mat(iv,r,m,nc,1) = mult(iv)
	      end do
	      do q = nc + 1, nblk
		do iv = 1, nsys
		  mat(iv,r,m,q,1) = mat(iv,r,m,q,1) + mult(iv) * mat(iv,i,nc,q,3)
		end do
	      end do
	      do q = 1, nblk
		do iv = 1, nsys
		  mat(iv,r,m,q,2) = mat(iv,r,m,q,2) + mult(iv) * mat(iv,i,nc,q,4)
		  mat(iv,r,m,q,3) = mat(iv,r,m,q,3) + mult(iv) * mat(iv,i,nc,q,5)
		  mat(iv,r,m,q,4) = mat(iv,r,m,q,4) + mult(iv) * mat(iv,i,nc,q,1)
		end do
	      end do
	      
	    end do

	  end do	! end of loop on nc

!.... do the interior rows

 	do i = 3, n - 5

!.... do the first four columns

	  do nc = 1, nblk - 1
	  
	    do m = nc, nblk - 1
	      l = m + 1
	      do iv = 1, nsys
		mult(iv) = -mat(iv,i,l,nc,3) / mat(iv,i,nc,nc,3)
		mat(iv,i,l,nc,3) = mult(iv)
	      end do
	      do q = nc + 1, nblk
		do iv = 1, nsys
		  mat(iv,i,l,q,3) = mat(iv,i,l,q,3) + mult(iv) * mat(iv,i,nc,q,3)
		end do
	      end do
	      do q = 1, nblk
		do iv = 1, nsys
		  mat(iv,i,l,q,4) = mat(iv,i,l,q,4) + mult(iv) * mat(iv,i,nc,q,4)
		  mat(iv,i,l,q,5) = mat(iv,i,l,q,5) + mult(iv) * mat(iv,i,nc,q,5)
		end do
	      end do
	    end do	      

	  end do	! end of loop on nc

!.... for all columns
	  
	  do nc = 1, nblk

	    p = i + 1
	    r = i + 2
	    
	    do iv = 1, nsys
	      fact(iv) = one / mat(iv,i,nc,nc,3)
	    end do
	    
	    do m = 1, nblk
	      do iv = 1, nsys
		mult(iv) = -mat(iv,p,m,nc,2) * fact(iv)
		mat(iv,p,m,nc,2) = mult(iv)
	      end do
	      do q = nc + 1, nblk
		do iv = 1, nsys
		  mat(iv,p,m,q,2) = mat(iv,p,m,q,2) + mult(iv) * mat(iv,i,nc,q,3)
		end do
	      end do
	      do q = 1, nblk
		do iv = 1, nsys
		  mat(iv,p,m,q,3) = mat(iv,p,m,q,3) + mult(iv) * mat(iv,i,nc,q,4)
		  mat(iv,p,m,q,4) = mat(iv,p,m,q,4) + mult(iv) * mat(iv,i,nc,q,5)
		end do
	      end do
	      
	      do iv = 1, nsys
		mult(iv) = -mat(iv,r,m,nc,1) * fact(iv)
		mat(iv,r,m,nc,1) = mult(iv)
	      end do
	      do q = nc + 1, nblk
		do iv = 1, nsys
		  mat(iv,r,m,q,1) = mat(iv,r,m,q,1) + mult(iv) * mat(iv,i,nc,q,3)
		end do
	      end do
	      do q = 1, nblk
		do iv = 1, nsys
		  mat(iv,r,m,q,2) = mat(iv,r,m,q,2) + mult(iv) * mat(iv,i,nc,q,4)
		  mat(iv,r,m,q,3) = mat(iv,r,m,q,3) + mult(iv) * mat(iv,i,nc,q,5)
		end do
	      end do
	      
	    end do

	  end do	! end of loop on nc
	  
	end do		! end of loop on i

!.... do the (n-4)th row

	i = n - 4

!.... do the first four columns

	  do nc = 1, nblk - 1
	  
	    do m = nc, nblk - 1
	      l = m + 1
	      do iv = 1, nsys
		mult(iv) = -mat(iv,i,l,nc,3) / mat(iv,i,nc,nc,3)
		mat(iv,i,l,nc,3) = mult(iv)
	      end do
	      do q = nc + 1, nblk
		do iv = 1, nsys
		  mat(iv,i,l,q,3) = mat(iv,i,l,q,3) + mult(iv) * mat(iv,i,nc,q,3)
		end do
	      end do
	      do q = 1, nblk
		do iv = 1, nsys
		  mat(iv,i,l,q,4) = mat(iv,i,l,q,4) + mult(iv) * mat(iv,i,nc,q,4)
		  mat(iv,i,l,q,5) = mat(iv,i,l,q,5) + mult(iv) * mat(iv,i,nc,q,5)
		end do
	      end do
	    end do	      

	  end do	! end of loop on nc

!.... for all columns
	  
	  do nc = 1, nblk

	    p = i + 1
	    r = i + 2
	    s = i + 3
	    t = i + 4
	    
	    do iv = 1, nsys
	      fact(iv) = one / mat(iv,i,nc,nc,3)
	    end do
	    
	    do m = 1, nblk
	      do iv = 1, nsys
		mult(iv) = -mat(iv,p,m,nc,2) * fact(iv)
		mat(iv,p,m,nc,2) = mult(iv)
	      end do
	      do q = nc + 1, nblk
		do iv = 1, nsys
		  mat(iv,p,m,q,2) = mat(iv,p,m,q,2) + mult(iv) * mat(iv,i,nc,q,3)
		end do
	      end do
	      do q = 1, nblk
		do iv = 1, nsys
		  mat(iv,p,m,q,3) = mat(iv,p,m,q,3) + mult(iv) * mat(iv,i,nc,q,4)
		  mat(iv,p,m,q,4) = mat(iv,p,m,q,4) + mult(iv) * mat(iv,i,nc,q,5)
		end do
	      end do
	      
	      do iv = 1, nsys
		mult(iv) = -mat(iv,r,m,nc,1) * fact(iv)
		mat(iv,r,m,nc,1) = mult(iv)
	      end do
	      do q = nc + 1, nblk
		do iv = 1, nsys
		  mat(iv,r,m,q,1) = mat(iv,r,m,q,1) + mult(iv) * mat(iv,i,nc,q,3)
		end do
	      end do
	      do q = 1, nblk
		do iv = 1, nsys
		  mat(iv,r,m,q,2) = mat(iv,r,m,q,2) + mult(iv) * mat(iv,i,nc,q,4)
		  mat(iv,r,m,q,3) = mat(iv,r,m,q,3) + mult(iv) * mat(iv,i,nc,q,5)
		end do
	      end do
	      
	      do iv = 1, nsys
		mult(iv) = -mat(iv,s,m,nc,5) * fact(iv)
		mat(iv,s,m,nc,5) = mult(iv)
	      end do
	      do q = nc + 1, nblk
		do iv = 1, nsys
		  mat(iv,s,m,q,5) = mat(iv,s,m,q,5) + mult(iv) * mat(iv,i,nc,q,3)
		end do
	      end do
	      do q = 1, nblk
		do iv = 1, nsys
		  mat(iv,s,m,q,1) = mat(iv,s,m,q,1) + mult(iv) * mat(iv,i,nc,q,4)
		  mat(iv,s,m,q,2) = mat(iv,s,m,q,2) + mult(iv) * mat(iv,i,nc,q,5)
		end do
	      end do

	      do iv = 1, nsys
		mult(iv) = -mat(iv,t,m,nc,4) * fact(iv)
		mat(iv,t,m,nc,4) = mult(iv)
	      end do
	      do q = nc + 1, nblk
		do iv = 1, nsys
		  mat(iv,t,m,q,4) = mat(iv,t,m,q,4) + mult(iv) * mat(iv,i,nc,q,3)
		end do
	      end do
	      do q = 1, nblk
		do iv = 1, nsys
		  mat(iv,t,m,q,5) = mat(iv,t,m,q,5) + mult(iv) * mat(iv,i,nc,q,4)
		  mat(iv,t,m,q,1) = mat(iv,t,m,q,1) + mult(iv) * mat(iv,i,nc,q,5)
		end do
	      end do

	    end do

	  end do	! end of loop on nc

!.... do the (n-3)rd row

	i = n - 3

!.... do the first four columns

	  do nc = 1, nblk - 1
	  
	    do m = nc, nblk - 1
	      l = m + 1
	      do iv = 1, nsys
		mult(iv) = -mat(iv,i,l,nc,3) / mat(iv,i,nc,nc,3)
		mat(iv,i,l,nc,3) = mult(iv)
	      end do
	      do q = nc + 1, nblk
		do iv = 1, nsys
		  mat(iv,i,l,q,3) = mat(iv,i,l,q,3) + mult(iv) * mat(iv,i,nc,q,3)
		end do
	      end do
	      do q = 1, nblk
		do iv = 1, nsys
		  mat(iv,i,l,q,4) = mat(iv,i,l,q,4) + mult(iv) * mat(iv,i,nc,q,4)
		  mat(iv,i,l,q,5) = mat(iv,i,l,q,5) + mult(iv) * mat(iv,i,nc,q,5)
		end do
	      end do
	    end do	      

	  end do	! end of loop on nc

!.... for all columns
	  
	  do nc = 1, nblk

	    p = i + 1
	    r = i + 2
	    s = i + 3
	    
	    do iv = 1, nsys
	      fact(iv) = one / mat(iv,i,nc,nc,3)
	    end do
	    
	    do m = 1, nblk
	      do iv = 1, nsys
		mult(iv) = -mat(iv,p,m,nc,2) * fact(iv)
		mat(iv,p,m,nc,2) = mult(iv)
	      end do
	      do q = nc + 1, nblk
		do iv = 1, nsys
		  mat(iv,p,m,q,2) = mat(iv,p,m,q,2) + mult(iv) * mat(iv,i,nc,q,3)
		end do
	      end do
	      do q = 1, nblk
		do iv = 1, nsys
		  mat(iv,p,m,q,3) = mat(iv,p,m,q,3) + mult(iv) * mat(iv,i,nc,q,4)
		  mat(iv,p,m,q,4) = mat(iv,p,m,q,4) + mult(iv) * mat(iv,i,nc,q,5)
		end do
	      end do
	      
	      do iv = 1, nsys
		mult(iv) = -mat(iv,r,m,nc,1) * fact(iv)
		mat(iv,r,m,nc,1) = mult(iv)
	      end do
	      do q = nc + 1, nblk
		do iv = 1, nsys
		  mat(iv,r,m,q,1) = mat(iv,r,m,q,1) + mult(iv) * mat(iv,i,nc,q,3)
		end do
	      end do
	      do q = 1, nblk
		do iv = 1, nsys
		  mat(iv,r,m,q,2) = mat(iv,r,m,q,2) + mult(iv) * mat(iv,i,nc,q,4)
		  mat(iv,r,m,q,3) = mat(iv,r,m,q,3) + mult(iv) * mat(iv,i,nc,q,5)
		end do
	      end do
	      
	      do iv = 1, nsys
		mult(iv) = -mat(iv,s,m,nc,5) * fact(iv)
		mat(iv,s,m,nc,5) = mult(iv)
	      end do
	      do q = nc + 1, nblk
		do iv = 1, nsys
		  mat(iv,s,m,q,5) = mat(iv,s,m,q,5) + mult(iv) * mat(iv,i,nc,q,3)
		end do
	      end do
	      do q = 1, nblk
		do iv = 1, nsys
		  mat(iv,s,m,q,1) = mat(iv,s,m,q,1) + mult(iv) * mat(iv,i,nc,q,4)
		  mat(iv,s,m,q,2) = mat(iv,s,m,q,2) + mult(iv) * mat(iv,i,nc,q,5)
		end do
	      end do

	    end do

	  end do	! end of loop on nc

!.... do the (n-2)nd row

	i = n - 2

!.... do the first four columns

	  do nc = 1, nblk - 1
	  
	    do m = nc, nblk - 1
	      l = m + 1
	      do iv = 1, nsys
		mult(iv) = -mat(iv,i,l,nc,3) / mat(iv,i,nc,nc,3)
		mat(iv,i,l,nc,3) = mult(iv)
	      end do
	      do q = nc + 1, nblk
		do iv = 1, nsys
		  mat(iv,i,l,q,3) = mat(iv,i,l,q,3) + mult(iv) * mat(iv,i,nc,q,3)
		end do
	      end do
	      do q = 1, nblk
		do iv = 1, nsys
		  mat(iv,i,l,q,4) = mat(iv,i,l,q,4) + mult(iv) * mat(iv,i,nc,q,4)
		  mat(iv,i,l,q,5) = mat(iv,i,l,q,5) + mult(iv) * mat(iv,i,nc,q,5)
		end do
	      end do
	    end do	      

	  end do	! end of loop on nc

!.... for all columns
	  
	  do nc = 1, nblk

	    p = i + 1
	    r = i + 2
	    
	    do iv = 1, nsys
	      fact(iv) = one / mat(iv,i,nc,nc,3)
	    end do
	    
	    do m = 1, nblk
	      do iv = 1, nsys
		mult(iv) = -mat(iv,p,m,nc,2) * fact(iv)
		mat(iv,p,m,nc,2) = mult(iv)
	      end do
	      do q = nc + 1, nblk
		do iv = 1, nsys
		  mat(iv,p,m,q,2) = mat(iv,p,m,q,2) + mult(iv) * mat(iv,i,nc,q,3)
		end do
	      end do
	      do q = 1, nblk
		do iv = 1, nsys
		  mat(iv,p,m,q,3) = mat(iv,p,m,q,3) + mult(iv) * mat(iv,i,nc,q,4)
		  mat(iv,p,m,q,4) = mat(iv,p,m,q,4) + mult(iv) * mat(iv,i,nc,q,5)
		end do
	      end do
	      
	      do iv = 1, nsys
		mult(iv) = -mat(iv,r,m,nc,1) * fact(iv)
		mat(iv,r,m,nc,1) = mult(iv)
	      end do
	      do q = nc + 1, nblk
		do iv = 1, nsys
		  mat(iv,r,m,q,1) = mat(iv,r,m,q,1) + mult(iv) * mat(iv,i,nc,q,3)
		end do
	      end do
	      do q = 1, nblk
		do iv = 1, nsys
		  mat(iv,r,m,q,2) = mat(iv,r,m,q,2) + mult(iv) * mat(iv,i,nc,q,4)
		  mat(iv,r,m,q,3) = mat(iv,r,m,q,3) + mult(iv) * mat(iv,i,nc,q,5)
		end do
	      end do
	      
	    end do

	  end do	! end of loop on nc

!.... do the next-to-the-last row
	
	i = n - 1

!.... do the first four columns

	  do nc = 1, nblk - 1
	  
	    do m = nc, nblk - 1
	      l = m + 1
	      do iv = 1, nsys
		mult(iv) = -mat(iv,i,l,nc,3) / mat(iv,i,nc,nc,3)
		mat(iv,i,l,nc,3) = mult(iv)
	      end do
	      do q = nc + 1, nblk
		do iv = 1, nsys
		  mat(iv,i,l,q,3) = mat(iv,i,l,q,3) + mult(iv) * mat(iv,i,nc,q,3)
		end do
	      end do
	      do q = 1, nblk
		do iv = 1, nsys
		  mat(iv,i,l,q,4) = mat(iv,i,l,q,4) + mult(iv) * mat(iv,i,nc,q,4)
		end do
	      end do
	    end do	      

	  end do	! end of loop on nc

!.... for all columns
	  
	  do nc = 1, nblk

	    p = i + 1
	    
	    do iv = 1, nsys
	      fact(iv) = one / mat(iv,i,nc,nc,3)
	    end do
	    
	    do m = 1, nblk
	      do iv = 1, nsys
		mult(iv) = -mat(iv,p,m,nc,2) * fact(iv)
		mat(iv,p,m,nc,2) = mult(iv)
	      end do
	      do q = nc + 1, nblk
		do iv = 1, nsys
		  mat(iv,p,m,q,2) = mat(iv,p,m,q,2) + mult(iv) * mat(iv,i,nc,q,3)
		end do
	      end do
	      do q = 1, nblk
		do iv = 1, nsys
		  mat(iv,p,m,q,3) = mat(iv,p,m,q,3) + mult(iv) * mat(iv,i,nc,q,4)
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
	      mult(iv) = -mat(iv,i,l,nc,3) / mat(iv,i,nc,nc,3)
	      mat(iv,i,l,nc,3) = mult(iv)
	    end do
	    do q = nc + 1, nblk
	      do iv = 1, nsys
		mat(iv,i,l,q,3) = mat(iv,i,l,q,3) + mult(iv) * mat(iv,i,nc,q,3)
	      end do
	    end do
	  end do	      

	end do	! end of loop on nc

	end if
!=============================================================================!
!                     F O R W A R D   E L I M I N A T I O N 
!=============================================================================!
	if (code.eq.1 .or. code.eq.2) then
	
!.... do the first row

	i = 1
	
!.... do the first four columns

	  do nc = 1, nblk - 1
	  
	    do m = nc, nblk - 1
	      l = m + 1
	      do iv = 1, nsys
		mult(iv) = mat(iv,i,l,nc,3)
		rhs(iv,i,l) = rhs(iv,i,l) + mult(iv) * rhs(iv,i,nc)
	      end do
	    end do	      

	  end do	! end of loop on nc

!.... for all columns
	  
	  do nc = 1, nblk

	    p = i + 1
	    r = i + 2
	    
	    do m = 1, nblk
	      do iv = 1, nsys
		mult(iv) = mat(iv,p,m,nc,2)
		rhs(iv,p,m) = rhs(iv,p,m) + mult(iv) * rhs(iv,i,nc)
	      end do
	      
	      do iv = 1, nsys
		mult(iv) = mat(iv,r,m,nc,1)
		rhs(iv,r,m) = rhs(iv,r,m) + mult(iv) * rhs(iv,i,nc)
	      end do
	      
	    end do

	  end do	! end of loop on nc

!.... do the second row

	i = 2
	
!.... do the first four columns

	  do nc = 1, nblk - 1
	  
	    do m = nc, nblk - 1
	      l = m + 1
	      do iv = 1, nsys
	        mult(iv) = mat(iv,i,l,nc,3)
		rhs(iv,i,l) = rhs(iv,i,l) + mult(iv) * rhs(iv,i,nc)
	      end do
	    end do	      

	  end do	! end of loop on nc

!.... for all columns
	  
	  do nc = 1, nblk

	    p = i + 1
	    r = i + 2
	    
	    do m = 1, nblk
	      do iv = 1, nsys
		mult(iv) = mat(iv,p,m,nc,2)
		rhs(iv,p,m) = rhs(iv,p,m) + mult(iv) * rhs(iv,i,nc)
	      end do
	      
	      do iv = 1, nsys
		mult(iv) = mat(iv,r,m,nc,1)
		rhs(iv,r,m) = rhs(iv,r,m) + mult(iv) * rhs(iv,i,nc)
	      end do
	      
	    end do

	  end do	! end of loop on nc

!.... do the interior rows

 	do i = 3, n - 5

!.... do the first four columns

	  do nc = 1, nblk - 1
	  
	    do m = nc, nblk - 1
	      l = m + 1
	      do iv = 1, nsys
		mult(iv) = mat(iv,i,l,nc,3)
		rhs(iv,i,l) = rhs(iv,i,l) + mult(iv) * rhs(iv,i,nc)
	      end do
	    end do	      

	  end do	! end of loop on nc

!.... for all columns
	  
	  do nc = 1, nblk

	    p = i + 1
	    r = i + 2
	    
	    do m = 1, nblk
	      do iv = 1, nsys
		mult(iv) = mat(iv,p,m,nc,2)
		rhs(iv,p,m) = rhs(iv,p,m) + mult(iv) * rhs(iv,i,nc)
	      end do
	      
	      do iv = 1, nsys
		mult(iv) = mat(iv,r,m,nc,1)
		rhs(iv,r,m) = rhs(iv,r,m) + mult(iv) * rhs(iv,i,nc)
	      end do
	      
	    end do

	  end do	! end of loop on nc
	  
	end do		! end of loop on i

!.... do the (n-4)th row

	i = n - 4

!.... do the first four columns

	  do nc = 1, nblk - 1
	  
	    do m = nc, nblk - 1
	      l = m + 1
	      do iv = 1, nsys
		mult(iv) = mat(iv,i,l,nc,3)
		rhs(iv,i,l) = rhs(iv,i,l) + mult(iv) * rhs(iv,i,nc)
	      end do
	    end do	      

	  end do	! end of loop on nc

!.... for all columns
	  
	  do nc = 1, nblk

	    p = i + 1
	    r = i + 2
	    s = i + 3
	    t = i + 4
	    
	    do m = 1, nblk
	      do iv = 1, nsys
		mult(iv) = mat(iv,p,m,nc,2)
		rhs(iv,p,m) = rhs(iv,p,m) + mult(iv) * rhs(iv,i,nc)
	      end do
	      
	      do iv = 1, nsys
		mult(iv) = mat(iv,r,m,nc,1)
		rhs(iv,r,m) = rhs(iv,r,m) + mult(iv) * rhs(iv,i,nc)
	      end do
	      
	      do iv = 1, nsys
		mult(iv) = mat(iv,s,m,nc,5)
		rhs(iv,s,m) = rhs(iv,s,m) + mult(iv) * rhs(iv,i,nc)
	      end do

	      do iv = 1, nsys
		mult(iv) = mat(iv,t,m,nc,4)
		rhs(iv,t,m) = rhs(iv,t,m) + mult(iv) * rhs(iv,i,nc)
	      end do

	    end do

	  end do	! end of loop on nc

!.... do the (n-3)rd row

	i = n - 3

!.... do the first four columns

	  do nc = 1, nblk - 1
	  
	    do m = nc, nblk - 1
	      l = m + 1
	      do iv = 1, nsys
		mult(iv) = mat(iv,i,l,nc,3)
		rhs(iv,i,l) = rhs(iv,i,l) + mult(iv) * rhs(iv,i,nc)
	      end do
	    end do	      

	  end do	! end of loop on nc

!.... for all columns
	  
	  do nc = 1, nblk

	    p = i + 1
	    r = i + 2
	    s = i + 3
	    
	    do m = 1, nblk
	      do iv = 1, nsys
		mult(iv) = mat(iv,p,m,nc,2)
		rhs(iv,p,m) = rhs(iv,p,m) + mult(iv) * rhs(iv,i,nc)
	      end do
	      
	      do iv = 1, nsys
		mult(iv) = mat(iv,r,m,nc,1)
		rhs(iv,r,m) = rhs(iv,r,m) + mult(iv) * rhs(iv,i,nc)
	      end do
	      
	      do iv = 1, nsys
		mult(iv) = mat(iv,s,m,nc,5)
	        rhs(iv,s,m) = rhs(iv,s,m) + mult(iv) * rhs(iv,i,nc)
	      end do

	    end do

	  end do	! end of loop on nc

!.... do the (n-2)nd row

	i = n - 2

!.... do the first four columns

	  do nc = 1, nblk - 1
	  
	    do m = nc, nblk - 1
	      l = m + 1
	      do iv = 1, nsys
		mult(iv) = mat(iv,i,l,nc,3)
		rhs(iv,i,l) = rhs(iv,i,l) + mult(iv) * rhs(iv,i,nc)
	      end do
	    end do	      

	  end do	! end of loop on nc

!.... for all columns
	  
	  do nc = 1, nblk

	    p = i + 1
	    r = i + 2
	    
	    do m = 1, nblk
	      do iv = 1, nsys
		mult(iv) = mat(iv,p,m,nc,2)
		rhs(iv,p,m) = rhs(iv,p,m) + mult(iv) * rhs(iv,i,nc)
	      end do
	      
	      do iv = 1, nsys
		mult(iv) = mat(iv,r,m,nc,1)
		rhs(iv,r,m) = rhs(iv,r,m) + mult(iv) * rhs(iv,i,nc)
	      end do
	      
	    end do

	  end do	! end of loop on nc

!.... do the next-to-the-last row
	
	i = n - 1

!.... do the first four columns

	  do nc = 1, nblk - 1
	  
	    do m = nc, nblk - 1
	      l = m + 1
	      do iv = 1, nsys
		mult(iv) = mat(iv,i,l,nc,3)
		rhs(iv,i,l) = rhs(iv,i,l) + mult(iv) * rhs(iv,i,nc)
	      end do
	    end do	      

	  end do	! end of loop on nc

!.... for all columns
	  
	  do nc = 1, nblk

	    p = i + 1
	    
	    do m = 1, nblk
	      do iv = 1, nsys
		mult(iv) = mat(iv,p,m,nc,2)
		rhs(iv,p,m) = rhs(iv,p,m) + mult(iv) * rhs(iv,i,nc)
	      end do
	    end do

	  end do	! end of loop on nc

!.... do the last row

	i = n
	
	do nc = 1, nblk - 1
	
	  do m = nc, nblk - 1
	    l = m + 1
	    do iv = 1, nsys
	      mult(iv) = mat(iv,i,l,nc,3)
	      rhs(iv,i,l) = rhs(iv,i,l) + mult(iv) * rhs(iv,i,nc)
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
	      rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,i,q) * mat(iv,i,m,q,3)
	    end do
	  end do
	  do iv = 1, nsys
	    rhs(iv,i,m) = rhs(iv,i,m) / mat(iv,i,m,m,3)
	  end do
	end do

!.... do the next-to-the-last row

	i = n - 1
	p = i + 1
	
	do m = nblk, 1, -1
	  do q = m+1, nblk
	    do iv = 1, nsys
	      rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,i,q) * mat(iv,i,m,q,3)
	    end do
	  end do
	  do q = 1, nblk
	    do iv = 1, nsys
	      rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,p,q) * mat(iv,i,m,q,4)
	    end do
	  end do
	  do iv = 1, nsys
	    rhs(iv,i,m) = rhs(iv,i,m) / mat(iv,i,m,m,3)
	  end do
	end do

!.... now do the rest of the rows in reverse order

	do i = n - 2, 3, -1
	
	  p = i + 1
	  r = i + 2
	  
	  do m = nblk, 1, -1
	    do q = m+1, nblk
	      do iv = 1, nsys
		rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,i,q) * mat(iv,i,m,q,3)
	      end do
	    end do
	    do q = 1, nblk
	      do iv = 1, nsys
		rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,p,q) * mat(iv,i,m,q,4)
		rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,r,q) * mat(iv,i,m,q,5)
	      end do
	    end do
	    do iv = 1, nsys
	      rhs(iv,i,m) = rhs(iv,i,m) / mat(iv,i,m,m,3)
	    end do
	  end do
				  
	end do
	
!.... Do the last two boundary rows

	  i = 2
	  p = i + 1
	  r = i + 2
	  s = i + 3
	  
	  do m = nblk, 1, -1
	    do q = m+1, nblk
	      do iv = 1, nsys
		rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,i,q) * mat(iv,i,m,q,3)
	      end do
	    end do
	    do q = 1, nblk
	      do iv = 1, nsys
		rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,p,q) * mat(iv,i,m,q,4)
		rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,r,q) * mat(iv,i,m,q,5)
		rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,s,q) * mat(iv,i,m,q,1)
	      end do
	    end do
	    do iv = 1, nsys
	      rhs(iv,i,m) = rhs(iv,i,m) / mat(iv,i,m,m,3)
	    end do
	  end do

	  i = 1
	  p = i + 1
	  r = i + 2
	  s = i + 3
	  t = i + 4
	  
	  do m = nblk, 1, -1
	    do q = m+1, nblk
	      do iv = 1, nsys
		rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,i,q) * mat(iv,i,m,q,3)
	      end do
	    end do
	    do q = 1, nblk
	      do iv = 1, nsys
		rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,p,q) * mat(iv,i,m,q,4)
		rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,r,q) * mat(iv,i,m,q,5)
		rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,s,q) * mat(iv,i,m,q,1)
		rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,t,q) * mat(iv,i,m,q,2)
	      end do
	    end do
	    do iv = 1, nsys
	      rhs(iv,i,m) = rhs(iv,i,m) / mat(iv,i,m,m,3)
	    end do
	  end do
	
	end if
!=============================================================================!
	return
	end

