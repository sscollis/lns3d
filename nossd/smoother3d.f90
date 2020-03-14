!=============================================================================!
	subroutine smoother3D(rl, vl) 
!  
!       Fourth order explicit smoother for the NS equations.
!       From: Computational Fluid Mechanics and Heat Transfer
!             D.A. Anderson, J. C. Tannehill, Rl H. Pletcher
!             Page:  450, 495 
!
!	Written: 6-9-95
!
!       Revised: 7-9-96
!
!	Notes:   Need to account for periodicity smoothed fields.
!=============================================================================!
	use global
	use buff_mod
	use stencil
	implicit none

	complex :: rl(ny,nx,ndof), vl(ny,nx,ndof)
	real :: isign
	integer :: idof, i, j

!=============================================================================!

!.... fourth-order dissipation with reduced order near the boundaries

	if (eps_e .eq. zero) return

!.... \eta direction
	
	do idof = 1, ndof
	  do i = 1, nx

	    j = 2
	    rl(j,i,idof) = rl(j,i,idof) + eps_e * buff(j,i) * ( &
					   -1.0 * vl(j-1,i,idof) + &
					    2.0 * vl(j,i  ,idof) - &
					    1.0 * vl(j+1,i,idof) )
	    
	    j = 3
	    rl(j,i,idof) = rl(j,i,idof) + eps_e * buff(j,i) * ( &
					  vl(j-2,i,idof) - &
					  4.0 * vl(j-1,i,idof) + &
					  6.0 * vl(j  ,i,idof) - &
					  4.0 * vl(j+1,i,idof) + &
					  1.0 * vl(j+2,i,idof) )

	    do j = 4, ny-3
	      rl(j,i,idof) = rl(j,i,idof) + eps_e * buff(j,i) * ( &
					    fa1 * vl(j-3,i,idof) + &
					    fa2 * vl(j-2,i,idof) + &
					    fa3 * vl(j-1,i,idof) + &
					    fa4 * vl(j  ,i,idof) + &
					    fa5 * vl(j+1,i,idof) + &
					    fa6 * vl(j+2,i,idof) + &
					    fa7 * vl(j+3,i,idof) )
	    end do

	    j = ny-2
	    rl(j,i,idof) = rl(j,i,idof) + eps_e * buff(j,i) * ( &
					  vl(j-2,i,idof) - &
					  4.0 * vl(j-1,i,idof) + &
					  6.0 * vl(j  ,i,idof) - &
					  4.0 * vl(j+1,i,idof) + &
					  1.0 * vl(j+2,i,idof) )

	    j = ny-1
	    rl(j,i,idof) = rl(j,i,idof) + eps_e * buff(j,i) * ( &
					  -1.0 * vl(j-1,i,idof) + &
					    2.0 * vl(j  ,i,idof) - &
					    1.0 * vl(j+1,i,idof) )

	  end do
	end do
					  
!.... \xi direction

	if (lsym) then
	  do idof = 1, ndof
	    if (idof .eq. 3) then
	      isign = -one
	    else
	      isign = one
	    end if
	    do j = 1, ny
	      i = 1
	      rl(j,i,idof) = rl(j,i,idof) + eps_e * buff(j,i) * ( &
					    isign * fa1 * vl(j,4,idof) + &
					    isign * fa2 * vl(j,3,idof) + &
					    isign * fa3 * vl(j,2,idof) + &
					    fa4 * vl(j,1,idof) + &
					    fa5 * vl(j,2,idof) + &
					    fa6 * vl(j,3,idof) + &
					    fa7 * vl(j,4,idof) )
	      
	      i = 2
	      rl(j,i,idof) = rl(j,i,idof) + eps_e * buff(j,i) * ( &
					    isign * fa1 * vl(j,3,idof) + &
					    isign * fa2 * vl(j,2,idof) + &
					    fa3 * vl(j,1,idof) + &
					    fa4 * vl(j,2,idof) + &
					    fa5 * vl(j,3,idof) + &
					    fa6 * vl(j,4,idof) + &
					    fa7 * vl(j,5,idof) )

	      i = 3
	      rl(j,i,idof) = rl(j,i,idof) + eps_e * buff(j,i) * ( &
					    isign * fa1 * vl(j,2,idof) + &
					    fa2 * vl(j,1,idof) + &
					    fa3 * vl(j,2,idof) + &
					    fa4 * vl(j,3,idof) + &
					    fa5 * vl(j,4,idof) + &
					    fa6 * vl(j,5,idof) + &
					    fa7 * vl(j,6,idof) )
	    end do
	  end do
	else
	  do idof = 1, ndof
	    do j = 1, ny
	      i = 2
	      rl(j,i,idof) = rl(j,i,idof) + eps_e * buff(j,i) * ( &
					    -1.0 * vl(j,i-1,idof) + &
					    2.0 * vl(j,i  ,idof) - &
					    1.0 * vl(j,i+1,idof) )
	      
	      i = 3
	      rl(j,i,idof) = rl(j,i,idof) + eps_e * buff(j,i) * ( &
					    vl(j,i-2,idof) - &
					    4.0 * vl(j,i-1,idof) + &
					    6.0 * vl(j,i  ,idof) - &
					    4.0 * vl(j,i+1,idof) + &
					    1.0 * vl(j,i+2,idof) )
	    end do
	  end do
	end if

!.... Interior

	do idof = 1, ndof
	  do j = 1, ny
	    do i = 4, nx-3
	      rl(j,i,idof) = rl(j,i,idof) + &
				  eps_e * buff(j,i) * ( &
				  fa1 * vl(j,i-3,idof) + &
				  fa2 * vl(j,i-2,idof) + &
				  fa3 * vl(j,i-1,idof) + &
				  fa4 * vl(j,i  ,idof) + &
				  fa5 * vl(j,i+1,idof) + &
				  fa6 * vl(j,i+2,idof) + &
				  fa7 * vl(j,i+3,idof) )
	    end do
	  end do
	end do
	
	if (rsym) then
	  do idof = 1, ndof
	    if (idof .eq. 3) then
	      isign = -one
	    else
	      isign = one
	    end if
	    do j = 1, ny
	      i = nx-2
	      rl(j,i,idof) = rl(j,i,idof) + eps_e * buff(j,i) * ( &
					    fa1 * vl(j,nx-5,idof) + &
					    fa2 * vl(j,nx-4,idof) + &
					    fa3 * vl(j,nx-3,idof) + &
					    fa4 * vl(j,nx-2,idof) + &
					    fa5 * vl(j,nx-1,idof) + &
					    fa6 * vl(j,nx,idof) + &
					    isign * fa7 * vl(j,nx-1,idof) )
	      
	      i = nx-1
	      rl(j,i,idof) = rl(j,i,idof) + eps_e * buff(j,i) * ( &
					    fa1 * vl(j,nx-4,idof) + &
					    fa2 * vl(j,nx-3,idof) + &
					    fa3 * vl(j,nx-2,idof) + &
					    fa4 * vl(j,nx-1,idof) + &
					    fa5 * vl(j,nx,idof) + &
					    isign * fa6 * vl(j,nx-1,idof) + &
					    isign * fa7 * vl(j,nx-2,idof) )

	      i = nx
	      rl(j,i,idof) = rl(j,i,idof) + eps_e * buff(j,i) * ( &
					    fa1 * vl(j,nx-3,idof) + &
					    fa2 * vl(j,nx-2,idof) + &
					    fa3 * vl(j,nx-1,idof) + &
					    fa4 * vl(j,nx,idof) + &
					    isign * fa5 * vl(j,nx-1,idof) + &
					    isign * fa6 * vl(j,nx-2,idof) + &
					    isign * fa7 * vl(j,nx-3,idof) )
	    end do
	  end do
	else
	  i = nx-2
	  rl(j,i,idof) = rl(j,i,idof) + eps_e * buff(j,i) * ( &
					vl(j,i-2,idof) - &
					4.0 * vl(j,i-1,idof) + &
					6.0 * vl(j,i  ,idof) - &
					4.0 * vl(j,i+1,idof) + &
					1.0 * vl(j,i+2,idof) )
  
	  i = nx-1
	  rl(j,i,idof) = rl(j,i,idof) + eps_e * buff(j,i) * ( &
					-1.0 * vl(j,i-1,idof) + &
					2.0 * vl(j,i,idof) - &
					1.0 * vl(j,i+1,idof) )
	end if
	
	return
	end
