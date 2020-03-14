!=======================================================================================================!
	subroutine lhsbt2f3d( mat, Bh, Dh, Dhi, Vh, Bhi, spgl, spg2l, dtl, calcd )
!  
!  Correct the LHS for boundary treatment in the \eta direction
!
!  This version supports the 4th order LHS
!
!  Revised: 4-23-96
!=======================================================================================================!
	use global
	use stencil
	use buff_mod
	implicit none
	
	complex :: mat(ny,nx,ndof,ndof,5)
	real    :: Bh(ny,nx,ndof,ndof), Dh(ny,nx,ndof,ndof)
	real    :: Dhi(ny,nx,ndof,ndof), Bhi(ny,nx,6)
	real    :: spgl(ny,nx), spg2l(ny,nx), Vh(ny,nx,6)
	real    :: dtl(ny,nx)
	logical :: calcd

	real    :: a1, a2, a3, a4, a5
	complex :: c1, c2, c3, c4, c5

	real :: detainv, detasinv, rcalcd
	
	integer :: iv, j, idof, jdof
!=======================================================================================================!
	rcalcd = zero
	if (calcd) rcalcd = one

	if (yper) return

	detainv  = one / deta
	detasinv = one / deta**2
	
!=======================================================================================================!
!.... use higher-order tangent on body node
!=======================================================================================================!
	j = 1

	a1 = alfa * gc1 * detainv
	a2 = alfa * gc2 * detainv
	a3 = alfa * gc3 * detainv
	a4 = alfa * gc4 * detainv
	a5 = alfa * gc5 * detainv

        do idof = 1, ndof
	  do jdof = 1, ndof
	    do iv = 1, nx
	      mat(j,iv,idof,jdof,3) = dtl(j,iv) * ( a1 * Bh(j,iv,idof,jdof) + &
                                      rcalcd * alfa * (  Dh(j,iv,idof,jdof) + &
				                   im * Dhi(j,iv,idof,jdof) ) )
	      mat(j,iv,idof,jdof,4) = a2 * dtl(j,iv) * Bh(j,iv,idof,jdof)
	      mat(j,iv,idof,jdof,5) = a3 * dtl(j,iv) * Bh(j,iv,idof,jdof)
	      mat(j,iv,idof,jdof,1) = a4 * dtl(j,iv) * Bh(j,iv,idof,jdof)
	      mat(j,iv,idof,jdof,2) = a5 * dtl(j,iv) * Bh(j,iv,idof,jdof)
	    end do
	  end do
	end do
	
!.... \hat{B}_i term

	c1 = im * alfa * gc1 * detainv
	c2 = im * alfa * gc2 * detainv
	c3 = im * alfa * gc3 * detainv
	c4 = im * alfa * gc4 * detainv
	c5 = im * alfa * gc5 * detainv

	do iv = 1, nx

	mat(j,iv,2,4,3) = mat(j,iv,2,4,3) - c1 * dtl(j,iv) * Bhi(j,iv,3)
	mat(j,iv,3,4,3) = mat(j,iv,3,4,3) - c1 * dtl(j,iv) * Bhi(j,iv,4)
	mat(j,iv,4,2,3) = mat(j,iv,4,2,3) - c1 * dtl(j,iv) * Bhi(j,iv,3)
	mat(j,iv,4,3,3) = mat(j,iv,4,3,3) - c1 * dtl(j,iv) * Bhi(j,iv,4)

	mat(j,iv,2,4,4) = mat(j,iv,2,4,4) - c2 * dtl(j,iv) * Bhi(j,iv,3)
	mat(j,iv,3,4,4) = mat(j,iv,3,4,4) - c2 * dtl(j,iv) * Bhi(j,iv,4)
	mat(j,iv,4,2,4) = mat(j,iv,4,2,4) - c2 * dtl(j,iv) * Bhi(j,iv,3)
	mat(j,iv,4,3,4) = mat(j,iv,4,3,4) - c2 * dtl(j,iv) * Bhi(j,iv,4)

	mat(j,iv,2,4,5) = mat(j,iv,2,4,5) - c3 * dtl(j,iv) * Bhi(j,iv,3)
	mat(j,iv,3,4,5) = mat(j,iv,3,4,5) - c3 * dtl(j,iv) * Bhi(j,iv,4)
	mat(j,iv,4,2,5) = mat(j,iv,4,2,5) - c3 * dtl(j,iv) * Bhi(j,iv,3)
	mat(j,iv,4,3,5) = mat(j,iv,4,3,5) - c3 * dtl(j,iv) * Bhi(j,iv,4)

	mat(j,iv,2,4,1) = mat(j,iv,2,4,1) - c4 * dtl(j,iv) * Bhi(j,iv,3)
	mat(j,iv,3,4,1) = mat(j,iv,3,4,1) - c4 * dtl(j,iv) * Bhi(j,iv,4)
	mat(j,iv,4,2,1) = mat(j,iv,4,2,1) - c4 * dtl(j,iv) * Bhi(j,iv,3)
	mat(j,iv,4,3,1) = mat(j,iv,4,3,1) - c4 * dtl(j,iv) * Bhi(j,iv,4)

	mat(j,iv,2,4,2) = mat(j,iv,2,4,2) - c5 * dtl(j,iv) * Bhi(j,iv,3)
	mat(j,iv,3,4,2) = mat(j,iv,3,4,2) - c5 * dtl(j,iv) * Bhi(j,iv,4)
	mat(j,iv,4,2,2) = mat(j,iv,4,2,2) - c5 * dtl(j,iv) * Bhi(j,iv,3)
	mat(j,iv,4,3,2) = mat(j,iv,4,3,2) - c5 * dtl(j,iv) * Bhi(j,iv,4)
	  
	end do

!.... \hat{V}_{\eta\eta} term

	a1 = alfa * dd1 * detasinv
	a2 = alfa * dd2 * detasinv
	a3 = alfa * dd3 * detasinv
	a4 = alfa * dd4 * detasinv
	a5 = alfa * dd5 * detasinv

	do iv = 1, nx
	
	mat(j,iv,2,2,3) = mat(j,iv,2,2,3) - a1 * dtl(j,iv) * Vh(j,iv,1)
	mat(j,iv,2,3,3) = mat(j,iv,2,3,3) - a1 * dtl(j,iv) * Vh(j,iv,5)
	mat(j,iv,3,2,3) = mat(j,iv,3,2,3) - a1 * dtl(j,iv) * Vh(j,iv,6)
	mat(j,iv,3,3,3) = mat(j,iv,3,3,3) - a1 * dtl(j,iv) * Vh(j,iv,2)
	mat(j,iv,4,4,3) = mat(j,iv,4,4,3) - a1 * dtl(j,iv) * Vh(j,iv,3)
	mat(j,iv,5,5,3) = mat(j,iv,5,5,3) - a1 * dtl(j,iv) * Vh(j,iv,4)

	mat(j,iv,2,2,4) = mat(j,iv,2,2,4) - a2 * dtl(j,iv) * Vh(j,iv,1)
	mat(j,iv,2,3,4) = mat(j,iv,2,3,4) - a2 * dtl(j,iv) * Vh(j,iv,5)
	mat(j,iv,3,2,4) = mat(j,iv,3,2,4) - a2 * dtl(j,iv) * Vh(j,iv,6)
	mat(j,iv,3,3,4) = mat(j,iv,3,3,4) - a2 * dtl(j,iv) * Vh(j,iv,2)
	mat(j,iv,4,4,4) = mat(j,iv,4,4,4) - a2 * dtl(j,iv) * Vh(j,iv,3)
	mat(j,iv,5,5,4) = mat(j,iv,5,5,4) - a2 * dtl(j,iv) * Vh(j,iv,4)

	mat(j,iv,2,2,5) = mat(j,iv,2,2,5) - a3 * dtl(j,iv) * Vh(j,iv,1)
	mat(j,iv,2,3,5) = mat(j,iv,2,3,5) - a3 * dtl(j,iv) * Vh(j,iv,5)
	mat(j,iv,3,2,5) = mat(j,iv,3,2,5) - a3 * dtl(j,iv) * Vh(j,iv,6)
	mat(j,iv,3,3,5) = mat(j,iv,3,3,5) - a3 * dtl(j,iv) * Vh(j,iv,2)
	mat(j,iv,4,4,5) = mat(j,iv,4,4,5) - a3 * dtl(j,iv) * Vh(j,iv,3)
	mat(j,iv,5,5,5) = mat(j,iv,5,5,5) - a3 * dtl(j,iv) * Vh(j,iv,4)

	mat(j,iv,2,2,1) = mat(j,iv,2,2,1) - a4 * dtl(j,iv) * Vh(j,iv,1)
	mat(j,iv,2,3,1) = mat(j,iv,2,3,1) - a4 * dtl(j,iv) * Vh(j,iv,5)
	mat(j,iv,3,2,1) = mat(j,iv,3,2,1) - a4 * dtl(j,iv) * Vh(j,iv,6)
	mat(j,iv,3,3,1) = mat(j,iv,3,3,1) - a4 * dtl(j,iv) * Vh(j,iv,2)
	mat(j,iv,4,4,1) = mat(j,iv,4,4,1) - a4 * dtl(j,iv) * Vh(j,iv,3)
	mat(j,iv,5,5,1) = mat(j,iv,5,5,1) - a4 * dtl(j,iv) * Vh(j,iv,4)

	mat(j,iv,2,2,2) = mat(j,iv,2,2,2) - a5 * dtl(j,iv) * Vh(j,iv,1)
	mat(j,iv,2,3,2) = mat(j,iv,2,3,2) - a5 * dtl(j,iv) * Vh(j,iv,5)
	mat(j,iv,3,2,2) = mat(j,iv,3,2,2) - a5 * dtl(j,iv) * Vh(j,iv,6)
	mat(j,iv,3,3,2) = mat(j,iv,3,3,2) - a5 * dtl(j,iv) * Vh(j,iv,2)
	mat(j,iv,4,4,2) = mat(j,iv,4,4,2) - a5 * dtl(j,iv) * Vh(j,iv,3)
	mat(j,iv,5,5,2) = mat(j,iv,5,5,2) - a5 * dtl(j,iv) * Vh(j,iv,4)

!.... I term

	mat(j,iv,1,1,3) = mat(j,iv,1,1,3) + one
	mat(j,iv,2,2,3) = mat(j,iv,2,2,3) + one
	mat(j,iv,3,3,3) = mat(j,iv,3,3,3) + one
	mat(j,iv,4,4,3) = mat(j,iv,4,4,3) + one
	mat(j,iv,5,5,3) = mat(j,iv,5,5,3) + one

	end do

!.... sponge term

	if (.not. calcd) then

	if (ispg .eq. 1) then
	
	  do iv = 1, nx
	    mat(j,iv,1,1,3) = mat(j,iv,1,1,3) + alfa * dtl(j,iv) * spgl(j,iv)
	    mat(j,iv,2,2,3) = mat(j,iv,2,2,3) + alfa * dtl(j,iv) * spgl(j,iv)
	    mat(j,iv,3,3,3) = mat(j,iv,3,3,3) + alfa * dtl(j,iv) * spgl(j,iv)
	    mat(j,iv,4,4,3) = mat(j,iv,4,4,3) + alfa * dtl(j,iv) * spgl(j,iv)
	    mat(j,iv,5,5,3) = mat(j,iv,5,5,3) + alfa * dtl(j,iv) * spgl(j,iv)
	  end do

	else if (ispg .ge. 2) then
	
	  do iv = 1, nx
	    mat(j,iv,1,1,3) = mat(j,iv,1,1,3) + alfa * dtl(j,iv) * (spgl(j,iv) + spg2l(j,iv))
	    mat(j,iv,2,2,3) = mat(j,iv,2,2,3) + alfa * dtl(j,iv) * (spgl(j,iv) + spg2l(j,iv))
	    mat(j,iv,3,3,3) = mat(j,iv,3,3,3) + alfa * dtl(j,iv) * (spgl(j,iv) + spg2l(j,iv))
	    mat(j,iv,4,4,3) = mat(j,iv,4,4,3) + alfa * dtl(j,iv) * (spgl(j,iv) + spg2l(j,iv))
	    mat(j,iv,5,5,3) = mat(j,iv,5,5,3) + alfa * dtl(j,iv) * (spgl(j,iv) + spg2l(j,iv))
	  end do
	
	end if

	end if
!=======================================================================================================!
!.... use higher-order tangent on first node off the body
!=======================================================================================================!
	j = 2

	a1 = alfa * gb1 * detainv
	a2 = alfa * gb2 * detainv
	a3 = alfa * gb3 * detainv
	a4 = alfa * gb4 * detainv
	a5 = alfa * gb5 * detainv

        do idof = 1, ndof
	  do jdof = 1, ndof
	    do iv = 1, nx
	      mat(j,iv,idof,jdof,2) = a1 * dtl(j,iv) * Bh(j,iv,idof,jdof)
	      mat(j,iv,idof,jdof,3) = dtl(j,iv) * ( a2 * Bh(j,iv,idof,jdof) + &
                                      rcalcd * alfa * (  Dh(j,iv,idof,jdof) + &
				                   im * Dhi(j,iv,idof,jdof) ) )
	      mat(j,iv,idof,jdof,4) = a3 * dtl(j,iv) * Bh(j,iv,idof,jdof)
	      mat(j,iv,idof,jdof,5) = a4 * dtl(j,iv) * Bh(j,iv,idof,jdof)
	      mat(j,iv,idof,jdof,1) = a5 * dtl(j,iv) * Bh(j,iv,idof,jdof)
	    end do
	  end do
	end do
	
!.... \hat{B}_i term

	c1 = im * alfa * gb1 * detainv
	c2 = im * alfa * gb2 * detainv
	c3 = im * alfa * gb3 * detainv
	c4 = im * alfa * gb4 * detainv
	c5 = im * alfa * gb5 * detainv

	do iv = 1, nx

	mat(j,iv,2,4,2) = mat(j,iv,2,4,2) - c1 * dtl(j,iv) * Bhi(j,iv,3)
	mat(j,iv,3,4,2) = mat(j,iv,3,4,2) - c1 * dtl(j,iv) * Bhi(j,iv,4)
	mat(j,iv,4,2,2) = mat(j,iv,4,2,2) - c1 * dtl(j,iv) * Bhi(j,iv,3)
	mat(j,iv,4,3,2) = mat(j,iv,4,3,2) - c1 * dtl(j,iv) * Bhi(j,iv,4)

	mat(j,iv,2,4,3) = mat(j,iv,2,4,3) - c2 * dtl(j,iv) * Bhi(j,iv,3)
	mat(j,iv,3,4,3) = mat(j,iv,3,4,3) - c2 * dtl(j,iv) * Bhi(j,iv,4)
	mat(j,iv,4,2,3) = mat(j,iv,4,2,3) - c2 * dtl(j,iv) * Bhi(j,iv,3)
	mat(j,iv,4,3,3) = mat(j,iv,4,3,3) - c2 * dtl(j,iv) * Bhi(j,iv,4)

	mat(j,iv,2,4,4) = mat(j,iv,2,4,4) - c3 * dtl(j,iv) * Bhi(j,iv,3)
	mat(j,iv,3,4,4) = mat(j,iv,3,4,4) - c3 * dtl(j,iv) * Bhi(j,iv,4)
	mat(j,iv,4,2,4) = mat(j,iv,4,2,4) - c3 * dtl(j,iv) * Bhi(j,iv,3)
	mat(j,iv,4,3,4) = mat(j,iv,4,3,4) - c3 * dtl(j,iv) * Bhi(j,iv,4)

	mat(j,iv,2,4,5) = mat(j,iv,2,4,5) - c4 * dtl(j,iv) * Bhi(j,iv,3)
	mat(j,iv,3,4,5) = mat(j,iv,3,4,5) - c4 * dtl(j,iv) * Bhi(j,iv,4)
	mat(j,iv,4,2,5) = mat(j,iv,4,2,5) - c4 * dtl(j,iv) * Bhi(j,iv,3)
	mat(j,iv,4,3,5) = mat(j,iv,4,3,5) - c4 * dtl(j,iv) * Bhi(j,iv,4)

	mat(j,iv,2,4,1) = mat(j,iv,2,4,1) - c5 * dtl(j,iv) * Bhi(j,iv,3)
	mat(j,iv,3,4,1) = mat(j,iv,3,4,1) - c5 * dtl(j,iv) * Bhi(j,iv,4)
	mat(j,iv,4,2,1) = mat(j,iv,4,2,1) - c5 * dtl(j,iv) * Bhi(j,iv,3)
	mat(j,iv,4,3,1) = mat(j,iv,4,3,1) - c5 * dtl(j,iv) * Bhi(j,iv,4)
	  
	end do

!.... \hat{V}_{\eta\eta} term

	a1 = alfa * db1 * detasinv
	a2 = alfa * db2 * detasinv
	a3 = alfa * db3 * detasinv
	a4 = alfa * db4 * detasinv
	a5 = alfa * db5 * detasinv

	do iv = 1, nx
	
	mat(j,iv,2,2,2) = mat(j,iv,2,2,2) - a1 * dtl(j,iv) * Vh(j,iv,1)
	mat(j,iv,2,3,2) = mat(j,iv,2,3,2) - a1 * dtl(j,iv) * Vh(j,iv,5)
	mat(j,iv,3,2,2) = mat(j,iv,3,2,2) - a1 * dtl(j,iv) * Vh(j,iv,6)
	mat(j,iv,3,3,2) = mat(j,iv,3,3,2) - a1 * dtl(j,iv) * Vh(j,iv,2)
	mat(j,iv,4,4,2) = mat(j,iv,4,4,2) - a1 * dtl(j,iv) * Vh(j,iv,3)
	mat(j,iv,5,5,2) = mat(j,iv,5,5,2) - a1 * dtl(j,iv) * Vh(j,iv,4)

	mat(j,iv,2,2,3) = mat(j,iv,2,2,3) - a2 * dtl(j,iv) * Vh(j,iv,1)
	mat(j,iv,2,3,3) = mat(j,iv,2,3,3) - a2 * dtl(j,iv) * Vh(j,iv,5)
	mat(j,iv,3,2,3) = mat(j,iv,3,2,3) - a2 * dtl(j,iv) * Vh(j,iv,6)
	mat(j,iv,3,3,3) = mat(j,iv,3,3,3) - a2 * dtl(j,iv) * Vh(j,iv,2)
	mat(j,iv,4,4,3) = mat(j,iv,4,4,3) - a2 * dtl(j,iv) * Vh(j,iv,3)
	mat(j,iv,5,5,3) = mat(j,iv,5,5,3) - a2 * dtl(j,iv) * Vh(j,iv,4)

	mat(j,iv,2,2,4) = mat(j,iv,2,2,4) - a3 * dtl(j,iv) * Vh(j,iv,1)
	mat(j,iv,2,3,4) = mat(j,iv,2,3,4) - a3 * dtl(j,iv) * Vh(j,iv,5)
	mat(j,iv,3,2,4) = mat(j,iv,3,2,4) - a3 * dtl(j,iv) * Vh(j,iv,6)
	mat(j,iv,3,3,4) = mat(j,iv,3,3,4) - a3 * dtl(j,iv) * Vh(j,iv,2)
	mat(j,iv,4,4,4) = mat(j,iv,4,4,4) - a3 * dtl(j,iv) * Vh(j,iv,3)
	mat(j,iv,5,5,4) = mat(j,iv,5,5,4) - a3 * dtl(j,iv) * Vh(j,iv,4)

	mat(j,iv,2,2,5) = mat(j,iv,2,2,5) - a4 * dtl(j,iv) * Vh(j,iv,1)
	mat(j,iv,2,3,5) = mat(j,iv,2,3,5) - a4 * dtl(j,iv) * Vh(j,iv,5)
	mat(j,iv,3,2,5) = mat(j,iv,3,2,5) - a4 * dtl(j,iv) * Vh(j,iv,6)
	mat(j,iv,3,3,5) = mat(j,iv,3,3,5) - a4 * dtl(j,iv) * Vh(j,iv,2)
	mat(j,iv,4,4,5) = mat(j,iv,4,4,5) - a4 * dtl(j,iv) * Vh(j,iv,3)
	mat(j,iv,5,5,5) = mat(j,iv,5,5,5) - a4 * dtl(j,iv) * Vh(j,iv,4)

	mat(j,iv,2,2,1) = mat(j,iv,2,2,1) - a5 * dtl(j,iv) * Vh(j,iv,1)
	mat(j,iv,2,3,1) = mat(j,iv,2,3,1) - a5 * dtl(j,iv) * Vh(j,iv,5)
	mat(j,iv,3,2,1) = mat(j,iv,3,2,1) - a5 * dtl(j,iv) * Vh(j,iv,6)
	mat(j,iv,3,3,1) = mat(j,iv,3,3,1) - a5 * dtl(j,iv) * Vh(j,iv,2)
	mat(j,iv,4,4,1) = mat(j,iv,4,4,1) - a5 * dtl(j,iv) * Vh(j,iv,3)
	mat(j,iv,5,5,1) = mat(j,iv,5,5,1) - a5 * dtl(j,iv) * Vh(j,iv,4)

!.... I term

	mat(j,iv,1,1,3) = mat(j,iv,1,1,3) + one
	mat(j,iv,2,2,3) = mat(j,iv,2,2,3) + one
	mat(j,iv,3,3,3) = mat(j,iv,3,3,3) + one
	mat(j,iv,4,4,3) = mat(j,iv,4,4,3) + one
	mat(j,iv,5,5,3) = mat(j,iv,5,5,3) + one

	end do

!.... implicit damping term

	if (eps_e .ne. zero) then

	do iv = 1, nx
	  mat(j,iv,1,1,2) = mat(j,iv,1,1,2) + alfa * dtl(j,iv) * eps_e * buff(j,iv) * (-one)
	  mat(j,iv,2,2,2) = mat(j,iv,2,2,2) + alfa * dtl(j,iv) * eps_e * buff(j,iv) * (-one)
	  mat(j,iv,3,3,2) = mat(j,iv,3,3,2) + alfa * dtl(j,iv) * eps_e * buff(j,iv) * (-one)
	  mat(j,iv,4,4,2) = mat(j,iv,4,4,2) + alfa * dtl(j,iv) * eps_e * buff(j,iv) * (-one)
	  mat(j,iv,5,5,2) = mat(j,iv,5,5,2) + alfa * dtl(j,iv) * eps_e * buff(j,iv) * (-one)
  
	  mat(j,iv,1,1,3) = mat(j,iv,1,1,3) + alfa * dtl(j,iv) * eps_e * buff(j,iv) * two
	  mat(j,iv,2,2,3) = mat(j,iv,2,2,3) + alfa * dtl(j,iv) * eps_e * buff(j,iv) * two
	  mat(j,iv,3,3,3) = mat(j,iv,3,3,3) + alfa * dtl(j,iv) * eps_e * buff(j,iv) * two
	  mat(j,iv,4,4,3) = mat(j,iv,4,4,3) + alfa * dtl(j,iv) * eps_e * buff(j,iv) * two
	  mat(j,iv,5,5,3) = mat(j,iv,5,5,3) + alfa * dtl(j,iv) * eps_e * buff(j,iv) * two
  
	  mat(j,iv,1,1,4) = mat(j,iv,1,1,4) + alfa * dtl(j,iv) * eps_e * buff(j,iv) * (-one)
	  mat(j,iv,2,2,4) = mat(j,iv,2,2,4) + alfa * dtl(j,iv) * eps_e * buff(j,iv) * (-one)
	  mat(j,iv,3,3,4) = mat(j,iv,3,3,4) + alfa * dtl(j,iv) * eps_e * buff(j,iv) * (-one)
	  mat(j,iv,4,4,4) = mat(j,iv,4,4,4) + alfa * dtl(j,iv) * eps_e * buff(j,iv) * (-one)
	  mat(j,iv,5,5,4) = mat(j,iv,5,5,4) + alfa * dtl(j,iv) * eps_e * buff(j,iv) * (-one)
	end do
		
	end if

!.... sponge term

	if (.not. calcd) then
	
	if (ispg .eq. 1) then
	
	  do iv = 1, nx
	    mat(j,iv,1,1,3) = mat(j,iv,1,1,3) + alfa * dtl(j,iv) * spgl(j,iv)
	    mat(j,iv,2,2,3) = mat(j,iv,2,2,3) + alfa * dtl(j,iv) * spgl(j,iv)
	    mat(j,iv,3,3,3) = mat(j,iv,3,3,3) + alfa * dtl(j,iv) * spgl(j,iv)
	    mat(j,iv,4,4,3) = mat(j,iv,4,4,3) + alfa * dtl(j,iv) * spgl(j,iv)
	    mat(j,iv,5,5,3) = mat(j,iv,5,5,3) + alfa * dtl(j,iv) * spgl(j,iv)
	  end do
		
	else if (ispg .ge. 2) then
	
	  do iv = 1, nx
	    mat(j,iv,1,1,3) = mat(j,iv,1,1,3) + alfa * dtl(j,iv) * (spgl(j,iv) + spg2l(j,iv))
	    mat(j,iv,2,2,3) = mat(j,iv,2,2,3) + alfa * dtl(j,iv) * (spgl(j,iv) + spg2l(j,iv))
	    mat(j,iv,3,3,3) = mat(j,iv,3,3,3) + alfa * dtl(j,iv) * (spgl(j,iv) + spg2l(j,iv))
	    mat(j,iv,4,4,3) = mat(j,iv,4,4,3) + alfa * dtl(j,iv) * (spgl(j,iv) + spg2l(j,iv))
	    mat(j,iv,5,5,3) = mat(j,iv,5,5,3) + alfa * dtl(j,iv) * (spgl(j,iv) + spg2l(j,iv))
	  end do
		
	end if

	end if
!=======================================================================================================!
!.... use higher-order tangent on first node off the far-field boundary
!=======================================================================================================!
	j = ny-1

	a1 = -alfa * gb1 * detainv
	a2 = -alfa * gb2 * detainv
	a3 = -alfa * gb3 * detainv
	a4 = -alfa * gb4 * detainv
	a5 = -alfa * gb5 * detainv

        do idof = 1, ndof
	  do jdof = 1, ndof
	    do iv = 1, nx
	      mat(j,iv,idof,jdof,5) = a5 * dtl(j,iv) * Bh(j,iv,idof,jdof)
	      mat(j,iv,idof,jdof,1) = a4 * dtl(j,iv) * Bh(j,iv,idof,jdof)
	      mat(j,iv,idof,jdof,2) = a3 * dtl(j,iv) * Bh(j,iv,idof,jdof)
	      mat(j,iv,idof,jdof,3) = dtl(j,iv) * ( a2 * Bh(j,iv,idof,jdof) + &
                                      rcalcd * alfa * (  Dh(j,iv,idof,jdof) + &
				                   im * Dhi(j,iv,idof,jdof) ) )
	      mat(j,iv,idof,jdof,4) = a1 * dtl(j,iv) * Bh(j,iv,idof,jdof)
	    end do
	  end do
	end do
	
!.... \hat{B}_i term

	c1 = -im * alfa * gb1 * detainv
	c2 = -im * alfa * gb2 * detainv
	c3 = -im * alfa * gb3 * detainv
	c4 = -im * alfa * gb4 * detainv
	c5 = -im * alfa * gb5 * detainv

	do iv = 1, nx

	mat(j,iv,2,4,5) = mat(j,iv,2,4,5) - c5 * dtl(j,iv) * Bhi(j,iv,3)
	mat(j,iv,3,4,5) = mat(j,iv,3,4,5) - c5 * dtl(j,iv) * Bhi(j,iv,4)
	mat(j,iv,4,2,5) = mat(j,iv,4,2,5) - c5 * dtl(j,iv) * Bhi(j,iv,3)
	mat(j,iv,4,3,5) = mat(j,iv,4,3,5) - c5 * dtl(j,iv) * Bhi(j,iv,4)

	mat(j,iv,2,4,1) = mat(j,iv,2,4,1) - c4 * dtl(j,iv) * Bhi(j,iv,3)
	mat(j,iv,3,4,1) = mat(j,iv,3,4,1) - c4 * dtl(j,iv) * Bhi(j,iv,4)
	mat(j,iv,4,2,1) = mat(j,iv,4,2,1) - c4 * dtl(j,iv) * Bhi(j,iv,3)
	mat(j,iv,4,3,1) = mat(j,iv,4,3,1) - c4 * dtl(j,iv) * Bhi(j,iv,4)

	mat(j,iv,2,4,2) = mat(j,iv,2,4,2) - c3 * dtl(j,iv) * Bhi(j,iv,3)
	mat(j,iv,3,4,2) = mat(j,iv,3,4,2) - c3 * dtl(j,iv) * Bhi(j,iv,4)
	mat(j,iv,4,2,2) = mat(j,iv,4,2,2) - c3 * dtl(j,iv) * Bhi(j,iv,3)
	mat(j,iv,4,3,2) = mat(j,iv,4,3,2) - c3 * dtl(j,iv) * Bhi(j,iv,4)

	mat(j,iv,2,4,3) = mat(j,iv,2,4,3) - c2 * dtl(j,iv) * Bhi(j,iv,3)
	mat(j,iv,3,4,3) = mat(j,iv,3,4,3) - c2 * dtl(j,iv) * Bhi(j,iv,4)
	mat(j,iv,4,2,3) = mat(j,iv,4,2,3) - c2 * dtl(j,iv) * Bhi(j,iv,3)
	mat(j,iv,4,3,3) = mat(j,iv,4,3,3) - c2 * dtl(j,iv) * Bhi(j,iv,4)

	mat(j,iv,2,4,4) = mat(j,iv,2,4,4) - c1 * dtl(j,iv) * Bhi(j,iv,3)
	mat(j,iv,3,4,4) = mat(j,iv,3,4,4) - c1 * dtl(j,iv) * Bhi(j,iv,4)
	mat(j,iv,4,2,4) = mat(j,iv,4,2,4) - c1 * dtl(j,iv) * Bhi(j,iv,3)
	mat(j,iv,4,3,4) = mat(j,iv,4,3,4) - c1 * dtl(j,iv) * Bhi(j,iv,4)
	  
	end do

!.... \hat{V}_{\eta\eta} term

	a1 = alfa * db1 * detasinv
	a2 = alfa * db2 * detasinv
	a3 = alfa * db3 * detasinv
	a4 = alfa * db4 * detasinv
	a5 = alfa * db5 * detasinv

	do iv = 1, nx
	
	mat(j,iv,2,2,5) = mat(j,iv,2,2,5) - a5 * dtl(j,iv) * Vh(j,iv,1)
	mat(j,iv,2,3,5) = mat(j,iv,2,3,5) - a5 * dtl(j,iv) * Vh(j,iv,5)
	mat(j,iv,3,2,5) = mat(j,iv,3,2,5) - a5 * dtl(j,iv) * Vh(j,iv,6)
	mat(j,iv,3,3,5) = mat(j,iv,3,3,5) - a5 * dtl(j,iv) * Vh(j,iv,2)
	mat(j,iv,4,4,5) = mat(j,iv,4,4,5) - a5 * dtl(j,iv) * Vh(j,iv,3)
	mat(j,iv,5,5,5) = mat(j,iv,5,5,5) - a5 * dtl(j,iv) * Vh(j,iv,4)

	mat(j,iv,2,2,1) = mat(j,iv,2,2,1) - a4 * dtl(j,iv) * Vh(j,iv,1)
	mat(j,iv,2,3,1) = mat(j,iv,2,3,1) - a4 * dtl(j,iv) * Vh(j,iv,5)
	mat(j,iv,3,2,1) = mat(j,iv,3,2,1) - a4 * dtl(j,iv) * Vh(j,iv,6)
	mat(j,iv,3,3,1) = mat(j,iv,3,3,1) - a4 * dtl(j,iv) * Vh(j,iv,2)
	mat(j,iv,4,4,1) = mat(j,iv,4,4,1) - a4 * dtl(j,iv) * Vh(j,iv,3)
	mat(j,iv,5,5,1) = mat(j,iv,5,5,1) - a4 * dtl(j,iv) * Vh(j,iv,4)

	mat(j,iv,2,2,2) = mat(j,iv,2,2,2) - a3 * dtl(j,iv) * Vh(j,iv,1)
	mat(j,iv,2,3,2) = mat(j,iv,2,3,2) - a3 * dtl(j,iv) * Vh(j,iv,5)
	mat(j,iv,3,2,2) = mat(j,iv,3,2,2) - a3 * dtl(j,iv) * Vh(j,iv,6)
	mat(j,iv,3,3,2) = mat(j,iv,3,3,2) - a3 * dtl(j,iv) * Vh(j,iv,2)
	mat(j,iv,4,4,2) = mat(j,iv,4,4,2) - a3 * dtl(j,iv) * Vh(j,iv,3)
	mat(j,iv,5,5,2) = mat(j,iv,5,5,2) - a3 * dtl(j,iv) * Vh(j,iv,4)

	mat(j,iv,2,2,3) = mat(j,iv,2,2,3) - a2 * dtl(j,iv) * Vh(j,iv,1)
	mat(j,iv,2,3,3) = mat(j,iv,2,3,3) - a2 * dtl(j,iv) * Vh(j,iv,5)
	mat(j,iv,3,2,3) = mat(j,iv,3,2,3) - a2 * dtl(j,iv) * Vh(j,iv,6)
	mat(j,iv,3,3,3) = mat(j,iv,3,3,3) - a2 * dtl(j,iv) * Vh(j,iv,2)
	mat(j,iv,4,4,3) = mat(j,iv,4,4,3) - a2 * dtl(j,iv) * Vh(j,iv,3)
	mat(j,iv,5,5,3) = mat(j,iv,5,5,3) - a2 * dtl(j,iv) * Vh(j,iv,4)

	mat(j,iv,2,2,4) = mat(j,iv,2,2,4) - a1 * dtl(j,iv) * Vh(j,iv,1)
	mat(j,iv,2,3,4) = mat(j,iv,2,3,4) - a1 * dtl(j,iv) * Vh(j,iv,5)
	mat(j,iv,3,2,4) = mat(j,iv,3,2,4) - a1 * dtl(j,iv) * Vh(j,iv,6)
	mat(j,iv,3,3,4) = mat(j,iv,3,3,4) - a1 * dtl(j,iv) * Vh(j,iv,2)
	mat(j,iv,4,4,4) = mat(j,iv,4,4,4) - a1 * dtl(j,iv) * Vh(j,iv,3)
	mat(j,iv,5,5,4) = mat(j,iv,5,5,4) - a1 * dtl(j,iv) * Vh(j,iv,4)

!.... I term

	mat(j,iv,1,1,3) = mat(j,iv,1,1,3) + one
	mat(j,iv,2,2,3) = mat(j,iv,2,2,3) + one
	mat(j,iv,3,3,3) = mat(j,iv,3,3,3) + one
	mat(j,iv,4,4,3) = mat(j,iv,4,4,3) + one
	mat(j,iv,5,5,3) = mat(j,iv,5,5,3) + one

	end do
	
!.... implicit damping term

	if (eps_e .ne. zero) then

	do iv = 1, nx
	  mat(j,iv,1,1,2) = mat(j,iv,1,1,2) + alfa * dtl(j,iv) * eps_e * buff(j,iv) * (-one)
	  mat(j,iv,2,2,2) = mat(j,iv,2,2,2) + alfa * dtl(j,iv) * eps_e * buff(j,iv) * (-one)
	  mat(j,iv,3,3,2) = mat(j,iv,3,3,2) + alfa * dtl(j,iv) * eps_e * buff(j,iv) * (-one)
	  mat(j,iv,4,4,2) = mat(j,iv,4,4,2) + alfa * dtl(j,iv) * eps_e * buff(j,iv) * (-one)
	  mat(j,iv,5,5,2) = mat(j,iv,5,5,2) + alfa * dtl(j,iv) * eps_e * buff(j,iv) * (-one)
  
	  mat(j,iv,1,1,3) = mat(j,iv,1,1,3) + alfa * dtl(j,iv) * eps_e * buff(j,iv) * two
	  mat(j,iv,2,2,3) = mat(j,iv,2,2,3) + alfa * dtl(j,iv) * eps_e * buff(j,iv) * two
	  mat(j,iv,3,3,3) = mat(j,iv,3,3,3) + alfa * dtl(j,iv) * eps_e * buff(j,iv) * two
	  mat(j,iv,4,4,3) = mat(j,iv,4,4,3) + alfa * dtl(j,iv) * eps_e * buff(j,iv) * two
	  mat(j,iv,5,5,3) = mat(j,iv,5,5,3) + alfa * dtl(j,iv) * eps_e * buff(j,iv) * two
  
	  mat(j,iv,1,1,4) = mat(j,iv,1,1,4) + alfa * dtl(j,iv) * eps_e * buff(j,iv) * (-one)
	  mat(j,iv,2,2,4) = mat(j,iv,2,2,4) + alfa * dtl(j,iv) * eps_e * buff(j,iv) * (-one)
	  mat(j,iv,3,3,4) = mat(j,iv,3,3,4) + alfa * dtl(j,iv) * eps_e * buff(j,iv) * (-one)
	  mat(j,iv,4,4,4) = mat(j,iv,4,4,4) + alfa * dtl(j,iv) * eps_e * buff(j,iv) * (-one)
	  mat(j,iv,5,5,4) = mat(j,iv,5,5,4) + alfa * dtl(j,iv) * eps_e * buff(j,iv) * (-one)
	end do
	
	end if

!.... sponge term

	if (.not. calcd) then

	if (ispg .eq. 1) then

	  do iv = 1, nx
	    mat(j,iv,1,1,3) = mat(j,iv,1,1,3) + alfa * dtl(j,iv) * spgl(j,iv)
	    mat(j,iv,2,2,3) = mat(j,iv,2,2,3) + alfa * dtl(j,iv) * spgl(j,iv)
	    mat(j,iv,3,3,3) = mat(j,iv,3,3,3) + alfa * dtl(j,iv) * spgl(j,iv)
	    mat(j,iv,4,4,3) = mat(j,iv,4,4,3) + alfa * dtl(j,iv) * spgl(j,iv)
	    mat(j,iv,5,5,3) = mat(j,iv,5,5,3) + alfa * dtl(j,iv) * spgl(j,iv)
	  end do

	else if (ispg .ge. 2) then

	  do iv = 1, nx
	    mat(j,iv,1,1,3) = mat(j,iv,1,1,3) + alfa * dtl(j,iv) * (spgl(j,iv) + spg2l(j,iv))
	    mat(j,iv,2,2,3) = mat(j,iv,2,2,3) + alfa * dtl(j,iv) * (spgl(j,iv) + spg2l(j,iv))
	    mat(j,iv,3,3,3) = mat(j,iv,3,3,3) + alfa * dtl(j,iv) * (spgl(j,iv) + spg2l(j,iv))
	    mat(j,iv,4,4,3) = mat(j,iv,4,4,3) + alfa * dtl(j,iv) * (spgl(j,iv) + spg2l(j,iv))
	    mat(j,iv,5,5,3) = mat(j,iv,5,5,3) + alfa * dtl(j,iv) * (spgl(j,iv) + spg2l(j,iv))
	  end do

	end if

	end if
!=======================================================================================================!
!.... use higher-order tangent on the far-field boundary
!=======================================================================================================!
	j = ny

	a1 = -alfa * gc1 * detainv
	a2 = -alfa * gc2 * detainv
	a3 = -alfa * gc3 * detainv
	a4 = -alfa * gc4 * detainv
	a5 = -alfa * gc5 * detainv

        do idof = 1, ndof
	  do jdof = 1, ndof
	    do iv = 1, nx
	      mat(j,iv,idof,jdof,4) = a5 * dtl(j,iv) * Bh(j,iv,idof,jdof)
	      mat(j,iv,idof,jdof,5) = a4 * dtl(j,iv) * Bh(j,iv,idof,jdof)
	      mat(j,iv,idof,jdof,1) = a3 * dtl(j,iv) * Bh(j,iv,idof,jdof)
	      mat(j,iv,idof,jdof,2) = a2 * dtl(j,iv) * Bh(j,iv,idof,jdof)
	      mat(j,iv,idof,jdof,3) = dtl(j,iv) * ( a1 * Bh(j,iv,idof,jdof) + &
                                      rcalcd * alfa * (  Dh(j,iv,idof,jdof) + &
				                   im * Dhi(j,iv,idof,jdof) ) )
	    end do
	  end do
	end do
	
!.... \hat{B}_i term

	c1 = -im * alfa * gc1 * detainv
	c2 = -im * alfa * gc2 * detainv
	c3 = -im * alfa * gc3 * detainv
	c4 = -im * alfa * gc4 * detainv
	c5 = -im * alfa * gc5 * detainv

	do iv = 1, nx

	mat(j,iv,2,4,4) = mat(j,iv,2,4,4) - c5 * dtl(j,iv) * Bhi(j,iv,3)
	mat(j,iv,3,4,4) = mat(j,iv,3,4,4) - c5 * dtl(j,iv) * Bhi(j,iv,4)
	mat(j,iv,4,2,4) = mat(j,iv,4,2,4) - c5 * dtl(j,iv) * Bhi(j,iv,3)
	mat(j,iv,4,3,4) = mat(j,iv,4,3,4) - c5 * dtl(j,iv) * Bhi(j,iv,4)

	mat(j,iv,2,4,5) = mat(j,iv,2,4,5) - c4 * dtl(j,iv) * Bhi(j,iv,3)
	mat(j,iv,3,4,5) = mat(j,iv,3,4,5) - c4 * dtl(j,iv) * Bhi(j,iv,4)
	mat(j,iv,4,2,5) = mat(j,iv,4,2,5) - c4 * dtl(j,iv) * Bhi(j,iv,3)
	mat(j,iv,4,3,5) = mat(j,iv,4,3,5) - c4 * dtl(j,iv) * Bhi(j,iv,4)

	mat(j,iv,2,4,1) = mat(j,iv,2,4,1) - c3 * dtl(j,iv) * Bhi(j,iv,3)
	mat(j,iv,3,4,1) = mat(j,iv,3,4,1) - c3 * dtl(j,iv) * Bhi(j,iv,4)
	mat(j,iv,4,2,1) = mat(j,iv,4,2,1) - c3 * dtl(j,iv) * Bhi(j,iv,3)
	mat(j,iv,4,3,1) = mat(j,iv,4,3,1) - c3 * dtl(j,iv) * Bhi(j,iv,4)

	mat(j,iv,2,4,2) = mat(j,iv,2,4,2) - c2 * dtl(j,iv) * Bhi(j,iv,3)
	mat(j,iv,3,4,2) = mat(j,iv,3,4,2) - c2 * dtl(j,iv) * Bhi(j,iv,4)
	mat(j,iv,4,2,2) = mat(j,iv,4,2,2) - c2 * dtl(j,iv) * Bhi(j,iv,3)
	mat(j,iv,4,3,2) = mat(j,iv,4,3,2) - c2 * dtl(j,iv) * Bhi(j,iv,4)

	mat(j,iv,2,4,3) = mat(j,iv,2,4,3) - c1 * dtl(j,iv) * Bhi(j,iv,3)
	mat(j,iv,3,4,3) = mat(j,iv,3,4,3) - c1 * dtl(j,iv) * Bhi(j,iv,4)
	mat(j,iv,4,2,3) = mat(j,iv,4,2,3) - c1 * dtl(j,iv) * Bhi(j,iv,3)
	mat(j,iv,4,3,3) = mat(j,iv,4,3,3) - c1 * dtl(j,iv) * Bhi(j,iv,4)
	  
	end do

!.... \hat{V}_{\eta\eta} term

	a1 = alfa * dd1 * detasinv
	a2 = alfa * dd2 * detasinv
	a3 = alfa * dd3 * detasinv
	a4 = alfa * dd4 * detasinv
	a5 = alfa * dd5 * detasinv

	do iv = 1, nx
	
	mat(j,iv,2,2,4) = mat(j,iv,2,2,4) - a5 * dtl(j,iv) * Vh(j,iv,1)
	mat(j,iv,2,3,4) = mat(j,iv,2,3,4) - a5 * dtl(j,iv) * Vh(j,iv,5)
	mat(j,iv,3,2,4) = mat(j,iv,3,2,4) - a5 * dtl(j,iv) * Vh(j,iv,6)
	mat(j,iv,3,3,4) = mat(j,iv,3,3,4) - a5 * dtl(j,iv) * Vh(j,iv,2)
	mat(j,iv,4,4,4) = mat(j,iv,4,4,4) - a5 * dtl(j,iv) * Vh(j,iv,3)
	mat(j,iv,5,5,4) = mat(j,iv,5,5,4) - a5 * dtl(j,iv) * Vh(j,iv,4)

	mat(j,iv,2,2,5) = mat(j,iv,2,2,5) - a4 * dtl(j,iv) * Vh(j,iv,1)
	mat(j,iv,2,3,5) = mat(j,iv,2,3,5) - a4 * dtl(j,iv) * Vh(j,iv,5)
	mat(j,iv,3,2,5) = mat(j,iv,3,2,5) - a4 * dtl(j,iv) * Vh(j,iv,6)
	mat(j,iv,3,3,5) = mat(j,iv,3,3,5) - a4 * dtl(j,iv) * Vh(j,iv,2)
	mat(j,iv,4,4,5) = mat(j,iv,4,4,5) - a4 * dtl(j,iv) * Vh(j,iv,3)
	mat(j,iv,5,5,5) = mat(j,iv,5,5,5) - a4 * dtl(j,iv) * Vh(j,iv,4)

	mat(j,iv,2,2,1) = mat(j,iv,2,2,1) - a3 * dtl(j,iv) * Vh(j,iv,1)
	mat(j,iv,2,3,1) = mat(j,iv,2,3,1) - a3 * dtl(j,iv) * Vh(j,iv,5)
	mat(j,iv,3,2,1) = mat(j,iv,3,2,1) - a3 * dtl(j,iv) * Vh(j,iv,6)
	mat(j,iv,3,3,1) = mat(j,iv,3,3,1) - a3 * dtl(j,iv) * Vh(j,iv,2)
	mat(j,iv,4,4,1) = mat(j,iv,4,4,1) - a3 * dtl(j,iv) * Vh(j,iv,3)
	mat(j,iv,5,5,1) = mat(j,iv,5,5,1) - a3 * dtl(j,iv) * Vh(j,iv,4)

	mat(j,iv,2,2,2) = mat(j,iv,2,2,2) - a2 * dtl(j,iv) * Vh(j,iv,1)
	mat(j,iv,2,3,2) = mat(j,iv,2,3,2) - a2 * dtl(j,iv) * Vh(j,iv,5)
	mat(j,iv,3,2,2) = mat(j,iv,3,2,2) - a2 * dtl(j,iv) * Vh(j,iv,6)
	mat(j,iv,3,3,2) = mat(j,iv,3,3,2) - a2 * dtl(j,iv) * Vh(j,iv,2)
	mat(j,iv,4,4,2) = mat(j,iv,4,4,2) - a2 * dtl(j,iv) * Vh(j,iv,3)
	mat(j,iv,5,5,2) = mat(j,iv,5,5,2) - a2 * dtl(j,iv) * Vh(j,iv,4)

	mat(j,iv,2,2,3) = mat(j,iv,2,2,3) - a1 * dtl(j,iv) * Vh(j,iv,1)
	mat(j,iv,2,3,3) = mat(j,iv,2,3,3) - a1 * dtl(j,iv) * Vh(j,iv,5)
	mat(j,iv,3,2,3) = mat(j,iv,3,2,3) - a1 * dtl(j,iv) * Vh(j,iv,6)
	mat(j,iv,3,3,3) = mat(j,iv,3,3,3) - a1 * dtl(j,iv) * Vh(j,iv,2)
	mat(j,iv,4,4,3) = mat(j,iv,4,4,3) - a1 * dtl(j,iv) * Vh(j,iv,3)
	mat(j,iv,5,5,3) = mat(j,iv,5,5,3) - a1 * dtl(j,iv) * Vh(j,iv,4)

!.... I term

	mat(j,iv,1,1,3) = mat(j,iv,1,1,3) + one
	mat(j,iv,2,2,3) = mat(j,iv,2,2,3) + one
	mat(j,iv,3,3,3) = mat(j,iv,3,3,3) + one
	mat(j,iv,4,4,3) = mat(j,iv,4,4,3) + one
	mat(j,iv,5,5,3) = mat(j,iv,5,5,3) + one

	end do

!.... sponge term

	if (.not. calcd) then

	if (ispg .eq. 1) then
	
	  do iv = 1, nx
	    mat(j,iv,1,1,3) = mat(j,iv,1,1,3) + alfa * dtl(j,iv) * spgl(j,iv)
	    mat(j,iv,2,2,3) = mat(j,iv,2,2,3) + alfa * dtl(j,iv) * spgl(j,iv)
	    mat(j,iv,3,3,3) = mat(j,iv,3,3,3) + alfa * dtl(j,iv) * spgl(j,iv)
	    mat(j,iv,4,4,3) = mat(j,iv,4,4,3) + alfa * dtl(j,iv) * spgl(j,iv)
	    mat(j,iv,5,5,3) = mat(j,iv,5,5,3) + alfa * dtl(j,iv) * spgl(j,iv)
	  end do
      
	else if (ispg .ge. 2) then
	
	  do iv = 1, nx
	    mat(j,iv,1,1,3) = mat(j,iv,1,1,3) + alfa * dtl(j,iv) * (spgl(j,iv) + spg2l(j,iv))
	    mat(j,iv,2,2,3) = mat(j,iv,2,2,3) + alfa * dtl(j,iv) * (spgl(j,iv) + spg2l(j,iv))
	    mat(j,iv,3,3,3) = mat(j,iv,3,3,3) + alfa * dtl(j,iv) * (spgl(j,iv) + spg2l(j,iv))
	    mat(j,iv,4,4,3) = mat(j,iv,4,4,3) + alfa * dtl(j,iv) * (spgl(j,iv) + spg2l(j,iv))
	    mat(j,iv,5,5,3) = mat(j,iv,5,5,3) + alfa * dtl(j,iv) * (spgl(j,iv) + spg2l(j,iv))
	  end do
		
	end if

	end if

	return
	end
