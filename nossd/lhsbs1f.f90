!=======================================================================================================!
	subroutine lhsbs1f( mat, Ah, Dh, Vh, dtl, side )
!  
!  Correct the LHS for symmetry conditions in the \xi direction.
!
!  This version supports the 4th order LHS. 
!
!  Revised: 10-16-95
!
!=======================================================================================================!
	use global
	use stencil
	use buff_mod
	implicit none
	
	real :: mat(ny,nx,ndof,ndof,5), Ah(ny,nx,ndof,ndof)
	real :: Dh(ny,nx,ndof,ndof)
	real :: Vh(ny,nx,6)
	real :: dtl(ny,nx)
	integer :: side

	real :: dxiinv, dxisinv

	integer :: i, iv, idof, jdof

	real :: isign
!=======================================================================================================!
	dxiinv  = one / dxi
	dxisinv = one / dxi**2

	if (side .eq. 1) then

!=======================================================================================================!
!.... symmetric condition on left node
!=======================================================================================================!
        do idof = 1, ndof
	  do jdof = 1, ndof
	    if (jdof.eq.3) then
	      isign = -one
	    else
	      isign = one
	    end if
	    do iv = 1, ny
	      mat(iv,1,idof,jdof,3) = alfa * dtl(iv,1) * Dh(iv,1,idof,jdof)
	      mat(iv,1,idof,jdof,4) = alfa * dtl(iv,1) * (ga3+isign*ga2) * dxiinv * Ah(iv,1,idof,jdof)
	      mat(iv,1,idof,jdof,5) = alfa * dtl(iv,1) * (ga4+isign*ga1) * dxiinv * Ah(iv,1,idof,jdof)
	      mat(iv,1,idof,jdof,1) = zero
	      mat(iv,1,idof,jdof,2) = zero
	    end do
	  end do
	end do
	
!.... \hat{V}_{\eta\eta} term

	do iv = 1, ny
	
	mat(iv,1,2,2,3) = mat(iv,1,2,2,3) - alfa * dtl(iv,1) * da3 * dxisinv * Vh(iv,1,1)
	mat(iv,1,2,3,3) = mat(iv,1,2,3,3) - alfa * dtl(iv,1) * da3 * dxisinv * Vh(iv,1,5)
	mat(iv,1,3,2,3) = mat(iv,1,3,2,3) - alfa * dtl(iv,1) * da3 * dxisinv * Vh(iv,1,6)
	mat(iv,1,3,3,3) = mat(iv,1,3,3,3) - alfa * dtl(iv,1) * da3 * dxisinv * Vh(iv,1,2)
	mat(iv,1,4,4,3) = mat(iv,1,4,4,3) - alfa * dtl(iv,1) * da3 * dxisinv * Vh(iv,1,3)
	mat(iv,1,5,5,3) = mat(iv,1,5,5,3) - alfa * dtl(iv,1) * da3 * dxisinv * Vh(iv,1,4)

	mat(iv,1,2,2,4) = mat(iv,1,2,2,4) - alfa * dtl(iv,1) * (da4+da2) * dxisinv * Vh(iv,1,1)
	mat(iv,1,2,3,4) = mat(iv,1,2,3,4) - alfa * dtl(iv,1) * (da4-da2) * dxisinv * Vh(iv,1,5)
	mat(iv,1,3,2,4) = mat(iv,1,3,2,4) - alfa * dtl(iv,1) * (da4+da2) * dxisinv * Vh(iv,1,6)
	mat(iv,1,3,3,4) = mat(iv,1,3,3,4) - alfa * dtl(iv,1) * (da4-da2) * dxisinv * Vh(iv,1,2)
	mat(iv,1,4,4,4) = mat(iv,1,4,4,4) - alfa * dtl(iv,1) * (da4+da2) * dxisinv * Vh(iv,1,3)
	mat(iv,1,5,5,4) = mat(iv,1,5,5,4) - alfa * dtl(iv,1) * (da4+da2) * dxisinv * Vh(iv,1,4)

	mat(iv,1,2,2,5) = mat(iv,1,2,2,5) - alfa * dtl(iv,1) * (da5+da1) * dxisinv * Vh(iv,1,1)
	mat(iv,1,2,3,5) = mat(iv,1,2,3,5) - alfa * dtl(iv,1) * (da5-da1) * dxisinv * Vh(iv,1,5)
	mat(iv,1,3,2,5) = mat(iv,1,3,2,5) - alfa * dtl(iv,1) * (da5+da1) * dxisinv * Vh(iv,1,6)
	mat(iv,1,3,3,5) = mat(iv,1,3,3,5) - alfa * dtl(iv,1) * (da5-da1) * dxisinv * Vh(iv,1,2)
	mat(iv,1,4,4,5) = mat(iv,1,4,4,5) - alfa * dtl(iv,1) * (da5+da1) * dxisinv * Vh(iv,1,3)
	mat(iv,1,5,5,5) = mat(iv,1,5,5,5) - alfa * dtl(iv,1) * (da5+da1) * dxisinv * Vh(iv,1,4)

!.... I term

	mat(iv,1,1,1,3) = mat(iv,1,1,1,3) + one
	mat(iv,1,2,2,3) = mat(iv,1,2,2,3) + one
	mat(iv,1,3,3,3) = mat(iv,1,3,3,3) + one
	mat(iv,1,4,4,3) = mat(iv,1,4,4,3) + one
	mat(iv,1,5,5,3) = mat(iv,1,5,5,3) + one

!.... implicit damping term

	if (eps_e .ne. zero) then

	mat(iv,1,1,1,3) = mat(iv,1,1,1,3) + alfa * dtl(iv,1) * eps_e * buff(iv,1) * fb3
	mat(iv,1,2,2,3) = mat(iv,1,2,2,3) + alfa * dtl(iv,1) * eps_e * buff(iv,1) * fb3
	mat(iv,1,3,3,3) = mat(iv,1,3,3,3) + alfa * dtl(iv,1) * eps_e * buff(iv,1) * fb3
	mat(iv,1,4,4,3) = mat(iv,1,4,4,3) + alfa * dtl(iv,1) * eps_e * buff(iv,1) * fb3
	mat(iv,1,5,5,3) = mat(iv,1,5,5,3) + alfa * dtl(iv,1) * eps_e * buff(iv,1) * fb3

	mat(iv,1,1,1,4) = mat(iv,1,1,1,4) + alfa * dtl(iv,1) * eps_e * buff(iv,1) * (fb4+fb2)
	mat(iv,1,2,2,4) = mat(iv,1,2,2,4) + alfa * dtl(iv,1) * eps_e * buff(iv,1) * (fb4+fb2)
	mat(iv,1,3,3,4) = mat(iv,1,3,3,4) + alfa * dtl(iv,1) * eps_e * buff(iv,1) * (fb4-fb2)
	mat(iv,1,4,4,4) = mat(iv,1,4,4,4) + alfa * dtl(iv,1) * eps_e * buff(iv,1) * (fb4+fb2)
	mat(iv,1,5,5,4) = mat(iv,1,5,5,4) + alfa * dtl(iv,1) * eps_e * buff(iv,1) * (fb4+fb2)

	mat(iv,1,1,1,5) = mat(iv,1,1,1,5) + alfa * dtl(iv,1) * eps_e * buff(iv,1) * (fb5+fb1)
	mat(iv,1,2,2,5) = mat(iv,1,2,2,5) + alfa * dtl(iv,1) * eps_e * buff(iv,1) * (fb5+fb1)
	mat(iv,1,3,3,5) = mat(iv,1,3,3,5) + alfa * dtl(iv,1) * eps_e * buff(iv,1) * (fb5-fb1)
	mat(iv,1,4,4,5) = mat(iv,1,4,4,5) + alfa * dtl(iv,1) * eps_e * buff(iv,1) * (fb5+fb1)
	mat(iv,1,5,5,5) = mat(iv,1,5,5,5) + alfa * dtl(iv,1) * eps_e * buff(iv,1) * (fb5+fb1)

	end if

	end do
	
!=======================================================================================================!
!.... symmetric condition on first node off the left boundary
!=======================================================================================================!
        do idof = 1, ndof
	  do jdof = 1, ndof
	    if (jdof.eq.3) then
	      isign = -one
	    else
	      isign = one
	    end if
	    do iv = 1, ny
	      mat(iv,2,idof,jdof,2) = alfa * dtl(iv,2) * ga2 * dxiinv * Ah(iv,2,idof,jdof)
	      mat(iv,2,idof,jdof,3) = alfa * dtl(iv,2) * isign*ga1 * dxiinv * Ah(iv,2,idof,jdof) + &
	                              alfa * dtl(iv,2) * Dh(iv,2,idof,jdof)
	      mat(iv,2,idof,jdof,4) = alfa * dtl(iv,2) * ga3 * dxiinv * Ah(iv,2,idof,jdof)
	      mat(iv,2,idof,jdof,5) = alfa * dtl(iv,2) * ga4 * dxiinv * Ah(iv,2,idof,jdof)
	      mat(iv,2,idof,jdof,1) = zero
	    end do
	  end do
	end do
	
!.... \hat{V}_{\eta\eta} term

	do iv = 1, ny
	
	mat(iv,2,2,2,2) = mat(iv,2,2,2,2) - alfa * dtl(iv,2) * da2 * dxisinv * Vh(iv,2,1)
	mat(iv,2,2,3,2) = mat(iv,2,2,3,2) - alfa * dtl(iv,2) * da2 * dxisinv * Vh(iv,2,5)
	mat(iv,2,3,2,2) = mat(iv,2,3,2,2) - alfa * dtl(iv,2) * da2 * dxisinv * Vh(iv,2,6)
	mat(iv,2,3,3,2) = mat(iv,2,3,3,2) - alfa * dtl(iv,2) * da2 * dxisinv * Vh(iv,2,2)
	mat(iv,2,4,4,2) = mat(iv,2,4,4,2) - alfa * dtl(iv,2) * da2 * dxisinv * Vh(iv,2,3)
	mat(iv,2,5,5,2) = mat(iv,2,5,5,2) - alfa * dtl(iv,2) * da2 * dxisinv * Vh(iv,2,4)

	mat(iv,2,2,2,3) = mat(iv,2,2,2,3) - alfa * dtl(iv,2) * (da3+da1) * dxisinv * Vh(iv,2,1)
	mat(iv,2,2,3,3) = mat(iv,2,2,3,3) - alfa * dtl(iv,2) * (da3-da1) * dxisinv * Vh(iv,2,5)
	mat(iv,2,3,2,3) = mat(iv,2,3,2,3) - alfa * dtl(iv,2) * (da3+da1) * dxisinv * Vh(iv,2,6)
	mat(iv,2,3,3,3) = mat(iv,2,3,3,3) - alfa * dtl(iv,2) * (da3-da1) * dxisinv * Vh(iv,2,2)
	mat(iv,2,4,4,3) = mat(iv,2,4,4,3) - alfa * dtl(iv,2) * (da3+da1) * dxisinv * Vh(iv,2,3)
	mat(iv,2,5,5,3) = mat(iv,2,5,5,3) - alfa * dtl(iv,2) * (da3+da1) * dxisinv * Vh(iv,2,4)

	mat(iv,2,2,2,4) = mat(iv,2,2,2,4) - alfa * dtl(iv,2) * da4 * dxisinv * Vh(iv,2,1)
	mat(iv,2,2,3,4) = mat(iv,2,2,3,4) - alfa * dtl(iv,2) * da4 * dxisinv * Vh(iv,2,5)
	mat(iv,2,3,2,4) = mat(iv,2,3,2,4) - alfa * dtl(iv,2) * da4 * dxisinv * Vh(iv,2,6)
	mat(iv,2,3,3,4) = mat(iv,2,3,3,4) - alfa * dtl(iv,2) * da4 * dxisinv * Vh(iv,2,2)
	mat(iv,2,4,4,4) = mat(iv,2,4,4,4) - alfa * dtl(iv,2) * da4 * dxisinv * Vh(iv,2,3)
	mat(iv,2,5,5,4) = mat(iv,2,5,5,4) - alfa * dtl(iv,2) * da4 * dxisinv * Vh(iv,2,4)

	mat(iv,2,2,2,5) = mat(iv,2,2,2,5) - alfa * dtl(iv,2) * da5 * dxisinv * Vh(iv,2,1)
	mat(iv,2,2,3,5) = mat(iv,2,2,3,5) - alfa * dtl(iv,2) * da5 * dxisinv * Vh(iv,2,5)
	mat(iv,2,3,2,5) = mat(iv,2,3,2,5) - alfa * dtl(iv,2) * da5 * dxisinv * Vh(iv,2,6)
	mat(iv,2,3,3,5) = mat(iv,2,3,3,5) - alfa * dtl(iv,2) * da5 * dxisinv * Vh(iv,2,2)
	mat(iv,2,4,4,5) = mat(iv,2,4,4,5) - alfa * dtl(iv,2) * da5 * dxisinv * Vh(iv,2,3)
	mat(iv,2,5,5,5) = mat(iv,2,5,5,5) - alfa * dtl(iv,2) * da5 * dxisinv * Vh(iv,2,4)

!.... I term

	mat(iv,2,1,1,3) = mat(iv,2,1,1,3) + one
	mat(iv,2,2,2,3) = mat(iv,2,2,2,3) + one
	mat(iv,2,3,3,3) = mat(iv,2,3,3,3) + one
	mat(iv,2,4,4,3) = mat(iv,2,4,4,3) + one
	mat(iv,2,5,5,3) = mat(iv,2,5,5,3) + one

!.... implicit damping term

	if (eps_e .ne. zero) then

	mat(iv,2,1,1,2) = mat(iv,2,1,1,2) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * fb2
	mat(iv,2,2,2,2) = mat(iv,2,2,2,2) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * fb2
	mat(iv,2,3,3,2) = mat(iv,2,3,3,2) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * fb2
	mat(iv,2,4,4,2) = mat(iv,2,4,4,2) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * fb2
	mat(iv,2,5,5,2) = mat(iv,2,5,5,2) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * fb2

	mat(iv,2,1,1,3) = mat(iv,2,1,1,3) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * (fb3+fb1)
	mat(iv,2,2,2,3) = mat(iv,2,2,2,3) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * (fb3+fb1)
	mat(iv,2,3,3,3) = mat(iv,2,3,3,3) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * (fb3-fb1)
	mat(iv,2,4,4,3) = mat(iv,2,4,4,3) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * (fb3+fb1)
	mat(iv,2,5,5,3) = mat(iv,2,5,5,3) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * (fb3+fb1)

	mat(iv,2,1,1,4) = mat(iv,2,1,1,4) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * fb4
	mat(iv,2,2,2,4) = mat(iv,2,2,2,4) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * fb4
	mat(iv,2,3,3,4) = mat(iv,2,3,3,4) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * fb4
	mat(iv,2,4,4,4) = mat(iv,2,4,4,4) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * fb4
	mat(iv,2,5,5,4) = mat(iv,2,5,5,4) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * fb4

	mat(iv,2,1,1,5) = mat(iv,2,1,1,5) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * fb5
	mat(iv,2,2,2,5) = mat(iv,2,2,2,5) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * fb5
	mat(iv,2,3,3,5) = mat(iv,2,3,3,5) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * fb5
	mat(iv,2,4,4,5) = mat(iv,2,4,4,5) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * fb5
	mat(iv,2,5,5,5) = mat(iv,2,5,5,5) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * fb5

	end if

	end do

	else if (side .eq. 2) then

!=======================================================================================================!
!.... symmetric condition on first node off the right boundary
!=======================================================================================================!
	i = nx - 1

        do idof = 1, ndof
	  do jdof = 1, ndof
	    if (jdof.eq.3) then
	      isign = -one
	    else
	      isign = one
	    end if
	    do iv = 1, ny
	      mat(iv,i,idof,jdof,1) = alfa * dtl(iv,i) * ga1 * dxiinv * Ah(iv,i,idof,jdof)
	      mat(iv,i,idof,jdof,2) = alfa * dtl(iv,i) * ga2 * dxiinv * Ah(iv,i,idof,jdof)
	      mat(iv,i,idof,jdof,3) = alfa * dtl(iv,i) * isign*ga4 * dxiinv * Ah(iv,i,idof,jdof) + &
	                              alfa * dtl(iv,i) * Dh(iv,i,idof,jdof)
	      mat(iv,i,idof,jdof,4) = alfa * dtl(iv,i) * ga3 * dxiinv * Ah(iv,i,idof,jdof)
	      mat(iv,i,idof,jdof,5) = zero
	    end do
	  end do
	end do
	
!.... \hat{V}_{\eta\eta} term

	do iv = 1, ny
	
	mat(iv,i,2,2,1) = mat(iv,i,2,2,1) - alfa * dtl(iv,i) * da1 * dxisinv * Vh(iv,i,1)
	mat(iv,i,2,3,1) = mat(iv,i,2,3,1) - alfa * dtl(iv,i) * da1 * dxisinv * Vh(iv,i,5)
	mat(iv,i,3,2,1) = mat(iv,i,3,2,1) - alfa * dtl(iv,i) * da1 * dxisinv * Vh(iv,i,6)
	mat(iv,i,3,3,1) = mat(iv,i,3,3,1) - alfa * dtl(iv,i) * da1 * dxisinv * Vh(iv,i,2)
	mat(iv,i,4,4,1) = mat(iv,i,4,4,1) - alfa * dtl(iv,i) * da1 * dxisinv * Vh(iv,i,3)
	mat(iv,i,5,5,1) = mat(iv,i,5,5,1) - alfa * dtl(iv,i) * da1 * dxisinv * Vh(iv,i,4)

	mat(iv,i,2,2,2) = mat(iv,i,2,2,2) - alfa * dtl(iv,i) * da2 * dxisinv * Vh(iv,i,1)
	mat(iv,i,2,3,2) = mat(iv,i,2,3,2) - alfa * dtl(iv,i) * da2 * dxisinv * Vh(iv,i,5)
	mat(iv,i,3,2,2) = mat(iv,i,3,2,2) - alfa * dtl(iv,i) * da2 * dxisinv * Vh(iv,i,6)
	mat(iv,i,3,3,2) = mat(iv,i,3,3,2) - alfa * dtl(iv,i) * da2 * dxisinv * Vh(iv,i,2)
	mat(iv,i,4,4,2) = mat(iv,i,4,4,2) - alfa * dtl(iv,i) * da2 * dxisinv * Vh(iv,i,3)
	mat(iv,i,5,5,2) = mat(iv,i,5,5,2) - alfa * dtl(iv,i) * da2 * dxisinv * Vh(iv,i,4)

	mat(iv,i,2,2,3) = mat(iv,i,2,2,3) - alfa * dtl(iv,i) * (da3+da5) * dxisinv * Vh(iv,i,1)
	mat(iv,i,2,3,3) = mat(iv,i,2,3,3) - alfa * dtl(iv,i) * (da3-da5) * dxisinv * Vh(iv,i,5)
	mat(iv,i,3,2,3) = mat(iv,i,3,2,3) - alfa * dtl(iv,i) * (da3+da5) * dxisinv * Vh(iv,i,6)
	mat(iv,i,3,3,3) = mat(iv,i,3,3,3) - alfa * dtl(iv,i) * (da3-da5) * dxisinv * Vh(iv,i,2)
	mat(iv,i,4,4,3) = mat(iv,i,4,4,3) - alfa * dtl(iv,i) * (da3+da5) * dxisinv * Vh(iv,i,3)
	mat(iv,i,5,5,3) = mat(iv,i,5,5,3) - alfa * dtl(iv,i) * (da3+da5) * dxisinv * Vh(iv,i,4)

	mat(iv,i,2,2,4) = mat(iv,i,2,2,4) - alfa * dtl(iv,i) * da4 * dxisinv * Vh(iv,i,1)
	mat(iv,i,2,3,4) = mat(iv,i,2,3,4) - alfa * dtl(iv,i) * da4 * dxisinv * Vh(iv,i,5)
	mat(iv,i,3,2,4) = mat(iv,i,3,2,4) - alfa * dtl(iv,i) * da4 * dxisinv * Vh(iv,i,6)
	mat(iv,i,3,3,4) = mat(iv,i,3,3,4) - alfa * dtl(iv,i) * da4 * dxisinv * Vh(iv,i,2)
	mat(iv,i,4,4,4) = mat(iv,i,4,4,4) - alfa * dtl(iv,i) * da4 * dxisinv * Vh(iv,i,3)
	mat(iv,i,5,5,4) = mat(iv,i,5,5,4) - alfa * dtl(iv,i) * da4 * dxisinv * Vh(iv,i,4)

!.... I term

	mat(iv,i,1,1,3) = mat(iv,i,1,1,3) + one
	mat(iv,i,2,2,3) = mat(iv,i,2,2,3) + one
	mat(iv,i,3,3,3) = mat(iv,i,3,3,3) + one
	mat(iv,i,4,4,3) = mat(iv,i,4,4,3) + one
	mat(iv,i,5,5,3) = mat(iv,i,5,5,3) + one

!.... implicit damping term

	if (eps_e .ne. zero) then

	mat(iv,i,1,1,1) = mat(iv,i,1,1,1) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * fb1
	mat(iv,i,2,2,1) = mat(iv,i,2,2,1) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * fb1
	mat(iv,i,3,3,1) = mat(iv,i,3,3,1) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * fb1
	mat(iv,i,4,4,1) = mat(iv,i,4,4,1) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * fb1
	mat(iv,i,5,5,1) = mat(iv,i,5,5,1) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * fb1

	mat(iv,i,1,1,2) = mat(iv,i,1,1,2) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * fb2
	mat(iv,i,2,2,2) = mat(iv,i,2,2,2) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * fb2
	mat(iv,i,3,3,2) = mat(iv,i,3,3,2) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * fb2
	mat(iv,i,4,4,2) = mat(iv,i,4,4,2) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * fb2
	mat(iv,i,5,5,2) = mat(iv,i,5,5,2) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * fb2

	mat(iv,i,1,1,3) = mat(iv,i,1,1,3) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * (fb3+fb5)
	mat(iv,i,2,2,3) = mat(iv,i,2,2,3) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * (fb3+fb5)
	mat(iv,i,3,3,3) = mat(iv,i,3,3,3) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * (fb3-fb5)
	mat(iv,i,4,4,3) = mat(iv,i,4,4,3) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * (fb3+fb5)
	mat(iv,i,5,5,3) = mat(iv,i,5,5,3) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * (fb3+fb5)

	mat(iv,i,1,1,4) = mat(iv,i,1,1,4) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * fb4
	mat(iv,i,2,2,4) = mat(iv,i,2,2,4) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * fb4
	mat(iv,i,3,3,4) = mat(iv,i,3,3,4) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * fb4
	mat(iv,i,4,4,4) = mat(iv,i,4,4,4) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * fb4
	mat(iv,i,5,5,4) = mat(iv,i,5,5,4) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * fb4

	end if

	end do

!=======================================================================================================!
!.... symmetric condition on left node
!=======================================================================================================!
	i = nx

        do idof = 1, ndof
	  do jdof = 1, ndof
	    if (jdof.eq.3) then
	      isign = -one
	    else
	      isign = one
	    end if
	    do iv = 1, ny
	      mat(iv,nx,idof,jdof,1) = alfa * dtl(iv,nx) * (ga1+isign*ga4) * dxiinv * Ah(iv,nx,idof,jdof)
	      mat(iv,nx,idof,jdof,2) = alfa * dtl(iv,nx) * (ga2+isign*ga3) * dxiinv * Ah(iv,nx,idof,jdof)
	      mat(iv,nx,idof,jdof,3) = alfa * dtl(iv,nx) * Dh(iv,nx,idof,jdof)
	      mat(iv,nx,idof,jdof,4) = zero
	      mat(iv,nx,idof,jdof,5) = zero
	    end do
	  end do
	end do
	
!.... \hat{V}_{\eta\eta} term

	do iv = 1, ny
	
	mat(iv,nx,2,2,1) = mat(iv,nx,2,2,1) - alfa * dtl(iv,nx) * (da1+da5) * dxisinv * Vh(iv,nx,1)
	mat(iv,nx,2,3,1) = mat(iv,nx,2,3,1) - alfa * dtl(iv,nx) * (da1-da5) * dxisinv * Vh(iv,nx,5)
	mat(iv,nx,3,2,1) = mat(iv,nx,3,2,1) - alfa * dtl(iv,nx) * (da1+da5) * dxisinv * Vh(iv,nx,6)
	mat(iv,nx,3,3,1) = mat(iv,nx,3,3,1) - alfa * dtl(iv,nx) * (da1-da5) * dxisinv * Vh(iv,nx,2)
	mat(iv,nx,4,4,1) = mat(iv,nx,4,4,1) - alfa * dtl(iv,nx) * (da1+da5) * dxisinv * Vh(iv,nx,3)
	mat(iv,nx,5,5,1) = mat(iv,nx,5,5,1) - alfa * dtl(iv,nx) * (da1+da5) * dxisinv * Vh(iv,nx,4)

	mat(iv,nx,2,2,2) = mat(iv,nx,2,2,2) - alfa * dtl(iv,nx) * (da2+da4) * dxisinv * Vh(iv,nx,1)
	mat(iv,nx,2,3,2) = mat(iv,nx,2,3,2) - alfa * dtl(iv,nx) * (da2-da4) * dxisinv * Vh(iv,nx,5)
	mat(iv,nx,3,2,2) = mat(iv,nx,3,2,2) - alfa * dtl(iv,nx) * (da2+da4) * dxisinv * Vh(iv,nx,6)
	mat(iv,nx,3,3,2) = mat(iv,nx,3,3,2) - alfa * dtl(iv,nx) * (da2-da4) * dxisinv * Vh(iv,nx,2)
	mat(iv,nx,4,4,2) = mat(iv,nx,4,4,2) - alfa * dtl(iv,nx) * (da2+da4) * dxisinv * Vh(iv,nx,3)
	mat(iv,nx,5,5,2) = mat(iv,nx,5,5,2) - alfa * dtl(iv,nx) * (da2+da4) * dxisinv * Vh(iv,nx,4)

	mat(iv,nx,2,2,3) = mat(iv,nx,2,2,3) - alfa * dtl(iv,nx) * da3 * dxisinv * Vh(iv,nx,1)
	mat(iv,nx,2,3,3) = mat(iv,nx,2,3,3) - alfa * dtl(iv,nx) * da3 * dxisinv * Vh(iv,nx,5)
	mat(iv,nx,3,2,3) = mat(iv,nx,3,2,3) - alfa * dtl(iv,nx) * da3 * dxisinv * Vh(iv,nx,6)
	mat(iv,nx,3,3,3) = mat(iv,nx,3,3,3) - alfa * dtl(iv,nx) * da3 * dxisinv * Vh(iv,nx,2)
	mat(iv,nx,4,4,3) = mat(iv,nx,4,4,3) - alfa * dtl(iv,nx) * da3 * dxisinv * Vh(iv,nx,3)
	mat(iv,nx,5,5,3) = mat(iv,nx,5,5,3) - alfa * dtl(iv,nx) * da3 * dxisinv * Vh(iv,nx,4)

!.... I term

	mat(iv,nx,1,1,3) = mat(iv,nx,1,1,3) + one
	mat(iv,nx,2,2,3) = mat(iv,nx,2,2,3) + one
	mat(iv,nx,3,3,3) = mat(iv,nx,3,3,3) + one
	mat(iv,nx,4,4,3) = mat(iv,nx,4,4,3) + one
	mat(iv,nx,5,5,3) = mat(iv,nx,5,5,3) + one

!.... implicit damping term

	if (eps_e .ne. zero) then

	mat(iv,nx,1,1,1) = mat(iv,nx,1,1,1) + alfa * dtl(iv,nx) * eps_e * buff(iv,nx) * (fb1+fb5)
	mat(iv,nx,2,2,1) = mat(iv,nx,2,2,1) + alfa * dtl(iv,nx) * eps_e * buff(iv,nx) * (fb1+fb5)
	mat(iv,nx,3,3,1) = mat(iv,nx,3,3,1) + alfa * dtl(iv,nx) * eps_e * buff(iv,nx) * (fb1-fb5)
	mat(iv,nx,4,4,1) = mat(iv,nx,4,4,1) + alfa * dtl(iv,nx) * eps_e * buff(iv,nx) * (fb1+fb5)
	mat(iv,nx,5,5,1) = mat(iv,nx,5,5,1) + alfa * dtl(iv,nx) * eps_e * buff(iv,nx) * (fb1+fb5)

	mat(iv,nx,1,1,2) = mat(iv,nx,1,1,2) + alfa * dtl(iv,nx) * eps_e * buff(iv,nx) * (fb2+fb4)
	mat(iv,nx,2,2,2) = mat(iv,nx,2,2,2) + alfa * dtl(iv,nx) * eps_e * buff(iv,nx) * (fb2+fb4)
	mat(iv,nx,3,3,2) = mat(iv,nx,3,3,2) + alfa * dtl(iv,nx) * eps_e * buff(iv,nx) * (fb2-fb4)
	mat(iv,nx,4,4,2) = mat(iv,nx,4,4,2) + alfa * dtl(iv,nx) * eps_e * buff(iv,nx) * (fb2+fb4)
	mat(iv,nx,5,5,2) = mat(iv,nx,5,5,2) + alfa * dtl(iv,nx) * eps_e * buff(iv,nx) * (fb2+fb4)

	mat(iv,nx,1,1,3) = mat(iv,nx,1,1,3) + alfa * dtl(iv,nx) * eps_e * buff(iv,nx) * fb3
	mat(iv,nx,2,2,3) = mat(iv,nx,2,2,3) + alfa * dtl(iv,nx) * eps_e * buff(iv,nx) * fb3
	mat(iv,nx,3,3,3) = mat(iv,nx,3,3,3) + alfa * dtl(iv,nx) * eps_e * buff(iv,nx) * fb3
	mat(iv,nx,4,4,3) = mat(iv,nx,4,4,3) + alfa * dtl(iv,nx) * eps_e * buff(iv,nx) * fb3
	mat(iv,nx,5,5,3) = mat(iv,nx,5,5,3) + alfa * dtl(iv,nx) * eps_e * buff(iv,nx) * fb3

	end if

	end do

	end if        ! side

	return
	end
