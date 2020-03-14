!=============================================================================!
	subroutine genC(C, v)
!
!  Form the C matrix which multiplies U,z
!
!=============================================================================!
	use global
	use local
	implicit none
	
	real :: C(ny*nx,ndof,ndof), v(ny*nx,ndof)
	real :: fact1(ny*nx), fact2(ny*nx), gmsinv
!=============================================================================!

!.... Continuity equation

	C(:,1,1) = u3
	C(:,1,2) = zero
	C(:,1,3) = zero
	C(:,1,4) = rho
	C(:,1,5) = zero

!=============================================================================!

!.... Momentum equation -- x_1

	fact1  = rhoinv / Re
	gmsinv = one / (gamma * Ma**2)

	C(:,2,1) = zero
	C(:,2,2) = u3 - fact1 * g3mu
	C(:,2,3) = zero
	C(:,2,4) = -fact1 * g1lm
	C(:,2,5) = -fact1 * dmu * two * S(:,1,3)

!=============================================================================!

!.... Momentum equation -- x_2

	C(:,3,1) = zero
	C(:,3,2) = zero
	C(:,3,3) = u3 - fact1 * g3mu
	C(:,3,4) = -fact1 * g2lm
	C(:,3,5) = -fact1 * dmu * two * S(:,2,3)

!=============================================================================!

!.... Momentum equation -- x_3

	C(:,4,1) = rhoinv * t * gmsinv
	C(:,4,2) = -fact1 * g1mu
	C(:,4,3) = -fact1 * g2mu
	C(:,4,4) = u3 - fact1 * g3lm - fact1 * two * g3mu
	C(:,4,5) = gmsinv - fact1 * dlm * divu - &
                   fact1 * dmu * two * S(:,3,3)

!=============================================================================!

!.... Energy equation

	fact1 = gamma * rhoinv / (Pr * Re)
	fact2 = gamma * gamma1 * Ma**2 * rhoinv / Re

	C(:,5,1) = zero
	C(:,5,2) = -four * fact2 * mu * S(:,1,3)
	C(:,5,3) = -four * fact2 * mu * S(:,2,3)
	C(:,5,4) =  gamma1 * t - fact2 * two * lm * divu - &
		    four * fact2 * mu * S(:,3,3)
	C(:,5,5) =  u3 - fact1 * (g3con + dcon * gt(:,3))

	return
	end

!=============================================================================!
	subroutine bcC(C)
!
!  WARNING: not updated with new equations
!
!  Apply the viscous BC to outflow
!
!=============================================================================!
	use global
	use local
	implicit none
	
	real    :: C(ny*nx,ndof,ndof)
	real    :: fact1, fact2, fact3, gmsinv, ginv

	integer :: i, j, ji
!=============================================================================!
	fact1  = one / Re
	fact2  = one / Re
	gmsinv = one / (gamma * Ma**2)
!=============================================================================!
	write(*,*) 'WARNING: not updated with new equations'

	do j = 1, ny
	  i = 1
          ji = j + (i-1)*ny

	  C(ji,2,4) = zero

	  C(ji,4,2) = zero

	  i = nx
          ji = j + (i-1)*ny

	  C(ji,2,4) = zero

	  C(ji,4,2) = zero
	end do

	return
	end
