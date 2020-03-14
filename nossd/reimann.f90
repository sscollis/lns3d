!=============================================================================!
	subroutine Reimann( n, vl, vl1, vl2, bn, xl, yl, &
	                    rhob, ub, vb, wb, tb, pb, cb ) 
!=============================================================================!
!  
!  Enforce Reimann boundary condition:
!
!  Inputs:
!
!  n:         number of nodes on the boundary
!  vl:        primative variables on the boundary
!  vl1:       primative variables one node in from the boundary
!  vl2:       primative variables two nodes in from the boundary
!  bn:        unit boundary normal
!  xl:        x-coordinates on the boundary
!  yl:        y-coordinates on the boundary
!  rhob:      rho from the potential solution on the boundary
!  ub:        u from the potential solution on the boundary
!  vb:        v from the potential solution on the boundary
!  wb:        w from the potential solution on the boundary
!  tb:        t from the potential solution on the boundary
!  pb:        p from the potential solution on the boundary
!  cb:        c from the potential solution on the boundary
!
!  Outputs:
!
!  vl:        Corrected primative variables on the boundary
!
!  Notes:
!
!  1. I extrapolate primative variables and then form the outgoing
!     Reimann.  The other option is to directly extrapolate the outgoing 
!     Reimann but I haven't tried this.  (It might work better? who knows...)
!
!  2. This could be simplified since the boundary normal vectors are now
!     of unit length
!
!  3. You could also simplify the density term since, 
!     rho^\gamma / p = \gamma \M^2
!
!  Revised: 7-31-96
!=============================================================================!
	use global
	implicit none

	integer :: n
	real :: vl(n,ndof), vl1(n,ndof), vl2(n,ndof)
	real :: bn(n,2), xl(n), yl(n)

	real :: vn(n), vt(n), vint(n,ndof)
	real :: etal(n), xil(n)
	real :: ub(n), vb(n), wb(n), tb(n), rhob(n), pb(n), cb(n)

	real :: R1(n), R2(n), c(n), rho(n), p(n)
	integer :: i
!=============================================================================!

!.... compute Vn and Vt
	
	vn(:) = bn(:,1) * vl(:,2) + bn(:,2) * vl(:,3)
	vt(:) = bn(:,2) * vl(:,2) - bn(:,1) * vl(:,3)

!.... perform extrapolation

	vint(:,:) = two * vl1(:,:) - vl2(:,:)

!.... set the Riemann invariants

	R1(:)  = bn(:,1) * ub(:) + bn(:,2) * vb(:) - two * cb / gamma1
	c(:)   = sqrt( vint(:,5) ) / Ma
	R2(:)  = bn(:,1) * vint(:,2) + bn(:,2) * vint(:,3) + &
	         two * c(:) / gamma1

	do i = 1, n

	  if ( vn(i) .le. zero ) then		! inflow

	    rho(i) = rhob(i)
	    p(i)   = pb(i)

	    vn(i)  = pt5 * ( R1(i) + R2(i) )
	    vt(i)  = bn(i,2) * ub(i) - bn(i,1) * vb(i)

	    vl(i,1) = ( rho(i)**gamma / (gamma * p(i)) * &
	      (gamma1 / four)**2 * &
	      ( R2(i) - R1(i) )**2 ) ** (one/gamma1)
	    vl(i,2) = (-bn(i,1) * vn(i) - bn(i,2) * vt(i)) / &
	      (-bn(i,1)**2 - bn(i,2)**2)
	    vl(i,3) = (-bn(i,2) * vn(i) + bn(i,1) * vt(i)) / &
	      (-bn(i,1)**2 - bn(i,2)**2)
	    vl(i,4) = tan(alpha)
	    vl(i,5) = ( Ma * gamma1 / four * ( R2(i) - R1(i) ) )**2

	  else				        ! outflow

	    rho(i) = vint(i,1)
	    p(i)   = vint(i,1) * vint(i,5) / (gamma * Ma**2)

	    vn(i)  = pt5 * ( R1(i) + R2(i) )
	    vt(i)  = bn(i,2) * vint(i,2) - bn(i,1) * vint(i,3)

	    vl(i,1) = ( rho(i)**gamma / (gamma * p(i)) * &
	      (gamma1 / four)**2 * &
	      ( R2(i) - R1(i) )**2 ) ** (one/gamma1)
	    vl(i,2) = (-bn(i,1) * vn(i) - bn(i,2) * vt(i)) / &
	      (-bn(i,1)**2 - bn(i,2)**2)
	    vl(i,3) = (-bn(i,2) * vn(i) + bn(i,1) * vt(i)) / &
	      (-bn(i,1)**2 - bn(i,2)**2)
	    vl(i,4) = vint(i,4)
	    vl(i,5) = ( Ma * gamma1 / four * ( R2(i) - R1(i) ) )**2

	  end if

	end do

	return
	end

!=============================================================================!
	subroutine ReimannRHS( n, vl, vl1, vl2, bn, xl, yl, rl, &
	                       rhob, ub, vb, wb, tb, pb, cb ) 
!  
!  Enforce Reimann boundary condition on the RHS
!  
!=============================================================================!
	use global
	implicit none

	integer :: n
	real :: vl(n,ndof), vl1(n,ndof), vl2(n,ndof), rl(n,ndof)
	real :: bn(n,2), xl(n), yl(n)

	real :: vn(n), vt(n), vint(n,ndof)
	real :: etal(n), xil(n)
	real :: ub(n), vb(n), wb(n), tb(n), rhob(n), pb(n), cb(n)

	real :: R1(n), R2(n), c(n), rho(n), p(n)
	integer :: i
!=============================================================================!

!.... compute Vn and Vt
	
	vn(:) = bn(:,1) * vl(:,2) + bn(:,2) * vl(:,3)
	vt(:) = bn(:,2) * vl(:,2) - bn(:,1) * vl(:,3)

!.... perform extrapolation

	vint(:,:) = two * vl1(:,:) - vl2(:,:)

!.... set the Riemann invariants

	R1(:)  = bn(:,1) * ub(:) + bn(:,2) * vb(:) - two * cb / gamma1
	c(:)   = sqrt( vint(:,5) ) / Ma
	R2(:)  = bn(:,1) * vint(:,2) + bn(:,2) * vint(:,3) + &
	         two * c(:) / gamma1

	do i = 1, n

	  if ( vn(i) .le. zero ) then		! inflow

	    rho(i) = rhob(i)
	    p(i)   = pb(i)

	    vn(i)  = pt5 * ( R1(i) + R2(i) )
	    vt(i)  = bn(i,2) * ub(i) - bn(i,1) * vb(i)

	    rl(i,1) = vl(i,1) - ( rho(i)**gamma / (gamma * p(i)) * &
	      (gamma1 / four)**2 * &
	      ( R2(i) - R1(i) )**2 ) ** (one/gamma1)
	    rl(i,2) = vl(i,2) - (-bn(i,1) * vn(i) - bn(i,2) * vt(i)) / &
	      (-bn(i,1)**2 - bn(i,2)**2)
	    rl(i,3) = vl(i,3) - (-bn(i,2) * vn(i) + bn(i,1) * vt(i)) / &
	      (-bn(i,1)**2 - bn(i,2)**2)
	    rl(i,4) = vl(i,4) - tan(alpha)
	    rl(i,5) = vl(i,5) - ( Ma * gamma1 / four * ( R2(i) - R1(i) ) )**2

	  else				        ! outflow

	    rho(i) = vint(i,1)
	    p(i)   = vint(i,1) * vint(i,5) / (gamma * Ma**2)

	    vn(i)  = pt5 * ( R1(i) + R2(i) )
	    vt(i)  = bn(i,2) * vint(i,2) - bn(i,1) * vint(i,3)

	    rl(i,1) = vl(i,1) - ( rho(i)**gamma / (gamma * p(i)) * &
	      (gamma1 / four)**2 * &
	      ( R2(i) - R1(i) )**2 ) ** (one/gamma1)
	    rl(i,2) = vl(i,2) - (-bn(i,1) * vn(i) - bn(i,2) * vt(i)) / &
	      (-bn(i,1)**2 - bn(i,2)**2)
	    rl(i,3) = vl(i,3) - (-bn(i,2) * vn(i) + bn(i,1) * vt(i)) / &
	      (-bn(i,1)**2 - bn(i,2)**2)
	    rl(i,4) = vl(i,4) - vint(i,4)
	    rl(i,5) = vl(i,5) - ( Ma * gamma1 / four * ( R2(i) - R1(i) ) )**2

	  end if

	end do

	return
	end

!=============================================================================!
	subroutine ReimannLHS( n, vl, vl1, vl2, bn, xl, yl, &
                               mat1, mat2, mat3, &
	                       rhob, ub, vb, wb, tb, pb, cb )
!  
!  Enforce Reimann boundary condition on the LHS
!  
!=============================================================================!
	use global
	implicit none

	integer :: n
	real :: vl(n,ndof), vl1(n,ndof), vl2(n,ndof)
	real :: bn(n,2), xl(n), yl(n)
	real :: mat1(n,ndof,ndof), mat2(n,ndof,ndof), mat3(n,ndof,ndof)

	real :: vn(n), vt(n), vint(n,ndof)
	real :: etal(n), xil(n)
	real :: ub(n), vb(n), wb(n), tb(n), rhob(n), pb(n), cb(n)

	real :: R1(n), R2(n), c(n), rho(n), p(n)
	real :: fact(n), fact1(n), fact2(n)
	integer :: i
!=============================================================================!

!.... compute Vn and Vt
	
	vn(:) = bn(:,1) * vl(:,2) + bn(:,2) * vl(:,3)
	vt(:) = bn(:,2) * vl(:,2) - bn(:,1) * vl(:,3)

!.... perform extrapolation of primative variables

	vint(:,:) = two * vl1(:,:) - vl2(:,:)

!.... set the Riemann invariants

	R1(:)  = bn(:,1) * ub(:) + bn(:,2) * vb(:) - two * cb / gamma1
	c(:)   = sqrt( vint(:,5) ) / Ma
	R2(:)  = bn(:,1) * vint(:,2) + bn(:,2) * vint(:,3) + &
	         two * c(:) / gamma1

	do i = 1, n

	  if ( vn(i) .le. zero ) then		! inflow

	  rho(i) = rhob(i)
	  p(i)   = pb(i)

	  vn(i)  = pt5 * ( R1(i) + R2(i) )
	  vt(i)  = bn(i,2) * ub(i) - bn(i,1) * vb(i)

	  fact(i) = (rho(i)**gamma/(gamma*p(i))* &
                    (gamma1/four)**2)**(one/gamma1) * &
	            two/gamma1 * (R2(i) - R1(i))**((two-gamma1)/gamma1)
	  mat1(i,1,1) = -one
	  mat2(i,1,2) = two * (fact(i) * bn(i,1))
	  mat2(i,1,3) = two * (fact(i) * bn(i,2))
	  mat2(i,1,5) = two * (fact(i) / ( gamma1 * Ma * sqrt(vint(i,5)) ))
	  mat3(i,1,2) = -(fact(i) * bn(i,1))
	  mat3(i,1,3) = -(fact(i) * bn(i,2))
	  mat3(i,1,5) = -(fact(i) / ( gamma1 * Ma * sqrt(vint(i,5)) ))
	  
	  fact(i) = pt5 * bn(i,1) / (bn(i,1)**2 + bn(i,2)**2)
	  mat1(i,2,2) = -one
	  mat2(i,2,2) = two * (fact(i) * bn(i,1))
	  mat2(i,2,3) = two * (fact(i) * bn(i,2))
	  mat2(i,2,5) = two * (fact(i) / ( gamma1 * Ma * sqrt(vint(i,5)) ))
	  mat3(i,2,2) = -(fact(i) * bn(i,1))
	  mat3(i,2,3) = -(fact(i) * bn(i,2))
	  mat3(i,2,5) = -(fact(i) / ( gamma1 * Ma * sqrt(vint(i,5)) ))
  
	  fact(i) = pt5 * bn(i,2) / (bn(i,1)**2 + bn(i,2)**2)
	  mat1(i,3,3) = -one
	  mat2(i,3,2) = two * (fact(i) * bn(i,1))
	  mat2(i,3,3) = two * (fact(i) * bn(i,2))
	  mat2(i,3,5) = two * (fact(i) / ( gamma1 * Ma * sqrt(vint(i,5)) ))
	  mat3(i,3,2) = -(fact(i) * bn(i,1))
	  mat3(i,3,3) = -(fact(i) * bn(i,2))
	  mat3(i,3,5) = -(fact(i) / ( gamma1 * Ma * sqrt(vint(i,5)) ))
  
	  mat1(i,4,4) = -one
  
	  fact(i) = (Ma * gamma1 / four)**2 * two * ( R2(i) - R1(i) )
	  mat1(i,5,5) = -one
	  mat2(i,5,2) = two * (fact(i) * bn(i,1))
	  mat2(i,5,3) = two * (fact(i) * bn(i,2))
	  mat2(i,5,5) = two * (fact(i) / ( gamma1 * Ma * sqrt(vint(i,5)) ))
	  mat3(i,5,2) = -(fact(i) * bn(i,1))
	  mat3(i,5,3) = -(fact(i) * bn(i,2))
	  mat3(i,5,5) = -(fact(i) / ( gamma1 * Ma * sqrt(vint(i,5)) ))

	  else				        ! outflow

	  rho(i) = vint(i,1)
	  p(i)   = vint(i,1) * vint(i,5) / (gamma * Ma**2)

	  vn(i)  = pt5 * ( R1(i) + R2(i) )
	  vt(i)  = bn(i,2) * vint(i,2) - bn(i,1) * vint(i,3)

	  fact1(i) = (rho(i)**gamma/(gamma*p(i)))**(one/gamma1) * &
	             (gamma1/four)**(two/gamma1) * (two/gamma1) * &
		     (R2(i) - R1(i))**((two-gamma1)/gamma1)
	  fact2(i) = (gamma1/four * (R2(i) - R1(i)))**(two/gamma1) * &
	             one/(gamma*gamma1) * &
		     (rho(i)**gamma/(gamma*p(i)))**((two-gamma)/gamma1)
	  mat1(i,1,1) = -one
	  mat2(i,1,2) = two*(fact1(i) * bn(i,1) + &
	                gamma1 * rho(i)**gamma1 / p(i))
	  mat2(i,1,3) = two*(fact1(i) * bn(i,2))
	  mat2(i,1,5) = two*(fact1(i) / ( gamma1 * Ma * &
	                sqrt(vint(i,5)) ) - &
	                fact2(i) * rho(i)**(gamma+1)/(gamma*Ma**2*p(i)**2))
	  mat3(i,1,2) = -(fact1(i) * bn(i,1) + &
	                gamma1 * rho(i)**gamma1 / p(i))
	  mat3(i,1,3) = -(fact1(i) * bn(i,2))
	  mat3(i,1,5) = -(fact1(i) / ( gamma1 * Ma * &
	                sqrt(vint(i,5)) ) - &
	                fact2(i) * rho(i)**(gamma+1)/(gamma*Ma**2*p(i)**2))
	  
	  fact1(i) = pt5 * bn(i,1) / (bn(i,1)**2 + bn(i,2)**2)
	  fact2(i) =       bn(i,2) / (bn(i,1)**2 + bn(i,2)**2)
	  mat1(i,2,2) = -one
	  mat2(i,2,2) = two*(fact1(i) * bn(i,1) + fact2(i) * bn(i,2))
	  mat2(i,2,3) = two*(fact1(i) * bn(i,2) - fact2(i) * bn(i,1))
	  mat2(i,2,5) = two*(fact1(i) / ( gamma1 * Ma * sqrt(vint(i,5)) ))
	  mat3(i,2,2) = -(fact1(i) * bn(i,1) + fact2(i) * bn(i,2))
	  mat3(i,2,3) = -(fact1(i) * bn(i,2) - fact2(i) * bn(i,1))
	  mat3(i,2,5) = -(fact1(i) / ( gamma1 * Ma * sqrt(vint(i,5)) ))
  
	  fact1(i) = pt5 * bn(i,2) / (bn(i,1)**2 + bn(i,2)**2)
	  fact2(i) =       bn(i,1) / (bn(i,1)**2 + bn(i,2)**2)
	  mat1(i,3,3) = -one
	  mat2(i,3,2) = two*(fact1(i) * bn(i,1) - fact2(i) * bn(i,2))
	  mat2(i,3,3) = two*(fact1(i) * bn(i,2) + fact2(i) * bn(i,1))
	  mat2(i,3,5) = two*(fact1(i) / ( gamma1 * Ma * sqrt(vint(i,5)) ))
	  mat3(i,3,2) = -(fact1(i) * bn(i,1) - fact2(i) * bn(i,2))
	  mat3(i,3,3) = -(fact1(i) * bn(i,2) + fact2(i) * bn(i,1))
	  mat3(i,3,5) = -(fact1(i) / ( gamma1 * Ma * sqrt(vint(i,5)) ))
  
	  mat1(i,4,4) = -one
	  mat2(i,4,4) = two
	  mat3(i,4,4) = -one
  
	  fact(i) = (Ma * gamma1 / four)**2 * two * ( R2(i) - R1(i) )
	  mat1(i,5,5) = -one
	  mat2(i,5,2) = two*(fact(i) * bn(i,1))
	  mat2(i,5,3) = two*(fact(i) * bn(i,2))
	  mat2(i,5,5) = two*(fact(i) / ( gamma1 * Ma * sqrt(vint(i,5)) ))
	  mat3(i,5,2) = -(fact(i) * bn(i,1))
	  mat3(i,5,3) = -(fact(i) * bn(i,2))
	  mat3(i,5,5) = -(fact(i) / ( gamma1 * Ma * sqrt(vint(i,5)) ))

	  end if

	end do

	return
	end
