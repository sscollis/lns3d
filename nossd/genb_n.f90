!=======================================================================================================!
	subroutine genB_n(B, vl)
!
!  Form the B matrix which multiplies U,y
!
!  Revised: 9-19-95
!
!  This version has been parabolized in x and implements the Lele & Poinsot
!  boundary conditions for the nonlinear equations.
!=======================================================================================================!
	use global
	use local
	implicit none
	
	real :: B(ny*nx,ndof,ndof), vl(ny*nx,ndof)
	real :: fact1, fact2, fact3, fact4, gmsinv

	real :: pinf, cinf, Lx, rk, c
	real :: L1, L2, L3, L4, L5
	real :: d1, d2, d3, d4, d5, dc
	
	integer :: i, j, ji

	real :: bn1, bn2, bs1, bs2
	real :: un, us, gpn, gps, grhon, grhos, gtn, gts
	real :: gu1s, gu1n, gu2s, gu2n, gu3s, gu3n
	real :: gunn, gusn, guns, guss
!=======================================================================================================!

!.... compute some extra stuff needed for Lele's BC's

	pinf = one / (gamma * Ma**2)	! pressure at infinity
	cinf = one / Ma			! speed of sound at infinity
	Lx   = 100.0			! length of domain in x
	
	rk = sigma * ( one - Ma**2 ) * cinf / Lx
	
!=======================================================================================================!
!	S i d e   B o u n d a r i e s
!=======================================================================================================!
	if (left.eq.1 .or. right.eq.1) then

!.... do for both the left and right boundaries

	do i = 1, nx, nx-1
	  if ( (i.eq.1 .and. left.eq.1) .or. (i.eq.nx .and. right.eq.1) ) then	  
	  do j = 1, ny
	    ji = j + (i-1) * ny

!.... get the boundary normal and tangent vectors

	    if (i .eq. 1) then
	      bn1 = -m1(ji) / sqrt( m1(ji)**2 + m2(ji)**2 )
	      bn2 =  m2(ji) / sqrt( m1(ji)**2 + m2(ji)**2 )
	    else
	      bn1 =  m1(ji) / sqrt( m1(ji)**2 + m2(ji)**2 )
	      bn2 =  m2(ji) / sqrt( m1(ji)**2 + m2(ji)**2 )
	    end if

	    bs1 = -bn2
	    bs2 =  bn1

!.... compute Us and Un
	
	    un = bn1 * vl(ji,2) + bn2 * vl(ji,3)
	    us = bs1 * vl(ji,2) + bs2 * vl(ji,3)

	    c = sqrt( t(ji) ) / Ma
	    dc = pt5 / (Ma**2 * c)

!.... compute some derivatives in the boundary normal coordinate system

	    gpn = bn1 * gp(ji,1) + bn2 * gp(ji,2)
	    gps = bs1 * gp(ji,1) + bs2 * gp(ji,2)
	
	    grhon = bn1 * grho(ji,1) + bn2 * grho(ji,2)
	    grhos = bs1 * grho(ji,1) + bs2 * grho(ji,2)

	    gtn = bn1 * gt(ji,1) + bn2 * gt(ji,2)
	    gts = bs1 * gt(ji,1) + bs2 * gt(ji,2)

	    gu1n = bn1 * gu(ji,1,1) + bn2 * gu(ji,1,2)
	    gu1s = bs1 * gu(ji,1,1) + bs2 * gu(ji,1,2)

	    gu2n = bn1 * gu(ji,2,1) + bn2 * gu(ji,2,2)
	    gu2s = bs1 * gu(ji,2,1) + bs2 * gu(ji,2,2)

	    gu3n = bn1 * gu(ji,3,1) + bn2 * gu(ji,3,2)
	    gu3s = bs1 * gu(ji,3,1) + bs2 * gu(ji,3,2)

	    gunn = bn1*bn1 * gu(ji,1,1) + bn2*bn1 * gu(ji,2,1) + &
	           bn1*bn2 * gu(ji,1,2) + bn2*bn2 * gu(ji,2,2)

	    gusn = bs1*bn1 * gu(ji,1,1) + bs2*bn1 * gu(ji,2,1) + &
	           bs1*bn2 * gu(ji,1,2) + bs2*bn2 * gu(ji,2,2)

	    guns = bn1*bs1 * gu(ji,1,1) + bn2*bs1 * gu(ji,2,1) + &
	           bn1*bs2 * gu(ji,1,2) + bn2*bs2 * gu(ji,2,2)

	    guss = bs1*bs1 * gu(ji,1,1) + bs2*bs1 * gu(ji,2,1) + &
	           bs1*bs2 * gu(ji,1,2) + bs2*bs2 * gu(ji,2,2)
	    
	if (un .le. zero) then

!=======================================================================================================!
!	I n f l o w
!=======================================================================================================!

!.... Continuity equation

	gmsinv = one / ( gamma * Ma**2 )

	B(ji,1,1) = pt5/c**2*(un+c)*bn2*t(ji)*gmsinv + us*bs2
	B(ji,1,2) = pt5/c**2*(un+c)*rho(ji)*c*bn1*bn2 + rho(ji)*bs1*bs2
	B(ji,1,3) = pt5/c**2*(un+c)*rho(ji)*c*bn2*bn2 + rho(ji)*bs2*bs2
	B(ji,1,4) = zero
	B(ji,1,5) = pt5/c**2*(un+c)*bn2*rho(ji)*gmsinv

!=======================================================================================================!

!.... Momentum equation -- x_1

	fact1  = rhoinv(ji) / Re

	B(ji,2,1) = pt5*bn1*rhoinv(ji)/c*(un+c)*bn2*t(ji)*gmsinv +		&
                    bs1*rhoinv(ji)*gmsinv*t(ji)*bs2
	B(ji,2,2) = pt5*bn1*rhoinv(ji)/c*(un+c)*rho(ji)*c*bn1*bn2 + us*bs2 -	&
                    fact1 * g2mu(ji)
	B(ji,2,3) = pt5*bn1*rhoinv(ji)/c*(un+c)*rho(ji)*c*bn2*bn2 -		&
                    damp(ji) * fact1 * g1lm(ji)
	B(ji,2,4) = zero
	B(ji,2,5) = bn1/(two*rho(ji)*c)*(un+c)*bn2*rho(ji)*gmsinv +		&
		    bs1*gmsinv*bs2 -						&
                    fact1 * dmu(ji) * ( gu(ji,1,2) + damp(ji) * gu(ji,2,1) )

!=======================================================================================================!

!.... Momentum equation -- x_2

	B(ji,3,1) = pt5*bn2*rhoinv(ji)/c*(un+c)*bn2*t(ji)*gmsinv +	&
                    bs2 * rhoinv(ji) * gmsinv * t(ji) * bs2
	B(ji,3,2) = pt5*bn2*rhoinv(ji)/c*(un+c)*rho(ji)*c*bn1*bn2 -	&
                    damp(ji) * fact1 * g1mu(ji)
	B(ji,3,3) = pt5*bn2*rhoinv(ji)/c*(un+c)*rho(ji)*c*bn2*bn2 +	&
                    us*bs2 - fact1 * g2lm(ji) - fact1 * two * g2mu(ji)
	B(ji,3,4) = -fact1 * g3mu(ji)
	B(ji,3,5) = bn2/(two*rho(ji)*c)*(un+c)*bn2*rho(ji)*gmsinv +	&
                    bs2*gmsinv*bs2 - fact1 * dlm(ji) *			&
	            (damp(ji) * gu(ji,1,1) + gu(ji,2,2) + gu(ji,3,3)) - &
                    fact1 * dmu(ji) * ( gu(ji,2,2) + gu(ji,2,2) )

!=======================================================================================================!

!.... Momentum equation -- x_3

	B(ji,4,1) = zero
	B(ji,4,2) = zero
	B(ji,4,3) = -fact1 * g3lm(ji)
	B(ji,4,4) = bs2*us - fact1 * g2mu(ji)
	B(ji,4,5) = -fact1 * dmu(ji) * ( gu(ji,3,2) + gu(ji,2,3) )

!=======================================================================================================!

!.... Energy equation

	fact1 = gamma * rhoinv(ji) / (Pr * Re)
	fact2 = gamma * gamma1 * Ma**2 * rhoinv(ji) / Re

	B(ji,5,1) =  Ma**2*pt5*rhoinv(ji)*gamma1*(un+c)*bn2*gmsinv*t(ji)

	B(ji,5,2) =  Ma**2*pt5*rhoinv(ji)*gamma1*(un+c)*rho(ji)*c*bn1*bn2 +	&
                     gamma1*t(ji)*bs1*bs2 -					&
                     two*fact2*mu(ji)*( gu(ji,1,2) + damp(ji) * gu(ji,2,1) )

	B(ji,5,3) =  Ma**2*pt5*rhoinv(ji)*gamma1*(un+c)*rho(ji)*c*bn2*bn2 +	&
		     gamma1*t(ji)*bs2*bs2 - fact2 * two * lm(ji) *		&
	             (damp(ji) * gu(ji,1,1) + gu(ji,2,2) + gu(ji,3,3)) -	&
	             two * fact2 * mu(ji) * ( gu(ji,2,2) + gu(ji,2,2) )

	B(ji,5,4) = -two * fact2 * mu(ji) * ( gu(ji,3,2) + gu(ji,2,3) )

	B(ji,5,5) =  Ma**2*pt5*rhoinv(ji)*gamma1*(un+c)*bn2*rho(ji)*gmsinv +	&
		     us*bs2 - fact1 * (g2con(ji) + dcon(ji) * gt(ji,2))


	else

!=======================================================================================================!
!	O u t f l o w
!=======================================================================================================!

!.... Continuity equation

	gmsinv = one / ( gamma * Ma**2 )

	B(ji,1,1) = one/c**2*( un*(c**2*bn2-bn2*t(ji)*gmsinv) +		&
                    pt5*(un+c)*bn2*gmsinv*t(ji) ) + us*bs2
	B(ji,1,2) = one/c**2*pt5*(un+c)*rho(ji)*c*bn1*bn2 + rho(ji)*bs1*bs2
	B(ji,1,3) = one/c**2*pt5*(un+c)*rho(ji)*c*bn2*bn2 + rho(ji)*bs2*bs2
	B(ji,1,4) = zero
	B(ji,1,5) = one/c**2*( un*(-bn2*rho(ji)*gmsinv) +		&
                    pt5*(un+c)*bn2*rho(ji)*gmsinv )

!=======================================================================================================!

!.... Momentum equation -- x_1

	fact1  = rhoinv(ji) / Re

	B(ji,2,1) = pt5*bn1*rhoinv(ji)/c*(un+c)*bn2*t(ji)*gmsinv +		&
                    bs1*rhoinv(ji)*gmsinv*t(ji)*bs2
	B(ji,2,2) = pt5*bn1*rhoinv(ji)/c*(un+c)*rho(ji)*c*bn1*bn2 +		&
		    bs1*un*bs1*bn2 + us*bs2 -					&
                    fact1 * g2mu(ji)
	B(ji,2,3) = pt5*bn1*rhoinv(ji)/c*(un+c)*rho(ji)*c*bn2*bn2 +		&
		    bs1*un*bs2*bn2 -						&
                    damp(ji) * fact1 * g1lm(ji)
	B(ji,2,4) = zero
	B(ji,2,5) = bn1/(two*rho(ji)*c)*(un+c)*bn2*rho(ji)*gmsinv +		&
		    bs1*gmsinv*bs2 -						&
                    fact1 * dmu(ji) * ( gu(ji,1,2) + damp(ji) * gu(ji,2,1) )

!=======================================================================================================!

!.... Momentum equation -- x_2

	B(ji,3,1) = pt5*bn2*rhoinv(ji)/c*(un+c)*bn2*t(ji)*gmsinv +	&
                    bs2 * rhoinv(ji) * gmsinv * t(ji) * bs2
	B(ji,3,2) = pt5*bn2*rhoinv(ji)/c*(un+c)*rho(ji)*c*bn1*bn2 +	&
		    bs2*un*bs1*bn2 -					&
                    damp(ji) * fact1 * g1mu(ji)
	B(ji,3,3) = pt5*bn2*rhoinv(ji)/c*(un+c)*rho(ji)*c*bn2*bn2 +	&
                    us*bs2 + bs2*un*bs2*bn2 -				&
		    fact1 * g2lm(ji) - fact1 * two * g2mu(ji)
	B(ji,3,4) = -fact1 * g3mu(ji)
	B(ji,3,5) = bn2/(two*rho(ji)*c)*(un+c)*bn2*rho(ji)*gmsinv +	&
                    bs2*gmsinv*bs2 - fact1 * dlm(ji) *			&
	            (damp(ji) * gu(ji,1,1) + gu(ji,2,2) + gu(ji,3,3)) - &
                    fact1 * dmu(ji) * ( gu(ji,2,2) + gu(ji,2,2) )

!=======================================================================================================!

!.... Momentum equation -- x_3

	B(ji,4,1) = zero
	B(ji,4,2) = zero
	B(ji,4,3) = -fact1 * g3lm(ji)
	B(ji,4,4) = bn2*un + bs2*us - fact1 * g2mu(ji)
	B(ji,4,5) = -fact1 * dmu(ji) * ( gu(ji,3,2) + gu(ji,2,3) )

!=======================================================================================================!

!.... Energy equation

	fact1 = gamma * rhoinv(ji) / (Pr * Re)
	fact2 = gamma * gamma1 * Ma**2 * rhoinv(ji) / Re

	B(ji,5,1) = -Ma**2*rhoinv(ji)*( un*(c**2*bn2-bn2*t(ji)*gmsinv) -	&
		     pt5*gamma1*(un+c)*bn2*t(ji)*gmsinv ) 

	B(ji,5,2) =  Ma**2*rhoinv(ji)*pt5*gamma1*(un+c)*rho(ji)*c*bn1*bn2 +	&
                     gamma1*t(ji)*bs1*bs2 -					&
                     two*fact2*mu(ji)*( gu(ji,1,2) + damp(ji) * gu(ji,2,1) )

	B(ji,5,3) =  Ma**2*rhoinv(ji)*pt5*gamma1*(un+c)*rho(ji)*c*bn2*bn2 +	&
		     gamma1*t(ji)*bs2*bs2 - fact2 * two * lm(ji) *		&
	             (damp(ji) * gu(ji,1,1) + gu(ji,2,2) + gu(ji,3,3)) -	&
	             two * fact2 * mu(ji) * ( gu(ji,2,2) + gu(ji,2,2) )

	B(ji,5,4) = -two * fact2 * mu(ji) * ( gu(ji,3,2) + gu(ji,2,3) )

	B(ji,5,5) = -Ma**2*rhoinv(ji)*( un*(-bn2*rho(ji)*gmsinv) -		&
                     pt5*gamma1*(un+c)*bn2*rho(ji)*gmsinv ) +			&
		     us*bs2 - fact1 * (g2con(ji) + dcon(ji) * gt(ji,2))

	end if		! inflow/outflow

!.... old way

!.... Momentum equation -- x_2		( set tau21,1 = zero )

!	    B(ji,3,2) = zero

	  end do
	  end if
        end do
	end if
!=======================================================================================================!
!	T o p   B o u n d a r y
!=======================================================================================================!
	if (top.eq.1) then

	j = ny
	do i = 1, nx-1
	  ji = j + (i-1) * ny
    
!.... get the boundary normal and tangent vectors

	    bn1 = n1(ji) / sqrt( n1(ji)**2 + n2(ji)**2 )
	    bn2 = n2(ji) / sqrt( n1(ji)**2 + n2(ji)**2 )

	    bs1 = -bn2
	    bs2 =  bn1

!.... compute Us and Un
	
	    un = bn1 * vl(ji,2) + bn2 * vl(ji,3)
	    us = bs1 * vl(ji,2) + bs2 * vl(ji,3)

	    c = sqrt( t(ji) ) / Ma
	    dc = pt5 / (Ma**2 * c)

!.... compute some derivatives in the boundary normal coordinate system

	    gpn = bn1 * gp(ji,1) + bn2 * gp(ji,2)
	    gps = bs1 * gp(ji,1) + bs2 * gp(ji,2)
	
	    grhon = bn1 * grho(ji,1) + bn2 * grho(ji,2)
	    grhos = bs1 * grho(ji,1) + bs2 * grho(ji,2)

	    gtn = bn1 * gt(ji,1) + bn2 * gt(ji,2)
	    gts = bs1 * gt(ji,1) + bs2 * gt(ji,2)

	    gu1n = bn1 * gu(ji,1,1) + bn2 * gu(ji,1,2)
	    gu1s = bs1 * gu(ji,1,1) + bs2 * gu(ji,1,2)

	    gu2n = bn1 * gu(ji,2,1) + bn2 * gu(ji,2,2)
	    gu2s = bs1 * gu(ji,2,1) + bs2 * gu(ji,2,2)

	    gu3n = bn1 * gu(ji,3,1) + bn2 * gu(ji,3,2)
	    gu3s = bs1 * gu(ji,3,1) + bs2 * gu(ji,3,2)

	    gunn = bn1*bn1 * gu(ji,1,1) + bn2*bn1 * gu(ji,2,1) + &
	           bn1*bn2 * gu(ji,1,2) + bn2*bn2 * gu(ji,2,2)

	    gusn = bs1*bn1 * gu(ji,1,1) + bs2*bn1 * gu(ji,2,1) + &
	           bs1*bn2 * gu(ji,1,2) + bs2*bn2 * gu(ji,2,2)

	    guns = bn1*bs1 * gu(ji,1,1) + bn2*bs1 * gu(ji,2,1) + &
	           bn1*bs2 * gu(ji,1,2) + bn2*bs2 * gu(ji,2,2)

	    guss = bs1*bs1 * gu(ji,1,1) + bs2*bs1 * gu(ji,2,1) + &
	           bs1*bs2 * gu(ji,1,2) + bs2*bs2 * gu(ji,2,2)
	    
	if (un .le. zero) then

!=======================================================================================================!
!	I n f l o w
!=======================================================================================================!

!.... Continuity equation

	gmsinv = one / ( gamma * Ma**2 )

	B(ji,1,1) = pt5/c**2*(un+c)*bn2*t(ji)*gmsinv + us*bs2
	B(ji,1,2) = pt5/c**2*(un+c)*rho(ji)*c*bn1*bn2 + rho(ji)*bs1*bs2
	B(ji,1,3) = pt5/c**2*(un+c)*rho(ji)*c*bn2*bn2 + rho(ji)*bs2*bs2
	B(ji,1,4) = zero
	B(ji,1,5) = pt5/c**2*(un+c)*bn2*rho(ji)*gmsinv

!.... set u,v,w,T

!	B(ji,1,1) = gamma/c**2*(un+c)*bn2*t(ji)*gmsinv + us*bs2
!	B(ji,1,2) = gamma/c**2*(un+c)*rho(ji)*c*bn1*bn2 + rho(ji)*bs1*bs2
!	B(ji,1,3) = gamma/c**2*(un+c)*rho(ji)*c*bn2*bn2 + rho(ji)*bs2*bs2
!	B(ji,1,4) = zero
!	B(ji,1,5) = gamma/c**2*(un+c)*bn2*rho(ji)*gmsinv

!=======================================================================================================!

!.... Momentum equation -- x_1

	fact1  = rhoinv(ji) / Re

	B(ji,2,1) = pt5*bn1*rhoinv(ji)/c*(un+c)*bn2*t(ji)*gmsinv +		&
                    bs1*rhoinv(ji)*gmsinv*t(ji)*bs2
	B(ji,2,2) = pt5*bn1*rhoinv(ji)/c*(un+c)*rho(ji)*c*bn1*bn2 + us*bs2 -	&
                    fact1 * g2mu(ji)
	B(ji,2,3) = pt5*bn1*rhoinv(ji)/c*(un+c)*rho(ji)*c*bn2*bn2 -		&
                    damp(ji) * fact1 * g1lm(ji)
	B(ji,2,4) = zero
	B(ji,2,5) = bn1/(two*rho(ji)*c)*(un+c)*bn2*rho(ji)*gmsinv +		&
		    bs1*gmsinv*bs2 -						&
                    fact1 * dmu(ji) * ( gu(ji,1,2) + damp(ji) * gu(ji,2,1) )

!=======================================================================================================!

!.... Momentum equation -- x_2

	B(ji,3,1) = pt5*bn2*rhoinv(ji)/c*(un+c)*bn2*t(ji)*gmsinv +	&
                    bs2 * rhoinv(ji) * gmsinv * t(ji) * bs2
	B(ji,3,2) = pt5*bn2*rhoinv(ji)/c*(un+c)*rho(ji)*c*bn1*bn2 -	&
                    damp(ji) * fact1 * g1mu(ji)
	B(ji,3,3) = pt5*bn2*rhoinv(ji)/c*(un+c)*rho(ji)*c*bn2*bn2 +	&
                    us*bs2 - fact1 * g2lm(ji) - fact1 * two * g2mu(ji)
	B(ji,3,4) = -fact1 * g3mu(ji)
	B(ji,3,5) = bn2/(two*rho(ji)*c)*(un+c)*bn2*rho(ji)*gmsinv +	&
                    bs2*gmsinv*bs2 - fact1 * dlm(ji) *			&
	            (damp(ji) * gu(ji,1,1) + gu(ji,2,2) + gu(ji,3,3)) - &
                    fact1 * dmu(ji) * ( gu(ji,2,2) + gu(ji,2,2) )

!=======================================================================================================!

!.... Momentum equation -- x_3

	B(ji,4,1) = zero
	B(ji,4,2) = zero
	B(ji,4,3) = -fact1 * g3lm(ji)
	B(ji,4,4) = bs2*us - fact1 * g2mu(ji)
	B(ji,4,5) = -fact1 * dmu(ji) * ( gu(ji,3,2) + gu(ji,2,3) )

!=======================================================================================================!

!.... Energy equation

	fact1 = gamma * rhoinv(ji) / (Pr * Re)
	fact2 = gamma * gamma1 * Ma**2 * rhoinv(ji) / Re

	B(ji,5,1) =  Ma**2*pt5*rhoinv(ji)*gamma1*(un+c)*bn2*gmsinv*t(ji)

	B(ji,5,2) =  Ma**2*pt5*rhoinv(ji)*gamma1*(un+c)*rho(ji)*c*bn1*bn2 +	&
                     gamma1*t(ji)*bs1*bs2 -					&
                     two*fact2*mu(ji)*( gu(ji,1,2) + damp(ji) * gu(ji,2,1) )

	B(ji,5,3) =  Ma**2*pt5*rhoinv(ji)*gamma1*(un+c)*rho(ji)*c*bn2*bn2 +	&
		     gamma1*t(ji)*bs2*bs2 - fact2 * two * lm(ji) *		&
	             (damp(ji) * gu(ji,1,1) + gu(ji,2,2) + gu(ji,3,3)) -	&
	             two * fact2 * mu(ji) * ( gu(ji,2,2) + gu(ji,2,2) )

	B(ji,5,4) = -two * fact2 * mu(ji) * ( gu(ji,3,2) + gu(ji,2,3) )

	B(ji,5,5) =  Ma**2*pt5*rhoinv(ji)*gamma1*(un+c)*bn2*rho(ji)*gmsinv +	&
		     us*bs2 - fact1 * (g2con(ji) + dcon(ji) * gt(ji,2))


	else

!=======================================================================================================!
!	O u t f l o w
!=======================================================================================================!

!.... Continuity equation

	gmsinv = one / ( gamma * Ma**2 )

	B(ji,1,1) = one/c**2*( un*(c**2*bn2-bn2*t(ji)*gmsinv) +		&
                    pt5*(un+c)*bn2*gmsinv*t(ji) ) + us*bs2
	B(ji,1,2) = one/c**2*pt5*(un+c)*rho(ji)*c*bn1*bn2 + rho(ji)*bs1*bs2
	B(ji,1,3) = one/c**2*pt5*(un+c)*rho(ji)*c*bn2*bn2 + rho(ji)*bs2*bs2
	B(ji,1,4) = zero
	B(ji,1,5) = one/c**2*( un*(-bn2*rho(ji)*gmsinv) +		&
                    pt5*(un+c)*bn2*rho(ji)*gmsinv )

!=======================================================================================================!

!.... Momentum equation -- x_1

	fact1  = rhoinv(ji) / Re

	B(ji,2,1) = pt5*bn1*rhoinv(ji)/c*(un+c)*bn2*t(ji)*gmsinv +		&
                    bs1*rhoinv(ji)*gmsinv*t(ji)*bs2
	B(ji,2,2) = pt5*bn1*rhoinv(ji)/c*(un+c)*rho(ji)*c*bn1*bn2 +		&
		    bs1*un*bs1*bn2 + us*bs2 -					&
                    fact1 * g2mu(ji)
	B(ji,2,3) = pt5*bn1*rhoinv(ji)/c*(un+c)*rho(ji)*c*bn2*bn2 +		&
		    bs1*un*bs2*bn2 -						&
                    damp(ji) * fact1 * g1lm(ji)
	B(ji,2,4) = zero
	B(ji,2,5) = bn1/(two*rho(ji)*c)*(un+c)*bn2*rho(ji)*gmsinv +		&
		    bs1*gmsinv*bs2 -						&
                    fact1 * dmu(ji) * ( gu(ji,1,2) + damp(ji) * gu(ji,2,1) )

!=======================================================================================================!

!.... Momentum equation -- x_2

	B(ji,3,1) = pt5*bn2*rhoinv(ji)/c*(un+c)*bn2*t(ji)*gmsinv +	&
                    bs2 * rhoinv(ji) * gmsinv * t(ji) * bs2
	B(ji,3,2) = pt5*bn2*rhoinv(ji)/c*(un+c)*rho(ji)*c*bn1*bn2 +	&
		    bs2*un*bs1*bn2 -					&
                    damp(ji) * fact1 * g1mu(ji)
	B(ji,3,3) = pt5*bn2*rhoinv(ji)/c*(un+c)*rho(ji)*c*bn2*bn2 +	&
                    us*bs2 + bs2*un*bs2*bn2 -				&
		    fact1 * g2lm(ji) - fact1 * two * g2mu(ji)
	B(ji,3,4) = -fact1 * g3mu(ji)
	B(ji,3,5) = bn2/(two*rho(ji)*c)*(un+c)*bn2*rho(ji)*gmsinv +	&
                    bs2*gmsinv*bs2 - fact1 * dlm(ji) *			&
	            (damp(ji) * gu(ji,1,1) + gu(ji,2,2) + gu(ji,3,3)) - &
                    fact1 * dmu(ji) * ( gu(ji,2,2) + gu(ji,2,2) )

!=======================================================================================================!

!.... Momentum equation -- x_3

	B(ji,4,1) = zero
	B(ji,4,2) = zero
	B(ji,4,3) = -fact1 * g3lm(ji)
	B(ji,4,4) = bn2*un + bs2*us - fact1 * g2mu(ji)
	B(ji,4,5) = -fact1 * dmu(ji) * ( gu(ji,3,2) + gu(ji,2,3) )

!=======================================================================================================!

!.... Energy equation

	fact1 = gamma * rhoinv(ji) / (Pr * Re)
	fact2 = gamma * gamma1 * Ma**2 * rhoinv(ji) / Re

	B(ji,5,1) = -Ma**2*rhoinv(ji)*( un*(c**2*bn2-bn2*t(ji)*gmsinv) -	&
		     pt5*gamma1*(un+c)*bn2*t(ji)*gmsinv ) 

	B(ji,5,2) =  Ma**2*rhoinv(ji)*pt5*gamma1*(un+c)*rho(ji)*c*bn1*bn2 +	&
                     gamma1*t(ji)*bs1*bs2 -					&
                     two*fact2*mu(ji)*( gu(ji,1,2) + damp(ji) * gu(ji,2,1) )

	B(ji,5,3) =  Ma**2*rhoinv(ji)*pt5*gamma1*(un+c)*rho(ji)*c*bn2*bn2 +	&
		     gamma1*t(ji)*bs2*bs2 - fact2 * two * lm(ji) *		&
	             (damp(ji) * gu(ji,1,1) + gu(ji,2,2) + gu(ji,3,3)) -	&
	             two * fact2 * mu(ji) * ( gu(ji,2,2) + gu(ji,2,2) )

	B(ji,5,4) = -two * fact2 * mu(ji) * ( gu(ji,3,2) + gu(ji,2,3) )

	B(ji,5,5) = -Ma**2*rhoinv(ji)*( un*(-bn2*rho(ji)*gmsinv) -		&
                     pt5*gamma1*(un+c)*bn2*rho(ji)*gmsinv ) +			&
		     us*bs2 - fact1 * (g2con(ji) + dcon(ji) * gt(ji,2))

	end if		! inflow/outflow

	end do

	end if 		! top
	
!=======================================================================================================!
!	W a l l   B o u n d a r y
!=======================================================================================================!
	if (wall.eq.4) then

	j = 1
	do i = 1, nx
	  ji = j + (i-1) * ny
    
!.... get the boundary normal and tangent vectors

	    bn1 = n1(ji) / sqrt( n1(ji)**2 + n2(ji)**2 )
	    bn2 = n2(ji) / sqrt( n1(ji)**2 + n2(ji)**2 )

	    bs1 = -bn2
	    bs2 =  bn1

!.... compute Us and Un
	
	    un = bn1 * vl(ji,2) + bn2 * vl(ji,3)
	    us = bs1 * vl(ji,2) + bs2 * vl(ji,3)

	    c = sqrt( t(ji) ) / Ma
	    dc = pt5 / (Ma**2 * c)

!.... compute some derivatives in the boundary normal coordinate system

	    gpn = bn1 * gp(ji,1) + bn2 * gp(ji,2)
	    gps = bs1 * gp(ji,1) + bs2 * gp(ji,2)
	
	    grhon = bn1 * grho(ji,1) + bn2 * grho(ji,2)
	    grhos = bs1 * grho(ji,1) + bs2 * grho(ji,2)

	    gtn = bn1 * gt(ji,1) + bn2 * gt(ji,2)
	    gts = bs1 * gt(ji,1) + bs2 * gt(ji,2)

	    gu1n = bn1 * gu(ji,1,1) + bn2 * gu(ji,1,2)
	    gu1s = bs1 * gu(ji,1,1) + bs2 * gu(ji,1,2)

	    gu2n = bn1 * gu(ji,2,1) + bn2 * gu(ji,2,2)
	    gu2s = bs1 * gu(ji,2,1) + bs2 * gu(ji,2,2)

	    gu3n = bn1 * gu(ji,3,1) + bn2 * gu(ji,3,2)
	    gu3s = bs1 * gu(ji,3,1) + bs2 * gu(ji,3,2)

	    gunn = bn1*bn1 * gu(ji,1,1) + bn2*bn1 * gu(ji,2,1) + &
	           bn1*bn2 * gu(ji,1,2) + bn2*bn2 * gu(ji,2,2)

	    gusn = bs1*bn1 * gu(ji,1,1) + bs2*bn1 * gu(ji,2,1) + &
	           bs1*bn2 * gu(ji,1,2) + bs2*bn2 * gu(ji,2,2)

	    guns = bn1*bs1 * gu(ji,1,1) + bn2*bs1 * gu(ji,2,1) + &
	           bn1*bs2 * gu(ji,1,2) + bn2*bs2 * gu(ji,2,2)

	    guss = bs1*bs1 * gu(ji,1,1) + bs2*bs1 * gu(ji,2,1) + &
	           bs1*bs2 * gu(ji,1,2) + bs2*bs2 * gu(ji,2,2)
	    
!.... Continuity equation

	    gmsinv = one / ( gamma * Ma**2 )
    
	    B(ji,1,1) = one/c**2*(un-c)*bn2*t(ji)*gmsinv + us*bs2
	    B(ji,1,2) = -one/c**2*(un-c)*rho(ji)*c*bn1*bn2 + rho(ji)*bs1*bs2
	    B(ji,1,3) = -one/c**2*(un-c)*rho(ji)*c*bn2*bn2 + rho(ji)*bs2*bs2
	    B(ji,1,4) = zero
	    B(ji,1,5) = one/c**2*(un-c)*bn2*rho(ji)*gmsinv

	end do

	end if 		! wall

	return
	end
