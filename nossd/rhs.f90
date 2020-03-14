!=============================================================================!
	subroutine rhs(rl, vl) 
!  
!       Generates the RHS for the compressible 2D-3C, N-S equations 
!
!	Written: 5-10-95
!       Revised: 6-19-95
!
!=============================================================================!
	use global
	use local
	implicit none

	real :: rl(ny*nx,ndof), vl(ny*nx,ndof)
	real :: fact1, fact2
!=============================================================================!

!.... continuity

	rl(:,1) = rho(:) * divu(:) + grho(:,1) * u1 + grho(:,2) * u2

!.... momentum

	fact1 = one / Re

	rl(:,2) = u1 * gu(:,1,1) + u2 * gu(:,1,2) +			&
	          rhoinv * gp(:,1) -					&
	          fact1 * ( g1lm * divu + lm * g1divu +			&
		            two * (g1mu * S(:,1,1) + g2mu * S(:,1,2) +	&
			           mu * S1jj ) ) * rhoinv

	rl(:,3) = u1 * gu(:,2,1) + u2 * gu(:,2,2) +			&
	          rhoinv * gp(:,2) -					&
	          fact1 * ( g2lm * divu + lm * g2divu +			&
		            two * (g1mu * S(:,2,1) + g2mu * S(:,2,2) + &
			           mu * S2jj ) ) * rhoinv

	rl(:,4) = u1 * gu(:,3,1) + u2 * gu(:,3,2) -			&
	          fact1 * ( two * (g1mu * S(:,3,1) + g2mu * S(:,3,2) + &
			           mu * S3jj ) ) * rhoinv

!.... temperature

	fact1 = gamma / ( Pr * Re )
	fact2 = gamma * gamma1 * Ma**2 / Re

	rl(:,5) = u1 * gt(:,1) + u2 * gt(:,2) +				&
                  gamma1 * t * divu - rhoinv * (			&
		  fact1 * ( g1con * gt(:,1) + g2con * gt(:,2) +		&
	                    con * Lapt ) +				&
		  fact2 * ( lm * divu**2 + two * mu * (			&
			S(:,1,1)**2 + S(:,1,2)**2 + S(:,1,3)**2 +	&
                        S(:,2,1)**2 + S(:,2,2)**2 + S(:,2,3)**2 +	&
		        S(:,3,1)**2 + S(:,3,2)**2 ) ) )

!.... correct the outflow boundaries using Lele & Poinsot BC's

	call rhs_l(rl,vl)
	
!.... standard sponge

	if (ispg .gt. 0) then
	  if (ispg .eq. 1) then
	    call spg_it(rl,vl,spg)
	  else if (ispg .eq. 2) then
	    rl(:,1) = rl(:,1) + (spg(:) + spg2(:)) * ( vl(:,1) - one )
	    rl(:,2) = rl(:,2) + (spg(:) + spg2(:)) * ( vl(:,2) )
	    rl(:,3) = rl(:,3) + (spg(:) + spg2(:)) * ( vl(:,3) )
	    rl(:,4) = rl(:,4) + (spg(:) + spg2(:)) * ( vl(:,4) )
	    rl(:,5) = rl(:,5) + (spg(:) + spg2(:)) * ( vl(:,5) - one )
	  else
	    call error('rhs_p$','ispg > 2 is not supported$')
	  end if
	end if

!.... explicit smoother

	if (eps_e .ne. zero) call smoother( rl, vl )
	
	return
	end
