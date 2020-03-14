!=============================================================================!
	subroutine rk3d(vl, rl, vml, xl, yl, dtl)
!  
!  Computes the RK3 RHS for the linearized three-dimensional compressible 
!  Navier-Stokes solver.  Updated to use the hat matrices
!
!  NOTE: the sign is opposite from lrhs!
!
!  Revised: 4-22-96
!=============================================================================!
	use global
	use local
	use local3d
	implicit none
!=============================================================================!
	complex :: vl(ny*nx,ndof), rl(ny*nx,ndof)
	real    :: vml(ny*nx,ndof)
	real    :: xl(ny*nx), yl(ny*nx), dtl(ny*nx)

!.... local variables
	
	integer :: i, j, ij, idof, jdof
	
!.... mean flow variables

	real :: cm(ny*nx), um(ny*nx)
	real :: a, d, kk
	complex :: c3(ny*nx)

	character*80 name
	integer :: lrec, ier, istat
!=============================================================================!

!.... compute first derivatives
	
	call cgrad( ndof, nx, ny, vl, c1v, c2v, dxi, deta, optx, opty, &
	            xper, yper, lsym, rsym, bsym, tsym, .true.)
	
!.... compute second derivatives

	call cgrad2( ndof, nx, ny, vl, c1v, c11v, c12v, c22v, dxi, deta, &
	             optx, opty, xper, yper, lsym, rsym, bsym, tsym, .true.)

!.... enforce adiabatic wall temperature boundary condition

	call cgradbc( c1v, c2v, c11v, c12v, c22v )

!=============================================================================!
!.... form the RHS
!=============================================================================!

!.... U,xi term

	rl(:,1) =         - Ah(:,1,1) * c1v(:,1)	&
	                  - Ah(:,1,2) * c1v(:,2)	&
	                  - Ah(:,1,3) * c1v(:,3)	&
	                  - Ah(:,1,4) * c1v(:,4)	&
	                  - Ah(:,1,5) * c1v(:,5)	

	rl(:,2) =         - Ah(:,2,1) * c1v(:,1)	&
	                  - Ah(:,2,2) * c1v(:,2)	&
	                  - Ah(:,2,3) * c1v(:,3)	&
	                  - Ah(:,2,4) * c1v(:,4)	&
	                  - Ah(:,2,5) * c1v(:,5)

	rl(:,3) =         - Ah(:,3,1) * c1v(:,1)	&
	                  - Ah(:,3,2) * c1v(:,2)	&
	                  - Ah(:,3,3) * c1v(:,3)	&
	                  - Ah(:,3,4) * c1v(:,4)	&
	                  - Ah(:,3,5) * c1v(:,5)

	rl(:,4) =         - Ah(:,4,1) * c1v(:,1)	&
	                  - Ah(:,4,2) * c1v(:,2)	&
	                  - Ah(:,4,3) * c1v(:,3)	&
	                  - Ah(:,4,4) * c1v(:,4)	&
	                  - Ah(:,4,5) * c1v(:,5)

	rl(:,5) =         - Ah(:,5,1) * c1v(:,1)	&
	                  - Ah(:,5,2) * c1v(:,2)	&
	                  - Ah(:,5,3) * c1v(:,3)	&
	                  - Ah(:,5,4) * c1v(:,4)	&
	                  - Ah(:,5,5) * c1v(:,5)

!.... U,eta term

	rl(:,1) = rl(:,1) - Bh(:,1,1) * c2v(:,1)	&
	                  - Bh(:,1,2) * c2v(:,2)	&
	                  - Bh(:,1,3) * c2v(:,3)	&
	                  - Bh(:,1,4) * c2v(:,4)	&
	                  - Bh(:,1,5) * c2v(:,5)	

	rl(:,2) = rl(:,2) - Bh(:,2,1) * c2v(:,1)	&
	                  - Bh(:,2,2) * c2v(:,2)	&
	                  - Bh(:,2,3) * c2v(:,3)	&
	                  - Bh(:,2,4) * c2v(:,4)	&
	                  - Bh(:,2,5) * c2v(:,5)

	rl(:,3) = rl(:,3) - Bh(:,3,1) * c2v(:,1)	&
	                  - Bh(:,3,2) * c2v(:,2)	&
	                  - Bh(:,3,3) * c2v(:,3)	&
	                  - Bh(:,3,4) * c2v(:,4)	&
	                  - Bh(:,3,5) * c2v(:,5)

	rl(:,4) = rl(:,4) - Bh(:,4,1) * c2v(:,1)	&
	                  - Bh(:,4,2) * c2v(:,2)	&
	                  - Bh(:,4,3) * c2v(:,3)	&
	                  - Bh(:,4,4) * c2v(:,4)	&
	                  - Bh(:,4,5) * c2v(:,5)

	rl(:,5) = rl(:,5) - Bh(:,5,1) * c2v(:,1)	&
	                  - Bh(:,5,2) * c2v(:,2)	&
	                  - Bh(:,5,3) * c2v(:,3)	&
	                  - Bh(:,5,4) * c2v(:,4)	&
	                  - Bh(:,5,5) * c2v(:,5)

!.... U term

	rl(:,1) = rl(:,1) - Dh(:,1,1) * vl(:,1)	&
	                  - Dh(:,1,2) * vl(:,2)	&
	                  - Dh(:,1,3) * vl(:,3)	&
	                  - Dh(:,1,4) * vl(:,4)	&
	                  - Dh(:,1,5) * vl(:,5)

	rl(:,2) = rl(:,2) - Dh(:,2,1) * vl(:,1)	&
	                  - Dh(:,2,2) * vl(:,2)	&
	                  - Dh(:,2,3) * vl(:,3)	&
	                  - Dh(:,2,4) * vl(:,4)	&
	                  - Dh(:,2,5) * vl(:,5)

	rl(:,3) = rl(:,3) - Dh(:,3,1) * vl(:,1)	&
			  - Dh(:,3,2) * vl(:,2)	&
			  - Dh(:,3,3) * vl(:,3)	&
	                  - Dh(:,3,4) * vl(:,4)	&
	                  - Dh(:,3,5) * vl(:,5)

	rl(:,4) = rl(:,4) - Dh(:,4,1) * vl(:,1)	&
			  - Dh(:,4,2) * vl(:,2)	&
			  - Dh(:,4,3) * vl(:,3)	&
	                  - Dh(:,4,4) * vl(:,4)	&
	                  - Dh(:,4,5) * vl(:,5)

	rl(:,5) = rl(:,5) - Dh(:,5,1) * vl(:,1)	&
	                  - Dh(:,5,2) * vl(:,2)	&
	                  - Dh(:,5,3) * vl(:,3)	&
	                  - Dh(:,5,4) * vl(:,4)	&
	                  - Dh(:,5,5) * vl(:,5)

!.... U,\xi\xi term

	rl(:,2) = rl(:,2) + Vh11(:,1) * c11v(:,2) 	&
	                  + Vh11(:,5) * c11v(:,3)
	rl(:,3) = rl(:,3) + Vh11(:,6) * c11v(:,2)	&
			  + Vh11(:,2) * c11v(:,3)
	rl(:,4) = rl(:,4) + Vh11(:,3) * c11v(:,4)
	rl(:,5) = rl(:,5) + Vh11(:,4) * c11v(:,5)

!.... U,\xi\eta term

	rl(:,2) = rl(:,2) + Vh12(:,1) * c12v(:,2) 	&
	                  + Vh12(:,5) * c12v(:,3)
	rl(:,3) = rl(:,3) + Vh12(:,6) * c12v(:,2)	&
			  + Vh12(:,2) * c12v(:,3)
	rl(:,4) = rl(:,4) + Vh12(:,3) * c12v(:,4)
	rl(:,5) = rl(:,5) + Vh12(:,4) * c12v(:,5)

!.... U,\eta\eta term

	rl(:,2) = rl(:,2) + Vh22(:,1) * c22v(:,2) 	&
	                  + Vh22(:,5) * c22v(:,3)
	rl(:,3) = rl(:,3) + Vh22(:,6) * c22v(:,2)	&
			  + Vh22(:,2) * c22v(:,3)
	rl(:,4) = rl(:,4) + Vh22(:,3) * c22v(:,4)
	rl(:,5) = rl(:,5) + Vh22(:,4) * c22v(:,5)

!.... 3D terms

!.... \hat{A}_i, \hat{B}_i terms

	rl(:,2) = rl(:,2) + im * ( ABhi(:,1) * c1v(:,4) + ABhi(:,3) * c2v(:,4) )
	rl(:,3) = rl(:,3) + im * ( ABhi(:,2) * c1v(:,4) + ABhi(:,4) * c2v(:,4) )
	rl(:,4) = rl(:,4) + im * ( ABhi(:,1) * c1v(:,2) + ABhi(:,2) * c1v(:,3) &
	                         + ABhi(:,3) * c2v(:,2) + ABhi(:,4) * c2v(:,3) )
	
!.... \hat{D}_i term

	rl(:,1) = rl(:,1) - im * ( Dhi(:,1,1) * vl(:,1)	&
	                         + Dhi(:,1,2) * vl(:,2)	&
	                         + Dhi(:,1,3) * vl(:,3)	&
	                         + Dhi(:,1,4) * vl(:,4)	&
	                         + Dhi(:,1,5) * vl(:,5) )

	rl(:,2) = rl(:,2) - im * ( Dhi(:,2,1) * vl(:,1)	&
	                         + Dhi(:,2,2) * vl(:,2)	&
	                         + Dhi(:,2,3) * vl(:,3)	&
	                         + Dhi(:,2,4) * vl(:,4)	&
	                         + Dhi(:,2,5) * vl(:,5) )

	rl(:,3) = rl(:,3) - im * ( Dhi(:,3,1) * vl(:,1)	&
			         + Dhi(:,3,2) * vl(:,2)	&
			         + Dhi(:,3,3) * vl(:,3)	&
	                         + Dhi(:,3,4) * vl(:,4)	&
	                         + Dhi(:,3,5) * vl(:,5) )

	rl(:,4) = rl(:,4) - im * ( Dhi(:,4,1) * vl(:,1)	&
			         + Dhi(:,4,2) * vl(:,2)	&
			         + Dhi(:,4,3) * vl(:,3)	&
	                         + Dhi(:,4,4) * vl(:,4)	&
	                         + Dhi(:,4,5) * vl(:,5) )

	rl(:,5) = rl(:,5) - im * ( Dhi(:,5,1) * vl(:,1)	&
	                         + Dhi(:,5,2) * vl(:,2)	&
	                         + Dhi(:,5,3) * vl(:,3)	&
	                         + Dhi(:,5,4) * vl(:,4)	&
	                         + Dhi(:,5,5) * vl(:,5) )

!.... Sponge term

	if (ispg .eq. 1) then		!.... standard sponge

	  rl(:,1) = rl(:,1) - spg(:) * vl(:,1)
	  rl(:,2) = rl(:,2) - spg(:) * vl(:,2)
	  rl(:,3) = rl(:,3) - spg(:) * vl(:,3)
	  rl(:,4) = rl(:,4) - spg(:) * vl(:,4)
	  rl(:,5) = rl(:,5) - spg(:) * vl(:,5)
	
	else if (ispg .eq. 2) then		
	
	  rl(:,1) = rl(:,1) - (spg(:) + spg2(:)) * vl(:,1)
	  rl(:,2) = rl(:,2) - (spg(:) + spg2(:)) * vl(:,2)
	  rl(:,3) = rl(:,3) - (spg(:) + spg2(:)) * vl(:,3)
	  rl(:,4) = rl(:,4) - (spg(:) + spg2(:)) * vl(:,4)
	  rl(:,5) = rl(:,5) - (spg(:) + spg2(:)) * vl(:,5)

	else if (ispg .eq. 3) then	!.... double sponge

!.... now do the outflow sponge

	  rl(:,1) = rl(:,1) - spg2(:) * vl(:,1)
	  rl(:,2) = rl(:,2) - spg2(:) * vl(:,2)
	  rl(:,3) = rl(:,3) - spg2(:) * vl(:,3)
	  rl(:,4) = rl(:,4) - spg2(:) * vl(:,4)
	  rl(:,5) = rl(:,5) - spg2(:) * vl(:,5)	

!.... the inflow sponge

	  cm = sqrt( vml(:,5) ) / Ma
	  um = vml(:,2)

	  d  = pt5 * ( onept33 + gamma1 / Pr ) / Re

	  do i = 1, nx
	    do j = 1, ny
	      ij = j + (i-1)*ny
	      kk = omega / (cm(ij)+um(ij))
	      a  = omega**2 * d / (cm(ij)+um(ij))**3
	      c3(ij) = exp( -a * (xl(ij) - x0) ) * exp( im * kk * xl(ij) )
!	      c3(ij) = wamp(i) * exp( im * kk * xl(ij) )
	    end do
	  end do

	  do ij = 1, nx*ny
	    rl(ij,1) = rl(ij,1) - spg(ij) * ( vl(ij,1) - &
			pt5 * c3(ij) / cm(ij)**2 )
	    rl(ij,2) = rl(ij,2) - spg(ij) * ( vl(ij,2) - &
			c3(ij) * pt5 / ( vml(ij,1) * cm(ij) ) )
	    rl(ij,3) = rl(ij,3) - spg(ij) * ( vl(ij,3) )
	    rl(ij,4) = rl(ij,4) - spg(ij) * ( vl(ij,4) )
	    rl(ij,5) = rl(ij,5) - spg(ij) * ( vl(ij,5) - &
			(gamma*Ma**2 * c3(ij) * pt5 - &
			vml(ij,5) * pt5 * c3(ij) / cm(ij)**2) / vml(ij,1) )
	  end do

	else if (ispg .eq. 4) then		

	  call espg_it( rl, vl, spg, spg2 )

	end if

!.... explicit smoother

	if (eps_e .ne. zero) call smoother3D( rl, vl )

!=============================================================================!

!.... multiply by the time step

	do idof = 1, ndof
	  rl(:,idof) = dtl * rl(:,idof)
	end do
	
	call rkbc3D(rl)

!.... output RHS statistics

	if (mod(lstep,itout).eq.0 .and. iter.eq.4) then
	  call resstat3D(rl, dtl)
	end if
	
	return
	end
!=============================================================================!
	module eic

	integer :: ic_start=0
	complex, allocatable :: vic(:,:,:)
	
	end module eic

!=============================================================================!
        subroutine espg_it( rl, vl, spgl, spg2l )
!
!       This routine reads in the initial condition when first called, and
!       sponges to the initial condition.
!
!=============================================================================!
	use eic
	use global
	use local
	implicit none
	
	complex :: rl(ny,nx,ndof), vl(ny,nx,ndof)
	real    :: spgl(ny,nx), spg2l(ny,nx)
	integer :: i

	real :: rtmp
	integer :: itmp
!=============================================================================!
	if (ic_start .eq. 0) then
	  allocate( vic(ny,nx,ndof) )
	  ic_start = 1
	  open(10,file='output.R.0',form='unformatted',status='old')
	  read(10) itmp, rtmp, itmp, itmp, itmp, itmp, &
	           rtmp, rtmp, rtmp, rtmp, rtmp
	  read(10) vic
	  close(10)
	end if

	do i = 1, nx
	  rl(:,i,1) = rl(:,i,1) - (spgl(:,i) + spg2l(:,i)) * &
	              ( vl(:,i,1) - vic(:,i,1) )
	  rl(:,i,2) = rl(:,i,2) - (spgl(:,i) + spg2l(:,i)) * &
	              ( vl(:,i,2) - vic(:,i,2) )
	  rl(:,i,3) = rl(:,i,3) - (spgl(:,i) + spg2l(:,i)) * &
	              ( vl(:,i,3) - vic(:,i,3) )
	  rl(:,i,4) = rl(:,i,4) - (spgl(:,i) + spg2l(:,i)) * &
	              ( vl(:,i,4) - vic(:,i,4) )
	  rl(:,i,5) = rl(:,i,5) - (spgl(:,i) + spg2l(:,i)) * &
	              ( vl(:,i,5) - vic(:,i,5) )
	end do
	
	return
        end
