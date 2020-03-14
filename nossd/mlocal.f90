!=============================================================================!
	subroutine mlocal(opt)
!
!  allocate space for local variables
!
!=============================================================================!
	use global
	use local
	implicit none
	
	integer :: ier, mem

	character*10 :: opt
	character*80 :: code='mLocal$'
!=============================================================================!
	mem = 0

	if (opt .eq. 'allocate  ') then

!.... allocate memory for mean gradients

	mem = mem + ny*nx*ndof*5
	allocate( g1v(ny*nx,ndof),  g2v(ny*nx,ndof),  		&
	          g11v(ny*nx,ndof), g12v(ny*nx,ndof), 		&
		  g22v(ny*nx,ndof), STAT=ier )
	if (ier .ne. 0) call error(code,'Insufficient Memory for gradients$')

	mem = mem + ny*nx*5
	allocate( g1vl(ny*nx),      g2vl(ny*nx),		&
		  g11vl(ny*nx),     g12vl(ny*nx),		&
		  g22vl(ny*nx),	    STAT=ier )
	if (ier .ne. 0) call error(code, 			&
		'Insufficient Memory for local gradients$')

!.... stuff above this line is needed for linear and nonlinear analysis

	mem = mem + ny*nx*6 + ny*nx*nsd*(nsd+2)
	allocate( rho(ny*nx),   u1(ny*nx),    u2(ny*nx),	&
		  u3(ny*nx),    t(ny*nx),     p(ny*nx),		&
		  grho(ny*nx,nsd),    gu(ny*nx,nsd,nsd),	&
	          gt(ny*nx,nsd),      rhoinv(ny*nx),   STAT=ier )
	if (ier .ne. 0) call error(code,			&
	        'Insufficient Memory for local variables$')

	mem = mem + ny*nx*7
	allocate( g1divu(ny*nx), g2divu(ny*nx), g3divu(ny*nx), 	&
	          S1jj(ny*nx),   S2jj(ny*nx),   S3jj(ny*nx), 	&
		  Lapt(ny*nx),   STAT=ier )
	if (ier .ne. 0) call error(code,			&
	        'Insufficient Memory for local variables$')

	mem = mem + ny*nx*(1+nsd+nsd*nsd)
	allocate( divu(ny*nx), gp(ny*nx,nsd), S(ny*nx,nsd,nsd), STAT=ier )
	if (ier .ne. 0) call error(code,'Insufficient Memory$')

	mem = mem + ny*nx*27
	allocate( mu(ny*nx),     dmu(ny*nx),    d2mu(ny*nx),	&
		  lm(ny*nx),     dlm(ny*nx),    d2lm(ny*nx),	&
		  con(ny*nx),    dcon(ny*nx),   d2con(ny*nx),	&
		  g1mu(ny*nx),   g2mu(ny*nx),   g3mu(ny*nx),	&
		  g1lm(ny*nx),   g2lm(ny*nx),   g3lm(ny*nx),	&
		  g1con(ny*nx),  g2con(ny*nx),  g3con(ny*nx),	&
		  g1dmu(ny*nx),  g2dmu(ny*nx),  g3dmu(ny*nx),	&
		  g1dlm(ny*nx),  g2dlm(ny*nx),  g3dlm(ny*nx),	&
		  g1dcon(ny*nx), g2dcon(ny*nx), g3dcon(ny*nx),	&
		  STAT=ier )
	if (ier .ne. 0) call error(code,'Insufficient Memory$')

!.... stuff below this line is only needed for linear problems or 
!.... implicit problems

	if (linear.eq.1 .or. impl.ne.0) then
	  mem = mem + ny*nx*ndof*ndof*3
	  allocate( Ah(ny*nx,ndof,ndof), Bh(ny*nx,ndof,ndof), 	&
		    Dh(ny*nx,ndof,ndof), STAT=ier )
	  if (ier .ne. 0) call error(code,'Insufficient Memory$')
  
	  mem = mem + ny*nx*6*3
	  allocate( Vh11(ny*nx,6), Vh12(ny*nx,6), Vh22(ny*nx,6), STAT=ier )
	  if (ier .ne. 0) call error(code,'Insufficient Memory$')

!.... 3-d stuff

	  mem = mem + ny*nx*ndof*ndof
	  allocate( Dhi(ny*nx,ndof,ndof), STAT=ier )
	  if (ier .ne. 0) call error(code,'Insufficient Memory$')

	  mem = mem + ny*nx*6
	  allocate( ABhi(ny*nx,6), STAT=ier )
	  if (ier .ne. 0) call error(code,'Insufficient Memory$')

	end if

	write(*,"(' Mlocal allocated  ===> ',1pe13.6,' words')") float(mem)

	else if (opt .eq. 'deallocate') then

!.... deallocate stuff not needed for linear analysis

	  mem = ny*nx*7 + ny*nx*nsd*(nsd+2)
	  deallocate(rho, u1, u2, u3, t, p, rhoinv)
	  deallocate(gu, grho, gt)
	  mem = mem + ny*nx*7
	  deallocate(g1divu, g2divu, g3divu, S1jj, S2jj, S3jj, Lapt)
	  mem = mem + ny*nx*(1+nsd+nsd*nsd)
	  deallocate(divu, gp, S)
	  mem = mem + ny*nx*27
	  deallocate(mu, dmu, d2mu, 		&
	  	     lm, dlm, d2lm, 		&
		     con, dcon, d2con, 		&
	             g1mu, g2mu, g3mu, 		&
		     g1lm, g2lm, g3lm, 		&
		     g1con, g2con, g3con, 	&
		     g1dmu, g2dmu, g3dmu, 	&
		     g1dlm, g2dlm, g3dlm, 	&
		     g1dcon, g2dcon, g3dcon )
	  
	  write(*,"(' Mlocal deallocated ===> ',1pe13.6,' words')") float(mem)

	else

	  call error(code,'Illegal option$')
	  
	end if

!=============================================================================!
	return
	end
