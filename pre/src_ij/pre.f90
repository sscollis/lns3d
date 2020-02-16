!=============================================================================!
	program pre
!
!  This is a pre processor for the LNS code.  It takes a plot3d grid 
!  file as input and computes the metric file needed for LNS.
!
!  Revised:  10-12-00    Switched order of indices, removed B-splines
!
!=============================================================================!
	use const
	implicit none
	
	integer :: i, j, k, ier, ians
	integer :: nx, ny, nz
	
	real, allocatable :: xy(:,:,:),  g1xy(:,:,:), g2xy(:,:,:)
	real, allocatable :: g11xy(:,:,:), g12xy(:,:,:), g22xy(:,:,:)
	real, allocatable :: g1m(:,:,:), g2m(:,:,:)
	real, allocatable :: g1n(:,:,:), g2n(:,:,:)
	real, allocatable :: xi(:),      eta(:),      IdetJ(:,:)
	real, allocatable :: m(:,:,:),   m11(:,:),    m12(:,:),     m22(:,:)
	real, allocatable :: n(:,:,:),   n11(:,:),    n12(:,:),     n22(:,:)
	real, allocatable :: q(:,:,:)
	
	real, allocatable :: f(:,:)
	real, allocatable :: g1(:,:), g2(:,:), g1l(:,:), g2l(:,:)
	real, allocatable :: g11(:,:),  g12(:,:),  g22(:,:)
	real, allocatable :: g11l(:,:), g12l(:,:), g22l(:,:)
	
	real :: dxi, deta, tmp, kappa

	character*80 :: name, filen
	
        integer :: optx = -1, opty = -1
	logical :: xper = .false., yper = .false.
        logical :: lsym = .false., rsym = .false.
        logical :: tsym = .false., bsym = .false.
        logical :: carp = .false.
	
!=============================================================================!
	write(*,"('Hardwire periodicity (1,0) ==> ',$)") 
	read(*,*) ians
	if (ians.eq.1) then
	  xper = .true.
	else
	  xper = .false.
	end if

!.... read the physical grid file

 10     filen = 'grid.dat'
	write (*,"('Enter unformatted grid file name [',a,']? ',$)") &
     		    filen(1:index(filen,' '))
	read (*,"(a20)") name
	if (name(1:1) .ne. ' ') filen = name

!.... allocate memory for the physical grid

	open (unit=10, file=filen, form='unformatted', status='unknown')
	read(10) nx, ny, nz
	write(*,"('Nx = ',i4,', Ny = ',i4)") nx, ny
	allocate (xy(2,nx,ny), STAT=ier)
	if (ier .ne. 0) call error('pre$', &
	    'Insufficient Memory for phys. grid$')
	read(10) (((xy(1,i,j), i = 1, nx), j = 1, ny), k = 1, nz), &
	         (((xy(2,i,j), i = 1, nx), j = 1, ny), k = 1, nz), &
	         (((      tmp, i = 1, nx), j = 1, ny), k = 1, nz)
	close(10)

!.... allocate memory for the computational grid

	allocate (xi(nx), eta(ny), STAT=ier)
	if (ier .ne. 0) call error('pre$', &
	   'Insufficient Memory for comp. grid$')

!.... make the computational grid

	dxi  = one / float(nx-1)
	deta = one / float(ny-1)
	
	do i = 1, nx
	  xi(i) = (i-1) * dxi
	end do
	
	do j = 1, ny
	  eta(j) = (j-1) * deta
	end do
	
!.... allocate space for first derivatives

	allocate (g1xy(2,nx,ny), g2xy(2,nx,ny), STAT=ier)
	if (ier .ne. 0) call error('pre$','Insufficient Memory for gxy$')

!.... compute first derivative metrics using finite differences
	
	call grad(2, nx, ny, xy, g1xy, g2xy, dxi, deta, optx, opty, &
                  xper, yper, lsym, rsym, bsym, tsym, carp)

!.... form the determinate of the jacobian

	allocate (IdetJ(nx,ny), STAT=ier)
	if (ier .ne. 0) call error('pre$','Insufficient Memory for detJ$')

	IdetJ(:,:) = one / ( g1xy(1,:,:) * g2xy(2,:,:) - &
	                     g2xy(1,:,:) * g1xy(2,:,:) )
	
!.... allocate memory for the grid metrics

	allocate( m(2,nx,ny),  m11(nx,ny), m12(nx,ny), m22(nx,ny), &
	          n(2,nx,ny),  n11(nx,ny), n12(nx,ny), n22(nx,ny), &
		  STAT=ier)
	if (ier .ne. 0) call error('pre$','Insufficient Memory for metrics$')

!.... form the first derivative metrics
	
	m(1,:,:) =  g2xy(2,:,:) * IdetJ
	n(1,:,:) = -g1xy(2,:,:) * IdetJ
	m(2,:,:) = -g2xy(1,:,:) * IdetJ
	n(2,:,:) =  g1xy(1,:,:) * IdetJ

	allocate (g1m(2,nx,ny), g2m(2,nx,ny), &
	          g1n(2,nx,ny), g2n(2,nx,ny), STAT=ier)
	if (ier .ne. 0) call error('pre$','Insufficient Memory for gm and gn$')

!.... compute second derivatives of x using finite differences

	allocate (g11xy(2,nx,ny), g12xy(2,nx,ny), g22xy(2,nx,ny), STAT=ier)
	if (ier .ne. 0) call error('pre$','Insufficient Memory for gxy$')

	call grad(2, nx, ny, g1xy, g11xy, g12xy, dxi, deta, &
	          optx, opty, xper, yper, lsym, rsym, bsym, tsym, carp)
	call grad(2, nx, ny, g2xy, g12xy, g22xy, dxi, deta, &
	          optx, opty, xper, yper, lsym, rsym, bsym, tsym, carp)

!.... compute second derivative metrics using finite difference

	call grad(2, nx, ny, m, g1m, g2m, dxi, deta, &
	          optx, opty, xper, yper, lsym, rsym, bsym, tsym, carp)

	call grad(2, nx, ny, n, g1n, g2n, dxi, deta, &
	          optx, opty, xper, yper, lsym, rsym, bsym, tsym, carp)
	
 	m11 = m(1,:,:) * g1m(1,:,:) + n(1,:,:) * g2m(1,:,:)
 	n11 = m(1,:,:) * g1n(1,:,:) + n(1,:,:) * g2n(1,:,:)

 	m22 = m(2,:,:) * g1m(2,:,:) + n(2,:,:) * g2m(2,:,:)
 	n22 = m(2,:,:) * g1n(2,:,:) + n(2,:,:) * g2n(2,:,:)
	
	m12 = m(1,:,:) * g1m(2,:,:) + n(1,:,:) * g2m(2,:,:)	
	n12 = m(1,:,:) * g1n(2,:,:) + n(1,:,:) * g2n(2,:,:)
		
	j = 1
	do i = 1, nx
	  write(55,21) xy(1,i,j), m12(i,j), &
	              m(2,i,j)*g1m(1,i,j)+n(2,i,j)*g2m(1,i,j), &
		      m12(i,j) - (m(2,i,j)*g1m(1,i,j)+n(2,i,j)*g2m(1,i,j)), &
		      m(1,i,j)*g2m(1,i,j)+m(2,i,j)*g2m(2,i,j), &
		      n(1,i,j)*g1n(1,i,j)+n(2,i,j)*g1n(2,i,j)
	  write(56,21) xy(1,i,j), n(1,i,j), n(2,i,j), n11(i,j), n12(i,j), &
	               n22(i,j)
 21	  format(7(1pe14.7,1x))
	end do

	deallocate( g1m, g2m, g1n, g2n )

!.... write out the metric file

        filen = 'metric.new'
 	write (*,"('Enter grid metric file name [',a,']? ',$)") &
     		    filen(1:index(filen,' '))
	read (*,"(a20)") name
	if (name(1:1) .ne. ' ') filen = name
	open (unit=10, file=filen, form='unformatted', status='unknown')
	write(10) m(1,:,:), m(2,:,:), n(1,:,:), n(2,:,:), &
	          m11, m12, m22, n11, n12, n22
	close(10)

!.... write out plot3d files

	allocate (q(5,nx,ny), STAT=ier)
	if (ier .ne. 0) call error('pre$','Insufficient Memory for q$')

	q(1,:,:) = m(1,:,:)
	q(2,:,:) = n(1,:,:)
	q(3,:,:) = m(2,:,:)
	q(4,:,:) = n(2,:,:)
	q(5,:,:) = IdetJ
	
!.... diagnostic

	j = 1
	do i = 1, nx
	  write(54,21) xy(1,i,j), &
	             (-m(1,i,j)*m(2,i,j)*m11(i,j) + m(1,i,j)**2*m12(i,j) - &
	             m(2,i,j)**2*m12(i,j) + m(1,i,j)*m(2,i,j)*m22(i,j) ) * &
		     sqrt(n(1,i,j)**2 + n(2,i,j)**2) / 			&
		     sqrt(m(1,i,j)**2 + m(2,i,j)**2) /			&
		     (m(1,i,j)*n(2,i,j)-m(2,i,j)*n(1,i,j)), 		&
		     (n(1,i,j)*n(2,i,j)*n11(i,j) - n(1,i,j)**2*n12(i,j) + &
	             n(2,i,j)**2*n12(i,j) - n(1,i,j)*n(2,i,j)*n22(i,j) ) * &
		     sqrt(m(1,i,j)**2 + m(2,i,j)**2) / 			&
		     sqrt(n(1,i,j)**2 + n(2,i,j)**2) / 			&
		     (m(1,i,j)*n(2,i,j)-m(2,i,j)*n(1,i,j))
	end do

	name = 'm1.dat'
	call wdata(name, q, nx, ny, 1, 5, 0.0, 0.0, 0.0, 0.0)
		
	q(1,:,:) = m(1,:,:) - 2.0 * xy(2,:,:) / ( acos(-1.0) * &
	                      ( (xy(1,:,:)-0.1)**2 + xy(2,:,:)**2 ) )
	q(2,:,:) = n(1,:,:) - (xy(1,:,:) - 0.1) / sqrt( (xy(1,:,:)-0.1)**2 + &
	                       xy(2,:,:)**2 )
	q(3,:,:) = m(2,:,:) - 2.0 * (0.1 - xy(1,:,:)) / ( acos(-1.0) * &
	                      ( (xy(1,:,:)-0.1)**2 + xy(2,:,:)**2 ) )
	q(4,:,:) = n(2,:,:) - xy(2,:,:) / sqrt( (xy(1,:,:)-0.1)**2 + &
	                      xy(2,:,:)**2 )
	q(5,:,:) = IdetJ
	
	name = 'error.dat'
	call wdata(name, q, nx, ny, 1, 5, 0.0, 0.0, 0.0, 0.0)

	q(1,:,:) = m11
	q(2,:,:) = m12
	q(3,:,:) = m12
	q(4,:,:) = m22
	q(5,:,:) = IdetJ
	
	name = 'm11.dat'
	call wdata(name, q, nx, ny, 1, 5, 0.0, 0.0, 0.0, 0.0)

	q(1,:,:) = n11
	q(2,:,:) = n12
	q(3,:,:) = n12
	q(4,:,:) = n22
	q(5,:,:) = IdetJ
	
	name = 'n11.dat'
	call wdata(name, q, nx, ny, 1, 5, 0.0, 0.0, 0.0, 0.0)

	q(1,:,:) = g11xy(1,:,:)
	q(2,:,:) = g11xy(2,:,:)
	q(3,:,:) = g12xy(1,:,:)
	q(4,:,:) = g12xy(2,:,:)
	q(5,:,:) = g22xy(1,:,:)
	
	name = 'grad.dat'
	call wdata(name, q, nx, ny, 1, 5, 0.0, 0.0, 0.0, 0.0)

	deallocate( g1xy, g2xy)

!.... write out the derivatives on the body

	write(*,"('Enter j ==> ',$)") 
	read(*,*) j
	do i = 1, nx
	  write(50,20) xi(i), m(1,i,j), m(2,i,j)
	  write(51,20) xi(i), n(1,i,j), n(2,i,j)
	  write(52,20) xi(i), m11(i,j), n11(i,j), m22(i,j), m12(i,j), n12(i,j)
	  write(53,20) xi(i), IdetJ(i,j), one / Idetj(i,j)
	end do

!==============================================================================
!
! Test the derivative operators
!
!==============================================================================

!.... put a sine wave on the mesh

	allocate (f(nx,ny), g1l(nx,ny), g2l(nx,ny), g1(nx,ny), g2(nx,ny), &
             STAT=ier)
	if (ier .ne. 0) call error('pre$','Insufficient Memory for function$')

	kappa = two * pi / 2.5
	f(:,:) = sin( kappa * xy(1,:,:))

!.... compute the gradient with finite differences

	call grad( 1, nx, ny, f, g1l, g2l, dxi, deta, &
	           optx, opty, xper, yper, lsym, rsym, bsym, tsym, carp)

!.... transform to physical coordinates

	g1 = g1l(:,:) * m(1,:,:) + g2l(:,:) * n(1,:,:)
	g2 = g1l(:,:) * m(2,:,:) + g2l(:,:) * n(2,:,:)

	q(1,:,:) = f
	q(2,:,:) = g1 / kappa
	q(3,:,:) = g2
	q(4,:,:) = abs(g1 - kappa * cos( kappa * xy(1,:,:) ) )
	q(5,:,:) = abs(g2 - zero)
	
	name = 'f.dat'
	call wdata(name, q, nx, ny, 1, 5, 0.0, 0.0, 0.0, 0.0)

	deallocate( g1, g2)
	
!... allocate memory for second derivatives

	allocate (g11(nx,ny),  g12(nx,ny),  g22(nx,ny), &
	          g11l(nx,ny), g12l(nx,ny), g22l(nx,ny), STAT=ier)
	if (ier .ne. 0) call error('pre$','Insufficient Memory for g2$')

!.... compute second derivatives using finite differencs

	call grad2( 1, nx, ny, f, g1l, g11l, g12l, g22l, dxi, deta, &
	            optx, opty, xper, yper, lsym, rsym, bsym, tsym, carp)

	g11  =  g11l(:,:)        * m(1,:,:) * m(1,:,:)	+ &
		two * g12l(:,:)  * m(1,:,:) * n(1,:,:) 	+ &
		g22l(:,:)        * n(1,:,:) * n(1,:,:)	+ &
		g1l(:,:)         * m11(:,:)    		+ &
		g2l(:,:)         * n11(:,:)

	g12  =  g11l(:,:)        * m(1,:,:) * m(2,:,:)	+ &
		g12l(:,:)        * m(1,:,:) * n(2,:,:)	+ &
		g12l(:,:)        * m(2,:,:) * n(1,:,:) 	+ &
		g22l(:,:)        * n(1,:,:) * n(2,:,:)	+ &
		g1l(:,:)         * m12(:,:)		+ &
		g2l(:,:)         * n12(:,:)

	g22  =  g11l(:,:)        * m(2,:,:) * m(2,:,:)	+ &
		two * g12l(:,:)  * m(2,:,:) * n(2,:,:)	+ &
		g22l(:,:)        * n(2,:,:) * n(2,:,:)	+ &
		g1l(:,:)         * m22(:,:)		+ &
		g2l(:,:)         * n22(:,:)

	q(1,:,:) = g11 / kappa**2
	q(2,:,:) = g12
	q(3,:,:) = g22
	q(4,:,:) = abs(g11/kappa**2 + sin( kappa * xy(1,:,:) ) )
	q(5,:,:) = zero
	
	name = 'g2.dat'
	call wdata(name, q, nx, ny, 1, 5, 0.0, 0.0, 0.0, 0.0)

	stop
	
 20	format(8(1pe13.6,1x))

	end
