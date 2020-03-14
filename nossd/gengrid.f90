!-----------------------------------------------------------------------------
	subroutine gengrid 
!  
!  Generate the grid 
!  
!-----------------------------------------------------------------------------
	use global
	implicit none

	integer :: i, j, k, ier, mem = 0, ij
	real    :: tmp
	character*80 name
!-----------------------------------------------------------------------------

!.... read in the grid

	name = 'grid.dat'
	open (unit=10, file=name, form='unformatted', &
	      status='old', err=1000)
	read(10,err=1000,end=1000) nx, ny, nz
	write(*,"(' Reading grid for (',i4,',',i4,',',i4,')')") nx, ny, nz
	
!.... check to make sure that Ny is not a power of two for implicit

	mem = mem + ny*nx*2 + nz + nx + ny
	allocate (x(ny,nx), y(ny,nx), z(nz), xi(nx), eta(ny), STAT=ier)
	if (ier .ne. 0) call error('gengrid$','Insufficient Memory for grid$')
	read(10,err=1000,end=1000) &
	         (((x(j,i), i = 1, nx), j = 1, ny), k = 1, nz), &
	         (((y(j,i), i = 1, nx), j = 1, ny), k = 1, nz), &
	         (((   tmp, i = 1, nx), j = 1, ny), k = 1, nz)
	close(10)

!.... make the xi grid
	
	dxi = one / float(nx-1)
	
	do i = 1, nx
	  xi(i) = zero + (i-1) * dxi
	end do

!.... make the eta grid

	deta = one / float(ny-1)
	
	do j = 1, ny
	  eta(j) = float(j-1) * deta
	end do
	
!.... make z grid

	z = zero

!.... allocate storage for metrics

	mem = mem + ny*nx*10
	allocate (m1(ny*nx),  m2(ny*nx),  n1(ny*nx),  n2(ny*nx), &
	          m11(ny*nx), m12(ny*nx), m22(ny*nx),            &
	          n11(ny*nx), n12(ny*nx), n22(ny*nx), STAT=ier)
	if (ier .ne. 0) call error('gengrid$', &
	                'Insufficient Memory for mesh metrics$')

!.... read in the mesh metrics

	open (unit=10,file='metric.dat',form='unformatted', &
	      status='old', err=1100)
	read(10,err=1100,end=1000) &
	      m1, m2, n1, n2, m11, m12, m22, n11, n12, n22
	close(10)

!.... form the metric products

	mem = mem + ny*nx*3
	allocate (m1m1(ny*nx), m1m2(ny*nx), m2m2(ny*nx), STAT=ier)
	if (ier .ne. 0) call error('gengrid$','Insufficient Memory$')

	m1m1 = m1 * m1
	m1m2 = m1 * m2
	m2m2 = m2 * m2
	
	mem = mem + ny*nx*3
	allocate (n1n1(ny*nx), n1n2(ny*nx), n2n2(ny*nx), STAT=ier)
	if (ier .ne. 0) call error('gengrid$','Insufficient Memory$')

	n1n1 = n1 * n1
	n1n2 = n1 * n2
	n2n2 = n2 * n2
	
	mem = mem + ny*nx*3
	allocate (m1n1(ny*nx), m1n2(ny*nx), m2n1(ny*nx), m2n2(ny*nx), STAT=ier)
	if (ier .ne. 0) call error('gengrid$','Insufficient Memory$')

	m1n1 = m1 * n1
	m1n2 = m1 * n2
	m2n1 = m2 * n1
	m2n2 = m2 * n2

!.... compute the boundary normal vectors.  These are defined such that
!.... they allways point outside of the domain.

	allocate( bnl(ny,2), bnr(ny,2), bnb(nx,2), bnt(nx,2), STAT=ier )
	if (ier .ne. 0) call error('gengrid$','Insufficient Memory$')

	j = 1
	do i = 1, nx
	  ij = j + (i-1) * ny
	  bnb(i,1) = -n1(ij) / sqrt( n1(ij)**2 + n2(ij)**2 )
	  bnb(i,2) = -n2(ij) / sqrt( n1(ij)**2 + n2(ij)**2 )
	end do

	j = ny
	do i = 1, nx
	  ij = j + (i-1) * ny
	  bnt(i,1) = n1(ij) / sqrt( n1(ij)**2 + n2(ij)**2 )
	  bnt(i,2) = n2(ij) / sqrt( n1(ij)**2 + n2(ij)**2 )
	end do

	i = 1
	do j = 1, ny
	  ij = j + (i-1) * ny
	  bnl(j,1) = -m1(ij) / sqrt( m1(ij)**2 + m2(ij)**2 )
	  bnl(j,2) = -m2(ij) / sqrt( m1(ij)**2 + m2(ij)**2 )
	end do

	i = nx
	do j = 1, ny
	  ij = j + (i-1) * ny
	  bnr(j,1) = m1(ij) / sqrt( m1(ij)**2 + m2(ij)**2 )
	  bnr(j,2) = m2(ij) / sqrt( m1(ij)**2 + m2(ij)**2 )
	end do	  

	write(*,"(' GenGrid allocated ===> ',1pe13.6,' words')") float(mem)

	return
1000	call error('gengrid$','Error reading grid$')
1100	call error('gengrid$','Error reading metrics$')
	end
