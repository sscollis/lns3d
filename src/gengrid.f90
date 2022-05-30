!-----------------------------------------------------------------------------
        subroutine gengrid 
!  
!  Generate the grid 
!  
!-----------------------------------------------------------------------------
        use global
        implicit none

        integer :: i, j, k, ier, mem = 0
        real    :: tmp
        character(80) name
!-----------------------------------------------------------------------------

!.... read in the grid

        name = 'grid.dat'
        open (unit=10, file=name, form='unformatted', &
              status='old', err=1000)
        read(10,err=1000,end=1000) nx, ny, nz
        write(*,"(' Reading grid for (',i4,',',i4,',',i4,')')") nx, ny, nz
        
!.... check to make sure that Ny is not a power of two for implicit

        mem = mem + ny*nx*2 + nz + nx + ny
        allocate (x(nx,ny), y(nx,ny), z(nz), xi(nx), eta(ny), STAT=ier)
        if (ier .ne. 0) call error('gengrid$','Insufficient Memory for grid$')
        read(10,err=1000,end=1000) &
                 (((x(i,j), i = 1, nx), j = 1, ny), k = 1, nz), &
                 (((y(i,j), i = 1, nx), j = 1, ny), k = 1, nz), &
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

        mem = mem + nx*ny*10
        allocate (m1(nx,ny),  m2(nx,ny),  n1(nx,ny),  n2(nx,ny), &
                  m11(nx,ny), m12(nx,ny), m22(nx,ny),            &
                  n11(nx,ny), n12(nx,ny), n22(nx,ny), STAT=ier)
        if (ier .ne. 0) call error('gengrid$', &
                        'Insufficient Memory for mesh metrics$')

        !$doacross local(i)
        !$omp parallel do private(i)
        do j = 1, ny
          do i = 1, nx
            m1(i,j) = zero; m2(i,j) = zero; 
            n1(i,j) = zero; n2(i,j) = zero; 
            m11(i,j) = zero; m12(i,j) = zero; m22(i,j) = zero;
            n11(i,j) = zero; n12(i,j) = zero; n22(i,j) = zero;
          end do
        end do

!.... read in the mesh metrics

        open (unit=10,file='metric.dat',form='unformatted', &
              status='old', err=1100)
        read(10,err=1100,end=1000) &
              m1, m2, n1, n2, m11, m12, m22, n11, n12, n22
        close(10)

!.... form the metric products

        mem = mem + nx*ny*3
        allocate (m1m1(nx,ny), m1m2(nx,ny), m2m2(nx,ny), STAT=ier)
        if (ier .ne. 0) call error('gengrid$','Insufficient Memory$')

        !$doacross local(i)
        !$omp parallel do private(i)
        do j = 1, ny
          do i = 1, nx
            m1m1(i,j) = m1(i,j) * m1(i,j)
            m1m2(i,j) = m1(i,j) * m2(i,j)
            m2m2(i,j) = m2(i,j) * m2(i,j)
          end do
        end do
        
        mem = mem + nx*ny*3
        allocate (n1n1(nx,ny), n1n2(nx,ny), n2n2(nx,ny), STAT=ier)
        if (ier .ne. 0) call error('gengrid$','Insufficient Memory$')

        !$doacross local(i)
        !$omp parallel do private(i)
        do j = 1, ny
          do i = 1, nx
            n1n1(i,j) = n1(i,j) * n1(i,j)
            n1n2(i,j) = n1(i,j) * n2(i,j)
            n2n2(i,j) = n2(i,j) * n2(i,j)
          end do
        end do

        mem = mem + nx*ny*3
        allocate (m1n1(nx,ny), m1n2(nx,ny), m2n1(nx,ny), m2n2(nx,ny), STAT=ier)
        if (ier .ne. 0) call error('gengrid$','Insufficient Memory$')

        !$doacross
        !$omp parallel do
        do j = 1, ny
          do i = 1, nx
            m1n1(i,j) = m1(i,j) * n1(i,j)
            m1n2(i,j) = m1(i,j) * n2(i,j)
            m2n1(i,j) = m2(i,j) * n1(i,j)
            m2n2(i,j) = m2(i,j) * n2(i,j)
          end do
        end do

!.... compute the boundary normal vectors.  These are defined such that
!.... they allways point outside of the domain.

        allocate( bnl(ny,2), bnr(ny,2), bnb(nx,2), bnt(nx,2), STAT=ier )
        if (ier .ne. 0) call error('gengrid$','Insufficient Memory$')

        j = 1
        do i = 1, nx
          bnb(i,1) = -n1(i,j) / sqrt( n1(i,j)**2 + n2(i,j)**2 )
          bnb(i,2) = -n2(i,j) / sqrt( n1(i,j)**2 + n2(i,j)**2 )
        end do

        j = ny
        do i = 1, nx
          bnt(i,1) = n1(i,j) / sqrt( n1(i,j)**2 + n2(i,j)**2 )
          bnt(i,2) = n2(i,j) / sqrt( n1(i,j)**2 + n2(i,j)**2 )
        end do

        i = 1
        do j = 1, ny
          bnl(j,1) = -m1(i,j) / sqrt( m1(i,j)**2 + m2(i,j)**2 )
          bnl(j,2) = -m2(i,j) / sqrt( m1(i,j)**2 + m2(i,j)**2 )
        end do

        i = nx
        do j = 1, ny
          bnr(j,1) = m1(i,j) / sqrt( m1(i,j)**2 + m2(i,j)**2 )
          bnr(j,2) = m2(i,j) / sqrt( m1(i,j)**2 + m2(i,j)**2 )
        end do    

!       write(*,"(' GenGrid allocated ===> ',1pe13.6,' words')") float(mem)

        return
1000    call error('gengrid$','Error reading grid$')
1100    call error('gengrid$','Error reading metrics$')
        end
