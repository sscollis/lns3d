!=============================================================================!
        module local3d
!
!  Local variables for three-dimensional code
!
!=============================================================================!
        use global

        complex, allocatable :: c1v(:,:,:), c2v(:,:,:)
        complex, allocatable :: c11v(:,:,:), c12v(:,:,:), c22v(:,:,:)
        !$sgi distribute c1v(*,*,block), c2v(*,*,block)
        !$sgi distribute c11v(*,*,block), c12v(*,*,block), c22v(*,*,block)
        
        real, allocatable :: ABhi(:,:,:), Dhi(:,:,:,:)

        contains

!=============================================================================!
        subroutine mlocal3d
!
!  allocate space for 3d local variables
!
!=============================================================================!
        implicit none
        
        integer :: i, j, idof, ier
        character*80 :: code='mLocal3D$'
!=============================================================================!
        allocate( c1v(ndof,nx,ny),  c2v(ndof,nx,ny),            &
                  c11v(ndof,nx,ny), c12v(ndof,nx,ny),           &
                  c22v(ndof,nx,ny), STAT=ier )
        if (ier .ne. 0) call error(code,'Insufficient Memory for gradients$')
          
        !$doacross local(i,idof)
        !$omp parallel do private(i,idof)
        do j = 1, ny
          do i = 1, nx
            do idof = 1, ndof
              c1v(idof,i,j) = zero
              c2v(idof,i,j) = zero
              c11v(idof,i,j) = zero
              c12v(idof,i,j) = zero
              c22v(idof,i,j) = zero
            end do
          end do
        end do

!.... 3-d matrices

          allocate( Dhi(ndof,ndof,nx,ny), STAT=ier )
          if (ier .ne. 0) call error(code,'Insufficient Memory$')
          allocate( ABhi(6,nx,ny), STAT=ier )
          if (ier .ne. 0) call error(code,'Insufficient Memory$')

        end subroutine mlocal3d

        end module local3d
!=============================================================================!
