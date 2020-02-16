!=============================================================================!
        subroutine calcp( v, g2v, m1l, m2l, n1l, n2l, gpr )
!=============================================================================!
        use global
        use pot
        implicit none
        
        real :: v(ndof,nx,ny), g2v(ndof,nx,ny)
        real :: m1l(nx,ny), m2l(nx,ny), n1l(nx,ny), n2l(nx,ny)
        real :: gpr(2,nx,ny)

        !$sgi distribute v(*,*,block), g2v(*,*,block)
        !$sgi distribute m1l(*,block), m2l(*,block), n1l(*,block), n2l(*,block)
        !$sgi distribute gpr(*,*,block)

        integer :: i, j
!=============================================================================!
        !$doacross local(i)
        !$omp parallel do private(i)
        do j = 1, ny
          do i = 1, nx

            gpr(1,i,j) = m1l(i,j) * g1p(i,j) + n1l(i,j) / (gamma * Ma**2) * &
                       ( g2v(1,i,j) * v(5,i,j) + v(1,i,j) * g2v(5,1,j) )

            gpr(1,i,j) = m2l(i,j) * g1p(i,j) + n2l(i,j) / (gamma * Ma**2) * &
                       ( g2v(1,i,j) * v(5,i,j) + v(1,i,j) * g2v(5,i,j) )

          end do
        end do

        return
        end
