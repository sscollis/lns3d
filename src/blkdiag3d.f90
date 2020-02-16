!=============================================================================!
        subroutine blkdiag3D( mat, Dh, Dhi, Vh11, Vh22 )
!  
!  For the block diagonal matrix
!
!=============================================================================!
        use global
        use stencil
        implicit none

        complex :: mat(ndof,ndof,nx,ny)
        real    :: Dh(ndof,ndof,nx,ny), Dhi(ndof,ndof,nx,ny), 
        real    :: Vh11(6,nx,ny), Vh22(6,nx,ny)

        real    :: a3, b3

        integer :: i, j, idof, jdof
!=============================================================================!

!.... fourth-order stencil

        a3 = da3 / dxi**2
        b3 = da3 / deta**2

        !$omp parallel do private (i,idof,jdof)
        loop_j: do j = 1, ny
        loop_i: do i = 1, nx

!.... \hat{D} term

        do idof = 1, ndof
          do jdof = 1, ndof
            mat(idof,jdof,i,j) = Dh(idof,jdof,i,j) + im * Dhi(idof,jdof,i,j)
          end do
        end do

!.... \hat{V}_{\xi\xi} term

        mat(2,2,i,j) = mat(2,2,i,j) - a3 * Vh11(1,i,j)
        mat(2,3,i,j) = mat(2,3,i,j) - a3 * Vh11(5,i,j)
        mat(3,2,i,j) = mat(3,2,i,j) - a3 * Vh11(6,i,j)
        mat(3,3,i,j) = mat(3,3,i,j) - a3 * Vh11(2,i,j)
        mat(4,4,i,j) = mat(4,4,i,j) - a3 * Vh11(3,i,j)
        mat(5,5,i,j) = mat(5,5,i,j) - a3 * Vh11(4,i,j)

!.... \hat{v}_{\xi\eta} term

!.... \hat{V}_{\eta\eta} term

        mat(2,2,i,j) = mat(2,2,i,j) - b3 * Vh22(1,i,j)
        mat(2,3,i,j) = mat(2,3,i,j) - b3 * Vh22(5,i,j)
        mat(3,2,i,j) = mat(3,2,i,j) - b3 * Vh22(6,i,j)
        mat(3,3,i,j) = mat(3,3,i,j) - b3 * Vh22(2,i,j)
        mat(4,4,i,j) = mat(4,4,i,j) - b3 * Vh22(3,i,j)
        mat(5,5,i,j) = mat(5,5,i,j) - b3 * Vh22(4,i,j)

        end do loop_i
        end do loop_j

        end subroutine blkdiag3D
