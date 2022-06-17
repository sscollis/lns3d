!=============================================================================!
!
!  These routines generate the factored LHS for the Navier-Stokes Equations
!  using 2nd order differencing
!
!=============================================================================!
        subroutine lhs1( mat, Ah, Dh, Vh11, bc, dtl, vl)
!  
!  Form the LHS for the \xi direction
!
!=============================================================================!
        use global
        use buffer
        implicit none
        
        real :: mat(3,ndof,ndof,nx,ny), vl(ndof,nx,ny)
        real :: Ah(ndof,ndof,nx,ny), Dh(ndof,ndof,nx,ny), Vh11(6,nx,ny)
        real :: bc(14,ndof,ndof,ny)
        real :: dtl(nx,ny)

!.... first derivative operator

        real :: a1, a2, a3
        real :: b1, b2, b3
        real :: c1, c2, c3

        real :: eps_i

        integer :: i, j, idof, jdof
!=============================================================================!

!.... second-order stencil

        a1 = -alfa * pt5 / dxi
        a2 =  alfa * zero
        a3 =  alfa * pt5 / dxi

        b1 =  alfa / dxi**2
        b2 = -alfa * two / dxi**2
        b3 =  alfa / dxi**2
        
        eps_i = four * eps_e
        c1 = -alfa * eps_i
        c2 =  alfa * two * eps_i
        c3 = -alfa * eps_i

        !$omp parallel do private (i,idof,jdof)
        loop_j: do j = 1, ny
        loop_i: do i = 1, nx

!.... \hat{A} term

        do idof = 1, ndof
          do jdof = 1, ndof
            mat(1,idof,jdof,i,j) = a1 * dtl(i,j) * Ah(idof,jdof,i,j)    
            mat(2,idof,jdof,i,j) = a2 * dtl(i,j) * Ah(idof,jdof,i,j)    
            mat(3,idof,jdof,i,j) = a3 * dtl(i,j) * Ah(idof,jdof,i,j)
          end do
        end do

!.... \hat{D} term

        do idof = 1, ndof
          do jdof = 1, ndof
            mat(2,idof,jdof,i,j) = mat(2,idof,jdof,i,j) + &
                                   alfa * dtl(i,j) * Dh(idof,jdof,i,j)
          end do
        end do 

!.... \hat{V}_{\xi\xi} term

        mat(1,2,2,i,j) = mat(1,2,2,i,j) - b1 * dtl(i,j) * Vh11(1,i,j)
        mat(1,2,3,i,j) = mat(1,2,3,i,j) - b1 * dtl(i,j) * Vh11(5,i,j)
        mat(1,3,2,i,j) = mat(1,3,2,i,j) - b1 * dtl(i,j) * Vh11(6,i,j)
        mat(1,3,3,i,j) = mat(1,3,3,i,j) - b1 * dtl(i,j) * Vh11(2,i,j)
        mat(1,4,4,i,j) = mat(1,4,4,i,j) - b1 * dtl(i,j) * Vh11(3,i,j)
        mat(1,5,5,i,j) = mat(1,5,5,i,j) - b1 * dtl(i,j) * Vh11(4,i,j)

        mat(2,2,2,i,j) = mat(2,2,2,i,j) - b2 * dtl(i,j) * Vh11(1,i,j)
        mat(2,2,3,i,j) = mat(2,2,3,i,j) - b2 * dtl(i,j) * Vh11(5,i,j)
        mat(2,3,2,i,j) = mat(2,3,2,i,j) - b2 * dtl(i,j) * Vh11(6,i,j)
        mat(2,3,3,i,j) = mat(2,3,3,i,j) - b2 * dtl(i,j) * Vh11(2,i,j)
        mat(2,4,4,i,j) = mat(2,4,4,i,j) - b2 * dtl(i,j) * Vh11(3,i,j)
        mat(2,5,5,i,j) = mat(2,5,5,i,j) - b2 * dtl(i,j) * Vh11(4,i,j)

        mat(3,2,2,i,j) = mat(3,2,2,i,j) - b3 * dtl(i,j) * Vh11(1,i,j)
        mat(3,2,3,i,j) = mat(3,2,3,i,j) - b3 * dtl(i,j) * Vh11(5,i,j)
        mat(3,3,2,i,j) = mat(3,3,2,i,j) - b3 * dtl(i,j) * Vh11(6,i,j)
        mat(3,3,3,i,j) = mat(3,3,3,i,j) - b3 * dtl(i,j) * Vh11(2,i,j)
        mat(3,4,4,i,j) = mat(3,4,4,i,j) - b3 * dtl(i,j) * Vh11(3,i,j)
        mat(3,5,5,i,j) = mat(3,5,5,i,j) - b3 * dtl(i,j) * Vh11(4,i,j)

!.... I term

        mat(2,1,1,i,j) = mat(2,1,1,i,j) + one
        mat(2,2,2,i,j) = mat(2,2,2,i,j) + one
        mat(2,3,3,i,j) = mat(2,3,3,i,j) + one
        mat(2,4,4,i,j) = mat(2,4,4,i,j) + one
        mat(2,5,5,i,j) = mat(2,5,5,i,j) + one
                        
!.... implicit damping term

        if (eps_e .ne. zero) then

        mat(1,1,1,i,j) = mat(1,1,1,i,j) + c1 * dtl(i,j) * buff(i,j)
        mat(1,2,2,i,j) = mat(1,2,2,i,j) + c1 * dtl(i,j) * buff(i,j)
        mat(1,3,3,i,j) = mat(1,3,3,i,j) + c1 * dtl(i,j) * buff(i,j)
        mat(1,4,4,i,j) = mat(1,4,4,i,j) + c1 * dtl(i,j) * buff(i,j)
        mat(1,5,5,i,j) = mat(1,5,5,i,j) + c1 * dtl(i,j) * buff(i,j)

        mat(2,1,1,i,j) = mat(2,1,1,i,j) + c2 * dtl(i,j) * buff(i,j)
        mat(2,2,2,i,j) = mat(2,2,2,i,j) + c2 * dtl(i,j) * buff(i,j)
        mat(2,3,3,i,j) = mat(2,3,3,i,j) + c2 * dtl(i,j) * buff(i,j)
        mat(2,4,4,i,j) = mat(2,4,4,i,j) + c2 * dtl(i,j) * buff(i,j)
        mat(2,5,5,i,j) = mat(2,5,5,i,j) + c2 * dtl(i,j) * buff(i,j)

        mat(3,1,1,i,j) = mat(3,1,1,i,j) + c3 * dtl(i,j) * buff(i,j)
        mat(3,2,2,i,j) = mat(3,2,2,i,j) + c3 * dtl(i,j) * buff(i,j)
        mat(3,3,3,i,j) = mat(3,3,3,i,j) + c3 * dtl(i,j) * buff(i,j)
        mat(3,4,4,i,j) = mat(3,4,4,i,j) + c3 * dtl(i,j) * buff(i,j)
        mat(3,5,5,i,j) = mat(3,5,5,i,j) + c3 * dtl(i,j) * buff(i,j)

        end if

        end do loop_i
        end do loop_j

!.... apply boundary treatment to the LHS

        call lhsbt1( mat, Ah, Dh, Vh11, bc, dtl )

!.... apply boundary conditions to the LHS

        call lhsbc1( mat, bc, vl, vm )
        
        return
        end

!=============================================================================!
        subroutine lhs2( mat, Bh, Dh, Vh22, bc, dtl, vl)
!  
!  Form the LHS for the \eta direction
!
!=============================================================================!
        use global
        use buffer
        implicit none
        
        real :: mat(3,ndof,ndof,nx,ny), vl(ndof,nx,ny)
        real :: Bh(ndof,ndof,nx,ny), Dh(ndof,ndof,nx,ny), Vh22(6,nx,ny)
        real :: bc(14,ndof,ndof,nx)
        real :: dtl(nx,ny)

!.... first derivative operator

        real :: a1, a2, a3
        real :: b1, b2, b3
        real :: c1, c2, c3
        
        real :: eps_i
        
        integer :: i, j, idof, jdof
!=============================================================================!

!.... second-order stencil

        a1 = -alfa * pt5 / deta
        a2 =  alfa * zero
        a3 =  alfa * pt5 / deta

        b1 =  alfa / deta**2
        b2 = -alfa * two / deta**2
        b3 =  alfa / deta**2

        eps_i = four * eps_e
        c1 = -alfa * eps_i
        c2 =  alfa * two * eps_i
        c3 = -alfa * eps_i

        !$omp parallel do private (i,idof,jdof)
        loop_j: do j = 1, ny
        loop_i: do i = 1, nx

!.... \hat{B} term

        do idof = 1, ndof
          do jdof = 1, ndof
            mat(1,idof,jdof,i,j) = a1 * dtl(i,j) * Bh(idof,jdof,i,j)
            mat(2,idof,jdof,i,j) = a2 * dtl(i,j) * Bh(idof,jdof,i,j)
            mat(3,idof,jdof,i,j) = a3 * dtl(i,j) * Bh(idof,jdof,i,j)
          end do
        end do

!.... \hat{V}_{\eta\eta} term

        mat(1,2,2,i,j) = mat(1,2,2,i,j) - b1 * dtl(i,j) * Vh22(1,i,j)
        mat(1,2,3,i,j) = mat(1,2,3,i,j) - b1 * dtl(i,j) * Vh22(5,i,j)
        mat(1,3,2,i,j) = mat(1,3,2,i,j) - b1 * dtl(i,j) * Vh22(6,i,j)
        mat(1,3,3,i,j) = mat(1,3,3,i,j) - b1 * dtl(i,j) * Vh22(2,i,j)
        mat(1,4,4,i,j) = mat(1,4,4,i,j) - b1 * dtl(i,j) * Vh22(3,i,j)
        mat(1,5,5,i,j) = mat(1,5,5,i,j) - b1 * dtl(i,j) * Vh22(4,i,j)

        mat(2,2,2,i,j) = mat(2,2,2,i,j) - b2 * dtl(i,j) * Vh22(1,i,j)
        mat(2,2,3,i,j) = mat(2,2,3,i,j) - b2 * dtl(i,j) * Vh22(5,i,j)
        mat(2,3,2,i,j) = mat(2,3,2,i,j) - b2 * dtl(i,j) * Vh22(6,i,j)
        mat(2,3,3,i,j) = mat(2,3,3,i,j) - b2 * dtl(i,j) * Vh22(2,i,j)
        mat(2,4,4,i,j) = mat(2,4,4,i,j) - b2 * dtl(i,j) * Vh22(3,i,j)
        mat(2,5,5,i,j) = mat(2,5,5,i,j) - b2 * dtl(i,j) * Vh22(4,i,j)

        mat(3,2,2,i,j) = mat(3,2,2,i,j) - b3 * dtl(i,j) * Vh22(1,i,j)
        mat(3,2,3,i,j) = mat(3,2,3,i,j) - b3 * dtl(i,j) * Vh22(5,i,j)
        mat(3,3,2,i,j) = mat(3,3,2,i,j) - b3 * dtl(i,j) * Vh22(6,i,j)
        mat(3,3,3,i,j) = mat(3,3,3,i,j) - b3 * dtl(i,j) * Vh22(2,i,j)
        mat(3,4,4,i,j) = mat(3,4,4,i,j) - b3 * dtl(i,j) * Vh22(3,i,j)
        mat(3,5,5,i,j) = mat(3,5,5,i,j) - b3 * dtl(i,j) * Vh22(4,i,j)

!.... I term

        mat(2,1,1,i,j) = mat(2,1,1,i,j) + one
        mat(2,2,2,i,j) = mat(2,2,2,i,j) + one
        mat(2,3,3,i,j) = mat(2,3,3,i,j) + one
        mat(2,4,4,i,j) = mat(2,4,4,i,j) + one
        mat(2,5,5,i,j) = mat(2,5,5,i,j) + one

!.... implicit damping term

        if (eps_e .ne. zero) then

        mat(1,1,1,i,j) = mat(1,1,1,i,j) + c1 * dtl(i,j) * buff(i,j)
        mat(1,2,2,i,j) = mat(1,2,2,i,j) + c1 * dtl(i,j) * buff(i,j)
        mat(1,3,3,i,j) = mat(1,3,3,i,j) + c1 * dtl(i,j) * buff(i,j)
        mat(1,4,4,i,j) = mat(1,4,4,i,j) + c1 * dtl(i,j) * buff(i,j)
        mat(1,5,5,i,j) = mat(1,5,5,i,j) + c1 * dtl(i,j) * buff(i,j)

        mat(2,1,1,i,j) = mat(2,1,1,i,j) + c2 * dtl(i,j) * buff(i,j)
        mat(2,2,2,i,j) = mat(2,2,2,i,j) + c2 * dtl(i,j) * buff(i,j)
        mat(2,3,3,i,j) = mat(2,3,3,i,j) + c2 * dtl(i,j) * buff(i,j)
        mat(2,4,4,i,j) = mat(2,4,4,i,j) + c2 * dtl(i,j) * buff(i,j)
        mat(2,5,5,i,j) = mat(2,5,5,i,j) + c2 * dtl(i,j) * buff(i,j)

        mat(3,1,1,i,j) = mat(3,1,1,i,j) + c3 * dtl(i,j) * buff(i,j)
        mat(3,2,2,i,j) = mat(3,2,2,i,j) + c3 * dtl(i,j) * buff(i,j)
        mat(3,3,3,i,j) = mat(3,3,3,i,j) + c3 * dtl(i,j) * buff(i,j)
        mat(3,4,4,i,j) = mat(3,4,4,i,j) + c3 * dtl(i,j) * buff(i,j)
        mat(3,5,5,i,j) = mat(3,5,5,i,j) + c3 * dtl(i,j) * buff(i,j)

        end if
        
!.... sponge

        if (ispg.eq.1) then
          mat(2,1,1,i,j) = mat(2,1,1,i,j) + alfa * dtl(i,j) * spg(i,j)
          mat(2,2,2,i,j) = mat(2,2,2,i,j) + alfa * dtl(i,j) * spg(i,j)
          mat(2,3,3,i,j) = mat(2,3,3,i,j) + alfa * dtl(i,j) * spg(i,j)
          mat(2,4,4,i,j) = mat(2,4,4,i,j) + alfa * dtl(i,j) * spg(i,j)
          mat(2,5,5,i,j) = mat(2,5,5,i,j) + alfa * dtl(i,j) * spg(i,j)
        else if (ispg.ge.2) then
          mat(2,1,1,i,j) = mat(2,1,1,i,j) + alfa*dtl(i,j)*(spg(i,j)+spg2(i,j))
          mat(2,2,2,i,j) = mat(2,2,2,i,j) + alfa*dtl(i,j)*(spg(i,j)+spg2(i,j))
          mat(2,3,3,i,j) = mat(2,3,3,i,j) + alfa*dtl(i,j)*(spg(i,j)+spg2(i,j))
          mat(2,4,4,i,j) = mat(2,4,4,i,j) + alfa*dtl(i,j)*(spg(i,j)+spg2(i,j))
          mat(2,5,5,i,j) = mat(2,5,5,i,j) + alfa*dtl(i,j)*(spg(i,j)+spg2(i,j))
        end if

        end do loop_i
        end do loop_j

!.... apply boundary treatment to the LHS

        call lhsbt2( mat, Bh, Dh, Vh22, bc, spg, spg2, dtl )

!.... apply boundary conditions to the LHS

        call lhsbc2( mat, bc, vl )

        return
        end
