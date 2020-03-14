!=============================================================================!
!
!  These routines generate the factored LHS for the Navier-Stokes Equations
!  using 2nd order differencing
!
!=============================================================================!
        subroutine lhs1( mat, Ah, Dh, Vh11, bc, dtl, v)
!  
!  Form the LHS for the \xi direction
!
!=============================================================================!
        use global
        use buff_mod
        implicit none
        
        real :: mat(ny*nx,ndof,ndof,3), v(ny*nx,ndof)
        real :: Ah(ny*nx,ndof,ndof), Dh(ny*nx,ndof,ndof), Vh11(ny*nx,6)
        real :: bc(ny*ndof*ndof*14)
        real :: dtl(ny*nx)

!.... first derivative operator

        real :: a1, a2, a3

        integer :: lrec, ier, istat
        integer :: idof, jdof
!=============================================================================!

!.... second-order stencil

        a1 = -alfa * pt5 / dxi
        a2 =  alfa * zero
        a3 =  alfa * pt5 / dxi

!.... \hat{A} term

        do idof = 1, ndof
          do jdof = 1, ndof
            mat(:,idof,jdof,1) = a1 * dtl(:) * Ah(:,idof,jdof)  
            mat(:,idof,jdof,2) = a2 * dtl(:) * Ah(:,idof,jdof)  
            mat(:,idof,jdof,3) = a3 * dtl(:) * Ah(:,idof,jdof)
          end do
        end do

!.... \hat{D} term

        do idof = 1, ndof
          do jdof = 1, ndof
            mat(:,idof,jdof,2) = mat(:,idof,jdof,2) + &
                                 alfa * dtl(:) * Dh(:,idof,jdof)
          end do
        end do 

!.... \hat{V}_{\xi\xi} term and I term

        call lhs1l( mat, Vh11, dtl, buff )

!.... apply boundary treatment to the LHS

        call lhsbt1( mat, Ah, Dh, Vh11, bc, dtl )

!.... apply boundary conditions to the LHS

        call lhsbc1( mat, bc, v, vm )
        
        return
        end

!=============================================================================!
        subroutine lhs1l( mat, Vh11, dtl, buff )
!  
!  Form the LHS for the \xi direction
!
!=============================================================================!
        use global
        implicit none
        
        real :: mat(ny*nx,ndof,ndof,3), buff(ny*nx)
        real :: dtl(ny*nx)
        
        real :: Vh11(ny*nx,6)

!.... second derivative operator

        real :: b1, b2, b3

        real :: eps_i
!=============================================================================!

!.... second-order stencil

        b1 =  alfa / dxi**2
        b2 = -alfa * two / dxi**2
        b3 =  alfa / dxi**2
        
!.... \hat{V}_{\xi\xi} term

        mat(:,2,2,1) = mat(:,2,2,1) - b1 * dtl(:) * Vh11(:,1)
        mat(:,2,3,1) = mat(:,2,3,1) - b1 * dtl(:) * Vh11(:,5)
        mat(:,3,2,1) = mat(:,3,2,1) - b1 * dtl(:) * Vh11(:,6)
        mat(:,3,3,1) = mat(:,3,3,1) - b1 * dtl(:) * Vh11(:,2)
        mat(:,4,4,1) = mat(:,4,4,1) - b1 * dtl(:) * Vh11(:,3)
        mat(:,5,5,1) = mat(:,5,5,1) - b1 * dtl(:) * Vh11(:,4)

        mat(:,2,2,2) = mat(:,2,2,2) - b2 * dtl(:) * Vh11(:,1)
        mat(:,2,3,2) = mat(:,2,3,2) - b2 * dtl(:) * Vh11(:,5)
        mat(:,3,2,2) = mat(:,3,2,2) - b2 * dtl(:) * Vh11(:,6)
        mat(:,3,3,2) = mat(:,3,3,2) - b2 * dtl(:) * Vh11(:,2)
        mat(:,4,4,2) = mat(:,4,4,2) - b2 * dtl(:) * Vh11(:,3)
        mat(:,5,5,2) = mat(:,5,5,2) - b2 * dtl(:) * Vh11(:,4)

        mat(:,2,2,3) = mat(:,2,2,3) - b3 * dtl(:) * Vh11(:,1)
        mat(:,2,3,3) = mat(:,2,3,3) - b3 * dtl(:) * Vh11(:,5)
        mat(:,3,2,3) = mat(:,3,2,3) - b3 * dtl(:) * Vh11(:,6)
        mat(:,3,3,3) = mat(:,3,3,3) - b3 * dtl(:) * Vh11(:,2)
        mat(:,4,4,3) = mat(:,4,4,3) - b3 * dtl(:) * Vh11(:,3)
        mat(:,5,5,3) = mat(:,5,5,3) - b3 * dtl(:) * Vh11(:,4)

!.... I term

        mat(:,1,1,2) = mat(:,1,1,2) + one
        mat(:,2,2,2) = mat(:,2,2,2) + one
        mat(:,3,3,2) = mat(:,3,3,2) + one
        mat(:,4,4,2) = mat(:,4,4,2) + one
        mat(:,5,5,2) = mat(:,5,5,2) + one
                        
!.... implicit damping term

        if (eps_e .ne. zero) then

        eps_i = four * eps_e
        
        b1 = -alfa * eps_i
        b2 =  alfa * two * eps_i
        b3 = -alfa * eps_i

        mat(:,1,1,1) = mat(:,1,1,1) + b1 * dtl(:) * buff
        mat(:,2,2,1) = mat(:,2,2,1) + b1 * dtl(:) * buff
        mat(:,3,3,1) = mat(:,3,3,1) + b1 * dtl(:) * buff
        mat(:,4,4,1) = mat(:,4,4,1) + b1 * dtl(:) * buff
        mat(:,5,5,1) = mat(:,5,5,1) + b1 * dtl(:) * buff

        mat(:,1,1,2) = mat(:,1,1,2) + b2 * dtl(:) * buff
        mat(:,2,2,2) = mat(:,2,2,2) + b2 * dtl(:) * buff
        mat(:,3,3,2) = mat(:,3,3,2) + b2 * dtl(:) * buff
        mat(:,4,4,2) = mat(:,4,4,2) + b2 * dtl(:) * buff
        mat(:,5,5,2) = mat(:,5,5,2) + b2 * dtl(:) * buff

        mat(:,1,1,3) = mat(:,1,1,3) + b3 * dtl(:) * buff
        mat(:,2,2,3) = mat(:,2,2,3) + b3 * dtl(:) * buff
        mat(:,3,3,3) = mat(:,3,3,3) + b3 * dtl(:) * buff
        mat(:,4,4,3) = mat(:,4,4,3) + b3 * dtl(:) * buff
        mat(:,5,5,3) = mat(:,5,5,3) + b3 * dtl(:) * buff

        end if
        
        return
        end

!=============================================================================!
        subroutine lhs2( mat, Bh, Dh, Vh22, bc, dtl, v)
!  
!  Form the LHS for the \eta direction
!
!=============================================================================!
        use global
        use buff_mod
        implicit none
        
        real :: mat(ny*nx,ndof,ndof,3), v(ny*nx,ndof)
        real :: Bh(ny*nx,ndof,ndof), Dh(ny*nx,ndof,ndof), Vh22(ny*nx,6)
        real :: bc(nx*ndof*ndof*14)
        real :: dtl(ny*nx)

!.... first derivative operator

        real :: a1, a2, a3
        
        integer :: lrec, ier, istat
        integer :: idof, jdof
!=============================================================================!

!.... second-order stencil

        a1 = -alfa * pt5 / deta
        a2 =  alfa * zero
        a3 =  alfa * pt5 / deta

!.... \hat{B} term

        do idof = 1, ndof
          do jdof = 1, ndof
            mat(:,idof,jdof,1) = a1 * dtl(:) * Bh(:,idof,jdof)
            mat(:,idof,jdof,2) = a2 * dtl(:) * Bh(:,idof,jdof)
            mat(:,idof,jdof,3) = a3 * dtl(:) * Bh(:,idof,jdof)
          end do
        end do

!.... \hat{V}_{\eta\eta} and I terms

        call lhs2l( mat, Vh22, dtl, buff )

!.... apply boundary treatment to the LHS

        call lhsbt2( mat, Bh, Dh, Vh22, bc, spg, spg2, dtl )

!.... apply boundary conditions to the LHS

        call lhsbc2( mat, bc, v )

        return
        end

!=============================================================================!
        subroutine lhs2l( mat, Vh22, dtl, buff )
!  
!  Form the LHS for the \eta direction
!
!=============================================================================!
        use global
        implicit none
        
        real :: mat(ny*nx,ndof,ndof,3), buff(ny*nx)
        real :: dtl(ny*nx)
        
        real :: Vh22(ny*nx,6)

!.... second derivative operator

        real :: b1, b2, b3

        integer :: lrec, ier, istat

        real :: eps_i
!=============================================================================!

!.... second-order stencil

        b1 =  alfa / deta**2
        b2 = -alfa * two / deta**2
        b3 =  alfa / deta**2

!.... \hat{V}_{\eta\eta} term

        mat(:,2,2,1) = mat(:,2,2,1) - b1 * dtl(:) * Vh22(:,1)
        mat(:,2,3,1) = mat(:,2,3,1) - b1 * dtl(:) * Vh22(:,5)
        mat(:,3,2,1) = mat(:,3,2,1) - b1 * dtl(:) * Vh22(:,6)
        mat(:,3,3,1) = mat(:,3,3,1) - b1 * dtl(:) * Vh22(:,2)
        mat(:,4,4,1) = mat(:,4,4,1) - b1 * dtl(:) * Vh22(:,3)
        mat(:,5,5,1) = mat(:,5,5,1) - b1 * dtl(:) * Vh22(:,4)

        mat(:,2,2,2) = mat(:,2,2,2) - b2 * dtl(:) * Vh22(:,1)
        mat(:,2,3,2) = mat(:,2,3,2) - b2 * dtl(:) * Vh22(:,5)
        mat(:,3,2,2) = mat(:,3,2,2) - b2 * dtl(:) * Vh22(:,6)
        mat(:,3,3,2) = mat(:,3,3,2) - b2 * dtl(:) * Vh22(:,2)
        mat(:,4,4,2) = mat(:,4,4,2) - b2 * dtl(:) * Vh22(:,3)
        mat(:,5,5,2) = mat(:,5,5,2) - b2 * dtl(:) * Vh22(:,4)

        mat(:,2,2,3) = mat(:,2,2,3) - b3 * dtl(:) * Vh22(:,1)
        mat(:,2,3,3) = mat(:,2,3,3) - b3 * dtl(:) * Vh22(:,5)
        mat(:,3,2,3) = mat(:,3,2,3) - b3 * dtl(:) * Vh22(:,6)
        mat(:,3,3,3) = mat(:,3,3,3) - b3 * dtl(:) * Vh22(:,2)
        mat(:,4,4,3) = mat(:,4,4,3) - b3 * dtl(:) * Vh22(:,3)
        mat(:,5,5,3) = mat(:,5,5,3) - b3 * dtl(:) * Vh22(:,4)

!.... I term

        mat(:,1,1,2) = mat(:,1,1,2) + one
        mat(:,2,2,2) = mat(:,2,2,2) + one
        mat(:,3,3,2) = mat(:,3,3,2) + one
        mat(:,4,4,2) = mat(:,4,4,2) + one
        mat(:,5,5,2) = mat(:,5,5,2) + one

!.... implicit damping term

        if (eps_e .ne. zero) then

        eps_i = four * eps_e
        
        b1 = -alfa * eps_i
        b2 =  alfa * two * eps_i
        b3 = -alfa * eps_i

        mat(:,1,1,1) = mat(:,1,1,1) + b1 * dtl(:) * buff
        mat(:,2,2,1) = mat(:,2,2,1) + b1 * dtl(:) * buff
        mat(:,3,3,1) = mat(:,3,3,1) + b1 * dtl(:) * buff
        mat(:,4,4,1) = mat(:,4,4,1) + b1 * dtl(:) * buff
        mat(:,5,5,1) = mat(:,5,5,1) + b1 * dtl(:) * eps_i * buff

        mat(:,1,1,2) = mat(:,1,1,2) + b2 * dtl(:) * buff
        mat(:,2,2,2) = mat(:,2,2,2) + b2 * dtl(:) * buff
        mat(:,3,3,2) = mat(:,3,3,2) + b2 * dtl(:) * buff
        mat(:,4,4,2) = mat(:,4,4,2) + b2 * dtl(:) * buff
        mat(:,5,5,2) = mat(:,5,5,2) + b2 * dtl(:) * buff

        mat(:,1,1,3) = mat(:,1,1,3) + b3 * dtl(:) * buff
        mat(:,2,2,3) = mat(:,2,2,3) + b3 * dtl(:) * buff
        mat(:,3,3,3) = mat(:,3,3,3) + b3 * dtl(:) * buff
        mat(:,4,4,3) = mat(:,4,4,3) + b3 * dtl(:) * buff
        mat(:,5,5,3) = mat(:,5,5,3) + b3 * dtl(:) * buff

        end if
        
!.... sponge

        if (ispg .eq. 1) then
        
          mat(:,1,1,2) = mat(:,1,1,2) + alfa * dtl(:) * spg(:)
          mat(:,2,2,2) = mat(:,2,2,2) + alfa * dtl(:) * spg(:)
          mat(:,3,3,2) = mat(:,3,3,2) + alfa * dtl(:) * spg(:)
          mat(:,4,4,2) = mat(:,4,4,2) + alfa * dtl(:) * spg(:)
          mat(:,5,5,2) = mat(:,5,5,2) + alfa * dtl(:) * spg(:)

        else if (ispg .ge. 2) then
        
          mat(:,1,1,2) = mat(:,1,1,2) + alfa * dtl(:) * ( spg(:) + spg2(:) )
          mat(:,2,2,2) = mat(:,2,2,2) + alfa * dtl(:) * ( spg(:) + spg2(:) )
          mat(:,3,3,2) = mat(:,3,3,2) + alfa * dtl(:) * ( spg(:) + spg2(:) )
          mat(:,4,4,2) = mat(:,4,4,2) + alfa * dtl(:) * ( spg(:) + spg2(:) )
          mat(:,5,5,2) = mat(:,5,5,2) + alfa * dtl(:) * ( spg(:) + spg2(:) )

        end if

        return
        end
