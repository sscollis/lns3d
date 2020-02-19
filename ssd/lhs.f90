!=============================================================================!
!
!  These routines generate the factored LHS for the Navier-Stokes Equations
!  using 2nd order differencing
!
!=============================================================================!
        subroutine lhs1( mat, Q1, Q2, Vh1, bc, dtl, v)
!  
!  Form the LHS for the \xi direction
!
!=============================================================================!
        use stuff
        use buff_stuff
        implicit none
        
        real :: mat(ny*nx,ndof,ndof,3), v(ny*nx,ndof)
        real :: Q1(ny*nx,ndof,ndof), Q2(ny*nx,ndof,ndof), Vh1(ny*nx,6)
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

        lrec = ny * nx * ndof * ndof
        call readdr(MATRIX, Q1, lrec, 1, ier)           ! \hat{A}
        call waitdr(MATRIX, istat, ier)
        call readdr(MATRIX, Q2, lrec, 3, ier)           ! \hat{D}
        
!.... \hat{A} term

        do idof = 1, ndof
          do jdof = 1, ndof
            mat(:,idof,jdof,1) = a1 * dtl(:) * Q1(:,idof,jdof)  
            mat(:,idof,jdof,2) = a2 * dtl(:) * Q1(:,idof,jdof)  
            mat(:,idof,jdof,3) = a3 * dtl(:) * Q1(:,idof,jdof)
          end do
        end do

        call waitdr(MATRIX, istat, ier)
        lrec = ny * nx * 6
        call readdr(MATRIX, Vh1, lrec, 4, ier)          ! \hat{V}_{\xi\xi}

!.... \hat{D} term

        do idof = 1, ndof
          do jdof = 1, ndof
            mat(:,idof,jdof,2) = mat(:,idof,jdof,2) + &
                                 alfa * dtl(:) * Q2(:,idof,jdof)
          end do
        end do 

        call waitdr(MATRIX, istat, ier)

!.... \hat{V}_{\xi\xi} term and I term

        call lhs1l( mat, Vh1, dtl, buff )

!.... apply boundary treatment to the LHS

        call lhsbt1( mat, Q1, Q2, Vh1, bc, dtl )

!.... apply boundary conditions to the LHS

        call lhsbc1( mat, bc, v, vm )
        
        return
        end

!=============================================================================!
        subroutine lhs1l( mat, Vh1, dtl, buff )
!  
!  Form the LHS for the \xi direction
!
!=============================================================================!
        use stuff
        implicit none
        
        real :: mat(ny*nx,ndof,ndof,3), buff(ny*nx)
        real :: dtl(ny*nx)
        
        real :: Vh1(ny*nx,6)

!.... second derivative operator

        real :: b1, b2, b3

        real :: eps_i
!=============================================================================!

!.... second-order stencil

        b1 =  alfa / dxi**2
        b2 = -alfa * two / dxi**2
        b3 =  alfa / dxi**2
        
!.... \hat{V}_{\xi\xi} term

        mat(:,2,2,1) = mat(:,2,2,1) - b1 * dtl(:) * Vh1(:,1)
        mat(:,2,3,1) = mat(:,2,3,1) - b1 * dtl(:) * Vh1(:,5)
        mat(:,3,2,1) = mat(:,3,2,1) - b1 * dtl(:) * Vh1(:,6)
        mat(:,3,3,1) = mat(:,3,3,1) - b1 * dtl(:) * Vh1(:,2)
        mat(:,4,4,1) = mat(:,4,4,1) - b1 * dtl(:) * Vh1(:,3)
        mat(:,5,5,1) = mat(:,5,5,1) - b1 * dtl(:) * Vh1(:,4)

        mat(:,2,2,2) = mat(:,2,2,2) - b2 * dtl(:) * Vh1(:,1)
        mat(:,2,3,2) = mat(:,2,3,2) - b2 * dtl(:) * Vh1(:,5)
        mat(:,3,2,2) = mat(:,3,2,2) - b2 * dtl(:) * Vh1(:,6)
        mat(:,3,3,2) = mat(:,3,3,2) - b2 * dtl(:) * Vh1(:,2)
        mat(:,4,4,2) = mat(:,4,4,2) - b2 * dtl(:) * Vh1(:,3)
        mat(:,5,5,2) = mat(:,5,5,2) - b2 * dtl(:) * Vh1(:,4)

        mat(:,2,2,3) = mat(:,2,2,3) - b3 * dtl(:) * Vh1(:,1)
        mat(:,2,3,3) = mat(:,2,3,3) - b3 * dtl(:) * Vh1(:,5)
        mat(:,3,2,3) = mat(:,3,2,3) - b3 * dtl(:) * Vh1(:,6)
        mat(:,3,3,3) = mat(:,3,3,3) - b3 * dtl(:) * Vh1(:,2)
        mat(:,4,4,3) = mat(:,4,4,3) - b3 * dtl(:) * Vh1(:,3)
        mat(:,5,5,3) = mat(:,5,5,3) - b3 * dtl(:) * Vh1(:,4)

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
        subroutine lhs2( mat, Q1, Q2, Vh1, bc, dtl, v)
!  
!  Form the LHS for the \eta direction
!
!=============================================================================!
        use stuff
        use buff_stuff
        implicit none
        
        real :: mat(ny*nx,ndof,ndof,3), v(ny*nx,ndof)
        real :: Q1(ny*nx,ndof,ndof), Q2(ny*nx,ndof,ndof), Vh1(ny*nx,6)
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

        lrec = ny * nx * ndof * ndof
        call readdr(MATRIX, Q1, lrec, 2, ier)           ! \hat{B}
        call waitdr(MATRIX, istat, ier)
        lrec = ny * nx * 6
        call readdr(MATRIX, Vh1, lrec, 6, ier)          ! \hat{V}_{\eta\eta}
        
!.... \hat{B} term

        do idof = 1, ndof
          do jdof = 1, ndof
            mat(:,idof,jdof,1) = a1 * dtl(:) * Q1(:,idof,jdof)  
            mat(:,idof,jdof,2) = a2 * dtl(:) * Q1(:,idof,jdof)  
            mat(:,idof,jdof,3) = a3 * dtl(:) * Q1(:,idof,jdof)
          end do
        end do

        call waitdr(MATRIX, istat, ier)

!.... \hat{V}_{\eta\eta} and I terms

        call lhs2l( mat, Vh1, dtl, buff )

!.... apply boundary treatment to the LHS

        call lhsbt2( mat, Q1, Q2, Vh1, bc, spg, spg2, dtl )

!.... apply boundary conditions to the LHS

        call lhsbc2( mat, bc, v )

        return
        end

!=============================================================================!
        subroutine lhs2l( mat, Vh1, dtl, buff )
!  
!  Form the LHS for the \eta direction
!
!=============================================================================!
        use stuff
        implicit none
        
        real :: mat(ny*nx,ndof,ndof,3), buff(ny*nx)
        real :: dtl(ny*nx)
        
        real :: Vh1(ny*nx,6)

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

        mat(:,2,2,1) = mat(:,2,2,1) - b1 * dtl(:) * Vh1(:,1)
        mat(:,2,3,1) = mat(:,2,3,1) - b1 * dtl(:) * Vh1(:,5)
        mat(:,3,2,1) = mat(:,3,2,1) - b1 * dtl(:) * Vh1(:,6)
        mat(:,3,3,1) = mat(:,3,3,1) - b1 * dtl(:) * Vh1(:,2)
        mat(:,4,4,1) = mat(:,4,4,1) - b1 * dtl(:) * Vh1(:,3)
        mat(:,5,5,1) = mat(:,5,5,1) - b1 * dtl(:) * Vh1(:,4)

        mat(:,2,2,2) = mat(:,2,2,2) - b2 * dtl(:) * Vh1(:,1)
        mat(:,2,3,2) = mat(:,2,3,2) - b2 * dtl(:) * Vh1(:,5)
        mat(:,3,2,2) = mat(:,3,2,2) - b2 * dtl(:) * Vh1(:,6)
        mat(:,3,3,2) = mat(:,3,3,2) - b2 * dtl(:) * Vh1(:,2)
        mat(:,4,4,2) = mat(:,4,4,2) - b2 * dtl(:) * Vh1(:,3)
        mat(:,5,5,2) = mat(:,5,5,2) - b2 * dtl(:) * Vh1(:,4)

        mat(:,2,2,3) = mat(:,2,2,3) - b3 * dtl(:) * Vh1(:,1)
        mat(:,2,3,3) = mat(:,2,3,3) - b3 * dtl(:) * Vh1(:,5)
        mat(:,3,2,3) = mat(:,3,2,3) - b3 * dtl(:) * Vh1(:,6)
        mat(:,3,3,3) = mat(:,3,3,3) - b3 * dtl(:) * Vh1(:,2)
        mat(:,4,4,3) = mat(:,4,4,3) - b3 * dtl(:) * Vh1(:,3)
        mat(:,5,5,3) = mat(:,5,5,3) - b3 * dtl(:) * Vh1(:,4)

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
