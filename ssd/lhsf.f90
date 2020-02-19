!=============================================================================!
!
!  These routines generate the factored LHS for the Navier-Stokes Equations
!  using 4th order differencing.
!
!=============================================================================!
        subroutine lhs1f( mat, Q1, Q2, Vh1, dtl, v )
!  
!  Form the LHS for the \xi direction
!
!=============================================================================!
        use stuff
        use buff_stuff
        use diff
        implicit none
        
        real :: mat(ny*nx,ndof,ndof,5), v(ny*nx,ndof)
        real :: Q1(ny*nx,ndof,ndof), Q2(ny*nx,ndof,ndof), Vh1(ny*nx,6)
        real :: dtl(ny*nx)

!.... first derivative operator

        real :: a1, a2, a3, a4, a5

        integer :: lrec, ier, istat
        integer :: idof, jdof
!=============================================================================!

!.... fourth-order stencil

        a1 = alfa * ga1  / dxi
        a2 = alfa * ga2  / dxi
        a3 = alfa * zero / dxi
        a4 = alfa * ga3  / dxi
        a5 = alfa * ga4  / dxi

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
            mat(:,idof,jdof,4) = a4 * dtl(:) * Q1(:,idof,jdof)
            mat(:,idof,jdof,5) = a5 * dtl(:) * Q1(:,idof,jdof)
          end do
        end do

        call waitdr(MATRIX, istat, ier)
        lrec = ny * nx * 6
        call readdr(MATRIX, Vh1, lrec, 4, ier)          ! \hat{V}_{\xi\xi}
        
!.... \hat{D} term

        do idof = 1, ndof
          do jdof = 1, ndof
            mat(:,idof,jdof,3) = mat(:,idof,jdof,3) + &
                                 alfa * dtl(:) * Q2(:,idof,jdof)
          end do
        end do 

        call waitdr(MATRIX, istat, ier)

!.... \hat{V}_{\xi\xi} term and I term

        call lhs1lf( mat, Vh1, dtl, buff )

!.... apply boundary treatment to the LHS

        call lhsbt1f( mat, Q1, Q2, Vh1, dtl )

!.... apply boundary conditions to the LHS

        call lhsbc1f( mat, v, vm )
        
        return
        end

!=============================================================================!
        subroutine lhs1lf( mat, Vh1, dtl, buff )
!  
!  Form the LHS for the \xi direction
!
!=============================================================================!
        use stuff
        use diff
        implicit none
        
        real :: mat(ny*nx,ndof,ndof,5), buff(ny*nx)
        real :: dtl(ny*nx)
        
        real :: Vh1(ny*nx,6)

!.... second derivative operator

        real :: b1, b2, b3, b4, b5
!=============================================================================!
        
!.... fourth-order stencil

        b1 =  alfa * da1 / dxi**2
        b2 =  alfa * da2 / dxi**2
        b3 =  alfa * da3 / dxi**2
        b4 =  alfa * da4 / dxi**2
        b5 =  alfa * da5 / dxi**2
        
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

        mat(:,2,2,4) = mat(:,2,2,4) - b4 * dtl(:) * Vh1(:,1)
        mat(:,2,3,4) = mat(:,2,3,4) - b4 * dtl(:) * Vh1(:,5)
        mat(:,3,2,4) = mat(:,3,2,4) - b4 * dtl(:) * Vh1(:,6)
        mat(:,3,3,4) = mat(:,3,3,4) - b4 * dtl(:) * Vh1(:,2)
        mat(:,4,4,4) = mat(:,4,4,4) - b4 * dtl(:) * Vh1(:,3)
        mat(:,5,5,4) = mat(:,5,5,4) - b4 * dtl(:) * Vh1(:,4)

        mat(:,2,2,5) = mat(:,2,2,5) - b5 * dtl(:) * Vh1(:,1)
        mat(:,2,3,5) = mat(:,2,3,5) - b5 * dtl(:) * Vh1(:,5)
        mat(:,3,2,5) = mat(:,3,2,5) - b5 * dtl(:) * Vh1(:,6)
        mat(:,3,3,5) = mat(:,3,3,5) - b5 * dtl(:) * Vh1(:,2)
        mat(:,4,4,5) = mat(:,4,4,5) - b5 * dtl(:) * Vh1(:,3)
        mat(:,5,5,5) = mat(:,5,5,5) - b5 * dtl(:) * Vh1(:,4)

!.... I term

        mat(:,1,1,3) = mat(:,1,1,3) + one
        mat(:,2,2,3) = mat(:,2,2,3) + one
        mat(:,3,3,3) = mat(:,3,3,3) + one
        mat(:,4,4,3) = mat(:,4,4,3) + one
        mat(:,5,5,3) = mat(:,5,5,3) + one
                        
!.... implicit damping term updated to fourth order

        if (eps_e .ne. zero) then

        mat(:,1,1,1) = mat(:,1,1,1) + alfa * dtl(:) * eps_e * buff * fb1
        mat(:,2,2,1) = mat(:,2,2,1) + alfa * dtl(:) * eps_e * buff * fb1
        mat(:,3,3,1) = mat(:,3,3,1) + alfa * dtl(:) * eps_e * buff * fb1
        mat(:,4,4,1) = mat(:,4,4,1) + alfa * dtl(:) * eps_e * buff * fb1
        mat(:,5,5,1) = mat(:,5,5,1) + alfa * dtl(:) * eps_e * buff * fb1

        mat(:,1,1,2) = mat(:,1,1,2) + alfa * dtl(:) * eps_e * buff * fb2
        mat(:,2,2,2) = mat(:,2,2,2) + alfa * dtl(:) * eps_e * buff * fb2
        mat(:,3,3,2) = mat(:,3,3,2) + alfa * dtl(:) * eps_e * buff * fb2
        mat(:,4,4,2) = mat(:,4,4,2) + alfa * dtl(:) * eps_e * buff * fb2
        mat(:,5,5,2) = mat(:,5,5,2) + alfa * dtl(:) * eps_e * buff * fb2

        mat(:,1,1,3) = mat(:,1,1,3) + alfa * dtl(:) * eps_e * buff * fb3
        mat(:,2,2,3) = mat(:,2,2,3) + alfa * dtl(:) * eps_e * buff * fb3
        mat(:,3,3,3) = mat(:,3,3,3) + alfa * dtl(:) * eps_e * buff * fb3
        mat(:,4,4,3) = mat(:,4,4,3) + alfa * dtl(:) * eps_e * buff * fb3
        mat(:,5,5,3) = mat(:,5,5,3) + alfa * dtl(:) * eps_e * buff * fb3

        mat(:,1,1,4) = mat(:,1,1,4) + alfa * dtl(:) * eps_e * buff * fb4
        mat(:,2,2,4) = mat(:,2,2,4) + alfa * dtl(:) * eps_e * buff * fb4
        mat(:,3,3,4) = mat(:,3,3,4) + alfa * dtl(:) * eps_e * buff * fb4
        mat(:,4,4,4) = mat(:,4,4,4) + alfa * dtl(:) * eps_e * buff * fb4
        mat(:,5,5,4) = mat(:,5,5,4) + alfa * dtl(:) * eps_e * buff * fb4

        mat(:,1,1,5) = mat(:,1,1,5) + alfa * dtl(:) * eps_e * buff * fb5
        mat(:,2,2,5) = mat(:,2,2,5) + alfa * dtl(:) * eps_e * buff * fb5
        mat(:,3,3,5) = mat(:,3,3,5) + alfa * dtl(:) * eps_e * buff * fb5
        mat(:,4,4,5) = mat(:,4,4,5) + alfa * dtl(:) * eps_e * buff * fb5
        mat(:,5,5,5) = mat(:,5,5,5) + alfa * dtl(:) * eps_e * buff * fb5

        end if
        
        return
        end

!=============================================================================!
        subroutine lhs2f( mat, Q1, Q2, Vh1, dtl, v)
!  
!  Form the LHS for the \eta direction
!
!=============================================================================!
        use stuff
        use buff_stuff
        use diff
        implicit none
        
        real :: mat(ny*nx,ndof,ndof,5), v(ny*nx,ndof)
        real :: Q1(ny*nx,ndof,ndof), Q2(ny*nx,ndof,ndof), Vh1(ny*nx,6)
        real :: dtl(ny*nx)

!.... first derivative operator

        real :: a1, a2, a3, a4, a5
        
        integer :: lrec, ier, istat
        integer :: idof, jdof
!=============================================================================!

!.... fourth-order stencil

        a1 = alfa * ga1 / deta
        a2 = alfa * ga2 / deta
        a3 = alfa * zero / deta
        a4 = alfa * ga3 / deta
        a5 = alfa * ga4 / deta

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
            mat(:,idof,jdof,4) = a4 * dtl(:) * Q1(:,idof,jdof)
            mat(:,idof,jdof,5) = a5 * dtl(:) * Q1(:,idof,jdof)
          end do
        end do

        call waitdr(MATRIX, istat, ier)

!.... \hat{V}_{\eta\eta} and I terms

        call lhs2lf( mat, Vh1, dtl, buff )

!.... apply boundary treatment to the LHS

        call lhsbt2f( mat, Q1, Q2, Vh1, spg, spg2, dtl )

!.... apply boundary conditions to the LHS

        call lhsbc2f( mat, v )

        return
        end

!=============================================================================!
        subroutine lhs2lf( mat, Vh1, dtl, buff )
!  
!  Form the LHS for the \eta direction
!
!=============================================================================!
        use stuff
        use diff
        implicit none
        
        real :: mat(ny*nx,ndof,ndof,5), buff(ny*nx)
        real :: dtl(ny*nx)
        
        real :: Vh1(ny*nx,6)

!.... second derivative operator

        real :: b1, b2, b3, b4, b5

        integer :: lrec, ier, istat
!=============================================================================!

!.... fourth-order stencil

        b1 = alfa * da1 / deta**2
        b2 = alfa * da2 / deta**2
        b3 = alfa * da3 / deta**2
        b4 = alfa * da4 / deta**2
        b5 = alfa * da5 / deta**2

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

        mat(:,2,2,4) = mat(:,2,2,4) - b4 * dtl(:) * Vh1(:,1)
        mat(:,2,3,4) = mat(:,2,3,4) - b4 * dtl(:) * Vh1(:,5)
        mat(:,3,2,4) = mat(:,3,2,4) - b4 * dtl(:) * Vh1(:,6)
        mat(:,3,3,4) = mat(:,3,3,4) - b4 * dtl(:) * Vh1(:,2)
        mat(:,4,4,4) = mat(:,4,4,4) - b4 * dtl(:) * Vh1(:,3)
        mat(:,5,5,4) = mat(:,5,5,4) - b4 * dtl(:) * Vh1(:,4)

        mat(:,2,2,5) = mat(:,2,2,5) - b5 * dtl(:) * Vh1(:,1)
        mat(:,2,3,5) = mat(:,2,3,5) - b5 * dtl(:) * Vh1(:,5)
        mat(:,3,2,5) = mat(:,3,2,5) - b5 * dtl(:) * Vh1(:,6)
        mat(:,3,3,5) = mat(:,3,3,5) - b5 * dtl(:) * Vh1(:,2)
        mat(:,4,4,5) = mat(:,4,4,5) - b5 * dtl(:) * Vh1(:,3)
        mat(:,5,5,5) = mat(:,5,5,5) - b5 * dtl(:) * Vh1(:,4)

!.... I term

        mat(:,1,1,3) = mat(:,1,1,3) + one
        mat(:,2,2,3) = mat(:,2,2,3) + one
        mat(:,3,3,3) = mat(:,3,3,3) + one
        mat(:,4,4,3) = mat(:,4,4,3) + one
        mat(:,5,5,3) = mat(:,5,5,3) + one

!.... implicit damping term updated to fourth-order

        if (eps_e .ne. zero) then

        mat(:,1,1,1) = mat(:,1,1,1) + alfa * dtl(:) * eps_e * buff * fb1
        mat(:,2,2,1) = mat(:,2,2,1) + alfa * dtl(:) * eps_e * buff * fb1
        mat(:,3,3,1) = mat(:,3,3,1) + alfa * dtl(:) * eps_e * buff * fb1
        mat(:,4,4,1) = mat(:,4,4,1) + alfa * dtl(:) * eps_e * buff * fb1
        mat(:,5,5,1) = mat(:,5,5,1) + alfa * dtl(:) * eps_e * buff * fb1

        mat(:,1,1,2) = mat(:,1,1,2) + alfa * dtl(:) * eps_e * buff * fb2
        mat(:,2,2,2) = mat(:,2,2,2) + alfa * dtl(:) * eps_e * buff * fb2
        mat(:,3,3,2) = mat(:,3,3,2) + alfa * dtl(:) * eps_e * buff * fb2
        mat(:,4,4,2) = mat(:,4,4,2) + alfa * dtl(:) * eps_e * buff * fb2
        mat(:,5,5,2) = mat(:,5,5,2) + alfa * dtl(:) * eps_e * buff * fb2

        mat(:,1,1,3) = mat(:,1,1,3) + alfa * dtl(:) * eps_e * buff * fb3
        mat(:,2,2,3) = mat(:,2,2,3) + alfa * dtl(:) * eps_e * buff * fb3
        mat(:,3,3,3) = mat(:,3,3,3) + alfa * dtl(:) * eps_e * buff * fb3
        mat(:,4,4,3) = mat(:,4,4,3) + alfa * dtl(:) * eps_e * buff * fb3
        mat(:,5,5,3) = mat(:,5,5,3) + alfa * dtl(:) * eps_e * buff * fb3

        mat(:,1,1,4) = mat(:,1,1,4) + alfa * dtl(:) * eps_e * buff * fb4
        mat(:,2,2,4) = mat(:,2,2,4) + alfa * dtl(:) * eps_e * buff * fb4
        mat(:,3,3,4) = mat(:,3,3,4) + alfa * dtl(:) * eps_e * buff * fb4
        mat(:,4,4,4) = mat(:,4,4,4) + alfa * dtl(:) * eps_e * buff * fb4
        mat(:,5,5,4) = mat(:,5,5,4) + alfa * dtl(:) * eps_e * buff * fb4

        mat(:,1,1,5) = mat(:,1,1,5) + alfa * dtl(:) * eps_e * buff * fb5
        mat(:,2,2,5) = mat(:,2,2,5) + alfa * dtl(:) * eps_e * buff * fb5
        mat(:,3,3,5) = mat(:,3,3,5) + alfa * dtl(:) * eps_e * buff * fb5
        mat(:,4,4,5) = mat(:,4,4,5) + alfa * dtl(:) * eps_e * buff * fb5
        mat(:,5,5,5) = mat(:,5,5,5) + alfa * dtl(:) * eps_e * buff * fb5

        end if

!.... sponge

        if (ispg .eq. 1) then
        
          mat(:,1,1,3) = mat(:,1,1,3) + alfa * dtl(:) * spg(:)
          mat(:,2,2,3) = mat(:,2,2,3) + alfa * dtl(:) * spg(:)
          mat(:,3,3,3) = mat(:,3,3,3) + alfa * dtl(:) * spg(:)
          mat(:,4,4,3) = mat(:,4,4,3) + alfa * dtl(:) * spg(:)
          mat(:,5,5,3) = mat(:,5,5,3) + alfa * dtl(:) * spg(:)
        
        else if (ispg .ge. 2) then
        
          mat(:,1,1,3) = mat(:,1,1,3) + alfa * dtl(:) * ( spg(:) + spg2(:) )
          mat(:,2,2,3) = mat(:,2,2,3) + alfa * dtl(:) * ( spg(:) + spg2(:) )
          mat(:,3,3,3) = mat(:,3,3,3) + alfa * dtl(:) * ( spg(:) + spg2(:) )
          mat(:,4,4,3) = mat(:,4,4,3) + alfa * dtl(:) * ( spg(:) + spg2(:) )
          mat(:,5,5,3) = mat(:,5,5,3) + alfa * dtl(:) * ( spg(:) + spg2(:) )
        
        end if

        return
        end
