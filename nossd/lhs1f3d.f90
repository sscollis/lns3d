!=============================================================================!
        subroutine lhs1f3D( mat, Ah, Dh, Dhi, Vh11, ABhi, dtl, v )
!  
!  Form the LHS for the \xi direction for the three-dimensional
!  Navier-Stokes Equations using 4th order differencing.
!
!=============================================================================!
        use global
        use buff_mod
        use stencil
        implicit none
        
        complex :: mat(ny*nx,ndof,ndof,5), v(ny*nx,ndof)
        real    :: Ah(ny*nx,ndof,ndof), Dh(ny*nx,ndof,ndof)
        real    :: Dhi(ny*nx,ndof,ndof), Vh11(ny*nx,6), ABhi(ny*nx,6)
        real    :: dtl(ny*nx)

!.... first derivative operator

        real    :: a1, a2, a3, a4, a5
        complex :: c1, c2, c3, c4, c5
        
        integer :: lrec, ier, istat
        integer :: idof, jdof
        logical :: calcd

!=============================================================================!
!       if ( mod(real(istep-1+iter),two) .eq. zero ) then
          calcd = .true.
!       else
!         calcd = .false.
!       end if
!       write(88,*) 'lhs1f3D: ',istep, iter, calcd
        
!.... fourth-order stencil

        a1 = alfa * ga1  / dxi
        a2 = alfa * ga2  / dxi
        a3 = alfa * zero / dxi
        a4 = alfa * ga3  / dxi
        a5 = alfa * ga4  / dxi

!.... \hat{A} term

        do idof = 1, ndof
          do jdof = 1, ndof
            mat(:,idof,jdof,1) = a1 * dtl(:) * Ah(:,idof,jdof)  
            mat(:,idof,jdof,2) = a2 * dtl(:) * Ah(:,idof,jdof)  
            mat(:,idof,jdof,3) = a3 * dtl(:) * Ah(:,idof,jdof)
            mat(:,idof,jdof,4) = a4 * dtl(:) * Ah(:,idof,jdof)
            mat(:,idof,jdof,5) = a5 * dtl(:) * Ah(:,idof,jdof)
          end do
        end do

!.... fourth-order stencil

        c1 = im * alfa * ga1  / dxi
        c2 = im * alfa * ga2  / dxi
        c3 = im * alfa * zero / dxi
        c4 = im * alfa * ga3  / dxi
        c5 = im * alfa * ga4  / dxi

!.... \hat{A}_i term

        mat(:,2,4,1) = mat(:,2,4,1) - c1 * dtl(:) * ABhi(:,1)
        mat(:,3,4,1) = mat(:,3,4,1) - c1 * dtl(:) * ABhi(:,2)
        mat(:,4,2,1) = mat(:,4,2,1) - c1 * dtl(:) * ABhi(:,1)
        mat(:,4,3,1) = mat(:,4,3,1) - c1 * dtl(:) * ABhi(:,2)

        mat(:,2,4,2) = mat(:,2,4,2) - c2 * dtl(:) * ABhi(:,1)
        mat(:,3,4,2) = mat(:,3,4,2) - c2 * dtl(:) * ABhi(:,2)
        mat(:,4,2,2) = mat(:,4,2,2) - c2 * dtl(:) * ABhi(:,1)
        mat(:,4,3,2) = mat(:,4,3,2) - c2 * dtl(:) * ABhi(:,2)

        mat(:,2,4,3) = mat(:,2,4,3) - c3 * dtl(:) * ABhi(:,1)
        mat(:,3,4,3) = mat(:,3,4,3) - c3 * dtl(:) * ABhi(:,2)
        mat(:,4,2,3) = mat(:,4,2,3) - c3 * dtl(:) * ABhi(:,1)
        mat(:,4,3,3) = mat(:,4,3,3) - c3 * dtl(:) * ABhi(:,2)

        mat(:,2,4,4) = mat(:,2,4,4) - c4 * dtl(:) * ABhi(:,1)
        mat(:,3,4,4) = mat(:,3,4,4) - c4 * dtl(:) * ABhi(:,2)
        mat(:,4,2,4) = mat(:,4,2,4) - c4 * dtl(:) * ABhi(:,1)
        mat(:,4,3,4) = mat(:,4,3,4) - c4 * dtl(:) * ABhi(:,2)

        mat(:,2,4,5) = mat(:,2,4,5) - c5 * dtl(:) * ABhi(:,1)
        mat(:,3,4,5) = mat(:,3,4,5) - c5 * dtl(:) * ABhi(:,2)
        mat(:,4,2,5) = mat(:,4,2,5) - c5 * dtl(:) * ABhi(:,1)
        mat(:,4,3,5) = mat(:,4,3,5) - c5 * dtl(:) * ABhi(:,2)

        if (calcd) then

!.... \hat{D} term

          do idof = 1, ndof
            do jdof = 1, ndof
              mat(:,idof,jdof,3) = mat(:,idof,jdof,3) + &
                                   alfa * dtl(:) * Dh(:,idof,jdof)
            end do
          end do

!.... \hat{D}_i term

          do idof = 1, ndof
            do jdof = 1, ndof
              mat(:,idof,jdof,3) = mat(:,idof,jdof,3) + &
                                   im * alfa * dtl(:) * Dhi(:,idof,jdof)
            end do
          end do

        end if

!.... \hat{V}_{\xi\xi} term and I term

        call lhs1lf3D( mat, Vh11, dtl, buff, calcd )

!.... apply boundary treatment to the LHS

        call lhsbt1f3D( mat, Ah, Dh, Dhi, Vh11, ABhi, spg, spg2, dtl, calcd )

!.... apply boundary conditions to the LHS

        call lhsbc1f3D( mat, v, vm )
        
        return
        end

!=============================================================================!
        subroutine lhs1lf3D( mat, Vh11, dtl, buff, calcd )
!  
!  Form the LHS for the \xi direction
!
!=============================================================================!
        use global
        use stencil
        implicit none
        
        complex :: mat(ny*nx,ndof,ndof,5)
        real    :: Vh11(ny*nx,6), dtl(ny*nx), buff(ny*nx)
        logical :: calcd

!.... second derivative operator

        real :: b1, b2, b3, b4, b5
!=============================================================================!
        
!.... fourth-order stencil

        b1 = alfa * da1 / dxi**2
        b2 = alfa * da2 / dxi**2
        b3 = alfa * da3 / dxi**2
        b4 = alfa * da4 / dxi**2
        b5 = alfa * da5 / dxi**2
        
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

        mat(:,2,2,4) = mat(:,2,2,4) - b4 * dtl(:) * Vh11(:,1)
        mat(:,2,3,4) = mat(:,2,3,4) - b4 * dtl(:) * Vh11(:,5)
        mat(:,3,2,4) = mat(:,3,2,4) - b4 * dtl(:) * Vh11(:,6)
        mat(:,3,3,4) = mat(:,3,3,4) - b4 * dtl(:) * Vh11(:,2)
        mat(:,4,4,4) = mat(:,4,4,4) - b4 * dtl(:) * Vh11(:,3)
        mat(:,5,5,4) = mat(:,5,5,4) - b4 * dtl(:) * Vh11(:,4)

        mat(:,2,2,5) = mat(:,2,2,5) - b5 * dtl(:) * Vh11(:,1)
        mat(:,2,3,5) = mat(:,2,3,5) - b5 * dtl(:) * Vh11(:,5)
        mat(:,3,2,5) = mat(:,3,2,5) - b5 * dtl(:) * Vh11(:,6)
        mat(:,3,3,5) = mat(:,3,3,5) - b5 * dtl(:) * Vh11(:,2)
        mat(:,4,4,5) = mat(:,4,4,5) - b5 * dtl(:) * Vh11(:,3)
        mat(:,5,5,5) = mat(:,5,5,5) - b5 * dtl(:) * Vh11(:,4)

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
        
!.... sponge

        if (.not. calcd) then

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

        end if

        return
        end
