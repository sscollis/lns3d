!=======================================================================================================!
        subroutine lhsbt1( mat, Ah, Dh, Vh, bc, dtl )
!  
!  Correct the LHS for boundary treatment in the \xi direction.
!
!  This version supports the 4th order LHS. 
!
!  Revised: 10-16-95
!
!=======================================================================================================!
        use global
        use stencil
        use buff_mod
        implicit none
        
        real :: mat(ny,nx,ndof,ndof,3), Ah(ny,nx,ndof,ndof)
        real :: Dh(ny,nx,ndof,ndof), bc(ny,ndof,ndof,14)
        real :: Vh(ny,nx,6)
        real :: dtl(ny,nx)

        real :: dxiinv, dxisinv

        integer :: iv, idof, jdof
!=======================================================================================================!
        dxiinv  = one / dxi
        dxisinv = one / dxi**2

!.... initialize the bc array

        bc = zero

        if (xper) return

        if (lsym) then
        
          call lhsbs1( mat, Ah, Dh, Vh, bc, dtl, 1 )             ! apply a symmetry condition

        else            ! .not. lsym
                
!=======================================================================================================!
!.... use higher-order tangent on left node
!=======================================================================================================!
        do idof = 1, ndof
          do jdof = 1, ndof
            do iv = 1, ny
              mat(iv,1,idof,jdof,1) = zero
              mat(iv,1,idof,jdof,2) = alfa * dtl(iv,1) * gc1 * dxiinv * Ah(iv,1,idof,jdof) + &
                                      alfa * dtl(iv,1) * Dh(iv,1,idof,jdof)
              mat(iv,1,idof,jdof,3) = alfa * dtl(iv,1) * gc2 * dxiinv * Ah(iv,1,idof,jdof)
              bc(iv,idof,jdof,1)    = alfa * dtl(iv,1) * gc3 * dxiinv * Ah(iv,1,idof,jdof)
              bc(iv,idof,jdof,2)    = alfa * dtl(iv,1) * gc4 * dxiinv * Ah(iv,1,idof,jdof)
              bc(iv,idof,jdof,3)    = alfa * dtl(iv,1) * gc5 * dxiinv * Ah(iv,1,idof,jdof)
            end do
          end do
        end do
        
!.... \hat{V}_{\eta\eta} term

        do iv = 1, ny
        
        mat(iv,1,2,2,2) = mat(iv,1,2,2,2) - alfa * dtl(iv,1) * dd1 * dxisinv * Vh(iv,1,1)
        mat(iv,1,2,3,2) = mat(iv,1,2,3,2) - alfa * dtl(iv,1) * dd1 * dxisinv * Vh(iv,1,5)
        mat(iv,1,3,2,2) = mat(iv,1,3,2,2) - alfa * dtl(iv,1) * dd1 * dxisinv * Vh(iv,1,6)
        mat(iv,1,3,3,2) = mat(iv,1,3,3,2) - alfa * dtl(iv,1) * dd1 * dxisinv * Vh(iv,1,2)
        mat(iv,1,4,4,2) = mat(iv,1,4,4,2) - alfa * dtl(iv,1) * dd1 * dxisinv * Vh(iv,1,3)
        mat(iv,1,5,5,2) = mat(iv,1,5,5,2) - alfa * dtl(iv,1) * dd1 * dxisinv * Vh(iv,1,4)

        mat(iv,1,2,2,3) = mat(iv,1,2,2,3) - alfa * dtl(iv,1) * dd2 * dxisinv * Vh(iv,1,1)
        mat(iv,1,2,3,3) = mat(iv,1,2,3,3) - alfa * dtl(iv,1) * dd2 * dxisinv * Vh(iv,1,5)
        mat(iv,1,3,2,3) = mat(iv,1,3,2,3) - alfa * dtl(iv,1) * dd2 * dxisinv * Vh(iv,1,6)
        mat(iv,1,3,3,3) = mat(iv,1,3,3,3) - alfa * dtl(iv,1) * dd2 * dxisinv * Vh(iv,1,2)
        mat(iv,1,4,4,3) = mat(iv,1,4,4,3) - alfa * dtl(iv,1) * dd2 * dxisinv * Vh(iv,1,3)
        mat(iv,1,5,5,3) = mat(iv,1,5,5,3) - alfa * dtl(iv,1) * dd2 * dxisinv * Vh(iv,1,4)

        bc(iv,2,2,1)    = bc(iv,2,2,1)    - alfa * dtl(iv,1) * dd3 * dxisinv * Vh(iv,1,1)
        bc(iv,2,3,1)    = bc(iv,2,3,1)    - alfa * dtl(iv,1) * dd3 * dxisinv * Vh(iv,1,5)
        bc(iv,3,2,1)    = bc(iv,3,2,1)    - alfa * dtl(iv,1) * dd3 * dxisinv * Vh(iv,1,6)
        bc(iv,3,3,1)    = bc(iv,3,3,1)    - alfa * dtl(iv,1) * dd3 * dxisinv * Vh(iv,1,2)
        bc(iv,4,4,1)    = bc(iv,4,4,1)    - alfa * dtl(iv,1) * dd3 * dxisinv * Vh(iv,1,3)
        bc(iv,5,5,1)    = bc(iv,5,5,1)    - alfa * dtl(iv,1) * dd3 * dxisinv * Vh(iv,1,4)

        bc(iv,2,2,2)    = bc(iv,2,2,2)    - alfa * dtl(iv,1) * dd4 * dxisinv * Vh(iv,1,1)
        bc(iv,2,3,2)    = bc(iv,2,3,2)    - alfa * dtl(iv,1) * dd4 * dxisinv * Vh(iv,1,5)
        bc(iv,3,2,2)    = bc(iv,3,2,2)    - alfa * dtl(iv,1) * dd4 * dxisinv * Vh(iv,1,6)
        bc(iv,3,3,2)    = bc(iv,3,3,2)    - alfa * dtl(iv,1) * dd4 * dxisinv * Vh(iv,1,2)
        bc(iv,4,4,2)    = bc(iv,4,4,2)    - alfa * dtl(iv,1) * dd4 * dxisinv * Vh(iv,1,3)
        bc(iv,5,5,2)    = bc(iv,5,5,2)    - alfa * dtl(iv,1) * dd4 * dxisinv * Vh(iv,1,4)

        bc(iv,2,2,3)    = bc(iv,2,2,3)    - alfa * dtl(iv,1) * dd5 * dxisinv * Vh(iv,1,1)
        bc(iv,2,3,3)    = bc(iv,2,3,3)    - alfa * dtl(iv,1) * dd5 * dxisinv * Vh(iv,1,5)
        bc(iv,3,2,3)    = bc(iv,3,2,3)    - alfa * dtl(iv,1) * dd5 * dxisinv * Vh(iv,1,6)
        bc(iv,3,3,3)    = bc(iv,3,3,3)    - alfa * dtl(iv,1) * dd5 * dxisinv * Vh(iv,1,2)
        bc(iv,4,4,3)    = bc(iv,4,4,3)    - alfa * dtl(iv,1) * dd5 * dxisinv * Vh(iv,1,3)
        bc(iv,5,5,3)    = bc(iv,5,5,3)    - alfa * dtl(iv,1) * dd5 * dxisinv * Vh(iv,1,4)

!.... I term

        mat(iv,1,1,1,2) = mat(iv,1,1,1,2) + one
        mat(iv,1,2,2,2) = mat(iv,1,2,2,2) + one
        mat(iv,1,3,3,2) = mat(iv,1,3,3,2) + one
        mat(iv,1,4,4,2) = mat(iv,1,4,4,2) + one
        mat(iv,1,5,5,2) = mat(iv,1,5,5,2) + one

        end do

!=======================================================================================================!
!.... use higher-order tangent on first node off the left boundary
!=======================================================================================================!
        do idof = 1, ndof
          do jdof = 1, ndof
            do iv = 1, ny
              mat(iv,2,idof,jdof,1) = alfa * dtl(iv,2) * gb1 * dxiinv * Ah(iv,2,idof,jdof)
              mat(iv,2,idof,jdof,2) = alfa * dtl(iv,2) * gb2 * dxiinv * Ah(iv,2,idof,jdof) + &
                                      alfa * dtl(iv,2) * Dh(iv,2,idof,jdof)
              mat(iv,2,idof,jdof,3) = alfa * dtl(iv,2) * gb3 * dxiinv * Ah(iv,2,idof,jdof)
              bc(iv,idof,jdof,4)    = alfa * dtl(iv,2) * gb4 * dxiinv * Ah(iv,2,idof,jdof)
              bc(iv,idof,jdof,5)    = alfa * dtl(iv,2) * gb5 * dxiinv * Ah(iv,2,idof,jdof)
            end do
          end do
        end do
        
!.... \hat{V}_{\eta\eta} term

        do iv = 1, ny
        
        mat(iv,2,2,2,1) = mat(iv,2,2,2,1) - alfa * dtl(iv,2) * db1 * dxisinv * Vh(iv,2,1)
        mat(iv,2,2,3,1) = mat(iv,2,2,3,1) - alfa * dtl(iv,2) * db1 * dxisinv * Vh(iv,2,5)
        mat(iv,2,3,2,1) = mat(iv,2,3,2,1) - alfa * dtl(iv,2) * db1 * dxisinv * Vh(iv,2,6)
        mat(iv,2,3,3,1) = mat(iv,2,3,3,1) - alfa * dtl(iv,2) * db1 * dxisinv * Vh(iv,2,2)
        mat(iv,2,4,4,1) = mat(iv,2,4,4,1) - alfa * dtl(iv,2) * db1 * dxisinv * Vh(iv,2,3)
        mat(iv,2,5,5,1) = mat(iv,2,5,5,1) - alfa * dtl(iv,2) * db1 * dxisinv * Vh(iv,2,4)

        mat(iv,2,2,2,2) = mat(iv,2,2,2,2) - alfa * dtl(iv,2) * db2 * dxisinv * Vh(iv,2,1)
        mat(iv,2,2,3,2) = mat(iv,2,2,3,2) - alfa * dtl(iv,2) * db2 * dxisinv * Vh(iv,2,5)
        mat(iv,2,3,2,2) = mat(iv,2,3,2,2) - alfa * dtl(iv,2) * db2 * dxisinv * Vh(iv,2,6)
        mat(iv,2,3,3,2) = mat(iv,2,3,3,2) - alfa * dtl(iv,2) * db2 * dxisinv * Vh(iv,2,2)
        mat(iv,2,4,4,2) = mat(iv,2,4,4,2) - alfa * dtl(iv,2) * db2 * dxisinv * Vh(iv,2,3)
        mat(iv,2,5,5,2) = mat(iv,2,5,5,2) - alfa * dtl(iv,2) * db2 * dxisinv * Vh(iv,2,4)

        mat(iv,2,2,2,3) = mat(iv,2,2,2,3) - alfa * dtl(iv,2) * db3 * dxisinv * Vh(iv,2,1)
        mat(iv,2,2,3,3) = mat(iv,2,2,3,3) - alfa * dtl(iv,2) * db3 * dxisinv * Vh(iv,2,5)
        mat(iv,2,3,2,3) = mat(iv,2,3,2,3) - alfa * dtl(iv,2) * db3 * dxisinv * Vh(iv,2,6)
        mat(iv,2,3,3,3) = mat(iv,2,3,3,3) - alfa * dtl(iv,2) * db3 * dxisinv * Vh(iv,2,2)
        mat(iv,2,4,4,3) = mat(iv,2,4,4,3) - alfa * dtl(iv,2) * db3 * dxisinv * Vh(iv,2,3)
        mat(iv,2,5,5,3) = mat(iv,2,5,5,3) - alfa * dtl(iv,2) * db3 * dxisinv * Vh(iv,2,4)

        bc(iv,2,2,4)    = bc(iv,2,2,4)    - alfa * dtl(iv,2) * db4 * dxisinv * Vh(iv,2,1)
        bc(iv,2,3,4)    = bc(iv,2,3,4)    - alfa * dtl(iv,2) * db4 * dxisinv * Vh(iv,2,5)
        bc(iv,3,2,4)    = bc(iv,3,2,4)    - alfa * dtl(iv,2) * db4 * dxisinv * Vh(iv,2,6)
        bc(iv,3,3,4)    = bc(iv,3,3,4)    - alfa * dtl(iv,2) * db4 * dxisinv * Vh(iv,2,2)
        bc(iv,4,4,4)    = bc(iv,4,4,4)    - alfa * dtl(iv,2) * db4 * dxisinv * Vh(iv,2,3)
        bc(iv,5,5,4)    = bc(iv,5,5,4)    - alfa * dtl(iv,2) * db4 * dxisinv * Vh(iv,2,4)

        bc(iv,2,2,5)    = bc(iv,2,2,5)    - alfa * dtl(iv,2) * db5 * dxisinv * Vh(iv,2,1)
        bc(iv,2,3,5)    = bc(iv,2,3,5)    - alfa * dtl(iv,2) * db5 * dxisinv * Vh(iv,2,5)
        bc(iv,3,2,5)    = bc(iv,3,2,5)    - alfa * dtl(iv,2) * db5 * dxisinv * Vh(iv,2,6)
        bc(iv,3,3,5)    = bc(iv,3,3,5)    - alfa * dtl(iv,2) * db5 * dxisinv * Vh(iv,2,2)
        bc(iv,4,4,5)    = bc(iv,4,4,5)    - alfa * dtl(iv,2) * db5 * dxisinv * Vh(iv,2,3)
        bc(iv,5,5,5)    = bc(iv,5,5,5)    - alfa * dtl(iv,2) * db5 * dxisinv * Vh(iv,2,4)

!.... I term

        mat(iv,2,1,1,2) = mat(iv,2,1,1,2) + one
        mat(iv,2,2,2,2) = mat(iv,2,2,2,2) + one
        mat(iv,2,3,3,2) = mat(iv,2,3,3,2) + one
        mat(iv,2,4,4,2) = mat(iv,2,4,4,2) + one
        mat(iv,2,5,5,2) = mat(iv,2,5,5,2) + one

        end do

!.... implicit damping term

        if (eps_e .ne. zero) then

        do iv = 1, ny
        mat(iv,2,1,1,1) = mat(iv,2,1,1,1) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * (-one)
        mat(iv,2,2,2,1) = mat(iv,2,2,2,1) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * (-one)
        mat(iv,2,3,3,1) = mat(iv,2,3,3,1) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * (-one)
        mat(iv,2,4,4,1) = mat(iv,2,4,4,1) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * (-one)
        mat(iv,2,5,5,1) = mat(iv,2,5,5,1) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * (-one)

        mat(iv,2,1,1,2) = mat(iv,2,1,1,2) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * two
        mat(iv,2,2,2,2) = mat(iv,2,2,2,2) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * two
        mat(iv,2,3,3,2) = mat(iv,2,3,3,2) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * two
        mat(iv,2,4,4,2) = mat(iv,2,4,4,2) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * two
        mat(iv,2,5,5,2) = mat(iv,2,5,5,2) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * two

        mat(iv,2,1,1,3) = mat(iv,2,1,1,3) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * (-one)
        mat(iv,2,2,2,3) = mat(iv,2,2,2,3) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * (-one)
        mat(iv,2,3,3,3) = mat(iv,2,3,3,3) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * (-one)
        mat(iv,2,4,4,3) = mat(iv,2,4,4,3) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * (-one)
        mat(iv,2,5,5,3) = mat(iv,2,5,5,3) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * (-one)
        end do
                
        end if

        end if          ! lsym

!=======================================================================================================!
!.... use higher-order tangent on second node off the left boundary
!=======================================================================================================!
        do idof = 1, ndof
          do jdof = 1, ndof
            do iv = 1, ny
              bc(iv,idof,jdof,7)    = alfa * dtl(iv,3) * ga1 * dxiinv * Ah(iv,3,idof,jdof)
              mat(iv,3,idof,jdof,1) = alfa * dtl(iv,3) * ga2 * dxiinv * Ah(iv,3,idof,jdof)
              mat(iv,3,idof,jdof,2) = alfa * dtl(iv,3) * Dh(iv,3,idof,jdof)
              mat(iv,3,idof,jdof,3) = alfa * dtl(iv,3) * ga3 * dxiinv * Ah(iv,3,idof,jdof)
              bc(iv,idof,jdof,6)    = alfa * dtl(iv,3) * ga4 * dxiinv * Ah(iv,3,idof,jdof)
            end do
          end do
        end do
        
!.... \hat{V}_{\eta\eta} term

        do iv = 1, ny
        
        bc(iv,2,2,7)    = bc(iv,2,2,7)    - alfa * dtl(iv,3) * da1 * dxisinv * Vh(iv,3,1)
        bc(iv,2,3,7)    = bc(iv,2,3,7)    - alfa * dtl(iv,3) * da1 * dxisinv * Vh(iv,3,5)
        bc(iv,3,2,7)    = bc(iv,3,2,7)    - alfa * dtl(iv,3) * da1 * dxisinv * Vh(iv,3,6)
        bc(iv,3,3,7)    = bc(iv,3,3,7)    - alfa * dtl(iv,3) * da1 * dxisinv * Vh(iv,3,2)
        bc(iv,4,4,7)    = bc(iv,4,4,7)    - alfa * dtl(iv,3) * da1 * dxisinv * Vh(iv,3,3)
        bc(iv,5,5,7)    = bc(iv,5,5,7)    - alfa * dtl(iv,3) * da1 * dxisinv * Vh(iv,3,4)

        mat(iv,3,2,2,1) = mat(iv,3,2,2,1) - alfa * dtl(iv,3) * da2 * dxisinv * Vh(iv,3,1)
        mat(iv,3,2,3,1) = mat(iv,3,2,3,1) - alfa * dtl(iv,3) * da2 * dxisinv * Vh(iv,3,5)
        mat(iv,3,3,2,1) = mat(iv,3,3,2,1) - alfa * dtl(iv,3) * da2 * dxisinv * Vh(iv,3,6)
        mat(iv,3,3,3,1) = mat(iv,3,3,3,1) - alfa * dtl(iv,3) * da2 * dxisinv * Vh(iv,3,2)
        mat(iv,3,4,4,1) = mat(iv,3,4,4,1) - alfa * dtl(iv,3) * da2 * dxisinv * Vh(iv,3,3)
        mat(iv,3,5,5,1) = mat(iv,3,5,5,1) - alfa * dtl(iv,3) * da2 * dxisinv * Vh(iv,3,4)

        mat(iv,3,2,2,2) = mat(iv,3,2,2,2) - alfa * dtl(iv,3) * da3 * dxisinv * Vh(iv,3,1)
        mat(iv,3,2,3,2) = mat(iv,3,2,3,2) - alfa * dtl(iv,3) * da3 * dxisinv * Vh(iv,3,5)
        mat(iv,3,3,2,2) = mat(iv,3,3,2,2) - alfa * dtl(iv,3) * da3 * dxisinv * Vh(iv,3,6)
        mat(iv,3,3,3,2) = mat(iv,3,3,3,2) - alfa * dtl(iv,3) * da3 * dxisinv * Vh(iv,3,2)
        mat(iv,3,4,4,2) = mat(iv,3,4,4,2) - alfa * dtl(iv,3) * da3 * dxisinv * Vh(iv,3,3)
        mat(iv,3,5,5,2) = mat(iv,3,5,5,2) - alfa * dtl(iv,3) * da3 * dxisinv * Vh(iv,3,4)

        mat(iv,3,2,2,3) = mat(iv,3,2,2,3) - alfa * dtl(iv,3) * da4 * dxisinv * Vh(iv,3,1)
        mat(iv,3,2,3,3) = mat(iv,3,2,3,3) - alfa * dtl(iv,3) * da4 * dxisinv * Vh(iv,3,5)
        mat(iv,3,3,2,3) = mat(iv,3,3,2,3) - alfa * dtl(iv,3) * da4 * dxisinv * Vh(iv,3,6)
        mat(iv,3,3,3,3) = mat(iv,3,3,3,3) - alfa * dtl(iv,3) * da4 * dxisinv * Vh(iv,3,2)
        mat(iv,3,4,4,3) = mat(iv,3,4,4,3) - alfa * dtl(iv,3) * da4 * dxisinv * Vh(iv,3,3)
        mat(iv,3,5,5,3) = mat(iv,3,5,5,3) - alfa * dtl(iv,3) * da4 * dxisinv * Vh(iv,3,4)

        bc(iv,2,2,6)    = bc(iv,2,2,6)    - alfa * dtl(iv,3) * da5 * dxisinv * Vh(iv,3,1)
        bc(iv,2,3,6)    = bc(iv,2,3,6)    - alfa * dtl(iv,3) * da5 * dxisinv * Vh(iv,3,5)
        bc(iv,3,2,6)    = bc(iv,3,2,6)    - alfa * dtl(iv,3) * da5 * dxisinv * Vh(iv,3,6)
        bc(iv,3,3,6)    = bc(iv,3,3,6)    - alfa * dtl(iv,3) * da5 * dxisinv * Vh(iv,3,2)
        bc(iv,4,4,6)    = bc(iv,4,4,6)    - alfa * dtl(iv,3) * da5 * dxisinv * Vh(iv,3,3)
        bc(iv,5,5,6)    = bc(iv,5,5,6)    - alfa * dtl(iv,3) * da5 * dxisinv * Vh(iv,3,4)

!.... I term

        mat(iv,3,1,1,2) = mat(iv,3,1,1,2) + one
        mat(iv,3,2,2,2) = mat(iv,3,2,2,2) + one
        mat(iv,3,3,3,2) = mat(iv,3,3,3,2) + one
        mat(iv,3,4,4,2) = mat(iv,3,4,4,2) + one
        mat(iv,3,5,5,2) = mat(iv,3,5,5,2) + one

        end do

!.... implicit damping term updated to fourth order

        if (eps_e .ne. zero) then

        do iv = 1, ny

        bc(iv,1,1,7) = bc(iv,1,1,7) + alfa * dtl(iv,3) * eps_e * buff(iv,3) * fb1
        bc(iv,2,2,7) = bc(iv,2,2,7) + alfa * dtl(iv,3) * eps_e * buff(iv,3) * fb1
        bc(iv,3,3,7) = bc(iv,3,3,7) + alfa * dtl(iv,3) * eps_e * buff(iv,3) * fb1
        bc(iv,4,4,7) = bc(iv,4,4,7) + alfa * dtl(iv,3) * eps_e * buff(iv,3) * fb1
        bc(iv,5,5,7) = bc(iv,5,5,7) + alfa * dtl(iv,3) * eps_e * buff(iv,3) * fb1

        mat(iv,3,1,1,1) = mat(iv,3,1,1,1) + alfa * dtl(iv,3) * eps_e * buff(iv,3) * fb2
        mat(iv,3,2,2,1) = mat(iv,3,2,2,1) + alfa * dtl(iv,3) * eps_e * buff(iv,3) * fb2
        mat(iv,3,3,3,1) = mat(iv,3,3,3,1) + alfa * dtl(iv,3) * eps_e * buff(iv,3) * fb2
        mat(iv,3,4,4,1) = mat(iv,3,4,4,1) + alfa * dtl(iv,3) * eps_e * buff(iv,3) * fb2
        mat(iv,3,5,5,1) = mat(iv,3,5,5,1) + alfa * dtl(iv,3) * eps_e * buff(iv,3) * fb2

        mat(iv,3,1,1,2) = mat(iv,3,1,1,2) + alfa * dtl(iv,3) * eps_e * buff(iv,3) * fb3
        mat(iv,3,2,2,2) = mat(iv,3,2,2,2) + alfa * dtl(iv,3) * eps_e * buff(iv,3) * fb3
        mat(iv,3,3,3,2) = mat(iv,3,3,3,2) + alfa * dtl(iv,3) * eps_e * buff(iv,3) * fb3
        mat(iv,3,4,4,2) = mat(iv,3,4,4,2) + alfa * dtl(iv,3) * eps_e * buff(iv,3) * fb3
        mat(iv,3,5,5,2) = mat(iv,3,5,5,2) + alfa * dtl(iv,3) * eps_e * buff(iv,3) * fb3

        mat(iv,3,1,1,3) = mat(iv,3,1,1,3) + alfa * dtl(iv,3) * eps_e * buff(iv,3) * fb4
        mat(iv,3,2,2,3) = mat(iv,3,2,2,3) + alfa * dtl(iv,3) * eps_e * buff(iv,3) * fb4
        mat(iv,3,3,3,3) = mat(iv,3,3,3,3) + alfa * dtl(iv,3) * eps_e * buff(iv,3) * fb4
        mat(iv,3,4,4,3) = mat(iv,3,4,4,3) + alfa * dtl(iv,3) * eps_e * buff(iv,3) * fb4
        mat(iv,3,5,5,3) = mat(iv,3,5,5,3) + alfa * dtl(iv,3) * eps_e * buff(iv,3) * fb4

        bc(iv,1,1,6) = bc(iv,1,1,6) + alfa * dtl(iv,3) * eps_e * buff(iv,3) * fb5
        bc(iv,2,2,6) = bc(iv,2,2,6) + alfa * dtl(iv,3) * eps_e * buff(iv,3) * fb5
        bc(iv,3,3,6) = bc(iv,3,3,6) + alfa * dtl(iv,3) * eps_e * buff(iv,3) * fb5
        bc(iv,4,4,6) = bc(iv,4,4,6) + alfa * dtl(iv,3) * eps_e * buff(iv,3) * fb5
        bc(iv,5,5,6) = bc(iv,5,5,6) + alfa * dtl(iv,3) * eps_e * buff(iv,3) * fb5

        end do

        end if

!=======================================================================================================!
!.... use higher-order tangent on second node off the right boundary
!=======================================================================================================!
        do idof = 1, ndof
          do jdof = 1, ndof
            do iv = 1, ny
              bc(iv,idof,jdof,8)       = alfa * dtl(iv,nx-2) * ga1 * dxiinv * Ah(iv,nx-2,idof,jdof)
              mat(iv,nx-2,idof,jdof,1) = alfa * dtl(iv,nx-2) * ga2 * dxiinv * Ah(iv,nx-2,idof,jdof)
              mat(iv,nx-2,idof,jdof,2) = alfa * dtl(iv,nx-2) * Dh(iv,nx-2,idof,jdof)
              mat(iv,nx-2,idof,jdof,3) = alfa * dtl(iv,nx-2) * ga3 * dxiinv * Ah(iv,nx-2,idof,jdof)
              bc(iv,idof,jdof,9)       = alfa * dtl(iv,nx-2) * ga4 * dxiinv * Ah(iv,nx-2,idof,jdof)
            end do
          end do
        end do
                
!.... \hat{V}_{\eta\eta} term

        do iv = 1, ny
        
        bc(iv,2,2,8)       = bc(iv,2,2,8)       - alfa * dtl(iv,nx-2) * da1 * dxisinv * Vh(iv,nx-2,1)
        bc(iv,2,3,8)       = bc(iv,2,3,8)       - alfa * dtl(iv,nx-2) * da1 * dxisinv * Vh(iv,nx-2,5)
        bc(iv,3,2,8)       = bc(iv,3,2,8)       - alfa * dtl(iv,nx-2) * da1 * dxisinv * Vh(iv,nx-2,6)
        bc(iv,3,3,8)       = bc(iv,3,3,8)       - alfa * dtl(iv,nx-2) * da1 * dxisinv * Vh(iv,nx-2,2)
        bc(iv,4,4,8)       = bc(iv,4,4,8)       - alfa * dtl(iv,nx-2) * da1 * dxisinv * Vh(iv,nx-2,3)
        bc(iv,5,5,8)       = bc(iv,5,5,8)       - alfa * dtl(iv,nx-2) * da1 * dxisinv * Vh(iv,nx-2,4)

        mat(iv,nx-2,2,2,1) = mat(iv,nx-2,2,2,1) - alfa * dtl(iv,nx-2) * da2 * dxisinv * Vh(iv,nx-2,1)
        mat(iv,nx-2,2,3,1) = mat(iv,nx-2,2,3,1) - alfa * dtl(iv,nx-2) * da2 * dxisinv * Vh(iv,nx-2,5)
        mat(iv,nx-2,3,2,1) = mat(iv,nx-2,3,2,1) - alfa * dtl(iv,nx-2) * da2 * dxisinv * Vh(iv,nx-2,6)
        mat(iv,nx-2,3,3,1) = mat(iv,nx-2,3,3,1) - alfa * dtl(iv,nx-2) * da2 * dxisinv * Vh(iv,nx-2,2)
        mat(iv,nx-2,4,4,1) = mat(iv,nx-2,4,4,1) - alfa * dtl(iv,nx-2) * da2 * dxisinv * Vh(iv,nx-2,3)
        mat(iv,nx-2,5,5,1) = mat(iv,nx-2,5,5,1) - alfa * dtl(iv,nx-2) * da2 * dxisinv * Vh(iv,nx-2,4)

        mat(iv,nx-2,2,2,2) = mat(iv,nx-2,2,2,2) - alfa * dtl(iv,nx-2) * da3 * dxisinv * Vh(iv,nx-2,1)
        mat(iv,nx-2,2,3,2) = mat(iv,nx-2,2,3,2) - alfa * dtl(iv,nx-2) * da3 * dxisinv * Vh(iv,nx-2,5)
        mat(iv,nx-2,3,2,2) = mat(iv,nx-2,3,2,2) - alfa * dtl(iv,nx-2) * da3 * dxisinv * Vh(iv,nx-2,6)
        mat(iv,nx-2,3,3,2) = mat(iv,nx-2,3,3,2) - alfa * dtl(iv,nx-2) * da3 * dxisinv * Vh(iv,nx-2,2)
        mat(iv,nx-2,4,4,2) = mat(iv,nx-2,4,4,2) - alfa * dtl(iv,nx-2) * da3 * dxisinv * Vh(iv,nx-2,3)
        mat(iv,nx-2,5,5,2) = mat(iv,nx-2,5,5,2) - alfa * dtl(iv,nx-2) * da3 * dxisinv * Vh(iv,nx-2,4)

        mat(iv,nx-2,2,2,3) = mat(iv,nx-2,2,2,3) - alfa * dtl(iv,nx-2) * da4 * dxisinv * Vh(iv,nx-2,1)
        mat(iv,nx-2,2,3,3) = mat(iv,nx-2,2,3,3) - alfa * dtl(iv,nx-2) * da4 * dxisinv * Vh(iv,nx-2,5)
        mat(iv,nx-2,3,2,3) = mat(iv,nx-2,3,2,3) - alfa * dtl(iv,nx-2) * da4 * dxisinv * Vh(iv,nx-2,6)
        mat(iv,nx-2,3,3,3) = mat(iv,nx-2,3,3,3) - alfa * dtl(iv,nx-2) * da4 * dxisinv * Vh(iv,nx-2,2)
        mat(iv,nx-2,4,4,3) = mat(iv,nx-2,4,4,3) - alfa * dtl(iv,nx-2) * da4 * dxisinv * Vh(iv,nx-2,3)
        mat(iv,nx-2,5,5,3) = mat(iv,nx-2,5,5,3) - alfa * dtl(iv,nx-2) * da4 * dxisinv * Vh(iv,nx-2,4)

        bc(iv,2,2,9)       = bc(iv,2,2,9)       - alfa * dtl(iv,nx-2) * da5 * dxisinv * Vh(iv,nx-2,1)
        bc(iv,2,3,9)       = bc(iv,2,3,9)       - alfa * dtl(iv,nx-2) * da5 * dxisinv * Vh(iv,nx-2,5)
        bc(iv,3,2,9)       = bc(iv,3,2,9)       - alfa * dtl(iv,nx-2) * da5 * dxisinv * Vh(iv,nx-2,6)
        bc(iv,3,3,9)       = bc(iv,3,3,9)       - alfa * dtl(iv,nx-2) * da5 * dxisinv * Vh(iv,nx-2,2)
        bc(iv,4,4,9)       = bc(iv,4,4,9)       - alfa * dtl(iv,nx-2) * da5 * dxisinv * Vh(iv,nx-2,3)
        bc(iv,5,5,9)       = bc(iv,5,5,9)       - alfa * dtl(iv,nx-2) * da5 * dxisinv * Vh(iv,nx-2,4)

!.... I term

        mat(iv,nx-2,1,1,2) = mat(iv,nx-2,1,1,2) + one
        mat(iv,nx-2,2,2,2) = mat(iv,nx-2,2,2,2) + one
        mat(iv,nx-2,3,3,2) = mat(iv,nx-2,3,3,2) + one
        mat(iv,nx-2,4,4,2) = mat(iv,nx-2,4,4,2) + one
        mat(iv,nx-2,5,5,2) = mat(iv,nx-2,5,5,2) + one

        end do

!.... implicit damping term updated to fourth order

        if (eps_e .ne. zero) then

        do iv = 1, ny

        bc(iv,1,1,8) = bc(iv,1,1,8) + alfa * dtl(iv,nx-2) * eps_e * buff(iv,nx-2) * fb1
        bc(iv,2,2,8) = bc(iv,2,2,8) + alfa * dtl(iv,nx-2) * eps_e * buff(iv,nx-2) * fb1
        bc(iv,3,3,8) = bc(iv,3,3,8) + alfa * dtl(iv,nx-2) * eps_e * buff(iv,nx-2) * fb1
        bc(iv,4,4,8) = bc(iv,4,4,8) + alfa * dtl(iv,nx-2) * eps_e * buff(iv,nx-2) * fb1
        bc(iv,5,5,8) = bc(iv,5,5,8) + alfa * dtl(iv,nx-2) * eps_e * buff(iv,nx-2) * fb1

        mat(iv,nx-2,1,1,1) = mat(iv,nx-2,1,1,1) + alfa * dtl(iv,nx-2) * eps_e * buff(iv,nx-2) * fb2
        mat(iv,nx-2,2,2,1) = mat(iv,nx-2,2,2,1) + alfa * dtl(iv,nx-2) * eps_e * buff(iv,nx-2) * fb2
        mat(iv,nx-2,3,3,1) = mat(iv,nx-2,3,3,1) + alfa * dtl(iv,nx-2) * eps_e * buff(iv,nx-2) * fb2
        mat(iv,nx-2,4,4,1) = mat(iv,nx-2,4,4,1) + alfa * dtl(iv,nx-2) * eps_e * buff(iv,nx-2) * fb2
        mat(iv,nx-2,5,5,1) = mat(iv,nx-2,5,5,1) + alfa * dtl(iv,nx-2) * eps_e * buff(iv,nx-2) * fb2

        mat(iv,nx-2,1,1,2) = mat(iv,nx-2,1,1,2) + alfa * dtl(iv,nx-2) * eps_e * buff(iv,nx-2) * fb3
        mat(iv,nx-2,2,2,2) = mat(iv,nx-2,2,2,2) + alfa * dtl(iv,nx-2) * eps_e * buff(iv,nx-2) * fb3
        mat(iv,nx-2,3,3,2) = mat(iv,nx-2,3,3,2) + alfa * dtl(iv,nx-2) * eps_e * buff(iv,nx-2) * fb3
        mat(iv,nx-2,4,4,2) = mat(iv,nx-2,4,4,2) + alfa * dtl(iv,nx-2) * eps_e * buff(iv,nx-2) * fb3
        mat(iv,nx-2,5,5,2) = mat(iv,nx-2,5,5,2) + alfa * dtl(iv,nx-2) * eps_e * buff(iv,nx-2) * fb3

        mat(iv,nx-2,1,1,3) = mat(iv,nx-2,1,1,3) + alfa * dtl(iv,nx-2) * eps_e * buff(iv,nx-2) * fb4
        mat(iv,nx-2,2,2,3) = mat(iv,nx-2,2,2,3) + alfa * dtl(iv,nx-2) * eps_e * buff(iv,nx-2) * fb4
        mat(iv,nx-2,3,3,3) = mat(iv,nx-2,3,3,3) + alfa * dtl(iv,nx-2) * eps_e * buff(iv,nx-2) * fb4
        mat(iv,nx-2,4,4,3) = mat(iv,nx-2,4,4,3) + alfa * dtl(iv,nx-2) * eps_e * buff(iv,nx-2) * fb4
        mat(iv,nx-2,5,5,3) = mat(iv,nx-2,5,5,3) + alfa * dtl(iv,nx-2) * eps_e * buff(iv,nx-2) * fb4

        bc(iv,1,1,9) = bc(iv,1,1,9) + alfa * dtl(iv,nx-2) * eps_e * buff(iv,nx-2) * fb5
        bc(iv,2,2,9) = bc(iv,2,2,9) + alfa * dtl(iv,nx-2) * eps_e * buff(iv,nx-2) * fb5
        bc(iv,3,3,9) = bc(iv,3,3,9) + alfa * dtl(iv,nx-2) * eps_e * buff(iv,nx-2) * fb5
        bc(iv,4,4,9) = bc(iv,4,4,9) + alfa * dtl(iv,nx-2) * eps_e * buff(iv,nx-2) * fb5
        bc(iv,5,5,9) = bc(iv,5,5,9) + alfa * dtl(iv,nx-2) * eps_e * buff(iv,nx-2) * fb5

        end do

        end if

        if (rsym) then
        
          call lhsbs1( mat, Ah, Dh, Vh, bc, dtl, 2 )             ! apply a symmetry condition

        else
                
!=======================================================================================================!
!.... use higher-order tangent on first node off the right boundary
!=======================================================================================================!
        do idof = 1, ndof
          do jdof = 1, ndof
            do iv = 1, ny
              bc(iv,idof,jdof,10)      = -alfa * dtl(iv,nx-1) * gb5 * dxiinv * Ah(iv,nx-1,idof,jdof)
              bc(iv,idof,jdof,11)      = -alfa * dtl(iv,nx-1) * gb4 * dxiinv * Ah(iv,nx-1,idof,jdof)
              mat(iv,nx-1,idof,jdof,1) = -alfa * dtl(iv,nx-1) * gb3 * dxiinv * Ah(iv,nx-1,idof,jdof)
              mat(iv,nx-1,idof,jdof,2) = -alfa * dtl(iv,nx-1) * gb2 * dxiinv * Ah(iv,nx-1,idof,jdof) + &
                                          alfa * dtl(iv,nx-1) * Dh(iv,nx-1,idof,jdof)
              mat(iv,nx-1,idof,jdof,3) = -alfa * dtl(iv,nx-1) * gb1 * dxiinv * Ah(iv,nx-1,idof,jdof)
            end do
          end do
        end do
        
!.... \hat{V}_{\eta\eta} term

        do iv = 1, ny
        
        bc(iv,2,2,10)      = bc(iv,2,2,10)      - alfa * dtl(iv,nx-1) * db5 * dxisinv * Vh(iv,nx-1,1)
        bc(iv,2,3,10)      = bc(iv,2,3,10)      - alfa * dtl(iv,nx-1) * db5 * dxisinv * Vh(iv,nx-1,5)
        bc(iv,3,2,10)      = bc(iv,3,2,10)      - alfa * dtl(iv,nx-1) * db5 * dxisinv * Vh(iv,nx-1,6)
        bc(iv,3,3,10)      = bc(iv,3,3,10)      - alfa * dtl(iv,nx-1) * db5 * dxisinv * Vh(iv,nx-1,2)
        bc(iv,4,4,10)      = bc(iv,4,4,10)      - alfa * dtl(iv,nx-1) * db5 * dxisinv * Vh(iv,nx-1,3)
        bc(iv,5,5,10)      = bc(iv,5,5,10)      - alfa * dtl(iv,nx-1) * db5 * dxisinv * Vh(iv,nx-1,4)

        bc(iv,2,2,11)      = bc(iv,2,2,11)      - alfa * dtl(iv,nx-1) * db4 * dxisinv * Vh(iv,nx-1,1)
        bc(iv,2,3,11)      = bc(iv,2,3,11)      - alfa * dtl(iv,nx-1) * db4 * dxisinv * Vh(iv,nx-1,5)
        bc(iv,3,2,11)      = bc(iv,3,2,11)      - alfa * dtl(iv,nx-1) * db4 * dxisinv * Vh(iv,nx-1,6)
        bc(iv,3,3,11)      = bc(iv,3,3,11)      - alfa * dtl(iv,nx-1) * db4 * dxisinv * Vh(iv,nx-1,2)
        bc(iv,4,4,11)      = bc(iv,4,4,11)      - alfa * dtl(iv,nx-1) * db4 * dxisinv * Vh(iv,nx-1,3)
        bc(iv,5,5,11)      = bc(iv,5,5,11)      - alfa * dtl(iv,nx-1) * db4 * dxisinv * Vh(iv,nx-1,4)

        mat(iv,nx-1,2,2,1) = mat(iv,nx-1,2,2,1) - alfa * dtl(iv,nx-1) * db3 * dxisinv * Vh(iv,nx-1,1)
        mat(iv,nx-1,2,3,1) = mat(iv,nx-1,2,3,1) - alfa * dtl(iv,nx-1) * db3 * dxisinv * Vh(iv,nx-1,5)
        mat(iv,nx-1,3,2,1) = mat(iv,nx-1,3,2,1) - alfa * dtl(iv,nx-1) * db3 * dxisinv * Vh(iv,nx-1,6)
        mat(iv,nx-1,3,3,1) = mat(iv,nx-1,3,3,1) - alfa * dtl(iv,nx-1) * db3 * dxisinv * Vh(iv,nx-1,2)
        mat(iv,nx-1,4,4,1) = mat(iv,nx-1,4,4,1) - alfa * dtl(iv,nx-1) * db3 * dxisinv * Vh(iv,nx-1,3)
        mat(iv,nx-1,5,5,1) = mat(iv,nx-1,5,5,1) - alfa * dtl(iv,nx-1) * db3 * dxisinv * Vh(iv,nx-1,4)

        mat(iv,nx-1,2,2,2) = mat(iv,nx-1,2,2,2) - alfa * dtl(iv,nx-1) * db2 * dxisinv * Vh(iv,nx-1,1)
        mat(iv,nx-1,2,3,2) = mat(iv,nx-1,2,3,2) - alfa * dtl(iv,nx-1) * db2 * dxisinv * Vh(iv,nx-1,5)
        mat(iv,nx-1,3,2,2) = mat(iv,nx-1,3,2,2) - alfa * dtl(iv,nx-1) * db2 * dxisinv * Vh(iv,nx-1,6)
        mat(iv,nx-1,3,3,2) = mat(iv,nx-1,3,3,2) - alfa * dtl(iv,nx-1) * db2 * dxisinv * Vh(iv,nx-1,2)
        mat(iv,nx-1,4,4,2) = mat(iv,nx-1,4,4,2) - alfa * dtl(iv,nx-1) * db2 * dxisinv * Vh(iv,nx-1,3)
        mat(iv,nx-1,5,5,2) = mat(iv,nx-1,5,5,2) - alfa * dtl(iv,nx-1) * db2 * dxisinv * Vh(iv,nx-1,4)

        mat(iv,nx-1,2,2,3) = mat(iv,nx-1,2,2,3) - alfa * dtl(iv,nx-1) * db1 * dxisinv * Vh(iv,nx-1,1)
        mat(iv,nx-1,2,3,3) = mat(iv,nx-1,2,3,3) - alfa * dtl(iv,nx-1) * db1 * dxisinv * Vh(iv,nx-1,5)
        mat(iv,nx-1,3,2,3) = mat(iv,nx-1,3,2,3) - alfa * dtl(iv,nx-1) * db1 * dxisinv * Vh(iv,nx-1,6)
        mat(iv,nx-1,3,3,3) = mat(iv,nx-1,3,3,3) - alfa * dtl(iv,nx-1) * db1 * dxisinv * Vh(iv,nx-1,2)
        mat(iv,nx-1,4,4,3) = mat(iv,nx-1,4,4,3) - alfa * dtl(iv,nx-1) * db1 * dxisinv * Vh(iv,nx-1,3)
        mat(iv,nx-1,5,5,3) = mat(iv,nx-1,5,5,3) - alfa * dtl(iv,nx-1) * db1 * dxisinv * Vh(iv,nx-1,4)

!.... I term

        mat(iv,nx-1,1,1,2) = mat(iv,nx-1,1,1,2) + one
        mat(iv,nx-1,2,2,2) = mat(iv,nx-1,2,2,2) + one
        mat(iv,nx-1,3,3,2) = mat(iv,nx-1,3,3,2) + one
        mat(iv,nx-1,4,4,2) = mat(iv,nx-1,4,4,2) + one
        mat(iv,nx-1,5,5,2) = mat(iv,nx-1,5,5,2) + one

        end do

!.... implicit damping term

        if (eps_e .ne. zero) then

        do iv = 1, ny
        mat(iv,nx-1,1,1,1) = mat(iv,nx-1,1,1,1) + alfa * dtl(iv,nx-1) * eps_e * buff(iv,nx-1) * (-one)
        mat(iv,nx-1,2,2,1) = mat(iv,nx-1,2,2,1) + alfa * dtl(iv,nx-1) * eps_e * buff(iv,nx-1) * (-one)
        mat(iv,nx-1,3,3,1) = mat(iv,nx-1,3,3,1) + alfa * dtl(iv,nx-1) * eps_e * buff(iv,nx-1) * (-one)
        mat(iv,nx-1,4,4,1) = mat(iv,nx-1,4,4,1) + alfa * dtl(iv,nx-1) * eps_e * buff(iv,nx-1) * (-one)
        mat(iv,nx-1,5,5,1) = mat(iv,nx-1,5,5,1) + alfa * dtl(iv,nx-1) * eps_e * buff(iv,nx-1) * (-one)

        mat(iv,nx-1,1,1,2) = mat(iv,nx-1,1,1,2) + alfa * dtl(iv,nx-1) * eps_e * buff(iv,nx-1) * two
        mat(iv,nx-1,2,2,2) = mat(iv,nx-1,2,2,2) + alfa * dtl(iv,nx-1) * eps_e * buff(iv,nx-1) * two
        mat(iv,nx-1,3,3,2) = mat(iv,nx-1,3,3,2) + alfa * dtl(iv,nx-1) * eps_e * buff(iv,nx-1) * two
        mat(iv,nx-1,4,4,2) = mat(iv,nx-1,4,4,2) + alfa * dtl(iv,nx-1) * eps_e * buff(iv,nx-1) * two
        mat(iv,nx-1,5,5,2) = mat(iv,nx-1,5,5,2) + alfa * dtl(iv,nx-1) * eps_e * buff(iv,nx-1) * two

        mat(iv,nx-1,1,1,3) = mat(iv,nx-1,1,1,3) + alfa * dtl(iv,nx-1) * eps_e * buff(iv,nx-1) * (-one)
        mat(iv,nx-1,2,2,3) = mat(iv,nx-1,2,2,3) + alfa * dtl(iv,nx-1) * eps_e * buff(iv,nx-1) * (-one)
        mat(iv,nx-1,3,3,3) = mat(iv,nx-1,3,3,3) + alfa * dtl(iv,nx-1) * eps_e * buff(iv,nx-1) * (-one)
        mat(iv,nx-1,4,4,3) = mat(iv,nx-1,4,4,3) + alfa * dtl(iv,nx-1) * eps_e * buff(iv,nx-1) * (-one)
        mat(iv,nx-1,5,5,3) = mat(iv,nx-1,5,5,3) + alfa * dtl(iv,nx-1) * eps_e * buff(iv,nx-1) * (-one)
        end do
        
        end if

!=======================================================================================================!
!.... use higher-order tangent on the right boundary
!=======================================================================================================!
        do idof = 1, ndof
          do jdof = 1, ndof
            do iv = 1, ny
              bc(iv,idof,jdof,12)    = -alfa * dtl(iv,nx) * gc5 * dxiinv * Ah(iv,nx,idof,jdof)
              bc(iv,idof,jdof,13)    = -alfa * dtl(iv,nx) * gc4 * dxiinv * Ah(iv,nx,idof,jdof)
              bc(iv,idof,jdof,14)    = -alfa * dtl(iv,nx) * gc3 * dxiinv * Ah(iv,nx,idof,jdof)
              mat(iv,nx,idof,jdof,1) = -alfa * dtl(iv,nx) * gc2 * dxiinv * Ah(iv,nx,idof,jdof)
              mat(iv,nx,idof,jdof,2) = -alfa * dtl(iv,nx) * gc1 * dxiinv * Ah(iv,nx,idof,jdof) + &
                                        alfa * dtl(iv,nx) * Dh(iv,nx,idof,jdof)
              mat(iv,nx,idof,jdof,3) = zero
            end do
          end do
        end do
        
!.... \hat{V}_{\eta\eta} term

        do iv = 1, ny

        bc(iv,2,2,12) = bc(iv,2,2,12) - alfa * dtl(iv,nx) * dd5 * dxisinv * Vh(iv,nx,1)
        bc(iv,2,3,12) = bc(iv,2,3,12) - alfa * dtl(iv,nx) * dd5 * dxisinv * Vh(iv,nx,5)
        bc(iv,3,2,12) = bc(iv,3,2,12) - alfa * dtl(iv,nx) * dd5 * dxisinv * Vh(iv,nx,6)
        bc(iv,3,3,12) = bc(iv,3,3,12) - alfa * dtl(iv,nx) * dd5 * dxisinv * Vh(iv,nx,2)
        bc(iv,4,4,12) = bc(iv,4,4,12) - alfa * dtl(iv,nx) * dd5 * dxisinv * Vh(iv,nx,3)
        bc(iv,5,5,12) = bc(iv,5,5,12) - alfa * dtl(iv,nx) * dd5 * dxisinv * Vh(iv,nx,4)

        bc(iv,2,2,13) = bc(iv,2,2,13) - alfa * dtl(iv,nx) * dd4 * dxisinv * Vh(iv,nx,1)
        bc(iv,2,3,13) = bc(iv,2,3,13) - alfa * dtl(iv,nx) * dd4 * dxisinv * Vh(iv,nx,5)
        bc(iv,3,2,13) = bc(iv,3,2,13) - alfa * dtl(iv,nx) * dd4 * dxisinv * Vh(iv,nx,6)
        bc(iv,3,3,13) = bc(iv,3,3,13) - alfa * dtl(iv,nx) * dd4 * dxisinv * Vh(iv,nx,2)
        bc(iv,4,4,13) = bc(iv,4,4,13) - alfa * dtl(iv,nx) * dd4 * dxisinv * Vh(iv,nx,3)
        bc(iv,5,5,13) = bc(iv,5,5,13) - alfa * dtl(iv,nx) * dd4 * dxisinv * Vh(iv,nx,4)

        bc(iv,2,2,14) = bc(iv,2,2,14) - alfa * dtl(iv,nx) * dd3 * dxisinv * Vh(iv,nx,1)
        bc(iv,2,3,14) = bc(iv,2,3,14) - alfa * dtl(iv,nx) * dd3 * dxisinv * Vh(iv,nx,5)
        bc(iv,3,2,14) = bc(iv,3,2,14) - alfa * dtl(iv,nx) * dd3 * dxisinv * Vh(iv,nx,6)
        bc(iv,3,3,14) = bc(iv,3,3,14) - alfa * dtl(iv,nx) * dd3 * dxisinv * Vh(iv,nx,2)
        bc(iv,4,4,14) = bc(iv,4,4,14) - alfa * dtl(iv,nx) * dd3 * dxisinv * Vh(iv,nx,3)
        bc(iv,5,5,14) = bc(iv,5,5,14) - alfa * dtl(iv,nx) * dd3 * dxisinv * Vh(iv,nx,4)

        mat(iv,nx,2,2,1) = mat(iv,nx,2,2,1) - alfa * dtl(iv,nx) * dd2 * dxisinv * Vh(iv,nx,1)
        mat(iv,nx,2,3,1) = mat(iv,nx,2,3,1) - alfa * dtl(iv,nx) * dd2 * dxisinv * Vh(iv,nx,5)
        mat(iv,nx,3,2,1) = mat(iv,nx,3,2,1) - alfa * dtl(iv,nx) * dd2 * dxisinv * Vh(iv,nx,6)
        mat(iv,nx,3,3,1) = mat(iv,nx,3,3,1) - alfa * dtl(iv,nx) * dd2 * dxisinv * Vh(iv,nx,2)
        mat(iv,nx,4,4,1) = mat(iv,nx,4,4,1) - alfa * dtl(iv,nx) * dd2 * dxisinv * Vh(iv,nx,3)
        mat(iv,nx,5,5,1) = mat(iv,nx,5,5,1) - alfa * dtl(iv,nx) * dd2 * dxisinv * Vh(iv,nx,4)

        mat(iv,nx,2,2,2) = mat(iv,nx,2,2,2) - alfa * dtl(iv,nx) * dd1 * dxisinv * Vh(iv,nx,1)
        mat(iv,nx,2,3,2) = mat(iv,nx,2,3,2) - alfa * dtl(iv,nx) * dd1 * dxisinv * Vh(iv,nx,5)
        mat(iv,nx,3,2,2) = mat(iv,nx,3,2,2) - alfa * dtl(iv,nx) * dd1 * dxisinv * Vh(iv,nx,6)
        mat(iv,nx,3,3,2) = mat(iv,nx,3,3,2) - alfa * dtl(iv,nx) * dd1 * dxisinv * Vh(iv,nx,2)
        mat(iv,nx,4,4,2) = mat(iv,nx,4,4,2) - alfa * dtl(iv,nx) * dd1 * dxisinv * Vh(iv,nx,3)
        mat(iv,nx,5,5,2) = mat(iv,nx,5,5,2) - alfa * dtl(iv,nx) * dd1 * dxisinv * Vh(iv,nx,4)

!.... I term

        mat(iv,nx,1,1,2) = mat(iv,nx,1,1,2) + one
        mat(iv,nx,2,2,2) = mat(iv,nx,2,2,2) + one
        mat(iv,nx,3,3,2) = mat(iv,nx,3,3,2) + one
        mat(iv,nx,4,4,2) = mat(iv,nx,4,4,2) + one
        mat(iv,nx,5,5,2) = mat(iv,nx,5,5,2) + one

        end do

        end if                  ! rsym

        return
        end
