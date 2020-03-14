!=======================================================================================================!
        subroutine lhsbt2( mat, Bh, Dh, Vh, bc, spgl, spg2l, dtl )
!  
!  Correct the LHS for boundary treatment in the \eta direction
!
!  This version supports the 4th order LHS
!
!  Revised: 6-20-96
!
!=======================================================================================================!
        use global
        use stencil
        use buff_mod
        implicit none
        
        real :: mat(ny,nx,ndof,ndof,3), Bh(ny,nx,ndof,ndof)
        real :: Dh(ny,nx,ndof,ndof), bc(nx,ndof,ndof,14)
        real :: spgl(ny,nx), spg2l(ny,nx), Vh(ny,nx,6)
        real :: dtl(ny,nx)

        real :: detainv, detasinv
        
        integer :: iv, idof, jdof
!=======================================================================================================!
        detainv  = one / deta
        detasinv = one / deta**2
        
!.... initialize the bc array

        bc = zero

        if (.not. yper) then
        
!=======================================================================================================!
!.... use higher-order tangent on body node
!=======================================================================================================!
        do idof = 1, ndof
          do jdof = 1, ndof
            do iv = 1, nx
              mat(1,iv,idof,jdof,1) = zero
              mat(1,iv,idof,jdof,2) = alfa * dtl(1,iv) * gc1 * detainv * Bh(1,iv,idof,jdof)
              mat(1,iv,idof,jdof,3) = alfa * dtl(1,iv) * gc2 * detainv * Bh(1,iv,idof,jdof)
              bc(iv,idof,jdof,1)    = alfa * dtl(1,iv) * gc3 * detainv * Bh(1,iv,idof,jdof)
              bc(iv,idof,jdof,2)    = alfa * dtl(1,iv) * gc4 * detainv * Bh(1,iv,idof,jdof)
              bc(iv,idof,jdof,3)    = alfa * dtl(1,iv) * gc5 * detainv * Bh(1,iv,idof,jdof)
            end do
          end do
        end do
        
!.... \hat{V}_{\eta\eta} term

        do iv = 1, nx
        
        mat(1,iv,2,2,2) = mat(1,iv,2,2,2) - alfa * dtl(1,iv) * dd1 * detasinv * Vh(1,iv,1)
        mat(1,iv,2,3,2) = mat(1,iv,2,3,2) - alfa * dtl(1,iv) * dd1 * detasinv * Vh(1,iv,5)
        mat(1,iv,3,2,2) = mat(1,iv,3,2,2) - alfa * dtl(1,iv) * dd1 * detasinv * Vh(1,iv,6)
        mat(1,iv,3,3,2) = mat(1,iv,3,3,2) - alfa * dtl(1,iv) * dd1 * detasinv * Vh(1,iv,2)
        mat(1,iv,4,4,2) = mat(1,iv,4,4,2) - alfa * dtl(1,iv) * dd1 * detasinv * Vh(1,iv,3)
        mat(1,iv,5,5,2) = mat(1,iv,5,5,2) - alfa * dtl(1,iv) * dd1 * detasinv * Vh(1,iv,4)

        mat(1,iv,2,2,3) = mat(1,iv,2,2,3) - alfa * dtl(1,iv) * dd2 * detasinv * Vh(1,iv,1)
        mat(1,iv,2,3,3) = mat(1,iv,2,3,3) - alfa * dtl(1,iv) * dd2 * detasinv * Vh(1,iv,5)
        mat(1,iv,3,2,3) = mat(1,iv,3,2,3) - alfa * dtl(1,iv) * dd2 * detasinv * Vh(1,iv,6)
        mat(1,iv,3,3,3) = mat(1,iv,3,3,3) - alfa * dtl(1,iv) * dd2 * detasinv * Vh(1,iv,2)
        mat(1,iv,4,4,3) = mat(1,iv,4,4,3) - alfa * dtl(1,iv) * dd2 * detasinv * Vh(1,iv,3)
        mat(1,iv,5,5,3) = mat(1,iv,5,5,3) - alfa * dtl(1,iv) * dd2 * detasinv * Vh(1,iv,4)

        bc(iv,2,2,1) = bc(iv,2,2,1) - alfa * dtl(1,iv) * dd3 * detasinv * Vh(1,iv,1)
        bc(iv,2,3,1) = bc(iv,2,3,1) - alfa * dtl(1,iv) * dd3 * detasinv * Vh(1,iv,5)
        bc(iv,3,2,1) = bc(iv,3,2,1) - alfa * dtl(1,iv) * dd3 * detasinv * Vh(1,iv,6)
        bc(iv,3,3,1) = bc(iv,3,3,1) - alfa * dtl(1,iv) * dd3 * detasinv * Vh(1,iv,2)
        bc(iv,4,4,1) = bc(iv,4,4,1) - alfa * dtl(1,iv) * dd3 * detasinv * Vh(1,iv,3)
        bc(iv,5,5,1) = bc(iv,5,5,1) - alfa * dtl(1,iv) * dd3 * detasinv * Vh(1,iv,4)

        bc(iv,2,2,2) = bc(iv,2,2,2) - alfa * dtl(1,iv) * dd4 * detasinv * Vh(1,iv,1)
        bc(iv,2,3,2) = bc(iv,2,3,2) - alfa * dtl(1,iv) * dd4 * detasinv * Vh(1,iv,5)
        bc(iv,3,2,2) = bc(iv,3,2,2) - alfa * dtl(1,iv) * dd4 * detasinv * Vh(1,iv,6)
        bc(iv,3,3,2) = bc(iv,3,3,2) - alfa * dtl(1,iv) * dd4 * detasinv * Vh(1,iv,2)
        bc(iv,4,4,2) = bc(iv,4,4,2) - alfa * dtl(1,iv) * dd4 * detasinv * Vh(1,iv,3)
        bc(iv,5,5,2) = bc(iv,5,5,2) - alfa * dtl(1,iv) * dd4 * detasinv * Vh(1,iv,4)

        bc(iv,2,2,3) = bc(iv,2,2,3) - alfa * dtl(1,iv) * dd5 * detasinv * Vh(1,iv,1)
        bc(iv,2,3,3) = bc(iv,2,3,3) - alfa * dtl(1,iv) * dd5 * detasinv * Vh(1,iv,5)
        bc(iv,3,2,3) = bc(iv,3,2,3) - alfa * dtl(1,iv) * dd5 * detasinv * Vh(1,iv,6)
        bc(iv,3,3,3) = bc(iv,3,3,3) - alfa * dtl(1,iv) * dd5 * detasinv * Vh(1,iv,2)
        bc(iv,4,4,3) = bc(iv,4,4,3) - alfa * dtl(1,iv) * dd5 * detasinv * Vh(1,iv,3)
        bc(iv,5,5,3) = bc(iv,5,5,3) - alfa * dtl(1,iv) * dd5 * detasinv * Vh(1,iv,4)

!.... I term

        mat(1,iv,1,1,2) = mat(1,iv,1,1,2) + one
        mat(1,iv,2,2,2) = mat(1,iv,2,2,2) + one
        mat(1,iv,3,3,2) = mat(1,iv,3,3,2) + one
        mat(1,iv,4,4,2) = mat(1,iv,4,4,2) + one
        mat(1,iv,5,5,2) = mat(1,iv,5,5,2) + one

        end do
        
!.... sponge

        if (ispg .eq. 1) then
        
          do iv = 1, nx
            mat(1,iv,1,1,2) = mat(1,iv,1,1,2) + alfa * dtl(1,iv) * spgl(1,iv)
            mat(1,iv,2,2,2) = mat(1,iv,2,2,2) + alfa * dtl(1,iv) * spgl(1,iv)
            mat(1,iv,3,3,2) = mat(1,iv,3,3,2) + alfa * dtl(1,iv) * spgl(1,iv)
            mat(1,iv,4,4,2) = mat(1,iv,4,4,2) + alfa * dtl(1,iv) * spgl(1,iv)
            mat(1,iv,5,5,2) = mat(1,iv,5,5,2) + alfa * dtl(1,iv) * spgl(1,iv)
          end do

        else if (ispg .eq. 2) then
        
          do iv = 1, nx
            mat(1,iv,1,1,2) = mat(1,iv,1,1,2) + alfa * dtl(1,iv) * (spgl(1,iv) + spg2l(1,iv))
            mat(1,iv,2,2,2) = mat(1,iv,2,2,2) + alfa * dtl(1,iv) * (spgl(1,iv) + spg2l(1,iv))
            mat(1,iv,3,3,2) = mat(1,iv,3,3,2) + alfa * dtl(1,iv) * (spgl(1,iv) + spg2l(1,iv))
            mat(1,iv,4,4,2) = mat(1,iv,4,4,2) + alfa * dtl(1,iv) * (spgl(1,iv) + spg2l(1,iv))
            mat(1,iv,5,5,2) = mat(1,iv,5,5,2) + alfa * dtl(1,iv) * (spgl(1,iv) + spg2l(1,iv))
          end do
        
        end if

!=======================================================================================================!
!.... use higher-order tangent on first node off the body
!=======================================================================================================!
        do idof = 1, ndof
          do jdof = 1, ndof
            do iv = 1, nx
              mat(2,iv,idof,jdof,1) = alfa * dtl(2,iv) * gb1 * detainv * Bh(2,iv,idof,jdof)
              mat(2,iv,idof,jdof,2) = alfa * dtl(2,iv) * gb2 * detainv * Bh(2,iv,idof,jdof)
              mat(2,iv,idof,jdof,3) = alfa * dtl(2,iv) * gb3 * detainv * Bh(2,iv,idof,jdof)
              bc(iv,idof,jdof,4)    = alfa * dtl(2,iv) * gb4 * detainv * Bh(2,iv,idof,jdof)
              bc(iv,idof,jdof,5)    = alfa * dtl(2,iv) * gb5 * detainv * Bh(2,iv,idof,jdof)
            end do
          end do
        end do

!.... \hat{V}_{\eta\eta} term

        do iv = 1, nx

        mat(2,iv,2,2,1) = mat(2,iv,2,2,1) - alfa * dtl(2,iv) * db1 * detasinv * Vh(2,iv,1)
        mat(2,iv,2,3,1) = mat(2,iv,2,3,1) - alfa * dtl(2,iv) * db1 * detasinv * Vh(2,iv,5)
        mat(2,iv,3,2,1) = mat(2,iv,3,2,1) - alfa * dtl(2,iv) * db1 * detasinv * Vh(2,iv,6)
        mat(2,iv,3,3,1) = mat(2,iv,3,3,1) - alfa * dtl(2,iv) * db1 * detasinv * Vh(2,iv,2)
        mat(2,iv,4,4,1) = mat(2,iv,4,4,1) - alfa * dtl(2,iv) * db1 * detasinv * Vh(2,iv,3)
        mat(2,iv,5,5,1) = mat(2,iv,5,5,1) - alfa * dtl(2,iv) * db1 * detasinv * Vh(2,iv,4)

        mat(2,iv,2,2,2) = mat(2,iv,2,2,2) - alfa * dtl(2,iv) * db2 * detasinv * Vh(2,iv,1)
        mat(2,iv,2,3,2) = mat(2,iv,2,3,2) - alfa * dtl(2,iv) * db2 * detasinv * Vh(2,iv,5)
        mat(2,iv,3,2,2) = mat(2,iv,3,2,2) - alfa * dtl(2,iv) * db2 * detasinv * Vh(2,iv,6)
        mat(2,iv,3,3,2) = mat(2,iv,3,3,2) - alfa * dtl(2,iv) * db2 * detasinv * Vh(2,iv,2)
        mat(2,iv,4,4,2) = mat(2,iv,4,4,2) - alfa * dtl(2,iv) * db2 * detasinv * Vh(2,iv,3)
        mat(2,iv,5,5,2) = mat(2,iv,5,5,2) - alfa * dtl(2,iv) * db2 * detasinv * Vh(2,iv,4)

        mat(2,iv,2,2,3) = mat(2,iv,2,2,3) - alfa * dtl(2,iv) * db3 * detasinv * Vh(2,iv,1)
        mat(2,iv,2,3,3) = mat(2,iv,2,3,3) - alfa * dtl(2,iv) * db3 * detasinv * Vh(2,iv,5)
        mat(2,iv,3,2,3) = mat(2,iv,3,2,3) - alfa * dtl(2,iv) * db3 * detasinv * Vh(2,iv,6)
        mat(2,iv,3,3,3) = mat(2,iv,3,3,3) - alfa * dtl(2,iv) * db3 * detasinv * Vh(2,iv,2)
        mat(2,iv,4,4,3) = mat(2,iv,4,4,3) - alfa * dtl(2,iv) * db3 * detasinv * Vh(2,iv,3)
        mat(2,iv,5,5,3) = mat(2,iv,5,5,3) - alfa * dtl(2,iv) * db3 * detasinv * Vh(2,iv,4)

        bc(iv,2,2,4) = bc(iv,2,2,4) - alfa * dtl(2,iv) * db4 * detasinv * Vh(2,iv,1)
        bc(iv,2,3,4) = bc(iv,2,3,4) - alfa * dtl(2,iv) * db4 * detasinv * Vh(2,iv,5)
        bc(iv,3,2,4) = bc(iv,3,2,4) - alfa * dtl(2,iv) * db4 * detasinv * Vh(2,iv,6)
        bc(iv,3,3,4) = bc(iv,3,3,4) - alfa * dtl(2,iv) * db4 * detasinv * Vh(2,iv,2)
        bc(iv,4,4,4) = bc(iv,4,4,4) - alfa * dtl(2,iv) * db4 * detasinv * Vh(2,iv,3)
        bc(iv,5,5,4) = bc(iv,5,5,4) - alfa * dtl(2,iv) * db4 * detasinv * Vh(2,iv,4)

        bc(iv,2,2,5) = bc(iv,2,2,5) - alfa * dtl(2,iv) * db5 * detasinv * Vh(2,iv,1)
        bc(iv,2,3,5) = bc(iv,2,3,5) - alfa * dtl(2,iv) * db5 * detasinv * Vh(2,iv,5)
        bc(iv,3,2,5) = bc(iv,3,2,5) - alfa * dtl(2,iv) * db5 * detasinv * Vh(2,iv,6)
        bc(iv,3,3,5) = bc(iv,3,3,5) - alfa * dtl(2,iv) * db5 * detasinv * Vh(2,iv,2)
        bc(iv,4,4,5) = bc(iv,4,4,5) - alfa * dtl(2,iv) * db5 * detasinv * Vh(2,iv,3)
        bc(iv,5,5,5) = bc(iv,5,5,5) - alfa * dtl(2,iv) * db5 * detasinv * Vh(2,iv,4)

!.... I term

        mat(2,iv,1,1,2) = mat(2,iv,1,1,2) + one
        mat(2,iv,2,2,2) = mat(2,iv,2,2,2) + one
        mat(2,iv,3,3,2) = mat(2,iv,3,3,2) + one
        mat(2,iv,4,4,2) = mat(2,iv,4,4,2) + one
        mat(2,iv,5,5,2) = mat(2,iv,5,5,2) + one

        end do
        
!.... sponge

        if (ispg .eq. 1) then
        
          do iv = 1, nx
            mat(2,iv,1,1,2) = mat(2,iv,1,1,2) + alfa * dtl(2,iv) * spgl(2,iv)
            mat(2,iv,2,2,2) = mat(2,iv,2,2,2) + alfa * dtl(2,iv) * spgl(2,iv)
            mat(2,iv,3,3,2) = mat(2,iv,3,3,2) + alfa * dtl(2,iv) * spgl(2,iv)
            mat(2,iv,4,4,2) = mat(2,iv,4,4,2) + alfa * dtl(2,iv) * spgl(2,iv)
            mat(2,iv,5,5,2) = mat(2,iv,5,5,2) + alfa * dtl(2,iv) * spgl(2,iv)
          end do
                
        else if (ispg .ge. 2) then
        
          do iv = 1, nx
            mat(2,iv,1,1,2) = mat(2,iv,1,1,2) + alfa * dtl(2,iv) * (spgl(2,iv) + spg2l(2,iv))
            mat(2,iv,2,2,2) = mat(2,iv,2,2,2) + alfa * dtl(2,iv) * (spgl(2,iv) + spg2l(2,iv))
            mat(2,iv,3,3,2) = mat(2,iv,3,3,2) + alfa * dtl(2,iv) * (spgl(2,iv) + spg2l(2,iv))
            mat(2,iv,4,4,2) = mat(2,iv,4,4,2) + alfa * dtl(2,iv) * (spgl(2,iv) + spg2l(2,iv))
            mat(2,iv,5,5,2) = mat(2,iv,5,5,2) + alfa * dtl(2,iv) * (spgl(2,iv) + spg2l(2,iv))
          end do
                
        end if

!.... implicit damping term

        if (eps_e .ne. zero) then

        do iv = 1, nx
        mat(2,iv,1,1,1) = mat(2,iv,1,1,1) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * (-one)
        mat(2,iv,2,2,1) = mat(2,iv,2,2,1) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * (-one)
        mat(2,iv,3,3,1) = mat(2,iv,3,3,1) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * (-one)
        mat(2,iv,4,4,1) = mat(2,iv,4,4,1) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * (-one)
        mat(2,iv,5,5,1) = mat(2,iv,5,5,1) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * (-one)

        mat(2,iv,1,1,2) = mat(2,iv,1,1,2) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * two
        mat(2,iv,2,2,2) = mat(2,iv,2,2,2) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * two
        mat(2,iv,3,3,2) = mat(2,iv,3,3,2) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * two
        mat(2,iv,4,4,2) = mat(2,iv,4,4,2) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * two
        mat(2,iv,5,5,2) = mat(2,iv,5,5,2) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * two

        mat(2,iv,1,1,3) = mat(2,iv,1,1,3) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * (-one)
        mat(2,iv,2,2,3) = mat(2,iv,2,2,3) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * (-one)
        mat(2,iv,3,3,3) = mat(2,iv,3,3,3) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * (-one)
        mat(2,iv,4,4,3) = mat(2,iv,4,4,3) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * (-one)
        mat(2,iv,5,5,3) = mat(2,iv,5,5,3) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * (-one)
        end do
                
        end if

!=======================================================================================================!
!.... use higher-order tangent on second node off the body
!=======================================================================================================!
        do idof = 1, ndof
          do jdof = 1, ndof
            do iv = 1, nx
              bc(iv,idof,jdof,7)    = alfa * dtl(3,iv) * ga1 * detainv * Bh(3,iv,idof,jdof)
              mat(3,iv,idof,jdof,1) = alfa * dtl(3,iv) * ga2 * detainv * Bh(3,iv,idof,jdof)
              mat(3,iv,idof,jdof,2) = zero
              mat(3,iv,idof,jdof,3) = alfa * dtl(3,iv) * ga3 * detainv * Bh(3,iv,idof,jdof)
              bc(iv,idof,jdof,6)    = alfa * dtl(3,iv) * ga4 * detainv * Bh(3,iv,idof,jdof)
            end do
          end do
        end do
        
!.... \hat{V}_{\eta\eta} term

        do iv = 1, nx

        bc(iv,2,2,7)    = bc(iv,2,2,7)    - alfa * dtl(3,iv) * da1 * detasinv * Vh(3,iv,1)
        bc(iv,2,3,7)    = bc(iv,2,3,7)    - alfa * dtl(3,iv) * da1 * detasinv * Vh(3,iv,5)
        bc(iv,3,2,7)    = bc(iv,3,2,7)    - alfa * dtl(3,iv) * da1 * detasinv * Vh(3,iv,6)
        bc(iv,3,3,7)    = bc(iv,3,3,7)    - alfa * dtl(3,iv) * da1 * detasinv * Vh(3,iv,2)
        bc(iv,4,4,7)    = bc(iv,4,4,7)    - alfa * dtl(3,iv) * da1 * detasinv * Vh(3,iv,3)
        bc(iv,5,5,7)    = bc(iv,5,5,7)    - alfa * dtl(3,iv) * da1 * detasinv * Vh(3,iv,4)

        mat(3,iv,2,2,1) = mat(3,iv,2,2,1) - alfa * dtl(3,iv) * da2 * detasinv * Vh(3,iv,1)
        mat(3,iv,2,3,1) = mat(3,iv,2,3,1) - alfa * dtl(3,iv) * da2 * detasinv * Vh(3,iv,5)
        mat(3,iv,3,2,1) = mat(3,iv,3,2,1) - alfa * dtl(3,iv) * da2 * detasinv * Vh(3,iv,6)
        mat(3,iv,3,3,1) = mat(3,iv,3,3,1) - alfa * dtl(3,iv) * da2 * detasinv * Vh(3,iv,2)
        mat(3,iv,4,4,1) = mat(3,iv,4,4,1) - alfa * dtl(3,iv) * da2 * detasinv * Vh(3,iv,3)
        mat(3,iv,5,5,1) = mat(3,iv,5,5,1) - alfa * dtl(3,iv) * da2 * detasinv * Vh(3,iv,4)

        mat(3,iv,2,2,2) = mat(3,iv,2,2,2) - alfa * dtl(3,iv) * da3 * detasinv * Vh(3,iv,1)
        mat(3,iv,2,3,2) = mat(3,iv,2,3,2) - alfa * dtl(3,iv) * da3 * detasinv * Vh(3,iv,5)
        mat(3,iv,3,2,2) = mat(3,iv,3,2,2) - alfa * dtl(3,iv) * da3 * detasinv * Vh(3,iv,6)
        mat(3,iv,3,3,2) = mat(3,iv,3,3,2) - alfa * dtl(3,iv) * da3 * detasinv * Vh(3,iv,2)
        mat(3,iv,4,4,2) = mat(3,iv,4,4,2) - alfa * dtl(3,iv) * da3 * detasinv * Vh(3,iv,3)
        mat(3,iv,5,5,2) = mat(3,iv,5,5,2) - alfa * dtl(3,iv) * da3 * detasinv * Vh(3,iv,4)

        mat(3,iv,2,2,3) = mat(3,iv,2,2,3) - alfa * dtl(3,iv) * da4 * detasinv * Vh(3,iv,1)
        mat(3,iv,2,3,3) = mat(3,iv,2,3,3) - alfa * dtl(3,iv) * da4 * detasinv * Vh(3,iv,5)
        mat(3,iv,3,2,3) = mat(3,iv,3,2,3) - alfa * dtl(3,iv) * da4 * detasinv * Vh(3,iv,6)
        mat(3,iv,3,3,3) = mat(3,iv,3,3,3) - alfa * dtl(3,iv) * da4 * detasinv * Vh(3,iv,2)
        mat(3,iv,4,4,3) = mat(3,iv,4,4,3) - alfa * dtl(3,iv) * da4 * detasinv * Vh(3,iv,3)
        mat(3,iv,5,5,3) = mat(3,iv,5,5,3) - alfa * dtl(3,iv) * da4 * detasinv * Vh(3,iv,4)

        bc(iv,2,2,6)    = bc(iv,2,2,6)    - alfa * dtl(3,iv) * da5 * detasinv * Vh(3,iv,1)
        bc(iv,2,3,6)    = bc(iv,2,3,6)    - alfa * dtl(3,iv) * da5 * detasinv * Vh(3,iv,5)
        bc(iv,3,2,6)    = bc(iv,3,2,6)    - alfa * dtl(3,iv) * da5 * detasinv * Vh(3,iv,6)
        bc(iv,3,3,6)    = bc(iv,3,3,6)    - alfa * dtl(3,iv) * da5 * detasinv * Vh(3,iv,2)
        bc(iv,4,4,6)    = bc(iv,4,4,6)    - alfa * dtl(3,iv) * da5 * detasinv * Vh(3,iv,3)
        bc(iv,5,5,6)    = bc(iv,5,5,6)    - alfa * dtl(3,iv) * da5 * detasinv * Vh(3,iv,4)

!.... I term

        mat(3,iv,1,1,2) = mat(3,iv,1,1,2) + one
        mat(3,iv,2,2,2) = mat(3,iv,2,2,2) + one
        mat(3,iv,3,3,2) = mat(3,iv,3,3,2) + one
        mat(3,iv,4,4,2) = mat(3,iv,4,4,2) + one
        mat(3,iv,5,5,2) = mat(3,iv,5,5,2) + one

        end do
        
!.... sponge

        if (ispg .eq. 1) then
        
          do iv = 1, nx
            mat(3,iv,1,1,2) = mat(3,iv,1,1,2) + alfa * dtl(3,iv) * spgl(3,iv)
            mat(3,iv,2,2,2) = mat(3,iv,2,2,2) + alfa * dtl(3,iv) * spgl(3,iv)
            mat(3,iv,3,3,2) = mat(3,iv,3,3,2) + alfa * dtl(3,iv) * spgl(3,iv)
            mat(3,iv,4,4,2) = mat(3,iv,4,4,2) + alfa * dtl(3,iv) * spgl(3,iv)
            mat(3,iv,5,5,2) = mat(3,iv,5,5,2) + alfa * dtl(3,iv) * spgl(3,iv)
          end do
                
        else if (ispg .eq. 2) then
        
          do iv = 1, nx
            mat(3,iv,1,1,2) = mat(3,iv,1,1,2) + alfa * dtl(3,iv) * (spgl(3,iv) + spg2l(3,iv))
            mat(3,iv,2,2,2) = mat(3,iv,2,2,2) + alfa * dtl(3,iv) * (spgl(3,iv) + spg2l(3,iv))
            mat(3,iv,3,3,2) = mat(3,iv,3,3,2) + alfa * dtl(3,iv) * (spgl(3,iv) + spg2l(3,iv))
            mat(3,iv,4,4,2) = mat(3,iv,4,4,2) + alfa * dtl(3,iv) * (spgl(3,iv) + spg2l(3,iv))
            mat(3,iv,5,5,2) = mat(3,iv,5,5,2) + alfa * dtl(3,iv) * (spgl(3,iv) + spg2l(3,iv))
          end do
                
        end if

!.... implicit damping term updated to fourth order

        if (eps_e .ne. zero) then

        do iv = 1, nx

        bc(iv,1,1,7) = bc(iv,1,1,7) + alfa * dtl(3,iv) * eps_e * buff(3,iv) * fb1
        bc(iv,2,2,7) = bc(iv,2,2,7) + alfa * dtl(3,iv) * eps_e * buff(3,iv) * fb1
        bc(iv,3,3,7) = bc(iv,3,3,7) + alfa * dtl(3,iv) * eps_e * buff(3,iv) * fb1
        bc(iv,4,4,7) = bc(iv,4,4,7) + alfa * dtl(3,iv) * eps_e * buff(3,iv) * fb1
        bc(iv,5,5,7) = bc(iv,5,5,7) + alfa * dtl(3,iv) * eps_e * buff(3,iv) * fb1

        mat(3,iv,1,1,1) = mat(3,iv,1,1,1) + alfa * dtl(3,iv) * eps_e * buff(3,iv) * fb2
        mat(3,iv,2,2,1) = mat(3,iv,2,2,1) + alfa * dtl(3,iv) * eps_e * buff(3,iv) * fb2
        mat(3,iv,3,3,1) = mat(3,iv,3,3,1) + alfa * dtl(3,iv) * eps_e * buff(3,iv) * fb2
        mat(3,iv,4,4,1) = mat(3,iv,4,4,1) + alfa * dtl(3,iv) * eps_e * buff(3,iv) * fb2
        mat(3,iv,5,5,1) = mat(3,iv,5,5,1) + alfa * dtl(3,iv) * eps_e * buff(3,iv) * fb2

        mat(3,iv,1,1,2) = mat(3,iv,1,1,2) + alfa * dtl(3,iv) * eps_e * buff(3,iv) * fb3
        mat(3,iv,2,2,2) = mat(3,iv,2,2,2) + alfa * dtl(3,iv) * eps_e * buff(3,iv) * fb3
        mat(3,iv,3,3,2) = mat(3,iv,3,3,2) + alfa * dtl(3,iv) * eps_e * buff(3,iv) * fb3
        mat(3,iv,4,4,2) = mat(3,iv,4,4,2) + alfa * dtl(3,iv) * eps_e * buff(3,iv) * fb3
        mat(3,iv,5,5,2) = mat(3,iv,5,5,2) + alfa * dtl(3,iv) * eps_e * buff(3,iv) * fb3

        mat(3,iv,1,1,3) = mat(3,iv,1,1,3) + alfa * dtl(3,iv) * eps_e * buff(3,iv) * fb4
        mat(3,iv,2,2,3) = mat(3,iv,2,2,3) + alfa * dtl(3,iv) * eps_e * buff(3,iv) * fb4
        mat(3,iv,3,3,3) = mat(3,iv,3,3,3) + alfa * dtl(3,iv) * eps_e * buff(3,iv) * fb4
        mat(3,iv,4,4,3) = mat(3,iv,4,4,3) + alfa * dtl(3,iv) * eps_e * buff(3,iv) * fb4
        mat(3,iv,5,5,3) = mat(3,iv,5,5,3) + alfa * dtl(3,iv) * eps_e * buff(3,iv) * fb4

        bc(iv,1,1,6) = bc(iv,1,1,6) + alfa * dtl(3,iv) * eps_e * buff(3,iv) * fb5
        bc(iv,2,2,6) = bc(iv,2,2,6) + alfa * dtl(3,iv) * eps_e * buff(3,iv) * fb5
        bc(iv,3,3,6) = bc(iv,3,3,6) + alfa * dtl(3,iv) * eps_e * buff(3,iv) * fb5
        bc(iv,4,4,6) = bc(iv,4,4,6) + alfa * dtl(3,iv) * eps_e * buff(3,iv) * fb5
        bc(iv,5,5,6) = bc(iv,5,5,6) + alfa * dtl(3,iv) * eps_e * buff(3,iv) * fb5

        end do

        end if

!=======================================================================================================!
!.... use higher-order tangent on second node off the far-field
!=======================================================================================================!
        do idof = 1, ndof
          do jdof = 1, ndof
            do iv = 1, nx
              bc(iv,idof,jdof,8)       = alfa * dtl(ny-2,iv) * ga1 * detainv * Bh(ny-2,iv,idof,jdof)
              mat(ny-2,iv,idof,jdof,1) = alfa * dtl(ny-2,iv) * ga2 * detainv * Bh(ny-2,iv,idof,jdof)
              mat(ny-2,iv,idof,jdof,2) = zero
              mat(ny-2,iv,idof,jdof,3) = alfa * dtl(ny-2,iv) * ga3 * detainv * Bh(ny-2,iv,idof,jdof)
              bc(iv,idof,jdof,9)       = alfa * dtl(ny-2,iv) * ga4 * detainv * Bh(ny-2,iv,idof,jdof)
            end do
          end do
        end do
        
!.... \hat{V}_{\eta\eta} term

        do iv = 1, nx

        bc(iv,2,2,8)       = bc(iv,2,2,8)       - alfa * dtl(ny-2,iv) * da1 * detasinv * Vh(ny-2,iv,1)
        bc(iv,2,3,8)       = bc(iv,2,3,8)       - alfa * dtl(ny-2,iv) * da1 * detasinv * Vh(ny-2,iv,5)
        bc(iv,3,2,8)       = bc(iv,3,2,8)       - alfa * dtl(ny-2,iv) * da1 * detasinv * Vh(ny-2,iv,6)
        bc(iv,3,3,8)       = bc(iv,3,3,8)       - alfa * dtl(ny-2,iv) * da1 * detasinv * Vh(ny-2,iv,2)
        bc(iv,4,4,8)       = bc(iv,4,4,8)       - alfa * dtl(ny-2,iv) * da1 * detasinv * Vh(ny-2,iv,3)
        bc(iv,5,5,8)       = bc(iv,5,5,8)       - alfa * dtl(ny-2,iv) * da1 * detasinv * Vh(ny-2,iv,4)
 
        mat(ny-2,iv,2,2,1) = mat(ny-2,iv,2,2,1) - alfa * dtl(ny-2,iv) * da2 * detasinv * Vh(ny-2,iv,1)
        mat(ny-2,iv,2,3,1) = mat(ny-2,iv,2,3,1) - alfa * dtl(ny-2,iv) * da2 * detasinv * Vh(ny-2,iv,5)
        mat(ny-2,iv,3,2,1) = mat(ny-2,iv,3,2,1) - alfa * dtl(ny-2,iv) * da2 * detasinv * Vh(ny-2,iv,6)
        mat(ny-2,iv,3,3,1) = mat(ny-2,iv,3,3,1) - alfa * dtl(ny-2,iv) * da2 * detasinv * Vh(ny-2,iv,2)
        mat(ny-2,iv,4,4,1) = mat(ny-2,iv,4,4,1) - alfa * dtl(ny-2,iv) * da2 * detasinv * Vh(ny-2,iv,3)
        mat(ny-2,iv,5,5,1) = mat(ny-2,iv,5,5,1) - alfa * dtl(ny-2,iv) * da2 * detasinv * Vh(ny-2,iv,4)

        mat(ny-2,iv,2,2,2) = mat(ny-2,iv,2,2,2) - alfa * dtl(ny-2,iv) * da3 * detasinv * Vh(ny-2,iv,1)
        mat(ny-2,iv,2,3,2) = mat(ny-2,iv,2,3,2) - alfa * dtl(ny-2,iv) * da3 * detasinv * Vh(ny-2,iv,5)
        mat(ny-2,iv,3,2,2) = mat(ny-2,iv,3,2,2) - alfa * dtl(ny-2,iv) * da3 * detasinv * Vh(ny-2,iv,6)
        mat(ny-2,iv,3,3,2) = mat(ny-2,iv,3,3,2) - alfa * dtl(ny-2,iv) * da3 * detasinv * Vh(ny-2,iv,2)
        mat(ny-2,iv,4,4,2) = mat(ny-2,iv,4,4,2) - alfa * dtl(ny-2,iv) * da3 * detasinv * Vh(ny-2,iv,3)
        mat(ny-2,iv,5,5,2) = mat(ny-2,iv,5,5,2) - alfa * dtl(ny-2,iv) * da3 * detasinv * Vh(ny-2,iv,4)

        mat(ny-2,iv,2,2,3) = mat(ny-2,iv,2,2,3) - alfa * dtl(ny-2,iv) * da4 * detasinv * Vh(ny-2,iv,1)
        mat(ny-2,iv,2,3,3) = mat(ny-2,iv,2,3,3) - alfa * dtl(ny-2,iv) * da4 * detasinv * Vh(ny-2,iv,5)
        mat(ny-2,iv,3,2,3) = mat(ny-2,iv,3,2,3) - alfa * dtl(ny-2,iv) * da4 * detasinv * Vh(ny-2,iv,6)
        mat(ny-2,iv,3,3,3) = mat(ny-2,iv,3,3,3) - alfa * dtl(ny-2,iv) * da4 * detasinv * Vh(ny-2,iv,2)
        mat(ny-2,iv,4,4,3) = mat(ny-2,iv,4,4,3) - alfa * dtl(ny-2,iv) * da4 * detasinv * Vh(ny-2,iv,3)
        mat(ny-2,iv,5,5,3) = mat(ny-2,iv,5,5,3) - alfa * dtl(ny-2,iv) * da4 * detasinv * Vh(ny-2,iv,4)

        bc(iv,2,2,9)       = bc(iv,2,2,9)       - alfa * dtl(ny-2,iv) * da5 * detasinv * Vh(ny-2,iv,1)
        bc(iv,2,3,9)       = bc(iv,2,3,9)       - alfa * dtl(ny-2,iv) * da5 * detasinv * Vh(ny-2,iv,5)
        bc(iv,3,2,9)       = bc(iv,3,2,9)       - alfa * dtl(ny-2,iv) * da5 * detasinv * Vh(ny-2,iv,6)
        bc(iv,3,3,9)       = bc(iv,3,3,9)       - alfa * dtl(ny-2,iv) * da5 * detasinv * Vh(ny-2,iv,2)
        bc(iv,4,4,9)       = bc(iv,4,4,9)       - alfa * dtl(ny-2,iv) * da5 * detasinv * Vh(ny-2,iv,3)
        bc(iv,5,5,9)       = bc(iv,5,5,9)       - alfa * dtl(ny-2,iv) * da5 * detasinv * Vh(ny-2,iv,4)

!.... I term

        mat(ny-2,iv,1,1,2) = mat(ny-2,iv,1,1,2) + one
        mat(ny-2,iv,2,2,2) = mat(ny-2,iv,2,2,2) + one
        mat(ny-2,iv,3,3,2) = mat(ny-2,iv,3,3,2) + one
        mat(ny-2,iv,4,4,2) = mat(ny-2,iv,4,4,2) + one
        mat(ny-2,iv,5,5,2) = mat(ny-2,iv,5,5,2) + one

        end do
        
!.... sponge

        if (ispg .eq. 1) then
        
        do iv = 1, nx
          mat(ny-2,iv,1,1,2) = mat(ny-2,iv,1,1,2) + alfa * dtl(ny-2,iv) * spgl(ny-2,iv)
          mat(ny-2,iv,2,2,2) = mat(ny-2,iv,2,2,2) + alfa * dtl(ny-2,iv) * spgl(ny-2,iv)
          mat(ny-2,iv,3,3,2) = mat(ny-2,iv,3,3,2) + alfa * dtl(ny-2,iv) * spgl(ny-2,iv)
          mat(ny-2,iv,4,4,2) = mat(ny-2,iv,4,4,2) + alfa * dtl(ny-2,iv) * spgl(ny-2,iv)
          mat(ny-2,iv,5,5,2) = mat(ny-2,iv,5,5,2) + alfa * dtl(ny-2,iv) * spgl(ny-2,iv)
        end do
                
        else if (ispg .eq. 2) then
        
        do iv = 1, nx
        mat(ny-2,iv,1,1,2) = mat(ny-2,iv,1,1,2) + alfa * dtl(ny-2,iv) * (spgl(ny-2,iv) + spg2l(ny-2,iv))
        mat(ny-2,iv,2,2,2) = mat(ny-2,iv,2,2,2) + alfa * dtl(ny-2,iv) * (spgl(ny-2,iv) + spg2l(ny-2,iv))
        mat(ny-2,iv,3,3,2) = mat(ny-2,iv,3,3,2) + alfa * dtl(ny-2,iv) * (spgl(ny-2,iv) + spg2l(ny-2,iv))
        mat(ny-2,iv,4,4,2) = mat(ny-2,iv,4,4,2) + alfa * dtl(ny-2,iv) * (spgl(ny-2,iv) + spg2l(ny-2,iv))
        mat(ny-2,iv,5,5,2) = mat(ny-2,iv,5,5,2) + alfa * dtl(ny-2,iv) * (spgl(ny-2,iv) + spg2l(ny-2,iv))
        end do
                
        end if

!.... implicit damping term updated to fourth order

        if (eps_e .ne. zero) then

        do iv = 1, nx

        bc(iv,1,1,8) = bc(iv,1,1,8) + alfa * dtl(ny-2,iv) * eps_e * buff(ny-2,iv) * fb1
        bc(iv,2,2,8) = bc(iv,2,2,8) + alfa * dtl(ny-2,iv) * eps_e * buff(ny-2,iv) * fb1
        bc(iv,3,3,8) = bc(iv,3,3,8) + alfa * dtl(ny-2,iv) * eps_e * buff(ny-2,iv) * fb1
        bc(iv,4,4,8) = bc(iv,4,4,8) + alfa * dtl(ny-2,iv) * eps_e * buff(ny-2,iv) * fb1
        bc(iv,5,5,8) = bc(iv,5,5,8) + alfa * dtl(ny-2,iv) * eps_e * buff(ny-2,iv) * fb1

        mat(ny-2,iv,1,1,1) = mat(ny-2,iv,1,1,1) + alfa * dtl(ny-2,iv) * eps_e * buff(ny-2,iv) * fb2
        mat(ny-2,iv,2,2,1) = mat(ny-2,iv,2,2,1) + alfa * dtl(ny-2,iv) * eps_e * buff(ny-2,iv) * fb2
        mat(ny-2,iv,3,3,1) = mat(ny-2,iv,3,3,1) + alfa * dtl(ny-2,iv) * eps_e * buff(ny-2,iv) * fb2
        mat(ny-2,iv,4,4,1) = mat(ny-2,iv,4,4,1) + alfa * dtl(ny-2,iv) * eps_e * buff(ny-2,iv) * fb2
        mat(ny-2,iv,5,5,1) = mat(ny-2,iv,5,5,1) + alfa * dtl(ny-2,iv) * eps_e * buff(ny-2,iv) * fb2

        mat(ny-2,iv,1,1,2) = mat(ny-2,iv,1,1,2) + alfa * dtl(ny-2,iv) * eps_e * buff(ny-2,iv) * fb3
        mat(ny-2,iv,2,2,2) = mat(ny-2,iv,2,2,2) + alfa * dtl(ny-2,iv) * eps_e * buff(ny-2,iv) * fb3
        mat(ny-2,iv,3,3,2) = mat(ny-2,iv,3,3,2) + alfa * dtl(ny-2,iv) * eps_e * buff(ny-2,iv) * fb3
        mat(ny-2,iv,4,4,2) = mat(ny-2,iv,4,4,2) + alfa * dtl(ny-2,iv) * eps_e * buff(ny-2,iv) * fb3
        mat(ny-2,iv,5,5,2) = mat(ny-2,iv,5,5,2) + alfa * dtl(ny-2,iv) * eps_e * buff(ny-2,iv) * fb3

        mat(ny-2,iv,1,1,3) = mat(ny-2,iv,1,1,3) + alfa * dtl(ny-2,iv) * eps_e * buff(ny-2,iv) * fb4
        mat(ny-2,iv,2,2,3) = mat(ny-2,iv,2,2,3) + alfa * dtl(ny-2,iv) * eps_e * buff(ny-2,iv) * fb4
        mat(ny-2,iv,3,3,3) = mat(ny-2,iv,3,3,3) + alfa * dtl(ny-2,iv) * eps_e * buff(ny-2,iv) * fb4
        mat(ny-2,iv,4,4,3) = mat(ny-2,iv,4,4,3) + alfa * dtl(ny-2,iv) * eps_e * buff(ny-2,iv) * fb4
        mat(ny-2,iv,5,5,3) = mat(ny-2,iv,5,5,3) + alfa * dtl(ny-2,iv) * eps_e * buff(ny-2,iv) * fb4

        bc(iv,1,1,9) = bc(iv,1,1,9) + alfa * dtl(ny-2,iv) * eps_e * buff(ny-2,iv) * fb5
        bc(iv,2,2,9) = bc(iv,2,2,9) + alfa * dtl(ny-2,iv) * eps_e * buff(ny-2,iv) * fb5
        bc(iv,3,3,9) = bc(iv,3,3,9) + alfa * dtl(ny-2,iv) * eps_e * buff(ny-2,iv) * fb5
        bc(iv,4,4,9) = bc(iv,4,4,9) + alfa * dtl(ny-2,iv) * eps_e * buff(ny-2,iv) * fb5
        bc(iv,5,5,9) = bc(iv,5,5,9) + alfa * dtl(ny-2,iv) * eps_e * buff(ny-2,iv) * fb5

        end do

        end if

!=======================================================================================================!
!.... use higher-order tangent on first node off the far-field boundary
!=======================================================================================================!
        do idof = 1, ndof
          do jdof = 1, ndof
            do iv = 1, nx
              bc(iv,idof,jdof,10)      = -alfa * dtl(ny-1,iv) * gb5 * detainv * Bh(ny-1,iv,idof,jdof)
              bc(iv,idof,jdof,11)      = -alfa * dtl(ny-1,iv) * gb4 * detainv * Bh(ny-1,iv,idof,jdof)
              mat(ny-1,iv,idof,jdof,1) = -alfa * dtl(ny-1,iv) * gb3 * detainv * Bh(ny-1,iv,idof,jdof)
              mat(ny-1,iv,idof,jdof,2) = -alfa * dtl(ny-1,iv) * gb2 * detainv * Bh(ny-1,iv,idof,jdof)
              mat(ny-1,iv,idof,jdof,3) = -alfa * dtl(ny-1,iv) * gb1 * detainv * Bh(ny-1,iv,idof,jdof)
            end do
          end do
        end do
        
!.... \hat{V}_{\eta\eta} term

        do iv = 1, nx

        bc(iv,2,2,10) = bc(iv,2,2,10) - alfa * dtl(ny-1,iv) * db5 * detasinv * Vh(ny-1,iv,1)
        bc(iv,2,3,10) = bc(iv,2,3,10) - alfa * dtl(ny-1,iv) * db5 * detasinv * Vh(ny-1,iv,5)
        bc(iv,3,2,10) = bc(iv,3,2,10) - alfa * dtl(ny-1,iv) * db5 * detasinv * Vh(ny-1,iv,6)
        bc(iv,3,3,10) = bc(iv,3,3,10) - alfa * dtl(ny-1,iv) * db5 * detasinv * Vh(ny-1,iv,2)
        bc(iv,4,4,10) = bc(iv,4,4,10) - alfa * dtl(ny-1,iv) * db5 * detasinv * Vh(ny-1,iv,3)
        bc(iv,5,5,10) = bc(iv,5,5,10) - alfa * dtl(ny-1,iv) * db5 * detasinv * Vh(ny-1,iv,4)

        bc(iv,2,2,11) = bc(iv,2,2,11) - alfa * dtl(ny-1,iv) * db4 * detasinv * Vh(ny-1,iv,1)
        bc(iv,2,3,11) = bc(iv,2,3,11) - alfa * dtl(ny-1,iv) * db4 * detasinv * Vh(ny-1,iv,5)
        bc(iv,3,2,11) = bc(iv,3,2,11) - alfa * dtl(ny-1,iv) * db4 * detasinv * Vh(ny-1,iv,6)
        bc(iv,3,3,11) = bc(iv,3,3,11) - alfa * dtl(ny-1,iv) * db4 * detasinv * Vh(ny-1,iv,2)
        bc(iv,4,4,11) = bc(iv,4,4,11) - alfa * dtl(ny-1,iv) * db4 * detasinv * Vh(ny-1,iv,3)
        bc(iv,5,5,11) = bc(iv,5,5,11) - alfa * dtl(ny-1,iv) * db4 * detasinv * Vh(ny-1,iv,4)

        mat(ny-1,iv,2,2,1) = mat(ny-1,iv,2,2,1) - alfa * dtl(ny-1,iv) * db3 * detasinv * Vh(ny-1,iv,1)
        mat(ny-1,iv,2,3,1) = mat(ny-1,iv,2,3,1) - alfa * dtl(ny-1,iv) * db3 * detasinv * Vh(ny-1,iv,5)
        mat(ny-1,iv,3,2,1) = mat(ny-1,iv,3,2,1) - alfa * dtl(ny-1,iv) * db3 * detasinv * Vh(ny-1,iv,6)
        mat(ny-1,iv,3,3,1) = mat(ny-1,iv,3,3,1) - alfa * dtl(ny-1,iv) * db3 * detasinv * Vh(ny-1,iv,2)
        mat(ny-1,iv,4,4,1) = mat(ny-1,iv,4,4,1) - alfa * dtl(ny-1,iv) * db3 * detasinv * Vh(ny-1,iv,3)
        mat(ny-1,iv,5,5,1) = mat(ny-1,iv,5,5,1) - alfa * dtl(ny-1,iv) * db3 * detasinv * Vh(ny-1,iv,4)

        mat(ny-1,iv,2,2,2) = mat(ny-1,iv,2,2,2) - alfa * dtl(ny-1,iv) * db2 * detasinv * Vh(ny-1,iv,1)
        mat(ny-1,iv,2,3,2) = mat(ny-1,iv,2,3,2) - alfa * dtl(ny-1,iv) * db2 * detasinv * Vh(ny-1,iv,5)
        mat(ny-1,iv,3,2,2) = mat(ny-1,iv,3,2,2) - alfa * dtl(ny-1,iv) * db2 * detasinv * Vh(ny-1,iv,6)
        mat(ny-1,iv,3,3,2) = mat(ny-1,iv,3,3,2) - alfa * dtl(ny-1,iv) * db2 * detasinv * Vh(ny-1,iv,2)
        mat(ny-1,iv,4,4,2) = mat(ny-1,iv,4,4,2) - alfa * dtl(ny-1,iv) * db2 * detasinv * Vh(ny-1,iv,3)
        mat(ny-1,iv,5,5,2) = mat(ny-1,iv,5,5,2) - alfa * dtl(ny-1,iv) * db2 * detasinv * Vh(ny-1,iv,4)

        mat(ny-1,iv,2,2,3) = mat(ny-1,iv,2,2,3) - alfa * dtl(ny-1,iv) * db1 * detasinv * Vh(ny-1,iv,1)
        mat(ny-1,iv,2,3,3) = mat(ny-1,iv,2,3,3) - alfa * dtl(ny-1,iv) * db1 * detasinv * Vh(ny-1,iv,5)
        mat(ny-1,iv,3,2,3) = mat(ny-1,iv,3,2,3) - alfa * dtl(ny-1,iv) * db1 * detasinv * Vh(ny-1,iv,6)
        mat(ny-1,iv,3,3,3) = mat(ny-1,iv,3,3,3) - alfa * dtl(ny-1,iv) * db1 * detasinv * Vh(ny-1,iv,2)
        mat(ny-1,iv,4,4,3) = mat(ny-1,iv,4,4,3) - alfa * dtl(ny-1,iv) * db1 * detasinv * Vh(ny-1,iv,3)
        mat(ny-1,iv,5,5,3) = mat(ny-1,iv,5,5,3) - alfa * dtl(ny-1,iv) * db1 * detasinv * Vh(ny-1,iv,4)

!.... I term

        mat(ny-1,iv,1,1,2) = mat(ny-1,iv,1,1,2) + one
        mat(ny-1,iv,2,2,2) = mat(ny-1,iv,2,2,2) + one
        mat(ny-1,iv,3,3,2) = mat(ny-1,iv,3,3,2) + one
        mat(ny-1,iv,4,4,2) = mat(ny-1,iv,4,4,2) + one
        mat(ny-1,iv,5,5,2) = mat(ny-1,iv,5,5,2) + one

        end do
        
!.... sponge

        if (ispg .eq. 1) then

        do iv = 1, nx
          mat(ny-1,iv,1,1,2) = mat(ny-1,iv,1,1,2) + alfa * dtl(ny-1,iv) * spgl(ny-1,iv)
          mat(ny-1,iv,2,2,2) = mat(ny-1,iv,2,2,2) + alfa * dtl(ny-1,iv) * spgl(ny-1,iv)
          mat(ny-1,iv,3,3,2) = mat(ny-1,iv,3,3,2) + alfa * dtl(ny-1,iv) * spgl(ny-1,iv)
          mat(ny-1,iv,4,4,2) = mat(ny-1,iv,4,4,2) + alfa * dtl(ny-1,iv) * spgl(ny-1,iv)
          mat(ny-1,iv,5,5,2) = mat(ny-1,iv,5,5,2) + alfa * dtl(ny-1,iv) * spgl(ny-1,iv)
        end do

        else if (ispg .eq. 2) then

        do iv = 1, nx
        mat(ny-1,iv,1,1,2) = mat(ny-1,iv,1,1,2) + alfa * dtl(ny-1,iv) * (spgl(ny-1,iv) + spg2l(ny-1,iv))
        mat(ny-1,iv,2,2,2) = mat(ny-1,iv,2,2,2) + alfa * dtl(ny-1,iv) * (spgl(ny-1,iv) + spg2l(ny-1,iv))
        mat(ny-1,iv,3,3,2) = mat(ny-1,iv,3,3,2) + alfa * dtl(ny-1,iv) * (spgl(ny-1,iv) + spg2l(ny-1,iv))
        mat(ny-1,iv,4,4,2) = mat(ny-1,iv,4,4,2) + alfa * dtl(ny-1,iv) * (spgl(ny-1,iv) + spg2l(ny-1,iv))
        mat(ny-1,iv,5,5,2) = mat(ny-1,iv,5,5,2) + alfa * dtl(ny-1,iv) * (spgl(ny-1,iv) + spg2l(ny-1,iv))
        end do

        end if

!.... implicit damping term

        if (eps_e .ne. zero) then

        do iv = 1, nx
        mat(ny-1,iv,1,1,1) = mat(ny-1,iv,1,1,1) + alfa * dtl(ny-1,iv) * eps_e * buff(ny-1,iv) * (-one)
        mat(ny-1,iv,2,2,1) = mat(ny-1,iv,2,2,1) + alfa * dtl(ny-1,iv) * eps_e * buff(ny-1,iv) * (-one)
        mat(ny-1,iv,3,3,1) = mat(ny-1,iv,3,3,1) + alfa * dtl(ny-1,iv) * eps_e * buff(ny-1,iv) * (-one)
        mat(ny-1,iv,4,4,1) = mat(ny-1,iv,4,4,1) + alfa * dtl(ny-1,iv) * eps_e * buff(ny-1,iv) * (-one)
        mat(ny-1,iv,5,5,1) = mat(ny-1,iv,5,5,1) + alfa * dtl(ny-1,iv) * eps_e * buff(ny-1,iv) * (-one)

        mat(ny-1,iv,1,1,2) = mat(ny-1,iv,1,1,2) + alfa * dtl(ny-1,iv) * eps_e * buff(ny-1,iv) * two
        mat(ny-1,iv,2,2,2) = mat(ny-1,iv,2,2,2) + alfa * dtl(ny-1,iv) * eps_e * buff(ny-1,iv) * two
        mat(ny-1,iv,3,3,2) = mat(ny-1,iv,3,3,2) + alfa * dtl(ny-1,iv) * eps_e * buff(ny-1,iv) * two
        mat(ny-1,iv,4,4,2) = mat(ny-1,iv,4,4,2) + alfa * dtl(ny-1,iv) * eps_e * buff(ny-1,iv) * two
        mat(ny-1,iv,5,5,2) = mat(ny-1,iv,5,5,2) + alfa * dtl(ny-1,iv) * eps_e * buff(ny-1,iv) * two

        mat(ny-1,iv,1,1,3) = mat(ny-1,iv,1,1,3) + alfa * dtl(ny-1,iv) * eps_e * buff(ny-1,iv) * (-one)
        mat(ny-1,iv,2,2,3) = mat(ny-1,iv,2,2,3) + alfa * dtl(ny-1,iv) * eps_e * buff(ny-1,iv) * (-one)
        mat(ny-1,iv,3,3,3) = mat(ny-1,iv,3,3,3) + alfa * dtl(ny-1,iv) * eps_e * buff(ny-1,iv) * (-one)
        mat(ny-1,iv,4,4,3) = mat(ny-1,iv,4,4,3) + alfa * dtl(ny-1,iv) * eps_e * buff(ny-1,iv) * (-one)
        mat(ny-1,iv,5,5,3) = mat(ny-1,iv,5,5,3) + alfa * dtl(ny-1,iv) * eps_e * buff(ny-1,iv) * (-one)
        end do
                
        end if

!=======================================================================================================!
!.... use higher-order tangent on the far-field boundary
!=======================================================================================================!
        do idof = 1, ndof
          do jdof = 1, ndof
            do iv = 1, nx
              bc(iv,idof,jdof,12)    = -alfa * dtl(ny,iv) * gc5 * detainv * Bh(ny,iv,idof,jdof)
              bc(iv,idof,jdof,13)    = -alfa * dtl(ny,iv) * gc4 * detainv * Bh(ny,iv,idof,jdof)
              bc(iv,idof,jdof,14)    = -alfa * dtl(ny,iv) * gc3 * detainv * Bh(ny,iv,idof,jdof)
              mat(ny,iv,idof,jdof,1) = -alfa * dtl(ny,iv) * gc2 * detainv * Bh(ny,iv,idof,jdof)
              mat(ny,iv,idof,jdof,2) = -alfa * dtl(ny,iv) * gc1 * detainv * Bh(ny,iv,idof,jdof)
              mat(ny,iv,idof,jdof,3) = zero
            end do
          end do
        end do
        
!.... \hat{V}_{\eta\eta} term

        do iv = 1, nx

        bc(iv,2,2,12) = bc(iv,2,2,12) - alfa * dtl(ny,iv) * dd5 * detasinv * Vh(ny,iv,1)
        bc(iv,2,3,12) = bc(iv,2,3,12) - alfa * dtl(ny,iv) * dd5 * detasinv * Vh(ny,iv,5)
        bc(iv,3,2,12) = bc(iv,3,2,12) - alfa * dtl(ny,iv) * dd5 * detasinv * Vh(ny,iv,6)
        bc(iv,3,3,12) = bc(iv,3,3,12) - alfa * dtl(ny,iv) * dd5 * detasinv * Vh(ny,iv,2)
        bc(iv,4,4,12) = bc(iv,4,4,12) - alfa * dtl(ny,iv) * dd5 * detasinv * Vh(ny,iv,3)
        bc(iv,5,5,12) = bc(iv,5,5,12) - alfa * dtl(ny,iv) * dd5 * detasinv * Vh(ny,iv,4)

        bc(iv,2,2,13) = bc(iv,2,2,13) - alfa * dtl(ny,iv) * dd4 * detasinv * Vh(ny,iv,1)
        bc(iv,2,3,13) = bc(iv,2,3,13) - alfa * dtl(ny,iv) * dd4 * detasinv * Vh(ny,iv,5)
        bc(iv,3,2,13) = bc(iv,3,2,13) - alfa * dtl(ny,iv) * dd4 * detasinv * Vh(ny,iv,6)
        bc(iv,3,3,13) = bc(iv,3,3,13) - alfa * dtl(ny,iv) * dd4 * detasinv * Vh(ny,iv,2)
        bc(iv,4,4,13) = bc(iv,4,4,13) - alfa * dtl(ny,iv) * dd4 * detasinv * Vh(ny,iv,3)
        bc(iv,5,5,13) = bc(iv,5,5,13) - alfa * dtl(ny,iv) * dd4 * detasinv * Vh(ny,iv,4)

        bc(iv,2,2,14) = bc(iv,2,2,14) - alfa * dtl(ny,iv) * dd3 * detasinv * Vh(ny,iv,1)
        bc(iv,2,3,14) = bc(iv,2,3,14) - alfa * dtl(ny,iv) * dd3 * detasinv * Vh(ny,iv,5)
        bc(iv,3,2,14) = bc(iv,3,2,14) - alfa * dtl(ny,iv) * dd3 * detasinv * Vh(ny,iv,6)
        bc(iv,3,3,14) = bc(iv,3,3,14) - alfa * dtl(ny,iv) * dd3 * detasinv * Vh(ny,iv,2)
        bc(iv,4,4,14) = bc(iv,4,4,14) - alfa * dtl(ny,iv) * dd3 * detasinv * Vh(ny,iv,3)
        bc(iv,5,5,14) = bc(iv,5,5,14) - alfa * dtl(ny,iv) * dd3 * detasinv * Vh(ny,iv,4)

        mat(ny,iv,2,2,1) = mat(ny,iv,2,2,1) - alfa * dtl(ny,iv) * dd2 * detasinv * Vh(ny,iv,1)
        mat(ny,iv,2,3,1) = mat(ny,iv,2,3,1) - alfa * dtl(ny,iv) * dd2 * detasinv * Vh(ny,iv,5)
        mat(ny,iv,3,2,1) = mat(ny,iv,3,2,1) - alfa * dtl(ny,iv) * dd2 * detasinv * Vh(ny,iv,6)
        mat(ny,iv,3,3,1) = mat(ny,iv,3,3,1) - alfa * dtl(ny,iv) * dd2 * detasinv * Vh(ny,iv,2)
        mat(ny,iv,4,4,1) = mat(ny,iv,4,4,1) - alfa * dtl(ny,iv) * dd2 * detasinv * Vh(ny,iv,3)
        mat(ny,iv,5,5,1) = mat(ny,iv,5,5,1) - alfa * dtl(ny,iv) * dd2 * detasinv * Vh(ny,iv,4)

        mat(ny,iv,2,2,2) = mat(ny,iv,2,2,2) - alfa * dtl(ny,iv) * dd1 * detasinv * Vh(ny,iv,1)
        mat(ny,iv,2,3,2) = mat(ny,iv,2,3,2) - alfa * dtl(ny,iv) * dd1 * detasinv * Vh(ny,iv,5)
        mat(ny,iv,3,2,2) = mat(ny,iv,3,2,2) - alfa * dtl(ny,iv) * dd1 * detasinv * Vh(ny,iv,6)
        mat(ny,iv,3,3,2) = mat(ny,iv,3,3,2) - alfa * dtl(ny,iv) * dd1 * detasinv * Vh(ny,iv,2)
        mat(ny,iv,4,4,2) = mat(ny,iv,4,4,2) - alfa * dtl(ny,iv) * dd1 * detasinv * Vh(ny,iv,3)
        mat(ny,iv,5,5,2) = mat(ny,iv,5,5,2) - alfa * dtl(ny,iv) * dd1 * detasinv * Vh(ny,iv,4)

!.... I term

        mat(ny,iv,1,1,2) = mat(ny,iv,1,1,2) + one
        mat(ny,iv,2,2,2) = mat(ny,iv,2,2,2) + one
        mat(ny,iv,3,3,2) = mat(ny,iv,3,3,2) + one
        mat(ny,iv,4,4,2) = mat(ny,iv,4,4,2) + one
        mat(ny,iv,5,5,2) = mat(ny,iv,5,5,2) + one

        end do
        
!.... sponge

        if (ispg .eq. 1) then
        
          do iv = 1, nx
            mat(ny,iv,1,1,2) = mat(ny,iv,1,1,2) + alfa * dtl(ny,iv) * spgl(ny,iv)
            mat(ny,iv,2,2,2) = mat(ny,iv,2,2,2) + alfa * dtl(ny,iv) * spgl(ny,iv)
            mat(ny,iv,3,3,2) = mat(ny,iv,3,3,2) + alfa * dtl(ny,iv) * spgl(ny,iv)
            mat(ny,iv,4,4,2) = mat(ny,iv,4,4,2) + alfa * dtl(ny,iv) * spgl(ny,iv)
            mat(ny,iv,5,5,2) = mat(ny,iv,5,5,2) + alfa * dtl(ny,iv) * spgl(ny,iv)
          end do
      
        else if (ispg .eq. 2) then
        
          do iv = 1, nx
            mat(ny,iv,1,1,2) = mat(ny,iv,1,1,2) + alfa * dtl(ny,iv) * (spgl(ny,iv) + spg2l(ny,iv))
            mat(ny,iv,2,2,2) = mat(ny,iv,2,2,2) + alfa * dtl(ny,iv) * (spgl(ny,iv) + spg2l(ny,iv))
            mat(ny,iv,3,3,2) = mat(ny,iv,3,3,2) + alfa * dtl(ny,iv) * (spgl(ny,iv) + spg2l(ny,iv))
            mat(ny,iv,4,4,2) = mat(ny,iv,4,4,2) + alfa * dtl(ny,iv) * (spgl(ny,iv) + spg2l(ny,iv))
            mat(ny,iv,5,5,2) = mat(ny,iv,5,5,2) + alfa * dtl(ny,iv) * (spgl(ny,iv) + spg2l(ny,iv))
          end do
                
        end if

        end if          ! yper

        return
        end
