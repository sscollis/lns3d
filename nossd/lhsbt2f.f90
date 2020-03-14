!=======================================================================================================!
        subroutine lhsbt2f( mat, Bh, Dh, Vh, spgl, spg2l, dtl )
!  
!  Correct the LHS for boundary treatment in the \eta direction
!
!  This version supports the 4th order LHS
!
!  Revised: 10-17-95
!
!=======================================================================================================!
        use global
        use stencil
        use buff_mod
        implicit none
        
        real :: mat(ny,nx,ndof,ndof,5), Bh(ny,nx,ndof,ndof)
        real :: Dh(ny,nx,ndof,ndof)
        real :: spgl(ny,nx), spg2l(ny,nx), Vh(ny,nx,6)
        real :: dtl(ny,nx)

        real :: detainv, detasinv
        
        integer :: iv, idof, jdof
!=======================================================================================================!
        detainv  = one / deta
        detasinv = one / deta**2
        
        if (.not. yper) then
        
!=======================================================================================================!
!.... use higher-order tangent on body node
!=======================================================================================================!

        do idof = 1, ndof
          do jdof = 1, ndof
            do iv = 1, nx
              mat(1,iv,idof,jdof,3) = alfa * dtl(1,iv) * gc1 * detainv * Bh(1,iv,idof,jdof)
              mat(1,iv,idof,jdof,4) = alfa * dtl(1,iv) * gc2 * detainv * Bh(1,iv,idof,jdof)
              mat(1,iv,idof,jdof,5) = alfa * dtl(1,iv) * gc3 * detainv * Bh(1,iv,idof,jdof)
              mat(1,iv,idof,jdof,1) = alfa * dtl(1,iv) * gc4 * detainv * Bh(1,iv,idof,jdof)
              mat(1,iv,idof,jdof,2) = alfa * dtl(1,iv) * gc5 * detainv * Bh(1,iv,idof,jdof)
            end do
          end do
        end do
        
!.... \hat{V}_{\eta\eta} term

        do iv = 1, nx
        
        mat(1,iv,2,2,3) = mat(1,iv,2,2,3) - alfa * dtl(1,iv) * dd1 * detasinv * Vh(1,iv,1)
        mat(1,iv,2,3,3) = mat(1,iv,2,3,3) - alfa * dtl(1,iv) * dd1 * detasinv * Vh(1,iv,5)
        mat(1,iv,3,2,3) = mat(1,iv,3,2,3) - alfa * dtl(1,iv) * dd1 * detasinv * Vh(1,iv,6)
        mat(1,iv,3,3,3) = mat(1,iv,3,3,3) - alfa * dtl(1,iv) * dd1 * detasinv * Vh(1,iv,2)
        mat(1,iv,4,4,3) = mat(1,iv,4,4,3) - alfa * dtl(1,iv) * dd1 * detasinv * Vh(1,iv,3)
        mat(1,iv,5,5,3) = mat(1,iv,5,5,3) - alfa * dtl(1,iv) * dd1 * detasinv * Vh(1,iv,4)

        mat(1,iv,2,2,4) = mat(1,iv,2,2,4) - alfa * dtl(1,iv) * dd2 * detasinv * Vh(1,iv,1)
        mat(1,iv,2,3,4) = mat(1,iv,2,3,4) - alfa * dtl(1,iv) * dd2 * detasinv * Vh(1,iv,5)
        mat(1,iv,3,2,4) = mat(1,iv,3,2,4) - alfa * dtl(1,iv) * dd2 * detasinv * Vh(1,iv,6)
        mat(1,iv,3,3,4) = mat(1,iv,3,3,4) - alfa * dtl(1,iv) * dd2 * detasinv * Vh(1,iv,2)
        mat(1,iv,4,4,4) = mat(1,iv,4,4,4) - alfa * dtl(1,iv) * dd2 * detasinv * Vh(1,iv,3)
        mat(1,iv,5,5,4) = mat(1,iv,5,5,4) - alfa * dtl(1,iv) * dd2 * detasinv * Vh(1,iv,4)

        mat(1,iv,2,2,5) = mat(1,iv,2,2,5) - alfa * dtl(1,iv) * dd3 * detasinv * Vh(1,iv,1)
        mat(1,iv,2,3,5) = mat(1,iv,2,3,5) - alfa * dtl(1,iv) * dd3 * detasinv * Vh(1,iv,5)
        mat(1,iv,3,2,5) = mat(1,iv,3,2,5) - alfa * dtl(1,iv) * dd3 * detasinv * Vh(1,iv,6)
        mat(1,iv,3,3,5) = mat(1,iv,3,3,5) - alfa * dtl(1,iv) * dd3 * detasinv * Vh(1,iv,2)
        mat(1,iv,4,4,5) = mat(1,iv,4,4,5) - alfa * dtl(1,iv) * dd3 * detasinv * Vh(1,iv,3)
        mat(1,iv,5,5,5) = mat(1,iv,5,5,5) - alfa * dtl(1,iv) * dd3 * detasinv * Vh(1,iv,4)

        mat(1,iv,2,2,1) = mat(1,iv,2,2,1) - alfa * dtl(1,iv) * dd4 * detasinv * Vh(1,iv,1)
        mat(1,iv,2,3,1) = mat(1,iv,2,3,1) - alfa * dtl(1,iv) * dd4 * detasinv * Vh(1,iv,5)
        mat(1,iv,3,2,1) = mat(1,iv,3,2,1) - alfa * dtl(1,iv) * dd4 * detasinv * Vh(1,iv,6)
        mat(1,iv,3,3,1) = mat(1,iv,3,3,1) - alfa * dtl(1,iv) * dd4 * detasinv * Vh(1,iv,2)
        mat(1,iv,4,4,1) = mat(1,iv,4,4,1) - alfa * dtl(1,iv) * dd4 * detasinv * Vh(1,iv,3)
        mat(1,iv,5,5,1) = mat(1,iv,5,5,1) - alfa * dtl(1,iv) * dd4 * detasinv * Vh(1,iv,4)

        mat(1,iv,2,2,2) = mat(1,iv,2,2,2) - alfa * dtl(1,iv) * dd5 * detasinv * Vh(1,iv,1)
        mat(1,iv,2,3,2) = mat(1,iv,2,3,2) - alfa * dtl(1,iv) * dd5 * detasinv * Vh(1,iv,5)
        mat(1,iv,3,2,2) = mat(1,iv,3,2,2) - alfa * dtl(1,iv) * dd5 * detasinv * Vh(1,iv,6)
        mat(1,iv,3,3,2) = mat(1,iv,3,3,2) - alfa * dtl(1,iv) * dd5 * detasinv * Vh(1,iv,2)
        mat(1,iv,4,4,2) = mat(1,iv,4,4,2) - alfa * dtl(1,iv) * dd5 * detasinv * Vh(1,iv,3)
        mat(1,iv,5,5,2) = mat(1,iv,5,5,2) - alfa * dtl(1,iv) * dd5 * detasinv * Vh(1,iv,4)

!.... I term

        mat(1,iv,1,1,3) = mat(1,iv,1,1,3) + one
        mat(1,iv,2,2,3) = mat(1,iv,2,2,3) + one
        mat(1,iv,3,3,3) = mat(1,iv,3,3,3) + one
        mat(1,iv,4,4,3) = mat(1,iv,4,4,3) + one
        mat(1,iv,5,5,3) = mat(1,iv,5,5,3) + one

        end do
        
!.... sponge

        if (ispg .eq. 1) then
        
          do iv = 1, nx
            mat(1,iv,1,1,3) = mat(1,iv,1,1,3) + alfa * dtl(1,iv) * spgl(1,iv)
            mat(1,iv,2,2,3) = mat(1,iv,2,2,3) + alfa * dtl(1,iv) * spgl(1,iv)
            mat(1,iv,3,3,3) = mat(1,iv,3,3,3) + alfa * dtl(1,iv) * spgl(1,iv)
            mat(1,iv,4,4,3) = mat(1,iv,4,4,3) + alfa * dtl(1,iv) * spgl(1,iv)
            mat(1,iv,5,5,3) = mat(1,iv,5,5,3) + alfa * dtl(1,iv) * spgl(1,iv)
          end do

        else if (ispg .eq. 2) then
        
          do iv = 1, nx
            mat(1,iv,1,1,3) = mat(1,iv,1,1,3) + alfa * dtl(1,iv) * (spgl(1,iv) + spg2l(1,iv))
            mat(1,iv,2,2,3) = mat(1,iv,2,2,3) + alfa * dtl(1,iv) * (spgl(1,iv) + spg2l(1,iv))
            mat(1,iv,3,3,3) = mat(1,iv,3,3,3) + alfa * dtl(1,iv) * (spgl(1,iv) + spg2l(1,iv))
            mat(1,iv,4,4,3) = mat(1,iv,4,4,3) + alfa * dtl(1,iv) * (spgl(1,iv) + spg2l(1,iv))
            mat(1,iv,5,5,3) = mat(1,iv,5,5,3) + alfa * dtl(1,iv) * (spgl(1,iv) + spg2l(1,iv))
          end do
        
        end if

!=======================================================================================================!
!.... use higher-order tangent on first node off the body
!=======================================================================================================!

        do idof = 1, ndof
          do jdof = 1, ndof
            do iv = 1, nx
              mat(2,iv,idof,jdof,2) = alfa * dtl(2,iv) * gb1 * detainv * Bh(2,iv,idof,jdof)
              mat(2,iv,idof,jdof,3) = alfa * dtl(2,iv) * gb2 * detainv * Bh(2,iv,idof,jdof)
              mat(2,iv,idof,jdof,4) = alfa * dtl(2,iv) * gb3 * detainv * Bh(2,iv,idof,jdof)
              mat(2,iv,idof,jdof,5) = alfa * dtl(2,iv) * gb4 * detainv * Bh(2,iv,idof,jdof)
              mat(2,iv,idof,jdof,1) = alfa * dtl(2,iv) * gb5 * detainv * Bh(2,iv,idof,jdof)
            end do
          end do
        end do

!.... \hat{V}_{\eta\eta} term

        do iv = 1, nx

        mat(2,iv,2,2,2) = mat(2,iv,2,2,2) - alfa * dtl(2,iv) * db1 * detasinv * Vh(2,iv,1)
        mat(2,iv,2,3,2) = mat(2,iv,2,3,2) - alfa * dtl(2,iv) * db1 * detasinv * Vh(2,iv,5)
        mat(2,iv,3,2,2) = mat(2,iv,3,2,2) - alfa * dtl(2,iv) * db1 * detasinv * Vh(2,iv,6)
        mat(2,iv,3,3,2) = mat(2,iv,3,3,2) - alfa * dtl(2,iv) * db1 * detasinv * Vh(2,iv,2)
        mat(2,iv,4,4,2) = mat(2,iv,4,4,2) - alfa * dtl(2,iv) * db1 * detasinv * Vh(2,iv,3)
        mat(2,iv,5,5,2) = mat(2,iv,5,5,2) - alfa * dtl(2,iv) * db1 * detasinv * Vh(2,iv,4)

        mat(2,iv,2,2,3) = mat(2,iv,2,2,3) - alfa * dtl(2,iv) * db2 * detasinv * Vh(2,iv,1)
        mat(2,iv,2,3,3) = mat(2,iv,2,3,3) - alfa * dtl(2,iv) * db2 * detasinv * Vh(2,iv,5)
        mat(2,iv,3,2,3) = mat(2,iv,3,2,3) - alfa * dtl(2,iv) * db2 * detasinv * Vh(2,iv,6)
        mat(2,iv,3,3,3) = mat(2,iv,3,3,3) - alfa * dtl(2,iv) * db2 * detasinv * Vh(2,iv,2)
        mat(2,iv,4,4,3) = mat(2,iv,4,4,3) - alfa * dtl(2,iv) * db2 * detasinv * Vh(2,iv,3)
        mat(2,iv,5,5,3) = mat(2,iv,5,5,3) - alfa * dtl(2,iv) * db2 * detasinv * Vh(2,iv,4)

        mat(2,iv,2,2,4) = mat(2,iv,2,2,4) - alfa * dtl(2,iv) * db3 * detasinv * Vh(2,iv,1)
        mat(2,iv,2,3,4) = mat(2,iv,2,3,4) - alfa * dtl(2,iv) * db3 * detasinv * Vh(2,iv,5)
        mat(2,iv,3,2,4) = mat(2,iv,3,2,4) - alfa * dtl(2,iv) * db3 * detasinv * Vh(2,iv,6)
        mat(2,iv,3,3,4) = mat(2,iv,3,3,4) - alfa * dtl(2,iv) * db3 * detasinv * Vh(2,iv,2)
        mat(2,iv,4,4,4) = mat(2,iv,4,4,4) - alfa * dtl(2,iv) * db3 * detasinv * Vh(2,iv,3)
        mat(2,iv,5,5,4) = mat(2,iv,5,5,4) - alfa * dtl(2,iv) * db3 * detasinv * Vh(2,iv,4)

        mat(2,iv,2,2,5) = mat(2,iv,2,2,5) - alfa * dtl(2,iv) * db4 * detasinv * Vh(2,iv,1)
        mat(2,iv,2,3,5) = mat(2,iv,2,3,5) - alfa * dtl(2,iv) * db4 * detasinv * Vh(2,iv,5)
        mat(2,iv,3,2,5) = mat(2,iv,3,2,5) - alfa * dtl(2,iv) * db4 * detasinv * Vh(2,iv,6)
        mat(2,iv,3,3,5) = mat(2,iv,3,3,5) - alfa * dtl(2,iv) * db4 * detasinv * Vh(2,iv,2)
        mat(2,iv,4,4,5) = mat(2,iv,4,4,5) - alfa * dtl(2,iv) * db4 * detasinv * Vh(2,iv,3)
        mat(2,iv,5,5,5) = mat(2,iv,5,5,5) - alfa * dtl(2,iv) * db4 * detasinv * Vh(2,iv,4)

        mat(2,iv,2,2,1) = mat(2,iv,2,2,1) - alfa * dtl(2,iv) * db5 * detasinv * Vh(2,iv,1)
        mat(2,iv,2,3,1) = mat(2,iv,2,3,1) - alfa * dtl(2,iv) * db5 * detasinv * Vh(2,iv,5)
        mat(2,iv,3,2,1) = mat(2,iv,3,2,1) - alfa * dtl(2,iv) * db5 * detasinv * Vh(2,iv,6)
        mat(2,iv,3,3,1) = mat(2,iv,3,3,1) - alfa * dtl(2,iv) * db5 * detasinv * Vh(2,iv,2)
        mat(2,iv,4,4,1) = mat(2,iv,4,4,1) - alfa * dtl(2,iv) * db5 * detasinv * Vh(2,iv,3)
        mat(2,iv,5,5,1) = mat(2,iv,5,5,1) - alfa * dtl(2,iv) * db5 * detasinv * Vh(2,iv,4)

!.... I term

        mat(2,iv,1,1,3) = mat(2,iv,1,1,3) + one
        mat(2,iv,2,2,3) = mat(2,iv,2,2,3) + one
        mat(2,iv,3,3,3) = mat(2,iv,3,3,3) + one
        mat(2,iv,4,4,3) = mat(2,iv,4,4,3) + one
        mat(2,iv,5,5,3) = mat(2,iv,5,5,3) + one

        end do
        
!.... sponge

        if (ispg .eq. 1) then
        
          do iv = 1, nx
            mat(2,iv,1,1,3) = mat(2,iv,1,1,3) + alfa * dtl(2,iv) * spgl(2,iv)
            mat(2,iv,2,2,3) = mat(2,iv,2,2,3) + alfa * dtl(2,iv) * spgl(2,iv)
            mat(2,iv,3,3,3) = mat(2,iv,3,3,3) + alfa * dtl(2,iv) * spgl(2,iv)
            mat(2,iv,4,4,3) = mat(2,iv,4,4,3) + alfa * dtl(2,iv) * spgl(2,iv)
            mat(2,iv,5,5,3) = mat(2,iv,5,5,3) + alfa * dtl(2,iv) * spgl(2,iv)
          end do
                
        else if (ispg .ge. 2) then
        
          do iv = 1, nx
            mat(2,iv,1,1,3) = mat(2,iv,1,1,3) + alfa * dtl(2,iv) * (spgl(2,iv) + spg2l(2,iv))
            mat(2,iv,2,2,3) = mat(2,iv,2,2,3) + alfa * dtl(2,iv) * (spgl(2,iv) + spg2l(2,iv))
            mat(2,iv,3,3,3) = mat(2,iv,3,3,3) + alfa * dtl(2,iv) * (spgl(2,iv) + spg2l(2,iv))
            mat(2,iv,4,4,3) = mat(2,iv,4,4,3) + alfa * dtl(2,iv) * (spgl(2,iv) + spg2l(2,iv))
            mat(2,iv,5,5,3) = mat(2,iv,5,5,3) + alfa * dtl(2,iv) * (spgl(2,iv) + spg2l(2,iv))
          end do
                
        end if

!.... implicit damping term

        if (eps_e .ne. zero) then

        do iv = 1, nx
        mat(2,iv,1,1,2) = mat(2,iv,1,1,2) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * (-one)
        mat(2,iv,2,2,2) = mat(2,iv,2,2,2) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * (-one)
        mat(2,iv,3,3,2) = mat(2,iv,3,3,2) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * (-one)
        mat(2,iv,4,4,2) = mat(2,iv,4,4,2) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * (-one)
        mat(2,iv,5,5,2) = mat(2,iv,5,5,2) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * (-one)

        mat(2,iv,1,1,3) = mat(2,iv,1,1,3) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * two
        mat(2,iv,2,2,3) = mat(2,iv,2,2,3) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * two
        mat(2,iv,3,3,3) = mat(2,iv,3,3,3) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * two
        mat(2,iv,4,4,3) = mat(2,iv,4,4,3) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * two
        mat(2,iv,5,5,3) = mat(2,iv,5,5,3) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * two

        mat(2,iv,1,1,4) = mat(2,iv,1,1,4) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * (-one)
        mat(2,iv,2,2,4) = mat(2,iv,2,2,4) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * (-one)
        mat(2,iv,3,3,4) = mat(2,iv,3,3,4) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * (-one)
        mat(2,iv,4,4,4) = mat(2,iv,4,4,4) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * (-one)
        mat(2,iv,5,5,4) = mat(2,iv,5,5,4) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * (-one)
        end do
                
        end if

!=======================================================================================================!
!.... use higher-order tangent on first node off the far-field boundary
!=======================================================================================================!

        do idof = 1, ndof
          do jdof = 1, ndof
            do iv = 1, nx
              mat(ny-1,iv,idof,jdof,5) = -alfa * dtl(ny-1,iv) * gb5 * detainv * Bh(ny-1,iv,idof,jdof)
              mat(ny-1,iv,idof,jdof,1) = -alfa * dtl(ny-1,iv) * gb4 * detainv * Bh(ny-1,iv,idof,jdof)
              mat(ny-1,iv,idof,jdof,2) = -alfa * dtl(ny-1,iv) * gb3 * detainv * Bh(ny-1,iv,idof,jdof)
              mat(ny-1,iv,idof,jdof,3) = -alfa * dtl(ny-1,iv) * gb2 * detainv * Bh(ny-1,iv,idof,jdof)
              mat(ny-1,iv,idof,jdof,4) = -alfa * dtl(ny-1,iv) * gb1 * detainv * Bh(ny-1,iv,idof,jdof)
            end do
          end do
        end do
        
!.... \hat{V}_{\eta\eta} term

        do iv = 1, nx

        mat(ny-1,iv,2,2,5) = mat(ny-1,iv,2,2,5) - alfa * dtl(ny-1,iv) * db5 * detasinv * Vh(ny-1,iv,1)
        mat(ny-1,iv,2,3,5) = mat(ny-1,iv,2,3,5) - alfa * dtl(ny-1,iv) * db5 * detasinv * Vh(ny-1,iv,5)
        mat(ny-1,iv,3,2,5) = mat(ny-1,iv,3,2,5) - alfa * dtl(ny-1,iv) * db5 * detasinv * Vh(ny-1,iv,6)
        mat(ny-1,iv,3,3,5) = mat(ny-1,iv,3,3,5) - alfa * dtl(ny-1,iv) * db5 * detasinv * Vh(ny-1,iv,2)
        mat(ny-1,iv,4,4,5) = mat(ny-1,iv,4,4,5) - alfa * dtl(ny-1,iv) * db5 * detasinv * Vh(ny-1,iv,3)
        mat(ny-1,iv,5,5,5) = mat(ny-1,iv,5,5,5) - alfa * dtl(ny-1,iv) * db5 * detasinv * Vh(ny-1,iv,4)

        mat(ny-1,iv,2,2,1) = mat(ny-1,iv,2,2,1) - alfa * dtl(ny-1,iv) * db4 * detasinv * Vh(ny-1,iv,1)
        mat(ny-1,iv,2,3,1) = mat(ny-1,iv,2,3,1) - alfa * dtl(ny-1,iv) * db4 * detasinv * Vh(ny-1,iv,5)
        mat(ny-1,iv,3,2,1) = mat(ny-1,iv,3,2,1) - alfa * dtl(ny-1,iv) * db4 * detasinv * Vh(ny-1,iv,6)
        mat(ny-1,iv,3,3,1) = mat(ny-1,iv,3,3,1) - alfa * dtl(ny-1,iv) * db4 * detasinv * Vh(ny-1,iv,2)
        mat(ny-1,iv,4,4,1) = mat(ny-1,iv,4,4,1) - alfa * dtl(ny-1,iv) * db4 * detasinv * Vh(ny-1,iv,3)
        mat(ny-1,iv,5,5,1) = mat(ny-1,iv,5,5,1) - alfa * dtl(ny-1,iv) * db4 * detasinv * Vh(ny-1,iv,4)

        mat(ny-1,iv,2,2,2) = mat(ny-1,iv,2,2,2) - alfa * dtl(ny-1,iv) * db3 * detasinv * Vh(ny-1,iv,1)
        mat(ny-1,iv,2,3,2) = mat(ny-1,iv,2,3,2) - alfa * dtl(ny-1,iv) * db3 * detasinv * Vh(ny-1,iv,5)
        mat(ny-1,iv,3,2,2) = mat(ny-1,iv,3,2,2) - alfa * dtl(ny-1,iv) * db3 * detasinv * Vh(ny-1,iv,6)
        mat(ny-1,iv,3,3,2) = mat(ny-1,iv,3,3,2) - alfa * dtl(ny-1,iv) * db3 * detasinv * Vh(ny-1,iv,2)
        mat(ny-1,iv,4,4,2) = mat(ny-1,iv,4,4,2) - alfa * dtl(ny-1,iv) * db3 * detasinv * Vh(ny-1,iv,3)
        mat(ny-1,iv,5,5,2) = mat(ny-1,iv,5,5,2) - alfa * dtl(ny-1,iv) * db3 * detasinv * Vh(ny-1,iv,4)

        mat(ny-1,iv,2,2,3) = mat(ny-1,iv,2,2,3) - alfa * dtl(ny-1,iv) * db2 * detasinv * Vh(ny-1,iv,1)
        mat(ny-1,iv,2,3,3) = mat(ny-1,iv,2,3,3) - alfa * dtl(ny-1,iv) * db2 * detasinv * Vh(ny-1,iv,5)
        mat(ny-1,iv,3,2,3) = mat(ny-1,iv,3,2,3) - alfa * dtl(ny-1,iv) * db2 * detasinv * Vh(ny-1,iv,6)
        mat(ny-1,iv,3,3,3) = mat(ny-1,iv,3,3,3) - alfa * dtl(ny-1,iv) * db2 * detasinv * Vh(ny-1,iv,2)
        mat(ny-1,iv,4,4,3) = mat(ny-1,iv,4,4,3) - alfa * dtl(ny-1,iv) * db2 * detasinv * Vh(ny-1,iv,3)
        mat(ny-1,iv,5,5,3) = mat(ny-1,iv,5,5,3) - alfa * dtl(ny-1,iv) * db2 * detasinv * Vh(ny-1,iv,4)

        mat(ny-1,iv,2,2,4) = mat(ny-1,iv,2,2,4) - alfa * dtl(ny-1,iv) * db1 * detasinv * Vh(ny-1,iv,1)
        mat(ny-1,iv,2,3,4) = mat(ny-1,iv,2,3,4) - alfa * dtl(ny-1,iv) * db1 * detasinv * Vh(ny-1,iv,5)
        mat(ny-1,iv,3,2,4) = mat(ny-1,iv,3,2,4) - alfa * dtl(ny-1,iv) * db1 * detasinv * Vh(ny-1,iv,6)
        mat(ny-1,iv,3,3,4) = mat(ny-1,iv,3,3,4) - alfa * dtl(ny-1,iv) * db1 * detasinv * Vh(ny-1,iv,2)
        mat(ny-1,iv,4,4,4) = mat(ny-1,iv,4,4,4) - alfa * dtl(ny-1,iv) * db1 * detasinv * Vh(ny-1,iv,3)
        mat(ny-1,iv,5,5,4) = mat(ny-1,iv,5,5,4) - alfa * dtl(ny-1,iv) * db1 * detasinv * Vh(ny-1,iv,4)

!.... I term

        mat(ny-1,iv,1,1,3) = mat(ny-1,iv,1,1,3) + one
        mat(ny-1,iv,2,2,3) = mat(ny-1,iv,2,2,3) + one
        mat(ny-1,iv,3,3,3) = mat(ny-1,iv,3,3,3) + one
        mat(ny-1,iv,4,4,3) = mat(ny-1,iv,4,4,3) + one
        mat(ny-1,iv,5,5,3) = mat(ny-1,iv,5,5,3) + one

        end do
        
!.... sponge

        if (ispg .eq. 1) then

        do iv = 1, nx
          mat(ny-1,iv,1,1,3) = mat(ny-1,iv,1,1,3) + alfa * dtl(ny-1,iv) * spgl(ny-1,iv)
          mat(ny-1,iv,2,2,3) = mat(ny-1,iv,2,2,3) + alfa * dtl(ny-1,iv) * spgl(ny-1,iv)
          mat(ny-1,iv,3,3,3) = mat(ny-1,iv,3,3,3) + alfa * dtl(ny-1,iv) * spgl(ny-1,iv)
          mat(ny-1,iv,4,4,3) = mat(ny-1,iv,4,4,3) + alfa * dtl(ny-1,iv) * spgl(ny-1,iv)
          mat(ny-1,iv,5,5,3) = mat(ny-1,iv,5,5,3) + alfa * dtl(ny-1,iv) * spgl(ny-1,iv)
        end do

        else if (ispg .eq. 2) then

        do iv = 1, nx
          mat(ny-1,iv,1,1,3) = mat(ny-1,iv,1,1,3) + alfa * dtl(ny-1,iv) * (spgl(ny-1,iv) + spg2l(ny-1,iv))
          mat(ny-1,iv,2,2,3) = mat(ny-1,iv,2,2,3) + alfa * dtl(ny-1,iv) * (spgl(ny-1,iv) + spg2l(ny-1,iv))
          mat(ny-1,iv,3,3,3) = mat(ny-1,iv,3,3,3) + alfa * dtl(ny-1,iv) * (spgl(ny-1,iv) + spg2l(ny-1,iv))
          mat(ny-1,iv,4,4,3) = mat(ny-1,iv,4,4,3) + alfa * dtl(ny-1,iv) * (spgl(ny-1,iv) + spg2l(ny-1,iv))
          mat(ny-1,iv,5,5,3) = mat(ny-1,iv,5,5,3) + alfa * dtl(ny-1,iv) * (spgl(ny-1,iv) + spg2l(ny-1,iv))
        end do

        end if

!.... implicit damping term

        if (eps_e .ne. zero) then

        do iv = 1, nx
        mat(ny-1,iv,1,1,2) = mat(ny-1,iv,1,1,2) + alfa * dtl(ny-1,iv) * eps_e * buff(ny-1,iv) * (-one)
        mat(ny-1,iv,2,2,2) = mat(ny-1,iv,2,2,2) + alfa * dtl(ny-1,iv) * eps_e * buff(ny-1,iv) * (-one)
        mat(ny-1,iv,3,3,2) = mat(ny-1,iv,3,3,2) + alfa * dtl(ny-1,iv) * eps_e * buff(ny-1,iv) * (-one)
        mat(ny-1,iv,4,4,2) = mat(ny-1,iv,4,4,2) + alfa * dtl(ny-1,iv) * eps_e * buff(ny-1,iv) * (-one)
        mat(ny-1,iv,5,5,2) = mat(ny-1,iv,5,5,2) + alfa * dtl(ny-1,iv) * eps_e * buff(ny-1,iv) * (-one)

        mat(ny-1,iv,1,1,3) = mat(ny-1,iv,1,1,3) + alfa * dtl(ny-1,iv) * eps_e * buff(ny-1,iv) * two
        mat(ny-1,iv,2,2,3) = mat(ny-1,iv,2,2,3) + alfa * dtl(ny-1,iv) * eps_e * buff(ny-1,iv) * two
        mat(ny-1,iv,3,3,3) = mat(ny-1,iv,3,3,3) + alfa * dtl(ny-1,iv) * eps_e * buff(ny-1,iv) * two
        mat(ny-1,iv,4,4,3) = mat(ny-1,iv,4,4,3) + alfa * dtl(ny-1,iv) * eps_e * buff(ny-1,iv) * two
        mat(ny-1,iv,5,5,3) = mat(ny-1,iv,5,5,3) + alfa * dtl(ny-1,iv) * eps_e * buff(ny-1,iv) * two

        mat(ny-1,iv,1,1,4) = mat(ny-1,iv,1,1,4) + alfa * dtl(ny-1,iv) * eps_e * buff(ny-1,iv) * (-one)
        mat(ny-1,iv,2,2,4) = mat(ny-1,iv,2,2,4) + alfa * dtl(ny-1,iv) * eps_e * buff(ny-1,iv) * (-one)
        mat(ny-1,iv,3,3,4) = mat(ny-1,iv,3,3,4) + alfa * dtl(ny-1,iv) * eps_e * buff(ny-1,iv) * (-one)
        mat(ny-1,iv,4,4,4) = mat(ny-1,iv,4,4,4) + alfa * dtl(ny-1,iv) * eps_e * buff(ny-1,iv) * (-one)
        mat(ny-1,iv,5,5,4) = mat(ny-1,iv,5,5,4) + alfa * dtl(ny-1,iv) * eps_e * buff(ny-1,iv) * (-one)
        end do
                
        end if

!=======================================================================================================!
!.... use higher-order tangent on the far-field boundary
!=======================================================================================================!

        do idof = 1, ndof
          do jdof = 1, ndof
            do iv = 1, nx
              mat(ny,iv,idof,jdof,4) = -alfa * dtl(ny,iv) * gc5 * detainv * Bh(ny,iv,idof,jdof)
              mat(ny,iv,idof,jdof,5) = -alfa * dtl(ny,iv) * gc4 * detainv * Bh(ny,iv,idof,jdof)
              mat(ny,iv,idof,jdof,1) = -alfa * dtl(ny,iv) * gc3 * detainv * Bh(ny,iv,idof,jdof)
              mat(ny,iv,idof,jdof,2) = -alfa * dtl(ny,iv) * gc2 * detainv * Bh(ny,iv,idof,jdof)
              mat(ny,iv,idof,jdof,3) = -alfa * dtl(ny,iv) * gc1 * detainv * Bh(ny,iv,idof,jdof)
            end do
          end do
        end do
        
!.... \hat{V}_{\eta\eta} term

        do iv = 1, nx

        mat(ny,iv,2,2,4) = mat(ny,iv,2,2,4) - alfa * dtl(ny,iv) * dd5 * detasinv * Vh(ny,iv,1)
        mat(ny,iv,2,3,4) = mat(ny,iv,2,3,4) - alfa * dtl(ny,iv) * dd5 * detasinv * Vh(ny,iv,5)
        mat(ny,iv,3,2,4) = mat(ny,iv,3,2,4) - alfa * dtl(ny,iv) * dd5 * detasinv * Vh(ny,iv,6)
        mat(ny,iv,3,3,4) = mat(ny,iv,3,3,4) - alfa * dtl(ny,iv) * dd5 * detasinv * Vh(ny,iv,2)
        mat(ny,iv,4,4,4) = mat(ny,iv,4,4,4) - alfa * dtl(ny,iv) * dd5 * detasinv * Vh(ny,iv,3)
        mat(ny,iv,5,5,4) = mat(ny,iv,5,5,4) - alfa * dtl(ny,iv) * dd5 * detasinv * Vh(ny,iv,4)

        mat(ny,iv,2,2,5) = mat(ny,iv,2,2,5) - alfa * dtl(ny,iv) * dd4 * detasinv * Vh(ny,iv,1)
        mat(ny,iv,2,3,5) = mat(ny,iv,2,3,5) - alfa * dtl(ny,iv) * dd4 * detasinv * Vh(ny,iv,5)
        mat(ny,iv,3,2,5) = mat(ny,iv,3,2,5) - alfa * dtl(ny,iv) * dd4 * detasinv * Vh(ny,iv,6)
        mat(ny,iv,3,3,5) = mat(ny,iv,3,3,5) - alfa * dtl(ny,iv) * dd4 * detasinv * Vh(ny,iv,2)
        mat(ny,iv,4,4,5) = mat(ny,iv,4,4,5) - alfa * dtl(ny,iv) * dd4 * detasinv * Vh(ny,iv,3)
        mat(ny,iv,5,5,5) = mat(ny,iv,5,5,5) - alfa * dtl(ny,iv) * dd4 * detasinv * Vh(ny,iv,4)

        mat(ny,iv,2,2,1) = mat(ny,iv,2,2,1) - alfa * dtl(ny,iv) * dd3 * detasinv * Vh(ny,iv,1)
        mat(ny,iv,2,3,1) = mat(ny,iv,2,3,1) - alfa * dtl(ny,iv) * dd3 * detasinv * Vh(ny,iv,5)
        mat(ny,iv,3,2,1) = mat(ny,iv,3,2,1) - alfa * dtl(ny,iv) * dd3 * detasinv * Vh(ny,iv,6)
        mat(ny,iv,3,3,1) = mat(ny,iv,3,3,1) - alfa * dtl(ny,iv) * dd3 * detasinv * Vh(ny,iv,2)
        mat(ny,iv,4,4,1) = mat(ny,iv,4,4,1) - alfa * dtl(ny,iv) * dd3 * detasinv * Vh(ny,iv,3)
        mat(ny,iv,5,5,1) = mat(ny,iv,5,5,1) - alfa * dtl(ny,iv) * dd3 * detasinv * Vh(ny,iv,4)

        mat(ny,iv,2,2,2) = mat(ny,iv,2,2,2) - alfa * dtl(ny,iv) * dd2 * detasinv * Vh(ny,iv,1)
        mat(ny,iv,2,3,2) = mat(ny,iv,2,3,2) - alfa * dtl(ny,iv) * dd2 * detasinv * Vh(ny,iv,5)
        mat(ny,iv,3,2,2) = mat(ny,iv,3,2,2) - alfa * dtl(ny,iv) * dd2 * detasinv * Vh(ny,iv,6)
        mat(ny,iv,3,3,2) = mat(ny,iv,3,3,2) - alfa * dtl(ny,iv) * dd2 * detasinv * Vh(ny,iv,2)
        mat(ny,iv,4,4,2) = mat(ny,iv,4,4,2) - alfa * dtl(ny,iv) * dd2 * detasinv * Vh(ny,iv,3)
        mat(ny,iv,5,5,2) = mat(ny,iv,5,5,2) - alfa * dtl(ny,iv) * dd2 * detasinv * Vh(ny,iv,4)

        mat(ny,iv,2,2,3) = mat(ny,iv,2,2,3) - alfa * dtl(ny,iv) * dd1 * detasinv * Vh(ny,iv,1)
        mat(ny,iv,2,3,3) = mat(ny,iv,2,3,3) - alfa * dtl(ny,iv) * dd1 * detasinv * Vh(ny,iv,5)
        mat(ny,iv,3,2,3) = mat(ny,iv,3,2,3) - alfa * dtl(ny,iv) * dd1 * detasinv * Vh(ny,iv,6)
        mat(ny,iv,3,3,3) = mat(ny,iv,3,3,3) - alfa * dtl(ny,iv) * dd1 * detasinv * Vh(ny,iv,2)
        mat(ny,iv,4,4,3) = mat(ny,iv,4,4,3) - alfa * dtl(ny,iv) * dd1 * detasinv * Vh(ny,iv,3)
        mat(ny,iv,5,5,3) = mat(ny,iv,5,5,3) - alfa * dtl(ny,iv) * dd1 * detasinv * Vh(ny,iv,4)

!.... I term

        mat(ny,iv,1,1,3) = mat(ny,iv,1,1,3) + one
        mat(ny,iv,2,2,3) = mat(ny,iv,2,2,3) + one
        mat(ny,iv,3,3,3) = mat(ny,iv,3,3,3) + one
        mat(ny,iv,4,4,3) = mat(ny,iv,4,4,3) + one
        mat(ny,iv,5,5,3) = mat(ny,iv,5,5,3) + one

        end do
        
!.... sponge

        if (ispg .eq. 1) then
        
          do iv = 1, nx
            mat(ny,iv,1,1,3) = mat(ny,iv,1,1,3) + alfa * dtl(ny,iv) * spgl(ny,iv)
            mat(ny,iv,2,2,3) = mat(ny,iv,2,2,3) + alfa * dtl(ny,iv) * spgl(ny,iv)
            mat(ny,iv,3,3,3) = mat(ny,iv,3,3,3) + alfa * dtl(ny,iv) * spgl(ny,iv)
            mat(ny,iv,4,4,3) = mat(ny,iv,4,4,3) + alfa * dtl(ny,iv) * spgl(ny,iv)
            mat(ny,iv,5,5,3) = mat(ny,iv,5,5,3) + alfa * dtl(ny,iv) * spgl(ny,iv)
          end do
      
        else if (ispg .eq. 2) then
        
          do iv = 1, nx
            mat(ny,iv,1,1,3) = mat(ny,iv,1,1,3) + alfa * dtl(ny,iv) * (spgl(ny,iv) + spg2l(ny,iv))
            mat(ny,iv,2,2,3) = mat(ny,iv,2,2,3) + alfa * dtl(ny,iv) * (spgl(ny,iv) + spg2l(ny,iv))
            mat(ny,iv,3,3,3) = mat(ny,iv,3,3,3) + alfa * dtl(ny,iv) * (spgl(ny,iv) + spg2l(ny,iv))
            mat(ny,iv,4,4,3) = mat(ny,iv,4,4,3) + alfa * dtl(ny,iv) * (spgl(ny,iv) + spg2l(ny,iv))
            mat(ny,iv,5,5,3) = mat(ny,iv,5,5,3) + alfa * dtl(ny,iv) * (spgl(ny,iv) + spg2l(ny,iv))
          end do
                
        end if

        end if          ! yper

        return
        end
