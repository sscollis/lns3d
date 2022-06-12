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
        use buffer
        implicit none
        
        real :: mat(3,ndof,ndof,nx,ny), Bh(ndof,ndof,nx,ny)
        real :: Dh(ndof,ndof,nx,ny), bc(14,ndof,ndof,nx)
        real :: spgl(ny,nx), spg2l(ny,nx), Vh(6,nx,ny)
        real :: dtl(nx,ny)

        real :: detainv, detasinv
        
        integer :: iv, idof, jdof
!=======================================================================================================!
        detainv  = one / deta
        detasinv = one / deta**2
        
!.... initialize the bc array

        bc = zero

        if (.not. yper) then
        
        !$omp parallel do private(idof,jdof)
        do iv = 1, nx

!=======================================================================================================!
!.... use higher-order tangent on body node
!=======================================================================================================!
        do idof = 1, ndof
          do jdof = 1, ndof
            mat(1,idof,jdof,iv,1) = zero
            mat(2,idof,jdof,iv,1) = alfa * dtl(iv,1) * gc1 * detainv * Bh(idof,jdof,iv,1)
            mat(3,idof,jdof,iv,1) = alfa * dtl(iv,1) * gc2 * detainv * Bh(idof,jdof,iv,1)
            bc(1,idof,jdof,iv)    = alfa * dtl(iv,1) * gc3 * detainv * Bh(idof,jdof,iv,1)
            bc(2,idof,jdof,iv)    = alfa * dtl(iv,1) * gc4 * detainv * Bh(idof,jdof,iv,1)
            bc(3,idof,jdof,iv)    = alfa * dtl(iv,1) * gc5 * detainv * Bh(idof,jdof,iv,1)
          end do
        end do
        
!.... \hat{V}_{\eta\eta} term

        mat(2,2,2,iv,1) = mat(2,2,2,iv,1) - alfa * dtl(iv,1) * dd1 * detasinv * Vh(1,iv,1)
        mat(2,2,3,iv,1) = mat(2,2,3,iv,1) - alfa * dtl(iv,1) * dd1 * detasinv * Vh(5,iv,1)
        mat(2,3,2,iv,1) = mat(2,3,2,iv,1) - alfa * dtl(iv,1) * dd1 * detasinv * Vh(6,iv,1)
        mat(2,3,3,iv,1) = mat(2,3,3,iv,1) - alfa * dtl(iv,1) * dd1 * detasinv * Vh(2,iv,1)
        mat(2,4,4,iv,1) = mat(2,4,4,iv,1) - alfa * dtl(iv,1) * dd1 * detasinv * Vh(3,iv,1)
        mat(2,5,5,iv,1) = mat(2,5,5,iv,1) - alfa * dtl(iv,1) * dd1 * detasinv * Vh(4,iv,1)

        mat(3,2,2,iv,1) = mat(3,2,2,iv,1) - alfa * dtl(iv,1) * dd2 * detasinv * Vh(1,iv,1)
        mat(3,2,3,iv,1) = mat(3,2,3,iv,1) - alfa * dtl(iv,1) * dd2 * detasinv * Vh(5,iv,1)
        mat(3,3,2,iv,1) = mat(3,3,2,iv,1) - alfa * dtl(iv,1) * dd2 * detasinv * Vh(6,iv,1)
        mat(3,3,3,iv,1) = mat(3,3,3,iv,1) - alfa * dtl(iv,1) * dd2 * detasinv * Vh(2,iv,1)
        mat(3,4,4,iv,1) = mat(3,4,4,iv,1) - alfa * dtl(iv,1) * dd2 * detasinv * Vh(3,iv,1)
        mat(3,5,5,iv,1) = mat(3,5,5,iv,1) - alfa * dtl(iv,1) * dd2 * detasinv * Vh(4,iv,1)

        bc(1,2,2,iv) = bc(1,2,2,iv) - alfa * dtl(iv,1) * dd3 * detasinv * Vh(1,iv,1)
        bc(1,2,3,iv) = bc(1,2,3,iv) - alfa * dtl(iv,1) * dd3 * detasinv * Vh(5,iv,1)
        bc(1,3,2,iv) = bc(1,3,2,iv) - alfa * dtl(iv,1) * dd3 * detasinv * Vh(6,iv,1)
        bc(1,3,3,iv) = bc(1,3,3,iv) - alfa * dtl(iv,1) * dd3 * detasinv * Vh(2,iv,1)
        bc(1,4,4,iv) = bc(1,4,4,iv) - alfa * dtl(iv,1) * dd3 * detasinv * Vh(3,iv,1)
        bc(1,5,5,iv) = bc(1,5,5,iv) - alfa * dtl(iv,1) * dd3 * detasinv * Vh(4,iv,1)

        bc(2,2,2,iv) = bc(2,2,2,iv) - alfa * dtl(iv,1) * dd4 * detasinv * Vh(1,iv,1)
        bc(2,2,3,iv) = bc(2,2,3,iv) - alfa * dtl(iv,1) * dd4 * detasinv * Vh(5,iv,1)
        bc(2,3,2,iv) = bc(2,3,2,iv) - alfa * dtl(iv,1) * dd4 * detasinv * Vh(6,iv,1)
        bc(2,3,3,iv) = bc(2,3,3,iv) - alfa * dtl(iv,1) * dd4 * detasinv * Vh(2,iv,1)
        bc(2,4,4,iv) = bc(2,4,4,iv) - alfa * dtl(iv,1) * dd4 * detasinv * Vh(3,iv,1)
        bc(2,5,5,iv) = bc(2,5,5,iv) - alfa * dtl(iv,1) * dd4 * detasinv * Vh(4,iv,1)

        bc(3,2,2,iv) = bc(3,2,2,iv) - alfa * dtl(iv,1) * dd5 * detasinv * Vh(1,iv,1)
        bc(3,2,3,iv) = bc(3,2,3,iv) - alfa * dtl(iv,1) * dd5 * detasinv * Vh(5,iv,1)
        bc(3,3,2,iv) = bc(3,3,2,iv) - alfa * dtl(iv,1) * dd5 * detasinv * Vh(6,iv,1)
        bc(3,3,3,iv) = bc(3,3,3,iv) - alfa * dtl(iv,1) * dd5 * detasinv * Vh(2,iv,1)
        bc(3,4,4,iv) = bc(3,4,4,iv) - alfa * dtl(iv,1) * dd5 * detasinv * Vh(3,iv,1)
        bc(3,5,5,iv) = bc(3,5,5,iv) - alfa * dtl(iv,1) * dd5 * detasinv * Vh(4,iv,1)

!.... I term

        mat(2,1,1,iv,1) = mat(2,1,1,iv,1) + one
        mat(2,2,2,iv,1) = mat(2,2,2,iv,1) + one
        mat(2,3,3,iv,1) = mat(2,3,3,iv,1) + one
        mat(2,4,4,iv,1) = mat(2,4,4,iv,1) + one
        mat(2,5,5,iv,1) = mat(2,5,5,iv,1) + one

!.... sponge

        if (ispg.eq.1) then
          mat(2,1,1,iv,1) = mat(2,1,1,iv,1) + alfa * dtl(iv,1) * spgl(1,iv)
          mat(2,2,2,iv,1) = mat(2,2,2,iv,1) + alfa * dtl(iv,1) * spgl(1,iv)
          mat(2,3,3,iv,1) = mat(2,3,3,iv,1) + alfa * dtl(iv,1) * spgl(1,iv)
          mat(2,4,4,iv,1) = mat(2,4,4,iv,1) + alfa * dtl(iv,1) * spgl(1,iv)
          mat(2,5,5,iv,1) = mat(2,5,5,iv,1) + alfa * dtl(iv,1) * spgl(1,iv)
        else if (ispg.eq.2) then
          mat(2,1,1,iv,1) = mat(2,1,1,iv,1) + alfa * dtl(iv,1) * (spgl(1,iv) + spg2l(1,iv))
          mat(2,2,2,iv,1) = mat(2,2,2,iv,1) + alfa * dtl(iv,1) * (spgl(1,iv) + spg2l(1,iv))
          mat(2,3,3,iv,1) = mat(2,3,3,iv,1) + alfa * dtl(iv,1) * (spgl(1,iv) + spg2l(1,iv))
          mat(2,4,4,iv,1) = mat(2,4,4,iv,1) + alfa * dtl(iv,1) * (spgl(1,iv) + spg2l(1,iv))
          mat(2,5,5,iv,1) = mat(2,5,5,iv,1) + alfa * dtl(iv,1) * (spgl(1,iv) + spg2l(1,iv))
        end if

!=======================================================================================================!
!.... use higher-order tangent on first node off the body
!=======================================================================================================!
        do idof = 1, ndof
          do jdof = 1, ndof
            mat(1,idof,jdof,iv,2) = alfa * dtl(iv,2) * gb1 * detainv * Bh(idof,jdof,iv,2)
            mat(2,idof,jdof,iv,2) = alfa * dtl(iv,2) * gb2 * detainv * Bh(idof,jdof,iv,2)
            mat(3,idof,jdof,iv,2) = alfa * dtl(iv,2) * gb3 * detainv * Bh(idof,jdof,iv,2)
            bc(4,idof,jdof,iv)    = alfa * dtl(iv,2) * gb4 * detainv * Bh(idof,jdof,iv,2)
            bc(5,idof,jdof,iv)    = alfa * dtl(iv,2) * gb5 * detainv * Bh(idof,jdof,iv,2)
          end do
        end do

!.... \hat{V}_{\eta\eta} term

        mat(1,2,2,iv,2) = mat(1,2,2,iv,2) - alfa * dtl(iv,2) * db1 * detasinv * Vh(1,iv,2)
        mat(1,2,3,iv,2) = mat(1,2,3,iv,2) - alfa * dtl(iv,2) * db1 * detasinv * Vh(5,iv,2)
        mat(1,3,2,iv,2) = mat(1,3,2,iv,2) - alfa * dtl(iv,2) * db1 * detasinv * Vh(6,iv,2)
        mat(1,3,3,iv,2) = mat(1,3,3,iv,2) - alfa * dtl(iv,2) * db1 * detasinv * Vh(2,iv,2)
        mat(1,4,4,iv,2) = mat(1,4,4,iv,2) - alfa * dtl(iv,2) * db1 * detasinv * Vh(3,iv,2)
        mat(1,5,5,iv,2) = mat(1,5,5,iv,2) - alfa * dtl(iv,2) * db1 * detasinv * Vh(4,iv,2)

        mat(2,2,2,iv,2) = mat(2,2,2,iv,2) - alfa * dtl(iv,2) * db2 * detasinv * Vh(1,iv,2)
        mat(2,2,3,iv,2) = mat(2,2,3,iv,2) - alfa * dtl(iv,2) * db2 * detasinv * Vh(5,iv,2)
        mat(2,3,2,iv,2) = mat(2,3,2,iv,2) - alfa * dtl(iv,2) * db2 * detasinv * Vh(6,iv,2)
        mat(2,3,3,iv,2) = mat(2,3,3,iv,2) - alfa * dtl(iv,2) * db2 * detasinv * Vh(2,iv,2)
        mat(2,4,4,iv,2) = mat(2,4,4,iv,2) - alfa * dtl(iv,2) * db2 * detasinv * Vh(3,iv,2)
        mat(2,5,5,iv,2) = mat(2,5,5,iv,2) - alfa * dtl(iv,2) * db2 * detasinv * Vh(4,iv,2)

        mat(3,2,2,iv,2) = mat(3,2,2,iv,2) - alfa * dtl(iv,2) * db3 * detasinv * Vh(1,iv,2)
        mat(3,2,3,iv,2) = mat(3,2,3,iv,2) - alfa * dtl(iv,2) * db3 * detasinv * Vh(5,iv,2)
        mat(3,3,2,iv,2) = mat(3,3,2,iv,2) - alfa * dtl(iv,2) * db3 * detasinv * Vh(6,iv,2)
        mat(3,3,3,iv,2) = mat(3,3,3,iv,2) - alfa * dtl(iv,2) * db3 * detasinv * Vh(2,iv,2)
        mat(3,4,4,iv,2) = mat(3,4,4,iv,2) - alfa * dtl(iv,2) * db3 * detasinv * Vh(3,iv,2)
        mat(3,5,5,iv,2) = mat(3,5,5,iv,2) - alfa * dtl(iv,2) * db3 * detasinv * Vh(4,iv,2)

        bc(4,2,2,iv) = bc(4,2,2,iv) - alfa * dtl(iv,2) * db4 * detasinv * Vh(1,iv,2)
        bc(4,2,3,iv) = bc(4,2,3,iv) - alfa * dtl(iv,2) * db4 * detasinv * Vh(5,iv,2)
        bc(4,3,2,iv) = bc(4,3,2,iv) - alfa * dtl(iv,2) * db4 * detasinv * Vh(6,iv,2)
        bc(4,3,3,iv) = bc(4,3,3,iv) - alfa * dtl(iv,2) * db4 * detasinv * Vh(2,iv,2)
        bc(4,4,4,iv) = bc(4,4,4,iv) - alfa * dtl(iv,2) * db4 * detasinv * Vh(3,iv,2)
        bc(4,5,5,iv) = bc(4,5,5,iv) - alfa * dtl(iv,2) * db4 * detasinv * Vh(4,iv,2)

        bc(5,2,2,iv) = bc(5,2,2,iv) - alfa * dtl(iv,2) * db5 * detasinv * Vh(1,iv,2)
        bc(5,2,3,iv) = bc(5,2,3,iv) - alfa * dtl(iv,2) * db5 * detasinv * Vh(5,iv,2)
        bc(5,3,2,iv) = bc(5,3,2,iv) - alfa * dtl(iv,2) * db5 * detasinv * Vh(6,iv,2)
        bc(5,3,3,iv) = bc(5,3,3,iv) - alfa * dtl(iv,2) * db5 * detasinv * Vh(2,iv,2)
        bc(5,4,4,iv) = bc(5,4,4,iv) - alfa * dtl(iv,2) * db5 * detasinv * Vh(3,iv,2)
        bc(5,5,5,iv) = bc(5,5,5,iv) - alfa * dtl(iv,2) * db5 * detasinv * Vh(4,iv,2)

!.... I term

        mat(2,1,1,iv,2) = mat(2,1,1,iv,2) + one
        mat(2,2,2,iv,2) = mat(2,2,2,iv,2) + one
        mat(2,3,3,iv,2) = mat(2,3,3,iv,2) + one
        mat(2,4,4,iv,2) = mat(2,4,4,iv,2) + one
        mat(2,5,5,iv,2) = mat(2,5,5,iv,2) + one

!.... sponge

        if (ispg.eq.1) then
          mat(2,1,1,iv,2) = mat(2,1,1,iv,2) + alfa * dtl(iv,2) * spgl(2,iv)
          mat(2,2,2,iv,2) = mat(2,2,2,iv,2) + alfa * dtl(iv,2) * spgl(2,iv)
          mat(2,3,3,iv,2) = mat(2,3,3,iv,2) + alfa * dtl(iv,2) * spgl(2,iv)
          mat(2,4,4,iv,2) = mat(2,4,4,iv,2) + alfa * dtl(iv,2) * spgl(2,iv)
          mat(2,5,5,iv,2) = mat(2,5,5,iv,2) + alfa * dtl(iv,2) * spgl(2,iv)
        else if (ispg.ge.2) then
          mat(2,1,1,iv,2) = mat(2,1,1,iv,2) + alfa * dtl(iv,2) * (spgl(2,iv) + spg2l(2,iv))
          mat(2,2,2,iv,2) = mat(2,2,2,iv,2) + alfa * dtl(iv,2) * (spgl(2,iv) + spg2l(2,iv))
          mat(2,3,3,iv,2) = mat(2,3,3,iv,2) + alfa * dtl(iv,2) * (spgl(2,iv) + spg2l(2,iv))
          mat(2,4,4,iv,2) = mat(2,4,4,iv,2) + alfa * dtl(iv,2) * (spgl(2,iv) + spg2l(2,iv))
          mat(2,5,5,iv,2) = mat(2,5,5,iv,2) + alfa * dtl(iv,2) * (spgl(2,iv) + spg2l(2,iv))
        end if

!.... implicit damping term

        if (eps_e.ne.zero) then
          mat(1,1,1,iv,2) = mat(1,1,1,iv,2) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * (-one)
          mat(1,2,2,iv,2) = mat(1,2,2,iv,2) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * (-one)
          mat(1,3,3,iv,2) = mat(1,3,3,iv,2) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * (-one)
          mat(1,4,4,iv,2) = mat(1,4,4,iv,2) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * (-one)
          mat(1,5,5,iv,2) = mat(1,5,5,iv,2) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * (-one)
          
          mat(2,1,1,iv,2) = mat(2,1,1,iv,2) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * two
          mat(2,2,2,iv,2) = mat(2,2,2,iv,2) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * two
          mat(2,3,3,iv,2) = mat(2,3,3,iv,2) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * two
          mat(2,4,4,iv,2) = mat(2,4,4,iv,2) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * two
          mat(2,5,5,iv,2) = mat(2,5,5,iv,2) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * two
          
          mat(3,1,1,iv,2) = mat(3,1,1,iv,2) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * (-one)
          mat(3,2,2,iv,2) = mat(3,2,2,iv,2) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * (-one)
          mat(3,3,3,iv,2) = mat(3,3,3,iv,2) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * (-one)
          mat(3,4,4,iv,2) = mat(3,4,4,iv,2) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * (-one)
          mat(3,5,5,iv,2) = mat(3,5,5,iv,2) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * (-one)
        end if

!=======================================================================================================!
!.... use higher-order tangent on second node off the body
!=======================================================================================================!
        do idof = 1, ndof
          do jdof = 1, ndof
            bc(7,idof,jdof,iv)    = alfa * dtl(iv,3) * ga1 * detainv * Bh(idof,jdof,iv,3)
            mat(1,idof,jdof,iv,3) = alfa * dtl(iv,3) * ga2 * detainv * Bh(idof,jdof,iv,3)
            mat(2,idof,jdof,iv,3) = zero
            mat(3,idof,jdof,iv,3) = alfa * dtl(iv,3) * ga3 * detainv * Bh(idof,jdof,iv,3)
            bc(6,idof,jdof,iv)    = alfa * dtl(iv,3) * ga4 * detainv * Bh(idof,jdof,iv,3)
          end do
        end do
        
!.... \hat{V}_{\eta\eta} term

        bc(7,2,2,iv)    = bc(7,2,2,iv)    - alfa * dtl(iv,3) * da1 * detasinv * Vh(1,iv,3)
        bc(7,2,3,iv)    = bc(7,2,3,iv)    - alfa * dtl(iv,3) * da1 * detasinv * Vh(5,iv,3)
        bc(7,3,2,iv)    = bc(7,3,2,iv)    - alfa * dtl(iv,3) * da1 * detasinv * Vh(6,iv,3)
        bc(7,3,3,iv)    = bc(7,3,3,iv)    - alfa * dtl(iv,3) * da1 * detasinv * Vh(2,iv,3)
        bc(7,4,4,iv)    = bc(7,4,4,iv)    - alfa * dtl(iv,3) * da1 * detasinv * Vh(3,iv,3)
        bc(7,5,5,iv)    = bc(7,5,5,iv)    - alfa * dtl(iv,3) * da1 * detasinv * Vh(4,iv,3)

        mat(1,2,2,iv,3) = mat(1,2,2,iv,3) - alfa * dtl(iv,3) * da2 * detasinv * Vh(1,iv,3)
        mat(1,2,3,iv,3) = mat(1,2,3,iv,3) - alfa * dtl(iv,3) * da2 * detasinv * Vh(5,iv,3)
        mat(1,3,2,iv,3) = mat(1,3,2,iv,3) - alfa * dtl(iv,3) * da2 * detasinv * Vh(6,iv,3)
        mat(1,3,3,iv,3) = mat(1,3,3,iv,3) - alfa * dtl(iv,3) * da2 * detasinv * Vh(2,iv,3)
        mat(1,4,4,iv,3) = mat(1,4,4,iv,3) - alfa * dtl(iv,3) * da2 * detasinv * Vh(3,iv,3)
        mat(1,5,5,iv,3) = mat(1,5,5,iv,3) - alfa * dtl(iv,3) * da2 * detasinv * Vh(4,iv,3)

        mat(2,2,2,iv,3) = mat(2,2,2,iv,3) - alfa * dtl(iv,3) * da3 * detasinv * Vh(1,iv,3)
        mat(2,2,3,iv,3) = mat(2,2,3,iv,3) - alfa * dtl(iv,3) * da3 * detasinv * Vh(5,iv,3)
        mat(2,3,2,iv,3) = mat(2,3,2,iv,3) - alfa * dtl(iv,3) * da3 * detasinv * Vh(6,iv,3)
        mat(2,3,3,iv,3) = mat(2,3,3,iv,3) - alfa * dtl(iv,3) * da3 * detasinv * Vh(2,iv,3)
        mat(2,4,4,iv,3) = mat(2,4,4,iv,3) - alfa * dtl(iv,3) * da3 * detasinv * Vh(3,iv,3)
        mat(2,5,5,iv,3) = mat(2,5,5,iv,3) - alfa * dtl(iv,3) * da3 * detasinv * Vh(4,iv,3)

        mat(3,2,2,iv,3) = mat(3,2,2,iv,3) - alfa * dtl(iv,3) * da4 * detasinv * Vh(1,iv,3)
        mat(3,2,3,iv,3) = mat(3,2,3,iv,3) - alfa * dtl(iv,3) * da4 * detasinv * Vh(5,iv,3)
        mat(3,3,2,iv,3) = mat(3,3,2,iv,3) - alfa * dtl(iv,3) * da4 * detasinv * Vh(6,iv,3)
        mat(3,3,3,iv,3) = mat(3,3,3,iv,3) - alfa * dtl(iv,3) * da4 * detasinv * Vh(2,iv,3)
        mat(3,4,4,iv,3) = mat(3,4,4,iv,3) - alfa * dtl(iv,3) * da4 * detasinv * Vh(3,iv,3)
        mat(3,5,5,iv,3) = mat(3,5,5,iv,3) - alfa * dtl(iv,3) * da4 * detasinv * Vh(4,iv,3)

        bc(6,2,2,iv)    = bc(6,2,2,iv)    - alfa * dtl(iv,3) * da5 * detasinv * Vh(1,iv,3)
        bc(6,2,3,iv)    = bc(6,2,3,iv)    - alfa * dtl(iv,3) * da5 * detasinv * Vh(5,iv,3)
        bc(6,3,2,iv)    = bc(6,3,2,iv)    - alfa * dtl(iv,3) * da5 * detasinv * Vh(6,iv,3)
        bc(6,3,3,iv)    = bc(6,3,3,iv)    - alfa * dtl(iv,3) * da5 * detasinv * Vh(2,iv,3)
        bc(6,4,4,iv)    = bc(6,4,4,iv)    - alfa * dtl(iv,3) * da5 * detasinv * Vh(3,iv,3)
        bc(6,5,5,iv)    = bc(6,5,5,iv)    - alfa * dtl(iv,3) * da5 * detasinv * Vh(4,iv,3)

!.... I term

        mat(2,1,1,iv,3) = mat(2,1,1,iv,3) + one
        mat(2,2,2,iv,3) = mat(2,2,2,iv,3) + one
        mat(2,3,3,iv,3) = mat(2,3,3,iv,3) + one
        mat(2,4,4,iv,3) = mat(2,4,4,iv,3) + one
        mat(2,5,5,iv,3) = mat(2,5,5,iv,3) + one

!.... sponge

        if (ispg.eq.1) then
          mat(2,1,1,iv,3) = mat(2,1,1,iv,3) + alfa * dtl(iv,3) * spgl(3,iv)
          mat(2,2,2,iv,3) = mat(2,2,2,iv,3) + alfa * dtl(iv,3) * spgl(3,iv)
          mat(2,3,3,iv,3) = mat(2,3,3,iv,3) + alfa * dtl(iv,3) * spgl(3,iv)
          mat(2,4,4,iv,3) = mat(2,4,4,iv,3) + alfa * dtl(iv,3) * spgl(3,iv)
          mat(2,5,5,iv,3) = mat(2,5,5,iv,3) + alfa * dtl(iv,3) * spgl(3,iv)
        else if (ispg.eq.2) then
          mat(2,1,1,iv,3) = mat(2,1,1,iv,3) + alfa * dtl(iv,3) * (spgl(3,iv) + spg2l(3,iv))
          mat(2,2,2,iv,3) = mat(2,2,2,iv,3) + alfa * dtl(iv,3) * (spgl(3,iv) + spg2l(3,iv))
          mat(2,3,3,iv,3) = mat(2,3,3,iv,3) + alfa * dtl(iv,3) * (spgl(3,iv) + spg2l(3,iv))
          mat(2,4,4,iv,3) = mat(2,4,4,iv,3) + alfa * dtl(iv,3) * (spgl(3,iv) + spg2l(3,iv))
          mat(2,5,5,iv,3) = mat(2,5,5,iv,3) + alfa * dtl(iv,3) * (spgl(3,iv) + spg2l(3,iv))
        end if

!.... implicit damping term updated to fourth order

        if (eps_e.ne.zero) then
          bc(7,1,1,iv) = bc(7,1,1,iv) + alfa * dtl(iv,3) * eps_e * buff(iv,3) * fb1
          bc(7,2,2,iv) = bc(7,2,2,iv) + alfa * dtl(iv,3) * eps_e * buff(iv,3) * fb1
          bc(7,3,3,iv) = bc(7,3,3,iv) + alfa * dtl(iv,3) * eps_e * buff(iv,3) * fb1
          bc(7,4,4,iv) = bc(7,4,4,iv) + alfa * dtl(iv,3) * eps_e * buff(iv,3) * fb1
          bc(7,5,5,iv) = bc(7,5,5,iv) + alfa * dtl(iv,3) * eps_e * buff(iv,3) * fb1
          
          mat(1,1,1,iv,3) = mat(1,1,1,iv,3) + alfa * dtl(iv,3) * eps_e * buff(iv,3) * fb2
          mat(1,2,2,iv,3) = mat(1,2,2,iv,3) + alfa * dtl(iv,3) * eps_e * buff(iv,3) * fb2
          mat(1,3,3,iv,3) = mat(1,3,3,iv,3) + alfa * dtl(iv,3) * eps_e * buff(iv,3) * fb2
          mat(1,4,4,iv,3) = mat(1,4,4,iv,3) + alfa * dtl(iv,3) * eps_e * buff(iv,3) * fb2
          mat(1,5,5,iv,3) = mat(1,5,5,iv,3) + alfa * dtl(iv,3) * eps_e * buff(iv,3) * fb2
          
          mat(2,1,1,iv,3) = mat(2,1,1,iv,3) + alfa * dtl(iv,3) * eps_e * buff(iv,3) * fb3
          mat(2,2,2,iv,3) = mat(2,2,2,iv,3) + alfa * dtl(iv,3) * eps_e * buff(iv,3) * fb3
          mat(2,3,3,iv,3) = mat(2,3,3,iv,3) + alfa * dtl(iv,3) * eps_e * buff(iv,3) * fb3
          mat(2,4,4,iv,3) = mat(2,4,4,iv,3) + alfa * dtl(iv,3) * eps_e * buff(iv,3) * fb3
          mat(2,5,5,iv,3) = mat(2,5,5,iv,3) + alfa * dtl(iv,3) * eps_e * buff(iv,3) * fb3
          
          mat(3,1,1,iv,3) = mat(3,1,1,iv,3) + alfa * dtl(iv,3) * eps_e * buff(iv,3) * fb4
          mat(3,2,2,iv,3) = mat(3,2,2,iv,3) + alfa * dtl(iv,3) * eps_e * buff(iv,3) * fb4
          mat(3,3,3,iv,3) = mat(3,3,3,iv,3) + alfa * dtl(iv,3) * eps_e * buff(iv,3) * fb4
          mat(3,4,4,iv,3) = mat(3,4,4,iv,3) + alfa * dtl(iv,3) * eps_e * buff(iv,3) * fb4
          mat(3,5,5,iv,3) = mat(3,5,5,iv,3) + alfa * dtl(iv,3) * eps_e * buff(iv,3) * fb4
          
          bc(6,1,1,iv) = bc(6,1,1,iv) + alfa * dtl(iv,3) * eps_e * buff(iv,3) * fb5
          bc(6,2,2,iv) = bc(6,2,2,iv) + alfa * dtl(iv,3) * eps_e * buff(iv,3) * fb5
          bc(6,3,3,iv) = bc(6,3,3,iv) + alfa * dtl(iv,3) * eps_e * buff(iv,3) * fb5
          bc(6,4,4,iv) = bc(6,4,4,iv) + alfa * dtl(iv,3) * eps_e * buff(iv,3) * fb5
          bc(6,5,5,iv) = bc(6,5,5,iv) + alfa * dtl(iv,3) * eps_e * buff(iv,3) * fb5
        end if

!=======================================================================================================!
!.... use higher-order tangent on second node off the far-field
!=======================================================================================================!
        do idof = 1, ndof
          do jdof = 1, ndof
            bc(8,idof,jdof,iv)       = alfa * dtl(iv,ny-2) * ga1 * detainv * Bh(idof,jdof,iv,ny-2)
            mat(1,idof,jdof,iv,ny-2) = alfa * dtl(iv,ny-2) * ga2 * detainv * Bh(idof,jdof,iv,ny-2)
            mat(2,idof,jdof,iv,ny-2) = zero
            mat(3,idof,jdof,iv,ny-2) = alfa * dtl(iv,ny-2) * ga3 * detainv * Bh(idof,jdof,iv,ny-2)
            bc(9,idof,jdof,iv)       = alfa * dtl(iv,ny-2) * ga4 * detainv * Bh(idof,jdof,iv,ny-2)
          end do
        end do
        
!.... \hat{V}_{\eta\eta} term

        bc(8,2,2,iv)       = bc(8,2,2,iv)       - alfa * dtl(iv,ny-2) * da1 * detasinv * Vh(1,iv,ny-2)
        bc(8,2,3,iv)       = bc(8,2,3,iv)       - alfa * dtl(iv,ny-2) * da1 * detasinv * Vh(5,iv,ny-2)
        bc(8,3,2,iv)       = bc(8,3,2,iv)       - alfa * dtl(iv,ny-2) * da1 * detasinv * Vh(6,iv,ny-2)
        bc(8,3,3,iv)       = bc(8,3,3,iv)       - alfa * dtl(iv,ny-2) * da1 * detasinv * Vh(2,iv,ny-2)
        bc(8,4,4,iv)       = bc(8,4,4,iv)       - alfa * dtl(iv,ny-2) * da1 * detasinv * Vh(3,iv,ny-2)
        bc(8,5,5,iv)       = bc(8,5,5,iv)       - alfa * dtl(iv,ny-2) * da1 * detasinv * Vh(4,iv,ny-2)
 
        mat(1,2,2,iv,ny-2) = mat(1,2,2,iv,ny-2) - alfa * dtl(iv,ny-2) * da2 * detasinv * Vh(1,iv,ny-2)
        mat(1,2,3,iv,ny-2) = mat(1,2,3,iv,ny-2) - alfa * dtl(iv,ny-2) * da2 * detasinv * Vh(5,iv,ny-2)
        mat(1,3,2,iv,ny-2) = mat(1,3,2,iv,ny-2) - alfa * dtl(iv,ny-2) * da2 * detasinv * Vh(6,iv,ny-2)
        mat(1,3,3,iv,ny-2) = mat(1,3,3,iv,ny-2) - alfa * dtl(iv,ny-2) * da2 * detasinv * Vh(2,iv,ny-2)
        mat(1,4,4,iv,ny-2) = mat(1,4,4,iv,ny-2) - alfa * dtl(iv,ny-2) * da2 * detasinv * Vh(3,iv,ny-2)
        mat(1,5,5,iv,ny-2) = mat(1,5,5,iv,ny-2) - alfa * dtl(iv,ny-2) * da2 * detasinv * Vh(4,iv,ny-2)

        mat(2,2,2,iv,ny-2) = mat(2,2,2,iv,ny-2) - alfa * dtl(iv,ny-2) * da3 * detasinv * Vh(1,iv,ny-2)
        mat(2,2,3,iv,ny-2) = mat(2,2,3,iv,ny-2) - alfa * dtl(iv,ny-2) * da3 * detasinv * Vh(5,iv,ny-2)
        mat(2,3,2,iv,ny-2) = mat(2,3,2,iv,ny-2) - alfa * dtl(iv,ny-2) * da3 * detasinv * Vh(6,iv,ny-2)
        mat(2,3,3,iv,ny-2) = mat(2,3,3,iv,ny-2) - alfa * dtl(iv,ny-2) * da3 * detasinv * Vh(2,iv,ny-2)
        mat(2,4,4,iv,ny-2) = mat(2,4,4,iv,ny-2) - alfa * dtl(iv,ny-2) * da3 * detasinv * Vh(3,iv,ny-2)
        mat(2,5,5,iv,ny-2) = mat(2,5,5,iv,ny-2) - alfa * dtl(iv,ny-2) * da3 * detasinv * Vh(4,iv,ny-2)

        mat(3,2,2,iv,ny-2) = mat(3,2,2,iv,ny-2) - alfa * dtl(iv,ny-2) * da4 * detasinv * Vh(1,iv,ny-2)
        mat(3,2,3,iv,ny-2) = mat(3,2,3,iv,ny-2) - alfa * dtl(iv,ny-2) * da4 * detasinv * Vh(5,iv,ny-2)
        mat(3,3,2,iv,ny-2) = mat(3,3,2,iv,ny-2) - alfa * dtl(iv,ny-2) * da4 * detasinv * Vh(6,iv,ny-2)
        mat(3,3,3,iv,ny-2) = mat(3,3,3,iv,ny-2) - alfa * dtl(iv,ny-2) * da4 * detasinv * Vh(2,iv,ny-2)
        mat(3,4,4,iv,ny-2) = mat(3,4,4,iv,ny-2) - alfa * dtl(iv,ny-2) * da4 * detasinv * Vh(3,iv,ny-2)
        mat(3,5,5,iv,ny-2) = mat(3,5,5,iv,ny-2) - alfa * dtl(iv,ny-2) * da4 * detasinv * Vh(4,iv,ny-2)

        bc(9,2,2,iv)       = bc(9,2,2,iv)       - alfa * dtl(iv,ny-2) * da5 * detasinv * Vh(1,iv,ny-2)
        bc(9,2,3,iv)       = bc(9,2,3,iv)       - alfa * dtl(iv,ny-2) * da5 * detasinv * Vh(5,iv,ny-2)
        bc(9,3,2,iv)       = bc(9,3,2,iv)       - alfa * dtl(iv,ny-2) * da5 * detasinv * Vh(6,iv,ny-2)
        bc(9,3,3,iv)       = bc(9,3,3,iv)       - alfa * dtl(iv,ny-2) * da5 * detasinv * Vh(2,iv,ny-2)
        bc(9,4,4,iv)       = bc(9,4,4,iv)       - alfa * dtl(iv,ny-2) * da5 * detasinv * Vh(3,iv,ny-2)
        bc(9,5,5,iv)       = bc(9,5,5,iv)       - alfa * dtl(iv,ny-2) * da5 * detasinv * Vh(4,iv,ny-2)

!.... I term

        mat(2,1,1,iv,ny-2) = mat(2,1,1,iv,ny-2) + one
        mat(2,2,2,iv,ny-2) = mat(2,2,2,iv,ny-2) + one
        mat(2,3,3,iv,ny-2) = mat(2,3,3,iv,ny-2) + one
        mat(2,4,4,iv,ny-2) = mat(2,4,4,iv,ny-2) + one
        mat(2,5,5,iv,ny-2) = mat(2,5,5,iv,ny-2) + one

!.... sponge

        if (ispg.eq.1) then
          mat(2,1,1,iv,ny-2) = mat(2,1,1,iv,ny-2) + alfa * dtl(iv,ny-2) * spgl(ny-2,iv)
          mat(2,2,2,iv,ny-2) = mat(2,2,2,iv,ny-2) + alfa * dtl(iv,ny-2) * spgl(ny-2,iv)
          mat(2,3,3,iv,ny-2) = mat(2,3,3,iv,ny-2) + alfa * dtl(iv,ny-2) * spgl(ny-2,iv)
          mat(2,4,4,iv,ny-2) = mat(2,4,4,iv,ny-2) + alfa * dtl(iv,ny-2) * spgl(ny-2,iv)
          mat(2,5,5,iv,ny-2) = mat(2,5,5,iv,ny-2) + alfa * dtl(iv,ny-2) * spgl(ny-2,iv)
        else if (ispg.eq.2) then
          mat(2,1,1,iv,ny-2) = mat(2,1,1,iv,ny-2) + alfa * dtl(iv,ny-2) * (spgl(ny-2,iv) + spg2l(ny-2,iv))
          mat(2,2,2,iv,ny-2) = mat(2,2,2,iv,ny-2) + alfa * dtl(iv,ny-2) * (spgl(ny-2,iv) + spg2l(ny-2,iv))
          mat(2,3,3,iv,ny-2) = mat(2,3,3,iv,ny-2) + alfa * dtl(iv,ny-2) * (spgl(ny-2,iv) + spg2l(ny-2,iv))
          mat(2,4,4,iv,ny-2) = mat(2,4,4,iv,ny-2) + alfa * dtl(iv,ny-2) * (spgl(ny-2,iv) + spg2l(ny-2,iv))
          mat(2,5,5,iv,ny-2) = mat(2,5,5,iv,ny-2) + alfa * dtl(iv,ny-2) * (spgl(ny-2,iv) + spg2l(ny-2,iv))
        end if

!.... implicit damping term updated to fourth order

        if (eps_e.ne.zero) then
          bc(8,1,1,iv) = bc(8,1,1,iv) + alfa * dtl(iv,ny-2) * eps_e * buff(iv,ny-2) * fb1
          bc(8,2,2,iv) = bc(8,2,2,iv) + alfa * dtl(iv,ny-2) * eps_e * buff(iv,ny-2) * fb1
          bc(8,3,3,iv) = bc(8,3,3,iv) + alfa * dtl(iv,ny-2) * eps_e * buff(iv,ny-2) * fb1
          bc(8,4,4,iv) = bc(8,4,4,iv) + alfa * dtl(iv,ny-2) * eps_e * buff(iv,ny-2) * fb1
          bc(8,5,5,iv) = bc(8,5,5,iv) + alfa * dtl(iv,ny-2) * eps_e * buff(iv,ny-2) * fb1
          
          mat(1,1,1,iv,ny-2) = mat(1,1,1,iv,ny-2) + alfa * dtl(iv,ny-2) * eps_e * buff(iv,ny-2) * fb2
          mat(1,2,2,iv,ny-2) = mat(1,2,2,iv,ny-2) + alfa * dtl(iv,ny-2) * eps_e * buff(iv,ny-2) * fb2
          mat(1,3,3,iv,ny-2) = mat(1,3,3,iv,ny-2) + alfa * dtl(iv,ny-2) * eps_e * buff(iv,ny-2) * fb2
          mat(1,4,4,iv,ny-2) = mat(1,4,4,iv,ny-2) + alfa * dtl(iv,ny-2) * eps_e * buff(iv,ny-2) * fb2
          mat(1,5,5,iv,ny-2) = mat(1,5,5,iv,ny-2) + alfa * dtl(iv,ny-2) * eps_e * buff(iv,ny-2) * fb2
          
          mat(2,1,1,iv,ny-2) = mat(2,1,1,iv,ny-2) + alfa * dtl(iv,ny-2) * eps_e * buff(iv,ny-2) * fb3
          mat(2,2,2,iv,ny-2) = mat(2,2,2,iv,ny-2) + alfa * dtl(iv,ny-2) * eps_e * buff(iv,ny-2) * fb3
          mat(2,3,3,iv,ny-2) = mat(2,3,3,iv,ny-2) + alfa * dtl(iv,ny-2) * eps_e * buff(iv,ny-2) * fb3
          mat(2,4,4,iv,ny-2) = mat(2,4,4,iv,ny-2) + alfa * dtl(iv,ny-2) * eps_e * buff(iv,ny-2) * fb3
          mat(2,5,5,iv,ny-2) = mat(2,5,5,iv,ny-2) + alfa * dtl(iv,ny-2) * eps_e * buff(iv,ny-2) * fb3
          
          mat(3,1,1,iv,ny-2) = mat(3,1,1,iv,ny-2) + alfa * dtl(iv,ny-2) * eps_e * buff(iv,ny-2) * fb4
          mat(3,2,2,iv,ny-2) = mat(3,2,2,iv,ny-2) + alfa * dtl(iv,ny-2) * eps_e * buff(iv,ny-2) * fb4
          mat(3,3,3,iv,ny-2) = mat(3,3,3,iv,ny-2) + alfa * dtl(iv,ny-2) * eps_e * buff(iv,ny-2) * fb4
          mat(3,4,4,iv,ny-2) = mat(3,4,4,iv,ny-2) + alfa * dtl(iv,ny-2) * eps_e * buff(iv,ny-2) * fb4
          mat(3,5,5,iv,ny-2) = mat(3,5,5,iv,ny-2) + alfa * dtl(iv,ny-2) * eps_e * buff(iv,ny-2) * fb4
          
          bc(9,1,1,iv) = bc(9,1,1,iv) + alfa * dtl(iv,ny-2) * eps_e * buff(iv,ny-2) * fb5
          bc(9,2,2,iv) = bc(9,2,2,iv) + alfa * dtl(iv,ny-2) * eps_e * buff(iv,ny-2) * fb5
          bc(9,3,3,iv) = bc(9,3,3,iv) + alfa * dtl(iv,ny-2) * eps_e * buff(iv,ny-2) * fb5
          bc(9,4,4,iv) = bc(9,4,4,iv) + alfa * dtl(iv,ny-2) * eps_e * buff(iv,ny-2) * fb5
          bc(9,5,5,iv) = bc(9,5,5,iv) + alfa * dtl(iv,ny-2) * eps_e * buff(iv,ny-2) * fb5
        end if

!=======================================================================================================!
!.... use higher-order tangent on first node off the far-field boundary
!=======================================================================================================!
        do idof = 1, ndof
          do jdof = 1, ndof
            bc(10,idof,jdof,iv)      = -alfa * dtl(iv,ny-1) * gb5 * detainv * Bh(idof,jdof,iv,ny-1)
            bc(11,idof,jdof,iv)      = -alfa * dtl(iv,ny-1) * gb4 * detainv * Bh(idof,jdof,iv,ny-1)
            mat(1,idof,jdof,iv,ny-1) = -alfa * dtl(iv,ny-1) * gb3 * detainv * Bh(idof,jdof,iv,ny-1)
            mat(2,idof,jdof,iv,ny-1) = -alfa * dtl(iv,ny-1) * gb2 * detainv * Bh(idof,jdof,iv,ny-1)
            mat(3,idof,jdof,iv,ny-1) = -alfa * dtl(iv,ny-1) * gb1 * detainv * Bh(idof,jdof,iv,ny-1)
          end do
        end do
        
!.... \hat{V}_{\eta\eta} term

        bc(10,2,2,iv) = bc(10,2,2,iv) - alfa * dtl(iv,ny-1) * db5 * detasinv * Vh(1,iv,ny-1)
        bc(10,2,3,iv) = bc(10,2,3,iv) - alfa * dtl(iv,ny-1) * db5 * detasinv * Vh(5,iv,ny-1)
        bc(10,3,2,iv) = bc(10,3,2,iv) - alfa * dtl(iv,ny-1) * db5 * detasinv * Vh(6,iv,ny-1)
        bc(10,3,3,iv) = bc(10,3,3,iv) - alfa * dtl(iv,ny-1) * db5 * detasinv * Vh(2,iv,ny-1)
        bc(10,4,4,iv) = bc(10,4,4,iv) - alfa * dtl(iv,ny-1) * db5 * detasinv * Vh(3,iv,ny-1)
        bc(10,5,5,iv) = bc(10,5,5,iv) - alfa * dtl(iv,ny-1) * db5 * detasinv * Vh(4,iv,ny-1)

        bc(11,2,2,iv) = bc(11,2,2,iv) - alfa * dtl(iv,ny-1) * db4 * detasinv * Vh(1,iv,ny-1)
        bc(11,2,3,iv) = bc(11,2,3,iv) - alfa * dtl(iv,ny-1) * db4 * detasinv * Vh(5,iv,ny-1)
        bc(11,3,2,iv) = bc(11,3,2,iv) - alfa * dtl(iv,ny-1) * db4 * detasinv * Vh(6,iv,ny-1)
        bc(11,3,3,iv) = bc(11,3,3,iv) - alfa * dtl(iv,ny-1) * db4 * detasinv * Vh(2,iv,ny-1)
        bc(11,4,4,iv) = bc(11,4,4,iv) - alfa * dtl(iv,ny-1) * db4 * detasinv * Vh(3,iv,ny-1)
        bc(11,5,5,iv) = bc(11,5,5,iv) - alfa * dtl(iv,ny-1) * db4 * detasinv * Vh(4,iv,ny-1)

        mat(1,2,2,iv,ny-1) = mat(1,2,2,iv,ny-1) - alfa * dtl(iv,ny-1) * db3 * detasinv * Vh(1,iv,ny-1)
        mat(1,2,3,iv,ny-1) = mat(1,2,3,iv,ny-1) - alfa * dtl(iv,ny-1) * db3 * detasinv * Vh(5,iv,ny-1)
        mat(1,3,2,iv,ny-1) = mat(1,3,2,iv,ny-1) - alfa * dtl(iv,ny-1) * db3 * detasinv * Vh(6,iv,ny-1)
        mat(1,3,3,iv,ny-1) = mat(1,3,3,iv,ny-1) - alfa * dtl(iv,ny-1) * db3 * detasinv * Vh(2,iv,ny-1)
        mat(1,4,4,iv,ny-1) = mat(1,4,4,iv,ny-1) - alfa * dtl(iv,ny-1) * db3 * detasinv * Vh(3,iv,ny-1)
        mat(1,5,5,iv,ny-1) = mat(1,5,5,iv,ny-1) - alfa * dtl(iv,ny-1) * db3 * detasinv * Vh(4,iv,ny-1)

        mat(2,2,2,iv,ny-1) = mat(2,2,2,iv,ny-1) - alfa * dtl(iv,ny-1) * db2 * detasinv * Vh(1,iv,ny-1)
        mat(2,2,3,iv,ny-1) = mat(2,2,3,iv,ny-1) - alfa * dtl(iv,ny-1) * db2 * detasinv * Vh(5,iv,ny-1)
        mat(2,3,2,iv,ny-1) = mat(2,3,2,iv,ny-1) - alfa * dtl(iv,ny-1) * db2 * detasinv * Vh(6,iv,ny-1)
        mat(2,3,3,iv,ny-1) = mat(2,3,3,iv,ny-1) - alfa * dtl(iv,ny-1) * db2 * detasinv * Vh(2,iv,ny-1)
        mat(2,4,4,iv,ny-1) = mat(2,4,4,iv,ny-1) - alfa * dtl(iv,ny-1) * db2 * detasinv * Vh(3,iv,ny-1)
        mat(2,5,5,iv,ny-1) = mat(2,5,5,iv,ny-1) - alfa * dtl(iv,ny-1) * db2 * detasinv * Vh(4,iv,ny-1)

        mat(3,2,2,iv,ny-1) = mat(3,2,2,iv,ny-1) - alfa * dtl(iv,ny-1) * db1 * detasinv * Vh(1,iv,ny-1)
        mat(3,2,3,iv,ny-1) = mat(3,2,3,iv,ny-1) - alfa * dtl(iv,ny-1) * db1 * detasinv * Vh(5,iv,ny-1)
        mat(3,3,2,iv,ny-1) = mat(3,3,2,iv,ny-1) - alfa * dtl(iv,ny-1) * db1 * detasinv * Vh(6,iv,ny-1)
        mat(3,3,3,iv,ny-1) = mat(3,3,3,iv,ny-1) - alfa * dtl(iv,ny-1) * db1 * detasinv * Vh(2,iv,ny-1)
        mat(3,4,4,iv,ny-1) = mat(3,4,4,iv,ny-1) - alfa * dtl(iv,ny-1) * db1 * detasinv * Vh(3,iv,ny-1)
        mat(3,5,5,iv,ny-1) = mat(3,5,5,iv,ny-1) - alfa * dtl(iv,ny-1) * db1 * detasinv * Vh(4,iv,ny-1)

!.... I term

        mat(2,1,1,iv,ny-1) = mat(2,1,1,iv,ny-1) + one
        mat(2,2,2,iv,ny-1) = mat(2,2,2,iv,ny-1) + one
        mat(2,3,3,iv,ny-1) = mat(2,3,3,iv,ny-1) + one
        mat(2,4,4,iv,ny-1) = mat(2,4,4,iv,ny-1) + one
        mat(2,5,5,iv,ny-1) = mat(2,5,5,iv,ny-1) + one

!.... sponge

        if (ispg.eq.1) then
          mat(2,1,1,iv,ny-1) = mat(2,1,1,iv,ny-1) + alfa * dtl(iv,ny-1) * spgl(ny-1,iv)
          mat(2,2,2,iv,ny-1) = mat(2,2,2,iv,ny-1) + alfa * dtl(iv,ny-1) * spgl(ny-1,iv)
          mat(2,3,3,iv,ny-1) = mat(2,3,3,iv,ny-1) + alfa * dtl(iv,ny-1) * spgl(ny-1,iv)
          mat(2,4,4,iv,ny-1) = mat(2,4,4,iv,ny-1) + alfa * dtl(iv,ny-1) * spgl(ny-1,iv)
          mat(2,5,5,iv,ny-1) = mat(2,5,5,iv,ny-1) + alfa * dtl(iv,ny-1) * spgl(ny-1,iv)
        else if (ispg.eq.2) then
          mat(2,1,1,iv,ny-1) = mat(2,1,1,iv,ny-1) + alfa * dtl(iv,ny-1) * (spgl(ny-1,iv) + spg2l(ny-1,iv))
          mat(2,2,2,iv,ny-1) = mat(2,2,2,iv,ny-1) + alfa * dtl(iv,ny-1) * (spgl(ny-1,iv) + spg2l(ny-1,iv))
          mat(2,3,3,iv,ny-1) = mat(2,3,3,iv,ny-1) + alfa * dtl(iv,ny-1) * (spgl(ny-1,iv) + spg2l(ny-1,iv))
          mat(2,4,4,iv,ny-1) = mat(2,4,4,iv,ny-1) + alfa * dtl(iv,ny-1) * (spgl(ny-1,iv) + spg2l(ny-1,iv))
          mat(2,5,5,iv,ny-1) = mat(2,5,5,iv,ny-1) + alfa * dtl(iv,ny-1) * (spgl(ny-1,iv) + spg2l(ny-1,iv))
        end if

!.... implicit damping term

        if (eps_e.ne.zero) then
          mat(1,1,1,iv,ny-1) = mat(1,1,1,iv,ny-1) + alfa * dtl(iv,ny-1) * eps_e * buff(iv,ny-1) * (-one)
          mat(1,2,2,iv,ny-1) = mat(1,2,2,iv,ny-1) + alfa * dtl(iv,ny-1) * eps_e * buff(iv,ny-1) * (-one)
          mat(1,3,3,iv,ny-1) = mat(1,3,3,iv,ny-1) + alfa * dtl(iv,ny-1) * eps_e * buff(iv,ny-1) * (-one)
          mat(1,4,4,iv,ny-1) = mat(1,4,4,iv,ny-1) + alfa * dtl(iv,ny-1) * eps_e * buff(iv,ny-1) * (-one)
          mat(1,5,5,iv,ny-1) = mat(1,5,5,iv,ny-1) + alfa * dtl(iv,ny-1) * eps_e * buff(iv,ny-1) * (-one)
          
          mat(2,1,1,iv,ny-1) = mat(2,1,1,iv,ny-1) + alfa * dtl(iv,ny-1) * eps_e * buff(iv,ny-1) * two
          mat(2,2,2,iv,ny-1) = mat(2,2,2,iv,ny-1) + alfa * dtl(iv,ny-1) * eps_e * buff(iv,ny-1) * two
          mat(2,3,3,iv,ny-1) = mat(2,3,3,iv,ny-1) + alfa * dtl(iv,ny-1) * eps_e * buff(iv,ny-1) * two
          mat(2,4,4,iv,ny-1) = mat(2,4,4,iv,ny-1) + alfa * dtl(iv,ny-1) * eps_e * buff(iv,ny-1) * two
          mat(2,5,5,iv,ny-1) = mat(2,5,5,iv,ny-1) + alfa * dtl(iv,ny-1) * eps_e * buff(iv,ny-1) * two
          
          mat(3,1,1,iv,ny-1) = mat(3,1,1,iv,ny-1) + alfa * dtl(iv,ny-1) * eps_e * buff(iv,ny-1) * (-one)
          mat(3,2,2,iv,ny-1) = mat(3,2,2,iv,ny-1) + alfa * dtl(iv,ny-1) * eps_e * buff(iv,ny-1) * (-one)
          mat(3,3,3,iv,ny-1) = mat(3,3,3,iv,ny-1) + alfa * dtl(iv,ny-1) * eps_e * buff(iv,ny-1) * (-one)
          mat(3,4,4,iv,ny-1) = mat(3,4,4,iv,ny-1) + alfa * dtl(iv,ny-1) * eps_e * buff(iv,ny-1) * (-one)
          mat(3,5,5,iv,ny-1) = mat(3,5,5,iv,ny-1) + alfa * dtl(iv,ny-1) * eps_e * buff(iv,ny-1) * (-one)
        end if

!=======================================================================================================!
!.... use higher-order tangent on the far-field boundary
!=======================================================================================================!
        do idof = 1, ndof
          do jdof = 1, ndof
            bc(12,idof,jdof,iv)    = -alfa * dtl(iv,ny) * gc5 * detainv * Bh(idof,jdof,iv,ny)
            bc(13,idof,jdof,iv)    = -alfa * dtl(iv,ny) * gc4 * detainv * Bh(idof,jdof,iv,ny)
            bc(14,idof,jdof,iv)    = -alfa * dtl(iv,ny) * gc3 * detainv * Bh(idof,jdof,iv,ny)
            mat(1,idof,jdof,iv,ny) = -alfa * dtl(iv,ny) * gc2 * detainv * Bh(idof,jdof,iv,ny)
            mat(2,idof,jdof,iv,ny) = -alfa * dtl(iv,ny) * gc1 * detainv * Bh(idof,jdof,iv,ny)
            mat(3,idof,jdof,iv,ny) = zero
          end do
        end do
        
!.... \hat{V}_{\eta\eta} term

        bc(12,2,2,iv) = bc(12,2,2,iv) - alfa * dtl(iv,ny) * dd5 * detasinv * Vh(1,iv,ny)
        bc(12,2,3,iv) = bc(12,2,3,iv) - alfa * dtl(iv,ny) * dd5 * detasinv * Vh(5,iv,ny)
        bc(12,3,2,iv) = bc(12,3,2,iv) - alfa * dtl(iv,ny) * dd5 * detasinv * Vh(6,iv,ny)
        bc(12,3,3,iv) = bc(12,3,3,iv) - alfa * dtl(iv,ny) * dd5 * detasinv * Vh(2,iv,ny)
        bc(12,4,4,iv) = bc(12,4,4,iv) - alfa * dtl(iv,ny) * dd5 * detasinv * Vh(3,iv,ny)
        bc(12,5,5,iv) = bc(12,5,5,iv) - alfa * dtl(iv,ny) * dd5 * detasinv * Vh(4,iv,ny)

        bc(13,2,2,iv) = bc(13,2,2,iv) - alfa * dtl(iv,ny) * dd4 * detasinv * Vh(1,iv,ny)
        bc(13,2,3,iv) = bc(13,2,3,iv) - alfa * dtl(iv,ny) * dd4 * detasinv * Vh(5,iv,ny)
        bc(13,3,2,iv) = bc(13,3,2,iv) - alfa * dtl(iv,ny) * dd4 * detasinv * Vh(6,iv,ny)
        bc(13,3,3,iv) = bc(13,3,3,iv) - alfa * dtl(iv,ny) * dd4 * detasinv * Vh(2,iv,ny)
        bc(13,4,4,iv) = bc(13,4,4,iv) - alfa * dtl(iv,ny) * dd4 * detasinv * Vh(3,iv,ny)
        bc(13,5,5,iv) = bc(13,5,5,iv) - alfa * dtl(iv,ny) * dd4 * detasinv * Vh(4,iv,ny)

        bc(14,2,2,iv) = bc(14,2,2,iv) - alfa * dtl(iv,ny) * dd3 * detasinv * Vh(1,iv,ny)
        bc(14,2,3,iv) = bc(14,2,3,iv) - alfa * dtl(iv,ny) * dd3 * detasinv * Vh(5,iv,ny)
        bc(14,3,2,iv) = bc(14,3,2,iv) - alfa * dtl(iv,ny) * dd3 * detasinv * Vh(6,iv,ny)
        bc(14,3,3,iv) = bc(14,3,3,iv) - alfa * dtl(iv,ny) * dd3 * detasinv * Vh(2,iv,ny)
        bc(14,4,4,iv) = bc(14,4,4,iv) - alfa * dtl(iv,ny) * dd3 * detasinv * Vh(3,iv,ny)
        bc(14,5,5,iv) = bc(14,5,5,iv) - alfa * dtl(iv,ny) * dd3 * detasinv * Vh(4,iv,ny)

        mat(1,2,2,iv,ny) = mat(1,2,2,iv,ny) - alfa * dtl(iv,ny) * dd2 * detasinv * Vh(1,iv,ny)
        mat(1,2,3,iv,ny) = mat(1,2,3,iv,ny) - alfa * dtl(iv,ny) * dd2 * detasinv * Vh(5,iv,ny)
        mat(1,3,2,iv,ny) = mat(1,3,2,iv,ny) - alfa * dtl(iv,ny) * dd2 * detasinv * Vh(6,iv,ny)
        mat(1,3,3,iv,ny) = mat(1,3,3,iv,ny) - alfa * dtl(iv,ny) * dd2 * detasinv * Vh(2,iv,ny)
        mat(1,4,4,iv,ny) = mat(1,4,4,iv,ny) - alfa * dtl(iv,ny) * dd2 * detasinv * Vh(3,iv,ny)
        mat(1,5,5,iv,ny) = mat(1,5,5,iv,ny) - alfa * dtl(iv,ny) * dd2 * detasinv * Vh(4,iv,ny)

        mat(2,2,2,iv,ny) = mat(2,2,2,iv,ny) - alfa * dtl(iv,ny) * dd1 * detasinv * Vh(1,iv,ny)
        mat(2,2,3,iv,ny) = mat(2,2,3,iv,ny) - alfa * dtl(iv,ny) * dd1 * detasinv * Vh(5,iv,ny)
        mat(2,3,2,iv,ny) = mat(2,3,2,iv,ny) - alfa * dtl(iv,ny) * dd1 * detasinv * Vh(6,iv,ny)
        mat(2,3,3,iv,ny) = mat(2,3,3,iv,ny) - alfa * dtl(iv,ny) * dd1 * detasinv * Vh(2,iv,ny)
        mat(2,4,4,iv,ny) = mat(2,4,4,iv,ny) - alfa * dtl(iv,ny) * dd1 * detasinv * Vh(3,iv,ny)
        mat(2,5,5,iv,ny) = mat(2,5,5,iv,ny) - alfa * dtl(iv,ny) * dd1 * detasinv * Vh(4,iv,ny)

!.... I term

        mat(2,1,1,iv,ny) = mat(2,1,1,iv,ny) + one
        mat(2,2,2,iv,ny) = mat(2,2,2,iv,ny) + one
        mat(2,3,3,iv,ny) = mat(2,3,3,iv,ny) + one
        mat(2,4,4,iv,ny) = mat(2,4,4,iv,ny) + one
        mat(2,5,5,iv,ny) = mat(2,5,5,iv,ny) + one

!.... sponge

        if (ispg.eq.1) then
          mat(2,1,1,iv,ny) = mat(2,1,1,iv,ny) + alfa * dtl(iv,ny) * spgl(ny,iv)
          mat(2,2,2,iv,ny) = mat(2,2,2,iv,ny) + alfa * dtl(iv,ny) * spgl(ny,iv)
          mat(2,3,3,iv,ny) = mat(2,3,3,iv,ny) + alfa * dtl(iv,ny) * spgl(ny,iv)
          mat(2,4,4,iv,ny) = mat(2,4,4,iv,ny) + alfa * dtl(iv,ny) * spgl(ny,iv)
          mat(2,5,5,iv,ny) = mat(2,5,5,iv,ny) + alfa * dtl(iv,ny) * spgl(ny,iv)
        else if (ispg.eq.2) then
          mat(2,1,1,iv,ny) = mat(2,1,1,iv,ny) + alfa * dtl(iv,ny) * (spgl(ny,iv) + spg2l(ny,iv))
          mat(2,2,2,iv,ny) = mat(2,2,2,iv,ny) + alfa * dtl(iv,ny) * (spgl(ny,iv) + spg2l(ny,iv))
          mat(2,3,3,iv,ny) = mat(2,3,3,iv,ny) + alfa * dtl(iv,ny) * (spgl(ny,iv) + spg2l(ny,iv))
          mat(2,4,4,iv,ny) = mat(2,4,4,iv,ny) + alfa * dtl(iv,ny) * (spgl(ny,iv) + spg2l(ny,iv))
          mat(2,5,5,iv,ny) = mat(2,5,5,iv,ny) + alfa * dtl(iv,ny) * (spgl(ny,iv) + spg2l(ny,iv))
        end if

        end do

        end if          ! yper

        return
        end subroutine lhsbt2
