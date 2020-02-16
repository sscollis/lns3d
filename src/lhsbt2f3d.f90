!=======================================================================================================!
        subroutine lhsbt2f3d( mat, Bh, Dh, Dhi, Vh, Bhi, spgl, spg2l, dtl, calcd )
!  
!  Correct the LHS for boundary treatment in the \eta direction
!
!  This version supports the 4th order LHS
!
!  Revised: 4-23-96
!=======================================================================================================!
        use global
        use stencil
        use buffer
        implicit none
        
        complex :: mat(5,ndof,ndof,nx,ny)
        real    :: Bh(ndof,ndof,nx,ny), Dh(ndof,ndof,nx,ny)
        real    :: Dhi(ndof,ndof,nx,ny), Bhi(ny,nx,6)
        real    :: spgl(ny,nx), spg2l(ny,nx), Vh(6,nx,ny)
        real    :: dtl(nx,ny)
        logical :: calcd

        real    :: a1, a2, a3, a4, a5
        complex :: c1, c2, c3, c4, c5

        real :: detainv, detasinv, rcalcd
        
        integer :: iv, j, idof, jdof
!=======================================================================================================!
        rcalcd = zero
        if (calcd) rcalcd = one

        if (yper) return

        detainv  = one / deta
        detasinv = one / deta**2
        
        !$omp parallel do private (j,idof,jdof,a1,a2,a3,a4,a5,c1,c2,c3,c4,c5)
        do iv = 1, nx
!=======================================================================================================!
!.... use higher-order tangent on body node
!=======================================================================================================!
        j = 1

        a1 = alfa * gc1 * detainv
        a2 = alfa * gc2 * detainv
        a3 = alfa * gc3 * detainv
        a4 = alfa * gc4 * detainv
        a5 = alfa * gc5 * detainv

        do idof = 1, ndof
          do jdof = 1, ndof
              mat(3,idof,jdof,iv,j) = dtl(iv,j) * ( a1 * Bh(idof,jdof,iv,j) + &
                                      rcalcd * alfa * (  Dh(idof,jdof,iv,j) + &
                                                   im * Dhi(idof,jdof,iv,j) ) )
              mat(4,idof,jdof,iv,j) = a2 * dtl(iv,j) * Bh(idof,jdof,iv,j)
              mat(5,idof,jdof,iv,j) = a3 * dtl(iv,j) * Bh(idof,jdof,iv,j)
              mat(1,idof,jdof,iv,j) = a4 * dtl(iv,j) * Bh(idof,jdof,iv,j)
              mat(2,idof,jdof,iv,j) = a5 * dtl(iv,j) * Bh(idof,jdof,iv,j)
          end do
        end do
        
!.... \hat{B}_i term

        c1 = im * alfa * gc1 * detainv
        c2 = im * alfa * gc2 * detainv
        c3 = im * alfa * gc3 * detainv
        c4 = im * alfa * gc4 * detainv
        c5 = im * alfa * gc5 * detainv

        mat(3,2,4,iv,j) = mat(3,2,4,iv,j) - c1 * dtl(iv,j) * Bhi(j,iv,3)
        mat(3,3,4,iv,j) = mat(3,3,4,iv,j) - c1 * dtl(iv,j) * Bhi(j,iv,4)
        mat(3,4,2,iv,j) = mat(3,4,2,iv,j) - c1 * dtl(iv,j) * Bhi(j,iv,3)
        mat(3,4,3,iv,j) = mat(3,4,3,iv,j) - c1 * dtl(iv,j) * Bhi(j,iv,4)

        mat(4,2,4,iv,j) = mat(4,2,4,iv,j) - c2 * dtl(iv,j) * Bhi(j,iv,3)
        mat(4,3,4,iv,j) = mat(4,3,4,iv,j) - c2 * dtl(iv,j) * Bhi(j,iv,4)
        mat(4,4,2,iv,j) = mat(4,4,2,iv,j) - c2 * dtl(iv,j) * Bhi(j,iv,3)
        mat(4,4,3,iv,j) = mat(4,4,3,iv,j) - c2 * dtl(iv,j) * Bhi(j,iv,4)

        mat(5,2,4,iv,j) = mat(5,2,4,iv,j) - c3 * dtl(iv,j) * Bhi(j,iv,3)
        mat(5,3,4,iv,j) = mat(5,3,4,iv,j) - c3 * dtl(iv,j) * Bhi(j,iv,4)
        mat(5,4,2,iv,j) = mat(5,4,2,iv,j) - c3 * dtl(iv,j) * Bhi(j,iv,3)
        mat(5,4,3,iv,j) = mat(5,4,3,iv,j) - c3 * dtl(iv,j) * Bhi(j,iv,4)

        mat(1,2,4,iv,j) = mat(1,2,4,iv,j) - c4 * dtl(iv,j) * Bhi(j,iv,3)
        mat(1,3,4,iv,j) = mat(1,3,4,iv,j) - c4 * dtl(iv,j) * Bhi(j,iv,4)
        mat(1,4,2,iv,j) = mat(1,4,2,iv,j) - c4 * dtl(iv,j) * Bhi(j,iv,3)
        mat(1,4,3,iv,j) = mat(1,4,3,iv,j) - c4 * dtl(iv,j) * Bhi(j,iv,4)

        mat(2,2,4,iv,j) = mat(2,2,4,iv,j) - c5 * dtl(iv,j) * Bhi(j,iv,3)
        mat(2,3,4,iv,j) = mat(2,3,4,iv,j) - c5 * dtl(iv,j) * Bhi(j,iv,4)
        mat(2,4,2,iv,j) = mat(2,4,2,iv,j) - c5 * dtl(iv,j) * Bhi(j,iv,3)
        mat(2,4,3,iv,j) = mat(2,4,3,iv,j) - c5 * dtl(iv,j) * Bhi(j,iv,4)
          
!.... \hat{V}_{\eta\eta} term

        a1 = alfa * dd1 * detasinv
        a2 = alfa * dd2 * detasinv
        a3 = alfa * dd3 * detasinv
        a4 = alfa * dd4 * detasinv
        a5 = alfa * dd5 * detasinv

        mat(3,2,2,iv,j) = mat(3,2,2,iv,j) - a1 * dtl(iv,j) * Vh(1,iv,j)
        mat(3,2,3,iv,j) = mat(3,2,3,iv,j) - a1 * dtl(iv,j) * Vh(5,iv,j)
        mat(3,3,2,iv,j) = mat(3,3,2,iv,j) - a1 * dtl(iv,j) * Vh(6,iv,j)
        mat(3,3,3,iv,j) = mat(3,3,3,iv,j) - a1 * dtl(iv,j) * Vh(2,iv,j)
        mat(3,4,4,iv,j) = mat(3,4,4,iv,j) - a1 * dtl(iv,j) * Vh(3,iv,j)
        mat(3,5,5,iv,j) = mat(3,5,5,iv,j) - a1 * dtl(iv,j) * Vh(4,iv,j)

        mat(4,2,2,iv,j) = mat(4,2,2,iv,j) - a2 * dtl(iv,j) * Vh(1,iv,j)
        mat(4,2,3,iv,j) = mat(4,2,3,iv,j) - a2 * dtl(iv,j) * Vh(5,iv,j)
        mat(4,3,2,iv,j) = mat(4,3,2,iv,j) - a2 * dtl(iv,j) * Vh(6,iv,j)
        mat(4,3,3,iv,j) = mat(4,3,3,iv,j) - a2 * dtl(iv,j) * Vh(2,iv,j)
        mat(4,4,4,iv,j) = mat(4,4,4,iv,j) - a2 * dtl(iv,j) * Vh(3,iv,j)
        mat(4,5,5,iv,j) = mat(4,5,5,iv,j) - a2 * dtl(iv,j) * Vh(4,iv,j)

        mat(5,2,2,iv,j) = mat(5,2,2,iv,j) - a3 * dtl(iv,j) * Vh(1,iv,j)
        mat(5,2,3,iv,j) = mat(5,2,3,iv,j) - a3 * dtl(iv,j) * Vh(5,iv,j)
        mat(5,3,2,iv,j) = mat(5,3,2,iv,j) - a3 * dtl(iv,j) * Vh(6,iv,j)
        mat(5,3,3,iv,j) = mat(5,3,3,iv,j) - a3 * dtl(iv,j) * Vh(2,iv,j)
        mat(5,4,4,iv,j) = mat(5,4,4,iv,j) - a3 * dtl(iv,j) * Vh(3,iv,j)
        mat(5,5,5,iv,j) = mat(5,5,5,iv,j) - a3 * dtl(iv,j) * Vh(4,iv,j)

        mat(1,2,2,iv,j) = mat(1,2,2,iv,j) - a4 * dtl(iv,j) * Vh(1,iv,j)
        mat(1,2,3,iv,j) = mat(1,2,3,iv,j) - a4 * dtl(iv,j) * Vh(5,iv,j)
        mat(1,3,2,iv,j) = mat(1,3,2,iv,j) - a4 * dtl(iv,j) * Vh(6,iv,j)
        mat(1,3,3,iv,j) = mat(1,3,3,iv,j) - a4 * dtl(iv,j) * Vh(2,iv,j)
        mat(1,4,4,iv,j) = mat(1,4,4,iv,j) - a4 * dtl(iv,j) * Vh(3,iv,j)
        mat(1,5,5,iv,j) = mat(1,5,5,iv,j) - a4 * dtl(iv,j) * Vh(4,iv,j)

        mat(2,2,2,iv,j) = mat(2,2,2,iv,j) - a5 * dtl(iv,j) * Vh(1,iv,j)
        mat(2,2,3,iv,j) = mat(2,2,3,iv,j) - a5 * dtl(iv,j) * Vh(5,iv,j)
        mat(2,3,2,iv,j) = mat(2,3,2,iv,j) - a5 * dtl(iv,j) * Vh(6,iv,j)
        mat(2,3,3,iv,j) = mat(2,3,3,iv,j) - a5 * dtl(iv,j) * Vh(2,iv,j)
        mat(2,4,4,iv,j) = mat(2,4,4,iv,j) - a5 * dtl(iv,j) * Vh(3,iv,j)
        mat(2,5,5,iv,j) = mat(2,5,5,iv,j) - a5 * dtl(iv,j) * Vh(4,iv,j)

!.... I term

        mat(3,1,1,iv,j) = mat(3,1,1,iv,j) + one
        mat(3,2,2,iv,j) = mat(3,2,2,iv,j) + one
        mat(3,3,3,iv,j) = mat(3,3,3,iv,j) + one
        mat(3,4,4,iv,j) = mat(3,4,4,iv,j) + one
        mat(3,5,5,iv,j) = mat(3,5,5,iv,j) + one

!.... sponge term

        if (.not. calcd) then

        if (ispg .eq. 1) then
        
            mat(3,1,1,iv,j) = mat(3,1,1,iv,j) + alfa * dtl(iv,j) * spgl(j,iv)
            mat(3,2,2,iv,j) = mat(3,2,2,iv,j) + alfa * dtl(iv,j) * spgl(j,iv)
            mat(3,3,3,iv,j) = mat(3,3,3,iv,j) + alfa * dtl(iv,j) * spgl(j,iv)
            mat(3,4,4,iv,j) = mat(3,4,4,iv,j) + alfa * dtl(iv,j) * spgl(j,iv)
            mat(3,5,5,iv,j) = mat(3,5,5,iv,j) + alfa * dtl(iv,j) * spgl(j,iv)

        else if (ispg .ge. 2) then
        
            mat(3,1,1,iv,j) = mat(3,1,1,iv,j) + alfa * dtl(iv,j) * (spgl(j,iv) + spg2l(j,iv))
            mat(3,2,2,iv,j) = mat(3,2,2,iv,j) + alfa * dtl(iv,j) * (spgl(j,iv) + spg2l(j,iv))
            mat(3,3,3,iv,j) = mat(3,3,3,iv,j) + alfa * dtl(iv,j) * (spgl(j,iv) + spg2l(j,iv))
            mat(3,4,4,iv,j) = mat(3,4,4,iv,j) + alfa * dtl(iv,j) * (spgl(j,iv) + spg2l(j,iv))
            mat(3,5,5,iv,j) = mat(3,5,5,iv,j) + alfa * dtl(iv,j) * (spgl(j,iv) + spg2l(j,iv))
        
        end if

        end if
!=======================================================================================================!
!.... use higher-order tangent on first node off the body
!=======================================================================================================!
        j = 2

        a1 = alfa * gb1 * detainv
        a2 = alfa * gb2 * detainv
        a3 = alfa * gb3 * detainv
        a4 = alfa * gb4 * detainv
        a5 = alfa * gb5 * detainv

        do idof = 1, ndof
          do jdof = 1, ndof
              mat(2,idof,jdof,iv,j) = a1 * dtl(iv,j) * Bh(idof,jdof,iv,j)
              mat(3,idof,jdof,iv,j) = dtl(iv,j) * ( a2 * Bh(idof,jdof,iv,j) + &
                                      rcalcd * alfa * (  Dh(idof,jdof,iv,j) + &
                                                   im * Dhi(idof,jdof,iv,j) ) )
              mat(4,idof,jdof,iv,j) = a3 * dtl(iv,j) * Bh(idof,jdof,iv,j)
              mat(5,idof,jdof,iv,j) = a4 * dtl(iv,j) * Bh(idof,jdof,iv,j)
              mat(1,idof,jdof,iv,j) = a5 * dtl(iv,j) * Bh(idof,jdof,iv,j)
          end do
        end do
        
!.... \hat{B}_i term

        c1 = im * alfa * gb1 * detainv
        c2 = im * alfa * gb2 * detainv
        c3 = im * alfa * gb3 * detainv
        c4 = im * alfa * gb4 * detainv
        c5 = im * alfa * gb5 * detainv

        mat(2,2,4,iv,j) = mat(2,2,4,iv,j) - c1 * dtl(iv,j) * Bhi(j,iv,3)
        mat(2,3,4,iv,j) = mat(2,3,4,iv,j) - c1 * dtl(iv,j) * Bhi(j,iv,4)
        mat(2,4,2,iv,j) = mat(2,4,2,iv,j) - c1 * dtl(iv,j) * Bhi(j,iv,3)
        mat(2,4,3,iv,j) = mat(2,4,3,iv,j) - c1 * dtl(iv,j) * Bhi(j,iv,4)

        mat(3,2,4,iv,j) = mat(3,2,4,iv,j) - c2 * dtl(iv,j) * Bhi(j,iv,3)
        mat(3,3,4,iv,j) = mat(3,3,4,iv,j) - c2 * dtl(iv,j) * Bhi(j,iv,4)
        mat(3,4,2,iv,j) = mat(3,4,2,iv,j) - c2 * dtl(iv,j) * Bhi(j,iv,3)
        mat(3,4,3,iv,j) = mat(3,4,3,iv,j) - c2 * dtl(iv,j) * Bhi(j,iv,4)

        mat(4,2,4,iv,j) = mat(4,2,4,iv,j) - c3 * dtl(iv,j) * Bhi(j,iv,3)
        mat(4,3,4,iv,j) = mat(4,3,4,iv,j) - c3 * dtl(iv,j) * Bhi(j,iv,4)
        mat(4,4,2,iv,j) = mat(4,4,2,iv,j) - c3 * dtl(iv,j) * Bhi(j,iv,3)
        mat(4,4,3,iv,j) = mat(4,4,3,iv,j) - c3 * dtl(iv,j) * Bhi(j,iv,4)

        mat(5,2,4,iv,j) = mat(5,2,4,iv,j) - c4 * dtl(iv,j) * Bhi(j,iv,3)
        mat(5,3,4,iv,j) = mat(5,3,4,iv,j) - c4 * dtl(iv,j) * Bhi(j,iv,4)
        mat(5,4,2,iv,j) = mat(5,4,2,iv,j) - c4 * dtl(iv,j) * Bhi(j,iv,3)
        mat(5,4,3,iv,j) = mat(5,4,3,iv,j) - c4 * dtl(iv,j) * Bhi(j,iv,4)

        mat(1,2,4,iv,j) = mat(1,2,4,iv,j) - c5 * dtl(iv,j) * Bhi(j,iv,3)
        mat(1,3,4,iv,j) = mat(1,3,4,iv,j) - c5 * dtl(iv,j) * Bhi(j,iv,4)
        mat(1,4,2,iv,j) = mat(1,4,2,iv,j) - c5 * dtl(iv,j) * Bhi(j,iv,3)
        mat(1,4,3,iv,j) = mat(1,4,3,iv,j) - c5 * dtl(iv,j) * Bhi(j,iv,4)
          
!.... \hat{V}_{\eta\eta} term

        a1 = alfa * db1 * detasinv
        a2 = alfa * db2 * detasinv
        a3 = alfa * db3 * detasinv
        a4 = alfa * db4 * detasinv
        a5 = alfa * db5 * detasinv

        mat(2,2,2,iv,j) = mat(2,2,2,iv,j) - a1 * dtl(iv,j) * Vh(1,iv,j)
        mat(2,2,3,iv,j) = mat(2,2,3,iv,j) - a1 * dtl(iv,j) * Vh(5,iv,j)
        mat(2,3,2,iv,j) = mat(2,3,2,iv,j) - a1 * dtl(iv,j) * Vh(6,iv,j)
        mat(2,3,3,iv,j) = mat(2,3,3,iv,j) - a1 * dtl(iv,j) * Vh(2,iv,j)
        mat(2,4,4,iv,j) = mat(2,4,4,iv,j) - a1 * dtl(iv,j) * Vh(3,iv,j)
        mat(2,5,5,iv,j) = mat(2,5,5,iv,j) - a1 * dtl(iv,j) * Vh(4,iv,j)

        mat(3,2,2,iv,j) = mat(3,2,2,iv,j) - a2 * dtl(iv,j) * Vh(1,iv,j)
        mat(3,2,3,iv,j) = mat(3,2,3,iv,j) - a2 * dtl(iv,j) * Vh(5,iv,j)
        mat(3,3,2,iv,j) = mat(3,3,2,iv,j) - a2 * dtl(iv,j) * Vh(6,iv,j)
        mat(3,3,3,iv,j) = mat(3,3,3,iv,j) - a2 * dtl(iv,j) * Vh(2,iv,j)
        mat(3,4,4,iv,j) = mat(3,4,4,iv,j) - a2 * dtl(iv,j) * Vh(3,iv,j)
        mat(3,5,5,iv,j) = mat(3,5,5,iv,j) - a2 * dtl(iv,j) * Vh(4,iv,j)

        mat(4,2,2,iv,j) = mat(4,2,2,iv,j) - a3 * dtl(iv,j) * Vh(1,iv,j)
        mat(4,2,3,iv,j) = mat(4,2,3,iv,j) - a3 * dtl(iv,j) * Vh(5,iv,j)
        mat(4,3,2,iv,j) = mat(4,3,2,iv,j) - a3 * dtl(iv,j) * Vh(6,iv,j)
        mat(4,3,3,iv,j) = mat(4,3,3,iv,j) - a3 * dtl(iv,j) * Vh(2,iv,j)
        mat(4,4,4,iv,j) = mat(4,4,4,iv,j) - a3 * dtl(iv,j) * Vh(3,iv,j)
        mat(4,5,5,iv,j) = mat(4,5,5,iv,j) - a3 * dtl(iv,j) * Vh(4,iv,j)

        mat(5,2,2,iv,j) = mat(5,2,2,iv,j) - a4 * dtl(iv,j) * Vh(1,iv,j)
        mat(5,2,3,iv,j) = mat(5,2,3,iv,j) - a4 * dtl(iv,j) * Vh(5,iv,j)
        mat(5,3,2,iv,j) = mat(5,3,2,iv,j) - a4 * dtl(iv,j) * Vh(6,iv,j)
        mat(5,3,3,iv,j) = mat(5,3,3,iv,j) - a4 * dtl(iv,j) * Vh(2,iv,j)
        mat(5,4,4,iv,j) = mat(5,4,4,iv,j) - a4 * dtl(iv,j) * Vh(3,iv,j)
        mat(5,5,5,iv,j) = mat(5,5,5,iv,j) - a4 * dtl(iv,j) * Vh(4,iv,j)

        mat(1,2,2,iv,j) = mat(1,2,2,iv,j) - a5 * dtl(iv,j) * Vh(1,iv,j)
        mat(1,2,3,iv,j) = mat(1,2,3,iv,j) - a5 * dtl(iv,j) * Vh(5,iv,j)
        mat(1,3,2,iv,j) = mat(1,3,2,iv,j) - a5 * dtl(iv,j) * Vh(6,iv,j)
        mat(1,3,3,iv,j) = mat(1,3,3,iv,j) - a5 * dtl(iv,j) * Vh(2,iv,j)
        mat(1,4,4,iv,j) = mat(1,4,4,iv,j) - a5 * dtl(iv,j) * Vh(3,iv,j)
        mat(1,5,5,iv,j) = mat(1,5,5,iv,j) - a5 * dtl(iv,j) * Vh(4,iv,j)

!.... I term

        mat(3,1,1,iv,j) = mat(3,1,1,iv,j) + one
        mat(3,2,2,iv,j) = mat(3,2,2,iv,j) + one
        mat(3,3,3,iv,j) = mat(3,3,3,iv,j) + one
        mat(3,4,4,iv,j) = mat(3,4,4,iv,j) + one
        mat(3,5,5,iv,j) = mat(3,5,5,iv,j) + one

!.... implicit damping term

        if (eps_e .ne. zero) then

          mat(2,1,1,iv,j) = mat(2,1,1,iv,j) + alfa * dtl(iv,j) * eps_e * buff(iv,j) * (-one)
          mat(2,2,2,iv,j) = mat(2,2,2,iv,j) + alfa * dtl(iv,j) * eps_e * buff(iv,j) * (-one)
          mat(2,3,3,iv,j) = mat(2,3,3,iv,j) + alfa * dtl(iv,j) * eps_e * buff(iv,j) * (-one)
          mat(2,4,4,iv,j) = mat(2,4,4,iv,j) + alfa * dtl(iv,j) * eps_e * buff(iv,j) * (-one)
          mat(2,5,5,iv,j) = mat(2,5,5,iv,j) + alfa * dtl(iv,j) * eps_e * buff(iv,j) * (-one)
  
          mat(3,1,1,iv,j) = mat(3,1,1,iv,j) + alfa * dtl(iv,j) * eps_e * buff(iv,j) * two
          mat(3,2,2,iv,j) = mat(3,2,2,iv,j) + alfa * dtl(iv,j) * eps_e * buff(iv,j) * two
          mat(3,3,3,iv,j) = mat(3,3,3,iv,j) + alfa * dtl(iv,j) * eps_e * buff(iv,j) * two
          mat(3,4,4,iv,j) = mat(3,4,4,iv,j) + alfa * dtl(iv,j) * eps_e * buff(iv,j) * two
          mat(3,5,5,iv,j) = mat(3,5,5,iv,j) + alfa * dtl(iv,j) * eps_e * buff(iv,j) * two
  
          mat(4,1,1,iv,j) = mat(4,1,1,iv,j) + alfa * dtl(iv,j) * eps_e * buff(iv,j) * (-one)
          mat(4,2,2,iv,j) = mat(4,2,2,iv,j) + alfa * dtl(iv,j) * eps_e * buff(iv,j) * (-one)
          mat(4,3,3,iv,j) = mat(4,3,3,iv,j) + alfa * dtl(iv,j) * eps_e * buff(iv,j) * (-one)
          mat(4,4,4,iv,j) = mat(4,4,4,iv,j) + alfa * dtl(iv,j) * eps_e * buff(iv,j) * (-one)
          mat(4,5,5,iv,j) = mat(4,5,5,iv,j) + alfa * dtl(iv,j) * eps_e * buff(iv,j) * (-one)
                
        end if

!.... sponge term

        if (.not. calcd) then
        
        if (ispg .eq. 1) then
        
            mat(3,1,1,iv,j) = mat(3,1,1,iv,j) + alfa * dtl(iv,j) * spgl(j,iv)
            mat(3,2,2,iv,j) = mat(3,2,2,iv,j) + alfa * dtl(iv,j) * spgl(j,iv)
            mat(3,3,3,iv,j) = mat(3,3,3,iv,j) + alfa * dtl(iv,j) * spgl(j,iv)
            mat(3,4,4,iv,j) = mat(3,4,4,iv,j) + alfa * dtl(iv,j) * spgl(j,iv)
            mat(3,5,5,iv,j) = mat(3,5,5,iv,j) + alfa * dtl(iv,j) * spgl(j,iv)
                
        else if (ispg .ge. 2) then
        
            mat(3,1,1,iv,j) = mat(3,1,1,iv,j) + alfa * dtl(iv,j) * (spgl(j,iv) + spg2l(j,iv))
            mat(3,2,2,iv,j) = mat(3,2,2,iv,j) + alfa * dtl(iv,j) * (spgl(j,iv) + spg2l(j,iv))
            mat(3,3,3,iv,j) = mat(3,3,3,iv,j) + alfa * dtl(iv,j) * (spgl(j,iv) + spg2l(j,iv))
            mat(3,4,4,iv,j) = mat(3,4,4,iv,j) + alfa * dtl(iv,j) * (spgl(j,iv) + spg2l(j,iv))
            mat(3,5,5,iv,j) = mat(3,5,5,iv,j) + alfa * dtl(iv,j) * (spgl(j,iv) + spg2l(j,iv))
                
        end if

        end if
!=======================================================================================================!
!.... use higher-order tangent on first node off the far-field boundary
!=======================================================================================================!
        j = ny-1

        a1 = -alfa * gb1 * detainv
        a2 = -alfa * gb2 * detainv
        a3 = -alfa * gb3 * detainv
        a4 = -alfa * gb4 * detainv
        a5 = -alfa * gb5 * detainv

        do idof = 1, ndof
          do jdof = 1, ndof
              mat(5,idof,jdof,iv,j) = a5 * dtl(iv,j) * Bh(idof,jdof,iv,j)
              mat(1,idof,jdof,iv,j) = a4 * dtl(iv,j) * Bh(idof,jdof,iv,j)
              mat(2,idof,jdof,iv,j) = a3 * dtl(iv,j) * Bh(idof,jdof,iv,j)
              mat(3,idof,jdof,iv,j) = dtl(iv,j) * ( a2 * Bh(idof,jdof,iv,j) + &
                                      rcalcd * alfa * (  Dh(idof,jdof,iv,j) + &
                                                   im * Dhi(idof,jdof,iv,j) ) )
              mat(4,idof,jdof,iv,j) = a1 * dtl(iv,j) * Bh(idof,jdof,iv,j)
          end do
        end do
        
!.... \hat{B}_i term

        c1 = -im * alfa * gb1 * detainv
        c2 = -im * alfa * gb2 * detainv
        c3 = -im * alfa * gb3 * detainv
        c4 = -im * alfa * gb4 * detainv
        c5 = -im * alfa * gb5 * detainv

        mat(5,2,4,iv,j) = mat(5,2,4,iv,j) - c5 * dtl(iv,j) * Bhi(j,iv,3)
        mat(5,3,4,iv,j) = mat(5,3,4,iv,j) - c5 * dtl(iv,j) * Bhi(j,iv,4)
        mat(5,4,2,iv,j) = mat(5,4,2,iv,j) - c5 * dtl(iv,j) * Bhi(j,iv,3)
        mat(5,4,3,iv,j) = mat(5,4,3,iv,j) - c5 * dtl(iv,j) * Bhi(j,iv,4)

        mat(1,2,4,iv,j) = mat(1,2,4,iv,j) - c4 * dtl(iv,j) * Bhi(j,iv,3)
        mat(1,3,4,iv,j) = mat(1,3,4,iv,j) - c4 * dtl(iv,j) * Bhi(j,iv,4)
        mat(1,4,2,iv,j) = mat(1,4,2,iv,j) - c4 * dtl(iv,j) * Bhi(j,iv,3)
        mat(1,4,3,iv,j) = mat(1,4,3,iv,j) - c4 * dtl(iv,j) * Bhi(j,iv,4)

        mat(2,2,4,iv,j) = mat(2,2,4,iv,j) - c3 * dtl(iv,j) * Bhi(j,iv,3)
        mat(2,3,4,iv,j) = mat(2,3,4,iv,j) - c3 * dtl(iv,j) * Bhi(j,iv,4)
        mat(2,4,2,iv,j) = mat(2,4,2,iv,j) - c3 * dtl(iv,j) * Bhi(j,iv,3)
        mat(2,4,3,iv,j) = mat(2,4,3,iv,j) - c3 * dtl(iv,j) * Bhi(j,iv,4)

        mat(3,2,4,iv,j) = mat(3,2,4,iv,j) - c2 * dtl(iv,j) * Bhi(j,iv,3)
        mat(3,3,4,iv,j) = mat(3,3,4,iv,j) - c2 * dtl(iv,j) * Bhi(j,iv,4)
        mat(3,4,2,iv,j) = mat(3,4,2,iv,j) - c2 * dtl(iv,j) * Bhi(j,iv,3)
        mat(3,4,3,iv,j) = mat(3,4,3,iv,j) - c2 * dtl(iv,j) * Bhi(j,iv,4)

        mat(4,2,4,iv,j) = mat(4,2,4,iv,j) - c1 * dtl(iv,j) * Bhi(j,iv,3)
        mat(4,3,4,iv,j) = mat(4,3,4,iv,j) - c1 * dtl(iv,j) * Bhi(j,iv,4)
        mat(4,4,2,iv,j) = mat(4,4,2,iv,j) - c1 * dtl(iv,j) * Bhi(j,iv,3)
        mat(4,4,3,iv,j) = mat(4,4,3,iv,j) - c1 * dtl(iv,j) * Bhi(j,iv,4)
          
!.... \hat{V}_{\eta\eta} term

        a1 = alfa * db1 * detasinv
        a2 = alfa * db2 * detasinv
        a3 = alfa * db3 * detasinv
        a4 = alfa * db4 * detasinv
        a5 = alfa * db5 * detasinv

        mat(5,2,2,iv,j) = mat(5,2,2,iv,j) - a5 * dtl(iv,j) * Vh(1,iv,j)
        mat(5,2,3,iv,j) = mat(5,2,3,iv,j) - a5 * dtl(iv,j) * Vh(5,iv,j)
        mat(5,3,2,iv,j) = mat(5,3,2,iv,j) - a5 * dtl(iv,j) * Vh(6,iv,j)
        mat(5,3,3,iv,j) = mat(5,3,3,iv,j) - a5 * dtl(iv,j) * Vh(2,iv,j)
        mat(5,4,4,iv,j) = mat(5,4,4,iv,j) - a5 * dtl(iv,j) * Vh(3,iv,j)
        mat(5,5,5,iv,j) = mat(5,5,5,iv,j) - a5 * dtl(iv,j) * Vh(4,iv,j)

        mat(1,2,2,iv,j) = mat(1,2,2,iv,j) - a4 * dtl(iv,j) * Vh(1,iv,j)
        mat(1,2,3,iv,j) = mat(1,2,3,iv,j) - a4 * dtl(iv,j) * Vh(5,iv,j)
        mat(1,3,2,iv,j) = mat(1,3,2,iv,j) - a4 * dtl(iv,j) * Vh(6,iv,j)
        mat(1,3,3,iv,j) = mat(1,3,3,iv,j) - a4 * dtl(iv,j) * Vh(2,iv,j)
        mat(1,4,4,iv,j) = mat(1,4,4,iv,j) - a4 * dtl(iv,j) * Vh(3,iv,j)
        mat(1,5,5,iv,j) = mat(1,5,5,iv,j) - a4 * dtl(iv,j) * Vh(4,iv,j)

        mat(2,2,2,iv,j) = mat(2,2,2,iv,j) - a3 * dtl(iv,j) * Vh(1,iv,j)
        mat(2,2,3,iv,j) = mat(2,2,3,iv,j) - a3 * dtl(iv,j) * Vh(5,iv,j)
        mat(2,3,2,iv,j) = mat(2,3,2,iv,j) - a3 * dtl(iv,j) * Vh(6,iv,j)
        mat(2,3,3,iv,j) = mat(2,3,3,iv,j) - a3 * dtl(iv,j) * Vh(2,iv,j)
        mat(2,4,4,iv,j) = mat(2,4,4,iv,j) - a3 * dtl(iv,j) * Vh(3,iv,j)
        mat(2,5,5,iv,j) = mat(2,5,5,iv,j) - a3 * dtl(iv,j) * Vh(4,iv,j)

        mat(3,2,2,iv,j) = mat(3,2,2,iv,j) - a2 * dtl(iv,j) * Vh(1,iv,j)
        mat(3,2,3,iv,j) = mat(3,2,3,iv,j) - a2 * dtl(iv,j) * Vh(5,iv,j)
        mat(3,3,2,iv,j) = mat(3,3,2,iv,j) - a2 * dtl(iv,j) * Vh(6,iv,j)
        mat(3,3,3,iv,j) = mat(3,3,3,iv,j) - a2 * dtl(iv,j) * Vh(2,iv,j)
        mat(3,4,4,iv,j) = mat(3,4,4,iv,j) - a2 * dtl(iv,j) * Vh(3,iv,j)
        mat(3,5,5,iv,j) = mat(3,5,5,iv,j) - a2 * dtl(iv,j) * Vh(4,iv,j)

        mat(4,2,2,iv,j) = mat(4,2,2,iv,j) - a1 * dtl(iv,j) * Vh(1,iv,j)
        mat(4,2,3,iv,j) = mat(4,2,3,iv,j) - a1 * dtl(iv,j) * Vh(5,iv,j)
        mat(4,3,2,iv,j) = mat(4,3,2,iv,j) - a1 * dtl(iv,j) * Vh(6,iv,j)
        mat(4,3,3,iv,j) = mat(4,3,3,iv,j) - a1 * dtl(iv,j) * Vh(2,iv,j)
        mat(4,4,4,iv,j) = mat(4,4,4,iv,j) - a1 * dtl(iv,j) * Vh(3,iv,j)
        mat(4,5,5,iv,j) = mat(4,5,5,iv,j) - a1 * dtl(iv,j) * Vh(4,iv,j)

!.... I term

        mat(3,1,1,iv,j) = mat(3,1,1,iv,j) + one
        mat(3,2,2,iv,j) = mat(3,2,2,iv,j) + one
        mat(3,3,3,iv,j) = mat(3,3,3,iv,j) + one
        mat(3,4,4,iv,j) = mat(3,4,4,iv,j) + one
        mat(3,5,5,iv,j) = mat(3,5,5,iv,j) + one

!.... implicit damping term

        if (eps_e .ne. zero) then

          mat(2,1,1,iv,j) = mat(2,1,1,iv,j) + alfa * dtl(iv,j) * eps_e * buff(iv,j) * (-one)
          mat(2,2,2,iv,j) = mat(2,2,2,iv,j) + alfa * dtl(iv,j) * eps_e * buff(iv,j) * (-one)
          mat(2,3,3,iv,j) = mat(2,3,3,iv,j) + alfa * dtl(iv,j) * eps_e * buff(iv,j) * (-one)
          mat(2,4,4,iv,j) = mat(2,4,4,iv,j) + alfa * dtl(iv,j) * eps_e * buff(iv,j) * (-one)
          mat(2,5,5,iv,j) = mat(2,5,5,iv,j) + alfa * dtl(iv,j) * eps_e * buff(iv,j) * (-one)
  
          mat(3,1,1,iv,j) = mat(3,1,1,iv,j) + alfa * dtl(iv,j) * eps_e * buff(iv,j) * two
          mat(3,2,2,iv,j) = mat(3,2,2,iv,j) + alfa * dtl(iv,j) * eps_e * buff(iv,j) * two
          mat(3,3,3,iv,j) = mat(3,3,3,iv,j) + alfa * dtl(iv,j) * eps_e * buff(iv,j) * two
          mat(3,4,4,iv,j) = mat(3,4,4,iv,j) + alfa * dtl(iv,j) * eps_e * buff(iv,j) * two
          mat(3,5,5,iv,j) = mat(3,5,5,iv,j) + alfa * dtl(iv,j) * eps_e * buff(iv,j) * two
  
          mat(4,1,1,iv,j) = mat(4,1,1,iv,j) + alfa * dtl(iv,j) * eps_e * buff(iv,j) * (-one)
          mat(4,2,2,iv,j) = mat(4,2,2,iv,j) + alfa * dtl(iv,j) * eps_e * buff(iv,j) * (-one)
          mat(4,3,3,iv,j) = mat(4,3,3,iv,j) + alfa * dtl(iv,j) * eps_e * buff(iv,j) * (-one)
          mat(4,4,4,iv,j) = mat(4,4,4,iv,j) + alfa * dtl(iv,j) * eps_e * buff(iv,j) * (-one)
          mat(4,5,5,iv,j) = mat(4,5,5,iv,j) + alfa * dtl(iv,j) * eps_e * buff(iv,j) * (-one)
        
        end if

!.... sponge term

        if (.not. calcd) then

        if (ispg .eq. 1) then

            mat(3,1,1,iv,j) = mat(3,1,1,iv,j) + alfa * dtl(iv,j) * spgl(j,iv)
            mat(3,2,2,iv,j) = mat(3,2,2,iv,j) + alfa * dtl(iv,j) * spgl(j,iv)
            mat(3,3,3,iv,j) = mat(3,3,3,iv,j) + alfa * dtl(iv,j) * spgl(j,iv)
            mat(3,4,4,iv,j) = mat(3,4,4,iv,j) + alfa * dtl(iv,j) * spgl(j,iv)
            mat(3,5,5,iv,j) = mat(3,5,5,iv,j) + alfa * dtl(iv,j) * spgl(j,iv)

        else if (ispg .ge. 2) then

            mat(3,1,1,iv,j) = mat(3,1,1,iv,j) + alfa * dtl(iv,j) * (spgl(j,iv) + spg2l(j,iv))
            mat(3,2,2,iv,j) = mat(3,2,2,iv,j) + alfa * dtl(iv,j) * (spgl(j,iv) + spg2l(j,iv))
            mat(3,3,3,iv,j) = mat(3,3,3,iv,j) + alfa * dtl(iv,j) * (spgl(j,iv) + spg2l(j,iv))
            mat(3,4,4,iv,j) = mat(3,4,4,iv,j) + alfa * dtl(iv,j) * (spgl(j,iv) + spg2l(j,iv))
            mat(3,5,5,iv,j) = mat(3,5,5,iv,j) + alfa * dtl(iv,j) * (spgl(j,iv) + spg2l(j,iv))

        end if

        end if
!=======================================================================================================!
!.... use higher-order tangent on the far-field boundary
!=======================================================================================================!
        j = ny

        a1 = -alfa * gc1 * detainv
        a2 = -alfa * gc2 * detainv
        a3 = -alfa * gc3 * detainv
        a4 = -alfa * gc4 * detainv
        a5 = -alfa * gc5 * detainv

        do idof = 1, ndof
          do jdof = 1, ndof
              mat(4,idof,jdof,iv,j) = a5 * dtl(iv,j) * Bh(idof,jdof,iv,j)
              mat(5,idof,jdof,iv,j) = a4 * dtl(iv,j) * Bh(idof,jdof,iv,j)
              mat(1,idof,jdof,iv,j) = a3 * dtl(iv,j) * Bh(idof,jdof,iv,j)
              mat(2,idof,jdof,iv,j) = a2 * dtl(iv,j) * Bh(idof,jdof,iv,j)
              mat(3,idof,jdof,iv,j) = dtl(iv,j) * ( a1 * Bh(idof,jdof,iv,j) + &
                                      rcalcd * alfa * (  Dh(idof,jdof,iv,j) + &
                                                   im * Dhi(idof,jdof,iv,j) ) )
          end do
        end do
        
!.... \hat{B}_i term

        c1 = -im * alfa * gc1 * detainv
        c2 = -im * alfa * gc2 * detainv
        c3 = -im * alfa * gc3 * detainv
        c4 = -im * alfa * gc4 * detainv
        c5 = -im * alfa * gc5 * detainv

        mat(4,2,4,iv,j) = mat(4,2,4,iv,j) - c5 * dtl(iv,j) * Bhi(j,iv,3)
        mat(4,3,4,iv,j) = mat(4,3,4,iv,j) - c5 * dtl(iv,j) * Bhi(j,iv,4)
        mat(4,4,2,iv,j) = mat(4,4,2,iv,j) - c5 * dtl(iv,j) * Bhi(j,iv,3)
        mat(4,4,3,iv,j) = mat(4,4,3,iv,j) - c5 * dtl(iv,j) * Bhi(j,iv,4)

        mat(5,2,4,iv,j) = mat(5,2,4,iv,j) - c4 * dtl(iv,j) * Bhi(j,iv,3)
        mat(5,3,4,iv,j) = mat(5,3,4,iv,j) - c4 * dtl(iv,j) * Bhi(j,iv,4)
        mat(5,4,2,iv,j) = mat(5,4,2,iv,j) - c4 * dtl(iv,j) * Bhi(j,iv,3)
        mat(5,4,3,iv,j) = mat(5,4,3,iv,j) - c4 * dtl(iv,j) * Bhi(j,iv,4)

        mat(1,2,4,iv,j) = mat(1,2,4,iv,j) - c3 * dtl(iv,j) * Bhi(j,iv,3)
        mat(1,3,4,iv,j) = mat(1,3,4,iv,j) - c3 * dtl(iv,j) * Bhi(j,iv,4)
        mat(1,4,2,iv,j) = mat(1,4,2,iv,j) - c3 * dtl(iv,j) * Bhi(j,iv,3)
        mat(1,4,3,iv,j) = mat(1,4,3,iv,j) - c3 * dtl(iv,j) * Bhi(j,iv,4)

        mat(2,2,4,iv,j) = mat(2,2,4,iv,j) - c2 * dtl(iv,j) * Bhi(j,iv,3)
        mat(2,3,4,iv,j) = mat(2,3,4,iv,j) - c2 * dtl(iv,j) * Bhi(j,iv,4)
        mat(2,4,2,iv,j) = mat(2,4,2,iv,j) - c2 * dtl(iv,j) * Bhi(j,iv,3)
        mat(2,4,3,iv,j) = mat(2,4,3,iv,j) - c2 * dtl(iv,j) * Bhi(j,iv,4)

        mat(3,2,4,iv,j) = mat(3,2,4,iv,j) - c1 * dtl(iv,j) * Bhi(j,iv,3)
        mat(3,3,4,iv,j) = mat(3,3,4,iv,j) - c1 * dtl(iv,j) * Bhi(j,iv,4)
        mat(3,4,2,iv,j) = mat(3,4,2,iv,j) - c1 * dtl(iv,j) * Bhi(j,iv,3)
        mat(3,4,3,iv,j) = mat(3,4,3,iv,j) - c1 * dtl(iv,j) * Bhi(j,iv,4)
          
!.... \hat{V}_{\eta\eta} term

        a1 = alfa * dd1 * detasinv
        a2 = alfa * dd2 * detasinv
        a3 = alfa * dd3 * detasinv
        a4 = alfa * dd4 * detasinv
        a5 = alfa * dd5 * detasinv

        mat(4,2,2,iv,j) = mat(4,2,2,iv,j) - a5 * dtl(iv,j) * Vh(1,iv,j)
        mat(4,2,3,iv,j) = mat(4,2,3,iv,j) - a5 * dtl(iv,j) * Vh(5,iv,j)
        mat(4,3,2,iv,j) = mat(4,3,2,iv,j) - a5 * dtl(iv,j) * Vh(6,iv,j)
        mat(4,3,3,iv,j) = mat(4,3,3,iv,j) - a5 * dtl(iv,j) * Vh(2,iv,j)
        mat(4,4,4,iv,j) = mat(4,4,4,iv,j) - a5 * dtl(iv,j) * Vh(3,iv,j)
        mat(4,5,5,iv,j) = mat(4,5,5,iv,j) - a5 * dtl(iv,j) * Vh(4,iv,j)

        mat(5,2,2,iv,j) = mat(5,2,2,iv,j) - a4 * dtl(iv,j) * Vh(1,iv,j)
        mat(5,2,3,iv,j) = mat(5,2,3,iv,j) - a4 * dtl(iv,j) * Vh(5,iv,j)
        mat(5,3,2,iv,j) = mat(5,3,2,iv,j) - a4 * dtl(iv,j) * Vh(6,iv,j)
        mat(5,3,3,iv,j) = mat(5,3,3,iv,j) - a4 * dtl(iv,j) * Vh(2,iv,j)
        mat(5,4,4,iv,j) = mat(5,4,4,iv,j) - a4 * dtl(iv,j) * Vh(3,iv,j)
        mat(5,5,5,iv,j) = mat(5,5,5,iv,j) - a4 * dtl(iv,j) * Vh(4,iv,j)

        mat(1,2,2,iv,j) = mat(1,2,2,iv,j) - a3 * dtl(iv,j) * Vh(1,iv,j)
        mat(1,2,3,iv,j) = mat(1,2,3,iv,j) - a3 * dtl(iv,j) * Vh(5,iv,j)
        mat(1,3,2,iv,j) = mat(1,3,2,iv,j) - a3 * dtl(iv,j) * Vh(6,iv,j)
        mat(1,3,3,iv,j) = mat(1,3,3,iv,j) - a3 * dtl(iv,j) * Vh(2,iv,j)
        mat(1,4,4,iv,j) = mat(1,4,4,iv,j) - a3 * dtl(iv,j) * Vh(3,iv,j)
        mat(1,5,5,iv,j) = mat(1,5,5,iv,j) - a3 * dtl(iv,j) * Vh(4,iv,j)

        mat(2,2,2,iv,j) = mat(2,2,2,iv,j) - a2 * dtl(iv,j) * Vh(1,iv,j)
        mat(2,2,3,iv,j) = mat(2,2,3,iv,j) - a2 * dtl(iv,j) * Vh(5,iv,j)
        mat(2,3,2,iv,j) = mat(2,3,2,iv,j) - a2 * dtl(iv,j) * Vh(6,iv,j)
        mat(2,3,3,iv,j) = mat(2,3,3,iv,j) - a2 * dtl(iv,j) * Vh(2,iv,j)
        mat(2,4,4,iv,j) = mat(2,4,4,iv,j) - a2 * dtl(iv,j) * Vh(3,iv,j)
        mat(2,5,5,iv,j) = mat(2,5,5,iv,j) - a2 * dtl(iv,j) * Vh(4,iv,j)

        mat(3,2,2,iv,j) = mat(3,2,2,iv,j) - a1 * dtl(iv,j) * Vh(1,iv,j)
        mat(3,2,3,iv,j) = mat(3,2,3,iv,j) - a1 * dtl(iv,j) * Vh(5,iv,j)
        mat(3,3,2,iv,j) = mat(3,3,2,iv,j) - a1 * dtl(iv,j) * Vh(6,iv,j)
        mat(3,3,3,iv,j) = mat(3,3,3,iv,j) - a1 * dtl(iv,j) * Vh(2,iv,j)
        mat(3,4,4,iv,j) = mat(3,4,4,iv,j) - a1 * dtl(iv,j) * Vh(3,iv,j)
        mat(3,5,5,iv,j) = mat(3,5,5,iv,j) - a1 * dtl(iv,j) * Vh(4,iv,j)

!.... I term

        mat(3,1,1,iv,j) = mat(3,1,1,iv,j) + one
        mat(3,2,2,iv,j) = mat(3,2,2,iv,j) + one
        mat(3,3,3,iv,j) = mat(3,3,3,iv,j) + one
        mat(3,4,4,iv,j) = mat(3,4,4,iv,j) + one
        mat(3,5,5,iv,j) = mat(3,5,5,iv,j) + one

!.... sponge term

        if (.not. calcd) then

        if (ispg .eq. 1) then
        
            mat(3,1,1,iv,j) = mat(3,1,1,iv,j) + alfa * dtl(iv,j) * spgl(j,iv)
            mat(3,2,2,iv,j) = mat(3,2,2,iv,j) + alfa * dtl(iv,j) * spgl(j,iv)
            mat(3,3,3,iv,j) = mat(3,3,3,iv,j) + alfa * dtl(iv,j) * spgl(j,iv)
            mat(3,4,4,iv,j) = mat(3,4,4,iv,j) + alfa * dtl(iv,j) * spgl(j,iv)
            mat(3,5,5,iv,j) = mat(3,5,5,iv,j) + alfa * dtl(iv,j) * spgl(j,iv)
      
        else if (ispg .ge. 2) then
        
            mat(3,1,1,iv,j) = mat(3,1,1,iv,j) + alfa * dtl(iv,j) * (spgl(j,iv) + spg2l(j,iv))
            mat(3,2,2,iv,j) = mat(3,2,2,iv,j) + alfa * dtl(iv,j) * (spgl(j,iv) + spg2l(j,iv))
            mat(3,3,3,iv,j) = mat(3,3,3,iv,j) + alfa * dtl(iv,j) * (spgl(j,iv) + spg2l(j,iv))
            mat(3,4,4,iv,j) = mat(3,4,4,iv,j) + alfa * dtl(iv,j) * (spgl(j,iv) + spg2l(j,iv))
            mat(3,5,5,iv,j) = mat(3,5,5,iv,j) + alfa * dtl(iv,j) * (spgl(j,iv) + spg2l(j,iv))
                
        end if

        end if

        end do

        return
        end
