!=======================================================================================================!
        subroutine lhsbt1f3D( mat, Ah, Dh, Dhi, Vh, Ahi, spgl, spg2l, dtl, calcd )
!  
!  Correct the LHS for boundary treatment in the \xi direction.
!
!  This version supports the 4th order LHS and three-dimensional flow
!
!  Revised: 4-22-96
!=======================================================================================================!
        use global
        use stencil
        use buff_mod
        implicit none
        
        complex :: mat(ny,nx,ndof,ndof,5)
        real    :: Ah(ny,nx,ndof,ndof), Dh(ny,nx,ndof,ndof)
        real    :: Dhi(ny,nx,ndof,ndof), Ahi(ny,nx,6)
        real    :: spgl(ny,nx), spg2l(ny,nx), Vh(ny,nx,6), dtl(ny,nx)
        logical :: calcd

        real    :: a1, a2, a3, a4, a5
        complex :: c1, c2, c3, c4, c5

        real    :: dxiinv, dxisinv

        integer :: iv, i, idof, jdof
        real    :: rcalcd
!=======================================================================================================!
        rcalcd = zero
        if (calcd) rcalcd = one

        if (xper) return

        dxiinv  = one / dxi
        dxisinv = one / dxi**2

        if (lsym) then           !  apply a symmetry condition
        
          call lhsbs1f3D( mat, Ah, Dh, Dhi, Vh, Ahi, spgl, spg2l, dtl, calcd, 1 )       

        else
                
!=======================================================================================================!
!.... use higher-order tangent on left node
!=======================================================================================================!
        i = 1

        a1 = alfa * gc1 * dxiinv
        a2 = alfa * gc2 * dxiinv
        a3 = alfa * gc3 * dxiinv
        a4 = alfa * gc4 * dxiinv
        a5 = alfa * gc5 * dxiinv

        do idof = 1, ndof
          do jdof = 1, ndof
            do iv = 1, ny
              mat(iv,i,idof,jdof,3) = dtl(iv,i) * ( a1 * Ah(iv,i,idof,jdof) + &
                                      rcalcd * alfa * (  Dh(iv,i,idof,jdof) + &
                                                   im * Dhi(iv,i,idof,jdof) ) )
              mat(iv,i,idof,jdof,4) = a2 * dtl(iv,i) * Ah(iv,i,idof,jdof)
              mat(iv,i,idof,jdof,5) = a3 * dtl(iv,i) * Ah(iv,i,idof,jdof)
              mat(iv,i,idof,jdof,1) = a4 * dtl(iv,i) * Ah(iv,i,idof,jdof)
              mat(iv,i,idof,jdof,2) = a5 * dtl(iv,i) * Ah(iv,i,idof,jdof)
            end do
          end do
        end do
        
!.... \hat{A}_i term

        c1 = im * alfa * gc1 * dxiinv
        c2 = im * alfa * gc2 * dxiinv
        c3 = im * alfa * gc3 * dxiinv
        c4 = im * alfa * gc4 * dxiinv
        c5 = im * alfa * gc5 * dxiinv

        do iv = 1, ny

        mat(iv,i,2,4,3) = mat(iv,i,2,4,3) - c1 * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,3,4,3) = mat(iv,i,3,4,3) - c1 * dtl(iv,i) * Ahi(iv,i,2)
        mat(iv,i,4,2,3) = mat(iv,i,4,2,3) - c1 * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,4,3,3) = mat(iv,i,4,3,3) - c1 * dtl(iv,i) * Ahi(iv,i,2)

        mat(iv,i,2,4,4) = mat(iv,i,2,4,4) - c2 * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,3,4,4) = mat(iv,i,3,4,4) - c2 * dtl(iv,i) * Ahi(iv,i,2)
        mat(iv,i,4,2,4) = mat(iv,i,4,2,4) - c2 * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,4,3,4) = mat(iv,i,4,3,4) - c2 * dtl(iv,i) * Ahi(iv,i,2)

        mat(iv,i,2,4,5) = mat(iv,i,2,4,5) - c3 * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,3,4,5) = mat(iv,i,3,4,5) - c3 * dtl(iv,i) * Ahi(iv,i,2)
        mat(iv,i,4,2,5) = mat(iv,i,4,2,5) - c3 * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,4,3,5) = mat(iv,i,4,3,5) - c3 * dtl(iv,i) * Ahi(iv,i,2)

        mat(iv,i,2,4,1) = mat(iv,i,2,4,1) - c4 * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,3,4,1) = mat(iv,i,3,4,1) - c4 * dtl(iv,i) * Ahi(iv,i,2)
        mat(iv,i,4,2,1) = mat(iv,i,4,2,1) - c4 * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,4,3,1) = mat(iv,i,4,3,1) - c4 * dtl(iv,i) * Ahi(iv,i,2)

        mat(iv,i,2,4,2) = mat(iv,i,2,4,2) - c5 * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,3,4,2) = mat(iv,i,3,4,2) - c5 * dtl(iv,i) * Ahi(iv,i,2)
        mat(iv,i,4,2,2) = mat(iv,i,4,2,2) - c5 * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,4,3,2) = mat(iv,i,4,3,2) - c5 * dtl(iv,i) * Ahi(iv,i,2)
          
        end do

!.... \hat{V}_{\xi\xi} term

        a1 = alfa * dd1 * dxisinv
        a2 = alfa * dd2 * dxisinv
        a3 = alfa * dd3 * dxisinv
        a4 = alfa * dd4 * dxisinv
        a5 = alfa * dd5 * dxisinv

        do iv = 1, ny
        
        mat(iv,i,2,2,3) = mat(iv,i,2,2,3) - a1 * dtl(iv,i) * Vh(iv,i,1)
        mat(iv,i,2,3,3) = mat(iv,i,2,3,3) - a1 * dtl(iv,i) * Vh(iv,i,5)
        mat(iv,i,3,2,3) = mat(iv,i,3,2,3) - a1 * dtl(iv,i) * Vh(iv,i,6)
        mat(iv,i,3,3,3) = mat(iv,i,3,3,3) - a1 * dtl(iv,i) * Vh(iv,i,2)
        mat(iv,i,4,4,3) = mat(iv,i,4,4,3) - a1 * dtl(iv,i) * Vh(iv,i,3)
        mat(iv,i,5,5,3) = mat(iv,i,5,5,3) - a1 * dtl(iv,i) * Vh(iv,i,4)

        mat(iv,i,2,2,4) = mat(iv,i,2,2,4) - a2 * dtl(iv,i) * Vh(iv,i,1)
        mat(iv,i,2,3,4) = mat(iv,i,2,3,4) - a2 * dtl(iv,i) * Vh(iv,i,5)
        mat(iv,i,3,2,4) = mat(iv,i,3,2,4) - a2 * dtl(iv,i) * Vh(iv,i,6)
        mat(iv,i,3,3,4) = mat(iv,i,3,3,4) - a2 * dtl(iv,i) * Vh(iv,i,2)
        mat(iv,i,4,4,4) = mat(iv,i,4,4,4) - a2 * dtl(iv,i) * Vh(iv,i,3)
        mat(iv,i,5,5,4) = mat(iv,i,5,5,4) - a2 * dtl(iv,i) * Vh(iv,i,4)

        mat(iv,i,2,2,5) = mat(iv,i,2,2,5) - a3 * dtl(iv,i) * Vh(iv,i,1)
        mat(iv,i,2,3,5) = mat(iv,i,2,3,5) - a3 * dtl(iv,i) * Vh(iv,i,5)
        mat(iv,i,3,2,5) = mat(iv,i,3,2,5) - a3 * dtl(iv,i) * Vh(iv,i,6)
        mat(iv,i,3,3,5) = mat(iv,i,3,3,5) - a3 * dtl(iv,i) * Vh(iv,i,2)
        mat(iv,i,4,4,5) = mat(iv,i,4,4,5) - a3 * dtl(iv,i) * Vh(iv,i,3)
        mat(iv,i,5,5,5) = mat(iv,i,5,5,5) - a3 * dtl(iv,i) * Vh(iv,i,4)

        mat(iv,i,2,2,1) = mat(iv,i,2,2,1) - a4 * dtl(iv,i) * Vh(iv,i,1)
        mat(iv,i,2,3,1) = mat(iv,i,2,3,1) - a4 * dtl(iv,i) * Vh(iv,i,5)
        mat(iv,i,3,2,1) = mat(iv,i,3,2,1) - a4 * dtl(iv,i) * Vh(iv,i,6)
        mat(iv,i,3,3,1) = mat(iv,i,3,3,1) - a4 * dtl(iv,i) * Vh(iv,i,2)
        mat(iv,i,4,4,1) = mat(iv,i,4,4,1) - a4 * dtl(iv,i) * Vh(iv,i,3)
        mat(iv,i,5,5,1) = mat(iv,i,5,5,1) - a4 * dtl(iv,i) * Vh(iv,i,4)

        mat(iv,i,2,2,2) = mat(iv,i,2,2,2) - a5 * dtl(iv,i) * Vh(iv,i,1)
        mat(iv,i,2,3,2) = mat(iv,i,2,3,2) - a5 * dtl(iv,i) * Vh(iv,i,5)
        mat(iv,i,3,2,2) = mat(iv,i,3,2,2) - a5 * dtl(iv,i) * Vh(iv,i,6)
        mat(iv,i,3,3,2) = mat(iv,i,3,3,2) - a5 * dtl(iv,i) * Vh(iv,i,2)
        mat(iv,i,4,4,2) = mat(iv,i,4,4,2) - a5 * dtl(iv,i) * Vh(iv,i,3)
        mat(iv,i,5,5,2) = mat(iv,i,5,5,2) - a5 * dtl(iv,i) * Vh(iv,i,4)

!.... I term

        mat(iv,i,1,1,3) = mat(iv,i,1,1,3) + one
        mat(iv,i,2,2,3) = mat(iv,i,2,2,3) + one
        mat(iv,i,3,3,3) = mat(iv,i,3,3,3) + one
        mat(iv,i,4,4,3) = mat(iv,i,4,4,3) + one
        mat(iv,i,5,5,3) = mat(iv,i,5,5,3) + one

        end do

!.... sponge term

        if (.not. calcd) then

        if (ispg .eq. 1) then
        
          do iv = 1, ny
            mat(iv,i,1,1,3) = mat(iv,i,1,1,3) + alfa * dtl(iv,i) * spgl(iv,i)
            mat(iv,i,2,2,3) = mat(iv,i,2,2,3) + alfa * dtl(iv,i) * spgl(iv,i)
            mat(iv,i,3,3,3) = mat(iv,i,3,3,3) + alfa * dtl(iv,i) * spgl(iv,i)
            mat(iv,i,4,4,3) = mat(iv,i,4,4,3) + alfa * dtl(iv,i) * spgl(iv,i)
            mat(iv,i,5,5,3) = mat(iv,i,5,5,3) + alfa * dtl(iv,i) * spgl(iv,i)
          end do

        else if (ispg .ge. 2) then
        
          do iv = 1, ny
            mat(iv,i,1,1,3) = mat(iv,i,1,1,3) + alfa * dtl(iv,i) * (spgl(iv,i) + spg2l(iv,i))
            mat(iv,i,2,2,3) = mat(iv,i,2,2,3) + alfa * dtl(iv,i) * (spgl(iv,i) + spg2l(iv,i))
            mat(iv,i,3,3,3) = mat(iv,i,3,3,3) + alfa * dtl(iv,i) * (spgl(iv,i) + spg2l(iv,i))
            mat(iv,i,4,4,3) = mat(iv,i,4,4,3) + alfa * dtl(iv,i) * (spgl(iv,i) + spg2l(iv,i))
            mat(iv,i,5,5,3) = mat(iv,i,5,5,3) + alfa * dtl(iv,i) * (spgl(iv,i) + spg2l(iv,i))
          end do
        
        end if

        end if
!=======================================================================================================!
!.... use higher-order tangent on first node off the left boundary
!=======================================================================================================!
        i = 2

        a1 = alfa * gb1 * dxiinv
        a2 = alfa * gb2 * dxiinv
        a3 = alfa * gb3 * dxiinv
        a4 = alfa * gb4 * dxiinv
        a5 = alfa * gb5 * dxiinv

        do idof = 1, ndof
          do jdof = 1, ndof
            do iv = 1, ny
              mat(iv,i,idof,jdof,2) = a1 * dtl(iv,i) * Ah(iv,i,idof,jdof)
              mat(iv,i,idof,jdof,3) = dtl(iv,i) * ( a2 * Ah(iv,i,idof,jdof) + &
                                      rcalcd * alfa * (  Dh(iv,i,idof,jdof) + &
                                                   im * Dhi(iv,i,idof,jdof) ) )
              mat(iv,i,idof,jdof,4) = a3 * dtl(iv,i) * Ah(iv,i,idof,jdof)
              mat(iv,i,idof,jdof,5) = a4 * dtl(iv,i) * Ah(iv,i,idof,jdof)
              mat(iv,i,idof,jdof,1) = a5 * dtl(iv,i) * Ah(iv,i,idof,jdof)
            end do
          end do
        end do
        
!.... \hat{A}_i term

        c1 = im * alfa * gb1 * dxiinv
        c2 = im * alfa * gb2 * dxiinv
        c3 = im * alfa * gb3 * dxiinv
        c4 = im * alfa * gb4 * dxiinv
        c5 = im * alfa * gb5 * dxiinv

        do iv = 1, ny

        mat(iv,i,2,4,2) = mat(iv,i,2,4,2) - c1 * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,3,4,2) = mat(iv,i,3,4,2) - c1 * dtl(iv,i) * Ahi(iv,i,2)
        mat(iv,i,4,2,2) = mat(iv,i,4,2,2) - c1 * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,4,3,2) = mat(iv,i,4,3,2) - c1 * dtl(iv,i) * Ahi(iv,i,2)

        mat(iv,i,2,4,3) = mat(iv,i,2,4,3) - c2 * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,3,4,3) = mat(iv,i,3,4,3) - c2 * dtl(iv,i) * Ahi(iv,i,2)
        mat(iv,i,4,2,3) = mat(iv,i,4,2,3) - c2 * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,4,3,3) = mat(iv,i,4,3,3) - c2 * dtl(iv,i) * Ahi(iv,i,2)

        mat(iv,i,2,4,4) = mat(iv,i,2,4,4) - c3 * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,3,4,4) = mat(iv,i,3,4,4) - c3 * dtl(iv,i) * Ahi(iv,i,2)
        mat(iv,i,4,2,4) = mat(iv,i,4,2,4) - c3 * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,4,3,4) = mat(iv,i,4,3,4) - c3 * dtl(iv,i) * Ahi(iv,i,2)

        mat(iv,i,2,4,5) = mat(iv,i,2,4,5) - c4 * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,3,4,5) = mat(iv,i,3,4,5) - c4 * dtl(iv,i) * Ahi(iv,i,2)
        mat(iv,i,4,2,5) = mat(iv,i,4,2,5) - c4 * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,4,3,5) = mat(iv,i,4,3,5) - c4 * dtl(iv,i) * Ahi(iv,i,2)

        mat(iv,i,2,4,1) = mat(iv,i,2,4,1) - c5 * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,3,4,1) = mat(iv,i,3,4,1) - c5 * dtl(iv,i) * Ahi(iv,i,2)
        mat(iv,i,4,2,1) = mat(iv,i,4,2,1) - c5 * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,4,3,1) = mat(iv,i,4,3,1) - c5 * dtl(iv,i) * Ahi(iv,i,2)
          
        end do

!.... \hat{V}_{\xi\xi} term

        a1 = alfa * db1 * dxisinv
        a2 = alfa * db2 * dxisinv
        a3 = alfa * db3 * dxisinv
        a4 = alfa * db4 * dxisinv
        a5 = alfa * db5 * dxisinv

        do iv = 1, ny
        
        mat(iv,i,2,2,2) = mat(iv,i,2,2,2) - a1 * dtl(iv,i) * Vh(iv,i,1)
        mat(iv,i,2,3,2) = mat(iv,i,2,3,2) - a1 * dtl(iv,i) * Vh(iv,i,5)
        mat(iv,i,3,2,2) = mat(iv,i,3,2,2) - a1 * dtl(iv,i) * Vh(iv,i,6)
        mat(iv,i,3,3,2) = mat(iv,i,3,3,2) - a1 * dtl(iv,i) * Vh(iv,i,2)
        mat(iv,i,4,4,2) = mat(iv,i,4,4,2) - a1 * dtl(iv,i) * Vh(iv,i,3)
        mat(iv,i,5,5,2) = mat(iv,i,5,5,2) - a1 * dtl(iv,i) * Vh(iv,i,4)

        mat(iv,i,2,2,3) = mat(iv,i,2,2,3) - a2 * dtl(iv,i) * Vh(iv,i,1)
        mat(iv,i,2,3,3) = mat(iv,i,2,3,3) - a2 * dtl(iv,i) * Vh(iv,i,5)
        mat(iv,i,3,2,3) = mat(iv,i,3,2,3) - a2 * dtl(iv,i) * Vh(iv,i,6)
        mat(iv,i,3,3,3) = mat(iv,i,3,3,3) - a2 * dtl(iv,i) * Vh(iv,i,2)
        mat(iv,i,4,4,3) = mat(iv,i,4,4,3) - a2 * dtl(iv,i) * Vh(iv,i,3)
        mat(iv,i,5,5,3) = mat(iv,i,5,5,3) - a2 * dtl(iv,i) * Vh(iv,i,4)

        mat(iv,i,2,2,4) = mat(iv,i,2,2,4) - a3 * dtl(iv,i) * Vh(iv,i,1)
        mat(iv,i,2,3,4) = mat(iv,i,2,3,4) - a3 * dtl(iv,i) * Vh(iv,i,5)
        mat(iv,i,3,2,4) = mat(iv,i,3,2,4) - a3 * dtl(iv,i) * Vh(iv,i,6)
        mat(iv,i,3,3,4) = mat(iv,i,3,3,4) - a3 * dtl(iv,i) * Vh(iv,i,2)
        mat(iv,i,4,4,4) = mat(iv,i,4,4,4) - a3 * dtl(iv,i) * Vh(iv,i,3)
        mat(iv,i,5,5,4) = mat(iv,i,5,5,4) - a3 * dtl(iv,i) * Vh(iv,i,4)

        mat(iv,i,2,2,5) = mat(iv,i,2,2,5) - a4 * dtl(iv,i) * Vh(iv,i,1)
        mat(iv,i,2,3,5) = mat(iv,i,2,3,5) - a4 * dtl(iv,i) * Vh(iv,i,5)
        mat(iv,i,3,2,5) = mat(iv,i,3,2,5) - a4 * dtl(iv,i) * Vh(iv,i,6)
        mat(iv,i,3,3,5) = mat(iv,i,3,3,5) - a4 * dtl(iv,i) * Vh(iv,i,2)
        mat(iv,i,4,4,5) = mat(iv,i,4,4,5) - a4 * dtl(iv,i) * Vh(iv,i,3)
        mat(iv,i,5,5,5) = mat(iv,i,5,5,5) - a4 * dtl(iv,i) * Vh(iv,i,4)

        mat(iv,i,2,2,1) = mat(iv,i,2,2,1) - a5 * dtl(iv,i) * Vh(iv,i,1)
        mat(iv,i,2,3,1) = mat(iv,i,2,3,1) - a5 * dtl(iv,i) * Vh(iv,i,5)
        mat(iv,i,3,2,1) = mat(iv,i,3,2,1) - a5 * dtl(iv,i) * Vh(iv,i,6)
        mat(iv,i,3,3,1) = mat(iv,i,3,3,1) - a5 * dtl(iv,i) * Vh(iv,i,2)
        mat(iv,i,4,4,1) = mat(iv,i,4,4,1) - a5 * dtl(iv,i) * Vh(iv,i,3)
        mat(iv,i,5,5,1) = mat(iv,i,5,5,1) - a5 * dtl(iv,i) * Vh(iv,i,4)

!.... I term

        mat(iv,i,1,1,3) = mat(iv,i,1,1,3) + one
        mat(iv,i,2,2,3) = mat(iv,i,2,2,3) + one
        mat(iv,i,3,3,3) = mat(iv,i,3,3,3) + one
        mat(iv,i,4,4,3) = mat(iv,i,4,4,3) + one
        mat(iv,i,5,5,3) = mat(iv,i,5,5,3) + one

        end do

!.... implicit damping term

        if (eps_e .ne. zero) then

        do iv = 1, ny
        
        mat(iv,i,1,1,2) = mat(iv,i,1,1,2) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * (-one)
        mat(iv,i,2,2,2) = mat(iv,i,2,2,2) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * (-one)
        mat(iv,i,3,3,2) = mat(iv,i,3,3,2) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * (-one)
        mat(iv,i,4,4,2) = mat(iv,i,4,4,2) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * (-one)
        mat(iv,i,5,5,2) = mat(iv,i,5,5,2) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * (-one)

        mat(iv,i,1,1,3) = mat(iv,i,1,1,3) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * two
        mat(iv,i,2,2,3) = mat(iv,i,2,2,3) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * two
        mat(iv,i,3,3,3) = mat(iv,i,3,3,3) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * two
        mat(iv,i,4,4,3) = mat(iv,i,4,4,3) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * two
        mat(iv,i,5,5,3) = mat(iv,i,5,5,3) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * two

        mat(iv,i,1,1,4) = mat(iv,i,1,1,4) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * (-one)
        mat(iv,i,2,2,4) = mat(iv,i,2,2,4) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * (-one)
        mat(iv,i,3,3,4) = mat(iv,i,3,3,4) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * (-one)
        mat(iv,i,4,4,4) = mat(iv,i,4,4,4) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * (-one)
        mat(iv,i,5,5,4) = mat(iv,i,5,5,4) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * (-one)

        end do
                
        end if

!.... sponge term

        if (.not. calcd) then

        if (ispg .eq. 1) then
        
          do iv = 1, ny
            mat(iv,i,1,1,3) = mat(iv,i,1,1,3) + alfa * dtl(iv,i) * spgl(iv,i)
            mat(iv,i,2,2,3) = mat(iv,i,2,2,3) + alfa * dtl(iv,i) * spgl(iv,i)
            mat(iv,i,3,3,3) = mat(iv,i,3,3,3) + alfa * dtl(iv,i) * spgl(iv,i)
            mat(iv,i,4,4,3) = mat(iv,i,4,4,3) + alfa * dtl(iv,i) * spgl(iv,i)
            mat(iv,i,5,5,3) = mat(iv,i,5,5,3) + alfa * dtl(iv,i) * spgl(iv,i)
          end do

        else if (ispg .ge. 2) then
        
          do iv = 1, ny
            mat(iv,i,1,1,3) = mat(iv,i,1,1,3) + alfa * dtl(iv,i) * (spgl(iv,i) + spg2l(iv,i))
            mat(iv,i,2,2,3) = mat(iv,i,2,2,3) + alfa * dtl(iv,i) * (spgl(iv,i) + spg2l(iv,i))
            mat(iv,i,3,3,3) = mat(iv,i,3,3,3) + alfa * dtl(iv,i) * (spgl(iv,i) + spg2l(iv,i))
            mat(iv,i,4,4,3) = mat(iv,i,4,4,3) + alfa * dtl(iv,i) * (spgl(iv,i) + spg2l(iv,i))
            mat(iv,i,5,5,3) = mat(iv,i,5,5,3) + alfa * dtl(iv,i) * (spgl(iv,i) + spg2l(iv,i))
          end do
        
        end if

        end if

        end if          ! lsym

        if (rsym) then

          call lhsbs1f3D( mat, Ah, Dh, Dhi, Vh, Ahi, spgl, spg2l, dtl, calcd, 2 )       

        else

!=======================================================================================================!
!.... use higher-order tangent on first node off the right boundary
!=======================================================================================================!
        i = nx-1

        a1 = -alfa * gb1 * dxiinv
        a2 = -alfa * gb2 * dxiinv
        a3 = -alfa * gb3 * dxiinv
        a4 = -alfa * gb4 * dxiinv
        a5 = -alfa * gb5 * dxiinv

        do idof = 1, ndof
          do jdof = 1, ndof
            do iv = 1, ny
              mat(iv,i,idof,jdof,5) = a5 * dtl(iv,i) * Ah(iv,i,idof,jdof)
              mat(iv,i,idof,jdof,1) = a4 * dtl(iv,i) * Ah(iv,i,idof,jdof)
              mat(iv,i,idof,jdof,2) = a3 * dtl(iv,i) * Ah(iv,i,idof,jdof)
              mat(iv,i,idof,jdof,3) = dtl(iv,i) * ( a2 * Ah(iv,i,idof,jdof) + &
                                      rcalcd * alfa * (  Dh(iv,i,idof,jdof) + &
                                                   im * Dhi(iv,i,idof,jdof) ) )
              mat(iv,i,idof,jdof,4) = a1 * dtl(iv,i) * Ah(iv,i,idof,jdof)
            end do
          end do
        end do
        
!.... \hat{A}_i term

        c1 = -im * alfa * gb1 * dxiinv
        c2 = -im * alfa * gb2 * dxiinv
        c3 = -im * alfa * gb3 * dxiinv
        c4 = -im * alfa * gb4 * dxiinv
        c5 = -im * alfa * gb5 * dxiinv

        do iv = 1, ny

        mat(iv,i,2,4,5) = mat(iv,i,2,4,5) - c5 * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,3,4,5) = mat(iv,i,3,4,5) - c5 * dtl(iv,i) * Ahi(iv,i,2)
        mat(iv,i,4,2,5) = mat(iv,i,4,2,5) - c5 * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,4,3,5) = mat(iv,i,4,3,5) - c5 * dtl(iv,i) * Ahi(iv,i,2)

        mat(iv,i,2,4,1) = mat(iv,i,2,4,1) - c4 * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,3,4,1) = mat(iv,i,3,4,1) - c4 * dtl(iv,i) * Ahi(iv,i,2)
        mat(iv,i,4,2,1) = mat(iv,i,4,2,1) - c4 * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,4,3,1) = mat(iv,i,4,3,1) - c4 * dtl(iv,i) * Ahi(iv,i,2)

        mat(iv,i,2,4,2) = mat(iv,i,2,4,2) - c3 * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,3,4,2) = mat(iv,i,3,4,2) - c3 * dtl(iv,i) * Ahi(iv,i,2)
        mat(iv,i,4,2,2) = mat(iv,i,4,2,2) - c3 * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,4,3,2) = mat(iv,i,4,3,2) - c3 * dtl(iv,i) * Ahi(iv,i,2)

        mat(iv,i,2,4,3) = mat(iv,i,2,4,3) - c2 * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,3,4,3) = mat(iv,i,3,4,3) - c2 * dtl(iv,i) * Ahi(iv,i,2)
        mat(iv,i,4,2,3) = mat(iv,i,4,2,3) - c2 * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,4,3,3) = mat(iv,i,4,3,3) - c2 * dtl(iv,i) * Ahi(iv,i,2)

        mat(iv,i,2,4,4) = mat(iv,i,2,4,4) - c1 * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,3,4,4) = mat(iv,i,3,4,4) - c1 * dtl(iv,i) * Ahi(iv,i,2)
        mat(iv,i,4,2,4) = mat(iv,i,4,2,4) - c1 * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,4,3,4) = mat(iv,i,4,3,4) - c1 * dtl(iv,i) * Ahi(iv,i,2)
          
        end do

!.... \hat{V}_{\xi\xi} term

        a1 = alfa * db1 * dxisinv
        a2 = alfa * db2 * dxisinv
        a3 = alfa * db3 * dxisinv
        a4 = alfa * db4 * dxisinv
        a5 = alfa * db5 * dxisinv

        do iv = 1, ny
        
        mat(iv,i,2,2,5) = mat(iv,i,2,2,5) - a5 * dtl(iv,i) * Vh(iv,i,1)
        mat(iv,i,2,3,5) = mat(iv,i,2,3,5) - a5 * dtl(iv,i) * Vh(iv,i,5)
        mat(iv,i,3,2,5) = mat(iv,i,3,2,5) - a5 * dtl(iv,i) * Vh(iv,i,6)
        mat(iv,i,3,3,5) = mat(iv,i,3,3,5) - a5 * dtl(iv,i) * Vh(iv,i,2)
        mat(iv,i,4,4,5) = mat(iv,i,4,4,5) - a5 * dtl(iv,i) * Vh(iv,i,3)
        mat(iv,i,5,5,5) = mat(iv,i,5,5,5) - a5 * dtl(iv,i) * Vh(iv,i,4)

        mat(iv,i,2,2,1) = mat(iv,i,2,2,1) - a4 * dtl(iv,i) * Vh(iv,i,1)
        mat(iv,i,2,3,1) = mat(iv,i,2,3,1) - a4 * dtl(iv,i) * Vh(iv,i,5)
        mat(iv,i,3,2,1) = mat(iv,i,3,2,1) - a4 * dtl(iv,i) * Vh(iv,i,6)
        mat(iv,i,3,3,1) = mat(iv,i,3,3,1) - a4 * dtl(iv,i) * Vh(iv,i,2)
        mat(iv,i,4,4,1) = mat(iv,i,4,4,1) - a4 * dtl(iv,i) * Vh(iv,i,3)
        mat(iv,i,5,5,1) = mat(iv,i,5,5,1) - a4 * dtl(iv,i) * Vh(iv,i,4)

        mat(iv,i,2,2,2) = mat(iv,i,2,2,2) - a3 * dtl(iv,i) * Vh(iv,i,1)
        mat(iv,i,2,3,2) = mat(iv,i,2,3,2) - a3 * dtl(iv,i) * Vh(iv,i,5)
        mat(iv,i,3,2,2) = mat(iv,i,3,2,2) - a3 * dtl(iv,i) * Vh(iv,i,6)
        mat(iv,i,3,3,2) = mat(iv,i,3,3,2) - a3 * dtl(iv,i) * Vh(iv,i,2)
        mat(iv,i,4,4,2) = mat(iv,i,4,4,2) - a3 * dtl(iv,i) * Vh(iv,i,3)
        mat(iv,i,5,5,2) = mat(iv,i,5,5,2) - a3 * dtl(iv,i) * Vh(iv,i,4)

        mat(iv,i,2,2,3) = mat(iv,i,2,2,3) - a2 * dtl(iv,i) * Vh(iv,i,1)
        mat(iv,i,2,3,3) = mat(iv,i,2,3,3) - a2 * dtl(iv,i) * Vh(iv,i,5)
        mat(iv,i,3,2,3) = mat(iv,i,3,2,3) - a2 * dtl(iv,i) * Vh(iv,i,6)
        mat(iv,i,3,3,3) = mat(iv,i,3,3,3) - a2 * dtl(iv,i) * Vh(iv,i,2)
        mat(iv,i,4,4,3) = mat(iv,i,4,4,3) - a2 * dtl(iv,i) * Vh(iv,i,3)
        mat(iv,i,5,5,3) = mat(iv,i,5,5,3) - a2 * dtl(iv,i) * Vh(iv,i,4)

        mat(iv,i,2,2,4) = mat(iv,i,2,2,4) - a1 * dtl(iv,i) * Vh(iv,i,1)
        mat(iv,i,2,3,4) = mat(iv,i,2,3,4) - a1 * dtl(iv,i) * Vh(iv,i,5)
        mat(iv,i,3,2,4) = mat(iv,i,3,2,4) - a1 * dtl(iv,i) * Vh(iv,i,6)
        mat(iv,i,3,3,4) = mat(iv,i,3,3,4) - a1 * dtl(iv,i) * Vh(iv,i,2)
        mat(iv,i,4,4,4) = mat(iv,i,4,4,4) - a1 * dtl(iv,i) * Vh(iv,i,3)
        mat(iv,i,5,5,4) = mat(iv,i,5,5,4) - a1 * dtl(iv,i) * Vh(iv,i,4)

!.... I term

        mat(iv,i,1,1,3) = mat(iv,i,1,1,3) + one
        mat(iv,i,2,2,3) = mat(iv,i,2,2,3) + one
        mat(iv,i,3,3,3) = mat(iv,i,3,3,3) + one
        mat(iv,i,4,4,3) = mat(iv,i,4,4,3) + one
        mat(iv,i,5,5,3) = mat(iv,i,5,5,3) + one

        end do
        
!.... implicit damping term

        if (eps_e .ne. zero) then

        do iv = 1, ny
        
        mat(iv,i,1,1,2) = mat(iv,i,1,1,2) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * (-one)
        mat(iv,i,2,2,2) = mat(iv,i,2,2,2) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * (-one)
        mat(iv,i,3,3,2) = mat(iv,i,3,3,2) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * (-one)
        mat(iv,i,4,4,2) = mat(iv,i,4,4,2) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * (-one)
        mat(iv,i,5,5,2) = mat(iv,i,5,5,2) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * (-one)

        mat(iv,i,1,1,3) = mat(iv,i,1,1,3) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * two
        mat(iv,i,2,2,3) = mat(iv,i,2,2,3) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * two
        mat(iv,i,3,3,3) = mat(iv,i,3,3,3) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * two
        mat(iv,i,4,4,3) = mat(iv,i,4,4,3) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * two
        mat(iv,i,5,5,3) = mat(iv,i,5,5,3) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * two

        mat(iv,i,1,1,4) = mat(iv,i,1,1,4) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * (-one)
        mat(iv,i,2,2,4) = mat(iv,i,2,2,4) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * (-one)
        mat(iv,i,3,3,4) = mat(iv,i,3,3,4) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * (-one)
        mat(iv,i,4,4,4) = mat(iv,i,4,4,4) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * (-one)
        mat(iv,i,5,5,4) = mat(iv,i,5,5,4) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * (-one)

        end do
        
        end if

!.... sponge term

        if (.not. calcd) then

        if (ispg .eq. 1) then
        
          do iv = 1, ny
            mat(iv,i,1,1,3) = mat(iv,i,1,1,3) + alfa * dtl(iv,i) * spgl(iv,i)
            mat(iv,i,2,2,3) = mat(iv,i,2,2,3) + alfa * dtl(iv,i) * spgl(iv,i)
            mat(iv,i,3,3,3) = mat(iv,i,3,3,3) + alfa * dtl(iv,i) * spgl(iv,i)
            mat(iv,i,4,4,3) = mat(iv,i,4,4,3) + alfa * dtl(iv,i) * spgl(iv,i)
            mat(iv,i,5,5,3) = mat(iv,i,5,5,3) + alfa * dtl(iv,i) * spgl(iv,i)
          end do

        else if (ispg .ge. 2) then
        
          do iv = 1, ny
            mat(iv,i,1,1,3) = mat(iv,i,1,1,3) + alfa * dtl(iv,i) * (spgl(iv,i) + spg2l(iv,i))
            mat(iv,i,2,2,3) = mat(iv,i,2,2,3) + alfa * dtl(iv,i) * (spgl(iv,i) + spg2l(iv,i))
            mat(iv,i,3,3,3) = mat(iv,i,3,3,3) + alfa * dtl(iv,i) * (spgl(iv,i) + spg2l(iv,i))
            mat(iv,i,4,4,3) = mat(iv,i,4,4,3) + alfa * dtl(iv,i) * (spgl(iv,i) + spg2l(iv,i))
            mat(iv,i,5,5,3) = mat(iv,i,5,5,3) + alfa * dtl(iv,i) * (spgl(iv,i) + spg2l(iv,i))
          end do
        
        end if

        end if
!=======================================================================================================!
!.... use higher-order tangent on the right boundary
!=======================================================================================================!
        i = nx

        a1 = -alfa * gc1 * dxiinv
        a2 = -alfa * gc2 * dxiinv
        a3 = -alfa * gc3 * dxiinv
        a4 = -alfa * gc4 * dxiinv
        a5 = -alfa * gc5 * dxiinv

        do idof = 1, ndof
          do jdof = 1, ndof
            do iv = 1, ny
              mat(iv,i,idof,jdof,4) = a5 * dtl(iv,i) * Ah(iv,i,idof,jdof)
              mat(iv,i,idof,jdof,5) = a4 * dtl(iv,i) * Ah(iv,i,idof,jdof)
              mat(iv,i,idof,jdof,1) = a3 * dtl(iv,i) * Ah(iv,i,idof,jdof)
              mat(iv,i,idof,jdof,2) = a2 * dtl(iv,i) * Ah(iv,i,idof,jdof)
              mat(iv,i,idof,jdof,3) = dtl(iv,i) * ( a1 * Ah(iv,i,idof,jdof) + &
                                      rcalcd * alfa * (  Dh(iv,i,idof,jdof) + &
                                                   im * Dhi(iv,i,idof,jdof) ) )
            end do
          end do
        end do
        
!.... \hat{A}_i term

        c1 = -im * alfa * gc1 * dxiinv
        c2 = -im * alfa * gc2 * dxiinv
        c3 = -im * alfa * gc3 * dxiinv
        c4 = -im * alfa * gc4 * dxiinv
        c5 = -im * alfa * gc5 * dxiinv

        do iv = 1, ny

        mat(iv,i,2,4,4) = mat(iv,i,2,4,4) - c5 * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,3,4,4) = mat(iv,i,3,4,4) - c5 * dtl(iv,i) * Ahi(iv,i,2)
        mat(iv,i,4,2,4) = mat(iv,i,4,2,4) - c5 * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,4,3,4) = mat(iv,i,4,3,4) - c5 * dtl(iv,i) * Ahi(iv,i,2)

        mat(iv,i,2,4,5) = mat(iv,i,2,4,5) - c4 * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,3,4,5) = mat(iv,i,3,4,5) - c4 * dtl(iv,i) * Ahi(iv,i,2)
        mat(iv,i,4,2,5) = mat(iv,i,4,2,5) - c4 * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,4,3,5) = mat(iv,i,4,3,5) - c4 * dtl(iv,i) * Ahi(iv,i,2)

        mat(iv,i,2,4,1) = mat(iv,i,2,4,1) - c3 * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,3,4,1) = mat(iv,i,3,4,1) - c3 * dtl(iv,i) * Ahi(iv,i,2)
        mat(iv,i,4,2,1) = mat(iv,i,4,2,1) - c3 * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,4,3,1) = mat(iv,i,4,3,1) - c3 * dtl(iv,i) * Ahi(iv,i,2)

        mat(iv,i,2,4,2) = mat(iv,i,2,4,2) - c2 * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,3,4,2) = mat(iv,i,3,4,2) - c2 * dtl(iv,i) * Ahi(iv,i,2)
        mat(iv,i,4,2,2) = mat(iv,i,4,2,2) - c2 * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,4,3,2) = mat(iv,i,4,3,2) - c2 * dtl(iv,i) * Ahi(iv,i,2)

        mat(iv,i,2,4,3) = mat(iv,i,2,4,3) - c1 * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,3,4,3) = mat(iv,i,3,4,3) - c1 * dtl(iv,i) * Ahi(iv,i,2)
        mat(iv,i,4,2,3) = mat(iv,i,4,2,3) - c1 * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,4,3,3) = mat(iv,i,4,3,3) - c1 * dtl(iv,i) * Ahi(iv,i,2)
          
        end do

!.... \hat{V}_{\xi\xi} term

        a1 = alfa * dd1 * dxisinv
        a2 = alfa * dd2 * dxisinv
        a3 = alfa * dd3 * dxisinv
        a4 = alfa * dd4 * dxisinv
        a5 = alfa * dd5 * dxisinv

        do iv = 1, ny
        
        mat(iv,i,2,2,4) = mat(iv,i,2,2,4) - a5 * dtl(iv,i) * Vh(iv,i,1)
        mat(iv,i,2,3,4) = mat(iv,i,2,3,4) - a5 * dtl(iv,i) * Vh(iv,i,5)
        mat(iv,i,3,2,4) = mat(iv,i,3,2,4) - a5 * dtl(iv,i) * Vh(iv,i,6)
        mat(iv,i,3,3,4) = mat(iv,i,3,3,4) - a5 * dtl(iv,i) * Vh(iv,i,2)
        mat(iv,i,4,4,4) = mat(iv,i,4,4,4) - a5 * dtl(iv,i) * Vh(iv,i,3)
        mat(iv,i,5,5,4) = mat(iv,i,5,5,4) - a5 * dtl(iv,i) * Vh(iv,i,4)

        mat(iv,i,2,2,5) = mat(iv,i,2,2,5) - a4 * dtl(iv,i) * Vh(iv,i,1)
        mat(iv,i,2,3,5) = mat(iv,i,2,3,5) - a4 * dtl(iv,i) * Vh(iv,i,5)
        mat(iv,i,3,2,5) = mat(iv,i,3,2,5) - a4 * dtl(iv,i) * Vh(iv,i,6)
        mat(iv,i,3,3,5) = mat(iv,i,3,3,5) - a4 * dtl(iv,i) * Vh(iv,i,2)
        mat(iv,i,4,4,5) = mat(iv,i,4,4,5) - a4 * dtl(iv,i) * Vh(iv,i,3)
        mat(iv,i,5,5,5) = mat(iv,i,5,5,5) - a4 * dtl(iv,i) * Vh(iv,i,4)

        mat(iv,i,2,2,1) = mat(iv,i,2,2,1) - a3 * dtl(iv,i) * Vh(iv,i,1)
        mat(iv,i,2,3,1) = mat(iv,i,2,3,1) - a3 * dtl(iv,i) * Vh(iv,i,5)
        mat(iv,i,3,2,1) = mat(iv,i,3,2,1) - a3 * dtl(iv,i) * Vh(iv,i,6)
        mat(iv,i,3,3,1) = mat(iv,i,3,3,1) - a3 * dtl(iv,i) * Vh(iv,i,2)
        mat(iv,i,4,4,1) = mat(iv,i,4,4,1) - a3 * dtl(iv,i) * Vh(iv,i,3)
        mat(iv,i,5,5,1) = mat(iv,i,5,5,1) - a3 * dtl(iv,i) * Vh(iv,i,4)

        mat(iv,i,2,2,2) = mat(iv,i,2,2,2) - a2 * dtl(iv,i) * Vh(iv,i,1)
        mat(iv,i,2,3,2) = mat(iv,i,2,3,2) - a2 * dtl(iv,i) * Vh(iv,i,5)
        mat(iv,i,3,2,2) = mat(iv,i,3,2,2) - a2 * dtl(iv,i) * Vh(iv,i,6)
        mat(iv,i,3,3,2) = mat(iv,i,3,3,2) - a2 * dtl(iv,i) * Vh(iv,i,2)
        mat(iv,i,4,4,2) = mat(iv,i,4,4,2) - a2 * dtl(iv,i) * Vh(iv,i,3)
        mat(iv,i,5,5,2) = mat(iv,i,5,5,2) - a2 * dtl(iv,i) * Vh(iv,i,4)

        mat(iv,i,2,2,3) = mat(iv,i,2,2,3) - a1 * dtl(iv,i) * Vh(iv,i,1)
        mat(iv,i,2,3,3) = mat(iv,i,2,3,3) - a1 * dtl(iv,i) * Vh(iv,i,5)
        mat(iv,i,3,2,3) = mat(iv,i,3,2,3) - a1 * dtl(iv,i) * Vh(iv,i,6)
        mat(iv,i,3,3,3) = mat(iv,i,3,3,3) - a1 * dtl(iv,i) * Vh(iv,i,2)
        mat(iv,i,4,4,3) = mat(iv,i,4,4,3) - a1 * dtl(iv,i) * Vh(iv,i,3)
        mat(iv,i,5,5,3) = mat(iv,i,5,5,3) - a1 * dtl(iv,i) * Vh(iv,i,4)

!.... I term

        mat(iv,i,1,1,3) = mat(iv,i,1,1,3) + one
        mat(iv,i,2,2,3) = mat(iv,i,2,2,3) + one
        mat(iv,i,3,3,3) = mat(iv,i,3,3,3) + one
        mat(iv,i,4,4,3) = mat(iv,i,4,4,3) + one
        mat(iv,i,5,5,3) = mat(iv,i,5,5,3) + one

        end do

!.... sponge term

        if (.not. calcd) then

        if (ispg .eq. 1) then
        
          do iv = 1, ny
            mat(iv,i,1,1,3) = mat(iv,i,1,1,3) + alfa * dtl(iv,i) * spgl(iv,i)
            mat(iv,i,2,2,3) = mat(iv,i,2,2,3) + alfa * dtl(iv,i) * spgl(iv,i)
            mat(iv,i,3,3,3) = mat(iv,i,3,3,3) + alfa * dtl(iv,i) * spgl(iv,i)
            mat(iv,i,4,4,3) = mat(iv,i,4,4,3) + alfa * dtl(iv,i) * spgl(iv,i)
            mat(iv,i,5,5,3) = mat(iv,i,5,5,3) + alfa * dtl(iv,i) * spgl(iv,i)
          end do

        else if (ispg .ge. 2) then
        
          do iv = 1, ny
            mat(iv,i,1,1,3) = mat(iv,i,1,1,3) + alfa * dtl(iv,i) * (spgl(iv,i) + spg2l(iv,i))
            mat(iv,i,2,2,3) = mat(iv,i,2,2,3) + alfa * dtl(iv,i) * (spgl(iv,i) + spg2l(iv,i))
            mat(iv,i,3,3,3) = mat(iv,i,3,3,3) + alfa * dtl(iv,i) * (spgl(iv,i) + spg2l(iv,i))
            mat(iv,i,4,4,3) = mat(iv,i,4,4,3) + alfa * dtl(iv,i) * (spgl(iv,i) + spg2l(iv,i))
            mat(iv,i,5,5,3) = mat(iv,i,5,5,3) + alfa * dtl(iv,i) * (spgl(iv,i) + spg2l(iv,i))
          end do
        
        end if

        end if

        end if   ! rsym

        return
        end
