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
        use buffer
        implicit none
        
        complex :: mat(5,ndof,ndof,nx,ny)
        real    :: Ah(ndof,ndof,nx,ny), Dh(ndof,ndof,nx,ny)
        real    :: Dhi(ndof,ndof,nx,ny), Ahi(ny,nx,6)
        real    :: spgl(ny,nx), spg2l(ny,nx), Vh(6,nx,ny), dtl(nx,ny)
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
                
        !$omp parallel do private (i,idof,jdof,a1,a2,a3,a4,a5,c1,c2,c3,c4,c5)
        do iv = 1, ny
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
              mat(3,idof,jdof,i,iv) = dtl(i,iv) * ( a1 * Ah(idof,jdof,i,iv) + &
                                      rcalcd * alfa * (  Dh(idof,jdof,i,iv) + &
                                                   im * Dhi(idof,jdof,i,iv) ) )
              mat(4,idof,jdof,i,iv) = a2 * dtl(i,iv) * Ah(idof,jdof,i,iv)
              mat(5,idof,jdof,i,iv) = a3 * dtl(i,iv) * Ah(idof,jdof,i,iv)
              mat(1,idof,jdof,i,iv) = a4 * dtl(i,iv) * Ah(idof,jdof,i,iv)
              mat(2,idof,jdof,i,iv) = a5 * dtl(i,iv) * Ah(idof,jdof,i,iv)
          end do
        end do
        
!.... \hat{A}_i term

        c1 = im * alfa * gc1 * dxiinv
        c2 = im * alfa * gc2 * dxiinv
        c3 = im * alfa * gc3 * dxiinv
        c4 = im * alfa * gc4 * dxiinv
        c5 = im * alfa * gc5 * dxiinv

        mat(3,2,4,i,iv) = mat(3,2,4,i,iv) - c1 * dtl(i,iv) * Ahi(iv,i,1)
        mat(3,3,4,i,iv) = mat(3,3,4,i,iv) - c1 * dtl(i,iv) * Ahi(iv,i,2)
        mat(3,4,2,i,iv) = mat(3,4,2,i,iv) - c1 * dtl(i,iv) * Ahi(iv,i,1)
        mat(3,4,3,i,iv) = mat(3,4,3,i,iv) - c1 * dtl(i,iv) * Ahi(iv,i,2)

        mat(4,2,4,i,iv) = mat(4,2,4,i,iv) - c2 * dtl(i,iv) * Ahi(iv,i,1)
        mat(4,3,4,i,iv) = mat(4,3,4,i,iv) - c2 * dtl(i,iv) * Ahi(iv,i,2)
        mat(4,4,2,i,iv) = mat(4,4,2,i,iv) - c2 * dtl(i,iv) * Ahi(iv,i,1)
        mat(4,4,3,i,iv) = mat(4,4,3,i,iv) - c2 * dtl(i,iv) * Ahi(iv,i,2)

        mat(5,2,4,i,iv) = mat(5,2,4,i,iv) - c3 * dtl(i,iv) * Ahi(iv,i,1)
        mat(5,3,4,i,iv) = mat(5,3,4,i,iv) - c3 * dtl(i,iv) * Ahi(iv,i,2)
        mat(5,4,2,i,iv) = mat(5,4,2,i,iv) - c3 * dtl(i,iv) * Ahi(iv,i,1)
        mat(5,4,3,i,iv) = mat(5,4,3,i,iv) - c3 * dtl(i,iv) * Ahi(iv,i,2)

        mat(1,2,4,i,iv) = mat(1,2,4,i,iv) - c4 * dtl(i,iv) * Ahi(iv,i,1)
        mat(1,3,4,i,iv) = mat(1,3,4,i,iv) - c4 * dtl(i,iv) * Ahi(iv,i,2)
        mat(1,4,2,i,iv) = mat(1,4,2,i,iv) - c4 * dtl(i,iv) * Ahi(iv,i,1)
        mat(1,4,3,i,iv) = mat(1,4,3,i,iv) - c4 * dtl(i,iv) * Ahi(iv,i,2)

        mat(2,2,4,i,iv) = mat(2,2,4,i,iv) - c5 * dtl(i,iv) * Ahi(iv,i,1)
        mat(2,3,4,i,iv) = mat(2,3,4,i,iv) - c5 * dtl(i,iv) * Ahi(iv,i,2)
        mat(2,4,2,i,iv) = mat(2,4,2,i,iv) - c5 * dtl(i,iv) * Ahi(iv,i,1)
        mat(2,4,3,i,iv) = mat(2,4,3,i,iv) - c5 * dtl(i,iv) * Ahi(iv,i,2)
          
!.... \hat{V}_{\xi\xi} term

        a1 = alfa * dd1 * dxisinv
        a2 = alfa * dd2 * dxisinv
        a3 = alfa * dd3 * dxisinv
        a4 = alfa * dd4 * dxisinv
        a5 = alfa * dd5 * dxisinv

        mat(3,2,2,i,iv) = mat(3,2,2,i,iv) - a1 * dtl(i,iv) * Vh(1,i,iv)
        mat(3,2,3,i,iv) = mat(3,2,3,i,iv) - a1 * dtl(i,iv) * Vh(5,i,iv)
        mat(3,3,2,i,iv) = mat(3,3,2,i,iv) - a1 * dtl(i,iv) * Vh(6,i,iv)
        mat(3,3,3,i,iv) = mat(3,3,3,i,iv) - a1 * dtl(i,iv) * Vh(2,i,iv)
        mat(3,4,4,i,iv) = mat(3,4,4,i,iv) - a1 * dtl(i,iv) * Vh(3,i,iv)
        mat(3,5,5,i,iv) = mat(3,5,5,i,iv) - a1 * dtl(i,iv) * Vh(4,i,iv)

        mat(4,2,2,i,iv) = mat(4,2,2,i,iv) - a2 * dtl(i,iv) * Vh(1,i,iv)
        mat(4,2,3,i,iv) = mat(4,2,3,i,iv) - a2 * dtl(i,iv) * Vh(5,i,iv)
        mat(4,3,2,i,iv) = mat(4,3,2,i,iv) - a2 * dtl(i,iv) * Vh(6,i,iv)
        mat(4,3,3,i,iv) = mat(4,3,3,i,iv) - a2 * dtl(i,iv) * Vh(2,i,iv)
        mat(4,4,4,i,iv) = mat(4,4,4,i,iv) - a2 * dtl(i,iv) * Vh(3,i,iv)
        mat(4,5,5,i,iv) = mat(4,5,5,i,iv) - a2 * dtl(i,iv) * Vh(4,i,iv)

        mat(5,2,2,i,iv) = mat(5,2,2,i,iv) - a3 * dtl(i,iv) * Vh(1,i,iv)
        mat(5,2,3,i,iv) = mat(5,2,3,i,iv) - a3 * dtl(i,iv) * Vh(5,i,iv)
        mat(5,3,2,i,iv) = mat(5,3,2,i,iv) - a3 * dtl(i,iv) * Vh(6,i,iv)
        mat(5,3,3,i,iv) = mat(5,3,3,i,iv) - a3 * dtl(i,iv) * Vh(2,i,iv)
        mat(5,4,4,i,iv) = mat(5,4,4,i,iv) - a3 * dtl(i,iv) * Vh(3,i,iv)
        mat(5,5,5,i,iv) = mat(5,5,5,i,iv) - a3 * dtl(i,iv) * Vh(4,i,iv)

        mat(1,2,2,i,iv) = mat(1,2,2,i,iv) - a4 * dtl(i,iv) * Vh(1,i,iv)
        mat(1,2,3,i,iv) = mat(1,2,3,i,iv) - a4 * dtl(i,iv) * Vh(5,i,iv)
        mat(1,3,2,i,iv) = mat(1,3,2,i,iv) - a4 * dtl(i,iv) * Vh(6,i,iv)
        mat(1,3,3,i,iv) = mat(1,3,3,i,iv) - a4 * dtl(i,iv) * Vh(2,i,iv)
        mat(1,4,4,i,iv) = mat(1,4,4,i,iv) - a4 * dtl(i,iv) * Vh(3,i,iv)
        mat(1,5,5,i,iv) = mat(1,5,5,i,iv) - a4 * dtl(i,iv) * Vh(4,i,iv)

        mat(2,2,2,i,iv) = mat(2,2,2,i,iv) - a5 * dtl(i,iv) * Vh(1,i,iv)
        mat(2,2,3,i,iv) = mat(2,2,3,i,iv) - a5 * dtl(i,iv) * Vh(5,i,iv)
        mat(2,3,2,i,iv) = mat(2,3,2,i,iv) - a5 * dtl(i,iv) * Vh(6,i,iv)
        mat(2,3,3,i,iv) = mat(2,3,3,i,iv) - a5 * dtl(i,iv) * Vh(2,i,iv)
        mat(2,4,4,i,iv) = mat(2,4,4,i,iv) - a5 * dtl(i,iv) * Vh(3,i,iv)
        mat(2,5,5,i,iv) = mat(2,5,5,i,iv) - a5 * dtl(i,iv) * Vh(4,i,iv)

!.... I term

        mat(3,1,1,i,iv) = mat(3,1,1,i,iv) + one
        mat(3,2,2,i,iv) = mat(3,2,2,i,iv) + one
        mat(3,3,3,i,iv) = mat(3,3,3,i,iv) + one
        mat(3,4,4,i,iv) = mat(3,4,4,i,iv) + one
        mat(3,5,5,i,iv) = mat(3,5,5,i,iv) + one

!.... sponge term

        if (.not. calcd) then

        if (ispg .eq. 1) then
        
            mat(3,1,1,i,iv) = mat(3,1,1,i,iv) + alfa * dtl(i,iv) * spgl(iv,i)
            mat(3,2,2,i,iv) = mat(3,2,2,i,iv) + alfa * dtl(i,iv) * spgl(iv,i)
            mat(3,3,3,i,iv) = mat(3,3,3,i,iv) + alfa * dtl(i,iv) * spgl(iv,i)
            mat(3,4,4,i,iv) = mat(3,4,4,i,iv) + alfa * dtl(i,iv) * spgl(iv,i)
            mat(3,5,5,i,iv) = mat(3,5,5,i,iv) + alfa * dtl(i,iv) * spgl(iv,i)

        else if (ispg .ge. 2) then
        
            mat(3,1,1,i,iv) = mat(3,1,1,i,iv) + alfa * dtl(i,iv) * (spgl(iv,i) + spg2l(iv,i))
            mat(3,2,2,i,iv) = mat(3,2,2,i,iv) + alfa * dtl(i,iv) * (spgl(iv,i) + spg2l(iv,i))
            mat(3,3,3,i,iv) = mat(3,3,3,i,iv) + alfa * dtl(i,iv) * (spgl(iv,i) + spg2l(iv,i))
            mat(3,4,4,i,iv) = mat(3,4,4,i,iv) + alfa * dtl(i,iv) * (spgl(iv,i) + spg2l(iv,i))
            mat(3,5,5,i,iv) = mat(3,5,5,i,iv) + alfa * dtl(i,iv) * (spgl(iv,i) + spg2l(iv,i))
        
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
              mat(2,idof,jdof,i,iv) = a1 * dtl(i,iv) * Ah(idof,jdof,i,iv)
              mat(3,idof,jdof,i,iv) = dtl(i,iv) * ( a2 * Ah(idof,jdof,i,iv) + &
                                      rcalcd * alfa * (  Dh(idof,jdof,i,iv) + &
                                                   im * Dhi(idof,jdof,i,iv) ) )
              mat(4,idof,jdof,i,iv) = a3 * dtl(i,iv) * Ah(idof,jdof,i,iv)
              mat(5,idof,jdof,i,iv) = a4 * dtl(i,iv) * Ah(idof,jdof,i,iv)
              mat(1,idof,jdof,i,iv) = a5 * dtl(i,iv) * Ah(idof,jdof,i,iv)
          end do
        end do
        
!.... \hat{A}_i term

        c1 = im * alfa * gb1 * dxiinv
        c2 = im * alfa * gb2 * dxiinv
        c3 = im * alfa * gb3 * dxiinv
        c4 = im * alfa * gb4 * dxiinv
        c5 = im * alfa * gb5 * dxiinv

        mat(2,2,4,i,iv) = mat(2,2,4,i,iv) - c1 * dtl(i,iv) * Ahi(iv,i,1)
        mat(2,3,4,i,iv) = mat(2,3,4,i,iv) - c1 * dtl(i,iv) * Ahi(iv,i,2)
        mat(2,4,2,i,iv) = mat(2,4,2,i,iv) - c1 * dtl(i,iv) * Ahi(iv,i,1)
        mat(2,4,3,i,iv) = mat(2,4,3,i,iv) - c1 * dtl(i,iv) * Ahi(iv,i,2)

        mat(3,2,4,i,iv) = mat(3,2,4,i,iv) - c2 * dtl(i,iv) * Ahi(iv,i,1)
        mat(3,3,4,i,iv) = mat(3,3,4,i,iv) - c2 * dtl(i,iv) * Ahi(iv,i,2)
        mat(3,4,2,i,iv) = mat(3,4,2,i,iv) - c2 * dtl(i,iv) * Ahi(iv,i,1)
        mat(3,4,3,i,iv) = mat(3,4,3,i,iv) - c2 * dtl(i,iv) * Ahi(iv,i,2)

        mat(4,2,4,i,iv) = mat(4,2,4,i,iv) - c3 * dtl(i,iv) * Ahi(iv,i,1)
        mat(4,3,4,i,iv) = mat(4,3,4,i,iv) - c3 * dtl(i,iv) * Ahi(iv,i,2)
        mat(4,4,2,i,iv) = mat(4,4,2,i,iv) - c3 * dtl(i,iv) * Ahi(iv,i,1)
        mat(4,4,3,i,iv) = mat(4,4,3,i,iv) - c3 * dtl(i,iv) * Ahi(iv,i,2)

        mat(5,2,4,i,iv) = mat(5,2,4,i,iv) - c4 * dtl(i,iv) * Ahi(iv,i,1)
        mat(5,3,4,i,iv) = mat(5,3,4,i,iv) - c4 * dtl(i,iv) * Ahi(iv,i,2)
        mat(5,4,2,i,iv) = mat(5,4,2,i,iv) - c4 * dtl(i,iv) * Ahi(iv,i,1)
        mat(5,4,3,i,iv) = mat(5,4,3,i,iv) - c4 * dtl(i,iv) * Ahi(iv,i,2)

        mat(1,2,4,i,iv) = mat(1,2,4,i,iv) - c5 * dtl(i,iv) * Ahi(iv,i,1)
        mat(1,3,4,i,iv) = mat(1,3,4,i,iv) - c5 * dtl(i,iv) * Ahi(iv,i,2)
        mat(1,4,2,i,iv) = mat(1,4,2,i,iv) - c5 * dtl(i,iv) * Ahi(iv,i,1)
        mat(1,4,3,i,iv) = mat(1,4,3,i,iv) - c5 * dtl(i,iv) * Ahi(iv,i,2)
          
!.... \hat{V}_{\xi\xi} term

        a1 = alfa * db1 * dxisinv
        a2 = alfa * db2 * dxisinv
        a3 = alfa * db3 * dxisinv
        a4 = alfa * db4 * dxisinv
        a5 = alfa * db5 * dxisinv

        mat(2,2,2,i,iv) = mat(2,2,2,i,iv) - a1 * dtl(i,iv) * Vh(1,i,iv)
        mat(2,2,3,i,iv) = mat(2,2,3,i,iv) - a1 * dtl(i,iv) * Vh(5,i,iv)
        mat(2,3,2,i,iv) = mat(2,3,2,i,iv) - a1 * dtl(i,iv) * Vh(6,i,iv)
        mat(2,3,3,i,iv) = mat(2,3,3,i,iv) - a1 * dtl(i,iv) * Vh(2,i,iv)
        mat(2,4,4,i,iv) = mat(2,4,4,i,iv) - a1 * dtl(i,iv) * Vh(3,i,iv)
        mat(2,5,5,i,iv) = mat(2,5,5,i,iv) - a1 * dtl(i,iv) * Vh(4,i,iv)

        mat(3,2,2,i,iv) = mat(3,2,2,i,iv) - a2 * dtl(i,iv) * Vh(1,i,iv)
        mat(3,2,3,i,iv) = mat(3,2,3,i,iv) - a2 * dtl(i,iv) * Vh(5,i,iv)
        mat(3,3,2,i,iv) = mat(3,3,2,i,iv) - a2 * dtl(i,iv) * Vh(6,i,iv)
        mat(3,3,3,i,iv) = mat(3,3,3,i,iv) - a2 * dtl(i,iv) * Vh(2,i,iv)
        mat(3,4,4,i,iv) = mat(3,4,4,i,iv) - a2 * dtl(i,iv) * Vh(3,i,iv)
        mat(3,5,5,i,iv) = mat(3,5,5,i,iv) - a2 * dtl(i,iv) * Vh(4,i,iv)

        mat(4,2,2,i,iv) = mat(4,2,2,i,iv) - a3 * dtl(i,iv) * Vh(1,i,iv)
        mat(4,2,3,i,iv) = mat(4,2,3,i,iv) - a3 * dtl(i,iv) * Vh(5,i,iv)
        mat(4,3,2,i,iv) = mat(4,3,2,i,iv) - a3 * dtl(i,iv) * Vh(6,i,iv)
        mat(4,3,3,i,iv) = mat(4,3,3,i,iv) - a3 * dtl(i,iv) * Vh(2,i,iv)
        mat(4,4,4,i,iv) = mat(4,4,4,i,iv) - a3 * dtl(i,iv) * Vh(3,i,iv)
        mat(4,5,5,i,iv) = mat(4,5,5,i,iv) - a3 * dtl(i,iv) * Vh(4,i,iv)

        mat(5,2,2,i,iv) = mat(5,2,2,i,iv) - a4 * dtl(i,iv) * Vh(1,i,iv)
        mat(5,2,3,i,iv) = mat(5,2,3,i,iv) - a4 * dtl(i,iv) * Vh(5,i,iv)
        mat(5,3,2,i,iv) = mat(5,3,2,i,iv) - a4 * dtl(i,iv) * Vh(6,i,iv)
        mat(5,3,3,i,iv) = mat(5,3,3,i,iv) - a4 * dtl(i,iv) * Vh(2,i,iv)
        mat(5,4,4,i,iv) = mat(5,4,4,i,iv) - a4 * dtl(i,iv) * Vh(3,i,iv)
        mat(5,5,5,i,iv) = mat(5,5,5,i,iv) - a4 * dtl(i,iv) * Vh(4,i,iv)

        mat(1,2,2,i,iv) = mat(1,2,2,i,iv) - a5 * dtl(i,iv) * Vh(1,i,iv)
        mat(1,2,3,i,iv) = mat(1,2,3,i,iv) - a5 * dtl(i,iv) * Vh(5,i,iv)
        mat(1,3,2,i,iv) = mat(1,3,2,i,iv) - a5 * dtl(i,iv) * Vh(6,i,iv)
        mat(1,3,3,i,iv) = mat(1,3,3,i,iv) - a5 * dtl(i,iv) * Vh(2,i,iv)
        mat(1,4,4,i,iv) = mat(1,4,4,i,iv) - a5 * dtl(i,iv) * Vh(3,i,iv)
        mat(1,5,5,i,iv) = mat(1,5,5,i,iv) - a5 * dtl(i,iv) * Vh(4,i,iv)

!.... I term

        mat(3,1,1,i,iv) = mat(3,1,1,i,iv) + one
        mat(3,2,2,i,iv) = mat(3,2,2,i,iv) + one
        mat(3,3,3,i,iv) = mat(3,3,3,i,iv) + one
        mat(3,4,4,i,iv) = mat(3,4,4,i,iv) + one
        mat(3,5,5,i,iv) = mat(3,5,5,i,iv) + one

!.... implicit damping term

        if (eps_e .ne. zero) then

        mat(2,1,1,i,iv) = mat(2,1,1,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (-one)
        mat(2,2,2,i,iv) = mat(2,2,2,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (-one)
        mat(2,3,3,i,iv) = mat(2,3,3,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (-one)
        mat(2,4,4,i,iv) = mat(2,4,4,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (-one)
        mat(2,5,5,i,iv) = mat(2,5,5,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (-one)

        mat(3,1,1,i,iv) = mat(3,1,1,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * two
        mat(3,2,2,i,iv) = mat(3,2,2,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * two
        mat(3,3,3,i,iv) = mat(3,3,3,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * two
        mat(3,4,4,i,iv) = mat(3,4,4,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * two
        mat(3,5,5,i,iv) = mat(3,5,5,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * two

        mat(4,1,1,i,iv) = mat(4,1,1,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (-one)
        mat(4,2,2,i,iv) = mat(4,2,2,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (-one)
        mat(4,3,3,i,iv) = mat(4,3,3,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (-one)
        mat(4,4,4,i,iv) = mat(4,4,4,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (-one)
        mat(4,5,5,i,iv) = mat(4,5,5,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (-one)

        end if

!.... sponge term

        if (.not. calcd) then

        if (ispg .eq. 1) then
        
            mat(3,1,1,i,iv) = mat(3,1,1,i,iv) + alfa * dtl(i,iv) * spgl(iv,i)
            mat(3,2,2,i,iv) = mat(3,2,2,i,iv) + alfa * dtl(i,iv) * spgl(iv,i)
            mat(3,3,3,i,iv) = mat(3,3,3,i,iv) + alfa * dtl(i,iv) * spgl(iv,i)
            mat(3,4,4,i,iv) = mat(3,4,4,i,iv) + alfa * dtl(i,iv) * spgl(iv,i)
            mat(3,5,5,i,iv) = mat(3,5,5,i,iv) + alfa * dtl(i,iv) * spgl(iv,i)

        else if (ispg .ge. 2) then
        
            mat(3,1,1,i,iv) = mat(3,1,1,i,iv) + alfa * dtl(i,iv) * (spgl(iv,i) + spg2l(iv,i))
            mat(3,2,2,i,iv) = mat(3,2,2,i,iv) + alfa * dtl(i,iv) * (spgl(iv,i) + spg2l(iv,i))
            mat(3,3,3,i,iv) = mat(3,3,3,i,iv) + alfa * dtl(i,iv) * (spgl(iv,i) + spg2l(iv,i))
            mat(3,4,4,i,iv) = mat(3,4,4,i,iv) + alfa * dtl(i,iv) * (spgl(iv,i) + spg2l(iv,i))
            mat(3,5,5,i,iv) = mat(3,5,5,i,iv) + alfa * dtl(i,iv) * (spgl(iv,i) + spg2l(iv,i))
        
        end if

        end if

        end do

        end if          ! lsym

        if (rsym) then

          call lhsbs1f3D( mat, Ah, Dh, Dhi, Vh, Ahi, spgl, spg2l, dtl, calcd, 2 )       

        else

        !$omp parallel do private (i,idof,jdof,a1,a2,a3,a4,a5,c1,c2,c3,c4,c5)
        do iv = 1, ny
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
              mat(5,idof,jdof,i,iv) = a5 * dtl(i,iv) * Ah(idof,jdof,i,iv)
              mat(1,idof,jdof,i,iv) = a4 * dtl(i,iv) * Ah(idof,jdof,i,iv)
              mat(2,idof,jdof,i,iv) = a3 * dtl(i,iv) * Ah(idof,jdof,i,iv)
              mat(3,idof,jdof,i,iv) = dtl(i,iv) * ( a2 * Ah(idof,jdof,i,iv) + &
                                      rcalcd * alfa * (  Dh(idof,jdof,i,iv) + &
                                                   im * Dhi(idof,jdof,i,iv) ) )
              mat(4,idof,jdof,i,iv) = a1 * dtl(i,iv) * Ah(idof,jdof,i,iv)
          end do
        end do
        
!.... \hat{A}_i term

        c1 = -im * alfa * gb1 * dxiinv
        c2 = -im * alfa * gb2 * dxiinv
        c3 = -im * alfa * gb3 * dxiinv
        c4 = -im * alfa * gb4 * dxiinv
        c5 = -im * alfa * gb5 * dxiinv

        mat(5,2,4,i,iv) = mat(5,2,4,i,iv) - c5 * dtl(i,iv) * Ahi(iv,i,1)
        mat(5,3,4,i,iv) = mat(5,3,4,i,iv) - c5 * dtl(i,iv) * Ahi(iv,i,2)
        mat(5,4,2,i,iv) = mat(5,4,2,i,iv) - c5 * dtl(i,iv) * Ahi(iv,i,1)
        mat(5,4,3,i,iv) = mat(5,4,3,i,iv) - c5 * dtl(i,iv) * Ahi(iv,i,2)

        mat(1,2,4,i,iv) = mat(1,2,4,i,iv) - c4 * dtl(i,iv) * Ahi(iv,i,1)
        mat(1,3,4,i,iv) = mat(1,3,4,i,iv) - c4 * dtl(i,iv) * Ahi(iv,i,2)
        mat(1,4,2,i,iv) = mat(1,4,2,i,iv) - c4 * dtl(i,iv) * Ahi(iv,i,1)
        mat(1,4,3,i,iv) = mat(1,4,3,i,iv) - c4 * dtl(i,iv) * Ahi(iv,i,2)

        mat(2,2,4,i,iv) = mat(2,2,4,i,iv) - c3 * dtl(i,iv) * Ahi(iv,i,1)
        mat(2,3,4,i,iv) = mat(2,3,4,i,iv) - c3 * dtl(i,iv) * Ahi(iv,i,2)
        mat(2,4,2,i,iv) = mat(2,4,2,i,iv) - c3 * dtl(i,iv) * Ahi(iv,i,1)
        mat(2,4,3,i,iv) = mat(2,4,3,i,iv) - c3 * dtl(i,iv) * Ahi(iv,i,2)

        mat(3,2,4,i,iv) = mat(3,2,4,i,iv) - c2 * dtl(i,iv) * Ahi(iv,i,1)
        mat(3,3,4,i,iv) = mat(3,3,4,i,iv) - c2 * dtl(i,iv) * Ahi(iv,i,2)
        mat(3,4,2,i,iv) = mat(3,4,2,i,iv) - c2 * dtl(i,iv) * Ahi(iv,i,1)
        mat(3,4,3,i,iv) = mat(3,4,3,i,iv) - c2 * dtl(i,iv) * Ahi(iv,i,2)

        mat(4,2,4,i,iv) = mat(4,2,4,i,iv) - c1 * dtl(i,iv) * Ahi(iv,i,1)
        mat(4,3,4,i,iv) = mat(4,3,4,i,iv) - c1 * dtl(i,iv) * Ahi(iv,i,2)
        mat(4,4,2,i,iv) = mat(4,4,2,i,iv) - c1 * dtl(i,iv) * Ahi(iv,i,1)
        mat(4,4,3,i,iv) = mat(4,4,3,i,iv) - c1 * dtl(i,iv) * Ahi(iv,i,2)
          
!.... \hat{V}_{\xi\xi} term

        a1 = alfa * db1 * dxisinv
        a2 = alfa * db2 * dxisinv
        a3 = alfa * db3 * dxisinv
        a4 = alfa * db4 * dxisinv
        a5 = alfa * db5 * dxisinv

        mat(5,2,2,i,iv) = mat(5,2,2,i,iv) - a5 * dtl(i,iv) * Vh(1,i,iv)
        mat(5,2,3,i,iv) = mat(5,2,3,i,iv) - a5 * dtl(i,iv) * Vh(5,i,iv)
        mat(5,3,2,i,iv) = mat(5,3,2,i,iv) - a5 * dtl(i,iv) * Vh(6,i,iv)
        mat(5,3,3,i,iv) = mat(5,3,3,i,iv) - a5 * dtl(i,iv) * Vh(2,i,iv)
        mat(5,4,4,i,iv) = mat(5,4,4,i,iv) - a5 * dtl(i,iv) * Vh(3,i,iv)
        mat(5,5,5,i,iv) = mat(5,5,5,i,iv) - a5 * dtl(i,iv) * Vh(4,i,iv)

        mat(1,2,2,i,iv) = mat(1,2,2,i,iv) - a4 * dtl(i,iv) * Vh(1,i,iv)
        mat(1,2,3,i,iv) = mat(1,2,3,i,iv) - a4 * dtl(i,iv) * Vh(5,i,iv)
        mat(1,3,2,i,iv) = mat(1,3,2,i,iv) - a4 * dtl(i,iv) * Vh(6,i,iv)
        mat(1,3,3,i,iv) = mat(1,3,3,i,iv) - a4 * dtl(i,iv) * Vh(2,i,iv)
        mat(1,4,4,i,iv) = mat(1,4,4,i,iv) - a4 * dtl(i,iv) * Vh(3,i,iv)
        mat(1,5,5,i,iv) = mat(1,5,5,i,iv) - a4 * dtl(i,iv) * Vh(4,i,iv)

        mat(2,2,2,i,iv) = mat(2,2,2,i,iv) - a3 * dtl(i,iv) * Vh(1,i,iv)
        mat(2,2,3,i,iv) = mat(2,2,3,i,iv) - a3 * dtl(i,iv) * Vh(5,i,iv)
        mat(2,3,2,i,iv) = mat(2,3,2,i,iv) - a3 * dtl(i,iv) * Vh(6,i,iv)
        mat(2,3,3,i,iv) = mat(2,3,3,i,iv) - a3 * dtl(i,iv) * Vh(2,i,iv)
        mat(2,4,4,i,iv) = mat(2,4,4,i,iv) - a3 * dtl(i,iv) * Vh(3,i,iv)
        mat(2,5,5,i,iv) = mat(2,5,5,i,iv) - a3 * dtl(i,iv) * Vh(4,i,iv)

        mat(3,2,2,i,iv) = mat(3,2,2,i,iv) - a2 * dtl(i,iv) * Vh(1,i,iv)
        mat(3,2,3,i,iv) = mat(3,2,3,i,iv) - a2 * dtl(i,iv) * Vh(5,i,iv)
        mat(3,3,2,i,iv) = mat(3,3,2,i,iv) - a2 * dtl(i,iv) * Vh(6,i,iv)
        mat(3,3,3,i,iv) = mat(3,3,3,i,iv) - a2 * dtl(i,iv) * Vh(2,i,iv)
        mat(3,4,4,i,iv) = mat(3,4,4,i,iv) - a2 * dtl(i,iv) * Vh(3,i,iv)
        mat(3,5,5,i,iv) = mat(3,5,5,i,iv) - a2 * dtl(i,iv) * Vh(4,i,iv)

        mat(4,2,2,i,iv) = mat(4,2,2,i,iv) - a1 * dtl(i,iv) * Vh(1,i,iv)
        mat(4,2,3,i,iv) = mat(4,2,3,i,iv) - a1 * dtl(i,iv) * Vh(5,i,iv)
        mat(4,3,2,i,iv) = mat(4,3,2,i,iv) - a1 * dtl(i,iv) * Vh(6,i,iv)
        mat(4,3,3,i,iv) = mat(4,3,3,i,iv) - a1 * dtl(i,iv) * Vh(2,i,iv)
        mat(4,4,4,i,iv) = mat(4,4,4,i,iv) - a1 * dtl(i,iv) * Vh(3,i,iv)
        mat(4,5,5,i,iv) = mat(4,5,5,i,iv) - a1 * dtl(i,iv) * Vh(4,i,iv)

!.... I term

        mat(3,1,1,i,iv) = mat(3,1,1,i,iv) + one
        mat(3,2,2,i,iv) = mat(3,2,2,i,iv) + one
        mat(3,3,3,i,iv) = mat(3,3,3,i,iv) + one
        mat(3,4,4,i,iv) = mat(3,4,4,i,iv) + one
        mat(3,5,5,i,iv) = mat(3,5,5,i,iv) + one

!.... implicit damping term

        if (eps_e .ne. zero) then

        mat(2,1,1,i,iv) = mat(2,1,1,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (-one)
        mat(2,2,2,i,iv) = mat(2,2,2,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (-one)
        mat(2,3,3,i,iv) = mat(2,3,3,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (-one)
        mat(2,4,4,i,iv) = mat(2,4,4,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (-one)
        mat(2,5,5,i,iv) = mat(2,5,5,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (-one)

        mat(3,1,1,i,iv) = mat(3,1,1,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * two
        mat(3,2,2,i,iv) = mat(3,2,2,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * two
        mat(3,3,3,i,iv) = mat(3,3,3,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * two
        mat(3,4,4,i,iv) = mat(3,4,4,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * two
        mat(3,5,5,i,iv) = mat(3,5,5,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * two

        mat(4,1,1,i,iv) = mat(4,1,1,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (-one)
        mat(4,2,2,i,iv) = mat(4,2,2,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (-one)
        mat(4,3,3,i,iv) = mat(4,3,3,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (-one)
        mat(4,4,4,i,iv) = mat(4,4,4,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (-one)
        mat(4,5,5,i,iv) = mat(4,5,5,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (-one)

        end if

!.... sponge term

        if (.not. calcd) then

        if (ispg .eq. 1) then
        
            mat(3,1,1,i,iv) = mat(3,1,1,i,iv) + alfa * dtl(i,iv) * spgl(iv,i)
            mat(3,2,2,i,iv) = mat(3,2,2,i,iv) + alfa * dtl(i,iv) * spgl(iv,i)
            mat(3,3,3,i,iv) = mat(3,3,3,i,iv) + alfa * dtl(i,iv) * spgl(iv,i)
            mat(3,4,4,i,iv) = mat(3,4,4,i,iv) + alfa * dtl(i,iv) * spgl(iv,i)
            mat(3,5,5,i,iv) = mat(3,5,5,i,iv) + alfa * dtl(i,iv) * spgl(iv,i)

        else if (ispg .ge. 2) then
        
            mat(3,1,1,i,iv) = mat(3,1,1,i,iv) + alfa * dtl(i,iv) * (spgl(iv,i) + spg2l(iv,i))
            mat(3,2,2,i,iv) = mat(3,2,2,i,iv) + alfa * dtl(i,iv) * (spgl(iv,i) + spg2l(iv,i))
            mat(3,3,3,i,iv) = mat(3,3,3,i,iv) + alfa * dtl(i,iv) * (spgl(iv,i) + spg2l(iv,i))
            mat(3,4,4,i,iv) = mat(3,4,4,i,iv) + alfa * dtl(i,iv) * (spgl(iv,i) + spg2l(iv,i))
            mat(3,5,5,i,iv) = mat(3,5,5,i,iv) + alfa * dtl(i,iv) * (spgl(iv,i) + spg2l(iv,i))
        
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
              mat(4,idof,jdof,i,iv) = a5 * dtl(i,iv) * Ah(idof,jdof,i,iv)
              mat(5,idof,jdof,i,iv) = a4 * dtl(i,iv) * Ah(idof,jdof,i,iv)
              mat(1,idof,jdof,i,iv) = a3 * dtl(i,iv) * Ah(idof,jdof,i,iv)
              mat(2,idof,jdof,i,iv) = a2 * dtl(i,iv) * Ah(idof,jdof,i,iv)
              mat(3,idof,jdof,i,iv) = dtl(i,iv) * ( a1 * Ah(idof,jdof,i,iv) + &
                                      rcalcd * alfa * (  Dh(idof,jdof,i,iv) + &
                                                   im * Dhi(idof,jdof,i,iv) ) )
          end do
        end do
        
!.... \hat{A}_i term

        c1 = -im * alfa * gc1 * dxiinv
        c2 = -im * alfa * gc2 * dxiinv
        c3 = -im * alfa * gc3 * dxiinv
        c4 = -im * alfa * gc4 * dxiinv
        c5 = -im * alfa * gc5 * dxiinv

        mat(4,2,4,i,iv) = mat(4,2,4,i,iv) - c5 * dtl(i,iv) * Ahi(iv,i,1)
        mat(4,3,4,i,iv) = mat(4,3,4,i,iv) - c5 * dtl(i,iv) * Ahi(iv,i,2)
        mat(4,4,2,i,iv) = mat(4,4,2,i,iv) - c5 * dtl(i,iv) * Ahi(iv,i,1)
        mat(4,4,3,i,iv) = mat(4,4,3,i,iv) - c5 * dtl(i,iv) * Ahi(iv,i,2)

        mat(5,2,4,i,iv) = mat(5,2,4,i,iv) - c4 * dtl(i,iv) * Ahi(iv,i,1)
        mat(5,3,4,i,iv) = mat(5,3,4,i,iv) - c4 * dtl(i,iv) * Ahi(iv,i,2)
        mat(5,4,2,i,iv) = mat(5,4,2,i,iv) - c4 * dtl(i,iv) * Ahi(iv,i,1)
        mat(5,4,3,i,iv) = mat(5,4,3,i,iv) - c4 * dtl(i,iv) * Ahi(iv,i,2)

        mat(1,2,4,i,iv) = mat(1,2,4,i,iv) - c3 * dtl(i,iv) * Ahi(iv,i,1)
        mat(1,3,4,i,iv) = mat(1,3,4,i,iv) - c3 * dtl(i,iv) * Ahi(iv,i,2)
        mat(1,4,2,i,iv) = mat(1,4,2,i,iv) - c3 * dtl(i,iv) * Ahi(iv,i,1)
        mat(1,4,3,i,iv) = mat(1,4,3,i,iv) - c3 * dtl(i,iv) * Ahi(iv,i,2)

        mat(2,2,4,i,iv) = mat(2,2,4,i,iv) - c2 * dtl(i,iv) * Ahi(iv,i,1)
        mat(2,3,4,i,iv) = mat(2,3,4,i,iv) - c2 * dtl(i,iv) * Ahi(iv,i,2)
        mat(2,4,2,i,iv) = mat(2,4,2,i,iv) - c2 * dtl(i,iv) * Ahi(iv,i,1)
        mat(2,4,3,i,iv) = mat(2,4,3,i,iv) - c2 * dtl(i,iv) * Ahi(iv,i,2)

        mat(3,2,4,i,iv) = mat(3,2,4,i,iv) - c1 * dtl(i,iv) * Ahi(iv,i,1)
        mat(3,3,4,i,iv) = mat(3,3,4,i,iv) - c1 * dtl(i,iv) * Ahi(iv,i,2)
        mat(3,4,2,i,iv) = mat(3,4,2,i,iv) - c1 * dtl(i,iv) * Ahi(iv,i,1)
        mat(3,4,3,i,iv) = mat(3,4,3,i,iv) - c1 * dtl(i,iv) * Ahi(iv,i,2)

!.... \hat{V}_{\xi\xi} term

        a1 = alfa * dd1 * dxisinv
        a2 = alfa * dd2 * dxisinv
        a3 = alfa * dd3 * dxisinv
        a4 = alfa * dd4 * dxisinv
        a5 = alfa * dd5 * dxisinv

        mat(4,2,2,i,iv) = mat(4,2,2,i,iv) - a5 * dtl(i,iv) * Vh(1,i,iv)
        mat(4,2,3,i,iv) = mat(4,2,3,i,iv) - a5 * dtl(i,iv) * Vh(5,i,iv)
        mat(4,3,2,i,iv) = mat(4,3,2,i,iv) - a5 * dtl(i,iv) * Vh(6,i,iv)
        mat(4,3,3,i,iv) = mat(4,3,3,i,iv) - a5 * dtl(i,iv) * Vh(2,i,iv)
        mat(4,4,4,i,iv) = mat(4,4,4,i,iv) - a5 * dtl(i,iv) * Vh(3,i,iv)
        mat(4,5,5,i,iv) = mat(4,5,5,i,iv) - a5 * dtl(i,iv) * Vh(4,i,iv)

        mat(5,2,2,i,iv) = mat(5,2,2,i,iv) - a4 * dtl(i,iv) * Vh(1,i,iv)
        mat(5,2,3,i,iv) = mat(5,2,3,i,iv) - a4 * dtl(i,iv) * Vh(5,i,iv)
        mat(5,3,2,i,iv) = mat(5,3,2,i,iv) - a4 * dtl(i,iv) * Vh(6,i,iv)
        mat(5,3,3,i,iv) = mat(5,3,3,i,iv) - a4 * dtl(i,iv) * Vh(2,i,iv)
        mat(5,4,4,i,iv) = mat(5,4,4,i,iv) - a4 * dtl(i,iv) * Vh(3,i,iv)
        mat(5,5,5,i,iv) = mat(5,5,5,i,iv) - a4 * dtl(i,iv) * Vh(4,i,iv)

        mat(1,2,2,i,iv) = mat(1,2,2,i,iv) - a3 * dtl(i,iv) * Vh(1,i,iv)
        mat(1,2,3,i,iv) = mat(1,2,3,i,iv) - a3 * dtl(i,iv) * Vh(5,i,iv)
        mat(1,3,2,i,iv) = mat(1,3,2,i,iv) - a3 * dtl(i,iv) * Vh(6,i,iv)
        mat(1,3,3,i,iv) = mat(1,3,3,i,iv) - a3 * dtl(i,iv) * Vh(2,i,iv)
        mat(1,4,4,i,iv) = mat(1,4,4,i,iv) - a3 * dtl(i,iv) * Vh(3,i,iv)
        mat(1,5,5,i,iv) = mat(1,5,5,i,iv) - a3 * dtl(i,iv) * Vh(4,i,iv)

        mat(2,2,2,i,iv) = mat(2,2,2,i,iv) - a2 * dtl(i,iv) * Vh(1,i,iv)
        mat(2,2,3,i,iv) = mat(2,2,3,i,iv) - a2 * dtl(i,iv) * Vh(5,i,iv)
        mat(2,3,2,i,iv) = mat(2,3,2,i,iv) - a2 * dtl(i,iv) * Vh(6,i,iv)
        mat(2,3,3,i,iv) = mat(2,3,3,i,iv) - a2 * dtl(i,iv) * Vh(2,i,iv)
        mat(2,4,4,i,iv) = mat(2,4,4,i,iv) - a2 * dtl(i,iv) * Vh(3,i,iv)
        mat(2,5,5,i,iv) = mat(2,5,5,i,iv) - a2 * dtl(i,iv) * Vh(4,i,iv)

        mat(3,2,2,i,iv) = mat(3,2,2,i,iv) - a1 * dtl(i,iv) * Vh(1,i,iv)
        mat(3,2,3,i,iv) = mat(3,2,3,i,iv) - a1 * dtl(i,iv) * Vh(5,i,iv)
        mat(3,3,2,i,iv) = mat(3,3,2,i,iv) - a1 * dtl(i,iv) * Vh(6,i,iv)
        mat(3,3,3,i,iv) = mat(3,3,3,i,iv) - a1 * dtl(i,iv) * Vh(2,i,iv)
        mat(3,4,4,i,iv) = mat(3,4,4,i,iv) - a1 * dtl(i,iv) * Vh(3,i,iv)
        mat(3,5,5,i,iv) = mat(3,5,5,i,iv) - a1 * dtl(i,iv) * Vh(4,i,iv)

!.... I term

        mat(3,1,1,i,iv) = mat(3,1,1,i,iv) + one
        mat(3,2,2,i,iv) = mat(3,2,2,i,iv) + one
        mat(3,3,3,i,iv) = mat(3,3,3,i,iv) + one
        mat(3,4,4,i,iv) = mat(3,4,4,i,iv) + one
        mat(3,5,5,i,iv) = mat(3,5,5,i,iv) + one

!.... sponge term

        if (.not. calcd) then

        if (ispg .eq. 1) then
        
            mat(3,1,1,i,iv) = mat(3,1,1,i,iv) + alfa * dtl(i,iv) * spgl(iv,i)
            mat(3,2,2,i,iv) = mat(3,2,2,i,iv) + alfa * dtl(i,iv) * spgl(iv,i)
            mat(3,3,3,i,iv) = mat(3,3,3,i,iv) + alfa * dtl(i,iv) * spgl(iv,i)
            mat(3,4,4,i,iv) = mat(3,4,4,i,iv) + alfa * dtl(i,iv) * spgl(iv,i)
            mat(3,5,5,i,iv) = mat(3,5,5,i,iv) + alfa * dtl(i,iv) * spgl(iv,i)

        else if (ispg .ge. 2) then
        
            mat(3,1,1,i,iv) = mat(3,1,1,i,iv) + alfa * dtl(i,iv) * (spgl(iv,i) + spg2l(iv,i))
            mat(3,2,2,i,iv) = mat(3,2,2,i,iv) + alfa * dtl(i,iv) * (spgl(iv,i) + spg2l(iv,i))
            mat(3,3,3,i,iv) = mat(3,3,3,i,iv) + alfa * dtl(i,iv) * (spgl(iv,i) + spg2l(iv,i))
            mat(3,4,4,i,iv) = mat(3,4,4,i,iv) + alfa * dtl(i,iv) * (spgl(iv,i) + spg2l(iv,i))
            mat(3,5,5,i,iv) = mat(3,5,5,i,iv) + alfa * dtl(i,iv) * (spgl(iv,i) + spg2l(iv,i))
        
        end if

        end if

        end do

        end if   ! rsym

        return
        end
