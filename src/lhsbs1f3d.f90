!=======================================================================================================!
        subroutine lhsbs1f3D( mat, Ah, Dh, Dhi, Vh, Ahi, spgl, spg2l, dtl, calcd, side )
!  
!  Correct the LHS for symmetry conditions in the \xi direction.
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
        integer :: side

        real    :: a1, a2, a3, a4, a5
        complex :: c1, c2, c3, c4, c5

        real :: dxiinv, dxisinv

        integer :: iv, i, idof, jdof

        real :: isign, rcalcd
!=======================================================================================================!
        rcalcd = zero
        if (calcd) rcalcd = one

        dxiinv  = one / dxi
        dxisinv = one / dxi**2

        if (side .eq. 1) then

          write(*,*) js, je

        !$omp parallel do private (i,idof,jdof,a1,a2,a3,a4,a5,c1,c2,c3,c4,c5,isign)
        do iv = js, je
!=======================================================================================================!
!.... symmetric condition on left node
!=======================================================================================================!
        i = 1

        a1 = alfa * ga1 * dxiinv
        a2 = alfa * ga2 * dxiinv
        a3 = zero
        a4 = alfa * ga3 * dxiinv
        a5 = alfa * ga4 * dxiinv

        do idof = 1, ndof
          do jdof = 1, ndof
            if (jdof.eq.3) then
              isign = -one
            else
              isign = one
            end if
              mat(3,idof,jdof,i,iv) = rcalcd * alfa * dtl(i,iv) * ( Dh(idof,jdof,i,iv) + &
                                                              im * Dhi(idof,jdof,i,iv) )
              mat(4,idof,jdof,i,iv) = dtl(i,iv) * (a4+isign*a2) * Ah(idof,jdof,i,iv)
              mat(5,idof,jdof,i,iv) = dtl(i,iv) * (a5+isign*a1) * Ah(idof,jdof,i,iv)
              mat(1,idof,jdof,i,iv) = zero
              mat(2,idof,jdof,i,iv) = zero
          end do
        end do
        
!.... \hat{A}_i term

        c1 = im * alfa * ga1 * dxiinv
        c2 = im * alfa * ga2 * dxiinv
        c3 = zero
        c4 = im * alfa * ga3 * dxiinv
        c5 = im * alfa * ga4 * dxiinv

        mat(3,2,4,i,iv) = mat(3,2,4,i,iv) - c3 * dtl(i,iv) * Ahi(iv,i,1)
        mat(3,3,4,i,iv) = mat(3,3,4,i,iv) - c3 * dtl(i,iv) * Ahi(iv,i,2)
        mat(3,4,2,i,iv) = mat(3,4,2,i,iv) - c3 * dtl(i,iv) * Ahi(iv,i,1)
        mat(3,4,3,i,iv) = mat(3,4,3,i,iv) - c3 * dtl(i,iv) * Ahi(iv,i,2)

        mat(4,2,4,i,iv) = mat(4,2,4,i,iv) - (c4+c2) * dtl(i,iv) * Ahi(iv,i,1)
        mat(4,3,4,i,iv) = mat(4,3,4,i,iv) - (c4+c2) * dtl(i,iv) * Ahi(iv,i,2)
        mat(4,4,2,i,iv) = mat(4,4,2,i,iv) - (c4+c2) * dtl(i,iv) * Ahi(iv,i,1)
        mat(4,4,3,i,iv) = mat(4,4,3,i,iv) - (c4-c2) * dtl(i,iv) * Ahi(iv,i,2)

        mat(5,2,4,i,iv) = mat(5,2,4,i,iv) - (c5+c1) * dtl(i,iv) * Ahi(iv,i,1)
        mat(5,3,4,i,iv) = mat(5,3,4,i,iv) - (c5+c1) * dtl(i,iv) * Ahi(iv,i,2)
        mat(5,4,2,i,iv) = mat(5,4,2,i,iv) - (c5+c1) * dtl(i,iv) * Ahi(iv,i,1)
        mat(5,4,3,i,iv) = mat(5,4,3,i,iv) - (c5-c1) * dtl(i,iv) * Ahi(iv,i,2)
          
!.... \hat{V}_{\xi\xi} term

        a1 = alfa * da1 * dxisinv
        a2 = alfa * da2 * dxisinv
        a3 = alfa * da3 * dxisinv
        a4 = alfa * da4 * dxisinv
        a5 = alfa * da5 * dxisinv

        mat(3,2,2,i,iv) = mat(3,2,2,i,iv) - dtl(i,iv) * a3 * Vh(1,i,iv)
        mat(3,2,3,i,iv) = mat(3,2,3,i,iv) - dtl(i,iv) * a3 * Vh(5,i,iv)
        mat(3,3,2,i,iv) = mat(3,3,2,i,iv) - dtl(i,iv) * a3 * Vh(6,i,iv)
        mat(3,3,3,i,iv) = mat(3,3,3,i,iv) - dtl(i,iv) * a3 * Vh(2,i,iv)
        mat(3,4,4,i,iv) = mat(3,4,4,i,iv) - dtl(i,iv) * a3 * Vh(3,i,iv)
        mat(3,5,5,i,iv) = mat(3,5,5,i,iv) - dtl(i,iv) * a3 * Vh(4,i,iv)

        mat(4,2,2,i,iv) = mat(4,2,2,i,iv) - dtl(i,iv) * (a4+a2) * Vh(1,i,iv)
        mat(4,2,3,i,iv) = mat(4,2,3,i,iv) - dtl(i,iv) * (a4-a2) * Vh(5,i,iv)
        mat(4,3,2,i,iv) = mat(4,3,2,i,iv) - dtl(i,iv) * (a4+a2) * Vh(6,i,iv)
        mat(4,3,3,i,iv) = mat(4,3,3,i,iv) - dtl(i,iv) * (a4-a2) * Vh(2,i,iv)
        mat(4,4,4,i,iv) = mat(4,4,4,i,iv) - dtl(i,iv) * (a4+a2) * Vh(3,i,iv)
        mat(4,5,5,i,iv) = mat(4,5,5,i,iv) - dtl(i,iv) * (a4+a2) * Vh(4,i,iv)

        mat(5,2,2,i,iv) = mat(5,2,2,i,iv) - dtl(i,iv) * (a5+a1) * Vh(1,i,iv)
        mat(5,2,3,i,iv) = mat(5,2,3,i,iv) - dtl(i,iv) * (a5-a1) * Vh(5,i,iv)
        mat(5,3,2,i,iv) = mat(5,3,2,i,iv) - dtl(i,iv) * (a5+a1) * Vh(6,i,iv)
        mat(5,3,3,i,iv) = mat(5,3,3,i,iv) - dtl(i,iv) * (a5-a1) * Vh(2,i,iv)
        mat(5,4,4,i,iv) = mat(5,4,4,i,iv) - dtl(i,iv) * (a5+a1) * Vh(3,i,iv)
        mat(5,5,5,i,iv) = mat(5,5,5,i,iv) - dtl(i,iv) * (a5+a1) * Vh(4,i,iv)

!.... I term

        mat(3,1,1,i,iv) = mat(3,1,1,i,iv) + one
        mat(3,2,2,i,iv) = mat(3,2,2,i,iv) + one
        mat(3,3,3,i,iv) = mat(3,3,3,i,iv) + one
        mat(3,4,4,i,iv) = mat(3,4,4,i,iv) + one
        mat(3,5,5,i,iv) = mat(3,5,5,i,iv) + one

!.... implicit damping term

        if (eps_e .ne. zero) then

        mat(3,1,1,i,iv) = mat(3,1,1,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb3
        mat(3,2,2,i,iv) = mat(3,2,2,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb3
        mat(3,3,3,i,iv) = mat(3,3,3,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb3
        mat(3,4,4,i,iv) = mat(3,4,4,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb3
        mat(3,5,5,i,iv) = mat(3,5,5,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb3

        mat(4,1,1,i,iv) = mat(4,1,1,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (fb4+fb2)
        mat(4,2,2,i,iv) = mat(4,2,2,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (fb4+fb2)
        mat(4,3,3,i,iv) = mat(4,3,3,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (fb4-fb2)
        mat(4,4,4,i,iv) = mat(4,4,4,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (fb4+fb2)
        mat(4,5,5,i,iv) = mat(4,5,5,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (fb4+fb2)

        mat(5,1,1,i,iv) = mat(5,1,1,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (fb5+fb1)
        mat(5,2,2,i,iv) = mat(5,2,2,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (fb5+fb1)
        mat(5,3,3,i,iv) = mat(5,3,3,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (fb5-fb1)
        mat(5,4,4,i,iv) = mat(5,4,4,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (fb5+fb1)
        mat(5,5,5,i,iv) = mat(5,5,5,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (fb5+fb1)

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
!.... symmetric condition on first node off the left boundary
!=======================================================================================================!
        i = 2

        a1 = alfa * ga1 * dxiinv
        a2 = alfa * ga2 * dxiinv
        a3 = zero
        a4 = alfa * ga3 * dxiinv
        a5 = alfa * ga4 * dxiinv

        do idof = 1, ndof
          do jdof = 1, ndof
            if (jdof.eq.3) then
              isign = -one
            else
              isign = one
            end if
              mat(2,idof,jdof,i,iv) = dtl(i,iv) * a2 * Ah(idof,jdof,i,iv)
              mat(3,idof,jdof,i,iv) = dtl(i,iv) * isign*a1 * Ah(idof,jdof,i,iv) + &
                                      rcalcd * alfa * dtl(i,iv) * ( Dh(idof,jdof,i,iv) + &
                                                              im * Dhi(idof,jdof,i,iv) )
              mat(4,idof,jdof,i,iv) = dtl(i,iv) * a4 * Ah(idof,jdof,i,iv)
              mat(5,idof,jdof,i,iv) = dtl(i,iv) * a5 * Ah(idof,jdof,i,iv)
              mat(1,idof,jdof,i,iv) = zero
          end do
        end do
        
!.... \hat{A}_i term

        c1 = im * alfa * ga1 * dxiinv
        c2 = im * alfa * ga2 * dxiinv
        c3 = im * alfa * zero * dxiinv
        c4 = im * alfa * ga3 * dxiinv
        c5 = im * alfa * ga4 * dxiinv

        mat(2,2,4,i,iv) = mat(2,2,4,i,iv) - c2 * dtl(i,iv) * Ahi(iv,i,1)
        mat(2,3,4,i,iv) = mat(2,3,4,i,iv) - c2 * dtl(i,iv) * Ahi(iv,i,2)
        mat(2,4,2,i,iv) = mat(2,4,2,i,iv) - c2 * dtl(i,iv) * Ahi(iv,i,1)
        mat(2,4,3,i,iv) = mat(2,4,3,i,iv) - c2 * dtl(i,iv) * Ahi(iv,i,2)

        mat(3,2,4,i,iv) = mat(3,2,4,i,iv) - (c3+c1) * dtl(i,iv) * Ahi(iv,i,1)
        mat(3,3,4,i,iv) = mat(3,3,4,i,iv) - (c3+c1) * dtl(i,iv) * Ahi(iv,i,2)
        mat(3,4,2,i,iv) = mat(3,4,2,i,iv) - (c3+c1) * dtl(i,iv) * Ahi(iv,i,1)
        mat(3,4,3,i,iv) = mat(3,4,3,i,iv) - (c3-c1) * dtl(i,iv) * Ahi(iv,i,2)

        mat(4,2,4,i,iv) = mat(4,2,4,i,iv) - c4 * dtl(i,iv) * Ahi(iv,i,1)
        mat(4,3,4,i,iv) = mat(4,3,4,i,iv) - c4 * dtl(i,iv) * Ahi(iv,i,2)
        mat(4,4,2,i,iv) = mat(4,4,2,i,iv) - c4 * dtl(i,iv) * Ahi(iv,i,1)
        mat(4,4,3,i,iv) = mat(4,4,3,i,iv) - c4 * dtl(i,iv) * Ahi(iv,i,2)

        mat(5,2,4,i,iv) = mat(5,2,4,i,iv) - c5 * dtl(i,iv) * Ahi(iv,i,1)
        mat(5,3,4,i,iv) = mat(5,3,4,i,iv) - c5 * dtl(i,iv) * Ahi(iv,i,2)
        mat(5,4,2,i,iv) = mat(5,4,2,i,iv) - c5 * dtl(i,iv) * Ahi(iv,i,1)
        mat(5,4,3,i,iv) = mat(5,4,3,i,iv) - c5 * dtl(i,iv) * Ahi(iv,i,2)
          
!.... \hat{V}_{\xi\xi} term

        a1 = alfa * da1 * dxisinv
        a2 = alfa * da2 * dxisinv
        a3 = alfa * da3 * dxisinv
        a4 = alfa * da4 * dxisinv
        a5 = alfa * da5 * dxisinv

        mat(2,2,2,i,iv) = mat(2,2,2,i,iv) - dtl(i,iv) * a2 * Vh(1,i,iv)
        mat(2,2,3,i,iv) = mat(2,2,3,i,iv) - dtl(i,iv) * a2 * Vh(5,i,iv)
        mat(2,3,2,i,iv) = mat(2,3,2,i,iv) - dtl(i,iv) * a2 * Vh(6,i,iv)
        mat(2,3,3,i,iv) = mat(2,3,3,i,iv) - dtl(i,iv) * a2 * Vh(2,i,iv)
        mat(2,4,4,i,iv) = mat(2,4,4,i,iv) - dtl(i,iv) * a2 * Vh(3,i,iv)
        mat(2,5,5,i,iv) = mat(2,5,5,i,iv) - dtl(i,iv) * a2 * Vh(4,i,iv)

        mat(3,2,2,i,iv) = mat(3,2,2,i,iv) - dtl(i,iv) * (a3+a1) * Vh(1,i,iv)
        mat(3,2,3,i,iv) = mat(3,2,3,i,iv) - dtl(i,iv) * (a3-a1) * Vh(5,i,iv)
        mat(3,3,2,i,iv) = mat(3,3,2,i,iv) - dtl(i,iv) * (a3+a1) * Vh(6,i,iv)
        mat(3,3,3,i,iv) = mat(3,3,3,i,iv) - dtl(i,iv) * (a3-a1) * Vh(2,i,iv)
        mat(3,4,4,i,iv) = mat(3,4,4,i,iv) - dtl(i,iv) * (a3+a1) * Vh(3,i,iv)
        mat(3,5,5,i,iv) = mat(3,5,5,i,iv) - dtl(i,iv) * (a3+a1) * Vh(4,i,iv)

        mat(4,2,2,i,iv) = mat(4,2,2,i,iv) - dtl(i,iv) * a4 * Vh(1,i,iv)
        mat(4,2,3,i,iv) = mat(4,2,3,i,iv) - dtl(i,iv) * a4 * Vh(5,i,iv)
        mat(4,3,2,i,iv) = mat(4,3,2,i,iv) - dtl(i,iv) * a4 * Vh(6,i,iv)
        mat(4,3,3,i,iv) = mat(4,3,3,i,iv) - dtl(i,iv) * a4 * Vh(2,i,iv)
        mat(4,4,4,i,iv) = mat(4,4,4,i,iv) - dtl(i,iv) * a4 * Vh(3,i,iv)
        mat(4,5,5,i,iv) = mat(4,5,5,i,iv) - dtl(i,iv) * a4 * Vh(4,i,iv)

        mat(5,2,2,i,iv) = mat(5,2,2,i,iv) - dtl(i,iv) * a5 * Vh(1,i,iv)
        mat(5,2,3,i,iv) = mat(5,2,3,i,iv) - dtl(i,iv) * a5 * Vh(5,i,iv)
        mat(5,3,2,i,iv) = mat(5,3,2,i,iv) - dtl(i,iv) * a5 * Vh(6,i,iv)
        mat(5,3,3,i,iv) = mat(5,3,3,i,iv) - dtl(i,iv) * a5 * Vh(2,i,iv)
        mat(5,4,4,i,iv) = mat(5,4,4,i,iv) - dtl(i,iv) * a5 * Vh(3,i,iv)
        mat(5,5,5,i,iv) = mat(5,5,5,i,iv) - dtl(i,iv) * a5 * Vh(4,i,iv)

!.... I term

        mat(3,1,1,i,iv) = mat(3,1,1,i,iv) + one
        mat(3,2,2,i,iv) = mat(3,2,2,i,iv) + one
        mat(3,3,3,i,iv) = mat(3,3,3,i,iv) + one
        mat(3,4,4,i,iv) = mat(3,4,4,i,iv) + one
        mat(3,5,5,i,iv) = mat(3,5,5,i,iv) + one

!.... implicit damping term

        if (eps_e .ne. zero) then

        mat(2,1,1,i,iv) = mat(2,1,1,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb2
        mat(2,2,2,i,iv) = mat(2,2,2,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb2
        mat(2,3,3,i,iv) = mat(2,3,3,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb2
        mat(2,4,4,i,iv) = mat(2,4,4,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb2
        mat(2,5,5,i,iv) = mat(2,5,5,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb2

        mat(3,1,1,i,iv) = mat(3,1,1,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (fb3+fb1)
        mat(3,2,2,i,iv) = mat(3,2,2,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (fb3+fb1)
        mat(3,3,3,i,iv) = mat(3,3,3,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (fb3-fb1)
        mat(3,4,4,i,iv) = mat(3,4,4,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (fb3+fb1)
        mat(3,5,5,i,iv) = mat(3,5,5,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (fb3+fb1)

        mat(4,1,1,i,iv) = mat(4,1,1,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb4
        mat(4,2,2,i,iv) = mat(4,2,2,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb4
        mat(4,3,3,i,iv) = mat(4,3,3,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb4
        mat(4,4,4,i,iv) = mat(4,4,4,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb4
        mat(4,5,5,i,iv) = mat(4,5,5,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb4

        mat(5,1,1,i,iv) = mat(5,1,1,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb5
        mat(5,2,2,i,iv) = mat(5,2,2,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb5
        mat(5,3,3,i,iv) = mat(5,3,3,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb5
        mat(5,4,4,i,iv) = mat(5,4,4,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb5
        mat(5,5,5,i,iv) = mat(5,5,5,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb5

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

        else if (side .eq. 2) then

        !$omp parallel do private (i,idof,jdof,a1,a2,a3,a4,a5,c1,c2,c3,c4,c5,isign)
        do iv = js, je
!=======================================================================================================!
!.... symmetric condition on first node off the right boundary
!=======================================================================================================!
        i = nx - 1

        a1 = alfa * ga1 * dxiinv
        a2 = alfa * ga2 * dxiinv
        a3 = zero
        a4 = alfa * ga3 * dxiinv
        a5 = alfa * ga4 * dxiinv

        do idof = 1, ndof
          do jdof = 1, ndof
            if (jdof.eq.3) then
              isign = -one
            else
              isign = one
            end if
              mat(1,idof,jdof,i,iv) = dtl(i,iv) * a1 * Ah(idof,jdof,i,iv)
              mat(2,idof,jdof,i,iv) = dtl(i,iv) * a2 * Ah(idof,jdof,i,iv)
              mat(3,idof,jdof,i,iv) = dtl(i,iv) * isign*a5 * Ah(idof,jdof,i,iv) + &
                                      rcalcd * alfa * dtl(i,iv) * ( Dh(idof,jdof,i,iv) + &
                                                              im * Dhi(idof,jdof,i,iv) )
              mat(4,idof,jdof,i,iv) = dtl(i,iv) * a4 * Ah(idof,jdof,i,iv)
              mat(5,idof,jdof,i,iv) = zero
          end do
        end do
        
!.... \hat{A}_i term

        c1 = im * alfa * ga1 * dxiinv
        c2 = im * alfa * ga2 * dxiinv
        c3 = im * alfa * zero * dxiinv
        c4 = im * alfa * ga3 * dxiinv
        c5 = im * alfa * ga4 * dxiinv

        mat(1,2,4,i,iv) = mat(1,2,4,i,iv) - c1 * dtl(i,iv) * Ahi(iv,i,1)
        mat(1,3,4,i,iv) = mat(1,3,4,i,iv) - c1 * dtl(i,iv) * Ahi(iv,i,2)
        mat(1,4,2,i,iv) = mat(1,4,2,i,iv) - c1 * dtl(i,iv) * Ahi(iv,i,1)
        mat(1,4,3,i,iv) = mat(1,4,3,i,iv) - c1 * dtl(i,iv) * Ahi(iv,i,2)

        mat(2,2,4,i,iv) = mat(2,2,4,i,iv) - c2 * dtl(i,iv) * Ahi(iv,i,1)
        mat(2,3,4,i,iv) = mat(2,3,4,i,iv) - c2 * dtl(i,iv) * Ahi(iv,i,2)
        mat(2,4,2,i,iv) = mat(2,4,2,i,iv) - c2 * dtl(i,iv) * Ahi(iv,i,1)
        mat(2,4,3,i,iv) = mat(2,4,3,i,iv) - c2 * dtl(i,iv) * Ahi(iv,i,2)

        mat(3,2,4,i,iv) = mat(3,2,4,i,iv) - (c3+c5) * dtl(i,iv) * Ahi(iv,i,1)
        mat(3,3,4,i,iv) = mat(3,3,4,i,iv) - (c3+c5) * dtl(i,iv) * Ahi(iv,i,2)
        mat(3,4,2,i,iv) = mat(3,4,2,i,iv) - (c3+c5) * dtl(i,iv) * Ahi(iv,i,1)
        mat(3,4,3,i,iv) = mat(3,4,3,i,iv) - (c3-c5) * dtl(i,iv) * Ahi(iv,i,2)

        mat(4,2,4,i,iv) = mat(4,2,4,i,iv) - c4 * dtl(i,iv) * Ahi(iv,i,1)
        mat(4,3,4,i,iv) = mat(4,3,4,i,iv) - c4 * dtl(i,iv) * Ahi(iv,i,2)
        mat(4,4,2,i,iv) = mat(4,4,2,i,iv) - c4 * dtl(i,iv) * Ahi(iv,i,1)
        mat(4,4,3,i,iv) = mat(4,4,3,i,iv) - c4 * dtl(i,iv) * Ahi(iv,i,2)

!.... \hat{V}_{\xi\xi} term

        a1 = alfa * da1 * dxisinv
        a2 = alfa * da2 * dxisinv
        a3 = alfa * da3 * dxisinv
        a4 = alfa * da4 * dxisinv
        a5 = alfa * da5 * dxisinv

        mat(1,2,2,i,iv) = mat(1,2,2,i,iv) - dtl(i,iv) * a1 * Vh(1,i,iv)
        mat(1,2,3,i,iv) = mat(1,2,3,i,iv) - dtl(i,iv) * a1 * Vh(5,i,iv)
        mat(1,3,2,i,iv) = mat(1,3,2,i,iv) - dtl(i,iv) * a1 * Vh(6,i,iv)
        mat(1,3,3,i,iv) = mat(1,3,3,i,iv) - dtl(i,iv) * a1 * Vh(2,i,iv)
        mat(1,4,4,i,iv) = mat(1,4,4,i,iv) - dtl(i,iv) * a1 * Vh(3,i,iv)
        mat(1,5,5,i,iv) = mat(1,5,5,i,iv) - dtl(i,iv) * a1 * Vh(4,i,iv)

        mat(2,2,2,i,iv) = mat(2,2,2,i,iv) - dtl(i,iv) * a2 * Vh(1,i,iv)
        mat(2,2,3,i,iv) = mat(2,2,3,i,iv) - dtl(i,iv) * a2 * Vh(5,i,iv)
        mat(2,3,2,i,iv) = mat(2,3,2,i,iv) - dtl(i,iv) * a2 * Vh(6,i,iv)
        mat(2,3,3,i,iv) = mat(2,3,3,i,iv) - dtl(i,iv) * a2 * Vh(2,i,iv)
        mat(2,4,4,i,iv) = mat(2,4,4,i,iv) - dtl(i,iv) * a2 * Vh(3,i,iv)
        mat(2,5,5,i,iv) = mat(2,5,5,i,iv) - dtl(i,iv) * a2 * Vh(4,i,iv)

        mat(3,2,2,i,iv) = mat(3,2,2,i,iv) - dtl(i,iv) * (a3+a5) * Vh(1,i,iv)
        mat(3,2,3,i,iv) = mat(3,2,3,i,iv) - dtl(i,iv) * (a3-a5) * Vh(5,i,iv)
        mat(3,3,2,i,iv) = mat(3,3,2,i,iv) - dtl(i,iv) * (a3+a5) * Vh(6,i,iv)
        mat(3,3,3,i,iv) = mat(3,3,3,i,iv) - dtl(i,iv) * (a3-a5) * Vh(2,i,iv)
        mat(3,4,4,i,iv) = mat(3,4,4,i,iv) - dtl(i,iv) * (a3+a5) * Vh(3,i,iv)
        mat(3,5,5,i,iv) = mat(3,5,5,i,iv) - dtl(i,iv) * (a3+a5) * Vh(4,i,iv)

        mat(4,2,2,i,iv) = mat(4,2,2,i,iv) - dtl(i,iv) * a4 * Vh(1,i,iv)
        mat(4,2,3,i,iv) = mat(4,2,3,i,iv) - dtl(i,iv) * a4 * Vh(5,i,iv)
        mat(4,3,2,i,iv) = mat(4,3,2,i,iv) - dtl(i,iv) * a4 * Vh(6,i,iv)
        mat(4,3,3,i,iv) = mat(4,3,3,i,iv) - dtl(i,iv) * a4 * Vh(2,i,iv)
        mat(4,4,4,i,iv) = mat(4,4,4,i,iv) - dtl(i,iv) * a4 * Vh(3,i,iv)
        mat(4,5,5,i,iv) = mat(4,5,5,i,iv) - dtl(i,iv) * a4 * Vh(4,i,iv)

!.... I term

        mat(3,1,1,i,iv) = mat(3,1,1,i,iv) + one
        mat(3,2,2,i,iv) = mat(3,2,2,i,iv) + one
        mat(3,3,3,i,iv) = mat(3,3,3,i,iv) + one
        mat(3,4,4,i,iv) = mat(3,4,4,i,iv) + one
        mat(3,5,5,i,iv) = mat(3,5,5,i,iv) + one

!.... implicit damping term

        if (eps_e .ne. zero) then

        mat(1,1,1,i,iv) = mat(1,1,1,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb1
        mat(1,2,2,i,iv) = mat(1,2,2,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb1
        mat(1,3,3,i,iv) = mat(1,3,3,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb1
        mat(1,4,4,i,iv) = mat(1,4,4,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb1
        mat(1,5,5,i,iv) = mat(1,5,5,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb1

        mat(2,1,1,i,iv) = mat(2,1,1,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb2
        mat(2,2,2,i,iv) = mat(2,2,2,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb2
        mat(2,3,3,i,iv) = mat(2,3,3,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb2
        mat(2,4,4,i,iv) = mat(2,4,4,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb2
        mat(2,5,5,i,iv) = mat(2,5,5,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb2

        mat(3,1,1,i,iv) = mat(3,1,1,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (fb3+fb5)
        mat(3,2,2,i,iv) = mat(3,2,2,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (fb3+fb5)
        mat(3,3,3,i,iv) = mat(3,3,3,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (fb3-fb5)
        mat(3,4,4,i,iv) = mat(3,4,4,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (fb3+fb5)
        mat(3,5,5,i,iv) = mat(3,5,5,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (fb3+fb5)

        mat(4,1,1,i,iv) = mat(4,1,1,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb4
        mat(4,2,2,i,iv) = mat(4,2,2,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb4
        mat(4,3,3,i,iv) = mat(4,3,3,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb4
        mat(4,4,4,i,iv) = mat(4,4,4,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb4
        mat(4,5,5,i,iv) = mat(4,5,5,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb4

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
!.... symmetric condition on right node
!=======================================================================================================!
        i = nx

        a1 = alfa * ga1 * dxiinv
        a2 = alfa * ga2 * dxiinv
        a3 = zero
        a4 = alfa * ga3 * dxiinv
        a5 = alfa * ga4 * dxiinv

        do idof = 1, ndof
          do jdof = 1, ndof
            if (jdof.eq.3) then
              isign = -one
            else
              isign = one
            end if
              mat(1,idof,jdof,i,iv) = dtl(i,iv) * (a1+isign*a5) * Ah(idof,jdof,i,iv)
              mat(2,idof,jdof,i,iv) = dtl(i,iv) * (a2+isign*a4) * Ah(idof,jdof,i,iv)
              mat(3,idof,jdof,i,iv) = rcalcd * alfa * dtl(i,iv) * ( Dh(idof,jdof,i,iv) + &
                                                              im * Dhi(idof,jdof,i,iv) )
              mat(4,idof,jdof,i,iv) = zero
              mat(5,idof,jdof,i,iv) = zero
          end do
        end do
        
!.... \hat{A}_i term

        c1 = im * alfa * ga1 * dxiinv
        c2 = im * alfa * ga2 * dxiinv
        c3 = zero
        c4 = im * alfa * ga3 * dxiinv
        c5 = im * alfa * ga4 * dxiinv

        mat(1,2,4,i,iv) = mat(1,2,4,i,iv) - (c1+c5) * dtl(i,iv) * Ahi(iv,i,1)
        mat(1,3,4,i,iv) = mat(1,3,4,i,iv) - (c1+c5) * dtl(i,iv) * Ahi(iv,i,2)
        mat(1,4,2,i,iv) = mat(1,4,2,i,iv) - (c1+c5) * dtl(i,iv) * Ahi(iv,i,1)
        mat(1,4,3,i,iv) = mat(1,4,3,i,iv) - (c1-c5) * dtl(i,iv) * Ahi(iv,i,2)
          
        mat(2,2,4,i,iv) = mat(2,2,4,i,iv) - (c2+c4) * dtl(i,iv) * Ahi(iv,i,1)
        mat(2,3,4,i,iv) = mat(2,3,4,i,iv) - (c2+c4) * dtl(i,iv) * Ahi(iv,i,2)
        mat(2,4,2,i,iv) = mat(2,4,2,i,iv) - (c2+c4) * dtl(i,iv) * Ahi(iv,i,1)
        mat(2,4,3,i,iv) = mat(2,4,3,i,iv) - (c2-c4) * dtl(i,iv) * Ahi(iv,i,2)

        mat(3,2,4,i,iv) = mat(3,2,4,i,iv) - c3 * dtl(i,iv) * Ahi(iv,i,1)
        mat(3,3,4,i,iv) = mat(3,3,4,i,iv) - c3 * dtl(i,iv) * Ahi(iv,i,2)
        mat(3,4,2,i,iv) = mat(3,4,2,i,iv) - c3 * dtl(i,iv) * Ahi(iv,i,1)
        mat(3,4,3,i,iv) = mat(3,4,3,i,iv) - c3 * dtl(i,iv) * Ahi(iv,i,2)

!.... \hat{V}_{\xi\xi} term

        a1 = alfa * da1 * dxisinv
        a2 = alfa * da2 * dxisinv
        a3 = alfa * da3 * dxisinv
        a4 = alfa * da4 * dxisinv
        a5 = alfa * da5 * dxisinv

        mat(1,2,2,i,iv) = mat(1,2,2,i,iv) - dtl(i,iv) * (a1+a5) * Vh(1,i,iv)
        mat(1,2,3,i,iv) = mat(1,2,3,i,iv) - dtl(i,iv) * (a1-a5) * Vh(5,i,iv)
        mat(1,3,2,i,iv) = mat(1,3,2,i,iv) - dtl(i,iv) * (a1+a5) * Vh(6,i,iv)
        mat(1,3,3,i,iv) = mat(1,3,3,i,iv) - dtl(i,iv) * (a1-a5) * Vh(2,i,iv)
        mat(1,4,4,i,iv) = mat(1,4,4,i,iv) - dtl(i,iv) * (a1+a5) * Vh(3,i,iv)
        mat(1,5,5,i,iv) = mat(1,5,5,i,iv) - dtl(i,iv) * (a1+a5) * Vh(4,i,iv)

        mat(2,2,2,i,iv) = mat(2,2,2,i,iv) - dtl(i,iv) * (a2+a4) * Vh(1,i,iv)
        mat(2,2,3,i,iv) = mat(2,2,3,i,iv) - dtl(i,iv) * (a2-a4) * Vh(5,i,iv)
        mat(2,3,2,i,iv) = mat(2,3,2,i,iv) - dtl(i,iv) * (a2+a4) * Vh(6,i,iv)
        mat(2,3,3,i,iv) = mat(2,3,3,i,iv) - dtl(i,iv) * (a2-a4) * Vh(2,i,iv)
        mat(2,4,4,i,iv) = mat(2,4,4,i,iv) - dtl(i,iv) * (a2+a4) * Vh(3,i,iv)
        mat(2,5,5,i,iv) = mat(2,5,5,i,iv) - dtl(i,iv) * (a2+a4) * Vh(4,i,iv)

        mat(3,2,2,i,iv) = mat(3,2,2,i,iv) - dtl(i,iv) * a3 * Vh(1,i,iv)
        mat(3,2,3,i,iv) = mat(3,2,3,i,iv) - dtl(i,iv) * a3 * Vh(5,i,iv)
        mat(3,3,2,i,iv) = mat(3,3,2,i,iv) - dtl(i,iv) * a3 * Vh(6,i,iv)
        mat(3,3,3,i,iv) = mat(3,3,3,i,iv) - dtl(i,iv) * a3 * Vh(2,i,iv)
        mat(3,4,4,i,iv) = mat(3,4,4,i,iv) - dtl(i,iv) * a3 * Vh(3,i,iv)
        mat(3,5,5,i,iv) = mat(3,5,5,i,iv) - dtl(i,iv) * a3 * Vh(4,i,iv)

!.... I term

        mat(3,1,1,i,iv) = mat(3,1,1,i,iv) + one
        mat(3,2,2,i,iv) = mat(3,2,2,i,iv) + one
        mat(3,3,3,i,iv) = mat(3,3,3,i,iv) + one
        mat(3,4,4,i,iv) = mat(3,4,4,i,iv) + one
        mat(3,5,5,i,iv) = mat(3,5,5,i,iv) + one

!.... implicit damping term

        if (eps_e .ne. zero) then

        mat(1,1,1,i,iv) = mat(1,1,1,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (fb1+fb5)
        mat(1,2,2,i,iv) = mat(1,2,2,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (fb1+fb5)
        mat(1,3,3,i,iv) = mat(1,3,3,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (fb1-fb5)
        mat(1,4,4,i,iv) = mat(1,4,4,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (fb1+fb5)
        mat(1,5,5,i,iv) = mat(1,5,5,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (fb1+fb5)

        mat(2,1,1,i,iv) = mat(2,1,1,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (fb2+fb4)
        mat(2,2,2,i,iv) = mat(2,2,2,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (fb2+fb4)
        mat(2,3,3,i,iv) = mat(2,3,3,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (fb2-fb4)
        mat(2,4,4,i,iv) = mat(2,4,4,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (fb2+fb4)
        mat(2,5,5,i,iv) = mat(2,5,5,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (fb2+fb4)

        mat(3,1,1,i,iv) = mat(3,1,1,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb3
        mat(3,2,2,i,iv) = mat(3,2,2,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb3
        mat(3,3,3,i,iv) = mat(3,3,3,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb3
        mat(3,4,4,i,iv) = mat(3,4,4,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb3
        mat(3,5,5,i,iv) = mat(3,5,5,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb3

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

        end if         ! side  

        return
        end
