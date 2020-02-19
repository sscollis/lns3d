!=======================================================================================================!
        subroutine lhsbs1f3D( mat, Ah, Dh, Dhi, Vh, Ahi, spgl, spg2l, dtl, calcd, side )
!  
!  Correct the LHS for symmetry conditions in the \xi direction.
!
!  This version supports the 4th order LHS and three-dimensional flow
!
!  Revised: 4-22-96
!=======================================================================================================!
        use stuff
        use diff
        use buff_stuff
        implicit none
        
        complex :: mat(ny,nx,ndof,ndof,5)
        real    :: Ah(ny,nx,ndof,ndof), Dh(ny,nx,ndof,ndof)
        real    :: Dhi(ny,nx,ndof,ndof), Ahi(ny,nx,6)
        real    :: spgl(ny,nx), spg2l(ny,nx), Vh(ny,nx,6), dtl(ny,nx)
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
            do iv = 1, ny
              mat(iv,i,idof,jdof,3) = rcalcd * alfa * dtl(iv,i) * ( Dh(iv,i,idof,jdof) + &
                                                              im * Dhi(iv,i,idof,jdof) )
              mat(iv,i,idof,jdof,4) = dtl(iv,i) * (a4+isign*a2) * Ah(iv,i,idof,jdof)
              mat(iv,i,idof,jdof,5) = dtl(iv,i) * (a5+isign*a1) * Ah(iv,i,idof,jdof)
              mat(iv,i,idof,jdof,1) = zero
              mat(iv,i,idof,jdof,2) = zero
            end do
          end do
        end do
        
!.... \hat{A}_i term

        c1 = im * alfa * ga1 * dxiinv
        c2 = im * alfa * ga2 * dxiinv
        c3 = zero
        c4 = im * alfa * ga3 * dxiinv
        c5 = im * alfa * ga4 * dxiinv

        do iv = 1, ny

        mat(iv,i,2,4,3) = mat(iv,i,2,4,3) - c3 * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,3,4,3) = mat(iv,i,3,4,3) - c3 * dtl(iv,i) * Ahi(iv,i,2)
        mat(iv,i,4,2,3) = mat(iv,i,4,2,3) - c3 * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,4,3,3) = mat(iv,i,4,3,3) - c3 * dtl(iv,i) * Ahi(iv,i,2)

        mat(iv,i,2,4,4) = mat(iv,i,2,4,4) - (c4+c2) * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,3,4,4) = mat(iv,i,3,4,4) - (c4+c2) * dtl(iv,i) * Ahi(iv,i,2)
        mat(iv,i,4,2,4) = mat(iv,i,4,2,4) - (c4+c2) * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,4,3,4) = mat(iv,i,4,3,4) - (c4-c2) * dtl(iv,i) * Ahi(iv,i,2)

        mat(iv,i,2,4,5) = mat(iv,i,2,4,5) - (c5+c1) * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,3,4,5) = mat(iv,i,3,4,5) - (c5+c1) * dtl(iv,i) * Ahi(iv,i,2)
        mat(iv,i,4,2,5) = mat(iv,i,4,2,5) - (c5+c1) * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,4,3,5) = mat(iv,i,4,3,5) - (c5-c1) * dtl(iv,i) * Ahi(iv,i,2)
          
        end do

!.... \hat{V}_{\xi\xi} term

        a1 = alfa * da1 * dxisinv
        a2 = alfa * da2 * dxisinv
        a3 = alfa * da3 * dxisinv
        a4 = alfa * da4 * dxisinv
        a5 = alfa * da5 * dxisinv

        do iv = 1, ny
        
        mat(iv,i,2,2,3) = mat(iv,i,2,2,3) - dtl(iv,i) * a3 * Vh(iv,i,1)
        mat(iv,i,2,3,3) = mat(iv,i,2,3,3) - dtl(iv,i) * a3 * Vh(iv,i,5)
        mat(iv,i,3,2,3) = mat(iv,i,3,2,3) - dtl(iv,i) * a3 * Vh(iv,i,6)
        mat(iv,i,3,3,3) = mat(iv,i,3,3,3) - dtl(iv,i) * a3 * Vh(iv,i,2)
        mat(iv,i,4,4,3) = mat(iv,i,4,4,3) - dtl(iv,i) * a3 * Vh(iv,i,3)
        mat(iv,i,5,5,3) = mat(iv,i,5,5,3) - dtl(iv,i) * a3 * Vh(iv,i,4)

        mat(iv,i,2,2,4) = mat(iv,i,2,2,4) - dtl(iv,i) * (a4+a2) * Vh(iv,i,1)
        mat(iv,i,2,3,4) = mat(iv,i,2,3,4) - dtl(iv,i) * (a4-a2) * Vh(iv,i,5)
        mat(iv,i,3,2,4) = mat(iv,i,3,2,4) - dtl(iv,i) * (a4+a2) * Vh(iv,i,6)
        mat(iv,i,3,3,4) = mat(iv,i,3,3,4) - dtl(iv,i) * (a4-a2) * Vh(iv,i,2)
        mat(iv,i,4,4,4) = mat(iv,i,4,4,4) - dtl(iv,i) * (a4+a2) * Vh(iv,i,3)
        mat(iv,i,5,5,4) = mat(iv,i,5,5,4) - dtl(iv,i) * (a4+a2) * Vh(iv,i,4)

        mat(iv,i,2,2,5) = mat(iv,i,2,2,5) - dtl(iv,i) * (a5+a1) * Vh(iv,i,1)
        mat(iv,i,2,3,5) = mat(iv,i,2,3,5) - dtl(iv,i) * (a5-a1) * Vh(iv,i,5)
        mat(iv,i,3,2,5) = mat(iv,i,3,2,5) - dtl(iv,i) * (a5+a1) * Vh(iv,i,6)
        mat(iv,i,3,3,5) = mat(iv,i,3,3,5) - dtl(iv,i) * (a5-a1) * Vh(iv,i,2)
        mat(iv,i,4,4,5) = mat(iv,i,4,4,5) - dtl(iv,i) * (a5+a1) * Vh(iv,i,3)
        mat(iv,i,5,5,5) = mat(iv,i,5,5,5) - dtl(iv,i) * (a5+a1) * Vh(iv,i,4)

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
        
        mat(iv,i,1,1,3) = mat(iv,i,1,1,3) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * fb3
        mat(iv,i,2,2,3) = mat(iv,i,2,2,3) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * fb3
        mat(iv,i,3,3,3) = mat(iv,i,3,3,3) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * fb3
        mat(iv,i,4,4,3) = mat(iv,i,4,4,3) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * fb3
        mat(iv,i,5,5,3) = mat(iv,i,5,5,3) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * fb3

        mat(iv,i,1,1,4) = mat(iv,i,1,1,4) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * (fb4+fb2)
        mat(iv,i,2,2,4) = mat(iv,i,2,2,4) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * (fb4+fb2)
        mat(iv,i,3,3,4) = mat(iv,i,3,3,4) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * (fb4-fb2)
        mat(iv,i,4,4,4) = mat(iv,i,4,4,4) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * (fb4+fb2)
        mat(iv,i,5,5,4) = mat(iv,i,5,5,4) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * (fb4+fb2)

        mat(iv,i,1,1,5) = mat(iv,i,1,1,5) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * (fb5+fb1)
        mat(iv,i,2,2,5) = mat(iv,i,2,2,5) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * (fb5+fb1)
        mat(iv,i,3,3,5) = mat(iv,i,3,3,5) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * (fb5-fb1)
        mat(iv,i,4,4,5) = mat(iv,i,4,4,5) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * (fb5+fb1)
        mat(iv,i,5,5,5) = mat(iv,i,5,5,5) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * (fb5+fb1)

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
            do iv = 1, ny
              mat(iv,i,idof,jdof,2) = dtl(iv,i) * a2 * Ah(iv,i,idof,jdof)
              mat(iv,i,idof,jdof,3) = dtl(iv,i) * isign*a1 * Ah(iv,i,idof,jdof) + &
                                      rcalcd * alfa * dtl(iv,i) * ( Dh(iv,i,idof,jdof) + &
                                                              im * Dhi(iv,i,idof,jdof) )
              mat(iv,i,idof,jdof,4) = dtl(iv,i) * a4 * Ah(iv,i,idof,jdof)
              mat(iv,i,idof,jdof,5) = dtl(iv,i) * a5 * Ah(iv,i,idof,jdof)
              mat(iv,i,idof,jdof,1) = zero
            end do
          end do
        end do
        
!.... \hat{A}_i term

        c1 = im * alfa * ga1 * dxiinv
        c2 = im * alfa * ga2 * dxiinv
        c3 = im * alfa * zero * dxiinv
        c4 = im * alfa * ga3 * dxiinv
        c5 = im * alfa * ga4 * dxiinv

        do iv = 1, ny

        mat(iv,i,2,4,2) = mat(iv,i,2,4,2) - c2 * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,3,4,2) = mat(iv,i,3,4,2) - c2 * dtl(iv,i) * Ahi(iv,i,2)
        mat(iv,i,4,2,2) = mat(iv,i,4,2,2) - c2 * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,4,3,2) = mat(iv,i,4,3,2) - c2 * dtl(iv,i) * Ahi(iv,i,2)

        mat(iv,i,2,4,3) = mat(iv,i,2,4,3) - (c3+c1) * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,3,4,3) = mat(iv,i,3,4,3) - (c3+c1) * dtl(iv,i) * Ahi(iv,i,2)
        mat(iv,i,4,2,3) = mat(iv,i,4,2,3) - (c3+c1) * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,4,3,3) = mat(iv,i,4,3,3) - (c3-c1) * dtl(iv,i) * Ahi(iv,i,2)

        mat(iv,i,2,4,4) = mat(iv,i,2,4,4) - c4 * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,3,4,4) = mat(iv,i,3,4,4) - c4 * dtl(iv,i) * Ahi(iv,i,2)
        mat(iv,i,4,2,4) = mat(iv,i,4,2,4) - c4 * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,4,3,4) = mat(iv,i,4,3,4) - c4 * dtl(iv,i) * Ahi(iv,i,2)

        mat(iv,i,2,4,5) = mat(iv,i,2,4,5) - c5 * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,3,4,5) = mat(iv,i,3,4,5) - c5 * dtl(iv,i) * Ahi(iv,i,2)
        mat(iv,i,4,2,5) = mat(iv,i,4,2,5) - c5 * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,4,3,5) = mat(iv,i,4,3,5) - c5 * dtl(iv,i) * Ahi(iv,i,2)
          
        end do

!.... \hat{V}_{\xi\xi} term

        a1 = alfa * da1 * dxisinv
        a2 = alfa * da2 * dxisinv
        a3 = alfa * da3 * dxisinv
        a4 = alfa * da4 * dxisinv
        a5 = alfa * da5 * dxisinv

        do iv = 1, ny
        
        mat(iv,i,2,2,2) = mat(iv,i,2,2,2) - dtl(iv,i) * a2 * Vh(iv,i,1)
        mat(iv,i,2,3,2) = mat(iv,i,2,3,2) - dtl(iv,i) * a2 * Vh(iv,i,5)
        mat(iv,i,3,2,2) = mat(iv,i,3,2,2) - dtl(iv,i) * a2 * Vh(iv,i,6)
        mat(iv,i,3,3,2) = mat(iv,i,3,3,2) - dtl(iv,i) * a2 * Vh(iv,i,2)
        mat(iv,i,4,4,2) = mat(iv,i,4,4,2) - dtl(iv,i) * a2 * Vh(iv,i,3)
        mat(iv,i,5,5,2) = mat(iv,i,5,5,2) - dtl(iv,i) * a2 * Vh(iv,i,4)

        mat(iv,i,2,2,3) = mat(iv,i,2,2,3) - dtl(iv,i) * (a3+a1) * Vh(iv,i,1)
        mat(iv,i,2,3,3) = mat(iv,i,2,3,3) - dtl(iv,i) * (a3-a1) * Vh(iv,i,5)
        mat(iv,i,3,2,3) = mat(iv,i,3,2,3) - dtl(iv,i) * (a3+a1) * Vh(iv,i,6)
        mat(iv,i,3,3,3) = mat(iv,i,3,3,3) - dtl(iv,i) * (a3-a1) * Vh(iv,i,2)
        mat(iv,i,4,4,3) = mat(iv,i,4,4,3) - dtl(iv,i) * (a3+a1) * Vh(iv,i,3)
        mat(iv,i,5,5,3) = mat(iv,i,5,5,3) - dtl(iv,i) * (a3+a1) * Vh(iv,i,4)

        mat(iv,i,2,2,4) = mat(iv,i,2,2,4) - dtl(iv,i) * a4 * Vh(iv,i,1)
        mat(iv,i,2,3,4) = mat(iv,i,2,3,4) - dtl(iv,i) * a4 * Vh(iv,i,5)
        mat(iv,i,3,2,4) = mat(iv,i,3,2,4) - dtl(iv,i) * a4 * Vh(iv,i,6)
        mat(iv,i,3,3,4) = mat(iv,i,3,3,4) - dtl(iv,i) * a4 * Vh(iv,i,2)
        mat(iv,i,4,4,4) = mat(iv,i,4,4,4) - dtl(iv,i) * a4 * Vh(iv,i,3)
        mat(iv,i,5,5,4) = mat(iv,i,5,5,4) - dtl(iv,i) * a4 * Vh(iv,i,4)

        mat(iv,i,2,2,5) = mat(iv,i,2,2,5) - dtl(iv,i) * a5 * Vh(iv,i,1)
        mat(iv,i,2,3,5) = mat(iv,i,2,3,5) - dtl(iv,i) * a5 * Vh(iv,i,5)
        mat(iv,i,3,2,5) = mat(iv,i,3,2,5) - dtl(iv,i) * a5 * Vh(iv,i,6)
        mat(iv,i,3,3,5) = mat(iv,i,3,3,5) - dtl(iv,i) * a5 * Vh(iv,i,2)
        mat(iv,i,4,4,5) = mat(iv,i,4,4,5) - dtl(iv,i) * a5 * Vh(iv,i,3)
        mat(iv,i,5,5,5) = mat(iv,i,5,5,5) - dtl(iv,i) * a5 * Vh(iv,i,4)

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
        
        mat(iv,i,1,1,2) = mat(iv,i,1,1,2) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * fb2
        mat(iv,i,2,2,2) = mat(iv,i,2,2,2) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * fb2
        mat(iv,i,3,3,2) = mat(iv,i,3,3,2) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * fb2
        mat(iv,i,4,4,2) = mat(iv,i,4,4,2) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * fb2
        mat(iv,i,5,5,2) = mat(iv,i,5,5,2) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * fb2

        mat(iv,i,1,1,3) = mat(iv,i,1,1,3) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * (fb3+fb1)
        mat(iv,i,2,2,3) = mat(iv,i,2,2,3) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * (fb3+fb1)
        mat(iv,i,3,3,3) = mat(iv,i,3,3,3) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * (fb3-fb1)
        mat(iv,i,4,4,3) = mat(iv,i,4,4,3) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * (fb3+fb1)
        mat(iv,i,5,5,3) = mat(iv,i,5,5,3) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * (fb3+fb1)

        mat(iv,i,1,1,4) = mat(iv,i,1,1,4) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * fb4
        mat(iv,i,2,2,4) = mat(iv,i,2,2,4) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * fb4
        mat(iv,i,3,3,4) = mat(iv,i,3,3,4) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * fb4
        mat(iv,i,4,4,4) = mat(iv,i,4,4,4) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * fb4
        mat(iv,i,5,5,4) = mat(iv,i,5,5,4) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * fb4

        mat(iv,i,1,1,5) = mat(iv,i,1,1,5) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * fb5
        mat(iv,i,2,2,5) = mat(iv,i,2,2,5) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * fb5
        mat(iv,i,3,3,5) = mat(iv,i,3,3,5) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * fb5
        mat(iv,i,4,4,5) = mat(iv,i,4,4,5) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * fb5
        mat(iv,i,5,5,5) = mat(iv,i,5,5,5) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * fb5

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

        else if (side .eq. 2) then

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
            do iv = 1, ny
              mat(iv,i,idof,jdof,1) = dtl(iv,i) * a1 * Ah(iv,i,idof,jdof)
              mat(iv,i,idof,jdof,2) = dtl(iv,i) * a2 * Ah(iv,i,idof,jdof)
              mat(iv,i,idof,jdof,3) = dtl(iv,i) * isign*a5 * Ah(iv,i,idof,jdof) + &
                                      rcalcd * alfa * dtl(iv,i) * ( Dh(iv,i,idof,jdof) + &
                                                              im * Dhi(iv,i,idof,jdof) )
              mat(iv,i,idof,jdof,4) = dtl(iv,i) * a4 * Ah(iv,i,idof,jdof)
              mat(iv,i,idof,jdof,5) = zero
            end do
          end do
        end do
        
!.... \hat{A}_i term

        c1 = im * alfa * ga1 * dxiinv
        c2 = im * alfa * ga2 * dxiinv
        c3 = im * alfa * zero * dxiinv
        c4 = im * alfa * ga3 * dxiinv
        c5 = im * alfa * ga4 * dxiinv

        do iv = 1, ny

        mat(iv,i,2,4,1) = mat(iv,i,2,4,1) - c1 * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,3,4,1) = mat(iv,i,3,4,1) - c1 * dtl(iv,i) * Ahi(iv,i,2)
        mat(iv,i,4,2,1) = mat(iv,i,4,2,1) - c1 * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,4,3,1) = mat(iv,i,4,3,1) - c1 * dtl(iv,i) * Ahi(iv,i,2)

        mat(iv,i,2,4,2) = mat(iv,i,2,4,2) - c2 * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,3,4,2) = mat(iv,i,3,4,2) - c2 * dtl(iv,i) * Ahi(iv,i,2)
        mat(iv,i,4,2,2) = mat(iv,i,4,2,2) - c2 * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,4,3,2) = mat(iv,i,4,3,2) - c2 * dtl(iv,i) * Ahi(iv,i,2)

        mat(iv,i,2,4,3) = mat(iv,i,2,4,3) - (c3+c5) * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,3,4,3) = mat(iv,i,3,4,3) - (c3+c5) * dtl(iv,i) * Ahi(iv,i,2)
        mat(iv,i,4,2,3) = mat(iv,i,4,2,3) - (c3+c5) * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,4,3,3) = mat(iv,i,4,3,3) - (c3-c5) * dtl(iv,i) * Ahi(iv,i,2)

        mat(iv,i,2,4,4) = mat(iv,i,2,4,4) - c4 * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,3,4,4) = mat(iv,i,3,4,4) - c4 * dtl(iv,i) * Ahi(iv,i,2)
        mat(iv,i,4,2,4) = mat(iv,i,4,2,4) - c4 * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,4,3,4) = mat(iv,i,4,3,4) - c4 * dtl(iv,i) * Ahi(iv,i,2)

        end do

!.... \hat{V}_{\xi\xi} term

        a1 = alfa * da1 * dxisinv
        a2 = alfa * da2 * dxisinv
        a3 = alfa * da3 * dxisinv
        a4 = alfa * da4 * dxisinv
        a5 = alfa * da5 * dxisinv

        do iv = 1, ny
        
        mat(iv,i,2,2,1) = mat(iv,i,2,2,1) - dtl(iv,i) * a1 * Vh(iv,i,1)
        mat(iv,i,2,3,1) = mat(iv,i,2,3,1) - dtl(iv,i) * a1 * Vh(iv,i,5)
        mat(iv,i,3,2,1) = mat(iv,i,3,2,1) - dtl(iv,i) * a1 * Vh(iv,i,6)
        mat(iv,i,3,3,1) = mat(iv,i,3,3,1) - dtl(iv,i) * a1 * Vh(iv,i,2)
        mat(iv,i,4,4,1) = mat(iv,i,4,4,1) - dtl(iv,i) * a1 * Vh(iv,i,3)
        mat(iv,i,5,5,1) = mat(iv,i,5,5,1) - dtl(iv,i) * a1 * Vh(iv,i,4)

        mat(iv,i,2,2,2) = mat(iv,i,2,2,2) - dtl(iv,i) * a2 * Vh(iv,i,1)
        mat(iv,i,2,3,2) = mat(iv,i,2,3,2) - dtl(iv,i) * a2 * Vh(iv,i,5)
        mat(iv,i,3,2,2) = mat(iv,i,3,2,2) - dtl(iv,i) * a2 * Vh(iv,i,6)
        mat(iv,i,3,3,2) = mat(iv,i,3,3,2) - dtl(iv,i) * a2 * Vh(iv,i,2)
        mat(iv,i,4,4,2) = mat(iv,i,4,4,2) - dtl(iv,i) * a2 * Vh(iv,i,3)
        mat(iv,i,5,5,2) = mat(iv,i,5,5,2) - dtl(iv,i) * a2 * Vh(iv,i,4)

        mat(iv,i,2,2,3) = mat(iv,i,2,2,3) - dtl(iv,i) * (a3+a5) * Vh(iv,i,1)
        mat(iv,i,2,3,3) = mat(iv,i,2,3,3) - dtl(iv,i) * (a3-a5) * Vh(iv,i,5)
        mat(iv,i,3,2,3) = mat(iv,i,3,2,3) - dtl(iv,i) * (a3+a5) * Vh(iv,i,6)
        mat(iv,i,3,3,3) = mat(iv,i,3,3,3) - dtl(iv,i) * (a3-a5) * Vh(iv,i,2)
        mat(iv,i,4,4,3) = mat(iv,i,4,4,3) - dtl(iv,i) * (a3+a5) * Vh(iv,i,3)
        mat(iv,i,5,5,3) = mat(iv,i,5,5,3) - dtl(iv,i) * (a3+a5) * Vh(iv,i,4)

        mat(iv,i,2,2,4) = mat(iv,i,2,2,4) - dtl(iv,i) * a4 * Vh(iv,i,1)
        mat(iv,i,2,3,4) = mat(iv,i,2,3,4) - dtl(iv,i) * a4 * Vh(iv,i,5)
        mat(iv,i,3,2,4) = mat(iv,i,3,2,4) - dtl(iv,i) * a4 * Vh(iv,i,6)
        mat(iv,i,3,3,4) = mat(iv,i,3,3,4) - dtl(iv,i) * a4 * Vh(iv,i,2)
        mat(iv,i,4,4,4) = mat(iv,i,4,4,4) - dtl(iv,i) * a4 * Vh(iv,i,3)
        mat(iv,i,5,5,4) = mat(iv,i,5,5,4) - dtl(iv,i) * a4 * Vh(iv,i,4)

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
        
        mat(iv,i,1,1,1) = mat(iv,i,1,1,1) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * fb1
        mat(iv,i,2,2,1) = mat(iv,i,2,2,1) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * fb1
        mat(iv,i,3,3,1) = mat(iv,i,3,3,1) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * fb1
        mat(iv,i,4,4,1) = mat(iv,i,4,4,1) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * fb1
        mat(iv,i,5,5,1) = mat(iv,i,5,5,1) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * fb1

        mat(iv,i,1,1,2) = mat(iv,i,1,1,2) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * fb2
        mat(iv,i,2,2,2) = mat(iv,i,2,2,2) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * fb2
        mat(iv,i,3,3,2) = mat(iv,i,3,3,2) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * fb2
        mat(iv,i,4,4,2) = mat(iv,i,4,4,2) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * fb2
        mat(iv,i,5,5,2) = mat(iv,i,5,5,2) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * fb2

        mat(iv,i,1,1,3) = mat(iv,i,1,1,3) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * (fb3+fb5)
        mat(iv,i,2,2,3) = mat(iv,i,2,2,3) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * (fb3+fb5)
        mat(iv,i,3,3,3) = mat(iv,i,3,3,3) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * (fb3-fb5)
        mat(iv,i,4,4,3) = mat(iv,i,4,4,3) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * (fb3+fb5)
        mat(iv,i,5,5,3) = mat(iv,i,5,5,3) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * (fb3+fb5)

        mat(iv,i,1,1,4) = mat(iv,i,1,1,4) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * fb4
        mat(iv,i,2,2,4) = mat(iv,i,2,2,4) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * fb4
        mat(iv,i,3,3,4) = mat(iv,i,3,3,4) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * fb4
        mat(iv,i,4,4,4) = mat(iv,i,4,4,4) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * fb4
        mat(iv,i,5,5,4) = mat(iv,i,5,5,4) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * fb4

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
            do iv = 1, ny
              mat(iv,i,idof,jdof,1) = dtl(iv,i) * (a1+isign*a5) * Ah(iv,i,idof,jdof)
              mat(iv,i,idof,jdof,2) = dtl(iv,i) * (a2+isign*a4) * Ah(iv,i,idof,jdof)
              mat(iv,i,idof,jdof,3) = rcalcd * alfa * dtl(iv,i) * ( Dh(iv,i,idof,jdof) + &
                                                              im * Dhi(iv,i,idof,jdof) )
              mat(iv,i,idof,jdof,4) = zero
              mat(iv,i,idof,jdof,5) = zero
            end do
          end do
        end do
        
!.... \hat{A}_i term

        c1 = im * alfa * ga1 * dxiinv
        c2 = im * alfa * ga2 * dxiinv
        c3 = zero
        c4 = im * alfa * ga3 * dxiinv
        c5 = im * alfa * ga4 * dxiinv

        do iv = 1, ny

        mat(iv,i,2,4,1) = mat(iv,i,2,4,1) - (c1+c5) * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,3,4,1) = mat(iv,i,3,4,1) - (c1+c5) * dtl(iv,i) * Ahi(iv,i,2)
        mat(iv,i,4,2,1) = mat(iv,i,4,2,1) - (c1+c5) * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,4,3,1) = mat(iv,i,4,3,1) - (c1-c5) * dtl(iv,i) * Ahi(iv,i,2)
          
        mat(iv,i,2,4,2) = mat(iv,i,2,4,2) - (c2+c4) * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,3,4,2) = mat(iv,i,3,4,2) - (c2+c4) * dtl(iv,i) * Ahi(iv,i,2)
        mat(iv,i,4,2,2) = mat(iv,i,4,2,2) - (c2+c4) * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,4,3,2) = mat(iv,i,4,3,2) - (c2-c4) * dtl(iv,i) * Ahi(iv,i,2)

        mat(iv,i,2,4,3) = mat(iv,i,2,4,3) - c3 * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,3,4,3) = mat(iv,i,3,4,3) - c3 * dtl(iv,i) * Ahi(iv,i,2)
        mat(iv,i,4,2,3) = mat(iv,i,4,2,3) - c3 * dtl(iv,i) * Ahi(iv,i,1)
        mat(iv,i,4,3,3) = mat(iv,i,4,3,3) - c3 * dtl(iv,i) * Ahi(iv,i,2)

        end do

!.... \hat{V}_{\xi\xi} term

        a1 = alfa * da1 * dxisinv
        a2 = alfa * da2 * dxisinv
        a3 = alfa * da3 * dxisinv
        a4 = alfa * da4 * dxisinv
        a5 = alfa * da5 * dxisinv

        do iv = 1, ny
        
        mat(iv,i,2,2,1) = mat(iv,i,2,2,1) - dtl(iv,i) * (a1+a5) * Vh(iv,i,1)
        mat(iv,i,2,3,1) = mat(iv,i,2,3,1) - dtl(iv,i) * (a1-a5) * Vh(iv,i,5)
        mat(iv,i,3,2,1) = mat(iv,i,3,2,1) - dtl(iv,i) * (a1+a5) * Vh(iv,i,6)
        mat(iv,i,3,3,1) = mat(iv,i,3,3,1) - dtl(iv,i) * (a1-a5) * Vh(iv,i,2)
        mat(iv,i,4,4,1) = mat(iv,i,4,4,1) - dtl(iv,i) * (a1+a5) * Vh(iv,i,3)
        mat(iv,i,5,5,1) = mat(iv,i,5,5,1) - dtl(iv,i) * (a1+a5) * Vh(iv,i,4)

        mat(iv,i,2,2,2) = mat(iv,i,2,2,2) - dtl(iv,i) * (a2+a4) * Vh(iv,i,1)
        mat(iv,i,2,3,2) = mat(iv,i,2,3,2) - dtl(iv,i) * (a2-a4) * Vh(iv,i,5)
        mat(iv,i,3,2,2) = mat(iv,i,3,2,2) - dtl(iv,i) * (a2+a4) * Vh(iv,i,6)
        mat(iv,i,3,3,2) = mat(iv,i,3,3,2) - dtl(iv,i) * (a2-a4) * Vh(iv,i,2)
        mat(iv,i,4,4,2) = mat(iv,i,4,4,2) - dtl(iv,i) * (a2+a4) * Vh(iv,i,3)
        mat(iv,i,5,5,2) = mat(iv,i,5,5,2) - dtl(iv,i) * (a2+a4) * Vh(iv,i,4)

        mat(iv,i,2,2,3) = mat(iv,i,2,2,3) - dtl(iv,i) * a3 * Vh(iv,i,1)
        mat(iv,i,2,3,3) = mat(iv,i,2,3,3) - dtl(iv,i) * a3 * Vh(iv,i,5)
        mat(iv,i,3,2,3) = mat(iv,i,3,2,3) - dtl(iv,i) * a3 * Vh(iv,i,6)
        mat(iv,i,3,3,3) = mat(iv,i,3,3,3) - dtl(iv,i) * a3 * Vh(iv,i,2)
        mat(iv,i,4,4,3) = mat(iv,i,4,4,3) - dtl(iv,i) * a3 * Vh(iv,i,3)
        mat(iv,i,5,5,3) = mat(iv,i,5,5,3) - dtl(iv,i) * a3 * Vh(iv,i,4)

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
        
        mat(iv,i,1,1,1) = mat(iv,i,1,1,1) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * (fb1+fb5)
        mat(iv,i,2,2,1) = mat(iv,i,2,2,1) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * (fb1+fb5)
        mat(iv,i,3,3,1) = mat(iv,i,3,3,1) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * (fb1-fb5)
        mat(iv,i,4,4,1) = mat(iv,i,4,4,1) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * (fb1+fb5)
        mat(iv,i,5,5,1) = mat(iv,i,5,5,1) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * (fb1+fb5)

        mat(iv,i,1,1,2) = mat(iv,i,1,1,2) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * (fb2+fb4)
        mat(iv,i,2,2,2) = mat(iv,i,2,2,2) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * (fb2+fb4)
        mat(iv,i,3,3,2) = mat(iv,i,3,3,2) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * (fb2-fb4)
        mat(iv,i,4,4,2) = mat(iv,i,4,4,2) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * (fb2+fb4)
        mat(iv,i,5,5,2) = mat(iv,i,5,5,2) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * (fb2+fb4)

        mat(iv,i,1,1,3) = mat(iv,i,1,1,3) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * fb3
        mat(iv,i,2,2,3) = mat(iv,i,2,2,3) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * fb3
        mat(iv,i,3,3,3) = mat(iv,i,3,3,3) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * fb3
        mat(iv,i,4,4,3) = mat(iv,i,4,4,3) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * fb3
        mat(iv,i,5,5,3) = mat(iv,i,5,5,3) + alfa * dtl(iv,i) * eps_e * buff(iv,i) * fb3

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

        end if         ! side  

        return
        end
