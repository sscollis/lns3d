!=======================================================================================================!
        subroutine lhsbs1( mat, Ah, Dh, Vh, bc, dtl, side )
!  
!  Correct the LHS for symmetry conditions in the \xi direction.
!
!  This version supports the 4th order LHS. 
!
!  Revised: 10-16-95
!
!=======================================================================================================!
        use global
        use stencil
        use buffer
        implicit none
        
        real :: mat(3,ndof,ndof,nx,ny), Ah(ndof,ndof,nx,ny)
        real :: Dh(ndof,ndof,nx,ny), bc(14,ndof,ndof,ny)
        real :: Vh(6,nx,ny)
        real :: dtl(nx,ny)
        integer :: side

        real :: dxiinv, dxisinv

        integer :: i, iv, idof, jdof

        real :: isign
!=======================================================================================================!
        dxiinv  = one / dxi
        dxisinv = one / dxi**2

        if (side .eq. 1) then

        !$omp parallel do private(idof,jdof,isign)
        do iv = 1, ny

!=======================================================================================================!
!.... symmetric condition on left node
!=======================================================================================================!
        do idof = 1, ndof
          do jdof = 1, ndof
            if (jdof.eq.3) then
              isign = -one
            else
              isign = one
            end if
            mat(2,idof,jdof,1,iv) = alfa * dtl(1,iv) * Dh(idof,jdof,1,iv)
            mat(3,idof,jdof,1,iv) = alfa * dtl(1,iv) * (ga3+isign*ga2) * dxiinv * Ah(idof,jdof,1,iv)
            bc(1,idof,jdof,iv)    = alfa * dtl(1,iv) * (ga4+isign*ga1) * dxiinv * Ah(idof,jdof,1,iv)
            bc(2,idof,jdof,iv)    = zero
            bc(3,idof,jdof,iv)    = zero
          end do
        end do
        
!.... \hat{V}_{\eta\eta} term

        mat(2,2,2,1,iv) = mat(2,2,2,1,iv) - alfa * dtl(1,iv) * da3 * dxisinv * Vh(1,1,iv)
        mat(2,2,3,1,iv) = mat(2,2,3,1,iv) - alfa * dtl(1,iv) * da3 * dxisinv * Vh(5,1,iv)
        mat(2,3,2,1,iv) = mat(2,3,2,1,iv) - alfa * dtl(1,iv) * da3 * dxisinv * Vh(6,1,iv)
        mat(2,3,3,1,iv) = mat(2,3,3,1,iv) - alfa * dtl(1,iv) * da3 * dxisinv * Vh(2,1,iv)
        mat(2,4,4,1,iv) = mat(2,4,4,1,iv) - alfa * dtl(1,iv) * da3 * dxisinv * Vh(3,1,iv)
        mat(2,5,5,1,iv) = mat(2,5,5,1,iv) - alfa * dtl(1,iv) * da3 * dxisinv * Vh(4,1,iv)

        mat(3,2,2,1,iv) = mat(3,2,2,1,iv) - alfa * dtl(1,iv) * (da4+da2) * dxisinv * Vh(1,1,iv)
        mat(3,2,3,1,iv) = mat(3,2,3,1,iv) - alfa * dtl(1,iv) * (da4-da2) * dxisinv * Vh(5,1,iv)
        mat(3,3,2,1,iv) = mat(3,3,2,1,iv) - alfa * dtl(1,iv) * (da4+da2) * dxisinv * Vh(6,1,iv)
        mat(3,3,3,1,iv) = mat(3,3,3,1,iv) - alfa * dtl(1,iv) * (da4-da2) * dxisinv * Vh(2,1,iv)
        mat(3,4,4,1,iv) = mat(3,4,4,1,iv) - alfa * dtl(1,iv) * (da4+da2) * dxisinv * Vh(3,1,iv)
        mat(3,5,5,1,iv) = mat(3,5,5,1,iv) - alfa * dtl(1,iv) * (da4+da2) * dxisinv * Vh(4,1,iv)

        bc(1,2,2,iv) = bc(1,2,2,iv) - alfa * dtl(1,iv) * (da5+da1) * dxisinv * Vh(1,1,iv)
        bc(1,2,3,iv) = bc(1,2,3,iv) - alfa * dtl(1,iv) * (da5-da1) * dxisinv * Vh(5,1,iv)
        bc(1,3,2,iv) = bc(1,3,2,iv) - alfa * dtl(1,iv) * (da5+da1) * dxisinv * Vh(6,1,iv)
        bc(1,3,3,iv) = bc(1,3,3,iv) - alfa * dtl(1,iv) * (da5-da1) * dxisinv * Vh(2,1,iv)
        bc(1,4,4,iv) = bc(1,4,4,iv) - alfa * dtl(1,iv) * (da5+da1) * dxisinv * Vh(3,1,iv)
        bc(1,5,5,iv) = bc(1,5,5,iv) - alfa * dtl(1,iv) * (da5+da1) * dxisinv * Vh(4,1,iv)

!.... I term

        mat(2,1,1,1,iv) = mat(2,1,1,1,iv) + one
        mat(2,2,2,1,iv) = mat(2,2,2,1,iv) + one
        mat(2,3,3,1,iv) = mat(2,3,3,1,iv) + one
        mat(2,4,4,1,iv) = mat(2,4,4,1,iv) + one
        mat(2,5,5,1,iv) = mat(2,5,5,1,iv) + one

!.... implicit damping term

        if (eps_e .ne. zero) then

        mat(2,1,1,1,iv) = mat(2,1,1,1,iv) + alfa * dtl(1,iv) * eps_e * buff(1,iv) * fb3
        mat(2,2,2,1,iv) = mat(2,2,2,1,iv) + alfa * dtl(1,iv) * eps_e * buff(1,iv) * fb3
        mat(2,3,3,1,iv) = mat(2,3,3,1,iv) + alfa * dtl(1,iv) * eps_e * buff(1,iv) * fb3
        mat(2,4,4,1,iv) = mat(2,4,4,1,iv) + alfa * dtl(1,iv) * eps_e * buff(1,iv) * fb3
        mat(2,5,5,1,iv) = mat(2,5,5,1,iv) + alfa * dtl(1,iv) * eps_e * buff(1,iv) * fb3

        mat(3,1,1,1,iv) = mat(3,1,1,1,iv) + alfa * dtl(1,iv) * eps_e * buff(1,iv) * (fb4+fb2)
        mat(3,2,2,1,iv) = mat(3,2,2,1,iv) + alfa * dtl(1,iv) * eps_e * buff(1,iv) * (fb4+fb2)
        mat(3,3,3,1,iv) = mat(3,3,3,1,iv) + alfa * dtl(1,iv) * eps_e * buff(1,iv) * (fb4-fb2)
        mat(3,4,4,1,iv) = mat(3,4,4,1,iv) + alfa * dtl(1,iv) * eps_e * buff(1,iv) * (fb4+fb2)
        mat(3,5,5,1,iv) = mat(3,5,5,1,iv) + alfa * dtl(1,iv) * eps_e * buff(1,iv) * (fb4+fb2)

        bc(1,1,1,iv) = bc(1,1,1,iv) + alfa * dtl(1,iv) * eps_e * buff(1,iv) * (fb5+fb1)
        bc(1,2,2,iv) = bc(1,2,2,iv) + alfa * dtl(1,iv) * eps_e * buff(1,iv) * (fb5+fb1)
        bc(1,3,3,iv) = bc(1,3,3,iv) + alfa * dtl(1,iv) * eps_e * buff(1,iv) * (fb5-fb1)
        bc(1,4,4,iv) = bc(1,4,4,iv) + alfa * dtl(1,iv) * eps_e * buff(1,iv) * (fb5+fb1)
        bc(1,5,5,iv) = bc(1,5,5,iv) + alfa * dtl(1,iv) * eps_e * buff(1,iv) * (fb5+fb1)

        end if

!=======================================================================================================!
!.... symmetric condition on first node off the left boundary
!=======================================================================================================!
        do idof = 1, ndof
          do jdof = 1, ndof
            if (jdof.eq.3) then
              isign = -one
            else
              isign = one
            end if
            mat(1,idof,jdof,2,iv) = alfa * dtl(2,iv) * ga2 * dxiinv * Ah(idof,jdof,2,iv)
            mat(2,idof,jdof,2,iv) = alfa * dtl(2,iv) * isign*ga1 * dxiinv * Ah(idof,jdof,2,iv) + &
                                    alfa * dtl(2,iv) * Dh(idof,jdof,2,iv)
            mat(3,idof,jdof,2,iv) = alfa * dtl(2,iv) * ga3 * dxiinv * Ah(idof,jdof,2,iv)
            bc(4,idof,jdof,iv)    = alfa * dtl(2,iv) * ga4 * dxiinv * Ah(idof,jdof,2,iv)
            bc(5,idof,jdof,iv)    = zero
          end do
        end do
        
!.... \hat{V}_{\eta\eta} term

        mat(1,2,2,2,iv) = mat(1,2,2,2,iv) - alfa * dtl(2,iv) * da2 * dxisinv * Vh(1,2,iv)
        mat(1,2,3,2,iv) = mat(1,2,3,2,iv) - alfa * dtl(2,iv) * da2 * dxisinv * Vh(5,2,iv)
        mat(1,3,2,2,iv) = mat(1,3,2,2,iv) - alfa * dtl(2,iv) * da2 * dxisinv * Vh(6,2,iv)
        mat(1,3,3,2,iv) = mat(1,3,3,2,iv) - alfa * dtl(2,iv) * da2 * dxisinv * Vh(2,2,iv)
        mat(1,4,4,2,iv) = mat(1,4,4,2,iv) - alfa * dtl(2,iv) * da2 * dxisinv * Vh(3,2,iv)
        mat(1,5,5,2,iv) = mat(1,5,5,2,iv) - alfa * dtl(2,iv) * da2 * dxisinv * Vh(4,2,iv)

        mat(2,2,2,2,iv) = mat(2,2,2,2,iv) - alfa * dtl(2,iv) * (da3+da1) * dxisinv * Vh(1,2,iv)
        mat(2,2,3,2,iv) = mat(2,2,3,2,iv) - alfa * dtl(2,iv) * (da3-da1) * dxisinv * Vh(5,2,iv)
        mat(2,3,2,2,iv) = mat(2,3,2,2,iv) - alfa * dtl(2,iv) * (da3+da1) * dxisinv * Vh(6,2,iv)
        mat(2,3,3,2,iv) = mat(2,3,3,2,iv) - alfa * dtl(2,iv) * (da3-da1) * dxisinv * Vh(2,2,iv)
        mat(2,4,4,2,iv) = mat(2,4,4,2,iv) - alfa * dtl(2,iv) * (da3+da1) * dxisinv * Vh(3,2,iv)
        mat(2,5,5,2,iv) = mat(2,5,5,2,iv) - alfa * dtl(2,iv) * (da3+da1) * dxisinv * Vh(4,2,iv)

        mat(3,2,2,2,iv) = mat(3,2,2,2,iv) - alfa * dtl(2,iv) * da4 * dxisinv * Vh(1,2,iv)
        mat(3,2,3,2,iv) = mat(3,2,3,2,iv) - alfa * dtl(2,iv) * da4 * dxisinv * Vh(5,2,iv)
        mat(3,3,2,2,iv) = mat(3,3,2,2,iv) - alfa * dtl(2,iv) * da4 * dxisinv * Vh(6,2,iv)
        mat(3,3,3,2,iv) = mat(3,3,3,2,iv) - alfa * dtl(2,iv) * da4 * dxisinv * Vh(2,2,iv)
        mat(3,4,4,2,iv) = mat(3,4,4,2,iv) - alfa * dtl(2,iv) * da4 * dxisinv * Vh(3,2,iv)
        mat(3,5,5,2,iv) = mat(3,5,5,2,iv) - alfa * dtl(2,iv) * da4 * dxisinv * Vh(4,2,iv)

        bc(4,2,2,iv) = bc(4,2,2,iv) - alfa * dtl(2,iv) * da5 * dxisinv * Vh(1,2,iv)
        bc(4,2,3,iv) = bc(4,2,3,iv) - alfa * dtl(2,iv) * da5 * dxisinv * Vh(5,2,iv)
        bc(4,3,2,iv) = bc(4,3,2,iv) - alfa * dtl(2,iv) * da5 * dxisinv * Vh(6,2,iv)
        bc(4,3,3,iv) = bc(4,3,3,iv) - alfa * dtl(2,iv) * da5 * dxisinv * Vh(2,2,iv)
        bc(4,4,4,iv) = bc(4,4,4,iv) - alfa * dtl(2,iv) * da5 * dxisinv * Vh(3,2,iv)
        bc(4,5,5,iv) = bc(4,5,5,iv) - alfa * dtl(2,iv) * da5 * dxisinv * Vh(4,2,iv)

!.... I term

        mat(2,1,1,2,iv) = mat(2,1,1,2,iv) + one
        mat(2,2,2,2,iv) = mat(2,2,2,2,iv) + one
        mat(2,3,3,2,iv) = mat(2,3,3,2,iv) + one
        mat(2,4,4,2,iv) = mat(2,4,4,2,iv) + one
        mat(2,5,5,2,iv) = mat(2,5,5,2,iv) + one

!.... implicit damping term

        if (eps_e .ne. zero) then

        mat(1,1,1,2,iv) = mat(1,1,1,2,iv) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * fb2
        mat(1,2,2,2,iv) = mat(1,2,2,2,iv) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * fb2
        mat(1,3,3,2,iv) = mat(1,3,3,2,iv) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * fb2
        mat(1,4,4,2,iv) = mat(1,4,4,2,iv) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * fb2
        mat(1,5,5,2,iv) = mat(1,5,5,2,iv) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * fb2

        mat(2,1,1,2,iv) = mat(2,1,1,2,iv) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * (fb3+fb1)
        mat(2,2,2,2,iv) = mat(2,2,2,2,iv) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * (fb3+fb1)
        mat(2,3,3,2,iv) = mat(2,3,3,2,iv) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * (fb3-fb1)
        mat(2,4,4,2,iv) = mat(2,4,4,2,iv) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * (fb3+fb1)
        mat(2,5,5,2,iv) = mat(2,5,5,2,iv) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * (fb3+fb1)

        mat(3,1,1,2,iv) = mat(3,1,1,2,iv) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * fb4
        mat(3,2,2,2,iv) = mat(3,2,2,2,iv) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * fb4
        mat(3,3,3,2,iv) = mat(3,3,3,2,iv) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * fb4
        mat(3,4,4,2,iv) = mat(3,4,4,2,iv) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * fb4
        mat(3,5,5,2,iv) = mat(3,5,5,2,iv) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * fb4

        bc(4,1,1,iv) = bc(4,1,1,iv) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * fb5
        bc(4,2,2,iv) = bc(4,2,2,iv) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * fb5
        bc(4,3,3,iv) = bc(4,3,3,iv) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * fb5
        bc(4,4,4,iv) = bc(4,4,4,iv) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * fb5
        bc(4,5,5,iv) = bc(4,5,5,iv) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * fb5

        end if

        end do

        else if (side .eq. 2) then

        !$omp parallel do private(idof,jdof,isign,i)
        do iv = 1, ny

!=======================================================================================================!
!.... symmetric condition on first node off the right boundary
!=======================================================================================================!
        i = nx - 1

        do idof = 1, ndof
          do jdof = 1, ndof
            if (jdof.eq.3) then
              isign = -one
            else
              isign = one
            end if
            bc(10,idof,jdof,iv)   = zero
            bc(11,idof,jdof,iv)   = alfa * dtl(i,iv) * ga1 * dxiinv * Ah(idof,jdof,i,iv)
            mat(1,idof,jdof,i,iv) = alfa * dtl(i,iv) * ga2 * dxiinv * Ah(idof,jdof,i,iv)
            mat(2,idof,jdof,i,iv) = alfa * dtl(i,iv) * isign*ga4 * dxiinv * Ah(idof,jdof,i,iv) + &
                                    alfa * dtl(i,iv) * Dh(idof,jdof,i,iv)
            mat(3,idof,jdof,i,iv) = alfa * dtl(i,iv) * ga3 * dxiinv * Ah(idof,jdof,i,iv)
          end do
        end do
        
!.... \hat{V}_{\eta\eta} term

        bc(11,2,2,iv) = bc(11,2,2,iv) - alfa * dtl(i,iv) * da1 * dxisinv * Vh(1,i,iv)
        bc(11,2,3,iv) = bc(11,2,3,iv) - alfa * dtl(i,iv) * da1 * dxisinv * Vh(5,i,iv)
        bc(11,3,2,iv) = bc(11,3,2,iv) - alfa * dtl(i,iv) * da1 * dxisinv * Vh(6,i,iv)
        bc(11,3,3,iv) = bc(11,3,3,iv) - alfa * dtl(i,iv) * da1 * dxisinv * Vh(2,i,iv)
        bc(11,4,4,iv) = bc(11,4,4,iv) - alfa * dtl(i,iv) * da1 * dxisinv * Vh(3,i,iv)
        bc(11,5,5,iv) = bc(11,5,5,iv) - alfa * dtl(i,iv) * da1 * dxisinv * Vh(4,i,iv)

        mat(1,2,2,i,iv) = mat(1,2,2,i,iv) - alfa * dtl(i,iv) * da2 * dxisinv * Vh(1,i,iv)
        mat(1,2,3,i,iv) = mat(1,2,3,i,iv) - alfa * dtl(i,iv) * da2 * dxisinv * Vh(5,i,iv)
        mat(1,3,2,i,iv) = mat(1,3,2,i,iv) - alfa * dtl(i,iv) * da2 * dxisinv * Vh(6,i,iv)
        mat(1,3,3,i,iv) = mat(1,3,3,i,iv) - alfa * dtl(i,iv) * da2 * dxisinv * Vh(2,i,iv)
        mat(1,4,4,i,iv) = mat(1,4,4,i,iv) - alfa * dtl(i,iv) * da2 * dxisinv * Vh(3,i,iv)
        mat(1,5,5,i,iv) = mat(1,5,5,i,iv) - alfa * dtl(i,iv) * da2 * dxisinv * Vh(4,i,iv)

        mat(2,2,2,i,iv) = mat(2,2,2,i,iv) - alfa * dtl(i,iv) * (da3+da5) * dxisinv * Vh(1,i,iv)
        mat(2,2,3,i,iv) = mat(2,2,3,i,iv) - alfa * dtl(i,iv) * (da3-da5) * dxisinv * Vh(5,i,iv)
        mat(2,3,2,i,iv) = mat(2,3,2,i,iv) - alfa * dtl(i,iv) * (da3+da5) * dxisinv * Vh(6,i,iv)
        mat(2,3,3,i,iv) = mat(2,3,3,i,iv) - alfa * dtl(i,iv) * (da3-da5) * dxisinv * Vh(2,i,iv)
        mat(2,4,4,i,iv) = mat(2,4,4,i,iv) - alfa * dtl(i,iv) * (da3+da5) * dxisinv * Vh(3,i,iv)
        mat(2,5,5,i,iv) = mat(2,5,5,i,iv) - alfa * dtl(i,iv) * (da3+da5) * dxisinv * Vh(4,i,iv)

        mat(3,2,2,i,iv) = mat(3,2,2,i,iv) - alfa * dtl(i,iv) * da4 * dxisinv * Vh(1,i,iv)
        mat(3,2,3,i,iv) = mat(3,2,3,i,iv) - alfa * dtl(i,iv) * da4 * dxisinv * Vh(5,i,iv)
        mat(3,3,2,i,iv) = mat(3,3,2,i,iv) - alfa * dtl(i,iv) * da4 * dxisinv * Vh(6,i,iv)
        mat(3,3,3,i,iv) = mat(3,3,3,i,iv) - alfa * dtl(i,iv) * da4 * dxisinv * Vh(2,i,iv)
        mat(3,4,4,i,iv) = mat(3,4,4,i,iv) - alfa * dtl(i,iv) * da4 * dxisinv * Vh(3,i,iv)
        mat(3,5,5,i,iv) = mat(3,5,5,i,iv) - alfa * dtl(i,iv) * da4 * dxisinv * Vh(4,i,iv)

!.... I term

        mat(2,1,1,i,iv) = mat(2,1,1,i,iv) + one
        mat(2,2,2,i,iv) = mat(2,2,2,i,iv) + one
        mat(2,3,3,i,iv) = mat(2,3,3,i,iv) + one
        mat(2,4,4,i,iv) = mat(2,4,4,i,iv) + one
        mat(2,5,5,i,iv) = mat(2,5,5,i,iv) + one

!.... implicit damping term

        if (eps_e .ne. zero) then

        bc(11,1,1,iv) = bc(11,1,1,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb1
        bc(11,2,2,iv) = bc(11,2,2,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb1
        bc(11,3,3,iv) = bc(11,3,3,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb1
        bc(11,4,4,iv) = bc(11,4,4,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb1
        bc(11,5,5,iv) = bc(11,5,5,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb1

        mat(1,1,1,i,iv) = mat(1,1,1,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb2
        mat(1,2,2,i,iv) = mat(1,2,2,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb2
        mat(1,3,3,i,iv) = mat(1,3,3,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb2
        mat(1,4,4,i,iv) = mat(1,4,4,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb2
        mat(1,5,5,i,iv) = mat(1,5,5,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb2

        mat(2,1,1,i,iv) = mat(2,1,1,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (fb3+fb5)
        mat(2,2,2,i,iv) = mat(2,2,2,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (fb3+fb5)
        mat(2,3,3,i,iv) = mat(2,3,3,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (fb3-fb5)
        mat(2,4,4,i,iv) = mat(2,4,4,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (fb3+fb5)
        mat(2,5,5,i,iv) = mat(2,5,5,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (fb3+fb5)

        mat(3,1,1,i,iv) = mat(3,1,1,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb4
        mat(3,2,2,i,iv) = mat(3,2,2,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb4
        mat(3,3,3,i,iv) = mat(3,3,3,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb4
        mat(3,4,4,i,iv) = mat(3,4,4,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb4
        mat(3,5,5,i,iv) = mat(3,5,5,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb4

        end if

!=======================================================================================================!
!.... symmetric condition on left node
!=======================================================================================================!
        i = nx

        do idof = 1, ndof
          do jdof = 1, ndof
            if (jdof.eq.3) then
              isign = -one
            else
              isign = one
            end if
            bc(12,idof,jdof,iv)   = zero
            bc(13,idof,jdof,iv)   = zero
            bc(14,idof,jdof,iv)   = alfa * dtl(i,iv) * (ga1+isign*ga4) * dxiinv * Ah(idof,jdof,i,iv)
            mat(1,idof,jdof,i,iv) = alfa * dtl(i,iv) * (ga2+isign*ga3) * dxiinv * Ah(idof,jdof,i,iv)
            mat(2,idof,jdof,i,iv) = alfa * dtl(i,iv) * Dh(idof,jdof,i,iv)
          end do
        end do
        
!.... \hat{V}_{\eta\eta} term

        bc(14,2,2,iv) = bc(14,2,2,iv) - alfa * dtl(i,iv) * (da1+da5) * dxisinv * Vh(1,i,iv)
        bc(14,2,3,iv) = bc(14,2,3,iv) - alfa * dtl(i,iv) * (da1-da5) * dxisinv * Vh(5,i,iv)
        bc(14,3,2,iv) = bc(14,3,2,iv) - alfa * dtl(i,iv) * (da1+da5) * dxisinv * Vh(6,i,iv)
        bc(14,3,3,iv) = bc(14,3,3,iv) - alfa * dtl(i,iv) * (da1-da5) * dxisinv * Vh(2,i,iv)
        bc(14,4,4,iv) = bc(14,4,4,iv) - alfa * dtl(i,iv) * (da1+da5) * dxisinv * Vh(3,i,iv)
        bc(14,5,5,iv) = bc(14,5,5,iv) - alfa * dtl(i,iv) * (da1+da5) * dxisinv * Vh(4,i,iv)

        mat(1,2,2,i,iv) = mat(1,2,2,i,iv) - alfa * dtl(i,iv) * (da2+da4) * dxisinv * Vh(1,i,iv)
        mat(1,2,3,i,iv) = mat(1,2,3,i,iv) - alfa * dtl(i,iv) * (da2-da4) * dxisinv * Vh(5,i,iv)
        mat(1,3,2,i,iv) = mat(1,3,2,i,iv) - alfa * dtl(i,iv) * (da2+da4) * dxisinv * Vh(6,i,iv)
        mat(1,3,3,i,iv) = mat(1,3,3,i,iv) - alfa * dtl(i,iv) * (da2-da4) * dxisinv * Vh(2,i,iv)
        mat(1,4,4,i,iv) = mat(1,4,4,i,iv) - alfa * dtl(i,iv) * (da2+da4) * dxisinv * Vh(3,i,iv)
        mat(1,5,5,i,iv) = mat(1,5,5,i,iv) - alfa * dtl(i,iv) * (da2+da4) * dxisinv * Vh(4,i,iv)

        mat(2,2,2,i,iv) = mat(2,2,2,i,iv) - alfa * dtl(i,iv) * da3 * dxisinv * Vh(1,i,iv)
        mat(2,2,3,i,iv) = mat(2,2,3,i,iv) - alfa * dtl(i,iv) * da3 * dxisinv * Vh(5,i,iv)
        mat(2,3,2,i,iv) = mat(2,3,2,i,iv) - alfa * dtl(i,iv) * da3 * dxisinv * Vh(6,i,iv)
        mat(2,3,3,i,iv) = mat(2,3,3,i,iv) - alfa * dtl(i,iv) * da3 * dxisinv * Vh(2,i,iv)
        mat(2,4,4,i,iv) = mat(2,4,4,i,iv) - alfa * dtl(i,iv) * da3 * dxisinv * Vh(3,i,iv)
        mat(2,5,5,i,iv) = mat(2,5,5,i,iv) - alfa * dtl(i,iv) * da3 * dxisinv * Vh(4,i,iv)

!.... I term

        mat(2,1,1,i,iv) = mat(2,1,1,i,iv) + one
        mat(2,2,2,i,iv) = mat(2,2,2,i,iv) + one
        mat(2,3,3,i,iv) = mat(2,3,3,i,iv) + one
        mat(2,4,4,i,iv) = mat(2,4,4,i,iv) + one
        mat(2,5,5,i,iv) = mat(2,5,5,i,iv) + one

!.... implicit damping term

        if (eps_e .ne. zero) then

        bc(14,1,1,iv) = bc(14,1,1,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (fb1+fb5)
        bc(14,2,2,iv) = bc(14,2,2,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (fb1+fb5)
        bc(14,3,3,iv) = bc(14,3,3,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (fb1-fb5)
        bc(14,4,4,iv) = bc(14,4,4,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (fb1+fb5)
        bc(14,5,5,iv) = bc(14,5,5,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (fb1+fb5)

        mat(1,1,1,i,iv) = mat(1,1,1,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (fb2+fb4)
        mat(1,2,2,i,iv) = mat(1,2,2,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (fb2+fb4)
        mat(1,3,3,i,iv) = mat(1,3,3,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (fb2-fb4)
        mat(1,4,4,i,iv) = mat(1,4,4,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (fb2+fb4)
        mat(1,5,5,i,iv) = mat(1,5,5,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * (fb2+fb4)

        mat(2,1,1,i,iv) = mat(2,1,1,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb3
        mat(2,2,2,i,iv) = mat(2,2,2,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb3
        mat(2,3,3,i,iv) = mat(2,3,3,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb3
        mat(2,4,4,i,iv) = mat(2,4,4,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb3
        mat(2,5,5,i,iv) = mat(2,5,5,i,iv) + alfa * dtl(i,iv) * eps_e * buff(i,iv) * fb3

        end if

        end do

        end if

        return
        end subroutine lhsbs1
