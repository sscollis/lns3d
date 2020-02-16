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
        use buffer
        implicit none
        
        real :: mat(3,ndof,ndof,nx,ny), Ah(ndof,ndof,nx,ny)
        real :: Dh(ndof,ndof,nx,ny), bc(14,ndof,ndof,ny)
        real :: Vh(6,nx,ny)
        real :: dtl(nx,ny)

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
                
        !$omp parallel do private(idof,jdof)
        do iv = 1, ny

!=======================================================================================================!
!.... use higher-order tangent on left node
!=======================================================================================================!

          do idof = 1, ndof
            do jdof = 1, ndof
              mat(1,idof,jdof,1,iv) = zero
              mat(2,idof,jdof,1,iv) = alfa * dtl(1,iv) * gc1 * dxiinv * Ah(idof,jdof,1,iv) + &
                                      alfa * dtl(1,iv) * Dh(idof,jdof,1,iv)
              mat(3,idof,jdof,1,iv) = alfa * dtl(1,iv) * gc2 * dxiinv * Ah(idof,jdof,1,iv)
              bc(1,idof,jdof,iv)    = alfa * dtl(1,iv) * gc3 * dxiinv * Ah(idof,jdof,1,iv)
              bc(2,idof,jdof,iv)    = alfa * dtl(1,iv) * gc4 * dxiinv * Ah(idof,jdof,1,iv)
              bc(3,idof,jdof,iv)    = alfa * dtl(1,iv) * gc5 * dxiinv * Ah(idof,jdof,1,iv)
            end do
          end do
        
!.... \hat{V}_{\eta\eta} term

          mat(2,2,2,1,iv) = mat(2,2,2,1,iv) - alfa * dtl(1,iv) * dd1 * dxisinv * Vh(1,1,iv)
          mat(2,2,3,1,iv) = mat(2,2,3,1,iv) - alfa * dtl(1,iv) * dd1 * dxisinv * Vh(5,1,iv)
          mat(2,3,2,1,iv) = mat(2,3,2,1,iv) - alfa * dtl(1,iv) * dd1 * dxisinv * Vh(6,1,iv)
          mat(2,3,3,1,iv) = mat(2,3,3,1,iv) - alfa * dtl(1,iv) * dd1 * dxisinv * Vh(2,1,iv)
          mat(2,4,4,1,iv) = mat(2,4,4,1,iv) - alfa * dtl(1,iv) * dd1 * dxisinv * Vh(3,1,iv)
          mat(2,5,5,1,iv) = mat(2,5,5,1,iv) - alfa * dtl(1,iv) * dd1 * dxisinv * Vh(4,1,iv)
          
          mat(3,2,2,1,iv) = mat(3,2,2,1,iv) - alfa * dtl(1,iv) * dd2 * dxisinv * Vh(1,1,iv)
          mat(3,2,3,1,iv) = mat(3,2,3,1,iv) - alfa * dtl(1,iv) * dd2 * dxisinv * Vh(5,1,iv)
          mat(3,3,2,1,iv) = mat(3,3,2,1,iv) - alfa * dtl(1,iv) * dd2 * dxisinv * Vh(6,1,iv)
          mat(3,3,3,1,iv) = mat(3,3,3,1,iv) - alfa * dtl(1,iv) * dd2 * dxisinv * Vh(2,1,iv)
          mat(3,4,4,1,iv) = mat(3,4,4,1,iv) - alfa * dtl(1,iv) * dd2 * dxisinv * Vh(3,1,iv)
          mat(3,5,5,1,iv) = mat(3,5,5,1,iv) - alfa * dtl(1,iv) * dd2 * dxisinv * Vh(4,1,iv)
          
          bc(1,2,2,iv)    = bc(1,2,2,iv)    - alfa * dtl(1,iv) * dd3 * dxisinv * Vh(1,1,iv)
          bc(1,2,3,iv)    = bc(1,2,3,iv)    - alfa * dtl(1,iv) * dd3 * dxisinv * Vh(5,1,iv)
          bc(1,3,2,iv)    = bc(1,3,2,iv)    - alfa * dtl(1,iv) * dd3 * dxisinv * Vh(6,1,iv)
          bc(1,3,3,iv)    = bc(1,3,3,iv)    - alfa * dtl(1,iv) * dd3 * dxisinv * Vh(2,1,iv)
          bc(1,4,4,iv)    = bc(1,4,4,iv)    - alfa * dtl(1,iv) * dd3 * dxisinv * Vh(3,1,iv)
          bc(1,5,5,iv)    = bc(1,5,5,iv)    - alfa * dtl(1,iv) * dd3 * dxisinv * Vh(4,1,iv)
          
          bc(2,2,2,iv)    = bc(2,2,2,iv)    - alfa * dtl(1,iv) * dd4 * dxisinv * Vh(1,1,iv)
          bc(2,2,3,iv)    = bc(2,2,3,iv)    - alfa * dtl(1,iv) * dd4 * dxisinv * Vh(5,1,iv)
          bc(2,3,2,iv)    = bc(2,3,2,iv)    - alfa * dtl(1,iv) * dd4 * dxisinv * Vh(6,1,iv)
          bc(2,3,3,iv)    = bc(2,3,3,iv)    - alfa * dtl(1,iv) * dd4 * dxisinv * Vh(2,1,iv)
          bc(2,4,4,iv)    = bc(2,4,4,iv)    - alfa * dtl(1,iv) * dd4 * dxisinv * Vh(3,1,iv)
          bc(2,5,5,iv)    = bc(2,5,5,iv)    - alfa * dtl(1,iv) * dd4 * dxisinv * Vh(4,1,iv)
          
          bc(3,2,2,iv)    = bc(3,2,2,iv)    - alfa * dtl(1,iv) * dd5 * dxisinv * Vh(1,1,iv)
          bc(3,2,3,iv)    = bc(3,2,3,iv)    - alfa * dtl(1,iv) * dd5 * dxisinv * Vh(5,1,iv)
          bc(3,3,2,iv)    = bc(3,3,2,iv)    - alfa * dtl(1,iv) * dd5 * dxisinv * Vh(6,1,iv)
          bc(3,3,3,iv)    = bc(3,3,3,iv)    - alfa * dtl(1,iv) * dd5 * dxisinv * Vh(2,1,iv)
          bc(3,4,4,iv)    = bc(3,4,4,iv)    - alfa * dtl(1,iv) * dd5 * dxisinv * Vh(3,1,iv)
          bc(3,5,5,iv)    = bc(3,5,5,iv)    - alfa * dtl(1,iv) * dd5 * dxisinv * Vh(4,1,iv)

!.... I term

          mat(2,1,1,1,iv) = mat(2,1,1,1,iv) + one
          mat(2,2,2,1,iv) = mat(2,2,2,1,iv) + one
          mat(2,3,3,1,iv) = mat(2,3,3,1,iv) + one
          mat(2,4,4,1,iv) = mat(2,4,4,1,iv) + one
          mat(2,5,5,1,iv) = mat(2,5,5,1,iv) + one

!=======================================================================================================!
!.... use higher-order tangent on first node off the left boundary
!=======================================================================================================!
          do idof = 1, ndof
            do jdof = 1, ndof
              mat(1,idof,jdof,2,iv) = alfa * dtl(2,iv) * gb1 * dxiinv * Ah(idof,jdof,2,iv)
              mat(2,idof,jdof,2,iv) = alfa * dtl(2,iv) * gb2 * dxiinv * Ah(idof,jdof,2,iv) + &
                                      alfa * dtl(2,iv) * Dh(idof,jdof,2,iv)
              mat(3,idof,jdof,2,iv) = alfa * dtl(2,iv) * gb3 * dxiinv * Ah(idof,jdof,2,iv)
              bc(4,idof,jdof,iv)    = alfa * dtl(2,iv) * gb4 * dxiinv * Ah(idof,jdof,2,iv)
              bc(5,idof,jdof,iv)    = alfa * dtl(2,iv) * gb5 * dxiinv * Ah(idof,jdof,2,iv)
            end do
          end do
        
!.... \hat{V}_{\eta\eta} term

          mat(1,2,2,2,iv) = mat(1,2,2,2,iv) - alfa * dtl(2,iv) * db1 * dxisinv * Vh(1,2,iv)
          mat(1,2,3,2,iv) = mat(1,2,3,2,iv) - alfa * dtl(2,iv) * db1 * dxisinv * Vh(5,2,iv)
          mat(1,3,2,2,iv) = mat(1,3,2,2,iv) - alfa * dtl(2,iv) * db1 * dxisinv * Vh(6,2,iv)
          mat(1,3,3,2,iv) = mat(1,3,3,2,iv) - alfa * dtl(2,iv) * db1 * dxisinv * Vh(2,2,iv)
          mat(1,4,4,2,iv) = mat(1,4,4,2,iv) - alfa * dtl(2,iv) * db1 * dxisinv * Vh(3,2,iv)
          mat(1,5,5,2,iv) = mat(1,5,5,2,iv) - alfa * dtl(2,iv) * db1 * dxisinv * Vh(4,2,iv)
          
          mat(2,2,2,2,iv) = mat(2,2,2,2,iv) - alfa * dtl(2,iv) * db2 * dxisinv * Vh(1,2,iv)
          mat(2,2,3,2,iv) = mat(2,2,3,2,iv) - alfa * dtl(2,iv) * db2 * dxisinv * Vh(5,2,iv)
          mat(2,3,2,2,iv) = mat(2,3,2,2,iv) - alfa * dtl(2,iv) * db2 * dxisinv * Vh(6,2,iv)
          mat(2,3,3,2,iv) = mat(2,3,3,2,iv) - alfa * dtl(2,iv) * db2 * dxisinv * Vh(2,2,iv)
          mat(2,4,4,2,iv) = mat(2,4,4,2,iv) - alfa * dtl(2,iv) * db2 * dxisinv * Vh(3,2,iv)
          mat(2,5,5,2,iv) = mat(2,5,5,2,iv) - alfa * dtl(2,iv) * db2 * dxisinv * Vh(4,2,iv)
          
          mat(3,2,2,2,iv) = mat(3,2,2,2,iv) - alfa * dtl(2,iv) * db3 * dxisinv * Vh(1,2,iv)
          mat(3,2,3,2,iv) = mat(3,2,3,2,iv) - alfa * dtl(2,iv) * db3 * dxisinv * Vh(5,2,iv)
          mat(3,3,2,2,iv) = mat(3,3,2,2,iv) - alfa * dtl(2,iv) * db3 * dxisinv * Vh(6,2,iv)
          mat(3,3,3,2,iv) = mat(3,3,3,2,iv) - alfa * dtl(2,iv) * db3 * dxisinv * Vh(2,2,iv)
          mat(3,4,4,2,iv) = mat(3,4,4,2,iv) - alfa * dtl(2,iv) * db3 * dxisinv * Vh(3,2,iv)
          mat(3,5,5,2,iv) = mat(3,5,5,2,iv) - alfa * dtl(2,iv) * db3 * dxisinv * Vh(4,2,iv)
          
          bc(4,2,2,iv)    = bc(4,2,2,iv)    - alfa * dtl(2,iv) * db4 * dxisinv * Vh(1,2,iv)
          bc(4,2,3,iv)    = bc(4,2,3,iv)    - alfa * dtl(2,iv) * db4 * dxisinv * Vh(5,2,iv)
          bc(4,3,2,iv)    = bc(4,3,2,iv)    - alfa * dtl(2,iv) * db4 * dxisinv * Vh(6,2,iv)
          bc(4,3,3,iv)    = bc(4,3,3,iv)    - alfa * dtl(2,iv) * db4 * dxisinv * Vh(2,2,iv)
          bc(4,4,4,iv)    = bc(4,4,4,iv)    - alfa * dtl(2,iv) * db4 * dxisinv * Vh(3,2,iv)
          bc(4,5,5,iv)    = bc(4,5,5,iv)    - alfa * dtl(2,iv) * db4 * dxisinv * Vh(4,2,iv)
          
          bc(5,2,2,iv)    = bc(5,2,2,iv)    - alfa * dtl(2,iv) * db5 * dxisinv * Vh(1,2,iv)
          bc(5,2,3,iv)    = bc(5,2,3,iv)    - alfa * dtl(2,iv) * db5 * dxisinv * Vh(5,2,iv)
          bc(5,3,2,iv)    = bc(5,3,2,iv)    - alfa * dtl(2,iv) * db5 * dxisinv * Vh(6,2,iv)
          bc(5,3,3,iv)    = bc(5,3,3,iv)    - alfa * dtl(2,iv) * db5 * dxisinv * Vh(2,2,iv)
          bc(5,4,4,iv)    = bc(5,4,4,iv)    - alfa * dtl(2,iv) * db5 * dxisinv * Vh(3,2,iv)
          bc(5,5,5,iv)    = bc(5,5,5,iv)    - alfa * dtl(2,iv) * db5 * dxisinv * Vh(4,2,iv)

!.... I term

          mat(2,1,1,2,iv) = mat(2,1,1,2,iv) + one
          mat(2,2,2,2,iv) = mat(2,2,2,2,iv) + one
          mat(2,3,3,2,iv) = mat(2,3,3,2,iv) + one
          mat(2,4,4,2,iv) = mat(2,4,4,2,iv) + one
          mat(2,5,5,2,iv) = mat(2,5,5,2,iv) + one

!.... implicit damping term

          if (eps_e .ne. zero) then

            mat(1,1,1,2,iv) = mat(1,1,1,2,iv) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * (-one)
            mat(1,2,2,2,iv) = mat(1,2,2,2,iv) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * (-one)
            mat(1,3,3,2,iv) = mat(1,3,3,2,iv) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * (-one)
            mat(1,4,4,2,iv) = mat(1,4,4,2,iv) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * (-one)
            mat(1,5,5,2,iv) = mat(1,5,5,2,iv) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * (-one)
            
            mat(2,1,1,2,iv) = mat(2,1,1,2,iv) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * two
            mat(2,2,2,2,iv) = mat(2,2,2,2,iv) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * two
            mat(2,3,3,2,iv) = mat(2,3,3,2,iv) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * two
            mat(2,4,4,2,iv) = mat(2,4,4,2,iv) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * two
            mat(2,5,5,2,iv) = mat(2,5,5,2,iv) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * two
            
            mat(3,1,1,2,iv) = mat(3,1,1,2,iv) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * (-one)
            mat(3,2,2,2,iv) = mat(3,2,2,2,iv) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * (-one)
            mat(3,3,3,2,iv) = mat(3,3,3,2,iv) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * (-one)
            mat(3,4,4,2,iv) = mat(3,4,4,2,iv) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * (-one)
            mat(3,5,5,2,iv) = mat(3,5,5,2,iv) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * (-one)
                
          end if

        end do

        end if          ! lsym

!=======================================================================================================!
!.... use higher-order tangent on second node off the left boundary
!=======================================================================================================!
        !$omp parallel do private(idof,jdof)
        do iv = 1, ny

          do idof = 1, ndof
            do jdof = 1, ndof
              bc(7,idof,jdof,iv)    = alfa * dtl(3,iv) * ga1 * dxiinv * Ah(idof,jdof,3,iv)
              mat(1,idof,jdof,3,iv) = alfa * dtl(3,iv) * ga2 * dxiinv * Ah(idof,jdof,3,iv)
              mat(2,idof,jdof,3,iv) = alfa * dtl(3,iv) * Dh(idof,jdof,3,iv)
              mat(3,idof,jdof,3,iv) = alfa * dtl(3,iv) * ga3 * dxiinv * Ah(idof,jdof,3,iv)
              bc(6,idof,jdof,iv)    = alfa * dtl(3,iv) * ga4 * dxiinv * Ah(idof,jdof,3,iv)
            end do
          end do
        
!.... \hat{V}_{\eta\eta} term

          bc(7,2,2,iv)    = bc(7,2,2,iv)    - alfa * dtl(3,iv) * da1 * dxisinv * Vh(1,3,iv)
          bc(7,2,3,iv)    = bc(7,2,3,iv)    - alfa * dtl(3,iv) * da1 * dxisinv * Vh(5,3,iv)
          bc(7,3,2,iv)    = bc(7,3,2,iv)    - alfa * dtl(3,iv) * da1 * dxisinv * Vh(6,3,iv)
          bc(7,3,3,iv)    = bc(7,3,3,iv)    - alfa * dtl(3,iv) * da1 * dxisinv * Vh(2,3,iv)
          bc(7,4,4,iv)    = bc(7,4,4,iv)    - alfa * dtl(3,iv) * da1 * dxisinv * Vh(3,3,iv)
          bc(7,5,5,iv)    = bc(7,5,5,iv)    - alfa * dtl(3,iv) * da1 * dxisinv * Vh(4,3,iv)
          
          mat(1,2,2,3,iv) = mat(1,2,2,3,iv) - alfa * dtl(3,iv) * da2 * dxisinv * Vh(1,3,iv)
          mat(1,2,3,3,iv) = mat(1,2,3,3,iv) - alfa * dtl(3,iv) * da2 * dxisinv * Vh(5,3,iv)
          mat(1,3,2,3,iv) = mat(1,3,2,3,iv) - alfa * dtl(3,iv) * da2 * dxisinv * Vh(6,3,iv)
          mat(1,3,3,3,iv) = mat(1,3,3,3,iv) - alfa * dtl(3,iv) * da2 * dxisinv * Vh(2,3,iv)
          mat(1,4,4,3,iv) = mat(1,4,4,3,iv) - alfa * dtl(3,iv) * da2 * dxisinv * Vh(3,3,iv)
          mat(1,5,5,3,iv) = mat(1,5,5,3,iv) - alfa * dtl(3,iv) * da2 * dxisinv * Vh(4,3,iv)
          
          mat(2,2,2,3,iv) = mat(2,2,2,3,iv) - alfa * dtl(3,iv) * da3 * dxisinv * Vh(1,3,iv)
          mat(2,2,3,3,iv) = mat(2,2,3,3,iv) - alfa * dtl(3,iv) * da3 * dxisinv * Vh(5,3,iv)
          mat(2,3,2,3,iv) = mat(2,3,2,3,iv) - alfa * dtl(3,iv) * da3 * dxisinv * Vh(6,3,iv)
          mat(2,3,3,3,iv) = mat(2,3,3,3,iv) - alfa * dtl(3,iv) * da3 * dxisinv * Vh(2,3,iv)
          mat(2,4,4,3,iv) = mat(2,4,4,3,iv) - alfa * dtl(3,iv) * da3 * dxisinv * Vh(3,3,iv)
          mat(2,5,5,3,iv) = mat(2,5,5,3,iv) - alfa * dtl(3,iv) * da3 * dxisinv * Vh(4,3,iv)
          
          mat(3,2,2,3,iv) = mat(3,2,2,3,iv) - alfa * dtl(3,iv) * da4 * dxisinv * Vh(1,3,iv)
          mat(3,2,3,3,iv) = mat(3,2,3,3,iv) - alfa * dtl(3,iv) * da4 * dxisinv * Vh(5,3,iv)
          mat(3,3,2,3,iv) = mat(3,3,2,3,iv) - alfa * dtl(3,iv) * da4 * dxisinv * Vh(6,3,iv)
          mat(3,3,3,3,iv) = mat(3,3,3,3,iv) - alfa * dtl(3,iv) * da4 * dxisinv * Vh(2,3,iv)
          mat(3,4,4,3,iv) = mat(3,4,4,3,iv) - alfa * dtl(3,iv) * da4 * dxisinv * Vh(3,3,iv)
          mat(3,5,5,3,iv) = mat(3,5,5,3,iv) - alfa * dtl(3,iv) * da4 * dxisinv * Vh(4,3,iv)
          
          bc(6,2,2,iv)    = bc(6,2,2,iv)    - alfa * dtl(3,iv) * da5 * dxisinv * Vh(1,3,iv)
          bc(6,2,3,iv)    = bc(6,2,3,iv)    - alfa * dtl(3,iv) * da5 * dxisinv * Vh(5,3,iv)
          bc(6,3,2,iv)    = bc(6,3,2,iv)    - alfa * dtl(3,iv) * da5 * dxisinv * Vh(6,3,iv)
          bc(6,3,3,iv)    = bc(6,3,3,iv)    - alfa * dtl(3,iv) * da5 * dxisinv * Vh(2,3,iv)
          bc(6,4,4,iv)    = bc(6,4,4,iv)    - alfa * dtl(3,iv) * da5 * dxisinv * Vh(3,3,iv)
          bc(6,5,5,iv)    = bc(6,5,5,iv)    - alfa * dtl(3,iv) * da5 * dxisinv * Vh(4,3,iv)

!.... I term

          mat(2,1,1,3,iv) = mat(2,1,1,3,iv) + one
          mat(2,2,2,3,iv) = mat(2,2,2,3,iv) + one
          mat(2,3,3,3,iv) = mat(2,3,3,3,iv) + one
          mat(2,4,4,3,iv) = mat(2,4,4,3,iv) + one
          mat(2,5,5,3,iv) = mat(2,5,5,3,iv) + one

!.... implicit damping term updated to fourth order

          if (eps_e .ne. zero) then

            bc(7,1,1,iv) = bc(7,1,1,iv) + alfa * dtl(3,iv) * eps_e * buff(3,iv) * fb1
            bc(7,2,2,iv) = bc(7,2,2,iv) + alfa * dtl(3,iv) * eps_e * buff(3,iv) * fb1
            bc(7,3,3,iv) = bc(7,3,3,iv) + alfa * dtl(3,iv) * eps_e * buff(3,iv) * fb1
            bc(7,4,4,iv) = bc(7,4,4,iv) + alfa * dtl(3,iv) * eps_e * buff(3,iv) * fb1
            bc(7,5,5,iv) = bc(7,5,5,iv) + alfa * dtl(3,iv) * eps_e * buff(3,iv) * fb1
            
            mat(1,1,1,3,iv) = mat(1,1,1,3,iv) + alfa * dtl(3,iv) * eps_e * buff(3,iv) * fb2
            mat(1,2,2,3,iv) = mat(1,2,2,3,iv) + alfa * dtl(3,iv) * eps_e * buff(3,iv) * fb2
            mat(1,3,3,3,iv) = mat(1,3,3,3,iv) + alfa * dtl(3,iv) * eps_e * buff(3,iv) * fb2
            mat(1,4,4,3,iv) = mat(1,4,4,3,iv) + alfa * dtl(3,iv) * eps_e * buff(3,iv) * fb2
            mat(1,5,5,3,iv) = mat(1,5,5,3,iv) + alfa * dtl(3,iv) * eps_e * buff(3,iv) * fb2
            
            mat(2,1,1,3,iv) = mat(2,1,1,3,iv) + alfa * dtl(3,iv) * eps_e * buff(3,iv) * fb3
            mat(2,2,2,3,iv) = mat(2,2,2,3,iv) + alfa * dtl(3,iv) * eps_e * buff(3,iv) * fb3
            mat(2,3,3,3,iv) = mat(2,3,3,3,iv) + alfa * dtl(3,iv) * eps_e * buff(3,iv) * fb3
            mat(2,4,4,3,iv) = mat(2,4,4,3,iv) + alfa * dtl(3,iv) * eps_e * buff(3,iv) * fb3
            mat(2,5,5,3,iv) = mat(2,5,5,3,iv) + alfa * dtl(3,iv) * eps_e * buff(3,iv) * fb3
            
            mat(3,1,1,3,iv) = mat(3,1,1,3,iv) + alfa * dtl(3,iv) * eps_e * buff(3,iv) * fb4
            mat(3,2,2,3,iv) = mat(3,2,2,3,iv) + alfa * dtl(3,iv) * eps_e * buff(3,iv) * fb4
            mat(3,3,3,3,iv) = mat(3,3,3,3,iv) + alfa * dtl(3,iv) * eps_e * buff(3,iv) * fb4
            mat(3,4,4,3,iv) = mat(3,4,4,3,iv) + alfa * dtl(3,iv) * eps_e * buff(3,iv) * fb4
            mat(3,5,5,3,iv) = mat(3,5,5,3,iv) + alfa * dtl(3,iv) * eps_e * buff(3,iv) * fb4
            
            bc(6,1,1,iv) = bc(6,1,1,iv) + alfa * dtl(3,iv) * eps_e * buff(3,iv) * fb5
            bc(6,2,2,iv) = bc(6,2,2,iv) + alfa * dtl(3,iv) * eps_e * buff(3,iv) * fb5
            bc(6,3,3,iv) = bc(6,3,3,iv) + alfa * dtl(3,iv) * eps_e * buff(3,iv) * fb5
            bc(6,4,4,iv) = bc(6,4,4,iv) + alfa * dtl(3,iv) * eps_e * buff(3,iv) * fb5
            bc(6,5,5,iv) = bc(6,5,5,iv) + alfa * dtl(3,iv) * eps_e * buff(3,iv) * fb5

          end if

!=======================================================================================================!
!.... use higher-order tangent on second node off the right boundary
!=======================================================================================================!
          do idof = 1, ndof
            do jdof = 1, ndof
              bc(8,idof,jdof,iv)       = alfa * dtl(nx-2,iv) * ga1 * dxiinv * Ah(idof,jdof,nx-2,iv)
              mat(1,idof,jdof,nx-2,iv) = alfa * dtl(nx-2,iv) * ga2 * dxiinv * Ah(idof,jdof,nx-2,iv)
              mat(2,idof,jdof,nx-2,iv) = alfa * dtl(nx-2,iv) * Dh(idof,jdof,nx-2,iv)
              mat(3,idof,jdof,nx-2,iv) = alfa * dtl(nx-2,iv) * ga3 * dxiinv * Ah(idof,jdof,nx-2,iv)
              bc(9,idof,jdof,iv)       = alfa * dtl(nx-2,iv) * ga4 * dxiinv * Ah(idof,jdof,nx-2,iv)
            end do
          end do
                
!.... \hat{V}_{\eta\eta} term

          bc(8,2,2,iv)       = bc(8,2,2,iv)       - alfa * dtl(nx-2,iv) * da1 * dxisinv * Vh(1,nx-2,iv)
          bc(8,2,3,iv)       = bc(8,2,3,iv)       - alfa * dtl(nx-2,iv) * da1 * dxisinv * Vh(5,nx-2,iv)
          bc(8,3,2,iv)       = bc(8,3,2,iv)       - alfa * dtl(nx-2,iv) * da1 * dxisinv * Vh(6,nx-2,iv)
          bc(8,3,3,iv)       = bc(8,3,3,iv)       - alfa * dtl(nx-2,iv) * da1 * dxisinv * Vh(2,nx-2,iv)
          bc(8,4,4,iv)       = bc(8,4,4,iv)       - alfa * dtl(nx-2,iv) * da1 * dxisinv * Vh(3,nx-2,iv)
          bc(8,5,5,iv)       = bc(8,5,5,iv)       - alfa * dtl(nx-2,iv) * da1 * dxisinv * Vh(4,nx-2,iv)
          
          mat(1,2,2,nx-2,iv) = mat(1,2,2,nx-2,iv) - alfa * dtl(nx-2,iv) * da2 * dxisinv * Vh(1,nx-2,iv)
          mat(1,2,3,nx-2,iv) = mat(1,2,3,nx-2,iv) - alfa * dtl(nx-2,iv) * da2 * dxisinv * Vh(5,nx-2,iv)
          mat(1,3,2,nx-2,iv) = mat(1,3,2,nx-2,iv) - alfa * dtl(nx-2,iv) * da2 * dxisinv * Vh(6,nx-2,iv)
          mat(1,3,3,nx-2,iv) = mat(1,3,3,nx-2,iv) - alfa * dtl(nx-2,iv) * da2 * dxisinv * Vh(2,nx-2,iv)
          mat(1,4,4,nx-2,iv) = mat(1,4,4,nx-2,iv) - alfa * dtl(nx-2,iv) * da2 * dxisinv * Vh(3,nx-2,iv)
          mat(1,5,5,nx-2,iv) = mat(1,5,5,nx-2,iv) - alfa * dtl(nx-2,iv) * da2 * dxisinv * Vh(4,nx-2,iv)
          
          mat(2,2,2,nx-2,iv) = mat(2,2,2,nx-2,iv) - alfa * dtl(nx-2,iv) * da3 * dxisinv * Vh(1,nx-2,iv)
          mat(2,2,3,nx-2,iv) = mat(2,2,3,nx-2,iv) - alfa * dtl(nx-2,iv) * da3 * dxisinv * Vh(5,nx-2,iv)
          mat(2,3,2,nx-2,iv) = mat(2,3,2,nx-2,iv) - alfa * dtl(nx-2,iv) * da3 * dxisinv * Vh(6,nx-2,iv)
          mat(2,3,3,nx-2,iv) = mat(2,3,3,nx-2,iv) - alfa * dtl(nx-2,iv) * da3 * dxisinv * Vh(2,nx-2,iv)
          mat(2,4,4,nx-2,iv) = mat(2,4,4,nx-2,iv) - alfa * dtl(nx-2,iv) * da3 * dxisinv * Vh(3,nx-2,iv)
          mat(2,5,5,nx-2,iv) = mat(2,5,5,nx-2,iv) - alfa * dtl(nx-2,iv) * da3 * dxisinv * Vh(4,nx-2,iv)
          
          mat(3,2,2,nx-2,iv) = mat(3,2,2,nx-2,iv) - alfa * dtl(nx-2,iv) * da4 * dxisinv * Vh(1,nx-2,iv)
          mat(3,2,3,nx-2,iv) = mat(3,2,3,nx-2,iv) - alfa * dtl(nx-2,iv) * da4 * dxisinv * Vh(5,nx-2,iv)
          mat(3,3,2,nx-2,iv) = mat(3,3,2,nx-2,iv) - alfa * dtl(nx-2,iv) * da4 * dxisinv * Vh(6,nx-2,iv)
          mat(3,3,3,nx-2,iv) = mat(3,3,3,nx-2,iv) - alfa * dtl(nx-2,iv) * da4 * dxisinv * Vh(2,nx-2,iv)
          mat(3,4,4,nx-2,iv) = mat(3,4,4,nx-2,iv) - alfa * dtl(nx-2,iv) * da4 * dxisinv * Vh(3,nx-2,iv)
          mat(3,5,5,nx-2,iv) = mat(3,5,5,nx-2,iv) - alfa * dtl(nx-2,iv) * da4 * dxisinv * Vh(4,nx-2,iv)
          
          bc(9,2,2,iv)       = bc(9,2,2,iv)       - alfa * dtl(nx-2,iv) * da5 * dxisinv * Vh(1,nx-2,iv)
          bc(9,2,3,iv)       = bc(9,2,3,iv)       - alfa * dtl(nx-2,iv) * da5 * dxisinv * Vh(5,nx-2,iv)
          bc(9,3,2,iv)       = bc(9,3,2,iv)       - alfa * dtl(nx-2,iv) * da5 * dxisinv * Vh(6,nx-2,iv)
          bc(9,3,3,iv)       = bc(9,3,3,iv)       - alfa * dtl(nx-2,iv) * da5 * dxisinv * Vh(2,nx-2,iv)
          bc(9,4,4,iv)       = bc(9,4,4,iv)       - alfa * dtl(nx-2,iv) * da5 * dxisinv * Vh(3,nx-2,iv)
          bc(9,5,5,iv)       = bc(9,5,5,iv)       - alfa * dtl(nx-2,iv) * da5 * dxisinv * Vh(4,nx-2,iv)

!.... I term

          mat(2,1,1,nx-2,iv) = mat(2,1,1,nx-2,iv) + one
          mat(2,2,2,nx-2,iv) = mat(2,2,2,nx-2,iv) + one
          mat(2,3,3,nx-2,iv) = mat(2,3,3,nx-2,iv) + one
          mat(2,4,4,nx-2,iv) = mat(2,4,4,nx-2,iv) + one
          mat(2,5,5,nx-2,iv) = mat(2,5,5,nx-2,iv) + one

!.... implicit damping term updated to fourth order

          if (eps_e .ne. zero) then

            bc(8,1,1,iv) = bc(8,1,1,iv) + alfa * dtl(nx-2,iv) * eps_e * buff(nx-2,iv) * fb1
            bc(8,2,2,iv) = bc(8,2,2,iv) + alfa * dtl(nx-2,iv) * eps_e * buff(nx-2,iv) * fb1
            bc(8,3,3,iv) = bc(8,3,3,iv) + alfa * dtl(nx-2,iv) * eps_e * buff(nx-2,iv) * fb1
            bc(8,4,4,iv) = bc(8,4,4,iv) + alfa * dtl(nx-2,iv) * eps_e * buff(nx-2,iv) * fb1
            bc(8,5,5,iv) = bc(8,5,5,iv) + alfa * dtl(nx-2,iv) * eps_e * buff(nx-2,iv) * fb1
            
            mat(1,1,1,nx-2,iv) = mat(1,1,1,nx-2,iv) + alfa * dtl(nx-2,iv) * eps_e * buff(nx-2,iv) * fb2
            mat(1,2,2,nx-2,iv) = mat(1,2,2,nx-2,iv) + alfa * dtl(nx-2,iv) * eps_e * buff(nx-2,iv) * fb2
            mat(1,3,3,nx-2,iv) = mat(1,3,3,nx-2,iv) + alfa * dtl(nx-2,iv) * eps_e * buff(nx-2,iv) * fb2
            mat(1,4,4,nx-2,iv) = mat(1,4,4,nx-2,iv) + alfa * dtl(nx-2,iv) * eps_e * buff(nx-2,iv) * fb2
            mat(1,5,5,nx-2,iv) = mat(1,5,5,nx-2,iv) + alfa * dtl(nx-2,iv) * eps_e * buff(nx-2,iv) * fb2
            
            mat(2,1,1,nx-2,iv) = mat(2,1,1,nx-2,iv) + alfa * dtl(nx-2,iv) * eps_e * buff(nx-2,iv) * fb3
            mat(2,2,2,nx-2,iv) = mat(2,2,2,nx-2,iv) + alfa * dtl(nx-2,iv) * eps_e * buff(nx-2,iv) * fb3
            mat(2,3,3,nx-2,iv) = mat(2,3,3,nx-2,iv) + alfa * dtl(nx-2,iv) * eps_e * buff(nx-2,iv) * fb3
            mat(2,4,4,nx-2,iv) = mat(2,4,4,nx-2,iv) + alfa * dtl(nx-2,iv) * eps_e * buff(nx-2,iv) * fb3
            mat(2,5,5,nx-2,iv) = mat(2,5,5,nx-2,iv) + alfa * dtl(nx-2,iv) * eps_e * buff(nx-2,iv) * fb3
            
            mat(3,1,1,nx-2,iv) = mat(3,1,1,nx-2,iv) + alfa * dtl(nx-2,iv) * eps_e * buff(nx-2,iv) * fb4
            mat(3,2,2,nx-2,iv) = mat(3,2,2,nx-2,iv) + alfa * dtl(nx-2,iv) * eps_e * buff(nx-2,iv) * fb4
            mat(3,3,3,nx-2,iv) = mat(3,3,3,nx-2,iv) + alfa * dtl(nx-2,iv) * eps_e * buff(nx-2,iv) * fb4
            mat(3,4,4,nx-2,iv) = mat(3,4,4,nx-2,iv) + alfa * dtl(nx-2,iv) * eps_e * buff(nx-2,iv) * fb4
            mat(3,5,5,nx-2,iv) = mat(3,5,5,nx-2,iv) + alfa * dtl(nx-2,iv) * eps_e * buff(nx-2,iv) * fb4
            
            bc(9,1,1,iv) = bc(9,1,1,iv) + alfa * dtl(nx-2,iv) * eps_e * buff(nx-2,iv) * fb5
            bc(9,2,2,iv) = bc(9,2,2,iv) + alfa * dtl(nx-2,iv) * eps_e * buff(nx-2,iv) * fb5
            bc(9,3,3,iv) = bc(9,3,3,iv) + alfa * dtl(nx-2,iv) * eps_e * buff(nx-2,iv) * fb5
            bc(9,4,4,iv) = bc(9,4,4,iv) + alfa * dtl(nx-2,iv) * eps_e * buff(nx-2,iv) * fb5
            bc(9,5,5,iv) = bc(9,5,5,iv) + alfa * dtl(nx-2,iv) * eps_e * buff(nx-2,iv) * fb5

          end if

        end do

        if (rsym) then
        
          call lhsbs1( mat, Ah, Dh, Vh, bc, dtl, 2 )             ! apply a symmetry condition

        else
                
        !$omp parallel do private(idof,jdof)
        do iv = 1, ny

!=======================================================================================================!
!.... use higher-order tangent on first node off the right boundary
!=======================================================================================================!
          do idof = 1, ndof
            do jdof = 1, ndof
              bc(10,idof,jdof,iv)      = -alfa * dtl(nx-1,iv) * gb5 * dxiinv * Ah(idof,jdof,nx-1,iv)
              bc(11,idof,jdof,iv)      = -alfa * dtl(nx-1,iv) * gb4 * dxiinv * Ah(idof,jdof,nx-1,iv)
              mat(1,idof,jdof,nx-1,iv) = -alfa * dtl(nx-1,iv) * gb3 * dxiinv * Ah(idof,jdof,nx-1,iv)
              mat(2,idof,jdof,nx-1,iv) = -alfa * dtl(nx-1,iv) * gb2 * dxiinv * Ah(idof,jdof,nx-1,iv) + &
                                          alfa * dtl(nx-1,iv) * Dh(idof,jdof,nx-1,iv)
              mat(3,idof,jdof,nx-1,iv) = -alfa * dtl(nx-1,iv) * gb1 * dxiinv * Ah(idof,jdof,nx-1,iv)
            end do
          end do
        
!.... \hat{V}_{\eta\eta} term

          bc(10,2,2,iv)      = bc(10,2,2,iv)      - alfa * dtl(nx-1,iv) * db5 * dxisinv * Vh(1,nx-1,iv)
          bc(10,2,3,iv)      = bc(10,2,3,iv)      - alfa * dtl(nx-1,iv) * db5 * dxisinv * Vh(5,nx-1,iv)
          bc(10,3,2,iv)      = bc(10,3,2,iv)      - alfa * dtl(nx-1,iv) * db5 * dxisinv * Vh(6,nx-1,iv)
          bc(10,3,3,iv)      = bc(10,3,3,iv)      - alfa * dtl(nx-1,iv) * db5 * dxisinv * Vh(2,nx-1,iv)
          bc(10,4,4,iv)      = bc(10,4,4,iv)      - alfa * dtl(nx-1,iv) * db5 * dxisinv * Vh(3,nx-1,iv)
          bc(10,5,5,iv)      = bc(10,5,5,iv)      - alfa * dtl(nx-1,iv) * db5 * dxisinv * Vh(4,nx-1,iv)
          
          bc(11,2,2,iv)      = bc(11,2,2,iv)      - alfa * dtl(nx-1,iv) * db4 * dxisinv * Vh(1,nx-1,iv)
          bc(11,2,3,iv)      = bc(11,2,3,iv)      - alfa * dtl(nx-1,iv) * db4 * dxisinv * Vh(5,nx-1,iv)
          bc(11,3,2,iv)      = bc(11,3,2,iv)      - alfa * dtl(nx-1,iv) * db4 * dxisinv * Vh(6,nx-1,iv)
          bc(11,3,3,iv)      = bc(11,3,3,iv)      - alfa * dtl(nx-1,iv) * db4 * dxisinv * Vh(2,nx-1,iv)
          bc(11,4,4,iv)      = bc(11,4,4,iv)      - alfa * dtl(nx-1,iv) * db4 * dxisinv * Vh(3,nx-1,iv)
          bc(11,5,5,iv)      = bc(11,5,5,iv)      - alfa * dtl(nx-1,iv) * db4 * dxisinv * Vh(4,nx-1,iv)
          
          mat(1,2,2,nx-1,iv) = mat(1,2,2,nx-1,iv) - alfa * dtl(nx-1,iv) * db3 * dxisinv * Vh(1,nx-1,iv)
          mat(1,2,3,nx-1,iv) = mat(1,2,3,nx-1,iv) - alfa * dtl(nx-1,iv) * db3 * dxisinv * Vh(5,nx-1,iv)
          mat(1,3,2,nx-1,iv) = mat(1,3,2,nx-1,iv) - alfa * dtl(nx-1,iv) * db3 * dxisinv * Vh(6,nx-1,iv)
          mat(1,3,3,nx-1,iv) = mat(1,3,3,nx-1,iv) - alfa * dtl(nx-1,iv) * db3 * dxisinv * Vh(2,nx-1,iv)
          mat(1,4,4,nx-1,iv) = mat(1,4,4,nx-1,iv) - alfa * dtl(nx-1,iv) * db3 * dxisinv * Vh(3,nx-1,iv)
          mat(1,5,5,nx-1,iv) = mat(1,5,5,nx-1,iv) - alfa * dtl(nx-1,iv) * db3 * dxisinv * Vh(4,nx-1,iv)
          
          mat(2,2,2,nx-1,iv) = mat(2,2,2,nx-1,iv) - alfa * dtl(nx-1,iv) * db2 * dxisinv * Vh(1,nx-1,iv)
          mat(2,2,3,nx-1,iv) = mat(2,2,3,nx-1,iv) - alfa * dtl(nx-1,iv) * db2 * dxisinv * Vh(5,nx-1,iv)
          mat(2,3,2,nx-1,iv) = mat(2,3,2,nx-1,iv) - alfa * dtl(nx-1,iv) * db2 * dxisinv * Vh(6,nx-1,iv)
          mat(2,3,3,nx-1,iv) = mat(2,3,3,nx-1,iv) - alfa * dtl(nx-1,iv) * db2 * dxisinv * Vh(2,nx-1,iv)
          mat(2,4,4,nx-1,iv) = mat(2,4,4,nx-1,iv) - alfa * dtl(nx-1,iv) * db2 * dxisinv * Vh(3,nx-1,iv)
          mat(2,5,5,nx-1,iv) = mat(2,5,5,nx-1,iv) - alfa * dtl(nx-1,iv) * db2 * dxisinv * Vh(4,nx-1,iv)
          
          mat(3,2,2,nx-1,iv) = mat(3,2,2,nx-1,iv) - alfa * dtl(nx-1,iv) * db1 * dxisinv * Vh(1,nx-1,iv)
          mat(3,2,3,nx-1,iv) = mat(3,2,3,nx-1,iv) - alfa * dtl(nx-1,iv) * db1 * dxisinv * Vh(5,nx-1,iv)
          mat(3,3,2,nx-1,iv) = mat(3,3,2,nx-1,iv) - alfa * dtl(nx-1,iv) * db1 * dxisinv * Vh(6,nx-1,iv)
          mat(3,3,3,nx-1,iv) = mat(3,3,3,nx-1,iv) - alfa * dtl(nx-1,iv) * db1 * dxisinv * Vh(2,nx-1,iv)
          mat(3,4,4,nx-1,iv) = mat(3,4,4,nx-1,iv) - alfa * dtl(nx-1,iv) * db1 * dxisinv * Vh(3,nx-1,iv)
          mat(3,5,5,nx-1,iv) = mat(3,5,5,nx-1,iv) - alfa * dtl(nx-1,iv) * db1 * dxisinv * Vh(4,nx-1,iv)
          
!.... I term

          mat(2,1,1,nx-1,iv) = mat(2,1,1,nx-1,iv) + one
          mat(2,2,2,nx-1,iv) = mat(2,2,2,nx-1,iv) + one
          mat(2,3,3,nx-1,iv) = mat(2,3,3,nx-1,iv) + one
          mat(2,4,4,nx-1,iv) = mat(2,4,4,nx-1,iv) + one
          mat(2,5,5,nx-1,iv) = mat(2,5,5,nx-1,iv) + one

!.... implicit damping term

          if (eps_e .ne. zero) then
            mat(1,1,1,nx-1,iv) = mat(1,1,1,nx-1,iv) + alfa * dtl(nx-1,iv) * eps_e * buff(nx-1,iv) * (-one)
            mat(1,2,2,nx-1,iv) = mat(1,2,2,nx-1,iv) + alfa * dtl(nx-1,iv) * eps_e * buff(nx-1,iv) * (-one)
            mat(1,3,3,nx-1,iv) = mat(1,3,3,nx-1,iv) + alfa * dtl(nx-1,iv) * eps_e * buff(nx-1,iv) * (-one)
            mat(1,4,4,nx-1,iv) = mat(1,4,4,nx-1,iv) + alfa * dtl(nx-1,iv) * eps_e * buff(nx-1,iv) * (-one)
            mat(1,5,5,nx-1,iv) = mat(1,5,5,nx-1,iv) + alfa * dtl(nx-1,iv) * eps_e * buff(nx-1,iv) * (-one)
            
            mat(2,1,1,nx-1,iv) = mat(2,1,1,nx-1,iv) + alfa * dtl(nx-1,iv) * eps_e * buff(nx-1,iv) * two
            mat(2,2,2,nx-1,iv) = mat(2,2,2,nx-1,iv) + alfa * dtl(nx-1,iv) * eps_e * buff(nx-1,iv) * two
            mat(2,3,3,nx-1,iv) = mat(2,3,3,nx-1,iv) + alfa * dtl(nx-1,iv) * eps_e * buff(nx-1,iv) * two
            mat(2,4,4,nx-1,iv) = mat(2,4,4,nx-1,iv) + alfa * dtl(nx-1,iv) * eps_e * buff(nx-1,iv) * two
            mat(2,5,5,nx-1,iv) = mat(2,5,5,nx-1,iv) + alfa * dtl(nx-1,iv) * eps_e * buff(nx-1,iv) * two
            
            mat(3,1,1,nx-1,iv) = mat(3,1,1,nx-1,iv) + alfa * dtl(nx-1,iv) * eps_e * buff(nx-1,iv) * (-one)
            mat(3,2,2,nx-1,iv) = mat(3,2,2,nx-1,iv) + alfa * dtl(nx-1,iv) * eps_e * buff(nx-1,iv) * (-one)
            mat(3,3,3,nx-1,iv) = mat(3,3,3,nx-1,iv) + alfa * dtl(nx-1,iv) * eps_e * buff(nx-1,iv) * (-one)
            mat(3,4,4,nx-1,iv) = mat(3,4,4,nx-1,iv) + alfa * dtl(nx-1,iv) * eps_e * buff(nx-1,iv) * (-one)
            mat(3,5,5,nx-1,iv) = mat(3,5,5,nx-1,iv) + alfa * dtl(nx-1,iv) * eps_e * buff(nx-1,iv) * (-one)
          end if

!=======================================================================================================!
!.... use higher-order tangent on the right boundary
!=======================================================================================================!
          do idof = 1, ndof
            do jdof = 1, ndof
              bc(12,idof,jdof,iv)    = -alfa * dtl(nx,iv) * gc5 * dxiinv * Ah(idof,jdof,nx,iv)
              bc(13,idof,jdof,iv)    = -alfa * dtl(nx,iv) * gc4 * dxiinv * Ah(idof,jdof,nx,iv)
              bc(14,idof,jdof,iv)    = -alfa * dtl(nx,iv) * gc3 * dxiinv * Ah(idof,jdof,nx,iv)
              mat(1,idof,jdof,nx,iv) = -alfa * dtl(nx,iv) * gc2 * dxiinv * Ah(idof,jdof,nx,iv)
              mat(2,idof,jdof,nx,iv) = -alfa * dtl(nx,iv) * gc1 * dxiinv * Ah(idof,jdof,nx,iv) + &
                                        alfa * dtl(nx,iv) * Dh(idof,jdof,nx,iv)
              mat(3,idof,jdof,nx,iv) = zero
            end do
          end do
        
!.... \hat{V}_{\eta\eta} term

          bc(12,2,2,iv) = bc(12,2,2,iv) - alfa * dtl(nx,iv) * dd5 * dxisinv * Vh(1,nx,iv)
          bc(12,2,3,iv) = bc(12,2,3,iv) - alfa * dtl(nx,iv) * dd5 * dxisinv * Vh(5,nx,iv)
          bc(12,3,2,iv) = bc(12,3,2,iv) - alfa * dtl(nx,iv) * dd5 * dxisinv * Vh(6,nx,iv)
          bc(12,3,3,iv) = bc(12,3,3,iv) - alfa * dtl(nx,iv) * dd5 * dxisinv * Vh(2,nx,iv)
          bc(12,4,4,iv) = bc(12,4,4,iv) - alfa * dtl(nx,iv) * dd5 * dxisinv * Vh(3,nx,iv)
          bc(12,5,5,iv) = bc(12,5,5,iv) - alfa * dtl(nx,iv) * dd5 * dxisinv * Vh(4,nx,iv)
          
          bc(13,2,2,iv) = bc(13,2,2,iv) - alfa * dtl(nx,iv) * dd4 * dxisinv * Vh(1,nx,iv)
          bc(13,2,3,iv) = bc(13,2,3,iv) - alfa * dtl(nx,iv) * dd4 * dxisinv * Vh(5,nx,iv)
          bc(13,3,2,iv) = bc(13,3,2,iv) - alfa * dtl(nx,iv) * dd4 * dxisinv * Vh(6,nx,iv)
          bc(13,3,3,iv) = bc(13,3,3,iv) - alfa * dtl(nx,iv) * dd4 * dxisinv * Vh(2,nx,iv)
          bc(13,4,4,iv) = bc(13,4,4,iv) - alfa * dtl(nx,iv) * dd4 * dxisinv * Vh(3,nx,iv)
          bc(13,5,5,iv) = bc(13,5,5,iv) - alfa * dtl(nx,iv) * dd4 * dxisinv * Vh(4,nx,iv)
          
          bc(14,2,2,iv) = bc(14,2,2,iv) - alfa * dtl(nx,iv) * dd3 * dxisinv * Vh(1,nx,iv)
          bc(14,2,3,iv) = bc(14,2,3,iv) - alfa * dtl(nx,iv) * dd3 * dxisinv * Vh(5,nx,iv)
          bc(14,3,2,iv) = bc(14,3,2,iv) - alfa * dtl(nx,iv) * dd3 * dxisinv * Vh(6,nx,iv)
          bc(14,3,3,iv) = bc(14,3,3,iv) - alfa * dtl(nx,iv) * dd3 * dxisinv * Vh(2,nx,iv)
          bc(14,4,4,iv) = bc(14,4,4,iv) - alfa * dtl(nx,iv) * dd3 * dxisinv * Vh(3,nx,iv)
          bc(14,5,5,iv) = bc(14,5,5,iv) - alfa * dtl(nx,iv) * dd3 * dxisinv * Vh(4,nx,iv)
          
          mat(1,2,2,nx,iv) = mat(1,2,2,nx,iv) - alfa * dtl(nx,iv) * dd2 * dxisinv * Vh(1,nx,iv)
          mat(1,2,3,nx,iv) = mat(1,2,3,nx,iv) - alfa * dtl(nx,iv) * dd2 * dxisinv * Vh(5,nx,iv)
          mat(1,3,2,nx,iv) = mat(1,3,2,nx,iv) - alfa * dtl(nx,iv) * dd2 * dxisinv * Vh(6,nx,iv)
          mat(1,3,3,nx,iv) = mat(1,3,3,nx,iv) - alfa * dtl(nx,iv) * dd2 * dxisinv * Vh(2,nx,iv)
          mat(1,4,4,nx,iv) = mat(1,4,4,nx,iv) - alfa * dtl(nx,iv) * dd2 * dxisinv * Vh(3,nx,iv)
          mat(1,5,5,nx,iv) = mat(1,5,5,nx,iv) - alfa * dtl(nx,iv) * dd2 * dxisinv * Vh(4,nx,iv)
          
          mat(2,2,2,nx,iv) = mat(2,2,2,nx,iv) - alfa * dtl(nx,iv) * dd1 * dxisinv * Vh(1,nx,iv)
          mat(2,2,3,nx,iv) = mat(2,2,3,nx,iv) - alfa * dtl(nx,iv) * dd1 * dxisinv * Vh(5,nx,iv)
          mat(2,3,2,nx,iv) = mat(2,3,2,nx,iv) - alfa * dtl(nx,iv) * dd1 * dxisinv * Vh(6,nx,iv)
          mat(2,3,3,nx,iv) = mat(2,3,3,nx,iv) - alfa * dtl(nx,iv) * dd1 * dxisinv * Vh(2,nx,iv)
          mat(2,4,4,nx,iv) = mat(2,4,4,nx,iv) - alfa * dtl(nx,iv) * dd1 * dxisinv * Vh(3,nx,iv)
          mat(2,5,5,nx,iv) = mat(2,5,5,nx,iv) - alfa * dtl(nx,iv) * dd1 * dxisinv * Vh(4,nx,iv)

!.... I term

          mat(2,1,1,nx,iv) = mat(2,1,1,nx,iv) + one
          mat(2,2,2,nx,iv) = mat(2,2,2,nx,iv) + one
          mat(2,3,3,nx,iv) = mat(2,3,3,nx,iv) + one
          mat(2,4,4,nx,iv) = mat(2,4,4,nx,iv) + one
          mat(2,5,5,nx,iv) = mat(2,5,5,nx,iv) + one

        end do

        end if                  ! rsym

        return
        end subroutine lhsbt1
