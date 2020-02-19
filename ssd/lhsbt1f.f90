!=======================================================================================================!
        subroutine lhsbt1f( mat, Ah, Dh, Vh, dtl )
!  
!  Correct the LHS for boundary treatment in the \xi direction.
!
!  This version supports the 4th order LHS. 
!
!  Revised: 10-16-95
!
!=======================================================================================================!
        use stuff
        use diff
        use buff_stuff
        implicit none
        
        real :: mat(ny,nx,ndof,ndof,5), Ah(ny,nx,ndof,ndof)
        real :: Dh(ny,nx,ndof,ndof)
        real :: Vh(ny,nx,6)
        real :: dtl(ny,nx)

        real :: dxiinv, dxisinv

        integer :: iv, idof, jdof
!=======================================================================================================!
        dxiinv  = one / dxi
        dxisinv = one / dxi**2

        if (xper) return

        if (lsym) then
        
          call lhsbs1f( mat, Ah, Dh, Vh, dtl, 1 )               ! applies a symmetry condition

        else
                
!=======================================================================================================!
!.... use higher-order tangent on left node
!=======================================================================================================!
        do idof = 1, ndof
          do jdof = 1, ndof
            do iv = 1, ny
              mat(iv,1,idof,jdof,3) = alfa * dtl(iv,1) * gc1 * dxiinv * Ah(iv,1,idof,jdof) + &
                                      alfa * dtl(iv,1) * Dh(iv,1,idof,jdof)
              mat(iv,1,idof,jdof,4) = alfa * dtl(iv,1) * gc2 * dxiinv * Ah(iv,1,idof,jdof)
              mat(iv,1,idof,jdof,5) = alfa * dtl(iv,1) * gc3 * dxiinv * Ah(iv,1,idof,jdof)
              mat(iv,1,idof,jdof,1) = alfa * dtl(iv,1) * gc4 * dxiinv * Ah(iv,1,idof,jdof)
              mat(iv,1,idof,jdof,2) = alfa * dtl(iv,1) * gc5 * dxiinv * Ah(iv,1,idof,jdof)
            end do
          end do
        end do
        
!.... \hat{V}_{\eta\eta} term

        do iv = 1, ny
        
        mat(iv,1,2,2,3) = mat(iv,1,2,2,3) - alfa * dtl(iv,1) * dd1 * dxisinv * Vh(iv,1,1)
        mat(iv,1,2,3,3) = mat(iv,1,2,3,3) - alfa * dtl(iv,1) * dd1 * dxisinv * Vh(iv,1,5)
        mat(iv,1,3,2,3) = mat(iv,1,3,2,3) - alfa * dtl(iv,1) * dd1 * dxisinv * Vh(iv,1,6)
        mat(iv,1,3,3,3) = mat(iv,1,3,3,3) - alfa * dtl(iv,1) * dd1 * dxisinv * Vh(iv,1,2)
        mat(iv,1,4,4,3) = mat(iv,1,4,4,3) - alfa * dtl(iv,1) * dd1 * dxisinv * Vh(iv,1,3)
        mat(iv,1,5,5,3) = mat(iv,1,5,5,3) - alfa * dtl(iv,1) * dd1 * dxisinv * Vh(iv,1,4)

        mat(iv,1,2,2,4) = mat(iv,1,2,2,4) - alfa * dtl(iv,1) * dd2 * dxisinv * Vh(iv,1,1)
        mat(iv,1,2,3,4) = mat(iv,1,2,3,4) - alfa * dtl(iv,1) * dd2 * dxisinv * Vh(iv,1,5)
        mat(iv,1,3,2,4) = mat(iv,1,3,2,4) - alfa * dtl(iv,1) * dd2 * dxisinv * Vh(iv,1,6)
        mat(iv,1,3,3,4) = mat(iv,1,3,3,4) - alfa * dtl(iv,1) * dd2 * dxisinv * Vh(iv,1,2)
        mat(iv,1,4,4,4) = mat(iv,1,4,4,4) - alfa * dtl(iv,1) * dd2 * dxisinv * Vh(iv,1,3)
        mat(iv,1,5,5,4) = mat(iv,1,5,5,4) - alfa * dtl(iv,1) * dd2 * dxisinv * Vh(iv,1,4)

        mat(iv,1,2,2,5) = mat(iv,1,2,2,5) - alfa * dtl(iv,1) * dd3 * dxisinv * Vh(iv,1,1)
        mat(iv,1,2,3,5) = mat(iv,1,2,3,5) - alfa * dtl(iv,1) * dd3 * dxisinv * Vh(iv,1,5)
        mat(iv,1,3,2,5) = mat(iv,1,3,2,5) - alfa * dtl(iv,1) * dd3 * dxisinv * Vh(iv,1,6)
        mat(iv,1,3,3,5) = mat(iv,1,3,3,5) - alfa * dtl(iv,1) * dd3 * dxisinv * Vh(iv,1,2)
        mat(iv,1,4,4,5) = mat(iv,1,4,4,5) - alfa * dtl(iv,1) * dd3 * dxisinv * Vh(iv,1,3)
        mat(iv,1,5,5,5) = mat(iv,1,5,5,5) - alfa * dtl(iv,1) * dd3 * dxisinv * Vh(iv,1,4)

        mat(iv,1,2,2,1) = mat(iv,1,2,2,1) - alfa * dtl(iv,1) * dd4 * dxisinv * Vh(iv,1,1)
        mat(iv,1,2,3,1) = mat(iv,1,2,3,1) - alfa * dtl(iv,1) * dd4 * dxisinv * Vh(iv,1,5)
        mat(iv,1,3,2,1) = mat(iv,1,3,2,1) - alfa * dtl(iv,1) * dd4 * dxisinv * Vh(iv,1,6)
        mat(iv,1,3,3,1) = mat(iv,1,3,3,1) - alfa * dtl(iv,1) * dd4 * dxisinv * Vh(iv,1,2)
        mat(iv,1,4,4,1) = mat(iv,1,4,4,1) - alfa * dtl(iv,1) * dd4 * dxisinv * Vh(iv,1,3)
        mat(iv,1,5,5,1) = mat(iv,1,5,5,1) - alfa * dtl(iv,1) * dd4 * dxisinv * Vh(iv,1,4)

        mat(iv,1,2,2,2) = mat(iv,1,2,2,2) - alfa * dtl(iv,1) * dd5 * dxisinv * Vh(iv,1,1)
        mat(iv,1,2,3,2) = mat(iv,1,2,3,2) - alfa * dtl(iv,1) * dd5 * dxisinv * Vh(iv,1,5)
        mat(iv,1,3,2,2) = mat(iv,1,3,2,2) - alfa * dtl(iv,1) * dd5 * dxisinv * Vh(iv,1,6)
        mat(iv,1,3,3,2) = mat(iv,1,3,3,2) - alfa * dtl(iv,1) * dd5 * dxisinv * Vh(iv,1,2)
        mat(iv,1,4,4,2) = mat(iv,1,4,4,2) - alfa * dtl(iv,1) * dd5 * dxisinv * Vh(iv,1,3)
        mat(iv,1,5,5,2) = mat(iv,1,5,5,2) - alfa * dtl(iv,1) * dd5 * dxisinv * Vh(iv,1,4)

!.... I term

        mat(iv,1,1,1,3) = mat(iv,1,1,1,3) + one
        mat(iv,1,2,2,3) = mat(iv,1,2,2,3) + one
        mat(iv,1,3,3,3) = mat(iv,1,3,3,3) + one
        mat(iv,1,4,4,3) = mat(iv,1,4,4,3) + one
        mat(iv,1,5,5,3) = mat(iv,1,5,5,3) + one

        end do

!=======================================================================================================!
!.... use higher-order tangent on first node off the left boundary
!=======================================================================================================!
        do idof = 1, ndof
          do jdof = 1, ndof
            do iv = 1, ny
              mat(iv,2,idof,jdof,2) = alfa * dtl(iv,2) * gb1 * dxiinv * Ah(iv,2,idof,jdof)
              mat(iv,2,idof,jdof,3) = alfa * dtl(iv,2) * gb2 * dxiinv * Ah(iv,2,idof,jdof) + &
                                      alfa * dtl(iv,2) * Dh(iv,2,idof,jdof)
              mat(iv,2,idof,jdof,4) = alfa * dtl(iv,2) * gb3 * dxiinv * Ah(iv,2,idof,jdof)
              mat(iv,2,idof,jdof,5) = alfa * dtl(iv,2) * gb4 * dxiinv * Ah(iv,2,idof,jdof)
              mat(iv,2,idof,jdof,1) = alfa * dtl(iv,2) * gb5 * dxiinv * Ah(iv,2,idof,jdof)
            end do
          end do
        end do
        
!.... \hat{V}_{\eta\eta} term

        do iv = 1, ny
        
        mat(iv,2,2,2,2) = mat(iv,2,2,2,2) - alfa * dtl(iv,2) * db1 * dxisinv * Vh(iv,2,1)
        mat(iv,2,2,3,2) = mat(iv,2,2,3,2) - alfa * dtl(iv,2) * db1 * dxisinv * Vh(iv,2,5)
        mat(iv,2,3,2,2) = mat(iv,2,3,2,2) - alfa * dtl(iv,2) * db1 * dxisinv * Vh(iv,2,6)
        mat(iv,2,3,3,2) = mat(iv,2,3,3,2) - alfa * dtl(iv,2) * db1 * dxisinv * Vh(iv,2,2)
        mat(iv,2,4,4,2) = mat(iv,2,4,4,2) - alfa * dtl(iv,2) * db1 * dxisinv * Vh(iv,2,3)
        mat(iv,2,5,5,2) = mat(iv,2,5,5,2) - alfa * dtl(iv,2) * db1 * dxisinv * Vh(iv,2,4)

        mat(iv,2,2,2,3) = mat(iv,2,2,2,3) - alfa * dtl(iv,2) * db2 * dxisinv * Vh(iv,2,1)
        mat(iv,2,2,3,3) = mat(iv,2,2,3,3) - alfa * dtl(iv,2) * db2 * dxisinv * Vh(iv,2,5)
        mat(iv,2,3,2,3) = mat(iv,2,3,2,3) - alfa * dtl(iv,2) * db2 * dxisinv * Vh(iv,2,6)
        mat(iv,2,3,3,3) = mat(iv,2,3,3,3) - alfa * dtl(iv,2) * db2 * dxisinv * Vh(iv,2,2)
        mat(iv,2,4,4,3) = mat(iv,2,4,4,3) - alfa * dtl(iv,2) * db2 * dxisinv * Vh(iv,2,3)
        mat(iv,2,5,5,3) = mat(iv,2,5,5,3) - alfa * dtl(iv,2) * db2 * dxisinv * Vh(iv,2,4)

        mat(iv,2,2,2,4) = mat(iv,2,2,2,4) - alfa * dtl(iv,2) * db3 * dxisinv * Vh(iv,2,1)
        mat(iv,2,2,3,4) = mat(iv,2,2,3,4) - alfa * dtl(iv,2) * db3 * dxisinv * Vh(iv,2,5)
        mat(iv,2,3,2,4) = mat(iv,2,3,2,4) - alfa * dtl(iv,2) * db3 * dxisinv * Vh(iv,2,6)
        mat(iv,2,3,3,4) = mat(iv,2,3,3,4) - alfa * dtl(iv,2) * db3 * dxisinv * Vh(iv,2,2)
        mat(iv,2,4,4,4) = mat(iv,2,4,4,4) - alfa * dtl(iv,2) * db3 * dxisinv * Vh(iv,2,3)
        mat(iv,2,5,5,4) = mat(iv,2,5,5,4) - alfa * dtl(iv,2) * db3 * dxisinv * Vh(iv,2,4)

        mat(iv,2,2,2,5) = mat(iv,2,2,2,5) - alfa * dtl(iv,2) * db4 * dxisinv * Vh(iv,2,1)
        mat(iv,2,2,3,5) = mat(iv,2,2,3,5) - alfa * dtl(iv,2) * db4 * dxisinv * Vh(iv,2,5)
        mat(iv,2,3,2,5) = mat(iv,2,3,2,5) - alfa * dtl(iv,2) * db4 * dxisinv * Vh(iv,2,6)
        mat(iv,2,3,3,5) = mat(iv,2,3,3,5) - alfa * dtl(iv,2) * db4 * dxisinv * Vh(iv,2,2)
        mat(iv,2,4,4,5) = mat(iv,2,4,4,5) - alfa * dtl(iv,2) * db4 * dxisinv * Vh(iv,2,3)
        mat(iv,2,5,5,5) = mat(iv,2,5,5,5) - alfa * dtl(iv,2) * db4 * dxisinv * Vh(iv,2,4)

        mat(iv,2,2,2,1) = mat(iv,2,2,2,1) - alfa * dtl(iv,2) * db5 * dxisinv * Vh(iv,2,1)
        mat(iv,2,2,3,1) = mat(iv,2,2,3,1) - alfa * dtl(iv,2) * db5 * dxisinv * Vh(iv,2,5)
        mat(iv,2,3,2,1) = mat(iv,2,3,2,1) - alfa * dtl(iv,2) * db5 * dxisinv * Vh(iv,2,6)
        mat(iv,2,3,3,1) = mat(iv,2,3,3,1) - alfa * dtl(iv,2) * db5 * dxisinv * Vh(iv,2,2)
        mat(iv,2,4,4,1) = mat(iv,2,4,4,1) - alfa * dtl(iv,2) * db5 * dxisinv * Vh(iv,2,3)
        mat(iv,2,5,5,1) = mat(iv,2,5,5,1) - alfa * dtl(iv,2) * db5 * dxisinv * Vh(iv,2,4)

!.... I term

        mat(iv,2,1,1,3) = mat(iv,2,1,1,3) + one
        mat(iv,2,2,2,3) = mat(iv,2,2,2,3) + one
        mat(iv,2,3,3,3) = mat(iv,2,3,3,3) + one
        mat(iv,2,4,4,3) = mat(iv,2,4,4,3) + one
        mat(iv,2,5,5,3) = mat(iv,2,5,5,3) + one

        end do
 
!.... implicit damping term

        if (eps_e .ne. zero) then

        do iv = 1, ny
        mat(iv,2,1,1,2) = mat(iv,2,1,1,2) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * (-one)
        mat(iv,2,2,2,2) = mat(iv,2,2,2,2) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * (-one)
        mat(iv,2,3,3,2) = mat(iv,2,3,3,2) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * (-one)
        mat(iv,2,4,4,2) = mat(iv,2,4,4,2) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * (-one)
        mat(iv,2,5,5,2) = mat(iv,2,5,5,2) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * (-one)

        mat(iv,2,1,1,3) = mat(iv,2,1,1,3) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * two
        mat(iv,2,2,2,3) = mat(iv,2,2,2,3) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * two
        mat(iv,2,3,3,3) = mat(iv,2,3,3,3) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * two
        mat(iv,2,4,4,3) = mat(iv,2,4,4,3) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * two
        mat(iv,2,5,5,3) = mat(iv,2,5,5,3) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * two

        mat(iv,2,1,1,4) = mat(iv,2,1,1,4) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * (-one)
        mat(iv,2,2,2,4) = mat(iv,2,2,2,4) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * (-one)
        mat(iv,2,3,3,4) = mat(iv,2,3,3,4) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * (-one)
        mat(iv,2,4,4,4) = mat(iv,2,4,4,4) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * (-one)
        mat(iv,2,5,5,4) = mat(iv,2,5,5,4) + alfa * dtl(iv,2) * eps_e * buff(iv,2) * (-one)
        end do
                
        end if

        end if          ! lsym


        if (rsym) then
        
          call lhsbs1f( mat, Ah, Dh, Vh, dtl, 2 )               ! applies a symmetry condition

        else
                
!=======================================================================================================!
!.... use higher-order tangent on first node off the right boundary
!=======================================================================================================!
        do idof = 1, ndof
          do jdof = 1, ndof
            do iv = 1, ny
              mat(iv,nx-1,idof,jdof,5) = -alfa * dtl(iv,nx-1) * gb5 * dxiinv * Ah(iv,nx-1,idof,jdof)
              mat(iv,nx-1,idof,jdof,1) = -alfa * dtl(iv,nx-1) * gb4 * dxiinv * Ah(iv,nx-1,idof,jdof)
              mat(iv,nx-1,idof,jdof,2) = -alfa * dtl(iv,nx-1) * gb3 * dxiinv * Ah(iv,nx-1,idof,jdof)
              mat(iv,nx-1,idof,jdof,3) = -alfa * dtl(iv,nx-1) * gb2 * dxiinv * Ah(iv,nx-1,idof,jdof) + &
                                          alfa * dtl(iv,nx-1) * Dh(iv,nx-1,idof,jdof)
              mat(iv,nx-1,idof,jdof,4) = -alfa * dtl(iv,nx-1) * gb1 * dxiinv * Ah(iv,nx-1,idof,jdof)
            end do
          end do
        end do
        
!.... \hat{V}_{\eta\eta} term

        do iv = 1, ny
        
        mat(iv,nx-1,2,2,5) = mat(iv,nx-1,2,2,5) - alfa * dtl(iv,nx-1) * db5 * dxisinv * Vh(iv,nx-1,1)
        mat(iv,nx-1,2,3,5) = mat(iv,nx-1,2,3,5) - alfa * dtl(iv,nx-1) * db5 * dxisinv * Vh(iv,nx-1,5)
        mat(iv,nx-1,3,2,5) = mat(iv,nx-1,3,2,5) - alfa * dtl(iv,nx-1) * db5 * dxisinv * Vh(iv,nx-1,6)
        mat(iv,nx-1,3,3,5) = mat(iv,nx-1,3,3,5) - alfa * dtl(iv,nx-1) * db5 * dxisinv * Vh(iv,nx-1,2)
        mat(iv,nx-1,4,4,5) = mat(iv,nx-1,4,4,5) - alfa * dtl(iv,nx-1) * db5 * dxisinv * Vh(iv,nx-1,3)
        mat(iv,nx-1,5,5,5) = mat(iv,nx-1,5,5,5) - alfa * dtl(iv,nx-1) * db5 * dxisinv * Vh(iv,nx-1,4)

        mat(iv,nx-1,2,2,1) = mat(iv,nx-1,2,2,1) - alfa * dtl(iv,nx-1) * db4 * dxisinv * Vh(iv,nx-1,1)
        mat(iv,nx-1,2,3,1) = mat(iv,nx-1,2,3,1) - alfa * dtl(iv,nx-1) * db4 * dxisinv * Vh(iv,nx-1,5)
        mat(iv,nx-1,3,2,1) = mat(iv,nx-1,3,2,1) - alfa * dtl(iv,nx-1) * db4 * dxisinv * Vh(iv,nx-1,6)
        mat(iv,nx-1,3,3,1) = mat(iv,nx-1,3,3,1) - alfa * dtl(iv,nx-1) * db4 * dxisinv * Vh(iv,nx-1,2)
        mat(iv,nx-1,4,4,1) = mat(iv,nx-1,4,4,1) - alfa * dtl(iv,nx-1) * db4 * dxisinv * Vh(iv,nx-1,3)
        mat(iv,nx-1,5,5,1) = mat(iv,nx-1,5,5,1) - alfa * dtl(iv,nx-1) * db4 * dxisinv * Vh(iv,nx-1,4)

        mat(iv,nx-1,2,2,2) = mat(iv,nx-1,2,2,2) - alfa * dtl(iv,nx-1) * db3 * dxisinv * Vh(iv,nx-1,1)
        mat(iv,nx-1,2,3,2) = mat(iv,nx-1,2,3,2) - alfa * dtl(iv,nx-1) * db3 * dxisinv * Vh(iv,nx-1,5)
        mat(iv,nx-1,3,2,2) = mat(iv,nx-1,3,2,2) - alfa * dtl(iv,nx-1) * db3 * dxisinv * Vh(iv,nx-1,6)
        mat(iv,nx-1,3,3,2) = mat(iv,nx-1,3,3,2) - alfa * dtl(iv,nx-1) * db3 * dxisinv * Vh(iv,nx-1,2)
        mat(iv,nx-1,4,4,2) = mat(iv,nx-1,4,4,2) - alfa * dtl(iv,nx-1) * db3 * dxisinv * Vh(iv,nx-1,3)
        mat(iv,nx-1,5,5,2) = mat(iv,nx-1,5,5,2) - alfa * dtl(iv,nx-1) * db3 * dxisinv * Vh(iv,nx-1,4)

        mat(iv,nx-1,2,2,3) = mat(iv,nx-1,2,2,3) - alfa * dtl(iv,nx-1) * db2 * dxisinv * Vh(iv,nx-1,1)
        mat(iv,nx-1,2,3,3) = mat(iv,nx-1,2,3,3) - alfa * dtl(iv,nx-1) * db2 * dxisinv * Vh(iv,nx-1,5)
        mat(iv,nx-1,3,2,3) = mat(iv,nx-1,3,2,3) - alfa * dtl(iv,nx-1) * db2 * dxisinv * Vh(iv,nx-1,6)
        mat(iv,nx-1,3,3,3) = mat(iv,nx-1,3,3,3) - alfa * dtl(iv,nx-1) * db2 * dxisinv * Vh(iv,nx-1,2)
        mat(iv,nx-1,4,4,3) = mat(iv,nx-1,4,4,3) - alfa * dtl(iv,nx-1) * db2 * dxisinv * Vh(iv,nx-1,3)
        mat(iv,nx-1,5,5,3) = mat(iv,nx-1,5,5,3) - alfa * dtl(iv,nx-1) * db2 * dxisinv * Vh(iv,nx-1,4)

        mat(iv,nx-1,2,2,4) = mat(iv,nx-1,2,2,4) - alfa * dtl(iv,nx-1) * db1 * dxisinv * Vh(iv,nx-1,1)
        mat(iv,nx-1,2,3,4) = mat(iv,nx-1,2,3,4) - alfa * dtl(iv,nx-1) * db1 * dxisinv * Vh(iv,nx-1,5)
        mat(iv,nx-1,3,2,4) = mat(iv,nx-1,3,2,4) - alfa * dtl(iv,nx-1) * db1 * dxisinv * Vh(iv,nx-1,6)
        mat(iv,nx-1,3,3,4) = mat(iv,nx-1,3,3,4) - alfa * dtl(iv,nx-1) * db1 * dxisinv * Vh(iv,nx-1,2)
        mat(iv,nx-1,4,4,4) = mat(iv,nx-1,4,4,4) - alfa * dtl(iv,nx-1) * db1 * dxisinv * Vh(iv,nx-1,3)
        mat(iv,nx-1,5,5,4) = mat(iv,nx-1,5,5,4) - alfa * dtl(iv,nx-1) * db1 * dxisinv * Vh(iv,nx-1,4)

!.... I term

        mat(iv,nx-1,1,1,3) = mat(iv,nx-1,1,1,3) + one
        mat(iv,nx-1,2,2,3) = mat(iv,nx-1,2,2,3) + one
        mat(iv,nx-1,3,3,3) = mat(iv,nx-1,3,3,3) + one
        mat(iv,nx-1,4,4,3) = mat(iv,nx-1,4,4,3) + one
        mat(iv,nx-1,5,5,3) = mat(iv,nx-1,5,5,3) + one

        end do
        
!.... implicit damping term

        if (eps_e .ne. zero) then

        do iv = 1, ny
        mat(iv,nx-1,1,1,2) = mat(iv,nx-1,1,1,2) + alfa * dtl(iv,nx-1) * eps_e * buff(iv,nx-1) * (-one)
        mat(iv,nx-1,2,2,2) = mat(iv,nx-1,2,2,2) + alfa * dtl(iv,nx-1) * eps_e * buff(iv,nx-1) * (-one)
        mat(iv,nx-1,3,3,2) = mat(iv,nx-1,3,3,2) + alfa * dtl(iv,nx-1) * eps_e * buff(iv,nx-1) * (-one)
        mat(iv,nx-1,4,4,2) = mat(iv,nx-1,4,4,2) + alfa * dtl(iv,nx-1) * eps_e * buff(iv,nx-1) * (-one)
        mat(iv,nx-1,5,5,2) = mat(iv,nx-1,5,5,2) + alfa * dtl(iv,nx-1) * eps_e * buff(iv,nx-1) * (-one)

        mat(iv,nx-1,1,1,3) = mat(iv,nx-1,1,1,3) + alfa * dtl(iv,nx-1) * eps_e * buff(iv,nx-1) * two
        mat(iv,nx-1,2,2,3) = mat(iv,nx-1,2,2,3) + alfa * dtl(iv,nx-1) * eps_e * buff(iv,nx-1) * two
        mat(iv,nx-1,3,3,3) = mat(iv,nx-1,3,3,3) + alfa * dtl(iv,nx-1) * eps_e * buff(iv,nx-1) * two
        mat(iv,nx-1,4,4,3) = mat(iv,nx-1,4,4,3) + alfa * dtl(iv,nx-1) * eps_e * buff(iv,nx-1) * two
        mat(iv,nx-1,5,5,3) = mat(iv,nx-1,5,5,3) + alfa * dtl(iv,nx-1) * eps_e * buff(iv,nx-1) * two

        mat(iv,nx-1,1,1,4) = mat(iv,nx-1,1,1,4) + alfa * dtl(iv,nx-1) * eps_e * buff(iv,nx-1) * (-one)
        mat(iv,nx-1,2,2,4) = mat(iv,nx-1,2,2,4) + alfa * dtl(iv,nx-1) * eps_e * buff(iv,nx-1) * (-one)
        mat(iv,nx-1,3,3,4) = mat(iv,nx-1,3,3,4) + alfa * dtl(iv,nx-1) * eps_e * buff(iv,nx-1) * (-one)
        mat(iv,nx-1,4,4,4) = mat(iv,nx-1,4,4,4) + alfa * dtl(iv,nx-1) * eps_e * buff(iv,nx-1) * (-one)
        mat(iv,nx-1,5,5,4) = mat(iv,nx-1,5,5,4) + alfa * dtl(iv,nx-1) * eps_e * buff(iv,nx-1) * (-one)
        end do
        
        end if

!=======================================================================================================!
!.... use higher-order tangent on the right boundary
!=======================================================================================================!
        do idof = 1, ndof
          do jdof = 1, ndof
            do iv = 1, ny
              mat(iv,nx,idof,jdof,4) = -alfa * dtl(iv,nx) * gc5 * dxiinv * Ah(iv,nx,idof,jdof)
              mat(iv,nx,idof,jdof,5) = -alfa * dtl(iv,nx) * gc4 * dxiinv * Ah(iv,nx,idof,jdof)
              mat(iv,nx,idof,jdof,1) = -alfa * dtl(iv,nx) * gc3 * dxiinv * Ah(iv,nx,idof,jdof)
              mat(iv,nx,idof,jdof,2) = -alfa * dtl(iv,nx) * gc2 * dxiinv * Ah(iv,nx,idof,jdof)
              mat(iv,nx,idof,jdof,3) = -alfa * dtl(iv,nx) * gc1 * dxiinv * Ah(iv,nx,idof,jdof) + &
                                        alfa * dtl(iv,nx) * Dh(iv,nx,idof,jdof)
            end do
          end do
        end do
        
!.... \hat{V}_{\eta\eta} term

        do iv = 1, ny

        mat(iv,nx,2,2,4) = mat(iv,nx,2,2,4) - alfa * dtl(iv,nx) * dd5 * dxisinv * Vh(iv,nx,1)
        mat(iv,nx,2,3,4) = mat(iv,nx,2,3,4) - alfa * dtl(iv,nx) * dd5 * dxisinv * Vh(iv,nx,5)
        mat(iv,nx,3,2,4) = mat(iv,nx,3,2,4) - alfa * dtl(iv,nx) * dd5 * dxisinv * Vh(iv,nx,6)
        mat(iv,nx,3,3,4) = mat(iv,nx,3,3,4) - alfa * dtl(iv,nx) * dd5 * dxisinv * Vh(iv,nx,2)
        mat(iv,nx,4,4,4) = mat(iv,nx,4,4,4) - alfa * dtl(iv,nx) * dd5 * dxisinv * Vh(iv,nx,3)
        mat(iv,nx,5,5,4) = mat(iv,nx,5,5,4) - alfa * dtl(iv,nx) * dd5 * dxisinv * Vh(iv,nx,4)

        mat(iv,nx,2,2,5) = mat(iv,nx,2,2,5) - alfa * dtl(iv,nx) * dd4 * dxisinv * Vh(iv,nx,1)
        mat(iv,nx,2,3,5) = mat(iv,nx,2,3,5) - alfa * dtl(iv,nx) * dd4 * dxisinv * Vh(iv,nx,5)
        mat(iv,nx,3,2,5) = mat(iv,nx,3,2,5) - alfa * dtl(iv,nx) * dd4 * dxisinv * Vh(iv,nx,6)
        mat(iv,nx,3,3,5) = mat(iv,nx,3,3,5) - alfa * dtl(iv,nx) * dd4 * dxisinv * Vh(iv,nx,2)
        mat(iv,nx,4,4,5) = mat(iv,nx,4,4,5) - alfa * dtl(iv,nx) * dd4 * dxisinv * Vh(iv,nx,3)
        mat(iv,nx,5,5,5) = mat(iv,nx,5,5,5) - alfa * dtl(iv,nx) * dd4 * dxisinv * Vh(iv,nx,4)

        mat(iv,nx,2,2,1) = mat(iv,nx,2,2,1) - alfa * dtl(iv,nx) * dd3 * dxisinv * Vh(iv,nx,1)
        mat(iv,nx,2,3,1) = mat(iv,nx,2,3,1) - alfa * dtl(iv,nx) * dd3 * dxisinv * Vh(iv,nx,5)
        mat(iv,nx,3,2,1) = mat(iv,nx,3,2,1) - alfa * dtl(iv,nx) * dd3 * dxisinv * Vh(iv,nx,6)
        mat(iv,nx,3,3,1) = mat(iv,nx,3,3,1) - alfa * dtl(iv,nx) * dd3 * dxisinv * Vh(iv,nx,2)
        mat(iv,nx,4,4,1) = mat(iv,nx,4,4,1) - alfa * dtl(iv,nx) * dd3 * dxisinv * Vh(iv,nx,3)
        mat(iv,nx,5,5,1) = mat(iv,nx,5,5,1) - alfa * dtl(iv,nx) * dd3 * dxisinv * Vh(iv,nx,4)

        mat(iv,nx,2,2,2) = mat(iv,nx,2,2,2) - alfa * dtl(iv,nx) * dd2 * dxisinv * Vh(iv,nx,1)
        mat(iv,nx,2,3,2) = mat(iv,nx,2,3,2) - alfa * dtl(iv,nx) * dd2 * dxisinv * Vh(iv,nx,5)
        mat(iv,nx,3,2,2) = mat(iv,nx,3,2,2) - alfa * dtl(iv,nx) * dd2 * dxisinv * Vh(iv,nx,6)
        mat(iv,nx,3,3,2) = mat(iv,nx,3,3,2) - alfa * dtl(iv,nx) * dd2 * dxisinv * Vh(iv,nx,2)
        mat(iv,nx,4,4,2) = mat(iv,nx,4,4,2) - alfa * dtl(iv,nx) * dd2 * dxisinv * Vh(iv,nx,3)
        mat(iv,nx,5,5,2) = mat(iv,nx,5,5,2) - alfa * dtl(iv,nx) * dd2 * dxisinv * Vh(iv,nx,4)

        mat(iv,nx,2,2,3) = mat(iv,nx,2,2,3) - alfa * dtl(iv,nx) * dd1 * dxisinv * Vh(iv,nx,1)
        mat(iv,nx,2,3,3) = mat(iv,nx,2,3,3) - alfa * dtl(iv,nx) * dd1 * dxisinv * Vh(iv,nx,5)
        mat(iv,nx,3,2,3) = mat(iv,nx,3,2,3) - alfa * dtl(iv,nx) * dd1 * dxisinv * Vh(iv,nx,6)
        mat(iv,nx,3,3,3) = mat(iv,nx,3,3,3) - alfa * dtl(iv,nx) * dd1 * dxisinv * Vh(iv,nx,2)
        mat(iv,nx,4,4,3) = mat(iv,nx,4,4,3) - alfa * dtl(iv,nx) * dd1 * dxisinv * Vh(iv,nx,3)
        mat(iv,nx,5,5,3) = mat(iv,nx,5,5,3) - alfa * dtl(iv,nx) * dd1 * dxisinv * Vh(iv,nx,4)

!.... I term

        mat(iv,nx,1,1,3) = mat(iv,nx,1,1,3) + one
        mat(iv,nx,2,2,3) = mat(iv,nx,2,2,3) + one
        mat(iv,nx,3,3,3) = mat(iv,nx,3,3,3) + one
        mat(iv,nx,4,4,3) = mat(iv,nx,4,4,3) + one
        mat(iv,nx,5,5,3) = mat(iv,nx,5,5,3) + one

        end do

        end if         ! rsym

        return
        end
