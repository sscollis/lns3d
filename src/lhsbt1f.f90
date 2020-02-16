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
        use global
        use stencil
        use buffer
        implicit none
        
        real :: mat(5,ndof,ndof,nx,ny), Ah(ndof,ndof,nx,ny)
        real :: Dh(ndof,ndof,nx,ny)
        real :: Vh(6,nx,ny)
        real :: dtl(nx,ny)

        real :: dxiinv, dxisinv

        integer :: iv, idof, jdof
!=======================================================================================================!
        dxiinv  = one / dxi
        dxisinv = one / dxi**2

        if (xper) return

        if (lsym) then
        
          call lhsbs1f( mat, Ah, Dh, Vh, dtl, 1 )               ! applies a symmetry condition

        else
                
        !$omp parallel do private (idof,jdof)
        do iv = 1, ny
!=======================================================================================================!
!.... use higher-order tangent on left node
!=======================================================================================================!
          do idof = 1, ndof
            do jdof = 1, ndof
              mat(3,idof,jdof,1,iv) = alfa * dtl(1,iv) * gc1 * dxiinv * Ah(idof,jdof,1,iv) + &
                                      alfa * dtl(1,iv) * Dh(idof,jdof,1,iv)
              mat(4,idof,jdof,1,iv) = alfa * dtl(1,iv) * gc2 * dxiinv * Ah(idof,jdof,1,iv)
              mat(5,idof,jdof,1,iv) = alfa * dtl(1,iv) * gc3 * dxiinv * Ah(idof,jdof,1,iv)
              mat(1,idof,jdof,1,iv) = alfa * dtl(1,iv) * gc4 * dxiinv * Ah(idof,jdof,1,iv)
              mat(2,idof,jdof,1,iv) = alfa * dtl(1,iv) * gc5 * dxiinv * Ah(idof,jdof,1,iv)
            end do
          end do
        
!.... \hat{V}_{\eta\eta} term

          mat(3,2,2,1,iv) = mat(3,2,2,1,iv) - alfa * dtl(1,iv) * dd1 * dxisinv * Vh(1,1,iv)
          mat(3,2,3,1,iv) = mat(3,2,3,1,iv) - alfa * dtl(1,iv) * dd1 * dxisinv * Vh(5,1,iv)
          mat(3,3,2,1,iv) = mat(3,3,2,1,iv) - alfa * dtl(1,iv) * dd1 * dxisinv * Vh(6,1,iv)
          mat(3,3,3,1,iv) = mat(3,3,3,1,iv) - alfa * dtl(1,iv) * dd1 * dxisinv * Vh(2,1,iv)
          mat(3,4,4,1,iv) = mat(3,4,4,1,iv) - alfa * dtl(1,iv) * dd1 * dxisinv * Vh(3,1,iv)
          mat(3,5,5,1,iv) = mat(3,5,5,1,iv) - alfa * dtl(1,iv) * dd1 * dxisinv * Vh(4,1,iv)
          
          mat(4,2,2,1,iv) = mat(4,2,2,1,iv) - alfa * dtl(1,iv) * dd2 * dxisinv * Vh(1,1,iv)
          mat(4,2,3,1,iv) = mat(4,2,3,1,iv) - alfa * dtl(1,iv) * dd2 * dxisinv * Vh(5,1,iv)
          mat(4,3,2,1,iv) = mat(4,3,2,1,iv) - alfa * dtl(1,iv) * dd2 * dxisinv * Vh(6,1,iv)
          mat(4,3,3,1,iv) = mat(4,3,3,1,iv) - alfa * dtl(1,iv) * dd2 * dxisinv * Vh(2,1,iv)
          mat(4,4,4,1,iv) = mat(4,4,4,1,iv) - alfa * dtl(1,iv) * dd2 * dxisinv * Vh(3,1,iv)
          mat(4,5,5,1,iv) = mat(4,5,5,1,iv) - alfa * dtl(1,iv) * dd2 * dxisinv * Vh(4,1,iv)
          
          mat(5,2,2,1,iv) = mat(5,2,2,1,iv) - alfa * dtl(1,iv) * dd3 * dxisinv * Vh(1,1,iv)
          mat(5,2,3,1,iv) = mat(5,2,3,1,iv) - alfa * dtl(1,iv) * dd3 * dxisinv * Vh(5,1,iv)
          mat(5,3,2,1,iv) = mat(5,3,2,1,iv) - alfa * dtl(1,iv) * dd3 * dxisinv * Vh(6,1,iv)
          mat(5,3,3,1,iv) = mat(5,3,3,1,iv) - alfa * dtl(1,iv) * dd3 * dxisinv * Vh(2,1,iv)
          mat(5,4,4,1,iv) = mat(5,4,4,1,iv) - alfa * dtl(1,iv) * dd3 * dxisinv * Vh(3,1,iv)
          mat(5,5,5,1,iv) = mat(5,5,5,1,iv) - alfa * dtl(1,iv) * dd3 * dxisinv * Vh(4,1,iv)
          
          mat(1,2,2,1,iv) = mat(1,2,2,1,iv) - alfa * dtl(1,iv) * dd4 * dxisinv * Vh(1,1,iv)
          mat(1,2,3,1,iv) = mat(1,2,3,1,iv) - alfa * dtl(1,iv) * dd4 * dxisinv * Vh(5,1,iv)
          mat(1,3,2,1,iv) = mat(1,3,2,1,iv) - alfa * dtl(1,iv) * dd4 * dxisinv * Vh(6,1,iv)
          mat(1,3,3,1,iv) = mat(1,3,3,1,iv) - alfa * dtl(1,iv) * dd4 * dxisinv * Vh(2,1,iv)
          mat(1,4,4,1,iv) = mat(1,4,4,1,iv) - alfa * dtl(1,iv) * dd4 * dxisinv * Vh(3,1,iv)
          mat(1,5,5,1,iv) = mat(1,5,5,1,iv) - alfa * dtl(1,iv) * dd4 * dxisinv * Vh(4,1,iv)
          
          mat(2,2,2,1,iv) = mat(2,2,2,1,iv) - alfa * dtl(1,iv) * dd5 * dxisinv * Vh(1,1,iv)
          mat(2,2,3,1,iv) = mat(2,2,3,1,iv) - alfa * dtl(1,iv) * dd5 * dxisinv * Vh(5,1,iv)
          mat(2,3,2,1,iv) = mat(2,3,2,1,iv) - alfa * dtl(1,iv) * dd5 * dxisinv * Vh(6,1,iv)
          mat(2,3,3,1,iv) = mat(2,3,3,1,iv) - alfa * dtl(1,iv) * dd5 * dxisinv * Vh(2,1,iv)
          mat(2,4,4,1,iv) = mat(2,4,4,1,iv) - alfa * dtl(1,iv) * dd5 * dxisinv * Vh(3,1,iv)
          mat(2,5,5,1,iv) = mat(2,5,5,1,iv) - alfa * dtl(1,iv) * dd5 * dxisinv * Vh(4,1,iv)
          
!.... I term

          mat(3,1,1,1,iv) = mat(3,1,1,1,iv) + one
          mat(3,2,2,1,iv) = mat(3,2,2,1,iv) + one
          mat(3,3,3,1,iv) = mat(3,3,3,1,iv) + one
          mat(3,4,4,1,iv) = mat(3,4,4,1,iv) + one
          mat(3,5,5,1,iv) = mat(3,5,5,1,iv) + one

!=======================================================================================================!
!.... use higher-order tangent on first node off the left boundary
!=======================================================================================================!
          do idof = 1, ndof
            do jdof = 1, ndof
              mat(2,idof,jdof,2,iv) = alfa * dtl(2,iv) * gb1 * dxiinv * Ah(idof,jdof,2,iv)
              mat(3,idof,jdof,2,iv) = alfa * dtl(2,iv) * gb2 * dxiinv * Ah(idof,jdof,2,iv) + &
                                      alfa * dtl(2,iv) * Dh(idof,jdof,2,iv)
              mat(4,idof,jdof,2,iv) = alfa * dtl(2,iv) * gb3 * dxiinv * Ah(idof,jdof,2,iv)
              mat(5,idof,jdof,2,iv) = alfa * dtl(2,iv) * gb4 * dxiinv * Ah(idof,jdof,2,iv)
              mat(1,idof,jdof,2,iv) = alfa * dtl(2,iv) * gb5 * dxiinv * Ah(idof,jdof,2,iv)
            end do
          end do
        
!.... \hat{V}_{\eta\eta} term

          mat(2,2,2,2,iv) = mat(2,2,2,2,iv) - alfa * dtl(2,iv) * db1 * dxisinv * Vh(1,2,iv)
          mat(2,2,3,2,iv) = mat(2,2,3,2,iv) - alfa * dtl(2,iv) * db1 * dxisinv * Vh(5,2,iv)
          mat(2,3,2,2,iv) = mat(2,3,2,2,iv) - alfa * dtl(2,iv) * db1 * dxisinv * Vh(6,2,iv)
          mat(2,3,3,2,iv) = mat(2,3,3,2,iv) - alfa * dtl(2,iv) * db1 * dxisinv * Vh(2,2,iv)
          mat(2,4,4,2,iv) = mat(2,4,4,2,iv) - alfa * dtl(2,iv) * db1 * dxisinv * Vh(3,2,iv)
          mat(2,5,5,2,iv) = mat(2,5,5,2,iv) - alfa * dtl(2,iv) * db1 * dxisinv * Vh(4,2,iv)
          
          mat(3,2,2,2,iv) = mat(3,2,2,2,iv) - alfa * dtl(2,iv) * db2 * dxisinv * Vh(1,2,iv)
          mat(3,2,3,2,iv) = mat(3,2,3,2,iv) - alfa * dtl(2,iv) * db2 * dxisinv * Vh(5,2,iv)
          mat(3,3,2,2,iv) = mat(3,3,2,2,iv) - alfa * dtl(2,iv) * db2 * dxisinv * Vh(6,2,iv)
          mat(3,3,3,2,iv) = mat(3,3,3,2,iv) - alfa * dtl(2,iv) * db2 * dxisinv * Vh(2,2,iv)
          mat(3,4,4,2,iv) = mat(3,4,4,2,iv) - alfa * dtl(2,iv) * db2 * dxisinv * Vh(3,2,iv)
          mat(3,5,5,2,iv) = mat(3,5,5,2,iv) - alfa * dtl(2,iv) * db2 * dxisinv * Vh(4,2,iv)
          
          mat(4,2,2,2,iv) = mat(4,2,2,2,iv) - alfa * dtl(2,iv) * db3 * dxisinv * Vh(1,2,iv)
          mat(4,2,3,2,iv) = mat(4,2,3,2,iv) - alfa * dtl(2,iv) * db3 * dxisinv * Vh(5,2,iv)
          mat(4,3,2,2,iv) = mat(4,3,2,2,iv) - alfa * dtl(2,iv) * db3 * dxisinv * Vh(6,2,iv)
          mat(4,3,3,2,iv) = mat(4,3,3,2,iv) - alfa * dtl(2,iv) * db3 * dxisinv * Vh(2,2,iv)
          mat(4,4,4,2,iv) = mat(4,4,4,2,iv) - alfa * dtl(2,iv) * db3 * dxisinv * Vh(3,2,iv)
          mat(4,5,5,2,iv) = mat(4,5,5,2,iv) - alfa * dtl(2,iv) * db3 * dxisinv * Vh(4,2,iv)
          
          mat(5,2,2,2,iv) = mat(5,2,2,2,iv) - alfa * dtl(2,iv) * db4 * dxisinv * Vh(1,2,iv)
          mat(5,2,3,2,iv) = mat(5,2,3,2,iv) - alfa * dtl(2,iv) * db4 * dxisinv * Vh(5,2,iv)
          mat(5,3,2,2,iv) = mat(5,3,2,2,iv) - alfa * dtl(2,iv) * db4 * dxisinv * Vh(6,2,iv)
          mat(5,3,3,2,iv) = mat(5,3,3,2,iv) - alfa * dtl(2,iv) * db4 * dxisinv * Vh(2,2,iv)
          mat(5,4,4,2,iv) = mat(5,4,4,2,iv) - alfa * dtl(2,iv) * db4 * dxisinv * Vh(3,2,iv)
          mat(5,5,5,2,iv) = mat(5,5,5,2,iv) - alfa * dtl(2,iv) * db4 * dxisinv * Vh(4,2,iv)
          
          mat(1,2,2,2,iv) = mat(1,2,2,2,iv) - alfa * dtl(2,iv) * db5 * dxisinv * Vh(1,2,iv)
          mat(1,2,3,2,iv) = mat(1,2,3,2,iv) - alfa * dtl(2,iv) * db5 * dxisinv * Vh(5,2,iv)
          mat(1,3,2,2,iv) = mat(1,3,2,2,iv) - alfa * dtl(2,iv) * db5 * dxisinv * Vh(6,2,iv)
          mat(1,3,3,2,iv) = mat(1,3,3,2,iv) - alfa * dtl(2,iv) * db5 * dxisinv * Vh(2,2,iv)
          mat(1,4,4,2,iv) = mat(1,4,4,2,iv) - alfa * dtl(2,iv) * db5 * dxisinv * Vh(3,2,iv)
          mat(1,5,5,2,iv) = mat(1,5,5,2,iv) - alfa * dtl(2,iv) * db5 * dxisinv * Vh(4,2,iv)

!.... I term

          mat(3,1,1,2,iv) = mat(3,1,1,2,iv) + one
          mat(3,2,2,2,iv) = mat(3,2,2,2,iv) + one
          mat(3,3,3,2,iv) = mat(3,3,3,2,iv) + one
          mat(3,4,4,2,iv) = mat(3,4,4,2,iv) + one
          mat(3,5,5,2,iv) = mat(3,5,5,2,iv) + one

!.... implicit damping term

          if (eps_e .ne. zero) then

            mat(2,1,1,2,iv) = mat(2,1,1,2,iv) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * (-one)
            mat(2,2,2,2,iv) = mat(2,2,2,2,iv) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * (-one)
            mat(2,3,3,2,iv) = mat(2,3,3,2,iv) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * (-one)
            mat(2,4,4,2,iv) = mat(2,4,4,2,iv) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * (-one)
            mat(2,5,5,2,iv) = mat(2,5,5,2,iv) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * (-one)
            
            mat(3,1,1,2,iv) = mat(3,1,1,2,iv) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * two
            mat(3,2,2,2,iv) = mat(3,2,2,2,iv) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * two
            mat(3,3,3,2,iv) = mat(3,3,3,2,iv) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * two
            mat(3,4,4,2,iv) = mat(3,4,4,2,iv) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * two
            mat(3,5,5,2,iv) = mat(3,5,5,2,iv) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * two
            
            mat(4,1,1,2,iv) = mat(4,1,1,2,iv) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * (-one)
            mat(4,2,2,2,iv) = mat(4,2,2,2,iv) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * (-one)
            mat(4,3,3,2,iv) = mat(4,3,3,2,iv) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * (-one)
            mat(4,4,4,2,iv) = mat(4,4,4,2,iv) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * (-one)
            mat(4,5,5,2,iv) = mat(4,5,5,2,iv) + alfa * dtl(2,iv) * eps_e * buff(2,iv) * (-one)
                
        end if

        end do

        end if          ! lsym


        if (rsym) then
        
          call lhsbs1f( mat, Ah, Dh, Vh, dtl, 2 )               ! applies a symmetry condition

        else

        !$omp parallel do private (idof,jdof)
        do iv = 1, ny

!=======================================================================================================!
!.... use higher-order tangent on first node off the right boundary
!=======================================================================================================!
          do idof = 1, ndof
            do jdof = 1, ndof
              mat(5,idof,jdof,nx-1,iv) = -alfa * dtl(nx-1,iv) * gb5 * dxiinv * Ah(idof,jdof,nx-1,iv)
              mat(1,idof,jdof,nx-1,iv) = -alfa * dtl(nx-1,iv) * gb4 * dxiinv * Ah(idof,jdof,nx-1,iv)
              mat(2,idof,jdof,nx-1,iv) = -alfa * dtl(nx-1,iv) * gb3 * dxiinv * Ah(idof,jdof,nx-1,iv)
              mat(3,idof,jdof,nx-1,iv) = -alfa * dtl(nx-1,iv) * gb2 * dxiinv * Ah(idof,jdof,nx-1,iv) + &
                                          alfa * dtl(nx-1,iv) * Dh(idof,jdof,nx-1,iv)
              mat(4,idof,jdof,nx-1,iv) = -alfa * dtl(nx-1,iv) * gb1 * dxiinv * Ah(idof,jdof,nx-1,iv)
            end do
          end do
        
!.... \hat{V}_{\eta\eta} term

          mat(5,2,2,nx-1,iv) = mat(5,2,2,nx-1,iv) - alfa * dtl(nx-1,iv) * db5 * dxisinv * Vh(1,nx-1,iv)
          mat(5,2,3,nx-1,iv) = mat(5,2,3,nx-1,iv) - alfa * dtl(nx-1,iv) * db5 * dxisinv * Vh(5,nx-1,iv)
          mat(5,3,2,nx-1,iv) = mat(5,3,2,nx-1,iv) - alfa * dtl(nx-1,iv) * db5 * dxisinv * Vh(6,nx-1,iv)
          mat(5,3,3,nx-1,iv) = mat(5,3,3,nx-1,iv) - alfa * dtl(nx-1,iv) * db5 * dxisinv * Vh(2,nx-1,iv)
          mat(5,4,4,nx-1,iv) = mat(5,4,4,nx-1,iv) - alfa * dtl(nx-1,iv) * db5 * dxisinv * Vh(3,nx-1,iv)
          mat(5,5,5,nx-1,iv) = mat(5,5,5,nx-1,iv) - alfa * dtl(nx-1,iv) * db5 * dxisinv * Vh(4,nx-1,iv)

          mat(1,2,2,nx-1,iv) = mat(1,2,2,nx-1,iv) - alfa * dtl(nx-1,iv) * db4 * dxisinv * Vh(1,nx-1,iv)
          mat(1,2,3,nx-1,iv) = mat(1,2,3,nx-1,iv) - alfa * dtl(nx-1,iv) * db4 * dxisinv * Vh(5,nx-1,iv)
          mat(1,3,2,nx-1,iv) = mat(1,3,2,nx-1,iv) - alfa * dtl(nx-1,iv) * db4 * dxisinv * Vh(6,nx-1,iv)
          mat(1,3,3,nx-1,iv) = mat(1,3,3,nx-1,iv) - alfa * dtl(nx-1,iv) * db4 * dxisinv * Vh(2,nx-1,iv)
          mat(1,4,4,nx-1,iv) = mat(1,4,4,nx-1,iv) - alfa * dtl(nx-1,iv) * db4 * dxisinv * Vh(3,nx-1,iv)
          mat(1,5,5,nx-1,iv) = mat(1,5,5,nx-1,iv) - alfa * dtl(nx-1,iv) * db4 * dxisinv * Vh(4,nx-1,iv)
          
          mat(2,2,2,nx-1,iv) = mat(2,2,2,nx-1,iv) - alfa * dtl(nx-1,iv) * db3 * dxisinv * Vh(1,nx-1,iv)
          mat(2,2,3,nx-1,iv) = mat(2,2,3,nx-1,iv) - alfa * dtl(nx-1,iv) * db3 * dxisinv * Vh(5,nx-1,iv)
          mat(2,3,2,nx-1,iv) = mat(2,3,2,nx-1,iv) - alfa * dtl(nx-1,iv) * db3 * dxisinv * Vh(6,nx-1,iv)
          mat(2,3,3,nx-1,iv) = mat(2,3,3,nx-1,iv) - alfa * dtl(nx-1,iv) * db3 * dxisinv * Vh(2,nx-1,iv)
          mat(2,4,4,nx-1,iv) = mat(2,4,4,nx-1,iv) - alfa * dtl(nx-1,iv) * db3 * dxisinv * Vh(3,nx-1,iv)
          mat(2,5,5,nx-1,iv) = mat(2,5,5,nx-1,iv) - alfa * dtl(nx-1,iv) * db3 * dxisinv * Vh(4,nx-1,iv)
          
          mat(3,2,2,nx-1,iv) = mat(3,2,2,nx-1,iv) - alfa * dtl(nx-1,iv) * db2 * dxisinv * Vh(1,nx-1,iv)
          mat(3,2,3,nx-1,iv) = mat(3,2,3,nx-1,iv) - alfa * dtl(nx-1,iv) * db2 * dxisinv * Vh(5,nx-1,iv)
          mat(3,3,2,nx-1,iv) = mat(3,3,2,nx-1,iv) - alfa * dtl(nx-1,iv) * db2 * dxisinv * Vh(6,nx-1,iv)
          mat(3,3,3,nx-1,iv) = mat(3,3,3,nx-1,iv) - alfa * dtl(nx-1,iv) * db2 * dxisinv * Vh(2,nx-1,iv)
          mat(3,4,4,nx-1,iv) = mat(3,4,4,nx-1,iv) - alfa * dtl(nx-1,iv) * db2 * dxisinv * Vh(3,nx-1,iv)
          mat(3,5,5,nx-1,iv) = mat(3,5,5,nx-1,iv) - alfa * dtl(nx-1,iv) * db2 * dxisinv * Vh(4,nx-1,iv)
          
          mat(4,2,2,nx-1,iv) = mat(4,2,2,nx-1,iv) - alfa * dtl(nx-1,iv) * db1 * dxisinv * Vh(1,nx-1,iv)
          mat(4,2,3,nx-1,iv) = mat(4,2,3,nx-1,iv) - alfa * dtl(nx-1,iv) * db1 * dxisinv * Vh(5,nx-1,iv)
          mat(4,3,2,nx-1,iv) = mat(4,3,2,nx-1,iv) - alfa * dtl(nx-1,iv) * db1 * dxisinv * Vh(6,nx-1,iv)
          mat(4,3,3,nx-1,iv) = mat(4,3,3,nx-1,iv) - alfa * dtl(nx-1,iv) * db1 * dxisinv * Vh(2,nx-1,iv)
          mat(4,4,4,nx-1,iv) = mat(4,4,4,nx-1,iv) - alfa * dtl(nx-1,iv) * db1 * dxisinv * Vh(3,nx-1,iv)
          mat(4,5,5,nx-1,iv) = mat(4,5,5,nx-1,iv) - alfa * dtl(nx-1,iv) * db1 * dxisinv * Vh(4,nx-1,iv)

!.... I term

          mat(3,1,1,nx-1,iv) = mat(3,1,1,nx-1,iv) + one
          mat(3,2,2,nx-1,iv) = mat(3,2,2,nx-1,iv) + one
          mat(3,3,3,nx-1,iv) = mat(3,3,3,nx-1,iv) + one
          mat(3,4,4,nx-1,iv) = mat(3,4,4,nx-1,iv) + one
          mat(3,5,5,nx-1,iv) = mat(3,5,5,nx-1,iv) + one

!.... implicit damping term

          if (eps_e .ne. zero) then

            mat(2,1,1,nx-1,iv) = mat(2,1,1,nx-1,iv) + alfa * dtl(nx-1,iv) * eps_e * buff(nx-1,iv) * (-one)
            mat(2,2,2,nx-1,iv) = mat(2,2,2,nx-1,iv) + alfa * dtl(nx-1,iv) * eps_e * buff(nx-1,iv) * (-one)
            mat(2,3,3,nx-1,iv) = mat(2,3,3,nx-1,iv) + alfa * dtl(nx-1,iv) * eps_e * buff(nx-1,iv) * (-one)
            mat(2,4,4,nx-1,iv) = mat(2,4,4,nx-1,iv) + alfa * dtl(nx-1,iv) * eps_e * buff(nx-1,iv) * (-one)
            mat(2,5,5,nx-1,iv) = mat(2,5,5,nx-1,iv) + alfa * dtl(nx-1,iv) * eps_e * buff(nx-1,iv) * (-one)

            mat(3,1,1,nx-1,iv) = mat(3,1,1,nx-1,iv) + alfa * dtl(nx-1,iv) * eps_e * buff(nx-1,iv) * two
            mat(3,2,2,nx-1,iv) = mat(3,2,2,nx-1,iv) + alfa * dtl(nx-1,iv) * eps_e * buff(nx-1,iv) * two
            mat(3,3,3,nx-1,iv) = mat(3,3,3,nx-1,iv) + alfa * dtl(nx-1,iv) * eps_e * buff(nx-1,iv) * two
            mat(3,4,4,nx-1,iv) = mat(3,4,4,nx-1,iv) + alfa * dtl(nx-1,iv) * eps_e * buff(nx-1,iv) * two
            mat(3,5,5,nx-1,iv) = mat(3,5,5,nx-1,iv) + alfa * dtl(nx-1,iv) * eps_e * buff(nx-1,iv) * two
            
            mat(4,1,1,nx-1,iv) = mat(4,1,1,nx-1,iv) + alfa * dtl(nx-1,iv) * eps_e * buff(nx-1,iv) * (-one)
            mat(4,2,2,nx-1,iv) = mat(4,2,2,nx-1,iv) + alfa * dtl(nx-1,iv) * eps_e * buff(nx-1,iv) * (-one)
            mat(4,3,3,nx-1,iv) = mat(4,3,3,nx-1,iv) + alfa * dtl(nx-1,iv) * eps_e * buff(nx-1,iv) * (-one)
            mat(4,4,4,nx-1,iv) = mat(4,4,4,nx-1,iv) + alfa * dtl(nx-1,iv) * eps_e * buff(nx-1,iv) * (-one)
            mat(4,5,5,nx-1,iv) = mat(4,5,5,nx-1,iv) + alfa * dtl(nx-1,iv) * eps_e * buff(nx-1,iv) * (-one)

          end if

!=======================================================================================================!
!.... use higher-order tangent on the right boundary
!=======================================================================================================!
          do idof = 1, ndof
            do jdof = 1, ndof
              mat(4,idof,jdof,nx,iv) = -alfa * dtl(nx,iv) * gc5 * dxiinv * Ah(idof,jdof,nx,iv)
              mat(5,idof,jdof,nx,iv) = -alfa * dtl(nx,iv) * gc4 * dxiinv * Ah(idof,jdof,nx,iv)
              mat(1,idof,jdof,nx,iv) = -alfa * dtl(nx,iv) * gc3 * dxiinv * Ah(idof,jdof,nx,iv)
              mat(2,idof,jdof,nx,iv) = -alfa * dtl(nx,iv) * gc2 * dxiinv * Ah(idof,jdof,nx,iv)
              mat(3,idof,jdof,nx,iv) = -alfa * dtl(nx,iv) * gc1 * dxiinv * Ah(idof,jdof,nx,iv) + &
                                        alfa * dtl(nx,iv) * Dh(idof,jdof,nx,iv)
            end do
          end do
        
!.... \hat{V}_{\eta\eta} term

          mat(4,2,2,nx,iv) = mat(4,2,2,nx,iv) - alfa * dtl(nx,iv) * dd5 * dxisinv * Vh(1,nx,iv)
          mat(4,2,3,nx,iv) = mat(4,2,3,nx,iv) - alfa * dtl(nx,iv) * dd5 * dxisinv * Vh(5,nx,iv)
          mat(4,3,2,nx,iv) = mat(4,3,2,nx,iv) - alfa * dtl(nx,iv) * dd5 * dxisinv * Vh(6,nx,iv)
          mat(4,3,3,nx,iv) = mat(4,3,3,nx,iv) - alfa * dtl(nx,iv) * dd5 * dxisinv * Vh(2,nx,iv)
          mat(4,4,4,nx,iv) = mat(4,4,4,nx,iv) - alfa * dtl(nx,iv) * dd5 * dxisinv * Vh(3,nx,iv)
          mat(4,5,5,nx,iv) = mat(4,5,5,nx,iv) - alfa * dtl(nx,iv) * dd5 * dxisinv * Vh(4,nx,iv)
          
          mat(5,2,2,nx,iv) = mat(5,2,2,nx,iv) - alfa * dtl(nx,iv) * dd4 * dxisinv * Vh(1,nx,iv)
          mat(5,2,3,nx,iv) = mat(5,2,3,nx,iv) - alfa * dtl(nx,iv) * dd4 * dxisinv * Vh(5,nx,iv)
          mat(5,3,2,nx,iv) = mat(5,3,2,nx,iv) - alfa * dtl(nx,iv) * dd4 * dxisinv * Vh(6,nx,iv)
          mat(5,3,3,nx,iv) = mat(5,3,3,nx,iv) - alfa * dtl(nx,iv) * dd4 * dxisinv * Vh(2,nx,iv)
          mat(5,4,4,nx,iv) = mat(5,4,4,nx,iv) - alfa * dtl(nx,iv) * dd4 * dxisinv * Vh(3,nx,iv)
          mat(5,5,5,nx,iv) = mat(5,5,5,nx,iv) - alfa * dtl(nx,iv) * dd4 * dxisinv * Vh(4,nx,iv)
          
          mat(1,2,2,nx,iv) = mat(1,2,2,nx,iv) - alfa * dtl(nx,iv) * dd3 * dxisinv * Vh(1,nx,iv)
          mat(1,2,3,nx,iv) = mat(1,2,3,nx,iv) - alfa * dtl(nx,iv) * dd3 * dxisinv * Vh(5,nx,iv)
          mat(1,3,2,nx,iv) = mat(1,3,2,nx,iv) - alfa * dtl(nx,iv) * dd3 * dxisinv * Vh(6,nx,iv)
          mat(1,3,3,nx,iv) = mat(1,3,3,nx,iv) - alfa * dtl(nx,iv) * dd3 * dxisinv * Vh(2,nx,iv)
          mat(1,4,4,nx,iv) = mat(1,4,4,nx,iv) - alfa * dtl(nx,iv) * dd3 * dxisinv * Vh(3,nx,iv)
          mat(1,5,5,nx,iv) = mat(1,5,5,nx,iv) - alfa * dtl(nx,iv) * dd3 * dxisinv * Vh(4,nx,iv)
          
          mat(2,2,2,nx,iv) = mat(2,2,2,nx,iv) - alfa * dtl(nx,iv) * dd2 * dxisinv * Vh(1,nx,iv)
          mat(2,2,3,nx,iv) = mat(2,2,3,nx,iv) - alfa * dtl(nx,iv) * dd2 * dxisinv * Vh(5,nx,iv)
          mat(2,3,2,nx,iv) = mat(2,3,2,nx,iv) - alfa * dtl(nx,iv) * dd2 * dxisinv * Vh(6,nx,iv)
          mat(2,3,3,nx,iv) = mat(2,3,3,nx,iv) - alfa * dtl(nx,iv) * dd2 * dxisinv * Vh(2,nx,iv)
          mat(2,4,4,nx,iv) = mat(2,4,4,nx,iv) - alfa * dtl(nx,iv) * dd2 * dxisinv * Vh(3,nx,iv)
          mat(2,5,5,nx,iv) = mat(2,5,5,nx,iv) - alfa * dtl(nx,iv) * dd2 * dxisinv * Vh(4,nx,iv)
          
          mat(3,2,2,nx,iv) = mat(3,2,2,nx,iv) - alfa * dtl(nx,iv) * dd1 * dxisinv * Vh(1,nx,iv)
          mat(3,2,3,nx,iv) = mat(3,2,3,nx,iv) - alfa * dtl(nx,iv) * dd1 * dxisinv * Vh(5,nx,iv)
          mat(3,3,2,nx,iv) = mat(3,3,2,nx,iv) - alfa * dtl(nx,iv) * dd1 * dxisinv * Vh(6,nx,iv)
          mat(3,3,3,nx,iv) = mat(3,3,3,nx,iv) - alfa * dtl(nx,iv) * dd1 * dxisinv * Vh(2,nx,iv)
          mat(3,4,4,nx,iv) = mat(3,4,4,nx,iv) - alfa * dtl(nx,iv) * dd1 * dxisinv * Vh(3,nx,iv)
          mat(3,5,5,nx,iv) = mat(3,5,5,nx,iv) - alfa * dtl(nx,iv) * dd1 * dxisinv * Vh(4,nx,iv)

!.... I term

          mat(3,1,1,nx,iv) = mat(3,1,1,nx,iv) + one
          mat(3,2,2,nx,iv) = mat(3,2,2,nx,iv) + one
          mat(3,3,3,nx,iv) = mat(3,3,3,nx,iv) + one
          mat(3,4,4,nx,iv) = mat(3,4,4,nx,iv) + one
          mat(3,5,5,nx,iv) = mat(3,5,5,nx,iv) + one

        end do

        end if         ! rsym

        return
        end subroutine lhsbt1f
