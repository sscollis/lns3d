!=============================================================================!
        subroutine lhs2f3D( mat, Bh, Dh, Dhi, Vh22, ABhi, dtl, vl)
!  
!  Form the LHS for the \eta direction for the three-dimensional
!  Navier-Stokes Equations using 4th order differencing.
!
!=============================================================================!
        use global
        use buffer
        use stencil
        implicit none
        
        complex :: mat(5,ndof,ndof,nx,ny), vl(ndof,nx,ny)
        real    :: Bh(ndof,ndof,nx,ny), Dh(ndof,ndof,nx,ny)
        real    :: Dhi(ndof,ndof,nx,ny), Vh22(6,nx,ny), ABhi(6,nx,ny)
        real    :: dtl(nx,ny)

!.... first derivative operator

        real    :: a1, a2, a3, a4, a5
        real    :: b1, b2, b3, b4, b5
        complex :: c1, c2, c3, c4, c5
        
        integer :: i, j, idof, jdof
        logical :: calcd
!=============================================================================!
        if (useCalcd) then
          if (mod(real(istep-1+iter),two) .eq. zero ) then
            calcd = .false.
          else
            calcd = .true.
          end if
          write(88,*) 'lhs2f3D: ',istep, iter, calcd
          !write(*,*) 'lhs2f3D: ',istep, iter, calcd
        else
          calcd = .false.
        endif 

!.... fourth-order stencil

        a1 = alfa * ga1  / deta
        a2 = alfa * ga2  / deta
        a3 = alfa * zero / deta
        a4 = alfa * ga3  / deta
        a5 = alfa * ga4  / deta

        b1 = alfa * da1 / deta**2
        b2 = alfa * da2 / deta**2
        b3 = alfa * da3 / deta**2
        b4 = alfa * da4 / deta**2
        b5 = alfa * da5 / deta**2

        c1 = im * alfa * ga1  / deta
        c2 = im * alfa * ga2  / deta
        c3 = im * alfa * zero / deta
        c4 = im * alfa * ga3  / deta
        c5 = im * alfa * ga4  / deta

        !$omp parallel do private (i,idof,jdof)
        loop_j: do j = 1, ny
        loop_i: do i = 1, nx

!.... \hat{B} term

        do idof = 1, ndof
          do jdof = 1, ndof
            mat(1,idof,jdof,i,j) = a1 * dtl(i,j) * Bh(idof,jdof,i,j)    
            mat(2,idof,jdof,i,j) = a2 * dtl(i,j) * Bh(idof,jdof,i,j)    
            mat(3,idof,jdof,i,j) = a3 * dtl(i,j) * Bh(idof,jdof,i,j)
            mat(4,idof,jdof,i,j) = a4 * dtl(i,j) * Bh(idof,jdof,i,j)
            mat(5,idof,jdof,i,j) = a5 * dtl(i,j) * Bh(idof,jdof,i,j)
          end do
        end do

!.... \hat{B}_i term

        mat(1,2,4,i,j) = mat(1,2,4,i,j) - c1 * dtl(i,j) * ABhi(3,i,j)
        mat(1,3,4,i,j) = mat(1,3,4,i,j) - c1 * dtl(i,j) * ABhi(4,i,j)
        mat(1,4,2,i,j) = mat(1,4,2,i,j) - c1 * dtl(i,j) * ABhi(3,i,j)
        mat(1,4,3,i,j) = mat(1,4,3,i,j) - c1 * dtl(i,j) * ABhi(4,i,j)

        mat(2,2,4,i,j) = mat(2,2,4,i,j) - c2 * dtl(i,j) * ABhi(3,i,j)
        mat(2,3,4,i,j) = mat(2,3,4,i,j) - c2 * dtl(i,j) * ABhi(4,i,j)
        mat(2,4,2,i,j) = mat(2,4,2,i,j) - c2 * dtl(i,j) * ABhi(3,i,j)
        mat(2,4,3,i,j) = mat(2,4,3,i,j) - c2 * dtl(i,j) * ABhi(4,i,j)

        mat(3,2,4,i,j) = mat(3,2,4,i,j) - c3 * dtl(i,j) * ABhi(3,i,j)
        mat(3,3,4,i,j) = mat(3,3,4,i,j) - c3 * dtl(i,j) * ABhi(4,i,j)
        mat(3,4,2,i,j) = mat(3,4,2,i,j) - c3 * dtl(i,j) * ABhi(3,i,j)
        mat(3,4,3,i,j) = mat(3,4,3,i,j) - c3 * dtl(i,j) * ABhi(4,i,j)

        mat(4,2,4,i,j) = mat(4,2,4,i,j) - c4 * dtl(i,j) * ABhi(3,i,j)
        mat(4,3,4,i,j) = mat(4,3,4,i,j) - c4 * dtl(i,j) * ABhi(4,i,j)
        mat(4,4,2,i,j) = mat(4,4,2,i,j) - c4 * dtl(i,j) * ABhi(3,i,j)
        mat(4,4,3,i,j) = mat(4,4,3,i,j) - c4 * dtl(i,j) * ABhi(4,i,j)

        mat(5,2,4,i,j) = mat(5,2,4,i,j) - c5 * dtl(i,j) * ABhi(3,i,j)
        mat(5,3,4,i,j) = mat(5,3,4,i,j) - c5 * dtl(i,j) * ABhi(4,i,j)
        mat(5,4,2,i,j) = mat(5,4,2,i,j) - c5 * dtl(i,j) * ABhi(3,i,j)
        mat(5,4,3,i,j) = mat(5,4,3,i,j) - c5 * dtl(i,j) * ABhi(4,i,j)

        if (calcd) then

!.... \hat{D} term

          do idof = 1, ndof
            do jdof = 1, ndof
              mat(3,idof,jdof,i,j) = mat(3,idof,jdof,i,j) + &
                                     alfa * dtl(i,j) * Dh(idof,jdof,i,j)
            end do
          end do

!.... \hat{D}_i term

          do idof = 1, ndof
            do jdof = 1, ndof
              mat(3,idof,jdof,i,j) = mat(3,idof,jdof,i,j) + &
                                     im * alfa * dtl(i,j) * Dhi(idof,jdof,i,j)
            end do
          end do

        end if

!.... \hat{V}_{\eta\eta} term

        mat(1,2,2,i,j) = mat(1,2,2,i,j) - b1 * dtl(i,j) * Vh22(1,i,j)
        mat(1,2,3,i,j) = mat(1,2,3,i,j) - b1 * dtl(i,j) * Vh22(5,i,j)
        mat(1,3,2,i,j) = mat(1,3,2,i,j) - b1 * dtl(i,j) * Vh22(6,i,j)
        mat(1,3,3,i,j) = mat(1,3,3,i,j) - b1 * dtl(i,j) * Vh22(2,i,j)
        mat(1,4,4,i,j) = mat(1,4,4,i,j) - b1 * dtl(i,j) * Vh22(3,i,j)
        mat(1,5,5,i,j) = mat(1,5,5,i,j) - b1 * dtl(i,j) * Vh22(4,i,j)

        mat(2,2,2,i,j) = mat(2,2,2,i,j) - b2 * dtl(i,j) * Vh22(1,i,j)
        mat(2,2,3,i,j) = mat(2,2,3,i,j) - b2 * dtl(i,j) * Vh22(5,i,j)
        mat(2,3,2,i,j) = mat(2,3,2,i,j) - b2 * dtl(i,j) * Vh22(6,i,j)
        mat(2,3,3,i,j) = mat(2,3,3,i,j) - b2 * dtl(i,j) * Vh22(2,i,j)
        mat(2,4,4,i,j) = mat(2,4,4,i,j) - b2 * dtl(i,j) * Vh22(3,i,j)
        mat(2,5,5,i,j) = mat(2,5,5,i,j) - b2 * dtl(i,j) * Vh22(4,i,j)

        mat(3,2,2,i,j) = mat(3,2,2,i,j) - b3 * dtl(i,j) * Vh22(1,i,j)
        mat(3,2,3,i,j) = mat(3,2,3,i,j) - b3 * dtl(i,j) * Vh22(5,i,j)
        mat(3,3,2,i,j) = mat(3,3,2,i,j) - b3 * dtl(i,j) * Vh22(6,i,j)
        mat(3,3,3,i,j) = mat(3,3,3,i,j) - b3 * dtl(i,j) * Vh22(2,i,j)
        mat(3,4,4,i,j) = mat(3,4,4,i,j) - b3 * dtl(i,j) * Vh22(3,i,j)
        mat(3,5,5,i,j) = mat(3,5,5,i,j) - b3 * dtl(i,j) * Vh22(4,i,j)

        mat(4,2,2,i,j) = mat(4,2,2,i,j) - b4 * dtl(i,j) * Vh22(1,i,j)
        mat(4,2,3,i,j) = mat(4,2,3,i,j) - b4 * dtl(i,j) * Vh22(5,i,j)
        mat(4,3,2,i,j) = mat(4,3,2,i,j) - b4 * dtl(i,j) * Vh22(6,i,j)
        mat(4,3,3,i,j) = mat(4,3,3,i,j) - b4 * dtl(i,j) * Vh22(2,i,j)
        mat(4,4,4,i,j) = mat(4,4,4,i,j) - b4 * dtl(i,j) * Vh22(3,i,j)
        mat(4,5,5,i,j) = mat(4,5,5,i,j) - b4 * dtl(i,j) * Vh22(4,i,j)

        mat(5,2,2,i,j) = mat(5,2,2,i,j) - b5 * dtl(i,j) * Vh22(1,i,j)
        mat(5,2,3,i,j) = mat(5,2,3,i,j) - b5 * dtl(i,j) * Vh22(5,i,j)
        mat(5,3,2,i,j) = mat(5,3,2,i,j) - b5 * dtl(i,j) * Vh22(6,i,j)
        mat(5,3,3,i,j) = mat(5,3,3,i,j) - b5 * dtl(i,j) * Vh22(2,i,j)
        mat(5,4,4,i,j) = mat(5,4,4,i,j) - b5 * dtl(i,j) * Vh22(3,i,j)
        mat(5,5,5,i,j) = mat(5,5,5,i,j) - b5 * dtl(i,j) * Vh22(4,i,j)

!.... I term

        mat(3,1,1,i,j) = mat(3,1,1,i,j) + one
        mat(3,2,2,i,j) = mat(3,2,2,i,j) + one
        mat(3,3,3,i,j) = mat(3,3,3,i,j) + one
        mat(3,4,4,i,j) = mat(3,4,4,i,j) + one
        mat(3,5,5,i,j) = mat(3,5,5,i,j) + one

!.... implicit damping term updated to fourth-order

        if (eps_e .ne. zero) then

        mat(1,1,1,i,j) = mat(1,1,1,i,j) + alfa * dtl(i,j)*eps_e*buff(i,j)*fb1
        mat(1,2,2,i,j) = mat(1,2,2,i,j) + alfa * dtl(i,j)*eps_e*buff(i,j)*fb1
        mat(1,3,3,i,j) = mat(1,3,3,i,j) + alfa * dtl(i,j)*eps_e*buff(i,j)*fb1
        mat(1,4,4,i,j) = mat(1,4,4,i,j) + alfa * dtl(i,j)*eps_e*buff(i,j)*fb1
        mat(1,5,5,i,j) = mat(1,5,5,i,j) + alfa * dtl(i,j)*eps_e*buff(i,j)*fb1

        mat(2,1,1,i,j) = mat(2,1,1,i,j) + alfa * dtl(i,j)*eps_e*buff(i,j)*fb2
        mat(2,2,2,i,j) = mat(2,2,2,i,j) + alfa * dtl(i,j)*eps_e*buff(i,j)*fb2
        mat(2,3,3,i,j) = mat(2,3,3,i,j) + alfa * dtl(i,j)*eps_e*buff(i,j)*fb2
        mat(2,4,4,i,j) = mat(2,4,4,i,j) + alfa * dtl(i,j)*eps_e*buff(i,j)*fb2
        mat(2,5,5,i,j) = mat(2,5,5,i,j) + alfa * dtl(i,j)*eps_e*buff(i,j)*fb2

        mat(3,1,1,i,j) = mat(3,1,1,i,j) + alfa * dtl(i,j)*eps_e*buff(i,j)*fb3
        mat(3,2,2,i,j) = mat(3,2,2,i,j) + alfa * dtl(i,j)*eps_e*buff(i,j)*fb3
        mat(3,3,3,i,j) = mat(3,3,3,i,j) + alfa * dtl(i,j)*eps_e*buff(i,j)*fb3
        mat(3,4,4,i,j) = mat(3,4,4,i,j) + alfa * dtl(i,j)*eps_e*buff(i,j)*fb3
        mat(3,5,5,i,j) = mat(3,5,5,i,j) + alfa * dtl(i,j)*eps_e*buff(i,j)*fb3

        mat(4,1,1,i,j) = mat(4,1,1,i,j) + alfa * dtl(i,j)*eps_e*buff(i,j)*fb4
        mat(4,2,2,i,j) = mat(4,2,2,i,j) + alfa * dtl(i,j)*eps_e*buff(i,j)*fb4
        mat(4,3,3,i,j) = mat(4,3,3,i,j) + alfa * dtl(i,j)*eps_e*buff(i,j)*fb4
        mat(4,4,4,i,j) = mat(4,4,4,i,j) + alfa * dtl(i,j)*eps_e*buff(i,j)*fb4
        mat(4,5,5,i,j) = mat(4,5,5,i,j) + alfa * dtl(i,j)*eps_e*buff(i,j)*fb4

        mat(5,1,1,i,j) = mat(5,1,1,i,j) + alfa * dtl(i,j)*eps_e*buff(i,j)*fb5
        mat(5,2,2,i,j) = mat(5,2,2,i,j) + alfa * dtl(i,j)*eps_e*buff(i,j)*fb5
        mat(5,3,3,i,j) = mat(5,3,3,i,j) + alfa * dtl(i,j)*eps_e*buff(i,j)*fb5
        mat(5,4,4,i,j) = mat(5,4,4,i,j) + alfa * dtl(i,j)*eps_e*buff(i,j)*fb5
        mat(5,5,5,i,j) = mat(5,5,5,i,j) + alfa * dtl(i,j)*eps_e*buff(i,j)*fb5

        end if

!.... sponge

        if (.not. calcd) then

        if (ispg .eq. 1) then
        
          mat(3,1,1,i,j) = mat(3,1,1,i,j) + alfa * dtl(i,j) * spg(i,j)
          mat(3,2,2,i,j) = mat(3,2,2,i,j) + alfa * dtl(i,j) * spg(i,j)
          mat(3,3,3,i,j) = mat(3,3,3,i,j) + alfa * dtl(i,j) * spg(i,j)
          mat(3,4,4,i,j) = mat(3,4,4,i,j) + alfa * dtl(i,j) * spg(i,j)
          mat(3,5,5,i,j) = mat(3,5,5,i,j) + alfa * dtl(i,j) * spg(i,j)
        
        else if (ispg .ge. 2) then
        
          mat(3,1,1,i,j) = mat(3,1,1,i,j) + alfa*dtl(i,j)*(spg(i,j)+spg2(i,j))
          mat(3,2,2,i,j) = mat(3,2,2,i,j) + alfa*dtl(i,j)*(spg(i,j)+spg2(i,j))
          mat(3,3,3,i,j) = mat(3,3,3,i,j) + alfa*dtl(i,j)*(spg(i,j)+spg2(i,j))
          mat(3,4,4,i,j) = mat(3,4,4,i,j) + alfa*dtl(i,j)*(spg(i,j)+spg2(i,j))
          mat(3,5,5,i,j) = mat(3,5,5,i,j) + alfa*dtl(i,j)*(spg(i,j)+spg2(i,j))
        
        end if

        end if

        end do loop_i
        end do loop_j

!.... apply boundary treatment to the LHS

        call lhsbt2f3D( mat, Bh, Dh, Dhi, Vh22, ABhi, spg, spg2, dtl, calcd )

!.... apply boundary conditions to the LHS

        call lhsbc2f3D( mat, vl, vm )

        return
        end
