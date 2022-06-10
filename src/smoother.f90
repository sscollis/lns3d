!=============================================================================!
        subroutine smoother(rl, vl) 
!  
!       Fourth order explicit smoother for the NS equations.
!       From: Computational Fluid Mechanics and Heat Transfer
!             D.A. Anderson, J. C. Tannehill, Rl H. Pletcher
!             Page:  450, 495
!
!       Note this routine can also do sixth-order and second-order smoothers
!       but they are currently disabled.
!
!       Written: 6-9-95
!
!       Revised: 7-9-96
!
!       Notes:   Need to account for periodicity in the smoothed fields.
!
!       Revised: 11-20-99 Switched indices
!=============================================================================!
        use global
        use buffer
        use stencil
        implicit none

        real :: rl(ndof,nx,ny), vl(ndof,nx,ny)
        !$sgi distribute rl(*,*,block), vl(*,*,block)

        real :: isign
        integer :: idof, i, j

!=============================================================================!

        if (eps_e .eq. zero) return

!===========================================================================
!     F O U R T H   O R D E R   T E R M 
!===========================================================================

!.... \eta direction
        
        !$omp parallel do private(idof, j)
        do i = 1, nx
          do idof = 1, ndof

            j = 2
            rl(idof,i,j) = rl(idof,i,j) + eps_e * buff(i,j) * ( &
                                           -1.0 * vl(idof,i,j-1) + &
                                            2.0 * vl(idof,i,j  ) - &
                                            1.0 * vl(idof,i,j+1) )
            
            j = 3
            rl(idof,i,j) = rl(idof,i,j) + eps_e * buff(i,j) * ( &
                                                vl(idof,i,j-2) - &
                                          4.0 * vl(idof,i,j-1) + &
                                          6.0 * vl(idof,i,j  ) - &
                                          4.0 * vl(idof,i,j+1) + &
                                          1.0 * vl(idof,i,j+2) )

            do j = 4, ny-3
              rl(idof,i,j) = rl(idof,i,j) + eps_e * buff(i,j) * ( &
                                            fa1 * vl(idof,i,j-3) + &
                                            fa2 * vl(idof,i,j-2) + &
                                            fa3 * vl(idof,i,j-1) + &
                                            fa4 * vl(idof,i,j  ) + &
                                            fa5 * vl(idof,i,j+1) + &
                                            fa6 * vl(idof,i,j+2) + &
                                            fa7 * vl(idof,i,j+3) )
            end do

            j = ny-2
            rl(idof,i,j) = rl(idof,i,j) + eps_e * buff(i,j) * ( &
                                                vl(idof,i,j-2) - &
                                          4.0 * vl(idof,i,j-1) + &
                                          6.0 * vl(idof,i,j  ) - &
                                          4.0 * vl(idof,i,j+1) + &
                                          1.0 * vl(idof,i,j+2) )

            j = ny-1
            rl(idof,i,j) = rl(idof,i,j) + eps_e * buff(i,j) * ( &
                                           -1.0 * vl(idof,i,j-1) + &
                                            2.0 * vl(idof,i,j  ) - &
                                            1.0 * vl(idof,i,j+1) )

          end do
        end do
                                          
!.... \xi direction

        if (lsym) then
          !$omp parallel do private(idof, isign, i)
          do j = 1, ny
            do idof = 1, ndof
              if (idof .eq. 3) then
                isign = -one
              else
                isign = one
              end if

              i = 1
              rl(idof,i,j) = rl(idof,i,j) + eps_e * buff(i,j) * ( &
                                            isign * fa1 * vl(idof,4,j) + &
                                            isign * fa2 * vl(idof,3,j) + &
                                            isign * fa3 * vl(idof,2,j) + &
                                            fa4 * vl(idof,1,j) + &
                                            fa5 * vl(idof,2,j) + &
                                            fa6 * vl(idof,3,j) + &
                                            fa7 * vl(idof,4,j) )
              
              i = 2
              rl(idof,i,j) = rl(idof,i,j) + eps_e * buff(i,j) * ( &
                                            isign * fa1 * vl(idof,3,j) + &
                                            isign * fa2 * vl(idof,2,j) + &
                                            fa3 * vl(idof,1,j) + &
                                            fa4 * vl(idof,2,j) + &
                                            fa5 * vl(idof,3,j) + &
                                            fa6 * vl(idof,4,j) + &
                                            fa7 * vl(idof,5,j) )

              i = 3
              rl(idof,i,j) = rl(idof,i,j) + eps_e * buff(i,j) * ( &
                                            isign * fa1 * vl(idof,2,j) + &
                                            fa2 * vl(idof,1,j) + &
                                            fa3 * vl(idof,2,j) + &
                                            fa4 * vl(idof,3,j) + &
                                            fa5 * vl(idof,4,j) + &
                                            fa6 * vl(idof,5,j) + &
                                            fa7 * vl(idof,6,j) )
            end do
          end do
        else
          !$omp parallel do private(idof, i)
          do j = 1, ny
            do idof = 1, ndof
              i = 2
              rl(idof,i,j) = rl(idof,i,j) + eps_e * buff(i,j) * ( &
                                           -1.0 * vl(idof,i-1,j) + &
                                            2.0 * vl(idof,i,j  ) - &
                                            1.0 * vl(idof,i+1,j) )
              
              i = 3
              rl(idof,i,j) = rl(idof,i,j) + eps_e * buff(i,j) * ( &
                                                  vl(idof,i-2,j) - &
                                            4.0 * vl(idof,i-1,j) + &
                                            6.0 * vl(idof,i,j  ) - &
                                            4.0 * vl(idof,i+1,j) + &
                                            1.0 * vl(idof,i+2,j) )
            end do
          end do
        end if

!.... Interior

        !$omp parallel do private(i, idof)
        do j = 1, ny
          do i = 4, nx-3
            do idof = 1, ndof
              rl(idof,i,j) = rl(idof,i,j) + &
                                  eps_e * buff(i,j) * ( &
                                  fa1 * vl(idof,i-3,j) + &
                                  fa2 * vl(idof,i-2,j) + &
                                  fa3 * vl(idof,i-1,j) + &
                                  fa4 * vl(idof,i,j  ) + &
                                  fa5 * vl(idof,i+1,j) + &
                                  fa6 * vl(idof,i+2,j) + &
                                  fa7 * vl(idof,i+3,j) )
            end do
          end do
        end do
        
        if (rsym) then
          !$omp parallel do private(idof, isign, i)
          do j = 1, ny
            do idof = 1, ndof
              if (idof .eq. 3) then
                isign = -one
              else
                isign = one
              end if
              i = nx-2
              rl(idof,i,j) = rl(idof,i,j) + eps_e * buff(i,j) * ( &
                                            fa1 * vl(idof,nx-5,j) + &
                                            fa2 * vl(idof,nx-4,j) + &
                                            fa3 * vl(idof,nx-3,j) + &
                                            fa4 * vl(idof,nx-2,j) + &
                                            fa5 * vl(idof,nx-1,j) + &
                                            fa6 * vl(idof,nx,j  ) + &
                                            isign * fa7 * vl(idof,nx-1,j) )
              
              i = nx-1
              rl(idof,i,j) = rl(idof,i,j) + eps_e * buff(i,j) * ( &
                                            fa1 * vl(idof,nx-4,j) + &
                                            fa2 * vl(idof,nx-3,j) + &
                                            fa3 * vl(idof,nx-2,j) + &
                                            fa4 * vl(idof,nx-1,j) + &
                                            fa5 * vl(idof,nx,j  ) + &
                                            isign * fa6 * vl(idof,nx-1,j) + &
                                            isign * fa7 * vl(idof,nx-2,j) )

              i = nx
              rl(idof,i,j) = rl(idof,i,j) + eps_e * buff(i,j) * ( &
                                            fa1 * vl(idof,nx-3,j) + &
                                            fa2 * vl(idof,nx-2,j) + &
                                            fa3 * vl(idof,nx-1,j) + &
                                            fa4 * vl(idof,nx,j) + &
                                            isign * fa5 * vl(idof,nx-1,j) + &
                                            isign * fa6 * vl(idof,nx-2,j) + &
                                            isign * fa7 * vl(idof,nx-3,j) )
            end do
          end do
        else
          !$omp parallel do private(idof, i)
          do j = 1, ny
            do idof = 1, ndof
              i = nx-2
              rl(idof,i,j) = rl(idof,i,j) + eps_e * buff(i,j) * ( &
                                                  vl(idof,i-2,j) - &
                                            4.0 * vl(idof,i-1,j) + &
                                            6.0 * vl(idof,i,j  ) - &
                                            4.0 * vl(idof,i+1,j) + &
                                            1.0 * vl(idof,i+2,j) )
  
              i = nx-1
              rl(idof,i,j) = rl(idof,i,j) + eps_e * buff(i,j) * ( &
                                           -1.0 * vl(idof,i-1,j) + &
                                            2.0 * vl(idof,i,j  ) - &
                                            1.0 * vl(idof,i+1,j) )
            end do
          end do
        end if

!===========================================================================
!     S I X T H   O R D E R   T E R M 
!===========================================================================
  
        if (.false.) then
        
          do idof = 1, ndof
            do i = 1, nx

              j = 2
              rl(idof,i,j) = rl(idof,i,j) + eps_e * buff(i,j) * ( &
                                                -1.0 * vl(idof,i,j-1) + &
                                                 2.0 * vl(idof,i,j  ) - &
                                                 1.0 * vl(idof,i,j+1) )
              
              do j = 3, 4
                rl(idof,i,j) = rl(idof,i,j) + eps_e * buff(i,j) * ( &
                                                        vl(idof,i,j-2) - &
                                                  4.0 * vl(idof,i,j-1) + &
                                                  6.0 * vl(idof,i,j  ) - &
                                                  4.0 * vl(idof,i,j+1) + &
                                                  1.0 * vl(idof,i,j+2) )
              end do

              do j = 5, ny-4
                rl(idof,i,j) = rl(idof,i,j) + &
                                    eps_e * buff(i,j) * ( &
                                              fc1 * vl(idof,i,j-4) + &
                                              fc2 * vl(idof,i,j-3) + &
                                              fc3 * vl(idof,i,j-2) + &
                                              fc4 * vl(idof,i,j-1) + &
                                              fc5 * vl(idof,i,j  ) + &
                                              fc6 * vl(idof,i,j+1) + &
                                              fc7 * vl(idof,i,j+2) + &
                                              fc8 * vl(idof,i,j+3) + &
                                              fc9 * vl(idof,i,j+4) )
              end do

              do j = ny-3, ny-2
                rl(idof,i,j) = rl(idof,i,j) + eps_e * buff(i,j) * ( &
                                                        vl(idof,i,j-2) - &
                                                  4.0 * vl(idof,i,j-1) + &
                                                  6.0 * vl(idof,i,j  ) - &
                                                  4.0 * vl(idof,i,j+1) + &
                                                  1.0 * vl(idof,i,j+2) )
              end do

              j = ny-1
              rl(idof,i,j) = rl(idof,i,j) + eps_e * buff(i,j) * ( &
                                            -1.0 * vl(idof,i,j-1) + &
                                             2.0 * vl(idof,i,j  ) - &
                                             1.0 * vl(idof,i,j+1) )

            end do
          end do
                                          
          do idof = 1, ndof
            do j = 1, ny

              i = 2
              rl(idof,i,j) = rl(idof,i,j) + eps_e * buff(i,j) * ( &
                                                -1.0 * vl(idof,i-1,j) + &
                                                 2.0 * vl(idof,i,j  ) - &
                                                 1.0 * vl(idof,i+1,j) )
              
              do i = 3, 4
                rl(idof,i,j) = rl(idof,i,j) + eps_e * buff(i,j) * ( &
                                                        vl(idof,i-2,j) - &
                                                  4.0 * vl(idof,i-1,j) + &
                                                  6.0 * vl(idof,i,j  ) - &
                                                  4.0 * vl(idof,i+1,j) + &
                                                  1.0 * vl(idof,i+2,j) )
              end do

              do i = 5, nx-4
                rl(idof,i,j) = rl(idof,i,j) + &
                                    eps_e * buff(i,j) * ( &
                                              fc1 * vl(idof,i-4,j) + &
                                              fc2 * vl(idof,i-3,j) + &
                                              fc3 * vl(idof,i-2,j) + &
                                              fc4 * vl(idof,i-1,j) + &
                                              fc5 * vl(idof,i,j  ) + &
                                              fc6 * vl(idof,i+1,j) + &
                                              fc7 * vl(idof,i+2,j) + &
                                              fc8 * vl(idof,i+3,j) + &
                                              fc9 * vl(idof,i+4,j) )
              end do

              do i = nx-3, nx-2
                rl(idof,i,j) = rl(idof,i,j) + eps_e * buff(i,j) * ( &
                                                        vl(idof,i-2,j) - &
                                                  4.0 * vl(idof,i-1,j) + &
                                                  6.0 * vl(idof,i,j  ) - &
                                                  4.0 * vl(idof,i+1,j) + &
                                                  1.0 * vl(idof,i+2,j) )
              end do

              i = nx-1
              rl(idof,i,j) = rl(idof,i,j) + eps_e * buff(i,j) * ( &
                                            -1.0 * vl(idof,i-1,j) + &
                                             2.0 * vl(idof,i,j  ) - &
                                             1.0 * vl(idof,i+1,j) )
            end do
          end do

        endif
        
!.... second order dissipation

        if (.false.) then

          do idof = 1, ndof
            do i = 1, nx
            
              j = 2
              rl(idof,i,j) = rl(idof,i,j) + eps_e * buff(i,j) * ( &
                                                -1.0 * vl(idof,i,j-1) + &
                                                 2.0 * vl(idof,i,j  ) - &
                                                 1.0 * vl(idof,i,j+1) )
              
              do j = 3, ny-2
                rl(idof,i,j) = rl(idof,i,j) + eps_e * buff(i,j) * ( &
                                                        vl(idof,i,j-2) - &
                                                  4.0 * vl(idof,i,j-1) + &
                                                  6.0 * vl(idof,i,j  ) - &
                                                  4.0 * vl(idof,i,j+1) + &
                                                  1.0 * vl(idof,i,j+2) )
              end do

              j = ny-1
              rl(idof,i,j) = rl(idof,i,j) + eps_e * buff(i,j) * ( &
                                            -1.0 * vl(idof,i,j-1) + &
                                             2.0 * vl(idof,i,j  ) - &
                                             1.0 * vl(idof,i,j+1) )

            end do
          end do
                                          
          do idof = 1, ndof
            do j = 1, ny

              i = 2
              rl(idof,i,j) = rl(idof,i,j) + eps_e * buff(i,j) * ( &
                                                -1.0 * vl(idof,i-1,j) + &
                                                 2.0 * vl(idof,i,j  ) - &
                                                 1.0 * vl(idof,i+1,j) )
              
              do i = 3, nx-2
                rl(idof,i,j) = rl(idof,i,j) + eps_e * buff(i,j) * ( &
                                                        vl(idof,i-2,j) - &
                                                  4.0 * vl(idof,i-1,j) + &
                                                  6.0 * vl(idof,i,j  ) - &
                                                  4.0 * vl(idof,i+1,j) + &
                                                  1.0 * vl(idof,i+2,j) )
              end do

              i = nx-1
              rl(idof,i,j) = rl(idof,i,j) + eps_e * buff(i,j) * ( &
                                            -1.0 * vl(idof,i-1,j) + &
                                             2.0 * vl(idof,i,j  ) - &
                                             1.0 * vl(idof,i+1,j) )

            end do
          end do
                
        end if
                                
        return
        end
