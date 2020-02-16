!=============================================================================!
        subroutine solve2bc( n, nsys, nblk, mat, rhs, bc )
!=============================================================================!
!  
!  Solves a block tridiagonal system of equations without pivoting. 
!  Vectorized over the first index.
!
!  This routine has been modified to allow for high-order boundary
!  conditions.
!
!  WARNING: Since pivoting is not done, this algorithm could fail
!
!  WARNING: This routine works for nblk = 5 only!
!
!  input:
!    nsys       : number of systems to vectorize over
!    n          : number of rows in the systems
!    nblk       : dimension of the blocks in each row
!
!  inout:
!    mat        : the block tridiagonal matrix to be solved
!    rhs        : the right hand side
!    bc         : the bc array
!=============================================================================!
        implicit none

        integer, intent(in) :: n, nblk, nsys
        real, intent(inout) :: mat(3,nblk,nblk,n,nsys), rhs(nblk,n,nsys)
        real, intent(inout) :: bc(14,nblk,nblk,nsys)

        real :: mult, fact
        
!.... useful constants

        real, parameter :: zero = 0.0, one = 1.0

!.... local variables
        
        integer :: nc, i, l, m, p, q, s, r, t, iv

        !$omp parallel do private(nc,i,l,m,p,q,r,s,t,iv,fact,mult)
        do iv = 1, nsys

!=============================================================================!
!                     F O R W A R D   E L I M I N A T I O N 
!=============================================================================!

!.... first 3 rows are special due to boundary conditions

        i = 1
        
        do nc = 1, nblk - 1

          do m = nc, nblk - 1
            l = m + 1
            mult = -mat(2,l,nc,i,iv) / mat(2,nc,nc,i,iv)
            do q = nc + 1, nblk
              mat(2,l,q,i,iv) = mat(2,l,q,i,iv) + mult * mat(2,nc,q,i,iv)
            end do
            do q = 1, nblk
              mat(3,l,q,i,iv) = mat(3,l,q,i,iv) + mult * mat(3,nc,q,i,iv)
              bc(1,l,q,iv) = bc(1,l,q,iv) + mult * bc(1,nc,q,iv)
              bc(2,l,q,iv) = bc(2,l,q,iv) + mult * bc(2,nc,q,iv)
              bc(3,l,q,iv) = bc(3,l,q,iv) + mult * bc(3,nc,q,iv)
            end do
            rhs(l,i,iv) = rhs(l,i,iv) + mult * rhs(nc,i,iv)
          end do              

        end do  ! end of loop on nc

!.... for all columns
          
        p = i + 1
        r = i + 2
          
        do nc = 1, nblk

          fact = one / mat(2,nc,nc,i,iv)
          
          do m = 1, nblk

!.... second row

            mult = -mat(1,m,nc,p,iv) * fact

            do q = nc + 1, nblk
              mat(1,m,q,p,iv) = mat(1,m,q,p,iv) + mult * mat(2,nc,q,i,iv)
            end do
            do q = 1, nblk
              mat(2,m,q,p,iv) = mat(2,m,q,p,iv) + mult * mat(3,nc,q,i,iv)
              mat(3,m,q,p,iv) = mat(3,m,q,p,iv) + mult * bc(1,nc,q,iv)
              bc(4,m,q,iv)    = bc(4,m,q,iv)    + mult * bc(2,nc,q,iv)
              bc(5,m,q,iv)    = bc(5,m,q,iv)    + mult * bc(3,nc,q,iv)
            end do  
            rhs(m,p,iv) = rhs(m,p,iv) + mult * rhs(nc,i,iv)

!.... third row

            mult = -bc(7,m,nc,iv) * fact
            do q = nc + 1, nblk
              bc(7,m,q,iv) = bc(7,m,q,iv) + mult * mat(2,nc,q,i,iv)
            end do
            do q = 1, nblk
              mat(1,m,q,r,iv) = mat(1,m,q,r,iv) + mult * mat(3,nc,q,i,iv)
              mat(2,m,q,r,iv) = mat(2,m,q,r,iv) + mult * bc(1,nc,q,iv)
              mat(3,m,q,r,iv) = mat(3,m,q,r,iv) + mult * bc(2,nc,q,iv)
              bc(6,m,q,iv)  = bc(6,m,q,iv)  + mult * bc(3,nc,q,iv)
            end do
            rhs(m,r,iv) = rhs(m,r,iv) + mult * rhs(nc,i,iv)
          end do

        end do  ! end of loop on nc

!.... second special row

        i = 2
        
        do nc = 1, nblk - 1

          do m = nc, nblk - 1
            l = m + 1
            mult = -mat(2,l,nc,i,iv) / mat(2,nc,nc,i,iv)
            do q = nc + 1, nblk
              mat(2,l,q,i,iv) = mat(2,l,q,i,iv) + mult * mat(2,nc,q,i,iv)
            end do
            do q = 1, nblk
              mat(3,l,q,i,iv) = mat(3,l,q,i,iv) + mult * mat(3,nc,q,i,iv)
              bc(4,l,q,iv) = bc(4,l,q,iv) + mult * bc(4,nc,q,iv)
              bc(5,l,q,iv) = bc(5,l,q,iv) + mult * bc(5,nc,q,iv)
            end do
            rhs(l,i,iv) = rhs(l,i,iv) + mult * rhs(nc,i,iv)
          end do              

        end do  ! end of loop on nc

!.... for all columns
          
        p  = i + 1
          
        do nc = 1, nblk

          fact = one / mat(2,nc,nc,i,iv)
          do m = 1, nblk
            mult = -mat(1,m,nc,p,iv) * fact
            do q = nc + 1, nblk
              mat(1,m,q,p,iv) = mat(1,m,q,p,iv) + mult * mat(2,nc,q,i,iv)
            end do
            do q = 1, nblk
              mat(2,m,q,p,iv) = mat(2,m,q,p,iv) + mult * mat(3,nc,q,i,iv)
              mat(3,m,q,p,iv) = mat(3,m,q,p,iv) + mult * bc(4,nc,q,iv)
              bc(6,m,q,iv)  = bc(6,m,q,iv)  + mult * bc(5,nc,q,iv)
            end do
            rhs(m,p,iv) = rhs(m,p,iv) + mult * rhs(nc,i,iv)
          end do

        end do  ! end of loop on nc

!.... the third and final special row

        i = 3
        
        do nc = 1, nblk - 1

          do m = nc, nblk - 1
            l = m + 1
            mult = -mat(2,l,nc,i,iv) / mat(2,nc,nc,i,iv)
            do q = nc + 1, nblk
              mat(2,l,q,i,iv) = mat(2,l,q,i,iv) + mult * mat(2,nc,q,i,iv)
            end do
            do q = 1, nblk
              mat(3,l,q,i,iv) = mat(3,l,q,i,iv) + mult * mat(3,nc,q,i,iv)
              bc(6,l,q,iv) = bc(6,l,q,iv) + mult * bc(6,nc,q,iv)
            end do
            rhs(l,i,iv) = rhs(l,i,iv) + mult * rhs(nc,i,iv)
          end do              

        end do  ! end of loop on nc

!.... for all columns
          
        p  = i + 1
          
        do nc = 1, nblk

          fact = one / mat(2,nc,nc,i,iv)
          
          do m = 1, nblk
            mult = -mat(1,m,nc,p,iv) * fact
            do q = nc + 1, nblk
              mat(1,m,q,p,iv) = mat(1,m,q,p,iv) + mult * mat(2,nc,q,i,iv)
            end do
            do q = 1, nblk
              mat(2,m,q,p,iv) = mat(2,m,q,p,iv) + mult * mat(3,nc,q,i,iv)
              mat(3,m,q,p,iv) = mat(3,m,q,p,iv) + mult * bc(6,nc,q,iv)
            end do
            rhs(m,p,iv) = rhs(m,p,iv) + mult * rhs(nc,i,iv)
          end do

        end do  ! end of loop on nc

!.... the rest of the matrix is a standard tridiagonal up to the other boundary

        do i = 4, n - 5

!.... do the first four columns

          do nc = 1, nblk - 1
          
            do m = nc, nblk - 1
              l = m + 1
              mult = -mat(2,l,nc,i,iv) / mat(2,nc,nc,i,iv)
              do q = nc + 1, nblk
                mat(2,l,q,i,iv) = mat(2,l,q,i,iv) + mult * mat(2,nc,q,i,iv)
              end do
              do q = 1, nblk
                mat(3,l,q,i,iv) = mat(3,l,q,i,iv) + mult * mat(3,nc,q,i,iv)
              end do
              rhs(l,i,iv) = rhs(l,i,iv) + mult * rhs(nc,i,iv)
            end do            

          end do        ! end of loop on nc

!.... for all columns
          
          p  = i + 1
            
          do nc = 1, nblk

            fact = one / mat(2,nc,nc,i,iv)
            
            do m = 1, nblk
              mult = -mat(1,m,nc,p,iv) * fact
              do q = nc + 1, nblk
                mat(1,m,q,p,iv) = mat(1,m,q,p,iv) + mult * mat(2,nc,q,i,iv)
              end do
              do q = 1, nblk
                mat(2,m,q,p,iv) = mat(2,m,q,p,iv) + mult * mat(3,nc,q,i,iv)
              end do
              rhs(m,p,iv) = rhs(m,p,iv) + mult * rhs(nc,i,iv)
            end do

          end do        ! end of loop on nc
          
        end do          ! end of loop on i

!.... do the FIFTH to the last row

        i = n - 4
        
!.... first four columns

          do nc = 1, nblk - 1
          
            do m = nc, nblk - 1
              l = m + 1
              mult = -mat(2,l,nc,i,iv) / mat(2,nc,nc,i,iv)
              do q = nc + 1, nblk
                mat(2,l,q,i,iv) = mat(2,l,q,i,iv) + mult * mat(2,nc,q,i,iv)
              end do
              do q = 1, nblk
                mat(3,l,q,i,iv) = mat(3,l,q,i,iv) + mult * mat(3,nc,q,i,iv)
              end do
              rhs(l,i,iv) = rhs(l,i,iv) + mult * rhs(nc,i,iv)
            end do            

          end do        ! end of loop on nc

!.... for all columns
          
          p = i + 1
          r = i + 2
          s = i + 3
          t = i + 4
                    
          do nc = 1, nblk

            fact = one / mat(2,nc,nc,i,iv)
            
            do m = 1, nblk
              mult = -mat(1,m,nc,p,iv) * fact
              do q = nc + 1, nblk
                mat(1,m,q,p,iv) = mat(1,m,q,p,iv) + mult * mat(2,nc,q,i,iv)
              end do
              do q = 1, nblk
                mat(2,m,q,p,iv) = mat(2,m,q,p,iv) + mult * mat(3,nc,q,i,iv)
              end do
              rhs(m,p,iv) = rhs(m,p,iv) + mult * rhs(nc,i,iv)

!.... eliminate the boundary conditions: bc(8)

              mult = -bc(8,m,nc,iv) * fact
              do q = nc + 1, nblk
                bc(8,m,q,iv) = bc(8,m,q,iv) + mult * mat(2,nc,q,i,iv)
              end do
              do q = 1, nblk
                mat(1,m,q,r,iv) = mat(1,m,q,r,iv) + mult * mat(3,nc,q,i,iv)
              end do
              rhs(m,r,iv) = rhs(m,r,iv) + mult * rhs(nc,i,iv)

!.... eliminate the boundary conditions: bc(10)

              mult = -bc(10,m,nc,iv) * fact
              do q = nc + 1, nblk
                bc(10,m,q,iv) = bc(10,m,q,iv) + mult * mat(2,nc,q,i,iv)
              end do
              do q = 1, nblk
                bc(11,m,q,iv) = bc(11,m,q,iv) + mult * mat(3,nc,q,i,iv)
              end do
              rhs(m,s,iv) = rhs(m,s,iv) + mult * rhs(nc,i,iv)
                
!.... eliminate the boundary conditions: bc(12)

              mult = -bc(12,m,nc,iv) * fact
              do q = nc + 1, nblk
                bc(12,m,q,iv) = bc(12,m,q,iv) + mult * mat(2,nc,q,i,iv)
              end do
              do q = 1, nblk
                bc(13,m,q,iv) = bc(13,m,q,iv) + mult * mat(3,nc,q,i,iv)
              end do
              rhs(m,t,iv) = rhs(m,t,iv) + mult * rhs(nc,i,iv)
            end do

          end do        ! end of loop on nc

!.... do the FOURTH to the last row

        i = n - 3
        
!.... first four columns

          do nc = 1, nblk - 1
          
            do m = nc, nblk - 1
              l = m + 1
              mult = -mat(2,l,nc,i,iv) / mat(2,nc,nc,i,iv)
              do q = nc + 1, nblk
                mat(2,l,q,i,iv) = mat(2,l,q,i,iv) + mult * mat(2,nc,q,i,iv)
              end do
              do q = 1, nblk
                mat(3,l,q,i,iv) = mat(3,l,q,i,iv) + mult * mat(3,nc,q,i,iv)
              end do
              rhs(l,i,iv) = rhs(l,i,iv) + mult * rhs(nc,i,iv)
            end do            

          end do        ! end of loop on nc

!.... for all columns
          
          p = i + 1
          r = i + 2
          s = i + 3
                    
          do nc = 1, nblk

            fact = one / mat(2,nc,nc,i,iv)
            
            do m = 1, nblk
              mult = -mat(1,m,nc,p,iv) * fact
              do q = nc + 1, nblk
                mat(1,m,q,p,iv) = mat(1,m,q,p,iv) + mult * mat(2,nc,q,i,iv)
              end do
              do q = 1, nblk
                mat(2,m,q,p,iv) = mat(2,m,q,p,iv) + mult * mat(3,nc,q,i,iv)
              end do
              rhs(m,p,iv) = rhs(m,p,iv) + mult * rhs(nc,i,iv)
              
!.... eliminate the boundary conditions: bc(11)

              mult = -bc(11,m,nc,iv) * fact
              do q = nc + 1, nblk
                bc(11,m,q,iv) = bc(11,m,q,iv) + mult * mat(2,nc,q,i,iv)
              end do
              do q = 1, nblk
                mat(1,m,q,r,iv) = mat(1,m,q,r,iv) + mult * mat(3,nc,q,i,iv)
              end do
              rhs(m,r,iv) = rhs(m,r,iv) + mult * rhs(nc,i,iv)
              
!.... eliminate the boundary conditions: bc(13)

              mult = -bc(13,m,nc,iv) * fact
              do q = nc + 1, nblk
                bc(13,m,q,iv) = bc(13,m,q,iv) + mult * mat(2,nc,q,i,iv)
              end do  
              do q = 1, nblk
                bc(14,m,q,iv) = bc(14,m,q,iv) + mult * mat(3,nc,q,i,iv)
              end do
              rhs(m,s,iv) = rhs(m,s,iv) + mult * rhs(nc,i,iv)
            end do

          end do        ! end of loop on nc

!.... do the THIRD to the last row

        i = n - 2
        
!.... first four columns

          do nc = 1, nblk - 1
          
            do m = nc, nblk - 1
              l = m + 1
              mult = -mat(2,l,nc,i,iv) / mat(2,nc,nc,i,iv)
              do q = nc + 1, nblk
                mat(2,l,q,i,iv) = mat(2,l,q,i,iv) + mult * mat(2,nc,q,i,iv)
              end do
              do q = 1, nblk
                mat(3,l,q,i,iv) = mat(3,l,q,i,iv) + mult * mat(3,nc,q,i,iv)
                bc(9,l,q,iv) = bc(9,l,q,iv) + mult * bc(9,nc,q,iv)
              end do
              rhs(l,i,iv) = rhs(l,i,iv) + mult * rhs(nc,i,iv)
            end do            

          end do        ! end of loop on nc

!.... for all columns
          
          p = i + 1
          r = i + 2
            
          do nc = 1, nblk

            fact = one / mat(2,nc,nc,i,iv)
            
            do m = 1, nblk
              mult = -mat(1,m,nc,p,iv) * fact
              do q = nc + 1, nblk
                mat(1,m,q,p,iv) = mat(1,m,q,p,iv) + mult * mat(2,nc,q,i,iv)
              end do
              do q = 1, nblk
                mat(2,m,q,p,iv) = mat(2,m,q,p,iv) + mult * mat(3,nc,q,i,iv)
                mat(3,m,q,p,iv) = mat(3,m,q,p,iv) + mult * bc(9,nc,q,iv)
              end do
              rhs(m,p,iv) = rhs(m,p,iv) + mult * rhs(nc,i,iv)
              
!.... eliminate the boundary condition:  bc(14)

              mult = -bc(14,m,nc,iv) * fact
              do q = nc + 1, nblk
                bc(14,m,q,iv) = bc(14,m,q,iv) + mult * mat(2,nc,q,i,iv)
              end do
              do q = 1, nblk
                mat(1,m,q,r,iv) = mat(1,m,q,r,iv) + mult * mat(3,nc,q,i,iv)
                mat(2,m,q,r,iv) = mat(2,m,q,r,iv) + mult * bc(9,nc,q,iv)
              end do
              rhs(m,r,iv) = rhs(m,r,iv) + mult * rhs(nc,i,iv)
            end do

          end do        ! end of loop on nc

!.... do the SECOND to the last row (which is normal tridiagonal)

        i = n - 1
        
!.... first four columns

          do nc = 1, nblk - 1
          
            do m = nc, nblk - 1
              l = m + 1
              mult = -mat(2,l,nc,i,iv) / mat(2,nc,nc,i,iv)
              do q = nc + 1, nblk
                mat(2,l,q,i,iv) = mat(2,l,q,i,iv) + mult * mat(2,nc,q,i,iv)
              end do
              do q = 1, nblk
                mat(3,l,q,i,iv) = mat(3,l,q,i,iv) + mult * mat(3,nc,q,i,iv)
              end do
              rhs(l,i,iv) = rhs(l,i,iv) + mult * rhs(nc,i,iv)
            end do            

          end do        ! end of loop on nc

!.... for all columns
          
          p = i + 1
            
          do nc = 1, nblk

            fact = one / mat(2,nc,nc,i,iv)
            
            do m = 1, nblk
              mult = -mat(1,m,nc,p,iv) * fact
              do q = nc + 1, nblk
                mat(1,m,q,p,iv) = mat(1,m,q,p,iv) + mult * mat(2,nc,q,i,iv)
              end do
              do q = 1, nblk
                mat(2,m,q,p,iv) = mat(2,m,q,p,iv) + mult * mat(3,nc,q,i,iv)
              end do
              rhs(m,p,iv) = rhs(m,p,iv) + mult * rhs(nc,i,iv)
            end do

          end do        ! end of loop on nc

!.... do the last row
          
        i = n
        
        do nc = 1, nblk - 1
        
          do m = nc, nblk - 1
            l = m + 1
            mult = -mat(2,l,nc,i,iv) / mat(2,nc,nc,i,iv)
            do q = nc + 1, nblk
              mat(2,l,q,i,iv) = mat(2,l,q,i,iv) + mult * mat(2,nc,q,i,iv)
            end do
            rhs(l,i,iv) = rhs(l,i,iv) + mult * rhs(nc,i,iv)
          end do              

        end do  ! end of loop on nc

!=============================================================================!
!                   B A C K W A R D   S U B S T I T U T I O N 
!=============================================================================!

!.... do the LAST row first

        i = n
        
        do m = nblk, 1, -1
          do q = m+1, nblk
            rhs(m,i,iv) = rhs(m,i,iv) - rhs(q,i,iv) * mat(2,m,q,i,iv)
          end do
          rhs(m,i,iv) = rhs(m,i,iv) / mat(2,m,m,i,iv)
        end do
        
!.... do the SECOND to the last row (normal row)

        i = n - 1
        p = i + 1
        
        do m = nblk, 1, -1
          do q = m+1, nblk
            rhs(m,i,iv) = rhs(m,i,iv) - rhs(q,i,iv) * mat(2,m,q,i,iv)
          end do
          do q = 1, nblk
            rhs(m,i,iv) = rhs(m,i,iv) - rhs(q,p,iv) * mat(3,m,q,i,iv)
          end do
          rhs(m,i,iv) = rhs(m,i,iv) / mat(2,m,m,i,iv)
        end do

!.... do the THIRD to the last row

        i = n - 2
        p = i + 1
        r = i + 2
        
        do m = nblk, 1, -1
          do q = m+1, nblk
            rhs(m,i,iv) = rhs(m,i,iv) - rhs(q,i,iv) * mat(2,m,q,i,iv)
          end do
          do q = 1, nblk
            rhs(m,i,iv) = rhs(m,i,iv) - rhs(q,p,iv) * mat(3,m,q,i,iv)
            rhs(m,i,iv) = rhs(m,i,iv) - rhs(q,r,iv) * bc(9,m,q,iv)
          end do
          rhs(m,i,iv) = rhs(m,i,iv) / mat(2,m,m,i,iv)
        end do

!.... now do the rest of the rows in reverse order upto the boundary

        do i = n - 3, 4, -1
        
          p = i + 1
          
          do m = nblk, 1, -1
            do q = m+1, nblk
              rhs(m,i,iv) = rhs(m,i,iv) - rhs(q,i,iv) * mat(2,m,q,i,iv)
            end do
            do q = 1, nblk
              rhs(m,i,iv) = rhs(m,i,iv) - rhs(q,p,iv) * mat(3,m,q,i,iv)
            end do
            rhs(m,i,iv) = rhs(m,i,iv) / mat(2,m,m,i,iv)
          end do
                                  
        end do
        
!.... now do the boundary rows

        i = 3
        
        p = i + 1
        r = i + 2
        
        do m = nblk, 1, -1
          do q = m+1, nblk
            rhs(m,i,iv) = rhs(m,i,iv) - rhs(q,i,iv) * mat(2,m,q,i,iv)
          end do
          do q = 1, nblk
            rhs(m,i,iv) = rhs(m,i,iv) - rhs(q,p,iv) * mat(3,m,q,i,iv)
            rhs(m,i,iv) = rhs(m,i,iv) - rhs(q,r,iv) * bc(6,m,q,iv)
          end do
          rhs(m,i,iv) = rhs(m,i,iv) / mat(2,m,m,i,iv)
        end do

        i = 2
        
        p = i + 1
        r = i + 2
        s = i + 3
        
        do m = nblk, 1, -1
          do q = m+1, nblk
            rhs(m,i,iv) = rhs(m,i,iv) - rhs(q,i,iv) * mat(2,m,q,i,iv)
          end do
          do q = 1, nblk
            rhs(m,i,iv) = rhs(m,i,iv) - rhs(q,p,iv) * mat(3,m,q,i,iv)
            rhs(m,i,iv) = rhs(m,i,iv) - rhs(q,r,iv) * bc(4,m,q,iv)
            rhs(m,i,iv) = rhs(m,i,iv) - rhs(q,s,iv) * bc(5,m,q,iv)
          end do
          rhs(m,i,iv) = rhs(m,i,iv) / mat(2,m,m,i,iv)
        end do

        i = 1
        
        p = i + 1
        r = i + 2
        s = i + 3
        t = i + 4
        
        do m = nblk, 1, -1
          do q = m+1, nblk
            rhs(m,i,iv) = rhs(m,i,iv) - rhs(q,i,iv) * mat(2,m,q,i,iv)
          end do
          do q = 1, nblk
            rhs(m,i,iv) = rhs(m,i,iv) - rhs(q,p,iv) * mat(3,m,q,i,iv)
            rhs(m,i,iv) = rhs(m,i,iv) - rhs(q,r,iv) * bc(1,m,q,iv)
            rhs(m,i,iv) = rhs(m,i,iv) - rhs(q,s,iv) * bc(2,m,q,iv)
            rhs(m,i,iv) = rhs(m,i,iv) - rhs(q,t,iv) * bc(3,m,q,iv)
          end do
          rhs(m,i,iv) = rhs(m,i,iv) / mat(2,m,m,i,iv)
        end do

        end do
!=============================================================================!
        return
        end
