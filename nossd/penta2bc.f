!=============================================================================!
        subroutine penta2bc( n, nsys, nblk, mat, rhs, mult, fact )
!=============================================================================!
!  
!  Solves a block pentadiagonal system of equations without pivoting.
!  Vectorized over the second index.
!
!  This version includes the Boundary Treatment
!
!  WARNING: Since pivoting is not done, this algorithm could fail
!
!  input:
!    nsys       : number of systems to vectorize over
!    n          : number of rows in the systems
!    nblk       : dimension of the blocks in each row
!
!  inout:
!    mat        : the block pentadiagonal matrix to be solved
!    rhs        : the right hand side
!
!  Revised:  7-17-95
!=============================================================================!
        implicit none

        integer n, nblk, nsys
        real mat(n,nsys,nblk,nblk,5), rhs(n,nsys,nblk)
        real mult(nsys), fact(nsys)
        
!.... useful constants

        real zero, one
        parameter (zero = 0.0, one = 1.0)

!.... local variables
        
        integer nc, i, j, l, m, p, q, r, s, t, iv
!=============================================================================!
!                     F O R W A R D   E L I M I N A T I O N 
!=============================================================================!

!.... do the first row

        i = 1
        
!.... do the first four columns

          do nc = 1, nblk - 1
          
            do m = nc, nblk - 1
              l = m + 1
              do iv = 1, nsys
                mult(iv) = -mat(i,iv,l,nc,3) / mat(i,iv,nc,nc,3)
              end do
              do q = nc + 1, nblk
                do iv = 1, nsys
                  mat(i,iv,l,q,3) = mat(i,iv,l,q,3) + mult(iv) * mat(i,iv,nc,q,3)
                end do
              end do
              do q = 1, nblk
                do iv = 1, nsys
                  mat(i,iv,l,q,4) = mat(i,iv,l,q,4) + mult(iv) * mat(i,iv,nc,q,4)
                  mat(i,iv,l,q,5) = mat(i,iv,l,q,5) + mult(iv) * mat(i,iv,nc,q,5)
                  mat(i,iv,l,q,1) = mat(i,iv,l,q,1) + mult(iv) * mat(i,iv,nc,q,1)
                  mat(i,iv,l,q,2) = mat(i,iv,l,q,2) + mult(iv) * mat(i,iv,nc,q,2)
                end do
              end do
              do iv = 1, nsys
                rhs(i,iv,l) = rhs(i,iv,l) + mult(iv) * rhs(i,iv,nc)
              end do
            end do            

          end do        ! end of loop on nc

!.... for all columns
          
          do nc = 1, nblk

            p = i + 1
            r = i + 2
            
            do iv = 1, nsys
              fact(iv) = one / mat(i,iv,nc,nc,3)
            end do
            
            do m = 1, nblk
              do iv = 1, nsys
                mult(iv) = -mat(p,iv,m,nc,2) * fact(iv)
              end do
              do q = nc + 1, nblk
                do iv = 1, nsys
                  mat(p,iv,m,q,2) = mat(p,iv,m,q,2) + mult(iv) * mat(i,iv,nc,q,3)
                end do
              end do
              do q = 1, nblk
                do iv = 1, nsys
                  mat(p,iv,m,q,3) = mat(p,iv,m,q,3) + mult(iv) * mat(i,iv,nc,q,4)
                  mat(p,iv,m,q,4) = mat(p,iv,m,q,4) + mult(iv) * mat(i,iv,nc,q,5)
                  mat(p,iv,m,q,5) = mat(p,iv,m,q,5) + mult(iv) * mat(i,iv,nc,q,1)
                  mat(p,iv,m,q,1) = mat(p,iv,m,q,1) + mult(iv) * mat(i,iv,nc,q,2)
                end do
              end do
              do iv = 1, nsys
                rhs(p,iv,m) = rhs(p,iv,m) + mult(iv) * rhs(i,iv,nc)
              end do
              
              do iv = 1, nsys
                mult(iv) = -mat(r,iv,m,nc,1) * fact(iv)
              end do
              do q = nc + 1, nblk
                do iv = 1, nsys
                  mat(r,iv,m,q,1) = mat(r,iv,m,q,1) + mult(iv) * mat(i,iv,nc,q,3)
                end do
              end do
              do q = 1, nblk
                do iv = 1, nsys
                  mat(r,iv,m,q,2) = mat(r,iv,m,q,2) + mult(iv) * mat(i,iv,nc,q,4)
                  mat(r,iv,m,q,3) = mat(r,iv,m,q,3) + mult(iv) * mat(i,iv,nc,q,5)
                  mat(r,iv,m,q,4) = mat(r,iv,m,q,4) + mult(iv) * mat(i,iv,nc,q,1)
                  mat(r,iv,m,q,5) = mat(r,iv,m,q,5) + mult(iv) * mat(i,iv,nc,q,2)
                end do
              end do
              do iv = 1, nsys
                rhs(r,iv,m) = rhs(r,iv,m) + mult(iv) * rhs(i,iv,nc)
              end do
              
            end do

          end do        ! end of loop on nc

!.... do the second row

        i = 2
        
!.... do the first four columns

          do nc = 1, nblk - 1
          
            do m = nc, nblk - 1
              l = m + 1
              do iv = 1, nsys
                mult(iv) = -mat(i,iv,l,nc,3) / mat(i,iv,nc,nc,3)
              end do
              do q = nc + 1, nblk
                do iv = 1, nsys
                  mat(i,iv,l,q,3) = mat(i,iv,l,q,3) + mult(iv) * mat(i,iv,nc,q,3)
                end do
              end do
              do q = 1, nblk
                do iv = 1, nsys
                  mat(i,iv,l,q,4) = mat(i,iv,l,q,4) + mult(iv) * mat(i,iv,nc,q,4)
                  mat(i,iv,l,q,5) = mat(i,iv,l,q,5) + mult(iv) * mat(i,iv,nc,q,5)
                  mat(i,iv,l,q,1) = mat(i,iv,l,q,1) + mult(iv) * mat(i,iv,nc,q,1)
                end do
              end do
              do iv = 1, nsys
                rhs(i,iv,l) = rhs(i,iv,l) + mult(iv) * rhs(i,iv,nc)
              end do
            end do            

          end do        ! end of loop on nc

!.... for all columns
          
          do nc = 1, nblk

            p = i + 1
            r = i + 2
            
            do iv = 1, nsys
              fact(iv) = one / mat(i,iv,nc,nc,3)
            end do
            
            do m = 1, nblk
              do iv = 1, nsys
                mult(iv) = -mat(p,iv,m,nc,2) * fact(iv)
              end do
              do q = nc + 1, nblk
                do iv = 1, nsys
                  mat(p,iv,m,q,2) = mat(p,iv,m,q,2) + mult(iv) * mat(i,iv,nc,q,3)
                end do
              end do
              do q = 1, nblk
                do iv = 1, nsys
                  mat(p,iv,m,q,3) = mat(p,iv,m,q,3) + mult(iv) * mat(i,iv,nc,q,4)
                  mat(p,iv,m,q,4) = mat(p,iv,m,q,4) + mult(iv) * mat(i,iv,nc,q,5)
                  mat(p,iv,m,q,5) = mat(p,iv,m,q,5) + mult(iv) * mat(i,iv,nc,q,1)
                end do
              end do
              do iv = 1, nsys
                rhs(p,iv,m) = rhs(p,iv,m) + mult(iv) * rhs(i,iv,nc)
              end do
              
              do iv = 1, nsys
                mult(iv) = -mat(r,iv,m,nc,1) * fact(iv)
              end do
              do q = nc + 1, nblk
                do iv = 1, nsys
                  mat(r,iv,m,q,1) = mat(r,iv,m,q,1) + mult(iv) * mat(i,iv,nc,q,3)
                end do
              end do
              do q = 1, nblk
                do iv = 1, nsys
                  mat(r,iv,m,q,2) = mat(r,iv,m,q,2) + mult(iv) * mat(i,iv,nc,q,4)
                  mat(r,iv,m,q,3) = mat(r,iv,m,q,3) + mult(iv) * mat(i,iv,nc,q,5)
                  mat(r,iv,m,q,4) = mat(r,iv,m,q,4) + mult(iv) * mat(i,iv,nc,q,1)
                end do
              end do
              do iv = 1, nsys
                rhs(r,iv,m) = rhs(r,iv,m) + mult(iv) * rhs(i,iv,nc)
              end do
              
            end do

          end do        ! end of loop on nc

!.... do the interior rows

        do i = 3, n - 5

!.... do the first four columns

          do nc = 1, nblk - 1
          
            do m = nc, nblk - 1
              l = m + 1
              do iv = 1, nsys
                mult(iv) = -mat(i,iv,l,nc,3) / mat(i,iv,nc,nc,3)
              end do
              do q = nc + 1, nblk
                do iv = 1, nsys
                  mat(i,iv,l,q,3) = mat(i,iv,l,q,3) + mult(iv) * mat(i,iv,nc,q,3)
                end do
              end do
              do q = 1, nblk
                do iv = 1, nsys
                  mat(i,iv,l,q,4) = mat(i,iv,l,q,4) + mult(iv) * mat(i,iv,nc,q,4)
                  mat(i,iv,l,q,5) = mat(i,iv,l,q,5) + mult(iv) * mat(i,iv,nc,q,5)
                end do
              end do
              do iv = 1, nsys
                rhs(i,iv,l) = rhs(i,iv,l) + mult(iv) * rhs(i,iv,nc)
              end do
            end do            

          end do        ! end of loop on nc

!.... for all columns
          
          do nc = 1, nblk

            p = i + 1
            r = i + 2
            
            do iv = 1, nsys
              fact(iv) = one / mat(i,iv,nc,nc,3)
            end do
            
            do m = 1, nblk
              do iv = 1, nsys
                mult(iv) = -mat(p,iv,m,nc,2) * fact(iv)
              end do
              do q = nc + 1, nblk
                do iv = 1, nsys
                  mat(p,iv,m,q,2) = mat(p,iv,m,q,2) + mult(iv) * mat(i,iv,nc,q,3)
                end do
              end do
              do q = 1, nblk
                do iv = 1, nsys
                  mat(p,iv,m,q,3) = mat(p,iv,m,q,3) + mult(iv) * mat(i,iv,nc,q,4)
                  mat(p,iv,m,q,4) = mat(p,iv,m,q,4) + mult(iv) * mat(i,iv,nc,q,5)
                end do
              end do
              do iv = 1, nsys
                rhs(p,iv,m) = rhs(p,iv,m) + mult(iv) * rhs(i,iv,nc)
              end do
              
              do iv = 1, nsys
                mult(iv) = -mat(r,iv,m,nc,1) * fact(iv)
              end do
              do q = nc + 1, nblk
                do iv = 1, nsys
                  mat(r,iv,m,q,1) = mat(r,iv,m,q,1) + mult(iv) * mat(i,iv,nc,q,3)
                end do
              end do
              do q = 1, nblk
                do iv = 1, nsys
                  mat(r,iv,m,q,2) = mat(r,iv,m,q,2) + mult(iv) * mat(i,iv,nc,q,4)
                  mat(r,iv,m,q,3) = mat(r,iv,m,q,3) + mult(iv) * mat(i,iv,nc,q,5)
                end do
              end do
              do iv = 1, nsys
                rhs(r,iv,m) = rhs(r,iv,m) + mult(iv) * rhs(i,iv,nc)
              end do
              
            end do

          end do        ! end of loop on nc
          
        end do          ! end of loop on i

!.... do the (n-4)th row

        i = n - 4

!.... do the first four columns

          do nc = 1, nblk - 1
          
            do m = nc, nblk - 1
              l = m + 1
              do iv = 1, nsys
                mult(iv) = -mat(i,iv,l,nc,3) / mat(i,iv,nc,nc,3)
              end do
              do q = nc + 1, nblk
                do iv = 1, nsys
                  mat(i,iv,l,q,3) = mat(i,iv,l,q,3) + mult(iv) * mat(i,iv,nc,q,3)
                end do
              end do
              do q = 1, nblk
                do iv = 1, nsys
                  mat(i,iv,l,q,4) = mat(i,iv,l,q,4) + mult(iv) * mat(i,iv,nc,q,4)
                  mat(i,iv,l,q,5) = mat(i,iv,l,q,5) + mult(iv) * mat(i,iv,nc,q,5)
                end do
              end do
              do iv = 1, nsys
                rhs(i,iv,l) = rhs(i,iv,l) + mult(iv) * rhs(i,iv,nc)
              end do
            end do            

          end do        ! end of loop on nc

!.... for all columns
          
          do nc = 1, nblk

            p = i + 1
            r = i + 2
            s = i + 3
            t = i + 4
            
            do iv = 1, nsys
              fact(iv) = one / mat(i,iv,nc,nc,3)
            end do
            
            do m = 1, nblk
              do iv = 1, nsys
                mult(iv) = -mat(p,iv,m,nc,2) * fact(iv)
              end do
              do q = nc + 1, nblk
                do iv = 1, nsys
                  mat(p,iv,m,q,2) = mat(p,iv,m,q,2) + mult(iv) * mat(i,iv,nc,q,3)
                end do
              end do
              do q = 1, nblk
                do iv = 1, nsys
                  mat(p,iv,m,q,3) = mat(p,iv,m,q,3) + mult(iv) * mat(i,iv,nc,q,4)
                  mat(p,iv,m,q,4) = mat(p,iv,m,q,4) + mult(iv) * mat(i,iv,nc,q,5)
                end do
              end do
              do iv = 1, nsys
                rhs(p,iv,m) = rhs(p,iv,m) + mult(iv) * rhs(i,iv,nc)
              end do
              
              do iv = 1, nsys
                mult(iv) = -mat(r,iv,m,nc,1) * fact(iv)
              end do
              do q = nc + 1, nblk
                do iv = 1, nsys
                  mat(r,iv,m,q,1) = mat(r,iv,m,q,1) + mult(iv) * mat(i,iv,nc,q,3)
                end do
              end do
              do q = 1, nblk
                do iv = 1, nsys
                  mat(r,iv,m,q,2) = mat(r,iv,m,q,2) + mult(iv) * mat(i,iv,nc,q,4)
                  mat(r,iv,m,q,3) = mat(r,iv,m,q,3) + mult(iv) * mat(i,iv,nc,q,5)
                end do
              end do
              do iv = 1, nsys
                rhs(r,iv,m) = rhs(r,iv,m) + mult(iv) * rhs(i,iv,nc)
              end do
              
              do iv = 1, nsys
                mult(iv) = -mat(s,iv,m,nc,5) * fact(iv)
              end do
              do q = nc + 1, nblk
                do iv = 1, nsys
                  mat(s,iv,m,q,5) = mat(s,iv,m,q,5) + mult(iv) * mat(i,iv,nc,q,3)
                end do
              end do
              do q = 1, nblk
                do iv = 1, nsys
                  mat(s,iv,m,q,1) = mat(s,iv,m,q,1) + mult(iv) * mat(i,iv,nc,q,4)
                  mat(s,iv,m,q,2) = mat(s,iv,m,q,2) + mult(iv) * mat(i,iv,nc,q,5)
                end do
              end do
              do iv = 1, nsys
                rhs(s,iv,m) = rhs(s,iv,m) + mult(iv) * rhs(i,iv,nc)
              end do

              do iv = 1, nsys
                mult(iv) = -mat(t,iv,m,nc,4) * fact(iv)
              end do
              do q = nc + 1, nblk
                do iv = 1, nsys
                  mat(t,iv,m,q,4) = mat(t,iv,m,q,4) + mult(iv) * mat(i,iv,nc,q,3)
                end do
              end do
              do q = 1, nblk
                do iv = 1, nsys
                  mat(t,iv,m,q,5) = mat(t,iv,m,q,5) + mult(iv) * mat(i,iv,nc,q,4)
                  mat(t,iv,m,q,1) = mat(t,iv,m,q,1) + mult(iv) * mat(i,iv,nc,q,5)
                end do
              end do
              do iv = 1, nsys
                rhs(t,iv,m) = rhs(t,iv,m) + mult(iv) * rhs(i,iv,nc)
              end do

            end do

          end do        ! end of loop on nc

!.... do the (n-3)rd row

        i = n - 3

!.... do the first four columns

          do nc = 1, nblk - 1
          
            do m = nc, nblk - 1
              l = m + 1
              do iv = 1, nsys
                mult(iv) = -mat(i,iv,l,nc,3) / mat(i,iv,nc,nc,3)
              end do
              do q = nc + 1, nblk
                do iv = 1, nsys
                  mat(i,iv,l,q,3) = mat(i,iv,l,q,3) + mult(iv) * mat(i,iv,nc,q,3)
                end do
              end do
              do q = 1, nblk
                do iv = 1, nsys
                  mat(i,iv,l,q,4) = mat(i,iv,l,q,4) + mult(iv) * mat(i,iv,nc,q,4)
                  mat(i,iv,l,q,5) = mat(i,iv,l,q,5) + mult(iv) * mat(i,iv,nc,q,5)
                end do
              end do
              do iv = 1, nsys
                rhs(i,iv,l) = rhs(i,iv,l) + mult(iv) * rhs(i,iv,nc)
              end do
            end do            

          end do        ! end of loop on nc

!.... for all columns
          
          do nc = 1, nblk

            p = i + 1
            r = i + 2
            s = i + 3
            
            do iv = 1, nsys
              fact(iv) = one / mat(i,iv,nc,nc,3)
            end do
            
            do m = 1, nblk
              do iv = 1, nsys
                mult(iv) = -mat(p,iv,m,nc,2) * fact(iv)
              end do
              do q = nc + 1, nblk
                do iv = 1, nsys
                  mat(p,iv,m,q,2) = mat(p,iv,m,q,2) + mult(iv) * mat(i,iv,nc,q,3)
                end do
              end do
              do q = 1, nblk
                do iv = 1, nsys
                  mat(p,iv,m,q,3) = mat(p,iv,m,q,3) + mult(iv) * mat(i,iv,nc,q,4)
                  mat(p,iv,m,q,4) = mat(p,iv,m,q,4) + mult(iv) * mat(i,iv,nc,q,5)
                end do
              end do
              do iv = 1, nsys
                rhs(p,iv,m) = rhs(p,iv,m) + mult(iv) * rhs(i,iv,nc)
              end do
              
              do iv = 1, nsys
                mult(iv) = -mat(r,iv,m,nc,1) * fact(iv)
              end do
              do q = nc + 1, nblk
                do iv = 1, nsys
                  mat(r,iv,m,q,1) = mat(r,iv,m,q,1) + mult(iv) * mat(i,iv,nc,q,3)
                end do
              end do
              do q = 1, nblk
                do iv = 1, nsys
                  mat(r,iv,m,q,2) = mat(r,iv,m,q,2) + mult(iv) * mat(i,iv,nc,q,4)
                  mat(r,iv,m,q,3) = mat(r,iv,m,q,3) + mult(iv) * mat(i,iv,nc,q,5)
                end do
              end do
              do iv = 1, nsys
                rhs(r,iv,m) = rhs(r,iv,m) + mult(iv) * rhs(i,iv,nc)
              end do
              
              do iv = 1, nsys
                mult(iv) = -mat(s,iv,m,nc,5) * fact(iv)
              end do
              do q = nc + 1, nblk
                do iv = 1, nsys
                  mat(s,iv,m,q,5) = mat(s,iv,m,q,5) + mult(iv) * mat(i,iv,nc,q,3)
                end do
              end do
              do q = 1, nblk
                do iv = 1, nsys
                  mat(s,iv,m,q,1) = mat(s,iv,m,q,1) + mult(iv) * mat(i,iv,nc,q,4)
                  mat(s,iv,m,q,2) = mat(s,iv,m,q,2) + mult(iv) * mat(i,iv,nc,q,5)
                end do
              end do
              do iv = 1, nsys
                rhs(s,iv,m) = rhs(s,iv,m) + mult(iv) * rhs(i,iv,nc)
              end do

            end do

          end do        ! end of loop on nc

!.... do the (n-2)nd row

        i = n - 2

!.... do the first four columns

          do nc = 1, nblk - 1
          
            do m = nc, nblk - 1
              l = m + 1
              do iv = 1, nsys
                mult(iv) = -mat(i,iv,l,nc,3) / mat(i,iv,nc,nc,3)
              end do
              do q = nc + 1, nblk
                do iv = 1, nsys
                  mat(i,iv,l,q,3) = mat(i,iv,l,q,3) + mult(iv) * mat(i,iv,nc,q,3)
                end do
              end do
              do q = 1, nblk
                do iv = 1, nsys
                  mat(i,iv,l,q,4) = mat(i,iv,l,q,4) + mult(iv) * mat(i,iv,nc,q,4)
                  mat(i,iv,l,q,5) = mat(i,iv,l,q,5) + mult(iv) * mat(i,iv,nc,q,5)
                end do
              end do
              do iv = 1, nsys
                rhs(i,iv,l) = rhs(i,iv,l) + mult(iv) * rhs(i,iv,nc)
              end do
            end do            

          end do        ! end of loop on nc

!.... for all columns
          
          do nc = 1, nblk

            p = i + 1
            r = i + 2
            
            do iv = 1, nsys
              fact(iv) = one / mat(i,iv,nc,nc,3)
            end do
            
            do m = 1, nblk
              do iv = 1, nsys
                mult(iv) = -mat(p,iv,m,nc,2) * fact(iv)
              end do
              do q = nc + 1, nblk
                do iv = 1, nsys
                  mat(p,iv,m,q,2) = mat(p,iv,m,q,2) + mult(iv) * mat(i,iv,nc,q,3)
                end do
              end do
              do q = 1, nblk
                do iv = 1, nsys
                  mat(p,iv,m,q,3) = mat(p,iv,m,q,3) + mult(iv) * mat(i,iv,nc,q,4)
                  mat(p,iv,m,q,4) = mat(p,iv,m,q,4) + mult(iv) * mat(i,iv,nc,q,5)
                end do
              end do
              do iv = 1, nsys
                rhs(p,iv,m) = rhs(p,iv,m) + mult(iv) * rhs(i,iv,nc)
              end do
              
              do iv = 1, nsys
                mult(iv) = -mat(r,iv,m,nc,1) * fact(iv)
              end do
              do q = nc + 1, nblk
                do iv = 1, nsys
                  mat(r,iv,m,q,1) = mat(r,iv,m,q,1) + mult(iv) * mat(i,iv,nc,q,3)
                end do
              end do
              do q = 1, nblk
                do iv = 1, nsys
                  mat(r,iv,m,q,2) = mat(r,iv,m,q,2) + mult(iv) * mat(i,iv,nc,q,4)
                  mat(r,iv,m,q,3) = mat(r,iv,m,q,3) + mult(iv) * mat(i,iv,nc,q,5)
                end do
              end do
              do iv = 1, nsys
                rhs(r,iv,m) = rhs(r,iv,m) + mult(iv) * rhs(i,iv,nc)
              end do
              
            end do

          end do        ! end of loop on nc

!.... do the next-to-the-last row
        
        i = n - 1

!.... do the first four columns

          do nc = 1, nblk - 1
          
            do m = nc, nblk - 1
              l = m + 1
              do iv = 1, nsys
                mult(iv) = -mat(i,iv,l,nc,3) / mat(i,iv,nc,nc,3)
              end do
              do q = nc + 1, nblk
                do iv = 1, nsys
                  mat(i,iv,l,q,3) = mat(i,iv,l,q,3) + mult(iv) * mat(i,iv,nc,q,3)
                end do
              end do
              do q = 1, nblk
                do iv = 1, nsys
                  mat(i,iv,l,q,4) = mat(i,iv,l,q,4) + mult(iv) * mat(i,iv,nc,q,4)
                end do
              end do
              do iv = 1, nsys
                rhs(i,iv,l) = rhs(i,iv,l) + mult(iv) * rhs(i,iv,nc)
              end do
            end do            

          end do        ! end of loop on nc

!.... for all columns
          
          do nc = 1, nblk

            p = i + 1
            
            do iv = 1, nsys
              fact(iv) = one / mat(i,iv,nc,nc,3)
            end do
            
            do m = 1, nblk
              do iv = 1, nsys
                mult(iv) = -mat(p,iv,m,nc,2) * fact(iv)
              end do
              do q = nc + 1, nblk
                do iv = 1, nsys
                  mat(p,iv,m,q,2) = mat(p,iv,m,q,2) + mult(iv) * mat(i,iv,nc,q,3)
                end do
              end do
              do q = 1, nblk
                do iv = 1, nsys
                  mat(p,iv,m,q,3) = mat(p,iv,m,q,3) + mult(iv) * mat(i,iv,nc,q,4)
                end do
              end do
              do iv = 1, nsys
                rhs(p,iv,m) = rhs(p,iv,m) + mult(iv) * rhs(i,iv,nc)
              end do
            end do

          end do        ! end of loop on nc

!.... do the last row

        i = n
        
        do nc = 1, nblk - 1
        
          do m = nc, nblk - 1
            l = m + 1
            do iv = 1, nsys
              mult(iv) = -mat(i,iv,l,nc,3) / mat(i,iv,nc,nc,3)
            end do
            do q = nc + 1, nblk
              do iv = 1, nsys
                mat(i,iv,l,q,3) = mat(i,iv,l,q,3) + mult(iv) * mat(i,iv,nc,q,3)
              end do
            end do
            do iv = 1, nsys
              rhs(i,iv,l) = rhs(i,iv,l) + mult(iv) * rhs(i,iv,nc)
            end do
          end do              

        end do  ! end of loop on nc

!=============================================================================!
!                   B A C K W A R D   S U B S T I T U T I O N 
!=============================================================================!

!.... do the last row first

        i = n
        
        do m = nblk, 1, -1
          do q = m+1, nblk
            do iv = 1, nsys
              rhs(i,iv,m) = rhs(i,iv,m) - rhs(i,iv,q) * mat(i,iv,m,q,3)
            end do
          end do
          do iv = 1, nsys
            rhs(i,iv,m) = rhs(i,iv,m) / mat(i,iv,m,m,3)
          end do
        end do

!.... do the next-to-the-last row

        i = n - 1
        p = i + 1
        
        do m = nblk, 1, -1
          do q = m+1, nblk
            do iv = 1, nsys
              rhs(i,iv,m) = rhs(i,iv,m) - rhs(i,iv,q) * mat(i,iv,m,q,3)
            end do
          end do
          do q = 1, nblk
            do iv = 1, nsys
              rhs(i,iv,m) = rhs(i,iv,m) - rhs(p,iv,q) * mat(i,iv,m,q,4)
            end do
          end do
          do iv = 1, nsys
            rhs(i,iv,m) = rhs(i,iv,m) / mat(i,iv,m,m,3)
          end do
        end do

!.... now do the rest of the rows in reverse order

        do i = n - 2, 3, -1
        
          p = i + 1
          r = i + 2
          
          do m = nblk, 1, -1
            do q = m+1, nblk
              do iv = 1, nsys
                rhs(i,iv,m) = rhs(i,iv,m) - rhs(i,iv,q) * mat(i,iv,m,q,3)
              end do
            end do
            do q = 1, nblk
              do iv = 1, nsys
                rhs(i,iv,m) = rhs(i,iv,m) - rhs(p,iv,q) * mat(i,iv,m,q,4)
                rhs(i,iv,m) = rhs(i,iv,m) - rhs(r,iv,q) * mat(i,iv,m,q,5)
              end do
            end do
            do iv = 1, nsys
              rhs(i,iv,m) = rhs(i,iv,m) / mat(i,iv,m,m,3)
            end do
          end do
                                  
        end do
        
!.... Do the last two boundary rows

          i = 2
          p = i + 1
          r = i + 2
          s = i + 3
          
          do m = nblk, 1, -1
            do q = m+1, nblk
              do iv = 1, nsys
                rhs(i,iv,m) = rhs(i,iv,m) - rhs(i,iv,q) * mat(i,iv,m,q,3)
              end do
            end do
            do q = 1, nblk
              do iv = 1, nsys
                rhs(i,iv,m) = rhs(i,iv,m) - rhs(p,iv,q) * mat(i,iv,m,q,4)
                rhs(i,iv,m) = rhs(i,iv,m) - rhs(r,iv,q) * mat(i,iv,m,q,5)
                rhs(i,iv,m) = rhs(i,iv,m) - rhs(s,iv,q) * mat(i,iv,m,q,1)
              end do
            end do
            do iv = 1, nsys
              rhs(i,iv,m) = rhs(i,iv,m) / mat(i,iv,m,m,3)
            end do
          end do

          i = 1
          p = i + 1
          r = i + 2
          s = i + 3
          t = i + 4
          
          do m = nblk, 1, -1
            do q = m+1, nblk
              do iv = 1, nsys
                rhs(i,iv,m) = rhs(i,iv,m) - rhs(i,iv,q) * mat(i,iv,m,q,3)
              end do
            end do
            do q = 1, nblk
              do iv = 1, nsys
                rhs(i,iv,m) = rhs(i,iv,m) - rhs(p,iv,q) * mat(i,iv,m,q,4)
                rhs(i,iv,m) = rhs(i,iv,m) - rhs(r,iv,q) * mat(i,iv,m,q,5)
                rhs(i,iv,m) = rhs(i,iv,m) - rhs(s,iv,q) * mat(i,iv,m,q,1)
                rhs(i,iv,m) = rhs(i,iv,m) - rhs(t,iv,q) * mat(i,iv,m,q,2)
              end do
            end do
            do iv = 1, nsys
              rhs(i,iv,m) = rhs(i,iv,m) / mat(i,iv,m,m,3)
            end do
          end do
!=============================================================================!
        return
        end

