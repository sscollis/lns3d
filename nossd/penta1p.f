!=============================================================================!
        subroutine penta1p( nsys, np, nblk, mat, rhs, per, per2, mult, fact )
!=============================================================================!
!  
!  Solves a block pentadiagonal system of equations without pivoting.
!  Vectorized over the first index.
!
!  This routine solves a PERIODIC block pentadiagonal system.
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
!=============================================================================!
        implicit none

        integer np, nblk, nsys
        real mat(nsys,np,nblk,nblk,5), rhs(nsys,np,nblk)
        real per(nsys,np,nblk,nblk,2), per2(nsys,nblk,nblk,6)
        real mult(nsys), fact(nsys)
        
!.... useful constants

        real zero, one
        parameter (zero = 0.0, one = 1.0)

!.... local variables
        
        integer nc, i, j, l, m, n, p, q, r, iv
!=============================================================================!

!.... don't solve for the redundant nodes

        n = np - 1

!.... setup the periodic fillin vector

        do q = 1, nblk
          do p = 1, nblk
            do iv = 1, nsys
              per(iv,1,p,q,1) = mat(iv,1,p,q,1)
              per(iv,1,p,q,2) = mat(iv,1,p,q,2)
              per(iv,2,p,q,1) = zero
              per(iv,2,p,q,2) = mat(iv,2,p,q,1)
              
              per2(iv,p,q,1)  = mat(iv,n-1,p,q,5)
              per2(iv,p,q,2)  = zero
              per2(iv,p,q,4)  = mat(iv,n,p,q,4)
              per2(iv,p,q,5)  = mat(iv,n,p,q,5)
            end do
          end do
        end do

!=============================================================================!
!                     F O R W A R D   E L I M I N A T I O N 
!=============================================================================!

        do i = 1, n - 6

!.... do the first four columns

          do nc = 1, nblk - 1
          
            do m = nc, nblk - 1
              l = m + 1
              do iv = 1, nsys
                mult(iv) = -mat(iv,i,l,nc,3) / mat(iv,i,nc,nc,3)
              end do
              do q = nc + 1, nblk
                do iv = 1, nsys
                  mat(iv,i,l,q,3) = mat(iv,i,l,q,3) + mult(iv) * mat(iv,i,nc,q,3)
                end do
              end do
              do q = 1, nblk
                do iv = 1, nsys
                  mat(iv,i,l,q,4) = mat(iv,i,l,q,4) + mult(iv) * mat(iv,i,nc,q,4)
                  mat(iv,i,l,q,5) = mat(iv,i,l,q,5) + mult(iv) * mat(iv,i,nc,q,5)
                  per(iv,i,l,q,1) = per(iv,i,l,q,1) + mult(iv) * per(iv,i,nc,q,1)
                  per(iv,i,l,q,2) = per(iv,i,l,q,2) + mult(iv) * per(iv,i,nc,q,2)
                end do
              end do
              do iv = 1, nsys
                rhs(iv,i,l) = rhs(iv,i,l) + mult(iv) * rhs(iv,i,nc)
              end do
            end do            

          end do        ! end of loop on nc

!.... for all columns
          
          nc = 1
            
          p = i + 1
          r = i + 2
          
          do iv = 1, nsys
            fact(iv) = one / mat(iv,i,nc,nc,3)
          end do
          do m = 1, nblk
            do iv = 1, nsys
              mult(iv) = -mat(iv,p,m,nc,2) * fact(iv)
            end do
            do q = nc + 1, nblk
              do iv = 1, nsys
                mat(iv,p,m,q,2) = mat(iv,p,m,q,2) + mult(iv) * mat(iv,i,nc,q,3)
              end do
            end do
            do q = 1, nblk
              do iv = 1, nsys
                mat(iv,p,m,q,3) = mat(iv,p,m,q,3) + mult(iv) * mat(iv,i,nc,q,4)
                mat(iv,p,m,q,4) = mat(iv,p,m,q,4) + mult(iv) * mat(iv,i,nc,q,5)
                per(iv,p,m,q,1) = per(iv,p,m,q,1) + mult(iv) * per(iv,i,nc,q,1)
                per(iv,p,m,q,2) = per(iv,p,m,q,2) + mult(iv) * per(iv,i,nc,q,2)
              end do
            end do
            do iv = 1, nsys
              rhs(iv,p,m) = rhs(iv,p,m) + mult(iv) * rhs(iv,i,nc)
            end do
            
            do iv = 1, nsys
              mult(iv) = -mat(iv,r,m,nc,1) * fact(iv)
            end do
            do q = nc + 1, nblk
              do iv = 1, nsys
                mat(iv,r,m,q,1) = mat(iv,r,m,q,1) + mult(iv) * mat(iv,i,nc,q,3)
              end do
            end do
            do q = 1, nblk
              do iv = 1, nsys
                mat(iv,r,m,q,2) = mat(iv,r,m,q,2) + mult(iv) * mat(iv,i,nc,q,4)
                mat(iv,r,m,q,3) = mat(iv,r,m,q,3) + mult(iv) * mat(iv,i,nc,q,5)
                per(iv,r,m,q,1) = mult(iv) * per(iv,i,nc,q,1)
                per(iv,r,m,q,2) = mult(iv) * per(iv,i,nc,q,2)
              end do
            end do
            do iv = 1, nsys
              rhs(iv,r,m) = rhs(iv,r,m) + mult(iv) * rhs(iv,i,nc)
            end do

!.... the periodic stuff on the next-to-the-last row

            do iv = 1, nsys
              mult(iv) = -per2(iv,m,nc,1) * fact(iv)
            end do
            do q = nc + 1, nblk
              do iv = 1, nsys
                per2(iv,m,q,1) = per2(iv,m,q,1) + mult(iv) * mat(iv,i,nc,q,3)
              end do
            end do
            do q = 1, nblk
              do iv = 1, nsys
                per2(iv,m,q,2) = per2(iv,m,q,2) + mult(iv) * mat(iv,i,nc,q,4)
                per2(iv,m,q,3) = mult(iv) * mat(iv,i,nc,q,5)
                mat(iv,n-1,m,q,3) = mat(iv,n-1,m,q,3) + mult(iv) * per(iv,i,nc,q,1)
                mat(iv,n-1,m,q,4) = mat(iv,n-1,m,q,4) + mult(iv) * per(iv,i,nc,q,2)
              end do
            end do
            do iv = 1, nsys
              rhs(iv,n-1,m) = rhs(iv,n-1,m) + mult(iv) * rhs(iv,i,nc)
            end do
      
!.... the periodic stuff on the last row

            do iv = 1, nsys
              mult(iv) = -per2(iv,m,nc,4) * fact(iv)
            end do
            do q = nc + 1, nblk
              do iv = 1, nsys
                per2(iv,m,q,4) = per2(iv,m,q,4) + mult(iv) * mat(iv,i,nc,q,3)
              end do
            end do
            do q = 1, nblk
              do iv = 1, nsys
                per2(iv,m,q,5) = per2(iv,m,q,5) + mult(iv) * mat(iv,i,nc,q,4)
                per2(iv,m,q,6) = mult(iv) * mat(iv,i,nc,q,5)
                mat(iv,n,m,q,2) = mat(iv,n,m,q,2) + mult(iv) * per(iv,i,nc,q,1)
                mat(iv,n,m,q,3) = mat(iv,n,m,q,3) + mult(iv) * per(iv,i,nc,q,2)
              end do
            end do
            do iv = 1, nsys
              rhs(iv,n,m) = rhs(iv,n,m) + mult(iv) * rhs(iv,i,nc)
            end do

          end do        ! end of loop on m

          do nc = 2, nblk

            p = i + 1
            r = i + 2
            
            do iv = 1, nsys
              fact(iv) = one / mat(iv,i,nc,nc,3)
            end do
            do m = 1, nblk
              do iv = 1, nsys
                mult(iv) = -mat(iv,p,m,nc,2) * fact(iv)
              end do
              do q = nc + 1, nblk
                do iv = 1, nsys
                  mat(iv,p,m,q,2) = mat(iv,p,m,q,2) + mult(iv) * mat(iv,i,nc,q,3)
                end do
              end do
              do q = 1, nblk
                do iv = 1, nsys
                  mat(iv,p,m,q,3) = mat(iv,p,m,q,3) + mult(iv) * mat(iv,i,nc,q,4)
                  mat(iv,p,m,q,4) = mat(iv,p,m,q,4) + mult(iv) * mat(iv,i,nc,q,5)
                  per(iv,p,m,q,1) = per(iv,p,m,q,1) + mult(iv) * per(iv,i,nc,q,1)
                  per(iv,p,m,q,2) = per(iv,p,m,q,2) + mult(iv) * per(iv,i,nc,q,2)
                end do
              end do
              do iv = 1, nsys
                rhs(iv,p,m) = rhs(iv,p,m) + mult(iv) * rhs(iv,i,nc)
              end do
              
              do iv = 1, nsys
                mult(iv) = -mat(iv,r,m,nc,1) * fact(iv)
              end do
              do q = nc + 1, nblk
                do iv = 1, nsys
                  mat(iv,r,m,q,1) = mat(iv,r,m,q,1) + mult(iv) * mat(iv,i,nc,q,3)
                end do
              end do
              do q = 1, nblk
                do iv = 1, nsys
                  mat(iv,r,m,q,2) = mat(iv,r,m,q,2) + mult(iv) * mat(iv,i,nc,q,4)
                  mat(iv,r,m,q,3) = mat(iv,r,m,q,3) + mult(iv) * mat(iv,i,nc,q,5)
                  per(iv,r,m,q,1) = per(iv,r,m,q,1) + mult(iv) * per(iv,i,nc,q,1)
                  per(iv,r,m,q,2) = per(iv,r,m,q,2) + mult(iv) * per(iv,i,nc,q,2)
                end do
              end do
              do iv = 1, nsys
                rhs(iv,r,m) = rhs(iv,r,m) + mult(iv) * rhs(iv,i,nc)
              end do

!.... the periodic stuff on the next-to-the-last row

              do iv = 1, nsys
                mult(iv) = -per2(iv,m,nc,1) * fact(iv)
              end do
              do q = nc + 1, nblk
                do iv = 1, nsys
                  per2(iv,m,q,1) = per2(iv,m,q,1) + mult(iv) * mat(iv,i,nc,q,3)
                end do
              end do
              do q = 1, nblk
                do iv = 1, nsys
                  per2(iv,m,q,2) = per2(iv,m,q,2) + mult(iv) * mat(iv,i,nc,q,4)
                  per2(iv,m,q,3) = per2(iv,m,q,3) + mult(iv) * mat(iv,i,nc,q,5)
                  mat(iv,n-1,m,q,3) = mat(iv,n-1,m,q,3) + mult(iv) * per(iv,i,nc,q,1)
                  mat(iv,n-1,m,q,4) = mat(iv,n-1,m,q,4) + mult(iv) * per(iv,i,nc,q,2)
                end do
              end do
              do iv = 1, nsys
                rhs(iv,n-1,m) = rhs(iv,n-1,m) + mult(iv) * rhs(iv,i,nc)
              end do
        
!.... the periodic stuff on the last row

              do iv = 1, nsys
                mult(iv) = -per2(iv,m,nc,4) * fact(iv)
              end do
              do q = nc + 1, nblk
                do iv = 1, nsys
                  per2(iv,m,q,4) = per2(iv,m,q,4) + mult(iv) * mat(iv,i,nc,q,3)
                end do
              end do
              do q = 1, nblk
                do iv = 1, nsys
                  per2(iv,m,q,5) = per2(iv,m,q,5) + mult(iv) * mat(iv,i,nc,q,4)
                  per2(iv,m,q,6) = per2(iv,m,q,6) + mult(iv) * mat(iv,i,nc,q,5)
                  mat(iv,n,m,q,2) = mat(iv,n,m,q,2) + mult(iv) * per(iv,i,nc,q,1)
                  mat(iv,n,m,q,3) = mat(iv,n,m,q,3) + mult(iv) * per(iv,i,nc,q,2)
                end do
              end do
              do iv = 1, nsys
                rhs(iv,n,m) = rhs(iv,n,m) + mult(iv) * rhs(iv,i,nc)
              end do

            end do      ! end of loop on m

          end do        ! end of loop on nc
          
          do q = 1, nblk
            do p = 1, nblk
              do iv = 1, nsys
                per2(iv,p,q,1)  = per2(iv,p,q,2)
                per2(iv,p,q,2)  = per2(iv,p,q,3)
                per2(iv,p,q,4)  = per2(iv,p,q,5)
                per2(iv,p,q,5)  = per2(iv,p,q,6)
              end do
            end do
          end do

        end do          ! end of loop on i

!====================================================================================
        i = n - 5

!.... do the first four columns

          do nc = 1, nblk - 1
          
            do m = nc, nblk - 1
              l = m + 1
              do iv = 1, nsys
                mult(iv) = -mat(iv,i,l,nc,3) / mat(iv,i,nc,nc,3)
              end do
              do q = nc + 1, nblk
                do iv = 1, nsys
                  mat(iv,i,l,q,3) = mat(iv,i,l,q,3) + mult(iv) * mat(iv,i,nc,q,3)
                end do
              end do
              do q = 1, nblk
                do iv = 1, nsys
                  mat(iv,i,l,q,4) = mat(iv,i,l,q,4) + mult(iv) * mat(iv,i,nc,q,4)
                  mat(iv,i,l,q,5) = mat(iv,i,l,q,5) + mult(iv) * mat(iv,i,nc,q,5)
                  per(iv,i,l,q,1) = per(iv,i,l,q,1) + mult(iv) * per(iv,i,nc,q,1)
                  per(iv,i,l,q,2) = per(iv,i,l,q,2) + mult(iv) * per(iv,i,nc,q,2)
                end do
              end do
              do iv = 1, nsys
                rhs(iv,i,l) = rhs(iv,i,l) + mult(iv) * rhs(iv,i,nc)
              end do
            end do            

          end do        ! end of loop on nc

!.... for all columns
          
          nc = 1
          
          p = i + 1
          r = i + 2
          
          do iv = 1, nsys
            fact(iv) = one / mat(iv,i,nc,nc,3)
          end do
          do m = 1, nblk
            do iv = 1, nsys
              mult(iv) = -mat(iv,p,m,nc,2) * fact(iv)
            end do
            do q = nc + 1, nblk
              do iv = 1, nsys
                mat(iv,p,m,q,2) = mat(iv,p,m,q,2) + mult(iv) * mat(iv,i,nc,q,3)
              end do
            end do
            do q = 1, nblk
              do iv = 1, nsys
                mat(iv,p,m,q,3) = mat(iv,p,m,q,3) + mult(iv) * mat(iv,i,nc,q,4)
                mat(iv,p,m,q,4) = mat(iv,p,m,q,4) + mult(iv) * mat(iv,i,nc,q,5)
                per(iv,p,m,q,1) = per(iv,p,m,q,1) + mult(iv) * per(iv,i,nc,q,1)
                per(iv,p,m,q,2) = per(iv,p,m,q,2) + mult(iv) * per(iv,i,nc,q,2)
              end do
            end do
            do iv = 1, nsys
              rhs(iv,p,m) = rhs(iv,p,m) + mult(iv) * rhs(iv,i,nc)
            end do
            
            do iv = 1, nsys
              mult(iv) = -mat(iv,r,m,nc,1) * fact(iv)
            end do
            do q = nc + 1, nblk
              do iv = 1, nsys
                mat(iv,r,m,q,1) = mat(iv,r,m,q,1) + mult(iv) * mat(iv,i,nc,q,3)
              end do
            end do
            do q = 1, nblk
              do iv = 1, nsys
                mat(iv,r,m,q,2) = mat(iv,r,m,q,2) + mult(iv) * mat(iv,i,nc,q,4)
                mat(iv,r,m,q,3) = mat(iv,r,m,q,3) + mult(iv) * mat(iv,i,nc,q,5)
                mat(iv,r,m,q,4) = mat(iv,r,m,q,4) + mult(iv) * per(iv,i,nc,q,1)
                per(iv,r,m,q,2) = mult(iv) * per(iv,i,nc,q,2)
              end do
            end do
            do iv = 1, nsys
              rhs(iv,r,m) = rhs(iv,r,m) + mult(iv) * rhs(iv,i,nc)
            end do

!.... the periodic stuff on the next-to-the-last row

            do iv = 1, nsys
              mult(iv) = -per2(iv,m,nc,1) * fact(iv)
            end do
            do q = nc + 1, nblk
              do iv = 1, nsys
                per2(iv,m,q,1) = per2(iv,m,q,1) + mult(iv) * mat(iv,i,nc,q,3)
              end do
            end do
            do q = 1, nblk
              do iv = 1, nsys
                per2(iv,m,q,2) = per2(iv,m,q,2) + mult(iv) * mat(iv,i,nc,q,4)
                mat(iv,n-1,m,q,1) = mat(iv,n-1,m,q,1) + mult(iv) * mat(iv,i,nc,q,5)
                mat(iv,n-1,m,q,3) = mat(iv,n-1,m,q,3) + mult(iv) * per(iv,i,nc,q,1)
                mat(iv,n-1,m,q,4) = mat(iv,n-1,m,q,4) + mult(iv) * per(iv,i,nc,q,2)
              end do
            end do
            do iv = 1, nsys
              rhs(iv,n-1,m) = rhs(iv,n-1,m) + mult(iv) * rhs(iv,i,nc)
            end do
      
!.... the periodic stuff on the last row

            do iv = 1, nsys
              mult(iv) = -per2(iv,m,nc,4) * fact(iv)
            end do
            do q = nc + 1, nblk
              do iv = 1, nsys
                per2(iv,m,q,4) = per2(iv,m,q,4) + mult(iv) * mat(iv,i,nc,q,3)
              end do
            end do
            do q = 1, nblk
              do iv = 1, nsys
                per2(iv,m,q,5) = per2(iv,m,q,5) + mult(iv) * mat(iv,i,nc,q,4)
                per2(iv,m,q,6) = mult(iv) * mat(iv,i,nc,q,5)
                mat(iv,n,m,q,2) = mat(iv,n,m,q,2) + mult(iv) * per(iv,i,nc,q,1)
                mat(iv,n,m,q,3) = mat(iv,n,m,q,3) + mult(iv) * per(iv,i,nc,q,2)
              end do
            end do
            do iv = 1, nsys
              rhs(iv,n,m) = rhs(iv,n,m) + mult(iv) * rhs(iv,i,nc)
            end do

          end do        ! end of loop on m

          do nc = 2, nblk

            p = i + 1
            r = i + 2
            
            do iv = 1, nsys
              fact(iv) = one / mat(iv,i,nc,nc,3)
            end do
            do m = 1, nblk
              do iv = 1, nsys
                mult(iv) = -mat(iv,p,m,nc,2) * fact(iv)
              end do
              do q = nc + 1, nblk
                do iv = 1, nsys
                  mat(iv,p,m,q,2) = mat(iv,p,m,q,2) + mult(iv) * mat(iv,i,nc,q,3)
                end do
              end do
              do q = 1, nblk
                do iv = 1, nsys
                  mat(iv,p,m,q,3) = mat(iv,p,m,q,3) + mult(iv) * mat(iv,i,nc,q,4)
                  mat(iv,p,m,q,4) = mat(iv,p,m,q,4) + mult(iv) * mat(iv,i,nc,q,5)
                  per(iv,p,m,q,1) = per(iv,p,m,q,1) + mult(iv) * per(iv,i,nc,q,1)
                  per(iv,p,m,q,2) = per(iv,p,m,q,2) + mult(iv) * per(iv,i,nc,q,2)
                end do
              end do
              do iv = 1, nsys
                rhs(iv,p,m) = rhs(iv,p,m) + mult(iv) * rhs(iv,i,nc)
              end do
              
              do iv = 1, nsys
                mult(iv) = -mat(iv,r,m,nc,1) * fact(iv)
              end do
              do q = nc + 1, nblk
                do iv = 1, nsys
                  mat(iv,r,m,q,1) = mat(iv,r,m,q,1) + mult(iv) * mat(iv,i,nc,q,3)
                end do
              end do
              do q = 1, nblk
                do iv = 1, nsys
                  mat(iv,r,m,q,2) = mat(iv,r,m,q,2) + mult(iv) * mat(iv,i,nc,q,4)
                  mat(iv,r,m,q,3) = mat(iv,r,m,q,3) + mult(iv) * mat(iv,i,nc,q,5)
                  mat(iv,r,m,q,4) = mat(iv,r,m,q,4) + mult(iv) * per(iv,i,nc,q,1)
                  per(iv,r,m,q,2) = per(iv,r,m,q,2) + mult(iv) * per(iv,i,nc,q,2)
                end do
              end do
              do iv = 1, nsys
                rhs(iv,r,m) = rhs(iv,r,m) + mult(iv) * rhs(iv,i,nc)
              end do

!.... the periodic stuff on the next-to-the-last row

              do iv = 1, nsys
                mult(iv) = -per2(iv,m,nc,1) * fact(iv)
              end do
              do q = nc + 1, nblk
                do iv = 1, nsys
                  per2(iv,m,q,1) = per2(iv,m,q,1) + mult(iv) * mat(iv,i,nc,q,3)
                end do
              end do
              do q = 1, nblk
                do iv = 1, nsys
                  per2(iv,m,q,2) = per2(iv,m,q,2) + mult(iv) * mat(iv,i,nc,q,4)
                  mat(iv,n-1,m,q,1) = mat(iv,n-1,m,q,1) + mult(iv) * mat(iv,i,nc,q,5)
                  mat(iv,n-1,m,q,3) = mat(iv,n-1,m,q,3) + mult(iv) * per(iv,i,nc,q,1)
                  mat(iv,n-1,m,q,4) = mat(iv,n-1,m,q,4) + mult(iv) * per(iv,i,nc,q,2)
                end do
              end do
              do iv = 1, nsys
                rhs(iv,n-1,m) = rhs(iv,n-1,m) + mult(iv) * rhs(iv,i,nc)
              end do
        
!.... the periodic stuff on the last row

              do iv = 1, nsys
                mult(iv) = -per2(iv,m,nc,4) * fact(iv)
              end do
              do q = nc + 1, nblk
                do iv = 1, nsys
                  per2(iv,m,q,4) = per2(iv,m,q,4) + mult(iv) * mat(iv,i,nc,q,3)
                end do
              end do
              do q = 1, nblk
                do iv = 1, nsys
                  per2(iv,m,q,5) = per2(iv,m,q,5) + mult(iv) * mat(iv,i,nc,q,4)
                  per2(iv,m,q,6) = per2(iv,m,q,6) + mult(iv) * mat(iv,i,nc,q,5)
                  mat(iv,n,m,q,2) = mat(iv,n,m,q,2) + mult(iv) * per(iv,i,nc,q,1)
                  mat(iv,n,m,q,3) = mat(iv,n,m,q,3) + mult(iv) * per(iv,i,nc,q,2)
                end do
              end do
              do iv = 1, nsys
                rhs(iv,n,m) = rhs(iv,n,m) + mult(iv) * rhs(iv,i,nc)
              end do

            end do      ! end of loop on m

          end do        ! end of loop on nc
          
          do q = 1, nblk
            do p = 1, nblk
              do iv = 1, nsys
                per2(iv,p,q,1)  = per2(iv,p,q,2)
                per2(iv,p,q,4)  = per2(iv,p,q,5)
                per2(iv,p,q,5)  = per2(iv,p,q,6)
              end do
            end do
          end do
        
!====================================================================================
        i = n - 4

!.... do the first four columns

          do nc = 1, nblk - 1
          
            do m = nc, nblk - 1
              l = m + 1
              do iv = 1, nsys
                mult(iv) = -mat(iv,i,l,nc,3) / mat(iv,i,nc,nc,3)
              end do
              do q = nc + 1, nblk
                do iv = 1, nsys
                  mat(iv,i,l,q,3) = mat(iv,i,l,q,3) + mult(iv) * mat(iv,i,nc,q,3)
                end do
              end do
              do q = 1, nblk
                do iv = 1, nsys
                  mat(iv,i,l,q,4) = mat(iv,i,l,q,4) + mult(iv) * mat(iv,i,nc,q,4)
                  mat(iv,i,l,q,5) = mat(iv,i,l,q,5) + mult(iv) * mat(iv,i,nc,q,5)
                  per(iv,i,l,q,1) = per(iv,i,l,q,1) + mult(iv) * per(iv,i,nc,q,1)
                  per(iv,i,l,q,2) = per(iv,i,l,q,2) + mult(iv) * per(iv,i,nc,q,2)
                end do
              end do
              do iv = 1, nsys
                rhs(iv,i,l) = rhs(iv,i,l) + mult(iv) * rhs(iv,i,nc)
              end do
            end do            

          end do        ! end of loop on nc

!.... for all columns
          
          do nc = 1, nblk

            p = i + 1
            r = i + 2
            
            do iv = 1, nsys
              fact(iv) = one / mat(iv,i,nc,nc,3)
            end do
            do m = 1, nblk
              do iv = 1, nsys
                mult(iv) = -mat(iv,p,m,nc,2) * fact(iv)
              end do
              do q = nc + 1, nblk
                do iv = 1, nsys
                  mat(iv,p,m,q,2) = mat(iv,p,m,q,2) + mult(iv) * mat(iv,i,nc,q,3)
                end do
              end do
              do q = 1, nblk
                do iv = 1, nsys
                  mat(iv,p,m,q,3) = mat(iv,p,m,q,3) + mult(iv) * mat(iv,i,nc,q,4)
                  mat(iv,p,m,q,4) = mat(iv,p,m,q,4) + mult(iv) * mat(iv,i,nc,q,5)
                  mat(iv,p,m,q,5) = mat(iv,p,m,q,5) + mult(iv) * per(iv,i,nc,q,1)
                  per(iv,p,m,q,2) = per(iv,p,m,q,2) + mult(iv) * per(iv,i,nc,q,2)
                end do
              end do
              do iv = 1, nsys
                rhs(iv,p,m) = rhs(iv,p,m) + mult(iv) * rhs(iv,i,nc)
              end do
              
              do iv = 1, nsys
                mult(iv) = -mat(iv,r,m,nc,1) * fact(iv)
              end do
              do q = nc + 1, nblk
                do iv = 1, nsys
                  mat(iv,r,m,q,1) = mat(iv,r,m,q,1) + mult(iv) * mat(iv,i,nc,q,3)
                end do
              end do
              do q = 1, nblk
                do iv = 1, nsys
                  mat(iv,r,m,q,2) = mat(iv,r,m,q,2) + mult(iv) * mat(iv,i,nc,q,4)
                  mat(iv,r,m,q,3) = mat(iv,r,m,q,3) + mult(iv) * mat(iv,i,nc,q,5)
                  mat(iv,r,m,q,4) = mat(iv,r,m,q,4) + mult(iv) * per(iv,i,nc,q,1)
                  mat(iv,r,m,q,5) = mat(iv,r,m,q,5) + mult(iv) * per(iv,i,nc,q,2)
                end do
              end do
              do iv = 1, nsys
                rhs(iv,r,m) = rhs(iv,r,m) + mult(iv) * rhs(iv,i,nc)
              end do

!.... the periodic stuff on the next-to-the-last row

              do iv = 1, nsys
                mult(iv) = -per2(iv,m,nc,1) * fact(iv)
              end do
              do q = nc + 1, nblk
                do iv = 1, nsys
                  per2(iv,m,q,1) = per2(iv,m,q,1) + mult(iv) * mat(iv,i,nc,q,3)
                end do
              end do
              do q = 1, nblk
                do iv = 1, nsys
                  mat(iv,n-1,m,q,1) = mat(iv,n-1,m,q,1) + mult(iv) * mat(iv,i,nc,q,4)
                  mat(iv,n-1,m,q,2) = mat(iv,n-1,m,q,2) + mult(iv) * mat(iv,i,nc,q,5)
                  mat(iv,n-1,m,q,3) = mat(iv,n-1,m,q,3) + mult(iv) * per(iv,i,nc,q,1)
                  mat(iv,n-1,m,q,4) = mat(iv,n-1,m,q,4) + mult(iv) * per(iv,i,nc,q,2)
                end do
              end do
              do iv = 1, nsys
                rhs(iv,n-1,m) = rhs(iv,n-1,m) + mult(iv) * rhs(iv,i,nc)
              end do
        
!.... the periodic stuff on the last row

              do iv = 1, nsys
                mult(iv) = -per2(iv,m,nc,4) * fact(iv)
              end do
              do q = nc + 1, nblk
                do iv = 1, nsys
                  per2(iv,m,q,4) = per2(iv,m,q,4) + mult(iv) * mat(iv,i,nc,q,3)
                end do
              end do
              do q = 1, nblk
                do iv = 1, nsys
                  per2(iv,m,q,5) = per2(iv,m,q,5) + mult(iv) * mat(iv,i,nc,q,4)
                  mat(iv,n,m,q,1) = mat(iv,n,m,q,1) + mult(iv) * mat(iv,i,nc,q,5)
                  mat(iv,n,m,q,2) = mat(iv,n,m,q,2) + mult(iv) * per(iv,i,nc,q,1)
                  mat(iv,n,m,q,3) = mat(iv,n,m,q,3) + mult(iv) * per(iv,i,nc,q,2)
                end do
              end do
              do iv = 1, nsys
                rhs(iv,n,m) = rhs(iv,n,m) + mult(iv) * rhs(iv,i,nc)
              end do

            end do      ! end of loop on m

          end do        ! end of loop on nc
          
          do q = 1, nblk
            do p = 1, nblk
              do iv = 1, nsys
                per2(iv,p,q,4)  = per2(iv,p,q,5)
              end do
            end do
          end do

!====================================================================================
        i = n - 3

!.... do the first four columns

          do nc = 1, nblk - 1
          
            do m = nc, nblk - 1
              l = m + 1
              do iv = 1, nsys
                mult(iv) = -mat(iv,i,l,nc,3) / mat(iv,i,nc,nc,3)
              end do
              do q = nc + 1, nblk
                do iv = 1, nsys
                  mat(iv,i,l,q,3) = mat(iv,i,l,q,3) + mult(iv) * mat(iv,i,nc,q,3)
                end do
              end do
              do q = 1, nblk
                do iv = 1, nsys
                  mat(iv,i,l,q,4) = mat(iv,i,l,q,4) + mult(iv) * mat(iv,i,nc,q,4)
                  mat(iv,i,l,q,5) = mat(iv,i,l,q,5) + mult(iv) * mat(iv,i,nc,q,5)
                  per(iv,i,l,q,2) = per(iv,i,l,q,2) + mult(iv) * per(iv,i,nc,q,2)
                end do
              end do
              do iv = 1, nsys
                rhs(iv,i,l) = rhs(iv,i,l) + mult(iv) * rhs(iv,i,nc)
              end do
            end do            

          end do        ! end of loop on nc

!.... for all columns
          
          do nc = 1, nblk

            p = i + 1
            r = i + 2
            
            do iv = 1, nsys
              fact(iv) = one / mat(iv,i,nc,nc,3)
            end do
            do m = 1, nblk
              do iv = 1, nsys
                mult(iv) = -mat(iv,p,m,nc,2) * fact(iv)
              end do
              do q = nc + 1, nblk
                do iv = 1, nsys
                  mat(iv,p,m,q,2) = mat(iv,p,m,q,2) + mult(iv) * mat(iv,i,nc,q,3)
                end do
              end do
              do q = 1, nblk
                do iv = 1, nsys
                  mat(iv,p,m,q,3) = mat(iv,p,m,q,3) + mult(iv) * mat(iv,i,nc,q,4)
                  mat(iv,p,m,q,4) = mat(iv,p,m,q,4) + mult(iv) * mat(iv,i,nc,q,5)
                  mat(iv,p,m,q,5) = mat(iv,p,m,q,5) + mult(iv) * per(iv,i,nc,q,2)
                end do
              end do
              do iv = 1, nsys
                rhs(iv,p,m) = rhs(iv,p,m) + mult(iv) * rhs(iv,i,nc)
              end do
              
              do iv = 1, nsys
                mult(iv) = -mat(iv,r,m,nc,1) * fact(iv)
              end do
              do q = nc + 1, nblk
                do iv = 1, nsys
                  mat(iv,r,m,q,1) = mat(iv,r,m,q,1) + mult(iv) * mat(iv,i,nc,q,3)
                end do
              end do
              do q = 1, nblk
                do iv = 1, nsys
                  mat(iv,r,m,q,2) = mat(iv,r,m,q,2) + mult(iv) * mat(iv,i,nc,q,4)
                  mat(iv,r,m,q,3) = mat(iv,r,m,q,3) + mult(iv) * mat(iv,i,nc,q,5)
                  mat(iv,r,m,q,4) = mat(iv,r,m,q,4) + mult(iv) * per(iv,i,nc,q,2)
                end do
              end do
              do iv = 1, nsys
                rhs(iv,r,m) = rhs(iv,r,m) + mult(iv) * rhs(iv,i,nc)
              end do

!.... the periodic stuff on the last row

              do iv = 1, nsys
                mult(iv) = -per2(iv,m,nc,4) * fact(iv)
              end do
              do q = nc + 1, nblk
                do iv = 1, nsys
                  per2(iv,m,q,4) = per2(iv,m,q,4) + mult(iv) * mat(iv,i,nc,q,3)
                end do
              end do
              do q = 1, nblk
                do iv = 1, nsys
                  mat(iv,n,m,q,1) = mat(iv,n,m,q,1) + mult(iv) * mat(iv,i,nc,q,4)
                  mat(iv,n,m,q,2) = mat(iv,n,m,q,2) + mult(iv) * mat(iv,i,nc,q,5)
                  mat(iv,n,m,q,3) = mat(iv,n,m,q,3) + mult(iv) * per(iv,i,nc,q,2)
                end do
              end do
              do iv = 1, nsys
                rhs(iv,n,m) = rhs(iv,n,m) + mult(iv) * rhs(iv,i,nc)
              end do

            end do      ! end of loop on m

          end do        ! end of loop on nc
          
!====================================================================================
        i = n - 2

!.... do the first four columns

          do nc = 1, nblk - 1
          
            do m = nc, nblk - 1
              l = m + 1
              do iv = 1, nsys
                mult(iv) = -mat(iv,i,l,nc,3) / mat(iv,i,nc,nc,3)
              end do
              do q = nc + 1, nblk
                do iv = 1, nsys
                  mat(iv,i,l,q,3) = mat(iv,i,l,q,3) + mult(iv) * mat(iv,i,nc,q,3)
                end do
              end do
              do q = 1, nblk
                do iv = 1, nsys
                  mat(iv,i,l,q,4) = mat(iv,i,l,q,4) + mult(iv) * mat(iv,i,nc,q,4)
                  mat(iv,i,l,q,5) = mat(iv,i,l,q,5) + mult(iv) * mat(iv,i,nc,q,5)
                end do
              end do
              do iv = 1, nsys
                rhs(iv,i,l) = rhs(iv,i,l) + mult(iv) * rhs(iv,i,nc)
              end do
            end do            

          end do        ! end of loop on nc

!.... for all columns
          
          do nc = 1, nblk

            p = i + 1
            r = i + 2
            
            do iv = 1, nsys
              fact(iv) = one / mat(iv,i,nc,nc,3)
            end do
            do m = 1, nblk
              do iv = 1, nsys
                mult(iv) = -mat(iv,p,m,nc,2) * fact(iv)
              end do
              do q = nc + 1, nblk
                do iv = 1, nsys
                  mat(iv,p,m,q,2) = mat(iv,p,m,q,2) + mult(iv) * mat(iv,i,nc,q,3)
                end do
              end do
              do q = 1, nblk
                do iv = 1, nsys
                  mat(iv,p,m,q,3) = mat(iv,p,m,q,3) + mult(iv) * mat(iv,i,nc,q,4)
                  mat(iv,p,m,q,4) = mat(iv,p,m,q,4) + mult(iv) * mat(iv,i,nc,q,5)
                end do
              end do
              do iv = 1, nsys
                rhs(iv,p,m) = rhs(iv,p,m) + mult(iv) * rhs(iv,i,nc)
              end do
              
              do iv = 1, nsys
                mult(iv) = -mat(iv,r,m,nc,1) * fact(iv)
              end do
              do q = nc + 1, nblk
                do iv = 1, nsys
                  mat(iv,r,m,q,1) = mat(iv,r,m,q,1) + mult(iv) * mat(iv,i,nc,q,3)
                end do
              end do
              do q = 1, nblk
                do iv = 1, nsys
                  mat(iv,r,m,q,2) = mat(iv,r,m,q,2) + mult(iv) * mat(iv,i,nc,q,4)
                  mat(iv,r,m,q,3) = mat(iv,r,m,q,3) + mult(iv) * mat(iv,i,nc,q,5)
                end do
              end do
              do iv = 1, nsys
                rhs(iv,r,m) = rhs(iv,r,m) + mult(iv) * rhs(iv,i,nc)
              end do

            end do      ! end of loop on m

          end do        ! end of loop on nc

!====================================================================================
        i = n - 1

!.... do the first four columns

          do nc = 1, nblk - 1
          
            do m = nc, nblk - 1
              l = m + 1
              do iv = 1, nsys
                mult(iv) = -mat(iv,i,l,nc,3) / mat(iv,i,nc,nc,3)
              end do
              do q = nc + 1, nblk
                do iv = 1, nsys
                  mat(iv,i,l,q,3) = mat(iv,i,l,q,3) + mult(iv) * mat(iv,i,nc,q,3)
                end do
              end do
              do q = 1, nblk
                do iv = 1, nsys
                  mat(iv,i,l,q,4) = mat(iv,i,l,q,4) + mult(iv) * mat(iv,i,nc,q,4)
                end do
              end do
              do iv = 1, nsys
                rhs(iv,i,l) = rhs(iv,i,l) + mult(iv) * rhs(iv,i,nc)
              end do
            end do            

          end do        ! end of loop on nc

!.... for all columns
          
          do nc = 1, nblk

            p = i + 1
            
            do iv = 1, nsys
              fact(iv) = one / mat(iv,i,nc,nc,3)
            end do
            
            do m = 1, nblk
              do iv = 1, nsys
                mult(iv) = -mat(iv,p,m,nc,2) * fact(iv)
              end do
              do q = nc + 1, nblk
                do iv = 1, nsys
                  mat(iv,p,m,q,2) = mat(iv,p,m,q,2) + mult(iv) * mat(iv,i,nc,q,3)
                end do
              end do
              do q = 1, nblk
                do iv = 1, nsys
                  mat(iv,p,m,q,3) = mat(iv,p,m,q,3) + mult(iv) * mat(iv,i,nc,q,4)
                end do
              end do
              do iv = 1, nsys
                rhs(iv,p,m) = rhs(iv,p,m) + mult(iv) * rhs(iv,i,nc)
              end do
            end do

          end do        ! end of loop on nc

!.... do the last row

        i = n
        
        do nc = 1, nblk - 1
        
          do m = nc, nblk - 1
            l = m + 1
            do iv = 1, nsys
              mult(iv) = -mat(iv,i,l,nc,3) / mat(iv,i,nc,nc,3)
            end do
            do q = nc + 1, nblk
              do iv = 1, nsys
                mat(iv,i,l,q,3) = mat(iv,i,l,q,3) + mult(iv) * mat(iv,i,nc,q,3)
              end do
            end do
            do iv = 1, nsys
              rhs(iv,i,l) = rhs(iv,i,l) + mult(iv) * rhs(iv,i,nc)
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
              rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,i,q) * mat(iv,i,m,q,3)
            end do
          end do
          do iv = 1, nsys
            rhs(iv,i,m) = rhs(iv,i,m) / mat(iv,i,m,m,3)
          end do
        end do

!.... do the next-to-the-last row

        i = n - 1
        p = i + 1
        
        do m = nblk, 1, -1
          do q = m+1, nblk
            do iv = 1, nsys
              rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,i,q) * mat(iv,i,m,q,3)
            end do
          end do
          do q = 1, nblk
            do iv = 1, nsys
              rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,p,q) * mat(iv,i,m,q,4)
            end do
          end do
          do iv = 1, nsys
            rhs(iv,i,m) = rhs(iv,i,m) / mat(iv,i,m,m,3)
          end do
        end do

        i = n - 2
        p = i + 1
        r = i + 2
        
        do m = nblk, 1, -1
          do q = m+1, nblk
            do iv = 1, nsys
              rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,i,q) * mat(iv,i,m,q,3)
            end do
          end do
          do q = 1, nblk
            do iv = 1, nsys
              rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,p,q) * mat(iv,i,m,q,4)
              rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,r,q) * mat(iv,i,m,q,5)
            end do
          end do
          do iv = 1, nsys
            rhs(iv,i,m) = rhs(iv,i,m) / mat(iv,i,m,m,3)
          end do
        end do

        i = n - 3
        p = i + 1
        r = i + 2
        
        do m = nblk, 1, -1
          do q = m+1, nblk
            do iv = 1, nsys
              rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,i,q) * mat(iv,i,m,q,3)
            end do
          end do
          do q = 1, nblk
            do iv = 1, nsys
              rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,p,q) * mat(iv,i,m,q,4)
              rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,r,q) * mat(iv,i,m,q,5)
              rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,n,q) * per(iv,i,m,q,2)
            end do
          end do
          do iv = 1, nsys
            rhs(iv,i,m) = rhs(iv,i,m) / mat(iv,i,m,m,3)
          end do
        end do

!.... now do the rest of the rows in reverse order

        do i = n - 4, 1, -1
        
          p = i + 1
          r = i + 2
          
          do m = nblk, 1, -1
            do q = m+1, nblk
              do iv = 1, nsys
                rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,i,q) * mat(iv,i,m,q,3)
              end do
            end do
            do q = 1, nblk
              do iv = 1, nsys
                rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,p,q) * mat(iv,i,m,q,4)
                rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,r,q) * mat(iv,i,m,q,5)
                rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,n-1,q) * per(iv,i,m,q,1)
                rhs(iv,i,m) = rhs(iv,i,m) - rhs(iv,n  ,q) * per(iv,i,m,q,2)
              end do
            end do
            do iv = 1, nsys
              rhs(iv,i,m) = rhs(iv,i,m) / mat(iv,i,m,m,3)
            end do
          end do
                                  
        end do

!.... account for periodicity

        do q = 1, nblk
          do iv = 1, nsys
            rhs(iv,np,q) = rhs(iv,1,q)
          end do
        end do
!=============================================================================!
        return
        end
