!=============================================================================!
        subroutine cpenta2bc( n, nsys, nblk, mat, rhs, code )
!=============================================================================!
!  
!  Solves a block pentadiagonal system of equations without pivoting.
!  Vectorized over the first index.
!
!  This version includes the Boundary Treatment, and complex math
!
!  WARNING: Since pivoting is not done, this algorithm could fail
!
!  input:
!    nsys       : number of systems to vectorize over
!    n          : number of rows in the systems
!    nblk       : dimension of the blocks in each row
!    code       : 0 = LU factorization only
!               : 1 = Forward and back solves given LU
!               : 2 = LU and solve
!
!  inout:
!    mat        : the block pentadiagonal matrix to be solved
!    rhs        : the right hand side
!
!  Revised:  8-6-96
!=============================================================================!
        implicit none

        integer n, nblk, nsys, code
        complex mat(5,nblk,nblk,n,nsys), rhs(nblk,n,nsys)
        complex mult, fact
        
!.... useful constants

        real zero, one
        parameter (zero = 0.0, one = 1.0)

!.... local variables
        
        integer nc, i, j, l, m, p, q, r, s, t, iv

!=============================================================================!
!                       L U   F A C T O R I Z A T I O N 
!=============================================================================!
        if (code .eq. 0 .or. code .eq. 2) then

!$omp parallel do private(i,l,m,p,q,r,s,t,iv,fact,mult)
        do iv = 1, nsys

!.... do the first row

        i = 1
        
!.... do the first four columns

          do nc = 1, nblk - 1
          
            do m = nc, nblk - 1
              l = m + 1
              mult = -mat(3,l,nc,i,iv) / mat(3,nc,nc,i,iv)
              mat(3,l,nc,i,iv) = mult
              do q = nc + 1, nblk
                mat(3,l,q,i,iv) = mat(3,l,q,i,iv) + mult * mat(3,nc,q,i,iv)
              end do
              do q = 1, nblk
                mat(4,l,q,i,iv) = mat(4,l,q,i,iv) + mult * mat(4,nc,q,i,iv)
                mat(5,l,q,i,iv) = mat(5,l,q,i,iv) + mult * mat(5,nc,q,i,iv)
                mat(1,l,q,i,iv) = mat(1,l,q,i,iv) + mult * mat(1,nc,q,i,iv)
                mat(2,l,q,i,iv) = mat(2,l,q,i,iv) + mult * mat(2,nc,q,i,iv)
              end do
            end do

          end do        ! end of loop on nc

!.... for all columns
          
          do nc = 1, nblk

            p = i + 1
            r = i + 2
            
            fact = one / mat(3,nc,nc,i,iv)
            do m = 1, nblk
              mult = -mat(2,m,nc,p,iv) * fact
              mat(2,m,nc,p,iv) = mult
              do q = nc + 1, nblk
                mat(2,m,q,p,iv) = mat(2,m,q,p,iv) + mult * mat(3,nc,q,i,iv)
              end do
              do q = 1, nblk
                mat(3,m,q,p,iv) = mat(3,m,q,p,iv) + mult * mat(4,nc,q,i,iv)
                mat(4,m,q,p,iv) = mat(4,m,q,p,iv) + mult * mat(5,nc,q,i,iv)
                mat(5,m,q,p,iv) = mat(5,m,q,p,iv) + mult * mat(1,nc,q,i,iv)
                mat(1,m,q,p,iv) = mat(1,m,q,p,iv) + mult * mat(2,nc,q,i,iv)
              end do
              
              mult = -mat(1,m,nc,r,iv) * fact
              mat(1,m,nc,r,iv) = mult
              do q = nc + 1, nblk
                mat(1,m,q,r,iv) = mat(1,m,q,r,iv) + mult * mat(3,nc,q,i,iv)
              end do
              do q = 1, nblk
                mat(2,m,q,r,iv) = mat(2,m,q,r,iv) + mult * mat(4,nc,q,i,iv)
                mat(3,m,q,r,iv) = mat(3,m,q,r,iv) + mult * mat(5,nc,q,i,iv)
                mat(4,m,q,r,iv) = mat(4,m,q,r,iv) + mult * mat(1,nc,q,i,iv)
                mat(5,m,q,r,iv) = mat(5,m,q,r,iv) + mult * mat(2,nc,q,i,iv)
              end do
              
            end do

          end do        ! end of loop on nc

!.... do the second row

        i = 2
        
!.... do the first four columns

          do nc = 1, nblk - 1
          
            do m = nc, nblk - 1
              l = m + 1
              mult = -mat(3,l,nc,i,iv) / mat(3,nc,nc,i,iv)
              mat(3,l,nc,i,iv) = mult
              do q = nc + 1, nblk
                mat(3,l,q,i,iv) = mat(3,l,q,i,iv) + mult * mat(3,nc,q,i,iv)
              end do
              do q = 1, nblk
                mat(4,l,q,i,iv) = mat(4,l,q,i,iv) + mult * mat(4,nc,q,i,iv)
                mat(5,l,q,i,iv) = mat(5,l,q,i,iv) + mult * mat(5,nc,q,i,iv)
                mat(1,l,q,i,iv) = mat(1,l,q,i,iv) + mult * mat(1,nc,q,i,iv)
              end do
            end do            

          end do        ! end of loop on nc

!.... for all columns
          
          do nc = 1, nblk

            p = i + 1
            r = i + 2
            
            fact = one / mat(3,nc,nc,i,iv)
            
            do m = 1, nblk
              mult = -mat(2,m,nc,p,iv) * fact
              mat(2,m,nc,p,iv) = mult
              do q = nc + 1, nblk
                mat(2,m,q,p,iv) = mat(2,m,q,p,iv) + mult * mat(3,nc,q,i,iv)
              end do
              do q = 1, nblk
                mat(3,m,q,p,iv) = mat(3,m,q,p,iv) + mult * mat(4,nc,q,i,iv)
                mat(4,m,q,p,iv) = mat(4,m,q,p,iv) + mult * mat(5,nc,q,i,iv)
                mat(5,m,q,p,iv) = mat(5,m,q,p,iv) + mult * mat(1,nc,q,i,iv)
              end do
              
              mult = -mat(1,m,nc,r,iv) * fact
              mat(1,m,nc,r,iv) = mult
              do q = nc + 1, nblk
                mat(1,m,q,r,iv) = mat(1,m,q,r,iv) + mult * mat(3,nc,q,i,iv)
              end do
              do q = 1, nblk
                mat(2,m,q,r,iv) = mat(2,m,q,r,iv) + mult * mat(4,nc,q,i,iv)
                mat(3,m,q,r,iv) = mat(3,m,q,r,iv) + mult * mat(5,nc,q,i,iv)
                mat(4,m,q,r,iv) = mat(4,m,q,r,iv) + mult * mat(1,nc,q,i,iv)
              end do
              
            end do

          end do        ! end of loop on nc

!.... do the interior rows

        do i = 3, n - 5

!.... do the first four columns

          do nc = 1, nblk - 1
          
            do m = nc, nblk - 1
              l = m + 1
              mult = -mat(3,l,nc,i,iv) / mat(3,nc,nc,i,iv)
              mat(3,l,nc,i,iv) = mult
              do q = nc + 1, nblk
                mat(3,l,q,i,iv) = mat(3,l,q,i,iv) + mult * mat(3,nc,q,i,iv)
              end do
              do q = 1, nblk
                mat(4,l,q,i,iv) = mat(4,l,q,i,iv) + mult * mat(4,nc,q,i,iv)
                mat(5,l,q,i,iv) = mat(5,l,q,i,iv) + mult * mat(5,nc,q,i,iv)
              end do
            end do            

          end do        ! end of loop on nc

!.... for all columns
          
          do nc = 1, nblk

            p = i + 1
            r = i + 2
            
            fact = one / mat(3,nc,nc,i,iv)
            
            do m = 1, nblk
              mult = -mat(2,m,nc,p,iv) * fact
              mat(2,m,nc,p,iv) = mult
              do q = nc + 1, nblk
                mat(2,m,q,p,iv) = mat(2,m,q,p,iv) + mult * mat(3,nc,q,i,iv)
              end do
              do q = 1, nblk
                mat(3,m,q,p,iv) = mat(3,m,q,p,iv) + mult * mat(4,nc,q,i,iv)
                mat(4,m,q,p,iv) = mat(4,m,q,p,iv) + mult * mat(5,nc,q,i,iv)
              end do
              
              mult = -mat(1,m,nc,r,iv) * fact
              mat(1,m,nc,r,iv) = mult
              do q = nc + 1, nblk
                mat(1,m,q,r,iv) = mat(1,m,q,r,iv) + mult * mat(3,nc,q,i,iv)
              end do
              do q = 1, nblk
                mat(2,m,q,r,iv) = mat(2,m,q,r,iv) + mult * mat(4,nc,q,i,iv)
                mat(3,m,q,r,iv) = mat(3,m,q,r,iv) + mult * mat(5,nc,q,i,iv)
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
              mult = -mat(3,l,nc,i,iv) / mat(3,nc,nc,i,iv)
              mat(3,l,nc,i,iv) = mult
              do q = nc + 1, nblk
                mat(3,l,q,i,iv) = mat(3,l,q,i,iv) + mult * mat(3,nc,q,i,iv)
              end do
              do q = 1, nblk
                mat(4,l,q,i,iv) = mat(4,l,q,i,iv) + mult * mat(4,nc,q,i,iv)
                mat(5,l,q,i,iv) = mat(5,l,q,i,iv) + mult * mat(5,nc,q,i,iv)
              end do
            end do            

          end do        ! end of loop on nc

!.... for all columns
          
          do nc = 1, nblk

            p = i + 1
            r = i + 2
            s = i + 3
            t = i + 4
            
            fact = one / mat(3,nc,nc,i,iv)
            
            do m = 1, nblk
              mult = -mat(2,m,nc,p,iv) * fact
              mat(2,m,nc,p,iv) = mult
              do q = nc + 1, nblk
                mat(2,m,q,p,iv) = mat(2,m,q,p,iv) + mult * mat(3,nc,q,i,iv)
              end do
              do q = 1, nblk
                mat(3,m,q,p,iv) = mat(3,m,q,p,iv) + mult * mat(4,nc,q,i,iv)
                mat(4,m,q,p,iv) = mat(4,m,q,p,iv) + mult * mat(5,nc,q,i,iv)
              end do
              
              mult = -mat(1,m,nc,r,iv) * fact
              mat(1,m,nc,r,iv) = mult
              do q = nc + 1, nblk
                mat(1,m,q,r,iv) = mat(1,m,q,r,iv) + mult * mat(3,nc,q,i,iv)
              end do
              do q = 1, nblk
                mat(2,m,q,r,iv) = mat(2,m,q,r,iv) + mult * mat(4,nc,q,i,iv)
                mat(3,m,q,r,iv) = mat(3,m,q,r,iv) + mult * mat(5,nc,q,i,iv)
              end do
              
              mult = -mat(5,m,nc,s,iv) * fact
              mat(5,m,nc,s,iv) = mult
              do q = nc + 1, nblk
                mat(5,m,q,s,iv) = mat(5,m,q,s,iv) + mult * mat(3,nc,q,i,iv)
              end do
              do q = 1, nblk
                mat(1,m,q,s,iv) = mat(1,m,q,s,iv) + mult * mat(4,nc,q,i,iv)
                mat(2,m,q,s,iv) = mat(2,m,q,s,iv) + mult * mat(5,nc,q,i,iv)
              end do
              
              mult = -mat(4,m,nc,t,iv) * fact
              mat(4,m,nc,t,iv) = mult
              do q = nc + 1, nblk
                mat(4,m,q,t,iv) = mat(4,m,q,t,iv) + mult * mat(3,nc,q,i,iv)
              end do
              do q = 1, nblk
                mat(5,m,q,t,iv) = mat(5,m,q,t,iv) + mult * mat(4,nc,q,i,iv)
                mat(1,m,q,t,iv) = mat(1,m,q,t,iv) + mult * mat(5,nc,q,i,iv)
              end do

            end do

          end do        ! end of loop on nc

!.... do the (n-3)rd row

        i = n - 3

!.... do the first four columns

          do nc = 1, nblk - 1
          
            do m = nc, nblk - 1
              l = m + 1
              mult = -mat(3,l,nc,i,iv) / mat(3,nc,nc,i,iv)
              mat(3,l,nc,i,iv) = mult
              do q = nc + 1, nblk
                mat(3,l,q,i,iv) = mat(3,l,q,i,iv) + mult * mat(3,nc,q,i,iv)
              end do
              do q = 1, nblk
                mat(4,l,q,i,iv) = mat(4,l,q,i,iv) + mult * mat(4,nc,q,i,iv)
                mat(5,l,q,i,iv) = mat(5,l,q,i,iv) + mult * mat(5,nc,q,i,iv)
              end do
            end do

          end do        ! end of loop on nc

!.... for all columns
          
          do nc = 1, nblk

            p = i + 1
            r = i + 2
            s = i + 3
            
            fact = one / mat(3,nc,nc,i,iv)
            
            do m = 1, nblk
              mult = -mat(2,m,nc,p,iv) * fact
              mat(2,m,nc,p,iv) = mult
              do q = nc + 1, nblk
                mat(2,m,q,p,iv) = mat(2,m,q,p,iv) + mult * mat(3,nc,q,i,iv)
              end do
              do q = 1, nblk
                mat(3,m,q,p,iv) = mat(3,m,q,p,iv) + mult * mat(4,nc,q,i,iv)
                mat(4,m,q,p,iv) = mat(4,m,q,p,iv) + mult * mat(5,nc,q,i,iv)
              end do
              
              mult = -mat(1,m,nc,r,iv) * fact
              mat(1,m,nc,r,iv) = mult
              do q = nc + 1, nblk
                mat(1,m,q,r,iv) = mat(1,m,q,r,iv) + mult * mat(3,nc,q,i,iv)
              end do
              do q = 1, nblk
                mat(2,m,q,r,iv) = mat(2,m,q,r,iv) + mult * mat(4,nc,q,i,iv)
                mat(3,m,q,r,iv) = mat(3,m,q,r,iv) + mult * mat(5,nc,q,i,iv)
              end do
              
              mult = -mat(5,m,nc,s,iv) * fact
              mat(5,m,nc,s,iv) = mult
              do q = nc + 1, nblk
                mat(5,m,q,s,iv) = mat(5,m,q,s,iv) + mult * mat(3,nc,q,i,iv)
              end do
              do q = 1, nblk
                mat(1,m,q,s,iv) = mat(1,m,q,s,iv) + mult * mat(4,nc,q,i,iv)
                mat(2,m,q,s,iv) = mat(2,m,q,s,iv) + mult * mat(5,nc,q,i,iv)
              end do

            end do

          end do        ! end of loop on nc
          
!.... do the (n-2)nd row

        i = n - 2

!.... do the first four columns

          do nc = 1, nblk - 1
          
            do m = nc, nblk - 1
              l = m + 1
              mult = -mat(3,l,nc,i,iv) / mat(3,nc,nc,i,iv)
              mat(3,l,nc,i,iv) = mult
              do q = nc + 1, nblk
                mat(3,l,q,i,iv) = mat(3,l,q,i,iv) + mult * mat(3,nc,q,i,iv)
              end do
              do q = 1, nblk
                mat(4,l,q,i,iv) = mat(4,l,q,i,iv) + mult * mat(4,nc,q,i,iv)
                mat(5,l,q,i,iv) = mat(5,l,q,i,iv) + mult * mat(5,nc,q,i,iv)
              end do
            end do

          end do        ! end of loop on nc

!.... for all columns
          
          do nc = 1, nblk

            p = i + 1
            r = i + 2
            
            fact = one / mat(3,nc,nc,i,iv)
            
            do m = 1, nblk
              mult = -mat(2,m,nc,p,iv) * fact
              mat(2,m,nc,p,iv) = mult
              do q = nc + 1, nblk
                mat(2,m,q,p,iv) = mat(2,m,q,p,iv) + mult * mat(3,nc,q,i,iv)
              end do
              do q = 1, nblk
                mat(3,m,q,p,iv) = mat(3,m,q,p,iv) + mult * mat(4,nc,q,i,iv)
                mat(4,m,q,p,iv) = mat(4,m,q,p,iv) + mult * mat(5,nc,q,i,iv)
              end do
              
              mult = -mat(1,m,nc,r,iv) * fact
              mat(1,m,nc,r,iv) = mult
              do q = nc + 1, nblk
                mat(1,m,q,r,iv) = mat(1,m,q,r,iv) + mult * mat(3,nc,q,i,iv)
              end do
              do q = 1, nblk
                mat(2,m,q,r,iv) = mat(2,m,q,r,iv) + mult * mat(4,nc,q,i,iv)
                mat(3,m,q,r,iv) = mat(3,m,q,r,iv) + mult * mat(5,nc,q,i,iv)
              end do
              
            end do

          end do        ! end of loop on nc

!.... do the next-to-the-last row
        
        i = n - 1

!.... do the first four columns

          do nc = 1, nblk - 1
          
            do m = nc, nblk - 1
              l = m + 1
              mult = -mat(3,l,nc,i,iv) / mat(3,nc,nc,i,iv)
              mat(3,l,nc,i,iv) = mult
              do q = nc + 1, nblk
                mat(3,l,q,i,iv) = mat(3,l,q,i,iv) + mult * mat(3,nc,q,i,iv)
              end do
              do q = 1, nblk
                mat(4,l,q,i,iv) = mat(4,l,q,i,iv) + mult * mat(4,nc,q,i,iv)
              end do
            end do
            
          end do        ! end of loop on nc

!.... for all columns
          
          do nc = 1, nblk

            p = i + 1
            
            fact = one / mat(3,nc,nc,i,iv)
            
            do m = 1, nblk
              mult = -mat(2,m,nc,p,iv) * fact
              mat(2,m,nc,p,iv) = mult
              do q = nc + 1, nblk
                mat(2,m,q,p,iv) = mat(2,m,q,p,iv) + mult * mat(3,nc,q,i,iv)
              end do
              do q = 1, nblk
                mat(3,m,q,p,iv) = mat(3,m,q,p,iv) + mult * mat(4,nc,q,i,iv)
              end do
            end do
            
          end do        ! end of loop on nc

!.... do the last row

        i = n
        
        do nc = 1, nblk - 1
        
          do m = nc, nblk - 1
            l = m + 1
            mult = -mat(3,l,nc,i,iv) / mat(3,nc,nc,i,iv)
            mat(3,l,nc,i,iv) = mult
            do q = nc + 1, nblk
              mat(3,l,q,i,iv) = mat(3,l,q,i,iv) + mult * mat(3,nc,q,i,iv)
            end do
          end do

        end do  ! end of loop on nc

        end do  ! loop on nsys

        end if
!=============================================================================!
!                     F O R W A R D   E L I M I N A T I O N 
!=============================================================================!
        if (code.eq.1 .or. code.eq.2) then
        
!$omp parallel do private(i,l,m,p,q,r,s,t,iv,fact,mult)
        do iv = 1, nsys

!.... do the first row

        i = 1
        
!.... do the first four columns

          do nc = 1, nblk - 1
          
            do m = nc, nblk - 1
              l = m + 1
              mult = mat(3,l,nc,i,iv)
              rhs(l,i,iv) = rhs(l,i,iv) + mult * rhs(nc,i,iv)
            end do

          end do        ! end of loop on nc

!.... for all columns
          
          do nc = 1, nblk

            p = i + 1
            r = i + 2
            
            do m = 1, nblk
              mult = mat(2,m,nc,p,iv)
              rhs(m,p,iv) = rhs(m,p,iv) + mult * rhs(nc,i,iv)
              
              mult = mat(1,m,nc,r,iv)
              rhs(m,r,iv) = rhs(m,r,iv) + mult * rhs(nc,i,iv)
            end do

          end do        ! end of loop on nc

!.... do the second row

        i = 2
        
!.... do the first four columns

          do nc = 1, nblk - 1
          
            do m = nc, nblk - 1
              l = m + 1
              mult = mat(3,l,nc,i,iv)
              rhs(l,i,iv) = rhs(l,i,iv) + mult * rhs(nc,i,iv)
            end do            

          end do        ! end of loop on nc

!.... for all columns
          
          do nc = 1, nblk

            p = i + 1
            r = i + 2
            
            do m = 1, nblk
              mult = mat(2,m,nc,p,iv)
              rhs(m,p,iv) = rhs(m,p,iv) + mult * rhs(nc,i,iv)
              
              mult = mat(1,m,nc,r,iv)
              rhs(m,r,iv) = rhs(m,r,iv) + mult * rhs(nc,i,iv)
            end do

          end do        ! end of loop on nc

!.... do the interior rows

        do i = 3, n - 5

!.... do the first four columns

          do nc = 1, nblk - 1
          
            do m = nc, nblk - 1
              l = m + 1
              mult = mat(3,l,nc,i,iv)
              rhs(l,i,iv) = rhs(l,i,iv) + mult * rhs(nc,i,iv)
            end do            

          end do        ! end of loop on nc

!.... for all columns
          
          do nc = 1, nblk

            p = i + 1
            r = i + 2
            
            do m = 1, nblk
              mult = mat(2,m,nc,p,iv)
              rhs(m,p,iv) = rhs(m,p,iv) + mult * rhs(nc,i,iv)

              mult = mat(1,m,nc,r,iv)
              rhs(m,r,iv) = rhs(m,r,iv) + mult * rhs(nc,i,iv)
            end do

          end do        ! end of loop on nc
          
        end do          ! end of loop on i

!.... do the (n-4)th row

        i = n - 4

!.... do the first four columns

          do nc = 1, nblk - 1
          
            do m = nc, nblk - 1
              l = m + 1
              mult = mat(3,l,nc,i,iv)
              rhs(l,i,iv) = rhs(l,i,iv) + mult * rhs(nc,i,iv)
            end do            

          end do        ! end of loop on nc

!.... for all columns
          
          do nc = 1, nblk

            p = i + 1
            r = i + 2
            s = i + 3
            t = i + 4
            
            do m = 1, nblk
              mult = mat(2,m,nc,p,iv)
              rhs(m,p,iv) = rhs(m,p,iv) + mult * rhs(nc,i,iv)

              mult = mat(1,m,nc,r,iv)
              rhs(m,r,iv) = rhs(m,r,iv) + mult * rhs(nc,i,iv)

              mult = mat(5,m,nc,s,iv)
              rhs(m,s,iv) = rhs(m,s,iv) + mult * rhs(nc,i,iv)

              mult = mat(4,m,nc,t,iv)
              rhs(m,t,iv) = rhs(m,t,iv) + mult * rhs(nc,i,iv)
            end do

          end do        ! end of loop on nc

!.... do the (n-3)rd row

        i = n - 3

!.... do the first four columns

          do nc = 1, nblk - 1
          
            do m = nc, nblk - 1
              l = m + 1
              mult = mat(3,l,nc,i,iv)
              rhs(l,i,iv) = rhs(l,i,iv) + mult * rhs(nc,i,iv)
            end do            

          end do        ! end of loop on nc

!.... for all columns
          
          do nc = 1, nblk

            p = i + 1
            r = i + 2
            s = i + 3
            
            do m = 1, nblk
              mult = mat(2,m,nc,p,iv)
              rhs(m,p,iv) = rhs(m,p,iv) + mult * rhs(nc,i,iv)
              
              mult = mat(1,m,nc,r,iv)
              rhs(iv,m,r) = rhs(iv,m,r) + mult * rhs(iv,nc,i)
              
              mult = mat(5,m,nc,s,iv)
              rhs(m,s,iv) = rhs(m,s,iv) + mult * rhs(nc,i,iv)
            end do

          end do        ! end of loop on nc

!.... do the (n-2)nd row

        i = n - 2

!.... do the first four columns

          do nc = 1, nblk - 1
          
            do m = nc, nblk - 1
              l = m + 1
              mult = mat(3,l,nc,i,iv)
              rhs(l,i,iv) = rhs(l,i,iv) + mult * rhs(nc,i,iv)
            end do            

          end do        ! end of loop on nc

!.... for all columns
          
          do nc = 1, nblk

            p = i + 1
            r = i + 2
            
            do m = 1, nblk
              mult = mat(2,m,nc,p,iv)
              rhs(m,p,iv) = rhs(m,p,iv) + mult * rhs(nc,i,iv)

              mult = mat(1,m,nc,r,iv)
              rhs(m,r,iv) = rhs(m,r,iv) + mult * rhs(nc,i,iv)
            end do

          end do        ! end of loop on nc

!.... do the next-to-the-last row
        
        i = n - 1

!.... do the first four columns

          do nc = 1, nblk - 1
          
            do m = nc, nblk - 1
              l = m + 1
              mult = mat(3,l,nc,i,iv)
              rhs(l,i,iv) = rhs(l,i,iv) + mult * rhs(nc,i,iv)
            end do            

          end do        ! end of loop on nc

!.... for all columns
          
          do nc = 1, nblk

            p = i + 1
            
            do m = 1, nblk
              mult = mat(2,m,nc,p,iv)
              rhs(m,p,iv) = rhs(m,p,iv) + mult * rhs(nc,i,iv)
            end do

          end do        ! end of loop on nc

!.... do the last row

        i = n
        
        do nc = 1, nblk - 1
        
          do m = nc, nblk - 1
            l = m + 1
            mult = mat(3,l,nc,i,iv)
            rhs(l,i,iv) = rhs(l,i,iv) + mult * rhs(nc,i,iv)
          end do              

        end do  ! end of loop on nc

!=============================================================================!
!                   B A C K W A R D   S U B S T I T U T I O N 
!=============================================================================!

!.... do the last row first

        i = n
        
        do m = nblk, 1, -1
          do q = m+1, nblk
            rhs(m,i,iv) = rhs(m,i,iv) - rhs(q,i,iv) * mat(3,m,q,i,iv)
          end do
          rhs(m,i,iv) = rhs(m,i,iv) / mat(3,m,m,i,iv)
        end do

!.... do the next-to-the-last row

        i = n - 1
        p = i + 1
        
        do m = nblk, 1, -1
          do q = m+1, nblk
            rhs(m,i,iv) = rhs(m,i,iv) - rhs(q,i,iv) * mat(3,m,q,i,iv)
          end do
          do q = 1, nblk
            rhs(m,i,iv) = rhs(m,i,iv) - rhs(q,p,iv) * mat(4,m,q,i,iv)
          end do
          rhs(m,i,iv) = rhs(m,i,iv) / mat(3,m,m,i,iv)
        end do

!.... now do the rest of the rows in reverse order

        do i = n - 2, 3, -1
        
          p = i + 1
          r = i + 2
          
          do m = nblk, 1, -1
            do q = m+1, nblk
              rhs(m,i,iv) = rhs(m,i,iv) - rhs(q,i,iv) * mat(3,m,q,i,iv)
            end do
            do q = 1, nblk
              rhs(m,i,iv) = rhs(m,i,iv) - rhs(q,p,iv) * mat(4,m,q,i,iv)
              rhs(m,i,iv) = rhs(m,i,iv) - rhs(q,r,iv) * mat(5,m,q,i,iv)
            end do
            rhs(m,i,iv) = rhs(m,i,iv) / mat(3,m,m,i,iv)
          end do
                                  
        end do
        
!.... Do the last two boundary rows

          i = 2
          p = i + 1
          r = i + 2
          s = i + 3
          
          do m = nblk, 1, -1
            do q = m+1, nblk
              rhs(m,i,iv) = rhs(m,i,iv) - rhs(q,i,iv) * mat(3,m,q,i,iv)
            end do
            do q = 1, nblk
              rhs(m,i,iv) = rhs(m,i,iv) - rhs(q,p,iv) * mat(4,m,q,i,iv)
              rhs(m,i,iv) = rhs(m,i,iv) - rhs(q,r,iv) * mat(5,m,q,i,iv)
              rhs(m,i,iv) = rhs(m,i,iv) - rhs(q,s,iv) * mat(1,m,q,i,iv)
            end do
            rhs(m,i,iv) = rhs(m,i,iv) / mat(3,m,m,i,iv)
          end do

          i = 1
          p = i + 1
          r = i + 2
          s = i + 3
          t = i + 4
          
          do m = nblk, 1, -1
            do q = m+1, nblk
              rhs(m,i,iv) = rhs(m,i,iv) - rhs(q,i,iv) * mat(3,m,q,i,iv)
            end do
            do q = 1, nblk
              rhs(m,i,iv) = rhs(m,i,iv) - rhs(q,p,iv) * mat(4,m,q,i,iv)
              rhs(m,i,iv) = rhs(m,i,iv) - rhs(q,r,iv) * mat(5,m,q,i,iv)
              rhs(m,i,iv) = rhs(m,i,iv) - rhs(q,s,iv) * mat(1,m,q,i,iv)
              rhs(m,i,iv) = rhs(m,i,iv) - rhs(q,t,iv) * mat(2,m,q,i,iv)
            end do
            rhs(m,i,iv) = rhs(m,i,iv) / mat(3,m,m,i,iv)
          end do
        
        end do

        end if
!=============================================================================!
        return
        end

