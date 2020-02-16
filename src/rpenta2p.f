!=============================================================================!
        subroutine rpenta2p( np, nsys, nblk, mat, rhs, per, per2, code )
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
!    code       : 0 = LU factorization only
!               : 1 = Forward and back solves given LU
!               : 2 = LU and solve
!
!  inout:
!    mat        : the block pentadiagonal matrix to be solved
!    rhs        : the right hand side
!=============================================================================!
        implicit none

        integer np, nblk, nsys, code
        real mat(5,nblk,nblk,np,nsys), rhs(nblk,np,nsys)
        real per(2,nblk,nblk,np,nsys), per2(np-2:np-1,nblk,nblk,np,nsys)
        real mult, fact
        
!.... useful constants

        real zero, one
        parameter (zero = 0.0, one = 1.0)

!.... local variables
        
        integer nc, i, j, l, m, n, p, q, r, iv
!=============================================================================!

!.... don't solve for the redundant nodes

        n = np - 1

!=============================================================================!
!                       L U   F A C T O R I Z A T I O N 
!=============================================================================!
        if (code .eq. 0 .or. code .eq. 2) then

!$omp parallel do private(nc,i,l,m,p,q,r,iv,fact,mult)
        do iv = 1, nsys

        do i = 1, np
           do q = 1, nblk
              do p = 1, nblk
                 per(1,p,q,i,iv) = zero
                 per(2,p,q,i,iv) = zero
                 per2(n-1,p,q,i,iv) = zero
                 per2(n,p,q,i,iv) = zero
              end do
           end do
        end do

!.... setup the periodic fillin vector

        do q = 1, nblk
          do p = 1, nblk
            per(1,p,q,1,iv) = mat(1,p,q,1,iv)
            per(2,p,q,1,iv) = mat(2,p,q,1,iv)
            per(1,p,q,2,iv) = zero
            per(2,p,q,2,iv) = mat(1,p,q,2,iv)
            
            per2(n-1,p,q,1,iv) = mat(5,p,q,n-1,iv)
            per2(n-1,p,q,2,iv) = zero
            per2(n,p,q,1,iv)   = mat(4,p,q,n,iv)
            per2(n,p,q,2,iv)   = mat(5,p,q,n,iv)
          end do
        end do

        do i = 1, n - 6

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
                per(1,l,q,i,iv) = per(1,l,q,i,iv) + mult * per(1,nc,q,i,iv)
                per(2,l,q,i,iv) = per(2,l,q,i,iv) + mult * per(2,nc,q,i,iv)
              end do
            end do

          end do        ! end of loop on nc

!.... for all columns
          
          nc = 1
            
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
              per(1,m,q,p,iv) = per(1,m,q,p,iv) + mult * per(1,nc,q,i,iv)
              per(2,m,q,p,iv) = per(2,m,q,p,iv) + mult * per(2,nc,q,i,iv)
            end do
            
            mult = -mat(1,m,nc,r,iv) * fact
            mat(1,m,nc,r,iv) = mult
            do q = nc + 1, nblk
              mat(1,m,q,r,iv) = mat(1,m,q,r,iv) + mult * mat(3,nc,q,i,iv)
            end do
            do q = 1, nblk
              mat(2,m,q,r,iv) = mat(2,m,q,r,iv) + mult * mat(4,nc,q,i,iv)
              mat(3,m,q,r,iv) = mat(3,m,q,r,iv) + mult * mat(5,nc,q,i,iv)
              per(1,m,q,r,iv) = mult * per(1,nc,q,i,iv)
              per(2,m,q,r,iv) = mult * per(2,nc,q,i,iv)
            end do

!.... the periodic stuff on the next-to-the-last row

            mult = -per2(n-1,m,nc,i,iv) * fact
            per2(n-1,m,nc,i,iv) = mult
            do q = nc + 1, nblk
              per2(n-1,m,q,i,iv) = per2(n-1,m,q,i,iv) + mult * mat(3,nc,q,i,iv)
            end do
            do q = 1, nblk
              per2(n-1,m,q,p,iv) = per2(n-1,m,q,p,iv) + mult * mat(4,nc,q,i,iv)
              per2(n-1,m,q,r,iv) = mult * mat(5,nc,q,i,iv)
              mat(3,m,q,n-1,iv) = mat(3,m,q,n-1,iv) + mult * per(1,nc,q,i,iv)
              mat(4,m,q,n-1,iv) = mat(4,m,q,n-1,iv) + mult * per(2,nc,q,i,iv)
            end do
      
!.... the periodic stuff on the last row

            mult = -per2(n,m,nc,i,iv) * fact
            per2(n,m,nc,i,iv) = mult
            do q = nc + 1, nblk
              per2(n,m,q,i,iv) = per2(n,m,q,i,iv) + mult * mat(3,nc,q,i,iv)
            end do
            do q = 1, nblk
              per2(n,m,q,p,iv) = per2(n,m,q,p,iv) + mult * mat(4,nc,q,i,iv)
              per2(n,m,q,r,iv) = mult * mat(5,nc,q,i,iv)
              mat(2,m,q,n,iv) = mat(2,m,q,n,iv) + mult * per(1,nc,q,i,iv)
              mat(3,m,q,n,iv) = mat(3,m,q,n,iv) + mult * per(2,nc,q,i,iv)
            end do

          end do        ! end of loop on m

          do nc = 2, nblk

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
                per(1,m,q,p,iv) = per(1,m,q,p,iv) + mult * per(1,nc,q,i,iv)
                per(2,m,q,p,iv) = per(2,m,q,p,iv) + mult * per(2,nc,q,i,iv)
              end do
              
              mult = -mat(1,m,nc,r,iv) * fact
              mat(1,m,nc,r,iv) = mult
              do q = nc + 1, nblk
                mat(1,m,q,r,iv) = mat(1,m,q,r,iv) + mult * mat(3,nc,q,i,iv)
              end do
              do q = 1, nblk
                mat(2,m,q,r,iv) = mat(2,m,q,r,iv) + mult * mat(4,nc,q,i,iv)
                mat(3,m,q,r,iv) = mat(3,m,q,r,iv) + mult * mat(5,nc,q,i,iv)
                per(1,m,q,r,iv) = per(1,m,q,r,iv) + mult * per(1,nc,q,i,iv)
                per(2,m,q,r,iv) = per(2,m,q,r,iv) + mult * per(2,nc,q,i,iv)
              end do

!.... the periodic stuff on the next-to-the-last row

              mult = -per2(n-1,m,nc,i,iv) * fact
              per2(n-1,m,nc,i,iv) = mult
              do q = nc + 1, nblk
                per2(n-1,m,q,i,iv) = per2(n-1,m,q,i,iv) + mult * mat(3,nc,q,i,iv)
              end do
              do q = 1, nblk
                per2(n-1,m,q,p,iv) = per2(n-1,m,q,p,iv) + mult * mat(4,nc,q,i,iv)
                per2(n-1,m,q,r,iv) = per2(n-1,m,q,r,iv) + mult * mat(5,nc,q,i,iv)
                mat(3,m,q,n-1,iv) = mat(3,m,q,n-1,iv) + mult * per(1,nc,q,i,iv)
                mat(4,m,q,n-1,iv) = mat(4,m,q,n-1,iv) + mult * per(2,nc,q,i,iv)
              end do
        
!.... the periodic stuff on the last row

              mult = -per2(n,m,nc,i,iv) * fact
              per2(n,m,nc,i,iv) = mult
              do q = nc + 1, nblk
                per2(n,m,q,i,iv) = per2(n,m,q,i,iv) + mult * mat(3,nc,q,i,iv)
              end do
              do q = 1, nblk
                per2(n,m,q,p,iv) = per2(n,m,q,p,iv) + mult * mat(4,nc,q,i,iv)
                per2(n,m,q,r,iv) = per2(n,m,q,r,iv) + mult * mat(5,nc,q,i,iv)
                mat(2,m,q,n,iv) = mat(2,m,q,n,iv) + mult * per(1,nc,q,i,iv)
                mat(3,m,q,n,iv) = mat(3,m,q,n,iv) + mult * per(2,nc,q,i,iv)
              end do

            end do      ! end of loop on m

          end do        ! end of loop on nc
          
        end do          ! end of loop on i

!====================================================================================
        i = n - 5

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
                per(1,l,q,i,iv) = per(1,l,q,i,iv) + mult * per(1,nc,q,i,iv)
                per(2,l,q,i,iv) = per(2,l,q,i,iv) + mult * per(2,nc,q,i,iv)
              end do
            end do

          end do        ! end of loop on nc

!.... for all columns
          
          nc = 1
          
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
              per(1,m,q,p,iv) = per(1,m,q,p,iv) + mult * per(1,nc,q,i,iv)
              per(2,m,q,p,iv) = per(2,m,q,p,iv) + mult * per(2,nc,q,i,iv)
            end do
            
            mult = -mat(1,m,nc,r,iv) * fact
            mat(1,m,nc,r,iv) = mult
            do q = nc + 1, nblk
              mat(1,m,q,r,iv) = mat(1,m,q,r,iv) + mult * mat(3,nc,q,i,iv)
            end do
            do q = 1, nblk
              mat(2,m,q,r,iv) = mat(2,m,q,r,iv) + mult * mat(4,nc,q,i,iv)
              mat(3,m,q,r,iv) = mat(3,m,q,r,iv) + mult * mat(5,nc,q,i,iv)
              mat(4,m,q,r,iv) = mat(4,m,q,r,iv) + mult * per(1,nc,q,i,iv)
              per(2,m,q,r,iv) = mult * per(2,nc,q,i,iv)
            end do

!.... the periodic stuff on the next-to-the-last row

            mult = -per2(n-1,m,nc,i,iv) * fact
            per2(n-1,m,nc,i,iv) = mult
            do q = nc + 1, nblk
              per2(n-1,m,q,i,iv) = per2(n-1,m,q,i,iv) + mult * mat(3,nc,q,i,iv)
            end do
            do q = 1, nblk
              per2(n-1,m,q,p,iv) = per2(n-1,m,q,p,iv) + mult * mat(4,nc,q,i,iv)
              mat(1,m,q,n-1,iv) = mat(1,m,q,n-1,iv) + mult * mat(5,nc,q,i,iv)
              mat(3,m,q,n-1,iv) = mat(3,m,q,n-1,iv) + mult * per(1,nc,q,i,iv)
              mat(4,m,q,n-1,iv) = mat(4,m,q,n-1,iv) + mult * per(2,nc,q,i,iv)
            end do
      
!.... the periodic stuff on the last row

            mult = -per2(n,m,nc,i,iv) * fact
            per2(n,m,nc,i,iv) = mult
            do q = nc + 1, nblk
              per2(n,m,q,i,iv) = per2(n,m,q,i,iv) + mult * mat(3,nc,q,i,iv)
            end do
            do q = 1, nblk
              per2(n,m,q,p,iv) = per2(n,m,q,p,iv) + mult * mat(4,nc,q,i,iv)
              per2(n,m,q,r,iv) = mult * mat(5,nc,q,i,iv)
              mat(2,m,q,n,iv) = mat(2,m,q,n,iv) + mult * per(1,nc,q,i,iv)
              mat(3,m,q,n,iv) = mat(3,m,q,n,iv) + mult * per(2,nc,q,i,iv)
            end do
            
          end do        ! end of loop on m

          do nc = 2, nblk

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
                per(1,m,q,p,iv) = per(1,m,q,p,iv) + mult * per(1,nc,q,i,iv)
                per(2,m,q,p,iv) = per(2,m,q,p,iv) + mult * per(2,nc,q,i,iv)
              end do
              
              mult = -mat(1,m,nc,r,iv) * fact
              mat(1,m,nc,r,iv) = mult
              do q = nc + 1, nblk
                mat(1,m,q,r,iv) = mat(1,m,q,r,iv) + mult * mat(3,nc,q,i,iv)
              end do
              do q = 1, nblk
                mat(2,m,q,r,iv) = mat(2,m,q,r,iv) + mult * mat(4,nc,q,i,iv)
                mat(3,m,q,r,iv) = mat(3,m,q,r,iv) + mult * mat(5,nc,q,i,iv)
                mat(4,m,q,r,iv) = mat(4,m,q,r,iv) + mult * per(1,nc,q,i,iv)
                per(2,m,q,r,iv) = per(2,m,q,r,iv) + mult * per(2,nc,q,i,iv)
              end do

!.... the periodic stuff on the next-to-the-last row

              mult = -per2(n-1,m,nc,i,iv) * fact
              per2(n-1,m,nc,i,iv) = mult
              do q = nc + 1, nblk
                per2(n-1,m,q,i,iv) = per2(n-1,m,q,i,iv) + mult * mat(3,nc,q,i,iv)
              end do
              do q = 1, nblk
                per2(n-1,m,q,p,iv) = per2(n-1,m,q,p,iv) + mult * mat(4,nc,q,i,iv)
                mat(1,m,q,n-1,iv) = mat(1,m,q,n-1,iv) + mult * mat(5,nc,q,i,iv)
                mat(3,m,q,n-1,iv) = mat(3,m,q,n-1,iv) + mult * per(1,nc,q,i,iv)
                mat(4,m,q,n-1,iv) = mat(4,m,q,n-1,iv) + mult * per(2,nc,q,i,iv)
              end do
        
!.... the periodic stuff on the last row

              mult = -per2(n,m,nc,i,iv) * fact
              per2(n,m,nc,i,iv) = mult
              do q = nc + 1, nblk
                per2(n,m,q,i,iv) = per2(n,m,q,i,iv) + mult * mat(3,nc,q,i,iv)
              end do
              do q = 1, nblk
                per2(n,m,q,p,iv) = per2(n,m,q,p,iv) + mult * mat(4,nc,q,i,iv)
                per2(n,m,q,r,iv) = per2(n,m,q,r,iv) + mult * mat(5,nc,q,i,iv)
                mat(2,m,q,n,iv) = mat(2,m,q,n,iv) + mult * per(1,nc,q,i,iv)
                mat(3,m,q,n,iv) = mat(3,m,q,n,iv) + mult * per(2,nc,q,i,iv)
              end do

            end do      ! end of loop on m

          end do        ! end of loop on nc
          
!====================================================================================
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
                per(1,l,q,i,iv) = per(1,l,q,i,iv) + mult * per(1,nc,q,i,iv)
                per(2,l,q,i,iv) = per(2,l,q,i,iv) + mult * per(2,nc,q,i,iv)
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
                mat(5,m,q,p,iv) = mat(5,m,q,p,iv) + mult * per(1,nc,q,i,iv)
                per(2,m,q,p,iv) = per(2,m,q,p,iv) + mult * per(2,nc,q,i,iv)
              end do
              
              mult = -mat(1,m,nc,r,iv) * fact
              mat(1,m,nc,r,iv) = mult
              do q = nc + 1, nblk
                mat(1,m,q,r,iv) = mat(1,m,q,r,iv) + mult * mat(3,nc,q,i,iv)
              end do
              do q = 1, nblk
                mat(2,m,q,r,iv) = mat(2,m,q,r,iv) + mult * mat(4,nc,q,i,iv)
                mat(3,m,q,r,iv) = mat(3,m,q,r,iv) + mult * mat(5,nc,q,i,iv)
                mat(4,m,q,r,iv) = mat(4,m,q,r,iv) + mult * per(1,nc,q,i,iv)
                mat(5,m,q,r,iv) = mat(5,m,q,r,iv) + mult * per(2,nc,q,i,iv)
              end do

!.... the periodic stuff on the next-to-the-last row

              mult = -per2(n-1,m,nc,i,iv) * fact
              per2(n-1,m,nc,i,iv) = mult
              do q = nc + 1, nblk
                per2(n-1,m,q,i,iv) = per2(n-1,m,q,i,iv) + mult * mat(3,nc,q,i,iv)
              end do
              do q = 1, nblk
                mat(1,m,q,n-1,iv) = mat(1,m,q,n-1,iv) + mult * mat(4,nc,q,i,iv)
                mat(2,m,q,n-1,iv) = mat(2,m,q,n-1,iv) + mult * mat(5,nc,q,i,iv)
                mat(3,m,q,n-1,iv) = mat(3,m,q,n-1,iv) + mult * per(1,nc,q,i,iv)
                mat(4,m,q,n-1,iv) = mat(4,m,q,n-1,iv) + mult * per(2,nc,q,i,iv)
              end do
        
!.... the periodic stuff on the last row

              mult = -per2(n,m,nc,i,iv) * fact
              per2(n,m,nc,i,iv) = mult
              do q = nc + 1, nblk
                per2(n,m,q,i,iv) = per2(n,m,q,i,iv) + mult * mat(3,nc,q,i,iv)
              end do
              do q = 1, nblk
                per2(n,m,q,p,iv) = per2(n,m,q,p,iv) + mult * mat(4,nc,q,i,iv)
                mat(1,m,q,n,iv) = mat(1,m,q,n,iv) + mult * mat(5,nc,q,i,iv)
                mat(2,m,q,n,iv) = mat(2,m,q,n,iv) + mult * per(1,nc,q,i,iv)
                mat(3,m,q,n,iv) = mat(3,m,q,n,iv) + mult * per(2,nc,q,i,iv)
              end do

            end do      ! end of loop on m

          end do        ! end of loop on nc
          
!====================================================================================
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
                per(2,l,q,i,iv) = per(2,l,q,i,iv) + mult * per(2,nc,q,i,iv)
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
                mat(5,m,q,p,iv) = mat(5,m,q,p,iv) + mult * per(2,nc,q,i,iv)
              end do
              
              mult = -mat(1,m,nc,r,iv) * fact
              mat(1,m,nc,r,iv) = mult
              do q = nc + 1, nblk
                mat(1,m,q,r,iv) = mat(1,m,q,r,iv) + mult * mat(3,nc,q,i,iv)
              end do
              do q = 1, nblk
                mat(2,m,q,r,iv) = mat(2,m,q,r,iv) + mult * mat(4,nc,q,i,iv)
                mat(3,m,q,r,iv) = mat(3,m,q,r,iv) + mult * mat(5,nc,q,i,iv)
                mat(4,m,q,r,iv) = mat(4,m,q,r,iv) + mult * per(2,nc,q,i,iv)
              end do

!.... the periodic stuff on the last row

              mult = -per2(n,m,nc,i,iv) * fact
              per2(n,m,nc,i,iv) = mult
              do q = nc + 1, nblk
                per2(n,m,q,i,iv) = per2(n,m,q,i,iv) + mult * mat(3,nc,q,i,iv)
              end do
              do q = 1, nblk
                mat(1,m,q,n,iv) = mat(1,m,q,n,iv) + mult * mat(4,nc,q,i,iv)
                mat(2,m,q,n,iv) = mat(2,m,q,n,iv) + mult * mat(5,nc,q,i,iv)
                mat(3,m,q,n,iv) = mat(3,m,q,n,iv) + mult * per(2,nc,q,i,iv)
              end do

            end do      ! end of loop on m

          end do        ! end of loop on nc
          
!====================================================================================
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

            end do      ! end of loop on m

          end do        ! end of loop on nc

!====================================================================================
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

        end do  ! loop over nsys

        end if
!=============================================================================!
!                     F O R W A R D   E L I M I N A T I O N 
!=============================================================================!
        if (code.eq.1 .or. code.eq.2) then
        
!$omp parallel do private(nc,i,l,m,p,q,r,iv,fact,mult)
        do iv = 1, nsys

        do i = 1, n - 6

!.... do the first four columns

          do nc = 1, nblk - 1
          
            do m = nc, nblk - 1
              l = m + 1
              mult = mat(3,l,nc,i,iv)
              rhs(l,i,iv) = rhs(l,i,iv) + mult * rhs(nc,i,iv)
            end do            

          end do        ! end of loop on nc

!.... for all columns
          
          nc = 1
            
          p = i + 1
          r = i + 2
          
          do m = 1, nblk
            mult = mat(2,m,nc,p,iv)
            rhs(m,p,iv) = rhs(m,p,iv) + mult * rhs(nc,i,iv)
            
            mult = mat(1,m,nc,r,iv)
            rhs(m,r,iv) = rhs(m,r,iv) + mult * rhs(nc,i,iv)

!.... the periodic stuff on the next-to-the-last row

            mult = per2(n-1,m,nc,i,iv)
            rhs(m,n-1,iv) = rhs(m,n-1,iv) + mult * rhs(nc,i,iv)
      
!.... the periodic stuff on the last row

            mult = per2(n,m,nc,i,iv)
            rhs(m,n,iv) = rhs(m,n,iv) + mult * rhs(nc,i,iv)

          end do        ! end of loop on m

          do nc = 2, nblk

            p = i + 1
            r = i + 2
            
            do m = 1, nblk
              mult = mat(2,m,nc,p,iv)
              rhs(m,p,iv) = rhs(m,p,iv) + mult * rhs(nc,i,iv)
              
              mult = mat(1,m,nc,r,iv)
              rhs(m,r,iv) = rhs(m,r,iv) + mult * rhs(nc,i,iv)

!.... the periodic stuff on the next-to-the-last row

              mult = per2(n-1,m,nc,i,iv)
              rhs(m,n-1,iv) = rhs(m,n-1,iv) + mult * rhs(nc,i,iv)
        
!.... the periodic stuff on the last row

              mult = per2(n,m,nc,i,iv)
              rhs(m,n,iv) = rhs(m,n,iv) + mult * rhs(nc,i,iv)

            end do      ! end of loop on m

          end do        ! end of loop on nc
          
        end do          ! end of loop on i

!====================================================================================
        i = n - 5

!.... do the first four columns

          do nc = 1, nblk - 1
          
            do m = nc, nblk - 1
              l = m + 1
              mult = mat(3,l,nc,i,iv)
              rhs(l,i,iv) = rhs(l,i,iv) + mult * rhs(nc,i,iv)
            end do            

          end do        ! end of loop on nc

!.... for all columns
          
          nc = 1
          
          p = i + 1
          r = i + 2
          
          do m = 1, nblk
            mult = mat(2,m,nc,p,iv)
            rhs(m,p,iv) = rhs(m,p,iv) + mult * rhs(nc,i,iv)
            
            mult = mat(1,m,nc,r,iv)
            rhs(m,r,iv) = rhs(m,r,iv) + mult * rhs(nc,i,iv)

!.... the periodic stuff on the next-to-the-last row

            mult = per2(n-1,m,nc,i,iv)
            rhs(m,n-1,iv) = rhs(m,n-1,iv) + mult * rhs(nc,i,iv)
      
!.... the periodic stuff on the last row

            mult = per2(n,m,nc,i,iv)
            rhs(m,n,iv) = rhs(m,n,iv) + mult * rhs(nc,i,iv)
          end do        ! end of loop on m

          do nc = 2, nblk

            p = i + 1
            r = i + 2
            
            do m = 1, nblk
              mult = mat(2,m,nc,p,iv)
              rhs(m,p,iv) = rhs(m,p,iv) + mult * rhs(nc,i,iv)
              
              mult = mat(1,m,nc,r,iv)
              rhs(m,r,iv) = rhs(m,r,iv) + mult * rhs(nc,i,iv)

!.... the periodic stuff on the next-to-the-last row

              mult = per2(n-1,m,nc,i,iv)
              rhs(m,n-1,iv) = rhs(m,n-1,iv) + mult * rhs(nc,i,iv)
        
!.... the periodic stuff on the last row

              mult = per2(n,m,nc,i,iv)
              rhs(m,n,iv) = rhs(m,n,iv) + mult * rhs(nc,i,iv)

            end do      ! end of loop on m

          end do        ! end of loop on nc
          
!====================================================================================
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
            
            do m = 1, nblk
              mult = mat(2,m,nc,p,iv)
              rhs(m,p,iv) = rhs(m,p,iv) + mult * rhs(nc,i,iv)
              
              mult = mat(1,m,nc,r,iv)
              rhs(m,r,iv) = rhs(m,r,iv) + mult * rhs(nc,i,iv)

!.... the periodic stuff on the next-to-the-last row

              mult = per2(n-1,m,nc,i,iv)
              rhs(m,n-1,iv) = rhs(m,n-1,iv) + mult * rhs(nc,i,iv)
        
!.... the periodic stuff on the last row

              mult = per2(n,m,nc,i,iv)
              rhs(m,n,iv) = rhs(m,n,iv) + mult * rhs(nc,i,iv)

            end do      ! end of loop on m

          end do        ! end of loop on nc
          
!====================================================================================
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
            
            do m = 1, nblk
              mult = mat(2,m,nc,p,iv)
              rhs(m,p,iv) = rhs(m,p,iv) + mult * rhs(nc,i,iv)
              
              mult = mat(1,m,nc,r,iv)
              rhs(m,r,iv) = rhs(m,r,iv) + mult * rhs(nc,i,iv)

!.... the periodic stuff on the last row

              mult = per2(n,m,nc,i,iv)
              rhs(m,n,iv) = rhs(m,n,iv) + mult * rhs(nc,i,iv)

            end do      ! end of loop on m

          end do        ! end of loop on nc
          
!====================================================================================
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
            end do      ! end of loop on m

          end do        ! end of loop on nc

!====================================================================================
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

        i = n - 2
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

        i = n - 3
        p = i + 1
        r = i + 2
        
        do m = nblk, 1, -1
          do q = m+1, nblk
            rhs(m,i,iv) = rhs(m,i,iv) - rhs(q,i,iv) * mat(3,m,q,i,iv)
          end do
          do q = 1, nblk
            rhs(m,i,iv) = rhs(m,i,iv) - rhs(q,p,iv) * mat(4,m,q,i,iv)
            rhs(m,i,iv) = rhs(m,i,iv) - rhs(q,r,iv) * mat(5,m,q,i,iv)
            rhs(m,i,iv) = rhs(m,i,iv) - rhs(q,n,iv) * per(2,m,q,i,iv)
          end do
          rhs(m,i,iv) = rhs(m,i,iv) / mat(3,m,m,i,iv)
        end do

!.... now do the rest of the rows in reverse order

        do i = n - 4, 1, -1
        
          p = i + 1
          r = i + 2
          
          do m = nblk, 1, -1
            do q = m+1, nblk
              rhs(m,i,iv) = rhs(m,i,iv) - rhs(q,i,iv) * mat(3,m,q,i,iv)
            end do
            do q = 1, nblk
              rhs(m,i,iv) = rhs(m,i,iv) - rhs(q,p,iv) * mat(4,m,q,i,iv)
              rhs(m,i,iv) = rhs(m,i,iv) - rhs(q,r,iv) * mat(5,m,q,i,iv)
              rhs(m,i,iv) = rhs(m,i,iv) - rhs(q,n-1,iv) * per(1,m,q,i,iv)
              rhs(m,i,iv) = rhs(m,i,iv) - rhs(q,n  ,iv) * per(2,m,q,i,iv)
            end do
            rhs(m,i,iv) = rhs(m,i,iv) / mat(3,m,m,i,iv)
          end do
                                  
        end do

!.... account for periodicity

        do q = 1, nblk
          rhs(q,np,iv) = rhs(q,1,iv)
        end do
        
        end do  ! loop over nsys

        end if
!=============================================================================!
        return
        end

     
        
