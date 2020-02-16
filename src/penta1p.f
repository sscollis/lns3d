!=============================================================================!
        subroutine penta1p( nsys, np, nblk, mat, rhs, per, per2 )
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
!    np         : number of rows in the systems
!    nblk       : dimension of the blocks in each row
!
!  inout:
!    mat        : the block pentadiagonal matrix to be solved
!    rhs        : the right hand side
!=============================================================================!
        implicit none

        integer np, nblk, nsys
        real mat(5,nblk,nblk,nsys,np), rhs(nblk,nsys,np)
        real per(2,nblk,nblk,nsys,np), per2(6,nblk,nblk,nsys)
        real mult, fact
        
!.... useful constants

        real zero, one
        parameter (zero = 0.0, one = 1.0)

!.... local variables
        
        integer nc, i, l, m, n, p, q, r, iv
!=============================================================================!

!.... don't solve for the redundant nodes

        n = np - 1

!$omp parallel do private(nc,i,l,m,p,q,r,iv,fact,mult)
        do iv = 1, nsys

!.... setup the periodic fillin vector

        do q = 1, nblk
           do p = 1, nblk
              per(1,p,q,iv,1) = mat(1,p,q,iv,1)
              per(2,p,q,iv,1) = mat(2,p,q,iv,1)
              per(1,p,q,iv,2) = zero
              per(2,p,q,iv,2) = mat(1,p,q,iv,2)
              
              per2(1,p,q,iv)  = mat(5,p,q,iv,n-1)
              per2(2,p,q,iv)  = zero
              per2(4,p,q,iv)  = mat(4,p,q,iv,n)
              per2(5,p,q,iv)  = mat(5,p,q,iv,n)
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
                 mult = -mat(3,l,nc,iv,i) / mat(3,nc,nc,iv,i)
                 do q = nc + 1, nblk
                    mat(3,l,q,iv,i) = mat(3,l,q,iv,i) + mult * mat(3,nc,q,iv,i)
                 end do
                 do q = 1, nblk
                    mat(4,l,q,iv,i) = mat(4,l,q,iv,i) + mult * mat(4,nc,q,iv,i)
                    mat(5,l,q,iv,i) = mat(5,l,q,iv,i) + mult * mat(5,nc,q,iv,i)
                    per(1,l,q,iv,i) = per(1,l,q,iv,i) + mult * per(1,nc,q,iv,i)
                    per(2,l,q,iv,i) = per(2,l,q,iv,i) + mult * per(2,nc,q,iv,i)
                 end do
                 rhs(l,iv,i) = rhs(l,iv,i) + mult * rhs(nc,iv,i)
              end do          

           end do               ! end of loop on nc

!.... for all columns
          
           nc = 1
           
           p = i + 1
           r = i + 2
           
           fact = one / mat(3,nc,nc,iv,i)
           do m = 1, nblk
              mult = -mat(2,m,nc,iv,p) * fact
              do q = nc + 1, nblk
                 mat(2,m,q,iv,p) = mat(2,m,q,iv,p) + mult * mat(3,nc,q,iv,i)
              end do
              do q = 1, nblk
                 mat(3,m,q,iv,p) = mat(3,m,q,iv,p) + mult * mat(4,nc,q,iv,i)
                 mat(4,m,q,iv,p) = mat(4,m,q,iv,p) + mult * mat(5,nc,q,iv,i)
                 per(1,m,q,iv,p) = per(1,m,q,iv,p) + mult * per(1,nc,q,iv,i)
                 per(2,m,q,iv,p) = per(2,m,q,iv,p) + mult * per(2,nc,q,iv,i)
              end do
              rhs(m,iv,p) = rhs(m,iv,p) + mult * rhs(nc,iv,i)
              
              mult = -mat(1,m,nc,iv,r) * fact
              do q = nc + 1, nblk
                 mat(1,m,q,iv,r) = mat(1,m,q,iv,r) + mult * mat(3,nc,q,iv,i)
              end do
              do q = 1, nblk
                 mat(2,m,q,iv,r) = mat(2,m,q,iv,r) + mult * mat(4,nc,q,iv,i)
                 mat(3,m,q,iv,r) = mat(3,m,q,iv,r) + mult * mat(5,nc,q,iv,i)
                 per(1,m,q,iv,r) = mult * per(1,nc,q,iv,i)
                 per(2,m,q,iv,r) = mult * per(2,nc,q,iv,i)
              end do
              rhs(m,iv,r) = rhs(m,iv,r) + mult * rhs(nc,iv,i)

!.... the periodic stuff on the next-to-the-last row

              mult = -per2(1,m,nc,iv) * fact
              do q = nc + 1, nblk
                 per2(1,m,q,iv) = per2(1,m,q,iv) + mult * mat(3,nc,q,iv,i)
              end do
              do q = 1, nblk
                 per2(2,m,q,iv) = per2(2,m,q,iv) + mult * mat(4,nc,q,iv,i)
                 per2(3,m,q,iv) = mult * mat(5,nc,q,iv,i)
                 mat(3,m,q,iv,n-1) = mat(3,m,q,iv,n-1) + mult * per(1,nc,q,iv,i)
                 mat(4,m,q,iv,n-1) = mat(4,m,q,iv,n-1) + mult * per(2,nc,q,iv,i)
              end do
              rhs(m,iv,n-1) = rhs(m,iv,n-1) + mult * rhs(nc,iv,i)
              
!.... the periodic stuff on the last row

              mult = -per2(4,m,nc,iv) * fact
              do q = nc + 1, nblk
                 per2(4,m,q,iv) = per2(4,m,q,iv) + mult * mat(3,nc,q,iv,i)
              end do
              do q = 1, nblk
                 per2(5,m,q,iv) = per2(5,m,q,iv) + mult * mat(4,nc,q,iv,i)
                 per2(6,m,q,iv) = mult * mat(5,nc,q,iv,i)
                 mat(2,m,q,iv,n) = mat(2,m,q,iv,n) + mult * per(1,nc,q,iv,i)
                 mat(3,m,q,iv,n) = mat(3,m,q,iv,n) + mult * per(2,nc,q,iv,i)
              end do
              rhs(m,iv,n) = rhs(m,iv,n) + mult * rhs(nc,iv,i)

           end do               ! end of loop on m

           do nc = 2, nblk

              p = i + 1
              r = i + 2
              
              fact = one / mat(3,nc,nc,iv,i)
              do m = 1, nblk
                 mult = -mat(2,m,nc,iv,p) * fact
                 do q = nc + 1, nblk
                    mat(2,m,q,iv,p) = mat(2,m,q,iv,p) + mult * mat(3,nc,q,iv,i)
                 end do
                 do q = 1, nblk
                    mat(3,m,q,iv,p) = mat(3,m,q,iv,p) + mult * mat(4,nc,q,iv,i)
                    mat(4,m,q,iv,p) = mat(4,m,q,iv,p) + mult * mat(5,nc,q,iv,i)
                    per(1,m,q,iv,p) = per(1,m,q,iv,p) + mult * per(1,nc,q,iv,i)
                    per(2,m,q,iv,p) = per(2,m,q,iv,p) + mult * per(2,nc,q,iv,i)
                 end do
                 rhs(m,iv,p) = rhs(m,iv,p) + mult * rhs(nc,iv,i)
                 
                 mult = -mat(1,m,nc,iv,r) * fact
                 do q = nc + 1, nblk
                    mat(1,m,q,iv,r) = mat(1,m,q,iv,r) + mult * mat(3,nc,q,iv,i)
                 end do
                 do q = 1, nblk
                    mat(2,m,q,iv,r) = mat(2,m,q,iv,r) + mult * mat(4,nc,q,iv,i)
                    mat(3,m,q,iv,r) = mat(3,m,q,iv,r) + mult * mat(5,nc,q,iv,i)
                    per(1,m,q,iv,r) = per(1,m,q,iv,r) + mult * per(1,nc,q,iv,i)
                    per(2,m,q,iv,r) = per(2,m,q,iv,r) + mult * per(2,nc,q,iv,i)
                 end do
                 rhs(m,iv,r) = rhs(m,iv,r) + mult * rhs(nc,iv,i)

!.... the periodic stuff on the next-to-the-last row

                 mult = -per2(1,m,nc,iv) * fact
                 do q = nc + 1, nblk
                    per2(1,m,q,iv) = per2(1,m,q,iv) + mult * mat(3,nc,q,iv,i)
                 end do
                 do q = 1, nblk
                    per2(2,m,q,iv) = per2(2,m,q,iv) + mult * mat(4,nc,q,iv,i)
                    per2(3,m,q,iv) = per2(3,m,q,iv) + mult * mat(5,nc,q,iv,i)
                    mat(3,m,q,iv,n-1) = mat(3,m,q,iv,n-1) + mult * per(1,nc,q,iv,i)
                    mat(4,m,q,iv,n-1) = mat(4,m,q,iv,n-1) + mult * per(2,nc,q,iv,i)
                 end do
                 rhs(m,iv,n-1) = rhs(m,iv,n-1) + mult * rhs(nc,iv,i)
                 
!.... the periodic stuff on the last row

                 mult = -per2(4,m,nc,iv) * fact
                 do q = nc + 1, nblk
                    per2(4,m,q,iv) = per2(4,m,q,iv) + mult * mat(3,nc,q,iv,i)
                 end do
                 do q = 1, nblk
                    per2(5,m,q,iv) = per2(5,m,q,iv) + mult * mat(4,nc,q,iv,i)
                    per2(6,m,q,iv) = per2(6,m,q,iv) + mult * mat(5,nc,q,iv,i)
                    mat(2,m,q,iv,n) = mat(2,m,q,iv,n) + mult * per(1,nc,q,iv,i)
                    mat(3,m,q,iv,n) = mat(3,m,q,iv,n) + mult * per(2,nc,q,iv,i)
                 end do
                 rhs(m,iv,n) = rhs(m,iv,n) + mult * rhs(nc,iv,i)

              end do            ! end of loop on m

           end do               ! end of loop on nc
           
           do q = 1, nblk
              do p = 1, nblk
                 per2(1,p,q,iv)  = per2(2,p,q,iv)
                 per2(2,p,q,iv)  = per2(3,p,q,iv)
                 per2(4,p,q,iv)  = per2(5,p,q,iv)
                 per2(5,p,q,iv)  = per2(6,p,q,iv)
              end do
           end do

        end do                  ! end of loop on i

!====================================================================================
        i = n - 5

!.... do the first four columns

        do nc = 1, nblk - 1
           
           do m = nc, nblk - 1
              l = m + 1
              mult = -mat(3,l,nc,iv,i) / mat(3,nc,nc,iv,i)
              do q = nc + 1, nblk
                 mat(3,l,q,iv,i) = mat(3,l,q,iv,i) + mult * mat(3,nc,q,iv,i)
              end do
              do q = 1, nblk
                 mat(4,l,q,iv,i) = mat(4,l,q,iv,i) + mult * mat(4,nc,q,iv,i)
                 mat(5,l,q,iv,i) = mat(5,l,q,iv,i) + mult * mat(5,nc,q,iv,i)
                 per(1,l,q,iv,i) = per(1,l,q,iv,i) + mult * per(1,nc,q,iv,i)
                 per(2,l,q,iv,i) = per(2,l,q,iv,i) + mult * per(2,nc,q,iv,i)
              end do
              rhs(l,iv,i) = rhs(l,iv,i) + mult * rhs(nc,iv,i)
           end do             

        end do                  ! end of loop on nc

!.... for all columns
          
        nc = 1
        
        p = i + 1
        r = i + 2
        
        fact = one / mat(3,nc,nc,iv,i)
        do m = 1, nblk
           mult = -mat(2,m,nc,iv,p) * fact
           do q = nc + 1, nblk
              mat(2,m,q,iv,p) = mat(2,m,q,iv,p) + mult * mat(3,nc,q,iv,i)
           end do
           do q = 1, nblk
              mat(3,m,q,iv,p) = mat(3,m,q,iv,p) + mult * mat(4,nc,q,iv,i)
              mat(4,m,q,iv,p) = mat(4,m,q,iv,p) + mult * mat(5,nc,q,iv,i)
              per(1,m,q,iv,p) = per(1,m,q,iv,p) + mult * per(1,nc,q,iv,i)
              per(2,m,q,iv,p) = per(2,m,q,iv,p) + mult * per(2,nc,q,iv,i)
           end do
           rhs(m,iv,p) = rhs(m,iv,p) + mult * rhs(nc,iv,i)
           
           mult = -mat(1,m,nc,iv,r) * fact
           do q = nc + 1, nblk
              mat(1,m,q,iv,r) = mat(1,m,q,iv,r) + mult * mat(3,nc,q,iv,i)
           end do
           do q = 1, nblk
              mat(2,m,q,iv,r) = mat(2,m,q,iv,r) + mult * mat(4,nc,q,iv,i)
              mat(3,m,q,iv,r) = mat(3,m,q,iv,r) + mult * mat(5,nc,q,iv,i)
              mat(4,m,q,iv,r) = mat(4,m,q,iv,r) + mult * per(1,nc,q,iv,i)
              per(2,m,q,iv,r) = mult * per(2,nc,q,iv,i)
           end do
           rhs(m,iv,r) = rhs(m,iv,r) + mult * rhs(nc,iv,i)

!.... the periodic stuff on the next-to-the-last row

           mult = -per2(1,m,nc,iv) * fact
           do q = nc + 1, nblk
              per2(1,m,q,iv) = per2(1,m,q,iv) + mult * mat(3,nc,q,iv,i)
           end do
           do q = 1, nblk
              per2(2,m,q,iv) = per2(2,m,q,iv) + mult * mat(4,nc,q,iv,i)
              mat(1,m,q,iv,n-1) = mat(1,m,q,iv,n-1) + mult * mat(5,nc,q,iv,i)
              mat(3,m,q,iv,n-1) = mat(3,m,q,iv,n-1) + mult * per(1,nc,q,iv,i)
              mat(4,m,q,iv,n-1) = mat(4,m,q,iv,n-1) + mult * per(2,nc,q,iv,i)
           end do
           rhs(m,iv,n-1) = rhs(m,iv,n-1) + mult * rhs(nc,iv,i)
           
!.... the periodic stuff on the last row

           mult = -per2(4,m,nc,iv) * fact
           do q = nc + 1, nblk
              per2(4,m,q,iv) = per2(4,m,q,iv) + mult * mat(3,nc,q,iv,i)
           end do
           do q = 1, nblk
              per2(5,m,q,iv) = per2(5,m,q,iv) + mult * mat(4,nc,q,iv,i)
              per2(6,m,q,iv) = mult * mat(5,nc,q,iv,i)
              mat(2,m,q,iv,n) = mat(2,m,q,iv,n) + mult * per(1,nc,q,iv,i)
              mat(3,m,q,iv,n) = mat(3,m,q,iv,n) + mult * per(2,nc,q,iv,i)
           end do
           rhs(m,iv,n) = rhs(m,iv,n) + mult * rhs(nc,iv,i)

        end do                  ! end of loop on m

        do nc = 2, nblk

           p = i + 1
           r = i + 2
           
           fact = one / mat(3,nc,nc,iv,i)
           do m = 1, nblk
              mult = -mat(2,m,nc,iv,p) * fact
              do q = nc + 1, nblk
                 mat(2,m,q,iv,p) = mat(2,m,q,iv,p) + mult * mat(3,nc,q,iv,i)
              end do
              do q = 1, nblk
                 mat(3,m,q,iv,p) = mat(3,m,q,iv,p) + mult * mat(4,nc,q,iv,i)
                 mat(4,m,q,iv,p) = mat(4,m,q,iv,p) + mult * mat(5,nc,q,iv,i)
                 per(1,m,q,iv,p) = per(1,m,q,iv,p) + mult * per(1,nc,q,iv,i)
                 per(2,m,q,iv,p) = per(2,m,q,iv,p) + mult * per(2,nc,q,iv,i)
              end do
              rhs(m,iv,p) = rhs(m,iv,p) + mult * rhs(nc,iv,i)
              
              mult = -mat(1,m,nc,iv,r) * fact
              do q = nc + 1, nblk
                 mat(1,m,q,iv,r) = mat(1,m,q,iv,r) + mult * mat(3,nc,q,iv,i)
              end do
              do q = 1, nblk
                 mat(2,m,q,iv,r) = mat(2,m,q,iv,r) + mult * mat(4,nc,q,iv,i)
                 mat(3,m,q,iv,r) = mat(3,m,q,iv,r) + mult * mat(5,nc,q,iv,i)
                 mat(4,m,q,iv,r) = mat(4,m,q,iv,r) + mult * per(1,nc,q,iv,i)
                 per(2,m,q,iv,r) = per(2,m,q,iv,r) + mult * per(2,nc,q,iv,i)
              end do
              rhs(m,iv,r) = rhs(m,iv,r) + mult * rhs(nc,iv,i)

!.... the periodic stuff on the next-to-the-last row

              mult = -per2(1,m,nc,iv) * fact
              do q = nc + 1, nblk
                 per2(1,m,q,iv) = per2(1,m,q,iv) + mult * mat(3,nc,q,iv,i)
              end do
              do q = 1, nblk
                 per2(2,m,q,iv) = per2(2,m,q,iv) + mult * mat(4,nc,q,iv,i)
                 mat(1,m,q,iv,n-1) = mat(1,m,q,iv,n-1) + mult * mat(5,nc,q,iv,i)
                 mat(3,m,q,iv,n-1) = mat(3,m,q,iv,n-1) + mult * per(1,nc,q,iv,i)
                 mat(4,m,q,iv,n-1) = mat(4,m,q,iv,n-1) + mult * per(2,nc,q,iv,i)
              end do
              rhs(m,iv,n-1) = rhs(m,iv,n-1) + mult * rhs(nc,iv,i)
              
!.... the periodic stuff on the last row

              mult = -per2(4,m,nc,iv) * fact
              do q = nc + 1, nblk
                 per2(4,m,q,iv) = per2(4,m,q,iv) + mult * mat(3,nc,q,iv,i)
              end do
              do q = 1, nblk
                 per2(5,m,q,iv) = per2(5,m,q,iv) + mult * mat(4,nc,q,iv,i)
                 per2(6,m,q,iv) = per2(6,m,q,iv) + mult * mat(5,nc,q,iv,i)
                 mat(2,m,q,iv,n) = mat(2,m,q,iv,n) + mult * per(1,nc,q,iv,i)
                 mat(3,m,q,iv,n) = mat(3,m,q,iv,n) + mult * per(2,nc,q,iv,i)
              end do
              rhs(m,iv,n) = rhs(m,iv,n) + mult * rhs(nc,iv,i)

           end do               ! end of loop on m

        end do                  ! end of loop on nc
        
        do q = 1, nblk
           do p = 1, nblk
              per2(1,p,q,iv)  = per2(2,p,q,iv)
              per2(4,p,q,iv)  = per2(5,p,q,iv)
              per2(5,p,q,iv)  = per2(6,p,q,iv)
           end do
        end do
        
!====================================================================================
        i = n - 4

!.... do the first four columns

        do nc = 1, nblk - 1
           
           do m = nc, nblk - 1
              l = m + 1
              mult = -mat(3,l,nc,iv,i) / mat(3,nc,nc,iv,i)
              do q = nc + 1, nblk
                 mat(3,l,q,iv,i) = mat(3,l,q,iv,i) + mult * mat(3,nc,q,iv,i)
              end do
              do q = 1, nblk
                 mat(4,l,q,iv,i) = mat(4,l,q,iv,i) + mult * mat(4,nc,q,iv,i)
                 mat(5,l,q,iv,i) = mat(5,l,q,iv,i) + mult * mat(5,nc,q,iv,i)
                 per(1,l,q,iv,i) = per(1,l,q,iv,i) + mult * per(1,nc,q,iv,i)
                 per(2,l,q,iv,i) = per(2,l,q,iv,i) + mult * per(2,nc,q,iv,i)
              end do
              rhs(l,iv,i) = rhs(l,iv,i) + mult * rhs(nc,iv,i)
           end do             

        end do                  ! end of loop on nc

!.... for all columns
          
        do nc = 1, nblk

           p = i + 1
           r = i + 2
           
           fact = one / mat(3,nc,nc,iv,i)
           do m = 1, nblk
              mult = -mat(2,m,nc,iv,p) * fact
              do q = nc + 1, nblk
                 mat(2,m,q,iv,p) = mat(2,m,q,iv,p) + mult * mat(3,nc,q,iv,i)
              end do
              do q = 1, nblk
                 mat(3,m,q,iv,p) = mat(3,m,q,iv,p) + mult * mat(4,nc,q,iv,i)
                 mat(4,m,q,iv,p) = mat(4,m,q,iv,p) + mult * mat(5,nc,q,iv,i)
                 mat(5,m,q,iv,p) = mat(5,m,q,iv,p) + mult * per(1,nc,q,iv,i)
                 per(2,m,q,iv,p) = per(2,m,q,iv,p) + mult * per(2,nc,q,iv,i)
              end do
              rhs(m,iv,p) = rhs(m,iv,p) + mult * rhs(nc,iv,i)
              
              mult = -mat(1,m,nc,iv,r) * fact
              do q = nc + 1, nblk
                 mat(1,m,q,iv,r) = mat(1,m,q,iv,r) + mult * mat(3,nc,q,iv,i)
              end do
              do q = 1, nblk
                 mat(2,m,q,iv,r) = mat(2,m,q,iv,r) + mult * mat(4,nc,q,iv,i)
                 mat(3,m,q,iv,r) = mat(3,m,q,iv,r) + mult * mat(5,nc,q,iv,i)
                 mat(4,m,q,iv,r) = mat(4,m,q,iv,r) + mult * per(1,nc,q,iv,i)
                 mat(5,m,q,iv,r) = mat(5,m,q,iv,r) + mult * per(2,nc,q,iv,i)
              end do
              rhs(m,iv,r) = rhs(m,iv,r) + mult * rhs(nc,iv,i)

!.... the periodic stuff on the next-to-the-last row

              mult = -per2(1,m,nc,iv) * fact
              do q = nc + 1, nblk
                 per2(1,m,q,iv) = per2(1,m,q,iv) + mult * mat(3,nc,q,iv,i)
              end do
              do q = 1, nblk
                 mat(1,m,q,iv,n-1) = mat(1,m,q,iv,n-1) + mult * mat(4,nc,q,iv,i)
                 mat(2,m,q,iv,n-1) = mat(2,m,q,iv,n-1) + mult * mat(5,nc,q,iv,i)
                 mat(3,m,q,iv,n-1) = mat(3,m,q,iv,n-1) + mult * per(1,nc,q,iv,i)
                 mat(4,m,q,iv,n-1) = mat(4,m,q,iv,n-1) + mult * per(2,nc,q,iv,i)
              end do
              rhs(m,iv,n-1) = rhs(m,iv,n-1) + mult * rhs(nc,iv,i)
              
!.... the periodic stuff on the last row

              mult = -per2(4,m,nc,iv) * fact
              do q = nc + 1, nblk
                 per2(4,m,q,iv) = per2(4,m,q,iv) + mult * mat(3,nc,q,iv,i)
              end do
              do q = 1, nblk
                 per2(5,m,q,iv) = per2(5,m,q,iv) + mult * mat(4,nc,q,iv,i)
                 mat(1,m,q,iv,n) = mat(1,m,q,iv,n) + mult * mat(5,nc,q,iv,i)
                 mat(2,m,q,iv,n) = mat(2,m,q,iv,n) + mult * per(1,nc,q,iv,i)
                 mat(3,m,q,iv,n) = mat(3,m,q,iv,n) + mult * per(2,nc,q,iv,i)
              end do
              rhs(m,iv,n) = rhs(m,iv,n) + mult * rhs(nc,iv,i)

           end do               ! end of loop on m

        end do                  ! end of loop on nc
        
        do q = 1, nblk
           do p = 1, nblk
              per2(4,p,q,iv)  = per2(5,p,q,iv)
           end do
        end do

!====================================================================================
        i = n - 3

!.... do the first four columns

        do nc = 1, nblk - 1
           
           do m = nc, nblk - 1
              l = m + 1
              mult = -mat(3,l,nc,iv,i) / mat(3,nc,nc,iv,i)
              do q = nc + 1, nblk
                 mat(3,l,q,iv,i) = mat(3,l,q,iv,i) + mult * mat(3,nc,q,iv,i)
              end do
              do q = 1, nblk
                 mat(4,l,q,iv,i) = mat(4,l,q,iv,i) + mult * mat(4,nc,q,iv,i)
                 mat(5,l,q,iv,i) = mat(5,l,q,iv,i) + mult * mat(5,nc,q,iv,i)
                 per(2,l,q,iv,i) = per(2,l,q,iv,i) + mult * per(2,nc,q,iv,i)
              end do
              rhs(l,iv,i) = rhs(l,iv,i) + mult * rhs(nc,iv,i)
           end do             

        end do                  ! end of loop on nc

!.... for all columns
          
        do nc = 1, nblk

           p = i + 1
           r = i + 2
           
           fact = one / mat(3,nc,nc,iv,i)
           do m = 1, nblk
              mult = -mat(2,m,nc,iv,p) * fact
              do q = nc + 1, nblk
                 mat(2,m,q,iv,p) = mat(2,m,q,iv,p) + mult * mat(3,nc,q,iv,i)
              end do
              do q = 1, nblk
                 mat(3,m,q,iv,p) = mat(3,m,q,iv,p) + mult * mat(4,nc,q,iv,i)
                 mat(4,m,q,iv,p) = mat(4,m,q,iv,p) + mult * mat(5,nc,q,iv,i)
                 mat(5,m,q,iv,p) = mat(5,m,q,iv,p) + mult * per(2,nc,q,iv,i)
              end do
              rhs(m,iv,p) = rhs(m,iv,p) + mult * rhs(nc,iv,i)
              
              mult = -mat(1,m,nc,iv,r) * fact
              do q = nc + 1, nblk
                 mat(1,m,q,iv,r) = mat(1,m,q,iv,r) + mult * mat(3,nc,q,iv,i)
              end do
              do q = 1, nblk
                 mat(2,m,q,iv,r) = mat(2,m,q,iv,r) + mult * mat(4,nc,q,iv,i)
                 mat(3,m,q,iv,r) = mat(3,m,q,iv,r) + mult * mat(5,nc,q,iv,i)
                 mat(4,m,q,iv,r) = mat(4,m,q,iv,r) + mult * per(2,nc,q,iv,i)
              end do
              rhs(m,iv,r) = rhs(m,iv,r) + mult * rhs(nc,iv,i)

!.... the periodic stuff on the last row

              mult = -per2(4,m,nc,iv) * fact
              do q = nc + 1, nblk
                 per2(4,m,q,iv) = per2(4,m,q,iv) + mult * mat(3,nc,q,iv,i)
              end do
              do q = 1, nblk
                 mat(1,m,q,iv,n) = mat(1,m,q,iv,n) + mult * mat(4,nc,q,iv,i)
                 mat(2,m,q,iv,n) = mat(2,m,q,iv,n) + mult * mat(5,nc,q,iv,i)
                 mat(3,m,q,iv,n) = mat(3,m,q,iv,n) + mult * per(2,nc,q,iv,i)
              end do
              rhs(m,iv,n) = rhs(m,iv,n) + mult * rhs(nc,iv,i)

           end do               ! end of loop on m

        end do                  ! end of loop on nc
        
!====================================================================================
        i = n - 2

!.... do the first four columns

        do nc = 1, nblk - 1
           
           do m = nc, nblk - 1
              l = m + 1
              mult = -mat(3,l,nc,iv,i) / mat(3,nc,nc,iv,i)
              do q = nc + 1, nblk
                 mat(3,l,q,iv,i) = mat(3,l,q,iv,i) + mult * mat(3,nc,q,iv,i)
              end do
              do q = 1, nblk
                 mat(4,l,q,iv,i) = mat(4,l,q,iv,i) + mult * mat(4,nc,q,iv,i)
                 mat(5,l,q,iv,i) = mat(5,l,q,iv,i) + mult * mat(5,nc,q,iv,i)
              end do
              rhs(l,iv,i) = rhs(l,iv,i) + mult * rhs(nc,iv,i)
           end do             

        end do                  ! end of loop on nc

!.... for all columns
          
        do nc = 1, nblk

           p = i + 1
           r = i + 2
           
           fact = one / mat(3,nc,nc,iv,i)
           do m = 1, nblk
              mult = -mat(2,m,nc,iv,p) * fact
              do q = nc + 1, nblk
                 mat(2,m,q,iv,p) = mat(2,m,q,iv,p) + mult * mat(3,nc,q,iv,i)
              end do
              do q = 1, nblk
                 mat(3,m,q,iv,p) = mat(3,m,q,iv,p) + mult * mat(4,nc,q,iv,i)
                 mat(4,m,q,iv,p) = mat(4,m,q,iv,p) + mult * mat(5,nc,q,iv,i)
              end do
              rhs(m,iv,p) = rhs(m,iv,p) + mult * rhs(nc,iv,i)
              
              mult = -mat(1,m,nc,iv,r) * fact
              do q = nc + 1, nblk
                 mat(1,m,q,iv,r) = mat(1,m,q,iv,r) + mult * mat(3,nc,q,iv,i)
              end do
              do q = 1, nblk
                 mat(2,m,q,iv,r) = mat(2,m,q,iv,r) + mult * mat(4,nc,q,iv,i)
                 mat(3,m,q,iv,r) = mat(3,m,q,iv,r) + mult * mat(5,nc,q,iv,i)
              end do
              rhs(m,iv,r) = rhs(m,iv,r) + mult * rhs(nc,iv,i)

           end do               ! end of loop on m

        end do                  ! end of loop on nc

!====================================================================================
        i = n - 1

!.... do the first four columns

        do nc = 1, nblk - 1
           
           do m = nc, nblk - 1
              l = m + 1
              mult = -mat(3,l,nc,iv,i) / mat(3,nc,nc,iv,i)
              do q = nc + 1, nblk
                 mat(3,l,q,iv,i) = mat(3,l,q,iv,i) + mult * mat(3,nc,q,iv,i)
              end do
              do q = 1, nblk
                 mat(4,l,q,iv,i) = mat(4,l,q,iv,i) + mult * mat(4,nc,q,iv,i)
              end do
              rhs(l,iv,i) = rhs(l,iv,i) + mult * rhs(nc,iv,i)
           end do             

        end do                  ! end of loop on nc

!.... for all columns
          
        do nc = 1, nblk

           p = i + 1
           
           fact = one / mat(3,nc,nc,iv,i)
           
           do m = 1, nblk
              mult = -mat(2,m,nc,iv,p) * fact
              do q = nc + 1, nblk
                 mat(2,m,q,iv,p) = mat(2,m,q,iv,p) + mult * mat(3,nc,q,iv,i)
              end do
              do q = 1, nblk
                 mat(3,m,q,iv,p) = mat(3,m,q,iv,p) + mult * mat(4,nc,q,iv,i)
              end do
              rhs(m,iv,p) = rhs(m,iv,p) + mult * rhs(nc,iv,i)
           end do

        end do                  ! end of loop on nc

!.... do the last row

        i = n
        
        do nc = 1, nblk - 1
           
           do m = nc, nblk - 1
              l = m + 1
              mult = -mat(3,l,nc,iv,i) / mat(3,nc,nc,iv,i)
              do q = nc + 1, nblk
                 mat(3,l,q,iv,i) = mat(3,l,q,iv,i) + mult * mat(3,nc,q,iv,i)
              end do
              rhs(l,iv,i) = rhs(l,iv,i) + mult * rhs(nc,iv,i)
           end do             

        end do                  ! end of loop on nc

!=============================================================================!
!                   B A C K W A R D   S U B S T I T U T I O N 
!=============================================================================!

!.... do the last row first

        i = n
        
        do m = nblk, 1, -1
           do q = m+1, nblk
              rhs(m,iv,i) = rhs(m,iv,i) - rhs(q,iv,i) * mat(3,m,q,iv,i)
           end do
           rhs(m,iv,i) = rhs(m,iv,i) / mat(3,m,m,iv,i)
        end do

!.... do the next-to-the-last row

        i = n - 1
        p = i + 1
        
        do m = nblk, 1, -1
           do q = m+1, nblk
              rhs(m,iv,i) = rhs(m,iv,i) - rhs(q,iv,i) * mat(3,m,q,iv,i)
           end do
           do q = 1, nblk
              rhs(m,iv,i) = rhs(m,iv,i) - rhs(q,iv,p) * mat(4,m,q,iv,i)
           end do
           rhs(m,iv,i) = rhs(m,iv,i) / mat(3,m,m,iv,i)
        end do

        i = n - 2
        p = i + 1
        r = i + 2
        
        do m = nblk, 1, -1
           do q = m+1, nblk
              rhs(m,iv,i) = rhs(m,iv,i) - rhs(q,iv,i) * mat(3,m,q,iv,i)
           end do
           do q = 1, nblk
              rhs(m,iv,i) = rhs(m,iv,i) - rhs(q,iv,p) * mat(4,m,q,iv,i)
              rhs(m,iv,i) = rhs(m,iv,i) - rhs(q,iv,r) * mat(5,m,q,iv,i)
           end do
           rhs(m,iv,i) = rhs(m,iv,i) / mat(3,m,m,iv,i)
        end do

        i = n - 3
        p = i + 1
        r = i + 2
        
        do m = nblk, 1, -1
           do q = m+1, nblk
              rhs(m,iv,i) = rhs(m,iv,i) - rhs(q,iv,i) * mat(3,m,q,iv,i)
           end do
           do q = 1, nblk
              rhs(m,iv,i) = rhs(m,iv,i) - rhs(q,iv,p) * mat(4,m,q,iv,i)
              rhs(m,iv,i) = rhs(m,iv,i) - rhs(q,iv,r) * mat(5,m,q,iv,i)
              rhs(m,iv,i) = rhs(m,iv,i) - rhs(q,iv,n) * per(2,m,q,iv,i)
           end do
           rhs(m,iv,i) = rhs(m,iv,i) / mat(3,m,m,iv,i)
        end do

!.... now do the rest of the rows in reverse order

        do i = n - 4, 1, -1
           
           p = i + 1
           r = i + 2
           
           do m = nblk, 1, -1
              do q = m+1, nblk
                 rhs(m,iv,i) = rhs(m,iv,i) - rhs(q,iv,i) * mat(3,m,q,iv,i)
              end do
              do q = 1, nblk
                 rhs(m,iv,i) = rhs(m,iv,i) - rhs(q,iv,p) * mat(4,m,q,iv,i)
                 rhs(m,iv,i) = rhs(m,iv,i) - rhs(q,iv,r) * mat(5,m,q,iv,i)
                 rhs(m,iv,i) = rhs(m,iv,i) - rhs(q,iv,n-1) * per(1,m,q,iv,i)
                 rhs(m,iv,i) = rhs(m,iv,i) - rhs(q,iv,n  ) * per(2,m,q,iv,i)
              end do
              rhs(m,iv,i) = rhs(m,iv,i) / mat(3,m,m,iv,i)
           end do
           
        end do

!.... account for periodicity

        do q = 1, nblk
           rhs(q,iv,np) = rhs(q,iv,1)
        end do

        end do

        return
        end
