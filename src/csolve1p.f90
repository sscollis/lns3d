!=============================================================================!
        subroutine csolve1p( nsys, np, nblk, mat, rhs, per, per2, mult, fact)
!=============================================================================!
!  
!  Solves a block tridiagonal system of equations without pivoting.
!  Vectorized over the second index.
!
!  This routine solves a PERIODIC block tridiagonal system.
!
!  WARNING: Since pivoting is not done, this algorithm could fail
!
!  input:
!    nsys       : number of systems to vectorize over
!    np         : number of rows in the systems including periodic node
!    nblk       : dimension of the blocks in each row
!
!  inout:
!    mat        : the block tridiagonal matrix to be solved
!    rhs        : the right hand side
!=============================================================================!
        implicit none

        integer, intent(in) :: np, nblk, nsys
        complex, intent(inout) :: mat(3,nblk,nblk,nsys,np), rhs(nblk,nsys,np)
        complex, intent(inout) :: per(nblk,nblk,nsys,np)
        complex, intent(inout) :: per2(2,nblk,nblk,nsys)
        complex, intent(inout) :: mult(nsys), fact(nsys)
        
!.... useful constants

        real, parameter :: zero = 0.0, one = 1.0

!.... local variables
        
        integer :: nc, nr, i, l, m, n, p, q, r, iv
!=============================================================================!

c!.... don't solve for the redundant nodes

        n = np - 1

        !$omp parallel do private(nc,i,l,m,p,q,r,iv)
        do iv = 1, nsys

!.... setup the periodic fillin vector

        per(:,1:nblk,iv,1)  = mat(1,:,1:nblk,iv,1)
        per2(1,:,1:nblk,iv) = mat(3,:,1:nblk,iv,n)
        
!=============================================================================!
!                     F O R W A R D   E L I M I N A T I O N 
!=============================================================================!
        do i = 1, n - 3

!.... do the first four columns

          do nc = 1, nblk - 1
          
            do m = nc, nblk - 1
              l = m + 1
              mult(iv) = -mat(2,l,nc,iv,i) / mat(2,nc,nc,iv,i)
              do q = nc + 1, nblk
                mat(2,l,q,iv,i) = mat(2,l,q,iv,i) + mult(iv) * mat(2,nc,q,iv,i)
              end do
              do q = 1, nblk
                mat(3,l,q,iv,i) = mat(3,l,q,iv,i) + mult(iv) * mat(3,nc,q,iv,i)
                per(l,q,iv,i) = per(l,q,iv,i) + mult(iv) * per(nc,q,iv,i)
              end do
              rhs(l,iv,i) = rhs(l,iv,i) + mult(iv) * rhs(nc,iv,i)
            end do            

          end do        ! end of loop on nc

!.... for all columns
          
          nc = 1
        
          p = i + 1
  
          fact(iv) = one / mat(2,nc,nc,iv,i)
          do m = 1, nblk
            mult(iv) = -mat(1,m,nc,iv,p) * fact(iv)
            do q = nc + 1, nblk
              mat(1,m,q,iv,p) = mat(1,m,q,iv,p) + mult(iv) * mat(2,nc,q,iv,i)
            end do
            do q = 1, nblk
              mat(2,m,q,iv,p) = mat(2,m,q,iv,p) + mult(iv) * mat(3,nc,q,iv,i)
              per(m,q,iv,p) = mult(iv) * per(nc,q,iv,i)
            end do
            rhs(m,iv,p) = rhs(m,iv,p) + mult(iv) * rhs(nc,iv,i)

!.... the periodic stuff on the last row

            mult(iv) = -per2(1,m,nc,iv) * fact(iv)
            do q = nc + 1, nblk
              per2(1,m,q,iv) = per2(1,m,q,iv) + mult(iv) * mat(2,nc,q,iv,i)
            end do
            do q = 1, nblk
              per2(2,m,q,iv) = mult(iv) * mat(3,nc,q,iv,i)
              mat(2,m,q,iv,n) = mat(2,m,q,iv,n) + mult(iv) * per(nc,q,iv,i)
            end do
            rhs(m,iv,n) = rhs(m,iv,n) + mult(iv) * rhs(nc,iv,i)
          end do

          do nc = 2, nblk

            fact(iv) = one / mat(2,nc,nc,iv,i)
            do m = 1, nblk
              mult(iv) = -mat(1,m,nc,iv,p) * fact(iv)
              do q = nc + 1, nblk
                mat(1,m,q,iv,p) = mat(1,m,q,iv,p) + mult(iv) * mat(2,nc,q,iv,i)
              end do
              do q = 1, nblk
                mat(2,m,q,iv,p) = mat(2,m,q,iv,p) + mult(iv) * mat(3,nc,q,iv,i)
                per(m,q,iv,p) = per(m,q,iv,p) + mult(iv) * per(nc,q,iv,i)
              end do
              rhs(m,iv,p) = rhs(m,iv,p) + mult(iv) * rhs(nc,iv,i)

!.... the periodic stuff on the last row

              mult(iv) = -per2(1,m,nc,iv) * fact(iv)
              do q = nc + 1, nblk
                per2(1,m,q,iv) = per2(1,m,q,iv) + mult(iv) * mat(2,nc,q,iv,i)
              end do
              do q = 1, nblk
                per2(2,m,q,iv) = per2(2,m,q,iv) + mult(iv) * mat(3,nc,q,iv,i)
                mat(2,m,q,iv,n) = mat(2,m,q,iv,n) + mult(iv) * per(nc,q,iv,i)
              end do
              rhs(m,iv,n) = rhs(m,iv,n) + mult(iv) * rhs(nc,iv,i)
            end do

          end do        ! end of loop on nc
          
          per2(1,:,:,iv) = per2(2,:,:,iv)

        end do          ! end of loop on i

!.... do the third to the last row

        i = n - 2
        
!.... do the first four columns

        do nc = 1, nblk - 1
        
          do m = nc, nblk - 1
            l = m + 1
            mult(iv) = -mat(2,l,nc,iv,i) / mat(2,nc,nc,iv,i)
            do q = nc + 1, nblk
              mat(2,l,q,iv,i) = mat(2,l,q,iv,i) + mult(iv) * mat(2,nc,q,iv,i)
            end do
            do q = 1, nblk
              mat(3,l,q,iv,i) = mat(3,l,q,iv,i) + mult(iv) * mat(3,nc,q,iv,i)
              per(l,q,iv,i) = per(l,q,iv,i) + mult(iv) * per(nc,q,iv,i)
            end do
            rhs(l,iv,i) = rhs(l,iv,i) + mult(iv) * rhs(nc,iv,i)
          end do              

        end do  ! end of loop on nc

!.... for all columns
          
        do nc = 1, nblk

          p  = i + 1
          
          fact(iv) = one / mat(2,nc,nc,iv,i)
          do m = 1, nblk
            mult(iv) = -mat(1,m,nc,iv,p) * fact(iv)
            do q = nc + 1, nblk
              mat(1,m,q,iv,p) = mat(1,m,q,iv,p) + mult(iv) * mat(2,nc,q,iv,i)
            end do
            do q = 1, nblk
              mat(2,m,q,iv,p) = mat(2,m,q,iv,p) + mult(iv) * mat(3,nc,q,iv,i)
              mat(3,m,q,iv,p) = mat(3,m,q,iv,p) + mult(iv) * per(nc,q,iv,i)
            end do
            rhs(m,iv,p) = rhs(m,iv,p) + mult(iv) * rhs(nc,iv,i)

!.... the periodic stuff on the last row

            mult(iv) = -per2(1,m,nc,iv) * fact(iv)
            do q = nc + 1, nblk
              per2(1,m,q,iv) = per2(1,m,q,iv) + mult(iv) * mat(2,nc,q,iv,i)
            end do
            do q = 1, nblk
              mat(1,m,q,iv,n) = mat(1,m,q,iv,n) + mult(iv) * mat(3,nc,q,iv,i)
              mat(2,m,q,iv,n) = mat(2,m,q,iv,n) + mult(iv) * per(nc,q,iv,i)
            end do
            rhs(m,iv,n) = rhs(m,iv,n) + mult(iv) * rhs(nc,iv,i)
          end do

        end do  ! end of loop on nc
        
!.... do the second to the last row

        i = n - 1
        
!.... do the first four columns

          do nc = 1, nblk - 1
          
            do m = nc, nblk - 1
              l = m + 1
              mult(iv) = -mat(2,l,nc,iv,i) / mat(2,nc,nc,iv,i)
              do q = nc + 1, nblk
                mat(2,l,q,iv,i) = mat(2,l,q,iv,i) + mult(iv) * mat(2,nc,q,iv,i)
              end do
              do q = 1, nblk
                mat(3,l,q,iv,i) = mat(3,l,q,iv,i) + mult(iv) * mat(3,nc,q,iv,i)
              end do
              rhs(l,iv,i) = rhs(l,iv,i) + mult(iv) * rhs(nc,iv,i)
            end do            

          end do        ! end of loop on nc

!.... for all columns
          
          do nc = 1, nblk

            p  = i + 1
            
            fact(iv) = one / mat(2,nc,nc,iv,i)
            do m = 1, nblk
              mult(iv) = -mat(1,m,nc,iv,p) * fact(iv)
              do q = nc + 1, nblk
                mat(1,m,q,iv,p) = mat(1,m,q,iv,p) + mult(iv) * mat(2,nc,q,iv,i)
              end do
              do q = 1, nblk
                mat(2,m,q,iv,p) = mat(2,m,q,iv,p) + mult(iv) * mat(3,nc,q,iv,i)
              end do
              rhs(m,iv,p) = rhs(m,iv,p) + mult(iv) * rhs(nc,iv,i)
            end do

          end do        ! end of loop on nc
        
!.... do the last row
          
        i = n
        
        do nc = 1, nblk - 1
        
          do m = nc, nblk - 1
            l = m + 1
            mult(iv) = -mat(2,l,nc,iv,i) / mat(2,nc,nc,iv,i)
            do q = nc + 1, nblk
              mat(2,l,q,iv,i) = mat(2,l,q,iv,i) + mult(iv) * mat(2,nc,q,iv,i)
            end do
            rhs(l,iv,i) = rhs(l,iv,i) + mult(iv) * rhs(nc,iv,i)
          end do              

        end do  ! end of loop on nc

!=============================================================================!
!                   B A C K W A R D   S U B S T I T U T I O N 
!=============================================================================!

!.... do the last row first

        i = n
        
        do m = nblk, 1, -1
          do q = m+1, nblk
            rhs(m,iv,i) = rhs(m,iv,i) - rhs(q,iv,i) * mat(2,m,q,iv,i)
          end do
          rhs(m,iv,i) = rhs(m,iv,i) / mat(2,m,m,iv,i)
        end do

!.... do the second to the last row

        i = n - 1
        
        p = i + 1
        
        do m = nblk, 1, -1
          do q = m+1, nblk
            rhs(m,iv,i) = rhs(m,iv,i) - rhs(q,iv,i) * mat(2,m,q,iv,i)
          end do
          do q= 1, nblk
            rhs(m,iv,i) = rhs(m,iv,i) - rhs(q,iv,p) * mat(3,m,q,iv,i)
          end do
          rhs(m,iv,i) = rhs(m,iv,i) / mat(2,m,m,iv,i)
        end do

!.... now do the rest of the rows in reverse order

        do i = n - 2, 1, -1
        
          p = i + 1
          
          do m = nblk, 1, -1
            do q = m+1, nblk
              rhs(m,iv,i) = rhs(m,iv,i) - rhs(q,iv,i) * mat(2,m,q,iv,i)
            end do
            do q = 1, nblk
              rhs(m,iv,i) = rhs(m,iv,i) - rhs(q,iv,p) * mat(3,m,q,iv,i)
              rhs(m,iv,i) = rhs(m,iv,i) - rhs(q,iv,n) * per(m,q,iv,i)
            end do
            rhs(m,iv,i) = rhs(m,iv,i) / mat(2,m,m,iv,i)
          end do
                                  
        end do

!.... account for periodicity

        rhs(1,iv,np) = rhs(1,iv,1)
        rhs(2,iv,np) = rhs(2,iv,1)
        rhs(3,iv,np) = rhs(3,iv,1)
        rhs(4,iv,np) = rhs(4,iv,1)
        rhs(5,iv,np) = rhs(5,iv,1)

        end do

        return
        end
