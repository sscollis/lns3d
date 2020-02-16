!=============================================================================!
        subroutine solve2p( np, nsys, nblk, mat, rhs, per, per2 )
!=============================================================================!
!  
!  Solves a block tridiagonal system of equations without pivoting.
!  Vectorized over the first index.
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
        real, intent(inout) :: mat(3,nblk,nblk,np,nsys), rhs(nblk,np,nsys)
        real, intent(inout) :: per(nblk,nblk,np,nsys), per2(2,nblk,nblk,nsys)

        real :: mult, fact
        
!.... useful constants

        real, parameter :: zero = 0.0, one = 1.0

!.... local variables
        
        integer :: nc, nr, i, l, m, n, p, q, r, iv
!=============================================================================!

!.... don't solve for the redundant nodes

        n = np - 1

        !$omp parallel do private(nc,i,l,m,p,q,r,iv,fact,mult)
        do iv = 1, nsys

!.... setup the periodic fillin vector

        per(:,1:nblk,1,iv) = mat(1,:,1:nblk,1,iv)
        per2(1,:,1:nblk,iv)  = mat(3,:,1:nblk,n,iv)
        
!=============================================================================!
!                     F O R W A R D   E L I M I N A T I O N 
!=============================================================================!
        do i = 1, n - 3

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
                per(l,q,i,iv) = per(l,q,i,iv) + mult * per(nc,q,i,iv)
              end do
              rhs(l,i,iv) = rhs(l,i,iv) + mult * rhs(nc,i,iv)
            end do            

          end do        ! end of loop on nc

!.... for all columns
          
          nc = 1
        
          p = i + 1
  
          fact = one / mat(2,nc,nc,i,iv)
          do m = 1, nblk
            mult = -mat(1,m,nc,p,iv) * fact
            do q = nc + 1, nblk
              mat(1,m,q,p,iv) = mat(1,m,q,p,iv) + mult * mat(2,nc,q,i,iv)
            end do
            do q = 1, nblk
              mat(2,m,q,p,iv) = mat(2,m,q,p,iv) + mult * mat(3,nc,q,i,iv)
              per(m,q,p,iv) = mult * per(nc,q,i,iv)
            end do
            rhs(m,p,iv) = rhs(m,p,iv) + mult * rhs(nc,i,iv)

!.... the periodic stuff on the last row

            mult = -per2(1,m,nc,iv) * fact
            do q = nc + 1, nblk
              per2(1,m,q,iv) = per2(1,m,q,iv) + mult * mat(2,nc,q,i,iv)
            end do
            do q = 1, nblk
              per2(2,m,q,iv) = mult * mat(3,nc,q,i,iv)
              mat(2,m,q,n,iv) = mat(2,m,q,n,iv) + mult * per(nc,q,i,iv)
            end do
            rhs(m,n,iv) = rhs(m,n,iv) + mult * rhs(nc,i,iv)
          end do

          do nc = 2, nblk

            fact = one / mat(2,nc,nc,i,iv)
            do m = 1, nblk
              mult = -mat(1,m,nc,p,iv) * fact
              do q = nc + 1, nblk
                mat(1,m,q,p,iv) = mat(1,m,q,p,iv) + mult * mat(2,nc,q,i,iv)
              end do
              do q = 1, nblk
                mat(2,m,q,p,iv) = mat(2,m,q,p,iv) + mult * mat(3,nc,q,i,iv)
                per(m,q,p,iv) = per(m,q,p,iv) + mult * per(nc,q,i,iv)
              end do
              rhs(m,p,iv) = rhs(m,p,iv) + mult * rhs(nc,i,iv)

!.... the periodic stuff on the last row

              mult = -per2(1,m,nc,iv) * fact
              do q = nc + 1, nblk
                per2(1,m,q,iv) = per2(1,m,q,iv) + mult * mat(2,nc,q,i,iv)
              end do
              do q = 1, nblk
                per2(2,m,q,iv) = per2(2,m,q,iv) + mult * mat(3,nc,q,i,iv)
                mat(2,m,q,n,iv) = mat(2,m,q,n,iv) + mult * per(nc,q,i,iv)
              end do
              rhs(m,n,iv) = rhs(m,n,iv) + mult * rhs(nc,i,iv)
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
            mult = -mat(2,l,nc,i,iv) / mat(2,nc,nc,i,iv)
            do q = nc + 1, nblk
              mat(2,l,q,i,iv) = mat(2,l,q,i,iv) + mult * mat(2,nc,q,i,iv)
            end do
            do q = 1, nblk
              mat(3,l,q,i,iv) = mat(3,l,q,i,iv) + mult * mat(3,nc,q,i,iv)
              per(l,q,i,iv) = per(l,q,i,iv) + mult * per(nc,q,i,iv)
            end do
            rhs(l,i,iv) = rhs(l,i,iv) + mult * rhs(nc,i,iv)
          end do              

        end do  ! end of loop on nc

!.... for all columns
          
        do nc = 1, nblk

          p  = i + 1
          
          fact = one / mat(2,nc,nc,i,iv)
          do m = 1, nblk
            mult = -mat(1,m,nc,p,iv) * fact
            do q = nc + 1, nblk
              mat(1,m,q,p,iv) = mat(1,m,q,p,iv) + mult * mat(2,nc,q,i,iv)
            end do
            do q = 1, nblk
              mat(2,m,q,p,iv) = mat(2,m,q,p,iv) + mult * mat(3,nc,q,i,iv)
              mat(3,m,q,p,iv) = mat(3,m,q,p,iv) + mult * per(nc,q,i,iv)
            end do
            rhs(m,p,iv) = rhs(m,p,iv) + mult * rhs(nc,i,iv)

!.... the periodic stuff on the last row

            mult = -per2(1,m,nc,iv) * fact
            do q = nc + 1, nblk
              per2(1,m,q,iv) = per2(1,m,q,iv) + mult * mat(2,nc,q,i,iv)
            end do
            do q = 1, nblk
              mat(1,m,q,n,iv) = mat(1,m,q,n,iv) + mult * mat(3,nc,q,i,iv)
              mat(2,m,q,n,iv) = mat(2,m,q,n,iv) + mult * per(nc,q,i,iv)
            end do
            rhs(m,n,iv) = rhs(m,n,iv) + mult * rhs(nc,i,iv)
          end do

        end do  ! end of loop on nc
        
!.... do the second to the last row

        i = n - 1
        
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
          
          do nc = 1, nblk

            p  = i + 1
            
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

!.... do the last row first

        i = n
        
        do m = nblk, 1, -1
          do q = m+1, nblk
            rhs(m,i,iv) = rhs(m,i,iv) - rhs(q,i,iv) * mat(2,m,q,i,iv)
          end do
          rhs(m,i,iv) = rhs(m,i,iv) / mat(2,m,m,i,iv)
        end do

!.... do the second to the last row

        i = n - 1
        
        p = i + 1
        
        do m = nblk, 1, -1
          do q = m+1, nblk
            rhs(m,i,iv) = rhs(m,i,iv) - rhs(q,i,iv) * mat(2,m,q,i,iv)
          end do
          do q= 1, nblk
            rhs(m,i,iv) = rhs(m,i,iv) - rhs(q,p,iv) * mat(3,m,q,i,iv)
          end do
          rhs(m,i,iv) = rhs(m,i,iv) / mat(2,m,m,i,iv)
        end do

!.... now do the rest of the rows in reverse order

        do i = n - 2, 1, -1
        
          p = i + 1
          
          do m = nblk, 1, -1
            do q = m+1, nblk
              rhs(m,i,iv) = rhs(m,i,iv) - rhs(q,i,iv) * mat(2,m,q,i,iv)
            end do
            do q = 1, nblk
              rhs(m,i,iv) = rhs(m,i,iv) - rhs(q,p,iv) * mat(3,m,q,i,iv)
              rhs(m,i,iv) = rhs(m,i,iv) - rhs(q,n,iv) * per(m,q,i,iv)
            end do
            rhs(m,i,iv) = rhs(m,i,iv) / mat(2,m,m,i,iv)
          end do
                                  
        end do

!.... account for periodicity

        rhs(1,np,iv) = rhs(1,1,iv)
        rhs(2,np,iv) = rhs(2,1,iv)
        rhs(3,np,iv) = rhs(3,1,iv)
        rhs(4,np,iv) = rhs(4,1,iv)
        rhs(5,np,iv) = rhs(5,1,iv)

        end do

!=============================================================================!
        return
        end
