!=============================================================================!
        subroutine csolve1p( nsys, np, nblk, mat, rhs, per, per2, mult, fact)
!=============================================================================!
!  
!  Solves a block tridiagonal system of equations without pivoting.
!  Vectorized over the first index.
!
!  This routine solves a PERIODIC block tridiagonal system.
!
!  WARNING: Since pivoting is not done, this algorithm could fail
!
!  WARNING: This routine works for nblk = 5 only!
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
        complex, intent(inout) :: mat(nsys,np,nblk,nblk*3), rhs(nsys,np,nblk)
        complex, intent(inout) :: per(nsys,np,nblk,nblk),per2(nsys,nblk,2*nblk)
        complex, intent(inout) :: mult(nsys), fact(nsys)
        
!.... useful constants

        real, parameter :: zero = 0.0, one = 1.0

!.... local variables
        
        integer :: nc, nr, i, j, l, m, n, p, q, r
!=============================================================================!

!.... don't solve for the redundant nodes

        n = np - 1

!.... setup the periodic fillin vector

        per(:,1,:,1:nblk) = mat(:,1,:,1:nblk)
        per2(:,:,1:nblk)  = mat(:,n,:,2*nblk+1:3*nblk)
        
!=============================================================================!
!                     F O R W A R D   E L I M I N A T I O N 
!=============================================================================!
        do i = 1, n - 3

!.... do the first four columns

          do nc = 1, nblk - 1
          
            j = nblk + nc

            do m = nc, nblk - 1
              l = m + 1
              mult(:) = -mat(:,i,l,j) / mat(:,i,nc,j)

              do q = nc + nblk + 1, 2 * nblk
                mat(:,i,l,q) = mat(:,i,l,q) + mult(:) * mat(:,i,nc,q)
              end do
              
              mat(:,i,l,11) = mat(:,i,l,11) + mult(:) * mat(:,i,nc,11)
              mat(:,i,l,12) = mat(:,i,l,12) + mult(:) * mat(:,i,nc,12)
              mat(:,i,l,13) = mat(:,i,l,13) + mult(:) * mat(:,i,nc,13)
              mat(:,i,l,14) = mat(:,i,l,14) + mult(:) * mat(:,i,nc,14)
              mat(:,i,l,15) = mat(:,i,l,15) + mult(:) * mat(:,i,nc,15)

              per(:,i,l,1) = per(:,i,l,1) + mult(:) * per(:,i,nc,1)
              per(:,i,l,2) = per(:,i,l,2) + mult(:) * per(:,i,nc,2)
              per(:,i,l,3) = per(:,i,l,3) + mult(:) * per(:,i,nc,3)
              per(:,i,l,4) = per(:,i,l,4) + mult(:) * per(:,i,nc,4)
              per(:,i,l,5) = per(:,i,l,5) + mult(:) * per(:,i,nc,5)

              rhs(:,i,l) = rhs(:,i,l) + mult(:) * rhs(:,i,nc)
            end do            

          end do        ! end of loop on nc

!.... for all columns
          
          p = i + 1
  
          nc = 1
        
          j = nc + nblk
          fact(:) = one / mat(:,i,nc,j)
          
          do m = 1, nblk

            mult(:) = -mat(:,p,m,nc) * fact(:)
            do q = nc + 1, nblk
              l = q + nblk
              mat(:,p,m,q) = mat(:,p,m,q) + mult(:) * mat(:,i,nc,l)
            end do

            mat(:,p,m, 6) = mat(:,p,m, 6) + mult(:) * mat(:,i,nc,11)
            mat(:,p,m, 7) = mat(:,p,m, 7) + mult(:) * mat(:,i,nc,12)
            mat(:,p,m, 8) = mat(:,p,m, 8) + mult(:) * mat(:,i,nc,13)
            mat(:,p,m, 9) = mat(:,p,m, 9) + mult(:) * mat(:,i,nc,14)
            mat(:,p,m,10) = mat(:,p,m,10) + mult(:) * mat(:,i,nc,15)

            per(:,p,m,1) = mult(:) * per(:,i,nc,1)
            per(:,p,m,2) = mult(:) * per(:,i,nc,2)
            per(:,p,m,3) = mult(:) * per(:,i,nc,3)
            per(:,p,m,4) = mult(:) * per(:,i,nc,4)
            per(:,p,m,5) = mult(:) * per(:,i,nc,5)

            rhs(:,p,m) = rhs(:,p,m) + mult(:) * rhs(:,i,nc)

!.... the periodic stuff on the last row

            mult(:) = -per2(:,m,nc) * fact(:)
            do q = nc + 1, nblk
              l = q + nblk
              per2(:,m,q) = per2(:,m,q) + mult(:) * mat(:,i,nc,l)
            end do

            per2(:,m, 6) = mult(:) * mat(:,i,nc,11)
            per2(:,m, 7) = mult(:) * mat(:,i,nc,12)
            per2(:,m, 8) = mult(:) * mat(:,i,nc,13)
            per2(:,m, 9) = mult(:) * mat(:,i,nc,14)
            per2(:,m,10) = mult(:) * mat(:,i,nc,15)
            
            mat(:,n,m, 6) = mat(:,n,m, 6) + mult(:) * per(:,i,nc,1)
            mat(:,n,m, 7) = mat(:,n,m, 7) + mult(:) * per(:,i,nc,2)
            mat(:,n,m, 8) = mat(:,n,m, 8) + mult(:) * per(:,i,nc,3)
            mat(:,n,m, 9) = mat(:,n,m, 9) + mult(:) * per(:,i,nc,4)
            mat(:,n,m,10) = mat(:,n,m,10) + mult(:) * per(:,i,nc,5)

            rhs(:,n,m) = rhs(:,n,m) + mult(:) * rhs(:,i,nc)
          end do

          do nc = 2, nblk

            j = nc + nblk
            fact(:) = one / mat(:,i,nc,j)
            
            do m = 1, nblk

              mult(:) = -mat(:,p,m,nc) * fact(:)
              do q = nc + 1, nblk
                l = q + nblk
                mat(:,p,m,q) = mat(:,p,m,q) + mult(:) * mat(:,i,nc,l)
              end do

              mat(:,p,m, 6) = mat(:,p,m, 6) + mult(:) * mat(:,i,nc,11)
              mat(:,p,m, 7) = mat(:,p,m, 7) + mult(:) * mat(:,i,nc,12)
              mat(:,p,m, 8) = mat(:,p,m, 8) + mult(:) * mat(:,i,nc,13)
              mat(:,p,m, 9) = mat(:,p,m, 9) + mult(:) * mat(:,i,nc,14)
              mat(:,p,m,10) = mat(:,p,m,10) + mult(:) * mat(:,i,nc,15)

              per(:,p,m,1) = per(:,p,m,1) + mult(:) * per(:,i,nc,1)
              per(:,p,m,2) = per(:,p,m,2) + mult(:) * per(:,i,nc,2)
              per(:,p,m,3) = per(:,p,m,3) + mult(:) * per(:,i,nc,3)
              per(:,p,m,4) = per(:,p,m,4) + mult(:) * per(:,i,nc,4)
              per(:,p,m,5) = per(:,p,m,5) + mult(:) * per(:,i,nc,5)

              rhs(:,p,m) = rhs(:,p,m) + mult(:) * rhs(:,i,nc)

!.... the periodic stuff on the last row

              mult(:) = -per2(:,m,nc) * fact(:)
              do q = nc + 1, nblk
                l = q + nblk
                per2(:,m,q) = per2(:,m,q) + mult(:) * mat(:,i,nc,l)
              end do

              per2(:,m, 6) = per2(:,m, 6) + mult(:) * mat(:,i,nc,11)
              per2(:,m, 7) = per2(:,m, 7) + mult(:) * mat(:,i,nc,12)
              per2(:,m, 8) = per2(:,m, 8) + mult(:) * mat(:,i,nc,13)
              per2(:,m, 9) = per2(:,m, 9) + mult(:) * mat(:,i,nc,14)
              per2(:,m,10) = per2(:,m,10) + mult(:) * mat(:,i,nc,15)
              
              mat(:,n,m, 6) = mat(:,n,m, 6) + mult(:) * per(:,i,nc,1)
              mat(:,n,m, 7) = mat(:,n,m, 7) + mult(:) * per(:,i,nc,2)
              mat(:,n,m, 8) = mat(:,n,m, 8) + mult(:) * per(:,i,nc,3)
              mat(:,n,m, 9) = mat(:,n,m, 9) + mult(:) * per(:,i,nc,4)
              mat(:,n,m,10) = mat(:,n,m,10) + mult(:) * per(:,i,nc,5)

              rhs(:,n,m) = rhs(:,n,m) + mult(:) * rhs(:,i,nc)
            end do

          end do        ! end of loop on nc
          
          per2(:,:,1:nblk) = per2(:,:,nblk+1:2*nblk)

        end do          ! end of loop on i

!.... do the third to the last row

        i = n - 2
        
!.... do the first four columns

        do nc = 1, nblk - 1
        
          j = nblk + nc

          do m = nc, nblk - 1
            l = m + 1
            mult(:) = -mat(:,i,l,j) / mat(:,i,nc,j)
            do q = nc + nblk + 1, nblk * 2
              mat(:,i,l,q) = mat(:,i,l,q) + mult(:) * mat(:,i,nc,q)
            end do

            mat(:,i,l,11) = mat(:,i,l,11) + mult(:) * mat(:,i,nc,11)
            mat(:,i,l,12) = mat(:,i,l,12) + mult(:) * mat(:,i,nc,12)
            mat(:,i,l,13) = mat(:,i,l,13) + mult(:) * mat(:,i,nc,13)
            mat(:,i,l,14) = mat(:,i,l,14) + mult(:) * mat(:,i,nc,14)
            mat(:,i,l,15) = mat(:,i,l,15) + mult(:) * mat(:,i,nc,15)

            per(:,i,l,1) = per(:,i,l,1) + mult(:) * per(:,i,nc,1)
            per(:,i,l,2) = per(:,i,l,2) + mult(:) * per(:,i,nc,2)
            per(:,i,l,3) = per(:,i,l,3) + mult(:) * per(:,i,nc,3)
            per(:,i,l,4) = per(:,i,l,4) + mult(:) * per(:,i,nc,4)
            per(:,i,l,5) = per(:,i,l,5) + mult(:) * per(:,i,nc,5)

            rhs(:,i,l) = rhs(:,i,l) + mult(:) * rhs(:,i,nc)
          end do              

        end do  ! end of loop on nc

!.... for all columns
          
        do nc = 1, nblk

          p  = i + 1
          
          j = nc + nblk
          fact(:) = one / mat(:,i,nc,j)
          
          do m = 1, nblk
            mult(:) = -mat(:,p,m,nc) * fact(:)
            do q = nc + 1, nblk
              l = q + nblk
              mat(:,p,m,q) = mat(:,p,m,q) + mult(:) * mat(:,i,nc,l)
            end do
            
            mat(:,p,m, 6) = mat(:,p,m, 6) + mult(:) * mat(:,i,nc,11)
            mat(:,p,m, 7) = mat(:,p,m, 7) + mult(:) * mat(:,i,nc,12)
            mat(:,p,m, 8) = mat(:,p,m, 8) + mult(:) * mat(:,i,nc,13)
            mat(:,p,m, 9) = mat(:,p,m, 9) + mult(:) * mat(:,i,nc,14)
            mat(:,p,m,10) = mat(:,p,m,10) + mult(:) * mat(:,i,nc,15)

            mat(:,p,m,11) = mat(:,p,m,11) + mult(:) * per(:,i,nc,1)
            mat(:,p,m,12) = mat(:,p,m,12) + mult(:) * per(:,i,nc,2)
            mat(:,p,m,13) = mat(:,p,m,13) + mult(:) * per(:,i,nc,3)
            mat(:,p,m,14) = mat(:,p,m,14) + mult(:) * per(:,i,nc,4)
            mat(:,p,m,15) = mat(:,p,m,15) + mult(:) * per(:,i,nc,5)

            rhs(:,p,m) = rhs(:,p,m) + mult(:) * rhs(:,i,nc)

!.... the periodic stuff on the last row

            mult(:) = -per2(:,m,nc) * fact(:)
            do q = nc + 1, nblk
              l = q + nblk
              per2(:,m,q) = per2(:,m,q) + mult(:) * mat(:,i,nc,l)
            end do

            mat(:,n,m,1) = mat(:,n,m,1) + mult(:) * mat(:,i,nc,11)
            mat(:,n,m,2) = mat(:,n,m,2) + mult(:) * mat(:,i,nc,12)
            mat(:,n,m,3) = mat(:,n,m,3) + mult(:) * mat(:,i,nc,13)
            mat(:,n,m,4) = mat(:,n,m,4) + mult(:) * mat(:,i,nc,14)
            mat(:,n,m,5) = mat(:,n,m,5) + mult(:) * mat(:,i,nc,15)

            mat(:,n,m, 6) = mat(:,n,m, 6) + mult(:) * per(:,i,nc,1)
            mat(:,n,m, 7) = mat(:,n,m, 7) + mult(:) * per(:,i,nc,2)
            mat(:,n,m, 8) = mat(:,n,m, 8) + mult(:) * per(:,i,nc,3)
            mat(:,n,m, 9) = mat(:,n,m, 9) + mult(:) * per(:,i,nc,4)
            mat(:,n,m,10) = mat(:,n,m,10) + mult(:) * per(:,i,nc,5)

            rhs(:,n,m) = rhs(:,n,m) + mult(:) * rhs(:,i,nc)
          end do

        end do  ! end of loop on nc
        
!.... do the second to the last row

        i = n - 1
        
!.... do the first four columns

          do nc = 1, nblk - 1
          
            j = nblk + nc

            do m = nc, nblk - 1
              l = m + 1
              mult(:) = -mat(:,i,l,j) / mat(:,i,nc,j)
              do q = nc + nblk + 1, nblk * 2
                mat(:,i,l,q) = mat(:,i,l,q) + mult(:) * mat(:,i,nc,q)
              end do
              
              mat(:,i,l,11) = mat(:,i,l,11) + mult(:) * mat(:,i,nc,11)
              mat(:,i,l,12) = mat(:,i,l,12) + mult(:) * mat(:,i,nc,12)
              mat(:,i,l,13) = mat(:,i,l,13) + mult(:) * mat(:,i,nc,13)
              mat(:,i,l,14) = mat(:,i,l,14) + mult(:) * mat(:,i,nc,14)
              mat(:,i,l,15) = mat(:,i,l,15) + mult(:) * mat(:,i,nc,15)

              rhs(:,i,l) = rhs(:,i,l) + mult(:) * rhs(:,i,nc)
            end do            

          end do        ! end of loop on nc

!.... for all columns
          
          do nc = 1, nblk

            p  = i + 1
            
            j = nc + nblk
            fact(:) = one / mat(:,i,nc,j)
            
            do m = 1, nblk
              mult(:) = -mat(:,p,m,nc) * fact(:)
              do q = nc + 1, nblk
                l = q + nblk
                mat(:,p,m,q) = mat(:,p,m,q) + mult(:) * mat(:,i,nc,l)
              end do
              
              mat(:,p,m, 6) = mat(:,p,m, 6) + mult(:) * mat(:,i,nc,11)
              mat(:,p,m, 7) = mat(:,p,m, 7) + mult(:) * mat(:,i,nc,12)
              mat(:,p,m, 8) = mat(:,p,m, 8) + mult(:) * mat(:,i,nc,13)
              mat(:,p,m, 9) = mat(:,p,m, 9) + mult(:) * mat(:,i,nc,14)
              mat(:,p,m,10) = mat(:,p,m,10) + mult(:) * mat(:,i,nc,15)

              rhs(:,p,m) = rhs(:,p,m) + mult(:) * rhs(:,i,nc)
            end do

          end do        ! end of loop on nc
        
!.... do the last row
          
        i = n
        
        do nc = 1, nblk - 1
        
          j = nblk + nc

          do m = nc, nblk - 1
            l = m + 1
            mult(:) = -mat(:,i,l,j) / mat(:,i,nc,j)
            do q = nc + nblk + 1, nblk * 2
              mat(:,i,l,q) = mat(:,i,l,q) + mult(:) * mat(:,i,nc,q)
            end do
            rhs(:,i,l) = rhs(:,i,l) + mult(:) * rhs(:,i,nc)
          end do              

        end do  ! end of loop on nc

!=============================================================================!
!                   B A C K W A R D   S U B S T I T U T I O N 
!=============================================================================!

!.... do the last row first

        i = n
        
        do m = nblk, 1, -1
          j = m + nblk
          do q = m+1, nblk
            l = q + nblk
            rhs(:,i,m) = rhs(:,i,m) - rhs(:,i,q) * mat(:,i,m,l)
          end do
          rhs(:,i,m) = rhs(:,i,m) / mat(:,i,m,j)
        end do

!.... do the second to the last row

        i = n - 1
        
        p = i + 1
        
        do m = nblk, 1, -1
          j = m + nblk
          do q = m+1, nblk
            l = q + nblk
            rhs(:,i,m) = rhs(:,i,m) - rhs(:,i,q) * mat(:,i,m,l)
          end do
          
          rhs(:,i,m) = rhs(:,i,m) - rhs(:,p,1) * mat(:,i,m,11)
          rhs(:,i,m) = rhs(:,i,m) - rhs(:,p,2) * mat(:,i,m,12)
          rhs(:,i,m) = rhs(:,i,m) - rhs(:,p,3) * mat(:,i,m,13)
          rhs(:,i,m) = rhs(:,i,m) - rhs(:,p,4) * mat(:,i,m,14)
          rhs(:,i,m) = rhs(:,i,m) - rhs(:,p,5) * mat(:,i,m,15)

          rhs(:,i,m) = rhs(:,i,m) / mat(:,i,m,j)
        end do

!.... now do the rest of the rows in reverse order

        do i = n - 2, 1, -1
        
          p = i + 1
          
          do m = nblk, 1, -1
            j = m + nblk
            do q = m+1, nblk
              l = q + nblk
              rhs(:,i,m) = rhs(:,i,m) - rhs(:,i,q) * mat(:,i,m,l)
            end do

            rhs(:,i,m) = rhs(:,i,m) - rhs(:,p,1) * mat(:,i,m,11)
            rhs(:,i,m) = rhs(:,i,m) - rhs(:,p,2) * mat(:,i,m,12)
            rhs(:,i,m) = rhs(:,i,m) - rhs(:,p,3) * mat(:,i,m,13)
            rhs(:,i,m) = rhs(:,i,m) - rhs(:,p,4) * mat(:,i,m,14)
            rhs(:,i,m) = rhs(:,i,m) - rhs(:,p,5) * mat(:,i,m,15)

            rhs(:,i,m) = rhs(:,i,m) - rhs(:,n,1) * per(:,i,m,1)
            rhs(:,i,m) = rhs(:,i,m) - rhs(:,n,2) * per(:,i,m,2)
            rhs(:,i,m) = rhs(:,i,m) - rhs(:,n,3) * per(:,i,m,3)
            rhs(:,i,m) = rhs(:,i,m) - rhs(:,n,4) * per(:,i,m,4)
            rhs(:,i,m) = rhs(:,i,m) - rhs(:,n,5) * per(:,i,m,5)

            rhs(:,i,m) = rhs(:,i,m) / mat(:,i,m,j)
          end do
                                  
        end do

!.... account for periodicity

        rhs(:,np,1) = rhs(:,1,1)
        rhs(:,np,2) = rhs(:,1,2)
        rhs(:,np,3) = rhs(:,1,3)
        rhs(:,np,4) = rhs(:,1,4)
        rhs(:,np,5) = rhs(:,1,5)
!=============================================================================!
        return
        end
