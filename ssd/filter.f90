!============================================================================!
        subroutine cfilter( ny, nx, ndof, v )
!  
!  Filter a complex 2-D field
!  
!============================================================================!
        implicit none
        
        integer :: nx, ny, ndof
        complex :: v(ny,nx,ndof), f(ny,nx), rhs(ny,nx)
        real    :: lhs(ny,nx,3)
        
        real, parameter :: alpha = 0.1
        real, parameter :: beta  = 0.0
        real, parameter :: d     = (6.0*beta-1.0)/16.0
        real, parameter :: a     = (2.0+3.0*alpha)/4.0
        real, parameter :: b     = (9.0 + 16.0*alpha + 10.0*beta)/16.0
        real, parameter :: c     = (alpha + 4.0*beta)/4.0
        
        integer :: iter, idof
!============================================================================!

        do idof = 1, ndof
        
        f = v(:,:,idof)
        
!.... setup RHS of filter in y

        rhs(1,:) = (15.0 * f(1,:) + 4.0 * f(2,:) - 6.0 * f(3,:) + &
                   4.0 * f(4,:) - f(5,:))/16.0
        rhs(2,:) = 3.0/4.0 * f(2,:) + (f(1,:) + 6.0 * f(3,:) - &
                   4.0 * f(4,:) + f(5,:))/16.0

        rhs(3,:) = 5.0/8.0 * f(3,:) + (-f(1,:) + 4.0 * f(2,:) + &
                   4.0 * f(4,:) - f(5,:))/16.0
        rhs(4:ny-3,:) = a * f(4:ny-3,:) + &
                        b/2.0 * ( f(5:ny-2,:) + f(3:ny-4,:) ) + &
                        c/2.0 * ( f(6:ny-1,:) + f(2:ny-5,:) ) + &
                        d/2.0 * ( f(7:ny,:)   + f(1:ny-6,:) )
        rhs(ny-2,:) = 5.0/8.0 * f(ny-2,:) + (-f(ny-4,:) + 4.0 * f(ny-3,:) + &
                      4.0 * f(ny-1,:) - f(ny,:))/16.0

        rhs(ny-1,:) = 3.0/4.0 * f(ny-1,:) + (f(ny,:) + 6.0 * f(ny-2,:) - &
                      4.0 * f(ny-3,:) + f(ny-4,:))/16.0
        rhs(ny,:) = (15.0 * f(ny,:) + 4.0 * f(ny-1,:) - 6.0 * f(ny-2,:) + &
                    4.0 * f(ny-3,:) - f(ny-4,:))/16.0

!.... setup LHS of filter in y

        if (beta.eq.0.0) then           !.... tridiagonal LHS

          lhs(:,:,1) = alpha
          lhs(:,:,2) = 1.0
          lhs(:,:,3) = alpha
          
!.... correct for the boundaries

          lhs(1,:,:) = 0.0
          lhs(1,:,2) = 1.0
          lhs(2,:,:) = 0.0
          lhs(2,:,2) = 1.0
          lhs(3,:,:) = 0.0
          lhs(3,:,2) = 1.0

          lhs(ny-2,:,:) = 0.0
          lhs(ny-2,:,2) = 1.0
          lhs(ny-1,:,:) = 0.0
          lhs(ny-1,:,2) = 1.0
          lhs(ny,:,:)   = 0.0
          lhs(ny,:,2)   = 1.0

          call csolve2( ny, nx, 1, lhs, rhs )
                
        else
          write(*,*) 'Error:  only tridiagonal LHS is currently supported'
          call exit(1)
        end if
        
        v(:,:,idof) = rhs(:,:)
        
        end do          ! idof
        
        return
        end

!=============================================================================!
        subroutine csolve2( n, nsys, nblk, mat, rhs )
!=============================================================================!
!  
!  Solves a block tridiagonal system of equations without pivoting. 
!  Vectorized over the second index. 
!
!  WARNING: Since pivoting is not done, this algorithm could fail
!
!  input:
!    nsys       : number of systems to vectorize over
!    n          : number of rows in the systems
!    nblk       : dimension of the blocks in each row
!
!  inout:
!    mat        : the block tridiagonal matrix to be solved
!    rhs        : the right hand side (complex)
!=============================================================================!
        implicit none

        integer, intent(in)    :: n, nblk, nsys
        real, intent(inout)    :: mat(n,nsys,nblk,nblk*3)
        complex, intent(inout) :: rhs(n,nsys,nblk)
        
!.... useful constants

        real, parameter :: zero = 0.0, one = 1.0

!.... local variables
        
        integer :: nc, nr, i, j, l, m, p, q
        real    :: mult(nsys), fact(nsys)

!=============================================================================!
!                     F O R W A R D   E L I M I N A T I O N 
!=============================================================================!
        do i = 1, n - 1

!.... do the first four columns

          do nc = 1, nblk - 1
          
            j = nblk + nc

            do m = nc, nblk - 1
              l = m + 1
              mult(:) = -mat(i,:,l,j) / mat(i,:,nc,j)
              do q = nc + nblk + 1, nblk * 3
                mat(i,:,l,q) = mat(i,:,l,q) + mult(:) * mat(i,:,nc,q)
              end do
              rhs(i,:,l) = rhs(i,:,l) + mult(:) * rhs(i,:,nc)
            end do            

          end do        ! end of loop on nc

!.... for all columns
          
          do nc = 1, nblk

            p  = i + 1
            
            j = nc + nblk
            fact(:) = one / mat(i,:,nc,j)
            
            do m = 1, nblk
              mult(:) = -mat(p,:,m,nc) * fact(:)
              do q = nc + 1, nblk * 2
                l = q + nblk
                mat(p,:,m,q) = mat(p,:,m,q) + mult(:) * mat(i,:,nc,l)
              end do
              rhs(p,:,m) = rhs(p,:,m) + mult(:) * rhs(i,:,nc)
            end do

          end do        ! end of loop on nc
          
        end do          ! end of loop on i

!.... do the last row
          
        i = n
        
        do nc = 1, nblk - 1
        
          j = nblk + nc

          do m = nc, nblk - 1
            l = m + 1
            mult(:) = -mat(i,:,l,j) / mat(i,:,nc,j)
            do q = nc + nblk + 1, nblk * 2
              mat(i,:,l,q) = mat(i,:,l,q) + mult(:) * mat(i,:,nc,q)
            end do
            rhs(i,:,l) = rhs(i,:,l) + mult(:) * rhs(i,:,nc)
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
            rhs(i,:,m) = rhs(i,:,m) - rhs(i,:,q) * mat(i,:,m,l)
          end do
          rhs(i,:,m) = rhs(i,:,m) / mat(i,:,m,j)
        end do

!.... now do the rest of the rows in reverse order

        do i = n - 1, 1, -1
        
          p = i + 1
          
          do m = nblk, 1, -1
            j = m + nblk
            do q = m+1, nblk
              l = q + nblk
              rhs(i,:,m) = rhs(i,:,m) - rhs(i,:,q) * mat(i,:,m,l)
            end do
            do q = 1, nblk
              l = q + nblk * 2
              rhs(i,:,m) = rhs(i,:,m) - rhs(p,:,q) * mat(i,:,m,l)
            end do
            rhs(i,:,m) = rhs(i,:,m) / mat(i,:,m,j)
          end do
                                  
        end do

        return
        end

!============================================================================!
        subroutine filter( ny, nx, f )
!  
!  Filter a REAL 2-D field
!  
!============================================================================!
        implicit none
        
        integer:: nx, ny
        real :: f(ny,nx), rhs(ny,nx)
        real, allocatable :: lhs(:,:,:)
        
        real, parameter :: alpha = 0.1
        real, parameter :: beta  = 0.0
        real, parameter :: d     = (6.0*beta-1.0)/16.0
        real, parameter :: a     = (2.0+3.0*alpha)/4.0
        real, parameter :: b     = (9.0 + 16.0*alpha + 10.0*beta)/16.0
        real, parameter :: c     = (alpha + 4.0*beta)/4.0
        
        integer :: iter
!============================================================================!

!.... setup RHS of filter in y

        rhs(1,:) = (15.0 * f(1,:) + 4.0 * f(2,:) - 6.0 * f(3,:) + &
                   4.0 * f(4,:) - f(5,:))/16.0
        rhs(2,:) = 3.0/4.0 * f(2,:) + (f(1,:) + 6.0 * f(3,:) - &
                   4.0 * f(4,:) + f(5,:))/16.0

        rhs(3,:) = 5.0/8.0 * f(3,:) + (-f(1,:) + 4.0 * f(2,:) + &
                   4.0 * f(4,:) - f(5,:))/16.0
        rhs(4:ny-3,:) = a * f(4:ny-3,:) + &
                        b/2.0 * ( f(5:ny-2,:) + f(3:ny-4,:) ) + &
                        c/2.0 * ( f(6:ny-1,:) + f(2:ny-5,:) ) + &
                        d/2.0 * ( f(7:ny,:)   + f(1:ny-6,:) )
        rhs(ny-2,:) = 5.0/8.0 * f(ny-2,:) + (-f(ny-4,:) + 4.0 * f(ny-3,:) + &
                      4.0 * f(ny-1,:) - f(ny,:))/16.0

        rhs(ny-1,:) = 3.0/4.0 * f(ny-1,:) + (f(ny,:) + 6.0 * f(ny-2,:) - &
                      4.0 * f(ny-3,:) + f(ny-4,:))/16.0
        rhs(ny,:) = (15.0 * f(ny,:) + 4.0 * f(ny-1,:) - 6.0 * f(ny-2,:) + &
                    4.0 * f(ny-3,:) - f(ny-4,:))/16.0

!.... setup LHS of filter in y

        if (beta.eq.0.0) then           !.... tridiagonal LHS

          allocate ( lhs(ny,nx,3) )
          
          lhs(:,:,1) = alpha
          lhs(:,:,2) = 1.0
          lhs(:,:,3) = alpha
          
!.... correct for the boundaries

          lhs(1,:,:) = 0.0
          lhs(1,:,2) = 1.0
          lhs(2,:,:) = 0.0
          lhs(2,:,2) = 1.0
          lhs(3,:,:) = 0.0
          lhs(3,:,2) = 1.0

          lhs(ny-2,:,:) = 0.0
          lhs(ny-2,:,2) = 1.0
          lhs(ny-1,:,:) = 0.0
          lhs(ny-1,:,2) = 1.0
          lhs(ny,:,:)   = 0.0
          lhs(ny,:,2)   = 1.0

          call solve2( ny, nx, 1, lhs, rhs )
          
          deallocate ( lhs )
        
        else
          write(*,*) 'Error:  only tridiagonal LHS is currently supported'
          call exit(1)
        end if
        
        f(:,:) = rhs(:,:)
        
        return
        end

!=============================================================================!
        subroutine solve2( n, nsys, nblk, mat, rhs )
!=============================================================================!
!  
!  Solves a block tridiagonal system of equations without pivoting. 
!  Vectorized over the second index. 
!
!  WARNING: Since pivoting is not done, this algorithm could fail
!
!  input:
!    nsys       : number of systems to vectorize over
!    n          : number of rows in the systems
!    nblk       : dimension of the blocks in each row
!
!  inout:
!    mat        : the block tridiagonal matrix to be solved
!    rhs        : the right hand side
!=============================================================================!
        implicit none

        integer, intent(in) :: n, nblk, nsys
        real, intent(inout) :: mat(n,nsys,nblk,nblk*3), rhs(n,nsys,nblk)
        
!.... useful constants

        real, parameter :: zero = 0.0, one = 1.0

!.... local variables
        
        integer :: nc, nr, i, j, l, m, p, q
        real    :: mult(nsys), fact(nsys)

!=============================================================================!
!                     F O R W A R D   E L I M I N A T I O N 
!=============================================================================!
        do i = 1, n - 1

!.... do the first four columns

          do nc = 1, nblk - 1
          
            j = nblk + nc

            do m = nc, nblk - 1
              l = m + 1
              mult(:) = -mat(i,:,l,j) / mat(i,:,nc,j)
              do q = nc + nblk + 1, nblk * 3
                mat(i,:,l,q) = mat(i,:,l,q) + mult(:) * mat(i,:,nc,q)
              end do
              rhs(i,:,l) = rhs(i,:,l) + mult(:) * rhs(i,:,nc)
            end do            

          end do        ! end of loop on nc

!.... for all columns
          
          do nc = 1, nblk

            p  = i + 1
            
            j = nc + nblk
            fact(:) = one / mat(i,:,nc,j)
            
            do m = 1, nblk
              mult(:) = -mat(p,:,m,nc) * fact(:)
              do q = nc + 1, nblk * 2
                l = q + nblk
                mat(p,:,m,q) = mat(p,:,m,q) + mult(:) * mat(i,:,nc,l)
              end do
              rhs(p,:,m) = rhs(p,:,m) + mult(:) * rhs(i,:,nc)
            end do

          end do        ! end of loop on nc
          
        end do          ! end of loop on i

!.... do the last row
          
        i = n
        
        do nc = 1, nblk - 1
        
          j = nblk + nc

          do m = nc, nblk - 1
            l = m + 1
            mult(:) = -mat(i,:,l,j) / mat(i,:,nc,j)
            do q = nc + nblk + 1, nblk * 2
              mat(i,:,l,q) = mat(i,:,l,q) + mult(:) * mat(i,:,nc,q)
            end do
            rhs(i,:,l) = rhs(i,:,l) + mult(:) * rhs(i,:,nc)
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
            rhs(i,:,m) = rhs(i,:,m) - rhs(i,:,q) * mat(i,:,m,l)
          end do
          rhs(i,:,m) = rhs(i,:,m) / mat(i,:,m,j)
        end do

!.... now do the rest of the rows in reverse order

        do i = n - 1, 1, -1
        
          p = i + 1
          
          do m = nblk, 1, -1
            j = m + nblk
            do q = m+1, nblk
              l = q + nblk
              rhs(i,:,m) = rhs(i,:,m) - rhs(i,:,q) * mat(i,:,m,l)
            end do
            do q = 1, nblk
              l = q + nblk * 2
              rhs(i,:,m) = rhs(i,:,m) - rhs(p,:,q) * mat(i,:,m,l)
            end do
            rhs(i,:,m) = rhs(i,:,m) / mat(i,:,m,j)
          end do
                                  
        end do
!=============================================================================!
        return
        end
