!-----------------------------------------------------------------------------
        subroutine wdata (name, q, nx, ny, nz, ndof, Ma, alpha, Re, time)
!  
!       write out a plot3d data file 
!  
!-----------------------------------------------------------------------------
        implicit none
        
        character*80 :: name, msg
        integer :: nx, ny, nz, ndof
        real :: q(ny,ndof,nx,nz), Ma, alpha, Re, time
        integer :: i, j, k, idof
        
        open (unit=10, file=name, form='unformatted', status='unknown', &
              err=1000)
        
        write(10,err=1000) nx, ny, nz
        write(10,err=1000) Ma, 0.0, Re, time
        write(10,err=1000) ((((q(j,idof,i,k), i = 1, nx), j = 1, ny),  &
                             k = 1, nz), idof = 1, ndof)
        close(10,err=1000)

        
        return

1000    write(msg,"('Error writing to file ',a,'\')") name
        call error('wdata\',msg)

        end
