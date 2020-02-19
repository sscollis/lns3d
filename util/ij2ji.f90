program ij2ji 
implicit none
real, allocatable :: q(:,:,:)
real :: time, Re, Ma, Pr, gamma, Cv
integer :: lstep, nx, ny, nz, ndof,idof,i,j

character*80 :: name

write(*,"('Enter file in IJ order ==> ',$)")
read(*,*) name

!.... read an LNS restart file in IJ order

open(20,file=name,form='unformatted',status='old')
read(20) lstep, time, nx, ny, nz, ndof, Re, Ma, Pr, gamma, Cv
write(*,*) lstep, time, nx, ny, nz, ndof, Re, Ma, Pr, gamma, Cv
allocate( q(ndof,nx,ny) )
read(20) (((q(idof,i,j), idof = 1, 5), i = 1, nx), j = 1, ny)
close(20)

write(*,"('Enter file in JI order ==> ',$)")
read(*,*) name

!.... write an LNS restart file in JI order

open(20,file=name,form='unformatted')
write(20) lstep, time, nx, ny, nz, ndof, Re, Ma, Pr, gamma, Cv
write(20) (((q(idof,i,j), j = 1, ny), i = 1, nx), idof = 1, ndof)
close(20)

end program ij2ji 
