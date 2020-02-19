program getevec

! Getevec for use with LNS3d eigendriver

! f90 -r8 -O getevec.f90 -o getevec ~nawaf/ARPACK/libarpack_SGI.a

implicit none

integer :: nx, ny, nz, nev, ndof
real :: time, Re, Ma, Pr, gamma, cv
complex, allocatable :: eval(:), evec(:,:,:,:)
real, allocatable :: reval(:,:), xy(:,:,:)
real :: tmp
complex :: scale

integer :: i, j, k, idof, which

!.... read in the grid.dat

open (10, file='grid.dat', form='unformatted', status='old')
read(10) nx, ny, nz
allocate(xy(2,nx,ny))
read(10) (((xy(1,i,j), i = 1, nx), j = 1, ny), k = 1, nz), &
         (((xy(2,i,j), i = 1, nx), j = 1, ny), k = 1, nz), &
         (((      tmp, i = 1, nx), j = 1, ny), k = 1, nz)
close(10)

!.... read in the eig.dat

open(10,file='eig.dat',form='unformatted', status='old')
read(10) nev, time, nx, ny, nz, ndof, Re, Ma, Pr, gamma, cv
allocate( eval(nev), reval(nev,2), evec(ndof,nx,ny,nev) )
read(10) eval
read(10) evec
close(10)

reval(:,1) = real(eval)
reval(:,2) = imag(eval)

call dmout(6, nev, 2, reval, nev, -6, 'Eigenvalues (Real, Imag)')

!.... write out restart file

101 continue
write (*,"(/,'Which eigenfunction ==> ',$)")
read (*,*) which
if ( which.lt.0 .or. which.gt.nev) goto 101
if (which .ne. 0) then

!.... normalized the eigenvector

  if (.false.) then
  scale = evec(1,1,1,which)
  do j = 1, ny
    do i = 1, nx
      do idof = 1, ndof
        evec(idof,i,j,which) = evec(idof,i,j,which)/scale
      end do
    end do
  end do
  end if

  open(11,file='evec.dat', form='unformatted', status='unknown')
  write(11) which, 0.0, nx, ny, nz, ndof, Re, Ma, Pr, gamma, cv
  write(11) evec(:,:,:,which)
  close(11)
  j = 1
  do i = 1, nx
    write(20,10) xy(1,i,j), &
         real(evec(1,i,j,which)), imag(evec(1,i,j,which)), &
         real(evec(2,i,j,which)), imag(evec(2,i,j,which)), &
         real(evec(3,i,j,which)), imag(evec(3,i,j,which)), &
         real(evec(4,i,j,which)), imag(evec(4,i,j,which)), &
         real(evec(5,i,j,which)), imag(evec(5,i,j,which))
  end do
  close(20)
  i = 1
  do j = 1, ny
    write(21,10) xy(2,i,j), &
         real(evec(1,i,j,which)), imag(evec(1,i,j,which)), &
         real(evec(2,i,j,which)), imag(evec(2,i,j,which)), &
         real(evec(3,i,j,which)), imag(evec(3,i,j,which)), &
         real(evec(4,i,j,which)), imag(evec(4,i,j,which)), &
         real(evec(5,i,j,which)), imag(evec(5,i,j,which))
  end do
  close(21)
  goto 101
end if

10  format(12(e12.5,1x))

end program getevec
