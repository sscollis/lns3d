program main
      integer :: nx, ny, nz, nsd=3, isd
      real*4, allocatable :: xyz(:,:,:,:)
      real*8, allocatable :: xyz2(:,:,:,:)

      open(10,file='grid.dat',form='unformatted')
      read(10) nx, ny, nz
      allocate( xyz(nx,ny,nz,nsd), xyz2(nx,ny,nz,nsd) )
      read(10) ((((xyz(i,j,k,isd),i=1,nx),j=1,ny),k=1,nz),isd=1,nsd)
      close(10)
      xyz2 = xyz
      open(10,file='grid.r8',form='unformatted')
      write(10) nx, ny, nz
      write(10) ((((xyz2(i,j,k,isd),i=1,nx),j=1,ny),k=1,nz),isd=1,nsd)
      close(10)
end program main
