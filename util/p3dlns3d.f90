!-----------------------------------------------------------------------
      program p3dlns3d
!-----------------------------------------------------------------------

!.... Convert a PLOT3D q-file to a restart file for LNS3D

      integer :: nx, ny, ndof=5

      real, allocatable :: v(:,:,:)

      character*80 :: fname
!-----------------------------------------------------------------------

      write(*,"(/,'Enter Plot3d data filename ==> ',$)")
      read(*,"(a)") fname
      open (10, file=fname, form='unformatted')
      read (10) nx, ny, nz
      write(*,"('Nx = ',i4,', Ny = ',i4)") nx, ny
      allocate (v(ny,nx,ndof))
      read (10) dum, dum, dum, dum
      read (10) (((v(j,i,k), i= 1,nx), j= 1,ny), k= 1,ndof)
      close(10)

      write(*,"(/,'Enter LNS3d data filename ==> ',$)")
      read(*,"(a)") fname
      open (10, file=fname, form='unformatted')
      write (10) 0, 0.0, nx, ny, nz, ndof, 1.68e6, 0.058, 1.0, 1.4, 715.0
      write (10) v
      close (10)

      call exit(0)
      end
