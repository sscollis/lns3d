!-----------------------------------------------------------------------
      program p3dlns3d
!-----------------------------------------------------------------------

!.... Convert a convert a HDIR field to Plot3d format 

      integer :: nx, ny, ndof=5

      real, allocatable :: v(:,:,:), x(:), y(:), cur(:)
      real, allocatable :: bx(:), by(:), bn1(:), bn2(:)

      character*80 :: fname
!-----------------------------------------------------------------------

      write(*,"(/,'Enter DIRST data filename ==> ',$)")
      read(*,"(a)") fname
      open (10, file=fname, form='unformatted')
      read (10) nx, ny
      allocate (v(ny,nx,ndof), x(nx), y(ny), cur(nx), &
                bx(nx), by(nx), bn1(nx), bn2(nx))
      v(:,:,1) = 1.0
      read (10) (((v(j,i,k+1), j= 1,ny), i= 1,nx), k= 1,4)
      read (10) (x(i), i= 1,nx)
      read (10) (y(j), j= 1,ny)
      read (10) (cur(i), i= 1,nx), &
                (bx(i), i= 1,nx),  (by(i), i=1, nx), &
                (bn1(i), i= 1,nx), (bn2(i), i= 1,nx)
      close (10)

      open(1, file='grid.dat', form='unformatted')
      write(1) nx, ny, 1
      write(1) ((bx(i)+bn1(i)*y(j), i= 1,nx), j= 1,ny), &
               ((by(i)+bn2(i)*y(j), i= 1,nx), j= 1,ny), &
               ((              0.0, i= 1,nx), j= 1,ny)
      close(1)
      
      open(1, file='comp.dat', form='unformatted')
      write(1) nx, ny, 1
      write(1) ((x(i), i= 1,nx), j= 1,ny), &
               ((y(j), i= 1,nx), j= 1,ny), &
               (( 0.0, i= 1,nx), j= 1,ny)
      close(1)

      write(*,"(/,'Enter Plot3d data filename ==> ',$)")
      read(*,"(a)") fname
      open (10, file=fname, form='unformatted')
      write (10) nx, ny, 1
      write(*,"('Nx = ',i4,', Ny = ',i4)") nx, ny
      write (10) 0.0, 0.0, 0.0, 0.0
      write (10) (((v(j,i,k), i= 1,nx), j= 1,ny), k= 1,ndof)
      close(10)

      call exit(0)
      end
