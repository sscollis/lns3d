!=============================================================================!
      program mkamp
!  
!  Determine the amplitude of a forced acoustic wave on the top boundary
!  accounting for viscous decay of the plane acoustic wave using the
!  formulae from Pierce.
!
!  Author:   S. Scott Collis
!
!  Written:  8-22-96
!=============================================================================!
      use const
      implicit none

      real, allocatable :: vm(:,:,:)

      real :: Ma, Re, Pr, gamma, cv, Rgas, time
      integer :: lstep, nx, ny, nz, ndof
      integer :: nxm, nym, nzm, ndofm
      integer :: i, j, k, idof, ix

      character*80 file, temp
      integer :: iloc, iend

      integer, external :: iargc

      real, allocatable :: rhom(:,:), um(:,:), tm(:,:), cm(:,:)
      real, allocatable :: x(:,:), y(:,:)

      real, allocatable :: wamp(:)
      logical :: lamp

      real :: tmp
      integer :: itmp

      real :: kk, alpha, delta, x0, r, rho0, c0, omega, th
!=============================================================================!

!.... read in the mean flow file

      open(unit=10, file='mean.dat', form='unformatted', &
           status='old', err=1000)
      read(10,err=1000) itmp, tmp, nx, ny, nz, ndof, &
                        Re, Ma, Pr, gamma, cv
      allocate( vm(ny,nx,ndof) )
      read(10,err=1000) vm
      close(10)

!.... read in the grid file

      open (unit=10, file='grid.dat', form='unformatted', status='unknown')
      read(10) nx, ny, nz
      allocate (x(ny,nx), y(ny,nx))
      read(10) (((x(j,i), i = 1, nx), j = 1, ny), k = 1, nz), &
               (((y(j,i), i = 1, nx), j = 1, ny), k = 1, nz), &
               (((   tmp, i = 1, nx), j = 1, ny), k = 1, nz)
      close(10)
       
!.... compute the acoustic wave

      allocate( rhom(ny,nx), um(ny,nx), tm(ny,nx), cm(ny,nx), wamp(nx) )

      rhom = vm(:,:,1)
      um   = vm(:,:,2)
      tm   = vm(:,:,5)
      cm   = sqrt( tm ) / Ma

      write(*,"('Correcting for Re = ',1pe13.6,', Pr ',1pe13.6)") Re, Pr
      rho0  = one
      
      if (Re .eq. zero) then
        delta = zero
      else
        delta = (one/Re) * pt5 / rho0 * ( onept33 + (gamma - one) / Pr )
      end if
        
      x0 = x(ny,1) 
      write(*,"('Origin, X0 = ',1pe13.6)") x0

      write(*,"('Enter the angular frequency ==> ',$)")
      read(*,*) omega

      open(10,file='amp.dat')
      j = ny
      do i = 1, nx
        kk     = omega / ( cm(j,i) + um(j,i) )
        alpha  = omega**2 * delta / ( cm(j,i) + um(j,i) )**3
        wamp(i) = exp( -alpha * (x(j,i) - x0) )
        write(10,"(2(1pe20.13,1x))") x(j,i), wamp(i)
      end do
      close(10)

      call exit(0)      

1000  write(*,"('>> Error reading mean field:  mean.dat')")
      call exit(1)

      end
