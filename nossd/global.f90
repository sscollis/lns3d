!=============================================================================!
        module global 
!
!  Global variables for LNS3D
!
!  Revised:  8-7-98
!=============================================================================!

!.... Useful constants

          real, parameter :: zero    = 0.0000000000000000000d+0
          real, parameter :: pt25    = 2.5000000000000000000d-1
          real, parameter :: pt33    = 3.3333333333333333333d-1
          real, parameter :: pt5     = 5.0000000000000000000d-1
          real, parameter :: pt66    = 6.6666666666666666666d-1
          real, parameter :: one     = 1.0000000000000000000d+0
          real, parameter :: onept25 = 1.2500000000000000000d+0   
          real, parameter :: onept33 = 1.3333333333333333333d+0
          real, parameter :: onept5  = 1.5000000000000000000d+0
          real, parameter :: two     = 2.0000000000000000000d+0
          real, parameter :: three   = 3.0000000000000000000d+0
          real, parameter :: four    = 4.0000000000000000000d+0
          real, parameter :: pi      = 3.1415926535897932385d+0
          complex, parameter :: im   = (0.0,1.0)

!.... Run parameters

          real    :: Ma, Re, Pr
          real    :: rhoref, Pref, Tref, muref
          
          real    :: alpha, beta, omega, eps_e, x0
          integer :: ires, linear, impl, loctime, iLHS

          logical :: tflag = .false.    ! time traces
          logical :: fflag = .false.    ! temporal statistics
          logical :: sflag = .false.    ! spatial statistics

!.... boundary conditions stuff
        
          integer :: top=-1, wall=-1, left=-1, right=-1, wallt=-1
          integer :: extrap=-1
          logical :: lsym=.false., rsym=.false., bsym=.false., tsym=.false.
          logical :: xper=.false., yper=.false.
          logical :: Navier=.true.
          real, allocatable :: bnl(:,:), bnr(:,:), bnb(:,:), bnt(:,:)

!.... Poinsot & Lele BC parameter

          real :: sigma = zero

!.... fluid properties

          real :: gamma  = 1.4
          real :: gamma1 = 0.4
          real :: cv     = 716.5
          real :: cp     = 1003.1
          real :: Rgas   = 286.6

          integer :: mattyp
          real    :: mu0, T0, S0   ! Sutherland's

!.... problem dimensions
          
          integer :: ny = 1, nx = 1, nz = 1, nsd = 3, ndof = 5, nbl
          real    :: xmin, xmax, ymin, ymax, zmin, zmax
          real    :: yi

!.... space grid and metrics

          real, allocatable :: x(:,:), y(:,:), z(:)
          real, allocatable :: xi(:),  eta(:)
          real, allocatable :: m1(:), m2(:), m11(:), m12(:), m22(:)
          real, allocatable :: n1(:), n2(:), n11(:), n12(:), n22(:)
          real, allocatable :: m1m1(:), m1m2(:), m2m2(:)
          real, allocatable :: n1n1(:), n1n2(:), n2n2(:)
          real, allocatable :: m1n1(:), m1n2(:), m2n1(:), m2n2(:)
          real :: dxi, deta, dz, kz
          
!.... differencing parameters

          integer :: optx = 0
          integer :: opty = 0

!.... time grid

          real :: Delt = zero, time = zero, alfa
          real :: dtmax, cflmax
          integer :: istep, nstep = 1, iter, niter = 1, lstep = 0

!.... mean flow

          real, allocatable :: vm(:)    ! (ny*nx*ndof)

!.... inflow disturbances

          real, allocatable :: rhor(:), rhoi(:), ur(:), ui(:)
          real, allocatable :: vr(:), vi(:), wr(:), wi(:), tr(:), ti(:)

!.... sponge variables

          integer :: ispg       
          real, allocatable :: spg(:)
          real, allocatable :: spg2(:)

!.... the boundary forcing amplitude (should be complex stupid)

          real, allocatable :: wamp(:)
          
!.... damping function
          
          integer :: idamp      
          real, allocatable  :: damp(:)

!.... unit numbers

          integer, parameter :: ifile=10, ihist=11

!.... file info

          character*80 :: base, filen
          integer, parameter :: lfile=80
          integer :: iver = 0, itout, ntout
          
        end module global 

!=============================================================================!
        module material
!       
!  Interfaces for fluid material routines
!
!  Revised:  4-16-96
!=============================================================================!
        interface getmat
          subroutine getmat(t, mu, lm, con, &
                            dmu, d2mu, dlm, d2lm, dcon, d2con)
            real t(:), mu(:), lm(:), con(:), dmu(:), d2mu(:)
            real dlm(:), d2lm(:), dcon(:), d2con(:)
          end
          subroutine sgetmat(t, mu, lm, con, &
                             dmu, d2mu, dlm, d2lm, dcon, d2con)
            real t, mu, lm, con, dmu, d2mu, dlm, d2lm, dcon, d2con
          end
        end interface
        
        end module material
