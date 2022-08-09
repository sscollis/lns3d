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
          
          real    :: theta, beta, omega, eps_e, x0
          integer :: ires, linear, impl, loctime, iLHS
          logical :: updateLHS = .false.
          logical :: useCalcd = .false.

          logical :: tflag = .false.    ! time traces
          logical :: fflag = .false.    ! temporal statistics
          logical :: sflag = .false.    ! spatial statistics

          logical :: complex_analysis = .false.  ! complex analysis flag

!.... boundary conditions
        
          integer :: top=-1, wall=-1, left=-1, right=-1, wallt=-1
          integer :: extrap=-1
          logical :: lsym=.false., rsym=.false., bsym=.false., tsym=.false.
          logical :: xper=.false., yper=.false.
          logical :: Navier=.true.
          real, allocatable :: bnl(:,:), bnr(:,:), bnb(:,:), bnt(:,:)

          integer :: is, ie, js, je

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
          
          integer :: ny=1, nx=1, nz=1, nsd=3, ndof=5, nbl=0
          logical :: update_nbl=.false., output_nbl=.true.
          real    :: xmin, xmax, ymin, ymax, zmin, zmax
          real    :: yi

!.... space grid and metrics

          real, allocatable :: x(:,:), y(:,:), z(:)
          real, allocatable :: xi(:),  eta(:)
          real, allocatable :: m1(:,:), m2(:,:), m11(:,:), m12(:,:), m22(:,:)
          !$sgi distribute m1(*,block), m2(*,block)
          !$sgi distribute m11(*,block), m12(*,block), m22(*,block)
          real, allocatable :: n1(:,:), n2(:,:), n11(:,:), n12(:,:), n22(:,:)
          !$sgi distribute n1(*,block), n2(*,block)
          !$sgi distribute n11(*,block), n12(*,block), n22(*,block)
          real, allocatable :: m1m1(:,:), m1m2(:,:), m2m2(:,:)
          !$sgi distribute m1m1(*,block), m1m2(*,block), m2m2(*,block)
          real, allocatable :: n1n1(:,:), n1n2(:,:), n2n2(:,:)
          !$sgi distribute n1n1(*,block), n1n2(*,block), n2n2(*,block)
          real, allocatable :: m1n1(:,:), m1n2(:,:), m2n1(:,:), m2n2(:,:)
          !$sgi distribute m1n1(*,block),m1n2(*,block),m2n1(*,block),m2n2(*,block)
          real :: dxi, deta, dz, kz
          
!.... differencing parameters

          integer :: optx = 0
          integer :: opty = 0
          logical :: carp = .false.
          logical :: general = .true.

!.... time grid

          real :: Delt = zero, time = zero, alfa
          real :: dtmax, cflmax, cfl
          integer :: istep, nstep = 1, iter, niter = 1, lstep = 0

!.... mean flow

          real, allocatable :: vm(:,:,:)

!.... inflow disturbances

#if 0
          ! Try out F90 types
          type inflow_type
            character(40) :: filename="inflow.dat"
            logical :: echo=.false.
            complex :: omega=0, alpha=0, beta=0
            !namelist /input/ filename, omega, alpha, beta, echo 
          end type
          type (inflow_type) :: infl
          namelist /inpt/ filename, echo, 
#endif
          character(40) :: inflow_file="inflow.dat"
          logical :: inflow_echo=.false.
          complex :: lomega, lalpha, lbeta
          namelist /inflow/ inflow_file, lomega, lalpha, lbeta, inflow_echo

          real, allocatable :: rhor(:), rhoi(:), ur(:), ui(:)
          real, allocatable :: vr(:), vi(:), wr(:), wi(:), tr(:), ti(:)

!.... sponge variables

          integer :: ispg       
          logical :: compSpg = .true., compSpg2 = .true.
          real, allocatable :: spg(:,:)
          real, allocatable :: spg2(:,:)

!.... the boundary forcing amplitude (should be complex)

          logical :: useAmp = .false.
          character(40) :: amp_file="amp.dat"
          real, allocatable :: wamp(:)
          
!.... damping function
          
          integer :: idamp      
          real, allocatable  :: damp(:,:)

!.... unit numbers

          integer, parameter :: ifile=10, ihist=11

!.... file info

          character(80) :: base, filen
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
