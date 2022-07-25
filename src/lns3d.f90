!=============================================================================!
        program lns3d
!  
!  Linearized 2-D and 3-D compressible Navier-Stokes solver and
!  Nonlinear 2-D compressible Navier-Stokes solver.
!
!  Equations are formulated in primative form with the state vector
!  V = {rho, u, v, w, T} in a generalized coordinate system {xi, eta, z}
!  
!  Spatial discretization is up-to sixth-order central in xi and eta 
!  with a single Fourier mode in z.
!
!  Time advancement is either: RK4, backward Euler, or two-step backward.
!  For the implicit schemes, approximate factorization is used to solve
!  the linear system of equations.
!
!  A variety of boundary conditions are supported.  See the routines: 
!  input.f90, itrbc.f90, itrbc3d.f90 for details.
!
!  When adding new boundary conditions, all routines with *bc*.f90 must be
!  updated.
!
!  Author:     S. Scott Collis
!
!  Revised:    12-18-00
!
!  Copyright:  S. Scott Collis
!              Department of Mechanical Engineering and Materials Science
!              Rice University, MS 321
!              Houston, TX 77005-1892
!              (713) 348-3617
!              collis@rice.edu
!=============================================================================!
        use global
        use local2d
        use local3d
        implicit none

        character(80) :: code='lns3D$'

        logical :: echo = .false.
        namelist /extras/ is, ie, js, je, echo, updateLHS, useCalcd, &
                          update_nbl, nbl, output_nbl, useAmp, &
                          amp_file, compSpg, compSpg2
        integer, parameter :: mfile=1
        integer :: narg, iarg, nfile=0, jfile(mfile)
        character(80) :: filename, temp, arg
        logical :: nlist=.false.

!$      integer, external :: omp_get_num_threads, omp_get_thread_num
!$      integer, external :: omp_get_num_procs
!=============================================================================!
        write(*,"(' LNS is initializing...',/)")
!$omp parallel
!$      if (omp_get_thread_num() == 0) then
!$        write(*,*) 'Running on ',omp_get_num_threads(),' processor(s)'
!$        write(*,*) 'There are  ',omp_get_num_procs(),' processor(s) available'
!$      end if
!$omp end parallel

        call init

!.... input parameters

        call input

!.... generate the grid

        call gengrid
        
!.... Boundary condition indices

        is = 1; ie = nx
        js = 1; je = ny

!.... Read filename from arguments

        narg = iargc()
        do iarg = 1, narg
          call getarg(iarg,arg)
          if (arg(1:1) .ne. '-') then
            nfile = nfile + 1
            if (nfile .gt. mfile) then
              write(*,*) '>> Error in argument list, too many files'
              call exit(1)
            end if
            jfile(nfile) = iarg
          endif
        end do
        if (nfile.gt.0) then
          call getarg(jfile(1),temp)
          filename = temp
          inquire(file=filename,exist=nlist)
          if (.not.nlist) then
            write(*,*) "Namelist file ", &
              filename(1:index(filename,' ')),"does not exist"
            filename = 'lns3d.nml'
          endif
        else
          filename = 'lns3d.nml'
        endif
        write(*,*) "Reading namelist input from ", &
          filename(1:index(filename,' '))
        inquire(file=filename,exist=nlist)
        if (nlist) then  
          open(10,file=filename,status='old',err=100)
          read(10,extras)
          close(10)
  100     continue
          if (echo) then
            write(*,extras)
          end if 
        else
          write(*,*) "Namelist ",filename(1:index(filename,' ')),&
                     "does not exist, using defaults"
        endif

!.... generate the mean flow

        call genmean
        
!.... generate the inflow disturbance

        call gendist

!.... generate the potential boundary values

        call potential
        
!.... generate the sponge term

        call init_sponge

!.... generate the buffer term

        call init_buffer

!.... generate the damping term

        call init_damper

!.... allocate local variables

        call mlocal2d
        if (complex_analysis) call mlocal3d

!.... build the matrix if using linear solver
        
        if (linear.eq.1) call genmtrx(vm,kz,vm)
                
!.... main Navier-Stokes driver

        if (complex_analysis) then
          if (linear.eq.1) then   ! linear
            if (impl.eq.0) then
              call expdrv3d
            else if (impl.eq.1 .or. impl.eq.2 .or. impl.eq.3) then
              call impdrv3d
#ifdef USE_ARPACK
            else if (impl.eq.4) then ! eigenvalue analysis
              call eigdrv3d
#endif
            else
              call error(code,'Invalid value for impl$')
            end if
          else
            call error(code,'3D code must be linear$')
          end if
        else                      ! not complex
          if (impl.eq.0) then
            call expdrv
          else if (impl.eq.1 .or. impl.eq.2 .or. impl.eq.3) then
            call impdrv
          else if (impl.eq.4) then
            call error(code,'impRK not currently working$')
!           call impRK
          else
            call error(code,'Invalid value for impl$')
          end if
        end if

!.... end program

        stop
        end
