!=============================================================================!
        program lns3d
!  
!  Linearized 2-D and 3-D compressible Navier-Stokes solver and
!  Nonlinear 2-D compressible Navier-Stokes solver.
!
!  Equations are formulated in primative form with the state vector
!  V = {rho, u, v, w, T} in a generalized coordinate system {xi,eta,z}
!  
!  Spatial discretization is fourth-order central in xi and eta with a single
!  Fourier mode in z.
!
!  Time advancement is either: RK4, backward Euler, or two-step backward.
!  For the implicit schemes, approximate factorization is used to solver
!  the linear system of equations.
!
!  A variety of boundary conditions are supported.  See the routines: 
!  input.f90, itrbc.f90, itrbc3d.f90 for details.
!
!  Author:     S. Scott Collis
!
!  Revised:    6-3-97
!
!  Copyright:  S. Scott Collis
!              Department of Mechanical Engineering and Materials Science
!              Rice University, MS 321
!              Houston, TX 77005-1892
!              (713) 527-8101 x3617
!=============================================================================!
        use stuff
        use local
        implicit none
        integer :: ier, mem = 0

        logical :: comp

!.... SSD stuff for use on the Cray

        integer      :: iloc
        character*80 :: name='FASTDIR', value
        logical      :: mexist, async=.true.
        integer      :: lrec, istat
        character*80 :: code='lns3D$'
!=============================================================================!

        write(*,"(' LNS is initializing...',/)")

!.... input parameters

        call input

!.... Determine whether to use complex or real analysis

        if ( kz .ne. zero .or. linear .eq. 2 ) then
          comp = .true.
        else
          comp = .false.
        end if
        if (linear .eq. 2) linear = 1

!.... generate the grid

        call gengrid
        
!.... generate the mean flow

        if (linear.eq.1) call genmean
        
!.... generate the inflow disturbance

        call gendist

!.... generate the potential boundary values

        call potential
        
!.... allocate memory for local variables

        call mlocal('allocate  ')

!.... open the ssd (if required)

        if (linear.eq.1 .or. impl.ne.0) then
          inquire(FILE=MATRIX, EXIST=mexist)
          if (mexist) then
            call getenv(name,value)
            iloc = index(value(2:len(value)),' ')
            write(*,*) 'rm '//value(1:iloc)//'/MATRIX'
            call ishell('rm '//value(1:iloc)//'/MATRIX')
          end if
          if (async) then
            call opendr(MATRIX, idxMTRX, 512, 2, ier)
            call asyncdr(MATRIX,ier)
          else
            call opendr(MATRIX, idxMTRX, 512, 1, ier)
          end if
        end if

!.... generate the sponge term (if required)

        if (ispg .gt. 0) then
          mem = ny*nx
          allocate( spg(ny*nx), STAT=ier)
          if (ispg .ge. 2) then
            mem = mem + ny*nx
            allocate( spg2(ny*nx), STAT=ier)
          end if
          if (ier .ne. 0) call error(code,'Insufficient Memory for spg$')
          call sponge( spg, spg2 )
        end if

!.... generate the buffer term

        if (eps_e .ne. zero) call buffer()

!.... generate the damping term

        mem = mem + ny*nx
        allocate( damp(ny*nx), STAT=ier)
        if (ier .ne. 0) call error(code,'Insufficient Memory for damp$')
        call damper(damp)

        write(*,"(' LNS3D allocated  ===> ',1pe13.6,' words')") float(mem)

!.... build the matrix if using linear solver
        
        if (linear.eq.1) call genmtrx(vm,kz,comp)
                
!.... main Navier-Stokes driver

!       if (kz.ne.zero .or. omega.ne.zero) then
!       if (kz.ne.zero) then
        if (comp) then
          call mlocal3d
          if (linear.eq.1) then
            if (impl.eq.0) then
              call expdrv3d
            else if (impl.eq.1 .or. impl.eq.2 .or. impl.eq.3) then
              call impdrv3d
            else
              call error(code,'Invalid value for impl$')
            end if
          else
            call error(code,'3D code must be linear$')
          end if
        else
          if (impl.eq.0) then
            call expdrv
          else if (impl.eq.1 .or. impl.eq.2 .or. impl.eq.3) then
            call impdrv
          else if (impl.eq.4) then
            call error(code,'impRK not currently working$')
            call impRK
          else
            call error(code,'Invalid value for impl$')
          end if
        end if

!.... close the ssd

        if (linear.eq.1 .or. impl.ne.0) then
          call closdr(MATRIX)
        end if

!.... end program

        stop
        end
