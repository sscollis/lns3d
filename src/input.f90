!=============================================================================!
        subroutine input
!=============================================================================!
        use global
        use material
        use sponge
        use buffer
        implicit none
        
        real :: Lz, Uref
        character(80) :: code='Input$'
!=============================================================================!

        Ma = 0.3
!       write (*,"('Enter Ma ==> ',$)")
        read (*,*) Ma

!.... Constant mu or Sutherland's law

        mattyp = 0
!       write (*,"(/,'(0) for constant Mu, (1) for Sutherland ==> ',$)")
        read (*,*) mattyp

        Re = 10000
!       write (*,"(/,'Enter Re ==> ',$)") 
        read (*,*) Re

!       write (*,"('Enter Pinf (Pa) , Tinf (K) ==> ',$)")
        read (*,*) Pref, Tref

        if (mattyp .eq. 0) then              ! Constant viscosity
          if (Re.eq.zero) then
            Navier = .false.
            Re = 1.0e99
          else
            Navier = .true.
          end if
        else if (mattyp .eq. 1) then         ! Sutherland's viscosity law 
          mu0    = 1.715336725523065e-05     ! (AIR) from White
          T0     = 273.0                     ! K
          S0     = 110.4                     ! K
          Rgas   = 286.6
          muref  = mu0*(Tref/T0)**onept5*(T0+S0)/(Tref+S0)
          rhoref = Pref / ( Rgas * Tref )
          Uref   = Re / rhoref * muref
          write(*,*) 'Uref * Lref= ', Uref
          S0     = S0 / Tref                 ! nondimensionalize
        else if (mattyp .eq. 2 .or. &        ! linear viscosity
                 mattyp .eq. 3 ) then        ! No second viscosity

!....     Note that right now I match Sutherland at the reference temperature
!....     So, if you want to match at the wall temperature, you need to make
!....     that the reference temperature.

          mu0    = 1.715336725523065e-05     ! (AIR) from White
          T0     = 273.0                     ! K
          S0     = 110.4                     ! K
          Rgas   = 286.6
          muref  = mu0*(Tref/T0)**onept5*(T0+S0)/(Tref+S0)
          rhoref = Pref / ( Rgas * Tref )
          Uref   = Re / rhoref * muref
          write(*,*) 'Uref * Lref= ', Uref
        else
          call error('input$','Illegal value for mattyp$')
        end if

        Pr = 0.7
!       write (*,"('Enter Pr ==> ',$)")
        read (*,*) Pr
!       write (*,"(/,'Enter dtmax, cflmax ==> ',$)")
        read (*,*) dtmax, cflmax
        nstep = 100
!       write (*,"(/,'Enter nstep ==> ',$)")
        read (*,*) nstep

!.... define the differencing scheme

        optx = 0
        opty = 0
!       write (*,"(/,'Enter optx, opty ==> ',$)")
        read (*,*) optx, opty
                
!.... define periodic directions

        xper = .false.
        yper = .false.
!       write (*,"(/,'Enter xper, yper ==> ',$)")
        read (*,*) xper, yper

!.... set restart file name

        base = 'output'
        base = base(1:index(base,' ')-1)//char(0)       ! null terminate
        filen = char(0)
        
!       write (*,"(/,'Enter itout, ntout ==> ',$)")
        read (*,*) itout, ntout
        
!.... residual output flag
!
!       0 = normal
!       1 = form exit residual and print norm
!       2 = 1 + save Plot3d file of residual on last step
!
!       write (*,"(/,'Enter ires ==> ',$)")
        read (*,*) ires
        
!.... sponge flag
!
!       0 = no sponge
!       1 = sponge on outflow
!       2 = sponge on top and outflow
!       3 = sponge on outflow, and acoustic wave forced on top
        
!       write(*,"(/,'Enter ispg ==> ',$)")
        read(*,*) ispg

!.... linear flag
!
!       0 = nonlinear
!       1 = linear
!       2 = linear / complex

!       write(*,"(/,'Enter linear ==> ',$)")
        read(*,*) linear
        
!.... implicit flag
!
!       0 = explicit (fourth-order Runge-Kutta)
!       1 = implicit (first-order backward Euler)
!       2 = implicit (second-order two-step backward)
!       3 = implicit (second-order mid-point rule)
!       4 = eigen-analysis
!
!.... iLHS   = order of LHS (0 = second-order, 1 = fourth-order)
!.... niter  = number of implicit iterations
!.... loctim = local time steping (1/0)

!       write(*,"(/,'Enter impl, iLHS, niter, loctim ==> ',$)")
        read(*,*) impl, iLHS, niter, loctime

!.... set Carpenter boundary difference flag

        carp = .false.
        if (impl.eq.0) carp = .true.

!=============================================================================!
!     L i n e a r   B o u n d a r y   C o n d i t i o n s
!=============================================================================!
!
!.... Wall BC flag
!
!       0 = full continuity at the wall
!       1 = drho/dn = zero
!       2 = Extrapolate to get rho at the wall
!       3 = wall bump bc, (3D version only)
!       4 = Poinsot & Lele modified continuity at the wall
!
!.... Wall temperature BC flag
!
!       0 = isothermal wall
!       1 = hard FD adiabatic
!       2 = soft adiabatic energy equation
!
!.... top BC flag
!
!       0 = zero disturbance BC
!       1 = Forced acoustic wave with characteristics
!
!.... Left and Right BC flags
!
!       0 = zero disturbance BC
!       1 = nonreflecting (c4=zero)
!       2 = symmetry plane, PL continuity [left side only]
!       3 = nonreflecting (characteristic filter, not in LHS)
!       4 = forced eigenfunction inflow disturbance [left side only]
!       5 = acoustic inflow disturbance [left side only]
!       6 = Poinsot & Lele outflow with viscous conditions
!       7 = exact symmetry [left side only]
!       8 = hold ic on the boundary [3D right boundary only]
!
!.... Sponge boundary treatments
!
!       0 = no sponge
!       1 = sponge with the forcing acoustic wave subtracted off
!       2 = sponge the disturbances to zero
!=============================================================================!
!     N o n L i n e a r   B o u n d a r y   C o n d i t i o n s
!=============================================================================!
!
!.... wall BC flag
!
!       0 = continuity is solved at the wall
!       1 = for each iteration drho/dn is set to zero and then corrected
!           by the wall normal momentum equation (nonlinear only)
!       2 = Wall normal momentum equation (nonlinear only)
!       3 = Extrapolation of density (third order)
!       4 = Poinsot & Lele continuity equation at the wall (not in LHS)
!       7 = Symmetry in finite difference stencil
!
!.... Wall temperature BC flag
!
!       0 = isothermal wall
!       1 = hard FD adiabatic
!       2 = soft adiabatic energy equation
!
!.... top BC flag
!
!       0 = First-order Riemann
!       1 = Poinsot & Lele (nonreflecting)
!       2 = Finite Wave
!       3 = Hard Freestream BC
!
!.... left and right BC
!
!       0 = First-order Riemann
!       1 = Poinsot & Lele (outflow)
!       2 = Symmetry plane (only on left boundary)
!       7 = Symmetry in finite difference stencil
!       8 = Hold IC [right boundary only]
!
!=============================================================================!
!       write(*,"(/,'Enter top, wall, left, right, wallt ==> ',$)")
        read(*,*) top, wall, left, right, wallt
        
        if (top.eq.7)   then
          if (impl.ne.0) call error('input$','top=7 on valid for impl=0$')
          tsym=.true.
        end if
          
        if (wall.eq.7) then
          if (impl.ne.0) call error('input$','wall=7 on valid for impl=0$')
          bsym=.true.
        end if
        if (left.eq.7)  lsym=.true.
        if (right.eq.7) rsym=.true.

        if (xper) then
          left = -1
          right = -1
          lsym = .false.
          rsym = .false.
        end if

        if (yper) then
          wall = -1
          wallt = -1
          top = -1
          bsym = .false.
          tsym = .false.
        end if
                
!.... set the explicit artificial dissipation coefficient

!       write(*,"(/,'Enter eps_e, ibuff ==> ',$)")
        read(*,*) eps_e, ibuff

!.... set the forcing frequency and virtual origin

!       write(*,"(/,'Enter omega, x0 ==> ',$)")
        read(*,*) omega, x0
        
!.... damping flag (parabolic flag)
!
!       0 = no damping
!       1 = damping
        
!       write(*,"(/,'Enter idamp ==> ',$)")
        read(*,*) idamp 
        if (idamp.ne.0 .and. linear.eq.1) then
          call error(code,'idamp <> 0 with linear = 1$')
        end if

!.... sponge parameters

!       write(*,"(/,'Enter As, Ns, xs, xt ==> ',$)")
        read(*,*) As, Ns, xs, xt
!       write(*,"(/,'Enter As2, Ns2, xs2, xt2 ==> ',$)")
        read(*,*) As2, Ns2, xs2, xt2
        
!.... extrapolation flag

!       write(*,"(/,'Enter extrap ==> ',$)")
        read(*,*) extrap

!.... Sweep angle

!       write(*,"(/,'Enter alpha (deg.) ==> ',$)")
        !read(*,*) alpha
        !alpha = alpha * pi / 180.0
!       write(*,"(/,'Enter theta (deg.) ==> ',$)")
        read(*,*) theta 

!.... traces and statistics

!       write(*,"(/,'Enter tflag, fflag, sflag ==> ',$)")
        read(*,*) tflag, fflag, sflag

!.... spanwise wave number

!       write(*,"(/,'Lz (0=2d) ==> ',$)")
        read(*,*) Lz
        
        if (Lz .ne. zero) then
          kz = two * pi / Lz      ! 3-d run
        else
          kz = zero               ! 2-d run
        end if
        
!.... Echo the run parameters

        open(10,file='echo.dat')
        write (10,*)
        write (10,*) 'R U N   P A R A M E T E R S'
        write (10,*)
        write (10,*) 'Mach    = ',Ma
        write (10,*) 'Re      = ',Re
        write (10,*) 'Pr      = ',Pr
        write (10,*) 'Tref    = ',Tref
        write (10,*) 'Delt    = ',Dtmax
        write (10,*) 'CFL     = ',cflmax
        write (10,*) 'niter   = ',niter
        write (10,*) 'nstep   = ',nstep
        write (10,*) 'optx    = ',optx
        write (10,*) 'opty    = ',opty
        write (10,*) 'xper    = ',xper
        write (10,*) 'yper    = ',yper
        write (10,*) 'ispg    = ',ispg
        write (10,*) 'linear  = ',linear
        write (10,*) 'ires    = ',ires
        write (10,*) 'impl    = ',impl
        write (10,*) 'loctime = ',loctime
        write (10,*) 'wall    = ',wall
        write (10,*) 'eps_e   = ',eps_e
        write (10,*) 'ibuff   = ',ibuff
        write (10,*) 'omega   = ',omega
        write (10,*) 'x0      = ',x0
        write (10,*) 'idamp   = ',idamp
        write (10,*) 'top     = ',top
        write (10,*) 'wall    = ',wall
        write (10,*) 'left    = ',left
        write (10,*) 'right   = ',right
        write (10,*) 'wallt   = ',wallt
        write (10,*) 'As      = ',As
        write (10,*) 'Ns      = ',Ns
        write (10,*) 'xs      = ',xs
        write (10,*) 'xt      = ',xt
        write (10,*) 'Extrap  = ',extrap
        !write (10,*) 'Alpha   = ',alpha
        write (10,*) 'Theta   = ',theta
        write (10,*) 'Lz      = ',Lz
        write (10,*) 'kz      = ',kz
        write (10,*)
        close (10)

!.... Sanity checks

        if (eps_e .ne. zero) then
          if (xper .or. yper) then
            call warning(code,'Periodic BC not supported with eps_e .ne. 0$')
          else if (tsym .or. bsym) then
            call warning(code,'tsym/bsym not supported with eps_e .ne. 0$')
          end if
        end if

!.... Determine whether to use complex or real analysis

        if ( kz .ne. zero .or. linear .eq. 2 ) then
          complex_analysis = .true.
        else
          complex_analysis = .false.
        end if
        if (linear .eq. 2) linear = 1

!.... Set generalized coordinates flag
        
        general = .true.

        return
        end
