!=============================================================================!
        subroutine init
!=============================================================================!
!
!  Initialize variables
!
!=============================================================================!
        use global
        use buff_mod
        use ic
        implicit none

!.... stuff

        tflag = .false.   ! time traces
        fflag = .false.   ! temporal statistics
        sflag = .false.   ! spatial statistics

        top    = -1; wall=-1; left=-1; right=-1; wallt=-1
        extrap = -1
        lsym   = .false.; rsym=.false.; bsym=.false.; tsym=.false.
        xper   =.false.; yper=.false.
        Navier = .true.

        sigma  = zero

        gamma  = 1.4
        gamma1 = 0.4
        cv     = 716.5
        cp     = 1003.1
        Rgas   = 286.6

        ny = 1; nx = 1; nz = 1; nsd = 3; ndof = 5

        optx = 0
        opty = 0

        Delt = zero; time = zero; nstep = 1; niter = 1; lstep = 0

        iver = 0

!.... buff_stuff

        ibuff = 0

!.... warning the IC modules aren't initialized

        return
        end
