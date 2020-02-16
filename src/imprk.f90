!=============================================================================!
        subroutine impRK()
!  
!  Driver for implicit R-K time advancement of the compressible Navier-Stokes
!  solver.  This routine could also be used as the backbone of a MF-GMRES 
!  driver.  However, the RK scheme doesn't really work well enough to justify
!  further work right now.  The R-K scheme is the scheme proposed by Jameson
!  that is supposed to be quite dissipative.  The coefficients come from
!  Charles Fenno's thesis from NC State [1993].
!
!  Author:   Scott Collis
!
!  Written:  6-5-95
!
!=============================================================================!
        use global
        use local
        implicit none
!=============================================================================!
        real, allocatable :: v(:), vold(:), v2old(:), r(:)
        real, allocatable :: vint(:)
        
        real :: rtime, cpu, dtau
        integer :: i, j, idof
        integer :: ier, mem = 0
        integer, external :: igetver
        real, external :: second
!=============================================================================!

!.... allocate the storage area for the disturbance field

        mem = mem + ny*nx*ndof
        allocate (v(ny*nx*ndof), vold(ny*nx*ndof), STAT=ier)
        if (ier .ne. 0) call error('imprk$','Insufficient Memory for v$')

!.... generate the initial condition

        call genini(v, vold)

!.... allocate memory for the RHS

        mem = mem + ny*nx*ndof
        allocate (r(ny*nx*ndof), STAT=ier)
        if (ier .ne. 0) call error('impRK$','Insufficient Memory for r$')

!.... allocate memory for the temporal v

        mem = mem + ny*nx*ndof
        allocate (vint(ny*nx*ndof), STAT=ier)
        if (ier .ne. 0) call error('impRK$','Insufficient Memory for vint$')

!.... allocate memory for the previous, previous disturbance field

        if (impl.eq.2) then
          mem = mem + ny*nx*ndof
          allocate (v2old(ny*nx*ndof), STAT=ier)
          if (ier .ne. 0) call error('impRK$', &
                                     'Insufficient Memory for v2old$')
        end if

        write(*,"(' ImpRK allocated ===> ',1pe13.6,' words')") float(mem)
!=============================================================================!

!.... dtau is the time step in pseudotime.  I don't have a rational way
!.... setting it right now.

        dtau = 0.5
        
!.... open the history file

        iver = igetver( base, 'h'//char(0) ) + 1
        call addver(base, 'h'//char(0), iver, filen, lfile)
        open(unit=ihist, file=filen, form='formatted', status='unknown')

!.... compute the CFL

        call set_time_step(vm)

!.... loop through the time steps

        do istep = 1, nstep
          
!.... -----------------------> predictor phase <-----------------------

          if (lstep .gt. 0) v2old = vold
          vold  = v
          time  = time + Delt
          call itrbc(v,vm)
          
!.... --------------------> multi-corrector phase <--------------------

        do iter = 1, niter
        
!.... Two-step backward

            if ( lstep .gt. 0 ) then

              call lrhs(v, r, vm, x, y)
              r = v - onept33 * vold + pt33 * v2old +           &
                  pt66 * delt * r
              call rkbc(r)
              call resstat(r)

              vint = v - 0.25 * dtau * r
              call itrbc(vint,vm)
              call lrhs(vint, r, vm, x, y)
              r = vint - onept33 * vold + pt33 * v2old +        &
                  pt66 * delt * r
              
              vint = v - 0.50 * dtau * r
              call itrbc(vint,vm)
              call lrhs(vint, r, vm, x, y)
              r = vint - onept33 * vold + pt33 * v2old +        &
                  pt66 * delt * r

              vint = v - 0.55 * dtau * r
              call itrbc(vint,vm)
              call lrhs(vint, r, vm, x, y)
              r = vint - onept33 * vold + pt33 * v2old +        &
                  pt66 * delt * r

            else

!.... Backward Euler

              call lrhs( v, r, vm, x, y)
              r = v - vold + delt * r
              call rkbc(r)
              call resstat(r)

              vint = v - 0.25 * dtau * r
              call itrbc(vint,vm)
              call lrhs( vint, r, vm, x, y)
              r = vint - vold + delt * r
              
              vint = v - 0.50 * dtau * r
              call itrbc(vint,vm)
              call lrhs( vint, r, vm, x, y)
              r = vint - vold + delt * r

              vint = v - 0.55 * dtau * r
              call itrbc(vint,vm)
              call lrhs( vint, r, vm, x, y)
              r = vint - vold + delt * r

            end if
                    
!.... update v

            v = v - dtau * r
              
!.... satisfy boundary conditions on v

            call itrbc(v,vm)
            
          end do                ! loop on iter

          lstep = lstep + 1

!.... quit if less than 10 seconds remain

#ifdef CRAY
          call tremain(rtime)
          if (rtime .le. 10) goto 5000
#endif

!.... compute temporal statistics

          if(.true. .and. xper .and. istep .ne. nstep .and. &
             mod(istep-1,itout) .eq. 0 ) then
            call fstat2d(.false., v, vm)        ! don't print profiles
          end if

!.... compute spatial statistics

          if (.false. .and. linear.eq.1) then
            call sstat(v, vm, x, y)
          end if

!.... write out restart and q files 

          if ( mod(lstep,ntout) .eq. 0 .and. istep .ne. nstep) then
            
            iver = igetver( base, 'R'//char(0) ) + 1
            call addver(base, 'R'//char(0), iver, filen, lfile)
            open(unit=ifile, file=filen, form='unformatted', &
                  status='unknown')
            write(ifile) lstep, time, nx, ny, nz, ndof
            write(ifile) v
            if (impl.eq.2) write(ifile) vold
            close(ifile)

            write(*,10) filen(1:index(filen,char(0))-1), lstep, time
            write(ihist,10) filen(1:index(filen,char(0))-1), lstep, time

          end if
          
        end do          ! loop on istep

!.... ---------------------->  Post Processing  <----------------------

5000    continue

!.... compute final statistics

        if (.true. .and. xper) then
          call fstat2d(.true., v, vm)   ! print profiles
        end if

!.... write the restart file

        iver = igetver( base, 'R'//char(0) ) + 1
        call addver(base, 'R'//char(0), iver, filen, lfile)
        open(unit=ifile, file=filen, form='unformatted', &
             status='unknown')
        write(ifile) lstep, time, nx, ny, nz, ndof
        write(ifile) v
        if (impl.eq.2) write(ifile) vold
        close(ifile)

        write(*,10) filen(1:index(filen,char(0))-1), lstep, time
        write(ihist,10) filen(1:index(filen,char(0))-1), lstep, time

!.... close the history file

        close(ihist)

        return
        
  10    format(' Wrote ',a,' at step = ',i7,' time = ',1pe13.6)

        end
