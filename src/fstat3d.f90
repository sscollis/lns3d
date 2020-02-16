!=============================================================================!
        subroutine fstat3d(eig, vl, vml)
!=============================================================================!
!
!  This routine computes FFT of unsteady boundary statistics
!
!  Need to implement FFT routines on non Cray platforms
!
!=============================================================================!
        use global
        implicit none

        logical :: eig
        complex :: vl(ndof,nx,ny), vml(ndof,nx,ny)

!.... allocate space for FFT routines

        complex       :: rt(nx-1), ut(nx-1)
        complex       :: vt(nx-1), wt(nx-1)
        complex       :: tt(nx-1)
        real, save, allocatable :: table(:)     ! 100 + 8*(nx-1)
        real          :: work(8*(nx-1))
        integer, save :: ifft = 0
        integer       :: i, j, idof
        real          :: dy1, dy2, dum
        real          :: esum0, esum1, esum2, esum3

!=============================================================================!
#ifndef CRAY
        call error('fstat3d$','FFT routines only on Cray platforms$')
#endif

!.... initialize the table on the first call

        if (ifft .eq. 0) then
           write(*,*) 'Initializing FFT table...'
           allocate( table( 100 + 8*(nx-1) ) )
#ifdef CRAY
           call CCFFT(0, nx-1, 1.0/float(nx-1), dum, dum, table, work, 0)
#endif
           open(unit=35, file='stat.dat', form='formatted', &
                status='unknown', position='append')
           ifft = 1
        end if

        if (eig) then
          open(unit=10, file='wave.dat', form='unformatted', &
               status='unknown')
        end if

!.... loop over y

        do j = 1, ny

!.... take the FFT of rho, u, v, w, and t

#ifdef CRAY
           call CCFFT(1, nx-1, 1.0/float(nx-1), vl(1,:,j), rt, table, work, 0)
           call CCFFT(1, nx-1, 1.0/float(nx-1), vl(2,:,j), ut, table, work, 0)
           call CCFFT(1, nx-1, 1.0/float(nx-1), vl(3,:,j), vt, table, work, 0)
           call CCFFT(1, nx-1, 1.0/float(nx-1), vl(4,:,j), wt, table, work, 0)
           call CCFFT(1, nx-1, 1.0/float(nx-1), vl(5,:,j), tt, table, work, 0)
#endif

!.... write out the eigenfunctions

           if (eig) then
              write(10) y(1,j), rt, ut, vt, wt, tt
           end if

!.... integrate in the wall normal direction

           if (j .eq. 1) then
              dy2 = y(1,j+1) - y(1,j)
              esum0 = pt5 * dy2 * vml(1,1,j) *                          &
                  ( vml(2,1,j)**2 + vml(3,1,j)**2 + vml(4,1,j)**2 )
              esum1 = pt5 * dy2 * vml(1,1,j) *                          &
                  ( ut(2) * CONJG(ut(2)) +                              &
                    vt(2) * CONJG(vt(2)) + wt(2) * CONJG(wt(2)) )
              esum2 = pt5 * dy2 * ( tt(2) * CONJG(tt(2)) )
              esum3 = pt5 * dy2 * vml(1,1,j) *                          &
                  ( ut(3) * CONJG(ut(3)) +                              &
                    vt(3) * CONJG(vt(3)) )
           else if ( j .eq. ny) then
              dy1 = y(1,j) - y(1,j-1)
              esum0 = esum0 + pt5 * dy1 * vml(1,1,j) *                  &
                  ( vml(2,1,j)**2 + vml(3,1,j)**2 + vml(4,1,j)**2 )
              esum1 = esum1 + pt5 * dy1 * vml(1,1,j) *                  &
                  ( ut(2) * CONJG(ut(2)) +                              &
                    vt(2) * CONJG(vt(2)) + wt(2) * CONJG(wt(2)) )
              esum2 = esum2 + pt5 * dy1 * ( tt(2) * CONJG(tt(2)) )
              esum3 = esum3 + pt5 * dy1 * vml(1,1,j) *                  &
                  ( ut(3) * CONJG(ut(3)) +                              &
                    vt(3) * CONJG(vt(3)) )
           else
              dy1 = y(1,j) - y(1,j-1)
              dy2 = y(1,j+1) - y(1,j)
              esum0 = esum0 + pt5 * (dy1 + dy2) * vml(1,1,j) *          &
                  ( vml(2,1,j)**2 + vml(3,1,j)**2 + vml(4,1,j)**2 )
              esum1 = esum1 + pt5 * (dy1 + dy2) * vml(1,1,j) *          &
                  ( ut(2) * CONJG(ut(2)) +                              &
                    vt(2) * CONJG(vt(2)) + wt(2) * CONJG(wt(2)) )
              esum2 = esum2 + pt5 * (dy1 + dy2) * ( tt(2) * CONJG(tt(2)) )
              esum3 = esum3 + pt5 * (dy1 + dy2) * vml(1,1,j) *          &
                  ( ut(3) * CONJG(ut(3)) +                              &
                    vt(3) * CONJG(vt(3)) )
           end if

        end do

!.... write out the statistics

        write(35,10) time, esum0, esum1, esum2, esum3
        call flush(35)
 10     format(5(1pe20.13,1x))

        if (eig) then
          close(10)
        end if
          
        return
        end
