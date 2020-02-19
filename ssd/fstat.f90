!=============================================================================!
        subroutine fstat(eig, vl, vml)
!=============================================================================!
!
!  This routine computes FFT of unsteady boundary statistics
!
!=============================================================================!
        use stuff
        implicit none

        logical :: eig
        real    :: vl(ny,nx,ndof), vml(ny,nx,ndof)

!.... allocate space for FFT routines

        complex       :: rt((nx-1)/2+1), ut((nx-1)/2+1)
        complex       :: vt((nx-1)/2+1), wt((nx-1)/2+1)
        complex       :: tt((nx-1)/2+1)
        real, save, allocatable :: table(:) ! 100 + 4*(nx-1))
        real          :: work(4 + 4*(nx-1))
        integer, save :: ifft = 0
        integer       :: i, j, idof
        real          :: dy1, dy2, dum
        real          :: esum0, esum1, esum2, esum3

!.... initialize the table on the first call

        if (ifft .eq. 0) then
           write(*,*) 'Initializing FFT table...'
           allocate( table(100+4*(nx-1)) )
           call SCFFT(0, nx-1, 1.0/float(nx-1), dum, dum, table, work, 0)
           open(unit=35, file='stat.dat', form='formatted', &
                status='unknown', position='append')
           ifft = 1
        end if

        if (eig) then
          open(unit=10, file='profile.dat', form='unformatted', &
               status='unknown')
        end if

!.... loop over y

        do j = 1, ny

!.... take the FFT of rho, u, v, w, and t

           call SCFFT(1, nx-1, 1.0/float(nx-1), vl(j,:,1), rt, table, work, 0)
           call SCFFT(1, nx-1, 1.0/float(nx-1), vl(j,:,2), ut, table, work, 0)
           call SCFFT(1, nx-1, 1.0/float(nx-1), vl(j,:,3), vt, table, work, 0)
           call SCFFT(1, nx-1, 1.0/float(nx-1), vl(j,:,4), wt, table, work, 0)
           call SCFFT(1, nx-1, 1.0/float(nx-1), vl(j,:,5), tt, table, work, 0)

!.... write out the eigenfunctions

           if (eig) then
              write(10) y(j,1), rt, ut, vt, wt, tt
           end if

!.... integrate in the wall normal direction

           if (j .eq. 1) then
              dy2 = y(j+1,1) - y(j,1)
              esum0 = pt5 * dy2 * vml(j,1,1) *                  &
                  ( vml(j,1,2)**2 + vml(j,1,3)**2 )
              esum1 = pt5 * dy2 * vml(j,1,1) *                  &
                  ( ut(2) * CONJG(ut(2)) +                      &
                    vt(2) * CONJG(vt(2)) )
              esum2 = pt5 * dy2 *                               &
                  ( tt(2) * CONJG(tt(2)) )
              esum3 = pt5 * dy2 * vml(j,1,1) *                  &
                  ( ut(3) * CONJG(ut(3)) +                      &
                    vt(3) * CONJG(vt(3)) )
           else if ( j .eq. ny) then
              dy1 = y(j,1) - y(j-1,1)
              esum0 = esum0 + pt5 * dy1 * vml(j,1,1) *          &
                  ( vml(j,1,2)**2 + vml(j,1,3)**2 )
              esum1 = esum1 + pt5 * dy1 * vml(j,1,1) *          &
                  ( ut(2) * CONJG(ut(2)) +                      &
                    vt(2) * CONJG(vt(2)) )
              esum2 = esum2 + pt5 * dy1 *                       &
                  ( tt(2) * CONJG(tt(2)) )
              esum3 = esum3 + pt5 * dy1 * vml(j,1,1) *          &
                  ( ut(3) * CONJG(ut(3)) +                      &
                    vt(3) * CONJG(vt(3)) )
           else
              dy1 = y(j,1) - y(j-1,1)
              dy2 = y(j+1,1) - y(j,1)
              esum0 = esum0 + pt5 * (dy1 + dy2) * vml(j,1,1) *  &
                  ( vml(j,1,2)**2 + vml(j,1,3)**2 )
              esum1 = esum1 + pt5 * (dy1 + dy2) * vml(j,1,1) *  &
                  ( ut(2) * CONJG(ut(2)) +                      &
                    vt(2) * CONJG(vt(2)) )
              esum2 = esum2 + pt5 * (dy1 + dy2) *               &
                  ( tt(2) * CONJG(tt(2)) )
              esum3 = esum3 + pt5 * (dy1 + dy2) * vml(j,1,1) *  &
                  ( ut(3) * CONJG(ut(3)) +                      &
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
