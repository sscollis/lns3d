!=============================================================================!
        subroutine resstat(res, dtl)
!  
!  This subroutine calculates the statistics of the residual.
!
!  Note that the norm of the residual is divided by dtl
!  
!=============================================================================!
        use global
        use material
        implicit none
        
        real :: res(ny,nx,ndof), dtl(ny,nx)

        real    :: resnod
        real    :: totres = 0.0d0, resfrt = 0.0d0, resmax = 0.0d0
        real    :: resmx(ndof)
        integer :: jtotrs, jresmx, numnp
        integer :: i, j, ixrmax, iyrmax
        
#ifdef CRAY
        real :: CPU
        real, external :: second
#else
        real*4 :: CPU
        real*4, external :: second
#endif
!=============================================================================!
        ixrmax = 0
        iyrmax = 0

        numnp = nx * ny
        do i = 1, nx
          do j = 1, ny
            resnod = ( res(j,i,1)**2 + &
                       res(j,i,2)**2 + &
                       res(j,i,3)**2 + &
                       res(j,i,4)**2 + &
                       res(j,i,5)**2 ) / dtl(j,i)**2
            totres = totres + resnod
            if (resnod .gt. resmax) then
              resmax = resnod
              ixrmax = i
              iyrmax = j
            end if
          end do
        end do

        if (ixrmax.ne.0 .and. iyrmax.ne.0) then
          resmx(:) = res(ixrmax,iyrmax,:) / dtl(iyrmax,ixrmax)
        else
          resmx(:) = zero
        end if

        totres = sqrt(totres / real(numnp))
        resmax = sqrt(resmax)

        if (resfrt .eq. zero) resfrt = totres
        if (totres .ne. zero) then
          jtotrs = int ( 10.d0 * log10 ( totres / resfrt ) )
          jresmx = int ( 10.d0 * log10 ( resmax / totres ) )
        else
          jtotrs = zero
          jresmx = zero
        end if

!.... calculate the CPU-time

        CPU = second()

!.... output the result

        write(*,1000) lstep+1, time, CPU, totres, jtotrs, &
                      ixrmax, iyrmax, jresmx

        write(ihist,1000) lstep+1, time, CPU, totres, jtotrs, &
                          ixrmax, iyrmax, jresmx
        call flush(ihist)

!.... this outputs the value of the maximum residual

!       write(60,1010) lstep+1,  ixrmax,   iyrmax,  resmx(1), resmx(2), &
!                      resmx(3), resmx(4), resmx(5)
                       
!.... reinitialize the variables

        totres = zero
        resmax = zero

        return
        
1000    format(1p,i6,e10.3,e10.3,e10.3,2x,'(',i4,')',&
               2x,'<',' (',i3,',',i3,') ','|',i4,'>')
1010    format(3(i6,1x),5(1pe13.6,1x))

        end

