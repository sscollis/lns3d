!=============================================================================!
        subroutine resstat3D(res, dtl)
!  
!  This subroutine calculates the statistics of the residual.
!
!  Note that the norm of the residual is divided by dtl
!  
!=============================================================================!
        use global
        use material
        implicit none
        
        complex :: res(ndof,nx,ny)
        real    :: dtl(nx,ny)
        !$sgi distribute res(*,*,block), dtl(*,block)

        real    :: resnod
        real    :: totres = 0.0d0, resfrt = 0.0d0, resmax = 0.0d0
        real    :: resmx(ndof)
        integer :: jtotrs, jresmx, numnp
        integer :: i, j, ixrmax=0, iyrmax=0
        real    :: CPUl
        
#ifdef CRAY
        real, save :: CPU
        real, external :: second
#else
        real*4, save :: CPU
        real*4, external :: second
#endif
!=============================================================================!
        ixrmax = 0
        iyrmax = 0

        numnp = nx * ny
!       !$omp parallel do private(i,resnod,ixrmax,iyrmax) &
!       !$omp& reduction(+: totres) reduction(max: resmax)
        do j = 1, ny
          do i = 1, nx
            resnod = ( res(1,i,j) * conjg(res(1,i,j)) + &
                       res(2,i,j) * conjg(res(2,i,j)) + &
                       res(3,i,j) * conjg(res(3,i,j)) + &
                       res(4,i,j) * conjg(res(4,i,j)) + &
                       res(5,i,j) * conjg(res(5,i,j)) ) / dtl(i,j)**2
            totres = totres + resnod
!           resmax = max(resnod,resmax)
            if (resnod .gt. resmax) then
              resmax = resnod
              ixrmax = i
              iyrmax = j
            end if
          end do
        end do

        if (ixrmax.ne.0 .and. iyrmax.ne.0) then
          resmx(:) = abs(res(:,ixrmax,iyrmax)) / dtl(ixrmax,iyrmax)
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

        if (istep.eq.0) then 
          CPUl = zero
        else
          CPUl = second() - CPU
        end if
        CPU = second()

!.... output the result

        write(*,1000) lstep+1, time, CPU, CPUl, totres, jtotrs, &
                      ixrmax, iyrmax, jresmx

        write(ihist,1000) lstep+1, time, CPU, CPUl, totres, jtotrs, &
                          ixrmax, iyrmax, jresmx
        call flush(ihist)

!.... this outputs the value of the maximum residual

!       write(60,1010) lstep+1,  ixrmax,   iyrmax,  resmx(1), resmx(2), &
!                      resmx(3), resmx(4), resmx(5)
                       
!.... reinitialize the variables

        totres = zero
        resmax = zero

        return
        
1000    format(1p,i6,e10.3,e10.3,e10.3,e10.3,2x,'(',i4,')',&
               2x,'<',' (',i3,',',i3,') ','|',i4,'>')
1010    format(3(i6,1x),5(1pe13.6,1x))

        end

!=============================================================================!
        subroutine vstat3d(vl, dtl)
!  
!  This subroutine calculates the statistics of the solution.
!
!=============================================================================!
        use global
        use material
        implicit none
        
        complex :: vl(ndof,nx,ny)
        real    :: dtl(nx,ny)
        !$sgi distribute vl(*,*,block), dtl(*,block)

        real    :: resnod
        real    :: totres = 0.0d0, resfrt = 0.0d0, resmax = 0.0d0
        real    :: resmx(ndof)
        integer :: jtotrs, jresmx, numnp
        integer :: i, j, ixrmax, iyrmax
        real    :: CPUl
        
#ifdef CRAY
        real, save :: CPU
        real, external :: second
#else
        real*4, save :: CPU
        real*4, external :: second
#endif
!=============================================================================!
        ixrmax = 0
        iyrmax = 0

        numnp = nx * ny
!       !$omp parallel do private(i,resnod,ixrmax,iyrmax) &
!       !$omp& reduction(+: totres) reduction(max: resmax)
        do j = 1, ny
          do i = 1, nx
            resnod = vl(1,i,j) * conjg(vl(1,i,j)) + &
                     vl(2,i,j) * conjg(vl(2,i,j)) + &
                     vl(3,i,j) * conjg(vl(3,i,j)) + &
                     vl(4,i,j) * conjg(vl(4,i,j)) + &
                     vl(5,i,j) * conjg(vl(5,i,j))
            totres = totres + resnod
!           resmax = max(resnod,resmax)
            if (resnod .gt. resmax) then
              resmax = resnod
              ixrmax = i
              iyrmax = j
            end if
          end do
        end do

        if (ixrmax.ne.0 .and. iyrmax.ne.0) then
          resmx(:) = abs(vl(:,ixrmax,iyrmax)) / dtl(ixrmax,iyrmax)
        else
          resmx(:) = zero
        end if

        totres = sqrt(totres / real(numnp))
        resmax = sqrt(resmax)

        if (resfrt .eq. zero) resfrt = totres
        if (totres .ne. zero) then
          jtotrs = int ( 10.0 * log10 ( totres / resfrt ) )
          jresmx = int ( 10.0 * log10 ( resmax / totres ) )
        else
          jtotrs = zero
          jresmx = zero
        end if

!.... calculate the CPU-time

        if (istep.eq.0) then 
          CPUl = zero
        else
          CPUl = second() - CPU
        end if
        CPU  = second()

!.... output the result

        write(*,1000) lstep, time, Delt, cfl, totres, jtotrs, &
                      jresmx, CPU, CPUl

        write(ihist,1000) lstep, time, Delt, cfl, totres, jtotrs, &
                          jresmx, CPU, CPUl

!!$     write(*,1000) lstep, time, Delt, cfl, totres, jtotrs, &
!!$                      ixrmax, iyrmax, jresmx, CPU, CPUl
!!$
!!$     write(ihist,1000) lstep, time, Delt, cfl, totres, jtotrs, &
!!$                          ixrmax, iyrmax, jresmx, CPU, CPUl
        call flush(ihist)

!.... this outputs the value of the maximum residual

!       write(60,1010) lstep,  ixrmax,   iyrmax,  resmx(1), resmx(2), &
!                      resmx(3), resmx(4), resmx(5)
                       
!.... reinitialize the variables

        totres = zero
        resmax = zero

        return
        
1000    format(1p,i6,e10.3,e10.3,e10.3,e10.3,1x,'(',i4,') ',i4,e10.3,e10.3)
!!$1000    format(1p,i6,e10.3,e10.3,e10.3,e10.3,2x,'(',i4,')',&
!!$               2x,'<',' (',i3,',',i3,') ','|',i4,'>',e10.3,e10.3)
1010    format(3(i6,1x),5(1pe13.6,1x))

        end
