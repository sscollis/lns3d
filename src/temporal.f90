        module temporal
        
        real, allocatable :: volinv(:,:), sxj(:,:), syj(:,:), sxi(:,:)
        real, allocatable :: acoustic(:,:), viscous(:,:), syi(:,:)
        
        !$sgi distribute volinv(*,block), sxj(*,block), syj(*,block)
        !$sgi distribute sxi(*,block), syi(*,block)
        !$sgi distribute acoustic(*,block), viscous(*,block)
        
        end module temporal
!=============================================================================!
        subroutine set_time_step(vl, dtl)
!  
!  Compute the CFL and time step
!  
!  Note: this isn't correct for three-dimensional flow
!
!  Revised:  6-11-96
!            11-10-99  SGI parallel
!=============================================================================!
        use global
        use temporal
        use material
        use local2d
        implicit none

        integer :: i, j
        real :: vl(ndof,nx,ny), dtl(nx,ny)
        !$sgi distribute vl(*,*,block), dtl(*,block)
        
        real :: dtlmax
        
        real    :: cflv, cfla, cflc
        integer :: iv, jv, ia, ja, ic, jc, nloc(2)
!=============================================================================!

        if (.false.) then  ! simple constant time-stepping

          Delt = dtmax
          !$doacross local(i)
          !$omp parallel do private(i)
          do j = 1, ny
            do i = 1, nx
              dtl(i,j) = Delt
            end do
          end do

        else

!.... first time setup

        if (istep.eq.0 .or. linear.eq.1) then
          allocate( volinv(nx,ny), sxj(nx,ny), syj(nx,ny), sxi(nx,ny), &
                    acoustic(nx,ny), viscous(nx,ny), syi(nx,ny) )

!.... temporarily one/detJ

          volinv = one / (m1 * n2 - m2 * n1)
          
          sxj =  n2 * dxi * volinv
          syj = -m2 * dxi * volinv
          
          sxi = -n1 * deta * volinv
          syi =  m1 * deta * volinv

          volinv = one / ( dxi * deta * volinv )
          
          acoustic = sqrt( sxi**2 + two * sxi * syi + syi**2 + &
                           sxj**2 + two * sxj * syj + syj**2 ) / Ma
                           
          viscous = two * ( sxi**2 + two * sxi * syi + syi**2 + &
                    two * (sxi * sxj + sxi * syj + syi * sxj + syi * syj) + &
                    sxj**2 + two * sxj * syj + syj**2 ) * volinv
        end if
        
!.... initialize

        !$omp parallel do private( i, rho, u1, u2, u3, t, rhoinv, mu, lm, &
        !$omp& con, dmu, d2mu, dlm, d2lm, dcon, d2con )
        do j = 1, ny
          do i = 1, nx
            rho = vl(1,i,j)
            u1  = vl(2,i,j)
            u2  = vl(3,i,j)
            u3  = vl(4,i,j)
            t   = vl(5,i,j)
            rhoinv = one / rho

!.... get the material properties

            call getmat(t, mu,   lm,    con,            &
                           dmu,  d2mu,  dlm,            &
                           d2lm, dcon,  d2con           )

            mu = max( mu / Re, (lm + mu) / Re, con / ( Re * Pr ) ) * rhoinv

!.... put sqrt(t) into t for use as temporary speed of sound

            t = sqrt(t)
        
!.... Total CFL

!           dtl(i,j) = (abs( u1 * sxi(i,j) + u2 * syi(i,j) ) + &
!                       abs( u1 * sxj(i,j) + u2 * syj(i,j) ) + &
!                       t * acoustic(i,j) + &
!                       mu * viscous(i,j) ) * volinv(i,j)

!.... convective CFL

!           dtl(i,j) = (abs( u1 * sxi(i,j) + u2 * syi(i,j) ) + &
!                       abs( u1 * sxj(i,j) + u2 * syj(i,j) ) ) * volinv(i,j)

!.... acoustic CFL

            dtl(i,j) = (abs( u1 * sxi(i,j) + u2 * syi(i,j) ) + &
                        abs( u1 * sxj(i,j) + u2 * syj(i,j) ) + &
                        t * acoustic(i,j) ) * volinv(i,j)

!.... These are other simplier expressions for 1-d flows

!           dtl(i,j) = ( abs( u1 ) +  t / Ma ) * dxi / m1(i,j)
!           dtl(i,j) = ( t / Ma ) * dxi / m1(i,j)

          end do
        end do

        cfla = maxval( dtl )
        nloc = maxloc( dtl )
        ia = nloc(1)
        ja = nloc(2)

        cfl = cfla
        
        if (cfl*dtmax .le. cflmax) then
          delt = dtmax
          !$doacross local(i)
          !$omp parallel do private(i)
          do j = 1, ny
            do i = 1, nx
              dtl(i,j) = dtmax
            end do
          end do
        else
          delt = cflmax / cfl
          if (loctime.eq.0) then
            !$doacross local(i)
            !$omp parallel do private(i)
            do j = 1, ny
              do i = 1, nx
                dtl(i,j) = delt
              end do
            end do
          else
            !$doacross local(i)
            !$omp parallel do private(i)
            do j = 1, ny
              do i = 1, nx
                dtl(i,j) = cflmax / dtl(i,j)
              end do
            end do
          end if
        end if

        cfl = cfl * delt

!.... WARNING:  When this is on, things can be crazy!
!.... enforce x symmetry (don't need this when there is no symmetry)

!       if (loctime.eq.1) then
!         do i = 1, nx/2
!           dtl(:,i) = dtl(:,nx-i+1)
!         end do
!       end if
                   
        if (istep.eq.0 .or. linear.eq.1) then
!         write(*,"(/,' Maximum Total CFL/dt      = ',1pe13.6, &
!               & ' at (',i3,',',i3,')')") cflv, iv, jv
          write(*,"  (' Maximum Acoustic CFL/dt   = ',1pe13.6, &
                & ' at (',i3,',',i3,')')") cfla, ia, ja
!         write(*,"  (' Maximum Convective CFL/dt = ',1pe13.6, &
!               & ' at (',i3,',',i3,')')") cflc, ic, jc
          write(*,"(/,' Delt = ',1pe13.6,',  CFLa = ',1pe13.6,/)") &
                delt, cfl
!         write(*,  "(' CFLa = ',1pe13.6,',  CFLc = ',1pe13.6,/)") &
!               cfla*delt, cflc*delt
        end if

!.... output the maximum dt

!       if (loctime.eq.1) then
!         dtlmax = maxval( dtl )
!         nloc = maxloc( dtl )
!         ia = (nloc(1)-1) / ny + 1
!         ja = nloc(1) - (ia-1)*ny
!         write(80,"(i5,1x,2(1pe13.6,1x,i5,1x,i5,1x))") istep, delt, &
!                 iv, jv, dtlmax, ia, ja
!       end if

!.... the linear solver doesn't need this stuff anymore

        if (linear.eq.1) then
          deallocate( volinv, sxj, syj, sxi, syi, acoustic, viscous)
        end if

        end if

        return
        end
