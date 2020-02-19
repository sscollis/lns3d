        module timemod
        
        real, allocatable :: volinv(:), sxj(:), syj(:), sxi(:), syi(:)
        real, allocatable :: acoustic(:), viscous(:)
        
        end module timemod
!=============================================================================!
        subroutine dtcfl(vl, dtl)
!  
!  Compute the CFL and time step
!  
!  Note: this isn't correct for three-dimensional flow
!
!  WARNING:  I think that the CFL is defined incorrectly!!!!!!!
!
!  Revised:  6-11-96
!=============================================================================!
        use stuff
        use timemod
        use material
        use local
        implicit none

        real :: vl(ny*nx,ndof), dtl(ny*nx)
        
        real :: cfl, dtlmax
        
        real    :: cflv, cfla, cflc
        integer :: iv, jv, ia, ja, ic, jc, nloc(1)
!=============================================================================!

!.... first time setup

        if (istep.eq.1 .or. linear.eq.1) then
          allocate( volinv(ny*nx), sxj(ny*nx), syj(ny*nx), sxi(ny*nx), &
                    acoustic(ny*nx), viscous(ny*nx), syi(ny*nx) )

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

        rho = vl(:,1)
        u1  = vl(:,2)
        u2  = vl(:,3)
        u3  = vl(:,4)
        t   = vl(:,5)
        rhoinv = one / rho

!.... get the material properties

        call getmat(t(:),                               &
                    mu(:),   lm(:),    con(:),          &
                    dmu(:),  d2mu(:),  dlm(:),          &
                    d2lm(:), dcon(:),  d2con(:)         )

        mu = max( mu / Re, (lm + mu) / Re, con / ( Re * Pr ) ) * rhoinv

!.... put sqrt(t) into t for use as temporary speed of sound

        t = sqrt(t)
        
!.... Total CFL

!       dtl = (abs( u1 * sxi + u2 * syi ) + &
!              abs( u1 * sxj + u2 * syj ) + &
!              t * acoustic + &
!              mu * viscous ) * volinv

!       cflv = maxval( dtl )
!       nloc = maxloc( dtl )
!       iv = (nloc(1)-1) / ny + 1
!       jv = nloc(1) - (iv-1)*ny

!.... convective CFL

!       dtl = (abs( u1 * sxi + u2 * syi ) + &
!              abs( u1 * sxj + u2 * syj ) ) * volinv

!       cflc = maxval( dtl )
!       nloc = maxloc( dtl )
!       ic = (nloc(1)-1) / ny + 1
!       jc = nloc(1) - (ic-1)*ny
                                
!.... acoustic CFL

        dtl = (abs( u1 * sxi + u2 * syi ) + &
               abs( u1 * sxj + u2 * syj ) + &
               t * acoustic ) * volinv

        cfla = maxval( dtl )
        nloc = maxloc( dtl )
        ia = (nloc(1)-1) / ny + 1
        ja = nloc(1) - (ia-1)*ny

        cfl = cfla              ! use the acoustic CFL for scaling
        
        if (cfl*dtmax .le. cflmax) then
          delt = dtmax
          dtl  = dtmax
        else
          delt = cflmax / cfl
          if (loctime.eq.0) then
            dtl = delt
          else
            dtl = cflmax / dtl
          end if
        end if

!.... WARNING:  When this is on, things can be crazy!
!.... enforce x symmetry (don't need this when there is no symmetry)

!       if (loctime.eq.1) then
!         do i = 1, nx/2
!           dtl(:,i) = dtl(:,nx-i+1)
!         end do
!       end if
                   
        if (istep.eq.1 .or. linear.eq.1) then
!         write(*,"(/,'Maximum Total CFL/dt      = ',1pe13.6, &
!               & ' at (',i3,',',i3,')')") cflv, iv, jv
          write(*,"  ('Maximum Acoustic CFL/dt   = ',1pe13.6, &
                & ' at (',i3,',',i3,')')") cfla, ia, ja
!         write(*,"  ('Maximum Convective CFL/dt = ',1pe13.6, &
!               & ' at (',i3,',',i3,')')") cflc, ic, jc
          write(*,"(/,'Delt = ',1pe13.6,'  CFLa = ',1pe13.6)") &
                delt, cfl*delt
!         write(*,  "('CFLa = ',1pe13.6,'  CFLc = ',1pe13.6,/)") &
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

        return
        end
