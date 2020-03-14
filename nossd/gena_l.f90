!=======================================================================================================!
        subroutine genA_l(A, vl)
!
!  Form the A matrix which multiplies U,x
!
!  Revised: 10-16-95
!
!  Implements the Lele & Poinsot boundary conditions for the linear equations.
!=======================================================================================================!
        use global
        use local
        implicit none
        
        real :: A(ny*nx,ndof,ndof), vl(ny*nx,ndof)
        real :: fact1, fact2, fact3, fact4, gmsinv

        real :: pinf, cinf, Lx, rk, c
        real :: L1, L2, L3, L4, L5
        real :: d1, d2, d3, d4, d5, dc
        
        integer :: i, j, ji

        real :: bn1, bn2, bs1, bs2
        real :: un, us, gpn, gps, grhon, grhos, gtn, gts
        real :: gu1s, gu1n, gu2s, gu2n, gu3s, gu3n
        real :: gunn, gusn, guns, guss
!=======================================================================================================!

!.... compute some extra stuff needed for Lele's BC's

        pinf = one / (gamma * Ma**2)    ! pressure at infinity
        cinf = one / Ma                 ! speed of sound at infinity
        Lx   = 100.0                    ! length of domain in x
        
        rk = sigma * ( one - Ma**2 ) * cinf / Lx

!=======================================================================================================!
!       W a l l   B o u n d a r y
!=======================================================================================================!
        if (wall.eq.4) then

        j = 1
        do i = 1, nx
          ji = j + (i-1) * ny
    
!.... get the boundary normal and tangent vectors

            bn1 = n1(ji) / sqrt( n1(ji)**2 + n2(ji)**2 )
            bn2 = n2(ji) / sqrt( n1(ji)**2 + n2(ji)**2 )

            bs1 = -bn2
            bs2 =  bn1

!.... compute Us and Un
        
            un = bn1 * vl(ji,2) + bn2 * vl(ji,3)
            us = bs1 * vl(ji,2) + bs2 * vl(ji,3)

            c = sqrt( t(ji) ) / Ma
            dc = pt5 / (Ma**2 * c)

!.... compute some derivatives in the boundary normal coordinate system

            gpn = bn1 * gp(ji,1) + bn2 * gp(ji,2)
            gps = bs1 * gp(ji,1) + bs2 * gp(ji,2)
        
            grhon = bn1 * grho(ji,1) + bn2 * grho(ji,2)
            grhos = bs1 * grho(ji,1) + bs2 * grho(ji,2)

            gtn = bn1 * gt(ji,1) + bn2 * gt(ji,2)
            gts = bs1 * gt(ji,1) + bs2 * gt(ji,2)

            gu1n = bn1 * gu(ji,1,1) + bn2 * gu(ji,1,2)
            gu1s = bs1 * gu(ji,1,1) + bs2 * gu(ji,1,2)

            gu2n = bn1 * gu(ji,2,1) + bn2 * gu(ji,2,2)
            gu2s = bs1 * gu(ji,2,1) + bs2 * gu(ji,2,2)

            gu3n = bn1 * gu(ji,3,1) + bn2 * gu(ji,3,2)
            gu3s = bs1 * gu(ji,3,1) + bs2 * gu(ji,3,2)

            gunn = bn1*bn1 * gu(ji,1,1) + bn2*bn1 * gu(ji,2,1) + &
                   bn1*bn2 * gu(ji,1,2) + bn2*bn2 * gu(ji,2,2)

            gusn = bs1*bn1 * gu(ji,1,1) + bs2*bn1 * gu(ji,2,1) + &
                   bs1*bn2 * gu(ji,1,2) + bs2*bn2 * gu(ji,2,2)

            guns = bn1*bs1 * gu(ji,1,1) + bn2*bs1 * gu(ji,2,1) + &
                   bn1*bs2 * gu(ji,1,2) + bn2*bs2 * gu(ji,2,2)

            guss = bs1*bs1 * gu(ji,1,1) + bs2*bs1 * gu(ji,2,1) + &
                   bs1*bs2 * gu(ji,1,2) + bs2*bs2 * gu(ji,2,2)
   
!.... Continuity equation  (only need continuity at the wall)

            gmsinv = one / (gamma * Ma**2)
    
            A(ji,1,1) = one/c**2*(un-c)*bn1*t(ji)*gmsinv + us*bs1
            A(ji,1,2) = -one/c**2*(un-c)*rho(ji)*c*bn1*bn1 + rho(ji)*bs1*bs1
            A(ji,1,3) = -one/c**2*(un-c)*rho(ji)*c*bn2*bn1 + rho(ji)*bs2*bs1
            A(ji,1,4) = zero
            A(ji,1,5) = one/c**2*(un-c)*bn1*rho(ji)*gmsinv

        end do
        
        end if          ! wall 

!=======================================================================================================!
!       S l i p   W a l l   B o u n d a r y
!=======================================================================================================!
        if (left.eq.2) then

        i = 1
        do j = 1, ny
          ji = j + (i-1) * ny
    
!.... get the boundary normal and tangent vectors

            bn1 = m1(ji) / sqrt( m1(ji)**2 + m2(ji)**2 )
            bn2 = m2(ji) / sqrt( m1(ji)**2 + m2(ji)**2 )

            bs1 = -bn2
            bs2 =  bn1

!.... compute Us and Un
        
            un = bn1 * vl(ji,2) + bn2 * vl(ji,3)
            us = bs1 * vl(ji,2) + bs2 * vl(ji,3)

            c = sqrt( t(ji) ) / Ma
            dc = pt5 / (Ma**2 * c)

!.... compute some derivatives in the boundary normal coordinate system

            gpn = bn1 * gp(ji,1) + bn2 * gp(ji,2)
            gps = bs1 * gp(ji,1) + bs2 * gp(ji,2)
        
            grhon = bn1 * grho(ji,1) + bn2 * grho(ji,2)
            grhos = bs1 * grho(ji,1) + bs2 * grho(ji,2)

            gtn = bn1 * gt(ji,1) + bn2 * gt(ji,2)
            gts = bs1 * gt(ji,1) + bs2 * gt(ji,2)

            gu1n = bn1 * gu(ji,1,1) + bn2 * gu(ji,1,2)
            gu1s = bs1 * gu(ji,1,1) + bs2 * gu(ji,1,2)

            gu2n = bn1 * gu(ji,2,1) + bn2 * gu(ji,2,2)
            gu2s = bs1 * gu(ji,2,1) + bs2 * gu(ji,2,2)

            gu3n = bn1 * gu(ji,3,1) + bn2 * gu(ji,3,2)
            gu3s = bs1 * gu(ji,3,1) + bs2 * gu(ji,3,2)

            gunn = bn1*bn1 * gu(ji,1,1) + bn2*bn1 * gu(ji,2,1) + &
                   bn1*bn2 * gu(ji,1,2) + bn2*bn2 * gu(ji,2,2)

            gusn = bs1*bn1 * gu(ji,1,1) + bs2*bn1 * gu(ji,2,1) + &
                   bs1*bn2 * gu(ji,1,2) + bs2*bn2 * gu(ji,2,2)

            guns = bn1*bs1 * gu(ji,1,1) + bn2*bs1 * gu(ji,2,1) + &
                   bn1*bs2 * gu(ji,1,2) + bn2*bs2 * gu(ji,2,2)

            guss = bs1*bs1 * gu(ji,1,1) + bs2*bs1 * gu(ji,2,1) + &
                   bs1*bs2 * gu(ji,1,2) + bs2*bs2 * gu(ji,2,2)
   
!.... Continuity equation (only need continuity at the wall)

            gmsinv = one / (gamma * Ma**2)
    
            A(ji,1,1) = one/c**2*(un-c)*bn1*t(ji)*gmsinv + us*bs1
            A(ji,1,2) = -one/c**2*(un-c)*rho(ji)*c*bn1*bn1 + rho(ji)*bs1*bs1
            A(ji,1,3) = -one/c**2*(un-c)*rho(ji)*c*bn2*bn1 + rho(ji)*bs2*bs1
            A(ji,1,4) = zero
            A(ji,1,5) = one/c**2*(un-c)*bn1*rho(ji)*gmsinv

        end do
        
        end if          ! slip wall

!=======================================================================================================!
!       S i d e   B o u n d a r i e s
!=======================================================================================================!
        if (left.eq.6 .or. right.eq.6) then

!.... do for both the left and right boundaries

        do i = 1, nx, nx-1
          if ( (i.eq.1 .and. left.eq.6) .or. (i.eq.nx .and. right.eq.6) ) then    
          do j = 1, ny
            ji = j + (i-1) * ny
    
            c = sqrt( t(ji) ) / Ma
            dc = pt5 / (Ma**2 * c)

!.... compute the Characteristic amplitudes
        
            L1 = ( u1(ji) - c ) * ( gp(ji,1) - rho(ji) * c * gu(ji,1,1) )
            L2 = u1(ji) * ( c**2 * grho(ji,1) - gp(ji,1) )
            L3 = u1(ji) * gu(ji,2,1)
            L4 = u1(ji) * gu(ji,3,1)
            L5 = ( u1(ji) + c ) * ( gp(ji,1) + rho(ji) * c * gu(ji,1,1) )
        
!.... nonreflecting outflow boundary condition
        
            L1 = rk * ( p(ji) - pinf )  ! nonreflecting condition

!.... compute the streamwise derivative terms

            d1 = one / c**2 * ( L2 + pt5 * ( L5 + L1 ) )
            d2 = pt5 * ( L5 + L1 )
            d3 = pt5 / ( rho(ji) * c ) * ( L5 - L1 )
            d4 = L3
            d5 = L4

!=======================================================================================================!

!.... Continuity equation

        gmsinv = one / (gamma * Ma**2)

        A(ji,1,1) = one/c**2*(u1(ji)*(c**2-t(ji)*gmsinv)+               &
                    pt5*(u1(ji)+c)*t(ji)*gmsinv)
        A(ji,1,2) = one/c**2*pt5*(u1(ji)+c)*rho(ji)*c
        A(ji,1,3) = zero
        A(ji,1,4) = zero
        A(ji,1,5) = one/c**2*( -u1(ji)*gmsinv*rho(ji) +                 &
                    pt5*(u1(ji)+c)*rho(ji)*gmsinv )

!=======================================================================================================!

!.... Momentum equation -- x_1

        fact1  = rhoinv(ji) / Re

        A(ji,2,1) = pt5*rhoinv(ji)/c*(u1(ji)+c)*t(ji)*gmsinv
        A(ji,2,2) = pt5*(u1(ji)+c) - fact1 * g1lm(ji) - fact1 * two * g1mu(ji)
        A(ji,2,3) = -fact1 * g2mu(ji)
        A(ji,2,4) = -fact1 * g3mu(ji)
        A(ji,2,5) = one/(two*rho(ji)*c)*(u1(ji)+c)*rho(ji)*gmsinv - &
                    fact1 * dlm(ji) * ( gu(ji,1,1) + gu(ji,2,2) + gu(ji,3,3) ) - &
                    fact1 * dmu(ji) * ( gu(ji,1,1) + gu(ji,1,1) )

!=======================================================================================================!

!.... Energy equation

        fact1 = gamma * rhoinv(ji) / (Pr * Re)
        fact2 = gamma * gamma1 * Ma**2 * rhoinv(ji) / Re

        A(ji,5,1) = -Ma**2 * rhoinv(ji) * ( u1(ji) * (c**2 - gmsinv * t(ji)) - &
                     pt5*gamma1*( (u1(ji)+c)*gmsinv*t(ji) ) )
        A(ji,5,2) =  Ma**2*rhoinv(ji)*(pt5*gamma1*(u1(ji)+c)*rho(ji)*c) - &
                     fact2 * two * lm(ji) * divu(ji) - &
                     two * fact2 * mu(ji) * ( gu(ji,1,1) + gu(ji,1,1) )
        A(ji,5,3) = -two * fact2 * mu(ji) * ( gu(ji,2,1) + gu(ji,1,2) )
        A(ji,5,4) = -two * fact2 * mu(ji) * ( gu(ji,3,1) + gu(ji,1,3) )
        A(ji,5,5) = -Ma**2*rhoinv(ji)*( (-u1(ji))*rho(ji)*gmsinv - &
                     pt5*gamma1*(u1(ji)+c)*rho(ji)*gmsinv ) - fact1 * (g1con(ji) + dcon(ji) * gt(ji,1))

          end do
        end if
        end do
        end if

        return
        end
