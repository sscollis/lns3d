!=======================================================================================================!
        subroutine genD_l(D, vl)
!
!  Form the D matrix which multiplies U
!
!  Revised: 10-15-95
!
!  This version has been parabolized in x and implements the Lele & Poinsot
!  boundary conditions for the linear equations.
!=======================================================================================================!
        use global
        use local
        implicit none
        
        real :: D(ny*nx,ndof,ndof), vl(ny*nx,ndof)
        real :: fact1, fact2, fact3, fact4, gmsinv

        real :: pinf, cinf, Lx, rk, c
        real :: L1, L2, L3, L4, L5
        real :: d1, d2, d3, d4, d5, dc
        
        integer :: i, j, ji

        real :: bn1, bn2, bs1, bs2
        real :: un, us, gpn, gps, grhon, grhos, gtn, gts
        real :: gu1s, gu1n, gu2s, gu2n, gu3s, gu3n
        real :: gunn, gusn, guns, guss

        real :: b, rbeta, ub1, ub2, gub11, gub12, gub21, gub22
        real :: ubn, ubs, gubnn, gubsn, tb, rhob, pb, cb, gtbn, grhobn, gpbn
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
            
!.... compute the Characteristic amplitudes using the interior values
        
            L1 = ( un - c ) * ( gpn - rho(ji) * c * gunn )
            L2 = un * ( c**2 * grhon - gpn )
            L3 = un * gusn
            L4 = un * gu3n
            L5 = ( un + c ) * ( gpn + rho(ji) * c * gunn )

!.... Boundary conditions

            L2 = zero
            L3 = zero
            L4 = zero
            L5 = L1

!.... Continuity equation

            gmsinv = one / (gamma * Ma**2)
    
            D(ji,1,1) = one/c**2*(un-c)*(gtn*gmsinv - c*gunn) + guss + gu(ji,3,3)
            D(ji,1,2) = one/c**2*bn1*(gpn - rho(ji)*c*gunn) + grhos*bs1
            D(ji,1,3) = one/c**2*bn2*(gpn - rho(ji)*c*gunn) + grhos*bs2
            D(ji,1,4) = grho(ji,3)
            D(ji,1,5) = -two/c**3*dc*(L2 + pt5*(L5 + L1)) +             &
                        one/c**2*(-dc*(gpn-rho(ji)*c*gunn) +            &
                        (un-c)*(grhon*gmsinv - rho(ji)*dc*gunn) )

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
            
!.... compute the Characteristic amplitudes using the interior values
        
            L1 = ( un - c ) * ( gpn - rho(ji) * c * gunn )
            L2 = un * ( c**2 * grhon - gpn )
            L3 = un * gusn
            L4 = un * gu3n
            L5 = ( un + c ) * ( gpn + rho(ji) * c * gunn )

!.... Boundary conditions

            L2 = zero
            L3 = zero
            L4 = zero
            L5 = L1

!.... Continuity equation

            gmsinv = one / (gamma * Ma**2)
    
            D(ji,1,1) = one/c**2*(un-c)*(gtn*gmsinv - c*gunn) + guss + gu(ji,3,3)
            D(ji,1,2) = one/c**2*bn1*(gpn - rho(ji)*c*gunn) + grhos*bs1
            D(ji,1,3) = one/c**2*bn2*(gpn - rho(ji)*c*gunn) + grhos*bs2
            D(ji,1,4) = grho(ji,3)
            D(ji,1,5) = -two/c**3*dc*(L2 + pt5*(L5 + L1)) +             &
                        one/c**2*(-dc*(gpn-rho(ji)*c*gunn) +            &
                        (un-c)*(grhon*gmsinv - rho(ji)*dc*gunn) )
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
    
!.... get the metrics along this boundary

            if (i .eq. 1) then
!             bn1 = -m1(ji) / sqrt( m1(ji)**2 + m2(ji)**2 )
!             bn2 =  m2(ji) / sqrt( m1(ji)**2 + m2(ji)**2 )
              bn1 = one
              bn2 = zero
            else
!             bn1 =  m1(ji) / sqrt( m1(ji)**2 + m2(ji)**2 )
!             bn2 =  m2(ji) / sqrt( m1(ji)**2 + m2(ji)**2 )
              bn1 = one
              bn2 = zero
            end if

            bs1 = -bn2
            bs2 =  bn1

!.... compute Us and Un
        
            un = bn1 * vl(ji,2) + bn2 * vl(ji,3)
            us = bs1 * vl(ji,2) + bs2 * vl(ji,3)

            c = sqrt( t(ji) ) / Ma
            dc = pt5 / (Ma**2 * c)

!.... compute the Characteristic amplitudes
        
            L1 = ( u1(ji) - c ) * ( gp(ji,1) - rho(ji) * c * gu(ji,1,1) )
            L2 = u1(ji) * ( c**2 * grho(ji,1) - gp(ji,1) )
            L3 = u1(ji) * gu(ji,2,1)
            L4 = u1(ji) * gu(ji,3,1)
            L5 = ( u1(ji) + c ) * ( gp(ji,1) + rho(ji) * c * gu(ji,1,1) )
        
            L1 = rk * ( p(ji) - pinf )
        
!.... compute the streamwise derivative terms

            d1 = one / c**2 * ( L2 + pt5 * ( L5 + L1 ) )
            d2 = pt5 * ( L5 + L1 )
            d3 = pt5 / ( rho(ji) * c ) * ( L5 - L1 )
            d4 = L3
            d5 = L4

!=======================================================================================================!

            gmsinv = one / (gamma * Ma**2)

!.... Continuity equation

            D(ji,1,1) = gu(ji,2,2) + gu(ji,3,3) +                                       &
                        one/c**2 * ( -u1(ji)*gt(ji,1)*gmsinv +                          &
                        pt5*(u1(ji)+c)*(gmsinv*gt(ji,1)+c*gu(ji,1,1)) +                 &
                        pt5*rk*t(ji)*gmsinv )
            D(ji,1,2) = one/c**2*( c**2*grho(ji,1) - gp(ji,1) +                         &
                        pt5*(gp(ji,1)+rho(ji)*c*gu(ji,1,1) ) )
            D(ji,1,3) = grho(ji,2)
            D(ji,1,4) = grho(ji,3)
            D(ji,1,5) = -two * dc / c**3 * ( L2 + pt5 * ( L5 + L1 ) ) +                 &
                        one/c**2 * ( u1(ji)*(two*c*dc*grho(ji,1) -                      &
                        grho(ji,1)*gmsinv ) +                                           &
                        pt5*( dc*(gp(ji,1)+rho(ji)*c*gu(ji,1,1)) +                      &
                        (u1(ji)+c)*(gmsinv*grho(ji,1) + rho(ji)*dc*gu(ji,1,1))) +       &
                        pt5*rk*rho(ji)*gmsinv )

!=======================================================================================================!

!.... Momentum equation -- x_1

        fact1 = rhoinv(ji) / Re

        D(ji,2,1) = -one/(two*rho(ji)**2*c)*(L5 - L1) +                                 &
            one/(two*rho(ji)*c)*( (u1(ji)+c)*(gmsinv*gt(ji,1)+c*gu(ji,1,1)) -           &
            rk*t(ji)*gmsinv ) +                                                         &
            rhoinv(ji) * (d3 + u2(ji) * gu(ji,1,2) + u3(ji) * gu(ji,1,3))
        D(ji,2,2) = one/(two*rho(ji)*c) * (gp(ji,1) + rho(ji)*c*gu(ji,1,1))
        D(ji,2,3) = gu(ji,1,2)
        D(ji,2,4) = gu(ji,1,3)
        D(ji,2,5) = -dc/(two*rho(ji)*c**2)*( L5 - L1 ) +                                &
                    one/(two*rho(ji)*c)*( dc*(gp(ji,1)+rho(ji)*c*gu(ji,1,1)) +          &
                    (u1(ji)+c)*( gmsinv*grho(ji,1) + rho(ji)*gu(ji,1,1)*dc ) -          &
                    rk*rho(ji)*gmsinv ) -                                               &
                   fact1 * ( g1dlm(ji) * divu(ji) + dlm(ji) * g1divu(ji) ) -    &
                   fact1 * ( g1dmu(ji) * ( gu(ji,1,1) + gu(ji,1,1) ) +  &
                             g2dmu(ji) * ( gu(ji,1,2) + gu(ji,2,1) ) +  &
                             g3dmu(ji) * ( gu(ji,1,3) + gu(ji,3,1) ) +  &
                               dmu(ji) * ( two * g11v(ji,2) +           &
                                           g22v(ji,2) + g12v(ji,3) ) )

!=======================================================================================================!

!.... Momentum equation -- x_2          ( set tau21,1 = zero )

        D(ji,3,5) = rhoinv(ji) * gmsinv * grho(ji,2)  -                                 &
                   fact1 * ( g2dlm(ji) * (gu(ji,1,1) +                  &
                                      gu(ji,2,2) + gu(ji,3,3)) +                        &
                             dlm(ji) * (g12v(ji,2) + g22v(ji,3)) ) -            &
                   fact1 * ( g2dmu(ji) * ( gu(ji,2,2) + gu(ji,2,2) ) +                  &
                             g3dmu(ji) * ( gu(ji,2,3) + gu(ji,3,2) ) +                  &
                               dmu(ji) * ( g22v(ji,3) + g22v(ji,3) ) )

!=======================================================================================================!

!.... Energy equation                   ( set q1,1 = zero )

        fact1 = gamma * rhoinv(ji) / (Pr * Re)
        fact2 = gamma * gamma1 * Ma**2 * rhoinv(ji) / Re

        D(ji,5,1) = Ma**2 * rhoinv(ji)**2 * ( L2 - pt5*gamma1*( L5 + L1 ) ) -   &
                    Ma**2 * rhoinv(ji) * ( u1(ji)*(-gt(ji,1))*gmsinv -          &
                    pt5*gamma1*((u1(ji)+c)*(gt(ji,1)*gmsinv+                    &
                    c*gu(ji,1,1))+rk*t(ji)*gmsinv) ) +                          &
                    rhoinv(ji) * ( gamma*Ma**2*rhoinv(ji)*d2 -                  &
                    t(ji)*rhoinv(ji)*d1 + u2(ji) * gt(ji,2) +                   &
                    u3(ji) * gt(ji,3) + gamma1*t(ji)*(gu(ji,2,2)+gu(ji,3,3)) )

        D(ji,5,2) = -Ma**2*rhoinv(ji)*( c**2*grho(ji,1)-gp(ji,1) -              &
                     pt5*gamma1*(gp(ji,1)+rho(ji)*c*gu(ji,1,1)) )
        D(ji,5,3) = gt(ji,2)
        D(ji,5,4) = gt(ji,3)
        D(ji,5,5) = -Ma**2*rhoinv(ji)*( u1(ji)*(two*c*dc*grho(ji,1) -           &
                    grho(ji,1)*gmsinv) - pt5*gamma1*(dc*(gp(ji,1) +             &
                    rho(ji)*c*gu(ji,1,1)) + (u1(ji)+c)*(gmsinv*grho(ji,1) +     &
                    rho(ji)*dc*gu(ji,1,1)) + rk*rho(ji)*gmsinv) ) -             &
                    gamma1 * rhoinv(ji) * ( d1 + u2(ji) * grho(ji,2) +          &
                    u3(ji) * grho(ji,3) ) -                                     &
                    fact1 * (g2dcon(ji) * gt(ji,2) +                            &
                             g3dcon(ji) * gt(ji,3) +                            &
                             dcon(ji) * (g22v(ji,5))) -                         &
                    fact2 * dlm(ji) * ( gu(ji,1,1) +            &
                            gu(ji,2,2) + gu(ji,3,3) )**2 -      &
                    two * fact2 * dmu(ji) * (                   &
                    ( pt5 * ( gu(ji,1,1) + gu(ji,1,1) ) )**2 +  &
                    ( pt5 * ( gu(ji,1,2) + gu(ji,2,1) ) )**2 +  &
                    ( pt5 * ( gu(ji,1,3) + gu(ji,3,1) ) )**2 +  &
                    ( pt5 * ( gu(ji,2,1) + gu(ji,1,2) ) )**2 +  &
                    ( pt5 * ( gu(ji,2,2) + gu(ji,2,2) ) )**2 +  &
                    ( pt5 * ( gu(ji,2,3) + gu(ji,2,3) ) )**2 +  &
                    ( pt5 * ( gu(ji,3,1) + gu(ji,1,3) ) )**2 +  &
                    ( pt5 * ( gu(ji,3,2) + gu(ji,2,3) ) )**2 +  &
                    ( pt5 * ( gu(ji,3,3) + gu(ji,3,3) ) )**2 )

          end do
        end if
        end do
        end if
        
        return
        end
