!=======================================================================================================!
        subroutine genA_n(A, vl, i, j)
!
!  Form the A matrix which multiplies U,x
!
!  Revised: 6-6-22 
!
!  This version has been parabolized in x and implements the Lele & Poinsot
!  boundary conditions for the nonlinear equations.
!=======================================================================================================!
        use global
        use local2d
        use local3d
        implicit none
        
        real :: A(ndof,ndof), vl(ndof)
        real :: fact1, fact2, fact3, fact4, gmsinv

        real :: pinf, cinf, Lx, rk, c
        real :: L1, L2, L3, L4, L5
        real :: d1, d2, d3, d4, d5, dc
        
        integer :: i, j

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
!       S i d e   B o u n d a r i e s
!=======================================================================================================!
        if (left.eq.1 .or. right.eq.1) then

!.... do for both the left and right boundaries

        if ( (i.eq.1 .and. left.eq.1) .or. (i.eq.nx .and. right.eq.1) ) then    
    
!.... get the boundary normal and tangent vectors

            if (i .eq. 1) then
              bn1 = -m1(i,j) / sqrt( m1(i,j)**2 + m2(i,j)**2 )
              bn2 =  m2(i,j) / sqrt( m1(i,j)**2 + m2(i,j)**2 )
            else
              bn1 =  m1(i,j) / sqrt( m1(i,j)**2 + m2(i,j)**2 )
              bn2 =  m2(i,j) / sqrt( m1(i,j)**2 + m2(i,j)**2 )
            end if

            bs1 = -bn2
            bs2 =  bn1

!.... compute Us and Un
        
            un = bn1 * vl(2) + bn2 * vl(3)
            us = bs1 * vl(2) + bs2 * vl(3)

            c = sqrt( t ) / Ma
            dc = pt5 / (Ma**2 * c)

!.... compute some derivatives in the boundary normal coordinate system

            gpn = bn1 * gp(1) + bn2 * gp(2)
            gps = bs1 * gp(1) + bs2 * gp(2)
        
            grhon = bn1 * grho(1) + bn2 * grho(2)
            grhos = bs1 * grho(1) + bs2 * grho(2)

            gtn = bn1 * gt(1) + bn2 * gt(2)
            gts = bs1 * gt(1) + bs2 * gt(2)

            gu1n = bn1 * gu(1,1) + bn2 * gu(1,2)
            gu1s = bs1 * gu(1,1) + bs2 * gu(1,2)

            gu2n = bn1 * gu(2,1) + bn2 * gu(2,2)
            gu2s = bs1 * gu(2,1) + bs2 * gu(2,2)

            gu3n = bn1 * gu(3,1) + bn2 * gu(3,2)
            gu3s = bs1 * gu(3,1) + bs2 * gu(3,2)

            gunn = bn1*bn1 * gu(1,1) + bn2*bn1 * gu(2,1) + &
                   bn1*bn2 * gu(1,2) + bn2*bn2 * gu(2,2)

            gusn = bs1*bn1 * gu(1,1) + bs2*bn1 * gu(2,1) + &
                   bs1*bn2 * gu(1,2) + bs2*bn2 * gu(2,2)

            guns = bn1*bs1 * gu(1,1) + bn2*bs1 * gu(2,1) + &
                   bn1*bs2 * gu(1,2) + bn2*bs2 * gu(2,2)

            guss = bs1*bs1 * gu(1,1) + bs2*bs1 * gu(2,1) + &
                   bs1*bs2 * gu(1,2) + bs2*bs2 * gu(2,2)
   
        if (un .le. zero) then

!=======================================================================================================!
!       I n f l o w
!=======================================================================================================!

!.... Continuity equation

        gmsinv = one / (gamma * Ma**2)

        A(1,1) = pt5/c**2*(un+c)*bn1*t*gmsinv + us*bs1
        A(1,2) = pt5/c**2*(un+c)*rho*c*bn1*bn1 + rho*bs1*bs1
        A(1,3) = pt5/c**2*(un+c)*rho*c*bn2*bn1 + rho*bs2*bs1
        A(1,4) = zero
        A(1,5) = pt5/c**2*(un+c)*bn1*rho*gmsinv

!=======================================================================================================!

!.... Momentum equation -- x_1

        fact1  = rhoinv / Re

        A(2,1) = pt5*bn1*rhoinv/c*(un+c)*bn1*t*gmsinv +              &
                    bs1*rhoinv*gmsinv*t*bs1
        A(2,2) = pt5*bn1*rhoinv/c*(un+c)*rho*c*bn1*bn1 + us*bs1 -    &
                    damp(i,j) * fact1 * g1lm -                       &
                    damp(i,j) * fact1 * two * g1mu
        A(2,3) = pt5*bn1*rhoinv/c*(un+c)*rho*c*bn2*bn1 -             &
                    damp(i,j) * fact1 * g2mu
        A(2,4) = -damp(i,j) * fact1 * g3mu
        A(2,5) = bn1/(two*rho*c)*(un+c)*bn1*rho*gmsinv +             &
                    bs1*gmsinv*bs1 -                                 &
                    damp(i,j) * fact1 * dlm *                        &
                    ( gu(1,1) + gu(2,2) + gu(3,3) ) -                &
                    damp(i,j) * fact1 * dmu * ( gu(1,1) + gu(1,1) )

!=======================================================================================================!

!.... Momentum equation -- x_2

        A(3,1) = pt5*bn2*rhoinv/c*(un+c)*bn1*t*gmsinv +      &
                    bs2*rhoinv*gmsinv*t*bs1
        A(3,2) = pt5*bn2*rhoinv/c*(un+c)*rho*c*bn1*bn1 -     &
                    damp(i,j) * fact1 * g2lm
        A(3,3) = pt5*bn2*rhoinv/c*(un+c)*rho*c*bn2*bn1 +     &
                    us*bs1 - damp(i,j) * fact1 * g1mu
        A(3,4) = zero
        A(3,5) = bn2/(two*rho*c)*(un+c)*bn1*rho*gmsinv +     &
                    bs2*gmsinv*bs1 -                         &
                    damp(i,j)*fact1*dmu*( gu(2,1) + gu(1,2) )

!=======================================================================================================!

!.... Momentum equation -- x_3

        A(4,1) = zero
        A(4,2) = -damp(i,j) * fact1 * g3lm
        A(4,3) = zero
        A(4,4) = bs1*us - damp(i,j) * fact1 * g1mu
        A(4,5) = -damp(i,j) * fact1 * dmu * ( gu(3,1) + gu(1,3) )

!=======================================================================================================!

!.... Energy equation

        fact1 = gamma * rhoinv / (Pr * Re)
        fact2 = gamma * gamma1 * Ma**2 * rhoinv / Re

        A(5,1) =  Ma**2*pt5*rhoinv*gamma1*(un+c)*bn1*gmsinv*t
        A(5,2) =  Ma**2*pt5*rhoinv*gamma1*(un+c)*rho*c*bn1*bn1 +     &
                     gamma1*t*bs1*bs1 -                                &
                     damp(i,j) * fact2 * two * lm * divu -             &
                     damp(i,j) * two * fact2 * mu * ( gu(1,1) + gu(1,1))
        A(5,3) =  Ma**2*pt5*rhoinv*gamma1*(un+c)*rho*c*bn2*bn1 +       &
                     gamma1*t*bs2*bs1 -                                &
                     damp(i,j) * two * fact2 * mu * ( gu(2,1) + gu(1,2))
        A(5,4) = -damp(i,j) * two * fact2 * mu * ( gu(3,1) + gu(1,3) )
        A(5,5) =  Ma**2*pt5*rhoinv*gamma1*(un+c)*bn1*rho*gmsinv + us*bs1 -   &
                     damp(i,j) * fact1 * (g1con + dcon * gt(1))

        else

!=======================================================================================================!
!       O u t f l o w
!=======================================================================================================!

!.... Continuity equation

        gmsinv = one / (gamma * Ma**2)

        A(1,1) = one/c**2*( un*(c**2*bn1-bn1*t*gmsinv) +         &
                    pt5*(un+c)*bn1*gmsinv*t ) + us*bs1
        A(1,2) = one/c**2*pt5*(un+c)*rho*c*bn1*bn1 +             &
                    rho*bs1*bs1
        A(1,3) = one/c**2*pt5*(un+c)*rho*c*bn2*bn1 +             &
                    rho*bs2*bs1
        A(1,4) = zero
        A(1,5) = one/c**2*( un*(-bn1*rho*gmsinv) +               &
                    pt5*(un+c)*bn1*rho*gmsinv )

!=======================================================================================================!

!.... Momentum equation -- x_1

        fact1  = rhoinv / Re

        A(2,1) = pt5*bn1*rhoinv/c*(un+c)*bn1*t*gmsinv +              &
                    bs1*rhoinv*gmsinv*t*bs1
        A(2,2) = pt5*bn1*rhoinv/c*(un+c)*rho*c*bn1*bn1 +             &
                    bs1*un*bs1*bn1 + us*bs1 -                        &
                    damp(i,j) * fact1 * g1lm -                       &
                    damp(i,j) * fact1 * two * g1mu
        A(2,3) = pt5*bn1*rhoinv/c*(un+c)*rho*c*bn2*bn1 +             &
                    bs1*un*bs2*bn1 -                                 &
                    damp(i,j) * fact1 * g2mu
        A(2,4) = -damp(i,j) * fact1 * g3mu
        A(2,5) = bn1/(two*rho*c)*(un+c)*bn1*rho*gmsinv +             &
                    bs1*gmsinv*bs1 -                                 &
                    damp(i,j) * fact1 * dlm *                        &
                    ( gu(1,1) + gu(2,2) + gu(3,3) ) -                &
                    damp(i,j) * fact1 * dmu * ( gu(1,1) + gu(1,1) )

!=======================================================================================================!

!.... Momentum equation -- x_2

        A(3,1) = pt5*bn2*rhoinv/c*(un+c)*bn1*t*gmsinv +      &
                    bs2*rhoinv*gmsinv*t*bs1
        A(3,2) = pt5*bn2*rhoinv/c*(un+c)*rho*c*bn1*bn1 +     &
                    bs2*un*bs1*bn1 -                         &
                    damp(i,j) * fact1 * g2lm
        A(3,3) = pt5*bn2*rhoinv/c*(un+c)*rho*c*bn2*bn1 +     &
                    bs2*un*bs2*bn1 + us*bs1 - damp(i,j) * fact1 * g1mu
        A(3,4) = zero
        A(3,5) = bn2/(two*rho*c)*(un+c)*bn1*rho*gmsinv +     &
                    bs2*gmsinv*bs1 -                         &
                    damp(i,j)*fact1*dmu*( gu(2,1) + gu(1,2) )

!=======================================================================================================!

!.... Momentum equation -- x_3

        A(4,1) = zero
        A(4,2) = -damp(i,j) * fact1 * g3lm
        A(4,3) = zero
        A(4,4) = bn1*un + bs1*us - damp(i,j) * fact1 * g1mu
        A(4,5) = -damp(i,j) * fact1 * dmu * ( gu(3,1) + gu(1,3) )

!=======================================================================================================!

!.... Energy equation

        fact1 = gamma * rhoinv / (Pr * Re)
        fact2 = gamma * gamma1 * Ma**2 * rhoinv / Re

        A(5,1) = -Ma**2*rhoinv*( un*(c**2*bn1-bn1*t*gmsinv) -        &
                     pt5*gamma1*(un+c)*bn1*t*gmsinv )

        A(5,2) =  Ma**2*rhoinv*pt5*gamma1*(un+c)*rho*c*bn1*bn1 +     &
                     gamma1*t*bs1*bs1 -                                     &
                     damp(i,j) * fact2 * two * lm * divu -               &
                     damp(i,j) * two * fact2 * mu * ( gu(1,1) + gu(1,1) )

        A(5,3) =  Ma**2*rhoinv*pt5*gamma1*(un+c)*rho*c*bn2*bn1 +     &
                     gamma1*t*bs2*bs1 -                                     &
                     damp(i,j) * two * fact2 * mu * ( gu(2,1) + gu(1,2) )

        A(5,4) = -damp(i,j) * two * fact2 * mu * ( gu(3,1) + gu(1,3) )

        A(5,5) = -Ma**2*rhoinv*( un*(-bn1*rho*gmsinv) -              &
                     pt5*gamma1*(un+c)*bn1*rho*gmsinv ) + us*bs1 -          &
                     damp(i,j) * fact1 * (g1con + dcon * gt(1))

        end if          ! inflow/outflow

        end if
        end if
                
!=======================================================================================================!
!       T o p   B o u n d a r y
!=======================================================================================================!
        if (top.eq.1 .and. j.eq.ny .and. i.lt.nx) then

!.... get the boundary normal and tangent vectors

            bn1 = n1(i,j) / sqrt( n1(i,j)**2 + n2(i,j)**2 )
            bn2 = n2(i,j) / sqrt( n1(i,j)**2 + n2(i,j)**2 )

            bs1 = -bn2
            bs2 =  bn1

!.... compute Us and Un
        
            un = bn1 * vl(2) + bn2 * vl(3)
            us = bs1 * vl(2) + bs2 * vl(3)

            c = sqrt( t ) / Ma
            dc = pt5 / (Ma**2 * c)

!.... compute some derivatives in the boundary normal coordinate system

            gpn = bn1 * gp(1) + bn2 * gp(2)
            gps = bs1 * gp(1) + bs2 * gp(2)
        
            grhon = bn1 * grho(1) + bn2 * grho(2)
            grhos = bs1 * grho(1) + bs2 * grho(2)

            gtn = bn1 * gt(1) + bn2 * gt(2)
            gts = bs1 * gt(1) + bs2 * gt(2)

            gu1n = bn1 * gu(1,1) + bn2 * gu(1,2)
            gu1s = bs1 * gu(1,1) + bs2 * gu(1,2)

            gu2n = bn1 * gu(2,1) + bn2 * gu(2,2)
            gu2s = bs1 * gu(2,1) + bs2 * gu(2,2)

            gu3n = bn1 * gu(3,1) + bn2 * gu(3,2)
            gu3s = bs1 * gu(3,1) + bs2 * gu(3,2)

            gunn = bn1*bn1 * gu(1,1) + bn2*bn1 * gu(2,1) + &
                   bn1*bn2 * gu(1,2) + bn2*bn2 * gu(2,2)

            gusn = bs1*bn1 * gu(1,1) + bs2*bn1 * gu(2,1) + &
                   bs1*bn2 * gu(1,2) + bs2*bn2 * gu(2,2)

            guns = bn1*bs1 * gu(1,1) + bn2*bs1 * gu(2,1) + &
                   bn1*bs2 * gu(1,2) + bn2*bs2 * gu(2,2)

            guss = bs1*bs1 * gu(1,1) + bs2*bs1 * gu(2,1) + &
                   bs1*bs2 * gu(1,2) + bs2*bs2 * gu(2,2)
   
        if (un .le. zero) then

!=======================================================================================================!
!       I n f l o w
!=======================================================================================================!

!.... Continuity equation

        gmsinv = one / (gamma * Ma**2)

        A(1,1) = pt5/c**2*(un+c)*bn1*t*gmsinv + us*bs1
        A(1,2) = pt5/c**2*(un+c)*rho*c*bn1*bn1 + rho*bs1*bs1
        A(1,3) = pt5/c**2*(un+c)*rho*c*bn2*bn1 + rho*bs2*bs1
        A(1,4) = zero
        A(1,5) = pt5/c**2*(un+c)*bn1*rho*gmsinv

!.... set u,v,w,T

!       A(1,1) = gamma/c**2*(un+c)*bn1*t*gmsinv + us*bs1
!       A(1,2) = gamma/c**2*(un+c)*rho*c*bn1*bn1 + rho*bs1*bs1
!       A(1,3) = gamma/c**2*(un+c)*rho*c*bn2*bn1 + rho*bs2*bs1
!       A(1,4) = zero
!       A(1,5) = gamma/c**2*(un+c)*bn1*rho*gmsinv

!=======================================================================================================!

!.... Momentum equation -- x_1

        fact1  = rhoinv / Re

        A(2,1) = pt5*bn1*rhoinv/c*(un+c)*bn1*t*gmsinv +              &
                    bs1*rhoinv*gmsinv*t*bs1
        A(2,2) = pt5*bn1*rhoinv/c*(un+c)*rho*c*bn1*bn1 + us*bs1 -    &
                    damp(i,j) * fact1 * g1lm -                               &
                    damp(i,j) * fact1 * two * g1mu
        A(2,3) = pt5*bn1*rhoinv/c*(un+c)*rho*c*bn2*bn1 -             &
                    damp(i,j) * fact1 * g2mu
        A(2,4) = -damp(i,j) * fact1 * g3mu
        A(2,5) = bn1/(two*rho*c)*(un+c)*bn1*rho*gmsinv +             &
                    bs1*gmsinv*bs1 -                                            &
                    damp(i,j) * fact1 * dlm *                                &
                    ( gu(1,1) + gu(2,2) + gu(3,3) ) -                  &
                    damp(i,j) * fact1 * dmu * ( gu(1,1) + gu(1,1) )

!=======================================================================================================!

!.... Momentum equation -- x_2

        A(3,1) = pt5*bn2*rhoinv/c*(un+c)*bn1*t*gmsinv +      &
                    bs2*rhoinv*gmsinv*t*bs1
        A(3,2) = pt5*bn2*rhoinv/c*(un+c)*rho*c*bn1*bn1 -     &
                    damp(i,j) * fact1 * g2lm
        A(3,3) = pt5*bn2*rhoinv/c*(un+c)*rho*c*bn2*bn1 +     &
                    us*bs1 - damp(i,j) * fact1 * g1mu
        A(3,4) = zero
        A(3,5) = bn2/(two*rho*c)*(un+c)*bn1*rho*gmsinv +     &
                    bs2*gmsinv*bs1 -                                    &
                    damp(i,j)*fact1*dmu*( gu(2,1) + gu(1,2) )

!=======================================================================================================!

!.... Momentum equation -- x_3

        A(4,1) = zero
        A(4,2) = -damp(i,j) * fact1 * g3lm
        A(4,3) = zero
        A(4,4) = bs1*us - damp(i,j) * fact1 * g1mu
        A(4,5) = -damp(i,j) * fact1 * dmu * ( gu(3,1) + gu(1,3) )

!=======================================================================================================!

!.... Energy equation

        fact1 = gamma * rhoinv / (Pr * Re)
        fact2 = gamma * gamma1 * Ma**2 * rhoinv / Re

        A(5,1) =  Ma**2*pt5*rhoinv*gamma1*(un+c)*bn1*gmsinv*t
        A(5,2) =  Ma**2*pt5*rhoinv*gamma1*(un+c)*rho*c*bn1*bn1 +     &
                     gamma1*t*bs1*bs1 -                                     &
                     damp(i,j) * fact2 * two * lm * divu -               &
                     damp(i,j) * two * fact2 * mu * ( gu(1,1) + gu(1,1) )
        A(5,3) =  Ma**2*pt5*rhoinv*gamma1*(un+c)*rho*c*bn2*bn1 +     &
                     gamma1*t*bs2*bs1 -                                     &
                     damp(i,j) * two * fact2 * mu * ( gu(2,1) + gu(1,2) )
        A(5,4) = -damp(i,j) * two * fact2 * mu * ( gu(3,1) + gu(1,3) )
        A(5,5) =  Ma**2*pt5*rhoinv*gamma1*(un+c)*bn1*rho*gmsinv + us*bs1 -   &
                     damp(i,j) * fact1 * (g1con + dcon * gt(1))

        else

!=======================================================================================================!
!       O u t f l o w
!=======================================================================================================!

!.... Continuity equation

        gmsinv = one / (gamma * Ma**2)

        A(1,1) = one/c**2*( un*(c**2*bn1-bn1*t*gmsinv) +         &
                    pt5*(un+c)*bn1*gmsinv*t ) + us*bs1
        A(1,2) = one/c**2*pt5*(un+c)*rho*c*bn1*bn1 +             &
                    rho*bs1*bs1
        A(1,3) = one/c**2*pt5*(un+c)*rho*c*bn2*bn1 +             &
                    rho*bs2*bs1
        A(1,4) = zero
        A(1,5) = one/c**2*( un*(-bn1*rho*gmsinv) +               &
                    pt5*(un+c)*bn1*rho*gmsinv )

!=======================================================================================================!

!.... Momentum equation -- x_1

        fact1  = rhoinv / Re

        A(2,1) = pt5*bn1*rhoinv/c*(un+c)*bn1*t*gmsinv +              &
                    bs1*rhoinv*gmsinv*t*bs1
        A(2,2) = pt5*bn1*rhoinv/c*(un+c)*rho*c*bn1*bn1 +             &
                    bs1*un*bs1*bn1 + us*bs1 -                                   &
                    damp(i,j) * fact1 * g1lm -                               &
                    damp(i,j) * fact1 * two * g1mu
        A(2,3) = pt5*bn1*rhoinv/c*(un+c)*rho*c*bn2*bn1 +             &
                    bs1*un*bs2*bn1 -                                            &
                    damp(i,j) * fact1 * g2mu
        A(2,4) = -damp(i,j) * fact1 * g3mu
        A(2,5) = bn1/(two*rho*c)*(un+c)*bn1*rho*gmsinv +             &
                    bs1*gmsinv*bs1 -                                            &
                    damp(i,j) * fact1 * dlm *                                &
                    ( gu(1,1) + gu(2,2) + gu(3,3) ) -                  &
                    damp(i,j) * fact1 * dmu * ( gu(1,1) + gu(1,1) )

!=======================================================================================================!

!.... Momentum equation -- x_2

        A(3,1) = pt5*bn2*rhoinv/c*(un+c)*bn1*t*gmsinv +      &
                    bs2*rhoinv*gmsinv*t*bs1
        A(3,2) = pt5*bn2*rhoinv/c*(un+c)*rho*c*bn1*bn1 +     &
                    bs2*un*bs1*bn1 -                                    &
                    damp(i,j) * fact1 * g2lm
        A(3,3) = pt5*bn2*rhoinv/c*(un+c)*rho*c*bn2*bn1 +     &
                    bs2*un*bs2*bn1 + us*bs1 - damp(i,j) * fact1 * g1mu
        A(3,4) = zero
        A(3,5) = bn2/(two*rho*c)*(un+c)*bn1*rho*gmsinv +     &
                    bs2*gmsinv*bs1 -                                    &
                    damp(i,j)*fact1*dmu*( gu(2,1) + gu(1,2) )

!=======================================================================================================!

!.... Momentum equation -- x_3

        A(4,1) = zero
        A(4,2) = -damp(i,j) * fact1 * g3lm
        A(4,3) = zero
        A(4,4) = bn1*un + bs1*us - damp(i,j) * fact1 * g1mu
        A(4,5) = -damp(i,j) * fact1 * dmu * ( gu(3,1) + gu(1,3) )

!=======================================================================================================!

!.... Energy equation

        fact1 = gamma * rhoinv / (Pr * Re)
        fact2 = gamma * gamma1 * Ma**2 * rhoinv / Re

        A(5,1) = -Ma**2*rhoinv*( un*(c**2*bn1-bn1*t*gmsinv) -        &
                     pt5*gamma1*(un+c)*bn1*t*gmsinv )

        A(5,2) =  Ma**2*rhoinv*pt5*gamma1*(un+c)*rho*c*bn1*bn1 +     &
                     gamma1*t*bs1*bs1 -                                     &
                     damp(i,j) * fact2 * two * lm * divu -               &
                     damp(i,j) * two * fact2 * mu * ( gu(1,1) + gu(1,1) )

        A(5,3) =  Ma**2*rhoinv*pt5*gamma1*(un+c)*rho*c*bn2*bn1 +     &
                     gamma1*t*bs2*bs1 -                                     &
                     damp(i,j) * two * fact2 * mu * ( gu(2,1) + gu(1,2) )

        A(5,4) = -damp(i,j) * two * fact2 * mu * ( gu(3,1) + gu(1,3) )

        A(5,5) = -Ma**2*rhoinv*( un*(-bn1*rho*gmsinv) -              &
                     pt5*gamma1*(un+c)*bn1*rho*gmsinv ) + us*bs1 -          &
                     damp(i,j) * fact1 * (g1con + dcon * gt(1))

        end if          ! inflow/outflow

        end if          ! top 
        
!=======================================================================================================!
!       W a l l   B o u n d a r y
!=======================================================================================================!
        if (wall.eq.4 .and. j.eq.1) then

!.... get the boundary normal and tangent vectors

          bn1 = n1(i,j) / sqrt( n1(i,j)**2 + n2(i,j)**2 )
          bn2 = n2(i,j) / sqrt( n1(i,j)**2 + n2(i,j)**2 )

          bs1 = -bn2
          bs2 =  bn1

!.... compute Us and Un
        
          un = bn1 * vl(2) + bn2 * vl(3)
          us = bs1 * vl(2) + bs2 * vl(3)

          c = sqrt( t ) / Ma
          dc = pt5 / (Ma**2 * c)

!.... compute some derivatives in the boundary normal coordinate system

          gpn = bn1 * gp(1) + bn2 * gp(2)
          gps = bs1 * gp(1) + bs2 * gp(2)
        
          grhon = bn1 * grho(1) + bn2 * grho(2)
          grhos = bs1 * grho(1) + bs2 * grho(2)

          gtn = bn1 * gt(1) + bn2 * gt(2)
          gts = bs1 * gt(1) + bs2 * gt(2)

          gu1n = bn1 * gu(1,1) + bn2 * gu(1,2)
          gu1s = bs1 * gu(1,1) + bs2 * gu(1,2)

          gu2n = bn1 * gu(2,1) + bn2 * gu(2,2)
          gu2s = bs1 * gu(2,1) + bs2 * gu(2,2)

          gu3n = bn1 * gu(3,1) + bn2 * gu(3,2)
          gu3s = bs1 * gu(3,1) + bs2 * gu(3,2)

          gunn = bn1*bn1 * gu(1,1) + bn2*bn1 * gu(2,1) + &
                 bn1*bn2 * gu(1,2) + bn2*bn2 * gu(2,2)

          gusn = bs1*bn1 * gu(1,1) + bs2*bn1 * gu(2,1) + &
                 bs1*bn2 * gu(1,2) + bs2*bn2 * gu(2,2)

          guns = bn1*bs1 * gu(1,1) + bn2*bs1 * gu(2,1) + &
                 bn1*bs2 * gu(1,2) + bn2*bs2 * gu(2,2)

          guss = bs1*bs1 * gu(1,1) + bs2*bs1 * gu(2,1) + &
                 bs1*bs2 * gu(1,2) + bs2*bs2 * gu(2,2)
   
!.... Continuity equation  (only need continuity at the wall)

          gmsinv = one / (gamma * Ma**2)
  
          A(1,1) = one/c**2*(un-c)*bn1*t*gmsinv + us*bs1
          A(1,2) = -one/c**2*(un-c)*rho*c*bn1*bn1 + rho*bs1*bs1
          A(1,3) = -one/c**2*(un-c)*rho*c*bn2*bn1 + rho*bs2*bs1
          A(1,4) = zero
          A(1,5) = one/c**2*(un-c)*bn1*rho*gmsinv

        end if          ! wall 

        return
        end
