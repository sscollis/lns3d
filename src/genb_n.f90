!=======================================================================================================!
        subroutine genB_n(B, vl, i, j)
!
!  Form the B matrix which multiplies U,y
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
        
        real :: B(ndof,ndof), vl(ndof)
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

       call error("genB_n$","Lele Poinsot BC's not updated to ij ordering$")

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

        gmsinv = one / ( gamma * Ma**2 )

        B(1,1) = pt5/c**2*(un+c)*bn2*t*gmsinv + us*bs2
        B(1,2) = pt5/c**2*(un+c)*rho*c*bn1*bn2 + rho*bs1*bs2
        B(1,3) = pt5/c**2*(un+c)*rho*c*bn2*bn2 + rho*bs2*bs2
        B(1,4) = zero
        B(1,5) = pt5/c**2*(un+c)*bn2*rho*gmsinv

!=======================================================================================================!

!.... Momentum equation -- x_1

        fact1  = rhoinv / Re

        B(2,1) = pt5*bn1*rhoinv/c*(un+c)*bn2*t*gmsinv +              &
                    bs1*rhoinv*gmsinv*t*bs2
        B(2,2) = pt5*bn1*rhoinv/c*(un+c)*rho*c*bn1*bn2 + us*bs2 -    &
                    fact1 * g2mu
        B(2,3) = pt5*bn1*rhoinv/c*(un+c)*rho*c*bn2*bn2 -             &
                    damp(i,j) * fact1 * g1lm
        B(2,4) = zero
        B(2,5) = bn1/(two*rho*c)*(un+c)*bn2*rho*gmsinv +             &
                    bs1*gmsinv*bs2 -                                 &
                    fact1 * dmu * ( gu(1,2) + damp(i,j) * gu(2,1) )

!=======================================================================================================!

!.... Momentum equation -- x_2

        B(3,1) = pt5*bn2*rhoinv/c*(un+c)*bn2*t*gmsinv +      &
                    bs2 * rhoinv * gmsinv * t * bs2
        B(3,2) = pt5*bn2*rhoinv/c*(un+c)*rho*c*bn1*bn2 -     &
                    damp(i,j) * fact1 * g1mu
        B(3,3) = pt5*bn2*rhoinv/c*(un+c)*rho*c*bn2*bn2 +     &
                    us*bs2 - fact1 * g2lm - fact1 * two * g2mu
        B(3,4) = -fact1 * g3mu
        B(3,5) = bn2/(two*rho*c)*(un+c)*bn2*rho*gmsinv +     &
                    bs2*gmsinv*bs2 - fact1 * dlm *           &
                    (damp(i,j) * gu(1,1) + gu(2,2) + gu(3,3)) -   &
                    fact1 * dmu * ( gu(2,2) + gu(2,2) )

!=======================================================================================================!

!.... Momentum equation -- x_3

        B(4,1) = zero
        B(4,2) = zero
        B(4,3) = -fact1 * g3lm
        B(4,4) = bs2*us - fact1 * g2mu
        B(4,5) = -fact1 * dmu * ( gu(3,2) + gu(2,3) )

!=======================================================================================================!

!.... Energy equation

        fact1 = gamma * rhoinv / (Pr * Re)
        fact2 = gamma * gamma1 * Ma**2 * rhoinv / Re

        B(5,1) =  Ma**2*pt5*rhoinv*gamma1*(un+c)*bn2*gmsinv*t

        B(5,2) =  Ma**2*pt5*rhoinv*gamma1*(un+c)*rho*c*bn1*bn2 +     &
                     gamma1*t*bs1*bs2 -                              &
                     two*fact2*mu*( gu(1,2) + damp(i,j) * gu(2,1) )

        B(5,3) =  Ma**2*pt5*rhoinv*gamma1*(un+c)*rho*c*bn2*bn2 +     &
                     gamma1*t*bs2*bs2 - fact2 * two * lm *           &
                     (damp(i,j) * gu(1,1) + gu(2,2) + gu(3,3)) -          &
                     two * fact2 * mu * ( gu(2,2) + gu(2,2) )

        B(5,4) = -two * fact2 * mu * ( gu(3,2) + gu(2,3) )

        B(5,5) =  Ma**2*pt5*rhoinv*gamma1*(un+c)*bn2*rho*gmsinv +    &
                     us*bs2 - fact1 * (g2con + dcon * gt(2))


        else

!=======================================================================================================!
!       O u t f l o w
!=======================================================================================================!

!.... Continuity equation

        gmsinv = one / ( gamma * Ma**2 )

        B(1,1) = one/c**2*( un*(c**2*bn2-bn2*t*gmsinv) +         &
                    pt5*(un+c)*bn2*gmsinv*t ) + us*bs2
        B(1,2) = one/c**2*pt5*(un+c)*rho*c*bn1*bn2 + rho*bs1*bs2
        B(1,3) = one/c**2*pt5*(un+c)*rho*c*bn2*bn2 + rho*bs2*bs2
        B(1,4) = zero
        B(1,5) = one/c**2*( un*(-bn2*rho*gmsinv) +               &
                    pt5*(un+c)*bn2*rho*gmsinv )

!=======================================================================================================!

!.... Momentum equation -- x_1

        fact1  = rhoinv / Re

        B(2,1) = pt5*bn1*rhoinv/c*(un+c)*bn2*t*gmsinv +              &
                    bs1*rhoinv*gmsinv*t*bs2
        B(2,2) = pt5*bn1*rhoinv/c*(un+c)*rho*c*bn1*bn2 +             &
                    bs1*un*bs1*bn2 + us*bs2 -                        &
                    fact1 * g2mu
        B(2,3) = pt5*bn1*rhoinv/c*(un+c)*rho*c*bn2*bn2 +             &
                    bs1*un*bs2*bn2 -                                 &
                    damp(i,j) * fact1 * g1lm
        B(2,4) = zero
        B(2,5) = bn1/(two*rho*c)*(un+c)*bn2*rho*gmsinv +             &
                    bs1*gmsinv*bs2 -                                 &
                    fact1 * dmu * ( gu(1,2) + damp(i,j) * gu(2,1) )

!=======================================================================================================!

!.... Momentum equation -- x_2

        B(3,1) = pt5*bn2*rhoinv/c*(un+c)*bn2*t*gmsinv +      &
                    bs2 * rhoinv * gmsinv * t * bs2
        B(3,2) = pt5*bn2*rhoinv/c*(un+c)*rho*c*bn1*bn2 +     &
                    bs2*un*bs1*bn2 -                         &
                    damp(i,j) * fact1 * g1mu
        B(3,3) = pt5*bn2*rhoinv/c*(un+c)*rho*c*bn2*bn2 +     &
                    us*bs2 + bs2*un*bs2*bn2 -                &
                    fact1 * g2lm - fact1 * two * g2mu
        B(3,4) = -fact1 * g3mu
        B(3,5) = bn2/(two*rho*c)*(un+c)*bn2*rho*gmsinv +     &
                    bs2*gmsinv*bs2 - fact1 * dlm *           &
                    (damp(i,j) * gu(1,1) + gu(2,2) + gu(3,3)) -   &
                    fact1 * dmu * ( gu(2,2) + gu(2,2) )

!=======================================================================================================!

!.... Momentum equation -- x_3

        B(4,1) = zero
        B(4,2) = zero
        B(4,3) = -fact1 * g3lm
        B(4,4) = bn2*un + bs2*us - fact1 * g2mu
        B(4,5) = -fact1 * dmu * ( gu(3,2) + gu(2,3) )

!=======================================================================================================!

!.... Energy equation

        fact1 = gamma * rhoinv / (Pr * Re)
        fact2 = gamma * gamma1 * Ma**2 * rhoinv / Re

        B(5,1) = -Ma**2*rhoinv*( un*(c**2*bn2-bn2*t*gmsinv) -        &
                     pt5*gamma1*(un+c)*bn2*t*gmsinv ) 

        B(5,2) =  Ma**2*rhoinv*pt5*gamma1*(un+c)*rho*c*bn1*bn2 +     &
                     gamma1*t*bs1*bs2 -                              &
                     two*fact2*mu*( gu(1,2) + damp(i,j) * gu(2,1) )

        B(5,3) =  Ma**2*rhoinv*pt5*gamma1*(un+c)*rho*c*bn2*bn2 +     &
                     gamma1*t*bs2*bs2 - fact2 * two * lm *           &
                     (damp(i,j) * gu(1,1) + gu(2,2) + gu(3,3)) -          &
                     two * fact2 * mu * ( gu(2,2) + gu(2,2) )

        B(5,4) = -two * fact2 * mu * ( gu(3,2) + gu(2,3) )

        B(5,5) = -Ma**2*rhoinv*( un*(-bn2*rho*gmsinv) -              &
                     pt5*gamma1*(un+c)*bn2*rho*gmsinv ) +            &
                     us*bs2 - fact1 * (g2con + dcon * gt(2))

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

        gmsinv = one / ( gamma * Ma**2 )

        B(1,1) = pt5/c**2*(un+c)*bn2*t*gmsinv + us*bs2
        B(1,2) = pt5/c**2*(un+c)*rho*c*bn1*bn2 + rho*bs1*bs2
        B(1,3) = pt5/c**2*(un+c)*rho*c*bn2*bn2 + rho*bs2*bs2
        B(1,4) = zero
        B(1,5) = pt5/c**2*(un+c)*bn2*rho*gmsinv

!.... set u,v,w,T

!       B(1,1) = gamma/c**2*(un+c)*bn2*t*gmsinv + us*bs2
!       B(1,2) = gamma/c**2*(un+c)*rho*c*bn1*bn2 + rho*bs1*bs2
!       B(1,3) = gamma/c**2*(un+c)*rho*c*bn2*bn2 + rho*bs2*bs2
!       B(1,4) = zero
!       B(1,5) = gamma/c**2*(un+c)*bn2*rho*gmsinv

!=======================================================================================================!

!.... Momentum equation -- x_1

        fact1  = rhoinv / Re

        B(2,1) = pt5*bn1*rhoinv/c*(un+c)*bn2*t*gmsinv +              &
                    bs1*rhoinv*gmsinv*t*bs2
        B(2,2) = pt5*bn1*rhoinv/c*(un+c)*rho*c*bn1*bn2 + us*bs2 -    &
                    fact1 * g2mu
        B(2,3) = pt5*bn1*rhoinv/c*(un+c)*rho*c*bn2*bn2 -             &
                    damp(i,j) * fact1 * g1lm
        B(2,4) = zero
        B(2,5) = bn1/(two*rho*c)*(un+c)*bn2*rho*gmsinv +             &
                    bs1*gmsinv*bs2 -                                 &
                    fact1 * dmu * ( gu(1,2) + damp(i,j) * gu(2,1) )

!=======================================================================================================!

!.... Momentum equation -- x_2

        B(3,1) = pt5*bn2*rhoinv/c*(un+c)*bn2*t*gmsinv +      &
                    bs2 * rhoinv * gmsinv * t * bs2
        B(3,2) = pt5*bn2*rhoinv/c*(un+c)*rho*c*bn1*bn2 -     &
                    damp(i,j) * fact1 * g1mu
        B(3,3) = pt5*bn2*rhoinv/c*(un+c)*rho*c*bn2*bn2 +     &
                    us*bs2 - fact1 * g2lm - fact1 * two * g2mu
        B(3,4) = -fact1 * g3mu
        B(3,5) = bn2/(two*rho*c)*(un+c)*bn2*rho*gmsinv +     &
                    bs2*gmsinv*bs2 - fact1 * dlm *           &
                    (damp(i,j) * gu(1,1) + gu(2,2) + gu(3,3)) -   &
                    fact1 * dmu * ( gu(2,2) + gu(2,2) )

!=======================================================================================================!

!.... Momentum equation -- x_3

        B(4,1) = zero
        B(4,2) = zero
        B(4,3) = -fact1 * g3lm
        B(4,4) = bs2*us - fact1 * g2mu
        B(4,5) = -fact1 * dmu * ( gu(3,2) + gu(2,3) )

!=======================================================================================================!

!.... Energy equation

        fact1 = gamma * rhoinv / (Pr * Re)
        fact2 = gamma * gamma1 * Ma**2 * rhoinv / Re

        B(5,1) =  Ma**2*pt5*rhoinv*gamma1*(un+c)*bn2*gmsinv*t

        B(5,2) =  Ma**2*pt5*rhoinv*gamma1*(un+c)*rho*c*bn1*bn2 +     &
                     gamma1*t*bs1*bs2 -                              &
                     two*fact2*mu*( gu(1,2) + damp(i,j) * gu(2,1) )

        B(5,3) =  Ma**2*pt5*rhoinv*gamma1*(un+c)*rho*c*bn2*bn2 +     &
                     gamma1*t*bs2*bs2 - fact2 * two * lm *           &
                     (damp(i,j) * gu(1,1) + gu(2,2) + gu(3,3)) -          &
                     two * fact2 * mu * ( gu(2,2) + gu(2,2) )

        B(5,4) = -two * fact2 * mu * ( gu(3,2) + gu(2,3) )

        B(5,5) =  Ma**2*pt5*rhoinv*gamma1*(un+c)*bn2*rho*gmsinv +    &
                     us*bs2 - fact1 * (g2con + dcon * gt(2))


        else

!=======================================================================================================!
!       O u t f l o w
!=======================================================================================================!

!.... Continuity equation

        gmsinv = one / ( gamma * Ma**2 )

        B(1,1) = one/c**2*( un*(c**2*bn2-bn2*t*gmsinv) +         &
                    pt5*(un+c)*bn2*gmsinv*t ) + us*bs2
        B(1,2) = one/c**2*pt5*(un+c)*rho*c*bn1*bn2 + rho*bs1*bs2
        B(1,3) = one/c**2*pt5*(un+c)*rho*c*bn2*bn2 + rho*bs2*bs2
        B(1,4) = zero
        B(1,5) = one/c**2*( un*(-bn2*rho*gmsinv) +               &
                    pt5*(un+c)*bn2*rho*gmsinv )

!=======================================================================================================!

!.... Momentum equation -- x_1

        fact1  = rhoinv / Re

        B(2,1) = pt5*bn1*rhoinv/c*(un+c)*bn2*t*gmsinv +              &
                    bs1*rhoinv*gmsinv*t*bs2
        B(2,2) = pt5*bn1*rhoinv/c*(un+c)*rho*c*bn1*bn2 +             &
                    bs1*un*bs1*bn2 + us*bs2 -                        &
                    fact1 * g2mu
        B(2,3) = pt5*bn1*rhoinv/c*(un+c)*rho*c*bn2*bn2 +             &
                    bs1*un*bs2*bn2 -                                 &
                    damp(i,j) * fact1 * g1lm
        B(2,4) = zero
        B(2,5) = bn1/(two*rho*c)*(un+c)*bn2*rho*gmsinv +             &
                    bs1*gmsinv*bs2 -                                 &
                    fact1 * dmu * ( gu(1,2) + damp(i,j) * gu(2,1) )

!=======================================================================================================!

!.... Momentum equation -- x_2

        B(3,1) = pt5*bn2*rhoinv/c*(un+c)*bn2*t*gmsinv +      &
                    bs2 * rhoinv * gmsinv * t * bs2
        B(3,2) = pt5*bn2*rhoinv/c*(un+c)*rho*c*bn1*bn2 +     &
                    bs2*un*bs1*bn2 -                         & 
                    damp(i,j) * fact1 * g1mu
        B(3,3) = pt5*bn2*rhoinv/c*(un+c)*rho*c*bn2*bn2 +     &
                    us*bs2 + bs2*un*bs2*bn2 -                &
                    fact1 * g2lm - fact1 * two * g2mu
        B(3,4) = -fact1 * g3mu
        B(3,5) = bn2/(two*rho*c)*(un+c)*bn2*rho*gmsinv +     &
                    bs2*gmsinv*bs2 - fact1 * dlm *           &
                    (damp(i,j) * gu(1,1) + gu(2,2) + gu(3,3)) -   &
                    fact1 * dmu * ( gu(2,2) + gu(2,2) )

!=======================================================================================================!

!.... Momentum equation -- x_3

        B(4,1) = zero
        B(4,2) = zero
        B(4,3) = -fact1 * g3lm
        B(4,4) = bn2*un + bs2*us - fact1 * g2mu
        B(4,5) = -fact1 * dmu * ( gu(3,2) + gu(2,3) )

!=======================================================================================================!

!.... Energy equation

        fact1 = gamma * rhoinv / (Pr * Re)
        fact2 = gamma * gamma1 * Ma**2 * rhoinv / Re

        B(5,1) = -Ma**2*rhoinv*( un*(c**2*bn2-bn2*t*gmsinv) -        &
                     pt5*gamma1*(un+c)*bn2*t*gmsinv ) 

        B(5,2) =  Ma**2*rhoinv*pt5*gamma1*(un+c)*rho*c*bn1*bn2 +     &
                     gamma1*t*bs1*bs2 -                              &
                     two*fact2*mu*( gu(1,2) + damp(i,j) * gu(2,1) )

        B(5,3) =  Ma**2*rhoinv*pt5*gamma1*(un+c)*rho*c*bn2*bn2 +     &
                     gamma1*t*bs2*bs2 - fact2 * two * lm *           &
                     (damp(i,j) * gu(1,1) + gu(2,2) + gu(3,3)) -     &
                     two * fact2 * mu * ( gu(2,2) + gu(2,2) )

        B(5,4) = -two * fact2 * mu * ( gu(3,2) + gu(2,3) )

        B(5,5) = -Ma**2*rhoinv*( un*(-bn2*rho*gmsinv) -              &
                     pt5*gamma1*(un+c)*bn2*rho*gmsinv ) +            &
                     us*bs2 - fact1 * (g2con + dcon * gt(2))

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
            
!.... Continuity equation

            gmsinv = one / ( gamma * Ma**2 )
    
            B(1,1) = one/c**2*(un-c)*bn2*t*gmsinv + us*bs2
            B(1,2) = -one/c**2*(un-c)*rho*c*bn1*bn2 + rho*bs1*bs2
            B(1,3) = -one/c**2*(un-c)*rho*c*bn2*bn2 + rho*bs2*bs2
            B(1,4) = zero
            B(1,5) = one/c**2*(un-c)*bn2*rho*gmsinv

        end if          ! wall

        return
        end
