!=======================================================================================================!
        subroutine genD_n(D, vl)
!
!  Form the D matrix which multiplies U
!
!  Revised: 9-12-95
!
!  This version has been parabolized in x and implements the Lele & Poinsot
!  boundary conditions for the nonlinear equations.
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

       call error("genD_n$","Lele Poinsot BC's not updated to ij ordering$")

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

        do i = 1, nx, nx-1
          if ( (i.eq.1 .and. left.eq.1) .or. (i.eq.nx .and. right.eq.1) ) then    
          do j = 1, ny
            ji = j + (i-1) * ny
    
!.... get the boundary normal and tangent vectors

            if (i .eq. 1) then
              bn1 = -m1(ji) / sqrt( m1(ji)**2 + m2(ji)**2 )
              bn2 =  m2(ji) / sqrt( m1(ji)**2 + m2(ji)**2 )
            else
              bn1 =  m1(ji) / sqrt( m1(ji)**2 + m2(ji)**2 )
              bn2 =  m2(ji) / sqrt( m1(ji)**2 + m2(ji)**2 )
            end if

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
        
!.... compute the potential solution

            b = one / pi
            if (Ma.lt.one) then
              rbeta = sqrt(one - Ma**2)
            else
              rbeta = sqrt(Ma**2 - one)
            end if
            ub1 = one + (x(j,1)-b)/pi/rbeta**2/((x(j,i)-b)**2 + &
                         rbeta**2*y(j,i)**2)
            ub2 = rbeta*y(j,i)/pi/((x(j,i)-b)**2+rbeta**2*y(j,i)**2)

            gub11 = one/(pi*rbeta**2*((x(j,i)-b)**2+rbeta**2*y(j,i)**2)) - &
                    two*(x(j,i)-b)**2/(pi*rbeta**2*((x(j,i)-b)**2+rbeta**2*y(j,i)**2)**2)

            gub12 = -two*(x(j,i)-b)*rbeta**2*y(j,i)/(pi*rbeta**2*((x(j,i)-b)**2+ &
                     rbeta**2*y(j,i)**2)**2)

            gub21 = -rbeta*y(j,i)*two*(x(j,i)-b)/(pi*((x(j,i)-b)**2+ &
                     rbeta**2*y(j,i)**2)**2)

            gub22 = rbeta/(pi*((x(j,i)-b)**2+rbeta**2*y(j,i)**2)) - &
                    two*rbeta**3*y(j,i)**2/(pi*((x(j,i)-b)**2+ &
                    rbeta**2*y(j,i)**2)**2)

            ubn = bn1*ub1 + bn2*ub2
            ubs = bs1*ub1 + bs2*ub2

            gubnn = bn1*bn1 * gub11 + bn2*bn1 * gub21 + &
                    bn1*bn2 * gub12 + bn2*bn2 * gub22

            gubsn = bs1*bn1 * gub11 + bs2*bn1 * gub21 + &
                    bs1*bn2 * gub12 + bs2*bn2 * gub22

            tb   = one - pt5*gamma1*Ma**2*( ub1**2 + ub2**2 - one )
            rhob = tb**(one/gamma1)
            pb   = rhob * tb / ( gamma * Ma**2 )
            cb   = sqrt(tb)/Ma

            gtbn = bn1 * (-gamma1*Ma**2*( ub1*gub11 + ub2*gub21 ) ) + &
                   bn2 * (-gamma1*Ma**2*( ub1*gub12 + ub2*gub22 ) )
            
            grhobn = one/gamma1 * tb**((two-gamma)/gamma1) * gtbn
            gpbn   = one/(gamma*Ma**2) * ( tb*grhobn + rhob*gtbn )

!.... nonreflecting boundary conditions
        
            if (un .le. zero ) then
              L1 = (ubn - cb)*(gpbn - rhob*cb*gubnn)
              L2 = ubn*(cb**2*grhobn - gpbn)
              L3 = ubn*gubsn
              L4 = zero
            else
              L1 = rk * ( p(ji) - pinf ) + (ubn - cb) * (gpbn - rhob*cb*gubnn)  
            end if

        if (un .le. zero) then

!=======================================================================================================!
!       I n f l o w
!=======================================================================================================!

!.... Continuity equation

        gmsinv = one / (gamma * Ma**2)

        D(ji,1,1) = pt5/c**2*(un+c)*(gtn*gmsinv + c*gunn) + guss + gu(ji,3,3)
        D(ji,1,2) = pt5/c**2*bn1*(gpn + rho(ji)*c*gunn) + grhos*bs1
        D(ji,1,3) = pt5/c**2*bn2*(gpn + rho(ji)*c*gunn) + grhos*bs2
        D(ji,1,4) = grho(ji,3)
        D(ji,1,5) = -two/c**3*dc*(L2 + pt5*(L5 + L1)) +         &
                    pt5/c**2*( dc*(gpn+rho(ji)*c*gunn) +        &
                    (un+c)*(grhon*gmsinv + rho(ji)*dc*gunn) )

!=======================================================================================================!

!.... Momentum equation -- x_1

        fact1  = rhoinv(ji) / Re

        D(ji,2,1) =                                                                     &
           -bn1/(two*rho(ji)**2*c)*( L5 - L1 ) +                                        &
            bn1/(two*rho(ji)*c)*(un+c)*( gtn*gmsinv + c*gunn ) -                        &
            bs1*rhoinv(ji)**2*gps + bs1*rhoinv(ji)*gmsinv*gts +                         &
            rhoinv(ji) * (                                                              &
            fact1 * damp(ji) * ( g1lm(ji) * divu(ji) + lm(ji) * g1divu(ji) ) +          &
            fact1 * ( damp(ji) * g1mu(ji) * ( two * gu(ji,1,1) ) +                      &
                            g2mu(ji) * ( gu(ji,1,2) + damp(ji) * gu(ji,2,1) ) +         &
                            g3mu(ji) * ( gu(ji,1,3) + damp(ji) * gu(ji,3,1) ) +         &
                              mu(ji) * ( damp(ji) * two * g11v(ji,2) +                  &
                                      g22v(ji,2) + damp(ji) * g12v(ji,3) ) ) )

        D(ji,2,2) = bn1/(two*rho(ji)*c)*bn1*(gpn + rho(ji)*c*gunn) + bs1 * gu1s
        D(ji,2,3) = bn1/(two*rho(ji)*c)*bn2*(gpn + rho(ji)*c*gunn) + bs2 * gu1s
        D(ji,2,4) = gu(ji,1,3)
        D(ji,2,5) =                                                                     & 
                  -bn1*dc/(two*rho(ji)*c**2)*( L5 - L1 ) +                              &
                   bn1/(two*rho(ji)*c)*( dc*(gpn + rho(ji)*c*gunn) +                    &
                   (un+c)*(grhon*gmsinv + rho(ji)*dc*gunn) ) +                          &
                   bs1*rhoinv(ji)*gmsinv*grhos -                                        &
                   fact1 * damp(ji) * ( g1dlm(ji) * divu(ji) + dlm(ji) * g1divu(ji) ) - &
                   fact1 * ( damp(ji) * g1dmu(ji) * ( gu(ji,1,1) + gu(ji,1,1) ) +       &
                             g2dmu(ji) * ( gu(ji,1,2) + damp(ji) * gu(ji,2,1) ) +       &
                             g3dmu(ji) * ( gu(ji,1,3) + damp(ji) * gu(ji,3,1) ) +       &
                               dmu(ji) * ( damp(ji) * two * g11v(ji,2) +                &
                                       g22v(ji,2) + damp(ji) * g12v(ji,3) ) )

!=======================================================================================================!

!.... Momentum equation -- x_2

        D(ji,3,1) =                                                                     &
           -bn2/(two*rho(ji)**2*c)*( L5 - L1 ) +                                        &
            bn2/(two*rho(ji)*c)*(un+c)*(gtn*gmsinv + c*gunn) -                          &
            bs2*rhoinv(ji)**2*gps + bs2*rhoinv(ji)*gmsinv*gts +                         &
            rhoinv(ji) * (                                                              &
            fact1 * ( g2lm(ji) * (damp(ji) * gu(ji,1,1) + gu(ji,2,2) +                  &
                gu(ji,3,3)) + lm(ji) * (damp(ji) * g12v(ji,2) + g22v(ji,3)) ) +         &
            fact1 * ( damp(ji) * g1mu(ji) * (gu(ji,2,1) + gu(ji,1,2) ) +                &
                      g2mu(ji) * ( gu(ji,2,2) + gu(ji,2,2) ) +                          &
                      g3mu(ji) * ( gu(ji,2,3) + gu(ji,3,2) ) +                          &
                        mu(ji) * ( damp(ji) * g11v(ji,3) + damp(ji) * g12v(ji,2) +      &
                                g22v(ji,3) + g22v(ji,3) ) ) )

        D(ji,3,2) = bn2/(two*rho(ji)*c)*bn1*(gpn + rho(ji)*c*gunn) + bs1 * gu2s 
        D(ji,3,3) = bn2/(two*rho(ji)*c)*bn2*(gpn + rho(ji)*c*gunn) + bs2 * gu2s
        D(ji,3,4) = gu(ji,2,3)
        D(ji,3,5) =                                                                     & 
                  -bn2*dc/(two*rho(ji)*c**2)*( L5 - L1 ) +                              &
                   bn2/(two*rho(ji)*c)*( dc*(gpn + rho(ji)*c*gunn) +                    &
                   (un+c)*(grhon*gmsinv + rho(ji)*dc*gunn) ) +                          &
                   bs2*rhoinv(ji)*gmsinv*grhos -                                        &
                   fact1 * ( g2dlm(ji) * (damp(ji) * gu(ji,1,1) +                       &
                                      gu(ji,2,2) + gu(ji,3,3)) +                        &
                             dlm(ji) * (damp(ji) * g12v(ji,2) + g22v(ji,3)) ) -         &
                   fact1 * ( damp(ji) * g1dmu(ji) * ( gu(ji,2,1) + gu(ji,1,2) ) +       &
                             g2dmu(ji) * ( gu(ji,2,2) + gu(ji,2,2) ) +                  &
                             g3dmu(ji) * ( gu(ji,2,3) + gu(ji,3,2) ) +                  &
                               dmu(ji) * ( damp(ji) * g11v(ji,3) +                      &
                                       damp(ji) * g12v(ji,2) +                          &
                                       g22v(ji,3) + g22v(ji,3) ) )

!=======================================================================================================!

!.... Momentum equation -- x_3

        D(ji,4,1) = rhoinv(ji) * ( gmsinv * gt(ji,3) - rhoinv(ji) * gp(ji,3) +  &
            fact1 * ( g3lm(ji) * (damp(ji) * gu(ji,1,1) + gu(ji,2,2) +          &
                      gu(ji,3,3)) + lm(ji) * g3divu(ji) ) +                     &
            fact1 * ( damp(ji) * g1mu(ji) * (gu(ji,3,1) + gu(ji,1,3) ) +        &
                      g2mu(ji) * ( gu(ji,3,2) + gu(ji,2,3) ) +                  &
                      g3mu(ji) * ( gu(ji,3,3) + gu(ji,3,3) ) +                  &
                        mu(ji) * ( damp(ji) * g11v(ji,4) + g22v(ji,4) ) ) )

        D(ji,4,2) = bs1 * gu3s
        D(ji,4,3) = bs2 * gu3s
        D(ji,4,4) = gu(ji,3,3)
        D(ji,4,5) = rhoinv(ji) * gmsinv * grho(ji,3)  -                         &
                   fact1 * ( g3dlm(ji) * (damp(ji) * gu(ji,1,1) +               &
                                      gu(ji,2,2) + gu(ji,3,3)) +                &
                               dlm(ji) * g3divu(ji) ) -                         &
                   fact1 * ( damp(ji) * g1dmu(ji) * (gu(ji,3,1) + gu(ji,1,3) ) +&
                             g2dmu(ji) * ( gu(ji,3,2) + gu(ji,2,3) ) +          &
                             g3dmu(ji) * ( gu(ji,3,3) + gu(ji,3,3) ) +          &
                               dmu(ji) * ( damp(ji) * g11v(ji,4) + g22v(ji,4) ) )

!=======================================================================================================!

!.... Energy equation

        fact1 = gamma * rhoinv(ji) / (Pr * Re)
        fact2 = gamma * gamma1 * Ma**2 * rhoinv(ji) / Re

        D(ji,5,1) = Ma**2*rhoinv(ji)**2*( L2 - pt5*gamma1*(L5 + L1) ) +         &
                    Ma**2*pt5*rhoinv(ji)*gamma1*(un+c)*(gtn*gmsinv+c*gunn) +    &
                    fact1 * rhoinv(ji) * ( damp(ji) * g1con(ji) * gt(ji,1) +    &
                                      g2con(ji) * gt(ji,2) +                    &
                                      g3con(ji) * gt(ji,3) +                    &
                                      con(ji) * ( damp(ji) * g11v(ji,5) +       &
                                              g22v(ji,5) ) ) +                  &
                    fact2 * rhoinv(ji) * ( lm(ji) * ( damp(ji) * gu(ji,1,1) +   &
                                    gu(ji,2,2) + gu(ji,3,3) )**2 ) +            &
                    two * fact2 * rhoinv(ji) * mu(ji) * (                       &
                    ( pt5 * damp(ji) * ( gu(ji,1,1) + gu(ji,1,1) ) )**2 +       &
                    ( pt5 * ( gu(ji,1,2) + damp(ji) * gu(ji,2,1) ) )**2 +       &
                    ( pt5 * ( gu(ji,1,3) + damp(ji) * gu(ji,3,1) ) )**2 +       &
                    ( pt5 * ( damp(ji) * gu(ji,2,1) + gu(ji,1,2) ) )**2 +       &
                    ( pt5 * ( gu(ji,2,2) + gu(ji,2,2) ) )**2 +                  &
                    ( pt5 * ( gu(ji,2,3) + gu(ji,2,3) ) )**2 +                  &
                    ( pt5 * ( damp(ji) * gu(ji,3,1) + gu(ji,1,3) ) )**2 +       &
                    ( pt5 * ( gu(ji,3,2) + gu(ji,2,3) ) )**2 +                  &
                    ( pt5 * ( gu(ji,3,3) + gu(ji,3,3) ) )**2 )

        D(ji,5,2) = Ma**2*pt5*rhoinv(ji)*gamma1*bn1*(gpn+rho(ji)*c*gunn) + bs1*gts
        D(ji,5,3) = Ma**2*pt5*rhoinv(ji)*gamma1*bn2*(gpn+rho(ji)*c*gunn) + bs2*gts
        D(ji,5,4) = gt(ji,3)

        D(ji,5,5) = gamma1 * ( guss + gu(ji,3,3) ) +                            &
                    Ma**2*pt5*rhoinv(ji)*gamma1*( dc*(gpn+rho(ji)*c*gunn) +     &
                    (un+c)*(gmsinv*grhon + rho(ji)*dc*gunn) ) -                 &
                   fact1 * (damp(ji) * g1dcon(ji) * gt(ji,1) +                  &
                            g2dcon(ji) * gt(ji,2) +                             &
                            g3dcon(ji) * gt(ji,3) +                             &
                            dcon(ji) * (damp(ji) * g11v(ji,5) + g22v(ji,5)) ) - &
                   fact2 * dlm(ji) * ( damp(ji) * gu(ji,1,1) +                  &
                                   gu(ji,2,2) +                                 &
                                   gu(ji,3,3) )**2 -                            &
                   two * fact2 * dmu(ji) * (                                    &
                    ( pt5 * damp(ji) * ( gu(ji,1,1) + gu(ji,1,1) ) )**2 +       &
                    ( pt5 * ( gu(ji,1,2) + damp(ji) * gu(ji,2,1) ) )**2 +       &
                    ( pt5 * ( gu(ji,1,3) + damp(ji) * gu(ji,3,1) ) )**2 +       &
                    ( pt5 * ( damp(ji) * gu(ji,2,1) + gu(ji,1,2) ) )**2 +       &
                    ( pt5 * ( gu(ji,2,2) + gu(ji,2,2) ) )**2 +                  &
                    ( pt5 * ( gu(ji,2,3) + gu(ji,2,3) ) )**2 +                  &
                    ( pt5 * ( damp(ji) * gu(ji,3,1) + gu(ji,1,3) ) )**2 +       &
                    ( pt5 * ( gu(ji,3,2) + gu(ji,2,3) ) )**2 +                  &
                    ( pt5 * ( gu(ji,3,3) + gu(ji,3,3) ) )**2 )

        else

!=======================================================================================================!
!       O u t f l o w
!=======================================================================================================!

!.... Continuity equation

        gmsinv = one / (gamma * Ma**2)

        D(ji,1,1) = one/c**2*( un*(-gtn*gmsinv) +                               &
                    pt5*( (un+c)*(gtn*gmsinv+c*gunn) + rk*t(ji)*gmsinv ) ) +    &
                    guss + gu(ji,3,3)

        D(ji,1,2) = one/c**2*( bn1*(c**2*grhon-gpn) +                           &
                    pt5*bn1*(gpn+rho(ji)*c*gunn) ) + grhos*bs1

        D(ji,1,3) = one/c**2*( bn2*(c**2*grhon-gpn) +                           &
                    pt5*bn2*(gpn+rho(ji)*c*gunn) ) + grhos*bs2

        D(ji,1,4) = grho(ji,3)

        D(ji,1,5) = -two/c**3*dc*( L2 + pt5*( L5 + L1 ) ) +                     &
                    one/c**2*( un*(two*c*dc*grhon - grhon*gmsinv) +             &
                    pt5*( dc*(gpn+rho(ji)*c*gunn) + (un+c)*(grhon*gmsinv +      &
                    rho(ji)*dc*gunn) + rk*rho(ji)*gmsinv ) )

!=======================================================================================================!

!.... Momentum equation -- x_1

        fact1  = rhoinv(ji) / Re

        D(ji,2,1) =                                                                     &
           -bn1/(two*rho(ji)**2*c)*( L5 - L1 ) +                                        &
            bn1/(two*rho(ji)*c)*( (un+c)*( gtn*gmsinv + c*gunn ) - rk*t(ji)*gmsinv ) -  &
            bs1*rhoinv(ji)**2*gps + bs1*rhoinv(ji)*gmsinv*gts +                         &
            rhoinv(ji) * (                                                              &
            fact1 * damp(ji) * ( g1lm(ji) * divu(ji) + lm(ji) * g1divu(ji) ) +          &
            fact1 * ( damp(ji) * g1mu(ji) * ( two * gu(ji,1,1) ) +                      &
                            g2mu(ji) * ( gu(ji,1,2) + damp(ji) * gu(ji,2,1) ) +         &
                            g3mu(ji) * ( gu(ji,1,3) + damp(ji) * gu(ji,3,1) ) +         &
                              mu(ji) * ( damp(ji) * two * g11v(ji,2) +                  &
                                      g22v(ji,2) + damp(ji) * g12v(ji,3) ) ) )

        D(ji,2,2) = bn1/(two*rho(ji)*c)*bn1*(gpn + rho(ji)*c*gunn) + bs1*bn1*gusn + bs1*gu1s
        D(ji,2,3) = bn1/(two*rho(ji)*c)*bn2*(gpn + rho(ji)*c*gunn) + bs1*bn2*gusn + bs2*gu1s
        D(ji,2,4) = gu(ji,1,3)
        D(ji,2,5) =                                                                     & 
                  -bn1*dc/(two*rho(ji)*c**2)*( L5 - L1 ) +                              &
                   bn1/(two*rho(ji)*c)*( dc*(gpn + rho(ji)*c*gunn) +                    &
                   (un+c)*(grhon*gmsinv + rho(ji)*dc*gunn) - rk*rho(ji)*gmsinv ) +      &
                   bs1*rhoinv(ji)*gmsinv*grhos -                                        &
                   fact1 * damp(ji) * ( g1dlm(ji) * divu(ji) + dlm(ji) * g1divu(ji) ) - &
                   fact1 * ( damp(ji) * g1dmu(ji) * ( gu(ji,1,1) + gu(ji,1,1) ) +       &
                             g2dmu(ji) * ( gu(ji,1,2) + damp(ji) * gu(ji,2,1) ) +       &
                             g3dmu(ji) * ( gu(ji,1,3) + damp(ji) * gu(ji,3,1) ) +       &
                               dmu(ji) * ( damp(ji) * two * g11v(ji,2) +                &
                                       g22v(ji,2) + damp(ji) * g12v(ji,3) ) )

!=======================================================================================================!

!.... Momentum equation -- x_2

        D(ji,3,1) =                                                                     &
           -bn2/(two*rho(ji)**2*c)*( L5 - L1 ) +                                        &
            bn2/(two*rho(ji)*c)*( (un+c)*( gtn*gmsinv + c*gunn ) - rk*t(ji)*gmsinv ) -  &
            bs2*rhoinv(ji)**2*gps + bs2*rhoinv(ji)*gmsinv*gts +                         &
            rhoinv(ji) * (                                                              &
            fact1 * ( g2lm(ji) * (damp(ji) * gu(ji,1,1) + gu(ji,2,2) +                  &
                gu(ji,3,3)) + lm(ji) * (damp(ji) * g12v(ji,2) + g22v(ji,3)) ) +         &
            fact1 * ( damp(ji) * g1mu(ji) * (gu(ji,2,1) + gu(ji,1,2) ) +                &
                      g2mu(ji) * ( gu(ji,2,2) + gu(ji,2,2) ) +                          &
                      g3mu(ji) * ( gu(ji,2,3) + gu(ji,3,2) ) +                          &
                        mu(ji) * ( damp(ji) * g11v(ji,3) + damp(ji) * g12v(ji,2) +      &
                                g22v(ji,3) + g22v(ji,3) ) ) )

        D(ji,3,2) = bn2/(two*rho(ji)*c)*bn1*(gpn + rho(ji)*c*gunn) + bs2*bn1*gusn + bs1*gu2s 
        D(ji,3,3) = bn2/(two*rho(ji)*c)*bn2*(gpn + rho(ji)*c*gunn) + bs2*bn2*gusn + bs2*gu2s
        D(ji,3,4) = gu(ji,2,3)
        D(ji,3,5) =                                                                     & 
                  -bn2*dc/(two*rho(ji)*c**2)*( L5 - L1 ) +                              &
                   bn2/(two*rho(ji)*c)*( dc*(gpn + rho(ji)*c*gunn) +                    &
                   (un+c)*(grhon*gmsinv + rho(ji)*dc*gunn) - rk*rho(ji)*gmsinv ) +      &
                   bs2*rhoinv(ji)*gmsinv*grhos -                                        &
                   fact1 * ( g2dlm(ji) * (damp(ji) * gu(ji,1,1) +                       &
                                      gu(ji,2,2) + gu(ji,3,3)) +                        &
                             dlm(ji) * (damp(ji) * g12v(ji,2) + g22v(ji,3)) ) -         &
                   fact1 * ( damp(ji) * g1dmu(ji) * ( gu(ji,2,1) + gu(ji,1,2) ) +       &
                             g2dmu(ji) * ( gu(ji,2,2) + gu(ji,2,2) ) +                  &
                             g3dmu(ji) * ( gu(ji,2,3) + gu(ji,3,2) ) +                  &
                               dmu(ji) * ( damp(ji) * g11v(ji,3) +                      &
                                       damp(ji) * g12v(ji,2) +                          &
                                       g22v(ji,3) + g22v(ji,3) ) )

!=======================================================================================================!

!.... Momentum equation -- x_3

        D(ji,4,1) = rhoinv(ji) * ( gmsinv * gt(ji,3) - rhoinv(ji) * gp(ji,3) +  &
            fact1 * ( g3lm(ji) * (damp(ji) * gu(ji,1,1) + gu(ji,2,2) +          &
                      gu(ji,3,3)) + lm(ji) * g3divu(ji) ) +                     &
            fact1 * ( damp(ji) * g1mu(ji) * (gu(ji,3,1) + gu(ji,1,3) ) +        &
                      g2mu(ji) * ( gu(ji,3,2) + gu(ji,2,3) ) +                  &
                      g3mu(ji) * ( gu(ji,3,3) + gu(ji,3,3) ) +                  &
                        mu(ji) * ( damp(ji) * g11v(ji,4) + g22v(ji,4) ) ) )

        D(ji,4,2) = bn1 * gu3n + bs1 * gu3s
        D(ji,4,3) = bn2 * gu3n + bs2 * gu3s
        D(ji,4,4) = gu(ji,3,3)
        D(ji,4,5) = rhoinv(ji) * gmsinv * grho(ji,3)  -                         &
                   fact1 * ( g3dlm(ji) * (damp(ji) * gu(ji,1,1) +               &
                                      gu(ji,2,2) + gu(ji,3,3)) +                &
                               dlm(ji) * g3divu(ji) ) -                         &
                   fact1 * ( damp(ji) * g1dmu(ji) * (gu(ji,3,1) + gu(ji,1,3) ) +&
                             g2dmu(ji) * ( gu(ji,3,2) + gu(ji,2,3) ) +          &
                             g3dmu(ji) * ( gu(ji,3,3) + gu(ji,3,3) ) +          &
                               dmu(ji) * ( damp(ji) * g11v(ji,4) + g22v(ji,4) ) )

!=======================================================================================================!

!.... Energy equation

        fact1 = gamma * rhoinv(ji) / (Pr * Re)
        fact2 = gamma * gamma1 * Ma**2 * rhoinv(ji) / Re

        D(ji,5,1) = Ma**2*rhoinv(ji)**2*( L2 - pt5*gamma1*(L5 + L1) ) -         &
                    Ma**2*rhoinv(ji)*( un*(-gtn*gmsinv) - pt5*gamma1*           &
                    ((un+c)*(gtn*gmsinv+c*gunn) + rk*t(ji)*gmsinv) ) +          &
                    fact1 * rhoinv(ji) * ( damp(ji) * g1con(ji) * gt(ji,1) +    &
                                      g2con(ji) * gt(ji,2) +                    &
                                      g3con(ji) * gt(ji,3) +                    &
                                      con(ji) * ( damp(ji) * g11v(ji,5) +       &
                                              g22v(ji,5) ) ) +                  &
                    fact2 * rhoinv(ji) * ( lm(ji) * ( damp(ji) * gu(ji,1,1) +   &
                                    gu(ji,2,2) + gu(ji,3,3) )**2 ) +            &
                    two * fact2 * rhoinv(ji) * mu(ji) * (                       &
                    ( pt5 * damp(ji) * ( gu(ji,1,1) + gu(ji,1,1) ) )**2 +       &
                    ( pt5 * ( gu(ji,1,2) + damp(ji) * gu(ji,2,1) ) )**2 +       &
                    ( pt5 * ( gu(ji,1,3) + damp(ji) * gu(ji,3,1) ) )**2 +       &
                    ( pt5 * ( damp(ji) * gu(ji,2,1) + gu(ji,1,2) ) )**2 +       &
                    ( pt5 * ( gu(ji,2,2) + gu(ji,2,2) ) )**2 +                  &
                    ( pt5 * ( gu(ji,2,3) + gu(ji,2,3) ) )**2 +                  &
                    ( pt5 * ( damp(ji) * gu(ji,3,1) + gu(ji,1,3) ) )**2 +       &
                    ( pt5 * ( gu(ji,3,2) + gu(ji,2,3) ) )**2 +                  &
                    ( pt5 * ( gu(ji,3,3) + gu(ji,3,3) ) )**2 )

        D(ji,5,2) = -Ma**2*rhoinv(ji)*( bn1*(c**2*grhon-gpn) -                  &
                     pt5*gamma1*bn1*(gpn+rho(ji)*c*gunn) ) + bs1*gts

        D(ji,5,3) = -Ma**2*rhoinv(ji)*( bn2*(c**2*grhon-gpn) -                  &
                     pt5*gamma1*bn2*(gpn+rho(ji)*c*gunn) ) + bs2*gts

        D(ji,5,4) = gt(ji,3)

        D(ji,5,5) = gamma1 * ( guss + gu(ji,3,3) ) -                            &
                    Ma**2*rhoinv(ji)*( un*(two*c*dc*grhon-grhon*gmsinv) -       &
                    pt5*gamma1*( dc*(gpn+rho(ji)*c*gunn) +                      &
                    (un+c)*(grhon*gmsinv + rho(ji)*dc*gunn) +                   &
                    rk*rho(ji)*gmsinv ) ) -                                     &
                   fact1 * (damp(ji) * g1dcon(ji) * gt(ji,1) +                  &
                            g2dcon(ji) * gt(ji,2) +                             &
                            g3dcon(ji) * gt(ji,3) +                             &
                            dcon(ji) * (damp(ji) * g11v(ji,5) + g22v(ji,5)) ) - &
                   fact2 * dlm(ji) * ( damp(ji) * gu(ji,1,1) +                  &
                                   gu(ji,2,2) +                                 &
                                   gu(ji,3,3) )**2 -                            &
                   two * fact2 * dmu(ji) * (                                    &
                    ( pt5 * damp(ji) * ( gu(ji,1,1) + gu(ji,1,1) ) )**2 +       &
                    ( pt5 * ( gu(ji,1,2) + damp(ji) * gu(ji,2,1) ) )**2 +       &
                    ( pt5 * ( gu(ji,1,3) + damp(ji) * gu(ji,3,1) ) )**2 +       &
                    ( pt5 * ( damp(ji) * gu(ji,2,1) + gu(ji,1,2) ) )**2 +       &
                    ( pt5 * ( gu(ji,2,2) + gu(ji,2,2) ) )**2 +                  &
                    ( pt5 * ( gu(ji,2,3) + gu(ji,2,3) ) )**2 +                  &
                    ( pt5 * ( damp(ji) * gu(ji,3,1) + gu(ji,1,3) ) )**2 +       &
                    ( pt5 * ( gu(ji,3,2) + gu(ji,2,3) ) )**2 +                  &
                    ( pt5 * ( gu(ji,3,3) + gu(ji,3,3) ) )**2 )

        end if          ! inflow/outflow

          end do
          end if
        end do
        end if
        
!=======================================================================================================!
!       T o p   B o u n d a r y
!=======================================================================================================!
        if (top.eq.1) then

        j = ny
        do i = 1, nx-1
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
        
!.... compute the potential solution

            b = one / pi
            if (Ma.lt.one) then
              rbeta = sqrt(one - Ma**2)
            else
              rbeta = sqrt(Ma**2 - one)
            end if
            ub1 = one + (x(j,1)-b)/pi/rbeta**2/((x(j,i)-b)**2 + &
                         rbeta**2*y(j,i)**2)
            ub2 = rbeta*y(j,i)/pi/((x(j,i)-b)**2+rbeta**2*y(j,i)**2)

            gub11 = one/(pi*rbeta**2*((x(j,i)-b)**2+rbeta**2*y(j,i)**2)) - &
                    two*(x(j,i)-b)**2/(pi*rbeta**2*((x(j,i)-b)**2+rbeta**2*y(j,i)**2)**2)

            gub12 = -two*(x(j,i)-b)*rbeta**2*y(j,i)/(pi*rbeta**2*((x(j,i)-b)**2+ &
                     rbeta**2*y(j,i)**2)**2)

            gub21 = -rbeta*y(j,i)*two*(x(j,i)-b)/(pi*((x(j,i)-b)**2+ &
                     rbeta**2*y(j,i)**2)**2)

            gub22 = rbeta/(pi*((x(j,i)-b)**2+rbeta**2*y(j,i)**2)) - &
                    two*rbeta**3*y(j,i)**2/(pi*((x(j,i)-b)**2+ &
                    rbeta**2*y(j,i)**2)**2)

            ubn = bn1*ub1 + bn2*ub2
            ubs = bs1*ub1 + bs2*ub2

            gubnn = bn1*bn1 * gub11 + bn2*bn1 * gub21 + &
                    bn1*bn2 * gub12 + bn2*bn2 * gub22

            gubsn = bs1*bn1 * gub11 + bs2*bn1 * gub21 + &
                    bs1*bn2 * gub12 + bs2*bn2 * gub22

            tb   = one - pt5*gamma1*Ma**2*( ub1**2 + ub2**2 - one )
            rhob = tb**(one/gamma1)
            pb   = rhob * tb / ( gamma * Ma**2 )
            cb   = sqrt(tb)/Ma

            gtbn = bn1 * (-gamma1*Ma**2*( ub1*gub11 + ub2*gub21 ) ) + &
                   bn2 * (-gamma1*Ma**2*( ub1*gub12 + ub2*gub22 ) )
            
            grhobn = one/gamma1 * tb**((two-gamma)/gamma1) * gtbn
            gpbn   = one/(gamma*Ma**2) * ( tb*grhobn + rhob*gtbn )

!.... nonreflecting boundary conditions
        
            if (un .le. zero ) then
              L1 = (ubn - cb)*(gpbn - rhob*cb*gubnn)
              L2 = ubn*(cb**2*grhobn - gpbn)
              L3 = ubn*gubsn
              L4 = zero

!             L1 = L5
!             L2 = pt5*gamma1*(L5 + L1)
!             L3 = zero
!             L4 = zero
            else
              L1 = rk * ( p(ji) - pinf ) + (ubn - cb)*(gpbn - rhob*cb*gubnn)    
            end if

        if (un .le. zero) then

!=======================================================================================================!
!       I n f l o w
!=======================================================================================================!

!.... Continuity equation

        gmsinv = one / (gamma * Ma**2)

        D(ji,1,1) = pt5/c**2*(un+c)*(gtn*gmsinv + c*gunn) + guss + gu(ji,3,3)
        D(ji,1,2) = pt5/c**2*bn1*(gpn + rho(ji)*c*gunn) + grhos*bs1
        D(ji,1,3) = pt5/c**2*bn2*(gpn + rho(ji)*c*gunn) + grhos*bs2
        D(ji,1,4) = grho(ji,3)
        D(ji,1,5) = -two/c**3*dc*(L2 + pt5*(L5 + L1)) +         &
                    pt5/c**2*( dc*(gpn+rho(ji)*c*gunn) +        &
                    (un+c)*(grhon*gmsinv + rho(ji)*dc*gunn) )

!.... set u,v,w,T

!       D(ji,1,1) = gamma/c**2*(un+c)*(gtn*gmsinv + c*gunn) + guss + gu(ji,3,3)
!       D(ji,1,2) = gamma/c**2*bn1*(gpn + rho(ji)*c*gunn) + grhos*bs1
!       D(ji,1,3) = gamma/c**2*bn2*(gpn + rho(ji)*c*gunn) + grhos*bs2
!       D(ji,1,4) = grho(ji,3)
!       D(ji,1,5) = -two/c**3*dc*(L2 + pt5*(L5 + L1)) +         &
!                    gamma/c**2*( dc*(gpn+rho(ji)*c*gunn) +     &
!                    (un+c)*(grhon*gmsinv + rho(ji)*dc*gunn) )

!=======================================================================================================!

!.... Momentum equation -- x_1

        fact1  = rhoinv(ji) / Re

        D(ji,2,1) =                                                                     &
           -bn1/(two*rho(ji)**2*c)*( L5 - L1 ) +                                        &
            bn1/(two*rho(ji)*c)*(un+c)*( gtn*gmsinv + c*gunn ) -                        &
            bs1*rhoinv(ji)**2*gps + bs1*rhoinv(ji)*gmsinv*gts +                         &
            rhoinv(ji) * (                                                              &
            fact1 * damp(ji) * ( g1lm(ji) * divu(ji) + lm(ji) * g1divu(ji) ) +          &
            fact1 * ( damp(ji) * g1mu(ji) * ( two * gu(ji,1,1) ) +                      &
                            g2mu(ji) * ( gu(ji,1,2) + damp(ji) * gu(ji,2,1) ) +         &
                            g3mu(ji) * ( gu(ji,1,3) + damp(ji) * gu(ji,3,1) ) +         &
                              mu(ji) * ( damp(ji) * two * g11v(ji,2) +                  &
                                      g22v(ji,2) + damp(ji) * g12v(ji,3) ) ) )

        D(ji,2,2) = bn1/(two*rho(ji)*c)*bn1*(gpn + rho(ji)*c*gunn) + bs1 * gu1s
        D(ji,2,3) = bn1/(two*rho(ji)*c)*bn2*(gpn + rho(ji)*c*gunn) + bs2 * gu1s
        D(ji,2,4) = gu(ji,1,3)
        D(ji,2,5) =                                                                     & 
                  -bn1*dc/(two*rho(ji)*c**2)*( L5 - L1 ) +                              &
                   bn1/(two*rho(ji)*c)*( dc*(gpn + rho(ji)*c*gunn) +                    &
                   (un+c)*(grhon*gmsinv + rho(ji)*dc*gunn) ) +                          &
                   bs1*rhoinv(ji)*gmsinv*grhos -                                        &
                   fact1 * damp(ji) * ( g1dlm(ji) * divu(ji) + dlm(ji) * g1divu(ji) ) - &
                   fact1 * ( damp(ji) * g1dmu(ji) * ( gu(ji,1,1) + gu(ji,1,1) ) +       &
                             g2dmu(ji) * ( gu(ji,1,2) + damp(ji) * gu(ji,2,1) ) +       &
                             g3dmu(ji) * ( gu(ji,1,3) + damp(ji) * gu(ji,3,1) ) +       &
                               dmu(ji) * ( damp(ji) * two * g11v(ji,2) +                &
                                       g22v(ji,2) + damp(ji) * g12v(ji,3) ) )

!=======================================================================================================!

!.... Momentum equation -- x_2

        D(ji,3,1) =                                                                     &
           -bn2/(two*rho(ji)**2*c)*( L5 - L1 ) +                                        &
            bn2/(two*rho(ji)*c)*(un+c)*(gtn*gmsinv + c*gunn) -                          &
            bs2*rhoinv(ji)**2*gps + bs2*rhoinv(ji)*gmsinv*gts +                         &
            rhoinv(ji) * (                                                              &
            fact1 * ( g2lm(ji) * (damp(ji) * gu(ji,1,1) + gu(ji,2,2) +                  &
                gu(ji,3,3)) + lm(ji) * (damp(ji) * g12v(ji,2) + g22v(ji,3)) ) +         &
            fact1 * ( damp(ji) * g1mu(ji) * (gu(ji,2,1) + gu(ji,1,2) ) +                &
                      g2mu(ji) * ( gu(ji,2,2) + gu(ji,2,2) ) +                          &
                      g3mu(ji) * ( gu(ji,2,3) + gu(ji,3,2) ) +                          &
                        mu(ji) * ( damp(ji) * g11v(ji,3) + damp(ji) * g12v(ji,2) +      &
                                g22v(ji,3) + g22v(ji,3) ) ) )

        D(ji,3,2) = bn2/(two*rho(ji)*c)*bn1*(gpn + rho(ji)*c*gunn) + bs1 * gu2s 
        D(ji,3,3) = bn2/(two*rho(ji)*c)*bn2*(gpn + rho(ji)*c*gunn) + bs2 * gu2s
        D(ji,3,4) = gu(ji,2,3)
        D(ji,3,5) =                                                                     & 
                  -bn2*dc/(two*rho(ji)*c**2)*( L5 - L1 ) +                              &
                   bn2/(two*rho(ji)*c)*( dc*(gpn + rho(ji)*c*gunn) +                    &
                   (un+c)*(grhon*gmsinv + rho(ji)*dc*gunn) ) +                          &
                   bs2*rhoinv(ji)*gmsinv*grhos -                                        &
                   fact1 * ( g2dlm(ji) * (damp(ji) * gu(ji,1,1) +                       &
                                      gu(ji,2,2) + gu(ji,3,3)) +                        &
                             dlm(ji) * (damp(ji) * g12v(ji,2) + g22v(ji,3)) ) -         &
                   fact1 * ( damp(ji) * g1dmu(ji) * ( gu(ji,2,1) + gu(ji,1,2) ) +       &
                             g2dmu(ji) * ( gu(ji,2,2) + gu(ji,2,2) ) +                  &
                             g3dmu(ji) * ( gu(ji,2,3) + gu(ji,3,2) ) +                  &
                               dmu(ji) * ( damp(ji) * g11v(ji,3) +                      &
                                       damp(ji) * g12v(ji,2) +                          &
                                       g22v(ji,3) + g22v(ji,3) ) )

!=======================================================================================================!

!.... Momentum equation -- x_3

        D(ji,4,1) = rhoinv(ji) * ( gmsinv * gt(ji,3) - rhoinv(ji) * gp(ji,3) +  &
            fact1 * ( g3lm(ji) * (damp(ji) * gu(ji,1,1) + gu(ji,2,2) +          &
                      gu(ji,3,3)) + lm(ji) * g3divu(ji) ) +                     &
            fact1 * ( damp(ji) * g1mu(ji) * (gu(ji,3,1) + gu(ji,1,3) ) +        &
                      g2mu(ji) * ( gu(ji,3,2) + gu(ji,2,3) ) +                  &
                      g3mu(ji) * ( gu(ji,3,3) + gu(ji,3,3) ) +                  &
                        mu(ji) * ( damp(ji) * g11v(ji,4) + g22v(ji,4) ) ) )

        D(ji,4,2) = bs1 * gu3s
        D(ji,4,3) = bs2 * gu3s
        D(ji,4,4) = gu(ji,3,3)
        D(ji,4,5) = rhoinv(ji) * gmsinv * grho(ji,3)  -                         &
                   fact1 * ( g3dlm(ji) * (damp(ji) * gu(ji,1,1) +               &
                                      gu(ji,2,2) + gu(ji,3,3)) +                &
                               dlm(ji) * g3divu(ji) ) -                         &
                   fact1 * ( damp(ji) * g1dmu(ji) * (gu(ji,3,1) + gu(ji,1,3) ) +&
                             g2dmu(ji) * ( gu(ji,3,2) + gu(ji,2,3) ) +          &
                             g3dmu(ji) * ( gu(ji,3,3) + gu(ji,3,3) ) +          &
                               dmu(ji) * ( damp(ji) * g11v(ji,4) + g22v(ji,4) ) )

!=======================================================================================================!

!.... Energy equation

        fact1 = gamma * rhoinv(ji) / (Pr * Re)
        fact2 = gamma * gamma1 * Ma**2 * rhoinv(ji) / Re

        D(ji,5,1) = Ma**2*rhoinv(ji)**2*( L2 - pt5*gamma1*(L5 + L1) ) +         &
                    Ma**2*pt5*rhoinv(ji)*gamma1*(un+c)*(gtn*gmsinv+c*gunn) +    &
                    fact1 * rhoinv(ji) * ( damp(ji) * g1con(ji) * gt(ji,1) +    &
                                      g2con(ji) * gt(ji,2) +                    &
                                      g3con(ji) * gt(ji,3) +                    &
                                      con(ji) * ( damp(ji) * g11v(ji,5) +       &
                                              g22v(ji,5) ) ) +                  &
                    fact2 * rhoinv(ji) * ( lm(ji) * ( damp(ji) * gu(ji,1,1) +   &
                                    gu(ji,2,2) + gu(ji,3,3) )**2 ) +            &
                    two * fact2 * rhoinv(ji) * mu(ji) * (                       &
                    ( pt5 * damp(ji) * ( gu(ji,1,1) + gu(ji,1,1) ) )**2 +       &
                    ( pt5 * ( gu(ji,1,2) + damp(ji) * gu(ji,2,1) ) )**2 +       &
                    ( pt5 * ( gu(ji,1,3) + damp(ji) * gu(ji,3,1) ) )**2 +       &
                    ( pt5 * ( damp(ji) * gu(ji,2,1) + gu(ji,1,2) ) )**2 +       &
                    ( pt5 * ( gu(ji,2,2) + gu(ji,2,2) ) )**2 +                  &
                    ( pt5 * ( gu(ji,2,3) + gu(ji,2,3) ) )**2 +                  &
                    ( pt5 * ( damp(ji) * gu(ji,3,1) + gu(ji,1,3) ) )**2 +       &
                    ( pt5 * ( gu(ji,3,2) + gu(ji,2,3) ) )**2 +                  &
                    ( pt5 * ( gu(ji,3,3) + gu(ji,3,3) ) )**2 )

        D(ji,5,2) = Ma**2*pt5*rhoinv(ji)*gamma1*bn1*(gpn+rho(ji)*c*gunn) + bs1*gts
        D(ji,5,3) = Ma**2*pt5*rhoinv(ji)*gamma1*bn2*(gpn+rho(ji)*c*gunn) + bs2*gts
        D(ji,5,4) = gt(ji,3)

        D(ji,5,5) = gamma1 * ( guss + gu(ji,3,3) ) +                            &
                    Ma**2*pt5*rhoinv(ji)*gamma1*( dc*(gpn+rho(ji)*c*gunn) +     &
                    (un+c)*(gmsinv*grhon + rho(ji)*dc*gunn) ) -                 &
                   fact1 * (damp(ji) * g1dcon(ji) * gt(ji,1) +                  &
                            g2dcon(ji) * gt(ji,2) +                             &
                            g3dcon(ji) * gt(ji,3) +                             &
                            dcon(ji) * (damp(ji) * g11v(ji,5) + g22v(ji,5)) ) - &
                   fact2 * dlm(ji) * ( damp(ji) * gu(ji,1,1) +                  &
                                   gu(ji,2,2) +                                 &
                                   gu(ji,3,3) )**2 -                            &
                   two * fact2 * dmu(ji) * (                                    &
                    ( pt5 * damp(ji) * ( gu(ji,1,1) + gu(ji,1,1) ) )**2 +       &
                    ( pt5 * ( gu(ji,1,2) + damp(ji) * gu(ji,2,1) ) )**2 +       &
                    ( pt5 * ( gu(ji,1,3) + damp(ji) * gu(ji,3,1) ) )**2 +       &
                    ( pt5 * ( damp(ji) * gu(ji,2,1) + gu(ji,1,2) ) )**2 +       &
                    ( pt5 * ( gu(ji,2,2) + gu(ji,2,2) ) )**2 +                  &
                    ( pt5 * ( gu(ji,2,3) + gu(ji,2,3) ) )**2 +                  &
                    ( pt5 * ( damp(ji) * gu(ji,3,1) + gu(ji,1,3) ) )**2 +       &
                    ( pt5 * ( gu(ji,3,2) + gu(ji,2,3) ) )**2 +                  &
                    ( pt5 * ( gu(ji,3,3) + gu(ji,3,3) ) )**2 )

        else

!=======================================================================================================!
!       O u t f l o w
!=======================================================================================================!

!.... Continuity equation

        gmsinv = one / (gamma * Ma**2)

        D(ji,1,1) = one/c**2*( un*(-gtn*gmsinv) +                               &
                    pt5*( (un+c)*(gtn*gmsinv+c*gunn) + rk*t(ji)*gmsinv ) ) +    &
                    guss + gu(ji,3,3)

        D(ji,1,2) = one/c**2*( bn1*(c**2*grhon-gpn) +                           &
                    pt5*bn1*(gpn+rho(ji)*c*gunn) ) + grhos*bs1

        D(ji,1,3) = one/c**2*( bn2*(c**2*grhon-gpn) +                           &
                    pt5*bn2*(gpn+rho(ji)*c*gunn) ) + grhos*bs2

        D(ji,1,4) = grho(ji,3)

        D(ji,1,5) = -two/c**3*dc*( L2 + pt5*( L5 + L1 ) ) +                     &
                    one/c**2*( un*(two*c*dc*grhon - grhon*gmsinv) +             &
                    pt5*( dc*(gpn+rho(ji)*c*gunn) + (un+c)*(grhon*gmsinv +      &
                    rho(ji)*dc*gunn) + rk*rho(ji)*gmsinv ) )

!=======================================================================================================!

!.... Momentum equation -- x_1

        fact1  = rhoinv(ji) / Re

        D(ji,2,1) =                                                                     &
           -bn1/(two*rho(ji)**2*c)*( L5 - L1 ) +                                        &
            bn1/(two*rho(ji)*c)*( (un+c)*( gtn*gmsinv + c*gunn ) - rk*t(ji)*gmsinv ) -  &
            bs1*rhoinv(ji)**2*gps + bs1*rhoinv(ji)*gmsinv*gts +                         &
            rhoinv(ji) * (                                                              &
            fact1 * damp(ji) * ( g1lm(ji) * divu(ji) + lm(ji) * g1divu(ji) ) +          &
            fact1 * ( damp(ji) * g1mu(ji) * ( two * gu(ji,1,1) ) +                      &
                            g2mu(ji) * ( gu(ji,1,2) + damp(ji) * gu(ji,2,1) ) +         &
                            g3mu(ji) * ( gu(ji,1,3) + damp(ji) * gu(ji,3,1) ) +         &
                              mu(ji) * ( damp(ji) * two * g11v(ji,2) +                  &
                                      g22v(ji,2) + damp(ji) * g12v(ji,3) ) ) )

        D(ji,2,2) = bn1/(two*rho(ji)*c)*bn1*(gpn + rho(ji)*c*gunn) + bs1*bn1*gusn + bs1*gu1s
        D(ji,2,3) = bn1/(two*rho(ji)*c)*bn2*(gpn + rho(ji)*c*gunn) + bs1*bn2*gusn + bs2*gu1s
        D(ji,2,4) = gu(ji,1,3)
        D(ji,2,5) =                                                                     & 
                  -bn1*dc/(two*rho(ji)*c**2)*( L5 - L1 ) +                              &
                   bn1/(two*rho(ji)*c)*( dc*(gpn + rho(ji)*c*gunn) +                    &
                   (un+c)*(grhon*gmsinv + rho(ji)*dc*gunn) - rk*rho(ji)*gmsinv ) +      &
                   bs1*rhoinv(ji)*gmsinv*grhos -                                        &
                   fact1 * damp(ji) * ( g1dlm(ji) * divu(ji) + dlm(ji) * g1divu(ji) ) - &
                   fact1 * ( damp(ji) * g1dmu(ji) * ( gu(ji,1,1) + gu(ji,1,1) ) +       &
                             g2dmu(ji) * ( gu(ji,1,2) + damp(ji) * gu(ji,2,1) ) +       &
                             g3dmu(ji) * ( gu(ji,1,3) + damp(ji) * gu(ji,3,1) ) +       &
                               dmu(ji) * ( damp(ji) * two * g11v(ji,2) +                &
                                       g22v(ji,2) + damp(ji) * g12v(ji,3) ) )

!=======================================================================================================!

!.... Momentum equation -- x_2

        D(ji,3,1) =                                                                     &
           -bn2/(two*rho(ji)**2*c)*( L5 - L1 ) +                                        &
            bn2/(two*rho(ji)*c)*( (un+c)*( gtn*gmsinv + c*gunn ) - rk*t(ji)*gmsinv ) -  &
            bs2*rhoinv(ji)**2*gps + bs2*rhoinv(ji)*gmsinv*gts +                         &
            rhoinv(ji) * (                                                              &
            fact1 * ( g2lm(ji) * (damp(ji) * gu(ji,1,1) + gu(ji,2,2) +                  &
                gu(ji,3,3)) + lm(ji) * (damp(ji) * g12v(ji,2) + g22v(ji,3)) ) +         &
            fact1 * ( damp(ji) * g1mu(ji) * (gu(ji,2,1) + gu(ji,1,2) ) +                &
                      g2mu(ji) * ( gu(ji,2,2) + gu(ji,2,2) ) +                          &
                      g3mu(ji) * ( gu(ji,2,3) + gu(ji,3,2) ) +                          &
                        mu(ji) * ( damp(ji) * g11v(ji,3) + damp(ji) * g12v(ji,2) +      &
                                g22v(ji,3) + g22v(ji,3) ) ) )

        D(ji,3,2) = bn2/(two*rho(ji)*c)*bn1*(gpn + rho(ji)*c*gunn) + bs2*bn1*gusn + bs1*gu2s 
        D(ji,3,3) = bn2/(two*rho(ji)*c)*bn2*(gpn + rho(ji)*c*gunn) + bs2*bn2*gusn + bs2*gu2s
        D(ji,3,4) = gu(ji,2,3)
        D(ji,3,5) =                                                                     & 
                  -bn2*dc/(two*rho(ji)*c**2)*( L5 - L1 ) +                              &
                   bn2/(two*rho(ji)*c)*( dc*(gpn + rho(ji)*c*gunn) +                    &
                   (un+c)*(grhon*gmsinv + rho(ji)*dc*gunn) - rk*rho(ji)*gmsinv ) +      &
                   bs2*rhoinv(ji)*gmsinv*grhos -                                        &
                   fact1 * ( g2dlm(ji) * (damp(ji) * gu(ji,1,1) +                       &
                                      gu(ji,2,2) + gu(ji,3,3)) +                        &
                             dlm(ji) * (damp(ji) * g12v(ji,2) + g22v(ji,3)) ) -         &
                   fact1 * ( damp(ji) * g1dmu(ji) * ( gu(ji,2,1) + gu(ji,1,2) ) +       &
                             g2dmu(ji) * ( gu(ji,2,2) + gu(ji,2,2) ) +                  &
                             g3dmu(ji) * ( gu(ji,2,3) + gu(ji,3,2) ) +                  &
                               dmu(ji) * ( damp(ji) * g11v(ji,3) +                      &
                                       damp(ji) * g12v(ji,2) +                          &
                                       g22v(ji,3) + g22v(ji,3) ) )

!=======================================================================================================!

!.... Momentum equation -- x_3

        D(ji,4,1) = rhoinv(ji) * ( gmsinv * gt(ji,3) - rhoinv(ji) * gp(ji,3) +  &
            fact1 * ( g3lm(ji) * (damp(ji) * gu(ji,1,1) + gu(ji,2,2) +          &
                      gu(ji,3,3)) + lm(ji) * g3divu(ji) ) +                     &
            fact1 * ( damp(ji) * g1mu(ji) * (gu(ji,3,1) + gu(ji,1,3) ) +        &
                      g2mu(ji) * ( gu(ji,3,2) + gu(ji,2,3) ) +                  &
                      g3mu(ji) * ( gu(ji,3,3) + gu(ji,3,3) ) +                  &
                        mu(ji) * ( damp(ji) * g11v(ji,4) + g22v(ji,4) ) ) )

        D(ji,4,2) = bn1 * gu3n + bs1 * gu3s
        D(ji,4,3) = bn2 * gu3n + bs2 * gu3s
        D(ji,4,4) = gu(ji,3,3)
        D(ji,4,5) = rhoinv(ji) * gmsinv * grho(ji,3)  -                         &
                   fact1 * ( g3dlm(ji) * (damp(ji) * gu(ji,1,1) +               &
                                      gu(ji,2,2) + gu(ji,3,3)) +                &
                               dlm(ji) * g3divu(ji) ) -                         &
                   fact1 * ( damp(ji) * g1dmu(ji) * (gu(ji,3,1) + gu(ji,1,3) ) +&
                             g2dmu(ji) * ( gu(ji,3,2) + gu(ji,2,3) ) +          &
                             g3dmu(ji) * ( gu(ji,3,3) + gu(ji,3,3) ) +          &
                               dmu(ji) * ( damp(ji) * g11v(ji,4) + g22v(ji,4) ) )

!=======================================================================================================!

!.... Energy equation

        fact1 = gamma * rhoinv(ji) / (Pr * Re)
        fact2 = gamma * gamma1 * Ma**2 * rhoinv(ji) / Re

        D(ji,5,1) = Ma**2*rhoinv(ji)**2*( L2 - pt5*gamma1*(L5 + L1) ) -         &
                    Ma**2*rhoinv(ji)*( un*(-gtn*gmsinv) - pt5*gamma1*           &
                    ((un+c)*(gtn*gmsinv+c*gunn) + rk*t(ji)*gmsinv) ) +          &
                    fact1 * rhoinv(ji) * ( damp(ji) * g1con(ji) * gt(ji,1) +    &
                                      g2con(ji) * gt(ji,2) +                    &
                                      g3con(ji) * gt(ji,3) +                    &
                                      con(ji) * ( damp(ji) * g11v(ji,5) +       &
                                              g22v(ji,5) ) ) +                  &
                    fact2 * rhoinv(ji) * ( lm(ji) * ( damp(ji) * gu(ji,1,1) +   &
                                    gu(ji,2,2) + gu(ji,3,3) )**2 ) +            &
                    two * fact2 * rhoinv(ji) * mu(ji) * (                       &
                    ( pt5 * damp(ji) * ( gu(ji,1,1) + gu(ji,1,1) ) )**2 +       &
                    ( pt5 * ( gu(ji,1,2) + damp(ji) * gu(ji,2,1) ) )**2 +       &
                    ( pt5 * ( gu(ji,1,3) + damp(ji) * gu(ji,3,1) ) )**2 +       &
                    ( pt5 * ( damp(ji) * gu(ji,2,1) + gu(ji,1,2) ) )**2 +       &
                    ( pt5 * ( gu(ji,2,2) + gu(ji,2,2) ) )**2 +                  &
                    ( pt5 * ( gu(ji,2,3) + gu(ji,2,3) ) )**2 +                  &
                    ( pt5 * ( damp(ji) * gu(ji,3,1) + gu(ji,1,3) ) )**2 +       &
                    ( pt5 * ( gu(ji,3,2) + gu(ji,2,3) ) )**2 +                  &
                    ( pt5 * ( gu(ji,3,3) + gu(ji,3,3) ) )**2 )

        D(ji,5,2) = -Ma**2*rhoinv(ji)*( bn1*(c**2*grhon-gpn) -                  &
                     pt5*gamma1*bn1*(gpn+rho(ji)*c*gunn) ) + bs1*gts

        D(ji,5,3) = -Ma**2*rhoinv(ji)*( bn2*(c**2*grhon-gpn) -                  &
                     pt5*gamma1*bn2*(gpn+rho(ji)*c*gunn) ) + bs2*gts

        D(ji,5,4) = gt(ji,3)

        D(ji,5,5) = gamma1 * ( guss + gu(ji,3,3) ) -                            &
                    Ma**2*rhoinv(ji)*( un*(two*c*dc*grhon-grhon*gmsinv) -       &
                    pt5*gamma1*( dc*(gpn+rho(ji)*c*gunn) +                      &
                    (un+c)*(grhon*gmsinv + rho(ji)*dc*gunn) +                   &
                    rk*rho(ji)*gmsinv ) ) -                                     &
                   fact1 * (damp(ji) * g1dcon(ji) * gt(ji,1) +                  &
                            g2dcon(ji) * gt(ji,2) +                             &
                            g3dcon(ji) * gt(ji,3) +                             &
                            dcon(ji) * (damp(ji) * g11v(ji,5) + g22v(ji,5)) ) - &
                   fact2 * dlm(ji) * ( damp(ji) * gu(ji,1,1) +                  &
                                   gu(ji,2,2) +                                 &
                                   gu(ji,3,3) )**2 -                            &
                   two * fact2 * dmu(ji) * (                                    &
                    ( pt5 * damp(ji) * ( gu(ji,1,1) + gu(ji,1,1) ) )**2 +       &
                    ( pt5 * ( gu(ji,1,2) + damp(ji) * gu(ji,2,1) ) )**2 +       &
                    ( pt5 * ( gu(ji,1,3) + damp(ji) * gu(ji,3,1) ) )**2 +       &
                    ( pt5 * ( damp(ji) * gu(ji,2,1) + gu(ji,1,2) ) )**2 +       &
                    ( pt5 * ( gu(ji,2,2) + gu(ji,2,2) ) )**2 +                  &
                    ( pt5 * ( gu(ji,2,3) + gu(ji,2,3) ) )**2 +                  &
                    ( pt5 * ( damp(ji) * gu(ji,3,1) + gu(ji,1,3) ) )**2 +       &
                    ( pt5 * ( gu(ji,3,2) + gu(ji,2,3) ) )**2 +                  &
                    ( pt5 * ( gu(ji,3,3) + gu(ji,3,3) ) )**2 )

        end if          ! inflow/outflow

        end do

        end if          ! top
        
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

        return
        end
