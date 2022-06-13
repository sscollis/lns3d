!=============================================================================!
        subroutine rhs_l(rl, vl, i, j) 
!  
!       Generates the RHS for the compressible 2D-3C, N-S equations 
!
!       This version supports the buffer domain where the streamwise
!       viscous derivatives are smoothly removed from the equations.
!
!       This version is for computing the outflow boundary nodes
!       using the Lele & Poinsot formulation.
!
!       Written: 5-10-95
!       Revised: 9-08-95
!       Revised: 3-17-01
!       Revised: 6-06-22
!=============================================================================!
        use global
        use local2d
        use material
        implicit none

        real :: rl(ndof), vl(ndof)
        real :: fact1, fact2, fact3, fact4
        
        real :: pinf, cinf, Lx, rk, cc, c
        real :: L1, L2, L3, L4, L5
        real :: d1, d2, d3, d4, d5
        real :: e1, e2, e3, e4, e5
        
        integer :: i, j
        
        real :: g1vl(ndof), g2vl(ndof), g11vl(ndof), &
                g12vl(ndof), g22vl(ndof)

        real :: bn1, bn2, bs1, bs2, us, un, gps, gpn, grhos, grhon, &
                gts, gtn, gus, gvs, guss, gusn, guns, gunn, gu3s, gu3n
                
        real :: b, rbeta, ub1, ub2, gub11, gub12, gub21, gub22
        real :: ubn, ubs, gubnn, gubsn, tb, rhob, pb, cb, &
                gtbn, grhobn, gpbn
        
        real :: wrho(nx), pnorm(nx)
!=============================================================================!

!.... compute some extra stuff needed for Lele's BC's

        pinf = one / (gamma * Ma**2)    ! pressure at infinity
        cinf = one / Ma                 ! speed of sound at infinity
        Lx   = 100.0                    ! length of domain in x

        rk = sigma * ( one - Ma**2 ) * cinf / Lx

        fact1 = one / (gamma * Ma**2)
        fact2 = one / Re
        fact3 = gamma / ( Pr * Re )
        fact4 = gamma * gamma1 * Ma**2 / Re

!.... localize

        g1vl(:)  = g1v(:,i,j)
        g2vl(:)  = g2v(:,i,j)
        g11vl(:) = g11v(:,i,j)
        g12vl(:) = g12v(:,i,j)
        g22vl(:) = g22v(:,i,j)

!=============================================================================!
!       L e f t   a n d   R i g h t   B o u n d a r i e s
!=============================================================================!

        if (left.eq.1 .or. right.eq.1) then

!.... don't do on wall

        if ( ( (i.eq.1 .and. left.eq.1) .or. &
               (i.eq.nx .and. right.eq.1) ) .and. j.gt.2 ) then

           call error("rhs_l$","rhs_l not currently updated$")

!.... get the metrics along this boundary

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

!.... compute some derivatives in the boundary normal coordinate system

            gpn  = bn1 * gp(1) + bn2 * gp(2)
            gps  = bs1 * gp(1) + bs2 * gp(2)

            grhon = bn1 * grho(1) + bn2 * grho(2)
            grhos = bs1 * grho(1) + bs2 * grho(2)

            gtn  = bn1 * gt(1) + bn2 * gt(2)
            gts  = bs1 * gt(1) + bs2 * gt(2)

            gus  = bs1 * gu(1,1) + bs2 * gu(1,2)
            gvs  = bs1 * gu(2,1) + bs2 * gu(2,2)

            gunn = bn1*bn1 * gu(1,1) + bn1*bn2 * gu(2,1) + &
                   bn1*bn2 * gu(1,2) + bn2*bn2 * gu(2,2)

            gusn = bs1*bn1 * gu(1,1) + bs2*bn1 * gu(2,1) + &
                   bs1*bn2 * gu(1,2) + bs2*bn2 * gu(2,2)

            guns = bn1*bs1 * gu(1,1) + bn2*bs1 * gu(2,1) + &
                   bn1*bs2 * gu(1,2) + bn2*bs2 * gu(2,2)

            guss = bs1*bs1 * gu(1,1) + bs2*bs1 * gu(2,1) + &
                   bs1*bs2 * gu(1,2) + bs2*bs2 * gu(2,2)

            gu3n = bn1 * gu(3,1) + bn2 * gu(3,2)
            gu3s = bs1 * gu(3,1) + bs2 * gu(3,2)

!.... compute the Characteristic amplitudes from the interior

            L1 = ( un - c ) * ( gpn - rho * c * gunn )
            L2 = un * ( c**2 * grhon - gpn )
            L3 = un * gusn
            L4 = un * gu3n
            L5 = ( un + c ) * ( gpn + rho * c * gunn )

            write(*,*) rho, us, un, c, t, gpn, gps

!.... use the potential solution to compute the incomming characteristic

            if (.true.) then             ! turn off for acoustic wave test
              b = one / pi
              if (Ma.lt.one) then
                rbeta = sqrt(one - Ma**2)
              else
                rbeta = sqrt(Ma**2 - one)
              end if
              ub1 = one + (x(1,j)-b)/pi/rbeta**2/((x(i,j)-b)**2 + &
                            rbeta**2*y(i,j)**2)
              ub2 = rbeta*y(i,j)/pi/((x(i,j)-b)**2+rbeta**2*y(i,j)**2)

              gub11 = one/( pi*rbeta**2*((x(i,j)-b)**2+rbeta**2*y(i,j)**2) ) -&
                      two*(x(i,j)-b)**2/ &
                      ( pi*rbeta**2*((x(i,j)-b)**2+rbeta**2*y(i,j)**2)**2 )

              gub12 = -two*(x(i,j)-b)*rbeta**2*y(i,j)/ &
                      ( pi*rbeta**2*((x(i,j)-b)**2+ &
                      rbeta**2*y(i,j)**2)**2 )

              gub21 = -rbeta*y(i,j)*two*(x(i,j)-b)/( pi*((x(i,j)-b)**2+ &
                        rbeta**2*y(i,j)**2)**2 )

              gub22 = rbeta/(pi*((x(i,j)-b)**2+rbeta**2*y(i,j)**2)) - &
                      two*rbeta**3*y(i,j)**2/(pi*((x(i,j)-b)**2+ &
                      rbeta**2*y(i,j)**2)**2)

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

!.... nonreflecting outflow boundary condition
              L1 = (ubn - cb)*(gpbn - rhob*cb*gubnn) + rk * ( p - pinf )
            else
              L1 = rk * ( p - pinf )
            end if

!.... compute the boundary normal derivative terms

            d1 = one / c**2 * ( L2 + pt5 * ( L5 + L1 ) )
            d2 = pt5 * ( L5 + L1 )
            d3 = pt5 / ( rho * c ) * ( L5 - L1 )
            d4 = L3
            d5 = L4

!.... debug

!           write (69,"(8(1pe13.6,1x))") y(i,j), &
!             g1mu*two*S(2,1) + mu*(g11vl(3)+g12vl(2)), &
!             g1con*gt(1) + con*g11vl(5)

!.... continuity

            rl(1) = d1 + grhos * us + grho(3) * u3 + &
                       rho * ( guss + gu(3,3) )

!.... momentum

            rl(2) = bn1 * d3 + bs1 * d4 + us * gus + &
                       bs1 * rhoinv * gps - &
                       fact2 * ( damp(i,j) * (g1lm * divu + &
                       lm * g1divu) + &
                       two * (damp(i,j) * g1mu * S(1,1) + &
                              g2mu * pt5 * (gu(1,2) + &
                              damp(i,j) * gu(2,1)) ) + &
                       mu * (damp(i,j) * two * g11vl(2) +  &
                              g22vl(2) + damp(i,j) * g12vl(3) ) ) / rho

!           rl(3) = bn2 * d3 + bs2 * d4 + us * gvs + bs2 * rhoinv * gps -        &
!                 fact2 * ( g2lm * (damp(i,j) * gu(1,1) + gu(2,2) +            &
!                                   gu(3,3)) +                                       &
!                           lm * (damp(i,j) * g12vl(2) + g22vl(3)) +             &
!                           two * (damp(i,j) * g1mu * S(2,1) +                    &
!                                  g2mu * S(2,2)) +                              &
!                           mu * ( damp(i,j) * (g11vl(3) + g12vl(2)) +           &
!                                  two * g22vl(3) ) ) / rho

!.... set tau21,1 = zero

            rl(3) = bn2 * d3 + bs2 * d4 + us * gvs + bs2 * rhoinv * gps -        &
                  fact2 * ( g2lm * (damp(i,j) * gu(1,1) + gu(2,2) +            &
                                        gu(3,3)) +                                   &
                            lm * (damp(i,j) * g12vl(2) + g22vl(3)) +             &
                            two * ( g2mu * S(2,2) ) +                            &
                            mu * (two * g22vl(3) ) ) / rho

            rl(4) = d5 + us * gu3s + gp(3) * rhoinv -                         &
                  fact2 * ( two * (damp(i,j) * g1mu * S(3,1) +                    &
                                   g2mu * S(3,2)) +                              &
                            mu * ( damp(i,j) * g11vl(4) +                          &
                                   g22vl(4) ) ) / rho

!.... temperature

!           rl(5) = ( gamma * Ma**2 * d2 - t * d1 ) * rhoinv +               &
!                 us * gts + u3 * gt(3) +                                        &
!                 gamma1 * t * ( guss + gu(3,3) )  -                             &
!                 fact3 / rho * ( damp(i,j) * g1con * gt(1) +                 &
!                               g2con * gt(2) +                                  &
!                               con * (damp(i,j) * g11vl(5) +                      &
!                                      g22vl(5)) ) -                                  &
!                 fact4 / rho * ( lm * (damp(i,j) * gu(1,1) + gu(2,2) +    &
!                                     gu(3,3))**2 +                                  &
!                   two * mu * (                                                    &
!                   ( pt5 * (two * damp(i,j) * gu(1,1)) )**2 +                        &
!                   ( pt5 * (gu(1,2) + damp(i,j) * gu(2,1)) )**2 +                 &
!                   ( pt5 * (gu(1,3) + damp(i,j) * gu(3,1)) )**2 +                 &
!                   ( pt5 * (gu(1,2) + damp(i,j) * gu(2,1)) )**2 +                 &
!                   ( pt5 * (two * gu(2,2)) )**2 +                                   &
!                   ( pt5 * (gu(2,3) + gu(3,2)) )**2 +                            &
!                   ( pt5 * (gu(1,3) + damp(i,j) * gu(3,1)) )**2 +                 &
!                   ( pt5 * (gu(3,2) + gu(2,3)) )**2 +                            &
!                   ( pt5 * (two * gu(3,3)) )**2 ) )

!.... set q1,1 = 0

            rl(5) = ( gamma * Ma**2 * d2 - t * d1 ) * rhoinv +               &
                  us * gts + u3 * gt(3) +                                        &
                  gamma1 * t * ( guss + gu(3,3) )  -                             &
                  fact3 / rho * (g2con * gt(2) + con * g22vl(5) ) -    &
                  fact4 / rho * ( lm * (damp(i,j) * gu(1,1) + gu(2,2) +    &
                                      gu(3,3))**2 +                                  &
                    two * mu * (                                                    &
                    ( pt5 * (two * damp(i,j) * gu(1,1)) )**2 +                        &
                    ( pt5 * (gu(1,2) + damp(i,j) * gu(2,1)) )**2 +                 &
                    ( pt5 * (gu(1,3) + damp(i,j) * gu(3,1)) )**2 +                 &
                    ( pt5 * (gu(1,2) + damp(i,j) * gu(2,1)) )**2 +                 &
                    ( pt5 * (two * gu(2,2)) )**2 +                                   &
                    ( pt5 * (gu(2,3) + gu(3,2)) )**2 +                            &
                    ( pt5 * (gu(1,3) + damp(i,j) * gu(3,1)) )**2 +                 &
                    ( pt5 * (gu(3,2) + gu(2,3)) )**2 +                            &
                    ( pt5 * (two * gu(3,3)) )**2 ) )
          end if
        end if

!=============================================================================!
!       T o p
!=============================================================================!
!  Key to the generalized Poinsot & Lele boundary conditions.
!  Boundary conditions are required for all incomming 
!  characteristic amplitudes
!
!  n is the unit outward pointing normal vector 
!  (it points out from the boundary)
!
!  Subsonic inflow:     un < 0
!
!       L1 -> in
!       L2 -> in
!       L3 -> in
!       L4 -> in
!       L5 -> out
!
!  Subsonic outflow:    un > 0
!
!       L1 -> in
!       L2 -> in
!       L3 -> in
!       L4 -> in
!       L5 -> out
!=============================================================================!

        if (top.eq.1 .and. j.eq.ny) then

           call error("rhs_l$","rhs_l not currently updated$")

!.... get the metrics along this boundary

            bn1 = n1(i,j) / sqrt( n1(i,j)**2 + n2(i,j)**2 )
            bn2 = n2(i,j) / sqrt( n1(i,j)**2 + n2(i,j)**2 )

            bs1 = -bn2
            bs2 =  bn1

!.... compute Us and Un

            un = bn1 * vl(2) + bn2 * vl(3)
            us = bs1 * vl(2) + bs2 * vl(3)

            c = sqrt( t ) / Ma

!.... compute some derivatives in the boundary normal coordinate system

            gpn  = bn1 * gp(1) + bn2 * gp(2)
            gps  = bs1 * gp(1) + bs2 * gp(2)

            grhon = bn1 * grho(1) + bn2 * grho(2)
            grhos = bs1 * grho(1) + bs2 * grho(2)

            gtn  = bn1 * gt(1) + bn2 * gt(2)
            gts  = bs1 * gt(1) + bs2 * gt(2)

            gus  = bs1 * gu(1,1) + bs2 * gu(1,2)
            gvs  = bs1 * gu(2,1) + bs2 * gu(2,2)

            gunn = bn1*bn1 * gu(1,1) + bn1*bn2 * gu(2,1) + &
                   bn1*bn2 * gu(1,2) + bn2*bn2 * gu(2,2)

            gusn = bs1*bn1 * gu(1,1) + bs2*bn1 * gu(2,1) + &
                   bs1*bn2 * gu(1,2) + bs2*bn2 * gu(2,2)

            guns = bn1*bs1 * gu(1,1) + bn2*bs1 * gu(2,1) + &
                   bn1*bs2 * gu(1,2) + bn2*bs2 * gu(2,2)

            guss = bs1*bs1 * gu(1,1) + bs2*bs1 * gu(2,1) + &
                   bs1*bs2 * gu(1,2) + bs2*bs2 * gu(2,2)

            gu3n = bn1 * gu(3,1) + bn2 * gu(3,2)
            gu3s = bs1 * gu(3,1) + bs2 * gu(3,2)

!.... compute the Characteristic amplitudes from the interior

            L1 = ( un - c ) * ( gpn - rho * c * gunn )
            L2 = un * ( c**2 * grhon - gpn )
            L3 = un * gusn
            L4 = un * gu3n
            L5 = ( un + c ) * ( gpn + rho * c * gunn )

!.... use the potential solution

            b = one / pi
            if (Ma.lt.one) then
              rbeta = sqrt(one - Ma**2)
            else
              rbeta = sqrt(Ma**2 - one)
            end if
            ub1 = one + (x(1,j)-b)/pi/rbeta**2/((x(i,j)-b)**2 + &
                         rbeta**2*y(i,j)**2)
            ub2 = rbeta*y(i,j)/pi/((x(i,j)-b)**2+rbeta**2*y(i,j)**2)

            gub11 = one/( pi*rbeta**2*((x(i,j)-b)**2+rbeta**2*y(i,j)**2) ) - &
                    two*(x(i,j)-b)**2/( pi*rbeta**2*((x(i,j)-b)**2+rbeta**2*y(i,j)**2)**2 )

            gub12 = -two*(x(i,j)-b)*rbeta**2*y(i,j)/( pi*rbeta**2*((x(i,j)-b)**2+ &
                     rbeta**2*y(i,j)**2)**2 )

            gub21 = -rbeta*y(i,j)*two*(x(i,j)-b)/( pi*((x(i,j)-b)**2+ &
                     rbeta**2*y(i,j)**2)**2 )

            gub22 = rbeta/(pi*((x(i,j)-b)**2+rbeta**2*y(i,j)**2)) - &
                    two*rbeta**3*y(i,j)**2/(pi*((x(i,j)-b)**2+ &
                    rbeta**2*y(i,j)**2)**2)

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

            if ( un .le. zero ) then
              L1 = (ubn - cb)*(gpbn - rhob*cb*gubnn)
              L2 = ubn*(cb**2*grhobn - gpbn)
              L3 = ubn*gubsn
              L4 = zero
!             L1 = L5
!             L2 = pt5*gamma1*(L5 + L1)
!             L3 = zero
!             L4 = zero
            else
              L1 = rk * ( p - pinf ) + (ubn - cb)*(gpbn - rhob*cb*gubnn)    
            end if

!.... compute the boundary normal derivative terms

            d1 = one / c**2 * ( L2 + pt5 * ( L5 + L1 ) )
            d2 = pt5 * ( L5 + L1 )
            d3 = pt5 / ( rho * c ) * ( L5 - L1 )
            d4 = L3
            d5 = L4

!.... continuity

            rl(1) = d1 + grhos * us + grho(3) * u3 + &
                       rho * ( guss + gu(3,3) )

!.... momentum

            rl(2) = bn1 * d3 + bs1 * d4 + us * gus + bs1 * rhoinv * gps -        &
                  fact2 * ( damp(i,j) * (g1lm * divu + lm * g1divu) +    &
                            two * (damp(i,j) * g1mu * S(1,1) +                    &
                                   g2mu * pt5 * (gu(1,2) +                       &
                                   damp(i,j) * gu(2,1)) ) +                           &
                            mu * (damp(i,j) * two * g11vl(2) +                     &
                            g22vl(2) + damp(i,j) * g12vl(3) ) ) / rho

            rl(3) = bn2 * d3 + bs2 * d4 + us * gvs + bs2 * rhoinv * gps -        &
                  fact2 * ( g2lm * (damp(i,j) * gu(1,1) + gu(2,2) +            &
                                    gu(3,3)) +                                       &
                            lm * (damp(i,j) * g12vl(2) + g22vl(3)) +             &
                            two * (damp(i,j) * g1mu * S(2,1) +                    &
                                   g2mu * S(2,2)) +                              &
                            mu * ( damp(i,j) * (g11vl(3) + g12vl(2)) +           &
                                   two * g22vl(3) ) ) / rho

            rl(4) = d5 + us * gu3s + gp(3) * rhoinv -                         &
                  fact2 * ( two * (damp(i,j) * g1mu * S(3,1) +                    &
                                   g2mu * S(3,2)) +                              &
                            mu * ( damp(i,j) * g11vl(4) +                          &
                                   g22vl(4) ) ) / rho

!.... temperature

            rl(5) = ( gamma * Ma**2 * d2 - t * d1 ) * rhoinv +               &
                  us * gts + u3 * gt(3) +                                        &
                  gamma1 * t * ( guss + gu(3,3) )  -                             &
                  fact3 / rho * ( damp(i,j) * g1con * gt(1) +                 &
                                g2con * gt(2) +                                  &
                                con * (damp(i,j) * g11vl(5) +                      &
                                       g22vl(5)) ) -                                  &
                  fact4 / rho * ( lm * (damp(i,j) * gu(1,1) + gu(2,2) +    &
                                      gu(3,3))**2 +                                  &
                    two * mu * (                                                    &
                    ( pt5 * (two * damp(i,j) * gu(1,1)) )**2 +                        &
                    ( pt5 * (gu(1,2) + damp(i,j) * gu(2,1)) )**2 +                 &
                    ( pt5 * (gu(1,3) + damp(i,j) * gu(3,1)) )**2 +                 &
                    ( pt5 * (gu(1,2) + damp(i,j) * gu(2,1)) )**2 +                 &
                    ( pt5 * (two * gu(2,2)) )**2 +                                   &
                    ( pt5 * (gu(2,3) + gu(3,2)) )**2 +                            &
                    ( pt5 * (gu(1,3) + damp(i,j) * gu(3,1)) )**2 +                 &
                    ( pt5 * (gu(3,2) + gu(2,3)) )**2 +                            &
                    ( pt5 * (two * gu(3,3)) )**2 ) )

!           if ( un .le. zero ) then
!              rl(2) = zero
!              rl(3) = zero
!              rl(4) = zero
!              rl(5) = zero
!           end if

        end if

!=============================================================================!
!       W a l l
!=============================================================================!
!  Key to the generalized Poinsot & Lele boundary conditions.
!  Boundary conditions are required for all incomming 
!  characteristic amplitudes
!
!  n is the unit normal vector pointing INTO the fluid
!
!  Noslip wall:         un = 0
!
!       L1 -> out
!       L2 -> zero
!       L3 -> zero
!       L4 -> zero
!       L5 -> in
!
!  BC:  L5 = L1 
!=============================================================================!

        if (wall.eq.4 .and. j.eq.1) then
        
           call error("rhs_l$","rhs_l not currently updated$")

!.... compute the wall normal pressure gradient using the momentum equations

          call wallbc(vl,wrho,pnorm)

!.... get the metrics along this boundary

            bn1 = n1(i,j) / sqrt( n1(i,j)**2 + n2(i,j)**2 )
            bn2 = n2(i,j) / sqrt( n1(i,j)**2 + n2(i,j)**2 )

!.... compute Us and Un
        
            un =  bn1 * vl(2) + bn2 * vl(3)
            us = -bn2 * vl(2) + bn1 * vl(3)

            c = sqrt( t ) / Ma

!.... compute some derivatives in the boundary normal coordinate system

            gpn =  bn1 * gp(1) + bn2 * gp(2)
!           gpn =  sqrt( n1**2 + n2**2 ) * pnorm(i)
            
            gps = -bn2 * gp(1) + bn1 * gp(2)
        
            grhon =  bn1 * grho(1) + bn2 * grho(2)
            grhos = -bn2 * grho(1) + bn1 * grho(2)

            gtn =  bn1 * gt(1) + bn2 * gt(2)
            gts = -bn2 * gt(1) + bn1 * gt(2)

            gus  =  -bn2 * gu(1,1) + bn1 * gu(1,2)
            gvs  =  -bn2 * gu(2,1) + bn1 * gu(2,2)

            gunn =  bn1*bn1 * gu(1,1) + bn1*bn2 * gu(2,1) + &
                    bn1*bn2 * gu(1,2) + bn2*bn2 * gu(2,2)
            gusn = -bn1*bn2 * gu(1,1) + bn1*bn1 * gu(2,1) - &
                    bn2*bn2 * gu(1,2) + bn1*bn2 * gu(2,2)
            guns = -bn2*bn1 * gu(1,1) - bn2*bn2 * gu(2,1) + &
                    bn1*bn1 * gu(1,2) + bn1*bn2 * gu(2,2)
            guss =  bn2*bn2 * gu(1,1) - bn1*bn2 * gu(2,1) - &
                    bn1*bn2 * gu(1,2) + bn1*bn1 * gu(2,2)
            
            gu3n =  bn1 * gu(3,1) + bn2 * gu(3,2)
            gu3s = -bn2 * gu(3,1) + bn1 * gu(3,2)
            
!.... compute the Characteristic amplitudes
        
            L1 = ( un - c ) * ( gpn - rho * c * gunn )
            L2 = un * ( c**2 * grhon - gpn )
            L3 = un * gusn
            L4 = un * gu3n
            L5 = ( un + c ) * ( gpn + rho * c * gunn )
        
!.... adiabatic wall boundary conditions
        
            L2 = zero
            L3 = zero
            L4 = zero
            L5 = L1
                    
!.... compute the boundary normal derivative terms

            d1 = one / c**2 * ( L2 + pt5 * ( L5 + L1 ) )
            d2 = pt5 * ( L5 + L1 )
            d3 = pt5 / ( rho * c ) * ( L5 - L1 )
            d4 = L3
            d5 = L4

!.... continuity

            rl(1) = d1 + grhos * us + grho(3) * u3 + &
                       rho * ( guss + gu(3,3) )
    
!.... momentum

            rl(2) = bn1 * d3 - bn2 * d4 + us * gus - bn2 * rhoinv * gps -        &
                  fact2 * ( damp(i,j) * (g1lm * divu + lm * g1divu) +    &
                            two * (damp(i,j) * g1mu * S(1,1) +                    &
                                   g2mu * pt5 * (gu(1,2) +                       &
                                   damp(i,j) * gu(2,1)) ) +                           &
                            mu * (damp(i,j) * two * g11vl(2) +                     &
                            g22vl(2) + damp(i,j) * g12vl(3) ) ) / rho
    
            rl(3) = bn2 * d3 + bn1 * d4 + us * gvs + bn1 * rhoinv * gps -        &
                  fact2 * ( g2lm * (damp(i,j) * gu(1,1) + gu(2,2) +            &
                                    gu(3,3)) +                                       &
                            lm * (damp(i,j) * g12vl(2) + g22vl(3)) +             &
                            two * (damp(i,j) * g1mu * S(2,1) +                    &
                                   g2mu * S(2,2)) +                              &
                            mu * ( damp(i,j) * (g11vl(3) + g12vl(2)) +           &
                                   two * g22vl(3) ) ) / rho
    
            rl(4) = d5 + us * gu3s + gp(3) / rho -                            &
                  fact2 * ( two * (damp(i,j) * g1mu * S(3,1) +                    &
                                   g2mu * S(3,2)) +                              &
                            mu * ( damp(i,j) * g11vl(4) +                          &
                                   g22vl(4) ) ) / rho

!.... temperature

            rl(5) = ( gamma * Ma**2 * d2 - t * d1 ) * rhoinv +               &
                  us * gts + u3 * gt(3) +                                        &
                  gamma1 * t * ( guss + gu(3,3) )  -                             &
                  fact3 / rho * ( damp(i,j) * g1con * gt(1) +                 &
                                g2con * gt(2) +                                  &
                                con * (damp(i,j) * g11vl(5) +                      &
                                       g22vl(5)) ) -                                  &
                  fact4 / rho * ( lm * (damp(i,j) * gu(1,1) + gu(2,2) +    &
                                      gu(3,3))**2 +                                  &
                    two * mu * (                                                    &
                    ( pt5 * (two * damp(i,j) * gu(1,1)) )**2 +                        &
                    ( pt5 * (gu(1,2) + damp(i,j) * gu(2,1)) )**2 +                 &
                    ( pt5 * (gu(1,3) + damp(i,j) * gu(3,1)) )**2 +                 &
                    ( pt5 * (gu(1,2) + damp(i,j) * gu(2,1)) )**2 +                 &
                    ( pt5 * (two * gu(2,2)) )**2 +                                   &
                    ( pt5 * (gu(2,3) + gu(3,2)) )**2 +                            &
                    ( pt5 * (gu(1,3) + damp(i,j) * gu(3,1)) )**2 +                 &
                    ( pt5 * (gu(3,2) + gu(2,3)) )**2 +                            &
                    ( pt5 * (two * gu(3,3)) )**2 ) )
        end if

        return
        end subroutine rhs_l
