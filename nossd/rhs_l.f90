!=============================================================================!
        subroutine rhs_l(rl, vl) 
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
!       Revised: 9-8-95
!
!=============================================================================!
        use global
        use local
        implicit none

        real :: rl(ny*nx,ndof), vl(ny*nx,ndof)
        real :: fact1, fact2, fact3, fact4
        
        real :: pinf, cinf, Lx, rk, c
        real :: L1, L2, L3, L4, L5
        real :: d1, d2, d3, d4, d5
        real :: e1, e2, e3, e4, e5
        
        integer :: i, j, ji
        
        real :: bn1, bn2, bs1, bs2, us, un, gps, gpn, grhos, grhon, gts, gtn, &
                gus, gvs, guss, gusn, guns, gunn, gu3s, gu3n
                
        real :: b, rbeta, ub1, ub2, gub11, gub12, gub21, gub22
        real :: ubn, ubs, gubnn, gubsn, tb, rhob, pb, cb, gtbn, grhobn, gpbn
        
        real :: wrho(nx), pnorm(nx)
!=============================================================================!

!.... compute some extra stuff needed for Lele's BC's

        pinf = one / (gamma * Ma**2)    ! pressure at infinity
        cinf = one / Ma                 ! speed of sound at infinity
        Lx   = 100.0                    ! length of domain in x
        
        rk = sigma * ( one - Ma**2 ) * cinf / Lx
        
        fact1 = one / (gamma * Ma**2)
        fact2 = one / Re
        fact3 = gamma / ( Re * Pr )
        fact4 = gamma * gamma1 * Ma**2 / Re

!=============================================================================!
!       L e f t   a n d   R i g h t   B o u n d a r i e s
!=============================================================================!
        if (left.eq.1 .or. right.eq.1) then

        do i = 1, nx, nx-1
          if ( (i.eq.1 .and. left.eq.1) .or. &
               (i.eq.nx .and. right.eq.1) ) then          
          do j = 2, ny                  ! don't do on the wall
            ji = j + (i-1) * ny
    
!.... get the metrics along this boundary

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

!.... compute some derivatives in the boundary normal coordinate system

            gpn  = bn1 * gp(ji,1) + bn2 * gp(ji,2)
            gps  = bs1 * gp(ji,1) + bs2 * gp(ji,2)
        
            grhon = bn1 * grho(ji,1) + bn2 * grho(ji,2)
            grhos = bs1 * grho(ji,1) + bs2 * grho(ji,2)

            gtn  = bn1 * gt(ji,1) + bn2 * gt(ji,2)
            gts  = bs1 * gt(ji,1) + bs2 * gt(ji,2)

            gus  = bs1 * gu(ji,1,1) + bs2 * gu(ji,1,2)
            gvs  = bs1 * gu(ji,2,1) + bs2 * gu(ji,2,2)

            gunn = bn1*bn1 * gu(ji,1,1) + bn1*bn2 * gu(ji,2,1) + &
                   bn1*bn2 * gu(ji,1,2) + bn2*bn2 * gu(ji,2,2)

            gusn = bs1*bn1 * gu(ji,1,1) + bs2*bn1 * gu(ji,2,1) + &
                   bs1*bn2 * gu(ji,1,2) + bs2*bn2 * gu(ji,2,2)

            guns = bn1*bs1 * gu(ji,1,1) + bn2*bs1 * gu(ji,2,1) + &
                   bn1*bs2 * gu(ji,1,2) + bn2*bs2 * gu(ji,2,2)

            guss = bs1*bs1 * gu(ji,1,1) + bs2*bs1 * gu(ji,2,1) + &
                   bs1*bs2 * gu(ji,1,2) + bs2*bs2 * gu(ji,2,2)
            
            gu3n = bn1 * gu(ji,3,1) + bn2 * gu(ji,3,2)
            gu3s = bs1 * gu(ji,3,1) + bs2 * gu(ji,3,2)
            
!.... compute the Characteristic amplitudes from the interior
        
            L1 = ( un - c ) * ( gpn - rho(ji) * c * gunn )
            L2 = un * ( c**2 * grhon - gpn )
            L3 = un * gusn
            L4 = un * gu3n
            L5 = ( un + c ) * ( gpn + rho(ji) * c * gunn )

!.... use the potential solution to compute the incomming characteristic

            if (.true.) then             ! turn off for acoustic wave test
              b = one / pi
              if (Ma.lt.one) then
                rbeta = sqrt(one - Ma**2)
              else
                rbeta = sqrt(Ma**2 - one)
              end if
              ub1 = one + (x(j,1)-b)/pi/rbeta**2/((x(j,i)-b)**2 + &
                            rbeta**2*y(j,i)**2)
              ub2 = rbeta*y(j,i)/pi/((x(j,i)-b)**2+rbeta**2*y(j,i)**2)
  
              gub11 = one/( pi*rbeta**2*((x(j,i)-b)**2+rbeta**2*y(j,i)**2) ) -&
                      two*(x(j,i)-b)**2/ &
                      ( pi*rbeta**2*((x(j,i)-b)**2+rbeta**2*y(j,i)**2)**2 )
  
              gub12 = -two*(x(j,i)-b)*rbeta**2*y(j,i)/ &
                      ( pi*rbeta**2*((x(j,i)-b)**2+ &
                      rbeta**2*y(j,i)**2)**2 )
  
              gub21 = -rbeta*y(j,i)*two*(x(j,i)-b)/( pi*((x(j,i)-b)**2+ &
                        rbeta**2*y(j,i)**2)**2 )
  
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

!.... nonreflecting outflow boundary condition
        
              L1 = (ubn - cb)*(gpbn - rhob*cb*gubnn) + rk * ( p(ji) - pinf )
            else
              L1 = rk * ( p(ji) - pinf )
            end if
            
!.... compute the boundary normal derivative terms

            d1 = one / c**2 * ( L2 + pt5 * ( L5 + L1 ) )
            d2 = pt5 * ( L5 + L1 )
            d3 = pt5 / ( rho(ji) * c ) * ( L5 - L1 )
            d4 = L3
            d5 = L4

!.... debug
 
!           write (69,"(8(1pe13.6,1x))") y(j,i), &
!             g1mu(ji)*two*S(ji,2,1) + mu(ji)*(g11v(ji,3)+g12v(ji,2)), &
!             g1con(ji)*gt(ji,1) + con(ji)*g11v(ji,5)

!.... continuity

            rl(ji,1) = d1 + grhos * us + grho(ji,3) * u3(ji) + &
                       rho(ji) * ( guss + gu(ji,3,3) )
    
!.... momentum

            rl(ji,2) = bn1 * d3 + bs1 * d4 + us * gus + &
                       bs1 * rhoinv(ji) * gps - &
                       fact2 * ( damp(ji) * (g1lm(ji) * divu(ji) + &
                       lm(ji) * g1divu(ji)) + &
                       two * (damp(ji) * g1mu(ji) * S(ji,1,1) + &
                              g2mu(ji) * pt5 * (gu(ji,1,2) + &
                              damp(ji) * gu(ji,2,1)) ) + &
                       mu(ji) * (damp(ji) * two * g11v(ji,2) +  &
                                g22v(ji,2) + damp(ji) * g12v(ji,3) ) ) / rho(ji)
    
!           rl(ji,3) = bn2 * d3 + bs2 * d4 + us * gvs + bs2 * rhoinv(ji) * gps -        &
!                 fact2 * ( g2lm(ji) * (damp(ji) * gu(ji,1,1) + gu(ji,2,2) +            &
!                                   gu(ji,3,3)) +                                       &
!                           lm(ji) * (damp(ji) * g12v(ji,2) + g22v(ji,3)) +             &
!                           two * (damp(ji) * g1mu(ji) * S(ji,2,1) +                    &
!                                  g2mu(ji) * S(ji,2,2)) +                              &
!                           mu(ji) * ( damp(ji) * (g11v(ji,3) + g12v(ji,2)) +           &
!                                  two * g22v(ji,3) ) ) / rho(ji)

!.... set tau21,1 = zero

            rl(ji,3) = bn2 * d3 + bs2 * d4 + us * gvs + bs2 * rhoinv(ji) * gps -        &
                  fact2 * ( g2lm(ji) * (damp(ji) * gu(ji,1,1) + gu(ji,2,2) +            &
                                        gu(ji,3,3)) +                                   &
                            lm(ji) * (damp(ji) * g12v(ji,2) + g22v(ji,3)) +             &
                            two * ( g2mu(ji) * S(ji,2,2) ) +                            &
                            mu(ji) * (two * g22v(ji,3) ) ) / rho(ji)
        
            rl(ji,4) = d5 + us * gu3s + gp(ji,3) * rhoinv(ji) -                         &
                  fact2 * ( two * (damp(ji) * g1mu(ji) * S(ji,3,1) +                    &
                                   g2mu(ji) * S(ji,3,2)) +                              &
                            mu(ji) * ( damp(ji) * g11v(ji,4) +                          &
                                   g22v(ji,4) ) ) / rho(ji)

!.... temperature

!           rl(ji,5) = ( gamma * Ma**2 * d2 - t(ji) * d1 ) * rhoinv(ji) +               &
!                 us * gts + u3(ji) * gt(ji,3) +                                        &
!                 gamma1 * t(ji) * ( guss + gu(ji,3,3) )  -                             &
!                 fact3 / rho(ji) * ( damp(ji) * g1con(ji) * gt(ji,1) +                 &
!                               g2con(ji) * gt(ji,2) +                                  &
!                               con(ji) * (damp(ji) * g11v(ji,5) +                      &
!                                      g22v(ji,5)) ) -                                  &
!                 fact4 / rho(ji) * ( lm(ji) * (damp(ji) * gu(ji,1,1) + gu(ji,2,2) +    &
!                                     gu(ji,3,3))**2 +                                  &
!                   two * mu(ji) * (                                                    &
!                   ( pt5 * (two * damp(ji) * gu(ji,1,1)) )**2 +                        &
!                   ( pt5 * (gu(ji,1,2) + damp(ji) * gu(ji,2,1)) )**2 +                 &
!                   ( pt5 * (gu(ji,1,3) + damp(ji) * gu(ji,3,1)) )**2 +                 &
!                   ( pt5 * (gu(ji,1,2) + damp(ji) * gu(ji,2,1)) )**2 +                 &
!                   ( pt5 * (two * gu(ji,2,2)) )**2 +                                   &
!                   ( pt5 * (gu(ji,2,3) + gu(ji,3,2)) )**2 +                            &
!                   ( pt5 * (gu(ji,1,3) + damp(ji) * gu(ji,3,1)) )**2 +                 &
!                   ( pt5 * (gu(ji,3,2) + gu(ji,2,3)) )**2 +                            &
!                   ( pt5 * (two * gu(ji,3,3)) )**2 ) )

!.... set q1,1 = 0

            rl(ji,5) = ( gamma * Ma**2 * d2 - t(ji) * d1 ) * rhoinv(ji) +               &
                  us * gts + u3(ji) * gt(ji,3) +                                        &
                  gamma1 * t(ji) * ( guss + gu(ji,3,3) )  -                             &
                  fact3 / rho(ji) * (g2con(ji) * gt(ji,2) + con(ji) * g22v(ji,5) ) -    &
                  fact4 / rho(ji) * ( lm(ji) * (damp(ji) * gu(ji,1,1) + gu(ji,2,2) +    &
                                      gu(ji,3,3))**2 +                                  &
                    two * mu(ji) * (                                                    &
                    ( pt5 * (two * damp(ji) * gu(ji,1,1)) )**2 +                        &
                    ( pt5 * (gu(ji,1,2) + damp(ji) * gu(ji,2,1)) )**2 +                 &
                    ( pt5 * (gu(ji,1,3) + damp(ji) * gu(ji,3,1)) )**2 +                 &
                    ( pt5 * (gu(ji,1,2) + damp(ji) * gu(ji,2,1)) )**2 +                 &
                    ( pt5 * (two * gu(ji,2,2)) )**2 +                                   &
                    ( pt5 * (gu(ji,2,3) + gu(ji,3,2)) )**2 +                            &
                    ( pt5 * (gu(ji,1,3) + damp(ji) * gu(ji,3,1)) )**2 +                 &
                    ( pt5 * (gu(ji,3,2) + gu(ji,2,3)) )**2 +                            &
                    ( pt5 * (two * gu(ji,3,3)) )**2 ) )
          end do
          end if
        end do
        
        if (.false.) then               ! old way
        
        do i = 1, nx, nx-1
          do j = 2, ny                  ! don't do on the wall
            ji = j + (i-1) * ny
    
            c = sqrt( t(ji) ) / Ma

!.... compute the Characteristic amplitudes
        
!           L1 = ( u1(ji) - c ) * ( gp(ji,1) - rho(ji) * c * gu(ji,1,1) )
            L2 = u1(ji) * ( c**2 * grho(ji,1) - gp(ji,1) )
!           L3 = u1(ji) * gu(ji,2,1)
!           L4 = u1(ji) * gu(ji,3,1)
            L5 = ( u1(ji) + c ) * ( gp(ji,1) + rho(ji) * c * gu(ji,1,1) )
        
!.... nonreflecting outflow boundary condition
        
            L1 = rk * ( p(ji) - pinf )  ! nonreflecting condition

!.... compute the streamwise derivative terms

            d1 = one / c**2 * ( L2 + pt5 * ( L5 + L1 ) )
            d2 = pt5 * ( L5 + L1 )
            d3 = pt5 / ( rho(ji) * c ) * ( L5 - L1 )
!           d4 = L3
!           d5 = L4

!.... continuity

            rl(ji,1) = d1 + grho(ji,2) * u2(ji) + grho(ji,3) * u3(ji) +         &
                       rho(ji) * ( gu(ji,2,2) + gu(ji,3,3) )
    
!.... momentum

            rl(ji,2) = d3 + u2(ji) * gu(ji,1,2)  -                              &
                  fact2 * ( zero * damp(ji) * (g1lm(ji) * divu(ji) +            &
                            lm(ji) * g1divu(ji)) +                              &
                            two * (zero * damp(ji) * g1mu(ji) * S(ji,1,1) +     &
                                   g2mu(ji) * pt5 * (gu(ji,1,2) +               &
                                   damp(ji) * gu(ji,2,1)) ) +                   &
                            mu(ji) * (zero * damp(ji) * two * g11v(ji,2) +      &
                            g22v(ji,2) + damp(ji) * g12v(ji,3) ) ) / rho(ji)
    
            rl(ji,3) = u1(ji) * gu(ji,2,1) + u2(ji) * gu(ji,2,2) +              &
                       gp(ji,2) / rho(ji) -                                     &
                  fact2 * ( g2lm(ji) * (damp(ji) * gu(ji,1,1) + gu(ji,2,2) +    &
                                    gu(ji,3,3)) +                               &
                            lm(ji) * (damp(ji) * g12v(ji,2) + g22v(ji,3)) +     &
                            two * (damp(ji) * zero * g1mu(ji) * S(ji,2,1) +     &
                                   g2mu(ji) * S(ji,2,2)) +                      &
                            mu(ji) * ( damp(ji) * (zero * g11v(ji,3) +          &
                                       zero * g12v(ji,2)) +                     &
                                   two * g22v(ji,3) ) ) / rho(ji)
    
            rl(ji,4) = u1(ji) * gu(ji,3,1) + u2(ji) * gu(ji,3,2) +              &
                       gp(ji,3) / rho(ji) -                                     &
                  fact2 * ( two * (damp(ji) * g1mu(ji) * S(ji,3,1) +            &
                                   g2mu(ji) * S(ji,3,2)) +                      &
                            mu(ji) * ( damp(ji) * g11v(ji,4) +                  &
                                   g22v(ji,4) ) ) / rho(ji)

!.... temperature


            rl(ji,5) = ( gamma * Ma**2 * d2 - t(ji) * d1 ) / rho(ji) +          &
                  u2(ji) * gt(ji,2) + u3(ji) * gt(ji,3) +                       &
                  gamma1 * t(ji) * ( gu(ji,2,2) + gu(ji,3,3) )  -               &
                  fact3 / rho(ji) * ( damp(ji) * zero * g1con(ji) * gt(ji,1) +  &
                                g2con(ji) * gt(ji,2) +                          &
                                con(ji) * (damp(ji) * zero * g11v(ji,5) +       &
                                       g22v(ji,5)) ) -                          &
                  fact4 / rho(ji) * ( lm(ji) * (damp(ji) * gu(ji,1,1) +         &
                                      gu(ji,2,2) + gu(ji,3,3))**2 +             &
                    two * mu(ji) * (                                            &
                    ( pt5 * (two * damp(ji) * gu(ji,1,1)) )**2 +                &
                    ( pt5 * (gu(ji,1,2) + damp(ji) * gu(ji,2,1)) )**2 +         &
                    ( pt5 * (gu(ji,1,3) + damp(ji) * gu(ji,3,1)) )**2 +         &
                    ( pt5 * (gu(ji,1,2) + damp(ji) * gu(ji,2,1)) )**2 +         &
                    ( pt5 * (two * gu(ji,2,2)) )**2 +                           &
                    ( pt5 * (gu(ji,2,3) + gu(ji,3,2)) )**2 +                    &
                    ( pt5 * (gu(ji,1,3) + damp(ji) * gu(ji,3,1)) )**2 +         &
                    ( pt5 * (gu(ji,3,2) + gu(ji,2,3)) )**2 +                    &
                    ( pt5 * (two * gu(ji,3,3)) )**2 ) )

          end do
        end do
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
        if (top.eq.1) then

          j = ny
          do i = 1, nx-1
            ji = j + (i-1) * ny
    
!.... get the metrics along this boundary

            bn1 = n1(ji) / sqrt( n1(ji)**2 + n2(ji)**2 )
            bn2 = n2(ji) / sqrt( n1(ji)**2 + n2(ji)**2 )

            bs1 = -bn2
            bs2 =  bn1

!.... compute Us and Un
        
            un = bn1 * vl(ji,2) + bn2 * vl(ji,3)
            us = bs1 * vl(ji,2) + bs2 * vl(ji,3)

            c = sqrt( t(ji) ) / Ma

!.... compute some derivatives in the boundary normal coordinate system

            gpn  = bn1 * gp(ji,1) + bn2 * gp(ji,2)
            gps  = bs1 * gp(ji,1) + bs2 * gp(ji,2)
        
            grhon = bn1 * grho(ji,1) + bn2 * grho(ji,2)
            grhos = bs1 * grho(ji,1) + bs2 * grho(ji,2)

            gtn  = bn1 * gt(ji,1) + bn2 * gt(ji,2)
            gts  = bs1 * gt(ji,1) + bs2 * gt(ji,2)

            gus  = bs1 * gu(ji,1,1) + bs2 * gu(ji,1,2)
            gvs  = bs1 * gu(ji,2,1) + bs2 * gu(ji,2,2)

            gunn = bn1*bn1 * gu(ji,1,1) + bn1*bn2 * gu(ji,2,1) + &
                   bn1*bn2 * gu(ji,1,2) + bn2*bn2 * gu(ji,2,2)

            gusn = bs1*bn1 * gu(ji,1,1) + bs2*bn1 * gu(ji,2,1) + &
                   bs1*bn2 * gu(ji,1,2) + bs2*bn2 * gu(ji,2,2)

            guns = bn1*bs1 * gu(ji,1,1) + bn2*bs1 * gu(ji,2,1) + &
                   bn1*bs2 * gu(ji,1,2) + bn2*bs2 * gu(ji,2,2)

            guss = bs1*bs1 * gu(ji,1,1) + bs2*bs1 * gu(ji,2,1) + &
                   bs1*bs2 * gu(ji,1,2) + bs2*bs2 * gu(ji,2,2)
            
            gu3n = bn1 * gu(ji,3,1) + bn2 * gu(ji,3,2)
            gu3s = bs1 * gu(ji,3,1) + bs2 * gu(ji,3,2)
            
!.... compute the Characteristic amplitudes from the interior
        
            L1 = ( un - c ) * ( gpn - rho(ji) * c * gunn )
            L2 = un * ( c**2 * grhon - gpn )
            L3 = un * gusn
            L4 = un * gu3n
            L5 = ( un + c ) * ( gpn + rho(ji) * c * gunn )

!.... use the potential solution

            b = one / pi
            if (Ma.lt.one) then
              rbeta = sqrt(one - Ma**2)
            else
              rbeta = sqrt(Ma**2 - one)
            end if
            ub1 = one + (x(j,1)-b)/pi/rbeta**2/((x(j,i)-b)**2 + &
                         rbeta**2*y(j,i)**2)
            ub2 = rbeta*y(j,i)/pi/((x(j,i)-b)**2+rbeta**2*y(j,i)**2)

            gub11 = one/( pi*rbeta**2*((x(j,i)-b)**2+rbeta**2*y(j,i)**2) ) - &
                    two*(x(j,i)-b)**2/( pi*rbeta**2*((x(j,i)-b)**2+rbeta**2*y(j,i)**2)**2 )

            gub12 = -two*(x(j,i)-b)*rbeta**2*y(j,i)/( pi*rbeta**2*((x(j,i)-b)**2+ &
                     rbeta**2*y(j,i)**2)**2 )

            gub21 = -rbeta*y(j,i)*two*(x(j,i)-b)/( pi*((x(j,i)-b)**2+ &
                     rbeta**2*y(j,i)**2)**2 )

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
              L1 = rk * ( p(ji) - pinf ) + (ubn - cb)*(gpbn - rhob*cb*gubnn)    
            end if
            
!.... compute the boundary normal derivative terms

            d1 = one / c**2 * ( L2 + pt5 * ( L5 + L1 ) )
            d2 = pt5 * ( L5 + L1 )
            d3 = pt5 / ( rho(ji) * c ) * ( L5 - L1 )
            d4 = L3
            d5 = L4

!.... continuity

            rl(ji,1) = d1 + grhos * us + grho(ji,3) * u3(ji) + &
                       rho(ji) * ( guss + gu(ji,3,3) )
    
!.... momentum

            rl(ji,2) = bn1 * d3 + bs1 * d4 + us * gus + bs1 * rhoinv(ji) * gps -        &
                  fact2 * ( damp(ji) * (g1lm(ji) * divu(ji) + lm(ji) * g1divu(ji)) +    &
                            two * (damp(ji) * g1mu(ji) * S(ji,1,1) +                    &
                                   g2mu(ji) * pt5 * (gu(ji,1,2) +                       &
                                   damp(ji) * gu(ji,2,1)) ) +                           &
                            mu(ji) * (damp(ji) * two * g11v(ji,2) +                     &
                            g22v(ji,2) + damp(ji) * g12v(ji,3) ) ) / rho(ji)
    
            rl(ji,3) = bn2 * d3 + bs2 * d4 + us * gvs + bs2 * rhoinv(ji) * gps -        &
                  fact2 * ( g2lm(ji) * (damp(ji) * gu(ji,1,1) + gu(ji,2,2) +            &
                                    gu(ji,3,3)) +                                       &
                            lm(ji) * (damp(ji) * g12v(ji,2) + g22v(ji,3)) +             &
                            two * (damp(ji) * g1mu(ji) * S(ji,2,1) +                    &
                                   g2mu(ji) * S(ji,2,2)) +                              &
                            mu(ji) * ( damp(ji) * (g11v(ji,3) + g12v(ji,2)) +           &
                                   two * g22v(ji,3) ) ) / rho(ji)
    
            rl(ji,4) = d5 + us * gu3s + gp(ji,3) * rhoinv(ji) -                         &
                  fact2 * ( two * (damp(ji) * g1mu(ji) * S(ji,3,1) +                    &
                                   g2mu(ji) * S(ji,3,2)) +                              &
                            mu(ji) * ( damp(ji) * g11v(ji,4) +                          &
                                   g22v(ji,4) ) ) / rho(ji)

!.... temperature

            rl(ji,5) = ( gamma * Ma**2 * d2 - t(ji) * d1 ) * rhoinv(ji) +               &
                  us * gts + u3(ji) * gt(ji,3) +                                        &
                  gamma1 * t(ji) * ( guss + gu(ji,3,3) )  -                             &
                  fact3 / rho(ji) * ( damp(ji) * g1con(ji) * gt(ji,1) +                 &
                                g2con(ji) * gt(ji,2) +                                  &
                                con(ji) * (damp(ji) * g11v(ji,5) +                      &
                                       g22v(ji,5)) ) -                                  &
                  fact4 / rho(ji) * ( lm(ji) * (damp(ji) * gu(ji,1,1) + gu(ji,2,2) +    &
                                      gu(ji,3,3))**2 +                                  &
                    two * mu(ji) * (                                                    &
                    ( pt5 * (two * damp(ji) * gu(ji,1,1)) )**2 +                        &
                    ( pt5 * (gu(ji,1,2) + damp(ji) * gu(ji,2,1)) )**2 +                 &
                    ( pt5 * (gu(ji,1,3) + damp(ji) * gu(ji,3,1)) )**2 +                 &
                    ( pt5 * (gu(ji,1,2) + damp(ji) * gu(ji,2,1)) )**2 +                 &
                    ( pt5 * (two * gu(ji,2,2)) )**2 +                                   &
                    ( pt5 * (gu(ji,2,3) + gu(ji,3,2)) )**2 +                            &
                    ( pt5 * (gu(ji,1,3) + damp(ji) * gu(ji,3,1)) )**2 +                 &
                    ( pt5 * (gu(ji,3,2) + gu(ji,2,3)) )**2 +                            &
                    ( pt5 * (two * gu(ji,3,3)) )**2 ) )

!           if ( un .le. zero ) then
!              rl(ji,2) = zero
!              rl(ji,3) = zero
!              rl(ji,4) = zero
!              rl(ji,5) = zero
!           end if

          end do

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
        if (wall.eq.4) then
        
!.... compute the wall normal pressure gradient using the momentum equations

          call wallbc(vl,wrho,pnorm)

          j = 1
          do i = 1, nx  
            ji = j + (i-1) * ny
    
!.... get the metrics along this boundary

            bn1 = n1(ji) / sqrt( n1(ji)**2 + n2(ji)**2 )
            bn2 = n2(ji) / sqrt( n1(ji)**2 + n2(ji)**2 )

!.... compute Us and Un
        
            un =  bn1 * vl(ji,2) + bn2 * vl(ji,3)
            us = -bn2 * vl(ji,2) + bn1 * vl(ji,3)

            c = sqrt( t(ji) ) / Ma

!.... compute some derivatives in the boundary normal coordinate system

            gpn =  bn1 * gp(ji,1) + bn2 * gp(ji,2)
!           gpn =  sqrt( n1(ji)**2 + n2(ji)**2 ) * pnorm(i)
            
            gps = -bn2 * gp(ji,1) + bn1 * gp(ji,2)
        
            grhon =  bn1 * grho(ji,1) + bn2 * grho(ji,2)
            grhos = -bn2 * grho(ji,1) + bn1 * grho(ji,2)

            gtn =  bn1 * gt(ji,1) + bn2 * gt(ji,2)
            gts = -bn2 * gt(ji,1) + bn1 * gt(ji,2)

            gus  =  -bn2 * gu(ji,1,1) + bn1 * gu(ji,1,2)
            gvs  =  -bn2 * gu(ji,2,1) + bn1 * gu(ji,2,2)

            gunn =  bn1*bn1 * gu(ji,1,1) + bn1*bn2 * gu(ji,2,1) + &
                    bn1*bn2 * gu(ji,1,2) + bn2*bn2 * gu(ji,2,2)
            gusn = -bn1*bn2 * gu(ji,1,1) + bn1*bn1 * gu(ji,2,1) - &
                    bn2*bn2 * gu(ji,1,2) + bn1*bn2 * gu(ji,2,2)
            guns = -bn2*bn1 * gu(ji,1,1) - bn2*bn2 * gu(ji,2,1) + &
                    bn1*bn1 * gu(ji,1,2) + bn1*bn2 * gu(ji,2,2)
            guss =  bn2*bn2 * gu(ji,1,1) - bn1*bn2 * gu(ji,2,1) - &
                    bn1*bn2 * gu(ji,1,2) + bn1*bn1 * gu(ji,2,2)
            
            gu3n =  bn1 * gu(ji,3,1) + bn2 * gu(ji,3,2)
            gu3s = -bn2 * gu(ji,3,1) + bn1 * gu(ji,3,2)
            
!.... compute the Characteristic amplitudes
        
            L1 = ( un - c ) * ( gpn - rho(ji) * c * gunn )
            L2 = un * ( c**2 * grhon - gpn )
            L3 = un * gusn
            L4 = un * gu3n
            L5 = ( un + c ) * ( gpn + rho(ji) * c * gunn )
        
!.... adiabatic wall boundary conditions
        
            L2 = zero
            L3 = zero
            L4 = zero
            L5 = L1
                    
!.... compute the boundary normal derivative terms

            d1 = one / c**2 * ( L2 + pt5 * ( L5 + L1 ) )
            d2 = pt5 * ( L5 + L1 )
            d3 = pt5 / ( rho(ji) * c ) * ( L5 - L1 )
            d4 = L3
            d5 = L4

!.... continuity

            rl(ji,1) = d1 + grhos * us + grho(ji,3) * u3(ji) + &
                       rho(ji) * ( guss + gu(ji,3,3) )
    
!.... momentum

            rl(ji,2) = bn1 * d3 - bn2 * d4 + us * gus - bn2 * rhoinv(ji) * gps -        &
                  fact2 * ( damp(ji) * (g1lm(ji) * divu(ji) + lm(ji) * g1divu(ji)) +    &
                            two * (damp(ji) * g1mu(ji) * S(ji,1,1) +                    &
                                   g2mu(ji) * pt5 * (gu(ji,1,2) +                       &
                                   damp(ji) * gu(ji,2,1)) ) +                           &
                            mu(ji) * (damp(ji) * two * g11v(ji,2) +                     &
                            g22v(ji,2) + damp(ji) * g12v(ji,3) ) ) / rho(ji)
    
            rl(ji,3) = bn2 * d3 + bn1 * d4 + us * gvs + bn1 * rhoinv(ji) * gps -        &
                  fact2 * ( g2lm(ji) * (damp(ji) * gu(ji,1,1) + gu(ji,2,2) +            &
                                    gu(ji,3,3)) +                                       &
                            lm(ji) * (damp(ji) * g12v(ji,2) + g22v(ji,3)) +             &
                            two * (damp(ji) * g1mu(ji) * S(ji,2,1) +                    &
                                   g2mu(ji) * S(ji,2,2)) +                              &
                            mu(ji) * ( damp(ji) * (g11v(ji,3) + g12v(ji,2)) +           &
                                   two * g22v(ji,3) ) ) / rho(ji)
    
            rl(ji,4) = d5 + us * gu3s + gp(ji,3) / rho(ji) -                            &
                  fact2 * ( two * (damp(ji) * g1mu(ji) * S(ji,3,1) +                    &
                                   g2mu(ji) * S(ji,3,2)) +                              &
                            mu(ji) * ( damp(ji) * g11v(ji,4) +                          &
                                   g22v(ji,4) ) ) / rho(ji)

!.... temperature

            rl(ji,5) = ( gamma * Ma**2 * d2 - t(ji) * d1 ) * rhoinv(ji) +               &
                  us * gts + u3(ji) * gt(ji,3) +                                        &
                  gamma1 * t(ji) * ( guss + gu(ji,3,3) )  -                             &
                  fact3 / rho(ji) * ( damp(ji) * g1con(ji) * gt(ji,1) +                 &
                                g2con(ji) * gt(ji,2) +                                  &
                                con(ji) * (damp(ji) * g11v(ji,5) +                      &
                                       g22v(ji,5)) ) -                                  &
                  fact4 / rho(ji) * ( lm(ji) * (damp(ji) * gu(ji,1,1) + gu(ji,2,2) +    &
                                      gu(ji,3,3))**2 +                                  &
                    two * mu(ji) * (                                                    &
                    ( pt5 * (two * damp(ji) * gu(ji,1,1)) )**2 +                        &
                    ( pt5 * (gu(ji,1,2) + damp(ji) * gu(ji,2,1)) )**2 +                 &
                    ( pt5 * (gu(ji,1,3) + damp(ji) * gu(ji,3,1)) )**2 +                 &
                    ( pt5 * (gu(ji,1,2) + damp(ji) * gu(ji,2,1)) )**2 +                 &
                    ( pt5 * (two * gu(ji,2,2)) )**2 +                                   &
                    ( pt5 * (gu(ji,2,3) + gu(ji,3,2)) )**2 +                            &
                    ( pt5 * (gu(ji,1,3) + damp(ji) * gu(ji,3,1)) )**2 +                 &
                    ( pt5 * (gu(ji,3,2) + gu(ji,2,3)) )**2 +                            &
                    ( pt5 * (two * gu(ji,3,3)) )**2 ) )

          end do

        end if

!=============================================================================!
!       W a l l   C o r n e r s 
!=============================================================================!
        if (.false.) then
        
        j = 1
          i = nx
            ji = j + (i-1) * ny
    
            c = sqrt( t(ji) ) / Ma

!.... compute the streamwise Characteristic amplitudes
        
            L1 = ( zero - c ) * ( gp(ji,1) - rho(ji) * c * zero )
            L2 = zero
            L3 = zero
            L4 = zero
            L5 = ( zero + c ) * ( gp(ji,1) + rho(ji) * c * zero )
        
!.... nonreflecting outflow boundary condition
        
            L1 = rk * ( p(ji) - pinf )  ! nonreflecting condition

!.... compute the streamwise derivative terms

            d1 = one / c**2 * ( L2 + pt5 * ( L5 + L1 ) )
            d2 = pt5 * ( L5 + L1 )
            d3 = pt5 / ( rho(ji) * c ) * ( L5 - L1 )
            d4 = L3
            d5 = L4

!.... compute the wall-normal Characteristic amplitudes
        
            L1 = ( zero - c ) * ( gp(ji,2) - rho(ji) * c * gu(ji,2,2) )
            L2 = zero * ( c**2 * grho(ji,2) - gp(ji,2) )
            L3 = zero * gu(ji,2,2)
            L4 = zero * gu(ji,3,2)
            L5 = ( zero + c ) * ( gp(ji,2) + rho(ji) * c * gu(ji,2,2) )
        
!.... nonreflecting outflow boundary condition
        
            L5 = L1

!.... compute the wall-normal derivative terms

            e1 = one / c**2 * ( L2 + pt5 * ( L5 + L1 ) )
            e2 = pt5 * ( L5 + L1 )
            e3 = pt5 / ( rho(ji) * c ) * ( L5 - L1 )
            e4 = L3
            e5 = L4

!.... continuity

            rl(ji,1) = d1 + e1 + grho(ji,3) * u3(ji) + rho(ji) * gu(ji,3,3)
    
!.... momentum

            rl(ji,2) = d3 + e4  -                                               &
                  fact2 * ( zero * damp(ji) * (g1lm(ji) * divu(ji) +            &
                            lm(ji) * g1divu(ji)) +                              &
                            two * (zero * damp(ji) * g1mu(ji) * S(ji,1,1) +     &
                                   g2mu(ji) * pt5 * (gu(ji,1,2) +               &
                                   damp(ji) * gu(ji,2,1)) ) +                   &
                            mu(ji) * (zero * damp(ji) * two * g11v(ji,2) +      &
                            g22v(ji,2) + damp(ji) * g12v(ji,3) ) ) / rho(ji)
    
            rl(ji,3) = d4 + e3 +                                                &
                  fact2 * ( g2lm(ji) * (damp(ji) * gu(ji,1,1) + gu(ji,2,2) +    &
                                    gu(ji,3,3)) +                               &
                            lm(ji) * (damp(ji) * g12v(ji,2) + g22v(ji,3)) +     &
                            two * (damp(ji) * zero * g1mu(ji) * S(ji,2,1) +     &
                                   g2mu(ji) * S(ji,2,2)) +                      &
                            mu(ji) * ( damp(ji) * (zero * g11v(ji,3) +          &
                                       zero * g12v(ji,2)) +                     &
                                       two * g22v(ji,3) ) ) / rho(ji)
    
            rl(ji,4) = d5 + e5 + gp(ji,3) / rho(ji) -                           &
                  fact2 * ( two * (damp(ji) * g1mu(ji) * S(ji,3,1) +            &
                                   g2mu(ji) * S(ji,3,2)) +                      &
                            mu(ji) * ( damp(ji) * g11v(ji,4) +                  &
                                   g22v(ji,4) ) ) / rho(ji)

!.... temperature

            rl(ji,5) = ( gamma * Ma**2 * (d2 + e2) -                            &
                         t(ji) * (d1 + e1) ) / rho(ji) +                        &
                  u3(ji) * gt(ji,3) + gamma1 * t(ji) * gu(ji,3,3) -             &
                  fact3 / rho(ji) * ( damp(ji) * zero * g1con(ji) * gt(ji,1) +  &
                                g2con(ji) * gt(ji,2) +                          &
                                con(ji) * (damp(ji) * zero * g11v(ji,5) +       &
                                       g22v(ji,5)) ) -                          &
                  fact4 / rho(ji) * ( lm(ji) * (damp(ji) * gu(ji,1,1) +         &
                                      gu(ji,2,2) + gu(ji,3,3))**2 +             &
                    two * mu(ji) * (                                            &
                    ( pt5 * (two * damp(ji) * gu(ji,1,1)) )**2 +                &
                    ( pt5 * (gu(ji,1,2) + damp(ji) * gu(ji,2,1)) )**2 +         &
                    ( pt5 * (gu(ji,1,3) + damp(ji) * gu(ji,3,1)) )**2 +         &
                    ( pt5 * (gu(ji,1,2) + damp(ji) * gu(ji,2,1)) )**2 +         &
                    ( pt5 * (two * gu(ji,2,2)) )**2 +                           &
                    ( pt5 * (gu(ji,2,3) + gu(ji,3,2)) )**2 +                    &
                    ( pt5 * (gu(ji,1,3) + damp(ji) * gu(ji,3,1)) )**2 +         &
                    ( pt5 * (gu(ji,3,2) + gu(ji,2,3)) )**2 +                    &
                    ( pt5 * (two * gu(ji,3,3)) )**2 ) )

        end if

!=============================================================================!
!       S l i p  W a l l
!=============================================================================!
!  Key to the generalized Poinsot & Lele boundary conditions.
!  Boundary conditions are required for all incomming 
!  characteristic amplitudes
!
!  n is the unit normal vector pointing INTO the fluid
!
!  Slip wall:           un = 0
!
!       L1 -> out
!       L2 -> zero
!       L3 -> zero
!       L4 -> zero
!       L5 -> in
!
!  BC:  L5 = L1 
!=============================================================================!
        if (.false. .and. left.eq.2) then
        
          i = 1
          do j = 1, ny  
            ji = j + (i-1) * ny
    
!.... get the metrics along this boundary

            bn1 = m1(ji) / sqrt( m1(ji)**2 + m2(ji)**2 )
            bn2 = m2(ji) / sqrt( m1(ji)**2 + m2(ji)**2 )

!.... compute Us and Un
        
            un =  bn1 * vl(ji,2) + bn2 * vl(ji,3)
            us = -bn2 * vl(ji,2) + bn1 * vl(ji,3)

            c = sqrt( t(ji) ) / Ma

!.... compute some derivatives in the boundary normal coordinate system

            gpn =  bn1 * gp(ji,1) + bn2 * gp(ji,2)
            gps = -bn2 * gp(ji,1) + bn1 * gp(ji,2)
        
            grhon =  bn1 * grho(ji,1) + bn2 * grho(ji,2)
            grhos = -bn2 * grho(ji,1) + bn1 * grho(ji,2)

            gtn =  bn1 * gt(ji,1) + bn2 * gt(ji,2)
            gts = -bn2 * gt(ji,1) + bn1 * gt(ji,2)

            gus  =  -bn2 * gu(ji,1,1) + bn1 * gu(ji,1,2)
            gvs  =  -bn2 * gu(ji,2,1) + bn1 * gu(ji,2,2)

            gunn =  bn1*bn1 * gu(ji,1,1) + bn1*bn2 * gu(ji,2,1) + &
                    bn1*bn2 * gu(ji,1,2) + bn2*bn2 * gu(ji,2,2)
            gusn = -bn1*bn2 * gu(ji,1,1) + bn1*bn1 * gu(ji,2,1) - &
                    bn2*bn2 * gu(ji,1,2) + bn1*bn2 * gu(ji,2,2)
            guns = -bn2*bn1 * gu(ji,1,1) - bn2*bn2 * gu(ji,2,1) + &
                    bn1*bn1 * gu(ji,1,2) + bn1*bn2 * gu(ji,2,2)
            guss =  bn2*bn2 * gu(ji,1,1) - bn1*bn2 * gu(ji,2,1) - &
                    bn1*bn2 * gu(ji,1,2) + bn1*bn1 * gu(ji,2,2)
            
            gu3n =  bn1 * gu(ji,3,1) + bn2 * gu(ji,3,2)
            gu3s = -bn2 * gu(ji,3,1) + bn1 * gu(ji,3,2)
            
!.... compute the Characteristic amplitudes
        
            L1 = ( un - c ) * ( gpn - rho(ji) * c * gunn )
            L2 = un * ( c**2 * grhon - gpn )
            L3 = un * gusn
            L4 = un * gu3n
            L5 = ( un + c ) * ( gpn + rho(ji) * c * gunn )
        
!.... adiabatic slip-wall boundary conditions
        
            L2 = zero
            L3 = zero
            L4 = zero
            L5 = L1
                    
!.... compute the boundary normal derivative terms

            d1 = one / c**2 * ( L2 + pt5 * ( L5 + L1 ) )
            d2 = pt5 * ( L5 + L1 )
            d3 = pt5 / ( rho(ji) * c ) * ( L5 - L1 )
            d4 = L3
            d5 = L4

!.... continuity

            rl(ji,1) = d1 + grhos * us + grho(ji,3) * u3(ji) + &
                       rho(ji) * ( guss + gu(ji,3,3) )
    
!.... momentum

            rl(ji,2) = bn1 * d3 - bn2 * d4 + us * gus - bn2 * rhoinv(ji) * gps -        &
                  fact2 * ( damp(ji) * (g1lm(ji) * divu(ji) + lm(ji) * g1divu(ji)) +    &
                            two * (damp(ji) * g1mu(ji) * S(ji,1,1) +                    &
                                   g2mu(ji) * pt5 * (zero * gu(ji,1,2) +                &
                                   damp(ji) * gu(ji,2,1)) ) +                           &
                            mu(ji) * (damp(ji) * two * g11v(ji,2) +                     &
                            g22v(ji,2) + damp(ji) * g12v(ji,3) ) ) / rho(ji)
    
            rl(ji,3) = bn2 * d3 + bn1 * d4 + us * gvs + bn1 * rhoinv(ji) * gps -        &
                  fact2 * ( g2lm(ji) * (damp(ji) * gu(ji,1,1) + gu(ji,2,2) +            &
                                    gu(ji,3,3)) +                                       &
                            lm(ji) * (damp(ji) * g12v(ji,2) + g22v(ji,3)) +             &
                            two * (damp(ji) * g1mu(ji) * zero * S(ji,2,1) +             &
                                   g2mu(ji) * S(ji,2,2)) +                              &
                            mu(ji) * ( damp(ji) * zero * (g11v(ji,3) + g12v(ji,2)) +    &
                                   two * g22v(ji,3) ) ) / rho(ji)
    
            rl(ji,4) = d5 + us * gu3s + gp(ji,3) / rho(ji) -                            &
                  fact2 * ( two * (damp(ji) * g1mu(ji) * S(ji,3,1) +                    &
                                   g2mu(ji) * S(ji,3,2)) +                              &
                            mu(ji) * ( damp(ji) * g11v(ji,4) +                          &
                                   g22v(ji,4) ) ) / rho(ji)

!.... temperature

            rl(ji,5) = ( gamma * Ma**2 * d2 - t(ji) * d1 ) * rhoinv(ji) +               &
                  us * gts + u3(ji) * gt(ji,3) +                                        &
                  gamma1 * t(ji) * ( guss + gu(ji,3,3) )  -                             &
                  fact3 / rho(ji) * ( damp(ji) * g1con(ji) * gt(ji,1) +                 &
                                g2con(ji) * zero * gt(ji,2) +                           &
                                con(ji) * (damp(ji) * g11v(ji,5) +                      &
                                       g22v(ji,5)) ) -                                  &
                  fact4 / rho(ji) * ( lm(ji) * (damp(ji) * gu(ji,1,1) + gu(ji,2,2) +    &
                                      gu(ji,3,3))**2 +                                  &
                    two * mu(ji) * (                                                    &
                    ( pt5 * (two * damp(ji) * gu(ji,1,1)) )**2 +                        &
                    ( pt5 * (zero * gu(ji,1,2) + zero * damp(ji) * gu(ji,2,1)) )**2 +   &
                    ( pt5 * (gu(ji,1,3) + damp(ji) * gu(ji,3,1)) )**2 +                 &
                    ( pt5 * (damp(ji) * zero * gu(ji,2,1) + zero * gu(ji,1,2)) )**2 +   &
                    ( pt5 * (gu(ji,2,2) + gu(ji,2,2)) )**2 +                            &
                    ( pt5 * (zero * gu(ji,2,3) + zero * gu(ji,3,2)) )**2 +              &
                    ( pt5 * (gu(ji,1,3) + damp(ji) * gu(ji,3,1)) )**2 +                 &
                    ( pt5 * (gu(ji,3,2) + gu(ji,2,3)) )**2 +                            &
                    ( pt5 * (two * gu(ji,3,3)) )**2 ) )

          end do

        end if

        return
        end

