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
!       Revised: 3-17-01
!=============================================================================!
        use global
        use local2d
        use material
        implicit none

        real :: rl(ndof,nx,ny), vl(ndof,nx,ny)
        real :: fact1, fact2, fact3, fact4
        
        real :: pinf, cinf, Lx, rk, cc
        real :: L1, L2, L3, L4, L5
        real :: d1, d2, d3, d4, d5
        real :: e1, e2, e3, e4, e5
        
        integer :: i, j, idof
        
        real :: g1vl(ndof), g2vl(ndof), g11vl(ndof), g12vl(ndof), g22vl(ndof)

        real :: bn1, bn2, bs1, bs2, us, un, gps, gpn, grhos, grhon, gts, gtn, &
                gus, gvs, guss, gusn, guns, gunn, gu3s, gu3n
                
        real :: b, rbeta, ub1, ub2, gub11, gub12, gub21, gub22
        real :: ubn, ubs, gubnn, gubsn, tb, rhob, pb, cb, gtbn, grhobn, gpbn
        
        real :: wrho(nx), pnorm(nx)
!=============================================================================!

        fact1 = one / (gamma * Ma**2)
        fact2 = one / Re
        fact3 = gamma / ( Pr * Re )
        fact4 = gamma * gamma1 * Ma**2 / Re

!=============================================================================!
!       W a l l
!=============================================================================!
!  Poinsot & Lele boundary conditions.
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

        j = 1
        do i = 1, nx
        
        if (general) then

        do idof = 1, ndof
          g11vl(idof) = g11v(idof,i,j)       * m1m1(i,j)        + &
                        two * g12v(idof,i,j) * m1n1(i,j)        + &
                        g22v(idof,i,j)       * n1n1(i,j)        + &
                        g1v(idof,i,j)        * m11(i,j)         + &
                        g2v(idof,i,j)        * n11(i,j)

          g12vl(idof) = g11v(idof,i,j) * m1m2(i,j)                + &
                        g12v(idof,i,j) * (m1n2(i,j) + m2n1(i,j) ) + &
                        g22v(idof,i,j) * n1n2(i,j)                + &
                        g1v(idof,i,j)  * m12(i,j)                 + &
                        g2v(idof,i,j)  * n12(i,j)

          g22vl(idof) = g11v(idof,i,j)       * m2m2(i,j)        + &
                        two * g12v(idof,i,j) * m2n2(i,j)        + &
                        g22v(idof,i,j)       * n2n2(i,j)        + &
                        g1v(idof,i,j)        * m22(i,j)         + &
                        g2v(idof,i,j)        * n22(i,j)

          g1vl(idof) = g1v(idof,i,j)*m1(i,j) + g2v(idof,i,j)*n1(i,j)
          g2vl(idof) = g1v(idof,i,j)*m2(i,j) + g2v(idof,i,j)*n2(i,j)
        end do

        else  ! Cartesian

        do idof = 1, ndof
          g11vl(idof) = g11v(idof,i,j) * m1m1(i,j) + g1v(idof,i,j)  * m11(i,j)

          g12vl(idof) = g12v(idof,i,j) * m1n2(i,j)
                        
          g22vl(idof) = g22v(idof,i,j) * n2n2(i,j) + g2v(idof,i,j)  * n22(i,j)

          g1vl(idof) = g1v(idof,i,j) * m1(i,j)
          g2vl(idof) = g2v(idof,i,j) * n2(i,j)
        end do

        end if

!=============================================================================!
!.... setup local variables
!=============================================================================!

        rho     = vl(1,i,j)
        u1      = vl(2,i,j)
        u2      = vl(3,i,j)
        u3      = vl(4,i,j)
        t       = vl(5,i,j)
        rhoinv  = one / rho
        p       = fact1 * rho * t

        grho(1) = g1vl(1)
        grho(2) = g2vl(1)
        grho(3) = zero

        gu(1,1) = g1vl(2)
        gu(1,2) = g2vl(2)
        gu(1,3) = zero
        
        gu(2,1) = g1vl(3)
        gu(2,2) = g2vl(3)
        gu(2,3) = zero

        gu(3,1) = g1vl(4)
        gu(3,2) = g2vl(4)
        gu(3,3) = zero

        gt(1)   = g1vl(5)
        gt(2)   = g2vl(5)
        gt(3)   = zero

!.... compute the gradient of the divergence of um

        g1divu = g11vl(2) + g12vl(3)
        g2divu = g12vl(2) + g22vl(3)
        g3divu = zero

!.... compute strain rate tensor derivative for the viscous terms

        S1jj = pt5 * ( g11vl(2) + g11vl(2) + g22vl(2) + g12vl(3) )
                
        S2jj = pt5 * ( g11vl(3) + g12vl(2) + g22vl(3) + g22vl(3) )
        
        S3jj = pt5 * ( g11vl(4) + g22vl(4) )

!.... compute Laplacian of Temperature for the heat conduction term

        Lapt = g11vl(5) + g22vl(5)
        
!.... compute divergence
  
        divu = gu(1,1) + gu(2,2) + gu(3,3)
            
!.... compute gradient of mean pressure using chain-rule

        gp(1) = fact1 * ( grho(1) * t + rho * gt(1) )
        gp(2) = fact1 * ( grho(2) * t + rho * gt(2) )
        gp(3) = fact1 * ( grho(3) * t + rho * gt(3) )

!.... compute strain rate tensor for the viscous terms

        S(1,1) = gu(1,1)
        S(2,1) = pt5 * ( gu(2,1) + gu(1,2) )
        S(3,1) = pt5 * gu(3,1)

        S(1,2) = pt5 * ( gu(1,2) + gu(2,1) )
        S(2,2) = gu(2,2)
        S(3,2) = pt5 * gu(3,2)

        S(1,3) = pt5 * gu(3,1) 
        S(2,3) = pt5 * gu(3,2)

!.... compute the material properties

        call getmat(t, mu,   lm,    con,  &
                       dmu,  d2mu,  dlm,  &
                       d2lm, dcon,  d2con )

!.... compute gradients of viscosity using chain-rule

        g1mu = dmu * gt(1)
        g2mu = dmu * gt(2)
        g3mu = dmu * gt(3)

        g1dmu = d2mu * gt(1)
        g2dmu = d2mu * gt(2)
        g3dmu = d2mu * gt(3)

!.... compute gradients of conductivity using chain-rule

        g1con = dcon * gt(1)
        g2con = dcon * gt(2)
        g3con = dcon * gt(3)

        g1dcon = d2con * gt(1)
        g2dcon = d2con * gt(2)
        g3dcon = d2con * gt(3)

!.... compute gradients of bulk viscosity using chain-rule

        g1lm = dlm * gt(1)
        g2lm = dlm * gt(2)
        g3lm = dlm * gt(3)

        g1dlm = d2lm * gt(1)
        g2dlm = d2lm * gt(2)
        g3dlm = d2lm * gt(3)

!.... get outward pointing normal along this boundary

        bn1 = bnb(i,1)
        bn2 = bnb(i,2)

!.... compute Us and Un
        
        un =  bn1 * u1 + bn2 * u2
        us = -bn2 * u1 + bn1 * u2

        cc = sqrt( t ) / Ma

!.... compute some derivatives in the boundary normal coordinate system

        gpn =  bn1 * gp(1) + bn2 * gp(2)
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
        
        L1 = ( un - cc ) * ( gpn - rho * cc * gunn )
        L2 = un * ( cc**2 * grhon - gpn )
        L3 = un * gusn
        L4 = un * gu3n
        L5 = ( un + cc ) * ( gpn + rho * cc * gunn )
        
!.... adiabatic wall boundary conditions
        
        L1 = L5
        L2 = zero
        L3 = zero
        L4 = zero
                    
!.... compute the boundary normal derivative terms

        d1 = one / cc**2 * ( L2 + pt5 * ( L5 + L1 ) )
        d2 = pt5 * ( L5 + L1 )
        d3 = pt5 / ( rho * cc ) * ( L5 - L1 )
        d4 = L3
        d5 = L4

!.... continuity

        rl(1,i,j) = d1 + grhos * us + grho(3) * u3 + &
                    rho * ( guss + gu(3,3) )
    
!.... momentum

        rl(2,i,j) = bn1 * d3 - bn2 * d4 + us * gus - bn2 * rhoinv * gps -     &
                  fact2 * ( (g1lm * divu + lm * g1divu) +                     &
                  two * (g1mu * S(1,1) + g2mu * pt5 * (gu(1,2) + gu(2,1)) ) + &
                  mu * (two * g11vl(2) + g22vl(2) + g12vl(3) ) ) / rho
    
        rl(3,i,j) = bn2 * d3 + bn1 * d4 + us * gvs + bn1 * rhoinv * gps -     &
                    fact2 * ( g2lm * (gu(1,1) + gu(2,2) + gu(3,3) ) +         &
                            lm * (g12vl(2) + g22vl(3)) +                      &
                            two * (g1mu * S(2,1) + g2mu * S(2,2)) +           &
                     mu * ( (g11vl(3) + g12vl(2)) + two * g22vl(3) ) ) / rho
    
        rl(4,i,j) = d5 + us * gu3s + gp(3) / rho -                            &
                    fact2 * ( two * (g1mu * S(3,1) + g2mu * S(3,2)) +         &
                            mu * ( g11vl(4) + g22vl(4) ) ) / rho

!.... temperature

        rl(5,i,j) = ( gamma * Ma**2 * d2 - t * d1 ) * rhoinv +            &
                  us * gts + u3 * gt(3) +                                 &
                  gamma1 * t * ( guss + gu(3,3) )  -                      &
                  fact3 / rho * ( g1con * gt(1) + g2con * gt(2) +         &
                                  con * (g11vl(5) + g22vl(5)) ) -         &
                  fact4 / rho * ( lm * (gu(1,1) + gu(2,2) + gu(3,3))**2 + &
                    two * mu * (                                          &
                    ( pt5 * (two * gu(1,1)) )**2 +                        &
                    ( pt5 * (gu(1,2) + gu(2,1)) )**2 +                    &
                    ( pt5 * (gu(1,3) + gu(3,1)) )**2 +                    &
                    ( pt5 * (gu(1,2) + gu(2,1)) )**2 +                    &
                    ( pt5 * (two * gu(2,2)) )**2 +                        &
                    ( pt5 * (gu(2,3) + gu(3,2)) )**2 +                    &
                    ( pt5 * (gu(1,3) + gu(3,1)) )**2 +                    &
                    ( pt5 * (gu(3,2) + gu(2,3)) )**2 +                    &
                    ( pt5 * (two * gu(3,3)) )**2 ) )

        end do
        end if

        end subroutine rhs_l
