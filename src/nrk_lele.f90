!=============================================================================!
        subroutine nrk(vl, rl, dtl) 
!  
!  Computes the RK RHS for the nonlinear compressible Navier-Stokes solver. 
!
!  Written: 9-14-95
!
!  Revised: 6-11-96     Cleaned up extraneous calculations
!  Revised: 11-9-99     Added SGI parallel directives
!  Revised: 03-16-01    Try out the Lele-Poinsot form of the equations
!=============================================================================!
        use global
        use local2d
        use material
        implicit none

        real :: vl(ndof,nx,ny), rl(ndof,nx,ny), dtl(nx,ny)
        !$sgi distribute vl(*,*,block), rl(*,*,block), dtl(*,block)

        integer :: ier, i, j, idof
        real    :: fact1, fact2, fact3, fact4

        real :: gpr(2,nx,ny)
        real :: g1vl(ndof), g2vl(ndof), g11vl(ndof), g12vl(ndof), g22vl(ndof)

        real*4 :: cpul
        real*4, external :: second

        integer, external :: mp_my_threadnum

        real :: bn1, bn2, un, us, cc, gpn, gps, grhon, grhos, gus, gvs
        real :: gunn, gusn, guns, guss, gu3n, gu3s, l1, l2, l3, l4, l5
        real :: d1, d2, d3, d4, d5, gtn, gts
!=============================================================================!

!.... Compute first derivatives of field in the mapped space

        call grad(ndof, nx, ny, vl, g1v, g2v, dxi, deta, optx, opty, &
                  xper, yper, lsym, rsym, bsym, tsym, .true.)
        
!.... Compute second derivatives of field
        
        call grad2(ndof, nx, ny, vl, g1v, g11v, g12v, g22v, dxi, deta, &
                   optx, opty, xper, yper, lsym, rsym, bsym, tsym, .true.)

!.... enforce adiabatic wall BC on the field

        call gradbc( g1v, g2v, g11v, g12v, g22v )
        
!.... damp the second derivatives in the viscous terms and 
!.... compute the pressure gradient on the outflow plane

        fact1 = one / (gamma * Ma**2)
        fact2 = one / Re
        fact3 = gamma / ( Pr * Re )
        fact4 = gamma * gamma1 * Ma**2 / Re

        if (idamp.eq.1) then
          !$doacross local(i,idof)
          !$omp parallel do private(i,idof)
          do j = 1, ny
            do i = 1, nx
              do idof = 1, ndof
                g11v(idof,i,j) = damp(i,j) * g11v(idof,i,j)
              end do
            end do
          end do
          call calcp( vl, g2v, m1, m2, n1, n2, gpr )
        end if

        cpul = second()

!=============================================================================!
!                       N O N L I N E A R   R H S
!=============================================================================!

!!$doacross local(idof,rho,u1,u2,u3,t,rhoinv,p,g1divu,g2divu,s1jj,s2jj,s3jj, &
!!$& Lapt,divu,mu,lm,con,dmu,d2mu,dlm,d2lm,dcon,d2con,g1mu,g2mu,g1dmu, &
!!$& g2dmu,g1con,g2con,g1dcon,g2dcon,g1lm,g2lm,g1dlm,g2dlm,g11vl,g12vl,g22vl, &
!!$& g1vl, g2vl, grho, gt, gu, gp, S)

!$omp parallel do private &
!$omp (i,idof,rho,u1,u2,u3,t,rhoinv,p,g1divu,g2divu,s1jj,s2jj,&
!$omp s3jj,Lapt,divu,mu,lm,con,dmu,d2mu,dlm,d2lm,dcon,d2con,g1mu,g2mu,g1dmu,&
!$omp g2dmu,g1con,g2con,g1dcon,g2dcon,g1lm,g2lm,g1dlm,g2dlm,g11vl,g12vl,g22vl,&
!$omp g1vl, g2vl, grho, gt, gu, gp, S, bn1, bn2, un, us, cc, gpn, gps, grhon, &
!$omp grhos, gus, gvs, gunn, gusn, guns, guss, gu3n, gu3s, &
!$omp l1, l2, l3, l4, l5, d1, d2, d3, d4, d5, gtn, gts, &
!$omp g3mu, g3dmu, g3con, g3dcon, g3lm, g3dlm, g3divu)
loop_j: do j = 1, ny
loop_i: do i = 1, nx

!.... It requries alot of FLOPS to use a generalized coordinate mapping.
!.... The "general" flag gives you the option of using a simple stretched
!.... Cartesian system which is nearly twice as fast!

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

!.... continuity

        rl(1,i,j) = rho * divu + grho(1) * u1 + grho(2) * u2

!.... momentum

        rl(2,i,j) = u1 * gu(1,1) + u2 * gu(1,2) + rhoinv * gp(1) - &
                   fact2 * ( g1lm * divu + lm * g1divu + &
                   two * (g1mu * S(1,1) + g2mu * S(1,2) + &
                                   mu * S1jj ) ) * rhoinv

        rl(3,i,j) = u1 * gu(2,1) + u2 * gu(2,2) + rhoinv * gp(2) - &
                   fact2 * ( g2lm * divu + lm * g2divu + &
                   two * (g1mu * S(2,1) + g2mu * S(2,2) +  &
                                   mu * S2jj ) ) * rhoinv

        rl(4,i,j) = u1 * gu(3,1) + u2 * gu(3,2) -                 &
                   fact2 * ( two * (g1mu * S(3,1) +  g2mu * S(3,2) + &
                             mu * S3jj ) ) * rhoinv

!.... temperature

        rl(5,i,j) = u1 * gt(1) + u2 * gt(2) + &
                   gamma1 * t * divu - rhoinv * ( &
                   fact3 * ( g1con * gt(1) + g2con * gt(2) + con * Lapt ) +  &
                   fact4 * ( lm * divu**2 + two * mu * ( &
                             S(1,1)**2 + S(1,2)**2 + S(1,3)**2 + &
                             S(2,1)**2 + S(2,2)**2 + S(2,3)**2 + &
                             S(3,1)**2 + S(3,2)**2 ) ) )

!.... try Poinsot & Lele boundary conditions

        if (wall == 4 .and. j == 1) then    

!.... get the metrics along this boundary

          bn1 = -bnb(i,1)
          bn2 = -bnb(i,2)

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
        
          L2 = zero
          L3 = zero
          L4 = zero
          L5 = L1
                    
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

          rl(2,i,j) = bn1 * d3 - bn2 * d4 + us * gus - bn2 * rhoinv * gps - &
                  fact2 * ( (g1lm * divu + lm * g1divu) +                     &
                  two * (g1mu * S(1,1) + g2mu * pt5 * (gu(1,2) + gu(2,1)) ) + &
                  mu * (two * g11vl(2) + g22vl(2) + g12vl(3) ) ) / rho
    
          rl(3,i,j) = bn2 * d3 + bn1 * d4 + us * gvs + bn1 * rhoinv * gps - &
                  fact2 * ( g2lm * (gu(1,1) + gu(2,2) + gu(3,3) ) +         &
                            lm * (g12vl(2) + g22vl(3)) +                    &
                            two * (g1mu * S(2,1) + g2mu * S(2,2)) +         &
                     mu * ( (g11vl(3) + g12vl(2)) + two * g22vl(3) ) ) / rho
    
          rl(4,i,j) = d5 + us * gu3s + gp(3) / rho -                    &
                  fact2 * ( two * (g1mu * S(3,1) + g2mu * S(3,2)) +     &
                            mu * ( g11vl(4) + g22vl(4) ) ) / rho

!.... temperature

          rl(5,i,j) = ( gamma * Ma**2 * d2 - t * d1 ) * rhoinv +        &
                  us * gts + u3 * gt(3) +                               &
                  gamma1 * t * ( guss + gu(3,3) )  -                    &
                  fact3 / rho * ( g1con * gt(1) + g2con * gt(2) +       &
                                  con * (g11vl(5) + g22vl(5)) ) -       &
                  fact4 / rho * ( lm * (gu(1,1) + gu(2,2) + gu(3,3))**2 + &
                    two * mu * (                                &
                    ( pt5 * (two * gu(1,1)) )**2 +              &
                    ( pt5 * (gu(1,2) + gu(2,1)) )**2 +          &
                    ( pt5 * (gu(1,3) + gu(3,1)) )**2 +          &
                    ( pt5 * (gu(1,2) + gu(2,1)) )**2 +          &
                    ( pt5 * (two * gu(2,2)) )**2 +              &
                    ( pt5 * (gu(2,3) + gu(3,2)) )**2 +          &
                    ( pt5 * (gu(1,3) + gu(3,1)) )**2 +          &
                    ( pt5 * (gu(3,2) + gu(2,3)) )**2 +          &
                    ( pt5 * (two * gu(3,3)) )**2 ) )
        end if

        end do loop_i
        end do loop_j

!.... correct boundaries using Lele & Poinsot BC's

        call rhs_l(rl,vl)
        
!.... standard sponge

        if (ispg.gt.0) then
          if (ispg.eq.1) then
            call spg_it(rl,vl,spg)
          else if (ispg.eq.2) then
            !$doacross local(i)
            !$omp parallel do private(i)
            do j = 1, ny
              do i = 1, nx
                rl(1,i,j) = rl(1,i,j) + &
                            (spg(i,j) + spg2(i,j)) * ( vl(1,i,j) - one )
                rl(2,i,j) = rl(2,i,j) + &
                            (spg(i,j) + spg2(i,j)) * ( vl(2,i,j) )
                rl(3,i,j) = rl(3,i,j) + &
                            (spg(i,j) + spg2(i,j)) * ( vl(3,i,j) )
                rl(4,i,j) = rl(4,i,j) + &
                            (spg(i,j) + spg2(i,j)) * ( vl(4,i,j) )
                rl(5,i,j) = rl(5,i,j) + &
                            (spg(i,j) + spg2(i,j)) * ( vl(5,i,j) - one )
              end do
            end do
          else
            call error('nrk$','ispg > 2 is not supported$')
          end if
        end if

!.... explicit smoother

        if (eps_e.ne.zero) call smoother( rl, vl )

!       write(*,*) 'RHS ',second()-cpul
        cpul = second()

!.... multiply by the local time step.  
!.... Negate because RK is written as y' = f(y)

        !$doacross local(i,idof)
        !$omp parallel do private(i,idof)
        do j = 1, ny
          do i = 1, nx
            do idof = 1, ndof
              rl(idof,i,j) = -dtl(i,j) * rl(idof,i,j)
            end do
          end do
        end do

!.... satisfy the boundary conditions on the rhs

        call rkbc(rl)

        return
        end
        
