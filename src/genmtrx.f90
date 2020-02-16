!=============================================================================!
        subroutine genmtrx(vl, kzl, rl) 
!  
!       Generates the linearized coefficient matrices for the
!       compressible NS equations 
!
! Input:
!   vl:    The field for which the matrices are computed
!   kzl:   The local spanwise wave number
!
! Output:
!   matrices:  Stored on the SSD
!
!       Revised: 5-10-95
!
!       1) added compact storage for G and Vij
!       2) improved G.Q product
!       3) now compute gradients of p, mu, con, lm using chain rule
!          instead of calling grad individually
!       4) Jump around C matrix if 2D
!
!       Revised: 5-23-95
!
!       1) removed mean flow designation to make for general
!       2) now do allocation external to the routine
!       3) merged with conmtrx
!
!       Revised: 9-9-97
!
!       1) removed SSD for use on non-Cray computers
!
!       Revised: 12-14-00
!
!       1) Switched indices
!
!=============================================================================!
        use global
        use local2d
        use local3d
        use pot
        use material
        implicit none
        
        real    :: vl(ndof,nx,ny), rl(ndof,nx,ny), kzl
!=============================================================================!
        integer ier, i, j, idof
        character*80 :: name, code='genmtrx$'
        real :: fact

        integer :: lrec, istat

        real :: A(ndof,ndof), B(ndof,ndof), C(ndof,ndof), D(ndof,ndof), Vij(4)

        real :: gpr(2,nx,ny)
        real :: g1vl(ndof), g2vl(ndof), g11vl(ndof), g12vl(ndof), g22vl(ndof)

        real :: fact1, fact2, fact3, gmsinv

        real :: fact4, bn1, bn2, un, us, cc, gpn, gps, grhon, grhos, gus, gvs
        real :: gunn, gusn, guns, guss, gu3n, gu3s, l1, l2, l3, l4, l5
        real :: d1, d2, d3, d4, d5, gtn, gts
!=============================================================================!

!.... Compute first derivatives of field in the mapped space

        call grad(ndof, nx, ny, vl, g1v, g2v, dxi, deta, optx, opty, &
                  xper, yper, lsym, rsym, bsym, tsym, carp)
        
!.... Compute second derivatives of field
        
        call grad2(ndof, nx, ny, vl, g1v, g11v, g12v, g22v, dxi, deta, &
                   optx, opty, xper, yper, lsym, rsym, bsym, tsym, carp)

!.... enforce adiabatic wall BC on the field

        if (Navier) call gradbc( g1v, g2v, g11v, g12v, g22v )
        
!.... damp the second derivatives in the viscous terms and 
!.... compute the pressure gradient on the outflow plane

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
        
!$omp parallel do private &
!$omp (i,idof,rho,u1,u2,u3,t,rhoinv,p,g1divu,g2divu,s1jj,s2jj,&
!$omp s3jj,Lapt,divu,mu,lm,con,dmu,d2mu,dlm,d2lm,dcon,d2con,g1mu,g2mu,g1dmu,&
!$omp g2dmu,g1con,g2con,g1dcon,g2dcon,g1lm,g2lm,g1dlm,g2dlm,g11vl,g12vl,g22vl,&
!$omp g1vl, g2vl, grho, gt, gu, gp, S, A, B, C, D, Vij, fact, fact1, fact2,&
!$omp fact3, gmsinv,g3divu,g3mu,g3dmu,g3con,g3dcon,g3lm,g3dlm,&
!$omp fact4, bn1, bn2, un, us, cc, gpn, gps, grhon, grhos, gus, gvs, &
!$omp gunn, gusn, guns, guss, gu3n, gu3s, l1, l2, l3, l4, l5, &
!$omp d1, d2, d3, d4, d5, gtn, gts)
loop_j: do j = 1, ny
loop_i: do i = 1, nx

!.... transform the gradients to physical space

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

        fact = one / (gamma * Ma**2)
        fact1 = one / Re
        fact2 = gamma / ( Pr * Re )
        fact3 = gamma * gamma1 * Ma**2 / Re

        rho     = vl(1,i,j)
        u1      = vl(2,i,j)
        u2      = vl(3,i,j)
        u3      = vl(4,i,j)
        t       = vl(5,i,j)
        rhoinv  = one / rho

        p       = fact * rho * t

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

        gp(1) = fact * ( grho(1) * t + rho * gt(1) )
        gp(2) = fact * ( grho(2) * t + rho * gt(2) )
        gp(3) = fact * ( grho(3) * t + rho * gt(3) )

!.... compute strain rate tensor for the viscous terms

        S = zero

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

!.... correct the pressure near the outflow plane

        if (idamp.eq.1 .and. linear.eq.0) then
          gp(1) = damp(i,j) * gp(1) + (one-damp(i,j)) * gpr(1,i,j)
          gp(2) = damp(i,j) * gp(2) + (one-damp(i,j)) * gpr(2,i,j)
        end if
        
        if (linear.eq.1 .and. wall.eq.10) &
             call genbump( g1v, g2v, g11v, g12v, g22v )

!==========================================================================!
!                       N O N L I N E A R   R H S
!==========================================================================!

        if (linear.eq.0) then

!.... continuity

        rl(1,i,j) = rho * divu + grho(1) * u1 + grho(2) * u2

!.... momentum

        rl(2,i,j) = u1 * gu(1,1) + u2 * gu(1,2) + rhoinv * gp(1) - &
                   fact1 * ( g1lm * divu + lm * g1divu + &
                   two * (g1mu * S(1,1) + g2mu * S(1,2) + mu * S1jj) ) * rhoinv

        rl(3,i,j) = u1 * gu(2,1) + u2 * gu(2,2) + rhoinv * gp(2) - &
                   fact1 * ( g2lm * divu + lm * g2divu + &
                   two * (g1mu * S(2,1) + g2mu * S(2,2) + mu * S2jj) ) * rhoinv

        rl(4,i,j) = u1 * gu(3,1) + u2 * gu(3,2) - &
                   fact1 * ( two * (g1mu * S(3,1) + &
                            g2mu * S(3,2) + mu * S3jj ) ) * rhoinv

!.... temperature

        rl(5,i,j) = u1 * gt(1) + u2 * gt(2) + gamma1 * t * divu - rhoinv * ( &
                    fact2 * ( g1con * gt(1) + g2con * gt(2) + con * Lapt ) + &
                    fact3 * ( lm * divu**2 + two * mu * (       &
                        S(1,1)**2 + S(1,2)**2 + S(1,3)**2 + &
                        S(2,1)**2 + S(2,2)**2 + S(2,3)**2 + &
                        S(3,1)**2 + S(3,2)**2 ) ) )

!.... try Poinsot & Lele boundary conditions

        if (wall == 4 .and. j == 1) then    

          fact1 = one / (gamma * Ma**2)
          fact2 = one / Re
          fact3 = gamma / ( Re * Pr )
          fact4 = gamma * gamma1 * Ma**2 / Re

!.... get the metrics along this boundary

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
                  fact2 * ( (g1lm * divu + lm * g1divu) +                   &
                            two * (g1mu * S(1,1) +                          &
                                   g2mu * pt5 * (gu(1,2) + gu(2,1)) ) + &
                  mu * (two * g11vl(2) + g22vl(2) + g12vl(3) ) ) / rho
    
          rl(3,i,j) = bn2 * d3 + bn1 * d4 + us * gvs + bn1 * rhoinv * gps - &
                  fact2 * ( g2lm * (gu(1,1) + gu(2,2) + gu(3,3)) +          &
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

!.... standard sponge

        if (ispg .gt. 0) then
          if (ispg .eq. 1) then
            call error(code,'ispg .eq. 1 is not working$')
            call spg_it(rl,vl,spg)
          else if (ispg .eq. 2) then
              rl(1,i,j) = rl(1,i,j) + (spg(i,j)+spg2(i,j)) * (vl(1,i,j) - one)
              rl(2,i,j) = rl(2,i,j) + (spg(i,j)+spg2(i,j)) * (vl(2,i,j))
              rl(3,i,j) = rl(3,i,j) + (spg(i,j)+spg2(i,j)) * (vl(3,i,j))
              rl(4,i,j) = rl(4,i,j) + (spg(i,j)+spg2(i,j)) * (vl(4,i,j))
              rl(5,i,j) = rl(5,i,j) + (spg(i,j)+spg2(i,j)) * (vl(5,i,j) - one)
          else
            call error('rhs_p$','ispg > 2 is not supported$')
          end if
        end if

       end if

!==========================================================================!
!               L I N E A R I Z E D   M A T R I C E S
!==========================================================================!

!==========================================================================!
!.... compute the A matrix
!==========================================================================!

!.... Continuity equation

        A(1,1) = u1
        A(1,2) = rho
        A(1,3) = zero
        A(1,4) = zero
        A(1,5) = zero

!.... Momentum equation -- x_1

        fact1  = rhoinv / Re
        gmsinv = one / (gamma * Ma**2)

        A(2,1) = rhoinv * t * gmsinv
        A(2,2) = u1 - fact1 * g1lm - fact1 * two * g1mu
        A(2,3) = -fact1 * g2mu
        A(2,4) = -fact1 * g3mu
        A(2,5) = gmsinv - fact1 * dlm * divu - &
                 fact1 * dmu * two * S(1,1)

!.... Momentum equation -- x_2

        A(3,1) = zero
        A(3,2) = -fact1 * g2lm
        A(3,3) = u1 - fact1 * g1mu
        A(3,4) = zero
        A(3,5) = -fact1 * dmu * two * S(2,1)

!.... Momentum equation -- x_3

        A(4,1) = zero
        A(4,2) = -fact1 * g3lm
        A(4,3) = zero
        A(4,4) = u1 - fact1 * g1mu
        A(4,5) = -fact1 * dmu * two * S(3,1)

!.... Energy equation

        fact1 = gamma * rhoinv / (Pr * Re)
        fact2 = gamma * gamma1 * Ma**2 * rhoinv / Re

        A(5,1) =  zero
        A(5,2) =  gamma1 * t - fact2 * two * lm * divu - &
                  four * fact2 * mu * S(1,1)
        A(5,3) = -four * fact2 * mu * S(2,1)
        A(5,4) = -four * fact2 * mu * S(3,1)
        A(5,5) = u1 - fact1 * (g1con + dcon * gt(1))

!==========================================================================!
!.... compute the B matrix
!==========================================================================!

!.... Continuity equation

        B(1,1) = u2
        B(1,2) = zero
        B(1,3) = rho
        B(1,4) = zero
        B(1,5) = zero

!.... Momentum equation -- x_1

        fact1  = rhoinv / Re
        gmsinv = one / (gamma * Ma**2)

        B(2,1) = zero
        B(2,2) = u2 - fact1 * g2mu
        B(2,3) = -fact1 * g1lm
        B(2,4) = zero
        B(2,5) = -fact1 * dmu * two * S(1,2)

!.... Momentum equation -- x_2

        B(3,1) = rhoinv * t * gmsinv
        B(3,2) = -fact1 * g1mu
        B(3,3) = u2 - fact1 * g2lm - fact1 * two * g2mu
        B(3,4) = -fact1 * g3mu
        B(3,5) = gmsinv - fact1 * dlm * divu - &
                 fact1 * dmu * two * S(2,2)

!.... Momentum equation -- x_3 

        B(4,1) = zero
        B(4,2) = zero
        B(4,3) = -fact1 * g3lm
        B(4,4) = u2 - fact1 * g2mu
        B(4,5) = -fact1 * dmu * two * S(3,2)

!.... Energy equation

        fact1 = gamma * rhoinv / (Pr * Re)
        fact2 = gamma * gamma1 * Ma**2 * rhoinv / Re

        B(5,1) =  zero 
        B(5,2) = -four * fact2 * mu * S(1,2)
        B(5,3) =  gamma1 * t - fact2 * two * lm * divu - &
                  four * fact2 * mu * S(2,2)
        B(5,4) = -four * fact2 * mu * S(3,2)
        B(5,5) =  u2 - fact1 * (g2con + dcon * gt(2))

!==========================================================================!
!.... compute the Vij matrix in compact form
!==========================================================================!

        fact1 = one / Re
        fact2 = gamma / (Pr * Re)

        Vij(1) = fact1 * rhoinv * ( lm + two * mu )
        Vij(2) = fact1 * rhoinv * mu
        Vij(3) = fact2 * rhoinv * con
        Vij(4) = fact1 * rhoinv * ( lm + mu )

!==========================================================================!
!.... Compute \hat{A}
!==========================================================================!

!.... A, B, Vxx, Vyy, Vxy terms

        Ah(1,1,i,j) = A(1,1) * m1(i,j) + B(1,1) * m2(i,j)
        Ah(1,2,i,j) = A(1,2) * m1(i,j)
        Ah(1,3,i,j) =                    B(1,3) * m2(i,j)
        Ah(1,4,i,j) = zero
        Ah(1,5,i,j) = zero

        Ah(2,1,i,j) = A(2,1) * m1(i,j)
        Ah(2,2,i,j) = A(2,2) * m1(i,j) + B(2,2) * m2(i,j) - &
                      Vij(1) * m11(i,j) - Vij(2) * m22(i,j)
        Ah(2,3,i,j) = A(2,3) * m1(i,j) + B(2,3) * m2(i,j) - &
                      Vij(4) * m12(i,j)
        Ah(2,4,i,j) = A(2,4) * m1(i,j)
        Ah(2,5,i,j) = A(2,5) * m1(i,j) + B(2,5) * m2(i,j)

        Ah(3,1,i,j) =                    B(3,1) * m2(i,j)
        Ah(3,2,i,j) = A(3,2) * m1(i,j) + B(3,2) * m2(i,j) - &
                      Vij(4) * m12(i,j)
        Ah(3,3,i,j) = A(3,3) * m1(i,j) + B(3,3) * m2(i,j) - &
                      Vij(2) * m11(i,j) - Vij(1) * m22(i,j)
        Ah(3,4,i,j) =                    B(3,4) * m2(i,j)
        Ah(3,5,i,j) = A(3,5) * m1(i,j) + B(3,5) * m2(i,j)

        Ah(4,1,i,j) = zero
        Ah(4,2,i,j) = A(4,2) * m1(i,j)
        Ah(4,3,i,j) =                    B(4,3) * m2(i,j)
        Ah(4,4,i,j) = A(4,4) * m1(i,j) + B(4,4) * m2(i,j) - &
                      Vij(2) * m11(i,j) - Vij(2) * m22(i,j)
        Ah(4,5,i,j) = A(4,5) * m1(i,j) + B(4,5) * m2(i,j)

        Ah(5,1,i,j) = A(5,1) * m1(i,j) + B(5,1) * m2(i,j)
        Ah(5,2,i,j) = A(5,2) * m1(i,j) + B(5,2) * m2(i,j)
        Ah(5,3,i,j) = A(5,3) * m1(i,j) + B(5,3) * m2(i,j)
        Ah(5,4,i,j) = A(5,4) * m1(i,j) + B(5,4) * m2(i,j)
        Ah(5,5,i,j) = A(5,5) * m1(i,j) + B(5,5) * m2(i,j) - &
                      Vij(3) * m11(i,j) - Vij(3) * m22(i,j)

!.... correct the matrix for the parabolized outflow

        if (idamp.eq.1 .and. linear.eq.0) then
          Ah(2,1,i,j) = damp(i,j) * Ah(2,1,i,j)
          Ah(2,5,i,j) = Ah(2,5,i,j) - (one-damp(i,j))*m1(i,j)/(gamma*Ma**2)

          Ah(3,1,i,j) = damp(i,j) * Ah(3,1,i,j)
          Ah(3,5,i,j) = Ah(3,5,i,j) - (one-damp(i,j))*m2(i,j)/(gamma*Ma**2)
        end if

!==========================================================================!
!.... Compute \hat{B}
!==========================================================================!

!.... A, B, Vxx, Vyy, Vxy term

        Bh(1,1,i,j) = A(1,1) * n1(i,j) + B(1,1) * n2(i,j)
        Bh(1,2,i,j) = A(1,2) * n1(i,j)
        Bh(1,3,i,j) =                    B(1,3) * n2(i,j)
        Bh(1,4,i,j) = zero
        Bh(1,5,i,j) = zero

        Bh(2,1,i,j) = A(2,1) * n1(i,j)
        Bh(2,2,i,j) = A(2,2) * n1(i,j) + B(2,2) * n2(i,j) - &
                      Vij(1) * n11(i,j) - Vij(2) * n22(i,j)
        Bh(2,3,i,j) = A(2,3) * n1(i,j) + B(2,3) * n2(i,j) - &
                      Vij(4) * n12(i,j)
        Bh(2,4,i,j) = A(2,4) * n1(i,j)
        Bh(2,5,i,j) = A(2,5) * n1(i,j) + B(2,5) * n2(i,j)

        Bh(3,1,i,j) =                    B(3,1) * n2(i,j)
        Bh(3,2,i,j) = A(3,2) * n1(i,j) + B(3,2) * n2(i,j) - &
                      Vij(4) * n12(i,j)
        Bh(3,3,i,j) = A(3,3) * n1(i,j) + B(3,3) * n2(i,j) - &
                      Vij(2) * n11(i,j) - Vij(1) * n22(i,j)
        Bh(3,4,i,j) =                    B(3,4) * n2(i,j)
        Bh(3,5,i,j) = A(3,5) * n1(i,j) + B(3,5) * n2(i,j)

        Bh(4,1,i,j) = zero
        Bh(4,2,i,j) = A(4,2) * n1(i,j)
        Bh(4,3,i,j) =                    B(4,3) * n2(i,j)
        Bh(4,4,i,j) = A(4,4) * n1(i,j) + B(4,4) * n2(i,j) - &
                      Vij(2) * n11(i,j) - Vij(2) * n22(i,j)
        Bh(4,5,i,j) = A(4,5) * n1(i,j) + B(4,5) * n2(i,j)

        Bh(5,1,i,j) = A(5,1) * n1(i,j) + B(5,1) * n2(i,j)
        Bh(5,2,i,j) = A(5,2) * n1(i,j) + B(5,2) * n2(i,j)
        Bh(5,3,i,j) = A(5,3) * n1(i,j) + B(5,3) * n2(i,j)
        Bh(5,4,i,j) = A(5,4) * n1(i,j) + B(5,4) * n2(i,j)
        Bh(5,5,i,j) = A(5,5) * n1(i,j) + B(5,5) * n2(i,j) - &
                      Vij(3) * n11(i,j) - Vij(3) * n22(i,j)

!.... correct for the adiabatic condition at the wall:  dT/dn = 0

        if ( (.not. yper) .and. wallt.eq.2 .and. j.eq.1) then
          Bh(:,5,i,j) = zero
        end if

!==========================================================================!
!.... compute the D matrix and \hat{D}
!==========================================================================!

!.... Continuity equation

        D(1,1) = divu
        D(1,2) = grho(1)
        D(1,3) = grho(2)
        D(1,4) = grho(3)
        D(1,5) = zero

        fact1 = rhoinv / Re
        gmsinv = one / (gamma * Ma**2)

!.... Momentum equation -- x_1

        if (linear.eq.1) then
          D(2,1) = rhoinv * ( u1 * gu(1,1) + u2 * gu(1,2) + &
                              u3 * gu(1,3) + gmsinv * gt(1) )
        else
          D(2,1) = rhoinv * ( gmsinv * gt(1) - rhoinv * gp(1) + &
                   fact1 * ( g1lm * divu + lm * g1divu ) + &
                   fact1 * two * ( g1mu * S(1,1) + &
                                   g2mu * S(1,2) + &
                                   g3mu * S(1,3) + &
                                     mu * S1jj ) )
        end if  

        D(2,2) = gu(1,1)
        D(2,3) = gu(1,2)
        D(2,4) = gu(1,3)
        D(2,5) = rhoinv * gmsinv * grho(1) -            &
                   fact1 * ( g1dlm * divu + dlm * g1divu ) -    &
                   fact1 * two * ( g1dmu * S(1,1) + &
                                   g2dmu * S(1,2) + &
                                   g3dmu * S(1,3) + &
                                     dmu * S1jj )

!.... Momentum equation -- x_2

        if (linear.eq.1) then
          D(3,1) = rhoinv * ( u1 * gu(2,1) + u2 * gu(2,2) + &
                              u3 * gu(2,3) + gmsinv * gt(2) )
        else
          D(3,1) = rhoinv * ( gmsinv * gt(2) - rhoinv * gp(2) + &
                     fact1 * ( g2lm * divu + lm * g2divu ) + &
                     fact1 * two * ( g1mu * S(2,1) + &
                                     g2mu * S(2,2) + &
                                     g3mu * S(2,3) + &
                                       mu * S2jj  ) )
        end if  

        D(3,2) = gu(2,1)
        D(3,3) = gu(2,2)
        D(3,4) = gu(2,3)
        D(3,5) = rhoinv * gmsinv * grho(2)  - &
                   fact1 * ( g2dlm * divu + dlm * g2divu ) - &
                   fact1 * two * ( g1dmu * S(2,1) + &
                                   g2dmu * S(2,2) + &
                                   g3dmu * S(2,3) + &
                                     dmu * S2jj )

!.... Momentum equation -- x_3

        if (linear.eq.1) then
          D(4,1) = rhoinv * ( u1 * gu(3,1) + u2 * gu(3,2) + &
                                u3 * gu(3,3) + gmsinv * gt(3) )
        else
          D(4,1) = rhoinv * ( gmsinv * gt(3) - rhoinv * gp(3) + &
                     fact1 * ( g3lm * divu + lm * g3divu ) + &
                     fact1 * two * ( g1mu * S(3,1) + &
                                     g2mu * S(3,2) + &
                                     g3mu * S(3,3) + &
                                       mu * S3jj ) )
        end if  

        D(4,2) = gu(3,1)
        D(4,3) = gu(3,2)
        D(4,4) = gu(3,3)
        D(4,5) = rhoinv * gmsinv * grho(3)  - &
                   fact1 * ( g3dlm * divu + dlm * g3divu ) - &
                   fact1 * two * ( g1dmu * S(3,1) + &
                                   g2dmu * S(3,2) + &
                                   g3dmu * S(3,3) + &
                                     dmu * S3jj )

!.... Energy equation

        fact1 = gamma * rhoinv / (Pr * Re)
        fact2 = gamma * gamma1 * Ma**2 * rhoinv / Re

        if (linear .eq. 1) then
          D(5,1) = rhoinv * ( u1 * gt(1) + u2 * gt(2) + &
                                u3 * gt(3) + gamma1 * t * divu )
        else
          D(5,1) = fact1 * rhoinv * ( g1con * gt(1) + &
                                        g2con * gt(2) + &
                                        g3con * gt(3) + &
                                          con * Lapt ) + &
                     fact2 * rhoinv * lm * divu**2 + &
                     two * fact2 * rhoinv * mu * ( &
                        S(1,1)**2 + S(1,2)**2 + S(1,3)**2 + &
                        S(2,1)**2 + S(2,2)**2 + S(2,3)**2 + &
                        S(3,1)**2 + S(3,2)**2 + S(3,3)**2 )

        end if
        D(5,2) = gt(1)
        D(5,3) = gt(2)
        D(5,4) = gt(3)

!.... to get the linear verion, it is assumed that \rho,t = 0
!.... Note that it really isn't necessary, it would be okay to always
!.... use the nonlinear version.

        if (linear .eq. 1) then
          D(5,5) = -gamma1 * rhoinv * ( u1 * grho(1) + &
                      u2 * grho(2) + u3 * grho(3) )
        else
          D(5,5) = gamma1 * divu
        end if
        
        D(5,5) = D(5,5) - &
                   fact1 * ( g1dcon * gt(1) + &
                             g2dcon * gt(2) + &
                             g3dcon * gt(3) + &
                               dcon * Lapt ) - &
                   fact2 * dlm * divu**2 - &
                   two * fact2 * dmu * ( &
                        S(1,1)**2 + S(1,2)**2 + S(1,3)**2 + &
                        S(2,1)**2 + S(2,2)**2 + S(2,3)**2 + &
                        S(3,1)**2 + S(3,2)**2 + S(3,3)**2 )

!.... Put into global data structure

        Dh(:,:,i,j) = D(:,:)

!==========================================================================!

!.... Vzz term

        if (kzl .ne. zero) then
          Dh(2,2,i,j) = Dh(2,2,i,j) + kzl**2 * Vij(2)
          Dh(3,3,i,j) = Dh(3,3,i,j) + kzl**2 * Vij(2)
          Dh(4,4,i,j) = Dh(4,4,i,j) + kzl**2 * Vij(1)
          Dh(5,5,i,j) = Dh(5,5,i,j) + kzl**2 * Vij(3)
        end if

!==========================================================================!
!.... Compute \hat{V}_{\xi\xi} in compact form
!
!       Vh1(1) = V(2,2)
!       Vh1(2) = V(3,3)
!       Vh1(3) = V(4,4)
!       Vh1(4) = V(5,5)
!       Vh1(5) = V(2,3)
!       Vh1(6) = V(3,2)
!
!       Note:  Vh1(5) ussually equals Vh1(6) except on a boundary where
!              viscous BC's are applied.
!
!==========================================================================!

!.... Vxx term and Vyy term

        Vh11(1,i,j) = Vij(1) * m1m1(i,j) + Vij(2) * m2m2(i,j)
        Vh11(2,i,j) = Vij(2) * m1m1(i,j) + Vij(1) * m2m2(i,j)
        Vh11(3,i,j) = Vij(2) * m1m1(i,j) + Vij(2) * m2m2(i,j)
        Vh11(4,i,j) = Vij(3) * m1m1(i,j) + Vij(3) * m2m2(i,j)

!.... Vxy term

        Vh11(5,i,j) = Vij(4) * m1m2(i,j)
        Vh11(6,i,j) = Vij(4) * m1m2(i,j)

!.... damp the second derivatives in the viscous terms

        if (idamp.eq.1) then
          Vh11(:,i,j) = damp(i,j) * Vh11(:,i,j)
        end if

!==========================================================================!
!.... Compute \hat{V}_{\xi\eta}
!==========================================================================!

!.... Vxx term and Vyy

        Vh12(1,i,j) = two * ( Vij(1) * m1n1(i,j) + Vij(2) * m2n2(i,j) )
        Vh12(2,i,j) = two * ( Vij(2) * m1n1(i,j) + Vij(1) * m2n2(i,j) )
        Vh12(3,i,j) = two * ( Vij(2) * m1n1(i,j) + Vij(2) * m2n2(i,j) )
        Vh12(4,i,j) = two * ( Vij(3) * m1n1(i,j) + Vij(3) * m2n2(i,j) )

!.... Vxy term

        Vh12(5,i,j) = Vij(4) * ( m1n2(i,j) + m2n1(i,j) )
        Vh12(6,i,j) = Vij(4) * ( m1n2(i,j) + m2n1(i,j) )

!==========================================================================!
!.... Compute \hat{V}_{\eta\eta}
!==========================================================================!

!.... Vxx term and Vyy term

        Vh22(1,i,j) = Vij(1) * n1n1(i,j) + Vij(2) * n2n2(i,j)
        Vh22(2,i,j) = Vij(2) * n1n1(i,j) + Vij(1) * n2n2(i,j)
        Vh22(3,i,j) = Vij(2) * n1n1(i,j) + Vij(2) * n2n2(i,j)
        Vh22(4,i,j) = Vij(3) * n1n1(i,j) + Vij(3) * n2n2(i,j)

!.... Vxy term

        Vh22(5,i,j) = Vij(4) * n1n2(i,j)
        Vh22(6,i,j) = Vij(4) * n1n2(i,j)

!==========================================================================!
!.... Compute 3D terms
!==========================================================================!

        if (complex_analysis) then

!.... \hat{A}_i and \hat{B}_i

          ABhi(1,i,j) = kzl * Vij(4) * m1(i,j) 
          ABhi(2,i,j) = kzl * Vij(4) * m2(i,j)
          ABhi(3,i,j) = kzl * Vij(4) * n1(i,j)
          ABhi(4,i,j) = kzl * Vij(4) * n2(i,j)

!==========================================================================!
!.... compute the C matrix and \hat{D}_i
!==========================================================================!

!.... Continuity equation

          C(1,1) = u3
          C(1,2) = zero
          C(1,3) = zero
          C(1,4) = rho
          C(1,5) = zero

!.... Momentum equation -- x_1

          fact1  = rhoinv / Re
          gmsinv = one / (gamma * Ma**2)

          C(2,1) = zero
          C(2,2) = u3 - fact1 * g3mu
          C(2,3) = zero
          C(2,4) = -fact1 * g1lm
          C(2,5) = -fact1 * dmu * two * S(1,3)

!.... Momentum equation -- x_2

          C(3,1) = zero
          C(3,2) = zero
          C(3,3) = u3 - fact1 * g3mu
          C(3,4) = -fact1 * g2lm
          C(3,5) = -fact1 * dmu * two * S(2,3)

!.... Momentum equation -- x_3

          C(4,1) = rhoinv * t * gmsinv
          C(4,2) = -fact1 * g1mu
          C(4,3) = -fact1 * g2mu
          C(4,4) = u3 - fact1 * g3lm - fact1 * two * g3mu
          C(4,5) = gmsinv - fact1 * dlm * divu - &
                   fact1 * dmu * two * S(3,3)

!.... Energy equation

          fact1 = gamma * rhoinv / (Pr * Re)
          fact2 = gamma * gamma1 * Ma**2 * rhoinv / Re

          C(5,1) = zero
          C(5,2) = -four * fact2 * mu * S(1,3)
          C(5,3) = -four * fact2 * mu * S(2,3)
          C(5,4) =  gamma1 * t - fact2 * two * lm * divu - &
                    four * fact2 * mu * S(3,3)
          C(5,5) =  u3 - fact1 * (g3con + dcon * gt(3))


          Dhi(:,:,i,j) = kzl * C

!.... put in unsteady term

          Dhi(1,1,i,j) = Dhi(1,1,i,j) - omega
          Dhi(2,2,i,j) = Dhi(2,2,i,j) - omega
          Dhi(3,3,i,j) = Dhi(3,3,i,j) - omega
          Dhi(4,4,i,j) = Dhi(4,4,i,j) - omega
          Dhi(5,5,i,j) = Dhi(5,5,i,j) - omega

!.... put in imaginary omega term for shift and invert

!!$          Dh(1,1,i,j) = Dh(1,1,i,j) + omega_r
!!$          Dh(2,2,i,j) = Dh(2,2,i,j) + omega_r
!!$          Dh(3,3,i,j) = Dh(3,3,i,j) + omega_r
!!$          Dh(4,4,i,j) = Dh(4,4,i,j) + omega_r
!!$          Dh(5,5,i,j) = Dh(5,5,i,j) + omega_r

        end if

        end do loop_i
        end do loop_j

!.... explicit smoother

        if (linear.eq.0 .and. eps_e.ne.zero) call smoother( rl, vl )

        return
        end
