!=============================================================================!
        subroutine Reimann( n, vl, vl1, vl2, bn, xl, yl, &
                            rhob, ub, vb, wb, tb, pb, cb ) 
!=============================================================================!
!  
!  Enforce Reimann boundary condition:
!
!  Inputs:
!
!  n:         number of nodes on the boundary
!  vl:        primative variables on the boundary
!  vl1:       primative variables one node in from the boundary
!  vl2:       primative variables two nodes in from the boundary
!  bn:        unit boundary normal
!  xl:        x-coordinates on the boundary
!  yl:        y-coordinates on the boundary
!  rhob:      rho from the potential solution on the boundary
!  ub:        u from the potential solution on the boundary
!  vb:        v from the potential solution on the boundary
!  wb:        w from the potential solution on the boundary
!  tb:        t from the potential solution on the boundary
!  pb:        p from the potential solution on the boundary
!  cb:        c from the potential solution on the boundary
!
!  Outputs:
!
!  vl:        Corrected primative variables on the boundary
!
!  Notes:
!
!  1. I extrapolate primative variables and then form the outgoing
!     Reimann.  The other option is to directly extrapolate the outgoing 
!     Reimann but I haven't tried this.  (It might work better? who knows...)
!
!  2. This could be simplified since the boundary normal vectors are now
!     of unit length
!
!  3. You could also simplify the density term since, 
!     rho^\gamma / p = \gamma \M^2
!
!  Revised: 7-31-96
!=============================================================================!
        use global
        implicit none

        integer :: n
        real :: vl(ndof,n), vl1(ndof,n), vl2(ndof,n)
        real :: bn(n,2), xl(n), yl(n)

        real :: vn(n), vt(n), vint(ndof,n)
        real :: etal(n), xil(n)
        real :: ub(n), vb(n), wb(n), tb(n), rhob(n), pb(n), cb(n)

        real :: R1(n), R2(n), c(n), rho(n), p(n)
        integer :: i

!.... variables for smoothing

        real :: fact
        real :: rho_in, p_in, vn_in, vt_in, v_in(ndof)
        real :: rho_out, p_out, vn_out, vt_out, v_out(ndof)
!=============================================================================!

!.... compute Vn and Vt
        
        vn(:) = bn(:,1) * vl(2,:) + bn(:,2) * vl(3,:)
        vt(:) = bn(:,2) * vl(2,:) - bn(:,1) * vl(3,:)

!.... perform extrapolation

        vint(:,:) = two * vl1(:,:) - vl2(:,:)

!.... set the Riemann invariants

        R1(:)  = bn(:,1) * ub(:) + bn(:,2) * vb(:) - two * cb / gamma1
        c(:)   = sqrt( vint(5,:) ) / Ma
        R2(:)  = bn(:,1) * vint(2,:) + bn(:,2) * vint(3,:) + &
                 two * c(:) / gamma1

        do i = 1, n

!!$          rho_in = rhob(i)
!!$          p_in   = pb(i)
!!$
!!$          vn_in  = pt5 * ( R1(i) + R2(i) )
!!$          vt_in  = bn(i,2) * ub(i) - bn(i,1) * vb(i)
!!$
!!$          v_in(1) = ( rho_in**gamma / (gamma * p_in) * &
!!$               (gamma1 / four)**2 * &
!!$               ( R2(i) - R1(i) )**2 ) ** (one/gamma1)
!!$          v_in(2) = (-bn(i,1) * vn_in - bn(i,2) * vt_in) / &
!!$               (-bn(i,1)**2 - bn(i,2)**2)
!!$          v_in(3) = (-bn(i,2) * vn_in + bn(i,1) * vt_in) / &
!!$               (-bn(i,1)**2 - bn(i,2)**2)
!!$          v_in(4) = tan(theta)
!!$          v_in(5) = ( Ma * gamma1 / four * ( R2(i) - R1(i) ) )**2
!!$          
!!$          rho_out = vint(1,i)
!!$          p_out   = vint(1,i) * vint(5,i) / (gamma * Ma**2)
!!$          
!!$          vn_out = pt5 * ( R1(i) + R2(i) )
!!$          vt_out = bn(i,2) * vint(2,i) - bn(i,1) * vint(3,i)
!!$
!!$          v_out(1) = ( rho_out**gamma / (gamma * p_out) * &
!!$               (gamma1 / four)**2 * &
!!$               ( R2(i) - R1(i) )**2 ) ** (one/gamma1)
!!$          v_out(2) = (-bn(i,1) * vn_out - bn(i,2) * vt_out) / &
!!$               (-bn(i,1)**2 - bn(i,2)**2)
!!$          v_out(3) = (-bn(i,2) * vn_out + bn(i,1) * vt_out) / &
!!$               (-bn(i,1)**2 - bn(i,2)**2)
!!$          v_out(4) = vint(4,i)
!!$          v_out(5) = ( Ma * gamma1 / four * ( R2(i) - R1(i) ) )**2
!!$          
!!$          fact = pt5 * ( one - tanh( 1000.0 * vn(i) ) )
!!$
!!$          vl(:,i) = fact * v_in(:) + (one - fact) * v_out(:)
!!$
          if ( vn(i) .le. zero ) then           ! inflow

            rho(i) = rhob(i)
            p(i)   = pb(i)

            vn(i)  = pt5 * ( R1(i) + R2(i) )
            vt(i)  = bn(i,2) * ub(i) - bn(i,1) * vb(i)

            vl(1,i) = ( rho(i)**gamma / (gamma * p(i)) * &
              (gamma1 / four)**2 * &
              ( R2(i) - R1(i) )**2 ) ** (one/gamma1)
            vl(2,i) = (-bn(i,1) * vn(i) - bn(i,2) * vt(i)) / &
              (-bn(i,1)**2 - bn(i,2)**2)
            vl(3,i) = (-bn(i,2) * vn(i) + bn(i,1) * vt(i)) / &
              (-bn(i,1)**2 - bn(i,2)**2)
            vl(4,i) = tan(theta)
            vl(5,i) = ( Ma * gamma1 / four * ( R2(i) - R1(i) ) )**2

          else                                  ! outflow

            rho(i) = vint(1,i)
            p(i)   = vint(1,i) * vint(5,i) / (gamma * Ma**2)

            vn(i)  = pt5 * ( R1(i) + R2(i) )
            vt(i)  = bn(i,2) * vint(2,i) - bn(i,1) * vint(3,i)

            vl(1,i) = ( rho(i)**gamma / (gamma * p(i)) * &
              (gamma1 / four)**2 * &
              ( R2(i) - R1(i) )**2 ) ** (one/gamma1)
            vl(2,i) = (-bn(i,1) * vn(i) - bn(i,2) * vt(i)) / &
              (-bn(i,1)**2 - bn(i,2)**2)
            vl(3,i) = (-bn(i,2) * vn(i) + bn(i,1) * vt(i)) / &
              (-bn(i,1)**2 - bn(i,2)**2)
            vl(4,i) = vint(4,i)
            vl(5,i) = ( Ma * gamma1 / four * ( R2(i) - R1(i) ) )**2

          end if

        end do

        return
        end

!=============================================================================!
        subroutine ReimannRHS( n, vl, vl1, vl2, bn, xl, yl, rl, &
                               rhob, ub, vb, wb, tb, pb, cb ) 
!  
!  Enforce Reimann boundary condition on the RHS
!  
!=============================================================================!
        use global
        implicit none

        integer :: n
        real :: vl(ndof,n), vl1(ndof,n), vl2(ndof,n), rl(ndof,n)
        real :: bn(n,2), xl(n), yl(n)

        real :: vn(n), vt(n), vint(ndof,n)
        real :: etal(n), xil(n)
        real :: ub(n), vb(n), wb(n), tb(n), rhob(n), pb(n), cb(n)

        real :: R1(n), R2(n), c(n), rho(n), p(n)
        integer :: i
!=============================================================================!

!.... compute Vn and Vt
        
        vn(:) = bn(:,1) * vl(2,:) + bn(:,2) * vl(3,:)
        vt(:) = bn(:,2) * vl(2,:) - bn(:,1) * vl(3,:)

!.... perform extrapolation

        vint(:,:) = two * vl1(:,:) - vl2(:,:)

!.... set the Riemann invariants

        R1(:)  = bn(:,1) * ub(:) + bn(:,2) * vb(:) - two * cb / gamma1
        c(:)   = sqrt( vint(5,:) ) / Ma
        R2(:)  = bn(:,1) * vint(2,:) + bn(:,2) * vint(3,:) + &
                 two * c(:) / gamma1

        do i = 1, n

          if ( vn(i) .le. zero ) then           ! inflow

            rho(i) = rhob(i)
            p(i)   = pb(i)

            vn(i)  = pt5 * ( R1(i) + R2(i) )
            vt(i)  = bn(i,2) * ub(i) - bn(i,1) * vb(i)

            rl(1,i) = vl(1,i) - ( rho(i)**gamma / (gamma * p(i)) * &
              (gamma1 / four)**2 * &
              ( R2(i) - R1(i) )**2 ) ** (one/gamma1)
            rl(2,i) = vl(2,i) - (-bn(i,1) * vn(i) - bn(i,2) * vt(i)) / &
              (-bn(i,1)**2 - bn(i,2)**2)
            rl(3,i) = vl(3,i) - (-bn(i,2) * vn(i) + bn(i,1) * vt(i)) / &
              (-bn(i,1)**2 - bn(i,2)**2)
            rl(4,i) = vl(4,i) - tan(theta)
            rl(5,i) = vl(5,i) - ( Ma * gamma1 / four * ( R2(i) - R1(i) ) )**2

          else                                  ! outflow

            rho(i) = vint(1,i)
            p(i)   = vint(1,i) * vint(5,i) / (gamma * Ma**2)

            vn(i)  = pt5 * ( R1(i) + R2(i) )
            vt(i)  = bn(i,2) * vint(2,i) - bn(i,1) * vint(3,i)

            rl(1,i) = vl(1,i) - ( rho(i)**gamma / (gamma * p(i)) * &
              (gamma1 / four)**2 * &
              ( R2(i) - R1(i) )**2 ) ** (one/gamma1)
            rl(2,i) = vl(2,i) - (-bn(i,1) * vn(i) - bn(i,2) * vt(i)) / &
              (-bn(i,1)**2 - bn(i,2)**2)
            rl(3,i) = vl(3,i) - (-bn(i,2) * vn(i) + bn(i,1) * vt(i)) / &
              (-bn(i,1)**2 - bn(i,2)**2)
            rl(4,i) = vl(4,i) - vint(4,i)
            rl(5,i) = vl(5,i) - ( Ma * gamma1 / four * ( R2(i) - R1(i) ) )**2

          end if

        end do

        return
        end

!=============================================================================!
        subroutine ReimannLHS( n, vl, vl1, vl2, bn, xl, yl, &
                               mat1, mat2, mat3, &
                               rhob, ub, vb, wb, tb, pb, cb )
!  
!  Enforce Reimann boundary condition on the LHS
!  
!=============================================================================!
        use global
        implicit none

        integer :: n
        real :: vl(ndof,n), vl1(ndof,n), vl2(ndof,n)
        real :: bn(n,2), xl(n), yl(n)
        real :: mat1(ndof,ndof,n), mat2(ndof,ndof,n), mat3(ndof,ndof,n)

        real :: vn(n), vt(n), vint(ndof,n)
        real :: etal(n), xil(n)
        real :: ub(n), vb(n), wb(n), tb(n), rhob(n), pb(n), cb(n)

        real :: R1(n), R2(n), c(n), rho(n), p(n)
        real :: fact(n), fact1(n), fact2(n)
        integer :: i
!=============================================================================!

!.... compute Vn and Vt
        
        vn(:) = bn(:,1) * vl(2,:) + bn(:,2) * vl(3,:)
        vt(:) = bn(:,2) * vl(2,:) - bn(:,1) * vl(3,:)

!.... perform extrapolation of primative variables

        vint(:,:) = two * vl1(:,:) - vl2(:,:)

!.... set the Riemann invariants

        R1(:)  = bn(:,1) * ub(:) + bn(:,2) * vb(:) - two * cb / gamma1
        c(:)   = sqrt( vint(5,:) ) / Ma
        R2(:)  = bn(:,1) * vint(2,:) + bn(:,2) * vint(3,:) + &
                 two * c(:) / gamma1

        do i = 1, n

          if ( vn(i) .le. zero ) then           ! inflow

          rho(i) = rhob(i)
          p(i)   = pb(i)

          vn(i)  = pt5 * ( R1(i) + R2(i) )
          vt(i)  = bn(i,2) * ub(i) - bn(i,1) * vb(i)

          fact(i) = (rho(i)**gamma/(gamma*p(i))* &
                    (gamma1/four)**2)**(one/gamma1) * &
                    two/gamma1 * (R2(i) - R1(i))**((two-gamma1)/gamma1)
          mat1(1,1,i) = -one
          mat2(1,2,i) = two * (fact(i) * bn(i,1))
          mat2(1,3,i) = two * (fact(i) * bn(i,2))
          mat2(1,5,i) = two * (fact(i) / ( gamma1 * Ma * sqrt(vint(5,i)) ))
          mat3(1,2,i) = -(fact(i) * bn(i,1))
          mat3(1,3,i) = -(fact(i) * bn(i,2))
          mat3(1,5,i) = -(fact(i) / ( gamma1 * Ma * sqrt(vint(5,i)) ))
          
          fact(i) = pt5 * bn(i,1) / (bn(i,1)**2 + bn(i,2)**2)
          mat1(2,2,i) = -one
          mat2(2,2,i) = two * (fact(i) * bn(i,1))
          mat2(2,3,i) = two * (fact(i) * bn(i,2))
          mat2(2,5,i) = two * (fact(i) / ( gamma1 * Ma * sqrt(vint(5,i)) ))
          mat3(2,2,i) = -(fact(i) * bn(i,1))
          mat3(2,3,i) = -(fact(i) * bn(i,2))
          mat3(2,5,i) = -(fact(i) / ( gamma1 * Ma * sqrt(vint(5,i)) ))
  
          fact(i) = pt5 * bn(i,2) / (bn(i,1)**2 + bn(i,2)**2)
          mat1(3,3,i) = -one
          mat2(3,2,i) = two * (fact(i) * bn(i,1))
          mat2(3,3,i) = two * (fact(i) * bn(i,2))
          mat2(3,5,i) = two * (fact(i) / ( gamma1 * Ma * sqrt(vint(5,i)) ))
          mat3(3,2,i) = -(fact(i) * bn(i,1))
          mat3(3,3,i) = -(fact(i) * bn(i,2))
          mat3(3,5,i) = -(fact(i) / ( gamma1 * Ma * sqrt(vint(5,i)) ))
  
          mat1(4,4,i) = -one
  
          fact(i) = (Ma * gamma1 / four)**2 * two * ( R2(i) - R1(i) )
          mat1(5,5,i) = -one
          mat2(5,2,i) = two * (fact(i) * bn(i,1))
          mat2(5,3,i) = two * (fact(i) * bn(i,2))
          mat2(5,5,i) = two * (fact(i) / ( gamma1 * Ma * sqrt(vint(5,i)) ))
          mat3(5,2,i) = -(fact(i) * bn(i,1))
          mat3(5,3,i) = -(fact(i) * bn(i,2))
          mat3(5,5,i) = -(fact(i) / ( gamma1 * Ma * sqrt(vint(5,i)) ))

          else                                  ! outflow

          rho(i) = vint(1,i)
          p(i)   = vint(1,i) * vint(5,i) / (gamma * Ma**2)

          vn(i)  = pt5 * ( R1(i) + R2(i) )
          vt(i)  = bn(i,2) * vint(2,i) - bn(i,1) * vint(3,i)

          fact1(i) = (rho(i)**gamma/(gamma*p(i)))**(one/gamma1) * &
                     (gamma1/four)**(two/gamma1) * (two/gamma1) * &
                     (R2(i) - R1(i))**((two-gamma1)/gamma1)
          fact2(i) = (gamma1/four * (R2(i) - R1(i)))**(two/gamma1) * &
                     one/(gamma*gamma1) * &
                     (rho(i)**gamma/(gamma*p(i)))**((two-gamma)/gamma1)
          mat1(1,1,i) = -one
          mat2(1,2,i) = two*(fact1(i) * bn(i,1) + &
                        gamma1 * rho(i)**gamma1 / p(i))
          mat2(1,3,i) = two*(fact1(i) * bn(i,2))
          mat2(1,5,i) = two*(fact1(i) / ( gamma1 * Ma * &
                        sqrt(vint(5,i)) ) - &
                        fact2(i) * rho(i)**(gamma+1)/(gamma*Ma**2*p(i)**2))
          mat3(1,2,i) = -(fact1(i) * bn(i,1) + &
                        gamma1 * rho(i)**gamma1 / p(i))
          mat3(1,3,i) = -(fact1(i) * bn(i,2))
          mat3(1,5,i) = -(fact1(i) / ( gamma1 * Ma * &
                        sqrt(vint(5,i)) ) - &
                        fact2(i) * rho(i)**(gamma+1)/(gamma*Ma**2*p(i)**2))
          
          fact1(i) = pt5 * bn(i,1) / (bn(i,1)**2 + bn(i,2)**2)
          fact2(i) =       bn(i,2) / (bn(i,1)**2 + bn(i,2)**2)
          mat1(2,2,i) = -one
          mat2(2,2,i) = two*(fact1(i) * bn(i,1) + fact2(i) * bn(i,2))
          mat2(2,3,i) = two*(fact1(i) * bn(i,2) - fact2(i) * bn(i,1))
          mat2(2,5,i) = two*(fact1(i) / ( gamma1 * Ma * sqrt(vint(5,i)) ))
          mat3(2,2,i) = -(fact1(i) * bn(i,1) + fact2(i) * bn(i,2))
          mat3(2,3,i) = -(fact1(i) * bn(i,2) - fact2(i) * bn(i,1))
          mat3(2,5,i) = -(fact1(i) / ( gamma1 * Ma * sqrt(vint(5,i)) ))
  
          fact1(i) = pt5 * bn(i,2) / (bn(i,1)**2 + bn(i,2)**2)
          fact2(i) =       bn(i,1) / (bn(i,1)**2 + bn(i,2)**2)
          mat1(3,3,i) = -one
          mat2(3,2,i) = two*(fact1(i) * bn(i,1) - fact2(i) * bn(i,2))
          mat2(3,3,i) = two*(fact1(i) * bn(i,2) + fact2(i) * bn(i,1))
          mat2(3,5,i) = two*(fact1(i) / ( gamma1 * Ma * sqrt(vint(5,i)) ))
          mat3(3,2,i) = -(fact1(i) * bn(i,1) - fact2(i) * bn(i,2))
          mat3(3,3,i) = -(fact1(i) * bn(i,2) + fact2(i) * bn(i,1))
          mat3(3,5,i) = -(fact1(i) / ( gamma1 * Ma * sqrt(vint(5,i)) ))
  
          mat1(4,4,i) = -one
          mat2(4,4,i) = two
          mat3(4,4,i) = -one
  
          fact(i) = (Ma * gamma1 / four)**2 * two * ( R2(i) - R1(i) )
          mat1(5,5,i) = -one
          mat2(5,2,i) = two*(fact(i) * bn(i,1))
          mat2(5,3,i) = two*(fact(i) * bn(i,2))
          mat2(5,5,i) = two*(fact(i) / ( gamma1 * Ma * sqrt(vint(5,i)) ))
          mat3(5,2,i) = -(fact(i) * bn(i,1))
          mat3(5,3,i) = -(fact(i) * bn(i,2))
          mat3(5,5,i) = -(fact(i) / ( gamma1 * Ma * sqrt(vint(5,i)) ))

          end if

        end do

        return
        end
