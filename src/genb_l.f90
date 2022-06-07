!=======================================================================================================!
        subroutine genB_l(B, vl, i, j)
!
!  Form the B matrix which multiplies U,y
!
!  Revised: 6-6-22 
!
!  This version has been parabolized in x and implements the Lele & Poinsot
!  boundary conditions for the linear equations.
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

!.... compute some extra stuff needed for Lele's BC's

        pinf = one / (gamma * Ma**2)    ! pressure at infinity
        cinf = one / Ma                 ! speed of sound at infinity
        Lx   = 100.0                    ! length of domain in x
        
        rk = sigma * ( one - Ma**2 ) * cinf / Lx

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

!=======================================================================================================!
!       S l i p   W a l l   B o u n d a r y
!=======================================================================================================!
        if (left.eq.2 .and. i.eq.1) then
    
!.... get the boundary normal and tangent vectors

            bn1 = m1(i,j) / sqrt( m1(i,j)**2 + m2(i,j)**2 )
            bn2 = m2(i,j) / sqrt( m1(i,j)**2 + m2(i,j)**2 )

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

        end if        ! slip wall

        return
        end
