!=======================================================================================================!
        subroutine genB_l(B, vl)
!
!  Form the B matrix which multiplies U,y
!
!  Revised: 10-15-95
!
!  This version has been parabolized in x and implements the Lele & Poinsot
!  boundary conditions for the linear equations.
!=======================================================================================================!
        use stuff
        use local
        implicit none
        
        real :: B(ny*nx,ndof,ndof), vl(ny*nx,ndof)
        real :: fact1, fact2, fact3, fact4, gmsinv

        real :: pinf, cinf, Lx, rk, c
        real :: L1, L2, L3, L4, L5
        real :: d1, d2, d3, d4, d5, dc
        
        integer :: i, j, ji

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
            
!.... Continuity equation

            gmsinv = one / ( gamma * Ma**2 )
    
            B(ji,1,1) = one/c**2*(un-c)*bn2*t(ji)*gmsinv + us*bs2
            B(ji,1,2) = -one/c**2*(un-c)*rho(ji)*c*bn1*bn2 + rho(ji)*bs1*bs2
            B(ji,1,3) = -one/c**2*(un-c)*rho(ji)*c*bn2*bn2 + rho(ji)*bs2*bs2
            B(ji,1,4) = zero
            B(ji,1,5) = one/c**2*(un-c)*bn2*rho(ji)*gmsinv

        end do

        end if          ! wall

!=======================================================================================================!
!       S l i p   W a l l   B o u n d a r y
!=======================================================================================================!
        if (left.eq.2) then

        i = 1
        do j = 1, ny
          ji = j + (i-1) * ny
    
!.... get the boundary normal and tangent vectors

            bn1 = m1(ji) / sqrt( m1(ji)**2 + m2(ji)**2 )
            bn2 = m2(ji) / sqrt( m1(ji)**2 + m2(ji)**2 )

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
            
!.... Continuity equation

            gmsinv = one / ( gamma * Ma**2 )
    
            B(ji,1,1) = one/c**2*(un-c)*bn2*t(ji)*gmsinv + us*bs2
            B(ji,1,2) = -one/c**2*(un-c)*rho(ji)*c*bn1*bn2 + rho(ji)*bs1*bs2
            B(ji,1,3) = -one/c**2*(un-c)*rho(ji)*c*bn2*bn2 + rho(ji)*bs2*bs2
            B(ji,1,4) = zero
            B(ji,1,5) = one/c**2*(un-c)*bn2*rho(ji)*gmsinv

        end do

        end if        ! slip wall

        return
        end
