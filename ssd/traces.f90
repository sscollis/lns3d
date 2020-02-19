!=============================================================================!
        subroutine traces(vl)
!=============================================================================!
!  
!  This subroutine outputs time traces
!
!=============================================================================!
        use stuff
        use material
        implicit none

        real :: vl(ny,nx,ndof)

        integer :: i, j
        real :: rho, u1, u2, u3, p, t

        logical :: mflag = .false.
        integer :: itrace=50
        integer, external :: igetver
!=============================================================================!

!.... open the time-trace file

        if (istep.eq.1) then
          iver = igetver( base, 't'//char(0) ) + 1
          call addver(base, 't'//char(0), iver, filen, lfile)
          open(unit=itrace, file=filen, form='formatted', status='unknown')
        end if

!.... output a time trace at point (i,j)

        i = 160
        j = 30

        rho = vl(j,i,1)
        u1  = vl(j,i,2)
        u2  = vl(j,i,3)
        u3  = vl(j,i,4)
        t   = vl(j,i,5)
        p   = rho * t / (gamma * Ma**2)
        
        write(itrace,10) time, rho, u1, u2, p, t

!.... compute the mass flux through the boundaries

        if (mflag) call massflx(vl)

        return
10      format(6(1pe13.6,1x))
        end

!=============================================================================!
        subroutine massflx(vl)
!=============================================================================!
!  
!  This computes the mass flux crossing the computational boundaries
!
!=============================================================================!
        use stuff
        use material
        implicit none

        real :: vl(ny,nx,ndof)

        integer :: i, j, ij, ip1j, ijp1

        real, allocatable :: n1l(:), n2l(:), vn(:)
        real :: mass_top, mass_bot, mass_left, mass_right, mass_total, ds
!=============================================================================!

!.... integrate the mass flux over the top boundary
!.... This integration is not very good

        allocate( n1l(nx), n2l(nx), vn(nx) )

!.... get the metrics along this boundary
  
        j = ny
        do i = 1, nx
          ij = j + (i-1) * ny
          n1l(i) = n1(ij) / sqrt( n1(ij)**2 + n2(ij)**2 )
          n2l(i) = n2(ij) / sqrt( n1(ij)**2 + n2(ij)**2 )
        end do

!.... compute Vn and Vt
        
        vn(:) = n1l(:) * vl(j,:,2) + n2l(:) * vl(j,:,3)

        mass_top  = zero

        do i = 1, nx-1
          ij   = j + (i-1) * ny
          ip1j = j + i * ny
    
          mass_top = mass_top + pt5 * ( vl(j,i,1)   * vn(i)   / &
                     sqrt(m1(ij)**2   + m2(ij)**2  ) + &
                     vl(j,i+1,1) * vn(i+1) / &
                     sqrt(m1(ip1j)**2 + m2(ip1j)**2) ) * dxi
        end do

!.... get the metrics along this boundary
  
        j = 1
        do i = 1, nx
          ij = j + (i-1) * ny
          n1l(i) = -n1(ij) / sqrt( n1(ij)**2 + n2(ij)**2 )
          n2l(i) = -n2(ij) / sqrt( n1(ij)**2 + n2(ij)**2 )
        end do

!.... compute Vn and Vt
        
        vn(:) = n1l(:) * vl(j,:,2) + n2l(:) * vl(j,:,3)
      
        mass_bot  = zero
      
        do i = 1, nx-1
          ij   = j + (i-1) * ny
          ip1j = j + i * ny
    
          mass_bot = mass_bot + pt5 * ( vl(j,i,1)   * vn(i)   / &
                      sqrt(m1(ij)**2   + m2(ij)**2  ) + &
                      vl(j,i+1,1) * vn(i+1) / &
                      sqrt(m1(ip1j)**2 + m2(ip1j)**2) ) * dxi
        end do

        deallocate( n1l, n2l, vn )

!.... integrate the mass flux over the left boundary

        allocate( n1l(ny), n2l(ny), vn(ny) )

!.... get the metrics along this boundary
  
        i = 1
        do j = 1, ny
          ij = j + (i-1) * ny
          n1l(j)  = -m1(ij) / sqrt( m1(ij)**2 + m2(ij)**2 )
          n2l(j)  = -m2(ij) / sqrt( m1(ij)**2 + m2(ij)**2 )
        end do

!.... compute Vn and Vt
        
        vn(:) = n1l(:) * vl(:,i,2) + n2l(:) * vl(:,i,3)
      
        mass_left = zero
      
        do j = 1, ny-1
            ij   = j + (i-1) * ny
            ijp1 = j+1 + (i-1) * ny
      
            mass_left = mass_left + pt5 * ( vl(j,i,1)   * vn(j)   / &
                        sqrt(n1(ij)**2   + n2(ij)**2  ) + &
                        vl(j+1,i,1) * vn(j+1) / &
                        sqrt(n1(ijp1)**2 + n2(ijp1)**2) ) * deta
        end do

!.... get the metrics along this boundary
  
        i = nx
        do j = 1, ny
          ij = j + (i-1) * ny
          n1l(j)  = m1(ij) / sqrt( m1(ij)**2 + m2(ij)**2 )
          n2l(j)  = m2(ij) / sqrt( m1(ij)**2 + m2(ij)**2 )
        end do

!.... compute Vn and Vt
        
        vn(:) = n1l(:) * vl(:,i,2) + n2l(:) * vl(:,i,3)
      
        mass_right = zero
      
        do j = 1, ny-1
          ij   = j + (i-1) * ny
          ijp1 = j+1 + (i-1) * ny
    
          mass_right = mass_right + pt5 * ( vl(j,i,1)   * vn(j)   / &
                       sqrt(n1(ij)**2   + n2(ij)**2  ) + &
                       vl(j+1,i,1) * vn(j+1) / &
                       sqrt(n1(ijp1)**2 + n2(ijp1)**2) ) * deta
        end do

        deallocate( n1l, n2l, vn )
      
        mass_total = mass_top + mass_bot + mass_left + mass_right
      
        write(51,10) time, mass_top, mass_bot, mass_left, mass_right, mass_total

        return
10      format(6(1pe13.6,1x))
        end
