!=============================================================================!
        subroutine potential
!  
!  Determine the potential flow for use in far-field boundary conditions
! 
!  WARNING:  No error checking is done on the input files!
! 
!!  Revised: 5-24-95
!=============================================================================!
        use global
        use pot
        implicit none
        
        integer :: i, j
!=============================================================================!
        if (linear.eq.1) return
        
        if (top.eq.0) then

          allocate( rhobt(nx), ubt(nx), vbt(nx), wbt(nx), &
                    tbt(nx), pbt(nx), cbt(nx) ) 
  
          open(10,file='top.pot',status='old',err=100)
          do i = 1, nx
            read(10,*,err=100) rhobt(i), ubt(i), vbt(i), wbt(i), tbt(i)
          end do
          close(10)

          pbt = rhobt * tbt / ( gamma * Ma**2 )
          cbt = sqrt(tbt) / Ma
        
        end if
        
        if (right.eq.0) then
        
          allocate( rhobr(ny), ubr(ny), vbr(ny), wbr(ny), &
                    tbr(ny), pbr(ny), cbr(ny) ) 
          
          open(10,file='right.pot',status='old',err=110)
          do j = 1, ny
            read(10,*,err=110) rhobr(j), ubr(j), vbr(j), wbr(j), tbr(j)
          end do
          close(10)

!         rhobr = one
!         ubr   = zero
!         vbr   = zero
!         wbr   = zero
!         tbr   = one

          pbr   = rhobr * tbr / ( gamma * Ma**2 )
          cbr   = sqrt(tbr) / Ma
        end if
        
        if (left.eq.0) then
        
          allocate( rhobl(ny), ubl(ny), vbl(ny), wbl(ny), &
                    tbl(ny), pbl(ny), cbl(ny) ) 
          
          open(10,file='left.pot',status='old',err=130)
          do j = 1, ny
            read(10,*,err=130) rhobl(j), ubl(j), vbl(j), wbl(j), tbl(j)
          end do
          close(10)

!         rhobl = one
!         ubl   = zero
!         vbl   = zero
!         wbl   = zero
!         tbl   = one

          pbl   = rhobl * tbl / ( gamma * Ma**2 )
          cbl   = sqrt(tbl) / Ma

        end if

        if (idamp.eq.1) then
          allocate( g1p(ny,nx) )
          open(10,file='pg.pot',form='unformatted',status='old',err=120)
          read(10,err=120) g1p
          close(10)
        end if
        
        return
100     call error('potential$','Error in top.pot$')
110     call error('potential$','Error in right.pot$')
120     call error('potential$','Error in pg.pot$')
130     call error('potential$','Error in left.pot$')
        end

