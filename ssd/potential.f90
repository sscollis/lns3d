!=============================================================================!
        subroutine potential
!  
!  Determine the potential flow for use in far-field boundary conditions
! 
!  WARNING:  No error checking is done on the input files!
! 
!!  Revised: 5-24-95
!=============================================================================!
        use stuff
        use pot
        implicit none
        
        integer :: i, j
!=============================================================================!
        if (linear.eq.1) return
        
        if (top.eq.0) then

          allocate( rhobt(nx), ubt(nx), vbt(nx), wbt(nx), &
                    tbt(nx), pbt(nx), cbt(nx) ) 
  
          open(10,file='top.pot')
          do i = 1, nx
            read(10,*) rhobt(i), ubt(i), vbt(i), wbt(i), tbt(i)
          end do
          close(10)

!.... HACK: SSC

!         write(*,*) "****************  WARNING ******************"
!         write(*,*) "****************  WARNING ******************"
!         write(*,*) "****************  WARNING ******************"
!         write(*,*) "****************  WARNING ******************"
!         write(*,*) "****************  WARNING ******************"
!         write(*,*) "**** Hacked in uniform potential flow in potential.f90"
!         write(*,*)
!
!         rhobt(:) = one
!         ubt(:)   = one
!         vbt(:)   = zero
!         wbt(:)   = zero
!         tbt(:)   = one

          pbt = rhobt * tbt / ( gamma * Ma**2 )
          cbt = sqrt(tbt) / Ma
        
        end if
        
        if (right.eq.0) then
        
          allocate( rhobr(ny), ubr(ny), vbr(ny), wbr(ny), &
                    tbr(ny), pbr(ny), cbr(ny) ) 
          
          open(10,file='right.pot')
          do j = 1, ny
            read(10,*) rhobr(j), ubr(j), vbr(j), wbr(j), tbr(j)
          end do
          close(10)
          pbr = rhobr * tbr / ( gamma * Ma**2 )
          cbr = sqrt(tbr) / Ma

        end if
        
        if (idamp.eq.1) then
          allocate( g1p(ny,nx) )
          open(10,file='pg.pot',form='unformatted')
          read(10) g1p
          close(10)
        end if
        
        return
        end

