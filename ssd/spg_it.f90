!=============================================================================!
        module ic

        integer :: ic_start=0
        real, allocatable :: vic(:,:,:)
        
        end module ic

!=============================================================================!
        subroutine spg_it( rl, vl, spgl )
!
!       This routine reads in the initial condition when first called, and
!       sponges to the initial condition.
!
!=============================================================================!
        use ic
        use stuff
        use local
        implicit none
        
        real :: rl(ny,nx,ndof), vl(ny,nx,ndof), spgl(ny,nx)
        integer :: i

!.... local variables

        real :: etal(ny), xil(ny), u1l(ny), u2l(ny), tl(ny), rhol(ny)
        
        real :: rtmp
        integer :: itmp
!=============================================================================!
        if (ic_start .eq. 0) then
          allocate( vic(ny,nx,ndof) )
          ic_start = 1
          open(10,file='output.R.0',form='unformatted',status='old')
          read(10) itmp, rtmp, itmp, itmp, itmp, itmp, &
                   rtmp, rtmp, rtmp, rtmp, rtmp
          read(10) vic
          close(10)
        end if

        do i = 1, nx
          rl(:,i,1) = rl(:,i,1) + spgl(:,i) * ( vl(:,i,1) - vic(:,i,1) )
          rl(:,i,2) = rl(:,i,2) + spgl(:,i) * ( vl(:,i,2) - vic(:,i,2) )
          rl(:,i,3) = rl(:,i,3) + spgl(:,i) * ( vl(:,i,3) - vic(:,i,3) )
          rl(:,i,4) = rl(:,i,4) + spgl(:,i) * ( vl(:,i,4) - vic(:,i,4) )
          rl(:,i,5) = rl(:,i,5) + spgl(:,i) * ( vl(:,i,5) - vic(:,i,5) )
        end do
        
        return
        end
