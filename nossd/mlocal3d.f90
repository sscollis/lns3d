!=============================================================================!
        subroutine mlocal3d
!
!  allocate space for 3d local variables
!
!=============================================================================!
        use global
        use local3d
        implicit none
        
        integer :: ier, mem

        character*10 :: opt
        character*80 :: code='mLocal3D$'
!=============================================================================!
        mem = 0

        mem = mem + 2*ny*nx*ndof*5
        allocate( c1v(ny*nx,ndof),  c2v(ny*nx,ndof),            &
                  c11v(ny*nx,ndof), c12v(ny*nx,ndof),           &
                  c22v(ny*nx,ndof), STAT=ier )
        if (ier .ne. 0) call error(code,'Insufficient Memory for gradients$')
          
        write(*,"(' Mlocal3D allocated ===> ',1pe13.6,' words')") float(mem)

        return
        end
