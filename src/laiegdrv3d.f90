!=============================================================================!
        subroutine eigdrv3d()
!  
!  Driver for eigenvalue problem (LAPack)
!
!=============================================================================!
        use global
        use local2d
        use local3d
        implicit none
!=============================================================================!
        complex :: v(ndof,nx,ny)
        complex :: matx(5,ndof,ndof,nx,ny)
        complex :: mat(ndof,nx,ndof,nx)
        real :: dtl(nx,ny)
        
!.... LAPack variables

        integer :: lwork, info
        real, allocatable :: rwork(:)
        complex, allocatable :: work(:)
        complex :: eval(ndof*nx), evec(ndof,nx,ndof*nx)
        
        integer :: i, j, idof, jdof
        character*80 :: name, code='EigDrv$'
!=============================================================================!

        dtl = one

        call lhs1f3D( matx, Ah, Dh, Dhi, Vh11, ABhi, dtl, v )

!.... Remove identity matrix

        do j = 1, ny
          do i = 1, nx
            do idof = 1, ndof
              matx(3,idof,idof,i,j) = matx(3,idof,idof,i,j) - one
            end do
          end do
        end do

!.... Convert to LAPack general matrix storage

        j = 1

        i = 1
        do idof = 1, ndof
          do jdof = 1, ndof
            mat(idof,i,jdof,nx-2) = matx(1,idof,jdof,i,j)
            mat(idof,i,jdof,nx-1) = matx(2,idof,jdof,i,j)
            mat(idof,i,jdof,i  )  = matx(3,idof,jdof,i,j)
            mat(idof,i,jdof,i+1)  = matx(4,idof,jdof,i,j)
            mat(idof,i,jdof,i+2)  = matx(5,idof,jdof,i,j)
          end do
        end do
        
        i = 2
        do idof = 1, ndof
          do jdof = 1, ndof
            mat(idof,i,jdof,nx-1) = matx(1,idof,jdof,i,j)
            mat(idof,i,jdof,i-1)  = matx(2,idof,jdof,i,j)
            mat(idof,i,jdof,i  )  = matx(3,idof,jdof,i,j)
            mat(idof,i,jdof,i+1)  = matx(4,idof,jdof,i,j)
            mat(idof,i,jdof,i+2)  = matx(5,idof,jdof,i,j)
          end do
        end do

        do i = 3, nx-2
          do idof = 1, ndof
            do jdof = 1, ndof
              mat(idof,i,jdof,i-2) = matx(1,idof,jdof,i,j)
              mat(idof,i,jdof,i-1) = matx(2,idof,jdof,i,j)
              mat(idof,i,jdof,i  ) = matx(3,idof,jdof,i,j)
              mat(idof,i,jdof,i+1) = matx(4,idof,jdof,i,j)
              mat(idof,i,jdof,i+2) = matx(5,idof,jdof,i,j)
            end do
          end do
        end do

        i = nx-1
        do idof = 1, ndof
          do jdof = 1, ndof
            mat(idof,i,jdof,i-2) = matx(1,idof,jdof,i,j)
            mat(idof,i,jdof,i-1) = matx(2,idof,jdof,i,j)
            mat(idof,i,jdof,i  ) = matx(3,idof,jdof,i,j)
            mat(idof,i,jdof,i+1) = matx(4,idof,jdof,i,j)
            mat(idof,i,jdof,2)   = matx(5,idof,jdof,i,j)
          end do
        end do
 
        i = nx
        do idof = 1, ndof
          do jdof = 1, ndof
            mat(idof,i,jdof,i-2) = matx(1,idof,jdof,i,j)
            mat(idof,i,jdof,i-1) = matx(2,idof,jdof,i,j)
            mat(idof,i,jdof,i  ) = matx(3,idof,jdof,i,j)
            mat(idof,i,jdof,2)   = matx(4,idof,jdof,i,j)
            mat(idof,i,jdof,3)   = matx(5,idof,jdof,i,j)
          end do
        end do

        lwork = max(1,3*(ndof*nx))
        allocate( work(lwork), rwork(2*(ndof*nx)) )
        call ZGEGV ( 'N', 'V', ndof*nx, mat, ndof*nx, eval, &
                     evec, ndof*nx, evec, ndof*nx,          &
                     work, lwork, rwork, info )
        deallocate( work, rwork )

!.... Output the eigenvalues and eigenvectors

        open(11,file='eig.dat', form='unformatted', status='unknown')
        write(11) ndof*nx, 0.0, nx, 1, nz, ndof, &
             Re, Ma, Pr, gamma, cv
        write(11) eval
        write(11) evec
        close(11)

        end subroutine eigdrv3d
