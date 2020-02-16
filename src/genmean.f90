!=============================================================================!
        subroutine genmean 
!  
!  Allocate memory and generate the mean flow 
!  
!=============================================================================!
        use global
        use material
        implicit none

        integer :: ier, mem=0
        
        character*80 name
        
        integer :: mstep, mnx, mny, mnz, mndof
        real    :: mtime, mMa, mRe, mPr, mgamma, mcv, mRgas
!=============================================================================!
        if (linear.ne.1) return

!.... allocate the main storage area for the mean field

        mem = mem + ndof*nx*ny
        allocate (vm(ndof,nx,ny), STAT=ier)
        if (ier .ne. 0) call error('genmean$','Insufficient Memory$')

!.... read in the mean field

        open(unit=10, file='mean.dat', form='unformatted', status='old', &
             err=1020)
        read(10,end=1020,err=1020) mstep, mtime, mnx, mny, mnz, mndof, &
                                   mRe, mMa, mPr, mgamma, mcv
        read(10,end=1020,err=1020) vm
        close(10)
        
!       write(*,"(' GenMean allocated ===> ',1pe13.6,' words')") float(mem)

        return
1020    call error('genmean$','Error reading the mean field file$')
        end
