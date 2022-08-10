!=============================================================================!
        subroutine genbump( g1v, g2v, g11v, g12v, g22v )
!
!  Generate the wall bump
!
!  1. Changed to IJ ordering
!  2. In this version, the gradients coming in are in transform space so they
!     are first converted to physical space.
!
!  Revised: 8/8/22
!
!  Notes:
!    1. btype = 0 is bump relative to a flat plate
!    2. btype = 1 is bump relative to a parabolic cylinder
!    3. Bumps are all currently Gausian with center bmu and standard 
!       deviation of bsigma
!=============================================================================!
        use global
        use bump
        implicit none
        
        real    :: g1v(ndof,nx,ny), g2v(ndof,nx,ny)
        real    :: g11v(ndof,nx,ny), g12v(ndof,nx,ny), g22v(ndof,nx,ny)
        character (80) :: code='GenBump$'

        real :: g1vl(ndof), g2vl(ndof), g11vl(ndof), g12vl(ndof), g22vl(ndof)

        integer :: i, j, idof
        real    :: s(nx)
        
        real    :: bmu=0.7, bsigma=0.05
        integer :: btype=0
        logical :: bexist=.false., becho=.false., boutput=.true.
        character(80) :: bfile="bump.nml"

        namelist /bumpInput/ bmu, bsigma, btype, bfile, becho, boutput
!=============================================================================!
            
!.... read the mean and standard deviation from the bump.inp file

        inquire(file=bfile,exist=bexist)
        if (.not.bexist) then
          write(*,*) "Bump namelist file ", trim(bfile), &
                     "does not exist, trying input from bump.inp"
          bfile="bump.inp"
          open(30,file=bfile,status='old',err=1000)
          read(30,*,err=1010,end=1010) bmu, bsigma, btype
          close(30)
        else
          open(30,file=bfile,status='old',err=1100)
          read(30,bumpInput,err=1110)
          close(30)
          if (becho) write(*,bumpInput)
        endif

!.... specify the bump type

        if (btype .eq. 0) then
          write(*,*) "Bump is relative to flat plate..."
          s(:) = x(:,1)
        else if (btype .eq. 1) then
          write(*,*) "Bump is relative to parabolic cylinder..."
          s(:) = Sqrt(x(:,1) + two*x(:,1)**2)/Sqrt(two) + &
                 Log(one + four*x(:,1) + two**(onept5)*Sqrt(x(:,1) + &
                 two*x(:,1)**2)) / four
        else
          call error(code,'Illegal value of btype$')
        end if

!.... Allocate bump storage

        allocate( hw(nx) )
        allocate( u1p(nx), u2p(nx), u3p(nx), tp(nx), tpp(nx) )
        allocate( u1w(nx), u2w(nx), u3w(nx), tw(nx), twp(nx) )

        j = 1
        loop_i:  do i = 1, nx

!.... transform to physical space

        if (general) then

        do idof = 1, ndof
          g11vl(idof) = g11v(idof,i,j)       * m1m1(i,j)        + &
                        two * g12v(idof,i,j) * m1n1(i,j)        + &
                        g22v(idof,i,j)       * n1n1(i,j)        + &
                        g1v(idof,i,j)        * m11(i,j)         + &
                        g2v(idof,i,j)        * n11(i,j)

          g12vl(idof) = g11v(idof,i,j) * m1m2(i,j)                + &
                        g12v(idof,i,j) * (m1n2(i,j) + m2n1(i,j) ) + &
                        g22v(idof,i,j) * n1n2(i,j)                + &
                        g1v(idof,i,j)  * m12(i,j)                 + &
                        g2v(idof,i,j)  * n12(i,j)

          g22vl(idof) = g11v(idof,i,j)       * m2m2(i,j)        + &
                        two * g12v(idof,i,j) * m2n2(i,j)        + &
                        g22v(idof,i,j)       * n2n2(i,j)        + &
                        g1v(idof,i,j)        * m22(i,j)         + &
                        g2v(idof,i,j)        * n22(i,j)

          g1vl(idof) = g1v(idof,i,j)*m1(i,j) + g2v(idof,i,j)*n1(i,j)
          g2vl(idof) = g1v(idof,i,j)*m2(i,j) + g2v(idof,i,j)*n2(i,j)
        end do

        else  ! Cartesian

        do idof = 1, ndof
          g11vl(idof) = g11v(idof,i,j) * m1m1(i,j) + g1v(idof,i,j)  * m11(i,j)
          g12vl(idof) = g12v(idof,i,j) * m1n2(i,j)
          g22vl(idof) = g22v(idof,i,j) * n2n2(i,j) + g2v(idof,i,j)  * n22(i,j)
          g1vl(idof) = g1v(idof,i,j) * m1(i,j)
          g2vl(idof) = g2v(idof,i,j) * n2(i,j)
        end do

        end if

!.... Normalized Gaussian bump

        hw(i) = exp( -pt5 * ( ( s(i) - bmu ) / bsigma )**2 )
        
!.... note there is a negative sign since the normal points outside the domain

        u1p(i) = -( bnb(i,1) * g1vl(2) + bnb(i,2) * g2vl(2) )
        u2p(i) = -( bnb(i,1) * g1vl(3) + bnb(i,2) * g2vl(3) )
        u3p(i) = -( bnb(i,1) * g1vl(4) + bnb(i,2) * g2vl(4) )
        tp(i)  = -( bnb(i,1) * g1vl(5) + bnb(i,2) * g2vl(5) )
        tpp(i) = bnb(i,1) * (bnb(i,1) * g11vl(5) + bnb(i,2) * g12vl(5)) + &
                 bnb(i,2) * (bnb(i,1) * g12vl(5) + bnb(i,2) * g22vl(5))

        u1w(i) = -u1p(i) * hw(i)
        u2w(i) = -u2p(i) * hw(i)
        u3w(i) = -u3p(i) * hw(i)
        tw(i)  = -tp(i)  * hw(i)
        twp(i) = -tpp(i) * hw(i)

        end do loop_i

!.... save the bump profile 
        
        if (boutput) then
          write(*,*) 'Saving bump geometry to file bump.dat...'
          open(30,file='bump.dat')
          do i = 1, nx
            write(30,10) s(i), hw(i), real(u1w(i)), real(u2w(i)), &
                         real(u3w(i)), real(tw(i)), real(twp(i))
          end do
          close(30)
        endif
        
        return
  10    format( 8(1pe16.7e4,1x) )
1000    call error(code,'Error opening bump.inp$')
1010    call error(code,'Error reading bump.inp$')
1100    call error(code,'Error opening bump.nml$')
1110    call error(code,'Error reading bump.nml$')
        end subroutine genbump
