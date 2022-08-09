!=============================================================================!
        subroutine genbump( g1v, g2v, g11v, g12v, g22v )
!
!  Generate the wall bump
!
!  1. Changed to IJ ordering
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
        character(80) :: code='GenBump$'

        integer :: i
        real    :: s(nx)
        
        real    :: bmu=0.7, bsigma=0.05
        integer :: btype=0
        logical :: bexist=.false., becho=.false., boutput=.true.
        character(80) :: bfile="bump.nml"

        namelist /bumpInput/ bmu, bsigma, btype, bfile, becho
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
          s = x(:,1)
        else if (btype .eq. 1) then
          write(*,*) "Bump is relative to parabolic cylinder..."
          s = Sqrt(x(:,1) + two*x(:,1)**2)/Sqrt(two) + &
              Log(one + four*x(:,1) + two**(onept5)*Sqrt(x(:,1) + &
              two*x(:,1)**2)) / four
        else
          call error(code,'Illegal value of btype$')
        end if

!.... Normalized Gaussian bump

        allocate( hw(nx) )
        hw(:) = exp( -pt5 * ( ( s - bmu ) / bsigma )**2 )
        
        allocate( u1p(nx), u2p(nx), u3p(nx), tp(nx), tpp(nx) )
        allocate( u1w(nx), u2w(nx), u3w(nx), tw(nx), twp(nx) )

!.... note there is a negative sign since the normal points outside the domain

        u1p = -( bnb(:,1) * g1v(2,:,1) + bnb(:,2) * g2v(2,:,1) )
        u2p = -( bnb(:,1) * g1v(3,:,1) + bnb(:,2) * g2v(3,:,1) )
        u3p = -( bnb(:,1) * g1v(4,:,1) + bnb(:,2) * g2v(4,:,1) )
        tp  = -( bnb(:,1) * g1v(5,:,1) + bnb(:,2) * g2v(5,:,1) )
        tpp = bnb(:,1) * (bnb(:,1) * g11v(5,:,1) + bnb(:,2) * g12v(5,:,1)) + &
              bnb(:,2) * (bnb(:,1) * g12v(5,:,1) + bnb(:,2) * g22v(5,:,1))

        u1w = -u1p * hw
        u2w = -u2p * hw
        u3w = -u3p * hw
        tw  = -tp  * hw
        twp = -tpp * hw

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
