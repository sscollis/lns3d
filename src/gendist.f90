!=============================================================================!
        subroutine gendist
!  
!  Generate the inflow disturbance if required
!  
!=============================================================================!
        use global
        use material
        implicit none

        real :: dummy
        integer :: j, mem=0
!=============================================================================!

!.... return if eigenfunction BC is not used

        if (xper .or. left.ne.4) return
        
!.... allocate space for the eigenfunctions

        mem = mem + ny * 10
        allocate ( rhor(ny) )
        allocate ( rhoi(ny) )
        allocate ( ur(ny) )
        allocate ( ui(ny) )
        allocate ( vr(ny) )
        allocate ( vi(ny) )
        allocate ( wr(ny) )
        allocate ( wi(ny) )
        allocate ( tr(ny) )
        allocate ( ti(ny) )

        open(unit=100,file='inflow.nml')
        read(100,inflow)
        if (inflow_echo) then
          write(*,inflow)
        endif
        close(100)
        write (*,*) 'Reading inflow profile from ',inflow_file
        write (*,*) '  with Ny = ', ny
        open(unit=ifile,file=inflow_file,form='formatted',status='unknown')
        do j = 1, ny 
          read(ifile,*) dummy,rhor(j),rhoi(j),ur(j),ui(j),  &
                        vr(j),vi(j),wr(j),wi(j),tr(j),ti(j)
          if (inflow_echo) then
            write(*,10) dummy,rhor(j),rhoi(j),ur(j),ui(j),    &
                        vr(j),vi(j),wr(j),wi(j),tr(j),ti(j)
          endif 
        end do
        close(ifile)
        write(*,"(' GenDist allocated ===> ',1pe13.6,' words')") float(mem)

        return
  10    format(1p,11(e13.6,1x))
        end
