program unipot
      
      implicit none

      integer :: i, nx

      write(*,"('Enter n ==> ',$)")
      read(*,*) nx

      open(10,file='pot.out')
      do i = 1, nx
        write(10,"(5(1pe13.6,1x))") 1.0, 1.0, 0.0, 0.0, 1.0
      end do
      close(10)

      stop
end program unipot
      
