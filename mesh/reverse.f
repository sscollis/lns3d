        program reverse
        
        real*8 x(500), y(500)
        
        open(10,file='lev.bot')
        
        i = 0
        
 10     continue
          i = i + 1
          read(10,*,end=20) x(i), y(i)
          goto 10
 20     continue
        close(10)
        
        num = i - 1
        
        open(10,file='new.bot')
        
        do i = 1, num
          write(10,*) x(num-i+1), y(num-i+1)
        end do
        close(10)
        
        stop
        end
