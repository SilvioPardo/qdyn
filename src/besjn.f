      program test_besjn
      real(8) :: x = 32.
      y = bessel_j1(x)
      write(6,*) x,y
      end program test_besjn
