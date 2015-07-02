	program fp
	    implicit none
	    real*8 ans
            real*8 tmp
  	    real*8 eps
   	    real*8 x
    	    integer*8 count

    	    ans = 1
	    tmp = 1
	    eps = 1
	    x = -19.5
	    count = 0
   
    	    do while (eps > 1e-15)
     		count = count + 1
        	tmp = tmp * (x/count)
        	eps = abs(tmp)
        	ans = ans + tmp
                ! write(*,*) "count = ",count,"eps = ",eps
	    end do
    
    	    write(*,*) "exp(x) = ", exp(x)
  	    write(*,*) "ans = ", ans
	end program

