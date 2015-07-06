module testing
    use tools
    implicit none
    
contains
	subroutine my_testing1()
		use diffusion
		
		implicit none
		integer:: i, k
		real(kind=ieee_double):: answer, factor
		real(kind=ieee_double), dimension(8):: a = (/0.1,0.2,0.3,0.4,0.3,0.2,0.1,0.05/)
		real(kind=ieee_double), dimension(9, 9):: mat
		real(kind=ieee_double), dimension(9):: x, b
		
		answer = max_norm(8,a)
		print *,"max norm of a is ",answer
		i = 1
		do while (i<=9)
			mat(i,1) = 10-i
			mat(i,9) = i
			k = 2
			do while (k<=i)
				mat(i,k) = mat(i,k-1)+mat(i,1)
				k = k + 1
			end do
			k = 8
			do while (k>i)
				mat(i,k) = mat(i,k+1)+mat(i,9)
				k = k - 1
			end do
			i = i + 1
		end do
		call matrix_mul_number(9, mat, 1e-3)
		
		i = 1
		do while (i<=9)
			factor = i/10.
			b(i) = factor*(1-factor)*exp(factor)
			x(i) = 1
			i = i + 1
		end do		
		call seidel(9, mat, x, b, 1e-8)
	end subroutine my_testing1
	
	subroutine my_testing2()
        implicit none
        real(kind=ieee_double), dimension(10):: a=(/1,2,3,4,5,6,7,8,9,10/)
        
        print *, "sum of a = ", sum(a)
    end subroutine my_testing2
    
    subroutine my_testing3()
        implicit none
        real(kind=ieee_double):: a, b
        
        a = 1e-32
        b = 1e-33
        if (a > b) then
            print *, "a > b"
        end if
    end subroutine my_testing3
end module testing

