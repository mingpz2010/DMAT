! Fortran具有很多内置函数
! DOT_PRODUCT(A, B):	returns the dot product of A and B
! SUM(A):		returns the sum of the elements of A
! TRANSPOSE(A):	returns the transpose of the 2D array A
!
module tools
    implicit none
    
contains
    ! 向量之间的拷贝
    subroutine vector_copy(n, dest, src)
        implicit none
        integer:: n, i
        real, intent(out):: dest(:)
        real:: src(:)
        
        i = 1
        do while (i<=n)
            dest(i) = src(i)
            i = i + 1
        end do
    end subroutine
    ! 向量与矩阵的乘积
    subroutine matrix_mul_vector(n, a, x)
        implicit none
        integer:: n,i,j
        real:: a(:,:)
        real, intent(out):: x(:)
        real, dimension(n):: b
        real:: su
        
        i = 1
        do while (i<=n)
            b(i) = x(i)
            i = i + 1
        end do
        
        i = 1
        do while (i<=n)
            su = 0.
            j = 1
            do while (j<=n)
                su = su + a(i,j)*b(j)
                j = j + 1
            end do
            x(i) = su
            i = i + 1
        end do
    end subroutine
    
    ! 向量与矩阵的乘积(2)
    subroutine matrix_mul_vector2(n, y, a, x)
        implicit none
        integer:: n,i,j
        real:: a(:,:), x(:)
        real, intent(out):: y(:)
        real:: su
        
        i = 1
        do while (i<=n)
            su = 0.
            j = 1
            do while (j<=n)
                su = su + a(i,j)*x(j)
                j = j + 1
            end do
            y(i) = su
            i = i + 1
        end do
    end subroutine
        
    ! 矩阵乘以一个数字
    subroutine matrix_mul_number(n, a, x)
        implicit none
        integer:: n,i,j
        real, intent(out):: a(:,:)
        real:: x
        
        i = 1
        do while (i<=n)
            j = 1
            do while (j<=n)
                a(i,j) = a(i,j)*x
                j = j + 1
            end do
            i = i + 1
        end do
    end subroutine
    
    ! 求一个向量的最大范数
    real function max_norm(n, a)
        implicit none
        integer, intent(in):: n
        integer:: i
        real, dimension(n):: a
        
        max_norm = abs(a(1))
        i = 1
        do while (i <= n)
            if (abs(a(i)) > max_norm) then
                max_norm = abs(a(i))
            end if
            i = i+1
        end do
    end function max_norm
    
    ! 求两个向量差异的最大范数
    real function max_norm2(n, x1, x2)
        implicit none
        integer:: n
        integer:: i
        real, dimension(n):: x1, x2
        
        max_norm2 = abs((x1(1) - x2(1))/(x1(1)))
        i = 1
        do while (i <= n)
            if (abs((x1(1) - x2(1))/(x1(1))) > max_norm2) then
                max_norm2 = abs((x1(1) - x2(1))/(x1(1)))
            end if
            i = i+1
        end do
    end function max_norm2
    
    ! Gauss-Seidel Method
    subroutine seidel(n, a, x, b, eps)
        implicit none
        integer:: n,i,j,itr
        real, intent(out):: x(:)
        real:: a(:,:), b(:), eps
        real, dimension(n):: tmp, residual_vec
        real:: residual, factor, left, right
        
        ! print *,"Seidel method: n = ",n, " eps = ", eps
        residual_vec(1) = 0.1
        tmp(1) = 0.1
        itr = 1
        do while (itr<=100)
            i = 1
            do while (i<=n)
                factor = 1/a(i, i)
                left = 0.
                right = 0.
                j = 1
                do while (j<i)
                    left = left + a(i,j) * x(j)
                    j = j + 1
                end do
                j = i
                do while (j<=n)
                    right = right + a(i,j) * x(j)
                    j = j + 1
                end do
                x(i) = x(i) + factor * (b(i)-left-right)
                i = i+1
            end do
 
 !----------------------------------------------------           
            ! 收敛条件判断
            ! 第一种停机准则，判断残差的最大范数
            i = 1
            do while (i<=n)
                tmp(i) = x(i)
                i = i + 1
            end do
            call matrix_mul_vector(n,a,tmp)
            i = 1
            do while (i<=n)
                tmp(i) = b(i) - tmp(i)
                i = i + 1
            end do
            residual = max_norm(n, tmp)
 !----------------------------------------------------   
 !----------------------------------------------------              
            !i = 1
            !do while (i<=n)
            !    residual_vec(i) = x(i) - tmp(i)
            !    i = i + 1
            !end do
            ! 第二种停机准则，判断相邻两次迭代值差的最大范数
            ! residual = max_norm(n, residual_vec)
!----------------------------------------------------
 
!----------------------------------------------------   
            ! 第三种停机准则，判断相邻两次迭代值的相对误差            
            !residual = max_norm(n, residual_vec)/max_norm(n, tmp)
            !i = 1
            !do while (i<=n)
            !    tmp(i) = x(i)
            !    i = i + 1
            !end do
!----------------------------------------------------
            
            if (residual < eps) then
                print *, "convergence, itr = ", itr, ", residual = ", residual
                exit
            end if
            
            itr = itr + 1  ! 防止迭代一直不收敛
            ! 调试用
            ! print *, "itr = ", itr, ", residual = ", residual
        end do
    end subroutine
    
    subroutine chase_method(n, a, x, b, eps)
        implicit none
        integer:: n,i,j
        real, intent(out):: x(:)
        real:: a(:,:), b(:), eps
        real, dimension(n):: aa, dia, c, d, beta
        
        i = 1
        ! 追赶法中对中间向量赋值
        do while (i<=n)
            dia(i) = a(i,i)
            if (i>1) then
                aa(i) = a(i,i-1)
            end if
            if (i<n) then
                c(i) = a(i,i+1)
            end if
            i = i + 1
        end do
        
        beta(1) = dia(1)
        x(1) = b(1)
        i = 2
        do while (i<=n)
            d(i) = aa(i)/beta(i-1)
            beta(i) = dia(i) - d(i)*c(i-1)
            x(i) = b(i) - d(i)*x(i-1)
            i = i + 1
        end do
        x(n) = x(n) / beta(n)
        i = n-1
        do while (i>=1)
            x(i) = (x(i)-c(i)*x(i+1))/beta(i)
            i = i - 1
        end do
        
    end subroutine

end module tools

