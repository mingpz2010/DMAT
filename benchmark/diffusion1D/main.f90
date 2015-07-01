!----------------------------------------------------------------!
!  用途：项目二，用于完成作业一维两群堆芯数值计算                   !
!  重点：科学计算                                                 !
!  日期：2014.12                                                 !
!  作者：mingpz@mail.ustc.edu.cn                                 !
!----------------------------------------------------------------!
!
! Fortran arrays: first index accesses consecutive locations 
! (opposite in C)
!
! 本程序的目的是用多群扩散方程这种确定理论，采用数值计算的方法算出
! 燃料组件不同排列情况下一维两群堆芯的keff、中子通量分布及功率分布。
!
!

program main
    use diffusion
    use testing
    implicit none

    print *,"Start to perform 1D diffusion simulation"
    print *,"----------------------------------------"
    call my_testing2()
    call my_testing3()
    call diffusion_solver(1500)
end program main
