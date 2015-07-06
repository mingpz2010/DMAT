!----------------------------------------------------------------!
!  鐢ㄩ�锛氶」鐩簩锛岀敤浜庡畬鎴愪綔涓氫竴缁翠袱缇ゅ爢鑺暟鍊艰绠�                  !
!  閲嶇偣锛氱瀛﹁绠�                                                !
!  鏃ユ湡锛�014.12                                                 !
!  浣滆�锛歮ingpz@mail.ustc.edu.cn                                 !
!----------------------------------------------------------------!
!
! Fortran arrays: first index accesses consecutive locations 
! (opposite in C)
!
! 鏈▼搴忕殑鐩殑鏄敤澶氱兢鎵╂暎鏂圭▼杩欑纭畾鐞嗚锛岄噰鐢ㄦ暟鍊艰绠楃殑鏂规硶绠楀嚭
! 鐕冩枡缁勪欢涓嶅悓鎺掑垪鎯呭喌涓嬩竴缁翠袱缇ゅ爢鑺殑keff銆佷腑瀛愰�閲忓垎甯冨強鍔熺巼鍒嗗竷銆�!
!

program main
    use diffusion
    use testing

    print *,"Start to perform 1D diffusion simulation"
    print *,"----------------------------------------"
    call my_testing2()
    call my_testing3()
    call diffusion_solver(1200)
end program main
