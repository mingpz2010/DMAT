module diffusion
    implicit none
    
    type material
        real, dimension(2):: tr, a, vf, kf
        real:: tr1to2, d1, d2
        real:: exposure
    end type material
    
contains    
	subroutine diffusion_solver(single_mesh)
	    use tools
		implicit none
		
		integer:: single_mesh, mesh, i, j, k, idx, itr
		integer, dimension(13):: geom=(/3,3,3,3,3,3,3,3,3,3,3,3,3/)
		type(material):: mat1
		type(material):: mat2
		type(material):: mat3
		real:: chi1, chi2, mesh_space
		real, allocatable:: flux1(:), flux2(:), source1(:), source2(:), p(:)
		real, allocatable:: flux1_last(:), flux2_last(:)
		real, allocatable:: mat_of_coeff1(:,:), mat_of_coeff2(:,:)
		real:: keff, keff_last, cjk, cjf1, cjf2
		
		!-------------------------------
		! 材料参数1, 2, 3
		!-------------------------------
		mat1%exposure = 0.2
		mat1%tr(1)=2.368355e-1; mat1%tr(2)=9.082422e-1
		mat1%a(1)=8.603111e-3; mat1%a(2)=7.853449e-2
		mat1%vf(1)=6.160544e-3; mat1%vf(2)=1.207603e-1
		mat1%kf(1)=8.158099e-14; mat1%kf(2)=1.599168e-12
		mat1%tr1to2=1.708253e-2
		mat1%d1 = 1./(3*mat1%tr(1))  ! 快群扩散系数计算
		mat1%d2 = 1./(3*mat1%tr(2))  ! 满群扩散系数计算
		print *, mat1%d1, mat1%d2

		mat2%exposure = 8.11
		mat2%tr(1)=2.367121e-1; mat2%tr(2)=9.239672e-1
		mat2%a(1)=9.150915e-3; mat2%a(2)=8.516051e-2
		mat2%vf(1)=5.585696e-3; mat2%vf(2)=1.250261e-1
		mat2%kf(1)=7.144699e-14; mat2%kf(2)=1.599212e-12
		mat2%tr1to2=1.690581e-2
		mat2%d1 = 1./(3*mat2%tr(1))  ! 快群扩散系数计算
		mat2%d2 = 1./(3*mat2%tr(2))  ! 满群扩散系数计算
		print *, mat2%d1, mat2%d2

		mat3%exposure = 16.55
		mat3%tr(1)=2.366212e-1; mat3%tr(2)=9.308326e-1
		mat3%a(1)=9.668583e-3; mat3%a(2)=8.506164e-2
		mat3%vf(1)=5.050670e-3; mat3%vf(2)=1.188626e-1
		mat3%kf(1)=6.310672e-14; mat3%kf(2)=1.485153e-12
		mat3%tr1to2=1.675986e-2
		mat3%d1 = 1./(3*mat3%tr(1))  ! 快群扩散系数计算
		mat3%d2 = 1./(3*mat3%tr(2))  ! 满群扩散系数计算
		print *, mat3%d1, mat3%d2
		
		!-------------------------------
		! 能谱参数
		!-------------------------------
		chi1 = 1.0; chi2 = 0.
		!-------------------------------
		! 迭代参数
		!-------------------------------
		keff = 1.0
		mesh = single_mesh * 13; mesh_space = 20.0/single_mesh
		allocate(flux1(mesh), flux2(mesh), source1(mesh), &
		        source2(mesh), p(mesh), flux1_last(mesh), flux2_last(mesh))
		allocate(mat_of_coeff1(mesh, mesh), mat_of_coeff2(mesh, mesh))
		print *, "1D Reactor geometry is : (mesh = ", mesh, ")"
        print *,geom
		
		! 网格划分之后，两群扩散方程的系数矩阵初始化
label1:	do i=1,13
		    k = 1
		    idx = (i-1) * single_mesh  ! 起始序号
label2:		do while (k <= single_mesh)
		        mat_of_coeff1(idx+k, idx+k) = (2.*mat3%d1)/ &
		                (mesh_space*mesh_space)+mat3%a(1)+mat3%tr1to2
		        mat_of_coeff2(idx+k, idx+k) = (2.*mat3%d2)/ &
		                (mesh_space*mesh_space)+mat3%a(2)
		        k = k + 1
		    end do label2
		end do label1
		
		i = 1
		do while (i<=mesh)
		    if (i>=2) then
		        mat_of_coeff1(i,i-1) = -mat3%d1/(mesh_space*mesh_space)
		        mat_of_coeff2(i,i-1) = -mat3%d2/(mesh_space*mesh_space)
		    end if
		    if (i<mesh) then
		        mat_of_coeff1(i,i+1) = -mat3%d1/(mesh_space*mesh_space)
		        mat_of_coeff2(i,i+1) = -mat3%d2/(mesh_space*mesh_space)
		    end if
		    i = i + 1
		end do
		!i = 1
		!do while (i<=mesh)
		!    mat_of_coeff1(mesh,i) = 0
		!    mat_of_coeff2(mesh,i) = 0
		!    i = i + 1
		!end do
		mat_of_coeff1(mesh, mesh) = 1
		mat_of_coeff2(mesh, mesh) = 1

		!print *, "mat_of_coeff1(DIA) : "
		!print *, mat_of_coeff1(1,1)
		!print *, mat_of_coeff1(mesh/2,mesh/2)
		!print *, mat_of_coeff1(mesh,mesh)
		!print *, "mat_of_coeff1(left) : "
		!print *, mat_of_coeff1(2,1)
		!print *, mat_of_coeff1(mesh/2,mesh/2-1)
		!print *, mat_of_coeff1(mesh,mesh-1)
		!print *, "mat_of_coeff1(right) : "
		!print *, mat_of_coeff1(1,2)
		!print *, mat_of_coeff1(mesh/2,mesh/2+1)
		!print *, mat_of_coeff1(mesh-1,mesh)
		
		! 初始迭代参数和源项赋值
		keff = 0.9
		itr = 1
		i = 1
		do while (i <= mesh)
		    flux1(i) = keff/((mat3%vf(1)+mat3%vf(2))*(mesh*mesh_space))
		    flux2(i) = flux1(i)
		    p(i) = mat3%vf(1)*flux1(i)+mat3%vf(2)*flux2(i)
		    i = i + 1
		end do
		
		! 主要的inner iteration部分
		cjk = 1
		cjf1 = 1
		cjf2 = 1
label3:	do while (cjk >= 1e-16 .or. cjf1>=1e-16 .or. cjf2>=1e-16)
			i = 1
label4:		do while (i < mesh)
				source1(i) = (chi1*p(i))/keff
				i = i + 1
			end do label4
			
			call vector_copy(mesh, flux1_last, flux1)
			call vector_copy(mesh, flux2_last, flux2)
			source1(mesh) = 0
			
			!print *, "source1 : "
		   	!print *, source1(1)
		   	!print *, source1(mesh/2)
		   	!print *, source1(mesh)
			
			call chase_method(mesh, mat_of_coeff1, flux1, source1, 1e-16)
			!print *, "flux1 : "
		   	!print *, flux1(1)
		   	!print *, flux1(mesh/2)
		   	!print *, flux1(mesh)
			! call seidel(mesh, mat_of_coeff1, flux1, source1, 1e-16)
			i = 1
label5:		do while (i <= mesh)
				source2(i) = (chi2*p(i))/keff + mat3%tr1to2 * flux1(i)
				i = i + 1
			end do label5
			call chase_method(mesh, mat_of_coeff2, flux2, source2, 1e-16)
			!call seidel(mesh, mat_of_coeff2, flux2, source2, 1e-16)
			
			keff_last = keff
			keff = 0
			
			i = 1
label6:		do while (i <= mesh)
				p(i) = mat3%vf(1)*flux1(i)+mat3%vf(2)*flux2(i)
				keff = keff + p(i) * mesh_space
				i = i + 1
			end do label6
			
			cjk = abs(keff-keff_last)/keff_last
			cjf1 = max_norm2(mesh, flux1, flux1_last)
			cjf2 = max_norm2(mesh, flux2, flux2_last)
			itr = itr + 1
			!if (itr == 2) then
			!    stop
			!end if
			print *,"itr = ", itr-1, ": keff = ", keff, ", CJK = ", cjk &
			        , "CJF1 = ", cjf1, "CJF2 = ", cjf2
		end do label3
		
		print *, "keff = ", keff
		
		deallocate(mat_of_coeff1, mat_of_coeff2)
		deallocate(flux1, flux2, source1, source2, p)
		deallocate(flux1_last, flux2_last)
	end subroutine diffusion_solver

end module diffusion

