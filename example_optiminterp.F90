!! Copyright (C) 2005 Alexander Barth <barth.alexander@gmail.com>
!!
!! This program is free software; you can redistribute it and/or modify it under
!! the terms of the GNU General Public License as published by the Free Software
!! Foundation; either version 3 of the License, or (at your option) any later
!! version.
!!
!! This program is distributed in the hope that it will be useful, but WITHOUT
!! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!! FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
!! details.
!!
!! You should have received a copy of the GNU General Public License along with
!! this program; if not, see <http://www.gnu.org/licenses/>.

!! Example for using the n-dimentional optimal interpolation Fortran module

!!//tanslated into Chinese by Jiang Lizhi , 7 Jan 2019

program test_optiminterp
 use optimal_interpolation
 implicit none

 integer, parameter :: n=100, on=200, m=30
 ! n times n is the sized of the 2D background grid ! n: 2维背景场格点的大小，输出场
 ! on: number of observations                       ! on: 观测数目,输入场
 ! m: number of influential observations			! m: 影响观测的数目


 real(wp) :: xi(2,n*n), fi(1,n*n), vari(n*n)
 ! xi(1,:) and xi(2,:) represent the x and y coordindate of the grid of the interpolated field
 ! fi and vari are the interpolated field and its error variance resp.
 ! xi(1,:)和xi(2,:)表示插值格点的x和y坐标. fi和vari分别是插值场和它的误差方差
 

 real(wp) :: x(2,on), var(on), f(1,on)
 ! x(1,:) and x(2,:) represent the x and y coordindate of the observations
 ! f and var are observations and their error variance resp.
 ! x(1,:) 和 x(2,:)表示观测点的x和y坐标
 ! f和var是观测场和他的误差方差
 
 real(wp) :: param(2)
! param: inverse of the correlation length
! param：相关长度的倒数

 integer :: i,j,l

 ! we use a simple random generator instead of Fortran 90's random_number in the hope that the results will be platform independent
 ! 使用Fortran 90的 random_number 简单随机数生成器, 以希望结果会与平台无关

 integer(8), parameter :: A = 1664525_8, B = 1013904223_8, Mo = 4294967296_8
 integer(8) :: next = 0_8


 
 ! create a regular 2D grid
 ! 创建一个规则的2D格点
 l = 1
 do j=1,n
   do i=1,n
     xi(1,l) = (i-1.)/(n-1.)
     xi(2,l) = (j-1.)/(n-1.)
     l = l+1
   end do
 end do

 ! param is the inverse of the correlation length
 ! param 是相关长度的倒数
 param = 1./0.1

 
 ! the error variance of the observations
 ! 观测的误差方差
 var = 0.01

 ! location of observations
! 观测点的位置
 do j=1,on
   do i=1,2
     ! simple random generator
	   ! 简单随机数生成
     next = mod(A*next  + B,Mo)
     x(i,j) = real(next,8)/real(Mo,8)
   end do
 end do


 ! the underlying function to interpolate
! 要插值的函数
 f(1,:) = sin( x(1,:)*6 ) * cos( x(2,:)*6);

 
 ! fi is the interpolated function and vari its error variance
 ! fi 是插值函数, vari 是它的误差方差
 
 call optiminterp(x,f,var,param,m,xi,fi,vari)

 ! control value
 !//值的控制
 write(6,'(A,F10.6)') 'Expected:',2.2205936104348591E-002
 write(6,'(A,F10.6)') 'Obtained:',fi(1,1)

end program test_optiminterp

