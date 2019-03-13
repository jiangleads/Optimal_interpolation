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
	!use standard_compute 
	implicit none


	integer             :: i,j,is,ie,js,je,it,posi,l
	integer,parameter   :: nx=57,ny=57
	real,parameter      :: lon0=162.75-21,lat0=42.0-21!,lon1=114.5,lat1=38  !//起止经纬度!以上取定原始球坐标的计算范围
	real,parameter      :: dr=10.,clon=162.75,clat=42.0,rd=6.371E6 !(clon,clat):极坐标中心经纬度,须在加密后的网格点上
	real,parameter      :: resx=0.75,resy=0.75!  原始球坐标x和y向格距都是0.5度
	integer,parameter   :: mx=1,my=1  !//x,y方向加密倍数! 插值后将原始相邻两个格点分成的子区域个数(分辨率为dax/MX;resy/MY)
	integer,parameter   :: nu=(nx-1)*mx+1,nv=(ny-1)*my+1  !//加密后各方向格点数! 插值后的细网格x向(NU)和y(NV)向格点数
	real,parameter      :: abc = atan(1.0)/45.0!	abc是角度制与弧度制的换算比例:弧度/角度 角度→弧度

	!//转换极坐标参数
	integer,parameter   :: ntheta=36,ndd=2000,ibjx=10,ibjy=10
	!	ntheta=36,即36个可选半径，每个半径相隔5度
	!	ndd为以极坐标原点为中心覆盖的极径大小（单位 km）
	!	ibj为插值所选的计算半径(格点数)。ibjx，ibjy为x向和y向的格点计算范围



	integer,parameter      :: i_pg=int(ndd/dr) !//i_pg为沿极径方向的格点数（一个特征半径(搜索半径?)上的均分格点数，不包括极点）
	real,parameter      :: dlmt = resx*abc!	x向格距(单位:弧度)
	real,parameter      :: dsit = resy*abc!	y向格距(单位:弧度)
	real,parameter      :: dy =(rd*dsit)/my!//细网格的y向格距，单位：m 其中rd*dsit是y向步长，单位：m
	real			:: h(nu,nv),xlon(nu,nv),ylat(nu,nv),p_h(i_pg,ntheta),p_r(i_pg,ntheta),p_a(i_pg,ntheta),zjx,zjy,avdx,dx(nv),xrc(nu,nv)&
		,yrc(nu,nv),rr,ssum,hh(nx,ny),hxlon(nx),hylat_t(ny),hh_t(nx,ny),hxylat_t(nv),hxxlon(nu),xh_t(nu,nv),rrx(i_pg,ntheta),rry(i_pg,ntheta)
	!	h为细网格上的要插值的变量,xlon,ylat为细网格点(nu,nv)的经纬度
	!	P_H为插值后极坐标各个同心圆上的变量值
	!	P_R为对应的极坐标极径值；P_A为对应的极坐标极角值
	!	XRC、YRC为细网格相对于涡旋中心的坐标，单位：m
	

	real                :: x1 !//范围内格点的纬度序列,j越大纬度越高

	integer, parameter :: n=ntheta*i_pg, on=nx*ny, m=100 !输入输出场大小, m影响半径

	!integer,parameter :: wp=4
	real(wp) :: x(2,on), var(on), f(1,on),x2(2,on), var2(on), f2(1,on)
	! x(1,:) and x(2,:) represent the x and y coordindate of the observations
	! f and var are observations and their error variance resp.
	! x(1,:) 和 x(2,:)表示观测点的x和y坐标
	! f和var是观测场和他的误差方差

	real(wp) :: xi(2,n), fi(1,n), vari(n)
	! xi(1,:) and xi(2,:) represent the x and y coordindate of the grid of the interpolated field
	! fi and vari are the interpolated field and its error variance resp.
	! xi(1,:)和xi(2,:)表示插值格点的x和y坐标. fi和vari分别是插值场和它的误差方差


	real(wp) :: param(2)
	! param: inverse of the correlation length
	! param：相关长度的倒数
	
	real time_begin,time_end1

	call CPU_TIME(time_begin)

	! param is the inverse of the correlation length
	! param 是相关长度的倒数
	param = 1./(500*1000)

	! the error variance of the observations
	! 观测的误差方差
	var = 0.005







	!//打开文件，读取气旋地面数据
	open (10,file="G:\ec79-15\1986_DP48.grd",form='unformatted',access='stream',action='read')
	open (11,file="D:\Desktop\chazhi\temp.grd",form='unformatted',access='stream',action='write')
	it=4
	posi=4*(nx*ny*(8*22+4*1)*(it-1))+1
	read(10,Pos=posi)((hh(i,j),i=1,nx),j=1,ny)
	!write(11)hh
	close(11)
	close (10)

	!//计算网格点的经纬度坐标
	do j = 1,ny
		do i = 1,nx
			xlon(i,j) = real(lon0+(i-1)*(resx/mx))
			ylat(i,j) = real(lat0+(j-1)*(resy/my))
		end do
	end do

	!//网格加密

	!细网格点的(绝对)经纬度坐标转换成相对极坐标中心的(经纬度)坐标
	xlon = xlon - clon
	ylat = ylat - clat

	!//计算细网格x方向的平均格距，单位：m
	avdx = 0.
	do j = 1,nv
		x1 = lat0*abc + float(j-1)*dsit/my !范围内细网格格点的纬度序列,j越大纬度越高
		avdx=avdx+rd*cos(x1)*dlmt  !//范围内细网格的x向格距的累加，(单位:m)
		dx(j) = rd*cos(x1)*dlmt/mx     !细网格上每个纬度的x向格距，单位：m
	end do
	avdx = avdx/real(nv)  !计算细网格平均x向格距，单位：m
	print *,'meandx=',avdx  !//细网格平均x向格距，单位：m
	print *,'dy=',dy     !//细网格平均y向格距，单位：m


	do j = 1,nv
		do i = 1,nu
			!//计算细网格绝对坐标下,点(i,j)到涡旋中心(x、y方向)的距离。单位：m
			xrc(i,j) = rd*cos(lat0*abc+float(j-1)*dsit/my)*(xlon(i,j)*abc)
			yrc(i,j) = ylat(i,j)*abc*rd
		end do
	end do
	print *,'相对距离完毕!'




	!距离赋值
	x(1,:)=reshape(xrc,[nu*nv])
	x(2,:)=reshape(yrc,[nu*nv])

	! the underlying function to interpolate
	! 要插值的函数
	f(1,:) = reshape(hh,[nu*nv])











	!//建立极坐标系，极径增量dr*1000m，极角增量10°
	open(11,file='D:\Desktop\chazhi\test0.txt',Form='formatted')
	do j=1,on
		write(11,200)x(1,j),x(2,j),f(1,j)
	end do
200	format(3(1x,f15.3))




	open (20,File='D:\Desktop\chazhi\R0.txt',Form='formatted')
	open (30,File='D:\Desktop\chazhi\A0.txt',Form='formatted')
	open (50,File='D:\Desktop\chazhi\center-t.txt',Form='formatted')
	open (60,File='D:\Desktop\chazhi\rx.txt',Form='formatted')
	open (70,File='D:\Desktop\chazhi\ry.txt',Form='formatted')
	do i = 1,i_pg
		do j = 1,ntheta
			p_r(i,j) = i*dr*1000. !//极径的值 单位：m，从离极点最近的同心圆算起，随i增大
			p_a(i,j) = real(j-1)*(360/ntheta)*abc !//极角的36个值，以x正方向为0度，随j增大
			rrx(i,j) = cos(p_a(i,j))*p_r(i,j) !转换成相对涡旋中心x方向的距离 单位：m
			rry(i,j) = sin(p_a(i,j))*p_r(i,j) !转换成相对涡旋中心y方向的距离 单位：m

			write (20,*) p_r(i,j)    !//P_R为对应的极坐标极径值
			write (30,*) p_a(i,j)/abc !//P_A为对应的极坐标极角值,单位:度
			write (60,*) rrx(i,j)   !输出数据 单位为角度(?)
			write (70,*) rry(i,j)   !输出数据 单位为角度(?)
		end do
	end do


	!输出涡旋中心坐标以及变量的值
	write (50,*) 0.,0.,h((1+(clon-lon0)/(resx/mx)),(1+(clat-lat0)/(resy/my)))
	close (50)
	print *,'极坐标建立成功!'

	close (20)
	print *,'极径输出完毕!'

	close (30)
	print *,'极角输出完毕!'

	close (60)
	close (70)
	print *,'RX，RY输出完毕!'





	!call 
     ! optiminterp(ox,of,ovar,param,m,gx,gf,gvar)


	!==============================================================================================================
	!===========================处理完毕，输出数据=====================================

	!x2=x
	!call optiminterp(x,f,var,param,m,x2,f2,var2)
	!open (40,File='D:\Desktop\chazhi\test.txt',Form='formatted')
	!!open (41,File='D:\Desktop\chazhi\v2.txt',Form='formatted')
	!do i = 1,on
	!	!write (*,*) xi(1,i),xi(2,i),fi(1,i)
	!	write (40,*) x2(1,i),x2(2,i),f2(1,i)
	!	
	!end do
	!write (*,*) vari
	!do i=1,on,100
	!	write(*,*)x2(1,i),x2(2,i),f2(1,i)
	!end do
	!close (40)
	
	!坐标赋值
	xi(1,:)=reshape(rrx,[n])
	xi(2,:)=reshape(rry,[n])

	call optiminterp(x,f,var,param,m,xi,fi,vari)
	! fi is the interpolated function and vari its error variance
	! fi 是插值函数, vari 是它的误差方差

	open (40,File='D:\Desktop\chazhi\test.txt',Form='formatted')
	!open (41,File='D:\Desktop\chazhi\v2.txt',Form='formatted')
	do i = 1,i_pg*ntheta
		write (40,*) xi(1,i),xi(2,i),fi(1,i)
		
	end do
	close (40)
	!write (*,*) vari
	
	!do i=1,n,100
	!	write(*,*)xi(1,i),xi(2,i),fi(1,i)
	!end do
	
	print *,'极坐标上的内插数据输出完毕!'
    print *,"m=",m,"var=",var(1),"param=1/",1/param(1)

call CPU_TIME(time_end1)

write(*,*) "计算耗时：",time_end1-time_begin
	!==============================================================================================================

	end program test_optiminterp

