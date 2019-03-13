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
	real,parameter      :: lon0=162.75-21,lat0=42.0-21!,lon1=114.5,lat1=38  !//��ֹ��γ��!����ȡ��ԭʼ������ļ��㷶Χ
	real,parameter      :: dr=10.,clon=162.75,clat=42.0,rd=6.371E6 !(clon,clat):���������ľ�γ��,���ڼ��ܺ���������
	real,parameter      :: resx=0.75,resy=0.75!  ԭʼ������x��y���඼��0.5��
	integer,parameter   :: mx=1,my=1  !//x,y������ܱ���! ��ֵ��ԭʼ�����������ֳɵ����������(�ֱ���Ϊdax/MX;resy/MY)
	integer,parameter   :: nu=(nx-1)*mx+1,nv=(ny-1)*my+1  !//���ܺ����������! ��ֵ���ϸ����x��(NU)��y(NV)������
	real,parameter      :: abc = atan(1.0)/45.0!	abc�ǽǶ����뻡���ƵĻ������:����/�Ƕ� �Ƕȡ�����

	!//ת�����������
	integer,parameter   :: ntheta=36,ndd=2000,ibjx=10,ibjy=10
	!	ntheta=36,��36����ѡ�뾶��ÿ���뾶���5��
	!	nddΪ�Լ�����ԭ��Ϊ���ĸ��ǵļ�����С����λ km��
	!	ibjΪ��ֵ��ѡ�ļ���뾶(�����)��ibjx��ibjyΪx���y��ĸ����㷶Χ



	integer,parameter      :: i_pg=int(ndd/dr) !//i_pgΪ�ؼ�������ĸ������һ�������뾶(�����뾶?)�ϵľ��ָ���������������㣩
	real,parameter      :: dlmt = resx*abc!	x����(��λ:����)
	real,parameter      :: dsit = resy*abc!	y����(��λ:����)
	real,parameter      :: dy =(rd*dsit)/my!//ϸ�����y���࣬��λ��m ����rd*dsit��y�򲽳�����λ��m
	real			:: h(nu,nv),xlon(nu,nv),ylat(nu,nv),p_h(i_pg,ntheta),p_r(i_pg,ntheta),p_a(i_pg,ntheta),zjx,zjy,avdx,dx(nv),xrc(nu,nv)&
		,yrc(nu,nv),rr,ssum,hh(nx,ny),hxlon(nx),hylat_t(ny),hh_t(nx,ny),hxylat_t(nv),hxxlon(nu),xh_t(nu,nv),rrx(i_pg,ntheta),rry(i_pg,ntheta)
	!	hΪϸ�����ϵ�Ҫ��ֵ�ı���,xlon,ylatΪϸ�����(nu,nv)�ľ�γ��
	!	P_HΪ��ֵ���������ͬ��Բ�ϵı���ֵ
	!	P_RΪ��Ӧ�ļ����꼫��ֵ��P_AΪ��Ӧ�ļ����꼫��ֵ
	!	XRC��YRCΪϸ����������������ĵ����꣬��λ��m
	

	real                :: x1 !//��Χ�ڸ���γ������,jԽ��γ��Խ��

	integer, parameter :: n=ntheta*i_pg, on=nx*ny, m=100 !�����������С, mӰ��뾶

	!integer,parameter :: wp=4
	real(wp) :: x(2,on), var(on), f(1,on),x2(2,on), var2(on), f2(1,on)
	! x(1,:) and x(2,:) represent the x and y coordindate of the observations
	! f and var are observations and their error variance resp.
	! x(1,:) �� x(2,:)��ʾ�۲���x��y����
	! f��var�ǹ۲ⳡ����������

	real(wp) :: xi(2,n), fi(1,n), vari(n)
	! xi(1,:) and xi(2,:) represent the x and y coordindate of the grid of the interpolated field
	! fi and vari are the interpolated field and its error variance resp.
	! xi(1,:)��xi(2,:)��ʾ��ֵ����x��y����. fi��vari�ֱ��ǲ�ֵ������������


	real(wp) :: param(2)
	! param: inverse of the correlation length
	! param����س��ȵĵ���
	
	real time_begin,time_end1

	call CPU_TIME(time_begin)

	! param is the inverse of the correlation length
	! param ����س��ȵĵ���
	param = 1./(500*1000)

	! the error variance of the observations
	! �۲������
	var = 0.005







	!//���ļ�����ȡ������������
	open (10,file="G:\ec79-15\1986_DP48.grd",form='unformatted',access='stream',action='read')
	open (11,file="D:\Desktop\chazhi\temp.grd",form='unformatted',access='stream',action='write')
	it=4
	posi=4*(nx*ny*(8*22+4*1)*(it-1))+1
	read(10,Pos=posi)((hh(i,j),i=1,nx),j=1,ny)
	!write(11)hh
	close(11)
	close (10)

	!//���������ľ�γ������
	do j = 1,ny
		do i = 1,nx
			xlon(i,j) = real(lon0+(i-1)*(resx/mx))
			ylat(i,j) = real(lat0+(j-1)*(resy/my))
		end do
	end do

	!//�������

	!ϸ������(����)��γ������ת������Լ��������ĵ�(��γ��)����
	xlon = xlon - clon
	ylat = ylat - clat

	!//����ϸ����x�����ƽ����࣬��λ��m
	avdx = 0.
	do j = 1,nv
		x1 = lat0*abc + float(j-1)*dsit/my !��Χ��ϸ�������γ������,jԽ��γ��Խ��
		avdx=avdx+rd*cos(x1)*dlmt  !//��Χ��ϸ�����x������ۼӣ�(��λ:m)
		dx(j) = rd*cos(x1)*dlmt/mx     !ϸ������ÿ��γ�ȵ�x���࣬��λ��m
	end do
	avdx = avdx/real(nv)  !����ϸ����ƽ��x���࣬��λ��m
	print *,'meandx=',avdx  !//ϸ����ƽ��x���࣬��λ��m
	print *,'dy=',dy     !//ϸ����ƽ��y���࣬��λ��m


	do j = 1,nv
		do i = 1,nu
			!//����ϸ�������������,��(i,j)����������(x��y����)�ľ��롣��λ��m
			xrc(i,j) = rd*cos(lat0*abc+float(j-1)*dsit/my)*(xlon(i,j)*abc)
			yrc(i,j) = ylat(i,j)*abc*rd
		end do
	end do
	print *,'��Ծ������!'




	!���븳ֵ
	x(1,:)=reshape(xrc,[nu*nv])
	x(2,:)=reshape(yrc,[nu*nv])

	! the underlying function to interpolate
	! Ҫ��ֵ�ĺ���
	f(1,:) = reshape(hh,[nu*nv])











	!//����������ϵ����������dr*1000m����������10��
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
			p_r(i,j) = i*dr*1000. !//������ֵ ��λ��m�����뼫�������ͬ��Բ������i����
			p_a(i,j) = real(j-1)*(360/ntheta)*abc !//���ǵ�36��ֵ����x������Ϊ0�ȣ���j����
			rrx(i,j) = cos(p_a(i,j))*p_r(i,j) !ת���������������x����ľ��� ��λ��m
			rry(i,j) = sin(p_a(i,j))*p_r(i,j) !ת���������������y����ľ��� ��λ��m

			write (20,*) p_r(i,j)    !//P_RΪ��Ӧ�ļ����꼫��ֵ
			write (30,*) p_a(i,j)/abc !//P_AΪ��Ӧ�ļ����꼫��ֵ,��λ:��
			write (60,*) rrx(i,j)   !������� ��λΪ�Ƕ�(?)
			write (70,*) rry(i,j)   !������� ��λΪ�Ƕ�(?)
		end do
	end do


	!����������������Լ�������ֵ
	write (50,*) 0.,0.,h((1+(clon-lon0)/(resx/mx)),(1+(clat-lat0)/(resy/my)))
	close (50)
	print *,'�����꽨���ɹ�!'

	close (20)
	print *,'����������!'

	close (30)
	print *,'����������!'

	close (60)
	close (70)
	print *,'RX��RY������!'





	!call 
     ! optiminterp(ox,of,ovar,param,m,gx,gf,gvar)


	!==============================================================================================================
	!===========================������ϣ��������=====================================

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
	
	!���긳ֵ
	xi(1,:)=reshape(rrx,[n])
	xi(2,:)=reshape(rry,[n])

	call optiminterp(x,f,var,param,m,xi,fi,vari)
	! fi is the interpolated function and vari its error variance
	! fi �ǲ�ֵ����, vari ����������

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
	
	print *,'�������ϵ��ڲ�����������!'
    print *,"m=",m,"var=",var(1),"param=1/",1/param(1)

call CPU_TIME(time_end1)

write(*,*) "�����ʱ��",time_end1-time_begin
	!==============================================================================================================

	end program test_optiminterp

