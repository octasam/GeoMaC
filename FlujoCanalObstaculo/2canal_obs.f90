!
!	canal_obs.f90
!	Particle tracking in a channel with a circular obstacle
!
!	Created by Fernando J. Guerrero-Martínez on 12/26/11.
!	Copyright 2011 __MyCompanyName__. All rights reserved.
!
!

program channel_obst
implicit none

integer i, j, k, nx, ny, c, nn, ei, ej, it, itmax, max_iter1, max_iter2, dbkx, dbky, npx, npy
real dx, dy, dxw, dys, dxe, dyn, dt, x0, xl, y0, yl, delv, ue, uw, vn, vs, se, sw, sn, ss
real  Re, gamma, div, rad, xcoord, ycoord, residual, tolerance

real, allocatable, dimension(:,:) :: p, pp, u, v, u1, v1, ae, aw, as, an, ap, sp, uc, vc, du, dv, partic
real, allocatable, dimension(:) :: xc, yc, x, y
integer, allocatable, dimension(:,:) :: mark_cells
character*50 filename, text
open(1,file='partic.txt')  !Output file for particle tracking over time

!Domain and mesh definition (uniform mesh)
x0=0.0
xl=10.0
y0=-0.5
yl=0.5

nx=600
ny=60

dx=(xl-x0)/float(nx)
dy=(yl-y0)/float(ny)

!Maximum interation for algorithms and convergence criterion
max_iter1=500 !GAUSS-TDMA
max_iter2=100 !SIMPLEC

tolerance=5.0E-5



allocate (xc(0:nx+1),yc(0:ny+1),x(0:nx),y(0:ny))
allocate (p(0:nx+1,0:ny+1), pp(0:nx+1,0:ny+1))
allocate (u(0:nx,0:ny+1), v(0:nx+1,0:ny), u1(0:nx,0:ny+1), v1(0:nx+1,0:ny))
allocate (ap(nx,ny), ae(nx,ny), aw(nx,ny), as(nx,ny), an(nx,ny), sp(nx,ny))
allocate (uc(0:nx+1,0:ny+1), vc(0:nx+1,0:ny+1), du(0:nx+1,0:ny+1), dv(0:nx+1,0:ny+1), partic(23040,2))
allocate (mark_cells(0:nx+1,0:ny+1))

!Mesh building
call mesh1d(xc, x, x0, xl, nx)
call mesh1d(yc, y, y0, yl, ny)

!======================================================================================
!               Initial conditions - time step - governing parameters
u=0.0; v=0.0; p=0.0; ap=0.0; ae=0.0; sp=0.0; du=0; dv=0

!u,v store velocity for time iterations - u1,v1 store vel for SIMPLEC iterations
u(0,:)=0.2  !inlet velocity
u1=u
v1=v

Re=2000.0
gamma=1.0/Re

dt=0.005

!Total time steps
itmax=12000

!======================================================================================
! Obstacle definition - make sure velocities are zero in the obstacle and its boundaires
rad=0.15; xcoord=3.0; ycoord=0.0 !radius and center of circle
mark_cells=0
	do i=1, nx
		do j=1, ny
			if(sqrt((xc(i)-xcoord)**2+(yc(j)-ycoord)**2)<rad)then
				mark_cells(i,j)=1
			end if
		end do
	end do

!initialize array to store position of prticles for the particle tracking
partic=0.0


!Time loop
do it = 1, itmax

!initialize counter for SIMPLEC iteration and Divergence
c=1
Div=1.0

!SIMPLEC loop
do while ((c < max_iter2) .and. (Div > tolerance))
!====================================================
!					EQUATION U
!====================================================
ap=0.0;	aw=0.0;	ae=0.0;	as=0.0;	an=0.0;	sp=0.0

ei=ubound(u,1)-1
ej=ubound(u,2)-1

do i=1, ei
	do j=1, ej
        dxw=x(i)-x(i-1) !staggered gird: with respect to volume faces
		dxe=x(i+1)-x(i)

		dys=yc(j)-yc(j-1) !staggered gird: with respect to volume centers
		dyn=yc(j+1)-yc(j)

		se=y(j)-y(j-1); sw=se
		sn=xc(i+1)-xc(i); ss=sn

		delv=(xc(i+1)-xc(i))*(y(j)-y(j-1))
		
		ue=0.5*(u1(i,j)+u1(i+1,j))
		uw=0.5*(u1(i,j)+u1(i-1,j))
		vn=0.5*(v1(i,j)+v1(i+1,j))
		vs=0.5*(v1(i,j-1)+v1(i+1,j-1))
		
		ae(i,j)=gamma*(se/dxe)+max(0.0,-ue)*se
		aw(i,j)=gamma*(sw/dxw)+max(0.0,uw)*sw
		an(i,j)=gamma*(sn/dyn)+max(0.0,-vn)*sn
		as(i,j)=gamma*(ss/dys)+max(0.0,vs)*ss
		ap(i,j)=ae(i,j)+aw(i,j)+an(i,j)+as(i,j)+delv/dt
		sp(i,j)=u(i,j)*delv/dt-(p(i+1,j)-p(i,j))*se

        ! CORRECTIONS FOR THE OBSTACLE
		if(mark_cells(i,j) .eq. 0 .and. mark_cells(i+2,j) .eq. 1) then
		sp(i,j)=sp(i,j)+ae(i,j)*u1(i+1,j)
		ae(i,j)=0.0
		end if	
		
		if(mark_cells(i,j) .eq. 0 .and. mark_cells(i-1,j) .eq. 1) then
		sp(i,j)=sp(i,j)+aw(i,j)*u1(i-1,j)
		aw(i,j)=0.0
		end if			
		
		if(mark_cells(i,j) .eq. 0 .and. (mark_cells(i+1,j+1) .eq. 1 .or. mark_cells(i,j+1) .eq. 1)) then
		ap(i,j)=ap(i,j)+an(i,j)
		sp(i,j)=sp(i,j)+2.0*an(i,j)*u1(i,j+1)
		an(i,j)=0.0
		end if	
		
		if(mark_cells(i,j) .eq. 0 .and. (mark_cells(i+1,j-1) .eq. 1 .or. mark_cells(i,j-1) .eq. 1)) then
		ap(i,j)=ap(i,j)+as(i,j)
		sp(i,j)=sp(i,j)+2.0*as(i,j)*u1(i,j-1)
		as(i,j)=0.0
		end if	
		
		if(mark_cells(i,j) .eq. 1 .or. (mark_cells(i,j) .eq. 0 .and. mark_cells(i+1,j) .eq. 1)) then
		aw(i,j)=0.0; ae(i,j)=0.0; as(i,j)=0.0; an(i,j)=0.0; u1(i,j)=0.0
		ap(i,j)=delv/dt
		sp(i,j)=u1(i,j)*delv/dt
		end if
		
	end do
end do

!       BOUNDARY CONDITIONS U

!Este: Neumann - Fully developed flow
aP(ei,1:ej)=aP(ei,1:ej)-aE(ei,1:ej)
aE(ei,1:ej)=0.0

!Oeste: Dirichlet - Inlet velocity
sP(1,1:ej)=sP(1,1:ej)+aW(1,1:ej)*u1(0,1:ej)
aW(1,1:ej)=0.0

!Norte: Dirichlet - Non-slip condition
sP(1:ei,ej)=sP(1:ei,ej)+aN(1:ei,ej)*u1(1:ei,ej+1)
aN(1:ei,ej)=0.0

!Sur: Dirichlet - Non-slip condition
sP(1:ei,1)=sP(1:ei,1)+aS(1:ei,1)*u1(1:ei,0)
aS(1:ei,1)=0.0

do i=1, ei
	do j=1, ej
		se=(y(j)-y(j-1))
		du(i,j)=se/(ap(i,j)-ae(i,j)-aw(i,j)-an(i,j)-as(i,j))
	end do
end do

CALL Gauss_Seidel2D(u1,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny,max_iter1, tolerance,residual)

!====================================================
!                  EQUATION V
!====================================================
ap=0.0;	aw=0.0;	ae=0.0;	as=0.0;	an=0.0;	sp=0.0

ei=ubound(v,1)-1
ej=ubound(v,2)-1

do i=1, ei
	do j=1, ej
		dxw=xc(i)-xc(i-1)   !staggered gird: with respect to volume centers
		dxe=xc(i+1)-xc(i)
		dys=y(j)-y(j-1)     !staggered gird: with respect to volume faces
		dyn=y(j+1)-y(j)

		se=yc(j+1)-yc(j); sw=se
		sn=x(i)-x(i-1); ss=sn
	
		delv=(x(i)-x(i-1))*(yc(j+1)-yc(j))
	
		ue=0.5*(u1(i,j)+u1(i,j+1))
		uw=0.5*(u1(i-1,j)+u1(i-1,j+1))
		vn=0.5*(v1(i,j)+v1(i,j+1))
		vs=0.5*(v1(i,j)+v1(i,j-1))
		
		ae(i,j)=gamma*(se/dxe)+max(0.0,-ue)*se
		aw(i,j)=gamma*(sw/dxw)+max(0.0,uw)*sw
		an(i,j)=gamma*(sn/dyn)+max(0.0,-vn)*sn
		as(i,j)=gamma*(ss/dys)+max(0.0,vs)*ss
		ap(i,j)=ae(i,j)+aw(i,j)+an(i,j)+as(i,j)+delv/dt
		sp(i,j)=v(i,j)*delv/dt-(p(i,j+1)-P(i,j))*sn
		
		if((mark_cells(i,j) .eq. 0) .and. (mark_cells(i+1,j+1) .eq. 1 .or. mark_cells(i+1,j) .eq. 1)) then
		ap(i,j)=ap(i,j)+ae(i,j)
		sp(i,j)=sp(i,j)+2.0*ae(i,j)*v1(i+1,j)
		ae(i,j)=0.0
		end if	
		
		if((mark_cells(i,j) .eq. 0) .and. (mark_cells(i-1,j+1) .eq. 1 .or. mark_cells(i-1,j) .eq. 1)) then
		ap(i,j)=ap(i,j)+aw(i,j)
		sp(i,j)=sp(i,j)+2.0*aw(i,j)*v1(i-1,j)
		aw(i,j)=0.0
		end if			
		
		if(mark_cells(i,j) .eq. 0 .and. mark_cells(i,j+2) .eq. 1) then
		sp(i,j)=sp(i,j)+an(i,j)*v1(i,j+1)
		an(i,j)=0.0
		end if	
		
		if(mark_cells(i,j) .eq. 0 .and. mark_cells(i,j-1) .eq. 1) then
		sp(i,j)=sp(i,j)+as(i,j)*v1(i,j-1)
		as(i,j)=0.0
		end if	
		
		if(mark_cells(i,j) .eq. 1 .or. (mark_cells(i,j) .eq. 0 .and. mark_cells(i,j+1) .eq. 1)) then
		aw(i,j)=0.0; ae(i,j)=0.0; as(i,j)=0.0; an(i,j)=0.0; v1(i,j)=0.0
		ap(i,j)=delv/dt
		sp(i,j)=v1(i,j)*delv/dt
		end if
		
	end do
end do

!       BOUNDARY CONDITIONS V

!Este: Neumann
aP(ei,1:ej)=aP(ei,1:ej)-aE(ei,1:ej)
aE(ei,1:ej)=0.0

!Oeste
sP(1,1:ej)=sP(1,1:ej)+aW(1,1:ej)*v1(0,1:ej)
aW(1,1:ej)=0.0

!Norte
sP(1:ei,ej)=sP(1:ei,ej)+aN(1:ei,ej)*v1(1:ei,ej+1)
aN(1:ei,ej)=0.0

!Sur
sP(1:ei,1)=sP(1:ei,1)+aS(1:ei,1)*v1(1:ei,0)
aS(1:ei,1)=0.0

DO i=1, ei
	do j=1, ej
		sn=x(i)-x(i-1)
		dv(i,j)=sn/(ap(i,j)-ae(i,j)-aw(i,j)-an(i,j)-as(i,j))
	end do
end do

CALL Gauss_Seidel2D(v1,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny,max_iter1, tolerance,residual)

!====================================================
!                  EQUATION P
!====================================================
ap=0.0;	aw=0.0;	ae=0.0;	as=0.0;	an=0.0;	sp=0.0; pp=0.0
ei=ubound(pp,1)-1
ej=ubound(pp,2)-1

do i=1, ei
	do j=1, ej
	dxw=xc(i)-xc(i-1)
	dxe=xc(i+1)-xc(i)
	
	dys=yc(j)-yc(j-1)
	dyn=yc(j+1)-yc(j)
	
	se=y(j)-y(j-1); sw=se
	sn=x(i)-x(i-1); ss=sn
	
	delv=(x(i)-x(i-1))*(y(j)-y(j-1))
	
	ue=u1(i,j)
	uw=u1(i-1,j)
	vn=v1(i,j)
	vs=v1(i,j-1)		
	
	aE(i,j)=du(i,j)*Se
	aW(i,j)=du(i-1,j)*Sw
	aN(i,j)=dv(i,j)*Sn
	aS(i,j)=dv(i,j-1)*Ss	
	
	aP(i,j)=aE(i,j)+aW(i,j)+aN(i,j)+aS(i,j)
	sP(i,j)=-ue*Se+uw*Sw-vn*Sn+vs*Ss	
	end do
end do

call Gauss_Seidel2D(pp,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny,max_iter1,tolerance,residual)

!==============================================
!   Pressure and velocity corrections
p=p+pp

do i=1,nx-1
do j=1,ny
	if(mark_cells(i,j) .eq. 0 .and. mark_cells(i+1,j) .ne. 1)then
		u1(i,j)=u1(i,j)+du(i,j)*(pp(i,j)-pp(i+1,j))
	end if
enddo
enddo

do i=1,nx
do j=1,ny-1
	if(mark_cells(i,j) .eq. 0 .and. mark_cells(i,j+1) .ne. 1)then
		v1(i,j)=v1(i,j)+dv(i,j)*(pp(i,j)-pp(i,j+1))
	end if
enddo
enddo

!Divergence and SIMPLEC counter
Div=maxval(abs(sp(1:ei,1:ej)))
c=c+1

!Update SIMPLEC velocities
u1(nx,:)=u1(nx-1,:)
v1(nx+1,:)=v1(nx,:)

enddo   !End SIMPLEC loop

!Update previous time-step velocity
u=u1
v=v1

write(*,*) it, Div, c

!==================================================
!             PARTICLE TRACKING
call interpolateToNodesUs(uc,u,nx,ny)
call interpolateToNodesVs(vc,v,nx,ny)

i=1
do while (partic(i,1) .ne. 0.0)
call distance_to_mesh(x, nx, partic(i,1), npx, 1)
call distance_to_mesh(y, ny, partic(i,2), npy, 1)
partic(i,1)=partic(i,1)+uc(npx,npy)*dt
partic(i,2)=partic(i,2)+vc(npx,npy)*dt
if(partic(i,1) .gt. 10.0) partic(i,1)=xl

i=i+1
end do

!==================================================
!           WRITE OUTPUT FILES
nn=25
if(mod(it, nn) .eq. 0.0) then
!call WriteVectorField('uv',it,uc,vc,xc,yc,nx,ny)

ei=1+48*(it-nn)/nn
ej=48+48*(it-nn)/nn

j=1
do i=ei, ej
    partic(i,1)=2.0+0.5*dx
    partic(i,2)=(-0.2+0.5*dy)+0.5*float(j-1)*dy
j=j+1
end do

i=1
do while (partic(i,1) .ne. 0.0)

if(partic(i,1) .gt. 5.0)then
i=i+1
cycle
end if

write(1,'(F7.4, A, F7.4)') partic(i,1), ' ',partic(i,2)
i=i+1
end do

write(1,*) ''
write(1,*) ''

!write(2,*) 'p [2:5][-0.5:0.5] ''partic.txt'' index ', (it-nn)/nn,' t '' '
!write(3,*)  'p [2:5][-0.5:0.5] ''uv', it,'.txt'' u 1:2:($3*0.25):($4*0.25) w vec'
end if
!==================================================

enddo	!Fin del ciclo temporal

!call interpolateToNodesUs(uc,u,nx,ny)
!call interpolateToNodesVs(vc,v,nx,ny)
!call WriteVectorField('uv',0,uc,vc,xc,yc,nx,ny)


close(1)

END PROGRAM

!============================================================================
!                   ****      Subroutines       ****
!============================================================================

Subroutine Mesh1D(xc,x,x0,xl,nx)
integer i,j,nx
real*4 x0,xl,dx
real*4 x(0:nx),xc(0:nx+1)
dx=(1.0)/dfloat(nx)

do i=0,nx
x(i)=dfloat(i)*dx
x(i)=x0+(xl-x0)*x(i)
end do

xc(0)=x(0); xc(nx+1)=x(nx)
do i=1,nx
xc(i)=(x(i)+x(i-1))*0.5
end do

End Subroutine

!==================================================
!Returns a
!Variable control: 1) devuelve nodo; 2) cara.
subroutine distance_to_mesh(x, nx, xi, nvc, control)
implicit none
integer nx, nvc, control
real*4  x0, xl, xi
real*4  x(0:nx)

x0=x(0); xl=x(nx)

if(xi<x0 .or. xi>xl)then
write(*,*) 'distance_to_mesh:', xi,'is out of the domain: ',x0, xl
stop
end if

select case (control)

case (1)    !Returns node number
    nvc=0
    do while(x(nvc)<xi)
    nvc=nvc+1
    end do

    if(xi .eq. x(nx)) nvc=nvc+1

case (2)    !Returns face number
    nvc=0
    do while(x(nvc)<xi)
    nvc=nvc+1
    end do

    if((xi .ne. x(0)).and. (xi-x(nvc-1)<x(nvc)-xi)) nvc=nvc-1 ! Returns the closest face to xi

end select

end subroutine


!==================================================
subroutine interpolateToNodesUs(uc,us,nx,ny)

implicit none
integer nx,ny
real:: us(0:nx,0:ny+1),uc(0:nx+1,0:ny+1)
integer bi,ei,bj,ej,i,j
bi=1; bj=1;ej=ny; ei=nx

! Internal points
do i=bi,ei
do j=bj-1,ej+1
	uc(I,J) = ( us(I-1,J) + us(I,J) ) * 0.5
end do 
end do

uc(bi-1,:)=us(0,:)
uc(ei+1,:)=us(nx,:)

end subroutine

!==================================================
subroutine interpolateToNodesVs(vc,vs,nx,ny)
implicit none
integer nx,ny
real:: vs(0:nx+1,0:ny),vc(0:nx+1,0:ny+1)
integer bi,ei,bj,ej,i,j
bi=1; bj=1;ej=ny; ei=nx

do i=bi-1,ei+1
do j=bj,ej
	vc(I,J) = ( vs(I,J) + vs(I,J-1) ) * 0.5
end do 
end do

vc(:,bj-1)=vs(:,0)
vc(:,ej+1)=vs(:,ny)
end subroutine

!==================================================
Subroutine WriteVectorField(Name,kx,uc,vc,xc,yc,nx,ny)
integer i,j,nx,ny,kx
real*4 uc(0:nx+1,0:ny+1),vc(0:nx+1,0:ny+1),xc(0:nx+1),yc(0:ny+1)
character*(*)Name
character*50 txt,Filename
write(txt,'(i6)')kx
txt=ADJUSTL(txt)
Filename=name//txt(1:len_trim(txt))//".txt"

open(11,file=Filename(1:len_trim(Filename)))
	do i=0,nx+1
	do j=0,ny+1
	write(11,*)xc(i),yc(j),uc(i,j),vc(i,j)
	end do
	end do
close(11)
End Subroutine

!==================================================
real*4 function calcResidual(phi,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny)
implicit none
integer ei,ej,i,j,nx,ny
real*4 phi(0:ei+1,0:ej+1),aP(nx,ny),aE(nx,ny),aW(nx,ny),aN(nx,ny),aS(nx,ny),sP(nx,ny)
real*4 acum(ei,ej),residual,NINV
NINV = 1.0 / dfloat(ei*ej)
acum=0
do i=1,ei
do j=1,ej
acum(i,j) = aE(i,j) * phi(i+1,j) + aW(i,j) * phi(i-1,j) + aN(i,j) * phi(i,j+1) + &
aS(i,j) * phi(i,j-1) + sp(i,j) - aP(i,j) * phi(i,j)
end do
end do
residual = sqrt( NINV * sum(acum * acum) )
calcResidual=residual
end function

subroutine TDMA(x,a,b,c,d,N)

implicit none
integer N,k
real*4 a(n),b(n),c(n),d(n),x(n),m

do k=2,N
m=a(k)/b(k-1)
b(k)=b(k)-m*c(k-1)
d(k)=d(k)-m*d(k-1)
end do

x(n)=d(n)/b(n)

do k=n-1,1,-1
x(k)=(d(k)-c(k)*x(k+1))/b(k)
end do

end subroutine

!==================================================
subroutine Gauss_Seidel2D(phi,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny,max_iter,tolerance,residual)

implicit none

integer ei,ej,i,j,nx,ny,count_iter,max_iter
real*4 phi(0:ei+1,0:ej+1),aP(nx,ny),aE(nx,ny),aW(nx,ny),aN(nx,ny),aS(nx,ny),sP(nx,ny)
real*4 residual,tolerance

	interface
		real*4 function calcResidual(phi,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny)
		implicit none
		integer ei,ej,i,j,nx,ny
		real*4 phi(0:ei+1,0:ej+1),aP(nx,ny),aE(nx,ny),aW(nx,ny),aN(nx,ny),aS(nx,ny),sP(nx,ny)
		real*4 acum(ei,ej),residual,NINV
		end function
	end interface
	

count_iter=0;  residual=1.0


do while((count_iter <= max_iter).and.(residual > tolerance)) 
	do i=1,ei
	do j=1,ej
	phi(i,j)=(aE(i,j) * phi(i+1,j) + aW(i,j) * phi(i-1,j)+aN(i,j) * phi(i,j+1) + &
	aS(i,j) * phi(i,j-1)+sp(i,j))/aP(i,j)
	end do
	end do
residual = calcResidual(phi,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny)
count_iter=count_iter+1
end do

end subroutine

