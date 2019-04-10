!Two-dimensional modeling of free convection in sloping porous enclosures.
!Stream function - Fixed Point iteration
!Author: Fernando J. Guerrero-Martínez. 2014


program DarcyBoussinesq
implicit none

integer i, j, k, ii, jj, kk, nx, ny, itmax, max_iter1, max_iter2, it, c, n, n_mean
real*4  dx, dy, dxw, dys, dxe, dyn, dt, x0, xl, y0, yl, d1
real*4 delv, tolerance, ue, uw, vn, vs, residual, div, norm_inf
real*4 se, sw, sn, ss, lambda, Ra, Nu_g, pi, alpha
real*4 deltaTemp, ae1, aw1, an1, as1, ap1

real(4), allocatable, dimension(:,:) ::  T, TT, convt, f
real(4), allocatable, dimension(:,:) :: ae, aw, as, an, ap, sp, u, v, uc, vc
real(4), allocatable, dimension(:) :: xc, yc, x, y, steady_T, Nu

nx=300
ny=100

max_iter1=2000 !solver
max_iter2=500 !fixed point

x0=0.0
xl=3.0
y0=0.0
yl=1.0

tolerance=5.0e-6

pi=4.*ATAN(1.)
alpha=10.0*(pi/180.0)
Ra=100.0

!time step
dt=2.0E-4

n=1 !counter to define steady state
n_mean=22000 !number of norms infinite of T to be averaged

allocate (T(0:nx+1,0:ny+1), TT(0:nx+1,0:ny+1), convt(0:nx+1,0:ny+1))
allocate (f(0:nx,0:ny))
allocate (ap(nx,ny), ae(nx,ny), aw(nx,ny), as(nx,ny), an(nx,ny), sp(nx,ny))
allocate (u(0:nx,0:ny+1), v(0:nx+1,0:ny))
allocate (uc(0:nx+1,0:ny+1), vc(0:nx+1,0:ny+1))
allocate (xc(0:nx+1),yc(0:ny+1),x(0:nx),y(0:ny), Nu(nx), steady_T(n_mean))


call mesh1d(xc, x, x0, xl, nx)
call mesh1d(yc, y, y0, yl, ny)

dx=(xl-x0)/float(nx)
dy=(yl-y0)/float(ny)
delv=dx*dy

se=dy
sw=dy
sn=dx
ss=dx

lambda=0.9
itmax=100000

!================================================
!
norm_inf=1.0 !norm infinite to define steady state

!================================================
!		Initial conditions
u=0.0; v=0.0; uc=0.0; vc=0.0;
T=0.0; TT=0.0; steady_T=0.0

T(:,0)=1.0
T(:,ny+1)=0.0

!TT is the temperature at present fixed point iteration. T is the temperature at the previous time step
TT=T
!Convergence of TT
convt=1.0


!Stream function is zero at boundaries and as initial condition into the cavity
f=0.0
!================================================
!ciclo temporal
do it = 1, itmax
c=1

!FIXED POINT: Maxval of the difference between the last and the present iterations is stored in deltaTemp
deltaTemp=1.0
do while ((c < max_iter2) .and. (deltaTemp > 5.0e-6))
!================================================
!		    EQUATION T
!================================================
ap=0.0; aw=0.0; ae=0.0; as=0.0; an=0.0; sp=0.0; d1=0.0

do i=1, nx
	do j=1, ny
	dxe=xc(i+1)-xc(i); dxw=xc(i)-xc(i-1)
	dyn=yc(j+1)-yc(j); dys=yc(j)-yc(j-1)
	
	ue=u(i,j)
	uw=u(i-1,j)
	vn=v(i,j)
	vs=v(i,j-1)

	ae1=se/dxe-ue*se/2.0
	aw1=sw/dxw+uw*sw/2.0
	an1=sn/dyn-vn*sn/2.0
	as1=ss/dys+vs*ss/2.0
	
		!Neuman conditions
		if(i .eq. 1) aw1=0.0! West boundary, adiabatic
        if(i .eq. nx) ae1=0.0! East boundary, adiabatic
			
	ap1=ae1+aw1+an1+as1+delv/dt
	d1=lambda*(ap1*TT(i,j)-ae1*TT(i+1,j)-aw1*TT(i-1,j)-an1*TT(i,j+1)-as1*TT(i,j-1)-delv/dt*T(i,j))
!===
	ae(i,j)=se/dxe
	aw(i,j)=sw/dxw
	an(i,j)=sn/dyn
	as(i,j)=ss/dys
	
		!Neuman conditions
		if(i .eq. 1) aw(i,j)=0.0! West boundary, adiabatic
		if(i .eq. nx) ae(i,j)=0.0! East boundary, adiabatic
	
	ap(i,j)=ae(i,j)+aw(i,j)+an(i,j)+as(i,j)+delv/dt
	sp(i,j)=ap(i,j)*TT(i,j)-ae(i,j)*TT(i+1,j)-aw(i,j)*TT(i-1,j)-an(i,j)*TT(i,j+1)-as(i,j)*TT(i,j-1)-d1
	end do
end do

!BOUNDARY CONDITIONS

!Top boundary: dirichlet
sp(:,ny)=sp(:,ny)+an(:,ny)*T(1:nx,ny+1)
an(:,ny)=0.0

!Bottom boundary: dirichlet
sp(:,1)=sp(:,1)+as(:,1)*T(1:nx,0)
as(:,1)=0.0

call SOR2D(1.3,TT,nx,ny,aP,aE,aW,aN,aS,sP,nx,ny,max_iter1,tolerance,residual)
!call Gauss_Seidel2D(TT,nx,ny,aP,aE,aW,aN,aS,sP,nx,ny,max_iter1,tolerance,residual)

!================================================
!		    EQUATION f
!================================================
ap=0.0;	aw=0.0;	ae=0.0;	as=0.0;	an=0.0;	sp=0.0

do i=1, nx
	do j=1, ny-1
	dxe=x(i+1)-x(i); dxw=x(i)-x(i-1)
	dyn=y(j+1)-y(j); dys=y(j)-y(j-1)
	
	ae(i,j)=se/dxe
	aw(i,j)=sw/dxw
	an(i,j)=sn/dyn
	as(i,j)=ss/dys
	ap(i,j)=ae(i,j)+aw(i,j)+an(i,j)+as(i,j)
	
	sp(i,j)=Ra*(cos(alpha)*(0.5*(TT(i+1,j)+TT(i+1,j+1))-0.5*(TT(i,j)+TT(i,j+1)))*se &
	           -sin(alpha)*(0.5*(TT(i,j+1)+TT(i+1,j+1))-0.5*(TT(i,j)+TT(i+1,j)))*sn) !Tn*sn-Ts*ss
	end do
end do

!East
sp(nx-1,:)=sp(nx-1,:)+ae(nx-1,:)*f(nx,1:ny)
ae(nx-1,:)=0.0

!West
sp(1,:)=sp(1,:)+aw(1,:)*f(0,1:ny)
aw(1,:)=0.0

!North
sp(:,ny-1)=sp(:,ny-1)+an(:,ny-1)*f(1:nx,ny)
an(:,ny-1)=0.0

!South
sp(:,1)=sp(:,1)+as(:,1)*f(1:nx,0)
as(:,1)=0.0

call SOR2D(1.95,f,nx-1,ny-1,aP,aE,aW,aN,aS,sP,nx,ny,max_iter1,tolerance,residual)
!call Gauss_Seidel2D(f,nx-1,ny-1,aP,aE,aW,aN,aS,sP,nx,ny,max_iter1,tolerance,residual)

deltaTemp=maxval(abs(convt(1:nx,1:ny)-TT(1:nx,1:ny)))
convt(1:nx,1:ny)=TT(1:nx,1:ny)

!================================================
!       Velocity field
do i=1, nx-1
do j=1, ny
u(i,j)=(f(i,j)-f(i,j-1))/dy
end do
end do

do i=1, nx
do j=1, ny-1
v(i,j)=-(f(i,j)-f(i-1,j))/dx
end do
end do

!================================================

c=c+1
end do   !end of fixed point
!================================================
! Caculate mean norm infinite
if(n .eq. n_mean) n=1
steady_T(n)=maxval(abs(T(1:nx,1:ny)-TT(1:nx,1:ny)))
n=n+1

if(it .gt. n_mean)then
norm_inf=0.0
norm_inf=sum(steady_T)/float(n_mean)
end if
!================================================
!Update T 
T=TT

!================================================
! Nusselt Global on hot surface
nu=0.0
do i=1, nx
nu(i)=-2.0*(T(i,1)-T(i,0))/dy
end do

Nu_g=0.0
do i=1, nx
Nu_g=Nu_g+Nu(i)*dx
end do

!================================================
!Divergence
sp=0.0
do i=1, nx
	do j=1, ny
	sp(i,j)=u(i,j)-u(i-1,j)+v(i,j)-v(i,j-1)
	enddo
enddo

write(*,*) it, c, deltaTemp,  maxval(abs(sp)), norm_inf

!================================================
!	WRITE OUTPUT AND END
!if(mod(it, 1000) .eq. 0.0) then


!================================================

if(norm_inf .lt. 5.0e-6) exit
end do	!fin del ciclo temporal


!================================================
!	WRITE OUTPUT AND END

!	Local Nusselt number
open(3, file='nu2d_a10.txt')
do i=1, nx
write(3,*) xc(i), nu(i)
end do
close(3)

!open(3, file='nu_global_a40.txt')
!write(3,*)  Nu_g, dt*float(it)
!close(3)




T(0,:)=T(1,:); T(nx+1,:)=T(nx,:)
call WriteScalarField('T',0,T,xc,yc,nx,ny)
call WriteScalarField('f',0,f,xc,yc,nx-1,ny-1)

call interpolateToNodesUs(uc,u,nx,ny)
call interpolateToNodesVs(vc,v,nx,ny)
call WriteVectorField('uv',0,uc,vc,xc,yc,nx,ny)


!================================================

!do i=0, ny+1
!write(1,*) yc(i), T(10,i)
!end do

call computing_time
end program

!================================================
!================================================
!================================================

subroutine mesh1d(xc,x,x0,xl,nx)
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

end subroutine

!================================================

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

!================================================
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

!================================================

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

!================================================

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

!================================================

subroutine SOR2D(param,phi,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny,max_iter,tolerance,residual)

implicit none

integer ei,ej,i,j,nx,ny,count_iter,max_iter
real*4 phi(0:ei+1,0:ej+1),phi0(0:ei+1,0:ej+1),aP(nx,ny),aE(nx,ny),aW(nx,ny),aN(nx,ny),aS(nx,ny),sP(nx,ny)
real*4 residual,tolerance, param

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
phi0=phi
do i=1,ei
do j=1,ej
phi(i,j)=(1.0-param)*phi0(i,j)+param*(aE(i,j) * phi(i+1,j) + aW(i,j) * phi(i-1,j)+aN(i,j) * phi(i,j+1) + &
aS(i,j) * phi(i,j-1)+sp(i,j))/aP(i,j)
end do
end do
residual = calcResidual(phi,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny)
count_iter=count_iter+1
end do
!write(*,*) '---> Iterations', count_iter, residual

end subroutine

!================================================

Subroutine WriteScalarField(Name,kx,T,xc,yc,nx,ny)
integer i,j,nx,ny,kx
real*4 T(0:nx+1,0:ny+1),xc(0:nx+1),yc(0:ny+1)
character*(*)Name
character*50 txt,Filename
write(txt,'(i6)')kx
txt=ADJUSTL(txt)
Filename=name//txt(1:len_trim(txt))//".txt"

open(10,file=Filename(1:len_trim(Filename)))
do i=0,nx+1
do j=0,ny+1
write(10,*)xc(i),yc(j),T(i,j)
end do
write(10,*)" "
end do
close(10)
End Subroutine

!********************************************************************
!********************************************************************

Subroutine computing_time
real*4 cputime, min, hou, sec

call cpu_time(cputime)

if(cputime .lt. 60.0)then
write(*,*)'computation time:', cputime, 'seconds'
end if

if((cputime .ge. 60.0) .and. (cputime .lt. 3600.0))then
min=int(cputime/60.0)
sec=(cputime/60.0-int(cputime/60.0))*60.0
write(*,*)'computation time:', min, 'minutes', sec, 'seconds'
end if

if(cputime .ge. 3600.0)then
hou=int(cputime/3600.0)
min=int((cputime/3600.0-int(cputime/3600.0))*60.0)
sec=(((cputime/3600.0-int(cputime/3600.0))*60.0)-int((cputime/3600.0-int(cputime/3600.0))*60.0))*60.0
write(*,*)'computation time:', hou, 'hours,', min, 'minutes,', sec, 'seconds'
end if

end Subroutine
