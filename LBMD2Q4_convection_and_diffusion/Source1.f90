! LBM Code for 2-D, Advection-diffusion problems, D2Q4
parameter (n=100,m=100)
real f1(0:n,0:m),f2(0:n,0:m),f3(0:n,0:m),f4(0:n,0:m)
real feq1,feq2,feq3,feq4
real rho(0:n,0:m),x(0:n),y(0:m)
integer i
open(2,file="results.plt")
open(3,file="midplbm.plt")
!
u=0.1
v=0.2
dx=1.0
dy=dx
dt=1.0
x(0)=0
y(0)=0.0
do i=1,n
x(i)=x(i-1)+dx
end do
do j=1,m
y(j)=y(j-1)+dy
end do
ck=dx/dt
csq=ck*ck
alpha=1.0
omega=1.0/(2.*alpha/(dt*csq)+0.5)
mstep=200
do j=0,m
do i=0,n
rho(i,j)=0.0 ! initial values of the dependent variable
end do
end do
do j=0,m
do i=0,n
f1(i,j)=0.25*rho(i,j)
f2(i,j)=0.25*rho(i,j)
f3(i,j)=0.25*rho(i,j)
f4(i,j)=0.25*rho(i,j)
end do
end do
do kk=1,mstep
do j=0,m
do i=0,n
feq1=0.25*rho(i,j)*(1.0+2.0*u/ck)
feq2=0.25*rho(i,j)*(1.0-2.0*u/ck)
feq3=0.25*rho(i,j)*(1.0+2.0*v/ck)
feq4=0.25*rho(i,j)*(1.0-2.0*v/ck)
f1(i,j)=omega*feq1+(1.-omega)*f1(i,j)
f2(i,j)=omega*feq2+(1.-omega)*f2(i,j)
f3(i,j)=omega*feq3+(1.-omega)*f3(i,j)
f4(i,j)=omega*feq4+(1.-omega)*f4(i,j)
end do
end do
! Streaming
do j=0,m
do i=1,n-1
f1(n-i,j)=f1(n-i-1,j)
f2(i-1,j)=f2(i,j)
end do
end do
do i=0,n
do j=1,m-1
f3(i,m-j)=f3(i,m-j-1)
f4(i,j-1)=f4(i,j)
end do
end do
! Boundary conditions
do j=1,m
f1(0,j)=0.5-f2(0,j)
f3(0,j)=0.5-f4(0,j)
f1(n,j)=0.0
f2(n,j)=0.0
f3(n,j)=0.0
f4(n,j)=0.0
end do
do i=1,n
f1(i,m)=0.0
f2(i,m)=0.0
f3(i,m)=0.0
f4(i,m)=0.0
f1(i,0)=f1(i,1)
f2(i,0)=f2(i,1)
f3(i,0)=f3(i,1)
f4(i,0)=f4(i,1)
end do
do j=0,m
do i=0,n
rho(i,j)=f1(i,j)+f2(i,j)+f3(i,j)+f4(i,j)
end do
end do
end do ! end of the main loop
write(2,*)"VARIABLES =X,Y,T"
write(2,*)"ZONE ","I=",n+1,",","J=",m+1,",","F=BLOCK"
do j=0,m
write(2,*)(x(i),i=0,n)
end do
do j=0,m
write(2,*)(y(j),i=0,n)
end do
do j=0,m
write(2,*)(rho(i,j),i=0,n)
end do
do i=0,n
write(3,*)x(i),rho(i,m/2)
end do
stop
end