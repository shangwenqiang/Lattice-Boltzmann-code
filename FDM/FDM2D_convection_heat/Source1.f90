parameter (n=100,m=100)
real dens,fo(0:n,0:m),f(0:n,0:m),x(0:n),y(0:m)
integer i,j
open(2,file='fin2drs.plt')
open(3,file='finitial.plt')
open(4,file='qfin.plt')
open(5,file='ttfind.plt')
!
dx=1.0
time=0.0
dy=1.0
dt=0.20
u=0.1
v=0.2
alpha=1.0
mstep=1000
do j=0,m
do i=0,n
fo(i,j)=0.0
end do
end do
do j=0,n
fo(0,j)=1.0
f(0,j)=1.0
fo(n,j)=0.0
f(n,j)=0.0
end do
do i=0,n
fo(i,0)=0.0
f(i,0)=0.0
! adiabatic bottom boundary
fo(i,0)=fo(i,1)
f(i,0)=f(i,1)
! zero temp top boundary
fo(i,m)=0.0
f(i,m)=0.0
end do
do kk=1,mstep
do j=1,m-1
do i=1,n-1
termx=(fo(i+1,j)+fo(i-1,j))/(dx*dx)
termy=(fo(i,j+1)+fo(i,j-1))/(dy*dy)
dd=1./(dx*dx)+1./(dy*dy)
advx=u*(fo(i,j)-fo(i-1,j))/dx
advy=v*(fo(i,j)-fo(i,j-1))/dy
advt=dt*(advx+advy)
f(i,j)=fo(i,j)+dt*alpha*(termx+termy-2.0*fo(i,j)*dd)-advt
end do
end do
!
do j=1,m-1
do i=1,n-1
fo(i,j)=f(i,j)
end do
end do
do i=0,n
! adiabatic bottom boundary
! f(i,0)=f(i,1)
! fo(i,0)=f(i,1)
! zero temp. bottom boundary
f(i,0)=0.0
fo(i,0)=0.0
end do
! adiabatic right hand boundary
do j=0,m
f(n,j)=f(n-1,j)
fo(n,j)=fo(n-1,j)
end do
time=time+dt
write(5,*)time,f(5,m/2)
end do
x(0)=0.0
do i=1,n
x(i)=x(i-1)+dx
end do
y(0)=0.0
do j=1,m
y(j)=y(j-1)+dy
end do
write(2,*)"VARIABLES =X,Y,T"
write(2,*)"ZONE ","I=",n+1,",","J=",m+1,",","F=BLOCK"
do j=0,m
write(2,*)(x(i),i=0,n)
end do
do j=0,m
write(2,*)(y(j),i=0,n)
end do
do j=0,m
write(2,*)(f(i,j),i=0,n)
end do
do j=0,m
q=f(0,j)-f(1,j)
write(4,*)y(j),q
end do
stop
end