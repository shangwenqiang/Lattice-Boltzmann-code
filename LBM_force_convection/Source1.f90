parameter (n=1000,m=50)
real f(0:8,0:n,0:m)
real feq(0:8,0:n,0:m),rho(0:n,0:m)
real w(0:8), cx(0:8),cy(0:8)
real u(0:n,0:m), v(0:n,0:m)
real g(0:8,0:n,0:m), geq(0:8,0:n,0:m),th(0:n,0:m)
integer i
open(2,file='uvfield.plt')
open(3,file='uvely.plt')
open(4,file='vvelx.plt')
!
cx(:)=(/0.0,1.0,0.0,-1.0,0.0,1.0,-1.0,-1.0,1.0/)
cy(:)=(/0.0,0.0,1.0,0.0,-1.0,1.0,1.0,-1.0,-1.0/)
w(:)=(/4./9.,1./9.,1./9.,1./9.,1./9.,1./36.,1./36.,1./36.,1./36./)
uo=0.1
sumvelo=0.0
rhoo=5.00
dx=1.0
dy=dx
dt=1.0
tw=1.0
th=0.0
g=0.0
visco=0.05
pr=3.8
alpha=visco/pr
Re=uo*m/alpha
print *, "Re=", Re
omega=1.0/(3.*visco+0.5)
omegat=1.0/(3.*alpha+0.5)
mstep=10000
do j=0,m
do i=0,n
rho(i,j)=rhoo
u(i,j)=0.0
v(i,j)=0.0
end do
end do
do j=1,m-1
u(0,j)=uo
v(0,j)=0.0
end do
! main loop
do kk=1,mstep
call collesion(u,v,f,feq,rho,omega,w,cx,cy,n,m)
call streaming(f,n,m)
call sfbound(f,n,m,uo)
call rhouv(f,rho,u,v,cx,cy,n,m)
! ¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¨C
do j=0,m
do i=0,n
sum=0.0
do k=0,8
sum=sum+g(k,i,j)
th(i,j)=sum
end do
end do
end do
! collestion for scalar
call collt(u,v,g,geq,th,omegat,w,cx,cy,n,m)
! streaming for scalar
call streaming(g,n,m)
call gbound(g,tw,w,n,m)
print *, th(n/2,m/2),v(0,m/2),rho(0,m/2),u(n,m/2),v(n,m/2),rho(n,m/2)
END DO
! end of the main loop
call result(u,v,th,uo,n,m)
stop
end
! end of the main program
subroutine collesion(u,v,f,feq,rho,omega,w,cx,cy,n,m)
real f(0:8,0:n,0:m)
real feq(0:8,0:n,0:m),rho(0:n,0:m)
real w(0:8), cx(0:8),cy(0:8)
real u(0:n,0:m), v(0:n,0:m)
DO i=0,n
DO j=0,m
t1=u(i,j)*u(i,j)+v(i,j)*v(i,j)
DO k=0,8
t2=u(i,j)*cx(k)+v(i,j)*cy(k)
feq(k,i,j)=rho(i,j)*w(k)*(1.0+3.0*t2+4.50*t2*t2-1.50*t1)
f(k,i,j)=omega*feq(k,i,j)+(1.-omega)*f(k,i,j)
END DO
END DO
END DO
return
end
subroutine collt(u,v,g,geq,th,omegat,w,cx,cy,n,m)
real g(0:8,0:n,0:m),geq(0:8,0:n,0:m),th(0:n,0:m)
real w(0:8),cx(0:8),cy(0:8)
real u(0:n,0:m),v(0:n,0:m)
do i=0,n
do j=0,m
do k=0,8
geq(k,i,j)=th(i,j)*w(k)*(1.0+3.0*(u(i,j)*cx(k)+v(i,j)*cy(k)))
g(k,i,j)=omegat*geq(k,i,j)+(1.0-omegat)*g(k,i,j)
end do
end do
end do
return
end
subroutine streaming(f,n,m)
real f(0:8,0:n,0:m)
! streaming
DO j=0,m
DO i=n,1,-1 !RIGHT TO LEFT
f(1,i,j)=f(1,i-1,j)
END DO
DO i=0,n-1 !LEFT TO RIGHT
f(3,i,j)=f(3,i+1,j)
END DO
END DO
DO j=m,1,-1 !TOP TO BOTTOM
DO i=0,n
f(2,i,j)=f(2,i,j-1)
END DO
DO i=n,1,-1
f(5,i,j)=f(5,i-1,j-1)
END DO
DO i=0,n-1
f(6,i,j)=f(6,i+1,j-1)
END DO
END DO
DO j=0,m-1 !BOTTOM TO TOP
DO i=0,n
f(4,i,j)=f(4,i,j+1)
END DO
DO i=0,n-1
f(7,i,j)=f(7,i+1,j+1)
END DO
DO i=n,1,-1
f(8,i,j)=f(8,i-1,j+1)
END DO
END DO
return
end
subroutine sfbound(f,n,m,uo)
real f(0:8,0:n,0:m)
do j=0,m
! flow in on west boundary
rhow=(f(0,0,j)+f(2,0,j)+f(4,0,j)+2.*(f(3,0,j)+f(6,0,j)+f(7,0,j)))/(1.-uo)
f(1,0,j)=f(3,0,j)+2.*rhow*uo/3.
f(5,0,j)=f(7,0,j)+rhow*uo/6.
f(8,0,j)=f(6,0,j)+rhow*uo/6.
end do
! bounce back on south boundary
do i=0,n
f(2,i,0)=f(4,i,0)
f(5,i,0)=f(7,i,0)
f(6,i,0)=f(8,i,0)
end do
! bounce back, north boundary
do i=0,n
f(4,i,m)=f(2,i,m)
f(8,i,m)=f(6,i,m)
f(7,i,m)=f(5,i,m)
end do
! accont for open boundary condition at the outlet
do j=1,m
f(1,n,j)=2.*f(1,n-1,j)-f(1,n-2,j)
f(5,n,j)=2.*f(5,n-1,j)-f(5,n-2,j)
f(8,n,j)=2.*f(8,n-1,j)-f(8,n-2,j)
end do
return
end
subroutine gbound(g,tw,w,n,m)
real g(0:8,0:n,0:m)
real w(0:8)
! Boundary conditions
! Left boundary condition, the temperature is given, tw
do j=0,m
g(1,0,j)=-g(3,0,j)
g(5,0,j)=-g(7,0,j)
g(8,0,j)=-g(6,0,j)
end do
! Right hand boundary condition, open
do j=0,m
g(6,n,j)=2.*g(6,n-1,j)-g(6,n-2,j)
g(3,n,j)=2.*g(3,n-1,j)-g(3,n-2,j)
g(7,n,j)=2.*g(7,n-1,j)-g(7,n-2,j)
g(2,n,j)=2.*g(2,n-1,j)-g(2,n-2,j)
g(0,n,j)=2.*g(0,n-1,j)-g(0,n-2,j)
g(1,n,j)=2.*g(1,n-1,j)-g(1,n-2,j)
g(4,n,j)=2.*g(4,n-1,j)-g(4,n-2,j)
g(5,n,j)=2.*g(5,n-1,j)-g(5,n-2,j)
g(8,n,j)=2.*g(8,n-1,j)-g(8,n-2,j)
end do
! Top boundary conditions, T=0.0
do i=0,n
g(8,i,m)=tw*(w(8)+w(6))-g(6,i,m)
g(7,i,m)=tw*(w(7)+w(5))-g(5,i,m)
g(4,i,m)=tw*(w(4)+w(2))-g(2,i,m)
g(1,i,m)=tw*(w(1)+w(3))-g(3,i,m)
end do
!Bottom boundary conditions, Adiabatic
! g(1,i,0)=g(1,i,1)
! g(2,i,0)=g(2,i,1)
! g(3,i,0)=g(3,i,1)
! g(4,i,0)=g(4,i,1)
! g(5,i,0)=g(5,i,1)
! g(6,i,0)=g(6,i,1)
! g(7,i,0)=g(7,i,1)
! g(8,i,0)=g(8,i,1)
! T=0.0
do i=0,n
g(2,i,0)=tw*(w(2)+w(4))-g(4,i,0)
g(6,i,0)=tw*(w(6)+w(8))-g(8,i,0)
g(5,i,0)=tw*(w(5)+w(7))-g(7,i,0)
end do
return
end
subroutine tcalcu(g,th,n,m)
real g(0:8,0:n,0:m),th(0:n,0:m)
do j=1,m-1
do i=1,n-1
ssumt=0.0
do k=0,8
ssumt=ssumt+g(k,i,j)
end do
th(i,j)=ssumt
end do
end do
return
end
subroutine rhouv(f,rho,u,v,cx,cy,n,m)
real f(0:8,0:n,0:m),rho(0:n,0:m),u(0:n,0:m),v(0:n,0:m),cx(0:8),cy(0:8)
do j=0,m
do i=0,n
ssum=0.0
do k=0,8
ssum=ssum+f(k,i,j)
end do
rho(i,j)=ssum
end do
end do
DO i=1,n
DO j=1,m-1
usum=0.0
vsum=0.0
DO k=0,8
usum=usum+f(k,i,j)*cx(k)
vsum=vsum+f(k,i,j)*cy(k)
END DO
u(i,j)=usum/rho(i,j)
v(i,j)=vsum/rho(i,j)
END DO
END DO
do j=0,m
v(n,j)=0.0
end do
return
end
subroutine result(u,v,th,uo,n,m)
real u(0:n,0:m),v(0:n,0:m),th(0:n,0:m)
write(2,*)"VARIABLES =X, Y, U, V, T"
write(2,*)"ZONE ","I=",",",n+1,"J=",m+1,",","F=BLOCK"
do j=0,m
write(2,*)(i,i=0,n)
end do
do j=0,m
write(2,*)(j,i=0,n)
end do
do j=0,m
write(2,*)(u(i,j),i=0,n)
end do
do j=0,m
write(2,*)(v(i,j),i=0,n)
end do
do j=0,m
write(2,*)(th(i,j),i=0,n)
end do
do j=0,m
write(3,*)j/float(m),u(5,j)/uo,u(n/2,j)/uo,u(n-10,j)/uo
end do
do i=0,n
write(4,*) i/float(n),v(i,m/2)/uo
end do
return
end