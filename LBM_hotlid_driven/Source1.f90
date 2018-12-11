!complete computer code is given below for a heated moving lid cavity
!======================================================
parameter (n=100,m=100)
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
uo=0.2
sumvelo=0.0
rhoo=5.00
dx=1.0
dy=dx
dt=1.0
tw=1.0
th=0.0
g=0.0
visco=0.02
pr=0.71
alpha=visco/pr
re=uo*m/alpha
print *, "re=", re
omega=1.0/(2.*visco+0.5)
omegat=1.0/(2.*alpha+0.5)
mstep=20000
do j=0,m
do i=0,n
rho(i,j)=rhoo
u(i,j)=0.0
v(i,j)=0.0
end do
end do
do i=1,n-1
u(i,m)=uo
v(i,m)=0.0
end do
! main loop
do kk=1,mstep
call collesion(u,v,f,feq,rho,omega,w,cx,cy,n,m)
call streaming(f,n,m)
call sfbound(f,n,m,uo)
call rhouv(f,rho,u,v,cx,cy,n,m)
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
end do
! end of the main loop
call result(u,v,rho,th,uo,n,m)
stop
end
! end of the main program
subroutine collesion(u,v,f,feq,rho,omega,w,cx,cy,n,m)
real f(0:8,0:n,0:m)
real feq(0:8,0:n,0:m),rho(0:n,0:m)
real w(0:8), cx(0:8),cy(0:8)
real u(0:n,0:m), v(0:n,0:m)
do i=0,n
do j=0,m
t1=u(i,j)*u(i,j)+v(i,j)*v(i,j)
do k=0,8
t2=u(i,j)*cx(k)+v(i,j)*cy(k)
feq(k,i,j)=rho(i,j)*w(k)*(1.0+3.0*t2+4.50*t2*t2-1.50*t1)
f(k,i,j)=omega*feq(k,i,j)+(1.-omega)*f(k,i,j)
end do
end do
end do
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
do j=0,m
do i=n,1,-1 !right to left
f(1,i,j)=f(1,i-1,j)
end do
do i=0,n-1 !left to right
f(3,i,j)=f(3,i+1,j)
end do
end do
do j=m,1,-1 !top to bottom
do i=0,n
f(2,i,j)=f(2,i,j-1)
end do
do i=n,1,-1
f(5,i,j)=f(5,i-1,j-1)
end do
do i=0,n-1
f(6,i,j)=f(6,i+1,j-1)
end do
end do
do j=0,m-1 !bottom to top
do i=0,n
f(4,i,j)=f(4,i,j+1)
end do
do i=0,n-1
f(7,i,j)=f(7,i+1,j+1)
end do
do i=n,1,-1
f(8,i,j)=f(8,i-1,j+1)
end do
end do
return
end
subroutine sfbound(f,n,m,uo)
real f(0:8,0:n,0:m)
do j=0,m
! bounce back on west boundary
f(1,0,j)=f(3,0,j)
f(5,0,j)=f(7,0,j)
f(8,0,j)=f(6,0,j)
! bounce back on east boundary
f(3,n,j)=f(1,n,j)
f(7,n,j)=f(5,n,j)
f(6,n,j)=f(8,n,j)
end do
! bounce back on south boundary
do i=0,n
f(2,i,0)=f(4,i,0)
f(5,i,0)=f(7,i,0)
f(6,i,0)=f(8,i,0)
end do
! moving lid, north boundary
do i=1,n-1
rhon=f(0,i,m)+f(1,i,m)+f(3,i,m)+2.*(f(2,i,m)+f(6,i,m)+f(5,i,m))
f(4,i,m)=f(2,i,m)
f(8,i,m)=f(6,i,m)+rhon*uo/6.0
f(7,i,m)=f(5,i,m)-rhon*uo/6.0
end do
return
end
subroutine gbound(g,tw,w,n,m)
real g(0:8,0:n,0:m)
real w(0:8)
! boundary conditions
! west boundary condition, t=0.
do j=0,m
g(1,0,j)=-g(3,0,j)
g(5,0,j)=-g(7,0,j)
g(8,0,j)=-g(6,0,j)
end do
! east boundary condition, t=0.
do j=0,m
g(6,n,j)=-g(8,n,j)
g(3,n,j)=-g(1,n,j)
g(7,n,j)=-g(5,n,j)
g(2,n,j)=-g(4,n,j)
g(0,n,j)=0.0
end do
! top boundary conditions, t=tw=1.0
do i=0,n
g(8,i,m)=tw*(w(8)+w(6))-g(6,i,m)
g(7,i,m)=tw*(w(7)+w(5))-g(5,i,m)
g(4,i,m)=tw*(w(4)+w(2))-g(2,i,m)
g(1,i,m)=tw*(w(1)+w(3))-g(3,i,m)
end do
!bottom boundary conditions, adiabatic
do i=0,n
g(1,i,0)=g(1,i,1)
g(2,i,0)=g(2,i,1)
g(3,i,0)=g(3,i,1)
g(4,i,0)=g(4,i,1)
g(5,i,0)=g(5,i,1)
g(6,i,0)=g(6,i,1)
g(7,i,0)=g(7,i,1)
g(8,i,0)=g(8,i,1)
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
do i=1,n
do j=1,m-1
usum=0.0
vsum=0.0
do k=0,8
usum=usum+f(k,i,j)*cx(k)
vsum=vsum+f(k,i,j)*cy(k)
end do
u(i,j)=usum/rho(i,j)
v(i,j)=vsum/rho(i,j)
end do
end do
return
end
subroutine result(u,v,rho,th,uo,n,m)
real u(0:n,0:m),v(0:n,0:m),th(0:n,0:m)
real rho(0:n,0:m),strf(0:n,0:m)
open(5, file='streamf.plt')
open(7,file='tprof.plt')
! streamfunction calculations
strf(0,0)=0.0
do i=1,n
rhoav=0.5*(rho(i-1,0)+rho(i,0))
if(i.ne.0) strf(i,0)=strf(i-1,0)-rhoav*0.5*(v(i-1,0)+v(i,0))
end do
do i=0,n
do j=1,m
rhom=0.5*(rho(i,j)+rho(i,j-1))
strf(i,j)=strf(i,j-1)+rhom*0.5*(u(i,j-1)+u(i,j))
end do
end do
! ¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¨C
write(2,*)"VARIABLES =X, Y, U, V, T"
write(2,*)"ZONE ","I=",n+1,",","J=",m+1,",","F=BLOCK"
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
write(5,*)"VARIABLES=X,Y,ST"
write(5,*)"ZONE ","I=",n+1,",","J=",m+1,",","F=BLOCK"
do j=0,m
write(5,*)(i,i=0,n)
end do
do j=0,m
write(5,*)(j,i=0,n)
end do
do j=0,m
write(5,*)(strf(i,j),i=0,n)
end do
do j=0,m
write(3,*)j/float(m),u(5,j)/uo,u(n/2,j)/uo,u(n-10,j)/uo
write(7,*)j/float(m),th(n/4,j),th(n/2,j),th(3*n/4,j)
end do
do i=0,n
write(4,*) i/float(n),v(i,m/2)/uo
end do
return
end
!======================================================