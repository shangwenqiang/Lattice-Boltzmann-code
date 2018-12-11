parameter (n=100,m=100)
real f(0:8,0:n,0:m)
real feq(0:8,0:n,0:m),rho(0:n,0:m)
real w(0:8), cx(0:8),cy(0:8)
real u(0:n,0:m), v(0:n,0:m)
real g(0:8,0:n,0:m), geq(0:8,0:n,0:m),th(0:n,0:m)
integer i
open(2,file='out1/uvfield.plt')
open(3,file='out1/uvely.plt')
open(4,file='out1/tmx.plt')
!
cx(:)=(/0.0,1.0,0.0,-1.0,0.0,1.0,-1.0,-1.0,1.0/)
cy(:)=(/0.0,0.0,1.0,0.0,-1.0,1.0,1.0,-1.0,-1.0/)
w(:)=(/4./9.,1./9.,1./9.,1./9.,1./9.,1./36.,1./36.,1./36.,1./36./)
uo=0.0
sumvelo=0.0
rhoo=6.00
dx=1.0
dy=dx
dt=1.0
tw=1.0
th=0.0
ra=1.0e5
pr=0.71
visco=0.02
alpha=visco/pr
pr=visco/alpha
gbeta=ra*visco*alpha/(float(m*m*m))
!Re=uo*m/alpha
!print *, "Re=", Re
omega=1.0/(3.*visco+0.5)
omegat=1.0/(3.*alpha+0.5)
mstep=150000
do j=0,m
do i=0,n
rho(i,j)=rhoo
u(i,j)=0.0
v(i,j)=0.0
end do
end do
!do i=0,n
!u(i,m)=uo
!v(i,m)=0.0
!end do
! main loop
do kk=1,mstep
call collesion(u,v,f,feq,rho,omega,w,cx,cy,n,m,th,gbeta)
call streaming(f,n,m)
call bounceb(f,n,m)
call rhouv(f,rho,u,v,cx,cy,n,m)
! ¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¨C
! collestion for scalar
call collt(u,v,g,geq,th,omegat,w,cx,cy,n,m)
! streaming for scalar
call streaming(g,n,m)
call gbound(g,tw,w,n,m)
call tcalcu(g,th,n,m)
print *, th(n/2,m/2),rho(0,m/2),u(n/2,m/2),v(n/2,m/2),rho(n,m/2)
END DO
! end of the main loop
call result(u,v,rho,th,uo,n,m,ra)
stop
end
! end of the main program
subroutine collesion(u,v,f,feq,rho,omega,w,cx,cy,n,m,th,gbeta)
real f(0:8,0:n,0:m)
real feq(0:8,0:n,0:m),rho(0:n,0:m)
real w(0:8), cx(0:8),cy(0:8)
real u(0:n,0:m), v(0:n,0:m)
real th(0:n,0:m)
tref=0.50
DO i=0,n
DO j=0,m
t1=u(i,j)*u(i,j)+v(i,j)*v(i,j)
DO k=0,8
t2=u(i,j)*cx(k)+v(i,j)*cy(k)
force=3.*w(k)*gbeta*(th(i,j)-tref)*cy(k)*rho(i,j)
if(i.eq.0.or.i.eq.n) force =0.0
if(j.eq.0.or.j.eq.m) force =0.0
feq(k,i,j)=rho(i,j)*w(k)*(1.0+3.0*t2+4.50*t2*t2-1.50*t1)
f(k,i,j)=omega*feq(k,i,j)+(1.-omega)*f(k,i,j)+force
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
subroutine bounceb(f,n,m)
real f(0:8,0:n,0:m)
do j=0,m
!west boundary
f(1,0,j)=f(3,0,j)
f(5,0,j)=f(7,0,j)
f(8,0,j)=f(6,0,j)
!east boundary
f(3,n,j)=f(1,n,j)
f(7,n,j)=f(5,n,j)
f(6,n,j)=f(8,n,j)
end do
do i=0,n
!south boundary
f(2,i,0)=f(4,i,0)
f(5,i,0)=f(7,i,0)
f(6,i,0)=f(8,i,0)
!north boundary
f(4,i,m)=f(2,i,m)
f(8,i,m)=f(6,i,m)
f(7,i,m)=f(5,i,m)
end do
return
end
subroutine gbound(g,tw,w,n,m)
real g(0:8,0:n,0:m)
real w(0:8),tw
! Boundary conditions
! Bottom boundary condition, T=1.
do i=0,n
g(2,i,0)=tw*(w(2)+w(4))-g(4,i,0)
g(5,i,0)=tw*(w(5)+w(7))-g(7,i,0)
g(6,i,0)=tw*(w(6)+w(8))-g(8,i,0)
end do
! Top boundary condition, T=0.
do i=0,n
g(4,i,m)=-g(2,i,m)
g(7,i,m)=-g(5,i,m)
g(8,i,m)=-g(6,i,m)
end do
do j=0,m
    do k=0,8
g(k,n,j)=g(k,n-1,j) ! Top boundary conditions, Adiabatic
g(k,0,j)=g(k,1,j)   ! Bottom boundary conditions, Adiabatic
end do
end do
return
end
subroutine tcalcu(g,th,n,m)
real g(0:8,0:n,0:m),th(0:n,0:m)
do j=0,m
do i=0,n
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
DO i=0,n
DO j=0,m
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
return
end
subroutine result(u,v,rho,th,uo,n,m,ra)
real u(0:n,0:m),v(0:n,0:m),th(0:n,0:m)
real strf(0:n,0:m),rho(0:n,0:m)
open(5,file='out1/streamt.plt')
open(6,file='out1/nuav.plt')
! streamfunction calculations
strf(0,0)=0.
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
write(2,*)"VARIABLES =X, Y, U, V, RHO"
write(2,*)"ZONE ","I=",",",n+1,"J=",m+1,",","F=BLOCK"
do j=0,m
write(2,*)(i/float(m),i=0,n)
end do
do j=0,m
write(2,*)(j/float(m),i=0,n)
end do
do j=0,m
write(2,*)(u(i,j),i=0,n)
end do
do j=0,m
write(2,*)(v(i,j),i=0,n)
end do
do j=0,m
write(2,*)(rho(i,j),i=0,n)
end do
do j=0,m
write(3,*)j/float(m),u(5,j)/uo,u(n/2,j)/uo,u(n-10,j)/uo
end do
do i=0,n
write(4,*) i/float(n),th(i,m/2)
end do
write(5,*)"VARIABLES =X, Y, S, T"
write(5,*)"ZONE ","I=",",",n+1,"J=",m+1,",","F=BLOCK"
do j=0,m
write(5,*)(i/float(m),i=0,n)
end do
do j=0,m
write(5,*)(j/float(m),i=0,n)
end do
do j=0,m
write(5,*)(strf(i,j),i=0,n)
end do
do j=0,m
write(5,*)(th(i,j),i=0,n)
end do
! Nusselt number Calculation
snul=0.0
snur=0.0
do i=0,n
rnul=(th(i,0)-th(i,1))*float(m)
rnur=(th(i,m-1)-th(i,m))*float(m)
snul=snul+rnul
snur=snur+rnur
write(5,*)i/float(n),rnul,rnur
end do
avnl=snul/float(n)
avnr=snur/float(n)
write(6,*)ra,avnl,avnr
return
end