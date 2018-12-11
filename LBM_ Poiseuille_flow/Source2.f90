! computer code for planes
parameter (n=1000,m=40)
real f(0:8,0:n,0:m),ff(0:8,0:n,0:m)
real feq(0:8,0:n,0:m),rho(0:n,0:m)
real w(0:8), cx(0:8),cy(0:8)
real u(0:n,0:m), v(0:n,0:m)
!integer ::xc=500,yc=20,chang=100,kuang=30
integer i,j,k!,s(0:n,0:m)
!
cx(:)=(/0.0,1.0,0.0,-1.0,0.0,1.0,-1.0,-1.0,1.0/)
cy(:)=(/0.0,0.0,1.0,0.0,-1.0,1.0,1.0,-1.0,-1.0/)
w(:)=(/4./9.,1./9.,1./9.,1./9.,1./9.,1./36.,1./36.,1./36.,1./36./)
uo=0.1
rhoo=1.0
dx=1.0
dy=dx
dt=1.0
alpha=0.02
re=uo*m/alpha
print *, "re=", re
omega=1.0/(2.*alpha+0.5)
mstep=40000
do j=0,m
do i=0,n
rho(i,j)=rhoo
u(i,j)=0.0
v(i,j)=0.0
!s(i,j)=0
end do
    end do
!do i=0,n
!    do j=0,m
!        if((i.gt.(xc-int(chang/2))).or.(i.lt.(xc+int(chang/2))).or.(j.gt.(yc-int(kuang/2))).or.(j.lt.(yc+int(kuang/2)))) then
!            s(i,j)=1
!        end if
!        if(i.eq.0.or.i.eq.n.or.j.eq.0.or.j.eq.m) then
!            s(i,j)=1
!        end if
!    end do
!end do
do j=1,m-1
    u(1,j)=uo
end do
    
do i=0,n
do j=0,m
t1=u(i,j)*u(i,j)+v(i,j)*v(i,j)
do k=0,8
t2=u(i,j)*cx(k)+v(i,j)*cy(k)
feq(k,i,j)=rho(i,j)*w(k)*(1.0+3.0*t2+4.50*t2*t2-1.50*t1)
f(k,i,j)=feq(k,i,j)
end do
end do
end do
! main loop
do kk=1,mstep
!call collesion(u,v,f,ff,feq,rho,omega,w,cx,cy,n,m,s)
call collesion(u,v,f,ff,feq,rho,omega,w,cx,cy,n,m)
call streaming(f,ff,n,m,cx,cy)
! ¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¨C
!call sfbound(f,ff,n,m,uo,s)
call sfbound(f,ff,n,m,uo)
!call rhouv(f,rho,u,v,cx,cy,n,m,s)
call rhouv(f,rho,u,v,cx,cy,n,m)
print*, kk,rho(n/2,m/2),u(n/2,m/2),v(n/2,m/2)
end do
! end of the main loop
!call result(u,v,rho,n,m,s)
call result(u,v,rho,n,m)
stop
end
! end of the main program
subroutine collesion(u,v,f,ff,feq,rho,omega,w,cx,cy,n,m)
real f(0:8,0:n,0:m),ff(0:8,0:n,0:m)
real feq(0:8,0:n,0:m),rho(0:n,0:m)
real w(0:8), cx(0:8),cy(0:8)
real u(0:n,0:m), v(0:n,0:m)
real omega
!integer s(0:n,0:m)
do i=1,n-1
do j=1,m-1
    !if (s(i,j).eq.0) then
t1=u(i,j)*u(i,j)+v(i,j)*v(i,j)
do k=0,8
t2=u(i,j)*cx(k)+v(i,j)*cy(k)
feq(k,i,j)=rho(i,j)*w(k)*(1.0+3.0*t2+4.50*t2*t2-1.50*t1)
ff(k,i,j)=omega*feq(k,i,j)+(1.-omega)*f(k,i,j)
end do
!end if
end do
end do
return
end
subroutine streaming(f,ff,n,m,cx,cy)
real f(0:8,0:n,0:m),ff(0:8,0:n,0:m)
real cx(0:8),cy(0:8)
integer ip,jp
! streaming
do j=1,m-1
    do i=1,n-1
        do k=0,8
        ip=i-int(cx(k))
        jp=j-int(cy(k))
        f(k,i,j)=ff(k,ip,jp)
end do
end do
end do
return
end
subroutine sfbound(f,ff,n,m,uo)
real f(0:8,0:n,0:m),ff(0:8,0:n,0:m)
!integer s(0:n,0:m)
do j=1,m-1
! bounce back on west boundary
rhow=(f(0,1,j)+f(2,1,j)+f(4,1,j)+2*(f(3,1,j)+f(6,1,j)+f(7,1,j)))/(1.0-uo)
f(1,1,j)=f(3,1,j)+2*rhow*uo/3.0
f(5,1,j)=f(7,1,j)+rhow*uo/6.0!-0.5*(f(2,1,j)-f(4,1,j))
f(8,1,j)=f(6,1,j)+rhow*uo/6.0!+0.5*(f(2,1,j)-f(4,1,j))
f(3,n-1,j)=2*f(3,n-2,j)-f(3,n-3,j)
f(7,n-1,j)=2*f(7,n-2,j)-f(7,n-3,j)
f(6,n-1,j)=2*f(6,n-2,j)-f(6,n-3,j)
end do
! bounce back on south boundary
do i=1,n-1
    !do j=1,m-1
    !    if(s(i,j).eq.1) then
f(2,i,1)=f(4,i,1)
f(5,i,1)=f(7,i,1)
f(6,i,1)=f(8,i,1)
f(4,i,m-1)=f(2,i,m-1)
f(8,i,m-1)=f(6,i,m-1)
f(7,i,m-1)=f(5,i,m-1)
!ff(1,i,j)=f(3,i,j)
!ff(3,i,j)=f(1,i,j)
!ff(2,i,j)=f(4,i,j)
!ff(4,i,j)=f(2,i,j)
!ff(5,i,j)=f(7,i,j)
!ff(7,i,j)=f(5,i,j)
!ff(6,i,j)=f(8,i,j)
!ff(8,i,j)=f(6,i,j)
!        end if
!    end do
end do
return
end
subroutine rhouv(f,rho,u,v,cx,cy,n,m)
real f(0:8,0:n,0:m),rho(0:n,0:m),u(0:n,0:m),v(0:n,0:m),cx(0:8),cy(0:8)
!integer s(0:n,0:m)
do j=1,m-1
do i=1,n-1
ssum=0.0
usum=0.0
vsum=0.0
do k=0,8
ssum=ssum+f(k,i,j)
usum=usum+f(k,i,j)*cx(k)
vsum=vsum+f(k,i,j)*cy(k)
end do
rho(i,j)=ssum
u(i,j)=usum/rho(i,j)
v(i,j)=vsum/rho(i,j)
!if(s(i,j).eq.1)then 
!    u(i,j)=0
!    v(i,j)=0
!    rho(i,j)=rhoo
!    endif 
end do
end do
return
end
subroutine result(u,v,rho,n,m)
real u(0:n,0:m),v(0:n,0:m)
real rho(0:n,0:m)
integer i,j!,s(0:n,0:m)
open(2,file='out/uvfield.plt')
write(2,*)"VARIABLES=X,Y,U,V,RHO"
write(2,*)"ZONE I=",n-1,",J=",m-1,",","F=BLOCK"
!do j=1,m-1
!    do i=1,n-1
!        write(2,*) i,j,u(i,j),v(i,j),rho(i,j),s(i,j)
!end do
!end do
do j=1,m-1
write(2,*)(i,i=1,n-1)
end do
do j=1,m-1
write(2,*)(j,i=1,n-1)
end do
do j=1,m-1
write(2,*)(u(i,j),i=1,n-1)
end do
do j=1,m-1
write(2,*)(v(i,j),i=1,n-1)
end do
do j=1,m-1
write(2,*)(rho(i,j),i=1,n-1)
end do
!do j=1,m-1
!write(2,*)(s(i,j),i=1,n-1)
!end do
return
end
!============end of the program