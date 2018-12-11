!This code works fine and produce results for lid driven cavity,MRT code.
program main
parameter(n=100,m=100)
real f(0:8,0:n,0:m)
real feq(0:8,0:n,0:m),rho(0:n,0:m)
real w(0:8),cx(0:8),cy(0:8)
real u(0:n,0:m),v(0:n,0:m)
integer i
real tminv(0:8,0:8),sm(0:8),tm(0:8,0:8),stmiv(0:8,0:8)
real ev(0:8,0:8)
open(2,file='uvfield.plt')
open(3,file='uvely.plt')
open(4,file='vvelx.plt')
open(8,file='timeu.plt')
open(10,file='tmat.plt')
w(:)=(/4./9.,1./9.,1./9.,1./9.,1./9.,1./36.,1./36.,1./36.,1./36./)
cx(:)=(/0.,1.,0.,-1.,0.,1.,-1.,-1.,1./)
cy(:)=(/0.,0.,1.,0.,-1.,1.,1.,-1.,-1./)
tm(0,:)=(/1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0/)!¾ØÕóM
tm(1,:)=(/-4.,-1.,-1.,-1.,-1.,2.0,2.0,2.0,2.0/)
tm(2,:)=(/4.0,-2.,-2.,-2.,-2.,1.0,1.0,1.0,1.0/)
tm(3,:)=(/0.0,1.0,0.0,-1.,0.0,1.0,-1.,-1.,1.0/)
tm(4,:)=(/0.0,-2.,0.0,2.0,0.0,1.0,-1.,-1.,1.0/)
tm(5,:)=(/0.0,0.0,1.0,0.0,-1.,1.0,1.0,-1.,-1./)
tm(6,:)=(/0.0,0.0,-2.,0.0,2.0,1.0,1.0,-1.,-1./)
tm(7,:)=(/0.0,1.0,-1.,1.0,-1.,0.0,0.0,0.0,0.0/)
tm(8,:)=(/0.0,0.0,0.0,0.0,0.0,1.0,-1.,1.0,-1./)
a1=1./36.
tminv(0,:)=(/4.0*a1,-4.0*a1,4.0*a1,0.0,0.0,0.0,0.0,0.0,0.0/)!MµÄÄæ¾ØÕó
tminv(1,:)=(/4.0*a1,-a1,-2.*a1,6.*a1,-6.*a1,0.,0.,9.*a1,0./)
tminv(2,:)=(/4.*a1,-a1,-2.*a1,0.,0.,6.*a1,-6.*a1,-9.*a1,0./)
tminv(3,:)=(/4.*a1,-a1,-2.*a1,-6.*a1,6.*a1,0.0,0.,9.*a1,0./)
tminv(4,:)=(/4.*a1,-a1,-2.*a1,0.,0.,-6.*a1,6.*a1,-9.*a1,0./)
tminv(5,:)=(/4.*a1,2.*a1,a1,6.*a1,3.*a1,6.*a1,3.*a1,0.,9.*a1/)
tminv(6,:)=(/4.*a1,2.*a1,a1,-6.*a1,-3.*a1,6.*a1,3.*a1,0.,-9.*a1/)
tminv(7,:)=(/4.*a1,2.*a1,a1,-6.*a1,-3.*a1,-6.*a1,-3.*a1,0.,9.*a1/)
tminv(8,:)=(/4.*a1,2.*a1,a1,6.*a1,3.*a1,-6.*a1,-3.*a1,0.,-9.*a1/)
do i=0,8
do j=0,8
    sumcc=0.0
    do l=0,8
        sumcc=sumcc+tminv(i,l)*tm(l,j)
    end do
ev(i,j)=sumcc!µ¥Î»¾ØÕó
end do
    end do
do i=0,8
    print*,(ev(i,j),j=0,8)
end do
pause
    
uo=0.05
rhoo=1.00
dx=1.0
dy=dx
dt=1.0
alpha=0.001
Re=uo*m/alpha
print*,"Re=",Re
omega=1.0/(3.*alpha+0.5)
tau=1./omega

sm(:)=(/1.0,1.4,1.4,1.0,1.2,1.0,1.2,tau,tau/)!¶Ô½Ç¾ØÕóS
do i=0,8
do j=0,8
    stmiv(i,j)=tminv(i,j)*sm(j)
end do
    end do
    do i=0,8
        print*,(stmiv(i,j),j=0,8)
    end do

mstep=100000
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
    
!main loop
do kk=1,mstep
    call collesion(u,v,f,feq,rho,omega,w,cx,cy,n,m,tm,tminv,stmiv)
    call streaming(f,n,m)
!---------------------------
    call sfbound(f,n,m,uo)
    call rhouv(f,rho,u,v,cx,cy,n,m)
    print*,u(0,m/2),v(0,m/2),rho(0,m/2),u(n,m/2),v(n,m/2),rho(n,m/2)
    write(8,*)kk,u(n/2,m/2),v(n/2,m/2)
    end do
call result(u,v,rho,uo,n,m)
stop
    end
!end of the main program

subroutine collesion(u,v,f,feq,rho,omega,w,cx,cy,n,m,tm,tminv,stmiv)
real f(0:8,0:n,0:m)
real feq(0:8,0:n,0:m),rho(0:n,0:m)
real w(0:8),cx(0:8),cy(0:8)
real u(0:n,0:m),v(0:n,0:m)
real tm(0:8,0:8),tminv(0:8,0:8),stmiv(0:8,0:8)
real fmom(0:8,0:n,0:m),fmeq(0:8,0:n,0:m)
!Calculate equilibrium moments
do i=0,n
    do j=0,m
        fmeq(0,i,j)=rho(i,j)
        fmeq(1,i,j)=rho(i,j)*(-2.0+3.0*rho(i,j)*(u(i,j)*u(i,j)+v(i,j)*v(i,j)))
        fmeq(2,i,j)=rho(i,j)*(1.0-3.0*rho(i,j)*(u(i,j)*u(i,j)+v(i,j)*v(i,j)))
        fmeq(3,i,j)=rho(i,j)*u(i,j)
        fmeq(4,i,j)=-rho(i,j)*u(i,j)
        fmeq(5,i,j)=rho(i,j)*v(i,j)
        fmeq(6,i,j)=-rho(i,j)*v(i,j)
        fmeq(7,i,j)=rho(i,j)*(u(i,j)*u(i,j)-v(i,j)*v(i,j))
        fmeq(8,i,j)=rho(i,j)*u(i,j)*v(i,j)
    end do
end do
!Calculate Moments
do i=0,n
    do j=0,m
        do k=0,8
            suma=0.0
            do l=0,8
                suma=suma+tm(k,l)*f(l,i,j)
            end do
            fmom(k,i,j)=suma
        end do
    end do
end do
!Collesion in the moment space
do i=0,n
    do j=0,m
        do k=0,8
            sumb=0.0
            do l=0,8
                sumb=sumb+stmiv(k,l)*(fmom(l,i,j)-fmeq(l,i,j))
            end do
            f(k,i,j)=f(k,i,j)-sumb
        end do
    end do
end do
return
end
    
subroutine streaming(f,n,m)
real f(0:8,0:n,0:m)
!streaming
do j=0,m
    do i=n,1,-1      !right to left
        f(1,i,j)=f(1,i-1,j)
    end do
    do i=0,n-1       !left to right
        f(3,i,j)=f(3,i+1,j)
    end do
end do
do j=m,1,-1
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
do j=0,m-1
    do i=0,n        !bottom to top
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
f(8,i,m)=f(6,i,m)+rhon*uo/6.0!f(1,i,m)-f(3,i,m)=2/3*rhon*uo
f(7,i,m)=f(5,i,m)-rhon*uo/6.0
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
rho(i,m)=f(0,i,m)+f(1,i,m)+f(3,i,m)+2.*(f(2,i,m)+f(6,i,m)+f(5,i,m))
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
    
subroutine result(u,v,rho,uo,n,m)
real u(0:n,0:m),v(0:n,0:m)
real rho(0:n,0:m),strf(0:n,0:m)
open(5,file="streamf.plt")
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
write(2,*)"VARIABLES=X,Y,U,V,S"
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
write(2,*)(strf(i,j),i=0,n)
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
write(3,*)j/float(m),u(n/2,j)/uo,u(n/4,j)/uo,u(3*n/4,j)/uo
end do
do i=0,n
write(4,*) i/float(n),v(i,m/2)/uo
end do
return
    end