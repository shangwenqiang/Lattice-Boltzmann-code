!  Console1.f90 
!
!  FUNCTIONS:
!  Console1 - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: Console1
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************
program Console1
parameter(n=44,m=24,r=24)
real f(0:18,0:n,0:m,0:r),ff(0:18,0:n,0:m,0:r)
real feq(0:18,0:n,0:m,0:r),rho(0:n,0:m,0:r)
real w(0:18),cx(0:18),cy(0:18),cz(0:18)
real u(0:2,0:n,0:m,0:r)
integer i,j,z
cx(:)=(/0,1,-1,0,0,0,0,1,-1,1,-1,1,-1,1,-1,0,0,0,0/)
cy(:)=(/0,0,0,1,-1,0,0,1,1,-1,-1,0,0,0,0,1,-1,1,-1/)
cz(:)=(/0,0,0,0,0,1,-1,0,0,0,0,1,1,-1,-1,1,1,-1,-1/)
w(:)=(/12./36.,2./36.,2./36.,2./36.,2./36.,2./36.,2./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36./)
uo=0.1
sumvelo=0.0
rhoo=1.0
dx=1.0
dy=dx
dz=dx
dt=1.0
alpha=0.16667
Re=uo*n/alpha
print*,"Re=",Re
omega=1.0/(3.*alpha+0.5)
mstep=10000
do z=0,r
do j=0,m
do i=0,n
rho(i,j,z)=rhoo
u(0,i,j,z)=0.0
u(1,i,j,z)=0.0
u(2,i,j,z)=0.0
end do
end do
end do
do z=1,r-1
do j=1,m-1
u(0,1,j,z)=uo
end do
end do
DO i=0,n
DO j=0,m
do z=0,r
t1=u(0,i,j,z)*u(0,i,j,z)+u(1,i,j,z)*u(1,i,j,z)+u(2,i,j,z)*u(2,i,j,z)
DO k=0,18
t2=u(0,i,j,z)*cx(k)+u(1,i,j,z)*cy(k)+u(2,i,j,z)*cz(k)
feq(k,i,j,z)=rho(i,j,z)*w(k)*(1.0+3.0*t2+4.50*t2*t2-1.50*t1)
f(k,i,j,z)=feq(k,i,j,z)
end do
END DO
END DO
END DO 
!main loop
do kk=1,mstep
call collesion(u,f,ff,feq,rho,omega,w,cx,cy,cz,n,m,r)
call streaming(f,ff,n,m,r,cx,cy,cz)
call sfbound(f,n,m,r,uo)
call rhouv(f,rho,u,cx,cy,cz,n,m,r)
print*,kk,rho(n,m/2,r/2),u(0,n/2,m/2,r/2)
end do
!end of the main loop
call result(u,rho,n,m,r)
stop
    end

subroutine collesion(u,f,ff,feq,rho,omega,w,cx,cy,cz,n,m,r)   
real f(0:18,0:n,0:m,0:r),ff(0:18,0:n,0:m,0:r)
real feq(0:18,0:n,0:m,0:r),rho(0:n,0:m,0:r)
real w(0:18),cx(0:18),cy(0:18),cz(0:18)
real u(0:2,0:n,0:m,0:r)
real omega
do z=1,r-1
do j=1,m-1
do i=1,n-1
t1=u(0,i,j,z)*u(0,i,j,z)+u(1,i,j,z)*u(1,i,j,z)+u(2,i,j,z)*u(2,i,j,z)
do k=0,18
t2=u(0,i,j,z)*cx(k)+u(1,i,j,z)*cy(k)+u(2,i,j,z)*cz(k)
feq(k,i,j,z)=rho(i,j,z)*w(k)*(1.0+3.0*t2+4.50*t2*t2-1.50*t1)
ff(k,i,j,z)=omega*feq(k,i,j,z)+(1.-omega)*f(k,i,j,z)
end do
end do
end do
end do
return 
end
    

subroutine streaming(f,ff,n,m,r,cx,cy,cz)
real f(0:18,0:n,0:m,0:r),ff(0:18,0:n,0:m,0:r)
real w(0:18),cx(0:18),cy(0:18),cz(0:18)
integer ip,jp,zp
do z=1,r-1
do j=1,m-1
do i=1,n-1 
do k=0,18
ip=i-int(cx(k))
jp=j-int(cy(k))
zp=z-int(cz(k))
f(k,i,j,z)=ff(k,ip,jp,zp)
end do
end do
end do
end do
return 
end

subroutine sfbound(f,n,m,r,uo)
real f(0:18,0:n,0:m,0:r)
do z=1,r-1
do j=1,m-1
!left boundary
rhow=(f(0,1,j,z)+f(3,1,j,z)+f(4,1,j,z)+f(5,1,j,z)+f(6,1,j,z)+f(15,1,j,z)+f(16,1,j,z)+f(17,1,j,z)+f(18,1,j,z)+2*(f(2,1,j,z)+f(8,1,j,z)+f(10,1,j,z)+f(12,1,j,z)+f(14,1,j,z)))/(1-uo)
f(1,1,j,z)=f(2,1,j,z)+1./3.*rhow*uo
f(7,1,j,z)=f(10,1,j,z)+1./6.*rhow*uo+1./2.*(f(4,1,j,z)+f(16,1,j,z)+f(18,1,j,z)-f(3,1,j,z)-f(15,1,j,z)-f(17,1,j,z))
f(9,1,j,z)=f(8,1,j,z)+1./6.*rhow*uo-1./2.*(f(4,1,j,z)+f(16,1,j,z)+f(18,1,j,z)-f(3,1,j,z)-f(15,1,j,z)-f(17,1,j,z))
f(11,1,j,z)=f(14,1,j,z)+1./6.*rhow*uo+1./2.*(f(6,1,j,z)+f(17,1,j,z)+f(18,1,j,z)-f(5,1,j,z)-f(15,1,j,z)-f(16,1,j,z))
f(13,1,j,z)=f(12,1,j,z)+1./6.*rhow*uo-1./2.*(f(6,1,j,z)+f(17,1,j,z)+f(18,1,j,z)-f(5,1,j,z)-f(15,1,j,z)-f(16,1,j,z))
!right boundary
f(2,n-1,j,z)=2*f(2,n-2,j,z)-f(2,n-3,j,z)
f(8,n-1,j,z)=2*f(8,n-2,j,z)-f(8,n-3,j,z)
f(10,n-1,j,z)=2*f(10,n-2,j,z)-f(10,n-3,j,z)
f(12,n-1,j,z)=2*f(12,n-2,j,z)-f(12,n-3,j,z)
f(14,n-1,j,z)=2*f(14,n-2,j,z)-f(14,n-3,j,z)
end do
end do
!front boundry
do z=1,r-1
do i=1,n-1
f(3,i,1,z)=f(4,i,1,z)
f(7,i,1,z)=f(10,i,1,z)
f(8,i,1,z)=f(9,i,1,z)
f(15,i,1,z)=f(18,i,1,z)
f(17,i,1,z)=f(16,i,1,z)
!behind bounrdry
f(4,i,m-1,z)=f(3,i,m-1,z)
f(9,i,m-1,z)=f(8,i,m-1,z)
f(10,i,m-1,z)=f(7,i,m-1,z)
f(16,i,m-1,z)=f(17,i,m-1,z)
f(18,i,m-1,z)=f(15,i,m-1,z)
end do
end do
do j=1,m-1
do i=1,n-1
!down boundary
f(5,i,j,1)=f(6,i,j,1)
f(11,i,j,1)=f(14,i,j,1)
f(12,i,j,1)=f(13,i,j,1)
f(15,i,j,1)=f(18,i,j,1)
f(16,i,j,1)=f(17,i,j,1)
!up boundary
f(6,i,j,r-1)=f(5,i,j,r-1)
f(13,i,j,r-1)=f(12,i,j,r-1)
f(14,i,j,r-1)=f(11,i,j,r-1)
f(17,i,j,r-1)=f(16,i,j,r-1)
f(18,i,j,r-1)=f(15,i,j,r-1)
end do
end do
return
  end


subroutine rhouv(f,rho,u,cx,cy,cz,n,m,r)
real w(0:18),cx(0:18),cy(0:18),cz(0:18)
real u(0:2,0:n,0:m,0:r)
real f(0:18,0:n,0:m,0:r),rho(0:n,0:m,0:r)
!do i=i,n
!rho(i,m)=f(0,i,m)+f(1,i,m)+f(3,i,m)+2*(f(2,i,m)+f(6,i,m)+f(5,i,m))
!end do
do i=2,n-1 
do j=1,m-1
do z=1,r-1 
usum=0.0
vsum=0.0
wsum=0.0
ssum=0.0
DO k=0,18
ssum=ssum+f(k,i,j,z)
usum=usum+f(k,i,j,z)*cx(k)
vsum=vsum+f(k,i,j,z)*cy(k)
wsum=wsum+f(k,i,j,z)*cz(k)
END DO
rho(i,j,z)=ssum
u(0,i,j,z)=usum/rho(i,j,z)
u(1,i,j,z)=vsum/rho(i,j,z)
u(2,i,j,z)=wsum/rho(i,j,z)
!write(*,*) i,j,z,u(0,i,j,z)
END DO
END DO
end do
return
    end
    
    
subroutine  result(u,rho,n,m,r)
real u(0:2,0:n,0:m,0:r)
real rho(0:n,0:m,0:r)
integer i,j,z
open(2,file='out1/uvfield.dat')
write(2,*)"VARIABLES =X, Y, Z, U, V, W, RHO"
write(2,*)'ZONE',',I=',n-1,',J=',m-1,',K=',r-1,',T=POINTS'
do z=1,r-1
do j=1,m-1
do i=1,n-1
write(2,*) i,j,z,u(0,i,j,z),u(1,i,j,z),u(2,i,j,z),rho(i,j,z)    
end do
end do
end do
return
end
    



