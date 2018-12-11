! The LBM code for 3D
parameter (n=100,m=100,s=100)
real f(0:14,0:n,0:m,0:s),feq
real rho(0:n,0:m,0:s),x(0:n),y(0:m),z(0:s)
real w(0:14)
integer i,j,k,h
open(2,file='Qresu.plt')
!
dx=1.0
dy=dx
dz=dx
x(0)=0.0
y(0)=0.0
z(0)=0.0
do i=1,n
x(i)=x(i-1)+dx
end do
do j=1,m
y(j)=y(j-1)+dy
    end do
do h=1,s
z(h)=z(h-1)+dz
    end do
    
dt=1.0
tw=1.0
alpha=0.25
csq=(dx*dx)/(dt*dt)
omega=1.0/(3.*alpha/(csq*dt)+0.5)
mstep=400

!权重因子初始化
w(0)=16./72.
do i=1,6
w(i)=8./72.
end do
do i=7,14
w(i)=1./72.
    end do
    
do h=0,s 
do j=0,m
do i=0,n
rho(i,j,h)=0.0 ! initial field
end do
end do
    end do
    
do h=0,s
do j=0,m
do i=0,n
do k=0,14
f(k,i,j,h)=w(k)*rho(i,j,h)
if(j.eq.0) f(k,i,j,h)=w(k)*tw
end do
end do
end do
end do
    
!开始循环
do kk=1,mstep
do h=0,s
do j=0,m
do i=0,n
sum=0.0
do k=0,14
sum=sum+f(k,i,j,h)
end do
rho(i,j,h)=sum
end do
end do
end do
do h=0,s
do j=0,m
do i=0,n
do k=0,14
feq=w(k)*rho(i,j,h)
f(k,i,j,h)=omega*feq+(1.-omega)*f(k,i,j,h)
end do
end do
end do
end do

! streaming
do h=0,s
do j=m,1,-1
do i=0,n
f(1,i,j,h)=f(1,i,j-1,h)
end do
end do
end do

do h=s,1,-1
do j=0,m
do i=0,n
f(2,i,j,h)=f(2,i,j,h-1)
end do
end do
end do

do h=0,s
do j=0,m-1
do i=0,n
f(3,i,j,h)=f(3,i,j+1,h)
end do
end do
end do

do h=0,s-1
do j=0,m
do i=0,n
f(4,i,j,h)=f(4,i,j,h+1)
end do
end do
end do

do h=0,s
do j=0,m
do i=n,1,-1
f(5,i,j,h)=f(5,i-1,j,h)
end do
end do
end do

do h=0,s
do j=0,m
do i=0,n-1
f(6,i,j,h)=f(6,i+1,j,h)
end do
end do
end do

do h=s,1,-1
do j=m,1,-1
do i=n,1,-1
f(7,i,j,h)=f(7,i-1,j-1,h-1)
end do
end do
end do

do h=s,1,-1
do j=m,1,-1
do i=0,n-1
f(8,i,j,h)=f(8,i+1,j-1,h-1)
end do
end do
end do

do h=0,s-1
do j=m,1,-1
do i=0,n-1
f(9,i,j,h)=f(9,i+1,j-1,h+1)
end do
end do
end do

do h=0,s-1
do j=m,1,-1
do i=n,1,-1
f(10,i,j,h)=f(10,i-1,j-1,h+1)
end do
end do
end do

do h=s,1,-1
do j=0,m-1
do i=0,n-1
f(11,i,j,h)=f(11,i+1,j+1,h-1)
end do
end do
end do

do h=s,1,-1
do j=0,m-1
do i=n,1,-1
f(12,i,j,h)=f(12,i-1,j+1,h-1)
end do
end do
end do

do h=0,s-1
do j=0,m-1
do i=n,1,-1
f(13,i,j,h)=f(13,i-1,j+1,h+1)
end do
end do
end do

do h=0,s-1
do j=0,m-1
do i=0,n-1
f(14,i,j,h)=f(14,i+1,j+1,h+1)
end do
end do
end do

! Boundary conditions
!左右面
do h=0,s
do i=0,n
f(1,i,0,h)=w(1)*tw+w(3)*tw-f(3,i,0,h)
f(7,i,0,h)=w(7)*tw+w(14)*tw-f(14,i,0,h)
f(8,i,0,h)=w(8)*tw+w(13)*tw-f(13,i,0,h)
f(9,i,0,h)=w(9)*tw+w(12)*tw-f(12,i,0,h)
f(10,i,0,h)=w(10)*tw+w(11)*tw-f(11,i,0,h)
f(3,i,m,h)=-f(1,i,m,h)
f(14,i,m,h)=-f(7,i,m,h)
f(13,i,m,h)=-f(8,i,m,h)
f(12,i,m,h)=-f(9,i,m,h)
f(11,i,m,h)=-f(10,i,m,h)
end do
end do

!前面
do h=0,s
do j=0,m
f(6,n,j,h)=-f(5,n,j,h)
f(14,n,j,h)=-f(7,n,j,h)
f(11,n,j,h)=-f(10,n,j,h)
f(9,n,j,h)=-f(12,n,j,h)
f(8,n,j,h)=-f(13,n,j,h)

f(6,0,j,h)=-f(5,0,j,h)
f(14,0,j,h)=-f(7,0,j,h)
f(11,0,j,h)=-f(10,0,j,h)
f(9,0,j,h)=-f(12,0,j,h)
f(8,0,j,h)=-f(13,0,j,h)
end do
end do

!上面
do j=0,m
do i=0,n
f(4,i,j,s)=-f(2,i,j,s)
f(9,i,j,s)=-f(12,i,j,s)
f(10,i,j,s)=-f(11,i,j,s)
f(13,i,j,s)=-f(8,i,j,s)
f(14,i,j,s)=-f(7,i,j,s)

f(4,i,j,0)=-f(2,i,j,0)
f(9,i,j,0)=-f(12,i,j,0)
f(10,i,j,0)=-f(11,i,j,0)
f(13,i,j,0)=-f(8,i,j,0)
f(14,i,j,0)=-f(7,i,j,0)
end do
end do

    end do  
!循环结束

do h=0,s
do j=0,m
do i=0,n
sum=0.0
do k=0,14
sum=sum+f(k,i,j,h)
end do
rho(i,j,h)=sum
end do
end do
end do
    
print *, rho(n/2,0,s/2)
write(2,*)"VARIABLES =X,Y,Z,T"
write(2,*)"ZONE ","I=",n+1,",","J=",m+1,",","K=",s+1,",","F=BLOCK"
do h=0,s
    do j=0,m
        write(2,*)(x(i),i=0,n)
    end do
    end do
do h=0,s
    do j=0,m
        write(2,*)(y(j),i=0,n)
    end do
    end do
do h=0,s
    do j=0,m
        write(2,*)(z(h),i=0,n)
    end do 
    end do

do h=0,s
    do j=0,m
        write(2,*)(rho(i,j,h),i=0,n)
    end do
    end do
stop
end