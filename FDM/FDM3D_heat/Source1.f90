! Finite Difference Code for 3-D
parameter (n=100,m=100,l=100)
real fo(0:n,0:m,0:l),f(0:n,0:m,0:l),x(0:n),y(0:m),z(0:l)
integer i,j,k
open(2,file='fin2drs.plt')
!
dx=1.0
dy=dx
dz=dx
time=0.0
dt=0.20     
alpha=0.25  !扩散系数
mstep=2000  

!全场温度初始化为0
do k=0,l
do j=0,m
do i=0,n
fo(i,j,k)=0.0
end do
end do
end do
    
!设置温度边界条件
do k=0,l
    do i=0,n
fo(i,0,k)=1.0
f(i,0,k)=1.0
fo(i,m,k)=0.0
f(i,m,k)=0.0
! adiabatic boundary condition
fo(i,m-1,k)=fo(i,m,k)
f(i,m-1,k)=f(i,m,k)
    end do
end do
    
do k=0,l
    do j=0,m
fo(0,j,k)=0.0
f(0,j,k)=0.0
fo(n,j,k)=0.0
f(n,j,k)=0.0
! adiabatic boundary condition
!fo(1,j,k)=f(0,j,k)
!f(1,j,k)=f(0,j,k)
!fo(n-1,j,k)=f(n,j,k)
!f(n-1,j,k)=f(n,j,k)
    end do
    end do
    
do j=0,m
    do i=0,n
fo(i,j,0)=0.0
f(i,j,0)=0.0
fo(i,j,l)=0.0
f(i,j,l)=0.0
! adiabatic boundary condition
!fo(i,j,1)=f(i,j,0)
!f(i,j,1)=f(i,j,0)
!fo(i,j,l-1)=f(i,j,l)
!f(i,j,l-1)=f(i,j,l)
    end do
    end do

!开始循环    
do kk=1,mstep
do k=1,l-1
do j=1,m-1
do i=1,n-1
termx=(fo(i+1,j,k)+fo(i-1,j,k))/(dx*dx)
termy=(fo(i,j+1,k)+fo(i,j-1,k))/(dy*dy)
termz=(fo(i,j,k+1)+fo(i,j,k-1))/(dz*dz)
dd=1./(dx*dx)+1./(dy*dy)+1./(dz*dz)
f(i,j,k)=fo(i,j,k)+dt*alpha*(termx+termy+termz-2.0*fo(i,j,k)*dd) !公式推导
end do
end do
end do
!
do k=1,l-1
do j=1,m-1
do i=1,n-1
fo(i,j,k)=f(i,j,k)
end do
end do
end do

time=time+dt
end do
x(0)=0.0
do i=1,n
x(i)=x(i-1)+dx
end do
y(0)=0.0
do j=1,m
y(j)=y(j-1)+dy
    end do
z(0)=0.0
do k=1,l
z(k)=z(k-1)+dz
end do
write(2,*)"VARIABLES =X,Y,Z,T"
write(2,*)"ZONE ","I=",n+1,",","J=",m+1,",","K=",l+1,",","F=BLOCK"
do k=0,l
    do j=0,m
        write(2,*)(x(i),i=0,n)
    end do
    end do
do k=0,l
    do j=0,m
        write(2,*)(y(j),i=0,n)
    end do
    end do
do k=0,l
    do j=0,m
        write(2,*)(z(k),i=0,n)
    end do 
    end do

do k=0,l
    do j=0,m
        write(2,*)(f(i,j,k),i=0,n)
    end do
    end do
stop
end