!LBM Code for 1-D, diffusion problems,D1Q2
parameter (m=100) !m is the number of lattice nodes
real fo(0:m),f1(0:m),f2(0:m),rho(0:m),feq(0:m),x(0:m)
integer i
open(2,file='result1.plt')
dt=1.0
dx=1.0
x(0)=0.0
flux=100.0
tk=20.0
do i=1,m
x(i)=x(i-1)+dx
end do
csq=dx*dx/(dt*dt)
alpha=0.25
omega=1.0/(alpha/(dt*csq)+0.5)!红皮书，公式3.23推导而来，omega为松弛频率
mstep=200 ! The total number of time steps
twall=1.0 ! Left hand wall temperature
! Initial condition
do i=0,m
rho(i)=0.0 ! Initial value of the domain temperature
f1(i)=0.5*rho(i)
f2(i)=0.5*rho(i)
    end do
    
do kk=1,mstep
! main loop
! collision process:
do i=0,m
rho(i)=f1(i)+f2(i)
!source=qq/(rho(i)*cp)！内热源
feq(i)=0.5*rho(i)
! since k1=k2=0.5, then feq1=feq2=feq
!f1(i)=(1.-omega)*f1(i)+omega*feq(i)+dt*0.5*source
!f2(i)=(1.-omega)*f2(i)+omega*feq(i)+dt*0.5*source!含源项
f1(i)=(1.-omega)*f1(i)+omega*feq(i)!红皮书，公式3.18
f2(i)=(1.-omega)*f2(i)+omega*feq(i)
end do

!Streaming process:
do i=1,m-1
f1(m-i)=f1(m-i-1) ! f1 streaming
f2(i-1)=f2(i) ! f2 streaming
end do

! Boundary condition
f1(0)=twall-f2(0)! constant temperature boundary condition, x=0
!f1(0)=f1(1)+f2(1)-f2(0)+flux*dx/tk; !constant heat flux boundary condition,x=0
f1(m)=f1(m-1) ! adiabatic boundary condition, x=L
f2(m)=f2(m-1)! adiabatic boundary condition, x=L
    end do
! end of the main loop

do i=0,m
write(2,*)x(i), rho(i)
end do
stop
end