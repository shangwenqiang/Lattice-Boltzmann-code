! LBM for 1-D advection-diffusion
parameter (m=100) !m is the number of lattice nodes
real dens,fo(0:m),f1(0:m),f2(0:m),rho(0:m),feq1,x(0:m),feq2
integer i
open(2,file='result.plt')
u=0.1
dt=1.0
dx=1.0
x(0)=0.0
do i=1,m
x(i)=x(i-1)+dx
end do
ck=dx/dt
csq=ck*ck
alpha=0.25
omega=1.0/(alpha/(dt*csq)+0.5)
mstep=400 ! The total number of time steps
twall=1.0 ! Left hand wall temperature
! Initial condition
do i=0,m
rho(i)=0.0 ! Initial value of the domain temperature
f1(i)=0.5*rho(i)
f2(i)=0.5*rho(i)
end do
do kk=1,mstep ! main loop
! collision process:
do i=0,m
rho(i)=f1(i)+f2(i)
feq1=0.5*rho(i)*(1.0+u/ck) ! extra term added to simulate advection
feq2=0.5*rho(i)*(1.0-u/ck) ! w1=w2=0.5
f1(i)=(1.-omega)*f1(i)+omega*feq1
f2(i)=(1.-omega)*f2(i)+omega*feq2
end do
! Streaming process:
do i=1,m-1
f1(m-i)=f1(m-i-1) ! f1 streaming
f2(i-1)=f2(i) ! f2 streaming
end do
! Boundary condition
f1(0)=twall-f2(0) ! constant temperature boundary condition, x=0
f1(m)=f1(m-1) ! adiabatic boundary condition, x=L
f2(m)=f2(m-1) ! adiabatic boundary condition, x=L
end do ! end of the main loop
do i=0,m
write(2,*)x(i), rho(i)
end do
stop
end