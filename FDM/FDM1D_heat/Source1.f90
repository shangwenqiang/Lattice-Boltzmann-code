!The Finite Difference Code:FDM,Finite Difference Code for 1-D diffusion problems
parameter (m=100)
real dens,fo(0:m),f(0:m)
integer i
open(2,file='finitrs.plt')
dx=1.0
dt=0.500
alpha=0.25
mstep=400
do i=0,m
fo(i)=0.0
end do
fo(0)=1.0 !initial condition for old value of f at x=0.
f(0)=1.0 ! initial condition for updated value of f at x=0.
fo(m)=fo(m-1) ! initial condition for old value of f at x=L
f(m)=f(m-1) ! initial condition for updated value of f at x=L
do kk=1,mstep
! main loop
do i=1,m-1
f(i)=fo(i)+dt*alpha*(fo(i+1)-2.*fo(i)+fo(i-1))/(dx*dx)
end do
do i=1,m-1
fo(i)=f(i) ! updating
end do
fo(m)=f(m-1) ! updating the boundary condition at x=L
end do
! end of the main loop
x=0.0
do i=0,m
write(2,*)x,f(i)
x=x+dx
end do
stop
end