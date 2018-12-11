! Finite Difference code for 1-D advection-diffusion problem.
parameter (n=200)
real dens,fo(0:n),f(0:n)
integer i
open(2,file='finitrs.plt')
!
dx=0.500
dt=0.2500
u=0.10
alpha=0.25
mstep=1600
do i=0,n
fo(i)=0.0
end do
fo(0)=1.0
f(0)=1.0
fo(n)=0.0
f(n)=0.0
do kk=1,mstep
do i=1,n-1
adv=dt*u*(fo(i)-fo(i-1))/dx
f(i)=fo(i)+dt*alpha*(fo(i+1)-2.*fo(i)+fo(i-1))/(dx*dx)-adv
end do
!
do i=1,n-1
fo(i)=f(i)
end do
end do
x=0.0
do i=0,n
write(2,*)x,f(i)
x=x+dx
end do
stop
end