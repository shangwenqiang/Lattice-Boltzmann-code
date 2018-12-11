parameter (N=405,M=405)
real f(0:8,0:n,0:m),ff(0:8,0:n,0:m),feq(0:8,0:n,0:m)
real rho(0:n,0:m)
real w(0:8), cx(0:8),cy(0:8)
real uv(0:1,0:n,0:m),uuvv(0:1,0:n,0:m),velocity(0:n,0:m)
real g(0:8,0:n,0:m), gg(0:8,0:n,0:m),geq(0:8,0:n,0:m),th(0:n,0:m),strf(0:n,0:m)
integer s(0:n,0:m),bb(0:n,0:m)
integer i,j,step
real omega,alpha,visco,error1

cx(:)=(/0.0,1.0,0.0,-1.0,0.0,1.0,-1.0,-1.0,1.0/)
cy(:)=(/0.0,0.0,1.0,0.0,-1.0,1.0,1.0,-1.0,-1.0/)
w(:)=(/4./9.,1./9.,1./9.,1./9.,1./9.,1./36.,1./36.,1./36.,1./36./)
mstep=80000000
step=5000
uo=0.0005
sumvelo=0.0
rhoo=1.0
drho=0.1
rhoo1=rhoo+drho
rhoo2=rhoo
dx=1.0
dy=dx
dt=1.0
tw=1.0
th=0.0
g=0.0

visco=0.1666666666666667
pr=3.8286
alpha=visco/pr
!!Re=uo*m/alpha
!print *, "Re=", Re
!omega=1.0/(3*0.0004+0.5)
omega=1.0/(3.*visco+0.5)
omegat=1.0/(3.*1.*alpha+0.5)
omegat2=1.0/(3.*10.*alpha+0.5)


do j=0,m
do i=0,n
rho(i,j)=rhoo
uv(0,i,j)=0.0
uv(1,i,j)=0.0
end do
    end do
     
do j=1,m-1
rho(1,j)=rhoo1
rho(n-1,j)=rhoo2
end do 
 kk=0
 call outresult(kk,n,m,uv,rho,velocity,s,error1)   
    
do j=0,m
do i=0,n
    !============1阶
s(i,j)=0!fluid
if((j .eq. 0) .or. (j .eq. m) .or.(j .eq. 1) .or. (j .eq. m-1)) then
s(i,j)=1!solid
end if
bb(i,j)=0
!1阶固体
if (( (j .ge.0 .and. j .le. 135 ).and. (i .ge. 270  .and. i .le. n) ).or. ((j .ge.270.and. j .le. m).and.(i .ge. 0 .and. i .le. 135) )   ) then
s(i,j)=1!solid
end if
!==============================================================================
!2阶固体
if (          (j .ge.0 .and. j .le. 45) .and.(       (i .ge. 90 .and. i .le. 135).or.(i .ge. 225 .and. i .le. 270))         )then
s(i,j)=1!solid
end if
if (          (j .ge.90 .and. j .le. 135) .and.(       (i .ge. 0 .and. i .le. 45).or.(i .ge.135 .and. i .le. 180)          )         )then
s(i,j)=1!solid
end if
if (          (j .ge.135 .and. j .le. 180) .and.(       (i .ge. 90 .and. i .le. 135).or.(i .ge.225 .and. i .le. 270).or.(i .ge.360 .and. i .le. 405)          )         )then
s(i,j)=1!solid
end if
if (          (j .ge.225 .and. j .le. 270) .and.(       (i .ge. 0 .and. i .le. 45).or.(i .ge.135 .and. i .le. 180).or.(i .ge.270 .and. i .le. 315)          )         )then
s(i,j)=1!solid
end if
if (          (j .ge.270 .and. j .le. 315) .and.(       (i .ge. 225 .and. i .le. 270).or.(i .ge.360 .and. i .le. 405)          )         )then
s(i,j)=1!solid
end if
if (          (j .ge.360 .and. j .le. 405) .and.(       (i .ge. 270 .and. i .le. 315).or.(i .ge.135 .and. i .le. 180)          )         )then
s(i,j)=1!solid
end if
!==============================================================================
end do
 end do
   !区分界面形式1东,2西,3南,4北，边角
   do i=1,n-1
    do j=1,m-1
        !东
        if(  (s(i,j).eq.1)&
                .and.(s(i,j+1).eq.1).and.(s(i,j-1).eq.1).and.  (s(i-1,j).eq.1).and.(s(i+1,j).eq.1)&
                .and.(s(i+1,j+1).ne.1).and.(s(i-1,j-1).ne.1).and.  (s(i-1,j+1).eq.1).and.(s(i+1,j-1).eq.1) )then
               s(i+1,j+1)=1
               s(i-1,j-1)=1
                endif
      if( ( (s(i,j).eq.1).and.(s(i,j-1).eq.1).and.(s(i,j+1).eq.1).and.  (s(i+1,j).ne.1) ))then
          bb(i,j)=1!dong
      endif
      !xi西
         if(  (s(i,j).eq.1).and.(s(i,j-1).eq.1).and.(s(i,j+1).eq.1).and.  (s(i-1,j).ne.1)         )then
          bb(i,j)=2!xi
         endif
         !xi南
         if(  (s(i,j).eq.1).and.(s(i+1,j).eq.1).and.(s(i-1,j).eq.1).and.  (s(i,j-1).ne.1)         )then
          bb(i,j)=3!南
         endif
         if(  (s(i,j).eq.1).and.(s(i+1,j).eq.1).and.(s(i-1,j).eq.1).and.  (s(i,j+1).ne.1)         )then
          bb(i,j)=4!北
         endif
         !!左上
         if(  (s(i,j).eq.1).and.(s(i+1,j).eq.1).and.(s(i,j-1).eq.1).and.  (s(i+1,j-1).eq.1)&
            .and.(s(i-1,j+1).ne.1).and.(s(i,j+1).ne.1).and.(s(i-1,j).ne.1) )then
          bb(i,j)=5
            endif
            !左下
         if(  (s(i,j).eq.1).and.(s(i,j+1).eq.1).and.(s(i+1,j).eq.1).and.  (s(i+1,j+1).eq.1)&
            .and.(s(i-1,j).ne.1).and.(s(i,j-1).ne.1).and.(s(i-1,j-1).ne.1) )then
          bb(i,j)=7
            endif
            !右上
            if(  (s(i,j).eq.1).and.(s(i,j-1).eq.1).and.(s(i-1,j).eq.1).and.  (s(i-1,j-1).eq.1)&
            .and.(s(i+1,j).ne.1).and.(s(i,j+1).ne.1).and.(s(i+1,j+1).ne.1) )then
          bb(i,j)=6
            endif
            !右下
            if(  (s(i,j).eq.1).and.(s(i,j+1).eq.1).and.(s(i-1,j).eq.1).and.  (s(i-1,j+1).eq.1)&
            .and.(s(i+1,j).ne.1).and.(s(i,j-1).ne.1).and.(s(i+1,j-1).ne.1) )then
          bb(i,j)=8
            endif
            !内角点
            !左下
            if(  (s(i,j).eq.1)&
                .and.(s(i,j+1).eq.1).and.(s(i,j-1).eq.1).and.  (s(i-1,j).eq.1).and.(s(i+1,j).eq.1)&
                .and.(s(i+1,j+1).eq.1).and.(s(i-1,j-1).eq.1).and.  (s(i-1,j+1).eq.1).and.(s(i+1,j-1).ne.1) )then
          bb(i,j)=10
                endif
              
                  !右下
            if(  (s(i,j).eq.1)&
                .and.(s(i,j+1).eq.1).and.(s(i,j-1).eq.1).and.  (s(i-1,j).eq.1).and.(s(i+1,j).eq.1)&
                .and.(s(i+1,j+1).eq.1).and.(s(i-1,j-1).ne.1).and.  (s(i-1,j+1).eq.1).and.(s(i+1,j-1).eq.1) )then
          bb(i,j)=9
                endif
                   !youshang
            if(  (s(i,j).eq.1)&
                .and.(s(i,j+1).eq.1).and.(s(i,j-1).eq.1).and.  (s(i-1,j).eq.1).and.(s(i+1,j).eq.1)&
                .and.(s(i+1,j+1).eq.1).and.(s(i-1,j-1).eq.1).and.  (s(i-1,j+1).ne.1).and.(s(i+1,j-1).eq.1) )then
          bb(i,j)=12
                endif
                      !左上
            if(  (s(i,j).eq.1)&
                .and.(s(i,j+1).eq.1).and.(s(i,j-1).eq.1).and.  (s(i-1,j).eq.1).and.(s(i+1,j).eq.1)&
                .and.(s(i+1,j+1).ne.1).and.(s(i-1,j-1).eq.1).and.  (s(i-1,j+1).eq.1).and.(s(i+1,j-1).eq.1) )then
          bb(i,j)=11
                endif     
    end do
end do        
!do j=0,m-1
!u(0,j)=uo
!v(0,j)=0.0
!end do
    
DO i=0,n
DO j=0,m
t1=uv(0,i,j)*uv(0,i,j)+uv(1,i,j)*uv(1,i,j)
DO k=0,8
t2=uv(0,i,j)*cx(k)+uv(1,i,j)*cy(k)
feq(k,i,j)=rho(i,j)*w(k)*(1.0+3.0*t2+4.50*t2*t2-1.50*t1)
f(k,i,j)=feq(k,i,j)
END DO
END DO
END DO 
PRINT*,"粘度、热扩散率、LX、LY=",visco,alpha,404,404

!===========  
    do kk=1,mstep
   do i=1,n-1
       do j=1,m-1
          uuvv(0,i,j)=uv(0,i,j)
           uuvv(1,i,j)=uv(1,i,j)
       end do
   end do
   !!!======输出
!==========evo
  call collesion(uv,f,ff,feq,rho,omega,w,cx,cy,n,m,s)
  call streaming(f,ff,n,m,rhoo1,cx,cy)
  call sfbound(f,ff,n,m,rhoo1,rhoo2,s)
  call getrhouv(f,rho,uv,cx,cy,n,m,velocity,alpha,visco,s)
  call  Error(uv,uuvv,n,m,kk,error1)
  if (mod(kk,step).eq.0) then
    call outresult(kk,n,m,uv,rho,velocity,s,error1)
  end if
  if (error1.le.0.000001)  go to 1

  !=====主循环 
    end do 
    !====后输出
1 call outresult(kk,n,m,uv,rho,velocity,s,error1)   
 end    
subroutine collesion(uv,f,ff,feq,rho,omega,w,cx,cy,n,m,s)
real f(0:8,0:n,0:m),ff(0:8,0:n,0:m)
real feq(0:8,0:n,0:m),rho(0:n,0:m)
real w(0:8),cx(0:8),cy(0:8)
real omega
integer s(0:n,0:m)
real uv(0:1,0:n,0:m)
DO i=1,n-1
DO j=1,m-1
    if (s(i,j).eq.0) then
t1=uv(0,i,j)*uv(0,i,j)+uv(1,i,j)*uv(1,i,j)
DO k=0,8
t2=uv(0,i,j)*cx(k)+uv(1,i,j)*cy(k)
feq(k,i,j)=rho(i,j)*w(k)*(1.0+3.0*t2+4.50*t2*t2-1.50*t1)
ff(k,i,j)=omega*feq(k,i,j)+(1.-omega)*f(k,i,j)
END DO
endif 
END DO
END DO
return
    end
 subroutine streaming(f,ff,n,m,rhoo1,cx,cy)
real f(0:8,0:n,0:m),ff(0:8,0:n,0:m)
real cx(0:8),cy(0:8),w(0:8)
!cx(:)=(/0.0,1.0,0.0,-1.0,0.0,1.0,-1.0,-1.0,1.0/)
!cy(:)=(/0.0,0.0,1.0,0.0,-1.0,1.0,1.0,-1.0,-1.0/)
integer ip,jp
DO J=1,m-1
DO I=1,n-1
 
DO K=0,8
 ip=I-int(CX(K))
 jp=J-int(CY(K))
f(K,I,J)=ff(K,ip,jp)
ENDDO

ENDDO 
ENDDO   
return
end   
 subroutine sfbound(f,ff,n,m,rhoo1,rhoo2,s)
real f(0:8,0:n,0:m),ff(0:8,0:n,0:m)
integer s(0:n,0:m)
do j=1,m-1
uw=1-(f(0,1,j)+f(2,1,j)+f(4,1,j)+2.*(f(3,1,j)+f(6,1,j)+f(7,1,j)))/rhoo1
f(1,1,j)=f(3,1,j)+2.*rhoo1*uw/3.
f(5,1,j)=f(7,1,j)+rhoo1*uw/6.-0.5*(f(2,1,j)-f(4,1,j))
f(8,1,j)=f(6,1,j)+rhoo1*uw/6.+0.5*(f(2,1,j)-f(4,1,j))
end do
!!(1,1)下角点
!ff(1,1,1)=f(3,1,1)
!ff(2,1,1)=f(4,1,1)
!ff(5,1,1)=f(7,1,1)
!!(1,404)下角点
!ff(1,1,1)=f(3,1,1)
!ff(4,1,1)=f(2,1,1)
!ff(8,1,1)=f(6,1,1)
! bounce back on south boundary
do j=1,m-1
ue=-1+(f(0,n-1,j)+f(2,n-1,j)+f(4,n-1,j)+2.*(f(1,n-1,j)+f(5,n-1,j)+f(8,n-1,j)))/rhoo2
f(3,n-1,j)=f(1,n-1,j)-2.*rhoo2*ue/3.
f(7,n-1,j)=f(5,n-1,j)-ue*rhoo2/6.+0.5*(f(2,n-1,j)-f(4,n-1,j))
f(6,n-1,j)=f(8,n-1,j)-ue*rhoo2/6.-0.5*(f(2,n-1,j)-f(4,n-1,j)) 
end do
!chukou
!!(404,1)下角点
!ff(3,404,1)=f(1,404,1)
!ff(2,404,1)=f(4,404,1)
!ff(6,404,1)=f(8,404,1)
!!(1,404)下角点
!ff(3,404,404)=f(1,404,404)
!ff(4,404,404)=f(2,404,404)
!ff(7,404,404)=f(5,404,404)
do i=1,n-1
    do j=1,m-1
if(s(i,j).eq.1) then
ff(1,i,j)=f(3,i,j)
ff(3,i,j)=f(1,i,j)
ff(2,i,j)=f(4,i,j)
ff(4,i,j)=f(2,i,j)
ff(5,i,j)=f(7,i,j)
ff(7,i,j)=f(5,i,j)
ff(6,i,j)=f(8,i,j)
ff(8,i,j)=f(6,i,j)
endif
ENDDO
ENDDO
return
    end
subroutine getrhouv(f,rho,uv,cx,cy,n,m,velocity,alpha,visco,s)
real f(0:8,0:n,0:m),rho(0:n,0:m),uv(0:1,0:n,0:m),cx(0:8),cy(0:8),velocity(0:n,0:m)
integer s(0:n,0:m)
real usum,vsum,ssum,alpha,visco
cx(:)=(/0.0,1.0,0.0,-1.0,0.0,1.0,-1.0,-1.0,1.0/)
cy(:)=(/0.0,0.0,1.0,0.0,-1.0,1.0,1.0,-1.0,-1.0/)
do j=1,m-1
do i=2,n-1!入口处不回归宏观量
usum=0.0
vsum=0.0
ssum=0.0
do k=0,8
ssum=ssum+f(k,i,j)
usum=usum+f(k,i,j)*cx(k)
vsum=vsum+f(k,i,j)*cy(k)
end do
rho(i,j)=ssum
uv(0,i,j)=usum/rho(i,j)
uv(1,i,j)=vsum/rho(i,j)
if(s(i,j).eq.1)then 
    uv(0,i,j)=0
    uv(1,i,j)=0
    rho(i,j)=1.0
    endif 
end do
end do
uin=0.0
do j=0,m
     uin=uin+uv(0,5,j)
end do

!uin2=uin/n
!Re=uin2*m/visco
!print*, uin2,Re


!===0=====uv=0
do i=1,n-1
do j=1,m-1
!if(s(i,j).eq.1)then
!u(i,j)=0.0
!v(i,j)=0.0
!end if
velocity(i,j)=sqrt(uv(0,i,j)*uv(0,i,j)+uv(1,i,j)*uv(1,i,j))
end do
end do


!==========================end the fiest velocuty


return
end
subroutine outresult(kk,n,m,uv,rho,velocity,s,error1)
integer kk,n,m
CHARACTER*20::NAME,FNAME
real uv(0:1,0:n,0:m),rho(0:n,0:m),velocity(0:n,0:m)
integer s(0:n,0:m)
   print*,"时间=",kk,"error=",error1
   write(name,'(i7)') INT(kk)
fname='out/file'//TRIM(ADJUSTL(NAME))//'.plt'
OPEN(4,FILE=fname)
!real u(0:n,0:m),v(0:n,0:m),th(0:n,0:m)
!open(5, file='streamf1.0.dat')
!open(7,file='tprof1.0.dat')
write(4,*)'VARIABLES ="X" "Y" "U" "V" "RHO" "Velocity" "s"'
write(4,*)'ZONE solutiontime=',kk,',I=',n-1,',J=',m-1,'T=POINTS'
do j=1,m-1
    DO I=1,N-1
     write(4,*) I,J, uv(0,i,j),uv(1,i,j),rho(i,j),velocity(i,j),s(i,j)    
    ENDDO
ENDDO
 CLOSE(4)
     return
    end
 subroutine Error(uv,uuvv,n,m,kk,error1)
real uv(0:1,0:n,0:m),uuvv(0:1,0:n,0:m),error1
 temp1=0.
temp2=0.
error1=0
do i=1,n-1
do j=1,m-1
    temp1=temp1+  (uv(0,i,j)-uuvv(0,i,j))  *   (uv(0,i,j)-uuvv(0,i,j))+(uv(1,i,j)-uuvv(1,i,j))*(uv(1,i,j)-uuvv(1,i,j))
    temp2=temp2+(uv(0,i,j)*uv(0,i,j)+uv(1,i,j)*uv(1,i,j))
end do
end do

temp1=sqrt(temp1)
temp2=sqrt(temp2)
error1=temp1/(temp2+1e-30)
print *,"时间=",kk ,"误差=",error1

 
 return
 end