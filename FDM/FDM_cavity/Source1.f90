!有限差分法源程序：  
! this program is written to solve a 2-d driven flow problem in a square  
! cavity using the vorticity-stream function method. 
! l1   x方向的节点数 
! m1   y方向的节点数 
! dx   x方向的步长 
! dy   y方向的步长
!phai   流函数 
!omega  涡量 
!re     雷诺数 
!its    迭代次数  
!err1，err2    前后迭代步各点涡量之差   
program  main   
implicit none 
real dx,dy  
real err1,err2,err3  
integer,parameter ::l1=81,m1=81 
real::re=100.0 
integer i,j,k  
integer its             !计算涡量时的跌代数 
integer itss            !计算压力时的跌代数 
integer n               !内迭代次数 
integer ::itsmax=5000   !最大跌代数 
real::eps=1e-6          !迭代精度 
real::alpha1=0.3        !松弛因子 
real::alpha2=0.3  
real u(l1,m1),v(l1,m1)  !速度分量 
real p(l1,m1)           !压力分布 
real xl(l1,m1),ym(l1,m1)!存放每个点的横坐标和纵坐标 
real phai(l1,m1) 
real omega(l1,m1)  
real aw(l1,m1),ae(l1,m1),as(l1,m1),an(l1,m1) !涡量方程中的各项系数 
real ome_k(l1,m1)       !存放每一迭代步得到的涡量值 
real phai_k(l1,m1)      !存放每一迭代步得到的流函数值  
real p_k(l1,m1)         !存放压力泊松方程迭代时上一步的值 
double precision time_start,time_end   
! *...............................................................................
call cpu_time(time_start)  
call meshing(l1,m1,dx,dy,phai,omega,u,v,p,xl,ym) 
err1=1.0 
err2=1.0 
its=0 
do while(err1>=eps.or.err2>=eps.or.its<itsmax) 
    do i=1,l1                               
        do j=1,m1  
            ome_k(i,j)=omega(i,j)        !用ome_k储存前一迭代步的计算结果     
            phai_k(i,j)=phai(i,j)        !用phai_k储存前一迭代步的计算结果 
        end do 
    end do  
    do i=2,l1-1                       !涡量方程系数的更新 
        do j=2,m1-1     
            aw(i,j)=1.0+re*u(i,j)*dx/2.0    
            ae(i,j)=1.0-re*u(i,j)*dx/2.0    
            as(i,j)=1.0+re*v(i,j)*dy/2.0    
            an(i,j)=1.0-re*v(i,j)*dy/2.0  
        end do 
    end do    
    do n=1,200                     
        do j=m1-1,2,-1          !解涡量方程 
            do i=2,l1-1         
                omega(i,j)=(aw(i,j)*omega(i-1,j)+ae(i,j)*omega(i+1,j)&                  
                +as(i,j)*omega(i,j-1)+an(i,j)*omega(i,j+1))/4.0 
                omega(i,j)=ome_k(i,j)+alpha1*(omega(i,j)-ome_k(i,j))   
            end do 
        end do  
    end do  
    write(*,*) 'omega(3,11)=',omega(3,11) 
    do n=1,200                   
        do j=m1-1,2,-1           !解流函数方程  
            do i=2,l1-1          
                phai(i,j)=(phai(i-1,j)+phai(i+1,j)&            
                +phai(i,j-1)+phai(i,j+1)-omega(i,j)*dx**2)/4.0  
                phai(i,j)=phai_k(i,j)+alpha2*(phai(i,j)-phai_k(i,j))      
            end do 
        end do 
    end do  
    write(*,*) 'phai(3,11)=',phai(3,11)    
    do i=2,l1-1                    !由流函数计算速度分布      
        do j=2,m1-1       
            u(i,j)=(phai(i,j+1)-phai(i,j))/dy+(2.0*phai_k(i,j)-phai_k(i,j+1)-phai_k(i,j-1))/(2.0*dy)       
            v(i,j)=(phai(i-1,j)-phai(i,j))/dy+(2.0*phai_k(i,j)-phai_k(i+1,j)-phai_k(i-1,j))/(2.0*dx)      
        end do     
    end do  
    do j=1,m1                    !左、右边界的涡量值 
        omega(1,j)=2.0*(phai(2,j)-phai(1,j))/dx**2       
        omega(l1,j)=2.0*(phai(l1-1,j)-phai(l1,j))/dx**2 
    end do  
    do i=2,l1-1                   !上、下边界的涡量值      
        omega(i,m1)=2.0/dy              
        omega(i,1)=2.0*(phai(i,2)-phai(i,1))/dy**2 
    end do   
    omega(1,1)=(omega(1,2)+omega(2,1))/2.0   
    omega(1,m1)=(omega(1,m1-1)+omega(2,m1))/2.0   
    omega(l1,1)=(omega(l1-1,1)+omega(2,l1))/2.0   
    omega(l1,m1)=(omega(l1,m1-1)+omega(l1-1,m1))/2.0 
    its=its+1   
    err1=maxval(abs(ome_k-omega))   
    err2=maxval(abs(phai_k-phai))  
    write(*,*) 'err1=',err1  
    write(*,*) 'err2=',err2  
    write(*,*) 'its=',its           
end do  
!  *......................................................................................................  
err3=1.0   
itss=0                           !用计算得到的速度求解内点压力场分布                  
do while(err3>=eps.and.itss<=itsmax) 
    do i=1,l1    
        do j=1,m1      
            p_k(i,j)=p(i,j)    
        end do  
    end do   
    do i=2,l1-1                        
        do j=2,m1-1    
            p(i,j)=(p(i-1,j)+p(i+1,j)+p(i,j+1)+p(i,j-1)-0.25*((u(i+1,j)-u(i-1,j))*(v(i,j+1)-v(i,j-1))-& 
            (v(i+1,j)-v(i-1,j))*(u(i,j+1)-u(i,j-1))))/4.0          
            p(i,j)=p_k(i,j)+0.3*(p(i,j)-p_k(i,j))   
        end do  
    end do  
    !  *............................................................................................................    
    do j=2,m1-1         !压力场的左、右边界值      
        p(1,j)=4.0/3.0*p(2,j)-1.0/3.0*p(3,j)-(2.0/(3.0*re))*(-2.0*u(2,j)+u(3,j))/dx     
        p(l1,j)=4.0/3.0*p(l1-1,j)-1.0/3.0*p(l1-2,j)+(2.0/(3.0*re))*(-2.0*u(l1-1,j)+u(l1-2,j))/dx     
    end do    
    do i=2,l1-1                    
        !压力场的上、下边界值     
        p(i,m1)=(4.0/3.0)*p(i,m1-1)-(1.0/3.0)*p(i,m1-2)+(2.0/(3.0*re*dy))*(-2.0*v(i,m1-1)+v(i,m1-2))  
        p(i,1)=4.0/3.0*p(i,m1-1)-1.0/3.0*p(i,m1-2)-(2.0/(3.0*re))*(-2.0*v(i,2)+v(i,3))/dx    
    end do     
    !角点上的压力进行加权平均   
    p(1,1)=0.5*p(1,2)+0.5*p(2,1)   
    p(1,m1)=0.5*p(1,m1-1)+0.5*p(2,m1)   
    p(l1,1)=0.5*p(l1,2)+0.5*p(l1-1,1)   
    p(l1,m1)=0.5*p(l1,m1-1)+0.5*p(l1-1,m1)
    itss=itss+1  
    err3=maxval(abs(p_k-p))  
end do  
call cpu_time(time_end)  
write(*,*)'总共运行时间(s)',time_end-time_start  
!  *..........................................................................................................  
open(1,file='phai.plt')  
open(2,file='omega.plt')  
open(3,file='v.plt')  
open(4,file='u.plt')  
open(5,file='pressure.plt')  
write(1,*) 'title="phai"'  
write(1,*) 'variable="x","y","phai"'  
write(1,*) 'zone t="s"', 'i=',l1, 'j=',m1, 'c=black' 
write(2,*) 'title="omega"'  
write(2,*) 'variable="x","y","omega"'  
write(2,*) 'zone t="s"', 'i=',l1, 'j=',m1, 'c=black'  
write(3,*) 'title="v"'  
write(3,*) 'variable="x","y","v"'  
write(3,*) 'zone t="s"', 'i=',l1, 'j=',m1, 'c=black'  
write(4,*) 'title="u"'  
write(4,*) 'variable="x","y","u"'  
write(4,*) 'zone t="s"', 'i=',l1, 'j=',m1, 'c=black'  
write(5,*) 'title="pressure"'  
write(5,*) 'variable="x","y","pressure"'  
write(5,*) 'zone t="s"', 'i=',l1, 'j=',m1, 'c=black'     
do i=1,l1       
    do j=1,m1        
        write(1,*) xl(i,j),ym(i,j), phai(i,j)       
        write(2,*) xl(i,j),ym(i,j), omega(i,j)       
        write(3,*) xl(i,j),ym(i,j), v(i,j)      
        write(4,*) xl(i,j),ym(i,j), u(i,j)       
        write(5,*) xl(i,j),ym(i,j), p(i,j)       
    end do     
end do  
close (1)  
close (2)  
close (3)  
close (4)  
close (5)  
open(6,file='v_x.txt')  
open(7,file='u_y.txt')  
do i=1,l1       
    write(6,*) v(i,(m1+1)/2.0)  
end do   
do j=1,m1    
    write(7,*) u((l1+1)/2.0,j)  
end do  
close (6)  
close (7)  
!  *............................................................................................................  
end program
!网格划分并且设定初场  
subroutine meshing(l,m,x,y,phi,ome,u_x,v_y,p0,xl,ym)    
integer l,m   
real x,y,phi(l,m),ome(l,m),u_x(l,m),v_y(l,m),p0(l,m),xl(l,m),ym(l,m)  
x=1.0/(l-1)  
y=1.0/(m-1)  
xl(1,:)=0  
xl(l,:)=1  
ym(:,1)=0  
ym(:,m)=1  
do i=1,l    
    do j=1,m   
        xl(i,j)=xl(1,j)+(i-1)*x   
        ym(i,j)=ym(i,1)+(j-1)*y    
    end do  
end do   
do i=1,l    !流函数、压力初场    
    do j=1,m      
        phi(i,j)=0   
        p0(i,j)=1.0     
    end do  
end do  
do i=1,l     !涡量初场    
    do j=1,m-1    
        ome(i,j)=0    
    end do  
end do  
do i=1,l   
ome(i,m)=2.0/y  
end do              
do i=1,l         !速度初场    
    u_x(i,1)=0   
    v_y(i,1)=0   
    u_x(i,m)=1.0    
    v_y(i,m)=0  
end do  
do j=2,m-1    
    do i=1,l    
        u_x(i,j)=0    
        v_y(i,j)=0   
    end do  
end do  
end subroutine  