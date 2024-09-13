program aa
implicit none
integer i,imax
parameter(imax=1000)
real t,dt,x,dx,y,dy,gam
real, dimension(imax) :: xt


!system condition
dt=1.e-2
gam=0.1


!initial condition
open(1,file='aa')
t=0.
x=0.
y=1.


do i=1,imax
write(1,*) t,x
xt(i)=x
dx=y
dy=-gam*y-x
x=x+dx*dt
y=y+dy+dt
t=i*dt
enddo
end
