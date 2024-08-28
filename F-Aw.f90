program aa
implicit none
integer i,imax,j,k
real t,x,y,dt,dx,dx_dash,dy,u1,u2,z1,z2,sumi,sumr,om,x_dash
parameter(imax=20000) !!20000 // 500000
real, dimension(imax) :: noise,xt


!make file
open(1,file='xt-g01-t01')

!open(1,file='xt-g01-t100-w10')

!make noise
call random_seed()  

do i = 1, imax / 2
call random_number(u1)
call random_number(u2)

z1 = sqrt(-2. * log(u1)) * cos(2. * 3.14 * u2)
z2 = sqrt(-2. * log(u1)) * sin(2. * 3.14 * u2)

noise(2 * i - 1) = z1
noise(2 * i) = z2
end do

do k=0, 2
!con.dition
t = 0.
x = 13.6
x_dash = 0.
y = 0.
dt = 1.e-2
!k = 0.
j = 1.


!calculate trajectory
do i=1, imax
write(1,*) t, x
xt(i) = x

t = i * dt
!om = 0.9 +k * 0.2
!om = 0.1 !!om=0.1,1,10
om=0.1*10**k
!!gam=0.1
dy = -.1 * y -x +sign(1., x - x_dash) +(2.e-2) * noise(i) +0.01 * j * sin(om * t)
dx = y
dx_dash = (x - x_dash) / .1   !!tau = 0.1 // 100

y = y +dy*dt
x = x +dx*dt
x_dash  = x_dash + dx_dash * dt
end do
write(1,*) ''
enddo
end
