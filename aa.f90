program aa
implicit none

integer i,imax
real t,v,vs,y,dt,dv,dvs,dy,gam,taua,u1,u2,z1,z2
parameter(imax=500000)
real, dimension(imax) :: noise
open(1,file='xx')


!make noise
call random_seed()  

do i = 1, imax / 2
call random_number(u1)
call random_number(u2)

z1 = sqrt(-2. * log(u1)) * cos(2.*3.14 * u2)
z2 = sqrt(-2. * log(u1)) * sin(2.*3.14 * u2)

noise(2 * i - 1) = z1
noise(2 * i) = z2
end do


gam = 0.04 
v = 1. ! v0 = 450 mV
y = 0.
t = 0. ! t0 = 10 ms
vs = 0.

dt = 1.e-2
taua = 1000. ! taua/t0

do i=0,imax

write(1,*) t,v

dy = -gam*y -v +sign(1.,v-vs)+2.e-2*noise(i)
dv = y
dvs = (v-vs)/taua

y = y +dy*dt
v = v +dv*dt
vs=vs+dvs*dt
t=t+dt

end do
end
