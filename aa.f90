program aa
implicit none

integer i,imax,j
real t,v,vs,y,dt,dv,dvs,dy,gam,taua,u1,u2,z1,z2,sumi,sumr,om,tmax
parameter(tmax=1500.,dt=1.e-2, taua=0.1)
parameter(imax=int(tmax/dt))
real, dimension(imax) :: noise,xt



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


gam = 0.1 ! gam/w0 
v = 1. ! v0 = 450 mV
y = 0. !dv/dt = y
t = 0. ! t0 = 10 ms
vs = 0.


open(1,file='xx')
do i=0,imax

write(1,*) t,v

!dy = -gam*y -v +sign(1.,v-vs) +2.e-2*noise(i)
dy = -gam*y -v +2.e-2*noise(i) !forceoff
dv = y
dvs = (v-vs)/taua

y = y +dy*dt
v = v +dv*dt
vs=vs+dvs*dt
t=t+dt

end do
open(2,file='wx')
do j=1,50
om=6.28*j/(imax*dt)

!reset
t=0.
sumr=0.
sumi=0.
!sum
do i=1,imax
t=t+dt
sumr=sumr+xt(i)*cos(om*t)*dt
sumi=sumi-xt(i)*sin(om*t)*dt
enddo
write(2,*) om, sqrt(sumi**2+sumr**2)
enddo
end
