program aa
implicit none
integer i,imax,j,k,l
real t,v,vs,y,dt,dv,dvs,dy,gam,taua,u1,u2,z1,z2
real sumi,sumr,om,tmax,w,f
parameter(tmax=1000.,dt=1.e-2, taua=0.1,gam=0.1)
parameter(imax=int(tmax/dt))
real, dimension(imax) :: noise,xt

!!make file
open(1,file='xx')
open(2,file='ww')
open(3,file='ff')


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


!intial condital
w=0. !omega/omega0

do l=0,3
do k=10,200
f=k*5.e-3 !f is ext force magnitude

!initial condition
v = 1.   ! v/v0 
y = 0.   ! dv/dt = y
t = 0.   ! t/t0
vs = 0.  ! adaptation displacement voltage


!calculate trajectory
do i=1,imax

write(1,*) t,v
xt(i)=v
!dy = -gam*y -v +sign(1.,v-vs) +(2.e-2)*noise(i)	!original
!dy = -gam*y -v +(2.e-2)*noise(i) !adaptation force off
dy = -gam*y -v +sign(1.,v-vs) +(2.e-2)*noise(i)+f*sin(w*t)  !add ext signal

dv = y
dvs = (v-vs)/taua

y = y +dy*dt
v = v +dv*dt
vs=vs+dvs*dt
t=t+dt

end do


!reset
t=0.
sumr=0.
sumi=0.

!DFT to sensitivity - F
do i=1,imax
t=t+dt
sumr=sumr+xt(i)*cos(w*t)*dt
sumi=sumi-xt(i)*sin(w*t)*dt
enddo
write(3,*) f,2*sqrt(sumi**2+sumr**2)/f
enddo
write(3,*) ''
w=0.97+l*0.005
enddo






!DFT to Xw - w
do j=1,200
!om=6.28*j/tmax
om=j/tmax

!reset
t=0.
sumr=0.
sumi=0.


do i=1,imax
t=t+dt
sumr=sumr+xt(i)*cos(om*t)*dt
sumi=sumi-xt(i)*sin(om*t)*dt
enddo

write(2,*) om, sqrt(sumi**2+sumr**2)

enddo
end
