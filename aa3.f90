program aa
implicit none

integer i,imax,j,k,l
real t,v,vs,y,dt,dv,dvs,dy,gam,taua,u1,u2,z1,z2
real sumi,sumr,om,tmax,w,f
parameter(tmax=500.,dt=1.e-2, taua=0.1,gam=0.1)
parameter(imax=int(tmax/dt))
real, dimension(imax) :: noise,xt





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!make file
open(2,file='ww')


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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





f=0.01
w=sqrt(1.-0.01/4)


do k=1,5

!initial condition
v = 1.   ! v/v0 
y = 0.   ! dv/dt = y
t = 0.   ! t/t0
vs = 0.  ! adaptation displacement voltage


!calculate trajectory
do i=1,imax
xt(i)=v
dy = -gam*y -v +sign(1.,v-vs) +(2.e-2)*noise(i)+f*sin(w*t)  !add ext signal
dv = y
dvs = (v-vs)/taua
y = y +dy*dt
v = v +dv*dt
vs=vs+dvs*dt
t=t+dt
end do



!DFT to Xw - w
do j=1,imax/50
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
write(2,*) ''
f=f+0.02
enddo
end
