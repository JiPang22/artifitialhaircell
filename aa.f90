program aa
integer i,j,jmax
integer, parameter :: imax=1000000
real :: t,x,y,dx,dy,sumi,sumr,om,z1,z2,u1,u2,F,xa,dxa,eta,omext,dom
real, dimension(imax) :: xt,noise
real, parameter :: dt=1.e-2,tau=.1,gam=0.1

!!!!!!!!!!!!!!!!!!!!!!!!노이즈 생성 !!!!!!!!!!!!!!!!!!
call random_seed()  
do i = 1, imax / 2
call random_number(u1)
call random_number(u2)
z1 = sqrt(-2. * log(u1)) * cos(2.*3.14 * u2)
z2 = sqrt(-2. * log(u1)) * sin(2.*3.14 * u2)
noise(2 * i - 1) = z1
noise(2 * i) = z2
end do
!!!!!!!!!!!!!!!!!!!!!노이즈 생성 끝!!!!!!!!!!!!!!!!!!

t=0.
x=1.
y=0.
xa=0.
eta=1.
F=0.01*eta
omext=0.5

do i=1,imax
xt(i)=x
dy=-x+(1./2.)*eta*sign(1.,x-xa)-gam*y+(0.001)*noise(i)+F*sin(omext*t)
dx=y
dxa=(x-xa)/tau
y=y+dy*dt
x=x+dx*dt
xa=dxa*dt
enddo

open(1,file='bb')
sumr=0.
sumi=0.
dom=6.28/(imax*dt)
jmax=int(5/dom)
do j=1,jmax
om=6.28*j/(imax*dt)
t=0.
write(1,*) om,sqrt(sumi**2+sumr**2)
sumr=0.
sumi=0.

do i=1,imax
t=i*dt
sumi=sumi+xt(i)*sin(om*t)*dt
sumr=sumr+xt(i)*cos(om*t)*dt
enddo
enddo
end
