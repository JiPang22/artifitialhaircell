program aaa
IMPLICIT NONE
integer i,j,jmax,k,jmin,numdata
integer, parameter :: imax=10000,kmax=1000
real :: t,x,y,dx,dy,sumi,sumr,om,z1,z2,u1,u2,F,xa,dxa,eta,omext,dom,x_dash,dx_dash,A,tau_a,gam,dt,dom_ext,ommax,ommin
real, dimension(imax) :: xt,noise,A_omext,noise_tilda


!conditons
!>> w_ext == [0.9, 1.1]

numdata=1000 !>> number of data point

dom_ext = (ommax-ommin)/numdata  !>> delta om_ext
jfirst = 0.9/dom_ext





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


!>> t0 == 1/om0
!>> x_tlida == x/l 
!>> x_tilda^dash == x^dash/l
!>> eta_tilda == eta/(m * om0^2 * l)
!>> A_tilda == (A/eta) * eta_tilda
!>> tau_a^tilda == om0 * tau_a


!>> condations
tau_a = 0.1
eta = 1.
gam = 0.14
dt = 1.e-2


do j = jfirst,jmax
!special conditon
omext = j*dom_ext

!>> initial conditions
t = 0.
x = 1.
x_dash = 0.
y = 0.


do i=1,imax
xt(i) = x
noise_tilda(i) = 1.e-2*noise(i)
dy = -gam*y -x -noise_tilda(i) +(1./2.)*eta*sign(1.,x-x_dash) +A*sin(omext*t)
dx = y
dx_dash = (x-x_dash)/tau_a

y = y + dy*dt
x = x + dx*dt
x_dash = x_dash +dx_dash*dt 
end do !>> time simul end // !>> fixed omext >> recode xt >> i end


!do j=1,jmax
!om=6.28*j/(imax*dt)
!write(1,*) omext,sqrt(sumi**2+sumr**2)
A_omext(j) = 2.*sqrt(sumi**2 + sumr**2) !>> recode A_omext

sumr = 0.
sumi = 0.
t = 0.

do i = 1,imax
t = i*dt
sumi = sumi+xt(i)*sin(omext*t)*dt
sumr = sumr+xt(i)*cos(omext*t)*dt
enddo !>> time integral end >> i end

!enddo
!write(1,*) ''
!enddo

end do !>> j end >> fixed omext one special condition end


open(2,file = 'Aw')  ! >> recode Aw_ext vs w_ext
do j = 1,50
omext = 0.1 * j
write(2,*) omext, A_omext(j) !  >> fixed ext Force.
end do
end program