!>> <Dimensionless process>
!>> t0 == 1/om0
!>> x_tlida == x/l 
!>> x_tilda^dash == x^dash/l
!>> eta_tilda == eta/(m * om0^2 * l)
!>> A_tilda == (A/eta) * eta_tilda
!>> tau_a^tilda == om0 * tau_a


program aaa
IMPLICIT NONE
integer i,j,jmax,jmin,numdata,k
integer, parameter :: imax=10000,kmax=1000
real :: t,x,y,dx,dy,sumi,sumr,om,z1,z2,u1,u2,F,xa,dxa,eta,omext,dom,x_dash,dx_dash,A,tau_a,gam,dt,dom_ext,ommax,ommin
real, dimension(imax+1) :: xt,noise,noise_tilda


!>>conditons
parameter(ommax = 5.)  !>>  om_ext / om_0 == om_tilda
parameter(numdata = 10000) !>> number of data point
parameter(jmax=1+numdata) !>> j is index of A_om_ext 
real, dimension(jmax+1) :: A_omext
dom_ext=1.1/jmax

open(2,file = 'Aw')  ! >> recode Aw_ext vs w_ext

!>> condations22
tau_a = 0.1
eta = 1.
gam = 0.14
dt = 1.e-2

!>> make noise
call random_seed()  
do i = 1, imax / 2
call random_number(u1)
call random_number(u2)
z1 = sqrt(-2. * log(u1)) * cos(2.*3.14 * u2)
z2 = sqrt(-2. * log(u1)) * sin(2.*3.14 * u2)
noise(2 * i - 1) = z1
noise(2 * i) = z2
end do !>> end noise make




do k=1,10 !>> k is index of ext Force

F=0.01*eta*k !>> grow ext Force magnitude


do j = 1,jmax !>> j is index of om_ext

!special conditon
omext = j*dom_ext

!>> initial conditions
t = 0.
x = 1.
x_dash = 0.
y = 0.



do i=1,imax !>> i is time index

xt(i) = x  !>>> recode x(t)
noise_tilda(i) = 1.e-2*noise(i) !>>  lowing noie magnitude
dy = -gam*y -x +noise_tilda(i) +(1./2.)*eta*sign(1.,x-x_dash) +F*sin(omext*t)
dx = y
dx_dash = (x-x_dash)/tau_a

y = y + dy*dt
x = x + dx*dt
x_dash = x_dash +dx_dash*dt 

end do !>> i end // !>> fixed omext, recode xt

A_omext(j) = 2.*sqrt(sumi**2 + sumr**2) !>> recode A_omext

sumr = 0.
sumi = 0.
t = 0.

do i = 1,imax 

t = i*dt
sumi = sumi+xt(i)*sin(omext*t)*dt
sumr = sumr+xt(i)*cos(omext*t)*dt

end do !>> i end


end do !>> j end >> fixed omext one special condition end



do j = 1,jmax

omext = j*dom_ext
write(2,*) omext, A_omext(j) !  >> fixed ext Force.
write(*,*) omext
end do!>> j end
write(2,*) '' !>> overlap the plot
enddo!>> k end

end program
