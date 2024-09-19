program aaa
IMPLICIT NONE
integer i, j, jmax, jmin, k
integer, parameter :: imax=314159
real  t, x, y, dt, dx, dy, sumi, sumr, z1, z2, u1, u2, F, eta, om_ext
real x_dash, dx_dash, tau_a, gam, dom_ext, om_max, om_min
real, dimension(imax + 1) :: xt, noise, noise_tilda
parameter(om_max = 1.1)
parameter(om_min = 0.9)
parameter(dt = 1.e-2)
parameter(dom_ext = 6.28 / (imax * dt))
parameter(jmin = int(om_min / dom_ext))
parameter(jmax = int(om_max / dom_ext))
real, dimension(jmax+1) :: A_om_ext
parameter(tau_a = 0.1)
parameter(eta = 1.)
parameter(gam = 0.2)
open(2, file = 'Aw') !>> make file!


     !>>  make noise
call random_seed()  
do i = 1, imax / 2
call random_number(u1)
call random_number(u2)
z1 = sqrt(-2. * log(u1)) * cos(2. * 3.14 * u2)
z2 = sqrt(-2. * log(u1)) * sin(2. * 3.14 * u2)
noise(2 * i - 1) = z1
noise(2 * i) = z2
end do  !>>     end noise make


do k = 1, 5 !>> k is index of ext Force
F = 0.02 * eta * k       !>> grow ext Force magnitude

do j = jmin, jmax !>> j is index of om_ext
!special conditon
om_ext = j * dom_ext
!>> initial conditions
t = 0.
x = 1.
x_dash = 0.
y = 0.
do i = 1, imax !>> i is time index
xt(i) = x  !>>> recode x(t)
noise_tilda(i) = 1.e-3 * noise(i)   !>>  lowing noize magnitude
dy = - gam * y - x + noise_tilda(i) + (1. / 2.) * eta * sign(1., x - x_dash) + F * sin(om_ext * t)
dx = y
dx_dash = (x - x_dash) / tau_a
y = y + dy * dt
x = x + dx * dt
x_dash = x_dash + dx_dash * dt 
t = t + dt
end do   !>> i end // !>> fixed om_ext, recode xt
sumr = 0.
sumi = 0.
t = 0.
do i = 1, imax 
t = i * dt
sumi = sumi + xt(i) * sin(om_ext * t) * dt
sumr = sumr + xt(i) * cos(om_ext * t) * dt
end do  !>> i end
A_om_ext(j) = 2. * sqrt(sumi**2 + sumr**2) !>> recode A_om_ext
end do  !>> j end >> fixed om_ext one special condition end


    !>> recoding amp - ext om
do j = jmin, jmax
om_ext = j * dom_ext
write(2, *) om_ext, A_om_ext(j) !  >> fixed ext Force.
end do  !>>   j end
write(2, *) '' !>> overlap the plot
end do  !>>   k end
end program