program aa
IMPLICIT NONE
integer i,j,jmax,k
integer, parameter :: imax=1000,kmax=1000
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
