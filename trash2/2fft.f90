program read_csv
    implicit none
    integer, parameter :: n = 10002, m = 2 ! n: rows, m: columns
    real :: col1(n), col2(n)
    real :: temp(m)
    integer :: i
    character(len=200) :: line
    character(len=100) :: filename
real tt(n),xt(n), om,dt,sumr,sumi
integer j
    filename = 'data.csv'

    open(unit=10, file=filename, status='old', action='read')

    do i = 1, n
        read(10, '(A)') line  ! 파일에서 한 줄을 문자열로 읽음
        read(line,*) temp     ! 문자열을 실수형 배열로 변환

        col1(i) = temp(1)
        col2(i) = temp(2)
    end do

    close(10)

    ! 배열 내용 출력
  !  print *, 'First few entries of Column 1:', col1(1:10)
   ! print *, 'First few entries of Column 2:', col2(1:10)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!DFT!!!!!!!!!!!!!!!!!!!!!!!!!!
open(2,file='ww')
dt=0.1
do j=1,n
om=j/(n*dt)
!reset
sumr=0.
sumi=0.
!sum
do i=1,n
tt(i)=col1(i)
xt(i)=col2(i)
sumr=sumr+xt(i)*cos(om*tt(i))*dt
sumi=sumi-xt(i)*sin(om*tt(i))*dt
enddo
write(2,*) om, sqrt(sumi**2+sumr**2)
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!DFT 끝!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end
