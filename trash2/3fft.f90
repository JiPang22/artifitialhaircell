program read_csv
implicit none
integer, parameter :: n = 10002 ! n: rows, m: columns
integer :: i, j
real :: t(n), x(n), om, dt, sumr, sumi


open(1, file='data.csv', status='old')

do i = 1, n
read(1,*) t(i),x(i)
end do

close(1)

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!DFT!!!!!!!!!!!!!!!!!!!!!!!!!!
    open(2, file='ww')
    dt = 0.1
    do j = 1, n
        om = j / (n * dt)
        !reset
        sumr = 0.0
        sumi = 0.0
        !sum
        do i = 1, n
            sumr = sumr + x(i) * cos(om * t(i)) * dt
            sumi = sumi - x(i) * sin(om * t(i)) * dt
        end do
        write(2, *) om, sqrt(sumi**2 + sumr**2)
        write(*, *) j
    end do
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!DFT 끝!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program read_csv

