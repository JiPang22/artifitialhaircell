program dft_example
    implicit none
    integer, parameter :: dp = kind(1.0d0)
    integer, parameter :: n = 10002  ! 데이터 포인트 수, 필요에 따라 조정

    real(dp) :: t(n), x(n)
    real(dp) :: om, sumr, sumi, dt
    integer :: i, j
    character(len=256) :: line
    integer :: unit_in, unit_out, ios

    ! 파일 입출력을 위한 변수
    character(len=10) :: file_name_in, file_name_out
    file_name_in = 'data.csv'
    file_name_out = 'ww.dat'

    ! 파일 열기 및 데이터 읽기
    open(unit=unit_in, file=file_name_in, status='old', action='read')
    do i = 1, n
        read(unit_in, '(A)', iostat=ios) line
        if (ios /= 0) exit
        read(line, *) t(i), x(i)
    end do
    close(unit_in)

    ! 시간 간격 계산
    dt = t(2) - t(1)

    ! 결과 저장을 위한 파일 열기
    open(unit=unit_out, file=file_name_out, status='replace', action='write')

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!DFT!!!!!!!!!!!!!!!!!!!!!!!!!!
    do j = 1, n
        om = j / (n * dt)
        !reset
        sumr = 0.0_dp
        sumi = 0.0_dp
        !sum
        do i = 1, n
            sumr = sumr + x(i) * cos(om * t(i)) * dt
            sumi = sumi - x(i) * sin(om * t(i)) * dt
        end do
        write(unit_out, '(F10.4, 1X, F10.4)') om, sqrt(sumr**2 + sumi**2)
        write(*, *) 'Frequency:', om, 'Amplitude:', sqrt(sumr**2 + sumi**2)
    end do
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!DFT 끝!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    close(unit_out)

end program dft_example

