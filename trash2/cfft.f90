program fft_example
    implicit none
    integer, parameter :: dp = kind(1.0d0)
    integer, parameter :: n = 10000  ! 데이터 포인트 수, 필요에 따라 조정

    real(dp) :: time(n), displacement(n)
    real(dp) :: freq(n), amplitude(n)
    integer :: i

    ! FFT 라이브러리를 포함합니다.
    include 'fftw3.f03'

    ! FFTW 계획과 배열 선언
    complex(dp), allocatable :: in(:), out(:)
    type(fftw_plan) :: plan

    ! 파일 입출력을 위한 변수
    character(len=256) :: line
    character(len=10) :: file_name
    integer :: unit_in, unit_out, ios

    ! 파일 열기 및 데이터 읽기
    open(unit=unit_in, file='data.csv', status='old', action='read')
    do i = 1, n
        read(unit_in, '(A)', iostat=ios) line
        if (ios /= 0) exit
        read(line,*) time(i), displacement(i)
    end do
    close(unit_in)

    ! FFT 수행을 위한 메모리 할당
    allocate(in(n), out(n))

    ! 입력 데이터 설정
    do i = 1, n
        in(i) = cmplx(displacement(i), 0.0_dp)
    end do

    ! FFT 계획 생성 및 실행
    plan = fftw_plan_dft_1d(n, in, out, FFTW_FORWARD, FFTW_ESTIMATE)
    call fftw_execute(plan)

    ! 결과 처리
    do i = 1, n
        freq(i) = (i-1) / (n * (time(2) - time(1)))  ! 주파수 계산
        amplitude(i) = abs(out(i)) / n               ! 진폭 계산
    end do

    ! 결과 저장
    file_name = 'ww.dat'
    open(unit=unit_out, file=file_name, status='replace', action='write')
    do i = 1, n/2
        write(unit_out, '(F10.4, 1X, F10.4)') freq(i), amplitude(i)
    end do
    close(unit_out)

    ! FFT 계획 해제 및 메모리 해제
    call fftw_destroy_plan(plan)
    deallocate(in, out)

end program fft_example

