program fft_example
    use, intrinsic :: iso_c_binding, only: c_double
    implicit none

    ! FFTW 관련 변수
    integer, parameter :: dp = selected_real_kind(15)
    integer :: n, i
    real(dp), allocatable :: data(:)
    complex(dp), allocatable :: fft_result(:)
    integer :: plan

    ! CSV 파일 읽기 관련 변수
    character(len=100) :: filename
    integer :: unit_number, ios

    ! 파일 이름 및 크기 설정
    filename = '44.csv'
    n = 500 ! 데이터 포인트의 수를 설정 (오실로스코프 데이터 크기에 맞춰 설정)

    ! 메모리 할당
    allocate(data(n))
    allocate(fft_result(n))

    ! CSV 파일 열기
    open(unit=unit_number, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) then
        print *, '파일을 열 수 없습니다: ', filename
        stop
    end if

    ! CSV 파일에서 데이터 읽기
    do i = 1, n
        read(unit_number, *, iostat=ios) data(i)
        if (ios /= 0) exit
    end do

    ! 파일 닫기
    close(unit_number)

    ! FFTW 계획 생성
    plan = fftw_plan_dft_r2c_1d(n, data, fft_result, FFTW_ESTIMATE)

    ! 푸리에 변환 수행
    call fftw_execute_dft_r2c(plan, data, fft_result)

    ! 결과 출력
    do i = 1, n
        print *, 'FFT Result ', i, ': ', fft_result(i)
    end do

    ! 메모리 해제
    call fftw_destroy_plan(plan)
    deallocate(data)
    deallocate(fft_result)

contains
    ! FFTW Fortran 인터페이스 함수 선언
    interface
        subroutine fftw_plan_dft_r2c_1d(plan, n, in, out, flags) bind(C, name='fftw_plan_dft_r2c_1d')
            import :: c_double
            integer(c_int) :: plan, n, flags
            real(c_double) :: in(*)
            complex(c_double) :: out(*)
        end subroutine fftw_plan_dft_r2c_1d

        subroutine fftw_execute_dft_r2c(plan, in, out) bind(C, name='fftw_execute_dft_r2c')
            import :: c_double
            integer(c_int) :: plan
            real(c_double) :: in(*)
            complex(c_double) :: out(*)
        end subroutine fftw_execute_dft_r2c

        subroutine fftw_destroy_plan(plan) bind(C, name='fftw_destroy_plan')
            integer(c_int) :: plan
        end subroutine fftw_destroy_plan
    end interface
end program fft_example

