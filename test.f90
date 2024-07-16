program fft_example
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    implicit none

    integer(c_int), parameter :: n = 1024
    real(c_double) :: data(n)
    complex(c_double) :: fft_result(n)
    type(c_ptr) :: plan
    integer :: i, ios
    character(len=*), parameter :: filename = '44.csv'
    integer :: unit_number

    open(unit=unit_number, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) stop 'cant_open: ', filename
    do i = 1, n
        read(unit_number, *, iostat=ios) data(i)
        if (ios /= 0) exit
    end do
    close(unit_number)

    plan = fftw_plan_dft_r2c_1d(n, data, fft_result, FFTW_ESTIMATE)
    call fftw_execute_dft_r2c(plan, data, fft_result)

    do i = 1, n
        print *, 'FFT Result ', i, ': ', fft_result(i)
    end do

    call fftw_destroy_plan(plan)

contains
    interface
        function fftw_plan_dft_r2c_1d(n, in, out, flags) bind(C, name='fftw_plan_dft_r2c_1d')
            import :: c_double, c_int, c_ptr
            integer(c_int) :: n, flags
            real(c_double) :: in(*)
            complex(c_double) :: out(*)
            type(c_ptr) :: fftw_plan_dft_r2c_1d
        end function fftw_plan_dft_r2c_1d

        subroutine fftw_execute_dft_r2c(plan, in, out) bind(C, name='fftw_execute_dft_r2c')
            import :: c_double, c_ptr
            type(c_ptr), value :: plan
            real(c_double) :: in(*)
            complex(c_double) :: out(*)
        end subroutine fftw_execute_dft_r2c

        subroutine fftw_destroy_plan(plan) bind(C, name='fftw_destroy_plan')
            import :: c_ptr
            type(c_ptr), value :: plan
        end subroutine fftw_destroy_plan
    end interface
end program fft_example

