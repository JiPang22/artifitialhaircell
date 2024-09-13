program aa
    use omp_lib  ! OpenMP 라이브러리 사용
    integer :: i, j, jmax, k, kmax
    integer, parameter :: imax = 100000
    real :: t, x, y, dx, dy, sumi, sumr, om, z1, z2, u1, u2, F, xa, dxa, eta, omext, dom
    real, dimension(imax) :: xt, noise
    real, parameter :: dt = 1.e-2, tau = .1, gam = 0.1

    !!!!!!!!!!!!!!!!!!!!노이즈 생성 !!!!!!!!!!!!!!!!!!
    call random_seed()
    do i = 1, imax / 2
        call random_number(u1)
        call random_number(u2)
        z1 = sqrt(-2. * log(u1)) * cos(2. * 3.14 * u2)
        z2 = sqrt(-2. * log(u1)) * sin(2. * 3.14 * u2)
        noise(2 * i - 1) = z1
        noise(2 * i) = z2
    end do
    !!!!!!!!!!!!!!!!!!!!!노이즈 생성 끝!!!!!!!!!!!!!!!!!!
    
    open(1, file = 'bb', status = 'unknown')
    kmax = 1000  ! kmax 값을 적절하게 설정하십시오

    do k = 0, kmax
        t = 0.
        x = 1.
        y = 0.
        xa = 0.
        eta = 1.
        F = 0.1 * eta
        ! F = 0.01 * eta + 0.01 * k
        ! F = 1. * eta + 1. * k
        omext = 0.01 * k

        ! OpenMP 병렬화 시작
        !$omp parallel do private(i, dx, dy, dxa) shared(x, y, xa, t, xt, noise, F, omext)
        do i = 1, imax
            xt(i) = x
            dy = -x + (1. / 2.) * eta * sign(1., x - xa) - gam * y + (0.001) * noise(i) + F * sin(omext * t)
            dx = y
            dxa = (x - xa) / tau
            y = y + dy * dt
            x = x + dx * dt
            xa = xa + dxa * dt
            t = t + dt
        end do
        !$omp end parallel do

        sumr = 0.
        sumi = 0.
        dom = 6.28 / (imax * dt)
        jmax = int(5 / dom)

        ! OpenMP 병렬화 시작
        !$omp parallel do private(j, om, t, i) shared(xt) reduction(+:sumi, sumr)
        do j = 1, jmax
            om = 6.28 * j / (imax * dt)
            t = 0.
            sumi = 0.
            sumr = 0.

            do i = 1, imax
                t = i * dt
                sumi = sumi + xt(i) * sin(om * t) * dt
                sumr = sumr + xt(i) * cos(om * t) * dt
            end do

            write(1, *) om, sqrt(sumi ** 2 + sumr ** 2)
        end do
        !$omp end parallel do

        ! 빈 줄 추가하여 결과 구분
        write(1, *) ''

    end do

    close(1)
end program aa

