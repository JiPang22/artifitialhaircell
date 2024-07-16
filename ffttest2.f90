program read_csv
    implicit none
    integer, parameter :: n = 10002, m = 2 ! n: rows, m: columns
    real :: col1(n), col2(n)
    real :: temp(m)
    integer :: i
    character(len=200) :: line
    character(len=100) :: filename

    filename = '44.csv'

    open(unit=10, file=filename, status='old', action='read')

    do i = 1, n
        read(10, '(A)') line  ! 파일에서 한 줄을 문자열로 읽음
        read(line,*) temp     ! 문자열을 실수형 배열로 변환

        col1(i) = temp(1)
        col2(i) = temp(2)
    end do

    close(10)

    ! 배열 내용 출력
    print *, 'First few entries of Column 1:', col1(1:10)
    print *, 'First few entries of Column 2:', col2(1:10)

end program read_csv
