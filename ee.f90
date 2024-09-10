program aa
implicit none
integer i,imax, num_peaks, peak_search_interval
parameter(imax=10000)
real t,dt,x,dx,y,dy,gam, avg_log_decrement, damping_ratio
real, dimension(imax) :: xt
real, dimension(imax) :: peaks
integer, dimension(imax) :: peak_indices

! 피크 탐색 간격 설정 (사용자가 원하는 값으로 설정)
peak_search_interval = 2  ! 피크 탐색 간격 (예: 2로 설정하여 두 포인트마다 탐색)

! 시스템 조건 설정
dt = 2.e-3
gam = 0.1

! 초기 조건 설정
open(1, file='aa')
t = 0.
x = 0.
y = 1.

! 시스템 통합 (시뮬레이션 진행)
do i = 1,imax
write(1,*) t, x
xt(i) = x
dx = y
dy = -gam*y - x
x = x + dx*dt
y = y + dy*dt
t = i*dt
enddo

! 피크 찾기 (간격 조정 가능)
num_peaks = 0
do i = 2, imax-1, peak_search_interval
if (xt(i) > xt(i-1) .and. xt(i) > xt(i+1)) then
num_peaks = num_peaks + 1
peaks(num_peaks) = xt(i)
peak_indices(num_peaks) = i
end if
end do

! 대수 감쇠율 계산
avg_log_decrement = 0.0
if (num_peaks > 1) then
do i = 1, num_peaks-1
avg_log_decrement = avg_log_decrement + log(peaks(i) / peaks(i+1))
end do
avg_log_decrement = avg_log_decrement / (num_peaks - 1)

! 감쇠비 계산
damping_ratio = avg_log_decrement / (2.0 * 3.14159265)

print *, '평균 대수 감쇠율: ', avg_log_decrement
print *, '감쇠비: ', damping_ratio
else
print *, '감쇠비 계산에 필요한 충분한 피크를 찾지 못했음.'
end if

end program aa

