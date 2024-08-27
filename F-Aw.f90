program aa
    implicit none
    
    integer i,imax,j,k
    real t,x,y,dt,dx,dx_dash,dy,u1,u2,z1,z2,sumi,sumr
    parameter(imax=1000)
    real, dimension(imax) :: noise,xt
    
    !make file
    open(1,file='xx')
    open(2,file='ww')
    open(3,file='ff')
    
    
    !make noise
    call random_seed()  
    
    do i = 1, imax / 2
    call random_number(u1)
    call random_number(u2)
    
    z1 = sqrt(-2. * log(u1)) * cos(2.*3.14 * u2)
    z2 = sqrt(-2. * log(u1)) * sin(2.*3.14 * u2)
    
    noise(2 * i - 1) = z1
    noise(2 * i) = z2
    end do

!calculate trajectory
    do i=1,imax

        write(1,*) t,x
        xt(i)=x

        
        dy = -gam*y -v +sign(1.,v-vs) +(2.e-2)*noise(i)+f*sin(w*t)  !add ext signal
        
        dv = y
        dvs = (v-vs)/taua
        
        y = y +dy*dt
        v = v +dv*dt
        vs=vs+dvs*dt
        t=t+dt
        
 end do