program solvingwithRungeKutta
    implicit none
    real(8)::a,b,h,y_0,t

    write(*,*)"Enter the interval a,b, the value of the step-size h and the value of y_0"
    read(*,*)a,b,h,y_0
    open(unit=1, file="resultadoRungeKutta.dat",status='old')
    call RungeKutta(y_0,a,b,h)
    close(1)

    contains
    subroutine RungeKutta(y_0,a,b,h)
        implicit none
        real(8), intent(inout)::a,h,b, y_0
        real(8):: y,t, F_1, F_2, F_3, F_4
        t=a
        y=y_0
        do while(t<=b)
            F_1=h*f(y,t)
            F_2=h*f(y+(1d0/2d0)*F_1, t+h*(1d0/2d0))
            F_3=h*f(y+(1d0/2d0)*F_2, t+h*(1d0/2d0))
            F_4=h*f(y+F_3,t+h)

            y=y+(F_1+2d0*F_2+2d0*F_3+F_4)/6d0
            write(1,*)t, y 
            t=t+h
        end do
    end subroutine
    function f(y,t)
        implicit none
        real(8)::f,y,t,pi
        pi=acos(-1d0)
        f=-y+sin(pi*2d0*t)
    end function
end program
