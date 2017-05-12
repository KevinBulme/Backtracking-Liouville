program backTracking
    
    !...........................!
    !        Declarations       !
    !...........................!
    
    type Particle
        real :: x_pos
        real :: y_pos
        real :: z_pos
        real :: vx
        real :: vy
        real :: vz
        real :: q_m !Mass-to-charge ratio q/m (質量電荷比)
    end type Particle

    !The unitary volume is a cubic space with a magnetic and electric field.
    
    type UnitaryVolume
        real ::  E_field_x
        real ::  E_field_y
        real ::  E_field_z
        real ::  B_field_x
        real ::  B_field_y
        real ::  B_field_z
    end type UnitaryVolume
    
    integer :: x_max, y_max, z_max
    integer :: x_min, y_min, z_min
    integer :: x, y, z, n
    integer :: n_particles !the total number of particles to be simulated
    
    type(UnitaryVolume), dimension(:,:,:), allocatable :: studySpace
    type(Particle), dimension(:), allocatable :: particles

    type(UnitaryVolume) :: cell

    !------------------Runge Kutta declarations------------------

    ! n_equations:number of equations, nsteps:number of steps, tstep:length of steps
    !         y_rk(1): initial position, y_rk(2):initial velocity
    !         E_factor: Electric field component of the differential equation, B_factor : Magnetic field component of the differential equation
    real t_rk, tstep, y_rk(5), E_factor, B_factor
    integer    n_equations, i, j, nsteps
    

    !..........................!
    !        Parameters        !
    !..........................!

    !Runge-Kutta global parameters
    n_equations=2
    nsteps=300
    tstep=0.1

    !Let's create a 11x11x11 volume

    x_max = 10
    y_max = 10
    z_max = 10
    
    x_min = 0
    y_min = 0
    z_min = 0

    n_particles = 1;
    
    
    !....................!
    !        Init        !
    !....................!
    
    allocate ( studySpace(x_min:x_max, y_min:y_max, z_min:z_max) )
    DO x = x_min, x_max
        DO y = y_min, y_max
            DO z = z_min, z_max
                studySpace(x,y,z)%E_field_x = 1
                studySpace(x,y,z)%E_field_y = 0
                studySpace(x,y,z)%E_field_z = 0
                studySpace(x,y,z)%B_field_x = 0
                studySpace(x,y,z)%B_field_y = 0
                studySpace(x,y,z)%B_field_z = 10
            END DO 
        END DO 
    END DO

    !Let's give our particles an initial position and velocity

    allocate (particles(1:n_particles))
    DO n = 1, n_particles
        particles(n)%x_pos = 10
        particles(n)%y_pos = 10
        particles(n)%z_pos = 10
        particles(n)%vx = 1
        particles(n)%vy = 2
        particles(n)%vz = 0
        particles(n)%q_m = 1
    END DO
    
    
    Print *, "BackTracking Initialized!"


    !...........................!
    !        Runge-Kutta        !
    !...........................!

    !    Runge Kutta for set of first order differential equations

    DO n = 1, n_particles
        cell = studySpace(FLOOR(particles(n)%x_pos),FLOOR(particles(n)%y_pos),FLOOR(particles(n)%z_pos))
        !For the x direction
        y_rk(1)=particles(n)%x_pos
        y_rk(2)=particles(n)%vx
        E_factor = particles(n)%q_m * cell%E_field_x
        B_factor = particles(n)%q_m * cell%B_field_z * particles(n)%vy
        call runge_kutta(n_equations, nsteps, tstep, y_rk, E_factor, B_factor)
    END DO
    Print *, "Runge Kutta Finished!"
    
    STOP

end program BackTracking




!------------------------end of main program------------------------

    !........................!
    !        Routines        !
    !........................!

!Runge-Kutta main routine
SUBROUTINE runge_kutta(n_equations, nsteps, tstep, y_rk, E_factor, B_factor)
    IMPLICIT none

    !       declarations
    real t_rk, tstep, y_rk(5)
    integer    n_equations, i, j, nsteps
    real E_factor, B_factor
    !              open file
    OPEN(90, FILE='rungef.dat')
    WRITE (90,*) 0, y_rk(1), y_rk(2)

    !         do loop n steps of Runga-Kutta algorithm
    DO 60 j = 1, nsteps
    t_rk=j*tstep
    call rk4(t_rk, y_rk, tstep, n_equations, E_factor, B_factor)
    WRITE (90,*) t_rk, y_rk(1), y_rk(2)
    60    CONTINUE
    CLOSE(90)

END

!         fourth-order Runge-Kutta subroutine 
SUBROUTINE rk4(t_rk, y_rk, tstep, n_equations, E_factor, B_factor)
    IMPLICIT none

    !        declarations
    real DERIV, h, t_rk, tstep, y_rk(5) 
    real k1(5), k2(5),k3(5), k4(5), temp1(5), temp2(5), temp3(5)
    real E_factor, B_factor
    integer i, n_equations
    h=tstep/2.0

    DO 10 i = 1,n_equations
    k1(i) = tstep * DERIV(t_rk, y_rk, i, E_factor, B_factor)
    temp1(i) = y_rk(i) + 0.5*k1(i)
    10    CONTINUE

    DO 20 i = 1,n_equations
    k2(i) = tstep * DERIV(t_rk+h, temp1, i, E_factor, B_factor)
    temp2(i) = y_rk(i) + 0.5*k2(i)
    20    CONTINUE

    DO 30 i = 1,n_equations
    k3(i) = tstep * DERIV(t_rk+h, temp2, i, E_factor, B_factor)
    temp3(i) = y_rk(i) + k3(i)
    30    CONTINUE

    DO 40 i = 1,n_equations
    k4(i) = tstep * DERIV(t_rk+tstep, temp3, i, E_factor, B_factor)
    y_rk(i) = y_rk(i) + (k1(i) + (2.*(k2(i) + k3(i))) + k4(i))/6.0
    40    CONTINUE

    RETURN
END

!        function which returns the derivatives
FUNCTION DERIV(t_rk, temp, i, E_factor, B_factor)
    IMPLICIT none

    !        declarations
    real DERIV, t_rk, temp(5)
    real E_factor, B_factor
    integer i

    IF (i .EQ. 1) DERIV=temp(2) !f1
    IF (i .EQ. 2) DERIV=E_factor + B_factor !f2 = q/m x E   Constant speed => DERIV=0
    RETURN
END
