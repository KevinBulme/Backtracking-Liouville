program backTracking
    
    !...........................!
    !        Declarations        !
    !...........................!
    
    !The unitary volume is a cubic sapce containing one particule.
    
    type UnitaryVolume
        !integer :: x_pos
        !integer :: y_pos
        !integer :: z_pos
        !integer :: time
        real ::  E_field
        real ::  B_field
    end type UnitaryVolume
    
    integer :: x_max, y_max, z_max, t_max
    integer :: x_min, y_min, z_min, t_min
    integer :: x, y, z, t
    
    type(UnitaryVolume), dimension(:,:,:,:), allocatable :: studySpace
    
    !..........................!
    !        Parameters        !
    !..........................!

    !Let's create a 101x101x101 volume observed during 11 times

    x_max = 100
    y_max = 100
    z_max = 100
    t_max = 10
    
    x_min = 0
    y_min = 0
    z_min = 0
    t_min = 0
    
    
    !....................!
    !        Init        !
    !....................!
    
    allocate ( studySpace(x_min:x_max, y_min:y_max, z_min:z_max, t_min:t_max) )
    DO x = x_min, x_max
        !Print *, x
        DO y = y_min, y_max
            DO z = z_min, z_max
                DO t = t_min, t_max
                    studySpace(x,y,z,t)%E_field = 0
                    studySpace(x,y,z,t)%B_field = 0
                END DO
            END DO 
        END DO 
    END DO
    
    
    Print *, "BackTracking!"
end program BackTracking
