program main
    use mo_integral
    implicit none

    real(8) :: integration, comparison

    associate(                            &
        low_x     => set_integral%low_x,  &
        up_x      => set_integral%up_x,   &
        ntot      => set_integral%ntot    &
     )
    
    low_x = 4.d0
    up_x  = 9.d0
    ntot  = 1000

    abstract_func_h => func_sqrtx 
    ! appoint the integrand function as sqrt(), one could also give the name of the function directly: integration = integral( func_sqrtx, low_x, up_x ).
    ! Adjusting bin: integration = integral( func_sqrtx, low_x, up_x, ntot).  
    integration = integral( abstract_func_h, low_x, up_x )

    comparison = 38.d0 / 3.d0

print*, integration, comparison
        
    end associate 

end program    






