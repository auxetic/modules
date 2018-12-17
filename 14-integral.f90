module mo_integral

    type tpset_integral
        real(8) ::  low_x, up_x   ! the limit of integration 
        integer ::  ntot = 10000  ! the bin of integration 
    end type
    type(tpset_integral) :: set_integral

    abstract interface
            pure function abstract_func( tx ) result( tfunc ) ! the integrand function
            real(8), intent(in)  :: tx
            real(8)              :: tfunc
        end function
    end interface

    procedure(abstract_func), pointer :: abstract_func_h => null()

contains

    function integral( abstract_func, tlow_x, tup_x, tntot ) result( integration )
        !
        !para list
        real(8), external    :: abstract_func
        real(8), intent(in)  :: tlow_x, tup_x
        integer, intent(in), optional :: tntot
        !result
        real(8) :: integration        

        !local
        integer :: tnn, ntot 

        if( present(tntot) ) then
            ntot = tntot
        else
            ntot = set_integral%ntot 
        end if
        tnn = ntot

        call get_integral( abstract_func, tlow_x, tup_x, tnn, ntot, integration )
                            
    end function 

    recursive subroutine get_integral( abstract_func, tlow_x, tup_x, tnn, tntot, return_func ) 
        implicit none
        !
        !para list
        real(8), external      :: abstract_func
        real(8), intent(in)    :: tlow_x, tup_x
        integer, intent(in)    :: tntot
        integer, intent(inout) :: tnn
        !result
        real(8), intent(inout) :: return_func

        !local
        real(8), save :: x, x2, dx, dx2
        real(8), save :: tdfunc
        
        if ( tnn == tntot ) then
            return_func = 0.d0
            x     = tlow_x
            dx    = ( tup_x - tlow_x ) / dble(tntot)
            dx2   = dx * 0.5d0        
        end if

        if ( tnn > 0 ) then
            x   = x  + dx
            x2  = x  - dx2
            tdfunc = abstract_func( x2 )
            return_func  = return_func + tdfunc * dx
            tnn = tnn - 1
            call get_integral( abstract_func, tlow_x, tup_x, tnn, tntot, return_func ) 
        end if

        if ( tnn == 0 ) then
            return          
        end if

    end subroutine

    pure function func_xbessel( tx ) result( tfunc )
        implicit none
        !
        !para list
        real(8), intent(in) :: tx
        real(8)             :: tfunc
        
        tfunc = tx * bessel_j0( tx )   
 
    end function

    pure function func_sqrtx( tx ) result( tfunc )
        implicit none
        !
        !para list
        real(8), intent(in) :: tx
        real(8)             :: tfunc
    
        tfunc = dsqrt( tx )   
 
    end function

end module
