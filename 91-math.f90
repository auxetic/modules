module mo_math
    implicit none

contains

    function randperm(n) result( redata )
        implicit none
        integer, intent(in) :: n
        integer, dimension(n) :: redata
        integer :: i, itemp
        real(8) :: rtemp
        redata(1) = 1
        do i=2, n
            call random_number(rtemp)
            itemp = floor(rtemp*i) + 1
            if ( i /= itemp ) then
                redata(i) = redata(itemp)
                redata(itemp) = i
            else
                redata(i) = i
            end if
        end do
    end function randperm

    subroutine swapr( x, y )
        implicit none

        real(8), intent(inout) :: x, y
        real(8) :: tmp

        tmp = x
        x   = y
        y   = tmp
    end subroutine

    subroutine swapi( x, y )
        implicit none

        integer, intent(inout) :: x, y
        integer :: tmp

        tmp = x
        x   = y
        y   = tmp
    end subroutine

    subroutine init_rand(seed)
        implicit none
            
        integer :: seed
        integer :: n
        integer, allocatable, dimension(:) :: seed_array

        call random_seed(size=n)
        allocate( seed_array(n) )
        seed_array = seed
        call random_seed(put=seed_array)
    end subroutine

end module

