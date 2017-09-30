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

    pure function mean(a) result(re)
        implicit none

        real(8), intent(in), dimension(:) :: a
        real(8) :: re

        re = sum(a) / size(a)
    end function

    pure function std(a) result(re)
        implicit none

        real(8), intent(in), dimension(:) :: a
        integer :: i, ilen
        real(8) :: mean_of_a, re

        ilen = size(a)

        mean_of_a = sum(a) / ilen 

        re = sum( (a-mean_of_a)**2 ) / ilen
        re = sqrt(re)
    end function

    pure function corr(a, b) result(re)
        implicit none

        real(8), intent(in), dimension(:) :: a, b
        real(8) :: re
        integer :: ilen

        real(8) :: sigma_of_a, sigma_of_b

        ilen = size(a)

        sigma_of_a = std(a)
        sigma_of_b = std(b)

        re = sum(a*b)/ilen - mean(a)*mean(b)
        re = re / std(a) / std(b)
    end function

end module
