module mo_math
    implicit none

contains

    function randperm(n) result(redata)
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
    end function

    function randuvec(n) result(uvec)
        implicit none

        integer, intent(in) :: n
        real(8) :: uvec(n)

        call random_number(uvec)
        uvec = 2 * uvec - 1.d0
        do while ( norm2(uvec) > 1.d0 )
            call random_number(uvec)
            uvec = 2 * uvec - 1.d0
        end do
        
        uvec = uvec / norm2(uvec)
    end function

    subroutine swapr(x, y)
        implicit none

        real(8), intent(inout) :: x, y
        real(8) :: tmp

        tmp = x
        x   = y
        y   = tmp
    end subroutine

    subroutine swapi(x, y)
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
        integer :: ilen
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

    pure function dot(a,b) result(re)
        implicit none

        real(8), intent(in), dimension(:) :: a, b
        real(8) :: re

        re = sum(a*b)
    end function

    pure function times2(a,b) result(re)
        implicit none

        real(8), intent(in)  :: a(2), b(2)
        real(8) :: re

        re = a(1) * b(2) - a(2) * b(1)
    end function

    pure function times3(a,b) result(re)
        implicit none

        real(8), intent(in)  :: a(3), b(3)
        real(8) :: re(3)

        re(1) = a(2) * b(3) - a(3) * b(2)
        re(2) = a(3) * b(1) - a(1) * b(3)
        re(3) = a(1) * b(2) - a(2) * b(1)
    end function

    function unitv(vector) result(uvector)
        implicit none

        real(8), dimension(:) :: vector
        real(8), allocatable, dimension(:) :: uvector

        uvector = vector / norm2(vector)
    end function

end module
