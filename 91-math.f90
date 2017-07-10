module mo_math
    implicit none

contains

    function randperm(n) result( redata )
        use ifport
        implicit none
        integer, intent(in) :: n
        integer, dimension(n) :: redata
        integer :: i, itemp
        redata(1) = 1
        do i=2, n
            itemp = floor(rand(0)*i) + 1
            if ( i /= itemp ) then
                redata(i) = redata(itemp)
                redata(itemp) = i
            else
                redata(i) = i
            end if
        end do
    end function randperm

end module

