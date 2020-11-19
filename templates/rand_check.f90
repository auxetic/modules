program main
    use mo_math
    implicit none

    integer :: seed0 = 100
    integer :: i, j
    integer, parameter :: nsample = 10000
    real(8), dimension(8,nsample) :: samples
    real(8), dimension(8) :: xy, xy1
    real(8) :: r, r1

    do i=1, nsample
        if ( mod(i,1000) == 0 ) print*, i
        call init_rand(seed0+i)
        call random_number(xy)
        samples(:,i) = xy
    end do

    ! check
    do i=1, nsample
        if ( mod(i,1000) == 0 ) print*, i
        xy = samples(:,i)
        do j=i+1, nsample
            xy1 = samples(:,j)
            r = hypot( xy1(1)-xy(1), xy1(2)-xy(2) )
            if ( r<=1e-4 ) then
                r1 = hypot( xy1(3)-xy(3), xy1(4)-xy(4) )
                if ( r1 < 1e-2 ) then
                    print*, i, j, r1
                end if
            end if
        end do
    end do


contains

end program main
