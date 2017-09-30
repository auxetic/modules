module mo_corr
    implicit none

    type tpcorr
        integer :: nbins
        real(8) :: wbin
        real(8) :: vmin, vmax
        real(8), allocatable, dimension(:) :: cum_n
        real(8), allocatable, dimension(:) :: cum_corr
    contains
        procedure :: x    => calc_xi
        procedure :: corr => calc_corr
    end type

    type tpcorrab
        real(8), allocatable, dimension(:) :: a, b
        real(8) :: corrab
    contains
        procedure :: corr => calc_corrab
    end type

    type(tpcorr)   :: kvcorr
    type(tpcorrab) :: corrab

contains

    subroutine init_corr( tcorr, nbins, vmin, vmax )
        implicit none

        type(tpcorr) :: tcorr
        integer :: nbins
        real(8) :: vmin, vmax

        tcorr%nbins = nbins
        tcorr%vmin  = vmin
        tcorr%vmax  = vmax
        tcorr%wbin  = (vmax-vmin) / nbins
        allocate( tcorr%cum_corr(nbins), tcorr%cum_n(nbins) );
    end subroutine

    pure function calc_xi( this, i ) result(re)
        implicit none

        class(tpcorr), intent(in) :: this
        integer,       intent(in) :: i
        real(8) :: re

        re = this%vmin + ( i - 0.5 ) * this%wbin
    end function

    pure function calc_corr( this, i ) result(re)
        implicit none

        class(tpcorr), intent(in) :: this
        integer,       intent(in) :: i
        real(8) :: x

        re = 0
        if ( this%cum_n(i) /= 0 ) then
            re = this%cum_corr(i) / this%cum_n(i)
        end if
    end function

    pure function calc_corrab( this ) result( re )
        use mo_math, only: corr
        implicit none

        class(tpcorrab), intent(in) :: this
        real(8) :: re

        re = corr(this%a, this%b)
    end function

end module
