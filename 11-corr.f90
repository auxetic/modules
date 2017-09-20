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

    Pure function calc_xi( this, i ) result(x)
        implicit none

        class(tpcorr), intent(in) :: this
        integer,       intent(in) :: i
        real(8) :: x

        x = this%vmin + ( i - 0.5 ) * this%wbin
    end function

    Pure function calc_corr( this, i ) result( x )
        implicit none

        class(tpcorr), intent(in) :: this
        integer,       intent(in) :: i
        real(8) :: x

        x = 0
        if ( this%cum_n(i) /= 0 ) then
            x = this%cum_corr(i) / this%cum_n(i)
        end if
    end function

    Pure function calc_corrab( this ) result( x )
        implicit none

        class(tpcorrab), intent(in) :: this
        real(8) :: x

        integer :: n
        real(8) :: ava, avb, sigmaa, sigmab

        n = size( this%a )

        ava = sum(this%a) / n
        avb = sum(this%b) / n

        sigmaa = sqrt( sum( (this%a - ava)**2 ) / n )
        sigmab = sqrt( sum( (this%b - avb)**2 ) / n )

        x = sum( (this%a-ava)*(this%b-avb) ) / n / sigmaa / sigmab
    end function

end module
