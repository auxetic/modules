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

    type(tpcorr) :: kvcorr

contains

    subroutine init_corr( tcorr, nbins, vmin, vmax )
        implicit none

        type(tpcorr) :: tcorr
        integer :: nbins
        real(8) :: vmin, vmax

        tcorr%nbins = nbins
        tcorr%vmin = vmin
        tcorr%vmax = vmax
        tcorr%wbin = (vmax-vmin) / nbins
        allocate( tcorr%cum_corr(nbins), tcorr%cum_n(nbins) );
        
    end subroutine

    Pure function calc_xi( this, i ) result(x)
        implicit none

        class(tpcorr), intent(in) :: this
        integer,      intent(in) :: i
        real(8) :: x

        x = this%vmin + ( i - 0.5 ) * this%wbin
    end function

    Pure function calc_corr( this, i ) result( x )
        implicit none

        class(tpcorr), intent(in) :: this
        integer,      intent(in) :: i
        real(8) :: x

        x = 0
        if ( this%cum_n(i) /= 0 ) then
            x = this%cum_corr(i) / this%cum_n(i)
        end if
    end function

end module
