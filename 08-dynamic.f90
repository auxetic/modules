module mo_dynamic
    use mo_syst
    use mo_config
    implicit none

    type msd_t
        integer :: natom                                ! Number of particles
        integer :: nhistmax                             ! length of histdata
        integer :: ndt                                  ! sample frequency
        real(8) :: kx
        integer :: cumstep = -1
        integer :: histidx = 0
        logical :: calc_flag = .false.
        integer, allocatable, dimension(:)   :: msdcount
        real(8), allocatable, dimension(:)   :: msd     ! msd
        real(8), allocatable, dimension(:)   :: alpha2  ! non-gaussian parameter
        real(8), allocatable, dimension(:)   :: fkt     ! intemediate scatter function
        real(8), allocatable, dimension(:,:) :: histdata ! array used to stor history config
    end type

    type vcorr_t
        integer :: natom
        integer :: nhistmax
        integer :: ndt
        integer :: cumstep = -1
        integer :: histidx = 0
        logical :: calc_flag = .false.
        integer, allocatable, dimension(:)   :: vcorrcount
        real(8), allocatable, dimension(:)   :: vcorr
        real(8), allocatable, dimension(:,:) :: histdata
    end type

    type(msd_t) :: msd1
    type(vcorr_t) :: vcorr

contains

    subroutine init_msd( tcon, tmsd, tndt, tnhistmax )
        implicit none

        ! para list
        type(con_t), intent(in)    :: tcon
        type(msd_t), intent(inout) :: tmsd
        integer, intent(in)        :: tndt
        integer, intent(in)        :: tnhistmax

        associate(                    &
            natom    => tmsd%natom,   &
            ndt      => tmsd%ndt,     &
            nhistmax => tmsd%nhistmax &
            )

            natom = tcon%natom
            ndt   = tndt
            nhistmax = tnhistmax

            allocate(                         &
                tmsd%msd(nhistmax),           &
                tmsd%msdcount(nhistmax),      &
                tmsd%alpha2(nhistmax),        &
                tmsd%histdata(natom,nhistmax) &
                )

            tmsd%msd = 0.d0
            tmsd%msdcount = 0
            tmsd%alpha2 = 0.d0
            tmsd%histdata = 0.d0


        end associate
    end subroutine

    subroutine calc_msd( tcon, tmsd )
        implicit none

        ! para list
        type(con_t), intent(in)    :: tcon
        type(msd_t), intent(inout) :: tmsd

        ! local
        integer :: idx
        integer :: i

        associate(                       &
            ra        => tcon%ra,        &
            cumstep   => tmsd%cumstep,   &
            histidx   => tmsd%histidx,   &
            calc_flag => tmsd%calc_flag, &
            ndt       => tmsd%ndt,       &
            nhistmax  => tmsd%nhistmax,  &
            msd       => tmsd%msd,       &
            msdcount  => tmsd%msdcount,  &
            alpha2    => tmsd%alpha2,    &
            histdata  => tmsd%histdata   &
            )

            cumstep = cumstep + 1
            if ( cumstep >= ndt ) cumstep = 0

            if ( cumstep /= 0 ) return

            ! main v

            histidx = histidx + 1

            ! if histidx less then nhistmax, then just store configure, not calc msd
            if ( .not. calc_flag .and. histidx == nhistmax ) calc_flag = .true.

            ! if histidx greater then max dim of histdata, reset to one
            if ( histidx == nhistmax + 1 ) histidx = 1

            ! store config to "history data"
            histdata(:,histidx) = ra(1,:)

            if ( calc_flag ) then

                do i=1, nhistmax

                    ! if histidx > i; idx = histidx - i + 1
                    ! if histidx < i; idx = nhistmax + histidx - i + 1
                    !if ( i == histidx ) cycle
                    idx = histidx-i+1
                    if ( idx <= 0 ) idx = idx + nhistmax

                    msdcount(idx) = msdcount(idx) + 1
                    msd(idx)      = msd(idx)      + sum( ( histdata(:,i) - histdata(:,histidx) )**2 )
                    alpha2(idx)   = alpha2(idx)   + sum( ( histdata(:,i) - histdata(:,histidx) )**4 )

                end do

            end if

        end associate
    end subroutine

    subroutine endof_msd( tmsd )
        implicit none

        ! para list
        type(msd_t), intent(inout) :: tmsd

        associate(                    &
            natom    => tmsd%natom,   &
            msdcount => tmsd%msdcount &
            )

            tmsd%msd    = tmsd%msd    / ( msdcount * natom )
            tmsd%alpha2 = tmsd%alpha2 / ( msdcount * natom )
            tmsd%alpha2 = tmsd%alpha2 / tmsd%msd**2 / 3.d0 -1.d0

            ! todo fkt, ...

        end associate
    end subroutine

    subroutine init_vcorr( tcon, tvcorr, tndt, tnhistmax )
        implicit none

        ! para list
        type(con_t),   intent(in)    :: tcon
        type(vcorr_t), intent(inout) :: tvcorr
        integer,       intent(in)    :: tndt
        integer,       intent(in)    :: tnhistmax

        associate(                    &
            natom    => tvcorr%natom,   &
            ndt      => tvcorr%ndt,     &
            nhistmax => tvcorr%nhistmax &
            )

            natom = tcon%natom
            ndt   = tndt
            nhistmax = tnhistmax

            allocate(                           &
                tvcorr%vcorr(nhistmax),         &
                tvcorr%vcorrcount(nhistmax),    &
                tvcorr%histdata(natom,nhistmax) &
                )

            tvcorr%vcorr      = 0.d0
            tvcorr%vcorrcount = 0
            tvcorr%histdata   = 0.d0

        end associate
    end subroutine

    subroutine calc_vcorr( tcon, tvcorr )
        implicit none

        ! para list
        type(con_t),   intent(in)    :: tcon
        type(vcorr_t), intent(inout) :: tvcorr

        ! local
        integer :: idx
        integer :: i

        associate(                           &
            va         => tcon%va,           &
            cumstep    => tvcorr%cumstep,    &
            histidx    => tvcorr%histidx,    &
            calc_flag  => tvcorr%calc_flag,  &
            ndt        => tvcorr%ndt,        &
            nhistmax   => tvcorr%nhistmax,   &
            vcorr      => tvcorr%vcorr,      &
            vcorrcount => tvcorr%vcorrcount, &
            histdata   => tvcorr%histdata    &
            )

            cumstep = cumstep + 1
            if ( cumstep >= ndt ) cumstep = 0

            if ( cumstep /= 0 ) return

            ! main v

            histidx = histidx + 1

            ! if histidx less then nhistmax, then just store configure, not calc vcorr
            if ( .not. calc_flag .and. histidx == nhistmax ) calc_flag = .true.

            ! if histidx greater then max dim of histdata, reset to one
            if ( histidx == nhistmax + 1 ) histidx = 1

            ! store config to "history data"
            histdata(:,histidx) = va(1,:)

            if ( calc_flag ) then

                do i=1, nhistmax

                    ! if histidx > i; idx = histidx - i + 1
                    ! if histidx < i; idx = nhistmax + histidx - i + 1
                    !if ( i == histidx ) cycle
                    idx = histidx-i+1
                    if ( idx <= 0 ) idx = idx + nhistmax

                    vcorrcount(idx) = vcorrcount(idx) + 1
                    vcorr(idx)      = vcorr(idx)      + dot_product( histdata(:,i), histdata(:,histidx) )

                end do

            end if

        end associate
    end subroutine

    subroutine endof_vcorr( tvcorr )
        implicit none

        ! para list
        type(vcorr_t), intent(inout) :: tvcorr

        associate(                          &
            vcorrcount => tvcorr%vcorrcount &
            )

            tvcorr%vcorr = tvcorr%vcorr / vcorrcount
            tvcorr%vcorr = tvcorr%vcorr / tvcorr%vcorr(1)

        end associate
    end subroutine

end module
