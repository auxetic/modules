module mo_list
    use mo_syst
    use mo_config
    implicit none

    integer, private, parameter :: listmax = 16
    real(8), private, parameter :: nlcut = 0.35d0

    type tplistone
        integer    :: nbsum
        integer    :: nblist(listmax)
        integer(1) :: iround(free,listmax)
        integer(1) :: cory(listmax)
        real(4)    :: con0(free)
    end type

    type tplist
        integer                                    :: natom
        type(tplistone), allocatable, dimension(:) :: list
        integer, allocatable, dimension(:)         :: nbi
        integer, allocatable, dimension(:)         :: rattlerflag
    end type

    type(tplist) :: nb

contains

    subroutine init_list(tnb, tcon)
        implicit none

        ! var list
        type(tplist), intent(inout) :: tnb
        type(tpcon),  intent(in)    :: tcon

        ! local
        integer :: tnatom

        tnatom    = tcon%natom
        tnb%natom = tnatom

        allocate( tnb%list(tnatom) )

    end subroutine init_list

    subroutine make_list(tnb, tcon)
        implicit none

        ! var list
        type(tpcon),  intent(in)    :: tcon
        type(tplist), intent(inout) :: tnb

        ! local
        real(8) :: dra(free), rai(free), raj(free), ri, rj, rij2, dij
        integer :: cory, iround(free)
        integer :: i, j, k, itemp

        associate(                 &
            natom  => tcon%natom,  &
            ra     => tcon%ra,     &
            r      => tcon%r,      &
            la     => tcon%la,     &
            lainv  => tcon%lainv,  &
            strain => tcon%strain, &
            list   => tnb%list     &
            )

            ! set nbsum to zero
            list(:)%nbsum = 0

            do i=1, natom

                list(i)%con0 = ra(:,i)
                rai          = ra(:,i)
                ri           = r(i)

                do j=i+1, natom

                    raj = ra(:,j)
                    rj  = r(j)

                    dra = raj - rai

                    cory   = nint( dra(free) * lainv(free) )
                    dra(1) = dra(1) - strain * la(free) * cory

                    do k=1, free-1
                        iround(k) = nint( dra(k) * lainv(k) )
                    end do
                    iround(free) = cory

                    do k=1, free
                        dra(k) = dra(k) - iround(k) * la(k)
                    end do

                    rij2 = sum( dra**2 )
                    dij  = ri + rj

                    if ( rij2 > ( dij+nlcut )**2 ) cycle

                    if ( list(i)%nbsum < listmax ) then
                        itemp                   = list(i)%nbsum
                        itemp                   = itemp + 1
                        list(i)%nbsum           = itemp
                        list(i)%nblist(itemp)   = j
                        list(i)%iround(:,itemp) = iround
                        list(i)%cory(itemp)     = cory
                    end if

                end do
            end do

        end associate

    end subroutine make_list

    subroutine calc_z( tnb, tcon )
        implicit none

        ! var list
        type(tpcon),  intent(in)    :: tcon
        type(tplist), intent(inout) :: tnb

        ! local
        real(8) :: dra(free), rai(free), raj(free), ri, rj, rij2, dij
        integer :: cory, iround(free)
        integer :: i, j, k, jj, itemp

        if ( allocated( tnb%nbi ) .and. size(tnb%nbi) /= tcon%natom ) then
            deallocate( tnb%nbi )
        end if

        if ( .not. allocated( tnb%nbi ) ) then
            allocate( tnb%nbi(tcon%natom) )
        end if

        associate(                 &
            natom  => tcon%natom,  &
            ra     => tcon%ra,     &
            r      => tcon%r,      &
            la     => tcon%la,     &
            lainv  => tcon%lainv,  &
            strain => tcon%strain, &
            list   => tnb%list,    &
            nbi    => tnb%nbi      &
            )

            ! set nbsum to zero
            nbi = 0

            do i=1, natom

                rai          = ra(:,i)
                ri           = r(i)

                do jj=1, list(i)%nbsum

                    j = list(i)%nblist(jj)

                    raj = ra(:,j)
                    rj  = r(j)

                    dra = raj - rai

                    cory = list(i)%cory(jj)
                    iround = list(i)%iround(:,jj)

                    dra(1) = dra(1) - strain * la(free) * cory

                    do k=1, free
                        dra(k) = dra(k) - iround(k) * la(k)
                    end do

                    rij2 = sum( dra**2 )
                    dij  = ri + rj

                    if ( rij2 > ( dij )**2 ) cycle

                    nbi(i) = nbi(i) + 1
                    nbi(j) = nbi(j) + 1

                end do
            end do

        end associate

    end subroutine

    function check_list( tnb, tcon ) result(flag)
        implicit none

        ! var list
        type(tplist), intent(in) :: tnb
        type(tpcon),  intent(in) :: tcon
        logical                  :: flag

        ! local
        real(8) :: maxdis, dra(free), dr2
        integer :: i, j

        associate(               &
            natom => tcon%natom, &
            ra    => tcon%ra     &
            )

            maxdis = 0.d0
            do i=1, tcon%natom
                dra = tcon%ra(:,i) - tnb%list(i)%con0
                dr2 = sum( dra**2 )
                if ( maxdis < dr2 ) maxdis = dr2
            end do

        end associate

        flag = .false.
        if ( maxdis > nlcut**2 ) flag = .true.

    end function check_list

    ! ToDo
    ! subroutine check_rattler

end module
