module mo_network
    use mo_syst
    use mo_config
    implicit none

    type tpspring
        integer :: i, j
        integer :: cory, iround(free)
        real(8) :: l0, lvec(free)
        real(8) :: ks
    end type

    type tpnetwork
        integer :: nsps, max_of_springs
        integer :: natom
        type(tpspring), allocatable, dimension(:) :: sps
    end type

    type(tpnetwork) :: net

contains

    subroutine init_network( tnetwork, tcon )
        implicit none

        ! para list
        type(tpcon),     intent(in)    :: tcon
        type(tpnetwork), intent(inout) :: tnetwork

        associate(                                    &
            natom          => tnetwork.natom,         &
            max_of_springs => tnetwork.max_of_springs &
            )

            natom = tcon.natom
            max_of_springs = natom * 10
            allocate( tnetwork.sps(max_of_springs) )

        end associate
    end subroutine init_network

    subroutine make_network( tnetwork, tcon )
        implicit none

        ! para list
        type(tpnetwork), intent(inout) :: tnetwork
        type(tpcon),     intent(in)    :: tcon

        ! local
        integer :: i, j, k
        real(8) :: dra(free), rij2, dij
        integer :: cory, iround(free)

        associate(                   &
            natom  => tcon.natom,    &
            ra     => tcon.ra,       &
            r      => tcon.r,        &
            la     => tcon.la,       &
            lainv  => tcon.lainv,    &
            strain => tcon.strain,   &
            nsps   => tnetwork.nsps, &
            sps    => tnetwork.sps   &
            )

            nsps = 0
            do i=1, natom
                do j=i+1, natom

                    dra = ra(:,j) - ra(:,i)

                    cory = nint( dra(free) * lainv(free) )
                    dra(1) = dra(1) - cory * strain * la(free)

                    do k=1, free-1
                        iround(k) = nint( dra(k) * lainv(k) )
                    end do
                    iround(free) = cory

                    do k=1, free
                        dra(k) = dra(k) - iround(k) * la(k)
                    end do

                    rij2 = sum( dra**2 )
                    dij = r(i) + r(j)

                    if ( rij2 > dij**2 ) cycle

                    nsps             = nsps + 1
                    sps(nsps).i      = i
                    sps(nsps).j      = j
                    sps(nsps).cory   = cory
                    sps(nsps).iround = iround
                    sps(nsps).lvec   = dra
                    sps(nsps).l0     = sqrt(rij2)
                    sps(nsps).ks     = 1.d0

                end do
            end do

        end associate

        ! reallocate array of sps
        tnetwork.sps = tnetwork.sps(1:tnetwork.nsps)

    end subroutine make_network

    subroutine remake_network( tnetwork, tcon )
        implicit none

        ! para list
        type(tpnetwork), intent(inout) :: tnetwork
        type(tpcon),     intent(in)    :: tcon

        ! local
        integer :: i, j, k, ii
        real(8) :: dra(free), rij2, dij
        integer :: cory, iround(free)

        associate(                   &
            natom  => tcon.natom,    &
            ra     => tcon.ra,       &
            r      => tcon.r,        &
            la     => tcon.la,       &
            lainv  => tcon.lainv,    &
            strain => tcon.strain,   &
            nsps   => tnetwork.nsps, &
            sps    => tnetwork.sps   &
            )

            do ii=1, nsps

                i = sps(ii).i
                j = sps(ii).j

                dra = ra(:,j) - ra(:,i)

                cory = nint( dra(free) * lainv(free) )
                dra(1) = dra(1) - cory * strain * la(free)

                do k=1, free-1
                    iround(k) = nint( dra(k) * lainv(k) )
                end do
                iround(free) = cory

                do k=1, free
                    dra(k) = dra(k) - iround(k) * la(k)
                end do

                rij2 = sum( dra**2 )

                sps(ii).cory   = cory
                sps(ii).iround = iround
                sps(ii).lvec   = dra
                sps(ii).l0     = sqrt(rij2)

            end do

        end associate

    end subroutine remake_network

    subroutine change_k_spring( tnetwork, tcon, opktan, oph )
        use ifport
        implicit none

        ! para list
        type(tpnetwork), intent(inout) :: tnetwork
        type(tpcon),     intent(in)    :: tcon
        real(8),         optional      :: opktan
        real(8),         optional      :: oph

        ! local
        integer :: i, itemp
        real(8) :: ktan, h


        ! 1. ks = ls
!       do i=1, tnetwork.nsps
!           tnetwork.sps(i).ks = calc_len( tcon, i, tnetwork )
!       end do

!       do i=1, tnetwork.nsps
!           itemp = floor( rand(0) * tnetwork.nsps ) + 1
!           tnetwork.sps(i).ks = calc_len( tcon, itemp, tnetwork )
!       end do

        ktan = 1.d0
        if ( present( opktan ) ) ktan = opktan
!        do i=1, tnetwork.nsps
!            tnetwork.sps(i).ks = 1.0 + tanh( ktan * ( calc_len( tcon, i, tnetwork ) - 1.d0 ) )
!!           tnetwork.sps(i).ks = 1.0 + tanh( ktan * ( tnetwork.sps(i).l0 - 1.d0 ) )
!        end do
        do i=1, tnetwork.nsps
            itemp = floor( rand(0)*tnetwork.nsps ) + 1
            tnetwork.sps(i).ks = 1.0 + tanh( ktan * ( calc_len( tcon, itemp, tnetwork ) - 1.d0 ) )
!           tnetwork.sps(i).ks = 1.0 + tanh( ktan * ( tnetwork.sps(i).l0 - 1.d0 ) )
        end do

!       h = 0.1d0
!       if ( present( oph ) ) h = oph
!       do i=1, tnetwork.nsps
!           if (  calc_len( tcon, i, tnetwork ) > 1.d0 ) then
!               tnetwork.sps(i).ks = 2.d0 - h
!           else
!               tnetwork.sps(i).ks = h
!           end if
!       end do

    end subroutine change_k_spring

    subroutine change_k_spring_2( tnetwork, tcon, opktan, oph )
        use ifport
        implicit none

        ! para list
        type(tpnetwork), intent(inout) :: tnetwork
        type(tpcon),     intent(in)    :: tcon
        real(8),         optional      :: opktan
        real(8),         optional      :: oph

        ! local
        integer :: i
        real(8) :: h, ktan

        h = 0.1
        if ( present( oph ) ) h = oph

!       do i=1, tnetwork.nsps
!           if ( rand(0) > 0.5 ) then
!               tnetwork.sps(i).ks = 2.d0 - h
!           else
!               tnetwork.sps(i).ks = h
!           end if
!       end do

!       do i=1, tnetwork.nsps
!           if ( calc_len( tcon, i, tnetwork ) > 1.d0 ) then
!               tnetwork.sps(i).ks = 2.d0 - h
!           else
!               tnetwork.sps(i).ks = h
!           end if
!       end do

        ktan = 1.d0
        if ( present( opktan ) ) ktan = opktan
        do i=1, tnetwork.nsps
            tnetwork.sps(i).ks = 1.d0 + tanh( ktan * ( calc_len( tcon, i, tnetwork ) - 1.d0 ) )
        end do


    end subroutine change_k_spring_2

    function calc_len( tcon, tibond, tnetwork  ) result(tl)
        implicit none

        ! para list
        type(tpcon), intent(in)     :: tcon
        integer, intent(in)         :: tibond
        type(tpnetwork), intent(in) :: tnetwork
        real(8) :: tl

        ! local
        integer :: i, j, k, ii
        real(8) :: dra(free), rij2, dij
        integer :: cory, iround(free)

        associate(                   &
            natom  => tcon.natom,    &
            ra     => tcon.ra,       &
            r      => tcon.r,        &
            la     => tcon.la,       &
            lainv  => tcon.lainv,    &
            strain => tcon.strain,   &
            nsps   => tnetwork.nsps, &
            sps    => tnetwork.sps   &
            )


            i = sps(tibond).i
            j = sps(tibond).j

            dra = ra(:,j) - ra(:,i)

            cory = nint( dra(free) * lainv(free) )
            dra(1) = dra(1) - cory * strain * la(free)

            do k=1, free-1
                iround(k) = nint( dra(k) * lainv(k) )
            end do
            iround(free) = cory

            do k=1, free
                dra(k) = dra(k) - iround(k) * la(k)
            end do

            tl = norm2( dra )

        end associate

    end function calc_len


end module
