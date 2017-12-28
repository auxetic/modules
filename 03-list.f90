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

    type tpvoro_one
        integer :: nbsum
        integer :: nblist(listmax)
        integer :: vid1(listmax), vid2(listmax)
    end type

    type tpvoro
        integer :: natom
        type(tpvoro_one), allocatable, dimension(:) :: list
        real(8), allocatable, dimension(:,:) :: center, vertex
        real(8) :: lx, ly, strain
    contains
        procedure :: init      => init_voro
        procedure :: decompose => calc_voro
    end type

    type(tpvoro) :: voro

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

        if ( allocated(tnb%list) ) then
            if ( size(tnb%list) /= tnatom ) then
                deallocate( tnb%list )
                allocate( tnb%list(tnatom) )
            end if
        else
            allocate( tnb%list(tnatom) )
        end if
    end subroutine

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
    end subroutine

    subroutine calc_z( tnb, tcon )
        implicit none

        ! var list
        type(tpcon),  intent(in)    :: tcon
        type(tplist), intent(inout) :: tnb

        ! local
        real(8) :: dra(free), rai(free), raj(free), ri, rj, rij2, dij
        integer :: cory, iround(free)
        integer :: i, j, k, jj

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
        integer :: i

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
        if ( maxdis > 0.25 * nlcut**2 ) flag = .true.
    end function

    ! ToDo
    ! subroutine check_rattler
    function calc_rattler( tcon, tnblist ) result(flag)
        implicit none

        type(tpcon),  intent(in)           :: tcon
        type(tplist), intent(in), optional :: tnblist
        integer, dimension(tcon%natom)     :: flag

        type(tplist) :: lclist
        integer      :: i

        if ( present(tnblist) ) then
            lclist = tnblist
        else
            call init_list( lclist, tcon )
            call make_list( lclist, tcon )
        end if

        call calc_z( lclist, tcon )

        flag = 0
        do i=1, lclist%natom
            if ( lclist%nbi(i) == 0 ) then
                flag(i) = 1
            elseif ( lclist%nbi(i) <= (free-1) ) then
                print*, "There exist unstable particle(s)"
                stop
            end if
        end do
    end function

    ! voro
    subroutine init_voro( this, tcon )
        implicit none

        class(tpvoro), intent(inout) :: this
        type(tpcon),  intent(in)     :: tcon

        this%natom  = tcon%natom
        this%center = tcon%ra
        this%lx     = tcon%lx
        this%ly     = tcon%ly
        this%strain = tcon%strain

        allocate( this%list(this%natom) )
    end subroutine

    subroutine calc_voro( this )
        implicit none

        class(tpvoro), intent(inout) :: this

        integer, parameter :: maxcan = 200
        integer, dimension(maxcan) :: verts
        real(8), dimension(maxcan) :: px, py, ps
        integer, dimension(maxcan) :: tag
        real(8), parameter :: rcut = 3.5

        real(8) :: rai(free), raij(free), rijsq
        integer :: i, j, k, cory, ncan

        associate( natom  => this%natom,  &
                   list   => this%list,   &
                   con    => this%center, &
                   lx     => this%lx,     &
                   ly     => this%ly,     &
                   strain => this%strain  &
                   )

        list(:)%nbsum = 0

        ! main loop
        do i=1, natom
            rai = con(:,i)
            k   = 0
            do j=1, natom
                if ( i == j ) cycle

                raij = con(:,j) - rai
                cory = nint( raij(2) / ly )
                raij(1) = raij(1) - cory * strain * ly
                raij(1) = raij(1) - anint( raij(1) / lx ) * lx
                raij(2) = raij(2) - anint( raij(2) / ly ) * ly

                rijsq = sum(raij**2)

                if ( rijsq < rcut**2 ) then
                    if ( k == maxcan ) then
                        print*, k
                        cycle
                    end if
                    k      = k + 1
                    px(k)  = raij(1)
                    py(k)  = raij(2)
                    ps(k)  = rijsq
                    tag(k) = j
                end if
            end do

            ncan = k

            call sortps( px, py, ps, tag, ncan )

            call calc_voro_single( px, py, ps, maxcan, ncan, verts )

            do k=1, ncan
                if ( verts(k) /= 0 ) then
                    list(i)%nbsum = list(i)%nbsum + 1
                    list(i)%nblist(list(i)%nbsum) = tag(k)
                end if
            end do

        end do

        end associate
    end subroutine

    subroutine sortps( px, py, ps, tag, ncan )
        use mo_math, only: swap
        implicit none

        integer :: ncan
        real(8), dimension(ncan) :: px, py, ps
        integer, dimension(ncan) :: tag

        integer :: i, imin

        do i=1, ncan-1

            imin = minloc( ps(i:ncan), 1 ) + i - 1
            if ( i == imin ) cycle

            call swap( px(i),  px(imin)  )
            call swap( py(i),  py(imin)  )
            call swap( ps(i),  ps(imin)  )
            call swap( tag(i), tag(imin) )

        end do
    end subroutine

    subroutine calc_voro_single( px, py, ps, maxcan, ncan, verts )
        implicit none

        integer :: maxcan, ncan, nv, ne
        real(8), dimension(maxcan) :: px, py, ps
        integer, dimension(maxcan) :: verts
        integer, parameter :: maxv = 200
        real(8), dimension(maxv) :: vx, vy
        integer, dimension(maxv) :: iv, jv

        logical :: flag
        real(8) :: ai,bi,ci, aj,bj,cj, det, detinv
        real(8) :: vxij, vyij
        real(8), parameter :: tol=1.d-8

        integer :: i, j, l, v

        !-- check
        if ( ncan < 3 ) then
            write(*,'('' less than 3 points given to work '',i5)') ncan
            stop
        end if

        v = 0
        do i=1, ncan-1

            ai =  px(i)
            bi =  py(i)
            ci = -ps(i)

            do j=i+1, ncan

                aj =  px(j)
                bj =  py(j)
                cj = -ps(j)

                det = ai*bj - aj*bi

                if( abs(det) >= tol ) then

                    detinv = 1.d0 / det
                    vxij = ( bi * cj - bj * ci ) * detinv
                    vyij = ( aj * ci - ai * cj ) * detinv
                    flag = .true.
                    l  = 1
                    do while( flag .and. l < ncan )
                        if ( l /= i .and. l /= j ) then
                            flag = ( ( px(l) * vxij + py(l) * vyij ) .le. ps(l) )
                        end if
                        l = l + 1
                    end do

                    if ( flag ) then
                        v = v + 1
                        if ( v > maxv ) stop 'too many vertices'
                        iv(v)  = i
                        jv(v)  = j
                        vx(v) = 0.5 * vxij
                        vy(v) = 0.5 * vyij
                    end if
                end if
            end do
        end do

        nv = v                      ! total vertex found
        if ( nv < 3 ) then
            write(*,'('' less than 3 vertices found in work '',i5)') nv
            stop
        end if

        verts(1:ncan) = 0

        do i = 1, nv
            verts(iv(i)) = verts(iv(i)) + 1
            verts(jv(i)) = verts(jv(i)) + 1
        end do

        flag = .true.
        ne   = 0

        do i = 1, ncan
            if ( verts(i) > 0 ) then
                ne = ne + 1
                if ( verts(i) /= 2 ) then    ! it is supposed that every neighbor contribute 2 vertex
                    flag = .false.
                endif
            endif
        end do

        if ( flag .eqv. .false. ) then
            write (*,'('' **** vertex error: degeneracy ? **** '')')
        endif

        if ( ne /= nv ) then
            write(*,'('' **** edge   error: degeneracy ? ****  '')')
        endif
    end subroutine

end module
