module mo_mode
    use mo_syst
    use mo_config
    implicit none

    type tpmatrix
        real(8), allocatable, dimension(:,:) :: dymatrix, dymatrix0
        real(8), allocatable, dimension(:)   :: egdymatrix, pw
        real(8), allocatable, dimension(:,:,:) :: trimatrix
        real(8), allocatable, dimension(:)   :: modez
        integer :: ndim, mdim, natom
        logical :: boxflag   = .false.      ! can box length change
        logical :: xyflag    = .false.      ! can box length change seperately
        logical :: shearflag = .false.      ! can box be changed by shear
    contains
        procedure :: solve => solve_mode
        procedure :: calc_pw => calc_pw
    end type

    type(tpmatrix) :: mode

contains

    subroutine init_mode( tmode, tcon, opboxflag, opxyflag, opshearflag )
        implicit none

        ! para list
        type(tpmatrix), intent(inout) :: tmode
        type(tpcon), intent(in) ::  tcon
        logical, intent(in), optional :: opboxflag, opxyflag, opshearflag

        ! local

        associate(                       &
            natom     => tmode.natom,    &
            mdim      => tmode.mdim,     &
            ndim      => tmode.ndim,     &
            boxflag   => tmode.boxflag,  &
            xyflag    => tmode.xyflag,   &
            shearflag => tmode.shearflag &
            )

            if ( present(opboxflag) .and. present(opxyflag) ) stop "boxflag and xyflag can not exist simultaneously"
            if ( present(opboxflag) )   boxflag   = .true.
            if ( present(opxyflag) )    xyflag    = .true.
            if ( present(opshearflag) ) shearflag = .true.

            natom = tcon.natom
            ndim = free * tcon.natom
            mdim = ndim

            if ( boxflag   ) mdim = mdim + 1
            if ( xyflag    ) mdim = mdim + free
            if ( shearflag ) mdim = mdim + 1

            allocate( tmode.dymatrix(mdim,mdim), tmode.egdymatrix(mdim) )
           !allocate( tmode.dymatrix0(mdim,mdim) )
            allocate( tmode.trimatrix(mdim,mdim,mdim) )
            allocate( tmode.modez(mdim) )

        end associate

    end subroutine

    subroutine make_dymatrix( tmode, tcon )
        implicit none

        ! para list
        type(tpmatrix), intent(inout) :: tmode
        type(tpcon), intent(in) :: tcon

        ! local
        integer :: i, j
        real(8) :: dra(free), rij, rij2, rij3, dij, xj, yj
        real(8) :: vr, vrr, rx, ry, rxx, ryy, rxy
        real(8) :: mij(free,free)   ! \p^2v / ( \pxi * \pxj ) ....

        associate(                      &
            dymatrix => tmode.dymatrix, &
            radius   => tcon.r,         &
            natom    => tcon.natom      &
            )

            dymatrix = 0.d0
        
            do i=1, natom-1
                do j=i+1, natom

                    dra = calc_dra( tcon, i, j )
                    rij2 = sum(dra**2)

                    dij = radius(i) + radius(j)

                    if ( rij2 > dij**2 ) cycle

                    rij  = sqrt(rij2)
                    rij3 = rij2*rij

                    vr = - ( 1.d0 - rij/dij )**(alpha-1) / dij
                    vrr = ( alpha-1) * ( 1.d0 - rij/dij )**(alpha-2) / dij**2

                    xj = dra(1); yj = dra(2)

                    rx = xj / rij
                    ry = yj / rij

                    rxx = 1.d0/rij - xj**2/rij3
                    ryy = 1.d0/rij - yj**2/rij3
                    rxy =          - xj*yj/rij3

                    mij(1,1) = - ( vrr*rx**2 + vr*rxx )
                    mij(2,2) = - ( vrr*ry**2 + vr*ryy )
                    mij(1,2) = - ( vrr*rx*ry + vr*rxy )
                    mij(2,1) = mij(1,2)

                    dymatrix(free*(i-1)+1:free*i, free*(i-1)+1:free*i ) = dymatrix(free*(i-1)+1:free*i, free*(i-1)+1:free*i ) - mij
                    dymatrix(free*(j-1)+1:free*j, free*(j-1)+1:free*j ) = dymatrix(free*(j-1)+1:free*j, free*(j-1)+1:free*j ) - mij
                    dymatrix(free*(i-1)+1:free*i, free*(j-1)+1:free*j ) = mij
                    dymatrix(free*(j-1)+1:free*j, free*(i-1)+1:free*i ) = mij

                end do
            end do

        end associate

    end subroutine

    subroutine make_trimatrix( tmode, tcon )
        implicit none

        ! para list
        type(tpmatrix), intent(inout) :: tmode
        type(tpcon), intent(in) :: tcon

        ! local
        integer :: i, j
        real(8) :: dra(free), rij, rij2, rij3, rij5, dij, xi, xj, yi, yj
        real(8) :: vr, vrr, vrrr
        real(8) :: rx, ry 
        real(8) :: rxx, ryy, rxy
        real(8) :: rxxx, rxxy, rxyy, ryyy
        real(8) :: mijk(free,free,free)   ! \p^2v / ( \pxi * \pxj ) ....

        associate(                        &
            trimatrix => tmode.trimatrix, &
            radius    => tcon.r,          &
            natom     => tcon.natom       &
            )

            trimatrix = 0.d0
        
            do i=1, natom-1
                do j=i+1, natom

                    dra = calc_dra( tcon, i, j )
                    rij2 = sum(dra**2)

                    dij = radius(i) + radius(j)

                    if ( rij2 > dij**2 ) cycle

                    rij  = sqrt(rij2)
                    rij3 = rij2*rij
                    rij5 = rij3*rij2

                    !vr = - ( 1.d0 - rij/dij )**(alpha-1) / dij
                    !vrr = ( alpha-1) * ( 1.d0 - rij/dij )**(alpha-2) / dij**2
                    ! alpha = 2
                    vr   = - ( 1.d0 - rij/dij ) / dij
                    vrr  = 1.d0 / dij**2
                    !vrrr = 0.d0

                    xj = dra(1); yj = dra(2)
                    xi =-dra(1); yi =-dra(2)

                    rx = xj / rij
                    ry = yj / rij

                    rxx = 1.d0/rij - rx**2/rij3
                    ryy = 1.d0/rij - ry**2/rij3
                    rxy =          - rx*ry/rij3

                    ! \pr / ( \pxj \pxj \pxj )
                    rxxx = 3.d0*xj/rij3 + xj**3/rij5
                    ryyy = 3.d0*yj/rij3 + yj**3/rij5
                    rxxy = -yj/rij3 + 3.d0*xj**2*yj/rij5
                    rxyy = -xj/rij3 + 3.d0*xj*yj**2/rij5

                    ! \pv / ( \pxi \pxj \pxj )
                    mijk(1,1,1) = -( vrr*3.d0*rxx*rx - vr*rxxx )
                    mijk(1,1,2) = -( vrr*(2.d0*rxx*ry+rx*rxy) - vr*rxxy )
                    mijk(1,2,1) = mijk(1,1,2)
                    mijk(1,2,2) = -( vrr*(2.d0*ryy*rx+ry*rxy) - vr*rxyy )
                    mijk(2,1,1) = -mijk(1,1,2)
                    mijk(2,1,2) = -mijk(1,2,2)
                    mijk(2,2,1) = -mijk(1,2,2)
                    mijk(2,2,2) = -( vrr*3.d0*ryy*ry - vr*ryyy )

                    ! i
                    trimatrix(free*(i-1)+1:free*i, free*(j-1)+1:free*j, free*(j-1)+1:free*j ) =  mijk
                    trimatrix(free*(i-1)+1:free*i, free*(j-1)+1:free*j, free*(i-1)+1:free*i ) = -mijk
                    trimatrix(free*(i-1)+1:free*i, free*(i-1)+1:free*i, free*(j-1)+1:free*j ) = -mijk
                    ! j
                    trimatrix(free*(j-1)+1:free*j, free*(j-1)+1:free*j, free*(i-1)+1:free*i ) =  mijk
                    trimatrix(free*(j-1)+1:free*j, free*(i-1)+1:free*i, free*(i-1)+1:free*i ) = -mijk
                    trimatrix(free*(j-1)+1:free*j, free*(i-1)+1:free*i, free*(j-1)+1:free*j ) =  mijk
                    ! iii, jjj
                    trimatrix(free*(i-1)+1:free*i, free*(i-1)+1:free*i, free*(i-1)+1:free*i ) = &
                        trimatrix(free*(i-1)+1:free*i, free*(i-1)+1:free*i, free*(i-1)+1:free*i ) + mijk
                    trimatrix(free*(j-1)+1:free*j, free*(j-1)+1:free*j, free*(j-1)+1:free*j ) = &
                        trimatrix(free*(j-1)+1:free*j, free*(j-1)+1:free*j, free*(j-1)+1:free*j ) - mijk

                end do
            end do

        end associate
        
        
    end subroutine

    subroutine clac_modez( tmode, tb )
        implicit none

        ! para list
        type(tpmatrix), intent(inout) :: tmode
        real(8), intent(out) :: tb

        ! local
        integer :: i, j, k
        real(8) :: kappa, tau
        real(8) :: pkpz(tmode.mdim), ptpz(tmode.mdim), pbpz(tmode.mdim)

        associate( &
            mdim => tmode.mdim, &
            dymatrix => tmode.dymatrix0, &
            trimatrix => tmode.trimatrix, &
            modez => tmode.modez &
            )
        
            kappa = 0.d0
            do i=1, mdim
                do j=1, mdim
                    kappa = kappa + dymatrix(j,i) * modez(i) * modez(j)
                end do
            end do

            tau = 0.d0
            do i=1, mdim
                do j=1, mdim
                    do k=1, mdim
                        tau = tau + trimatrix(k,j,i) * modez(i) * modez(j) * modez(k)
                    end do
                end do
            end do

            pkpz = 0.d0
            do i=1, mdim
                do j=1, mdim
                    pkpz(i) = pkpz(i) + dymatrix(j,i) * modez(j)
                end do
            end do

            ptpz = 0.d0
            do i=1, mdim
                do j=1, mdim
                    do k=1, mdim
                        ptpz(i) = ptpz(i) + trimatrix(k,j,i) * modez(k) * modez(j)
                    end do
                end do
            end do

            pbpz = 0.5d0*kappa**2/tau**2 * pkpz - 4.d0*kappa**3/3.d0/tau**3 * ptpz

            modez = modez - pbpz*1.d-2

            modez = modez / norm2(modez)

            tb = 2.d0*kappa**3/3.d0/tau**2
            print*, norm2(pbpz), tb, dot_product(pbpz,modez)/norm2(pbpz)/norm2(modez)

        end associate

    end subroutine

    subroutine solve_mode( this, oprange )
        implicit none

        class(tpmatrix) :: this
        integer, optional :: oprange

        ! local
        integer :: rangevar

        rangevar = 0
        if ( present( oprange ) ) rangevar = oprange

        associate(                        &
            mdim       => this.mdim,      &
            dymatrix   => this.dymatrix,  &
            egdymatrix => this.egdymatrix &
            )

            call solve_matrix( dymatrix, mdim, egdymatrix, rangevar )

        end associate

    end subroutine

    subroutine calc_pw( this )
        implicit none

        class(tpmatrix) :: this

        ! local
        integer :: i, j 
        real(8) :: sum4

        if ( .not. allocated(this.pw) ) allocate( this.pw(this.mdim) )

        associate(                &
            mdim   => this.mdim,  &
            natom  => this.natom, &
            pw     => this.pw,    &
            dymatrix => this.dymatrix &
            )
            
            do i=1, mdim
                sum4 = 0.d0
                do j=1, natom
                    sum4 = sum4 + sum( dymatrix(free*(j-1)+1:free*j,i)**2 )**2
                end do
                pw(i) = 1.d0 / dble(natom) / sum4
            end do

        end associate
        
    end subroutine

    subroutine solve_matrix(a,order,b, rangevar)
        implicit none
        
        !-- variables input and output    
        integer :: order
        real(8),dimension(1:order,1:order) :: a
        real(8),dimension(1:order) :: b
        integer :: rangevar

        !-- local variables    
        character :: jobz = 'V'
        character :: range = 'A'
        character :: uplo = 'U'
        
        integer :: n
    !   real, dimension  :: a
        integer :: lda
        integer :: vl, vu, il, iu   ! will not referenced when range=A 
        real(8) :: abstol           ! important variable
        integer :: m
    !   real, dimension :: w        ! use b above ! output eigenvalues
        real(8), allocatable, dimension(:,:) :: z
        integer :: ldz
        integer, allocatable, dimension(:) :: isuppz
        real(8), allocatable, dimension(:) :: work
        integer :: lwork
        integer, allocatable, dimension(:) :: iwork
        integer :: liwork
        integer :: info
       
        n = order
        lda = order
        ldz= order

        if( rangevar == 0 ) then
            range = 'A'
        elseif( rangevar > 0 .and. rangevar <= order ) then
            range = 'I'
            il = 1
            iu = rangevar
        else
            print*, "error rangevar set"
            stop
        end if
      
        !-- 
        abstol = -1 
    !   abstol = dlamch('S')        ! high precision
       
    !   allocate(a(order,order)); allocate(b(order))
        allocate(z(order,order))
        allocate(isuppz(2*order))
        

        !- query the optimal workspace    
        allocate(work(1000000)); allocate(iwork(1000000))
        lwork  = -1
        liwork = -1
        
        call dsyevr(jobz,range,uplo,n,a,lda,vl,vu,il,iu,abstol,m,b,z,       &
                   & ldz, isuppz, work, lwork, iwork, liwork, info)

        lwork  = int(work(1))
        liwork = int(iwork(1))
        
        deallocate(work); deallocate(iwork)
        
        allocate(work(lwork)); allocate(iwork(liwork))
        
        !vvv- main

        call dsyevr(jobz,range,uplo,n,a,lda,vl,vu,il,iu,abstol,m,b,z,       &
                   & ldz, isuppz, work, lwork, iwork, liwork, info)

        if(info .ne. 0) then
            write(*,*) '** Find error in subroutine diago'
            stop
        end if
        !^^^- done    
        
        deallocate(work); deallocate(iwork)
        
        !- output results via a    
        a = z
        
    !    deallocate(a); deallocate(b)
        deallocate(z); deallocate(isuppz)
    end subroutine

end module
