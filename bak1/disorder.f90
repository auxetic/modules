module mo_disorder
    use mo_syst
    use mo_config
    implicit none
        
    real(8) :: eta, eta_spring
    real(8), allocatable, dimension(:,:) :: radisorder

contains

    subroutine init_disorder( tradisorder, tcon, tnp )
        use ifport
        implicit none
           
        type(tpcon) :: tcon
        real(8), allocatable, dimension(:,:) :: tradisorder
        integer :: tnp

        real(8) :: temp
        integer :: i

        temp = rand( tnp )
        tnp = 0

        associate( natom => tcon.natom )

            allocate( radisorder(free,natom) )
            do i=1, natom
                temp = rand(0) * 2.d0 * pi
                radisorder(1,i) = cos(temp)
                radisorder(2,i) = sin(temp)
    !           radisorder(1,i) = 2.0 * rand(0) - 1.0
    !           radisorder(2,i) = 2.0 * rand(0) - 1.0
            end do

        end associate

    end subroutine init_disorder

    subroutine make_disorder( tcon, tradisorder, teta )
        implicit none
    
        type(tpcon) :: tcon
        real(8), allocatable, dimension(:,:) :: tradisorder
        real(8) :: teta

        integer :: i

        associate( ra => tcon.ra )

            ra = ra + teta * tradisorder

        end associate

    end subroutine

end module
