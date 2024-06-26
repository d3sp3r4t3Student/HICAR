module grid_interface

    use icar_constants 
    use mod_wrf_constants, only : epsilon
    use mpi_f08
    implicit none

    private
    public :: grid_t

    type grid_t
        integer :: yimg,    ximg
        integer :: yimages, ximages
        integer :: ims, ime
        integer :: jms, jme
        integer :: kms, kme
        integer :: ns_halo_nx, ew_halo_ny, halo_nz, halo_size
        integer :: nx_e, ny_e
        integer :: nx_global, ny_global
        integer :: nx, ny, nz

        integer ::  ids,ide, jds,jde, kds,kde, & ! for the entire model domain    (d)
                    its,ite, jts,jte, kts,kte    ! for the data tile to process   (t)

        logical :: is2d, is3d
        character(len=kMAX_DIM_LENGTH), allocatable :: dimensions(:)

        type(MPI_Datatype) :: NS_halo
        type(MPI_Datatype) :: NS_win_halo

        type(MPI_Datatype) :: EW_halo
        type(MPI_Datatype) :: EW_win_halo

        type(MPI_Datatype) :: corner_halo
        type(MPI_Datatype) :: corner_NS_win_halo
        type(MPI_Datatype) :: corner_EW_win_halo
    contains
        procedure :: get_dims
        procedure :: domain_decomposition
        procedure :: set_grid_dimensions

    end type

interface
    module function get_dims(this) result(dims)
        implicit none
        class(grid_t), intent(in) :: this
        integer, allocatable :: dims(:)
    end function

    module subroutine domain_decomposition(this, nx, ny, nimages, image, ratio)
        class(grid_t),  intent(inout) :: this
        integer,        intent(in)    :: nx, ny, nimages
        integer,        intent(in)    :: image
        real,           intent(in), optional :: ratio
    end subroutine

    module subroutine set_grid_dimensions(this, nx, ny, nz, image, comms, global_nz, adv_order, nx_extra, ny_extra)
        implicit none
        class(grid_t),   intent(inout) :: this
        integer,         intent(in)    :: nx, ny, nz, image
        type(MPI_Comm), optional,  intent(in)    :: comms
        integer,         intent(in), optional :: global_nz, adv_order, nx_extra, ny_extra

    end subroutine
end interface
end module
