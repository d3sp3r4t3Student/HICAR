module variable_interface
    use icar_constants,          only : kMAX_DIM_LENGTH, kMAX_STRING_LENGTH, kMAX_NAME_LENGTH, kINTEGER, kREAL, kDOUBLE
    use grid_interface,          only : grid_t
    use meta_data_interface,     only : meta_data_t
    use iso_fortran_env,         only : real64

    implicit none

    ! defines a variable type that can store data and attributes
    ! have to think about how to handle multiple variable types (int, 2d, etc)
    ! could add multiple "local" variables or create multiple variable types...
    type, extends(meta_data_t) :: variable_t
        real, pointer :: data_3d(:,:,:) => null()
        real, pointer :: data_2d(:,:)   => null()
        real(kind=real64), pointer :: data_2dd(:,:) => null()
        integer, pointer ::           data_2di(:,:) => null()
        real, pointer :: dqdt_3d(:,:,:) => null()   ! Note these have to be pointers so they get referenced when variable_t is passed around(?)
        real, pointer :: dqdt_2d(:,:)   => null()   ! Note these have to be pointers so they get referenced when variable_t is passed around(?)

        logical                         :: unlimited_dim = .False.
        logical                         :: three_d = .False.
        logical                         :: two_d = .False.
        logical                         :: force_boundaries = .True.
        logical                         :: computed = .False.
        character(len=kMAX_NAME_LENGTH) :: forcing_var = ""

        integer :: n_dimensions
        integer :: dtype = kREAL
        integer,                        allocatable :: dim_len(:)
        integer,                        allocatable :: global_dim_len(:)
        type(grid_t)                                :: grid
        character(len=kMAX_DIM_LENGTH), allocatable :: dimensions(:)

        ! note these are used for netcdf output
        integer, allocatable    :: dim_ids(:)
        integer                 :: var_id = -1
        integer :: xstag = 0
        integer :: ystag = 0
    contains
        procedure, public  :: init_grid
        procedure, public  :: init_dims
        generic,   public  :: initialize => init_grid
        generic,   public  :: initialize => init_dims

    ! inherited from meta_data
    !     procedure, public : add_attribute

    end type

    interface

        module subroutine init_grid(this, grid, forcing_var, force_boundaries, dtype)
            implicit none
            class(variable_t),  intent(inout) :: this
            type(grid_t),       intent(in)    :: grid
            character(len=*),   intent(in), optional :: forcing_var
            logical,            intent(in), optional :: force_boundaries
            integer,            intent(in), optional :: dtype

        end subroutine

        module subroutine init_dims(this, dims, forcing_var, force_boundaries)
            implicit none
            class(variable_t),  intent(inout) :: this
            integer,            intent(in)    :: dims(:)
            character(len=*),   intent(in), optional :: forcing_var
            logical,            intent(in), optional :: force_boundaries
        end subroutine

    end interface

end module
