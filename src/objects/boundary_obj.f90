!>------------------------------------------------------------
!!  Handles reading boundary conditions from the forcing file(s)
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!------------------------------------------------------------
submodule(boundary_interface) boundary_implementation

    use io_routines,            only : io_getdims, io_read, io_maxDims, io_variable_is_present
    use time_io,                only : read_times, find_timestep_in_file
    use co_util,                only : broadcast
    use string,                 only : str
    use mod_atm_utilities,      only : rh_to_mr

    implicit none
contains

    !>------------------------------------------------------------
    !! Set default component values
    !! Reads initial conditions from the forcing file into image 1
    !!
    !! Distributes initial conditions to all other images
    !!
    !!------------------------------------------------------------
    module subroutine init(this, options)
        class(boundary_t), intent(inout) :: this
        class(options_t),  intent(inout) :: options

        character(len=kMAX_NAME_LENGTH), allocatable :: vars_to_read(:)
        integer,                         allocatable :: var_dimensions(:)

        this%file_list = options%parameters%boundary_files

        ! the parameters option type can't contain allocatable arrays because it is a coarray
        ! so we need to allocate the vars_to_read and var_dimensions outside of the options type
        call setup_variable_lists(options%parameters%vars_to_read, options%parameters%dim_list, vars_to_read, var_dimensions)

        ! Read through forcing variable names stored in "options"
        ! needs to read each one to find the grid information for it
        ! then create grid and initialize a variable...
        ! also need to explicitly save lat and lon data
        if (this_image() == 1) then
            call this%init_local(this%file_list,     &
                                 vars_to_read, var_dimensions,          &
                                 options%parameters%start_time,         &
                                 options%parameters%latvar,             &
                                 options%parameters%lonvar,             &
                                 options%parameters%zvar,               &
                                 options%parameters%time_var)

        endif

        call this%distribute_initial_conditions()

    end subroutine

    !>------------------------------------------------------------
    !! Set default component values
    !! Reads initial conditions from the forcing file
    !!
    !!------------------------------------------------------------
    module subroutine init_local(this, file_list, var_list, dim_list, start_time, &
                                 lat_var, lon_var, z_var, time_var)
        class(boundary_t),               intent(inout)  :: this
        character(len=kMAX_NAME_LENGTH), intent(in)     :: file_list(:)
        character(len=kMAX_NAME_LENGTH), intent(in)     :: var_list (:)
        integer,                         intent(in)     :: dim_list (:)
        type(Time_type),                 intent(in)     :: start_time
        character(len=kMAX_NAME_LENGTH), intent(in)     :: lat_var
        character(len=kMAX_NAME_LENGTH), intent(in)     :: lon_var
        character(len=kMAX_NAME_LENGTH), intent(in)     :: z_var
        character(len=kMAX_NAME_LENGTH), intent(in)     :: time_var

        type(variable_t)  :: test_variable


        integer :: i

        ! figure out while file and timestep contains the requested start_time
        call set_curfile_curstep(this, start_time, file_list, time_var)

        !  read in latitude and longitude coordinate data
        call io_read(file_list(this%curfile), lat_var, this%lat, this%curstep)
        call io_read(file_list(this%curfile), lon_var, this%lon, this%curstep)

        ! read in the height coordinate of the input data
        call io_read(file_list(this%curfile), z_var,   this%z,   this%curstep)

        ! call assert(size(var_list) == size(dim_list), "list of variable dimensions must match list of variables")

        do i=1, size(var_list)
            print*, trim(file_list(this%curfile)), "  ", trim(var_list(i))

            call add_var_to_dict(this%variables, file_list(this%curfile), var_list(i), dim_list(i), this%curstep)

        end do

        call setup_boundary_geo(this)

    end subroutine

    subroutine setup_boundary_geo(this)
        implicit none
        type(boundary_t), intent(inout) :: this

        if (allocated(this%geo%lat)) deallocate(this%geo%lat)
        allocate( this%geo%lat, source=this%lat)

        if (allocated(this%geo%lon)) deallocate(this%geo%lon)
        allocate( this%geo%lon, source=this%lon)

        if (allocated(this%geo%z)) deallocate(this%geo%z)
        allocate( this%geo%z, source=this%z)

    end subroutine



    !>------------------------------------------------------------
    !! Reads and adds a variable into the variable dictionary
    !!
    !! Given a filename, varname, number of dimensions (2,3) and timestep
    !! Read the timestep of varname from filename and stores the result in a variable structure
    !!
    !! Variable is then added to a master variable dictionary
    !!
    !!------------------------------------------------------------
    subroutine add_var_to_dict(var_dict, file_name, var_name, ndims, timestep)
        implicit none
        type(var_dict_t), intent(inout) :: var_dict
        character(len=*), intent(in)    :: file_name
        character(len=*), intent(in)    :: var_name
        integer,          intent(in)    :: ndims
        integer,          intent(in)    :: timestep

        real, allocatable :: temp_2d_data(:,:)
        real, allocatable :: temp_3d_data(:,:,:)
        type(variable_t)  :: new_variable

        if (ndims==2) then
            call io_read(file_name, var_name, temp_2d_data, timestep)

            call new_variable%initialize( shape( temp_2d_data ) )
            new_variable%data_2d = temp_2d_data

            call var_dict%add_var(var_name, new_variable)

            ! do not deallocate data arrays because they are pointed to inside the var_dict now
            ! deallocate(new_variable%data_2d)

        elseif (ndims==3) then
            call io_read(file_name, var_name, temp_3d_data, timestep)

            call new_variable%initialize( shape( temp_3d_data ) )
            new_variable%data_3d = temp_3d_data

            call var_dict%add_var(var_name, new_variable)

            ! do not deallocate data arrays because they are pointed to inside the var_dict now
            ! deallocate(new_variable%data_3d)
        endif

    end subroutine

    !>------------------------------------------------------------
    !! Reads a new set of forcing data for the next time step
    !!
    !!------------------------------------------------------------
    module subroutine update_forcing(this, options)
        class(boundary_t), intent(inout) :: this
        class(options_t),  intent(inout) :: options

        real, allocatable :: data3d(:,:,:), data2d(:,:)
        type(variable_t)  :: var
        character(len=kMAX_NAME_LENGTH) :: name


        if (this_image()==1) then
            call update_forcing_step(this, this%file_list, options%parameters%time_var)
            associate(list => this%variables)

            call list%reset_iterator()
            do while (list%has_more_elements())
                ! get the next variable in the structure
                var = list%next(name)

                ! because the data arrays are pointers, this should update the data stored in this%variables
                if (var%three_d) then
                    call io_read(this%file_list(this%curfile), name, data3d, this%curstep)
                    var%data_3d(:,:,:) = data3d(:,:,:)

                else if (var%two_d) then
                    call io_read(this%file_list(this%curfile), name, data2d, this%curstep)
                    var%data_2d(:,:) = data2d(:,:)
                endif

            end do

            end associate
        endif

        call this%distribute_update()

    end subroutine

    !>------------------------------------------------------------
    !! Sends the udpated forcing data to all other images
    !!
    !!------------------------------------------------------------
    module subroutine distribute_update(this)
        class(boundary_t), intent(inout) :: this

        type(variable_t)  :: var

        associate(list => this%variables)

        call list%reset_iterator()
        do while (list%has_more_elements())
            ! get the next variable in the structure
            var = list%next()

            if (var%three_d) then
                call broadcast(var%data_3d, 1, 1, num_images(), create_co_array = .True.)

            else if (var%two_d) then
                call broadcast(var%data_2d, 1, 1, num_images(), create_co_array = .True.)
            endif

        enddo

        end associate

    end subroutine


    subroutine distribute_2d_array(arr, source, first, last)
        implicit none
        real, allocatable, intent(inout) :: arr(:,:)
        integer,           intent(in)    :: source, first, last

        integer :: dims(2)

        ! because this array is probably not allocated, and may not be the correct shape if it is, we have to
        ! broadcast the shape of the array first.
        if (this_image()==source) dims = shape(arr)
        call broadcast(dims, source, first, last, create_co_array=.True.)

        if (allocated(arr)) then
            if (any(dims /= shape(arr))) deallocate(arr)
        endif
        if (.not.allocated(arr)) allocate( arr( dims(1), dims(2) ) )
        call broadcast(arr, source, first, last, create_co_array=.True.)

    end subroutine

    subroutine distribute_3d_array(arr, source, first, last)
        implicit none
        real, allocatable, intent(inout) :: arr(:,:,:)
        integer,           intent(in)    :: source, first, last

        integer :: dims(3)

        ! because this array is probably not allocated, and may not be the correct shape if it is, we have to
        ! broadcast the shape of the array first.
        if (this_image()==source) dims = shape(arr)
        call broadcast(dims, source, first, last, create_co_array=.True.)

        if (allocated(arr)) then
            if (any(dims /= shape(arr))) deallocate(arr)
        endif
        if (.not.allocated(arr)) allocate( arr( dims(1), dims(2), dims(3) ) )
        call broadcast(arr, source, first, last, create_co_array=.True.)

    end subroutine

    !>------------------------------------------------------------
    !! Distribute the initial conditions from image 1 into all other images
    !!
    !!------------------------------------------------------------
    module subroutine distribute_initial_conditions(this)
      class(boundary_t), intent(inout) :: this

      integer                         :: number_of_variables
      type(variable_t)                :: temp_variable
      integer                         :: i
      character(len=kMAX_NAME_LENGTH) :: name

      ! needs to distribute lat, lon, time and all vars in var_dict
      ! broadcast(variable, from:image, to: first_image, last_image)
      call distribute_2d_array(this%lat, 1, 1, num_images())
      call distribute_2d_array(this%lon, 1, 1, num_images())
      call distribute_3d_array(this%z,   1, 1, num_images())

      if (this_image() == 1) then
          call this%variables%reset_iterator()
          number_of_variables = this%variables%n_vars
      endif

      call broadcast(number_of_variables, 1, 1, num_images(), create_co_array=.True.)

      do i=1, number_of_variables
          ! if we are the sending image,
          if (this_image() == 1) then
              temp_variable = this%variables%next(name)
          endif

          ! broadcasting is handled within the variable_t object (though it calls broadcast on its components)
          ! we can't make a coarray of a variable_t because it has a pointer to data in it.
          call temp_variable%broadcast(1, 1, num_images())
          call broadcast(name, 1, 1, num_images(), create_co_array=.True.)

          if (this_image() /= 1) then
              call this%variables%add_var(name, temp_variable, save_state=.True.)
          endif
      enddo

    end subroutine


    ! count the number of variables specified, then allocate and store those variables in a list just their size.
    subroutine setup_variable_lists(master_var_list, master_dim_list, vars_to_read, var_dimensions)
        implicit none
        character(len=kMAX_NAME_LENGTH), intent(in)                 :: master_var_list(:)
        integer,                         intent(in)                 :: master_dim_list(:)
        character(len=kMAX_NAME_LENGTH), intent(inout), allocatable :: vars_to_read(:)
        integer,                         intent(inout), allocatable :: var_dimensions(:)

        integer :: n_valid_vars
        integer :: i, curvar, err

        do i=1, size(master_var_list)
            if (trim(master_var_list(i)) /= '') then
                n_valid_vars = n_valid_vars + 1
            endif
        enddo

        allocate(vars_to_read(  n_valid_vars), stat=err)
        if (err /= 0) stop "vars_to_read: Allocation request denied"

        allocate(var_dimensions(  n_valid_vars), stat=err)
        if (err /= 0) stop "var_dimensions: Allocation request denied"

        curvar = 1
        do i=1, size(master_var_list)
            if (trim(master_var_list(i)) /= '') then
                vars_to_read(curvar) = master_var_list(i)
                var_dimensions(curvar) = master_dim_list(i)
                curvar = curvar + 1
            endif
        enddo
    end subroutine

    !>------------------------------------------------------------
    !! Find the time step in the input forcing to start the model on
    !!
    !! The model start date (start_time) may not be the same as the first forcing
    !! date (initial_time).  Convert the difference between the two into forcing
    !! steps by dividing by the time delta between forcing steps (in_dt) after
    !! converting in_dt from seconds to days.
    !!
    !! @param  options  model options structure
    !! @retval step     integer number of steps into the forcing sequence
    !!
    !!------------------------------------------------------------
    ! function bc_find_step(options) result(step)
    !     implicit none
    !     type(options_t), intent(in) :: options
    !     integer :: step
    !     type(time_delta_t) :: dt
    !
    !     dt = (options%parameters%start_time - options%parameters%initial_time)
    !     step = dt%seconds() / options%parameters%in_dt + 1
    !
    !     if (options%parameters%debug) write(*,*) "bc_find_step: First forcing time step = ",trim(str(step))
    !
    ! end function bc_find_step


    function get_n_timesteps(filename, varname, var_space_dims) result(steps_in_file)
        implicit none
        character(len=*), intent(in) :: filename, varname
        integer,          intent(in), optional :: var_space_dims
        integer :: steps_in_file

        integer :: dims(io_maxDims)
        integer :: space_dims

        space_dims=3
        if (present(var_space_dims)) space_dims = var_space_dims

        call io_getdims(filename, varname, dims)

        if (dims(1) == space_dims) then
            steps_in_file = 1
        else
            steps_in_file = dims(dims(1)+1)
        endif

    end function get_n_timesteps


    subroutine set_curfile_curstep(bc, time, file_list, time_var)
        implicit none
        type(boundary_t),   intent(inout) :: bc
        type(Time_type),    intent(in) :: time
        character(len=*),   intent(in) :: file_list(:)
        character(len=*),   intent(in) :: time_var


        integer :: error

        ! these are module variables that should be correctly set when the subroutine returns
        bc%curfile = 1
        bc%curstep = 1
        error=1
        bc%curfile=0
        do while ( (error/=0) .and. (bc%curfile <= size(file_list)) )
            bc%curfile = bc%curfile + 1
            bc%curstep = find_timestep_in_file(file_list(bc%curfile), time_var, time, error=error)
        enddo

        if (error==1) then
            stop "Ran out of files to process while searching for matching time variable!"
        endif

    end subroutine set_curfile_curstep


    subroutine update_forcing_step(bc, file_list, time_var)
        implicit none
        type(boundary_t),   intent(inout) :: bc
        character(len=*),   intent(in) :: file_list(:)
        character(len=*),   intent(in) :: time_var


        integer :: error, steps_in_file

        error=1
        bc%curstep = bc%curstep + 1 ! this may be all we have to do most of the time

        ! check that we haven't stepped passed the end of the current file
        steps_in_file = get_n_timesteps(file_list(bc%curfile), time_var, 0)
        if (steps_in_file < bc%curstep) then
            ! if we have, use the next file
            bc%curfile = bc%curfile + 1
            ! and the first timestep in the next file
            bc%curstep = 1

            ! if we have run out of input files, stop with an error message
            if (bc%curfile <= size(file_list)) then
                stop "Ran out of files to process while searching for matching time variable!"
            endif

        endif
    end subroutine update_forcing_step



    !>------------------------------------------------------------
    !!  Generic routine to read a forcing variable (varname) from a netcdf file (filename) at a time step (curstep)
    !!
    !!  Data are interpolated to the high res grid either in 3D or at the boundaries only (boundary_only)
    !!  Applies modifications specificaly for U,V,T, and P variables
    !!
    !! @param highres       Allocated input array to store the result.
    !! @param filename      Name of the NetCDF file to read.
    !! @param varname       Name of the variable to read from <filename>.
    !! @param geolut        Geographic Look Up Table that defines the 2D horizontal interpolation.
    !! @param vlut          Vertical Look Up table that defines the vertical interpolation.
    !! @param curstep       The time step in <filename> to read.
    !! @param boundary_only Logical:interpolate to all points or just to boundaries
    !! @param options       Model options structure
    !! @param z_lo          Low-resolution (input) 3D vertical coordinate for pressure adjustment. [meters]
    !! @param z_hi          High-resolution (output) 3D vertical coordinate for pressure adjustment. [meters]
    !! @param time_varying_zlut Logical: if true, first adjust the low-resolution data to a common vertical coordinate
    !! @param interp_vertical Logical: if true, perform vertical interpolation (false for pressure)
    !!
    !!------------------------------------------------------------
    ! subroutine read_var(highres, filename, varname, geolut, vlut, curstep, boundary_only, options, &
    !                     z_lo, z_hi, time_varying_zlut, interp_vertical)
    !     implicit none
    !     real,dimension(:,:,:),   intent(inout):: highres
    !     character(len=*),        intent(in)   :: filename,varname
    !     type(geo_look_up_table), intent(in)   :: geolut
    !     type(vert_look_up_table),intent(in)   :: vlut
    !     integer,                 intent(in)   :: curstep
    !     logical,                 intent(in)   :: boundary_only
    !     type(options_t),      intent(in)   :: options
    !     real,dimension(:,:,:),   intent(in), optional :: z_lo, z_hi
    !     type(vert_look_up_table),intent(in), optional :: time_varying_zlut
    !     logical,                 intent(in), optional :: interp_vertical
    !
    !     ! local variables
    !     real,dimension(:,:,:),allocatable :: inputdata,extra_data
    !     integer :: nx,ny,nz,i
    !     logical :: apply_vertical_interpolation
    !
    !     ! the special case is to skip vertical interpolation (pressure and one pass of temperature)
    !     if (present(interp_vertical)) then
    !         apply_vertical_interpolation=interp_vertical
    !     else
    !         ! so default to applying interpolation
    !         apply_vertical_interpolation=.True.
    !     endif
    !
    !     ! Read the data in, should be relatively fast because we are reading a low resolution forcing file
    !     call io_read(filename,varname,inputdata,curstep)
    !
    !     nx=size(inputdata,1)
    !     ny=size(inputdata,2)
    !     nz=size(inputdata,3)
    !
    !     ! Variable specific options
    !     ! For wind variables run a first pass of smoothing over the low res data
    !     if (((varname==options%parameters%vvar).or.(varname==options%parameters%uvar)).and.(.not.options%parameters%ideal)) then
    !         call smooth_array(inputdata,1,2)
    !
    !     ! For Temperature, we may need to add an offset
    !     elseif ((varname==options%parameters%tvar).and.(options%parameters%t_offset/=0)) then
    !         inputdata=inputdata+options%parameters%t_offset
    !
    !     ! For pressure, we may need to add a base pressure offset read from pbvar
    !     else if ((varname==options%parameters%pvar).and.(options%parameters%pbvar/='')) then
    !         call io_read(filename,options%parameters%pbvar,extra_data,curstep)
    !         inputdata=inputdata+extra_data
    !         deallocate(extra_data)
    !     endif
    !
    !     ! if the z-axis of the input data varies over time, we need to first interpolate
    !     ! to the "standard" z-axis so that the hi-res vlut doesn't need to change
    !     if ((options%parameters%time_varying_z).and.present(time_varying_zlut)) then
    !         allocate(extra_data(nx,ny,nz))
    !         extra_data=inputdata
    !         call vinterp(inputdata, extra_data, time_varying_zlut, axis=3)
    !         deallocate(extra_data)
    !     endif
    !
    !
    !     ! just read the low res version with out interpolating for e.g. external wind data
    !     if ( (nx==size(highres,1)).and.(ny==size(highres,3)).and.(nz==size(highres,2)) ) then
    !         highres=reshape(inputdata,[nx,nz,ny],order=[1,3,2])
    !         deallocate(inputdata)
    !     else
    !         ! interpolate data onto the high resolution grid after re-arranging the dimensions.
    !         allocate(extra_data(nx,nz,ny))
    !         extra_data=reshape(inputdata,[nx,nz,ny],order=[1,3,2])
    !
    !         ! first interpolate to a high res grid (temporarily stored in inputdata)
    !         deallocate(inputdata)
    !         allocate(inputdata(size(highres,1),nz,size(highres,3)))
    !         call geo_interp(inputdata, &
    !                         extra_data, &
    !                         geolut,boundary_only)
    !         ! Then apply vertical interpolation on that grid
    !         if (apply_vertical_interpolation) then
    !             call vinterp(highres, inputdata, &
    !                          vlut,boundary_only)
    !         else
    !             ! if we aren't interpolating, just copy over the output
    !             highres=inputdata(:,:size(highres,2),:)
    !         endif
    !         deallocate(extra_data)
    !         deallocate(inputdata)
    !     endif
    !
    !     ! highres is the useful output of the subroutine
    ! end subroutine read_var

    !>------------------------------------------------------------
    !!  Same as read_var but for 2-dimensional data
    !!
    !!  Data are read and horizontally interpolated. This routine is simpler because it is used For
    !!  a more limited set of variables, and they are simpler 2D instead of 3D.
    !!  Primarily used for surface variables: Sensible and latent heat fluxes, PBL height, skin temperature, radiation
    !!
    !! @param highres       Allocated input array to store the result.
    !! @param filename      Name of the NetCDF file to read.
    !! @param varname       Name of the variable to read from <filename>.
    !! @param geolut        Geographic Look Up Table that defines the 2D horizontal interpolation.
    !! @param curstep       The time step in <filename> to read.
    !!
    !!------------------------------------------------------------
!     subroutine read_2dvar(highres,filename,varname,geolut,curstep)
!         implicit none
!         real,dimension(:,:),intent(inout)::highres
!         character(len=*),intent(in) :: filename,varname
!         type(geo_look_up_table),intent(in) :: geolut
!         integer,intent(in)::curstep
!
!         real,dimension(:,:),allocatable :: inputdata
!
! !       Read the data in
!         call io_read(filename,varname,inputdata,curstep)
! !       interpolate data onto the high resolution grid
!         call geo_interp2d(highres,inputdata,geolut)
!         deallocate(inputdata)
!
!     end subroutine read_2dvar

    !>------------------------------------------------------------
    !!  Read in the time step from a boundary conditions file if available
    !!
    !!  if not time_var is specified, nothing happens
    !!
    !! Should update this to save the times so they don't all have to be read again every timestep
    !!
    !! @param model_time    Double Scalar to hold time data
    !! @param filename      Name of the NetCDF file to read.
    !! @param varname       Name of the time variable to read from <filename>.
    !! @param curstep       The time step in <filename> to read.
    !!
    !!------------------------------------------------------------
    subroutine read_bc_time(model_time, filename, time_var, curstep)
        implicit none
        type(Time_type),    intent(inout) :: model_time
        character(len=*),   intent(in)    :: filename, time_var
        integer,            intent(in)    :: curstep

        type(Time_type), dimension(:), allocatable :: times

        if (time_var/="") then
            call read_times(filename, time_var, times)
            model_time = times(curstep)
            deallocate(times)
        endif
    end subroutine read_bc_time



    !>------------------------------------------------------------
    !! Check that two model grids have the same shape
    !!
    !! If the size of all dimensions in data1 and data2 are not exactly the same the model will stop
    !!
    !! @param data1     First 3D array to check dimension sizes
    !! @param data2     Second 3D array to check dimension sizes
    !!
    !!------------------------------------------------------------
    subroutine check_shapes_3d(data1,data2)
        implicit none
        real,dimension(:,:,:),intent(in)::data1,data2
        integer :: i
        do i=1,3
            if (size(data1,i).ne.size(data2,i)) then
                write(*,*) "Restart file 3D dimensions don't match domain"
                write(*,*) shape(data1)
                write(*,*) shape(data2)
                stop
            endif
        enddo
    end subroutine check_shapes_3d

    !>------------------------------------------------------------
    !! Swap the last two dimensions of an array
    !!
    !! Call reshape after finding the nx,ny,nz values
    !!
    !! @param data     3D array to be reshaped
    !!
    !!------------------------------------------------------------
    subroutine swap_y_z_dimensions(data)
        implicit none
        real,dimension(:,:,:),intent(inout),allocatable :: data
        real,dimension(:,:,:), allocatable :: temporary_data
        integer :: nx,ny,nz

        nx=size(data,1)
        ny=size(data,2)
        nz=size(data,3)
        allocate(temporary_data(nx,nz,ny))
        temporary_data = reshape(data, [nx,nz,ny], order=[1,3,2])

        deallocate(data)
        allocate(data(nx,nz,ny))
        data=temporary_data

    end subroutine swap_y_z_dimensions


    !>------------------------------------------------------------
    !!  Same as update_dxdt but only for the edges of the domains
    !!
    !!  This is used for fields that are calculated/updated internally
    !!  by the model physics (e.g. temperature and moisture)
    !!  In the output dxdt variable, the first dimension is the z axis,
    !!  The second dimension is either X or Y (which ever is specified)
    !!  And the third dimension specifies the boundary it applies to
    !!  1=left, 2=right, 3=bottom, 4=top
    !!
    !! @param dx_dt Change in variable X between time periods d1 and d2
    !! @param d1    Value of field X at time period 1
    !! @param d2    Value of field X at time period 2
    !! @retval dx_dt This field is updated along the boundaries to be (d1-d2)
    !!
    !!------------------------------------------------------------
    subroutine update_edges(dx_dt, d1, d2)
        implicit none
        real,dimension(:,:,:), intent(inout) :: dx_dt
        real,dimension(:,:,:), intent(in)    :: d1, d2
        integer :: nx, nz, ny, i

        nx = size(d1, 1)
        nz = size(d1, 2)
        ny = size(d1, 3)
        do i=1,nz
            dx_dt(i,:ny,1) = d1(1,i,:)  - d2(1,i,:)
            dx_dt(i,:ny,2) = d1(nx,i,:) - d2(nx,i,:)
            dx_dt(i,:nx,3) = d1(:,i,1)  - d2(:,i,1)
            dx_dt(i,:nx,4) = d1(:,i,ny) - d2(:,i,ny)
        enddo
    end subroutine update_edges



    !>------------------------------------------------------------
    !!  Adjust the pressure field for the vertical shift between the low and high-res domains
    !!
    !!  Ideally this should include temperature... but it isn't entirely clear
    !!  what it would mean to do that, what temperature do you use? Current even though you are adjusting future P?
    !!  Alternatively, could adjust input pressure to SLP with future T then adjust back to elevation with current T?
    !!  Currently if T is supplied, it uses the mean of the high and low-res T to split the difference.
    !!  Equations from : http://www.wmo.int/pages/prog/www/IMOP/meetings/SI/ET-Stand-1/Doc-10_Pressure-red.pdf
    !!  excerpt from CIMO Guide, Part I, Chapter 3 (Edition 2008, Updated in 2010) equation 3.2
    !!  http://www.meteormetrics.com/correctiontosealevel.htm
    !!
    !! @param pressure  The pressure field to be adjusted
    !! @param z_lo      The 3D vertical coordinate of the input pressures
    !! @param z_hi      The 3D vertical coordinate of the computed/adjusted pressures
    !! @param lowresT   OPTIONAL 3D temperature field of the input pressures
    !! @param lowresT   OPTIONAL 3D temperature field of the computed/adjusted pressures
    !! @retval pressure The pressure field after adjustment
    !!
    !!------------------------------------------------------------
    subroutine update_pressure(pressure,z_lo,z_hi, lowresT, hiresT)
        implicit none
        real,dimension(:,:,:), intent(inout) :: pressure
        real,dimension(:,:,:), intent(in) :: z_lo,z_hi
        real,dimension(:,:,:), intent(in), optional :: lowresT, hiresT
        ! local variables
        real,dimension(:),allocatable::slp !sea level pressure [Pa]
        ! vapor pressure, change in height, change in temperature with height and mean temperature
        real,dimension(:),allocatable:: dz, tmean !, e, dTdz
        integer :: nx,ny,nz,i,j
        nx=size(pressure,1)
        nz=size(pressure,2)
        ny=size(pressure,3)

        if (present(lowresT)) then
            !$omp parallel shared(pressure, z_lo,z_hi, lowresT, hiresT) &
            !$omp private(i,j, dz, tmean) firstprivate(nx,ny,nz)  !! private(e, dTdz)
            allocate(dz(nx))
            allocate(tmean(nx))
            ! allocate(e(nx))
            ! allocate(dTdz(nx))
            !$omp do
            do j=1,ny
                ! is an additional loop over z more cache friendly?
                do i=1,nz
                    ! vapor pressure
!                     e = qv(:,:,j) * pressure(:,:,j) / (0.62197+qv(:,:,j))
                    ! change in elevation (note reverse direction from "expected" because the formula is an SLP reduction)
                    dz   = (z_lo(:,i,j) - z_hi(:,i,j))
                    ! lapse rate (not sure if this should be positive or negative)
                    ! dTdz = (loresT(:,:,j) - hiresT(:,:,j)) / dz
                    ! mean temperature between levels
                    if (present(hiresT)) then
                        tmean= (hiresT(:,i,j) + lowresT(:,i,j)) / 2
                    else
                        tmean= lowresT(:,i,j)
                    endif
                    ! slp= ps*np.exp(((g/R)*Hp) / (ts - a*Hp/2.0 + e*Ch))
                    pressure(:,i,j) = pressure(:,i,j) * exp( ((gravity/Rd) * dz) / tmean )   !&
!                                         (tmean + (e * 0.12) ) )
                    ! alternative formulation M=0.029, R=8.314?
                    ! p= p0*(t0/(t0+dtdz*z))**((g*M)/(R*dtdz))
                    ! do i=1,nz
                    !     pressure(:,i,j) = pressure(:,i,j)*(t0/(tmean(:,i)+dTdz(:,i)*z))**((g*M)/(R*dtdz))
                    ! enddo
                enddo
            enddo
            !$omp end do
            deallocate(dz, tmean)
            ! deallocate(e, dTdz)
            !$omp end parallel
        else
            ! this is pretty foolish to convert to sea level pressure and back... should be done in one step
            ! just need to test that the relationship can work that way h=(zhi-zlo)
            ! this doesn't get used much though (only from bc_init) so it doesnt seem worth the time...
            !$omp parallel shared(pressure, z_lo,z_hi) &
            !$omp private(slp,i,j) firstprivate(nx,ny,nz)
            allocate(slp(nx))
            !$omp do
            do j=1,ny
                do i=1,nz
                    ! slp = pressure(:,i,j) / (1 - 2.25577E-5 * z_lo(:,i,j))**5.25588
                    pressure(:,i,j) = pressure(:,i,j) * (1 - 2.25577e-5 * (z_hi(:,i,j)-z_lo(:,i,j)))**5.25588
                enddo
            enddo
            !$omp end do
            deallocate(slp)
            !$omp end parallel
        endif
    end subroutine update_pressure

end submodule