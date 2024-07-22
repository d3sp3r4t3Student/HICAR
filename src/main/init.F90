module initialization
    use options_interface,  only : options_t
    use domain_interface,   only : domain_t
    use boundary_interface, only : boundary_t
    use microphysics,               only : mp_init
    use advection,                  only : adv_init
    use radiation,                  only : radiation_init
    use convection,                 only : init_convection
    use planetary_boundary_layer,   only : pbl_init
    use land_surface,               only : lsm_init
    use surface_layer,              only : sfc_init
    use io_routines,                only : io_read, io_write
    use mod_atm_utilities,          only : init_atm_utilities
    use wind,                       only : update_winds, init_winds
    use icar_constants!,             only : kITERATIVE_WINDS, kWIND_LINEAR
    use ioserver_interface,         only : ioserver_t
    use ioclient_interface,         only : ioclient_t
    use iso_fortran_env
    use mpi_f08

    implicit none
    private
    public :: split_processes, welcome_message, init_options, init_model, init_physics, init_model_state

contains

    subroutine unique(array, unique_array, count)
        integer, intent(in) :: array(:)
        integer, allocatable, intent(out) :: unique_array(:)
        integer, intent(out) :: count
        integer :: i, j, n, found

        n = size(array)
        count = 0

        if (n == 0) then
            ! Handle the case of an empty input array
            allocate(unique_array(0))
            return
        end if

        ! Initialize unique_array with maximum possible size
        allocate(unique_array(n))

        do i = 1, n
            found = 0
            do j = 1, count
                if (array(i) == unique_array(j)) then
                    found = 1
                    exit
                end if
            end do
            if (found == 0) then
                count = count + 1
                unique_array(count) = array(i)
            end if
        end do

        ! Resize unique_array to the exact number of unique elements
        if (count < n) then
            call resize(unique_array, count)
        end if
    end subroutine unique

    subroutine resize(array, new_size)
        integer, allocatable, intent(inout) :: array(:)
        integer, intent(in) :: new_size
        integer, allocatable :: temp(:)

        if (new_size == size(array)) return ! No resizing needed

        allocate(temp(new_size))
        temp = array(1:new_size)
        deallocate(array)
        array = temp
    end subroutine resize

    subroutine split_processes(exec_team, domain, ioserver, ioclient)
        implicit none
        integer, intent(inout) :: exec_team
        type(domain_t), intent(inout) :: domain
        type(ioserver_t), intent(inout) :: ioserver
        type(ioclient_t), intent(inout) :: ioclient
        
        integer :: n, k, name_len, color, ierr, node_name_i, num_PE
        integer :: unique_node_count, procs_per_node, remainder, node_id
        integer, allocatable :: node_names(:), unique_nodes(:)
        logical, allocatable :: node_used(:)
        character(len=MPI_MAX_PROCESSOR_NAME) :: node_name
        type(MPI_Comm) :: globalComm, splitComm
        integer :: num_threads
        
        num_threads = 1
        
        call MPI_Comm_size(MPI_COMM_WORLD, num_PE, ierr)
        
        ! Ensure dom_img_indx is allocated correctly
        allocate(dom_img_indx(num_PE))
        dom_img_indx = 0
        
        allocate(node_names(num_PE))
        allocate(node_used(num_PE))
        node_names = 0
        node_used = .false.
        node_name_i = 0
        
        ! First, determine information about node/CPU configuration
        call MPI_Get_processor_name(node_name, name_len, ierr)
        do n = 1, name_len
            node_name_i = node_name_i + ichar(node_name(n:n)) * n * 10
        end do
        
        node_names(PE_RANK_GLOBAL + 1) = node_name_i
        
        ! Get list of node names on all processes
        call MPI_Allgather(node_name_i, 1, MPI_INTEGER, node_names, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
        
        ! Determine the number of unique nodes
        call unique(node_names, unique_nodes, unique_node_count)
        
        ! Calculate the number of processes per node
        procs_per_node = num_PE / unique_node_count
        remainder = mod(num_PE, unique_node_count)
        
        ! Assign processes to nodes
        node_id = 0
        do n = 1, num_PE
            if (mod(n - 1, procs_per_node) == 0 .and. n /= 1) then
                node_id = node_id + 1
            end if
            if (node_id >= unique_node_count) then
                node_id = unique_node_count - 1
            end if
            if (node_names(n) == unique_nodes(node_id + 1)) then
                exec_team = kIO_TEAM
                color = 1
            else
                exec_team = kCOMPUTE_TEAM
                color = 0
            end if
            
            ! Ensure index is within bounds
            if (n >= 1 .and. n <= size(dom_img_indx)) then
                dom_img_indx(n) = exec_team
            else
                write(*,*) 'Error: Index out of bounds for dom_img_indx'
            end if
            
        end do
        
        ! Adjust for remainder processes
        if (remainder /= 0) then
            do n = 1, remainder
                ! Ensure index is within bounds
                if (n >= 1 .and. n <= size(dom_img_indx)) then
                    dom_img_indx(n) = kCOMPUTE_TEAM
                else
                    write(*,*) 'Error: Index out of bounds for dom_img_indx'
                end if
            end do
        end if
        
        ! Assign Compute and IO processes
        call MPI_Comm_dup(MPI_COMM_WORLD, globalComm, ierr)
        call MPI_Comm_split(globalComm, color, PE_RANK_GLOBAL, splitComm, ierr)
        
        select case (exec_team)
        case (kCOMPUTE_TEAM)
            call MPI_Comm_dup(splitComm, domain%compute_comms, ierr)
            ioserver%IO_comms = MPI_COMM_NULL
        case (kIO_TEAM)
            call MPI_Comm_dup(splitComm, ioserver%IO_comms, ierr)
            domain%compute_comms = MPI_COMM_NULL
        end select
        
        ! Group server and client processes
        call MPI_Comm_dup(MPI_COMM_WORLD, globalComm, ierr)
        call MPI_Comm_split(globalComm, node_names(PE_RANK_GLOBAL + 1), PE_RANK_GLOBAL, splitComm, ierr)
        
        select case (exec_team)
        case (kCOMPUTE_TEAM)
            call MPI_Comm_dup(splitComm, ioclient%parent_comms, ierr)
            ioserver%client_comms = MPI_COMM_NULL
        case (kIO_TEAM)
            call MPI_Comm_dup(splitComm, ioserver%client_comms, ierr)
            ioclient%parent_comms = MPI_COMM_NULL
        end select
    end subroutine split_processes
    
        

    subroutine init_options(options, namelist_file, info_only, gen_nml, only_namelist_check)
        implicit none
        type(options_t), intent(inout) :: options
        character(len=*), intent(in) :: namelist_file
        logical, optional, intent(in) :: info_only, gen_nml, only_namelist_check

        ! read in options file
        call options%init(namelist_file, info_only=info_only, gen_nml=gen_nml, only_namelist_check=only_namelist_check)

    end subroutine init_options

    subroutine init_model(options, domain, boundary)
        implicit none
        type(options_t), intent(inout) :: options
        type(domain_t),  intent(inout) :: domain
        type(boundary_t), intent(inout) :: boundary ! forcing file for init conditions

        if (STD_OUT_PE) write(*,*) "Initializing Domain"
        if (STD_OUT_PE) flush(output_unit)
        call domain%init(options)

        if (STD_OUT_PE) write(*,*) "Initializing boundary condition data structure"
        if (STD_OUT_PE) flush(output_unit)
        call boundary%init(options, domain%latitude%data_2d, domain%longitude%data_2d, domain%variables_to_force)

        if (STD_OUT_PE) write(*,*) "Initializing atmospheric utilities"
        if (STD_OUT_PE) flush(output_unit)
        ! initialize the atmospheric helper utilities
        call init_atm_utilities(options)

        if (STD_OUT_PE) write(*,'(/ A)') "Finished basic initialization"
        if (STD_OUT_PE) write(*,'(A /)') "---------------------------------------"

    end subroutine init_model

    subroutine init_physics(options, domain, forcing)
        implicit none
        type(options_t), intent(inout) :: options
        type(domain_t),  intent(inout) :: domain
        type(boundary_t), intent(in) :: forcing

        if (STD_OUT_PE) write(*,*) "Init initial winds"
        if (STD_OUT_PE) flush(output_unit)
        call init_winds(domain, options)

        if (STD_OUT_PE) write(*,*) "Updating initial winds"
        if (STD_OUT_PE) flush(output_unit)
        call update_winds(domain, forcing, options)

        ! initialize microphysics code (e.g. compute look up tables in Thompson et al)
        call mp_init(options) !this could easily be moved to init_model...
        if (STD_OUT_PE) flush(output_unit)
        call init_convection(domain, options)
        if (STD_OUT_PE) flush(output_unit)

        call pbl_init(domain, options)
        if (STD_OUT_PE) flush(output_unit)

        call radiation_init(domain, options)
        if (STD_OUT_PE) flush(output_unit)

        call lsm_init(domain, options)
        if (STD_OUT_PE) flush(output_unit)

        call sfc_init(domain, options)
        if (STD_OUT_PE) flush(output_unit)

        call adv_init(domain, options)
        if (STD_OUT_PE) flush(output_unit)

    end subroutine init_physics

    subroutine init_model_state(options, domain, boundary)
        implicit none
        type(options_t), intent(inout) :: options
        type(domain_t),  intent(inout) :: domain
        type(boundary_t), intent(inout) :: boundary ! forcing file for init conditions

        if (STD_OUT_PE) write(*,*) "Reading Initial conditions from boundary dataset"
        call domain%get_initial_conditions(boundary, options)

    end subroutine init_model_state

    subroutine welcome_message()
        implicit none

        write(*,*) ""
        write(*,*) "============================================================"
        write(*,*) "|                                                          |"
        write(*,*) "|  The Intermediate Complexity Atmospheric Research Model  |"
        write(*,*) "|                          (ICAR)                          |"
        write(*,*) "|                                                          |"
        write(*,*) "|   Developed at NCAR:                                     |"
        write(*,*) "|     The National Center for Atmospheric Research         |"
        write(*,*) "|     NCAR is sponsored by the National Science Foundation |"
        write(*,*) "|                                                          |"
        write(*,*) "|   Version: ", kVERSION_STRING, "                         |"
        write(*,*) "|                                                          |"
        write(*,*) "============================================================"
        write(*,*) ""

    end subroutine welcome_message

end module
