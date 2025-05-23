!>------------------------------------------------------------
!! Module to pre-compute and apply near-surface corrections to 
!! the domain wind field. This includes topographic 
!! exposure/sheltering via the Sx parameter (Winstral et al. 2017) 
!! and surface roughness.
!!
!!  @author
!!  Dylan Reynolds (dylan.reynolds@slf.ch)
!!
!!------------------------------------------------------------
module wind_surf
    use domain_interface,  only : domain_t
    use options_interface, only : options_t
    use io_routines,       only : io_write
    use mod_atm_utilities,     only : calc_thresh_ang
    use icar_constants,    only : kVARS, STD_OUT_PE

    implicit none
    private
    public :: calc_Sx, calc_TPI, apply_Sx
    real, parameter :: deg2rad=0.017453293 !2*pi/360
    real, parameter :: rad2deg=57.2957779371
    integer, save :: TPI_k_max, Sx_k_max = 0
    real, save    :: TPI_scale, Sx_scale_ang
    ! real, parameter :: SX_SCALE_ANG = 30.0 !This means that for a (SX_SCALE_ANG/2)-degree difference between threshold and Sx, 
                            !reversal starts. Can be thought of as minimum terrain slope necesarry for flow reversal
    ! real, parameter :: TPI_SCALE = 200.0
    real, parameter :: SX_Z_MAX = 1500.0
    real, parameter :: TPI_Z_MAX = 200.0
                                              
    
contains

    subroutine calc_TPI(domain, options, neighbor_TPI)
        implicit none
        class(domain_t),  intent(inout) :: domain
        class(options_t), intent(in)    :: options
        real, allocatable, intent(out)  :: neighbor_TPI(:,:)

        integer           :: Sx_search_max, TPI_search_max, i, j, i_s, j_s, i_start_buff, i_end_buff, j_start_buff, j_end_buff
        integer           :: ips, ipe, jps, jpe, TPI_num
        real              :: search_height, TPI_sum, dist, TPI_d_max, Sx_d_max
        
        TPI_d_max = options%wind%TPI_dmax
        TPI_search_max = min(floor(max(1.0,TPI_d_max/domain%dx)),domain%neighborhood_max)

        Sx_d_max = options%wind%Sx_dmax
        Sx_search_max = min(floor(max(1.0,Sx_d_max/domain%dx)),domain%neighborhood_max)

        ips = max(domain%ims - Sx_search_max,domain%ihs); ipe = min(domain%ime + Sx_search_max,domain%ihe);
        jps = max(domain%jms - Sx_search_max,domain%jhs); jpe = min(domain%jme + Sx_search_max,domain%jhe);

        allocate(neighbor_TPI(ips:ipe,jps:jpe))
              
        neighbor_TPI = 0
        
        !Now calc TPI
        do i=ips, ipe
            do j=jps, jpe
                TPI_num = 0
                TPI_sum = 0
                    
                ! Check to use buffers to avoid searching out of grid
                i_start_buff = max(domain%ihs,i-TPI_search_max)
                i_end_buff = min(domain%ihe,i+TPI_search_max)
                
                j_start_buff = max(domain%jhs,j-TPI_search_max)
                j_end_buff = min(domain%jhe,j+TPI_search_max)
                
                do i_s = i_start_buff, i_end_buff
                    do j_s = j_start_buff, j_end_buff
                    
                        !Determine distance of search point
                        dist = domain%dx*sqrt(real((i-i_s)**2+(j-j_s)**2))
                        
                        if (dist <= TPI_d_max .and. .not.(dist == 0) ) then
                            search_height = domain%vars_2d(domain%var_indx(kVARS%neighbor_terrain)%v)%data_2d(i_s,j_s)
                        
                            TPI_sum = TPI_sum + search_height
                            TPI_num = TPI_num + 1
                        end if                               
                    end do
                end do
                if (TPI_num > 0) neighbor_TPI(i,j) = domain%vars_2d(domain%var_indx(kVARS%neighbor_terrain)%v)%data_2d(i,j) - TPI_sum/TPI_num
            end do
        end do
        
        domain%vars_2d(domain%var_indx(kVARS%TPI)%v)%data_2d = neighbor_TPI(domain%grid2d%ims:domain%grid2d%ime,domain%grid2d%jms:domain%grid2d%jme)
        
        
        ! if ( STD_OUT_PE ) then
        !     write (*,*) "Saving *_TPI.nc"
        !     !Save file
        !     call io_write("neighbor_TPI.nc", "TPI", domain%vars_2d(domain%var_indx(kVARS%neighbor_TPI)%v)%data_2d(:,:) ) 
        ! endif
             
    end subroutine calc_TPI
    
    subroutine calc_Sx(domain, options)
        implicit none
        class(domain_t), intent(inout) :: domain
        type(options_t),intent(in)    :: options
        
        real, allocatable    :: Sx_array_temp(:,:,:,:), sheltering_TPI(:,:,:,:), temp_sheltering_TPI(:,:,:,:), neighbor_TPI(:,:)
        integer           :: search_max, i, j, k, ang, i_s, j_s, i_start_buff, i_end_buff, j_start_buff, j_end_buff, TPI_i_s, TPI_i_e, TPI_j_s, TPI_j_e
        integer           :: rear_ang, fore_ang, test_ang, rear_ang_diff, fore_ang_diff, ang_diff, k_max, window_rear, window_fore, maxSxLoc, window_width
        integer           :: azm_index
        real              :: search_height, pt_height, h_diff, Sx_temp, maxSxVal, TPI_Shelter_temp, exposed_TPI, z_mean, dist, azm, d_max

        
        if (Sx_k_max==0 .or. TPI_k_max==0) then
            do k = domain%grid%kms,domain%grid%kme
                z_mean = SUM(options%domain%dz_levels(1:k))
                if (z_mean > SX_Z_MAX .and. Sx_k_max==0) Sx_k_max = max(2,k-1)
                if (z_mean > TPI_Z_MAX .and. TPI_k_max==0) TPI_k_max = max(2,k-1)
            enddo
        endif

        Sx_k_max = max(Sx_k_max,TPI_k_max)

        call calc_TPI(domain, options,neighbor_TPI)

        TPI_i_s = lbound(neighbor_TPI,1)
        TPI_i_e = ubound(neighbor_TPI,1)
        TPI_j_s = lbound(neighbor_TPI,2)
        TPI_j_e = ubound(neighbor_TPI,2)

        d_max = options%wind%Sx_dmax
        TPI_scale = options%wind%TPI_scale
        Sx_scale_ang = options%wind%Sx_scale_ang

        search_max = min(floor(max(1.0,d_max/domain%dx)),domain%neighborhood_max)
               
        exposed_TPI = 10.0
        
        ! Initialize Sx and set to value of -90 (absolute minimum possible) for search algorithm
        ! allocate(domain%vars_4d(domain%var_indx(kVARS%Sx)%v)%data_4d( 1:72, domain%grid2d%ims:domain%grid2d%ime, 1:Sx_k_max, domain%grid2d%jms:domain%grid2d%jme ))
        allocate(Sx_array_temp(domain%grid2d%ims:domain%grid2d%ime, 1:Sx_k_max, domain%grid2d%jms:domain%grid2d%jme,  1:72))
        
        Sx_array_temp = -90.0
        domain%vars_4d(domain%var_indx(kVARS%Sx)%v)%data_4d = -90.0
        
        allocate(temp_sheltering_TPI( domain%grid2d%ims:domain%grid2d%ime, 1:Sx_k_max, domain%grid2d%jms:domain%grid2d%jme, 1:72 ))
        allocate(sheltering_TPI( domain%grid2d%ims:domain%grid2d%ime, 1:Sx_k_max, domain%grid2d%jms:domain%grid2d%jme, 1:72 ))

        temp_sheltering_TPI = 0.0
        sheltering_TPI = 0.0
                     
                
        !call unique_sort_ind(azm_indices,valid_ks)
        
        do i=domain%grid2d%ims, domain%grid2d%ime
            do j=domain%grid2d%jms, domain%grid2d%jme
                do k = 1, Sx_k_max

                    if (k == 1) then
                        pt_height = domain%vars_2d(domain%var_indx(kVARS%terrain)%v)%data_2d(i,j)
                    else if (k > 1) then
                        pt_height = pt_height + domain%vars_3d(domain%var_indx(kVARS%dz_interface)%v)%data_3d(i,k,j)
                    end if
                    
                    ! Check to use buffers to avoid searching out of grid
                    i_start_buff = max(TPI_i_s,i-search_max)
                    i_end_buff = min(TPI_i_e,i+search_max)

                    j_start_buff = max(TPI_j_s,j-search_max)
                    j_end_buff = min(TPI_j_e,j+search_max)
                
                    do i_s = i_start_buff, i_end_buff
                        do j_s = j_start_buff, j_end_buff
                            
                            !Determine distance of search point
                            dist = domain%dx*sqrt(real((i-i_s)**2+(j-j_s)**2))
                            
                            !Since dist is for a grid, and we want a search RADIUS, check that we are within d_max
                            if (dist <= d_max .and. .not.(dist == 0) ) then
                            
                                !Compute azimuth ind of point
                                azm = atan2(1.0*(i_s-i),1.0*(j_s-j))*rad2deg
                                if(azm < 0) then
                                    azm = 360+azm
                                else if(azm >= 360.0) then
                                    azm=0.0
                                endif
                                azm_index = int(azm/5)+1

                                !Calculate height difference
                                search_height = domain%vars_2d(domain%var_indx(kVARS%neighbor_terrain)%v)%data_2d(i_s,j_s)
                                h_diff = search_height - pt_height

                                !Calculate Sx slope to search-cell
                                Sx_temp = atan(h_diff/dist)*rad2deg
                                TPI_Shelter_temp = neighbor_TPI(i_s,j_s)
                            
                                ! If new Sx is greater than existing Sx for a given search angle, replace
                                if (Sx_temp > Sx_array_temp(i,k,j,azm_index) ) then
                                    !If we have found a sheltering Sx, and it is "exposed" (TPI >= 100.0), use this Sx
                                    !If we have found a sheltering Sx, but it is not "exposed", then this cell gets neither sheltering nor exposure (Sx = 0)
                                    !Sx_array_temp(azm_isndices(i_s,j_s),i,k,j) = Sx_temp
                                    if ( (Sx_temp > 0.0) .and. (TPI_Shelter_temp >= exposed_TPI) ) then ! .and. ( TPI_Shelter_temp >= temp_sheltering_TPI(azm_index,i,k,j)) ) then
                                        Sx_array_temp(i,k,j,azm_index) = Sx_temp
                                        temp_sheltering_TPI(i,k,j,azm_index) = TPI_Shelter_temp
                                    else if ( (Sx_temp <= 0.0) ) then !  .and. (neighbor_TPI(i,j) >= exposed_TPI) ) then
                                        Sx_array_temp(i,k,j,azm_index) = Sx_temp
                                    end if
                                end if
                            end if
                        end do
                    end do
                                    
                    !After finding Sx in each absolute direction around grid cell, 
                    !Pick max for each 30º window and perform interpolation to other directions if necesarry
                    
                    rear_ang = 1 
                    fore_ang = 1
                    
                    if (.not.( all((Sx_array_temp(i,k,j,:) <= -90.0)) )) then
                    
                        !Perform 40º window max search
                        window_width = 1
                        do ang = 1, 72
                            window_rear = ang-window_width
                            window_fore = ang+window_width
                        
                            if (ang <= window_width) then
                                window_rear = 72-(window_width-ang)
                                
                                maxSxVal = maxval(Sx_array_temp(i,k,j,window_rear:72))
                                maxSxLoc = maxloc(Sx_array_temp(i,k,j,window_rear:72),dim=1)
                                maxSxLoc = maxSxLoc + (window_rear - 1)
                                
                                !if (maxSxVal == -90.0) write(*,*) ang, i, k, j
                                
                                if (maxval(Sx_array_temp(i,k,j,1:window_fore)) > maxSxVal) then
                                    maxSxVal = maxval(Sx_array_temp(i,k,j,1:window_fore))
                                    maxSxLoc = maxloc(Sx_array_temp(i,k,j,1:window_fore),dim=1)
                                end if
                            else if ( ang >= (72-(window_width-1)) ) then
                                window_fore = window_width-(72-ang)
                                
                                maxSxVal = maxval(Sx_array_temp(i,k,j,window_rear:72))
                                maxSxLoc = maxloc(Sx_array_temp(i,k,j,window_rear:72),dim=1)
                                maxSxLoc = maxSxLoc + (window_rear - 1)
                                
                                if (maxval(Sx_array_temp(i,k,j,1:window_fore)) > maxSxVal) then
                                    maxSxVal = maxval(Sx_array_temp(i,k,j,1:window_fore))
                                    maxSxLoc = maxloc(Sx_array_temp(i,k,j,1:window_fore),dim=1)
                                end if
                            else
                                maxSxVal = maxval(Sx_array_temp(i,k,j,window_rear:window_fore))
                                maxSxLoc = maxloc(Sx_array_temp(i,k,j,window_rear:window_fore),dim=1)
                                maxSxLoc = maxSxLoc + (window_rear - 1)
                            end if
                            domain%vars_4d(domain%var_indx(kVARS%Sx)%v)%data_4d(i,k,j,ang) = maxSxVal
                            sheltering_TPI(i,k,j,ang) = temp_sheltering_TPI(i,k,j,maxSxLoc)
                        end do                    
                    
                        do ang = 1, 72
                            !Determine indices for interpolation
                            if ( (ang==fore_ang) ) then
                                !Update indices for interpolated Sx's
                                rear_ang = ang
                            
                                fore_ang = ang+1
                                if (fore_ang > 72) fore_ang = 1
                                
                                do while (domain%vars_4d(domain%var_indx(kVARS%Sx)%v)%data_4d(i,k,j,fore_ang) <= -90.0)
                                    fore_ang = fore_ang+1
                                    if (fore_ang > 72) fore_ang = 1
                                end do
                            
                            end if
                            
                            if (ang==1) then
                                rear_ang = 72
                                do while(domain%vars_4d(domain%var_indx(kVARS%Sx)%v)%data_4d(i,k,j,rear_ang) <= -90)
                                    rear_ang = rear_ang-1
                                end do
                            end if
                    
                            !If we did not calculate Sx for a given direction
                            if (domain%vars_4d(domain%var_indx(kVARS%Sx)%v)%data_4d(i,k,j,ang) == -90.0) then
                                !Weight the two surrounding Sx values based on our angular-distance to them
                                rear_ang_diff = ang-rear_ang
                                fore_ang_diff = fore_ang-ang
                                ang_diff = fore_ang-rear_ang
                        
                                !Handle wrap-around case
                                if (ang > fore_ang) then
                                    fore_ang_diff = fore_ang+(72-ang)
                                    ang_diff = fore_ang+(72-rear_ang)
                                end if
                        
                                !Interpolation, linearly-weighted by angular-distance from values
                                domain%vars_4d(domain%var_indx(kVARS%Sx)%v)%data_4d(i,k,j,ang) = (domain%vars_4d(domain%var_indx(kVARS%Sx)%v)%data_4d(i,k,j,rear_ang)*fore_ang_diff + &
                                                    domain%vars_4d(domain%var_indx(kVARS%Sx)%v)%data_4d(i,k,j,fore_ang)*rear_ang_diff)/ang_diff

                                !if (domain%Sx(ang,i,k,j) > 0) sheltering_TPI(ang,i,k,j) = (sheltering_TPI(rear_ang,i,k,j)*fore_ang_diff + &
                                !                    sheltering_TPI(fore_ang,i,k,j)*rear_ang_diff)/ang_diff

                            end if
                        end do

                    else
                        !IF we only have -90 for all entries, set to 0
                        domain%vars_4d(domain%var_indx(kVARS%Sx)%v)%data_4d(i,k,j,:) = 0.0
                    end if
                end do
            end do
        end do
                
        !Finally, mask out positive (sheltered) Sx where the sheltering element is not "prominent"        
        !       , mask out negative (exposed) Sx shere the element is not "porminent"
        do i=domain%grid2d%ims, domain%grid2d%ime
            do j=domain%grid2d%jms, domain%grid2d%jme
                do k = 1, Sx_k_max
                    do ang = 1, 72
                        if( (domain%vars_4d(domain%var_indx(kVARS%Sx)%v)%data_4d(i,k,j,ang) < 0.0) .and. (domain%vars_2d(domain%var_indx(kVARS%TPI)%v)%data_2d(i,j) <= exposed_TPI) ) domain%vars_4d(domain%var_indx(kVARS%Sx)%v)%data_4d(i,k,j,ang) = 0.0
                        if( (domain%vars_4d(domain%var_indx(kVARS%Sx)%v)%data_4d(i,k,j,ang) > 0.0) .and. (sheltering_TPI(i,k,j,ang) <= exposed_TPI) ) domain%vars_4d(domain%var_indx(kVARS%Sx)%v)%data_4d(i,k,j,ang) = 0.0
                    end do
                end do
            end do
        end do
        
        ! if ( STD_OUT_PE ) then
        !     write (*,*) "Saving *_Sx.nc"
        !     !Save file
        !     call io_write(filename, "Sx", domain%Sx(:,:,:,:) ) 
        !     call io_write("TPI_out.nc", "TPI", domain%vars_2d(domain%var_indx(kVARS%neighbor_TPI)%v)%data_2d(:,:) ) 
        !     call io_write("sheltering_TPI.nc", "Sx_shelter", sheltering_TPI(:,:,:,:) ) 
        ! endif
        
        deallocate(Sx_array_temp)
        
        !Sync images befoer exiting so that we don't try to read Sx while it is being written
        !sync images ([DOM_IMG_INDX])

    end subroutine calc_Sx
    
    subroutine unique_sort_ind(arr,us_arr)
        implicit none
        integer, allocatable, intent(in)    :: arr(:,:)
        integer, allocatable, intent(inout) :: us_arr(:)
        real, allocatable                :: unique(:)
        integer   :: i
        real      :: min_arr, max_arr
        
        allocate(unique(1:size(arr)))
        
        i = 0
        min_arr = minval(arr)-1
        max_arr = maxval(arr)
        
        do while (min_arr < max_arr)
            i = i+1
            min_arr = minval(arr, mask=(arr>min_arr))
            unique(i) = min_arr
        end do
        
        if (allocated(us_arr)) deallocate(us_arr)
        allocate(us_arr(1:i))
        us_arr = unique(1:i)
        
        deallocate(unique)
       
        
    end subroutine unique_sort_ind
    
    

    subroutine apply_Sx(Sx, TPI, u, v, Ri, dzdx, dzdy, ims, ime, kms, kme, jms, jme, its, ite, jts, jte)
        implicit none
        real, intent(in)                       :: Sx(ims:ime,kms:kme,jms:jme,1:72), TPI(ims:ime,jms:jme), Ri(ims:ime,kms:kme,jms:jme), dzdx(ims:ime,kms:kme,jms:jme), dzdy(ims:ime,kms:kme,jms:jme)
        real, intent(inout) :: u(ims:ime+1,kms:kme,jms:jme), v(ims:ime,kms:kme,jms:jme+1)
        integer, intent(in)                   :: its, ite, jts, jte, ims, ime, kms, kme, jms, jme
        
        real, allocatable, dimension(:,:)   :: x_norm, y_norm, thresh_ang
        real, allocatable, dimension(:,:,:) :: Sx_U_corr, Sx_V_corr, Sx_curr, Sx_corr, TPI_corr

        integer ::  i, j, k
        real    ::  Ri_num, WS, max_spd
                        
        
        allocate(x_norm(ims:ime,jms:jme))
        allocate(y_norm(ims:ime,jms:jme))
        allocate(thresh_ang(ims:ime,jms:jme))
       
        allocate(Sx_curr(ims:ime,kms:Sx_k_max,jms:jme))
        allocate(Sx_corr(ims:ime,kms:Sx_k_max,jms:jme))
        allocate(TPI_corr(ims:ime,kms:Sx_k_max,jms:jme))
        allocate(Sx_U_corr(ims:ime,kms:Sx_k_max,jms:jme))
        allocate(Sx_V_corr(ims:ime,kms:Sx_k_max,jms:jme))

        !Initialize Sx_curr. This will keep the border values from being updated
        Sx_curr = 0
        
        thresh_ang = 0.
        max_spd = 0.
        
        !Calculate horizontal normal vector components of terrain
        thresh_ang = sqrt( dzdx(ims:ime,kms,jms:jme)**2+dzdy(ims:ime,kms,jms:jme)**2)
        where(.not.(thresh_ang == 0.)) x_norm = -dzdx(ims:ime,kms,jms:jme)/thresh_ang
        where(.not.(thresh_ang == 0.)) y_norm = -dzdy(ims:ime,kms,jms:jme)/thresh_ang
        thresh_ang = 0.
        
        !Pick appropriate Sx for wind direction
        do k = kms, Sx_k_max
            call pick_Sx(Sx(:,k,:,:), Sx_curr(:,k,:), u(:,k,:), v(:,k,:),ims,ime,jms,jme)
        end do
        
        !Initialize Sx and TPI corr
        Sx_corr = 0
        TPI_corr = 0
        
        !Loop through i,j
        do i = ims, ime
            do j = jms, jme
                !Loop through vertical column
                do k = kms, Sx_k_max
                    !Compute threshold separation angle from atmospheric conditions
                    !If we are potentially sheltered (on a lee-slope)
                    if (Sx_curr(i,1,j) > 0)  then
                
                        Ri_num = Ri(i,1,j)
                        WS = sqrt( u(i,1,j)**2 + v(i,1,j)**2 )
                        thresh_ang(i,j) = calc_thresh_ang(Ri_num,WS)
                        
                        !if (Ri_num <= 0.) then
                        !    max_spd = 0.3*Ri_num+0.3
                        !else if (Ri_num < 0.25) then
                        !    max_spd = 1/(13.9*(Ri_num+0.15)) - 0.18
                        !else
                        !    max_spd = 0.0
                        !endif
                        !max_spd = max(max_spd,0.0)
              
                        !Consider sheltering corrections
                        !If surface was sheltered and we are still sheltered
                        if ((Sx_curr(i,kms,j) > thresh_ang(i,j)) .and. (Sx_curr(i,k,j) > thresh_ang(i,j))) then
                            !Sheltered correction
                            Sx_corr(i,k,j) = (Sx_curr(i,k,j)-thresh_ang(i,j))/Sx_scale_ang
                        end if
                    endif

                
                    !Compute TPI_corr, only for negative TPI's (valleys)
                    if ((k <= TPI_k_max) .and. (TPI(i,j) < 0)) then
                        !Scale TPI correction and exposure with height so we smoothly merge with forcing data
                        TPI_corr(i,k,j) = (TPI(i,j)/TPI_scale) * (1.0*(TPI_k_max+1-k)/TPI_k_max) !This is setup such that valleys <= -100
                    end if
                end do
            end do
        end do
        !Dummy bounding of corrections, for safety
        TPI_corr  = min(max(TPI_corr,-0.5),0.0)
        Sx_corr = min(max(Sx_corr,0.0),1.0)
    
        do k=kms,Sx_k_max
            Sx_corr(:,k,:) = 2*Sx_corr(:,k,:)*(x_norm*( (u(ims+1:ime+1,k,:)+u(ims:ime,k,:))/2 ) + &
                                 y_norm*( (v(:,k,jms+1:jme+1)+v(:,k,jms:jme))/2 ))
            Sx_U_corr(:,k,:) = Sx_corr(:,k,:)*x_norm
            Sx_V_corr(:,k,:) = Sx_corr(:,k,:)*y_norm
        enddo
        !Finally, apply TPI and Sx corrections, staggering corrections to U/V grids
        u(ims+1:ime,kms:Sx_k_max,:) = u(ims+1:ime,kms:Sx_k_max,:) - ( (Sx_U_corr(ims+1:ime,:,:) + Sx_U_corr(ims:ime-1,:,:))/2 )
        u(ims,kms:Sx_k_max,:) = u(ims,kms:Sx_k_max,:) - Sx_U_corr(ims,:,:)
        u(ime+1,kms:Sx_k_max,:) = u(ime+1,kms:Sx_k_max,:) - Sx_U_corr(ime,:,:)
        
        v(:,kms:Sx_k_max,jms+1:jme) = v(:,kms:Sx_k_max,jms+1:jme) - ( (Sx_V_corr(:,:,jms+1:jme) + Sx_V_corr(:,:,jms:jme-1))/2 )
        v(:,kms:Sx_k_max,jms) = v(:,kms:Sx_k_max,jms) - Sx_V_corr(:,:,jms)
        v(:,kms:Sx_k_max,jme+1) = v(:,kms:Sx_k_max,jme+1) - Sx_V_corr(:,:,jme)

        u(ims+1:ime,kms:Sx_k_max,:) = u(ims+1:ime,kms:Sx_k_max,:) * (1 + (TPI_corr(ims+1:ime,:,:)+TPI_corr(ims:ime-1,:,:))/2 )
        u(ims,kms:Sx_k_max,:) = u(ims,kms:Sx_k_max,:) * (1+TPI_corr(ims,:,:))
        u(ime+1,kms:Sx_k_max,:) = u(ime+1,kms:Sx_k_max,:) * (1+TPI_corr(ime,:,:))
        
        v(:,kms:Sx_k_max,jms+1:jme) = v(:,kms:Sx_k_max,jms+1:jme) * (1 + (TPI_corr(:,:,jms+1:jme)+TPI_corr(:,:,jms:jme-1))/2 )
        v(:,kms:Sx_k_max,jms) = v(:,kms:Sx_k_max,jms) * (1+TPI_corr(:,:,jms))
        v(:,kms:Sx_k_max,jme+1) = v(:,kms:Sx_k_max,jme+1) * (1+TPI_corr(:,:,jme))

    end subroutine apply_Sx
    
    subroutine pick_Sx(Sx,Sx_curr,u,v,ims,ime,jms,jme)
        implicit none
        real, intent(in)       :: Sx(ims:ime,jms:jme,1:72)
        real, intent(inout)    :: Sx_curr(ims:ime,jms:jme)
        real, intent(in)       :: u(ims:ime+1,jms:jme), v(ims:ime,jms:jme+1)
        integer, intent(in)    :: ims, ime, jms, jme

        real, allocatable, dimension(:,:)   :: winddir, u_m, v_m
        integer, allocatable                ::  dir_indices(:,:)
        integer ::  i, j
        
                
        allocate(winddir(ims:ime,jms:jme))
        allocate(u_m(ims:ime,jms:jme))
        allocate(v_m(ims:ime,jms:jme))
        allocate(dir_indices(ims:ime,jms:jme))

        u_m = (u(ims:ime,:) + u(ims+1:ime+1,:))/2
        v_m = (v(:,jms:jme) + v(:,jms+1:jme+1))/2
        
        !Compute wind direction for each cell on mass grid
        winddir = atan2(-u_m,-v_m)*rad2deg
        where(winddir < 0.0) winddir = winddir+360
        where(winddir == 360.0) winddir = 0.0
        dir_indices = int(winddir/5)+1
                
        !Build grid of Sx values based on wind direction at that cell
        do i = ims, ime
            do j = jms, jme
                Sx_curr(i,j) = Sx(i,j,dir_indices(i,j))
            end do
        end do
        

        
    end subroutine pick_Sx

    ! function calc_Sx_dim_extent(options) result(Sx_k_extent)
    !     implicit none
    !     class(options_t), intent(in) :: options
    !     integer :: Sx_k_extent

    !     if (Sx_k_max==0 .or. TPI_k_max==0) then
    !         do k = domain%grid%kms,domain%grid%kme
    !             z_mean = SUM(options%domain%dz_levels(1:k))
    !             if (z_mean > SX_Z_MAX .and. Sx_k_max==0) Sx_k_max = max(2,k-1)
    !             if (z_mean > TPI_Z_MAX .and. TPI_k_max==0) TPI_k_max = max(2,k-1)
    !         enddo
    !     endif

    !     Sx_k_max = max(Sx_k_max,TPI_k_max)

    ! end function calc_Sx_dim_extent

end module wind_surf
