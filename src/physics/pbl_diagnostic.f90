!>----------------------------------------------------------
!!  Simple PBL diffusion package for ICAR
!!
!!  Local-K diffusion type PBL as in Louis (1979) as documented in Hong and Pan (1996) = HP96
!!  Hong and Pan used this for their free atmosphere diffusion, but noted differences
!!  used in the "current operational model" notably the asymptotic length scale lambda
!!
!! <pre>
!! HP96 = Hong,S.-Y. and H.-L. Pan (1996) Monthly Weather Review v127 p2322
!!       Nonlocal Boundary Layer Vertical Diffusion in a Medium Range Forecast Model
!!
!! Implemented with K,shear,stability... on half levels
!!  rho on half levels for f=k*rho*dq/dz*dt
!!  rho on full levels for q=q+f/rho
!!   q,U,V on full levels
!! </pre>
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!----------------------------------------------------------
module pbl_diagnostic
    use data_structures
    use domain_interface,   only : domain_t
    use options_interface,  only : options_t
    use io_routines,          only : io_write

    private
    public :: diagnostic_pbl, finalize_diagnostic_pbl, init_diagnostic_pbl

!   NOTE *_m indicates module level variables
!   these variables are declared as module level variables so that they do not need to be allocated
!   deallocated, and re-allocated all the time.

!   gradient in the virtual potential temperature
    real, allocatable, dimension(:,:,:) :: BVF
!   gradient in the richardson number
    real, allocatable, dimension(:,:,:) :: rig_m
!   gradient in the richardson number
    real, allocatable, dimension(:,:,:) :: ri_flux
!   vertical wind shear (dwind / dz)
    real, allocatable, dimension(:,:,:) :: shear_m
!   atmospheric stability
    real, allocatable, dimension(:,:,:) :: stability_m
!   atmospheric stability
    real, allocatable, dimension(:,:,:) :: stability_h
!   length scale that asymptotes from karman*z to l0 (250m below)
    real, allocatable, dimension(:,:,:) :: l_m
!   diffusion term for scalars (K/prandtl)
    real, allocatable, dimension(:,:,:) :: Kq_m
!   input qv field to use in calculateing the qv_pbl_tendency
    real, allocatable, dimension(:,:,:) :: lastqv_m

    integer :: ims, ime, jms, jme, kms, kme
    integer :: ids, ide, jds, jde, kds, kde


    real, parameter :: asymp_length_scale = 500.0 !m from COSMO doccumentation (Part 2 section 3) 
    real, parameter :: N_substeps=100. ! number of substeps to allow (puts a cap on K to match CFL)


contains
    subroutine diagnostic_pbl(th, qv, cloud, ice, qrain, qsnow, lh, sh, um, vm, pii, rho, z, adv_dz, dz, jaco, jaco_w, terrain, its, ite, jts, jte, kts, kte_in, dt, tend_qv_pbl)
        real,   intent(inout), dimension(ims:ime, kms:kme, jms:jme) :: th            ! potential temperature [K]
        real,   intent(inout), dimension(ims:ime, kms:kme, jms:jme) :: qv            ! water vapor mixing ratio [kg/kg]
        real,   intent(inout), dimension(ims:ime, kms:kme, jms:jme) :: cloud         ! cloud water mixing ratio [kg/kg]
        real,   intent(inout), dimension(ims:ime, kms:kme, jms:jme) :: ice           ! cloud ice mixing ratio [kg/kg]
        real,   intent(inout), dimension(ims:ime, kms:kme, jms:jme) :: qrain         ! rain water mixing ratio [kg/kg]
        real,   intent(inout), dimension(ims:ime, kms:kme, jms:jme) :: qsnow         ! snow mixing ratio [kg/kg]
        real,   intent(in), dimension(ims:ime, jms:jme) :: lh         ! snow mixing ratio [kg/kg]
        real,   intent(in), dimension(ims:ime, jms:jme) :: sh         ! snow mixing ratio [kg/kg]
        real,   intent(in),    dimension(ims:ime, kms:kme, jms:jme) :: um            ! east-west wind speed on mass grid [m/s]
        real,   intent(in),    dimension(ims:ime, kms:kme, jms:jme) :: vm            ! north south wind speed on mass grid [m/s]
        real,   intent(in),    dimension(ims:ime, kms:kme, jms:jme) :: pii           ! exner function
        real,   intent(in),    dimension(ims:ime, kms:kme, jms:jme) :: rho           ! air density [kg / m^3]
        real,   intent(in),    dimension(ims:ime, kms:kme, jms:jme) :: z             ! model level heights [m]
        real,   intent(in),    dimension(ims:ime, kms:kme, jms:jme) :: adv_dz        ! computational model level thickness [m]
        real,   intent(in),    dimension(ims:ime, kms:kme, jms:jme) :: dz            ! physical model level thickness [m]
        real,   intent(in),    dimension(ims:ime, kms:kme, jms:jme) :: jaco          ! mass-centered jacobian
        real,   intent(in),    dimension(ims:ime, kms:kme, jms:jme) :: jaco_w        ! w-staggered jacobian
        real,   intent(in),    dimension(ims:ime, jms:jme)          :: terrain       ! terrain height above sea level [m]
        integer,intent(in) :: its, ite, jts, jte, kts, kte_in
        real,   intent(in) :: dt                                 ! time step [s]
        real,   intent(inout), dimension(ims:ime,kms:kme,jms:jme), optional :: tend_qv_pbl   ! output water vapor tendency [kg/kg/s] (for use in other physics)

        ! local
        real, dimension(ims:ime,kms:kme-1,jms:jme) :: dz_mass_i
        real, dimension(ims:ime)                   :: qv_flux, th_flux

        integer :: i,j,k, kte

        ! don't process the top model level regardless, there shouldn't be any real pbl diffusion occuring there
        kte = min(kme, kte_in)
        Kq_m = 1.0
        
        !For stability calculations, use physical dz
        dz_mass_i(:,kms:kme-1,:) = ((dz(:,kms:kme-1,:) + dz(:,kms+1:kme,:))/2)
        
        call calc_shear_sq(um, vm, dz_mass_i,kts,kte)

        call calc_BVF(th, qv, dz_mass_i, cloud, ice, qrain, qsnow, kts, kte)

        call calc_pbl_stability_function()
        !For diffusion calculations, use advection dz
        dz_mass_i(:,kms:kme-1,:) = ((adv_dz(:,kms:kme-1,:) + adv_dz(:,kms+1:kme,:))/2)

        do j = jts, jte
            do k = kts, kte
                lastqv_m(:,k,j) = qv(:,k,j)

                ! from eqn 12 in HP96
                !l_m(its:ite,k,j) = 1 / (1/(karman*(z(its:ite,k,j) - terrain(its:ite,j))) + asymp_length_scale)
                l_m(its:ite,k,j) = karman * (z(its:ite,k,j) - terrain(its:ite,j))/(1+(z(its:ite,k,j) - terrain(its:ite,j))/asymp_length_scale)
                where (stability_m(its:ite,k,j) < 0)
                    Kq_m(its:ite,k,j) = (l_m(its:ite,k,j)**2)  * sqrt(shear_m(its:ite,k,j)) * 0.007
                else where (stability_m(its:ite,k,j) >= 0)
                    Kq_m(its:ite,k,j) = l_m(its:ite,k,j)**2 * stability_m(its:ite,k,j)**(3./2.) * &
                                    sqrt(max(0.0, (shear_m(its:ite,k,j)  - stability_h(its:ite,k,j) * BVF(its:ite,k,j)) ))
                    ! diffusion for scalars
                    Kq_m(its:ite,k,j) = Kq_m(its:ite,k,j) * stability_h(its:ite,k,j)
                end where

            enddo
            
            th_flux = -(sh(ims:ime,j) * dt/cp)  &
             / (pii(ims:ime,kts,j))              
            qv_flux = -(lh(ims:ime,j) * dt / LH_vaporization )

            !Stagger Kq_m to vertical faces
            where ( Kq_m(its:ite,kts:kte,j) < 1.0) Kq_m(its:ite,kts:kte,j) = 1.0
            Kq_m(its:ite,kts:kte-1,j) = (Kq_m(its:ite,kts+1:kte,j) + Kq_m(its:ite,kts:kte-1,j))/2
            Kq_m(its:ite,kts:kte,j) = Kq_m(its:ite,kts:kte,j) * dt
            where ( Kq_m(its:ite,kts:kte,j) > 1000.0) Kq_m(its:ite,kts:kte,j) = 1000.0

            
            call pbl_diffusion(qv, th, cloud, ice, qrain, qsnow, qv_flux, th_flux, &
                               rho, adv_dz, dz_mass_i, its, ite, kts, kte, j, jaco, jaco_w)

        enddo
        !call io_write("shear_m.nc", "shear_m", shear_m(:,:,:) )
        !call io_write("rig_m.nc", "rig_m", rig_m(:,:,:) )
        !call io_write("l_m.nc", "l_m", l_m(:,:,:) )
        !call io_write("stability_h.nc", "stability_h", stability_h(:,:,:) )
        !call io_write("stability_m.nc", "stability_m", stability_m(:,:,:) )
        !call io_write("BVF.nc", "BVF", BVF(:,:,:) )
        !call io_write("Kq_m.nc", "Kq_m", Kq_m(:,:,:) )
    end subroutine diagnostic_pbl

    subroutine diffuse_variable(q, rho, rho_stag, dz, dz_mass_i, its, ite, kts, kte, j, jaco, jaco_w, bot_flux)
        real,   intent(inout),  dimension(ims:ime,kms:kme,jms:jme) :: q
        real,   intent(in),     dimension(ims:ime,kms:kme,jms:jme) :: jaco, rho, dz, jaco_w
        real,   intent(in),     dimension(ims:ime,kms:kme-1,jms:jme) :: dz_mass_i
        real,   intent(in),     dimension(ims:ime,kms:kme) :: rho_stag
        real,   intent(in),  optional, dimension(ims:ime) :: bot_flux


        integer,intent(in) :: its, ite, kts, kte, j

        real, dimension(its:ite,kts:kte+1) :: fluxes
        integer :: k
        
        fluxes(its:ite,kts) = 0
        if (present(bot_flux)) fluxes(its:ite,kts) = bot_flux(its:ite)

        do k = kts+1, kte
            ! Eventually this should be made into an implicit solution to avoid substepping
            ! if gradient is downward (qv[z+1]>qv[z]) then flux is negative
            fluxes(its:ite,k) = -Kq_m(its:ite, k-1, j)*rho_stag(its:ite, k-1)*(q(its:ite, k, j) - q(its:ite, k-1, j)) / &
                                    (dz_mass_i(its:ite, k-1, j)*jaco_w(its:ite, k-1, j))
        enddo
        fluxes(its:ite,kte+1) = 0 !fluxes(its:ite,kte)

        q(its:ite, kts:kte,j) = q(its:ite, kts:kte, j) - (fluxes(its:ite,kts+1:kte+1) - fluxes(its:ite,kts:kte)) / &
                                    (rho(its:ite,kts:kte,j)*jaco(its:ite,kts:kte,j)*dz(its:ite,kts:kte,j))  

        
    end subroutine diffuse_variable

    subroutine pbl_diffusion(qv, th, cloud, ice, qrain, qsnow, qv_flux, th_flux, rho, dz, dz_mass_i, its, ite, kts, kte, j, jaco, jaco_w)
        real, intent(inout), dimension(ims:ime, kms:kme, jms:jme) :: th            ! potential temperature [K]
        real, intent(inout), dimension(ims:ime, kms:kme, jms:jme) :: qv            ! water vapor mixing ratio [kg/kg]
        real, intent(inout), dimension(ims:ime, kms:kme, jms:jme) :: cloud         ! cloud water mixing ratio [kg/kg]
        real, intent(inout), dimension(ims:ime, kms:kme, jms:jme) :: ice           ! cloud ice mixing ratio [kg/kg]
        real, intent(inout), dimension(ims:ime, kms:kme, jms:jme) :: qrain         ! rain water mixing ratio [kg/kg]
        real, intent(inout), dimension(ims:ime, kms:kme, jms:jme) :: qsnow         ! snow mixing ratio [kg/kg]
        real, intent(in),    dimension(ims:ime) :: qv_flux         ! snow mixing ratio [kg/kg]
        real, intent(in),    dimension(ims:ime) :: th_flux         ! snow mixing ratio [kg/kg]
        real, intent(in),    dimension(ims:ime, kms:kme, jms:jme) :: rho           ! air density [kg / m^3]
        real, intent(in),    dimension(ims:ime, kms:kme, jms:jme) :: dz            ! model level thickness
        real, intent(in),    dimension(ims:ime, kms:kme-1, jms:jme) :: dz_mass_i   ! model level thickness from mass-point to mass-point [m]
        real, intent(in),    dimension(ims:ime, kms:kme, jms:jme) :: jaco          ! mass-centered jacobian 
        real, intent(in),    dimension(ims:ime, kms:kme, jms:jme) :: jaco_w        ! w-staggered jacobian 

        integer,intent(in) :: its, ite, kts, kte, j

        ! locals
        integer :: k, nsubsteps, t
        real, dimension(ims:ime, kms:kme) :: rho_stag

        do k = kts, kte
            if (k < kte) then
                rho_stag(ims:ime, k) = (rho(ims:ime, k, j) + rho(ims:ime, k+1, j)) / 2
            else
                rho_stag(ims:ime, k) = rho(ims:ime, k, j)
            endif
        enddo

        !where((Kq_m(its:ite, kts:kte,j)) > N_substeps*dz(its:ite, kts:kte,j))   &
        !      Kq_m(its:ite, kts:kte, j) = dz(its:ite, kts:kte, j) * N_substeps

        !nsubsteps = ceiling( 2 * maxval(Kq_m(its:ite, kts:kte, j) / dz(its:ite, kts:kte, j)))
        !Kq_m(its:ite, kts:kte, j) = Kq_m(its:ite, kts:kte, j) / nsubsteps
        do t = 1, 2
            ! First water vapor
            call diffuse_variable(qv, rho, rho_stag, dz, dz_mass_i, its, ite, kts, kte, j, jaco, jaco_w, bot_flux=qv_flux)
            ! and cloud water
            call diffuse_variable(cloud, rho, rho_stag, dz, dz_mass_i, its, ite, kts, kte, j, jaco, jaco_w)
            ! and cloud ice
            call diffuse_variable(ice, rho, rho_stag, dz, dz_mass_i, its, ite, kts, kte, j, jaco, jaco_w)
            ! and snow
            call diffuse_variable(qsnow, rho, rho_stag, dz, dz_mass_i, its, ite, kts, kte, j, jaco, jaco_w)
            ! and rain
            call diffuse_variable(qrain, rho, rho_stag, dz, dz_mass_i, its, ite, kts, kte, j, jaco, jaco_w)
            ! ditto for potential temperature
            call diffuse_variable(th, rho, rho_stag, dz, dz_mass_i, its, ite, kts, kte, j, jaco, jaco_w, bot_flux=th_flux)
            ! don't bother with graupel assuming they are falling fast *enough* not entirely fair...
        enddo
    end subroutine pbl_diffusion

    subroutine calc_shear_sq(u_mass, v_mass, dz_mass_i,kts,kte)
        real, intent(in),    dimension(ims:ime, kms:kme, jms:jme) :: u_mass     ! east-west wind speed [m/s]
        real, intent(in),    dimension(ims:ime, kms:kme, jms:jme) :: v_mass     ! north south wind speed [m/s]
        real, intent(in),    dimension(ims:ime, kms:kme-1, jms:jme) :: dz_mass_i  ! model level thickness [m]
        integer,intent(in) :: kts, kte

        shear_m(:, kts, :) = ((-3*u_mass(:, kts, :) + 4*u_mass(:, kts+1, :) - u_mass(:, kts+2, :)) / (dz_mass_i(:, kts+1, :) + dz_mass_i(:, kts,:)))**2  + &
                             ((-3*v_mass(:, kts, :) + 4*v_mass(:, kts+1, :) - v_mass(:, kts+2, :)) / (dz_mass_i(:, kts+1, :) + dz_mass_i(:, kts,:)))**2 
        shear_m(:, kts+1:kte-1, :) = ((u_mass(:, kts+2:kte, :) - u_mass(:, kts:kte-2, :)) / (dz_mass_i(:, kts+1:kte-1, :) + dz_mass_i(:, kts:kte-2,:)))**2  + &
                                     ((v_mass(:, kts+2:kte, :) - v_mass(:, kts:kte-2, :)) / (dz_mass_i(:, kts+1:kte-1, :) + dz_mass_i(:, kts:kte-2,:)))**2  

        shear_m(:, kte, :) = ((3*u_mass(:, kte, :) - 4*u_mass(:, kte-1, :) + u_mass(:, kte-2, :)) / (dz_mass_i(:, kte-1, :) + dz_mass_i(:, kte-2,:)))**2  + &
                             ((3*v_mass(:, kte, :) - 4*v_mass(:, kte-1, :) + v_mass(:, kte-2, :)) / (dz_mass_i(:, kte-1, :) + dz_mass_i(:, kte-2,:)))**2 
                             
        where(shear_m<1e-5) shear_m = 1e-5
    end subroutine calc_shear_sq

!   calculate the vertical gradient in virtual potential temperature
    subroutine calc_BVF(th, qv, dz_mass_i, cloud, ice, qrain, qsnow, kts, kte)
        real, intent(inout), dimension(ims:ime, kms:kme, jms:jme) :: th            ! potential temperature [K]
        real, intent(inout), dimension(ims:ime, kms:kme, jms:jme) :: qv            ! water vapor mixing ratio [kg/kg]
        real, intent(inout), dimension(ims:ime, kms:kme, jms:jme) :: cloud         ! cloud water mixing ratio [kg/kg]
        real, intent(inout), dimension(ims:ime, kms:kme, jms:jme) :: ice           ! cloud ice mixing ratio [kg/kg]
        real, intent(inout), dimension(ims:ime, kms:kme, jms:jme) :: qrain         ! rain water mixing ratio [kg/kg]
        real, intent(inout), dimension(ims:ime, kms:kme, jms:jme) :: qsnow         ! snow mixing ratio [kg/kg]
        real, intent(in),    dimension(ims:ime, kms:kme-1, jms:jme) :: dz_mass_i     ! model level thickness [m]
        integer,intent(in) :: kts, kte

        real, allocatable, dimension(:,:,:) :: vth                        ! virtual potential temperature
        real, allocatable, dimension(:,:,:) :: virt_pot_temp_zgradient_m  ! virtual potential z-gradient

        allocate(vth(ims:ime, kms:kme, jms:jme))
        allocate(virt_pot_temp_zgradient_m(ims:ime, kms:kme, jms:jme))
        ! first calculate the virtual potential temperature
        ! vth=th*(1+0.61*qv-(qc+qi+qr+qs))
        vth = th * (1 + 0.61 * qv - (cloud + ice + qrain + qsnow))
                             
                                             
        virt_pot_temp_zgradient_m(:, kts,:) = (-3*vth(:, kts,:) + 4*vth(:, kts+1,:) - vth(:, kts+2,:)) / (dz_mass_i(:,kts,:)+dz_mass_i(:,kts+1,:))
        virt_pot_temp_zgradient_m(:, kts+1:kte-1,:) = (vth(:, kts+2:kte,:) - vth(:, kts:kte-2,:)) / (dz_mass_i(:,kts+1:kte-1,:)+dz_mass_i(:,kts:kte-2,:))                 
        virt_pot_temp_zgradient_m(:, kte,:) = ( 3*vth(:, kte,:) - 4*vth(:, kte-1,:) + vth(:, kte-2,:)) / (dz_mass_i(:,kte-1,:)+dz_mass_i(:,kte-2,:))
                                             
        BVF = (gravity/vth)*virt_pot_temp_zgradient_m

    end subroutine calc_BVF

    ! calculate the stability function based on HP96
    subroutine calc_pbl_stability_function()
        real, allocatable, dimension(:,:,:) :: gamma 

        allocate(gamma(ims:ime, kms:kme, jms:jme))

        rig_m = BVF/shear_m
        
        where(rig_m > 0)  ri_flux = 0.8333*(rig_m + 0.2805 - sqrt(max(0.0, (rig_m**2 - 0.1122*rig_m + 0.2805**2) )))
        where(rig_m <= 0) ri_flux = 1.285 *(rig_m + 0.2305 - sqrt(max(0.0, (rig_m**2 + 0.1023*rig_m + 0.2305**2) )))

        gamma = ri_flux/(1-ri_flux)
        
        where(rig_m > 0)  stability_h = (1-2.5648*gamma)/(1-1.1388*gamma)
        where(rig_m <= 0) stability_h = (1-3.337*gamma) /(1-0.688*gamma)
        
        where(rig_m > 0)  stability_m = 1-3.7000*gamma/stability_h
        where(rig_m <= 0) stability_m = 1-4.025*gamma /stability_h

    end subroutine calc_pbl_stability_function



! memory allocation and deallocation
! Can/Should also add parameter definition from options%pbl (or something like that)
    subroutine init_diagnostic_pbl(domain,options)
        type(domain_t), intent(in) :: domain
        type(options_t),intent(in) :: options

        ! module level variables assumed constant throughout the run...
        ims = domain%ims
        ime = domain%ime
        jms = domain%jms
        jme = domain%jme
        kms = domain%kms
        kme = domain%kme
        ids = domain%ids
        ide = domain%ide
        jds = domain%jds
        jde = domain%jde
        kds = domain%kds
        kde = domain%kde

        allocate(BVF                      (ims:ime,kms:kme,jms:jme))
        allocate(rig_m                    (ims:ime,kms:kme,jms:jme))
        allocate(ri_flux                  (ims:ime,kms:kme,jms:jme))
        allocate(stability_m              (ims:ime,kms:kme,jms:jme))
        allocate(stability_h              (ims:ime,kms:kme,jms:jme))
        allocate(shear_m                  (ims:ime,kms:kme,jms:jme))
        allocate(Kq_m                     (ims:ime,kms:kme,jms:jme))
        allocate(l_m                      (ims:ime,kms:kme,jms:jme))
        allocate(lastqv_m                 (ims:ime,kms:kme,jms:jme))
    end subroutine init_diagnostic_pbl

!   deallocate memory if requested
    subroutine finalize_diagnostic_pbl()
        if (allocated(BVF)) then
            deallocate(BVF)
        endif
        if (allocated(rig_m)) then
            deallocate(rig_m)
        endif
        if (allocated(ri_flux)) then
            deallocate(ri_flux)
        endif
        if (allocated(stability_m)) then
            deallocate(stability_m)
        endif
        if (allocated(stability_h)) then
            deallocate(stability_h)
        endif
        if (allocated(shear_m)) then
            deallocate(shear_m)
        endif
        if (allocated(Kq_m)) then
            deallocate(Kq_m)
        endif
        if (allocated(l_m)) then
            deallocate(l_m)
        endif
        if (allocated(lastqv_m)) then
            deallocate(lastqv_m)
        endif
    end subroutine finalize_diagnostic_pbl
end module pbl_diagnostic
