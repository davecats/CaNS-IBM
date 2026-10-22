! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
module mod_utils
  implicit none
  private
  public bulk_mean,f_sizeof,swap,calc_mean_flow,calc_mean_flow_easy,chk_ibm
  public bulk_mean_ibm,bulk_forcing_ibm
#if defined(_OPENACC)
  public device_memory_footprint
#endif
contains
  subroutine bulk_mean(n,grid_vol_ratio,p,mean)
    !
    ! compute the mean value of an observable over the entire domain
    !
    use mpi
    use mod_types
    implicit none
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(0:) :: grid_vol_ratio
    real(rp), intent(in), dimension(0:,0:,0:) :: p
    real(rp), intent(out) :: mean
    integer :: i,j,k
    integer :: ierr
    mean = 0.
    !$acc data copy(mean) async(1)
    !$acc parallel loop collapse(3) default(present) reduction(+:mean) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)  REDUCTION(+:mean)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          mean = mean + p(i,j,k)*grid_vol_ratio(k)
        end do
      end do
    end do
    !$acc end data
    !$acc wait(1)
    call MPI_ALLREDUCE(MPI_IN_PLACE,mean,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
  end subroutine bulk_mean
  subroutine bulk_mean_ibm(n,grid_vol_ratio,mask,a,p,is_intrinsic,mean,gain)
    !
    ! IBM-aware bulk mean of p over the non-solid points (mask = .false.), and the gain,
    ! i.e. the mean response to a unit impulse added as 1/a (see `bulk_forcing_ibm`);
    ! both are normalized by the box volume (superficial) or the fluid volume (intrinsic)
    !
    use mpi
    use mod_types
    implicit none
    integer , intent(in), dimension(3)        :: n
    real(rp), intent(in), dimension(0:)       :: grid_vol_ratio
    logical , intent(in), dimension(0:,0:,0:) :: mask
    real(rp), intent(in), dimension(0:,0:,0:) :: a,p
    logical , intent(in)                      :: is_intrinsic
    real(rp), intent(out)                     :: mean,gain
    real(rp) :: vol,buf(3)
    integer  :: i,j,k,ierr
    mean = 0.; gain = 0.; vol = 0.
    !$acc data copy(mean,gain,vol) async(1)
    !$acc parallel loop collapse(3) default(present) reduction(+:mean,gain,vol) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)  REDUCTION(+:mean,gain,vol)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          if(.not.mask(i,j,k)) then
            mean = mean + p(i,j,k)*grid_vol_ratio(k)
            gain = gain + grid_vol_ratio(k)/a(i,j,k)
            vol  = vol  + grid_vol_ratio(k)
          end if
        end do
      end do
    end do
    !$acc end data
    !$acc wait(1)
    buf(:) = [mean,gain,vol]
    call MPI_ALLREDUCE(MPI_IN_PLACE,buf,3,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    mean = buf(1); gain = buf(2); vol = buf(3)
    if(is_intrinsic) then
      mean = mean/vol
      gain = gain/vol
    end if
  end subroutine bulk_mean_ibm
  subroutine bulk_forcing_ibm(n,mask,a,ff,p)
    !
    ! adds the bulk forcing impulse ff to p as the exact IB scheme would integrate it,
    ! ff/a, at the non-solid points only
    !
    use mod_types
    implicit none
    integer , intent(in   ), dimension(3)        :: n
    logical , intent(in   ), dimension(0:,0:,0:) :: mask
    real(rp), intent(in   ), dimension(0:,0:,0:) :: a
    real(rp), intent(in   )                      :: ff
    real(rp), intent(inout), dimension(0:,0:,0:) :: p
    integer :: i,j,k
    !$acc parallel loop collapse(3) default(present) async(1)
    !$OMP parallel do   collapse(3) DEFAULT(shared)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          if(.not.mask(i,j,k)) p(i,j,k) = p(i,j,k) + ff/a(i,j,k)
        end do
      end do
    end do
  end subroutine bulk_forcing_ibm
  subroutine calc_mean_flow(grid_vol_ratio,mask_id,field_id,mean,n)
    !
    ! compute the mean value of an observable over the entire domain
    !
    use mpi
    use mod_types
    implicit none
    integer , intent(in), dimension(3)        :: n
    logical,intent(in)                        :: mask_id(0:,0:,0:) 
    real(rp), intent(in), dimension(0:)       :: grid_vol_ratio
    real(rp), intent(in), dimension(0:,0:,0:) :: field_id
    real(rp), intent(out)                     :: mean
    real(rp)                                  :: vol  
    integer                                   :: i,j,k
    integer                                   :: ierr
    mean = 0.
    vol  = 0._rp
    !$acc data copy(mean,vol) async(1)
    !$acc parallel loop collapse(3) default(present) reduction(+:mean,vol) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)  REDUCTION(+:mean,vol)
    do k = 1,n(3)
      do j = 1,n(2)
        do i = 1,n(1)

          if(.not.mask_id(i,j,k))then
            mean=mean+field_id(i,j,k)*grid_vol_ratio(k)
            vol=vol+grid_vol_ratio(k) 
          endif
        end do
      end do
    end do
    !$acc end data
    !$acc wait(1)
    call MPI_ALLREDUCE(MPI_IN_PLACE,mean,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,vol,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    ! we do this bcs now whole domain isnt covered w/ fluid cells. 
    ! some of them are solid cells and we need take them out of calculation
    ! since CaNS grid_vol_ratio calcultes over the whole domain which includes
    ! the solids cells to therefore we find the fluid cells too. 
    mean=mean/vol
  end subroutine calc_mean_flow
  subroutine calc_mean_flow_easy(mask_id,field_id,q,n,dl,dir)
    use mpi
    use mod_types
    implicit none
    logical,intent(in)                        :: mask_id(0:,0:,0:) 
    integer , intent(in)                      :: dir
    real(rp), intent(in), dimension(3)        :: dl
    integer , intent(in), dimension(3)        :: n
    real(rp), intent(in), dimension(0:,0:,0:) :: field_id
    real(rp), intent(out)                     :: q
    integer                                   :: i,j,k
    real(rp)                                  :: aface
    integer                                   :: ierr
    select case(dir)
      case(1)
        aface=dl(2)*dl(3)
      case(2)
        aface=dl(3)*dl(1)
      case(3)
        aface=dl(2)*dl(1)
      case default
        print*, "invalid dir"
    end select
    q = 0._rp
    !$acc data copy(q) async(1)
    !$acc parallel loop collapse(3) default(present) reduction(+:q) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)  REDUCTION(+:q)
    do k = 1,n(3)
      do j = 1,n(2)
        do i = 1,n(1)
          if(.not.mask_id(i,j,k))then
            q = q + field_id(i,j,k)*aface
          endif
        end do
      end do
    end do
    !$acc end data
    !$acc wait(1)
    call MPI_ALLREDUCE(MPI_IN_PLACE,q,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    ! we do this bcs now whole domain isnt covered w/ fluid cells. 
    ! some of them are solid cells and we need take them out of calculation
    ! since CaNS grid_vol_ratio calcultes over the whole domain which includes
    ! the solids cells to therefore we find the fluid cells too. 
    q=real(q/n(dir),kind=rp)
  end subroutine calc_mean_flow_easy
  subroutine chk_ibm(n,l,dli,dzfi,grid_vol_ratio_f,mask_u,mask_v,mask_w,u,v,w, &
                     velmax_solid,divtot_m,divmax_m,ubulk_fl)
    !
    ! IBM diagnostics:
    !   velmax_solid(1:3): max|u|, max|v|, max|w| over solid faces
    !   divtot_m, divmax_m: divergence norms (as in chkdiv) with the solid faces zeroed
    !   ubulk_fl: mean of u over the fluid u-faces only
    !
    use mpi
    use mod_types
    implicit none
    integer , intent(in), dimension(3)        :: n
    real(rp), intent(in), dimension(3)        :: l,dli
    real(rp), intent(in), dimension(0:)       :: dzfi,grid_vol_ratio_f
    logical , intent(in), dimension(0:,0:,0:) :: mask_u,mask_v,mask_w
    real(rp), intent(in), dimension(0:,0:,0:) :: u,v,w
    real(rp), intent(out)                     :: velmax_solid(3),divtot_m,divmax_m,ubulk_fl
    real(rp) :: umax,vmax,wmax,div,vol,dxi,dyi
    integer  :: i,j,k,ierr
    dxi = dli(1)
    dyi = dli(2)
    umax = 0.; vmax = 0.; wmax = 0.
    divtot_m = 0.; divmax_m = 0.
    ubulk_fl = 0.; vol = 0.
    !$acc data copy(umax,vmax,wmax,divtot_m,divmax_m,ubulk_fl,vol) async(1)
    !$acc parallel loop collapse(3) default(present) private(div) &
    !$acc reduction(max:umax,vmax,wmax,divmax_m) reduction(+:divtot_m,ubulk_fl,vol) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)  PRIVATE(div) &
    !$OMP REDUCTION(max:umax,vmax,wmax,divmax_m) REDUCTION(+:divtot_m,ubulk_fl,vol)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          if(mask_u(i,j,k)) then
            umax = max(umax,abs(u(i,j,k)))
          else
            ubulk_fl = ubulk_fl + u(i,j,k)*grid_vol_ratio_f(k)
            vol      = vol      +          grid_vol_ratio_f(k)
          end if
          if(mask_v(i,j,k)) vmax = max(vmax,abs(v(i,j,k)))
          if(mask_w(i,j,k)) wmax = max(wmax,abs(w(i,j,k)))
          div = (merge(0._rp,w(i,j,k),mask_w(i,j,k))-merge(0._rp,w(i,j,k-1),mask_w(i,j,k-1)))*dzfi(k) + &
                (merge(0._rp,v(i,j,k),mask_v(i,j,k))-merge(0._rp,v(i,j-1,k),mask_v(i,j-1,k)))*dyi     + &
                (merge(0._rp,u(i,j,k),mask_u(i,j,k))-merge(0._rp,u(i-1,j,k),mask_u(i-1,j,k)))*dxi
          divmax_m = max(divmax_m,abs(div))
          divtot_m = divtot_m + abs(div)/(dxi*dyi*dzfi(k))
        end do
      end do
    end do
    !$acc end data
    !$acc wait(1)
    velmax_solid(:) = [umax,vmax,wmax]
    call MPI_ALLREDUCE(MPI_IN_PLACE,velmax_solid,3,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,divmax_m,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,divtot_m,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,ubulk_fl,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,vol,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    divtot_m = divtot_m/product(l(:))
    ubulk_fl = ubulk_fl/vol
  end subroutine chk_ibm
  pure integer function f_sizeof(val) result(isize)
    !
    ! returns storage size of the scalar argument val in bytes
    !
    implicit none
    class(*), intent(in) :: val
    isize = storage_size(val)/8
  end function f_sizeof
  subroutine swap(arr1,arr2)
    use mod_types, only: rp
    implicit none
    real(rp), intent(inout), pointer, contiguous, dimension(:,:,:) :: arr1,arr2
    real(rp),                pointer, contiguous, dimension(:,:,:) :: tmp
    tmp  => arr1
    arr1 => arr2
    arr2 => tmp
  end subroutine swap
#if defined(_OPENACC)
  function device_memory_footprint(n,n_z,nscal) result(itotal)
    !
    ! estimate GPU memory footprint, assuming one MPI task <-> one GPU
    !
    use mod_types, only: i8,rp
    use mod_param, only: is_impdiff,is_impdiff_1d,is_poisson_pcr_tdma
    integer, intent(in), dimension(3) :: n,n_z
    integer, intent(in) :: nscal
    integer :: nh(3)
    integer(i8) :: itotal,itemp,rp_size
    rp_size = f_sizeof(1._rp)
    itotal = 0
    !
    ! 1. 'main' arrays: u,v,w,p,pp,scalars
    !
    nh(:) = 1
    itotal = itotal + product(int(n(:),i8)+2*int(nh(:),i8))*rp_size*(5+nscal)
    !
    ! 2. grids arrays: zc,zf,dzc,dzf,dzci,dzfi,grid_vol_ratio_c,grid_vol_ratio_f (tiny footprint)
    !
    nh(:) = 1
    itotal = itotal + (n(3)+2*nh(3))*rp_size*8
    !
    ! 3. solver eigenvalues and Gauss elimination coefficient arrays (small footprint)
    !    rhs?%[x,y,z] arrays, lambdaxy? arrays, and a?,b?,c? arrays
    !
    block
      integer(i8) :: itemp1,itemp1_(3),itemp2,itemp3
      itemp1_(:) = [int(n_z(2),i8)*int(n_z(3),i8)*2, &
                    int(n_z(1),i8)*int(n_z(3),i8)*2, &
                    int(n_z(1),i8)*int(n_z(2),i8)*2]
      itemp1 = sum(itemp1_(:))   ! rhs
      itemp2 = product(int(n_z(1:2),i8)) ! lambdaxy
      itemp3 = int(n_z(3),i8)*3          ! a,b,c
      if(.not.is_impdiff) then
        !
        ! rhsbp, lambdaxyp, ap,bp,cp
        !
        itotal = itotal + itemp1*rp_size                        + itemp2*rp_size           + itemp3*rp_size
      else if(is_impdiff_1d) then
        !
        ! rhsbp,rhsb[u,v,w,scalars,buf]%z, lambdaxyp, a?,b?,c? [p,u,v,w,scalars,buf]
        !
        itotal = itotal + (itemp1+itemp1_(3)*(4+nscal))*rp_size + itemp2*rp_size           + itemp3*rp_size*(5+nscal)
      else
        !
        ! rhsbp,rhsb[u,v,w,scalars,buf]%[x,y,z], lambdaxy[p,u,v,w,scalars], (a?,b?,c?)[p,u,v,w,scalars,buf]
        !
        itotal = itotal + itemp1*rp_size*(1+4+nscal)            + itemp2*rp_size*(5+nscal) + itemp3*rp_size*(5+nscal)
      end if
    end block
    !
    ! 4. prediction velocity arrays arrays d[u,v,w]dtrk_t, d[u,v,w]dtrko_t + scalars equivalent
    !
    itemp  = product(int(n(:),i8))*rp_size
    itotal = itotal + itemp*(2*(3+nscal))
    if(is_impdiff) then
      itotal = itotal + itemp*(1*(3+nscal))
    end if
    !
    ! 5. transpose & FFT buffer arrays, halo buffer arrays, and solver arrays
    !    taken directly from `mod_common_cudecomp`
    !
    block
      use mod_common_cudecomp, only: work,work_halo,work_ptdma,solver_buf_0,solver_buf_1,pz_aux_1,pz_aux_2
      itemp = size(work        ,kind=i8) + size(work_halo   ,kind=i8) + &
              size(solver_buf_0,kind=i8) + size(solver_buf_1,kind=i8)
      if(allocated(pz_aux_1)) itemp = itemp + size(pz_aux_1,kind=i8)
      if(allocated(pz_aux_2)) itemp = itemp + size(pz_aux_2,kind=i8)
      if(is_poisson_pcr_tdma) itemp = itemp + size(work_ptdma,kind=i8)
      itotal = itotal + itemp*rp_size
    end block
  end function device_memory_footprint
#endif
end module mod_utils
