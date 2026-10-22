! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
module mod_correc
  use mod_types
  implicit none
  private
  public correc,correc_ibm
  contains
  subroutine correc(n,dli,dzci,dt,p,u,v,w)
    !
    ! corrects the velocity so that it is divergence free
    !
    implicit none
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(3 ) :: dli
    real(rp), intent(in), dimension(0:) :: dzci
    real(rp), intent(in) :: dt
    real(rp), intent(in   ), dimension(0:,0:,0:) :: p
    real(rp), intent(inout), dimension(0:,0:,0:) :: u,v,w
    real(rp) :: factori,factorj
    !real(rp), dimension(0:n(3)+1) :: factork
    integer :: i,j,k
    !
    !factor = rkcoeffab(rkiter)*dt
    !
    factori = dt*dli(1)
    factorj = dt*dli(2)
    !factork = dt*dzci!dli(3)
    !$acc parallel loop collapse(3) default(present) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
    do k=0,n(3)+1
      do j=0,n(2)+1
        do i=0,n(1)
          u(i,j,k) = u(i,j,k) - factori*(   p(i+1,j,k)-p(i,j,k))
        end do
      end do
    end do
    !$acc parallel loop collapse(3) default(present) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
    do k=0,n(3)+1
      do j=0,n(2)
        do i=0,n(1)+1
          v(i,j,k) = v(i,j,k) - factorj*(   p(i,j+1,k)-p(i,j,k))
        end do
      end do
    end do
    !$acc parallel loop collapse(3) default(present) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
    do k=0,n(3)
      do j=0,n(2)+1
        do i=0,n(1)+1
          w(i,j,k) = w(i,j,k) - dt*dzci(k)*(p(i,j,k+1)-p(i,j,k))
        end do
      end do
    end do
  end subroutine correc
  subroutine correc_ibm(n,dli,dzci,dt,p,mask_u,mask_v,mask_w,a_u,a_v,a_w,u,v,w)
    !
    ! IB-aware velocity correction (`ibm_proj = 'damped'`): as `correc`, but the correction
    ! is divided by the IB coefficient a at the non-solid points (the exact IB scheme
    ! integrates the pressure-gradient impulse as 1/a), and the solid points are not corrected;
    ! the velocity is then divergence free only up to a residual in the near-wall cells
    !
    implicit none
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(3 ) :: dli
    real(rp), intent(in), dimension(0:) :: dzci
    real(rp), intent(in) :: dt
    real(rp), intent(in   ), dimension(0:,0:,0:) :: p
    logical , intent(in   ), dimension(0:,0:,0:) :: mask_u,mask_v,mask_w
    real(rp), intent(in   ), dimension(0:,0:,0:) :: a_u,a_v,a_w
    real(rp), intent(inout), dimension(0:,0:,0:) :: u,v,w
    real(rp) :: factori,factorj
    integer :: i,j,k
    !
    factori = dt*dli(1)
    factorj = dt*dli(2)
    !$acc parallel loop collapse(3) default(present) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
    do k=0,n(3)+1
      do j=0,n(2)+1
        do i=0,n(1)
          if(.not.mask_u(i,j,k)) u(i,j,k) = u(i,j,k) - factori*(   p(i+1,j,k)-p(i,j,k))/a_u(i,j,k)
        end do
      end do
    end do
    !$acc parallel loop collapse(3) default(present) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
    do k=0,n(3)+1
      do j=0,n(2)
        do i=0,n(1)+1
          if(.not.mask_v(i,j,k)) v(i,j,k) = v(i,j,k) - factorj*(   p(i,j+1,k)-p(i,j,k))/a_v(i,j,k)
        end do
      end do
    end do
    !$acc parallel loop collapse(3) default(present) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)
    do k=0,n(3)
      do j=0,n(2)+1
        do i=0,n(1)+1
          if(.not.mask_w(i,j,k)) w(i,j,k) = w(i,j,k) - dt*dzci(k)*(p(i,j,k+1)-p(i,j,k))/a_w(i,j,k)
        end do
      end do
    end do
  end subroutine correc_ibm
end module mod_correc
