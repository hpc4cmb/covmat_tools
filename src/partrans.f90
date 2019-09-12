! Parallel transport module by Hans Kristian Eriksen, November 2013
! Adaptations for covmat_tools by RK

module partrans
  use healpix_types
  use pix_tools
  implicit none  

contains


  ! Parallel transport polarization angles from a given map to centers
  ! of a lower resolution map at nside_out. 
  ! Note: Maps must be in nested format
  subroutine qu_transport_map(nside_out, map)
    implicit none

    integer(i4b),                   intent(in)    :: nside_out
    real(dp),     dimension(0:,1:), intent(inout) :: map

    integer(i4b) :: i, j, q, ierr
    real(dp)     :: cos2psi, sin2psi, m(2)
    integer(i4b), save :: nside, npix=-1, nmaps=-1
    real(dp), allocatable, dimension(:,:,:), save :: vec

    if (npix /= size(map,1) .or. nmaps /= size(map,2)) then
       if (allocated(vec)) deallocate(vec)
    end if

    if (nside_out < 1) return

    ! Precompute pixel vectors
    if (.not. allocated(vec)) then
       npix  = size(map,1)
       nmaps = size(map,2)
       nside = npix2nside(npix)
       q     = (nside/nside_out)**2
       if (nmaps /= 3) then
          write(*,*) 'partrans -- nmaps must be equal to 3'
          stop
       end if
       allocate(vec(0:npix-1,3,2), stat=ierr)
       if (ierr /= 0) stop 'No room for parallel transport workspace'
       !$OMP PARALLEL DO SCHEDULE(STATIC,1000) PRIVATE(i)
       do i = 0, npix-1
          call pix2vec_nest(nside_out, i/q, vec(i,:,1)) ! Low-res pixel centers
          call pix2vec_nest(nside,     i,   vec(i,:,2)) ! High-res pixel centers
       end do
       !$OMP END PARALLEL DO
    end if

    if (nside == nside_out) return

    ! Perform the parallel transport
    !$OMP PARALLEL DO SCHEDULE(STATIC,1000) PRIVATE(i,m,cos2psi,sin2psi)
    do i = 0, npix-1
       call qu_transport_rot(vec(i,:,1), vec(i,:,2), cos2psi, sin2psi)
       m        = map(i,2:3)
       map(i,2) =  cos2psi*m(1) - sin2psi*m(2)
       map(i,3) =  sin2psi*m(1) + cos2psi*m(2)
    end do
    !$OMP END PARALLEL DO

  end subroutine qu_transport_map

  ! Returns cos(2*psi) and sin(2*psi) of the rotation induced by
  ! parallel transport from vec1 to vec2.
  subroutine qu_transport_rot(vec1, vec2, cos2psi, sin2psi)
    implicit none

    real(dp), dimension(3), intent(in)  :: vec1, vec2
    real(dp),               intent(out) :: cos2psi, sin2psi

    integer(i4b) :: i, j, nfield
    real(dp) :: len_u, len_v, z, cos_theta, sgn, c1, s1, c2, s2
    real(dp), dimension(3) :: u, v

    z = sum(vec1*vec2)
    if (abs(z) >= 1.d0-1.d-8) then
       cos2psi = 1; sin2psi = 0
       return
    end if

    sgn    = 1.d0
    if (vec1(1)*vec2(2)-vec1(2)*vec2(1) < 0.d0) sgn = -1.d0

    ! Rotation from vec1 to vec 2
    u         = vec1(3) * vec1 
    u(3)      = u(3) - 1.d0
    v         = vec2 - z * vec1
    len_u     = sqrt(sum(u*u))
    len_v     = sqrt(sum(v*v))
    ! Local angle from vertical
    cos_theta = max(min((z * vec1(3) - vec2(3)) / (len_u*len_v), 1.d0), -1.d0)
    ! Double it and calculate cos and sin
    c1 = 2*cos_theta**2-1
    s1 = 2*sgn*sqrt(1-cos_theta**2)*cos_theta

    ! Rotation from vec2 to vec 1; sgn is opposite from 1->2
    u          = vec2(3) * vec2
    u(3)       = u(3) - 1.d0
    v          = vec1 - z * vec2
    len_u      = sqrt(sum(u*u))
    len_v      = sqrt(sum(v*v))
    cos_theta  = max(min((z * vec2(3) - vec1(3)) / (len_u*len_v),1.d0),-1.d0)
    c2 =  2*cos_theta**2-1
    s2 = -2*sgn*sqrt(1-cos_theta**2)*cos_theta

    cos2psi = c1*c2+s1*s2
    sin2psi = c1*s2-s1*c2
  end subroutine qu_transport_rot


end module partrans
