!> @file
!! This file holds only the downgrade module.

!> @brief This module contains routines to downgrade Healpix maps.
!! @author Reijo Keskitalo
!! @date 2009-2010
!! @version 1.0

! Revisions:
! 2013-09-04 : Missing pixel value is now ALWAYS 0, not HPX_BADVAL.
!              Depending on the client code, it may need to replace the values
! 2013-11-22 : - Removed the unused Plm argument.
!              - Added support for user-supplied smoothing window.
!              - Added parallel transport of polarization vectors by Hans Kristian Eriksen

MODULE downgrade

  USE fitstools
  USE head_fits
  USE pix_tools
  USE healpix_types
  USE extension
  USE udgrade_nr
  USE alm_tools
  USE partrans ! Parallel transport of polarization vectors

  USE covmat_util, ONLY : tic, toc

  INTERFACE noise_weight_downgrade
     MODULE PROCEDURE &
          noise_weight_downgrade_nobs_preloaded, &
          noise_weight_downgrade_nobs_not_loaded
  END INTERFACE



CONTAINS



  !> Inverse noise-weight high resolution pixels for
  !! optimal low resolution map
  SUBROUTINE noise_weight_downgrade_nobs_not_loaded(nobs_file, map_in, map_out, &
       nside_in, nside_out, nmaps, rcond_limit_in, invert_in, diag_only_in)

    IMPLICIT NONE

    CHARACTER(len=*), INTENT(IN) :: nobs_file      !< fits file with (inverse) white noise covariance matrices
    REAL(dp), POINTER, INTENT(IN) :: map_in(:,:)   !< array containing the map to downgrade
    REAL(dp), POINTER, INTENT(INOUT) :: map_out(:,:) !< output map
    INTEGER(i4b), INTENT(IN) :: nside_in           !< input resolution
    INTEGER(i4b), INTENT(IN) :: nside_out          !< output resolution
    INTEGER(i4b), INTENT(IN) :: nmaps              !< number of Stokes components
    REAL(dp), OPTIONAL, INTENT(IN) :: rcond_limit_in !< threshold for excluding pixels from output map based on
    LOGICAL, OPTIONAL, INTENT(in) :: invert_in !< nobs have been inverted
    LOGICAL, OPTIONAL, INTENT(in) :: diag_only_in !< nobs have been inverted
    !! reciprocal condition number of the associated, combined matrix
    REAL(sp), POINTER, SAVE :: cc(:,:)=>null()     !< preloaded Nobs matrices

    INTEGER(i4b) :: ierr, ncols, nmaps_cc, ordering_cc, nside_cc
    INTEGER(i8b) :: npix_in, npix_cc
    REAL(dp) :: rcond_limit
    LOGICAL :: invert, diag_only
    real(sp), allocatable :: cctemp(:)

    IF (nmaps < 1) THEN
       IF (ASSOCIATED(cc)) DEALLOCATE(cc)
       RETURN
    END IF

    invert = .FALSE.
    IF (PRESENT(invert_in)) invert = invert_in

    diag_only = .FALSE.
    IF (PRESENT(diag_only_in)) diag_only = diag_only_in

    ncols = (nmaps * (nmaps + 1)) / 2

    npix_in = nside2npix(nside_in)

    rcond_limit = 0.0_dp
    IF (PRESENT(rcond_limit_in)) rcond_limit = rcond_limit_in

    IF (.NOT. ASSOCIATED(cc)) THEN
       ALLOCATE(cc(0:npix_in-1, ncols), stat=ierr)
       IF (ierr /= 0) STOP 'No room for cc'

       npix_cc = getsize_fits(nobs_file, nmaps=nmaps_cc, ordering=ordering_cc, nside=nside_cc)

       IF (nmaps_cc /= ncols) STOP 'ERROR: wrong number of columns in nobs'
       IF (nside_cc /= nside_in) STOP 'ERROR: wrong nside in nobs'

       CALL input_tod(nobs_file, cc, npix_in, nCols)
       SELECT CASE (ordering_cc)
       CASE (0)
          WRITE (*,*) 'Unknown ordering, assuming NESTED'
       CASE (1)
          WRITE (*,*) 'Converting matrix to NESTED'
          CALL convert_ring2nest(nside_in, cc)
       CASE (2)
          WRITE (*,*) 'Matrix is in NESTED scheme'
       END SELECT

       if (any(cc(:, ncols) < 0) .and. .not. any(cc(:, ncols-1) < 0)) then
          ! The QU and UU columns are swapped
          write(*,'(a)') 'WARNING: QU and UU columns are reversed in ' // trim(nobs_file)
          allocate(cctemp(npix_in), stat=ierr)
          if (ierr /= 0) stop 'No room to swap CC columns'
          cctemp = cc(:, ncols)
          cc(:, ncols) = cc(:, ncols-1)
          cc(:, ncols-1) = cctemp
          deallocate(cctemp)
       end if
    END IF

    CALL noise_weight_downgrade_nobs_preloaded(cc, map_in, map_out, &
         nside_in, nside_out, nmaps, rcond_limit, invert, diag_only)

  END SUBROUTINE noise_weight_downgrade_nobs_not_loaded



  !> Inverse noise-weight high resolution pixels for
  !! optimal low resolution map
  SUBROUTINE noise_weight_downgrade_nobs_preloaded(cc, map_in, map_out, &
       nside_in, nside_out, nmaps, rcond_limit_in, invert_in, diag_only_in)

    IMPLICIT NONE

    REAL(sp), POINTER, INTENT(IN) :: cc(:,:)       !< preloaded Nobs matrices
    REAL(dp), POINTER, INTENT(IN) :: map_in(:,:)   !< array containing the map to downgrade
    REAL(dp), POINTER, INTENT(INOUT) :: map_out(:,:) !< output map
    INTEGER(i4b), INTENT(IN) :: nside_in           !< input resolution
    INTEGER(i4b), INTENT(IN) :: nside_out          !< output resolution
    INTEGER(i4b), INTENT(IN) :: nmaps              !< number of Stokes components
    REAL(dp), OPTIONAL, INTENT(IN) :: rcond_limit_in !< threshold for excluding pixels from output map based on
    !! reciprocal condition number of the associated, combined matrix
    LOGICAL, OPTIONAL, INTENT(in) :: invert_in
    LOGICAL, OPTIONAL, INTENT(in) :: diag_only_in

    INTEGER(i4b) :: ncols, irow, icol, ierr, i
    INTEGER(i4b), PARAMETER :: workspacelen=100
    REAL(dp), ALLOCATABLE :: workspace(:)
    INTEGER(i8b) :: ipix_in, ipix_out, npix2one, npix_in, npix_out
    REAL(dp) :: rcond_limit
    REAL(dp), ALLOCATABLE :: eigenvals(:), eigenvecs(:,:)
    REAL(dp), ALLOCATABLE :: cc_one(:,:), cc_sum(:,:)
    LOGICAL :: eigeninvert, invert, diag_only
    REAL(dp), POINTER :: maptemp(:,:)

    IF (nmaps < 1) RETURN

    IF (nmaps == 3) THEN
       ALLOCATE(maptemp(0:12*nside_in**2-1, nmaps), stat=ierr)
       IF (ierr /= 0) STOP 'No room for parallel transport'
       maptemp = map_in
       CALL qu_transport_map(nside_out, maptemp)
    ELSE
       maptemp => map_in
    END IF

    eigeninvert = .FALSE.
    rcond_limit = 0
    IF (PRESENT(rcond_limit_in)) rcond_limit = rcond_limit_in
    IF (rcond_limit > 0.0) eigeninvert = .TRUE.

    invert = .FALSE.
    IF (PRESENT(invert_in)) invert = invert_in

    diag_only = .FALSE.
    IF (PRESENT(diag_only_in)) diag_only = diag_only_in

    ncols = (nmaps * (nmaps + 1)) / 2

    npix_in = nside2npix(nside_in)
    npix_out = nside2npix(nside_out)
    IF (npix_in < 0 .OR. npix_out < 0) THEN
       WRITE (*,'(" Illegal nside: ",i0,", ",i0)') nside_in, nside_out
       RETURN
    END IF
    npix2one = (nside_in / nside_out)**2

    IF (.NOT. ASSOCIATED(map_out)) THEN
       ALLOCATE(map_out(0:npix_out-1, nmaps), stat=ierr)
       IF (ierr /= 0) STOP 'Failed to allocate map_out'
    END IF
    map_out = 0.0_dp

    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP   SHARED(nmaps, npix_out, npix2one, cc, maptemp, map_out,&
    !$OMP          diag_only, eigeninvert, invert, rcond_limit)&
    !$OMP   PRIVATE(ipix_out, cc_sum, ipix_in, i, irow, icol, cc_one,&
    !$OMP           eigenvecs, eigenvals, workspace)

    ALLOCATE(cc_one(nmaps, nmaps), cc_sum(nmaps, nmaps))
    ALLOCATE(eigenvecs(nmaps, nmaps), eigenvals(nmaps))
    ALLOCATE(workspace(workspacelen))

    !$OMP DO SCHEDULE(STATIC, 4)
    DO ipix_out = 0, npix_out-1
       ! Sum the Nobs matrices of the subpixels, map is NESTED
       cc_sum = 0
       DO ipix_in = ipix_out*npix2one, (ipix_out+1)*npix2one-1
          i = 0
          IF (diag_only) cc_one = 0 ! only populate the diagonal
          DO irow = 1, nmaps
             DO icol = irow, nmaps
                i = i + 1
                IF (irow == icol) THEN
                   cc_one(irow,irow) = cc(ipix_in,i)
                ELSE IF (.NOT. diag_only) THEN
                   cc_one(irow,icol) = cc(ipix_in,i)
                   cc_one(icol,irow) = cc(ipix_in,i)
                END IF
             END DO
          END DO

          IF (invert) THEN
             IF (diag_only) THEN
                DO irow = 1, nmaps
                   cc_one(irow,irow) = 1 / cc_one(irow,irow)
                END DO
             ELSE IF (eigeninvert) THEN
                CALL eigen_invert_cc(cc_one, nmaps, rcond_limit, eigenvecs, &
                     eigenvals, workspace, workspacelen)
             ELSE
                CALL invert_cc(cc_one, nmaps)
             END IF
          END IF

          IF (maptemp(ipix_in,1).NE.HPX_DBADVAL) THEN
             cc_sum = cc_sum + cc_one
             map_out(ipix_out,:) = map_out(ipix_out,:) + &
                  MATMUL(cc_one, maptemp(ipix_in,:))
          ENDIF
       END DO

       IF (nmaps == 1) THEN
          IF (cc_sum(1,1)>0) THEN
             map_out(ipix_out, 1) = map_out(ipix_out, 1) / cc_sum(1, 1)
          ELSE
             map_out(ipix_out,1) = 0
          ENDIF
       ELSE
          ! Invert the summed Nobs matrix and weight the signals
          ! accordingly
          IF (diag_only) THEN
             DO irow = 1, nmaps
                cc_sum(irow,irow) = 1 / cc_sum(irow,irow)
             END DO
          ELSE IF (eigeninvert) THEN
             CALL eigen_invert_cc(cc_sum, nmaps, rcond_limit, eigenvecs, eigenvals, workspace, workspacelen)
          ELSE
             CALL invert_cc(cc_sum, nmaps)
          END IF

          IF (cc_sum(1, 1) == 0) THEN
             map_out(ipix_out,:) = 0
          ELSE
             map_out(ipix_out,:) = MATMUL(cc_sum, map_out(ipix_out,:))
          ENDIF
       END IF

    END DO
    !$OMP END DO
    DEALLOCATE(cc_one, cc_sum, eigenvecs, eigenvals, workspace)
    !$OMP END PARALLEL

    IF (nmaps == 3) DEALLOCATE(maptemp)

  CONTAINS

    SUBROUTINE eigen_invert_cc(cc, nmaps, rcondlim, eigenvecs, eigenvals,&
                               workspace, workspacelen)
      REAL(dp), INTENT(inout) :: cc(:,:),workspace(:)
      REAL(dp), INTENT(inout) :: eigenvecs(:,:), eigenvals(:)
      INTEGER(i4b), INTENT(in) :: nmaps
      INTEGER(i4b), INTENT(in) :: workspacelen
      REAL(dp), INTENT(in) :: rcondlim
      INTEGER(i4b) :: info, ieigen

      IF (cc(1,1)==0) RETURN

      IF (nmaps == 1) THEN
         IF (cc(1,1) /= 0) cc(1,1) = 1 / cc(1,1)
         RETURN
      END IF

      eigenvecs = cc
      CALL dsyev('V', 'U', nmaps, eigenvecs, nmaps, eigenvals, &
           workspace, workspacelen, info)

      IF (MINVAL(eigenvals)<=0) THEN
         cc = 0
         RETURN
      END IF

      IF (MINVAL(eigenvals)/MAXVAL(eigenvals) < rcondlim .OR. info /= 0) THEN
         cc = 0
         RETURN
      END IF
      DO ieigen = 1,nmaps
         eigenvecs(:,ieigen) = &
              eigenvecs(:,ieigen) / SQRT(eigenvals(ieigen))
      END DO
      cc = MATMUL(eigenvecs, TRANSPOSE(eigenvecs))
    END SUBROUTINE eigen_invert_cc


    RECURSIVE SUBROUTINE invert_cc(cc, nmaps)
      REAL(dp), INTENT(inout) :: cc(:,:)
      INTEGER(i4b), INTENT(in) :: nmaps
!!$      real(dp) :: a, b, c, d, e, f, norm
      INTEGER(i4b) :: info, i, j

      IF (nmaps == 1) THEN
         IF (cc(1,1) /= 0) cc(1,1) = 1 / cc(1,1)
         RETURN
      END IF

      info = 0
      CALL dpotrf('U', nmaps, cc, nmaps, info)

      IF (info /= 0) THEN
         cc = 0
         RETURN
      END IF

      info = 0
      CALL dpotri('U', nmaps, cc, nmaps, info)

      IF (info /= 0) THEN
         cc = 0
         RETURN
      END IF

      ! copy to lower triangle

      DO i = 1, nmaps
         DO j = i+1, nmaps
            cc(j,i) = cc(i,j)
         END DO
      END DO

      RETURN

!!$      if (nmaps /= 3) stop 'invert_cc not implemented for arbitrary sizes'
!!$
!!$      a = cc(1,1)
!!$      b = cc(1,2)
!!$      c = cc(1,3)
!!$      d = cc(2,2)
!!$      e = cc(2,3)
!!$      f = cc(3,3)
!!$      norm = 2.0*b*c*e - c**2*d - b**2*f + a*d*f - a*e**2
!!$      if (norm==0) then
!!$         cc = 0.0
!!$         return
!!$      endif
!!$
!!$      cc(1, 1) = (d*f  - e**2)/norm
!!$      cc(1, 2) = (-b*f  + c*e)/norm
!!$      cc(1, 3) = (b*e  - c*d)/norm
!!$      cc(2, 1) = cc(1, 2)
!!$      cc(2, 2) = (-c**2 + a*f)/norm
!!$      cc(2, 3) = (-a*e  + c*b)/norm
!!$      cc(3, 1) = cc(1, 3)
!!$      cc(3, 2) = cc(2, 3)
!!$      cc(3, 3) = (-b**2 + a*d)/norm
    END SUBROUTINE invert_cc

  END SUBROUTINE noise_weight_downgrade_nobs_preloaded


  !> @brief Smooth and optionally downgrade map using the cosine window.
  !!
  !! Smooth and optionally downgrade map by convolving the spherical harmonic
  !! expansion of the map with the cosine window.
  !!
  !! @note The cosine window is defined by two threshold multipoles,
  !! \f$\ell_1\f$ and \f$\ell_2\f$:
  !! \f[
  !!     w(\ell) = \left[1+\cos{\left((\ell-\ell_1)*\pi/(\ell_2-\ell_1)\right)}\right]/2
  !! \f]
  !! And the window obviously is 1 before \f$\ell_1\f$ and 0 after \f$\ell_2\f$.
  !! Currently these thresholds are hard-coded to
  !! \f$\ell_1 = 2*Nside_{out}\f$ and \f$\ell_2 = 3*Nside_{out}\f$
  !!
  !! @see Benabed et al., http://arxiv.org/pdf/0901.4537
  !!
  SUBROUTINE apodized_downgrade( &
       map_in, map_out, nside_in, nside_out, nmaps, verbose, polsmooth, &
       force_wpix, deconv_wpix)

    IMPLICIT NONE

    REAL(dp), POINTER, INTENT(INOUT) :: map_in(:,:) !< input map
    REAL(dp), POINTER, INTENT(INOUT) :: map_out(:,:) !< output map
    INTEGER(i4b), INTENT(IN) :: nside_in !< input resolution paramater
    INTEGER(i4b), INTENT(IN) :: nside_out !< output resolution. If same is input, no pixel window is applied
    INTEGER(i4b), INTENT(IN) :: nmaps !< number of Stokes components of the maps, 1 or 3
    LOGICAL(lgt), OPTIONAL, INTENT(IN) :: verbose !< optional input controlling verbosity
    LOGICAL, OPTIONAL, INTENT(IN) :: polsmooth !< smooth polarized components
    LOGICAL, OPTIONAL, INTENT(IN) :: force_wpix !< force pixel window smoothing even when nside_in=nside_out
    LOGICAL, OPTIONAL, INTENT(IN) :: deconv_wpix !< deconvolve input resolution pixel window
    REAL(dp), POINTER, SAVE :: pixelwindow(:,:)=>null(), beam(:,:)=>null()
    INTEGER(i4b) :: lmax, mmax, ell, ell1, ell2, ierr
    LOGICAL(lgt) :: verbose_in, nopol, force_wpix_in, deconv_wpix_in

    IF (nmaps /= 1 .AND. nmaps /= 3) THEN
       IF (ASSOCIATED(beam))  DEALLOCATE(beam)
       RETURN
    END IF

    verbose_in = .FALSE.
    IF (PRESENT(verbose)) verbose_in = verbose

    IF (PRESENT(polsmooth)) THEN
       nopol = .NOT. polsmooth
       IF (.NOT. polsmooth .AND. nside_in /= nside_out) &
            STOP 'Cannot downgrade with polsmooth=F'
    ELSE
       nopol = .FALSE.
    END IF

    force_wpix_in = .FALSE.
    IF (PRESENT(force_wpix)) force_wpix_in = force_wpix

    deconv_wpix_in = .TRUE.
    IF (PRESENT(deconv_wpix)) deconv_wpix_in = deconv_wpix

    lmax = 3*nside_out
    mmax = lmax

    IF (.NOT. ASSOCIATED(beam)) THEN
       ALLOCATE(beam(0:lmax, nmaps), pixelwindow(0:lmax, nmaps), stat=ierr)
       IF (ierr /= 0) STOP 'no room for smoothing beam'

       ! build the apodizing window
       beam = 0.0_dp

       ! DX11 choices establised in the CTP meeting in Paris on March 7 2014
       !ell1 = nside_out
       !ell2 = 3 * nside_out
       ! Experimental values to remove ringing from 857GHz
       ell1 = 1
       ell2 = 3 * nside_out
       ! original values
       !ell1 = 2*nside_out
       !ell2 = lmax ! 3*nside_out

       beam(:ell1, :) = 1.0_dp
       DO ell = ell1+1, ell2
          beam(ell, :) = (1.0 + COS((ell-ell1)*pi/(ell2-ell1))) / 2.0
       END DO

       ! add pixel window

       IF (nside_out < nside_in .OR. force_wpix_in) THEN
          CALL pixel_window(pixelwindow, nside_out)
          beam = beam * pixelwindow
          IF (nside_out < nside_in .AND. deconv_wpix_in) THEN
             ! Remove the input resolution pixel window
             CALL pixel_window(pixelwindow, nside_in)
             WHERE (pixelwindow < 1e-10) pixelwindow = 1
             beam = beam / pixelwindow
          END IF
       ELSE
          WRITE (*,*) 'WARNING: not applying pixel window.'
       END IF

       DEALLOCATE(pixelwindow)
    END IF

    CALL beam_downgrade(beam, lmax, mmax, map_in, map_out, nside_in, &
         nside_out, nmaps, verbose_in, nopol)

  END SUBROUTINE apodized_downgrade



  !! Apply the user-specified window function while downgrading
  SUBROUTINE windowfile_downgrade( &
       windowfile, map_in, map_out, nside_in, nside_out, nmaps, verbose, &
       polsmooth, force_wpix, deconv_wpix)

    IMPLICIT NONE

    CHARACTER(len=*) :: windowfile
    REAL(dp), POINTER, INTENT(INOUT) :: map_in(:,:) !< input map
    REAL(dp), POINTER, INTENT(INOUT) :: map_out(:,:) !< output map
    INTEGER(i4b), INTENT(IN) :: nside_in !< input resolution paramater
    INTEGER(i4b), INTENT(IN) :: nside_out !< output resolution. If same is input, no pixel window is applied
    INTEGER(i4b), INTENT(IN) :: nmaps !< number of Stokes components of the maps, 1 or 3
    LOGICAL(lgt), OPTIONAL, INTENT(IN) :: verbose !< optional input controlling verbosity
    LOGICAL, OPTIONAL, INTENT(IN) :: polsmooth !< smooth polarized components
    LOGICAL, OPTIONAL, INTENT(IN) :: force_wpix !< force pixel window smoothing even when nside_in=nside_out
    LOGICAL, OPTIONAL, INTENT(IN) :: deconv_wpix !< deconvolve input resolution pixel window
    REAL(dp), POINTER, SAVE :: pixelwindow(:,:)=>null(), beam(:,:)=>null()
    INTEGER(i4b) :: lmax, mmax, ierr, icol, ncol
    INTEGER(i8b) :: nelem
    LOGICAL(lgt) :: verbose_in, nopol, force_wpix_in, deconv_wpix_in
    REAL(dp), ALLOCATABLE :: inbeam(:,:)
    CHARACTER(len=80) :: header(1000)

    IF (nmaps /= 1 .AND. nmaps /= 3) THEN
       IF (ASSOCIATED(beam))  DEALLOCATE(beam)
       RETURN
    END IF

    verbose_in = .FALSE.
    IF (PRESENT(verbose)) verbose_in = verbose

    IF (PRESENT(polsmooth)) THEN
       nopol = .NOT. polsmooth
       IF (.NOT. polsmooth .AND. nside_in /= nside_out) &
            STOP 'Cannot downgrade with polsmooth=F'
    ELSE
       nopol = .FALSE.
    END IF

    force_wpix_in = .FALSE.
    IF (PRESENT(force_wpix)) force_wpix_in = force_wpix

    deconv_wpix_in = .TRUE.
    IF (PRESENT(deconv_wpix)) deconv_wpix_in = deconv_wpix

    lmax = 3*nside_out
    mmax = lmax

    IF (.NOT. ASSOCIATED(beam)) THEN
       ALLOCATE(beam(0:lmax, nmaps), pixelwindow(0:lmax, nmaps), stat=ierr)
       IF (ierr /= 0) STOP 'no room for smoothing beam'

       ! read the smoothing beam, the lines are from Hans Kristian with minor modifications

       ! Seem to remember there is something weird going on with the WMAP beams when reading only
       ! one component at a time. Remove this wrapper if you feel more comfortable with that...
       nelem = getsize_fits(windowfile, nmaps=ncol)
       IF (ncol < 4) THEN
          ncol = 1
       ELSE
          ncol = 4
       END IF

       ALLOCATE(inbeam(0:lmax,ncol))

       CALL fits2cl(windowfile, inbeam, lmax, ncol, header)
       IF (ncol == 1) THEN
          DO icol = 1,3
             beam(:,icol) = inbeam(:,1)
          END DO
       ELSE
          beam(0:lmax,1:nmaps) = inbeam(0:lmax,1:nmaps)
       END IF
       DEALLOCATE(inbeam)

       IF (nmaps > 1) THEN
          IF (SUM(beam(:,2)) < 1.d0) beam(:,2) = beam(:,1)
          IF (SUM(beam(:,3)) < 1.d0) beam(:,3) = beam(:,2)
       END IF

       ! add pixel window

       IF (nside_out < nside_in .OR. force_wpix_in) THEN
          CALL pixel_window(pixelwindow, nside_out)
          beam = beam * pixelwindow
          IF (nside_out < nside_in .AND. deconv_wpix_in) THEN
             ! Remove the input resolution pixel window
             CALL pixel_window(pixelwindow, nside_in)
             WHERE (pixelwindow < 1e-10) pixelwindow = 1
             beam = beam / pixelwindow
          END IF
       ELSE
          WRITE (*,*) 'WARNING: not applying pixel window.'
       END IF

       DEALLOCATE(pixelwindow)
    END IF

    CALL beam_downgrade(beam, lmax, mmax, map_in, map_out, nside_in, nside_out, &
         nmaps, verbose_in, nopol)

  END SUBROUTINE windowfile_downgrade



  !! Downgrade map using harmonic de/recomposition while applying a precomputed beam
  SUBROUTINE beam_downgrade(beam, lmax, mmax, map_in, map_out, nside_in, &
       nside_out, nmaps, verbose, nopol)

    IMPLICIT NONE

    REAL(dp), POINTER, INTENT(IN) :: beam(:,:)
    INTEGER(i4b), INTENT(IN) :: lmax, mmax
    REAL(dp), POINTER, INTENT(INOUT) :: map_in(:,:) !< input map
    REAL(dp), POINTER, INTENT(INOUT) :: map_out(:,:) !< output map
    INTEGER(i4b), INTENT(IN) :: nside_in !< input resolution paramater
    INTEGER(i4b), INTENT(IN) :: nside_out !< output resolution. If same is input, no pixel window is applied
    INTEGER(i4b), INTENT(IN) :: nmaps !< number of Stokes components of the maps, 1 or 3
    LOGICAL(lgt), INTENT(IN) :: verbose !< optional input controlling verbosity
    LOGICAL, INTENT(IN) :: nopol !< smooth polarized components
    REAL(dp), POINTER, SAVE :: weights(:,:)=>null()
    REAL(dp), POINTER, SAVE :: map(:,:)=>null()
    REAL(dp) :: zbounds(2)
    COMPLEX(dpc), POINTER, SAVE :: alm(:,:,:)=>null()
    INTEGER(i4b) :: ierr

    IF (.NOT. ASSOCIATED(beam)) STOP 'Cannot deconvolve an unallocated beam'

    zbounds = (/ 0_dp, 0_dp /)

    ALLOCATE(weights(2*nside_in, nmaps), alm(nmaps, 0:lmax, 0:mmax), &
         map(0:nside2npix(nside_in)-1,nmaps), stat=ierr)
    IF (ierr /= 0) STOP 'no room for harmonic analysis'

    weights = 1.0_dp
    map = map_in

    ! Expand the input map

    IF (verbose) CALL tic(555)
    CALL convert_nest2ring(nside_in, map)
    IF (nMaps == 1 .OR. nopol) THEN
       CALL map2alm(nside_in, lmax, mmax, map(:,1), alm(1:1,:,:), zbounds, weights)
    ELSE
       CALL map2alm(nside_in, lmax, mmax, map, alm, zbounds, weights)
    END IF
    IF (verbose) CALL toc('harmonic analysis', 555)

    ! Apply the window

    CALL alter_alm(nside_out, lmax, mmax, 0.0_dp, alm, window=beam)
    IF (verbose) CALL toc('apply windows', 555)

    ! re-synthesize

    IF (.NOT. ASSOCIATED(map_out)) THEN
       ALLOCATE(map_out(0:nside2npix(nside_out)-1, nmaps), stat=ierr)
       IF (ierr /= 0) STOP 'no room for output map'
    END IF

    IF (verbose) CALL tic(555)
    IF (nmaps == 1 .OR. nopol) THEN
       CALL alm2map(nside_out, lmax, mmax, alm(1:1,:,:), map_out(:,1))
    ELSE
       CALL alm2map(nside_out, lmax, mmax, alm, map_out)
    END IF
    CALL convert_ring2nest(nside_out, map_out)
    IF (verbose) CALL toc('resynthesize', 555)

    ! optionally return the polarization maps to their original form

    IF (nopol .AND. nmaps==3) map_out(:,2:3) = map_in(:,2:3)

  END SUBROUTINE beam_downgrade



  !> Smooth and optionally downgrade a map using a gaussian beam.
  SUBROUTINE smooth_downgrade( &
       fwhm, map_in, map_out, nside_in, nside_out, nmaps, verbose, polsmooth, &
       force_wpix, deconv_wpix)
    ! expand high resolution map in spherical harmonics and apply
    ! smoothing beam, before synthesizing into smaller resolution

    IMPLICIT NONE

    REAL(dp), INTENT(IN)  :: fwhm                  !< Width of the Gaussian beam in arc minutes
    REAL(dp), POINTER, INTENT(INOUT) :: map_in(:,:)!< input map
    REAL(dp), POINTER, INTENT(INOUT) :: map_out(:,:) !< smoothed map
    INTEGER(i4b), INTENT(IN) :: nside_in           !< input resolution
    INTEGER(i4b), INTENT(IN) :: nside_out          !< output resolution
    INTEGER(i4b), INTENT(IN) :: nmaps              !< number of Stokes components, 1 or 3
    LOGICAL(lgt), OPTIONAL, INTENT(IN) :: verbose !< optional input controlling verbosity
    LOGICAL, OPTIONAL, INTENT(IN) :: polsmooth !< smooth polarized components
    LOGICAL, OPTIONAL, INTENT(IN) :: force_wpix !< force pixel window smoothing even when nside_in=nside_out
    LOGICAL, OPTIONAL, INTENT(IN) :: deconv_wpix !< deconvolve input resolution pixel window
    REAL(dp), POINTER, SAVE :: pixelwindow(:,:)=>null(), beam(:,:)=>null()
    INTEGER(i4b) :: lmax, mmax, ierr
    LOGICAL(lgt) :: verbose_in, nopol, force_wpix_in, deconv_wpix_in

    IF (nmaps /= 1 .AND. nmaps /= 3) THEN
       IF (ASSOCIATED(beam))  DEALLOCATE(beam)
       RETURN
    END IF

    verbose_in = .FALSE.
    IF (PRESENT(verbose)) verbose_in = verbose

    IF (PRESENT(polsmooth)) THEN
       nopol = .NOT. polsmooth
       IF (.NOT. polsmooth .AND. nside_in /= nside_out) &
            STOP 'Cannot downgrade with polsmooth=F'
    ELSE
       nopol = .FALSE.
    END IF

    force_wpix_in = .FALSE.
    IF (PRESENT(force_wpix)) force_wpix_in = force_wpix

    deconv_wpix_in = .TRUE.
    IF (PRESENT(deconv_wpix)) deconv_wpix_in = deconv_wpix

    lmax = 3*nside_out
    mmax = lmax

    IF (.NOT. ASSOCIATED(beam)) THEN
       ALLOCATE(beam(0:lmax, nmaps), pixelwindow(0:lmax, nmaps), stat=ierr)
       IF (ierr /= 0) STOP 'no room for smoothing beam'

       ! build the gaussian window

       CALL gaussbeam(fwhm, lmax, beam)

       ! add pixel window to the beam

       IF (nside_out < nside_in .OR. force_wpix_in) THEN
          ! Apply output resolution pixel window
          CALL pixel_window(pixelwindow, nside_out)
          beam = beam * pixelwindow
          IF (nside_out < nside_in .AND. deconv_wpix_in) THEN
             ! Remove the input resolution pixel window
             CALL pixel_window(pixelwindow, nside_in)
             WHERE (pixelwindow < 1e-10) pixelwindow = 1
             beam = beam / pixelwindow
          END IF
       ELSE
          WRITE (*,*) 'WARNING: not applying pixel window.'
       END IF

       DEALLOCATE(pixelwindow)
    END IF

    CALL beam_downgrade(beam, lmax, mmax, map_in, map_out, nside_in, nside_out, &
         nmaps, verbose_in, nopol)


  END SUBROUTINE smooth_downgrade



  SUBROUTINE pixelsmooth_downgrade(map_in, map_out, nside_in, nside_out, nmaps)
    ! apply a neighbourhood average pixel space filter

    IMPLICIT NONE

    REAL(dp), POINTER :: map_in(:,:), map_out(:,:), map(:,:)
    REAL(dp), POINTER :: maptemp(:,:)
    INTEGER(i4b) :: nside_in, nside_out, nmaps, ierr, isig, nlist
    INTEGER(i8b) :: list(8), ipixel, npix_out

    IF (nmaps == 3) THEN
       ALLOCATE(maptemp(0:12*nside_in**2-1, nmaps), stat=ierr)
       IF (ierr /= 0) STOP 'No room for parallel transport'
       maptemp = map_in
       CALL qu_transport_map(nside_out, maptemp)
    ELSE
       maptemp => map_in
    END IF

    npix_out = nside2npix(nside_out)

    ALLOCATE(map(0:npix_out-1,nmaps))

    IF (nside_out /= nside_in) THEN
       CALL udgrade_nest(maptemp, nside_in, map, nside_out)
    ELSE
       map = map_in
    END IF

    DO ipixel = 0, npix_out-1
       CALL neighbours_nest(nside_out, ipixel, list, nlist)
       DO isig = 1, nmaps
          map_out(ipixel, isig) = SUM(map(list(1:nlist),isig))/nlist
       END DO
    END DO

    IF (nmaps == 3) DEALLOCATE(maptemp)

  END SUBROUTINE pixelsmooth_downgrade



  SUBROUTINE average_downgrade(map_in, map_out, nside_in, nside_out, nmaps)
    ! simple averaging

    IMPLICIT NONE

    REAL(dp), POINTER :: map_in(:,:), map_out(:,:)
    REAL(dp), POINTER :: maptemp(:,:)
    INTEGER(i4b) :: nside_in, nside_out, nmaps, ierr

    IF (nmaps == 3) THEN
       ALLOCATE(maptemp(0:12*nside_in**2-1, nmaps), stat=ierr)
       IF (ierr /= 0) STOP 'No room for parallel transport'
       maptemp = map_in
       CALL qu_transport_map(nside_out, maptemp)
    ELSE
       maptemp => map_in
    END IF

    CALL udgrade_nest(maptemp, nside_in, map_out, nside_out, fmissval=0d0)

    IF (nmaps == 3) DEALLOCATE(maptemp)

  END SUBROUTINE average_downgrade



END MODULE downgrade
