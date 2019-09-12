! 03-08-2007 RK
! short code that checks how much memory a code can allocate
! The method works only with operating systems that free
! deallocated memory immediately upon request

MODULE memtest

  IMPLICIT NONE

  INTEGER, PARAMETER, PRIVATE :: idp = SELECTED_INT_KIND(10)
  INTEGER, PARAMETER, PRIVATE :: dp  = SELECTED_REAL_KIND(10,100)
  INTEGER(kind=idp), PRIVATE  :: lower_limit, upper_limit, middle_limit
  REAL(kind=dp), POINTER      :: dummytable(:)
  INTEGER, PRIVATE            :: ierr=0, memtest_verbosity = 0

CONTAINS

  FUNCTION maximum_allocatable_memory(verbosity)
    ! This function returns the amount of allocatable
    ! memory in MB
    !
    ! verbosity = 0 : no printing
    ! verbosity = 1 : print available memory
    ! verbosity > 1 : print every step
    IMPLICIT NONE

    INTEGER, OPTIONAL, INTENT(in) :: verbosity
    INTEGER :: maximum_allocatable_memory

    IF (PRESENT(verbosity)) memtest_verbosity = verbosity

    lower_limit = 2**20/8 ! one megabyte in double precision
    upper_limit = lower_limit

    IF (memtest_verbosity > 1) &
         WRITE (*,'(a,i8,a,/)') 'Starting from ', upper_limit*8/2**20, 'MB'

    DO
       ALLOCATE(dummytable(upper_limit),stat=ierr)
       IF (ierr==0) THEN
          DEALLOCATE(dummytable)
          IF (memtest_verbosity > 1) &
               WRITE (*,'(a,i8,a)') &
               'Succesfully allocated ', upper_limit*8/2**20, 'MB'
          lower_limit = upper_limit
          upper_limit = upper_limit*2
       ELSE
          EXIT
       END IF
    END DO

    IF (memtest_verbosity > 1) &
         WRITE (*,'(a,i8,a)') 'Failed to allocate ', upper_limit*8/2**20, 'MB'
    
    ! After the loop, maximum allocable memory is between the limits.
    ! Use binary search to fork the limit.
    DO
       
       IF (lower_limit == upper_limit) EXIT
       
       middle_limit = lower_limit + (upper_limit-lower_limit)/2
       ALLOCATE(dummytable(middle_limit),stat=ierr)
       IF (ierr==0) THEN
          DEALLOCATE(dummytable)
          IF (memtest_verbosity > 1) &
               WRITE (*,'(a,i8,a)') &
               'Succesfully allocated ', upper_limit*8/2**20, 'MB'
          lower_limit = middle_limit
       ELSE
          IF (memtest_verbosity > 1) &
               WRITE (*,'(a,i8,a)') &
               'Failed to allocate ', middle_limit*8/2**20, 'MB'
          upper_limit = middle_limit
       END IF
       
    END DO

    IF (memtest_verbosity > 0) &
         WRITE (*,'(/,a,i8,a)') &
         'Maximum allocatable memory ', lower_limit*8/2**20, 'MB'

    maximum_allocatable_memory = lower_limit*8/2**20

  END FUNCTION maximum_allocatable_memory

END MODULE memtest
