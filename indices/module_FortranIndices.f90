MODULE module_FortranIndices
! Module to compute 1D Indices using Fortran

!!!!!!!! Contents

  USE module_definitions
  USE module_basic

  CONTAINS

  SUBROUTINE index1D_from3D_d3(d1, d2, d3, indexn, matv, missv, indexv)
  ! Subroutine to compute 1D indices using 3rd dimension from 3D arrays producing 2D arrays
  
    IMPLICIT NONE
    
    INTEGER, INTENT(in)                                  :: d1, d2, d3
    CHARACTER(len=*), INTENT(in)                         :: indexn
    REAL, INTENT(in)                                     :: missv    
    REAL, DIMENSION(d1,d2,d3), INTENT(in)                :: matv
    REAL, DIMENSION(d1,d2), INTENT(out)                  :: indexv

    ! Local
    INTEGER                                              :: i, j, k
    INTEGER                                              :: Navailindexn
    INTEGER                                              :: Ia
    CHARACTER(len=Sm), DIMENSION(:), ALLOCATABLE         :: availindexn
    CHARACTER(len=Sm)                                    :: fname
    
!!!!!!! Variables
! d1, d2, d3:  shape of the 3D array
! indexn: name of the index to compute
! matv: 3D array with the values
! missv: missing value of matv

  fname = 'index1D_from3D_d3'
  
  Navailindexn = 1
  IF (ALLOCATED(availindexn)) DEALLOCATE(availindexn)
  ALLOCATE(availindexn(Navailindexn))
  
  availindexn(1) = 'mean'
  
  indexv = fillValueR
  
  SELECT CASE (TRIM(indexn))
  
    CASE ('mean')
    
      indexv = 0.
      ! Amount of values
      Ia = 0
    
      DO i=1, d1
        DO j=1, d2
          DO k=1, d3
            IF (matv(i,j,k) /= missv) THEN
              indexv(j,k) = indexv(j,k) + matv(i,j,k)
              Ia = Ia + 1
            END IF
          END DO
          indexv = indexv / Ia
        END DO
      END DO

    CASE DEFAULT
    
      msg = "1D index '" // TRIM(indexn) // "' not ready !!"

      CALL StopRunAvail(msg, fname, Navailindexn, availindexn)
    
    END SELECT

    DEALLOCATE(availindexn)

    RETURN

  END SUBROUTINE index1D_from3D_d3
END MODULE module_FortranIndices

