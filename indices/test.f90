PROGRAM test
! Program to test Fotran subroutines & functions of retrieval of 1D indices with Fortran
!   L. Fita, Centro de Investigaciones del Mar y la Atmosfera (CIMA), UBA-CONICET, CNRS IRL 3351 IFAECI, IRD, Argentina
!
  
  USE module_definitions
  USE module_basic
  USE module_FortranIndices

  IMPLICIT NONE
  
  CHARACTER(len=Sm)                                      :: Sma, Smb
  INTEGER                                                :: i, j ,k 
  INTEGER                                                :: dimx, dimy, dimz, dimt
  REAL                                                   :: Ra, Rb
  REAL, DIMENSION(:), ALLOCATABLE                        :: RD1a, RD1b
  REAL, DIMENSION(:,:), ALLOCATABLE                      :: RD2a, RD2b
  REAL, DIMENSION(:,:,:), ALLOCATABLE                    :: RD3a, RD3b
  LOGICAL                                                :: La, Lb

  ! Test Quantiles
  dimt = 100
  Ra = FillValueR
  ALLOCATE(RD1a(dimt))
  DO i=1, dimt
    RD1a(i) = i*1.
  END DO
  ALLOCATE(RD1b(Npercents))
  
  PRINT *,' rank ____'
  CALL percentilesR(dimt, RD1a, Ra, Npercents, RD1b)
  DO i=1, Npercents
    PRINT *, (i-1)*5, ' %: ', RD1b(i) 
  END DO
  
  PRINT *,' interp ____'
  Sma = 'C0'
  CALL percentilesR(dimt, RD1a, Ra, Npercents, RD1b, Sma)
  DO i=1, Npercents
    PRINT *, (i-1)*5, ' %: ', RD1b(i) 
  END DO

  RD1a(3:7) = FillValueR
  
  CALL percentilesR(dimt, RD1a, Ra, Npercents, RD1b)
  DO i=1, Npercents
    PRINT *, (i-1)*5, ' %: ', RD1b(i) 
  END DO

  
  PRINT *,' interp ____'
  CALL percentilesR(dimt, RD1a, Ra, Npercents, RD1b, Sma)
  DO i=1, Npercents
    PRINT *, (i-1)*5, ' %: ', RD1b(i) 
  END DO

  STOP

END PROGRAM test
