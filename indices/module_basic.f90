MODULE module_basic
! Module with the basic instructions

! Content

!!!!!!! Subroutines/Functions
! ItoS: Function to transform an integer to String
! StopRun: Subroutine to stop execution and provide an error message
! StopRunAvail: Subroutine to stop running and print a message and showing available options

  USE module_definitions
  
  CONTAINS
  
  CHARACTER(LEN=Sm) FUNCTION ItoS(Ival)
! ItoS: Function to transform an integer to String

    IMPLICIT NONE

    INTEGER, INTENT(IN)                                  :: Ival

! Local
    CHARACTER(LEN=Sm)                                    :: itoS0

    WRITE(ItoS0,'(I50)')Ival
    CALL  removeNONnum(ItoS0, ItoS)

  END FUNCTION ItoS

  SUBROUTINE StopRunAvail(msg, fname, Navail, avail)
! Subroutine to stop running and print a message and showing available options

    IMPLICIT NONE

    INTEGER, INTENT(in)                                    :: Navail
    CHARACTER(LEN=*), INTENT(IN)                           :: fname
    CHARACTER(LEN=*), INTENT(IN)                           :: msg
    CHARACTER(LEN=*), DIMENSION(Navail), INTENT(IN)        :: avail

! local
    INTEGER                                                :: i
    CHARACTER(LEN=50)                                      :: errmsg, warnmsg

    errmsg = 'ERROR -- error -- ERROR -- error'

    PRINT *, TRIM(errmsg)
    PRINT *, '  ' // TRIM(fname) // ': ' // TRIM(msg)
    PRINT *, '    available ones:', TRIM(avail(1)), (', ' // TRIM(avail(i)), i=2, Navail)
    STOP

  END SUBROUTINE stoprunAvail
  
  SUBROUTINE StopRun(msg, funcn, errN)
  ! Subroutine to stop execution and provide an error message

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(in)                           :: msg, funcn
    INTEGER, INTENT(in)                                    :: errN

    ! Local

!!!!!!! Variables
! msg: message to print with the error
! funcn: name of the funciton, section to localize the error
! errN: number of the error provided for a given FORTRAN function

    IF (errN /= 0) THEN
      PRINT *, errormsg
      PRINT *, '  ' // TRIM(funcn) // ': ' // TRIM(msg) // ' !!'
      PRINT *, '    error number:', errN
      STOP
    END IF

    RETURN

  END SUBROUTINE StopRun

END MODULE module_basic
