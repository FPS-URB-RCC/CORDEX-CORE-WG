MODULE module_basic
! Module with the basic instructions

! Content

!!!!!!! Subroutines/Functions
! FindMinimumR: Function returns the location of the real minimum in the section between Start and End
! ItoS: Function to transform an integer to String
! removeNONnum: Subroutine to remove non numeric characters from a string
! SortR: Subroutine receives an array x real and sorts it into ascending order
! StopRun: Subroutine to stop execution and provide an error message
! StopRunAvail: Subroutine to stop running and print a message and showing available options
! SwapR: Subroutine swaps the real values of its two formal arguments

  USE module_definitions
  
  CONTAINS
  
  INTEGER FUNCTION FindMinimumR(x, dsize, Startv, Endv, missv)
  ! Function returns the location of the real minimum in the section between Start and End

    IMPLICIT NONE

    INTEGER, INTENT(in)                                  :: dsize
    REAL, DIMENSION(dsize), INTENT(in)                   :: x
    INTEGER, INTENT(in)                                  :: Startv, Endv
    REAL, INTENT(in)                                     :: missv

    ! Local
    REAL                                                 :: Minimum
    INTEGER                                              :: Location
    INTEGER                                              :: i
    CHARACTER(len=Sm)                                    :: fname
    
    fname = 'FindMinimumR'

    Minimum  = x(Startv)                                 ! assume the first is the min
    Location = Startv                                    ! record its position
    DO i = Startv+1, Endv                                ! start with next elements
      IF (x(i) /= missv .AND. x(i) < Minimum) THEN       !   if x(i) less than the min?
        Minimum  = x(i)                                  !      Yes, a new minimum found
        Location = i                                     !      record its position
      END IF
    END DO

    FindMinimumR = Location                              ! return the position

  END FUNCTION FindMinimumR

  SUBROUTINE SwapR(a, b)
  ! Subroutine swaps the real values of its two formal arguments

    IMPLICIT NONE

    REAL, INTENT(inout)                                  :: a, b

    ! Local
    REAL                                                 :: Temp
    CHARACTER(len=Sm)                                    :: fname
    
    fname = 'SwapR'

    Temp = a
    a    = b
    b    = Temp

    RETURN

  END SUBROUTINE SwapR

  SUBROUTINE SortR(x, Nx, missv)
  ! Subroutine receives an array x real and sorts it into ascending order

    IMPLICIT NONE

    INTEGER, INTENT(in)                                  :: Nx
    REAL, INTENT(in)                                     :: missv
    REAL, DIMENSION(Nx), INTENT(inout)                   :: x

    ! Local
    INTEGER                                              :: i
    INTEGER                                              :: Location
    CHARACTER(len=Sm)                                    :: fname

!!!!!!! Variables
! x: values to sort
! Nx: amount of values to sort
! missv: value of missing value
    
    fname = 'SortR'

    DO i = 1, Nx-1                                       ! except for the last
      Location = FindMinimumR(x, Nx, i, Nx, missv)       ! find min from this to last
      CALL  SwapR(x(i), x(Location))                     ! swap this and the minimum
    END DO
     
    RETURN

  END SUBROUTINE SortR
   
  SUBROUTINE removeNONnum(String, newString)
! Subroutine to remove non numeric characters from a string

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)                         :: String
    CHARACTER(LEN=*), INTENT(OUT)                        :: newString

! Local
    INTEGER                                              :: ic, inc, Lstring
    CHARACTER(len=50)                                    :: fname

!!!!!!! Variables
! String: string to remove non-numeric characters
! newString: resultant string
    fname = 'removeNONnum'

    Lstring = LEN_TRIM(String)
    newString = ''
    inc = 1
    DO ic=1, Lstring
      IF (ICHAR(String(ic:ic)) >= ICHAR('0') .AND. ICHAR(String(ic:ic)) <= ICHAR('9')) THEN
        newString(inc:inc) = String(ic:ic)
        inc = inc + 1
      END IF
    END DO

  END SUBROUTINE removeNONnum
  
! This does not seems to like to f2py
!  CHARACTER(LEN=Sm) FUNCTION ItoS(Ival)
  CHARACTER(LEN=50) FUNCTION ItoS(Ival)
! ItoS: Function to transform an integer to String

    IMPLICIT NONE

    INTEGER, INTENT(in)                                  :: Ival

! Local
    CHARACTER(LEN=Sm)                                    :: itoS0

    WRITE(ItoS0,'(I50)')Ival
    CALL  removeNONnum(ItoS0, ItoS)

  END FUNCTION ItoS

  SUBROUTINE StopRunAvail(msg, fname, Navail, avail)
! Subroutine to stop running and print a message and showing available options

    IMPLICIT NONE

    INTEGER, INTENT(in)                                    :: Navail
    CHARACTER(LEN=*), INTENT(in)                           :: fname
    CHARACTER(LEN=*), INTENT(in)                           :: msg
    CHARACTER(LEN=*), DIMENSION(Navail), INTENT(in)        :: avail

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
