MODULE logModule

  IMPLICIT NONE

  PRIVATE
  PUBLIC logInit, logInfo

  INTEGER :: logUnit = 6

CONTAINS

  SUBROUTINE logInit(unit)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: unit

    logUnit = unit

  END SUBROUTINE logInit

  SUBROUTINE logInfo(message)

    IMPLICIT NONE

    CHARACTER (LEN=*), INTENT(IN) :: message

    WRITE (logUnit, '(A,A)') 'info: ', message

  END SUBROUTINE logInfo

END MODULE logModule
