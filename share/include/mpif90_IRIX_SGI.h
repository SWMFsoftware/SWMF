!DIR$ ID "@(#)mpi/include/mpif.h	31.6	05/27/98 13:24:49"

!
!	(C) COPYRIGHT SILICON GRAPHICS, INC.
!	UNPUBLISHED PROPRIETARY INFORMATION.
!	ALL RIGHTS RESERVED.
!

!
! ---------------------------------------
! This is mpif.h for IRIX and PVP systems
! ---------------------------------------
!
! Copyright Notice
!  + 1993 University of Chicago
!  + 1993 Mississippi State University


        character*100 MPI_HEADER_FILE
        parameter (MPI_HEADER_FILE = 'mpif90_IRIX.h')

	integer MPI_VERSION
	integer MPI_SUBVERSION

	parameter (MPI_VERSION		= 1)
	parameter (MPI_SUBVERSION	= 2)

! MPI_Status

	integer MPI_SOURCE
	integer MPI_TAG
	integer MPI_ERROR
	integer MPI_STATUS_SIZE

	parameter (MPI_SOURCE		= 1)
	parameter (MPI_TAG		= 2)
	parameter (MPI_ERROR		= 3)
	parameter (MPI_STATUS_SIZE	= 6)

! MPI_Comm

	integer MPI_COMM_NULL
	integer MPI_COMM_WORLD
	integer MPI_COMM_SELF

	parameter (MPI_COMM_NULL	= 0)
	parameter (MPI_COMM_WORLD	= 1)
	parameter (MPI_COMM_SELF	= 2)

! MPI_Errhandler

	integer MPI_ERRHANDLER_NULL
	integer MPI_ERRORS_ARE_FATAL
	integer MPI_ERRORS_RETURN

	parameter (MPI_ERRHANDLER_NULL	= 0)
	parameter (MPI_ERRORS_ARE_FATAL	= 1)
	parameter (MPI_ERRORS_RETURN	= 2)

! MPI_Group

	integer MPI_GROUP_NULL
	integer MPI_GROUP_EMPTY

	parameter (MPI_GROUP_NULL	= 0)
	parameter (MPI_GROUP_EMPTY	= 1)

! MPI_Request

	integer MPI_REQUEST_NULL

	parameter (MPI_REQUEST_NULL	= 0)

! MPI_Op

	integer MPI_OP_NULL
	integer MPI_MAX
	integer MPI_MIN
	integer MPI_SUM
	integer MPI_PROD
	integer MPI_LAND
	integer MPI_BAND
	integer MPI_LOR
	integer MPI_BOR
	integer MPI_LXOR
	integer MPI_BXOR
	integer MPI_MAXLOC
	integer MPI_MINLOC

	parameter (MPI_OP_NULL	= 0)
	parameter (MPI_MAX	= 1)
	parameter (MPI_MIN	= 2)
	parameter (MPI_SUM	= 3)
	parameter (MPI_PROD	= 4)
	parameter (MPI_LAND	= 5)
	parameter (MPI_BAND	= 6)
	parameter (MPI_LOR	= 7)
	parameter (MPI_BOR	= 8)
	parameter (MPI_LXOR	= 9)
	parameter (MPI_BXOR	= 10)
	parameter (MPI_MAXLOC	= 11)
	parameter (MPI_MINLOC	= 12)

! MPI_Datatype

	integer MPI_DATATYPE_NULL

	integer MPI_CHAR
	integer MPI_SHORT
	integer MPI_INT
	integer MPI_LONG
	integer MPI_UNSIGNED_CHAR
	integer MPI_UNSIGNED_SHORT
	integer MPI_UNSIGNED
	integer MPI_UNSIGNED_LONG
	integer MPI_FLOAT
	integer MPI_DOUBLE
	integer MPI_LONG_DOUBLE
	integer MPI_LONG_LONG
	integer MPI_LONG_LONG_INT

	integer MPI_INTEGER
	integer MPI_REAL
	integer MPI_DOUBLE_PRECISION
	integer MPI_COMPLEX
	integer MPI_DOUBLE_COMPLEX
	integer MPI_LOGICAL
	integer MPI_CHARACTER
	integer MPI_INTEGER1
	integer MPI_INTEGER2
	integer MPI_INTEGER4
	integer MPI_INTEGER8
	integer MPI_REAL4
	integer MPI_REAL8
	integer MPI_REAL16

	integer MPI_BYTE
	integer MPI_PACKED
	integer MPI_UB
	integer MPI_LB

	integer MPI_FLOAT_INT
	integer MPI_DOUBLE_INT
	integer MPI_LONG_INT
	integer MPI_2INT
	integer MPI_SHORT_INT
	integer MPI_LONG_DOUBLE_INT

	integer MPI_2REAL
	integer MPI_2DOUBLE_PRECISION
	integer MPI_2INTEGER

	parameter (MPI_DATATYPE_NULL	= 0)

	parameter (MPI_CHAR		= 1)
	parameter (MPI_SHORT		= 2)
	parameter (MPI_INT		= 3)
	parameter (MPI_LONG		= 4)
	parameter (MPI_UNSIGNED_CHAR	= 5)
	parameter (MPI_UNSIGNED_SHORT	= 6)
	parameter (MPI_UNSIGNED		= 7)
	parameter (MPI_UNSIGNED_LONG	= 8)
	parameter (MPI_FLOAT		= 9)
	parameter (MPI_DOUBLE		= 10)
	parameter (MPI_LONG_DOUBLE	= 11)
	parameter (MPI_LONG_LONG	= 12)
	parameter (MPI_LONG_LONG_INT	= 12)

	parameter (MPI_INTEGER		= 13)
	parameter (MPI_REAL		= 14+iRealPrec)
	parameter (MPI_DOUBLE_PRECISION	= 15)
	parameter (MPI_COMPLEX		= 16+iRealPrec)
	parameter (MPI_DOUBLE_COMPLEX	= 17)
	parameter (MPI_LOGICAL		= 18)
	parameter (MPI_CHARACTER	= 19)
	parameter (MPI_INTEGER1		= 20)
	parameter (MPI_INTEGER2		= 21)
	parameter (MPI_INTEGER4		= 22)
	parameter (MPI_INTEGER8		= 23)
	parameter (MPI_REAL4		= 24)
	parameter (MPI_REAL8		= 25)
	parameter (MPI_REAL16		= 26)

	parameter (MPI_BYTE		= 27)
	parameter (MPI_PACKED		= 28)
	parameter (MPI_UB		= 29)
	parameter (MPI_LB		= 30)

	parameter (MPI_FLOAT_INT	= 31)
	parameter (MPI_DOUBLE_INT	= 32)
	parameter (MPI_LONG_INT		= 33)
	parameter (MPI_2INT		= 34)
	parameter (MPI_SHORT_INT	= 35)
	parameter (MPI_LONG_DOUBLE_INT	= 36)

	parameter (MPI_2REAL		= 37)
	parameter (MPI_2DOUBLE_PRECISION= 38)
	parameter (MPI_2INTEGER		= 39)

! MPI-1 error codes and classes

	integer MPI_SUCCESS
	integer MPI_ERR_BUFFER
	integer MPI_ERR_COUNT
	integer MPI_ERR_TYPE
	integer MPI_ERR_TAG
	integer MPI_ERR_COMM
	integer MPI_ERR_RANK
	integer MPI_ERR_REQUEST
	integer MPI_ERR_ROOT
	integer MPI_ERR_GROUP
	integer MPI_ERR_OP
	integer MPI_ERR_TOPOLOGY
	integer MPI_ERR_DIMS
	integer MPI_ERR_ARG
	integer MPI_ERR_UNKNOWN
	integer MPI_ERR_TRUNCATE
	integer MPI_ERR_OTHER
	integer MPI_ERR_INTERN
	integer MPI_ERR_IN_STATUS
	integer MPI_ERR_PENDING

	parameter (MPI_SUCCESS		= 0)
	parameter (MPI_ERR_BUFFER	= 1)
	parameter (MPI_ERR_COUNT	= 2)
	parameter (MPI_ERR_TYPE		= 3)
	parameter (MPI_ERR_TAG		= 4) 
	parameter (MPI_ERR_COMM		= 5)
	parameter (MPI_ERR_RANK		= 6)
	parameter (MPI_ERR_REQUEST	= 7)
	parameter (MPI_ERR_ROOT		= 8)
	parameter (MPI_ERR_GROUP	= 9)
	parameter (MPI_ERR_OP		= 10)
	parameter (MPI_ERR_TOPOLOGY	= 11)
	parameter (MPI_ERR_DIMS		= 12)
	parameter (MPI_ERR_ARG		= 13)
	parameter (MPI_ERR_UNKNOWN	= 14)
	parameter (MPI_ERR_TRUNCATE	= 15)
	parameter (MPI_ERR_OTHER	= 16)
	parameter (MPI_ERR_INTERN	= 17)
	parameter (MPI_ERR_IN_STATUS	= 18)
	parameter (MPI_ERR_PENDING	= 19)

! MPI-2 error codes and classes

	integer MPI_ERR_ACCESS
	integer MPI_ERR_AMODE
	integer MPI_ERR_ASSERT
	integer MPI_ERR_BAD_FILE
	integer MPI_ERR_BASE
	integer MPI_ERR_CONVERSION
	integer MPI_ERR_DISP
	integer MPI_ERR_DUP_DATAREP
	integer MPI_ERR_FILE_EXISTS
	integer MPI_ERR_FILE_IN_USE
	integer MPI_ERR_FILE
	integer MPI_ERR_INFO_KEY
	integer MPI_ERR_INFO_NOKEY
	integer MPI_ERR_INFO_VALUE
	integer MPI_ERR_INFO
	integer MPI_ERR_IO
	integer MPI_ERR_KEYVAL
	integer MPI_ERR_LOCKTYPE
	integer MPI_ERR_NAME
	integer MPI_ERR_NO_MEM
	integer MPI_ERR_NOT_SAME
	integer MPI_ERR_NO_SPACE
	integer MPI_ERR_NO_SUCH_FILE
	integer MPI_ERR_PORT
	integer MPI_ERR_QUOTA
	integer MPI_ERR_READ_ONLY
	integer MPI_ERR_RMA_CONFLICT
	integer MPI_ERR_RMA_SYNC
	integer MPI_ERR_SERVICE
	integer MPI_ERR_SIZE
	integer MPI_ERR_SPAWN
	integer MPI_ERR_UNSUPPORTED_DATAREP
	integer MPI_ERR_UNSUPPORTED_OPERATION
	integer MPI_ERR_WIN

	parameter (MPI_ERR_ACCESS		= 28)
	parameter (MPI_ERR_AMODE		= 29)
	parameter (MPI_ERR_ASSERT		= 30)
	parameter (MPI_ERR_BAD_FILE		= 31)
	parameter (MPI_ERR_BASE			= 32)
	parameter (MPI_ERR_CONVERSION		= 33)
	parameter (MPI_ERR_DISP			= 34)
	parameter (MPI_ERR_DUP_DATAREP		= 35)
	parameter (MPI_ERR_FILE_EXISTS		= 36)
	parameter (MPI_ERR_FILE_IN_USE		= 37)
	parameter (MPI_ERR_FILE			= 38)
	parameter (MPI_ERR_INFO_KEY		= 39)
	parameter (MPI_ERR_INFO_NOKEY		= 40)
	parameter (MPI_ERR_INFO_VALUE		= 41)
	parameter (MPI_ERR_INFO			= 42)
	parameter (MPI_ERR_IO			= 43)
	parameter (MPI_ERR_KEYVAL		= 44)
	parameter (MPI_ERR_LOCKTYPE		= 45)
	parameter (MPI_ERR_NAME			= 46)
	parameter (MPI_ERR_NO_MEM		= 47)
	parameter (MPI_ERR_NOT_SAME		= 48)
	parameter (MPI_ERR_NO_SPACE		= 49)
	parameter (MPI_ERR_NO_SUCH_FILE		= 50)
	parameter (MPI_ERR_PORT			= 51)
	parameter (MPI_ERR_QUOTA		= 52)
	parameter (MPI_ERR_READ_ONLY		= 53)
	parameter (MPI_ERR_RMA_CONFLICT		= 54)
	parameter (MPI_ERR_RMA_SYNC		= 55)
	parameter (MPI_ERR_SERVICE		= 56)
	parameter (MPI_ERR_SIZE			= 57)
	parameter (MPI_ERR_SPAWN		= 58)
	parameter (MPI_ERR_UNSUPPORTED_DATAREP	= 59)
	parameter (MPI_ERR_UNSUPPORTED_OPERATION= 60)
	parameter (MPI_ERR_WIN			= 61)

	integer MPI_ERR_LASTCODE
	parameter (MPI_ERR_LASTCODE		= 100)


! Permanent keyvals

	integer MPI_KEYVAL_INVALID
	integer MPI_TAG_UB
	integer MPI_HOST
	integer MPI_IO
	integer MPI_WTIME_IS_GLOBAL

	parameter (MPI_KEYVAL_INVALID	= 0)
	parameter (MPI_TAG_UB		= 5)
	parameter (MPI_HOST		= 6)
	parameter (MPI_IO		= 7)
	parameter (MPI_WTIME_IS_GLOBAL	= 8)

! Results of the compare operations

	integer MPI_IDENT
	integer MPI_CONGRUENT
	integer MPI_SIMILAR
	integer MPI_UNEQUAL

	parameter (MPI_IDENT		= 0)
	parameter (MPI_CONGRUENT	= 1)
	parameter (MPI_SIMILAR		= 2)
	parameter (MPI_UNEQUAL		= 3)

! Topology types

	integer MPI_GRAPH
	integer MPI_CART

	parameter (MPI_GRAPH	= 1)
	parameter (MPI_CART	= 2)

! Misc constants

	integer MPI_MAX_PROCESSOR_NAME
	parameter (MPI_MAX_PROCESSOR_NAME = 256)

	integer MPI_MAX_ERROR_STRING
	parameter (MPI_MAX_ERROR_STRING = 256)

	integer MPI_BSEND_OVERHEAD
	parameter (MPI_BSEND_OVERHEAD = 32)

	integer MPI_UNDEFINED
	parameter (MPI_UNDEFINED = -3)

	integer MPI_ANY_SOURCE
	parameter (MPI_ANY_SOURCE = -2)

	integer MPI_PROC_NULL
	parameter (MPI_PROC_NULL = -1)

	integer MPI_ANY_TAG
	parameter (MPI_ANY_TAG = -1)

! Misc Fortran declarations

	integer MPI_BOTTOM
!	pointer (MPI_BOTTOM_PTR, MPI_BOTTOM)
!	data MPI_BOTTOM_PTR / 0 /

	external MPI_NULL_COPY_FN, MPI_NULL_DELETE_FN, MPI_DUP_FN

	double precision MPI_WTIME, MPI_WTICK
	external MPI_WTIME, MPI_WTICK

! MPI-2

	integer MPI_FUNDAMENTAL
	parameter (MPI_FUNDAMENTAL = -1)

! Kind values for MPI-2

	integer MPI_OFFSET_KIND
	parameter (MPI_OFFSET_KIND = 8)

