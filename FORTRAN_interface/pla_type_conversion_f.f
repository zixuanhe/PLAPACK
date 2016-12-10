      subroutine pla_type_conversion_init_f( type_list )
      
      implicit none

      include "mpif.h"

      integer type_list( * )

      type_list(1) = MPI_CHARACTER
      type_list(2) = MPI_INTEGER
      type_list(3) = MPI_REAL
      type_list(4) = MPI_DOUBLE_PRECISION
      type_list(5) = MPI_COMPLEX
      type_list(6) = MPI_DOUBLE_COMPLEX

      return
      end



      subroutine PLA_TYPE_CONVERSION_CRAY_INIT_F( type_list )
      
      implicit none

      include "mpif.h"

      integer type_list( * )

      type_list(1) = MPI_CHARACTER
      type_list(2) = MPI_INTEGER
      type_list(3) = -1
      type_list(4) = MPI_REAL
      type_list(5) = -1
      type_list(6) = MPI_COMPLEX

      return
      end

