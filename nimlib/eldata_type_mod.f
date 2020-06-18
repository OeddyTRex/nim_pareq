c-----------------------------------------------------------------------
c     file eldata_type_mod.f
c     this module holds type definitions and routines for allocating
c     data organized by 2D element.  the primary purpose is for plotting
c     discontinuous fields.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0.  eldata_type_mod
c     1.  r4eldata_alloc
c     2.  r4eldata_dealloc
c-----------------------------------------------------------------------
c     subprogram 0. eldata_type_mod defines data structures for
c     packaging element-based fields.
c-----------------------------------------------------------------------
      MODULE eldata_type_mod
      USE local
      IMPLICIT NONE

      TYPE :: r4eldata_type
        CHARACTER(8) :: eltype
        CHARACTER(2) :: coordlabels
        INTEGER(i4) :: nele,nptperel,nfields
        INTEGER(i4), DIMENSION(:,:,:), POINTER :: connect
        REAL(r4), DIMENSION(:,:,:), POINTER :: coords
        REAL(r4), DIMENSION(:,:,:), POINTER :: data
        CHARACTER(16), DIMENSION(:), POINTER :: datalabels
      END TYPE r4eldata_type
c-----------------------------------------------------------------------
c     subprogram name interfaces
c-----------------------------------------------------------------------
      INTERFACE eldata_alloc
        MODULE PROCEDURE r4eldata_alloc
      END INTERFACE

      INTERFACE eldata_dealloc
        MODULE PROCEDURE r4eldata_dealloc
      END INTERFACE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. r4eldata_alloc
c     allocate the arrays and define intra-element connectivity
c     according to the passed element type parameter and basis
c     specification.
c
c     note that coordinates are not loaded here, so they must be loaded
c     in the same order as defined by the connectivity.
c
c     points have unique labels starting from iptst input, which
c     defines the starting index for this block.  iptst output is
c     the starting index for the next block.  if renumel is true
c     each element uses its own separate point indices, so the running
c     total index is not kept.
c-----------------------------------------------------------------------
      SUBROUTINE r4eldata_alloc(databl,eltype,poly_deg,nfields,nele,
     $                          iptst,renumel)

      TYPE(r4eldata_type), INTENT(OUT) :: databl
      CHARACTER(*), INTENT(IN) :: eltype
      INTEGER(i4), INTENT(IN) :: poly_deg,nfields,nele
      INTEGER(i4), INTENT(INOUT) :: iptst
      LOGICAL, INTENT(IN) :: renumel

      INTEGER(i4) :: nsubel,ix,iy,iel,isub,ipll,ipul
c-----------------------------------------------------------------------
c     copy parameters.
c-----------------------------------------------------------------------
      databl%eltype=eltype
      databl%nele=nele
      databl%nfields=nfields
      ALLOCATE(databl%datalabels(nfields))
c-----------------------------------------------------------------------
c     select element type and allocate arrays.
c-----------------------------------------------------------------------
      SELECT CASE(eltype)
      CASE ("quad")  !  conventional quadrilateral
        nsubel=poly_deg**2
        databl%nptperel=(poly_deg+1)**2
        ALLOCATE(databl%coords(2,databl%nptperel,nele))
        ALLOCATE(databl%data(nfields,databl%nptperel,nele))
        ALLOCATE(databl%connect(4,nsubel,nele))
        DO iel=1,nele
          isub=1
          DO iy=1,poly_deg
            ipll=(iy-1)*(poly_deg+1)+iptst
            ipul=iy*(poly_deg+1)+iptst
            DO ix=1,poly_deg
              databl%connect(:,isub,iel)=(/ipll+ix-1_i4,ipll+ix,
     $                                     ipul+ix,ipul+ix-1_i4/)
              isub=isub+1
            ENDDO
          ENDDO
          IF (.NOT.renumel) iptst=iptst+databl%nptperel
        ENDDO
      CASE DEFAULT
        CALL nim_stop("R4eldata_alloc: "//eltype//
     $                " element is not supported.")
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE r4eldata_alloc
c-----------------------------------------------------------------------
c     subprogram 2. r4eldata_dealloc
c     deallocate the r4eldata_type arrays.
c-----------------------------------------------------------------------
      SUBROUTINE r4eldata_dealloc(databl)

      TYPE(r4eldata_type), INTENT(INOUT) :: databl

      DEALLOCATE(databl%datalabels)
      DEALLOCATE(databl%coords)
      DEALLOCATE(databl%data)
      DEALLOCATE(databl%connect)

      RETURN
      END SUBROUTINE r4eldata_dealloc
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE eldata_type_mod

