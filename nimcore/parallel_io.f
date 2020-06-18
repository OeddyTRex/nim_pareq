c-----------------------------------------------------------------------
c     file parallel_io.f
c
c     routines for sending rblocks and seams from one processor to another
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     1.  send_rblock
c     2.  send_tblock
c     3.  send_seam
c     4.  bcast_str (routines 4-6 are MPI bug-fixes,
c     5.  send_str     can eventually be deleted)
c     6.  recv_str
c     7.  gather_rblock
c     8.  copy_bicube
c     9.  copy_lagr_quad_2D
c     10. copy_modal_quad_2D
c     11. gather_lagr_quad
c     12. gather_modal_quad
c     13. trim_rblock
c     14. trim_lagr_quad
c     15. trim_modal_quad
c     16. gather_tblock
c     17. copy_tri_linear_2D
c     18. gather_tri_linear
c     19. trim_tblock
c     20. trim_tri_linear

c-----------------------------------------------------------------------
c     subprogram 1. send_rblock:  send a rblock from sendproc to recvproc
c                                 receiver allocates space as it comes in
c-----------------------------------------------------------------------
      SUBROUTINE send_rblock(sendproc,rbsend,recvproc,rbrecv)
      USE local
      USE rblock_type_mod
      USE mpi_nim
      USE pardata
      IMPLICIT NONE
      INTEGER(i4), INTENT(IN) :: sendproc,recvproc
      TYPE (rblock_type), INTENT(IN) :: rbsend
      TYPE (rblock_type), INTENT(INOUT) :: rbrecv
      integer(i4) :: status(mpi_status_size)

      INTEGER(i4) :: notify,ierror,i

c sender code
c wait for request from recvproc
c send all data in rbsend to recvproc

      IF (node == sendproc) then

        CALL mpi_recv(notify,1,mpi_nim_int,recvproc,0,
     $       mpi_comm_world,status,ierror)

        CALL send_str(rbsend%name,64,recvproc)
        CALL mpi_send(rbsend%id,1,mpi_nim_int,recvproc,0,
     $       mpi_comm_world,ierror)
        CALL mpi_send(rbsend%mx,1,mpi_nim_int,recvproc,0,
     $       mpi_comm_world,ierror)
        CALL mpi_send(rbsend%my,1,mpi_nim_int,recvproc,0,
     $       mpi_comm_world,ierror)
        CALL mpi_send(rbsend%degenerate,1,mpi_nim_logical,
     $       recvproc,0,mpi_comm_world,ierror)

        CALL send_lagr_quad_2D(rbsend%rz)
        CALL send_lagr_quad_2D(rbsend%be_eq)
        CALL send_lagr_quad(rbsend%be)
        CALL send_lagr_quad_2D(rbsend%ja_eq)
        CALL send_lagr_quad_2D(rbsend%ve_eq)
        CALL send_lagr_quad(rbsend%ve)
        CALL send_lagr_quad_2D(rbsend%pres_eq)
        CALL send_lagr_quad(rbsend%pres)
        CALL send_lagr_quad_2D(rbsend%prese_eq)
        CALL send_lagr_quad(rbsend%prese)
        CALL send_lagr_quad_2D(rbsend%nd_eq)
        CALL send_lagr_quad(rbsend%nd)
        CALL send_lagr_quad_2D(rbsend%diff_shape)
        CALL send_lagr_quad(rbsend%conc)
        CALL send_lagr_quad(rbsend%tele)
        CALL send_lagr_quad(rbsend%tion)

c       not needed for diffusive corrections.
c       CALL send_modal_quad(rbsend%auxb)
c       CALL send_modal_quad(rbsend%auxv)

c receiver code
c tell sendproc I'm ready
c recv all data into rbrecv, allocate space as it comes in

      ELSE IF (node == recvproc) THEN

        CALL mpi_send(notify,1,mpi_nim_int,sendproc,0,
     $       mpi_comm_world,ierror)

        CALL recv_str(rbrecv%name,64,sendproc)
        CALL mpi_recv(rbrecv%id,1,mpi_nim_int,sendproc,0,
     $       mpi_comm_world,status,ierror)
        CALL mpi_recv(rbrecv%mx,1,mpi_nim_int,sendproc,0,
     $       mpi_comm_world,status,ierror)
        CALL mpi_recv(rbrecv%my,1,mpi_nim_int,sendproc,0,
     $       mpi_comm_world,status,ierror)
        CALL mpi_recv(rbrecv%degenerate,1,mpi_nim_logical,
     $       sendproc,0,mpi_comm_world,status,ierror)

        CALL receive_lagr_quad_2D(rbrecv%rz)
        CALL receive_lagr_quad_2D(rbrecv%be_eq)
        CALL receive_lagr_quad(rbrecv%be)
        CALL receive_lagr_quad_2D(rbrecv%ja_eq)
        CALL receive_lagr_quad_2D(rbrecv%ve_eq)
        CALL receive_lagr_quad(rbrecv%ve)
        CALL receive_lagr_quad_2D(rbrecv%pres_eq)
        CALL receive_lagr_quad(rbrecv%pres)
        CALL receive_lagr_quad_2D(rbrecv%prese_eq)
        CALL receive_lagr_quad(rbrecv%prese)
        CALL receive_lagr_quad_2D(rbrecv%nd_eq)
        CALL receive_lagr_quad(rbrecv%nd)
        CALL receive_lagr_quad_2D(rbrecv%diff_shape)
        CALL receive_lagr_quad(rbrecv%conc)
        CALL receive_lagr_quad(rbrecv%tele)
        CALL receive_lagr_quad(rbrecv%tion)

c       not needed for diffusive corrections.
c       CALL receive_modal_quad(rbrecv%auxb)
c       CALL receive_modal_quad(rbrecv%auxv)

      ENDIF

      RETURN

c internal subroutines used to pass rblock data structures:

      CONTAINS  

c   send a bicube structure.

        SUBROUTINE send_bicube(bc)

        TYPE(bicube_type), INTENT(IN) :: bc

        INTEGER(i4) :: isize

        CALL send_str(bc%name,6,recvproc)
        CALL mpi_send(bc%nqty,1,mpi_nim_int,recvproc,0,
     $       mpi_comm_world,ierror)
        DO i = 1,bc%nqty
          CALL send_str(bc%title(i),6,recvproc)
        ENDDO
        isize = SIZE(bc%fs)
        CALL mpi_send(bc%xs,SIZE(bc%xs),
     $       mpi_nim_real,recvproc,0,mpi_comm_world,ierror)
        CALL mpi_send(bc%ys,SIZE(bc%ys),
     $       mpi_nim_real,recvproc,0,mpi_comm_world,ierror)
        CALL mpi_send(bc%fs,isize,
     $       mpi_nim_real,recvproc,0,mpi_comm_world,ierror)
        CALL mpi_send(bc%fsx,isize,
     $       mpi_nim_real,recvproc,0,mpi_comm_world,ierror)
        CALL mpi_send(bc%fsy,isize,
     $       mpi_nim_real,recvproc,0,mpi_comm_world,ierror)
        CALL mpi_send(bc%fsxy,isize,
     $       mpi_nim_real,recvproc,0,mpi_comm_world,ierror)

        RETURN
        END SUBROUTINE send_bicube

c   send a lagr_quad structure.  nqty now serves as a flag: if it is 0,
c   the structure is not active, and if it is <0, only interior data
c   is used for discontinuous fields.

        SUBROUTINE send_lagr_quad(laq)

        TYPE(lagr_quad_type), INTENT(IN) :: laq

        CALL send_str(laq%name,6,recvproc)
        IF (ALLOCATED(laq%fs)) THEN  !  continuous field
          CALL mpi_send(laq%nqty,1,mpi_nim_int,recvproc,0,
     $         mpi_comm_world,ierror)
        ELSE                          !  discontinuous field
          CALL mpi_send(-laq%nqty,1,mpi_nim_int,recvproc,0,
     $         mpi_comm_world,ierror)
        ENDIF
        IF (laq%nqty==0) RETURN
        CALL mpi_send(laq%nfour,1,mpi_nim_int,recvproc,0,
     $       mpi_comm_world,ierror)
        CALL mpi_send(laq%n_side,1,mpi_nim_int,recvproc,0,
     $       mpi_comm_world,ierror)
        CALL mpi_send(laq%n_int,1,mpi_nim_int,recvproc,0,
     $       mpi_comm_world,ierror)
        DO i = 1,laq%nqty
          CALL send_str(laq%title(i),6,recvproc)
        ENDDO
        IF (ALLOCATED(laq%fs))
     $    CALL mpi_send(laq%fs,SIZE(laq%fs),
     $         mpi_nim_comp,recvproc,0,mpi_comm_world,ierror)
        IF (laq%n_side>0) THEN
          CALL mpi_send(laq%fsh,SIZE(laq%fsh),
     $         mpi_nim_comp,recvproc,0,mpi_comm_world,ierror)
          CALL mpi_send(laq%fsv,SIZE(laq%fsv),
     $         mpi_nim_comp,recvproc,0,mpi_comm_world,ierror)
        ENDIF
        IF (laq%n_int>0) THEN
          CALL mpi_send(laq%fsi,SIZE(laq%fsi),
     $         mpi_nim_comp,recvproc,0,mpi_comm_world,ierror)
        ENDIF

        RETURN
        END SUBROUTINE send_lagr_quad

c   send a 2D lagr_quad structure.

        SUBROUTINE send_lagr_quad_2D(laq)

        TYPE(lagr_quad_2D_type), INTENT(IN) :: laq

        CALL send_str(laq%name,6,recvproc)
        IF (ALLOCATED(laq%fs)) THEN  !  continuous field
          CALL mpi_send(laq%nqty,1,mpi_nim_int,recvproc,0,
     $         mpi_comm_world,ierror)
        ELSE                          !  discontinuous field
          CALL mpi_send(-laq%nqty,1,mpi_nim_int,recvproc,0,
     $         mpi_comm_world,ierror)
        ENDIF
        IF (laq%nqty==0) RETURN
        CALL mpi_send(laq%n_side,1,mpi_nim_int,recvproc,0,
     $       mpi_comm_world,ierror)
        CALL mpi_send(laq%n_int,1,mpi_nim_int,recvproc,0,
     $       mpi_comm_world,ierror)
        DO i = 1,laq%nqty
          CALL send_str(laq%title(i),6,recvproc)
        ENDDO
        IF (ALLOCATED(laq%fs))
     $    CALL mpi_send(laq%fs,SIZE(laq%fs),
     $         mpi_nim_real,recvproc,0,mpi_comm_world,ierror)
        IF (laq%n_side>0) THEN
          CALL mpi_send(laq%fsh,SIZE(laq%fsh),
     $         mpi_nim_real,recvproc,0,mpi_comm_world,ierror)
          CALL mpi_send(laq%fsv,SIZE(laq%fsv),
     $         mpi_nim_real,recvproc,0,mpi_comm_world,ierror)
        ENDIF
        IF (laq%n_int>0) THEN
          CALL mpi_send(laq%fsi,SIZE(laq%fsi),
     $         mpi_nim_real,recvproc,0,mpi_comm_world,ierror)
        ENDIF

        RETURN
        END SUBROUTINE send_lagr_quad_2D

c   send a modal_quad structure.

        SUBROUTINE send_modal_quad(modq)

        TYPE(modal_quad_type), INTENT(IN) :: modq

        CALL send_str(modq%name,6,recvproc)
        CALL mpi_send(modq%nqty,1,mpi_nim_int,recvproc,0,
     $       mpi_comm_world,ierror)
        IF (modq%nqty==0) RETURN
        CALL mpi_send(modq%nfour,1,mpi_nim_int,recvproc,0,
     $       mpi_comm_world,ierror)
        CALL mpi_send(modq%pd,1,mpi_nim_int,recvproc,0,
     $       mpi_comm_world,ierror)
        CALL mpi_send(modq%pdmin,1,mpi_nim_int,recvproc,0,
     $       mpi_comm_world,ierror)
        CALL mpi_send(modq%pdmax,1,mpi_nim_int,recvproc,0,
     $       mpi_comm_world,ierror)
        CALL mpi_send(modq%n_int,1,mpi_nim_int,recvproc,0,
     $       mpi_comm_world,ierror)
        DO i = 1,modq%nqty
          CALL send_str(modq%title(i),6,recvproc)
        ENDDO
        CALL mpi_send(modq%fsi,SIZE(modq%fsi),
     $       mpi_nim_comp,recvproc,0,mpi_comm_world,ierror)

        RETURN
        END SUBROUTINE send_modal_quad

c   send a 2D modal_quad structure.

        SUBROUTINE send_modal_quad_2D(modq)

        TYPE(modal_quad_2D_type), INTENT(IN) :: modq

        CALL send_str(modq%name,6,recvproc)
        CALL mpi_send(modq%nqty,1,mpi_nim_int,recvproc,0,
     $       mpi_comm_world,ierror)
        IF (modq%nqty==0) RETURN
        CALL mpi_send(modq%pd,1,mpi_nim_int,recvproc,0,
     $       mpi_comm_world,ierror)
        CALL mpi_send(modq%pdmin,1,mpi_nim_int,recvproc,0,
     $       mpi_comm_world,ierror)
        CALL mpi_send(modq%pdmax,1,mpi_nim_int,recvproc,0,
     $       mpi_comm_world,ierror)
        CALL mpi_send(modq%n_int,1,mpi_nim_int,recvproc,0,
     $       mpi_comm_world,ierror)
        DO i = 1,modq%nqty
          CALL send_str(modq%title(i),6,recvproc)
        ENDDO
        CALL mpi_send(modq%fsi,SIZE(modq%fsi),
     $       mpi_nim_real,recvproc,0,mpi_comm_world,ierror)

        RETURN
        END SUBROUTINE send_modal_quad_2D

c   receive a bicube structure.

        SUBROUTINE receive_bicube(bc)

        TYPE(bicube_type), INTENT(INOUT) :: bc

        INTEGER(i4) :: isize,nq
        CHARACTER(6) :: name

        CALL recv_str(name,6,sendproc)
        CALL mpi_recv(nq,1,mpi_nim_int,sendproc,0,
     $       mpi_comm_world,status,ierror)
        CALL bicube_alloc(bc,rbrecv%mx,rbrecv%my,nq)
        bc%name=name
        DO i = 1,bc%nqty
          CALL recv_str(bc%title(i),6,sendproc)
        ENDDO
        isize = SIZE(bc%fs)
        CALL mpi_recv(bc%xs,SIZE(bc%xs),
     $       mpi_nim_real,sendproc,0,
     $       mpi_comm_world,status,ierror)
        CALL mpi_recv(bc%ys,SIZE(bc%ys),
     $       mpi_nim_real,sendproc,0,
     $       mpi_comm_world,status,ierror)
        CALL mpi_recv(bc%fs,isize,mpi_nim_real,sendproc,0,
     $       mpi_comm_world,status,ierror)
        CALL mpi_recv(bc%fsx,isize,mpi_nim_real,sendproc,0,
     $       mpi_comm_world,status,ierror)
        CALL mpi_recv(bc%fsy,isize,mpi_nim_real,sendproc,0,
     $       mpi_comm_world,status,ierror)
        CALL mpi_recv(bc%fsxy,isize,mpi_nim_real,sendproc,0,
     $       mpi_comm_world,status,ierror)

        RETURN
        END SUBROUTINE receive_bicube

c   receive a lagr_quad structure.  the communicated nq serves as
c   flag for whether this structure is used and whether the field
c   is continuous.

        SUBROUTINE receive_lagr_quad(laq)

        TYPE(lagr_quad_type), INTENT(INOUT) :: laq

        INTEGER(i4) :: nq,nf,ns,ni
        CHARACTER(6) :: name

        CALL recv_str(name,6,sendproc)
        laq%name=name
        CALL mpi_recv(nq,1,mpi_nim_int,sendproc,0,
     $       mpi_comm_world,status,ierror)
        IF (nq==0) RETURN
        CALL mpi_recv(nf,1,mpi_nim_int,sendproc,0,
     $       mpi_comm_world,status,ierror)
        CALL mpi_recv(ns,1,mpi_nim_int,sendproc,0,
     $       mpi_comm_world,status,ierror)
        CALL mpi_recv(ni,1,mpi_nim_int,sendproc,0,
     $       mpi_comm_world,status,ierror)
        IF (nq>0) THEN                       !  continuous field
          CALL lagr_quad_alloc(laq,rbrecv%mx,rbrecv%my,
     $         nq,nf,ns+1_i4)
        ELSE
          CALL lagr_disc_alloc(laq,rbrecv%mx,rbrecv%my,
     $         -nq,nf,NINT(SQRT(REAL(ni)))-1_i4)
        ENDIF
        DO i = 1,laq%nqty
          CALL recv_str(laq%title(i),6,sendproc)
        ENDDO
        IF (ALLOCATED(laq%fs))
     $    CALL mpi_recv(laq%fs,SIZE(laq%fs),
     $         mpi_nim_comp,sendproc,0,
     $         mpi_comm_world,status,ierror)
        IF (laq%n_side>0) THEN
          CALL mpi_recv(laq%fsh,SIZE(laq%fsh),
     $         mpi_nim_comp,sendproc,0,
     $         mpi_comm_world,status,ierror)
          CALL mpi_recv(laq%fsv,SIZE(laq%fsv),
     $         mpi_nim_comp,sendproc,0,
     $         mpi_comm_world,status,ierror)
        ENDIF
        IF (laq%n_int>0) THEN
          CALL mpi_recv(laq%fsi,SIZE(laq%fsi),
     $         mpi_nim_comp,sendproc,0,
     $         mpi_comm_world,status,ierror)
        ENDIF

        RETURN
        END SUBROUTINE receive_lagr_quad

c   receive a 2D lagr_quad structure.

        SUBROUTINE receive_lagr_quad_2D(laq)

        TYPE(lagr_quad_2D_type), INTENT(INOUT) :: laq

        INTEGER(i4) :: nq,ns,ni
        CHARACTER(6) :: name

        CALL recv_str(name,6,sendproc)
        laq%name=name
        CALL mpi_recv(nq,1,mpi_nim_int,sendproc,0,
     $       mpi_comm_world,status,ierror)
        IF (nq==0) RETURN
        CALL mpi_recv(ns,1,mpi_nim_int,sendproc,0,
     $       mpi_comm_world,status,ierror)
        CALL mpi_recv(ni,1,mpi_nim_int,sendproc,0,
     $       mpi_comm_world,status,ierror)
        IF (nq>0) THEN                       !  continuous field
          CALL lagr_quad_alloc(laq,rbrecv%mx,rbrecv%my,
     $         nq,ns+1_i4)
        ELSE
          CALL lagr_disc_alloc(laq,rbrecv%mx,rbrecv%my,
     $         -nq,NINT(SQRT(REAL(ni)))-1_i4)
        ENDIF
        DO i = 1,laq%nqty
          CALL recv_str(laq%title(i),6,sendproc)
        ENDDO
        IF (ALLOCATED(laq%fs))
     $    CALL mpi_recv(laq%fs,SIZE(laq%fs),
     $         mpi_nim_real,sendproc,0,
     $         mpi_comm_world,status,ierror)
        IF (laq%n_side>0) THEN
          CALL mpi_recv(laq%fsh,SIZE(laq%fsh),
     $         mpi_nim_real,sendproc,0,
     $         mpi_comm_world,status,ierror)
          CALL mpi_recv(laq%fsv,SIZE(laq%fsv),
     $         mpi_nim_real,sendproc,0,
     $         mpi_comm_world,status,ierror)
        ENDIF
        IF (laq%n_int>0) THEN
          CALL mpi_recv(laq%fsi,SIZE(laq%fsi),
     $         mpi_nim_real,sendproc,0,
     $         mpi_comm_world,status,ierror)
        ENDIF

        RETURN
        END SUBROUTINE receive_lagr_quad_2D

c   receive a modal_quad structure.

        SUBROUTINE receive_modal_quad(modq)

        TYPE(modal_quad_type), INTENT(INOUT) :: modq

        INTEGER(i4) :: nq,nf,pd,pdmin,pdmax,ni
        CHARACTER(6) :: name

        CALL recv_str(name,6,sendproc)
        modq%name=name
        CALL mpi_recv(nq,1,mpi_nim_int,sendproc,0,
     $       mpi_comm_world,status,ierror)
        IF (nq==0) RETURN
        CALL mpi_recv(nf,1,mpi_nim_int,sendproc,0,
     $       mpi_comm_world,status,ierror)
        CALL mpi_recv(pd,1,mpi_nim_int,sendproc,0,
     $       mpi_comm_world,status,ierror)
        CALL mpi_recv(pdmin,1,mpi_nim_int,sendproc,0,
     $       mpi_comm_world,status,ierror)
        CALL mpi_recv(pdmax,1,mpi_nim_int,sendproc,0,
     $       mpi_comm_world,status,ierror)
        CALL mpi_recv(ni,1,mpi_nim_int,sendproc,0,
     $       mpi_comm_world,status,ierror)
        CALL modal_disc_alloc(modq,rbrecv%mx,rbrecv%my,
     $       nq,nf,pd,pdmin,pdmax)
        DO i = 1,modq%nqty
          CALL recv_str(modq%title(i),6,sendproc)
        ENDDO
        CALL mpi_recv(modq%fsi,SIZE(modq%fsi),
     $       mpi_nim_comp,sendproc,0,
     $       mpi_comm_world,status,ierror)

        RETURN
        END SUBROUTINE receive_modal_quad

c   receive a 2D modal_quad structure.

        SUBROUTINE receive_modal_quad_2D(modq)

        TYPE(modal_quad_2D_type), INTENT(INOUT) :: modq

        INTEGER(i4) :: nq,pd,pdmin,pdmax,ni
        CHARACTER(6) :: name

        CALL recv_str(name,6,sendproc)
        modq%name=name
        CALL mpi_recv(nq,1,mpi_nim_int,sendproc,0,
     $       mpi_comm_world,status,ierror)
        IF (nq==0) RETURN
        CALL mpi_recv(pd,1,mpi_nim_int,sendproc,0,
     $       mpi_comm_world,status,ierror)
        CALL mpi_recv(pdmin,1,mpi_nim_int,sendproc,0,
     $       mpi_comm_world,status,ierror)
        CALL mpi_recv(pdmax,1,mpi_nim_int,sendproc,0,
     $       mpi_comm_world,status,ierror)
        CALL mpi_recv(ni,1,mpi_nim_int,sendproc,0,
     $       mpi_comm_world,status,ierror)
        CALL modal_disc_alloc(modq,rbrecv%mx,rbrecv%my,
     $       nq,pd,pdmin,pdmax)
        DO i = 1,modq%nqty
          CALL recv_str(modq%title(i),6,sendproc)
        ENDDO
        CALL mpi_recv(modq%fsi,SIZE(modq%fsi),
     $       mpi_nim_real,sendproc,0,
     $       mpi_comm_world,status,ierror)

        RETURN
        END SUBROUTINE receive_modal_quad_2D

c   end of rblock interal routines.

      END SUBROUTINE send_rblock


c-----------------------------------------------------------------------
c     subprogram 2. send_tblock:  send a tblock from sendproc to recvproc
c                                 receiver allocates space as it comes in
c-----------------------------------------------------------------------
      SUBROUTINE send_tblock(sendproc,tbsend,recvproc,tbrecv)
      USE local
      USE tblock_type_mod
      USE mpi_nim
      USE pardata
      IMPLICIT NONE
      INTEGER(i4), INTENT(IN) :: sendproc,recvproc
      TYPE (tblock_type), INTENT(IN) :: tbsend
      TYPE (tblock_type), INTENT(INOUT) :: tbrecv
      integer(i4) :: status(mpi_status_size)

      INTEGER(i4) :: notify,ierror,i,ivert,nsize

c sender code
c wait for request from recvproc
c send all data in tbsend to recvproc

      IF (node == sendproc) then

        CALL mpi_recv(notify,1,mpi_nim_int,recvproc,0,
     $       mpi_comm_world,status,ierror)

        CALL send_str(tbsend%name,64,recvproc)
        CALL mpi_send(tbsend%id,1,mpi_nim_int,recvproc,0,
     $       mpi_comm_world,ierror)
        CALL mpi_send(tbsend%tgeom%mvert,1,mpi_nim_int,recvproc,0,
     $       mpi_comm_world,ierror)
        CALL mpi_send(tbsend%tgeom%mcell,1,mpi_nim_int,recvproc,0,
     $       mpi_comm_world,ierror)

        CALL mpi_send(tbsend%tgeom%xs,SIZE(tbsend%tgeom%xs),
     $       mpi_nim_real,recvproc,0,mpi_comm_world,ierror)
        CALL mpi_send(tbsend%tgeom%ys,SIZE(tbsend%tgeom%ys),
     $       mpi_nim_real,recvproc,0,mpi_comm_world,ierror)

        CALL mpi_send(tbsend%tgeom%vertex,SIZE(tbsend%tgeom%vertex),
     $       mpi_nim_int,recvproc,0,mpi_comm_world,ierror)

        DO ivert=0,tbsend%tgeom%mvert
          CALL mpi_send(SIZE(tbsend%tgeom%neighbor(ivert)%vertex),1,
     $         mpi_nim_int,recvproc,0,mpi_comm_world,ierror)
          CALL mpi_send(tbsend%tgeom%neighbor(ivert)%vertex,
     $         SIZE(tbsend%tgeom%neighbor(ivert)%vertex),
     $         mpi_nim_int,recvproc,0,mpi_comm_world,ierror)
        ENDDO

        CALL send_tri_linear_2D(tbsend%be_eq)
        CALL send_tri_linear(tbsend%be)
        CALL send_tri_linear_2D(tbsend%ja_eq)
        CALL send_tri_linear_2D(tbsend%ve_eq)
        CALL send_tri_linear(tbsend%ve)
        CALL send_tri_linear_2D(tbsend%pres_eq)
        CALL send_tri_linear(tbsend%pres)
        CALL send_tri_linear_2D(tbsend%prese_eq)
        CALL send_tri_linear(tbsend%prese)
        CALL send_tri_linear_2D(tbsend%nd_eq)
        CALL send_tri_linear(tbsend%nd)
        CALL send_tri_linear_2D(tbsend%diff_shape)
        CALL send_tri_linear(tbsend%conc)
        CALL send_tri_linear(tbsend%tele)
        CALL send_tri_linear(tbsend%tion)

c receiver code
c tell sendproc I'm ready
c recv all data into tbrecv, allocate space as it comes in

      ELSE IF (node == recvproc) THEN

        CALL mpi_send(notify,1,mpi_nim_int,sendproc,0,
     $       mpi_comm_world,ierror)

        CALL recv_str(tbrecv%name,64,sendproc)
        CALL mpi_recv(tbrecv%id,1,mpi_nim_int,sendproc,0,
     $       mpi_comm_world,status,ierror)
        CALL mpi_recv(tbrecv%tgeom%mvert,1,mpi_nim_int,sendproc,0,
     $       mpi_comm_world,status,ierror)
        CALL mpi_recv(tbrecv%tgeom%mcell,1,mpi_nim_int,sendproc,0,
     $       mpi_comm_world,status,ierror)

        CALL tri_linear_geom_alloc(tbrecv%tgeom,
     $       tbrecv%tgeom%mvert,tbrecv%tgeom%mcell)

        CALL mpi_recv(tbrecv%tgeom%xs,SIZE(tbrecv%tgeom%xs),
     $       mpi_nim_real,sendproc,0,
     $       mpi_comm_world,status,ierror)
        CALL mpi_recv(tbrecv%tgeom%ys,SIZE(tbrecv%tgeom%ys),
     $       mpi_nim_real,sendproc,0,
     $       mpi_comm_world,status,ierror)

        CALL mpi_recv(tbrecv%tgeom%vertex,SIZE(tbrecv%tgeom%vertex),
     $       mpi_nim_int,sendproc,0,
     $       mpi_comm_world,status,ierror)

        DO ivert=0,tbrecv%tgeom%mvert
          CALL mpi_recv(nsize,1,mpi_nim_int,sendproc,0,
     $         mpi_comm_world,status,ierror)
          ALLOCATE(tbrecv%tgeom%neighbor(ivert)%vertex(0:nsize-1))
          CALL mpi_recv(tbrecv%tgeom%neighbor(ivert)%vertex,
     $         SIZE(tbrecv%tgeom%neighbor(ivert)%vertex),
     $         mpi_nim_int,sendproc,0,
     $         mpi_comm_world,status,ierror)
        ENDDO

        tbrecv%mvert=tbrecv%tgeom%mvert
        tbrecv%mcell=tbrecv%tgeom%mcell

        CALL receive_tri_linear_2D(tbrecv%be_eq)
        CALL receive_tri_linear(tbrecv%be)
        CALL receive_tri_linear_2D(tbrecv%ja_eq)
        CALL receive_tri_linear_2D(tbrecv%ve_eq)
        CALL receive_tri_linear(tbrecv%ve)
        CALL receive_tri_linear_2D(tbrecv%pres_eq)
        CALL receive_tri_linear(tbrecv%pres)
        CALL receive_tri_linear_2D(tbrecv%prese_eq)
        CALL receive_tri_linear(tbrecv%prese)
        CALL receive_tri_linear_2D(tbrecv%nd_eq)
        CALL receive_tri_linear(tbrecv%nd)
        CALL receive_tri_linear_2D(tbrecv%diff_shape)
        CALL receive_tri_linear(tbrecv%conc)
        CALL receive_tri_linear(tbrecv%tele)
        CALL receive_tri_linear(tbrecv%tion)

      ENDIF

      RETURN

c internal subroutines used to pass tblock data structures:

      CONTAINS  

c   send a tri_linear_2D (real) structure.

        SUBROUTINE send_tri_linear_2D(tl2)

        TYPE(tri_linear_2D_type), INTENT(IN) :: tl2

        CALL send_str(tl2%name,6,recvproc)
        CALL mpi_send(tl2%nqty,1,mpi_nim_int,recvproc,0,
     $       mpi_comm_world,ierror)
        DO i = 1,tl2%nqty
          CALL send_str(tl2%title(i),6,recvproc)
        ENDDO
        CALL mpi_send(tl2%fs,SIZE(tl2%fs),
     $       mpi_nim_real,recvproc,0,mpi_comm_world,ierror)

        RETURN
        END SUBROUTINE send_tri_linear_2D

c   send a tri_linear (3D complex) structure.

        SUBROUTINE send_tri_linear(tl)

        TYPE(tri_linear_type), INTENT(IN) :: tl

        CALL send_str(tl%name,6,recvproc)
        CALL mpi_send(tl%nqty,1,mpi_nim_int,recvproc,0,
     $       mpi_comm_world,ierror)
        CALL mpi_send(tl%nfour,1,mpi_nim_int,recvproc,0,
     $       mpi_comm_world,ierror)
        DO i = 1,tl%nqty
          CALL send_str(tl%title(i),6,recvproc)
        ENDDO
        CALL mpi_send(tl%fs,SIZE(tl%fs),
     $       mpi_nim_comp,recvproc,0,mpi_comm_world,ierror)

        RETURN
        END SUBROUTINE send_tri_linear

c   receive a tri_linear_2D (real) structure.

        SUBROUTINE receive_tri_linear_2D(tl2)

        TYPE(tri_linear_2D_type), INTENT(INOUT) :: tl2

        INTEGER(i4) :: nq
        CHARACTER(6) :: name

        CALL recv_str(name,6,sendproc)
        CALL mpi_recv(nq,1,mpi_nim_int,sendproc,0,
     $       mpi_comm_world,status,ierror)
        CALL tri_linear_alloc(tl2,tbrecv%mvert,nq)
        tl2%name=name
        DO i = 1,tl2%nqty
          CALL recv_str(tl2%title(i),6,sendproc)
        ENDDO
        CALL mpi_recv(tl2%fs,SIZE(tl2%fs),
     $       mpi_nim_real,sendproc,0,mpi_comm_world,status,ierror)

        RETURN
        END SUBROUTINE receive_tri_linear_2D

c   receive a tri_linear (3D complex) structure.

        SUBROUTINE receive_tri_linear(tl)

        TYPE(tri_linear_type), INTENT(INOUT) :: tl

        INTEGER(i4) :: nq,nf
        CHARACTER(6) :: name

        CALL recv_str(name,6,sendproc)
        CALL mpi_recv(nq,1,mpi_nim_int,sendproc,0,
     $       mpi_comm_world,status,ierror)
        CALL mpi_recv(nf,1,mpi_nim_int,sendproc,0,
     $       mpi_comm_world,status,ierror)
        CALL tri_linear_alloc(tl,tbrecv%mvert,nq,nf)
        tl%name=name
        DO i = 1,tl%nqty
          CALL recv_str(tl%title(i),6,sendproc)
        ENDDO
        CALL mpi_recv(tl%fs,SIZE(tl%fs),
     $       mpi_nim_comp,sendproc,0,mpi_comm_world,status,ierror)

        RETURN
        END SUBROUTINE receive_tri_linear

c   end of tblock internal routines.

      END SUBROUTINE send_tblock


c-----------------------------------------------------------------------
c     subprogram 3. send_seam:  send a seam from sendproc to recvproc
c                               receiver allocates space as it comes in
c-----------------------------------------------------------------------
      SUBROUTINE send_seam(sendproc,seamsend,recvproc,seamrecv)
      USE local
      USE edge_type_mod
      USE mpi_nim
      USE pardata
      IMPLICIT NONE
      INTEGER(i4), INTENT(IN) :: sendproc,recvproc
      TYPE (edge_type), INTENT(IN) :: seamsend
      TYPE (edge_type), INTENT(INOUT) :: seamrecv
      integer(i4) :: status(mpi_status_size)

      INTEGER(i4) :: notify,ierror,i,np,ip

c sender code
c wait for request from recvproc
c send all data in seamsend to recvproc

      IF (node == sendproc) THEN

        CALL mpi_recv(notify,1,mpi_nim_int,recvproc,0,
     $       mpi_comm_world,status,ierror)

        CALL send_str(seamsend%name,64,recvproc)
c        CALL mpi_send(seamsend%name,64,mpi_nim_char,recvproc,0,
c     $       mpi_comm_world,ierror)
        CALL mpi_send(seamsend%id,1,mpi_nim_int,recvproc,0,
     $       mpi_comm_world,ierror)
        CALL mpi_send(seamsend%nvert,1,mpi_nim_int,recvproc,0,
     $       mpi_comm_world,ierror)

        DO i = 1,seamsend%nvert
          np = SIZE(seamsend%vertex(i)%ptr,2)
          CALL mpi_send(np,1,mpi_nim_int,recvproc,0,
     $         mpi_comm_world,ierror)
          DO ip = 1,np
            CALL mpi_send(seamsend%vertex(i)%ptr(1,ip),2,
     $           mpi_nim_int,recvproc,0,mpi_comm_world,ierror)
          ENDDO
          IF (seamsend%id > 0)
     $         CALL mpi_send(seamsend%vertex(i)%intxy,2,
     $         mpi_nim_int,recvproc,0,mpi_comm_world,ierror)
        ENDDO

        IF (seamsend%id == 0)
     $       CALL mpi_send(seamsend%excorner,seamsend%nvert,
     $       mpi_nim_logical,recvproc,0,mpi_comm_world,ierror)

c receiver code
c tell sendproc I'm ready
c recv all data into seamrecv, allocate space as it comes in

      ELSE IF (node == recvproc) THEN

        CALL mpi_send(notify,1,mpi_nim_int,sendproc,0,
     $       mpi_comm_world,ierror)

        CALL recv_str(seamrecv%name,64,sendproc)
c        CALL mpi_recv(seamrecv%name,64,mpi_nim_char,sendproc,0,
c     $       mpi_comm_world,status,ierror)
        CALL mpi_recv(seamrecv%id,1,mpi_nim_int,sendproc,0,
     $       mpi_comm_world,status,ierror)
        CALL mpi_recv(seamrecv%nvert,1,mpi_nim_int,sendproc,0,
     $       mpi_comm_world,status,ierror)

        ALLOCATE(seamrecv%vertex(seamrecv%nvert))

        DO i = 1,seamrecv%nvert
          CALL mpi_recv(np,1,mpi_nim_int,sendproc,0,
     $         mpi_comm_world,status,ierror)
          ALLOCATE(seamrecv%vertex(i)%ptr(2,np))
          DO ip = 1,np
            CALL mpi_recv(seamrecv%vertex(i)%ptr(1,ip),2,
     $           mpi_nim_int,sendproc,0,mpi_comm_world,status,ierror)
          ENDDO
          IF (seamrecv%id > 0) THEN
            CALL mpi_recv(seamrecv%vertex(i)%intxy,2,
     $           mpi_nim_int,sendproc,0,mpi_comm_world,status,ierror)
          ENDIF
        ENDDO

        IF (seamrecv%id == 0) THEN
          ALLOCATE(seamrecv%excorner(seamrecv%nvert))
          CALL mpi_recv(seamrecv%excorner,seamrecv%nvert,
     $         mpi_nim_logical,sendproc,0,
     $         mpi_comm_world,status,ierror)
        ENDIF

      ENDIF

      RETURN
      END SUBROUTINE send_seam


c-----------------------------------------------------------------------
c     subprograms 4-6. string message passing routines: 
c                      broadcast, send, receive a string from one 
c                        processor to another by treating them as integers !!
c                      should not have to call these from 
c                        send_rblock and send_seam
c                        if I can eventually get MPI characters data type
c                        to work correctly !!
c-----------------------------------------------------------------------

      SUBROUTINE bcast_str(str,n)
      USE local
      USE mpi_nim
      USE pardata
      character*(*) str
      INTEGER(i4) :: n,i,value,ierror

      do i = 1,n
        if (node == 0) value = ichar(str(i:i))
        CALL mpi_bcast(value,1,mpi_nim_int,0,mpi_comm_world,ierror)
        str(i:i) = char(value)
      enddo

      return
      end


      SUBROUTINE send_str(str,n,proc)
      USE local
      USE mpi_nim
      USE pardata
      character*(*) str
      INTEGER(i4) :: n,i,value,proc,ierror

      do i = 1,n
        value = ichar(str(i:i))
        CALL mpi_send(value,1,mpi_nim_int,proc,0,
     $       mpi_comm_world,ierror)
      enddo

      return
      end


      SUBROUTINE recv_str(str,n,proc)
      USE local
      USE mpi_nim
      USE pardata
      character*(*) str
      INTEGER(i4) :: n,i,value,proc,ierror
      integer(i4) :: status(mpi_status_size)

      do i = 1,n
        CALL mpi_recv(value,1,mpi_nim_int,proc,0,
     $       mpi_comm_world,status,ierror)
        str(i:i) = char(value)
      enddo

      return
      end

c-----------------------------------------------------------------------
c     subprogram 7. gather_rblock:  create an rblock with
c     complete information from all Fourier layers.
c-----------------------------------------------------------------------
      SUBROUTINE gather_rblock(rb_all,rb_loc,nmodes,nmodes_total)
      USE local
      USE rblock_type_mod
      USE pardata
      IMPLICIT NONE

      TYPE(rblock_type), INTENT(OUT) :: rb_all
      TYPE(rblock_type), INTENT(IN) :: rb_loc
      INTEGER(i4), INTENT(IN) :: nmodes,nmodes_total
c-----------------------------------------------------------------------
c     create an rblock structure in layer 0 with array lengths for all
c     Fourier components for the fundamental variables.  first copy
c     geometry info and equilibrium fields.
c-----------------------------------------------------------------------
      IF (ilayer==0) THEN
        rb_all%mx=rb_loc%mx
        rb_all%my=rb_loc%my
        rb_all%id=rb_loc%id
        rb_all%degenerate=rb_loc%degenerate
        rb_all%name=rb_loc%name
        CALL copy_lagr_quad_2D(rb_loc%rz,rb_all%rz)
        CALL copy_lagr_quad_2D(rb_loc%be_eq,rb_all%be_eq)
        CALL copy_lagr_quad_2D(rb_loc%ja_eq,rb_all%ja_eq)
        CALL copy_lagr_quad_2D(rb_loc%ve_eq,rb_all%ve_eq)
        CALL copy_lagr_quad_2D(rb_loc%pres_eq,rb_all%pres_eq)
        CALL copy_lagr_quad_2D(rb_loc%prese_eq,rb_all%prese_eq)
        CALL copy_lagr_quad_2D(rb_loc%nd_eq,rb_all%nd_eq)
        CALL copy_lagr_quad_2D(rb_loc%diff_shape,rb_all%diff_shape)
      ENDIF
c-----------------------------------------------------------------------
c     now gather all Fourier components of the dependent variables.
c-----------------------------------------------------------------------
      CALL gather_lagr_quad(rb_loc%be,rb_all%be,nmodes,nmodes_total)
      CALL gather_lagr_quad(rb_loc%ve,rb_all%ve,nmodes,nmodes_total)
      CALL gather_lagr_quad(rb_loc%pres,rb_all%pres,nmodes,nmodes_total)
      CALL gather_lagr_quad(rb_loc%prese,rb_all%prese,nmodes,
     $                      nmodes_total)
      CALL gather_lagr_quad(rb_loc%nd,rb_all%nd,nmodes,nmodes_total)
      CALL gather_lagr_quad(rb_loc%conc,rb_all%conc,nmodes,nmodes_total)
      CALL gather_lagr_quad(rb_loc%tele,rb_all%tele,nmodes,nmodes_total)
      CALL gather_lagr_quad(rb_loc%tion,rb_all%tion,nmodes,nmodes_total)

c     not needed for diffusive corrections.
c     IF (rb_loc%auxb%nqty>0)
c    $  CALL gather_modal_quad(rb_loc%auxb,rb_all%auxb,nmodes,
c    $                        nmodes_total)
c     IF (rb_loc%auxv%nqty>0)
c    $  CALL gather_modal_quad(rb_loc%auxv,rb_all%auxv,nmodes,
c    $                         nmodes_total)

      RETURN
      END SUBROUTINE gather_rblock

c-----------------------------------------------------------------------
c     subprogram 8. copy_bicube:  allocate a bicube structure and copy.
c-----------------------------------------------------------------------
      SUBROUTINE copy_bicube(bc_old,bc_new)
      USE local
      USE bicube
      IMPLICIT NONE

      TYPE(bicube_type), INTENT(IN) :: bc_old
      TYPE(bicube_type), INTENT(OUT) :: bc_new

      CALL bicube_alloc(bc_new,bc_old%mx,bc_old%my,bc_old%nqty)
      bc_new%name=bc_old%name
      bc_new%title=bc_old%title
      bc_new%xs=bc_old%xs
      bc_new%ys=bc_old%ys
      bc_new%fs=bc_old%fs
      bc_new%fsx=bc_old%fsx
      bc_new%fsy=bc_old%fsy
      bc_new%fsxy=bc_old%fsxy

      RETURN
      END SUBROUTINE copy_bicube

c-----------------------------------------------------------------------
c     subprogram 9. copy_lagr_quad_2D:  allocate a 2D lagrange_quad
c     structure and copy.
c-----------------------------------------------------------------------
      SUBROUTINE copy_lagr_quad_2D(laq_old,laq_new)
      USE local
      USE lagr_quad_mod
      IMPLICIT NONE

      TYPE(lagr_quad_2D_type), INTENT(IN) :: laq_old
      TYPE(lagr_quad_2D_type), INTENT(OUT) :: laq_new

      CALL lagr_quad_alloc(laq_new,laq_old%mx,laq_old%my,laq_old%nqty,
     $                     laq_old%n_side+1_i4)
      laq_new%name=laq_old%name
      laq_new%title=laq_old%title
      laq_new%fs=laq_old%fs
      IF (laq_old%n_side>0) THEN
        laq_new%fsh=laq_old%fsh
        laq_new%fsv=laq_old%fsv
      ENDIF
      IF (laq_old%n_int>0) THEN
        laq_new%fsi=laq_old%fsi
      ENDIF

      RETURN
      END SUBROUTINE copy_lagr_quad_2D

c-----------------------------------------------------------------------
c     subprogram 10. copy_modal_quad_2D:  allocate a 2D modal_quad
c     structure and copy.
c-----------------------------------------------------------------------
      SUBROUTINE copy_modal_quad_2D(modq_old,modq_new)
      USE local
      USE modal_disc_mod
      IMPLICIT NONE

      TYPE(modal_quad_2D_type), INTENT(IN) :: modq_old
      TYPE(modal_quad_2D_type), INTENT(OUT) :: modq_new

      CALL modal_disc_alloc(modq_new,modq_old%mx,modq_old%my,
     $                      modq_old%nqty,modq_old%pd,modq_old%pdmin,
     $                      modq_old%pdmax)
      modq_new%name=modq_old%name
      modq_new%title=modq_old%title
      modq_new%fsi=modq_old%fsi

      RETURN
      END SUBROUTINE copy_modal_quad_2D

c-----------------------------------------------------------------------
c     subprogram 11. gather_lagr_quad:  allocate a lagr_quad structure,
c     copy Fourier components from the local processor, then gather
c     components from the other layers.
c-----------------------------------------------------------------------
      SUBROUTINE gather_lagr_quad(laq_old,laq_new,nmodes,nmodes_total)
      USE local
      USE lagr_quad_mod
      USE lagr_disc_mod
      USE mpi_nim
      USE pardata
      IMPLICIT NONE

      TYPE(lagr_quad_type), INTENT(IN) :: laq_old
      TYPE(lagr_quad_type), INTENT(OUT) :: laq_new
      INTEGER(i4), INTENT(IN) :: nmodes,nmodes_total

      INTEGER(i4) :: nq,ns,nsh,nsv,nsi,ierror,il,nl,rem,nm,last
      INTEGER(i4) :: status(mpi_status_size)
      COMPLEX(r8) :: dummy=0
c-----------------------------------------------------------------------
c     first copy lagr_quad boiler-plate, and the Fourier components
c     in layer 0.
c-----------------------------------------------------------------------
      nq=laq_old%nqty
      IF (ilayer==0) THEN
        IF (ALLOCATED(laq_old%fs)) THEN     !  continuous field
          CALL lagr_quad_alloc(laq_new,laq_old%mx,laq_old%my,
     $                         nq,nmodes_total,laq_old%n_side+1_i4)
        ELSE                                 !  discontinuous field
          CALL lagr_disc_alloc(laq_new,laq_old%mx,laq_old%my,
     $                         nq,nmodes_total,
     $                         NINT(SQRT(REAL(laq_old%n_int)))-1_i4)
        ENDIF
        laq_new%name=laq_old%name
        laq_new%title=laq_old%title
        IF (ALLOCATED(laq_new%fs))
     $    laq_new%fs(:,:,:,mode_lo:mode_hi)=laq_old%fs
        IF (ALLOCATED(laq_new%fsh))
     $    laq_new%fsh(:,:,:,:,mode_lo:mode_hi)=laq_old%fsh
        IF (ALLOCATED(laq_new%fsv))
     $    laq_new%fsv(:,:,:,:,mode_lo:mode_hi)=laq_old%fsv
        IF (ALLOCATED(laq_new%fsi))
     $    laq_new%fsi(:,:,:,:,mode_lo:mode_hi)=laq_old%fsi
      ENDIF
c-----------------------------------------------------------------------
c     send Fourier components to layer 0.  the extra communication helps
c     synchronize.
c-----------------------------------------------------------------------
      IF (nmodes/=nmodes_total) THEN
        last=mode_hi
        nl=nprocs/nprocs_layer
        rem=MODULO(nmodes_total,nl)
        DO il=1,nl-1
          nm=nmodes_total/nl
          IF (il<rem) nm=nm+1
          IF (ALLOCATED(laq_old%fs))
     $      ns=nq*nm*(laq_old%mx+1)*(laq_old%my+1)
          IF (ALLOCATED(laq_old%fsh))
     $      nsh=nq*nm*laq_old%mx*(laq_old%my+1)*laq_old%n_side
          IF (ALLOCATED(laq_old%fsv))
     $      nsv=nq*nm*(laq_old%mx+1)*laq_old%my*laq_old%n_side
          IF (ALLOCATED(laq_old%fsi))
     $      nsi=nq*nm*laq_old%mx*laq_old%my*laq_old%n_int
          IF (layer2proc(il)==node) THEN
            CALL mpi_recv(dummy,1,mpi_nim_comp,layer2proc(0),0,
     $           mpi_comm_world,status,ierror)
            IF (ALLOCATED(laq_old%fs))
     $        CALL mpi_send(laq_old%fs,ns,mpi_nim_comp,layer2proc(0),0,
     $           mpi_comm_world,ierror)
            IF (ALLOCATED(laq_old%fsh))
     $        CALL mpi_send(laq_old%fsh,nsh,
     $             mpi_nim_comp,layer2proc(0),0,
     $             mpi_comm_world,ierror)
            IF (ALLOCATED(laq_old%fsv))
     $        CALL mpi_send(laq_old%fsv,nsv,
     $             mpi_nim_comp,layer2proc(0),0,
     $             mpi_comm_world,ierror)
            IF (ALLOCATED(laq_old%fsi))
     $        CALL mpi_send(laq_old%fsi,nsi,
     $             mpi_nim_comp,layer2proc(0),0,
     $             mpi_comm_world,ierror)
          ELSE IF (layer2proc(0)==node) THEN
            CALL mpi_send(dummy,1,mpi_nim_comp,layer2proc(il),0,
     $           mpi_comm_world,ierror)
            IF (ALLOCATED(laq_old%fs))
     $        CALL mpi_recv(laq_new%fs(:,:,:,last+1:last+nm),ns,
     $             mpi_nim_comp,layer2proc(il),0,
     $           mpi_comm_world,status,ierror)
            IF (ALLOCATED(laq_new%fsh))
     $        CALL mpi_recv(laq_new%fsh(:,:,:,:,last+1:last+nm),nsh,
     $             mpi_nim_comp,layer2proc(il),0,
     $             mpi_comm_world,status,ierror)
            IF (ALLOCATED(laq_old%fsv))
     $        CALL mpi_recv(laq_new%fsv(:,:,:,:,last+1:last+nm),nsv,
     $             mpi_nim_comp,layer2proc(il),0,
     $             mpi_comm_world,status,ierror)
            IF (ALLOCATED(laq_old%fsi))
     $        CALL mpi_recv(laq_new%fsi(:,:,:,:,last+1:last+nm),nsi,
     $             mpi_nim_comp,layer2proc(il),0,
     $             mpi_comm_world,status,ierror)
            last=last+nm
          ENDIF
        ENDDO
      ENDIF

      RETURN
      END SUBROUTINE gather_lagr_quad

c-----------------------------------------------------------------------
c     subprogram 12. gather_modal_quad:  allocate a modal_quad
c     structure, copy Fourier components from the local processor,
c     then gather components from the other layers.
c-----------------------------------------------------------------------
      SUBROUTINE gather_modal_quad(modq_old,modq_new,nmodes,
     $                             nmodes_total)
      USE local
      USE modal_disc_mod
      USE mpi_nim
      USE pardata
      IMPLICIT NONE

      TYPE(modal_quad_type), INTENT(IN) :: modq_old
      TYPE(modal_quad_type), INTENT(OUT) :: modq_new
      INTEGER(i4), INTENT(IN) :: nmodes,nmodes_total

      INTEGER(i4) :: nq,pd,pdmin,pdmax,nsi,ierror,il,nl,rem,nm,last
      INTEGER(i4) :: status(mpi_status_size)
      COMPLEX(r8) :: dummy=0
c-----------------------------------------------------------------------
c     first copy lagr_quad boiler-plate, and the Fourier components
c     in layer 0.
c-----------------------------------------------------------------------
      nq=modq_old%nqty
      pd=modq_old%pd
      pdmin=modq_old%pdmin
      pdmax=modq_old%pdmax
      IF (ilayer==0) THEN
        CALL modal_disc_alloc(modq_new,modq_old%mx,modq_old%my,
     $                        nq,nmodes_total,pd,pdmin,pdmax)
        modq_new%name=modq_old%name
        modq_new%title=modq_old%title
        modq_new%fsi(:,:,:,:,mode_lo:mode_hi)=modq_old%fsi
      ENDIF
c-----------------------------------------------------------------------
c     send Fourier components to layer 0.  the extra communication helps
c     synchronize.
c-----------------------------------------------------------------------
      IF (nmodes/=nmodes_total) THEN
        last=mode_hi
        nl=nprocs/nprocs_layer
        rem=MODULO(nmodes_total,nl)
        DO il=1,nl-1
          nm=nmodes_total/nl
          IF (il<rem) nm=nm+1
          nsi=nq*nm*modq_old%mx*modq_old%my*modq_old%n_int
          IF (layer2proc(il)==node) THEN
            CALL mpi_recv(dummy,1,mpi_nim_comp,layer2proc(0),0,
     $           mpi_comm_world,status,ierror)
            CALL mpi_send(modq_old%fsi,nsi,
     $           mpi_nim_comp,layer2proc(0),0,
     $           mpi_comm_world,ierror)
          ELSE IF (layer2proc(0)==node) THEN
            CALL mpi_send(dummy,1,mpi_nim_comp,layer2proc(il),0,
     $           mpi_comm_world,ierror)
            CALL mpi_recv(modq_new%fsi(:,:,:,:,last+1:last+nm),nsi,
     $           mpi_nim_comp,layer2proc(il),0,
     $           mpi_comm_world,status,ierror)
            last=last+nm
          ENDIF
        ENDDO
      ENDIF

      RETURN
      END SUBROUTINE gather_modal_quad

c-----------------------------------------------------------------------
c     subprogram 13. trim_rblock:  reduce an rblock to the Fourier
c     components in the local layer.  this is part of a scatter
c     operation.
c-----------------------------------------------------------------------
      SUBROUTINE trim_rblock(rb,nmodes,nmodes_total)
      USE local
      USE rblock_type_mod
      IMPLICIT NONE

      TYPE(rblock_type), INTENT(INOUT) :: rb
      INTEGER(i4), INTENT(IN) :: nmodes,nmodes_total

      CALL trim_lagr_quad(rb%be,nmodes,nmodes_total)
      CALL trim_lagr_quad(rb%ve,nmodes,nmodes_total)
      CALL trim_lagr_quad(rb%pres,nmodes,nmodes_total)
      CALL trim_lagr_quad(rb%prese,nmodes,nmodes_total)
      CALL trim_lagr_quad(rb%nd,nmodes,nmodes_total)
      CALL trim_lagr_quad(rb%conc,nmodes,nmodes_total)
      CALL trim_lagr_quad(rb%tele,nmodes,nmodes_total)
      CALL trim_lagr_quad(rb%tion,nmodes,nmodes_total)

c     not needed for diffusive corrections.
c     IF (rb%auxb%nqty>0)
c    $  CALL trim_modal_quad(rb%auxb,nmodes,nmodes_total)
c     IF (rb%auxv%nqty>0)
c    $  CALL trim_modal_quad(rb%auxv,nmodes,nmodes_total)

      RETURN
      END SUBROUTINE trim_rblock

c-----------------------------------------------------------------------
c     subprogram 14. trim_lagr_quad:  reduce a lagr_quad structure to
c     the Fourier components on the local layer.
c-----------------------------------------------------------------------
      SUBROUTINE trim_lagr_quad(laq,nmodes,nmodes_total)
      USE local
      USE lagr_quad_mod
      USE mpi_nim
      USE pardata
      IMPLICIT NONE

      TYPE(lagr_quad_type), INTENT(INOUT) :: laq
      INTEGER(i4), INTENT(IN) :: nmodes,nmodes_total

      INTEGER(i4) :: nq
      COMPLEX(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: tmp
      COMPLEX(r8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: tmp2

      nq=laq%nqty
      laq%nfour=nmodes
      DEALLOCATE(laq%f,laq%fx,laq%fy)
      ALLOCATE(laq%f(nq,nmodes))
      ALLOCATE(laq%fx(nq,nmodes))
      ALLOCATE(laq%fy(nq,nmodes))

      IF (ALLOCATED(laq%fs)) THEN
        ALLOCATE(tmp(nq,0:laq%mx,0:laq%my,nmodes))
        tmp=laq%fs(:,:,:,mode_lo:mode_hi)
        DEALLOCATE(laq%fs)
        ALLOCATE(laq%fs(nq,0:laq%mx,0:laq%my,nmodes))
        laq%fs=tmp
        DEALLOCATE(tmp)
      ENDIF

      IF (ALLOCATED(laq%fsh)) THEN
        ALLOCATE(tmp2(nq,laq%n_side,laq%mx,0:laq%my,nmodes))
        tmp2=laq%fsh(:,:,:,:,mode_lo:mode_hi)
        DEALLOCATE(laq%fsh)
        ALLOCATE(laq%fsh(nq,laq%n_side,laq%mx,0:laq%my,nmodes))
        laq%fsh=tmp2
        DEALLOCATE(tmp2)
      ENDIF

      IF (ALLOCATED(laq%fsv)) THEN
        ALLOCATE(tmp2(nq,laq%n_side,0:laq%mx,laq%my,nmodes))
        tmp2=laq%fsv(:,:,:,:,mode_lo:mode_hi)
        DEALLOCATE(laq%fsv)
        ALLOCATE(laq%fsv(nq,laq%n_side,0:laq%mx,laq%my,nmodes))
        laq%fsv=tmp2
        DEALLOCATE(tmp2)
      ENDIF

      IF (ALLOCATED(laq%fsi)) THEN
        ALLOCATE(tmp2(nq,laq%n_int,laq%mx,laq%my,nmodes))
        tmp2=laq%fsi(:,:,:,:,mode_lo:mode_hi)
        DEALLOCATE(laq%fsi)
        ALLOCATE(laq%fsi(nq,laq%n_int,laq%mx,laq%my,nmodes))
        laq%fsi=tmp2
        DEALLOCATE(tmp2)
      ENDIF

      RETURN
      END SUBROUTINE trim_lagr_quad

c-----------------------------------------------------------------------
c     subprogram 15. trim_modal_quad:  reduce a modal_quad structure to
c     the Fourier components on the local layer.
c-----------------------------------------------------------------------
      SUBROUTINE trim_modal_quad(modq,nmodes,nmodes_total)
      USE local
      USE modal_disc_mod
      USE mpi_nim
      USE pardata
      IMPLICIT NONE

      TYPE(modal_quad_type), INTENT(INOUT) :: modq
      INTEGER(i4), INTENT(IN) :: nmodes,nmodes_total

      INTEGER(i4) :: nq
      COMPLEX(r8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: tmp2

      nq=modq%nqty
      modq%nfour=nmodes
      DEALLOCATE(modq%f,modq%fx,modq%fy)
      ALLOCATE(modq%f(nq,nmodes))
      ALLOCATE(modq%fx(nq,nmodes))
      ALLOCATE(modq%fy(nq,nmodes))

      ALLOCATE(tmp2(nq,modq%n_int,modq%mx,modq%my,nmodes))
      tmp2=modq%fsi(:,:,:,:,mode_lo:mode_hi)
      DEALLOCATE(modq%fsi)
      ALLOCATE(modq%fsi(nq,modq%n_int,modq%mx,modq%my,nmodes))
      modq%fsi=tmp2
      DEALLOCATE(tmp2)

      RETURN
      END SUBROUTINE trim_modal_quad

c-----------------------------------------------------------------------
c     subprogram 16. gather_tblock:  create a tblock with
c     complete information from all Fourier layers.
c-----------------------------------------------------------------------
      SUBROUTINE gather_tblock(tb_all,tb_loc,nmodes,nmodes_total)
      USE local
      USE tblock_type_mod
      USE pardata
      IMPLICIT NONE

      TYPE(tblock_type), INTENT(OUT) :: tb_all
      TYPE(tblock_type), INTENT(IN) :: tb_loc
      INTEGER(i4), INTENT(IN) :: nmodes,nmodes_total

      INTEGER(i4) :: iv,nnbr
c-----------------------------------------------------------------------
c     create a tblock structure on layer 0 with array lengths for all
c     Fourier components for the fundamental variables.  first copy
c     geometry info and equilibrium fields.
c-----------------------------------------------------------------------
      IF (ilayer==0) THEN
        tb_all%name=tb_loc%name
        tb_all%id=tb_loc%id
        tb_all%mvert=tb_loc%mvert
        tb_all%mcell=tb_loc%mcell
        tb_all%tgeom%mvert=tb_loc%mvert
        tb_all%tgeom%mcell=tb_loc%mcell
        CALL tri_linear_geom_alloc(tb_all%tgeom,tb_all%mvert,
     $    tb_all%mcell)
        tb_all%tgeom%xs=tb_loc%tgeom%xs
        tb_all%tgeom%ys=tb_loc%tgeom%ys
        tb_all%tgeom%vertex=tb_loc%tgeom%vertex
        DO iv=0,tb_all%tgeom%mvert
           nnbr=SIZE(tb_loc%tgeom%neighbor(iv)%vertex)-1
           ALLOCATE(tb_all%tgeom%neighbor(iv)%vertex(0:nnbr))
           tb_all%tgeom%neighbor(iv)%vertex=
     $        tb_loc%tgeom%neighbor(iv)%vertex
        ENDDO
        CALL copy_tri_linear_2D(tb_loc%be_eq,tb_all%be_eq)
        CALL copy_tri_linear_2D(tb_loc%ja_eq,tb_all%ja_eq)
        CALL copy_tri_linear_2D(tb_loc%ve_eq,tb_all%ve_eq)
        CALL copy_tri_linear_2D(tb_loc%pres_eq,tb_all%pres_eq)
        CALL copy_tri_linear_2D(tb_loc%prese_eq,tb_all%prese_eq)
        CALL copy_tri_linear_2D(tb_loc%nd_eq,tb_all%nd_eq)
        CALL copy_tri_linear_2D(tb_loc%diff_shape,tb_all%diff_shape)
      ENDIF
c-----------------------------------------------------------------------
c     now gather Fourier components.
c-----------------------------------------------------------------------
      CALL gather_tri_linear(tb_loc%be,tb_all%be,
     $                       'gather',nmodes,nmodes_total)
      CALL gather_tri_linear(tb_loc%ve,tb_all%ve,
     $                       'gather',nmodes,nmodes_total)
      CALL gather_tri_linear(tb_loc%pres,tb_all%pres,
     $                       'gather',nmodes,nmodes_total)
      CALL gather_tri_linear(tb_loc%prese,tb_all%prese,
     $                       'gather',nmodes,nmodes_total)
      CALL gather_tri_linear(tb_loc%nd,tb_all%nd,
     $                       'gather',nmodes,nmodes_total)
      CALL gather_tri_linear(tb_loc%conc,tb_all%conc,
     $                       'gather',nmodes,nmodes_total)
      CALL gather_tri_linear(tb_loc%tele,tb_all%tele,
     $                       'gather',nmodes,nmodes_total)
      CALL gather_tri_linear(tb_loc%tion,tb_all%tion,
     $                       'gather',nmodes,nmodes_total)

      RETURN
      END SUBROUTINE gather_tblock

c-----------------------------------------------------------------------
c     subprogram 17. copy_tri_linear_2D:  copy a tri_linear_2D
c     structure.
c-----------------------------------------------------------------------
      SUBROUTINE copy_tri_linear_2D(tl_old,tl_new)
      USE local
      USE tri_linear
      IMPLICIT NONE

      TYPE(tri_linear_2D_type), INTENT(IN) :: tl_old
      TYPE(tri_linear_2D_type), INTENT(OUT) :: tl_new

      CALL tri_linear_2D_alloc(tl_new,tl_old%mvert,tl_old%nqty)
      tl_new%name=tl_old%name
      tl_new%title=tl_old%title
      tl_new%fs=tl_old%fs

      RETURN
      END SUBROUTINE copy_tri_linear_2D

c-----------------------------------------------------------------------
c     subprogram 18. gather_tri_linear:  allocate a tri_linear
c     structure, copy Fourier components from the local processor, then
c     gather components from the other layers.
c-----------------------------------------------------------------------
      SUBROUTINE gather_tri_linear(tl_old,tl_new,op,nmodes,nmodes_total)
      USE local
      USE tri_linear
      USE mpi_nim
      USE pardata
      IMPLICIT NONE

      TYPE(tri_linear_type), INTENT(IN) :: tl_old
      TYPE(tri_linear_type), INTENT(OUT) :: tl_new
      CHARACTER(*), INTENT(IN) :: op
      INTEGER(i4), INTENT(IN) :: nmodes,nmodes_total

      INTEGER(i4) :: nq,nqt,ns,ierror,il,nl,rem,nm,last
      INTEGER(i4) :: status(mpi_status_size)
      COMPLEX(r8) :: dummy=0
c-----------------------------------------------------------------------
c     first copy tri_linear boiler-plate, and the Fourier components
c     in layer 0.
c-----------------------------------------------------------------------
      IF (ilayer==0) THEN
        IF (op=='gather') THEN
          nq=SIZE(tl_old%fs,1)
          CALL tri_linear_alloc(tl_new,tl_old%mvert,nq,nmodes_total)
          tl_new%name=tl_old%name
          tl_new%title=tl_old%title(1)
          tl_new%fs(:,:,:,mode_lo:mode_hi)=tl_old%fs
        ELSE
          nq=tl_old%nqty
          CALL tri_linear_alloc(tl_new,tl_old%mvert,nq,tl_old%nfour)
          tl_new%name=tl_old%name
          tl_new%title=tl_old%title
          tl_new%fs=tl_old%fs
        ENDIF
      ENDIF
      IF (op/='gather') RETURN
c-----------------------------------------------------------------------
c     send Fourier components to layer 0.  the extra communication helps
c     synchronize.
c-----------------------------------------------------------------------
      nq=tl_old%nqty
      IF (nmodes/=nmodes_total) THEN
        last=mode_hi
        nl=nprocs/nprocs_layer
        rem=MODULO(nmodes_total,nl)
        DO il=1,nl-1
          nm=nmodes_total/nl
          IF (il<rem) nm=nm+1
          ns=nq*nm*(tl_old%mvert+1)
          IF (layer2proc(il)==node) THEN
            CALL mpi_recv(dummy,1,mpi_nim_comp,layer2proc(0),0,
     $           mpi_comm_world,status,ierror)
            CALL mpi_send(tl_old%fs,ns,mpi_nim_comp,layer2proc(0),0,
     $           mpi_comm_world,ierror)
          ELSE IF (layer2proc(0)==node) THEN
            CALL mpi_send(dummy,1,mpi_nim_comp,layer2proc(il),0,
     $           mpi_comm_world,ierror)
            CALL mpi_recv(tl_new%fs(:,:,:,last+1:last+nm),ns,
     $           mpi_nim_comp,layer2proc(il),0,
     $           mpi_comm_world,status,ierror)
            last=last+nm
          ENDIF
        ENDDO
      ENDIF

      RETURN
      END SUBROUTINE gather_tri_linear

c-----------------------------------------------------------------------
c     subprogram 19. trim_tblock:  reduce an tblock to the Fourier
c     components in the local layer.  this is part of a scatter
c     operation.
c-----------------------------------------------------------------------
      SUBROUTINE trim_tblock(tb,nmodes,nmodes_total)
      USE local
      USE tblock_type_mod
      IMPLICIT NONE

      TYPE(tblock_type), INTENT(INOUT) :: tb
      INTEGER(i4), INTENT(IN) :: nmodes,nmodes_total

      CALL trim_tri_linear(tb%be,nmodes,nmodes_total)
      CALL trim_tri_linear(tb%ve,nmodes,nmodes_total)
      CALL trim_tri_linear(tb%pres,nmodes,nmodes_total)
      CALL trim_tri_linear(tb%prese,nmodes,nmodes_total)
      CALL trim_tri_linear(tb%nd,nmodes,nmodes_total)
      CALL trim_tri_linear(tb%conc,nmodes,nmodes_total)
      CALL trim_tri_linear(tb%tele,nmodes,nmodes_total)
      CALL trim_tri_linear(tb%tion,nmodes,nmodes_total)

      RETURN
      END SUBROUTINE trim_tblock

c-----------------------------------------------------------------------
c     subprogram 20. trim_tri_linear:  reduce a tri_linear structure to
c     the Fourier components on the local layer.
c-----------------------------------------------------------------------
      SUBROUTINE trim_tri_linear(tl,nmodes,nmodes_total)
      USE local
      USE tri_linear
      USE mpi_nim
      USE pardata
      IMPLICIT NONE

      TYPE(tri_linear_type), INTENT(INOUT) :: tl
      INTEGER(i4), INTENT(IN) :: nmodes,nmodes_total

      INTEGER(i4) :: nq
      COMPLEX(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: tmp

      nq=tl%nqty
      tl%nfour=nmodes
      DEALLOCATE(tl%f,tl%fx,tl%fy)
      ALLOCATE(tl%f(nq,nmodes))
      ALLOCATE(tl%fx(nq,nmodes))
      ALLOCATE(tl%fy(nq,nmodes))

      ALLOCATE(tmp(nq,0:tl%mvert,0:0,nmodes))
      tmp=tl%fs(:,:,:,mode_lo:mode_hi)
      DEALLOCATE(tl%fs)
      ALLOCATE(tl%fs(nq,0:tl%mvert,0:0,nmodes))
      tl%fs=tmp
      DEALLOCATE(tmp)

      RETURN
      END SUBROUTINE trim_tri_linear

