c-----------------------------------------------------------------------
c     file pardata.f
c
c     data structures for performing block decomposition and 
c       seam communication in parallel
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c     data structures.
c-----------------------------------------------------------------------

      MODULE pardata
      USE local
      IMPLICIT NONE

c nprocs = total # of processors (assigned to blocks and layers)
c node = processor id of me (0 to nprocs-1)

      integer(i4) :: nprocs
      integer(i4) :: node                   

c layer2proc(0:nlayers-1) = proc ids of all layers for the
c   local grid blocks.

      integer(i4), dimension(:), allocatable :: layer2proc

c nprocs_layer = # of procs assigned to each layer
c node_layer = I am this proc within the layer (0 to nprocs_layer-1)
c ilayer = which layer this proc belong to (0 to nlayers-1)
c mode_lo = my 1st mode corresponds to this global mode (1 to nmodes_total)
c mode_hi = my last mode corresponds to this global mode (1 to nmodes_total)
c mode2layer = maps each global mode (1 to nmodes_total) to a unique layer

      integer(i4) :: nprocs_layer
      integer(i4) :: node_layer
      integer(i4) :: ilayer
      integer(i4) :: mode_lo,mode_hi
      integer(i4), dimension(:), allocatable :: mode2layer

c MPI communicators
c comm_layer = within a layer, across blocks
c              all the procs owning blocks in the same layer
c comm_mode =  withih a block, across modes (layers)
c              all the procs owning modes for the same block(s)

      integer(i4) :: comm_layer,comm_mode

c block2proc(1:nbl_total)
c   block2proc(i) = proc id of owner of global block i (within my layer)
c global2local(1:nbl_total)
c   global2local(i) = local index (1:nbl) of global block i 
c                     on proc who owns it (within my layer)
c loc2glob(1:nbl)
c   loc2glob(i) = global index (1:nbl_total) of local block i
c block_sizes(3,1:nbl_total)
c   block_sizes(1,i) = nvert (# of seam pts) in global block i (r or tblocks)
c   block_sizes(2,i) = mx (# of x vertices) in global rblock i (0 for tblocks)
c   block_sizes(3,i) = my (# of y vertices) in global rblock i (0 for tblocks)

      integer(i4), dimension(:), allocatable :: block2proc
      integer(i4), dimension(:), allocatable :: global2local
      integer(i4), dimension(:), allocatable :: loc2glob
      integer(i4), dimension(:,:), allocatable :: block_sizes

c data stucture for message sending of SEAM DATA
c nsend = # of messages I will send
c allocate send(1:nsend) of send_type, 
c   one for each proc a message will be sent to
c proc = processor id to send message to
c count = number of datums to send, a "datum" = nqty values
c block(count) = which of my local blocks (1:nbl) the datum comes from
c vertex(count) = which seam vertex (1:nvert) the datum comes from
c image(count) = which vertex image (1:nimage) the datum corresponds to,
c data(nqty*count) = all values packed into one message,
c                    is a 1-d vector so that data can be packed contiguously
c                    regardless of what nqty is on a particular seaming call
c cdata(nqty*count) = a complex version of data

      type :: send_type
        integer(i4) :: proc
        integer(i4) :: count
        integer(i4), dimension(:), pointer :: block
        integer(i4), dimension(:), pointer :: vertex
        integer(i4), dimension(:), pointer :: image
        real(r8), dimension(:), pointer :: data
        complex(r8), dimension(:), pointer :: cdata
      end type send_type

      integer(i4) :: nsend
      type(send_type), dimension(:), pointer :: send

c data stucture for message receiving of SEAM DATA
c nrecv = # of messages I will receive
c allocate recv(1:nrecv) of recv_type, 
c   one for each proc a message will come from
c proc = processor id of who will send me the messaeg
c count = number of datums I will receive, a "datum" = nqty values
c block(count) = which of my local blocks (1:nbl) the datum gets summed to
c vertex(count) = which seam vertex in my block the datum corresponds to
c order(count) = unique ordering of this image in all (except degenerate
c		 points) representations so sums produce identical round-off
c data(nqty*count) = all values packed into one message,
c                    received data will be packed contiguously
c                    regardless of what nqty is on a particular seaming call
c cdata(nqty*count) = a complex version of data
c recv_request = array of requests for posting asynchronous receives
      
      type :: recv_type
        integer(i4) :: proc
        integer(i4) :: count
        integer(i4), dimension(:), pointer :: block
        integer(i4), dimension(:), pointer :: vertex
        integer(i4), dimension(:), pointer :: order
        real(r8), dimension(:), pointer :: data
        complex(r8), dimension(:), pointer :: cdata
      end type recv_type

      integer(i4) :: nrecv
      type(recv_type), dimension(:), pointer :: recv
      integer(i4), dimension(:), allocatable :: recv_request
      
c data stucture for performing SEAMING of values between blocks I own
c nself = # of seam points whose images are also owned by me
c not needed when nprocs = 1, since can just loop over all seams
c allocate self(1:nself) of self type, 
c   one for each point, note that a pair of 
c   self-referencing pts will be stored twice
c block_out = which of my local blocks (1:nbl) the datum gets summed to
c vertex_out = which seam vertex (1:nvert) the datum gets summed to
c block_in = which of my local blocks (1:nbl) the datum comes from
c vertex_in = which seam vertex (1:nvert) the datum comes from
c order = unique ordering of this image in all (except degenerate
c	  points) representations so sums produce identical round-off

      type :: self_type
        integer(i4) block_out,vertex_out
        integer(i4) block_in,vertex_in
        integer(i4) order
      end type self_type

      integer(i4) :: nself
      type(self_type), dimension(:), pointer :: self

c data stucture for message sending of SEGMENT DATA
c nsendseg = # of messages I will send
c allocate sendseg(1:nsendseg) of sendseg_type, 
c   one for each proc a message will be sent to
c proc = processor id to send message to
c count = number of datums to send, a "datum" = nqty x nqty values
c block(count) = which of my local blocks (1:nbl) the datum comes from
c segment(count) = which segment (1:nvert) the datum comes from
c data(nqty*count) = all values packed into one message,
c                    is a 1-d vector so that data can be packed contiguously
c cdata(nqty*count) = a complex version of data
c mat_data(nmat*nqty**2*count) = holds nmat off-diagonal matrix
c                                elements without symmetry assumptions
c cmat_data(nmat*nqty**2*count) = complex version of mat_data.

      type :: sendseg_type
        integer(i4) :: proc
        integer(i4) :: count
        integer(i4), dimension(:), pointer :: block
        integer(i4), dimension(:), pointer :: segment
        real(r8), dimension(:), pointer :: data
        complex(r8), dimension(:), pointer :: cdata
        real(r8), dimension(:), pointer :: mat_data
        complex(r8), dimension(:), pointer :: cmat_data
      end type sendseg_type

      integer(i4) :: nsendseg
      type(sendseg_type), dimension(:), pointer :: sendseg

c data stucture for message receiving of SEGMENT DATA
c nrecvseg = # of messages I will receive
c allocate recvseg(1:nrecvseg) of recvseg_type, 
c   one for each proc a message will come from
c proc = processor id of who will send me the messaeg
c count = number of datums I will receive, a "datum" = nqty x nqty values
c block(count) = which of my local blocks (1:nbl) the datum gets summed to
c segment(count) = which segment in my block the datum corresponds to
c data(nqty*count) = all values packed into one message,
c                    received data will be packed contiguously
c cdata(nqty*count) = a complex version of data
c mat_data(nmat*nqty**2*count) = holds off-diagonal matrix elements
c                                without symmetry assumptions
c cmat_data(nmat*nqty**2*count) = complex version of mat_data.
c recvseg_request = array of requests for posting asynchronous receives
      
      type :: recvseg_type
        integer(i4) :: proc
        integer(i4) :: count
        integer(i4), dimension(:), pointer :: block
        integer(i4), dimension(:), pointer :: segment
        real(r8), dimension(:), pointer :: data
        complex(r8), dimension(:), pointer :: cdata
        real(r8), dimension(:), pointer :: mat_data
        complex(r8), dimension(:), pointer :: cmat_data
      end type recvseg_type

      integer(i4) :: nrecvseg
      type(recvseg_type), dimension(:), pointer :: recvseg
      integer(i4), dimension(:), allocatable :: recvseg_request
      
c data stucture for seaming SEGMENTS of values between blocks I own
c nselfseg = # of segments whose images are also owned by me
c not needed when nprocs = 1, since can just loop over all segments
c allocate selfseg(1:nselfseg) of selfseg type, 
c   one for each segment, note that a pair of 
c   self-referencing segments will be stored twice
c block_out = which of my local blocks (1:nbl) the datum gets summed to
c segment_out = which segment (1:nvert) the datum gets summed to
c block_in = which of my local blocks (1:nbl) the datum comes from
c segment_in = which segment (1:nvert) the datum comes from

      type :: selfseg_type
        integer(i4) block_out,segment_out
        integer(i4) block_in,segment_in
      end type selfseg_type

      integer(i4) :: nselfseg
      type(selfseg_type), dimension(:), pointer :: selfseg

c the line structure contains data required to initialize communication
c from blocks to global lines for preconditioning:

c data stored by line:
c nlinks = # of blocks contributing to this line.
c perpst & perpen = start and end indices in direction perpendicular to
c   the line.
c parast & paraen = start and end indices in the parallel direction for
c   computations with the line as a whole.  parast=1 for periodic lines
c   and paraen=0 for nonperiodic lines.
c ncomm = # of off-processor communications in the line.
c glb_bl_recv = global block number for the block with the same index
c   as this line.
c bl2line_parast & bl2line_paraen = parallel start and end indices
c   within a block that are sent to a line.
c line2bl_parast & line2bl_paraen = parallel start and end indices
c   within a block that are received from a line.
c par_dir = normal rblock index that represents the parallel direction

c data stored by link or segment of a line.
c node = node from which this line receives this segment.
c glb_bl = global block index holding this segment.
c loc_bl = local block index holding this segment.
c mpara & mperp = cell dimensions of the block containing this segment
c   in the directions parallel and perpendicular to this line.
c bl2line_segst & bl2line_segen = line parallel indices filled by this
c   segment.
c line2bl_segst & line2bl_segen = line parallel indices sent to this
c   segment.
c bl2line_count & line2bl_count = # of block vertices collected
c   from this segment and # shipped back, respectively.
c recv_req & send_req = request arrays used for nonblocking mpi calls
c   (both dimensioned ncomm for this line).
c req_index = converts from segment index to request array index.
c bl_perpst & bl_perpen = start and end perpendicular indices for
c   a segment of the block with this line's index.

      TYPE :: line_type
        INTEGER(i4) :: nlinks,perpst,perpen,parast,paraen,ncomm,
     $                 glb_bl_recv,bl2line_parast,bl2line_paraen,
     $                 line2bl_parast,line2bl_paraen,par_dir
        INTEGER(i4), DIMENSION(:), POINTER :: node,glb_bl,loc_bl,
     $               mpara,mperp,bl2line_segst,bl2line_segen,
     $               bl2line_count,line2bl_segst,line2bl_segen,
     $               line2bl_count,recv_req,send_req,req_index,
     $               bl_perpst,bl_perpen
        LOGICAL :: periodic
      END TYPE line_type

      TYPE(line_type), DIMENSION(:), ALLOCATABLE :: linex,liney

c slugrid_handle = used to store the location of SuperLU_DIST's
c   processor grid, which is allocated in C routines.
c slu_nrowp & slu_ncolp = row and column dimensions of SLUD's
c   processor grid.

      INTEGER(i8) :: slugrid_handle
      INTEGER(i4) :: slu_nrowp,slu_ncolp

c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------

      END MODULE pardata
