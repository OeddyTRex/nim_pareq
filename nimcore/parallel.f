c-----------------------------------------------------------------------
c     file parallel.f
c
c     routines for performing block decomposition and 
c       seam/segment communication in parallel
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     1.  parallel_block_init.
c     2.  parallel_block_init2.
c     3.  parallel_seam_init.
c     4.  parallel_seam_comm.
c     5.  parallel_seam_comm_comp.
c     6.  parallel_seg_init.
c     7.  parallel_seg_comm.
c     8.  parallel_seg_comm_comp.
c     9.  parallel_matseg_comm.
c     10. parallel_matseg_comm_comp.
c     11. parallel_line_init.
c     12. parallel_line_comm_bl2line.
c     13. parallel_line_comm_line2bl.
c     14. parallel_line_comm_bl2line_comp.
c     15. parallel_line_comm_line2bl_comp.
c     16. qp0_bcast
c     17. qpc_bcast
c     18. parallel_block_dealloc.
c     19. parallel_seam_dealloc.
c     20. parallel_seg_dealloc.
c-----------------------------------------------------------------------
c     subprogram 1. parallel_block_init
c     assign blocks to processors
c     setup block2proc and global2local vectors
c-----------------------------------------------------------------------

      SUBROUTINE parallel_block_init(nmodes,nmodes_total,nlayers,
     $                               decompflag,rname)
      USE fields
      USE mpi_nim
      USE pardata
      IMPLICIT NONE

      INTEGER(i4), INTENT(OUT) :: nmodes
      INTEGER(i4), INTENT(IN) :: nmodes_total,nlayers,decompflag
      CHARACTER(*), INTENT(IN) :: rname

      INTEGER(i4) :: ib,ibl,remainder,blocklo,blockhi,ierror,itmp,jtmp,
     $               i,k,nvert,ioerr,block2proc_l
      integer(i4), dimension(:), allocatable :: workvec
      CHARACTER(120), PARAMETER :: fmt1="(/i3,' processor layers have ',
     $  i3,' Fourier comps, and ',i3,' have ',i3)",
     $  fmt2="(/i3,' layers have ',i3,' Fourier comps.')",
     $  fmt3="(' Number of procs (',i4,')',' > nlayers (',i4,')',
     $         ' * number of grid blocks (',i4,').')"
      CHARACTER(120) :: msg
      CHARACTER(32) :: partfile


c check input before dividing by nlayers

      IF (nlayers <= 0)
     $  CALL nim_stop('Nlayers must be > 0.')

c check that nlayers divides into nprocs evenly

      IF (mod(nprocs,nlayers) /= 0)
     $  CALL nim_stop('Number of procs is not divisible by nlayers.')

c check for too many processors

      IF (nprocs>nlayers*nbl_total)THEN
         WRITE(msg,fmt=fmt3)nprocs,nlayers,nbl_total
         CALL nim_stop(msg)
      ENDIF

c check for too many layers

      IF (nlayers>nmodes_total) THEN
         CALL nim_stop
     $     ('nlayers must be <= the number of Fourier components.')
      ENDIF


c assign each layer to a group of procs
c 1st few procs get layer 1, next few get layer 2, etc

      nprocs_layer = nprocs/nlayers
      node_layer = mod(node,nprocs_layer)
      ilayer = node / nprocs_layer

c compute nmodes = # of Fourier modes per processor

      nmodes = nmodes_total / nlayers
      remainder = mod(nmodes_total,nlayers)

c report the distribution of Fourier components.

      IF (node == 0 .AND. nprocs > 1) THEN
        IF (.NOT.out_opened) THEN
          OPEN(UNIT=out_unit,FILE='nimrod.out',STATUS='UNKNOWN',
     $         POSITION='APPEND')
          out_opened=.true.
        ENDIF
        IF (remainder/=0) THEN
          WRITE(out_unit,fmt1) remainder,nmodes+1,nlayers-remainder,
     $      nmodes
          WRITE(nim_wr,fmt1) remainder,nmodes+1,nlayers-remainder,nmodes
        ELSE
          WRITE(out_unit,fmt2) nlayers,nmodes
          WRITE(nim_wr,fmt2) nlayers,nmodes
        ENDIF
      ENDIF

c 1st proc (within a block) gets 1st few modes, next proc gets next few, etc
c all procs in layer < remainder get one extra mode

      if (ilayer < remainder) nmodes = nmodes + 1

      mode_lo = ilayer * nmodes + 1
      if (ilayer >= remainder) mode_lo = remainder*(nmodes+1) +
     $     (ilayer - remainder)*nmodes + 1
      mode_hi = mode_lo + nmodes - 1

c group procs into 2 MPI communicators
c comm_layer = within a layer, across blocks
c comm_mode = within a block, across modes (layers)

      call mpi_comm_split(mpi_comm_world,ilayer,0,
     $     comm_layer,ierror)
      call mpi_comm_split(mpi_comm_world,node_layer,0,
     $     comm_mode,ierror)

c create the mode to layer map

      allocate(mode2layer(nmodes_total))
      allocate(workvec(nmodes_total))

      mode2layer=0
      do ib=mode_lo,mode_hi
        mode2layer(ib)=ilayer
      enddo

      CALL mpi_allreduce(mode2layer,workvec,nmodes_total,
     $     mpi_nim_int,mpi_sum,comm_mode,ierror)
      mode2layer=workvec

      deallocate(workvec)

c allocate space for global block vectors

      allocate(block2proc(nbl_total))
      allocate(global2local(nbl_total))
      block2proc = 0
      global2local = 0

c compute nbl = # of rblocks+tblocks on this processor
c for decompflag=0 or 1, all procs < remainder get one
c extra block

      if (decompflag<2) then
        nbl = nbl_total / nprocs_layer
        remainder = mod(nbl_total,nprocs_layer)
        if (node_layer < remainder) nbl = nbl + 1
      else
        nbl=0
      endif

c assign rblocks & tblocks to procs (within a layer),
c all layers do these computations at the same time
c for decompflag = 0, assign in clumps
c for decompflag = 1, assign in strided fashion
c also compute nrbl = # of rblocks on this processor
c mark the global vectors according to block assignment

      ibl = 0
      nrbl = 0
      if (decompflag == 0 .AND. nprocs>1) then
        blocklo = node_layer * nbl + 1
        if (node_layer >= remainder) blocklo = remainder*(nbl+1) +
     $       (node_layer - remainder)*nbl + 1
        blockhi = blocklo + nbl - 1
        do ib = 1,nbl_total
          if (ib >= blocklo .and. ib <= blockhi) then
            ibl = ibl + 1
            block2proc(ib) = node
            global2local(ib) = ibl
            if (ib <= nrbl_total) nrbl = nrbl + 1
          endif
        enddo
      else if (decompflag == 1 .OR. nprocs==1) then
        do ib = 1,nbl_total
          if (mod(ib-1_i4,nprocs_layer) == node_layer) then
            ibl = ibl + 1
            block2proc(ib) = node
            global2local(ib) = ibl
            if (ib <= nrbl_total) nrbl = nrbl + 1
          endif
        enddo
      else if (decompflag==2) then
        partfile=trim('block_'//trim(rname)//'.sup.part')
        open(unit=grd_unit,file=trim(partfile),form='formatted',
     $       status='old')
        do ib = 1,nbl_total
          read(grd_unit,*) block2proc_l
          if(block2proc_l==node_layer) then
            ibl=ibl+1
            nbl=nbl+1
            block2proc(ib)=node
            global2local(ib)=ibl
            if (ib <= nrbl_total) nrbl = nrbl+1
          endif
        enddo
        close(unit=grd_unit)
      endif

c error checking within and across layers

      ib = 0
      IF (ibl /= nbl) THEN
        ib = 1
        write (nim_wr,*) 'Ibl /= Nbl on proc',node,ibl,nbl
      ENDIF

      CALL mpi_allreduce(ib,itmp,1,
     $     mpi_nim_int,mpi_max,mpi_comm_world,ierror)
      IF (itmp == 1) CALL nim_stop('Decomposition incorrect.')

      CALL mpi_allreduce(nbl,itmp,1,mpi_nim_int,mpi_sum,
     $     comm_layer,ierror)
      jtmp = 0
      IF (itmp /= nbl_total) jtmp = 1
      CALL mpi_allreduce(jtmp,itmp,1,mpi_nim_int,mpi_max,
     $     mpi_comm_world,ierror)
      if (itmp /= 0_i4) CALL nim_stop('Block total incorrect.')

      CALL mpi_allreduce(nrbl,itmp,1,mpi_nim_int,mpi_sum,
     $     comm_layer,ierror)
      jtmp = 0
      IF (itmp /= nrbl_total) jtmp = 1
      CALL mpi_allreduce(jtmp,itmp,1,mpi_nim_int,mpi_max,
     $     mpi_comm_world,ierror)
      IF (itmp /= 0_i4) CALL nim_stop('Rblock total incorrect.')

c merge 2 global block vectors within each layer

      allocate(workvec(nbl_total))

      CALL mpi_allreduce(block2proc,workvec,nbl_total,
     $     mpi_nim_int,mpi_sum,comm_layer,ierror)
      block2proc = workvec

      CALL mpi_allreduce(global2local,workvec,nbl_total,
     $     mpi_nim_int,mpi_sum,comm_layer,ierror)
      global2local = workvec

      deallocate(workvec)

c create the local to global index conversion for convenience.

      allocate(loc2glob(nbl))

      DO ib=1,nbl_total
        IF (block2proc(ib)==node) loc2glob(global2local(ib))=ib
      ENDDO

c create a layer vector for the local blocks.

      allocate(layer2proc(0:nlayers-1),workvec(nlayers))

      layer2proc = 0
      layer2proc(ilayer) = node

      CALL mpi_allreduce(layer2proc,workvec,nlayers,
     $     mpi_nim_int,mpi_sum,comm_mode,ierror)
      layer2proc = workvec

      deallocate(workvec)

      RETURN
      END SUBROUTINE parallel_block_init


c-----------------------------------------------------------------------
c     subprogram 2. parallel_block_init2
c     setup block_sizes global data structure
c-----------------------------------------------------------------------

      SUBROUTINE parallel_block_init2
      USE fields
      USE seam_storage_mod
      USE mpi_nim
      USE pardata
      IMPLICIT NONE
      INTEGER(i4) :: ibl,ierror
      INTEGER(i4), dimension(:,:), allocatable :: workvec

c-----------------------------------------------------------------------
c     setup block_sizes data structure so all procs know size of
c     all r and tblocks (nvert) and bounds of all rblocks (mx,my)
c-----------------------------------------------------------------------

      ALLOCATE(block_sizes(3,nbl_total))
      ALLOCATE(workvec(3,nbl_total))

      block_sizes = 0

      DO ibl = 1,nrbl
        block_sizes(1,rb(ibl)%id) = seam(ibl)%nvert
        block_sizes(2,rb(ibl)%id) = rb(ibl)%mx
        block_sizes(3,rb(ibl)%id) = rb(ibl)%my
      ENDDO

      DO ibl = nrbl+1,nbl
        block_sizes(1,tb(ibl)%id) = seam(ibl)%nvert
        block_sizes(2,tb(ibl)%id) = tb(ibl)%mvert
        block_sizes(3,tb(ibl)%id) = tb(ibl)%mcell
      ENDDO

      CALL mpi_allreduce(block_sizes,workvec,3*nbl_total,
     $     mpi_nim_int,mpi_sum,comm_layer,ierror)
      block_sizes = workvec

      DEALLOCATE(workvec)

      RETURN
      END SUBROUTINE parallel_block_init2


c-----------------------------------------------------------------------
c     subprogram 3. parallel_seam_init
c     setup seam communication data stuctures
c-----------------------------------------------------------------------

      SUBROUTINE parallel_seam_init(nqty)
      USE fields
      USE seam_storage_mod
      USE mpi_nim
      USE pardata
      IMPLICIT NONE
      integer(i4), intent(in) :: nqty
      integer(i4), dimension(:), allocatable :: procs,workvec
      integer(i4), dimension(:), allocatable :: extra_request
      integer(i4), dimension(:), allocatable :: extra_request2
      integer(i4), dimension(:,:), allocatable :: statuses
      integer(i4), dimension(:), allocatable :: data_block,data_vertex
      integer(i4), dimension(:), allocatable :: data_order
      integer(i4) :: i,imageproc,ibl,ivert,image,isend,irecv,count
      integer(i4) :: globalblock,ierror
      integer(i4) :: status(mpi_status_size)

c -------------------------
c setup send and self data structures
c -------------------------

c allocate a work vector of length nprocs and zero it

      allocate(procs(0:nprocs-1))
      procs = 0

c procs(n) = 1 if any of my seam data will be sent to proc n (include self)
c procs(n) = 0 otherwise
c  globalblock = global block # of an image point
c  imageproc = owner of the global block

      do ibl = 1,nbl
        do ivert = 1,seam(ibl)%nvert
          do image = 1,seam(ibl)%vertex(ivert)%nimage
            globalblock = seam(ibl)%vertex(ivert)%ptr2(1,image)
            imageproc = block2proc(globalblock)
            procs(imageproc) = 1
          enddo
        enddo
      enddo

c nsend = # of procs I will send seam-data to (exclude self)

      procs(node) = 0
      nsend = sum(procs)
      if (nsend > 0) allocate(send(nsend))

c procs(n) = # of seam-datums I will send to proc n (include self)

      procs = 0
      do ibl = 1,nbl
        do ivert = 1,seam(ibl)%nvert
          do image = 1,seam(ibl)%vertex(ivert)%nimage
            globalblock = seam(ibl)%vertex(ivert)%ptr2(1,image)
            imageproc = block2proc(globalblock)
            procs(imageproc) = procs(imageproc) + 1
          enddo
        enddo
      enddo

c use procs to allocate fields in send-data structure (to send to other procs)
c use procs(node) to allocate self-data structure
c  nself = # of seam-datums I must exchange between my own blocks
c  do not allocate self-data structure if nprocs = 1 or nself = 0
c after allocation, set send(isend)%count and nself to 0,
c  so can use as counters in next stage
c also after allocation, store which message (1:nsend) is sent to 
c  imageproc in procs(imageproc) for use in next stage

      isend = 0
      do imageproc = 0,nprocs-1
        if (node == imageproc) then
          if (nprocs > 1) then
            nself = procs(imageproc)
            if (nself > 0) allocate(self(nself))
            nself = 0
          endif
        else if (procs(imageproc) > 0) then
          isend = isend + 1
          send(isend)%proc = imageproc
          count = procs(imageproc)
          allocate(send(isend)%block(count))
          allocate(send(isend)%vertex(count))
          allocate(send(isend)%image(count))
          allocate(send(isend)%data(nqty*count))
          allocate(send(isend)%cdata(nqty*count))
          send(isend)%count = 0
          procs(imageproc) = isend
        endif
      enddo

c loop over all seam-datums and store ptr2 info in send-data structure
c datums being exchanged on-processor are stored in self-data structure
c  (don't bother if nprocs = 1)
c use send(isend)%count and nself as incremented ptrs to next available space,
c  when done they will have been reset to correct totals
c regarding order, _out is the receiving block, so out's order(image) is
c the location the _in data will hold in _out's seam_hold.

      do ibl = 1,nbl
        do ivert = 1,seam(ibl)%nvert
          do image = 1,seam(ibl)%vertex(ivert)%nimage
            globalblock = seam(ibl)%vertex(ivert)%ptr2(1,image)
            imageproc = block2proc(globalblock)
            if (node == imageproc) then
              if (nprocs > 1) then
                nself = nself + 1
                self(nself)%block_out = ibl
                self(nself)%vertex_out = ivert
                self(nself)%block_in = global2local(globalblock)
                self(nself)%vertex_in = 
     $               seam(ibl)%vertex(ivert)%ptr2(2,image)
                if (seam(ibl)%vertex(ivert)%nimage==1) then
                  self(nself)%order = 1
                else
                  self(nself)%order =
     $                 seam(ibl)%vertex(ivert)%order(image)
                endif
              endif
            else
              isend = procs(imageproc)
              count = send(isend)%count + 1
              send(isend)%block(count) = ibl
              send(isend)%vertex(count) = ivert
              send(isend)%image(count) = image
              send(isend)%count = count
            endif
          enddo
        enddo
      enddo

c -------------------------
c setup recv data structure
c -------------------------

c procs(n) = 1 if a seam-data message will be sent to proc n (exclude self)
c procs(n) = 0 otherwise

      procs = 0
      do isend = 1,nsend
        procs(send(isend)%proc) = 1
      enddo

c sum procs vector (in minimal way) across all procs
c nrecv = # of seam-data messages I will receive
c then are finished with procs

      allocate(workvec(nprocs))
      workvec = 1

      CALL mpi_reduce_scatter(procs,nrecv,workvec,
     $     mpi_nim_int,mpi_sum,mpi_comm_world,ierror)

      deallocate(workvec)
      deallocate(procs)

c allocate recv data structure

      if (nrecv > 0) allocate(recv(nrecv))

c send count of how many send-datums to each receiving proc

      do isend = 1,nsend
        CALL mpi_send(send(isend)%count,1,mpi_nim_int,
     $       send(isend)%proc,0,mpi_comm_world,ierror)
      enddo

c receive counts and allocate recv data structure for each incoming message

      do irecv = 1,nrecv
        CALL mpi_recv(count,1,mpi_nim_int,mpi_any_source,0,
     $       mpi_comm_world,status,ierror)
        recv(irecv)%proc = status(mpi_source)
        recv(irecv)%count = count
        allocate(recv(irecv)%block(count))
        allocate(recv(irecv)%vertex(count))
        allocate(recv(irecv)%order(count))
        allocate(recv(irecv)%data(nqty*count))
        allocate(recv(irecv)%cdata(nqty*count))
      enddo

c post receives for incoming block and vertex numbers
c need extra vector to store requests for 2nd message

      if (nrecv > 0) then
        allocate(recv_request(nrecv))
        allocate(extra_request(nrecv))
        allocate(extra_request2(nrecv))
      endif

      do irecv = 1,nrecv
        CALL mpi_irecv(recv(irecv)%block(1),recv(irecv)%count,
     $       mpi_nim_int,recv(irecv)%proc,0,mpi_comm_world,
     $       extra_request(irecv),ierror)
        CALL mpi_irecv(recv(irecv)%vertex(1),recv(irecv)%count,
     $       mpi_nim_int,recv(irecv)%proc,0,mpi_comm_world,
     $       extra_request2(irecv),ierror)
        CALL mpi_irecv(recv(irecv)%order(1),recv(irecv)%count,
     $       mpi_nim_int,recv(irecv)%proc,0,mpi_comm_world,
     $       recv_request(irecv),ierror)
      enddo

c synchonize to insure all receives have been posted

      CALL mpi_barrier(mpi_comm_world,ierror)

c send block and vertex #s for each send-datum to each receiving proc
c these are local block # and local vertex # for where the send-datum goes
c  on the receiving proc, computed by the sender
c order for each communicated datum is the same in for all representations.
c here, the sending block's order is in its own order(0), and that is
c communicated to the receiving block.

      do isend = 1,nsend
        count = send(isend)%count
        allocate(data_block(count))
        allocate(data_vertex(count))
        allocate(data_order(count))
        do i = 1,count
          ibl = send(isend)%block(i)
          ivert = send(isend)%vertex(i)
          image = send(isend)%image(i)
          globalblock = seam(ibl)%vertex(ivert)%ptr2(1,image)
          data_block(i) = global2local(globalblock)
          data_vertex(i) = seam(ibl)%vertex(ivert)%ptr2(2,image)
          data_order(i) = seam(ibl)%vertex(ivert)%order(0)
        enddo
        CALL mpi_send(data_block,count,mpi_nim_int,
     $       send(isend)%proc,0,mpi_comm_world,ierror)
        CALL mpi_send(data_vertex,count,mpi_nim_int,
     $       send(isend)%proc,0,mpi_comm_world,ierror)
        CALL mpi_send(data_order,count,mpi_nim_int,
     $       send(isend)%proc,0,mpi_comm_world,ierror)
        deallocate(data_block)
        deallocate(data_vertex)
        deallocate(data_order)
      enddo

c loop until all messages are received
c wait on 3rd message sent, since it will guarantee first two already arrived
c free extra message request vector

      if (nrecv > 0) then
        allocate(statuses(mpi_status_size,nrecv))
        call mpi_waitall(nrecv,recv_request,statuses,ierror)
        deallocate(statuses)
        deallocate(extra_request,extra_request2)
      endif

c if there are only two images at a vertex, have the communicated data
c placed in the first element of seam_hold.

      do irecv = 1,nrecv
        do i = 1,recv(irecv)%count
          ibl=recv(irecv)%block(i)
          ivert=recv(irecv)%vertex(i)
          if (seam(ibl)%vertex(ivert)%nimage==1) then
            recv(irecv)%order(i)=1
          endif
        enddo
      enddo

      RETURN
      END SUBROUTINE parallel_seam_init


c-----------------------------------------------------------------------
c     subprogram 4. parallel_seam_comm
c     knit seams together across procs
c-----------------------------------------------------------------------

      SUBROUTINE parallel_seam_comm(nqty)
      USE seam_storage_mod
      USE mpi_nim
      USE pardata
      IMPLICIT NONE
      integer(i4), intent(in) :: nqty
      integer(i4) :: ibl,jbl,nbl,ivert,irecv,isend
      integer(i4) :: i,j,jstart,jend,iself,jvert,ierror
      integer(i4) :: status(mpi_status_size)
      integer(i4) :: ni,image,ip

c use seam to determine the number of blocks instead of accessing fields

      nbl=SIZE(seam)

c post receives for incoming seam data

      do irecv = 1,nrecv
        CALL mpi_irecv(recv(irecv)%data(1),nqty*recv(irecv)%count,
     $       mpi_nim_real,recv(irecv)%proc,0,
     $       mpi_comm_world,recv_request(irecv),ierror)
      enddo

c synchronize to insure all receives have been posted
c this is not strictly necessary, could be faster to comment this out
c synchronization is done ONLY within a layer

c     CALL mpi_barrier(comm_layer,ierror)

c send message with each send-datum (from seam_in) to each receiving proc

      do isend = 1,nsend
        jstart = 1
        jend = jstart + nqty - 1
        do i = 1,send(isend)%count
          ibl = send(isend)%block(i)
          ivert = send(isend)%vertex(i)
          send(isend)%data(jstart:jend) =
     $         seam(ibl)%vertex(ivert)%seam_in(1:nqty)
          jstart = jstart + nqty
          jend = jend + nqty
        enddo
        CALL mpi_send(send(isend)%data,nqty*send(isend)%count,
     $       mpi_nim_real,send(isend)%proc,0,
     $       mpi_comm_world,ierror)
      enddo

c perform my own seam_in->seam_copies for pairs of pts I own
c  using self-data structure
c conceptually identical to edge_accumulate routine

      do iself = 1,nself
        ibl = self(iself)%block_out
        ivert = self(iself)%vertex_out
        jbl = self(iself)%block_in
        jvert = self(iself)%vertex_in
        image = self(iself)%order
        seam(ibl)%vertex(ivert)%seam_hold(image,1:nqty) =
     $       seam(jbl)%vertex(jvert)%seam_in(1:nqty)
      enddo

c loop until all messages are received
c as each seam-data message comes in, perform data->seam_hold copy

      do j = 1,nrecv
        CALL mpi_waitany(nrecv,recv_request,irecv,status,ierror)
        jstart = 1
        jend = jstart + nqty - 1
        do i = 1,recv(irecv)%count
          ibl = recv(irecv)%block(i)
          ivert = recv(irecv)%vertex(i)
          image = recv(irecv)%order(i)
          seam(ibl)%vertex(ivert)%seam_hold(image,1:nqty) =
     $         recv(irecv)%data(jstart:jend)
          jstart = jstart + nqty
          jend = jend + nqty
        enddo
      enddo

c now sum data from all images (including self) in a pre-determined order,
c so that round-off error will be the same for all representations.

      do ibl = 1,nbl
        do ivert = 1,seam(ibl)%nvert
          ni=seam(ibl)%vertex(ivert)%nimage
          select case(ni)
          case(0)
            seam(ibl)%vertex(ivert)%seam_out(1:nqty)=
     $         seam(ibl)%vertex(ivert)%seam_in(1:nqty)
          case(1)
            seam(ibl)%vertex(ivert)%seam_out(1:nqty)=
     $         seam(ibl)%vertex(ivert)%seam_in(1:nqty)
     $        +seam(ibl)%vertex(ivert)%seam_hold(1,1:nqty)
          case default
            image=seam(ibl)%vertex(ivert)%order(0)
            seam(ibl)%vertex(ivert)%seam_hold(image,1:nqty)=
     $        seam(ibl)%vertex(ivert)%seam_in(1:nqty)
            seam(ibl)%vertex(ivert)%seam_out(1:nqty)=
     $        sum(seam(ibl)%vertex(ivert)%seam_hold(:,1:nqty),1)
          end select
        enddo
      enddo

      RETURN
      END SUBROUTINE parallel_seam_comm


c-----------------------------------------------------------------------
c     subprogram 5. parallel_seam_comm_comp
c     knit seams together across procs for complex vertex quantities.
c-----------------------------------------------------------------------

      SUBROUTINE parallel_seam_comm_comp(nqty)
      USE seam_storage_mod
      USE mpi_nim
      USE pardata
      IMPLICIT NONE
      integer(i4), intent(in) :: nqty
      integer(i4) :: ibl,jbl,nbl,ivert,irecv,isend
      integer(i4) :: i,j,jstart,jend,iself,jvert,ierror
      integer(i4) :: status(mpi_status_size)
      integer(i4) :: ni,image,ip

c use seam to determine the number of blocks instead of accessing fields

      nbl=SIZE(seam)

c post receives for incoming seam data

      do irecv = 1,nrecv
        CALL mpi_irecv(recv(irecv)%cdata(1),nqty*recv(irecv)%count,
     $       mpi_nim_comp,recv(irecv)%proc,0,
     $       mpi_comm_world,recv_request(irecv),ierror)
      enddo

c synchronize to insure all receives have been posted
c this is not strictly necessary, could be faster to comment this out
c synchronization is done ONLY within a layer

c     CALL mpi_barrier(comm_layer,ierror)

c send message with each send-datum (from seam_in) to each receiving proc

      do isend = 1,nsend
        jstart = 1
        jend = jstart + nqty - 1
        do i = 1,send(isend)%count
          ibl = send(isend)%block(i)
          ivert = send(isend)%vertex(i)
          send(isend)%cdata(jstart:jend) =
     $         seam(ibl)%vertex(ivert)%seam_cin(1:nqty)
          jstart = jstart + nqty
          jend = jend + nqty
        enddo
        CALL mpi_send(send(isend)%cdata,nqty*send(isend)%count,
     $       mpi_nim_comp,send(isend)%proc,0,
     $       mpi_comm_world,ierror)
      enddo

c perform my own seam_in->seam_copies for pairs of pts I own
c  using self-data structure
c conceptually identical to edge_accumulate routine

      do iself = 1,nself
        ibl = self(iself)%block_out
        ivert = self(iself)%vertex_out
        jbl = self(iself)%block_in
        jvert = self(iself)%vertex_in
        image = self(iself)%order
        seam(ibl)%vertex(ivert)%seam_chold(image,1:nqty) =
     $       seam(jbl)%vertex(jvert)%seam_cin(1:nqty)
      enddo

c loop until all messages are received
c as each seam-data message comes in, perform data->seam_hold copy

      do j = 1,nrecv
        CALL mpi_waitany(nrecv,recv_request,irecv,status,ierror)
        jstart = 1
        jend = jstart + nqty - 1
        do i = 1,recv(irecv)%count
          ibl = recv(irecv)%block(i)
          ivert = recv(irecv)%vertex(i)
          image = recv(irecv)%order(i)
          seam(ibl)%vertex(ivert)%seam_chold(image,1:nqty) =
     $         recv(irecv)%cdata(jstart:jend)
          jstart = jstart + nqty
          jend = jend + nqty
        enddo
      enddo

c now sum data from all images (including self) in a pre-determined order,
c so that round-off error will be the same for all representations.

      do ibl = 1,nbl
        do ivert = 1,seam(ibl)%nvert
          ni=seam(ibl)%vertex(ivert)%nimage
          select case(ni)
          case(0)
            seam(ibl)%vertex(ivert)%seam_cout(1:nqty)=
     $         seam(ibl)%vertex(ivert)%seam_cin(1:nqty)
          case(1)
            seam(ibl)%vertex(ivert)%seam_cout(1:nqty)=
     $         seam(ibl)%vertex(ivert)%seam_cin(1:nqty)
     $        +seam(ibl)%vertex(ivert)%seam_chold(1,1:nqty)
          case default
            image=seam(ibl)%vertex(ivert)%order(0)
            seam(ibl)%vertex(ivert)%seam_chold(image,1:nqty)=
     $        seam(ibl)%vertex(ivert)%seam_cin(1:nqty)
            seam(ibl)%vertex(ivert)%seam_cout(1:nqty)=
     $        sum(seam(ibl)%vertex(ivert)%seam_chold(:,1:nqty),1)
          end select
        enddo
      enddo

      RETURN
      END SUBROUTINE parallel_seam_comm_comp


c-----------------------------------------------------------------------
c     subprogram 6. parallel_seg_init
c     setup segment communication data stuctures
c-----------------------------------------------------------------------

      SUBROUTINE parallel_seg_init(nqty,nqm,offmat)
      USE fields
      USE seam_storage_mod
      USE mpi_nim
      USE pardata
      IMPLICIT NONE
      integer(i4), intent(in) :: nqty,nqm,offmat
      integer(i4), dimension(:), allocatable :: procs,workvec
      integer(i4), dimension(:), allocatable :: extra_request
      integer(i4), dimension(:,:), allocatable :: statuses
      integer(i4), dimension(:), allocatable :: data_block
      integer(i4), dimension(:), allocatable :: data_segment
      integer(i4) :: i,imageproc,ibl,iseg,isend,irecv,count
      integer(i4) :: globalblock,ierror
      integer(i4) :: status(mpi_status_size)

c -------------------------
c setup sendseg and selfseg data structures
c -------------------------

c allocate a work vector of length nprocs and zero it

      allocate(procs(0:nprocs-1))
      procs = 0

c procs(n) = 1 if any of my seg data will be sent to proc n (include self)
c procs(n) = 0 otherwise
c  globalblock = global block # of an image segment
c  imageproc = owner of the global block
c no seg data if ptr(1) = 0

      do ibl = 1,nbl
        do iseg = 1,seam(ibl)%nvert
           globalblock = seam(ibl)%segment(iseg)%ptr(1)
           if (globalblock /= 0) then
             imageproc = block2proc(globalblock)
             procs(imageproc) = 1
           endif
        enddo
      enddo

c nsendseg = # of procs I will send segment-data to (exclude self)

      procs(node) = 0
      nsendseg = sum(procs)
      if (nsendseg > 0) allocate(sendseg(nsendseg))

c procs(n) = # of segment-datums I will send to proc n (include self)
c no seg data if ptr(1) = 0

      procs = 0
      do ibl = 1,nbl
        do iseg = 1,seam(ibl)%nvert
          globalblock = seam(ibl)%segment(iseg)%ptr(1)
          if (globalblock /= 0) then 
            imageproc = block2proc(globalblock)
            procs(imageproc) = procs(imageproc) + 1
         endif
        enddo
      enddo

c use procs to allocate fields in seg-data structure (to send to other procs)
c use procs(node) to allocate selfseg-data structure
c  nselfseg = # of seg-datums I must exchange between my own blocks
c  do not allocate selfseg-data structure if nprocs = 1 or nselfseg = 0
c after allocation, set sendseg(isend)%count and nselfseg to 0,
c  so can use as counters in next stage
c also after allocation, store which message (1:nsendseg) is sent to 
c  imageproc in procs(imageproc) for use in next stage

      isend = 0
      do imageproc = 0,nprocs-1
        if (node == imageproc) then
          if (nprocs > 1) then
            nselfseg = procs(imageproc)
            if (nselfseg > 0) allocate(selfseg(nselfseg))
            nselfseg = 0
          endif
        else if (procs(imageproc) > 0) then
          isend = isend + 1
          sendseg(isend)%proc = imageproc
          count = procs(imageproc)
          allocate(sendseg(isend)%block(count))
          allocate(sendseg(isend)%segment(count))
          allocate(sendseg(isend)%data(nqty*count))
          allocate(sendseg(isend)%cdata(nqty*count))
          allocate(sendseg(isend)%mat_data(offmat*nqm**2*count))
          allocate(sendseg(isend)%cmat_data(offmat*nqm**2*count))
          sendseg(isend)%count = 0
          procs(imageproc) = isend
        endif
      enddo

c loop over all seg-datums and store ptr info in seg-data structure
c datums being exchanged on-processor are stored in selfseg-data structure
c  (don't bother if nprocs = 1)
c use sendseg(isend)%count and nselfseg as incremented ptrs to next 
c  available space, when done they will have been reset to correct totals
c no seg data if ptr(1) = 0

      do ibl = 1,nbl
        do iseg = 1,seam(ibl)%nvert
          globalblock = seam(ibl)%segment(iseg)%ptr(1)
          if (globalblock /= 0) then 
            imageproc = block2proc(globalblock)
            if (node == imageproc) then
              if (nprocs > 1) then
                nselfseg = nselfseg + 1
                selfseg(nselfseg)%block_out = ibl
                selfseg(nselfseg)%segment_out = iseg
                selfseg(nselfseg)%block_in = global2local(globalblock)
                selfseg(nselfseg)%segment_in = 
     $               seam(ibl)%segment(iseg)%ptr(2)
              endif
            else
              isend = procs(imageproc)
              count = sendseg(isend)%count + 1
              sendseg(isend)%block(count) = ibl
              sendseg(isend)%segment(count) = iseg
              sendseg(isend)%count = count
            endif
          endif
        enddo
      enddo

c -------------------------
c setup recvseg data structure
c -------------------------

c procs(n) = 1 if a seg-data message will be sent to proc n (exclude self)
c procs(n) = 0 otherwise

      procs = 0
      do isend = 1,nsendseg
        procs(sendseg(isend)%proc) = 1
      enddo

c sum procs vector (in minimal way) across all procs
c nrecvseg = # of seg-data messages I will receive
c then are finished with procs

      allocate(workvec(nprocs))
      workvec = 1

      CALL mpi_reduce_scatter(procs,nrecvseg,workvec,
     $     mpi_nim_int,mpi_sum,mpi_comm_world,ierror)

      deallocate(workvec)
      deallocate(procs)

c allocate recvseg data structure

      if (nrecvseg > 0) allocate(recvseg(nrecvseg))

c send count of how many seg-datums to each receiving proc

      do isend = 1,nsendseg
        CALL mpi_send(sendseg(isend)%count,1,mpi_nim_int,
     $       sendseg(isend)%proc,0,mpi_comm_world,ierror)
      enddo

c receive counts and allocate recvseg data structure for each incoming message

      do irecv = 1,nrecvseg
        CALL mpi_recv(count,1,mpi_nim_int,mpi_any_source,0,
     $       mpi_comm_world,status,ierror)
        recvseg(irecv)%proc = status(mpi_source)
        recvseg(irecv)%count = count
        allocate(recvseg(irecv)%block(count))
        allocate(recvseg(irecv)%segment(count))
        allocate(recvseg(irecv)%data(nqty*count))
        allocate(recvseg(irecv)%cdata(nqty*count))
        allocate(recvseg(irecv)%mat_data(offmat*nqm**2*count))
        allocate(recvseg(irecv)%cmat_data(offmat*nqm**2*count))
      enddo

c post receives for incoming block and segment numbers
c need extra vector to store requests for 2nd message

      if (nrecvseg > 0) then
        allocate(recvseg_request(nrecvseg))
        allocate(extra_request(nrecvseg))
      endif

      do irecv = 1,nrecvseg
        CALL mpi_irecv(recvseg(irecv)%block(1),recvseg(irecv)%count,
     $       mpi_nim_int,recvseg(irecv)%proc,0,mpi_comm_world,
     $       extra_request(irecv),ierror)
        CALL mpi_irecv(recvseg(irecv)%segment(1),recvseg(irecv)%count,
     $       mpi_nim_int,recvseg(irecv)%proc,0,mpi_comm_world,
     $       recvseg_request(irecv),ierror)
      enddo

c synchonize to insure all receives have been posted

      CALL mpi_barrier(mpi_comm_world,ierror)

c send block and segment #s for each seg-datum to each receiving proc
c these are local block # and local segment # for where the seg-datum goes
c  on the receiving proc, computed by the sender

      do isend = 1,nsendseg
        count = sendseg(isend)%count
        allocate(data_block(count))
        allocate(data_segment(count))
        do i = 1,count
          ibl = sendseg(isend)%block(i)
          iseg = sendseg(isend)%segment(i)
          globalblock = seam(ibl)%segment(iseg)%ptr(1)
          data_block(i) = global2local(globalblock)
          data_segment(i) = seam(ibl)%segment(iseg)%ptr(2)
        enddo
        CALL mpi_send(data_block,count,mpi_nim_int,
     $       sendseg(isend)%proc,0,mpi_comm_world,ierror)
        CALL mpi_send(data_segment,count,mpi_nim_int,
     $       sendseg(isend)%proc,0,mpi_comm_world,ierror)
        deallocate(data_block)
        deallocate(data_segment)
      enddo

c loop until all messages are received
c wait on 2nd message sent, since it will guarantee first one already arrived
c free extra message request vector

      if (nrecvseg > 0) then
        allocate(statuses(mpi_status_size,nrecvseg))
        call mpi_waitall(nrecvseg,recvseg_request,statuses,ierror)
        deallocate(statuses)
        deallocate(extra_request)
      endif

      RETURN
      END SUBROUTINE parallel_seg_init


c-----------------------------------------------------------------------
c     subprogram 7. parallel_seg_comm
c     knit segments together across procs for normal edge-centered data.
c-----------------------------------------------------------------------

      SUBROUTINE parallel_seg_comm(nqty,nside)
      USE seam_storage_mod
      USE mpi_nim
      USE pardata
      IMPLICIT NONE
      integer(i4), intent(in) :: nqty,nside
      integer(i4) :: ibl,jbl,iseg,irecv,isend,nqns,iqi,iqo
      integer(i4) :: i,j,jb,jstart,jend,iself,jseg,ierror
      integer(i4) :: status(mpi_status_size)

c post receives for incoming segment data

      nqns=nqty*nside
      do irecv = 1,nrecvseg
        CALL mpi_irecv(recvseg(irecv)%data(1),
     $        nqns*recvseg(irecv)%count,
     $        mpi_nim_real,recvseg(irecv)%proc,0,
     $        mpi_comm_world,recvseg_request(irecv),ierror)
      enddo

c synchonize to insure all receives have been posted
c this is not strictly necessary, could be faster to comment this out
c synchronization is done ONLY within a layer

c     CALL mpi_barrier(comm_layer,ierror)

c send message with each seg-datum (from seam_in) to each receiving proc

      do isend = 1,nsendseg
        jstart = 1
        jend = jstart + nqns - 1
        do i = 1,sendseg(isend)%count
          ibl = sendseg(isend)%block(i)
          iseg = sendseg(isend)%segment(i)
          sendseg(isend)%data(jstart:jend) =
     $         seam(ibl)%segment(iseg)%seam_in(1:nqns)
          jstart = jstart + nqns
          jend = jend + nqns
        enddo
        CALL mpi_send(sendseg(isend)%data,nqns*sendseg(isend)%count,
     $       mpi_nim_real,sendseg(isend)%proc,0,
     $       mpi_comm_world,ierror)
      enddo

c perform my own seam_in->seam_out summations for pairs of segments I own
c  using selfseg-data structure
c conceptually identical to edge_seg_accumulate routine, now with
c  side-basis ordering reversed between adjacent blocks

      do iself = 1,nselfseg
        ibl = selfseg(iself)%block_out
        iseg = selfseg(iself)%segment_out
        jbl = selfseg(iself)%block_in
        jseg = selfseg(iself)%segment_in
        iqi = 1
        iqo = nqty*(nside-1) + 1
        do jb = 1,nside
          seam(ibl)%segment(iseg)%seam_out(iqo:iqo+nqty-1) =
     $         seam(ibl)%segment(iseg)%seam_out(iqo:iqo+nqty-1) +
     $         seam(jbl)%segment(jseg)%seam_in(iqi:iqi+nqty-1)
          iqi = iqi + nqty
          iqo = iqo - nqty
        enddo
      enddo

c loop until all messages are received
c as each seg-data message comes in, perform data->seam_out summation

      do j = 1,nrecvseg
        CALL mpi_waitany(nrecvseg,recvseg_request,irecv,status,ierror)
        jstart = 1
        jend = jstart + nqty - 1
        do i = 1,recvseg(irecv)%count
          ibl = recvseg(irecv)%block(i)
          iseg = recvseg(irecv)%segment(i)
          iqo = nqty*(nside-1) + 1
          do jb = 1,nside
            seam(ibl)%segment(iseg)%seam_out(iqo:iqo+nqty-1) =
     $           seam(ibl)%segment(iseg)%seam_out(iqo:iqo+nqty-1) +
     $           recvseg(irecv)%data(jstart:jend)
            iqo = iqo - nqty
            jstart = jstart + nqty
            jend = jend + nqty
          enddo
        enddo
      enddo

      RETURN
      END SUBROUTINE parallel_seg_comm


c-----------------------------------------------------------------------
c     subprogram 8. parallel_seg_comm_comp
c     knit segments together across procs for complex edge-centered
c     data.
c-----------------------------------------------------------------------

      SUBROUTINE parallel_seg_comm_comp(nqty,nside)
      USE seam_storage_mod
      USE mpi_nim
      USE pardata
      IMPLICIT NONE
      integer(i4), intent(in) :: nqty,nside
      integer(i4) :: ibl,jbl,iseg,irecv,isend,nqns,iqi,iqo
      integer(i4) :: i,j,jb,jstart,jend,iself,jseg,ierror
      integer(i4) :: status(mpi_status_size)

c post receives for incoming segment data

      nqns=nqty*nside
      do irecv = 1,nrecvseg
        CALL mpi_irecv(recvseg(irecv)%cdata(1),
     $        nqns*recvseg(irecv)%count,
     $        mpi_nim_comp,recvseg(irecv)%proc,0,
     $        mpi_comm_world,recvseg_request(irecv),ierror)
      enddo

c synchonize to insure all receives have been posted
c this is not strictly necessary, could be faster to comment this out
c synchronization is done ONLY within a layer

c     CALL mpi_barrier(comm_layer,ierror)

c send message with each seg-datum (from seam_in) to each receiving proc

      do isend = 1,nsendseg
        jstart = 1
        jend = jstart + nqns - 1
        do i = 1,sendseg(isend)%count
          ibl = sendseg(isend)%block(i)
          iseg = sendseg(isend)%segment(i)
          sendseg(isend)%cdata(jstart:jend) =
     $         seam(ibl)%segment(iseg)%seam_cin(1:nqns)
          jstart = jstart + nqns
          jend = jend + nqns
        enddo
        CALL mpi_send(sendseg(isend)%cdata,nqns*sendseg(isend)%count,
     $       mpi_nim_comp,sendseg(isend)%proc,0,
     $       mpi_comm_world,ierror)
      enddo

c perform my own seam_in->seam_out summations for pairs of segments I own
c  using selfseg-data structure
c conceptually identical to edge_seg_accumulate routine, now with
c  side-basis ordering reversed between adjacent blocks

      do iself = 1,nselfseg
        ibl = selfseg(iself)%block_out
        iseg = selfseg(iself)%segment_out
        jbl = selfseg(iself)%block_in
        jseg = selfseg(iself)%segment_in
        iqi = 1
        iqo = nqty*(nside-1) + 1
        do jb = 1,nside
          seam(ibl)%segment(iseg)%seam_cout(iqo:iqo+nqty-1) =
     $         seam(ibl)%segment(iseg)%seam_cout(iqo:iqo+nqty-1) +
     $         seam(jbl)%segment(jseg)%seam_cin(iqi:iqi+nqty-1)
          iqi = iqi + nqty
          iqo = iqo - nqty
        enddo
      enddo

c loop until all messages are received
c as each seg-data message comes in, perform data->seam_out summation

      do j = 1,nrecvseg
        CALL mpi_waitany(nrecvseg,recvseg_request,irecv,status,ierror)
        jstart = 1
        jend = jstart + nqty - 1
        do i = 1,recvseg(irecv)%count
          ibl = recvseg(irecv)%block(i)
          iseg = recvseg(irecv)%segment(i)
          iqo = nqty*(nside-1) + 1
          do jb = 1,nside
            seam(ibl)%segment(iseg)%seam_cout(iqo:iqo+nqty-1) =
     $           seam(ibl)%segment(iseg)%seam_cout(iqo:iqo+nqty-1) +
     $           recvseg(irecv)%cdata(jstart:jend)
            iqo = iqo - nqty
            jstart = jstart + nqty
            jend = jend + nqty
          enddo
        enddo
      enddo

      RETURN
      END SUBROUTINE parallel_seg_comm_comp


c-----------------------------------------------------------------------
c     subprogram 9. parallel_matseg_comm
c     knit segments together across procs for off-diagonal matrix
c     elements.
c-----------------------------------------------------------------------

      SUBROUTINE parallel_matseg_comm(nqty,nmat)
      USE fields
      USE seam_storage_mod
      USE mpi_nim
      USE pardata
      IMPLICIT NONE
      integer(i4), intent(in) :: nqty,nmat
      integer(i4) :: ibl,jbl,iseg,irecv,isend,nqcom
      integer(i4) :: i,j,jstart,jend,iself,jseg,ierror,im
      integer(i4) :: status(mpi_status_size)

      nqcom = nmat*nqty**2

c post receives for incoming segment data

      do irecv = 1,nrecvseg
        CALL mpi_irecv(recvseg(irecv)%mat_data(1),
     $        nqcom*recvseg(irecv)%count,
     $        mpi_nim_real,recvseg(irecv)%proc,0,
     $        mpi_comm_world,recvseg_request(irecv),ierror)
      enddo

c synchonize to insure all receives have been posted
c this is not strictly necessary, could be faster to comment this out
c synchronization is done ONLY within a layer

c     CALL mpi_barrier(comm_layer,ierror)

c send message with each seg-datum (from seam_mat_in) to each receiving
c proc

      do isend = 1,nsendseg
        jstart = 1
        jend = jstart + nqcom - 1
        do i = 1,sendseg(isend)%count
          ibl = sendseg(isend)%block(i)
          iseg = sendseg(isend)%segment(i)
          sendseg(isend)%mat_data(jstart:jend) = RESHAPE(
     $      seam(ibl)%segment(iseg)%seam_mat_in(1:nqty,1:nqty,1:nmat),
     $      (/nqcom/) )
          jstart = jstart + nqcom
          jend = jend + nqcom
        enddo
        CALL mpi_send(sendseg(isend)%mat_data,
     $       nqcom*sendseg(isend)%count,
     $       mpi_nim_real,sendseg(isend)%proc,0,
     $       mpi_comm_world,ierror)
      enddo

c perform my own seam_in->seam_out summations for pairs of segments I own
c  using selfseg-data structure
c conceptually identical to edge_seg_accumulate routine

      do iself = 1,nselfseg
        ibl = selfseg(iself)%block_out
        iseg = selfseg(iself)%segment_out
        jbl = selfseg(iself)%block_in
        jseg = selfseg(iself)%segment_in
        do im=1,nmat,2
          seam(ibl)%segment(iseg)%seam_mat_out(1:nqty,1:nqty,im) =
     $         seam(ibl)%segment(iseg)%seam_mat_out(1:nqty,1:nqty,im) +
     $         seam(jbl)%segment(jseg)%seam_mat_in(1:nqty,1:nqty,im+1)
          seam(ibl)%segment(iseg)%seam_mat_out(1:nqty,1:nqty,im+1) =
     $         seam(ibl)%segment(iseg)%seam_mat_out(1:nqty,1:nqty,im+1)+
     $         seam(jbl)%segment(jseg)%seam_mat_in(1:nqty,1:nqty,im)
        enddo
      enddo

c loop until all messages are received
c as each seg-data message comes in, perform data->seam_out summation

      do j = 1,nrecvseg
        CALL mpi_waitany(nrecvseg,recvseg_request,irecv,status,ierror)
        jstart = 1
        jend = jstart + nqty**2 - 1
        do i = 1,recvseg(irecv)%count
          ibl = recvseg(irecv)%block(i)
          iseg = recvseg(irecv)%segment(i)
          do im=1,nmat,2
            seam(ibl)%segment(iseg)%seam_mat_out(1:nqty,1:nqty,im+1) =
     $        seam(ibl)%segment(iseg)%seam_mat_out(1:nqty,1:nqty,im+1) +
     $        RESHAPE( recvseg(irecv)%mat_data(jstart:jend),
     $                 (/nqty,nqty/) )
            jstart = jstart + nqty**2
            jend = jend + nqty**2
            seam(ibl)%segment(iseg)%seam_mat_out(1:nqty,1:nqty,im) =
     $        seam(ibl)%segment(iseg)%seam_mat_out(1:nqty,1:nqty,im) +
     $        RESHAPE( recvseg(irecv)%mat_data(jstart:jend),
     $                 (/nqty,nqty/) )
            jstart = jstart + nqty**2
            jend = jend + nqty**2
          enddo
        enddo
      enddo

      RETURN
      END SUBROUTINE parallel_matseg_comm


c-----------------------------------------------------------------------
c     subprogram 10. parallel_matseg_comm_comp
c     knit segments together across procs for complex off-diagonal
c     matrix elements.
c-----------------------------------------------------------------------

      SUBROUTINE parallel_matseg_comm_comp(nqty,nmat)
      USE fields
      USE seam_storage_mod
      USE mpi_nim
      USE pardata
      IMPLICIT NONE
      integer(i4), intent(in) :: nqty,nmat
      integer(i4) :: ibl,jbl,iseg,irecv,isend,nqcom
      integer(i4) :: i,j,jstart,jend,iself,jseg,ierror,im
      integer(i4) :: status(mpi_status_size)

      nqcom = nmat*nqty**2

c post receives for incoming segment data

      do irecv = 1,nrecvseg
        CALL mpi_irecv(recvseg(irecv)%cmat_data(1),
     $        nqcom*recvseg(irecv)%count,
     $        mpi_nim_comp,recvseg(irecv)%proc,0,
     $        mpi_comm_world,recvseg_request(irecv),ierror)
      enddo

c synchonize to insure all receives have been posted
c this is not strictly necessary, could be faster to comment this out
c synchronization is done ONLY within a layer

c     CALL mpi_barrier(comm_layer,ierror)

c send message with each seg-datum (from seam_mat_cin) to each receiving
c proc

      do isend = 1,nsendseg
        jstart = 1
        jend = jstart + nqcom - 1
        do i = 1,sendseg(isend)%count
          ibl = sendseg(isend)%block(i)
          iseg = sendseg(isend)%segment(i)
          sendseg(isend)%cmat_data(jstart:jend) = RESHAPE(
     $      seam(ibl)%segment(iseg)%seam_mat_cin(1:nqty,1:nqty,1:nmat),
     $      (/nqcom/) )
          jstart = jstart + nqcom
          jend = jend + nqcom
        enddo
        CALL mpi_send(sendseg(isend)%cmat_data,
     $       nqcom*sendseg(isend)%count,
     $       mpi_nim_comp,sendseg(isend)%proc,0,
     $       mpi_comm_world,ierror)
      enddo

c perform my own seam_in->seam_out summations for pairs of segments I own
c  using selfseg-data structure
c conceptually identical to edge_seg_accumulate routine

      do iself = 1,nselfseg
        ibl = selfseg(iself)%block_out
        iseg = selfseg(iself)%segment_out
        jbl = selfseg(iself)%block_in
        jseg = selfseg(iself)%segment_in
        do im=1,nmat,2
          seam(ibl)%segment(iseg)%seam_mat_cout(1:nqty,1:nqty,im) =
     $       seam(ibl)%segment(iseg)%seam_mat_cout(1:nqty,1:nqty,im) +
     $       seam(jbl)%segment(jseg)%seam_mat_cin(1:nqty,1:nqty,im+1)
          seam(ibl)%segment(iseg)%seam_mat_cout(1:nqty,1:nqty,im+1) =
     $       seam(ibl)%segment(iseg)%seam_mat_cout(1:nqty,1:nqty,im+1) +
     $       seam(jbl)%segment(jseg)%seam_mat_cin(1:nqty,1:nqty,im)
        enddo
      enddo

c loop until all messages are received
c as each seg-data message comes in, perform data->seam_out summation

      do j = 1,nrecvseg
        CALL mpi_waitany(nrecvseg,recvseg_request,irecv,status,ierror)
        jstart = 1
        jend = jstart + nqty**2 - 1
        do i = 1,recvseg(irecv)%count
          ibl = recvseg(irecv)%block(i)
          iseg = recvseg(irecv)%segment(i)
          do im=1,nmat,2
            seam(ibl)%segment(iseg)%seam_mat_cout(1:nqty,1:nqty,im+1) =
     $        seam(ibl)%segment(iseg)%seam_mat_cout(1:nqty,1:nqty,im+1)+
     $        RESHAPE( recvseg(irecv)%cmat_data(jstart:jend),
     $                 (/nqty,nqty/) )
            jstart = jstart + nqty**2
            jend = jend + nqty**2
            seam(ibl)%segment(iseg)%seam_mat_cout(1:nqty,1:nqty,im) =
     $        seam(ibl)%segment(iseg)%seam_mat_cout(1:nqty,1:nqty,im) +
     $        RESHAPE( recvseg(irecv)%cmat_data(jstart:jend),
     $                 (/nqty,nqty/) )
            jstart = jstart + nqty**2
            jend = jend + nqty**2
          enddo
        enddo
      enddo

      RETURN
      END SUBROUTINE parallel_matseg_comm_comp


c-----------------------------------------------------------------------
c     subprogram 11. parallel_line_init.
c     set up communication information for global line preconditioning
c     via data decomposition swaps.  note that this is needed for
c     serial computations that use global line preconditioning, too.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE parallel_line_init(rb,nrbl,nbl,poly_deg)
      USE local
      USE seam_storage_mod
      USE edge
      USE rblock_type_mod
      USE mpi_nim
      USE pardata
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: nrbl,nbl,poly_deg
      TYPE(rblock_type), DIMENSION(nrbl), INTENT(IN) :: rb

      TYPE :: chain_type
        INTEGER(i4) :: nlinks
        TYPE(chain_link_type), POINTER :: first,last,current
      END TYPE chain_type

      TYPE :: chain_link_type
        INTEGER(i4), DIMENSION(5) :: blinfo
        TYPE(chain_link_type), POINTER :: next,prev
      END TYPE chain_link_type

      TYPE(chain_type), DIMENSION(nrbl) :: chainx,chainy
      INTEGER(i4) :: ibl,ibg,iv,ip,ir,ibmax,ibmin,flag,itmp,
     $               ppst,ppen,past,paen
      INTEGER(i4), DIMENSION(nrbl) :: mxb,myb,glb_bl
      INTEGER :: ierror
      CHARACTER(80) :: msg

c=======================================================================
c     phase I determines what blocks appear in the global grid lines
c     using linked lists and seam communication.
c=======================================================================
c-----------------------------------------------------------------------
c     find the global block index for each local rblock.
c-----------------------------------------------------------------------
      DO ibg=1,SIZE(block2proc)
        IF (block2proc(ibg)==node) THEN
          ibl=global2local(ibg)
          IF (ibl<=nrbl) glb_bl(ibl)=ibg
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     find the rblock dimensions.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        mxb(ibl)=rb(ibl)%mx
        myb(ibl)=rb(ibl)%my
      ENDDO
c-----------------------------------------------------------------------
c     load local links into the chains as if each rblock were a separate
c     chain.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        chainx(ibl)%nlinks=1
        ALLOCATE(chainx(ibl)%current)
        chainx(ibl)%current%blinfo(1:5)=
     $    (/node,glb_bl(ibl),ibl,mxb(ibl),poly_deg*myb(ibl)/)
        NULLIFY(chainx(ibl)%current%next)
        NULLIFY(chainx(ibl)%current%prev)
        chainx(ibl)%first=>chainx(ibl)%current
        chainx(ibl)%last=>chainx(ibl)%current
        chainy(ibl)%nlinks=1
        ALLOCATE(chainy(ibl)%current)
        chainy(ibl)%current%blinfo(1:5)=
     $    (/node,glb_bl(ibl),ibl,myb(ibl),poly_deg*mxb(ibl)/)
        NULLIFY(chainy(ibl)%current%next)
        NULLIFY(chainy(ibl)%current%prev)
        chainy(ibl)%first=>chainy(ibl)%current
        chainy(ibl)%last=>chainy(ibl)%current
      ENDDO
c-----------------------------------------------------------------------
c     pass block sizes and labels to the left until all encounter 
c     the left border.
c-----------------------------------------------------------------------
      left: DO
        DO ibl=1,nbl
          DO iv=1,seam(ibl)%nvert
            seam(ibl)%vertex(iv)%seam_in(1:5)=0
          ENDDO
        ENDDO
        DO ibl=1,nrbl
          DO iv=2*mxb(ibl)+myb(ibl),2*(mxb(ibl)+myb(ibl))-1
            seam(ibl)%vertex(iv)%seam_in(1:5)=chainx(ibl)%current%blinfo
          ENDDO
        ENDDO
        CALL edge_network(5_i4,0_i4,0_i4,.false.)
c-----------------------------------------------------------------------
c       check that the passed information is from one block, if so
c       and it's a new block, push the stack.
c-----------------------------------------------------------------------
        flag=0
        DO ibl=1,nrbl
          ibmax=-HUGE(ibmax)
          ibmin= HUGE(ibmax)   
          DO iv=mxb(ibl)+1,mxb(ibl)+myb(ibl)
            ibg=NINT(seam(ibl)%vertex(iv)%seam_out(2))
            ibmax=MAX(ibmax,ibg)
            ibmin=MIN(ibmin,ibg)
          ENDDO
          IF (ibmax/=ibmin) THEN
            WRITE(msg,'(3(a,i3))')
     $        'Parallel_line_init: linex with block',ibl,
     $        ' encounters blocks ',ibmax,' and ', ibmin
            CALL nim_stop(msg)
          ELSE IF (ibg==chainx(ibl)%first%blinfo(2)) THEN
          ELSE IF (ibg>0) THEN
            IF (ibg/=chainx(ibl)%current%blinfo(2)) THEN
              ALLOCATE(chainx(ibl)%current%next)
              chainx(ibl)%last=>chainx(ibl)%current%next
              chainx(ibl)%last%prev=>chainx(ibl)%current
              NULLIFY(chainx(ibl)%last%next)
              chainx(ibl)%last%blinfo=
     $          NINT(seam(ibl)%vertex(mxb(ibl)+1)%seam_out(1:5))
              chainx(ibl)%current=>chainx(ibl)%last
              chainx(ibl)%nlinks=chainx(ibl)%nlinks+1
              flag=1
            ENDIF
          ENDIF
        ENDDO
        IF (nprocs>1) THEN
          CALL mpi_allreduce(flag,itmp,1,mpi_nim_int,mpi_sum,
     $                       mpi_comm_world,ierror)
          flag=itmp
        ENDIF
        IF (flag==0) EXIT left
      ENDDO left
c-----------------------------------------------------------------------
c-TMP write info after left passes.  write separate nodes to sep. files
c-----------------------------------------------------------------------
c     CALL parallel_line_write(chainx,linex,'x','After left passes.')
c-----------------------------------------------------------------------
c     pass block sizes and labels to the right until all encounter 
c     the left border.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        chainx(ibl)%current=>chainx(ibl)%first
      ENDDO
      right: DO
        DO ibl=1,nbl
          DO iv=1,seam(ibl)%nvert
            seam(ibl)%vertex(iv)%seam_in(1:5)=0
          ENDDO
        ENDDO
        DO ibl=1,nrbl
          DO iv=mxb(ibl)+1,mxb(ibl)+myb(ibl)
            seam(ibl)%vertex(iv)%seam_in(1:5)=chainx(ibl)%current%blinfo
          ENDDO
        ENDDO
        CALL edge_network(5_i4,0_i4,0_i4,.false.)
c-----------------------------------------------------------------------
c       check that the passed information is from one block, if so
c       and it's a new block, push the stack.
c-----------------------------------------------------------------------
        flag=0
        DO ibl=1,nrbl
          ibmax=-HUGE(ibmax)
          ibmin= HUGE(ibmax)   
          DO iv=2*mxb(ibl)+myb(ibl),2*(mxb(ibl)+myb(ibl))-1
            ibg=NINT(seam(ibl)%vertex(iv)%seam_out(2))
            ibmax=MAX(ibmax,ibg)
            ibmin=MIN(ibmin,ibg)
          ENDDO
          IF (ibmax/=ibmin) THEN
            WRITE(msg,'(3(a,i3))')
     $        'Parallel_line_init: linex with block',ibl,
     $        ' encounters blocks ',ibmax,' and ', ibmin
            CALL nim_stop(msg)
          ELSE IF (ibg==chainx(ibl)%last%blinfo(2)) THEN
            chainx(ibl)%first%prev=>chainx(ibl)%last
          ELSE IF (rb(ibl)%degenerate) THEN
          ELSE IF (ibg>0) THEN
            IF (ibg/=chainx(ibl)%current%blinfo(2)) THEN
              ALLOCATE(chainx(ibl)%current%prev)
              chainx(ibl)%first=>chainx(ibl)%current%prev
              chainx(ibl)%first%next=>chainx(ibl)%current
              NULLIFY(chainx(ibl)%first%prev)
              chainx(ibl)%first%blinfo=
     $          NINT(seam(ibl)%vertex(seam(ibl)%nvert-1)%seam_out(1:5))
              chainx(ibl)%current=>chainx(ibl)%first
              chainx(ibl)%nlinks=chainx(ibl)%nlinks+1
              flag=1
            ENDIF
          ENDIF
        ENDDO
        IF (nprocs>1) THEN
          CALL mpi_allreduce(flag,itmp,1,mpi_nim_int,mpi_sum,
     $                       mpi_comm_world,ierror)
          flag=itmp
        ENDIF
        IF (flag==0) EXIT right
      ENDDO right
c-----------------------------------------------------------------------
c     check for periodic x-lines and check consistency of duplicate
c     chains on the same processor.
c-----------------------------------------------------------------------
      CALL parallel_line_perset(chainx)
      CALL parallel_line_check(chainx,'x')
c-----------------------------------------------------------------------
c-TMP write info after right passes.  write separate nodes to sep. files
c-----------------------------------------------------------------------
c     CALL parallel_line_write(chainx,linex,'x','After right passes.')
c-----------------------------------------------------------------------
c     pass block sizes and labels downward until all encounter 
c     the top border.
c-----------------------------------------------------------------------
      down: DO
        DO ibl=1,nbl
          DO iv=1,seam(ibl)%nvert
            seam(ibl)%vertex(iv)%seam_in(1:5)=0
          ENDDO
        ENDDO
        DO ibl=1,nrbl
          DO iv=1,mxb(ibl)
            seam(ibl)%vertex(iv)%seam_in(1:5)=chainy(ibl)%current%blinfo
          ENDDO
        ENDDO
        CALL edge_network(5_i4,0_i4,0_i4,.false.)
c-----------------------------------------------------------------------
c       check that the passed information is from one block, if so
c       and it's a new block, push the stack.
c-----------------------------------------------------------------------
        flag=0
        DO ibl=1,nrbl
          ibmax=-HUGE(ibmax)
          ibmin= HUGE(ibmax)   
          DO iv=mxb(ibl)+myb(ibl),2*mxb(ibl)+myb(ibl)-1
            ibg=NINT(seam(ibl)%vertex(iv)%seam_out(2))
            ibmax=MAX(ibmax,ibg)
            ibmin=MIN(ibmin,ibg)
          ENDDO
          IF (ibmax/=ibmin) THEN
            WRITE(msg,'(3(a,i3))')
     $        'Parallel_line_init: liney with block',ibl,
     $        ' encounters blocks ',ibmax,' and ', ibmin
            CALL nim_stop(msg)
          ELSE IF (ibg==chainy(ibl)%first%blinfo(2)) THEN
          ELSE IF (ibg>0) THEN
            IF (ibg/=chainy(ibl)%current%blinfo(2)) THEN
              ALLOCATE(chainy(ibl)%current%next)
              chainy(ibl)%last=>chainy(ibl)%current%next
              chainy(ibl)%last%prev=>chainy(ibl)%current
              NULLIFY(chainy(ibl)%last%next)
              chainy(ibl)%last%blinfo=
     $          NINT(seam(ibl)%vertex(mxb(ibl)+myb(ibl))%seam_out(1:5))
              chainy(ibl)%current=>chainy(ibl)%last
              chainy(ibl)%nlinks=chainy(ibl)%nlinks+1
              flag=1
            ENDIF
          ENDIF
        ENDDO
        IF (nprocs>1) THEN
          CALL mpi_allreduce(flag,itmp,1,mpi_nim_int,mpi_sum,
     $                       mpi_comm_world,ierror)
          flag=itmp
        ENDIF
        IF (flag==0) EXIT down
      ENDDO down
c-----------------------------------------------------------------------
c-TMP write info after downward passes.
c-----------------------------------------------------------------------
c     CALL parallel_line_write(chainy,liney,'y','After down passes.')
c-----------------------------------------------------------------------
c     pass block sizes and labels upward until all encounter 
c     the bottom border.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        chainy(ibl)%current=>chainy(ibl)%first
      ENDDO
      up: DO
        DO ibl=1,nbl
          DO iv=1,seam(ibl)%nvert
            seam(ibl)%vertex(iv)%seam_in(1:5)=0
          ENDDO
        ENDDO
        DO ibl=1,nrbl
          DO iv=mxb(ibl)+myb(ibl),2*mxb(ibl)+myb(ibl)-1
            seam(ibl)%vertex(iv)%seam_in(1:5)=chainy(ibl)%current%blinfo
          ENDDO
        ENDDO
        CALL edge_network(5_i4,0_i4,0_i4,.false.)
c-----------------------------------------------------------------------
c       check that the passed information is from one block, if so
c       and it's a new block, push the stack.
c-----------------------------------------------------------------------
        flag=0
        DO ibl=1,nrbl
          ibmax=-HUGE(ibmax)
          ibmin= HUGE(ibmax)   
          DO iv=1,mxb(ibl)
            ibg=NINT(seam(ibl)%vertex(iv)%seam_out(2))
            ibmax=MAX(ibmax,ibg)
            ibmin=MIN(ibmin,ibg)
          ENDDO
          IF (ibmax/=ibmin) THEN
            WRITE(msg,'(3(a,i3))')
     $        'Parallel_line_init: liney with block',ibl,
     $        ' encounters blocks ',ibmax,' and ', ibmin
            CALL nim_stop(msg)
          ELSE IF (ibg==chainy(ibl)%last%blinfo(2)) THEN
            chainy(ibl)%first%prev=>chainy(ibl)%last
          ELSE IF (ibg>0) THEN
            IF (ibg/=chainy(ibl)%current%blinfo(2)) THEN
              ALLOCATE(chainy(ibl)%current%prev)
              chainy(ibl)%first=>chainy(ibl)%current%prev
              chainy(ibl)%first%next=>chainy(ibl)%current
              NULLIFY(chainy(ibl)%first%prev)
              chainy(ibl)%first%blinfo=
     $          NINT(seam(ibl)%vertex(1)%seam_out(1:5))
              chainy(ibl)%current=>chainy(ibl)%first
              chainy(ibl)%nlinks=chainy(ibl)%nlinks+1
              flag=1
            ENDIF
          ENDIF
        ENDDO
        IF (nprocs>1) THEN
          CALL mpi_allreduce(flag,itmp,1,mpi_nim_int,mpi_sum,
     $                       mpi_comm_world,ierror)
          flag=itmp
        ENDIF
        IF (flag==0) EXIT up
      ENDDO up
c-----------------------------------------------------------------------
c     check for periodic y-lines and check consistency of duplicate
c     chains on the same processor.
c-----------------------------------------------------------------------
      CALL parallel_line_perset(chainy)
      CALL parallel_line_check(chainy,'y')
c-----------------------------------------------------------------------
c-TMP write info after upward passes.
c-----------------------------------------------------------------------
c     CALL parallel_line_write(chainy,liney,'y','After up passes.')

c=======================================================================
c     phase II establishes line array dimensions.
c=======================================================================
c-----------------------------------------------------------------------
c     trasfer linked list information to 1D arrays.
c-----------------------------------------------------------------------
      ALLOCATE(linex(nrbl),liney(nrbl))
      DO ibl=1,nrbl
        linex(ibl)%par_dir=1
        liney(ibl)%par_dir=2
      ENDDO
      CALL parallel_line_tx(chainx,linex)
      CALL parallel_line_tx(chainy,liney)
c-----------------------------------------------------------------------
c     x-direction line limits and vertex counts.
c-----------------------------------------------------------------------
      CALL parallel_line_dims(linex)
      CALL parallel_line_count(linex)
c-----------------------------------------------------------------------
c-TMP write x-line info
c-----------------------------------------------------------------------
c     CALL parallel_line_write(chainx,linex,'x','After x-dims set.',
c    $                         wr_dims='yes')
c-----------------------------------------------------------------------
c     y-direction lines and vertex counts.
c-----------------------------------------------------------------------
      CALL parallel_line_dims(liney)
      CALL parallel_line_count(liney)
c-----------------------------------------------------------------------
c-TMP write y-line info
c-----------------------------------------------------------------------
c     CALL parallel_line_write(chainy,liney,'y','After y-dims set.',
c    $                         wr_dims='yes')

      RETURN
c=======================================================================
c     repeated operations are placed in internal subroutines.
c=======================================================================
      CONTAINS

c-----------------------------------------------------------------------
c       write existing line information.
c-----------------------------------------------------------------------
        SUBROUTINE parallel_line_write(chain,line,dir,message,wr_dims)

          TYPE(chain_type), DIMENSION(:), INTENT(IN) :: chain
          TYPE(line_type), DIMENSION(:), INTENT(IN) :: line
          CHARACTER(*), INTENT(IN) :: dir,message
          CHARACTER(3), OPTIONAL, INTENT(IN) :: wr_dims

          TYPE(chain_link_type), POINTER :: current
          INTEGER(i4) :: i
          CHARACTER(8) :: file_name

          WRITE(file_name,'(a,i1)') 'nodeout',node
          OPEN(UNIT=temp_unit,FILE=file_name,STATUS='UNKNOWN',
     $         POSITION='APPEND')

          WRITE(temp_unit,*) message
          DO ibl=1,nrbl
            IF (PRESENT(wr_dims)) THEN
              IF (wr_dims=='yes') THEN
                WRITE(temp_unit,*) 'ibl ',ibl,' perpst ',
     $            line(ibl)%perpst,' perpen ',line(ibl)%perpen,
     $            ' parast ',line(ibl)%parast,' paraen ',
     $            line(ibl)%paraen
                DO i=1,line(ibl)%nlinks
                  WRITE(temp_unit,*) 'ibl ',ibl,' ',dir,'link # ',i,
     $              ' nlinks ',line(ibl)%nlinks
                  WRITE(temp_unit,'(5i4)') line(ibl)%node(i),
     $              line(ibl)%glb_bl(i),line(ibl)%loc_bl(i),
     $              line(ibl)%mpara(i),line(ibl)%mperp(i)
                  WRITE(temp_unit,*) 'parast ',
     $              line(ibl)%bl2line_segst(i),
     $              ' paraen ',line(ibl)%bl2line_segen(i),' count ',
     $                line(ibl)%bl2line_count(i)
                ENDDO
              ENDIF
            ELSE
              i=1
              current=>chain(ibl)%first
              DO
                WRITE(temp_unit,*) 'ibl ',ibl,' ',dir,'link # ',i,
     $            ' nlinks ',chain(ibl)%nlinks
                WRITE(temp_unit,'(5i4)') current%blinfo
                IF (ASSOCIATED(current%next)) THEN
                  current=>current%next
                ELSE
                  EXIT
                ENDIF
                i=i+1
              ENDDO
            ENDIF
          ENDDO

          CLOSE(temp_unit)

        END SUBROUTINE parallel_line_write
c-----------------------------------------------------------------------
c       for periodic lines, reset the first pointer to the lowest global
c       block index.
c-----------------------------------------------------------------------
        SUBROUTINE parallel_line_perset(chain)

          TYPE(chain_type), DIMENSION(:), INTENT(INOUT) :: chain

          DO ibl=1,nrbl
            IF (ASSOCIATED(chain(ibl)%first%prev)) THEN
              itmp=HUGE(itmp)
              chain(ibl)%current=>chain(ibl)%first
              DO
                itmp=MIN(itmp,chain(ibl)%current%blinfo(2))
                IF (ASSOCIATED(chain(ibl)%current%next)) THEN
                  chain(ibl)%current=>chain(ibl)%current%next
                ELSE
                  EXIT
                ENDIF
              ENDDO
              chain(ibl)%last%next=>chain(ibl)%first
              DO
                IF (itmp==chain(ibl)%current%blinfo(2)) THEN
                  chain(ibl)%first=>chain(ibl)%current
                  chain(ibl)%last=>chain(ibl)%first%prev
                  NULLIFY(chain(ibl)%last%next)
                  chain(ibl)%current=>chain(ibl)%last
                ENDIF
                IF (ASSOCIATED(chain(ibl)%current%next)) THEN
                  chain(ibl)%current=>chain(ibl)%current%next
                ELSE
                  EXIT
                ENDIF
              ENDDO
            ENDIF
          ENDDO

        END SUBROUTINE parallel_line_perset
c-----------------------------------------------------------------------
c       eliminate duplicate chains on the same processor, but check
c       consistency as well.
c-----------------------------------------------------------------------
        SUBROUTINE parallel_line_check(chain,dir)

          TYPE(chain_type), DIMENSION(:), INTENT(IN) :: chain
          CHARACTER(*), INTENT(IN) :: dir

          DO ibl=1,nrbl
            DO ibg=ibl+1,nrbl
              IF (chain(ibg)%first%blinfo(2)==
     $            chain(ibl)%first%blinfo(2)) THEN
                IF (chain(ibg)%last%blinfo(2)/=
     $              chain(ibl)%last%blinfo(2).OR.
     $              chain(ibg)%nlinks/=chain(ibl)%nlinks) THEN
                  WRITE(msg,'(2(3a,i3),a)')
     $              'Parallel_line_check: line',dir,' ',ibl,
     $              ' and line',dir,' ',ibg,' are mis-matched.'
                  CALL nim_stop(msg)
                ENDIF
              ENDIF
            ENDDO
          ENDDO

        END SUBROUTINE parallel_line_check
c-----------------------------------------------------------------------
c       transfer linked list information to arrays indexed according
c       the the corresponding link index.  then, deallocate the linked
c       list, since it's no longer needed.
c-----------------------------------------------------------------------
        SUBROUTINE parallel_line_tx(chain,line)

          TYPE(chain_type), DIMENSION(:), INTENT(INOUT) :: chain
          TYPE(line_type), DIMENSION(:), INTENT(OUT) :: line

          INTEGER(i4) :: i,j

          DO ibl=1,nrbl
            line(ibl)%nlinks=chain(ibl)%nlinks
            ALLOCATE(line(ibl)%node(line(ibl)%nlinks))
            ALLOCATE(line(ibl)%glb_bl(line(ibl)%nlinks))
            ALLOCATE(line(ibl)%loc_bl(line(ibl)%nlinks))
            ALLOCATE(line(ibl)%mpara(line(ibl)%nlinks))
            ALLOCATE(line(ibl)%mperp(line(ibl)%nlinks))
            i=1
            chain(ibl)%current=>chain(ibl)%first
            IF (ASSOCIATED(chain(ibl)%first%prev)) THEN
              line(ibl)%periodic=.true.
            ELSE
              line(ibl)%periodic=.false.
            ENDIF
            line(ibl)%ncomm=0
            DO
              line(ibl)%node(i)=chain(ibl)%current%blinfo(1)
              IF (line(ibl)%node(i)/=node)
     $          line(ibl)%ncomm=line(ibl)%ncomm+1
              line(ibl)%glb_bl(i)=chain(ibl)%current%blinfo(2)
              line(ibl)%loc_bl(i)=chain(ibl)%current%blinfo(3)
              line(ibl)%mpara(i)=chain(ibl)%current%blinfo(4)
              line(ibl)%mperp(i)=chain(ibl)%current%blinfo(5)
              IF (ASSOCIATED(chain(ibl)%current%next)) THEN
                chain(ibl)%current=>chain(ibl)%current%next
                DEALLOCATE(chain(ibl)%current%prev)
              ELSE
                DEALLOCATE(chain(ibl)%current)
                EXIT
              ENDIF
              i=i+1
            ENDDO
            ALLOCATE(line(ibl)%recv_req(line(ibl)%ncomm))
            ALLOCATE(line(ibl)%send_req(line(ibl)%ncomm))
            ALLOCATE(line(ibl)%req_index(line(ibl)%nlinks))
            line(ibl)%req_index=HUGE(ibl)
            j=0
            DO i=1,line(ibl)%nlinks
              IF (line(ibl)%node(i)/=node) THEN
                j=j+1
                line(ibl)%req_index(i)=j
              ENDIF
            ENDDO
          ENDDO

        END SUBROUTINE parallel_line_tx
c-----------------------------------------------------------------------
c       determine the line array limits and segment limits contributed
c       from each link of a chain.
c-----------------------------------------------------------------------
        SUBROUTINE parallel_line_dims(line)

          TYPE(line_type), DIMENSION(:), INTENT(INOUT) :: line

          INTEGER(i4) :: i,nl,nl_bl,rem

          DO ibl=1,nrbl
            nl=line(ibl)%mperp(1)+1
            nl_bl=nl/line(ibl)%nlinks
            rem=MODULO(nl,line(ibl)%nlinks)
            ppen=-1
            past=0
            IF (line(ibl)%periodic) past=1
            line(ibl)%parast=past
            paen=line(ibl)%mpara(1)
            ALLOCATE(line(ibl)%bl2line_segst(line(ibl)%nlinks))
            ALLOCATE(line(ibl)%bl2line_segen(line(ibl)%nlinks))
            ALLOCATE(line(ibl)%line2bl_segst(line(ibl)%nlinks))
            ALLOCATE(line(ibl)%line2bl_segen(line(ibl)%nlinks))
            ALLOCATE(line(ibl)%bl_perpst(line(ibl)%nlinks))
            ALLOCATE(line(ibl)%bl_perpen(line(ibl)%nlinks))
            line(ibl)%glb_bl_recv=-1
            DO i=1,line(ibl)%nlinks
              ppst=ppen+1
              ppen=ppst+nl_bl-1
              IF (i<=rem) ppen=ppen+1
              line(ibl)%bl_perpst(i)=ppst
              line(ibl)%bl_perpen(i)=ppen
              IF (line(ibl)%glb_bl(i)==glb_bl(ibl)) THEN
                line(ibl)%perpst=ppst
                line(ibl)%perpen=ppen
              ENDIF
              line(ibl)%bl2line_segst(i)=past
              line(ibl)%bl2line_segen(i)=paen
              line(ibl)%line2bl_segen(i)=paen
              line(ibl)%line2bl_segst(i)=paen-line(ibl)%mpara(i)
              IF (line(ibl)%node(i)==node.AND.
     $            line(ibl)%loc_bl(i)==ibl) THEN
                line(ibl)%bl2line_paraen=line(ibl)%mpara(i)
                line(ibl)%bl2line_parast=line(ibl)%mpara(i)+past-paen
                line(ibl)%line2bl_paraen=line(ibl)%mpara(i)
                line(ibl)%line2bl_parast=0
              ENDIF
              IF (i<line(ibl)%nlinks) THEN
                past=paen+1
                paen=paen+line(ibl)%mpara(i+1)
              ENDIF
              IF (line(ibl)%node(i)==node.AND.
     $            line(ibl)%loc_bl(i)==ibl) THEN
                line(ibl)%glb_bl_recv=line(ibl)%glb_bl(i)
              ENDIF
            ENDDO
            line(ibl)%paraen=paen
            IF (line(ibl)%glb_bl_recv<=0)
     $        CALL nim_stop('Parallel_line_dims: error finding block.')
          ENDDO

        END SUBROUTINE parallel_line_dims
c-----------------------------------------------------------------------
c       for each line, count the number of vertices sent to/from each
c       link of the chain.
c-----------------------------------------------------------------------
        SUBROUTINE parallel_line_count(line)

          TYPE(line_type), DIMENSION(:), INTENT(INOUT) :: line

          INTEGER(i4) :: i

          DO ibl=1,nrbl
            ppst=line(ibl)%perpst
            ppen=line(ibl)%perpen
            ALLOCATE(line(ibl)%bl2line_count(line(ibl)%nlinks))
            ALLOCATE(line(ibl)%line2bl_count(line(ibl)%nlinks))
            DO i=1,line(ibl)%nlinks
              line(ibl)%bl2line_count(i)=(ppen-ppst+1)*
     $          (line(ibl)%bl2line_segen(i)
     $          -line(ibl)%bl2line_segst(i)+1)
              line(ibl)%line2bl_count(i)=(ppen-ppst+1)*
     $          (line(ibl)%line2bl_segen(i)
     $          -line(ibl)%line2bl_segst(i)+1)
            ENDDO
          ENDDO

        END SUBROUTINE parallel_line_count
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      END SUBROUTINE parallel_line_init


c-----------------------------------------------------------------------
c     subprogram 12. parallel_line_comm_bl2line.
c     communicate rblock-based data to global lines stretching across
c     the rblocks.  this is used for both parallel and serial
c     communication of this type.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE parallel_line_comm_bl2line(bl_data,line_data,nrbl,nbl,
     $                                      nq,line)
      USE local
      USE mpi_nim
      USE pardata
      USE vector_type_mod
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: nrbl,nbl,nq
      TYPE(vector_type), DIMENSION(nbl), INTENT(IN) :: bl_data
      TYPE(vector_type), DIMENSION(nrbl), INTENT(INOUT) :: line_data
      TYPE(line_type), DIMENSION(nrbl), INTENT(IN) :: line

      TYPE :: buffer_type
        TYPE(buffer_link_type), DIMENSION(:), POINTER :: link
      END TYPE buffer_type

      TYPE :: buffer_link_type
        REAL(r8), DIMENSION(:,:,:), POINTER :: buff
      END TYPE buffer_link_type

      INTEGER(i4) :: ibl,ibg,itmp,ix,iy,ip,ir,iq,ppst,ppen,past,paen
      INTEGER(i4), DIMENSION(1) :: loc
      INTEGER :: ierror
      INTEGER, DIMENSION(mpi_status_size) :: status
      TYPE(buffer_type), DIMENSION(nrbl) :: buff_bl
c-----------------------------------------------------------------------
c     post receives for off-proc communication.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        IF (line(ibl)%ncomm==0) CYCLE
        ppst=line(ibl)%perpst
        DO ip=1,line(ibl)%nlinks
          IF (line(ibl)%node(ip)/=node) THEN
            past=line(ibl)%bl2line_segst(ip)
            ir=line(ibl)%req_index(ip)
            CALL mpi_irecv(line_data(ibl)%arr(1,1,past),
     $                     nq*line(ibl)%bl2line_count(ip),
     $                     mpi_nim_real,line(ibl)%node(ip),
     $                     line(ibl)%glb_bl_recv,
     $                     mpi_comm_world,line(ibl)%recv_req(ir),
     $                     ierror)
          ENDIF
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     off-processor sends.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        IF (line(ibl)%ncomm==0) CYCLE
        ALLOCATE(buff_bl(ibl)%link(line(ibl)%nlinks))
        past=line(ibl)%bl2line_parast
        paen=line(ibl)%bl2line_paraen
        DO ip=1,line(ibl)%nlinks
          IF (line(ibl)%node(ip)==node) CYCLE
          ppst=line(ibl)%bl_perpst(ip)
          ppen=line(ibl)%bl_perpen(ip)
          ALLOCATE(buff_bl(ibl)%link(ip)%buff(ppst:ppen,nq,past:paen))
          IF (line(ibl)%par_dir==1) THEN
            DO ix=past,paen
              DO iq=1,nq
                buff_bl(ibl)%link(ip)%buff(:,iq,ix)=
     $            bl_data(ibl)%arr(iq,ix,ppst:ppen)
              ENDDO
            ENDDO
          ELSE
            DO iq=1,nq
              buff_bl(ibl)%link(ip)%buff(:,iq,:)=
     $          bl_data(ibl)%arr(iq,ppst:ppen,past:paen)
            ENDDO
          ENDIF
          itmp=SIZE(buff_bl(ibl)%link(ip)%buff)
          ir=line(ibl)%req_index(ip)
          CALL mpi_isend(buff_bl(ibl)%link(ip)%buff(ppst,1,past),itmp,
     $                   mpi_nim_real,line(ibl)%node(ip),
     $                   line(ibl)%glb_bl(ip),mpi_comm_world,
     $                   line(ibl)%send_req(ir),ierror)
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     on-processor copies.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        ppst=line(ibl)%perpst
        ppen=line(ibl)%perpen
        DO ip=1,line(ibl)%nlinks
          IF (line(ibl)%node(ip)==node) THEN
            past=line(ibl)%bl2line_segst(ip)
            paen=line(ibl)%bl2line_segen(ip)
            ibg=line(ibl)%loc_bl(ip)
            IF (line(ibl)%par_dir==1) THEN
              DO ix=past,paen
                DO iq=1,nq
                  line_data(ibl)%arr(:,iq,ix)=
     $             bl_data(ibg)%arr(iq,ix+line(ibg)%bl2line_parast-past,
     $                              ppst:ppen)
                ENDDO
              ENDDO
            ELSE
              DO iq=1,nq
                line_data(ibl)%arr(:,iq,past:paen)=
     $            bl_data(ibg)%arr(iq,ppst:ppen,
     $                line(ibg)%bl2line_parast:line(ibg)%bl2line_paraen)
              ENDDO
            ENDIF
          ENDIF
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     check that data has arrived.  block loop?
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        DO itmp=1,line(ibl)%ncomm
          CALL mpi_waitany(line(ibl)%ncomm,line(ibl)%recv_req,ir,
     $                     status,ierror)
          loc=MINLOC(ABS(line(ibl)%req_index-ir))
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     deallocate send buffers.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        IF (line(ibl)%ncomm==0) CYCLE
        DO itmp=1,line(ibl)%ncomm
          CALL mpi_waitany(line(ibl)%ncomm,line(ibl)%send_req,ir,
     $                     status,ierror)
          loc=MINLOC(ABS(line(ibl)%req_index-ir))
          ip=loc(1)
          DEALLOCATE(buff_bl(ibl)%link(ip)%buff)
        ENDDO
        DEALLOCATE(buff_bl(ibl)%link)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE parallel_line_comm_bl2line


c-----------------------------------------------------------------------
c     subprogram 13. parallel_line_comm_line2bl.
c     communicate rblock-based data from global lines stretching across
c     the rblocks.  this is used for both parallel and serial
c     communication of this type.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE parallel_line_comm_line2bl(bl_data,line_data,nrbl,nbl,
     $                                      nq,line)
      USE local
      USE mpi_nim
      USE pardata
      USE vector_type_mod
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: nrbl,nbl,nq
      TYPE(vector_type), DIMENSION(nbl), INTENT(INOUT) :: bl_data
      TYPE(vector_type), DIMENSION(nrbl), INTENT(IN) :: line_data
      TYPE(line_type), DIMENSION(nrbl), INTENT(IN) :: line

      TYPE :: buffer_type
        TYPE(buffer_link_type), DIMENSION(:), POINTER :: link
      END TYPE buffer_type

      TYPE :: buffer_link_type
        REAL(r8), DIMENSION(:,:,:), POINTER :: buff
      END TYPE buffer_link_type

      INTEGER(i4) :: ibl,ibg,itmp,ix,iy,ip,ir,iq,ppst,ppen,past,paen
      INTEGER(i4), DIMENSION(1) :: loc
      INTEGER :: ierror
      INTEGER, DIMENSION(mpi_status_size) :: status
      TYPE(buffer_type), DIMENSION(nrbl) :: buff_bl
c-----------------------------------------------------------------------
c     post receives for off-proc communication.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        IF (line(ibl)%ncomm==0) CYCLE
        ALLOCATE(buff_bl(ibl)%link(line(ibl)%nlinks))
        past=line(ibl)%line2bl_parast
        paen=line(ibl)%line2bl_paraen
        DO ip=1,line(ibl)%nlinks
          IF (line(ibl)%node(ip)/=node) THEN
            ppst=line(ibl)%bl_perpst(ip)
            ppen=line(ibl)%bl_perpen(ip)
            ALLOCATE(buff_bl(ibl)%link(ip)%buff(ppst:ppen,nq,past:paen))
            itmp=SIZE(buff_bl(ibl)%link(ip)%buff)
            ir=line(ibl)%req_index(ip)
            CALL mpi_irecv(buff_bl(ibl)%link(ip)%buff(ppst,1,past),itmp,
     $                     mpi_nim_real,line(ibl)%node(ip),
     $                     line(ibl)%glb_bl(ip),
     $                     mpi_comm_world,line(ibl)%recv_req(ir),
     $                     ierror)
          ENDIF
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     off-processor sends.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        IF (line(ibl)%ncomm==0) CYCLE
        ppst=line(ibl)%perpst
        DO ip=1,line(ibl)%nlinks
          IF (line(ibl)%node(ip)==node) CYCLE
          past=line(ibl)%line2bl_segst(ip)
          ir=line(ibl)%req_index(ip)
          CALL mpi_isend(line_data(ibl)%arr(1,1,past),
     $                   nq*line(ibl)%line2bl_count(ip),
     $                   mpi_nim_real,line(ibl)%node(ip),
     $                   line(ibl)%glb_bl_recv,mpi_comm_world,
     $                   line(ibl)%send_req(ir),ierror)
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     on-processor copies.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        ppst=line(ibl)%perpst
        ppen=line(ibl)%perpen
        DO ip=1,line(ibl)%nlinks
          IF (line(ibl)%node(ip)==node) THEN
            past=line(ibl)%line2bl_segst(ip)
            paen=line(ibl)%line2bl_segen(ip)
            ibg=line(ibl)%loc_bl(ip)
            IF (line(ibl)%par_dir==1) THEN
              DO ix=past,paen
                DO iq=1,nq
                  bl_data(ibg)%arr(iq,ix+line(ibg)%line2bl_parast-past,
     $                             ppst:ppen)=
     $              line_data(ibl)%arr(:,iq,ix)
                ENDDO
              ENDDO
            ELSE
              DO iq=1,nq
                bl_data(ibg)%arr(iq,ppst:ppen,
     $            line(ibg)%line2bl_parast:line(ibg)%line2bl_paraen)=
     $            line_data(ibl)%arr(:,iq,past:paen)
              ENDDO
            ENDIF
          ENDIF
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     check that data has arrived.  block loop?
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        IF (line(ibl)%ncomm==0) CYCLE
        past=line(ibl)%line2bl_parast
        paen=line(ibl)%line2bl_paraen
        DO itmp=1,line(ibl)%ncomm
          CALL mpi_waitany(line(ibl)%ncomm,line(ibl)%recv_req,ir,
     $                     status,ierror)
          loc=MINLOC(ABS(line(ibl)%req_index-ir))
          ip=loc(1)
          ppst=line(ibl)%bl_perpst(ip)
          ppen=line(ibl)%bl_perpen(ip)
          IF (line(ibl)%par_dir==1) THEN
            DO ix=past,paen
              DO iq=1,nq
                bl_data(ibl)%arr(iq,ix,ppst:ppen)=
     $            buff_bl(ibl)%link(ip)%buff(:,iq,ix)
              ENDDO
            ENDDO
          ELSE
            DO iq=1,nq
              bl_data(ibl)%arr(iq,ppst:ppen,past:paen)=
     $          buff_bl(ibl)%link(ip)%buff(:,iq,:)
            ENDDO
          ENDIF
          DEALLOCATE(buff_bl(ibl)%link(ip)%buff)
        ENDDO
        DEALLOCATE(buff_bl(ibl)%link)
      ENDDO
c-----------------------------------------------------------------------
c     check sends.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        DO itmp=1,line(ibl)%ncomm
          CALL mpi_waitany(line(ibl)%ncomm,line(ibl)%send_req,ir,
     $                     status,ierror)
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE parallel_line_comm_line2bl
c-----------------------------------------------------------------------
c     subprogram 14. parallel_line_comm_bl2line_comp.
c     communicate rblock-based data to global lines stretching across
c     the rblocks.  this is used for both parallel and serial
c     communication of this type.  this version communicates complex
c     data.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE parallel_line_comm_bl2line_comp(bl_data,line_data,nrbl,
     $                                           nbl,nq,line)
      USE local
      USE mpi_nim
      USE pardata
      USE vector_type_mod
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: nrbl,nbl,nq
      TYPE(cvector_2D_type), DIMENSION(nbl), INTENT(IN) :: bl_data
      TYPE(cvector_2D_type), DIMENSION(nrbl), INTENT(INOUT) :: line_data
      TYPE(line_type), DIMENSION(nrbl), INTENT(IN) :: line

      TYPE :: buffer_type
        TYPE(buffer_link_type), DIMENSION(:), POINTER :: link
      END TYPE buffer_type

      TYPE :: buffer_link_type
        COMPLEX(r8), DIMENSION(:,:,:), POINTER :: buff
      END TYPE buffer_link_type

      INTEGER(i4) :: ibl,ibg,itmp,ix,iy,ip,ir,iq,ppst,ppen,past,paen
      INTEGER(i4), DIMENSION(1) :: loc
      INTEGER :: ierror
      INTEGER, DIMENSION(mpi_status_size) :: status
      TYPE(buffer_type), DIMENSION(nrbl) :: buff_bl
c-----------------------------------------------------------------------
c     post receives for off-proc communication.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        IF (line(ibl)%ncomm==0) CYCLE
        ppst=line(ibl)%perpst
        DO ip=1,line(ibl)%nlinks
          IF (line(ibl)%node(ip)/=node) THEN
            past=line(ibl)%bl2line_segst(ip)
            ir=line(ibl)%req_index(ip)
            CALL mpi_irecv(line_data(ibl)%arr(1,1,past),
     $                     nq*line(ibl)%bl2line_count(ip),
     $                     mpi_nim_comp,line(ibl)%node(ip),
     $                     line(ibl)%glb_bl_recv,
     $                     mpi_comm_world,line(ibl)%recv_req(ir),
     $                     ierror)
          ENDIF
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     off-processor sends.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        IF (line(ibl)%ncomm==0) CYCLE
        ALLOCATE(buff_bl(ibl)%link(line(ibl)%nlinks))
        past=line(ibl)%bl2line_parast
        paen=line(ibl)%bl2line_paraen
        DO ip=1,line(ibl)%nlinks
          IF (line(ibl)%node(ip)==node) CYCLE
          ppst=line(ibl)%bl_perpst(ip)
          ppen=line(ibl)%bl_perpen(ip)
          ALLOCATE(buff_bl(ibl)%link(ip)%buff(ppst:ppen,nq,past:paen))
          IF (line(ibl)%par_dir==1) THEN
            DO ix=past,paen
              DO iq=1,nq
                buff_bl(ibl)%link(ip)%buff(:,iq,ix)=
     $            bl_data(ibl)%arr(iq,ix,ppst:ppen)
              ENDDO
            ENDDO
          ELSE
            DO iq=1,nq
              buff_bl(ibl)%link(ip)%buff(:,iq,:)=
     $          bl_data(ibl)%arr(iq,ppst:ppen,past:paen)
            ENDDO
          ENDIF
          itmp=SIZE(buff_bl(ibl)%link(ip)%buff)
          ir=line(ibl)%req_index(ip)
          CALL mpi_isend(buff_bl(ibl)%link(ip)%buff(ppst,1,past),itmp,
     $                   mpi_nim_comp,line(ibl)%node(ip),
     $                   line(ibl)%glb_bl(ip),mpi_comm_world,
     $                   line(ibl)%send_req(ir),ierror)
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     on-processor copies.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        ppst=line(ibl)%perpst
        ppen=line(ibl)%perpen
        DO ip=1,line(ibl)%nlinks
          IF (line(ibl)%node(ip)==node) THEN
            past=line(ibl)%bl2line_segst(ip)
            paen=line(ibl)%bl2line_segen(ip)
            ibg=line(ibl)%loc_bl(ip)
            IF (line(ibl)%par_dir==1) THEN
              DO iq=1,nq
                DO ix=past,paen
                  line_data(ibl)%arr(:,iq,ix)=
     $             bl_data(ibg)%arr(iq,ix+line(ibg)%bl2line_parast-past,
     $                              ppst:ppen)
                ENDDO
              ENDDO
            ELSE
              DO iq=1,nq
                line_data(ibl)%arr(:,iq,past:paen)=
     $            bl_data(ibg)%arr(iq,ppst:ppen,
     $                line(ibg)%bl2line_parast:line(ibg)%bl2line_paraen)
              ENDDO
            ENDIF
          ENDIF
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     check that data has arrived.  block loop?
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        DO itmp=1,line(ibl)%ncomm
          CALL mpi_waitany(line(ibl)%ncomm,line(ibl)%recv_req,ir,
     $                     status,ierror)
          loc=MINLOC(ABS(line(ibl)%req_index-ir))
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     deallocate send buffers.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        IF (line(ibl)%ncomm==0) CYCLE
        DO itmp=1,line(ibl)%ncomm
          CALL mpi_waitany(line(ibl)%ncomm,line(ibl)%send_req,ir,
     $                     status,ierror)
          loc=MINLOC(ABS(line(ibl)%req_index-ir))
          ip=loc(1)
          DEALLOCATE(buff_bl(ibl)%link(ip)%buff)
        ENDDO
        DEALLOCATE(buff_bl(ibl)%link)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE parallel_line_comm_bl2line_comp


c-----------------------------------------------------------------------
c     subprogram 15. parallel_line_comm_line2bl_comp.
c     communicate rblock-based data from global lines stretching across
c     the rblocks.  this is used for both parallel and serial
c     communication of this type.  this version communicates complex
c     data.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE parallel_line_comm_line2bl_comp(bl_data,line_data,nrbl,
     $                                           nbl,nq,line)
      USE local
      USE mpi_nim
      USE pardata
      USE vector_type_mod
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: nrbl,nbl,nq
      TYPE(cvector_2D_type), DIMENSION(nbl), INTENT(INOUT) :: bl_data
      TYPE(cvector_2D_type), DIMENSION(nrbl), INTENT(IN) :: line_data
      TYPE(line_type), DIMENSION(nrbl), INTENT(IN) :: line

      TYPE :: buffer_type
        TYPE(buffer_link_type), DIMENSION(:), POINTER :: link
      END TYPE buffer_type

      TYPE :: buffer_link_type
        COMPLEX(r8), DIMENSION(:,:,:), POINTER :: buff
      END TYPE buffer_link_type

      INTEGER(i4) :: ibl,ibg,itmp,ix,iy,ip,ir,iq,ppst,ppen,past,paen
      INTEGER(i4), DIMENSION(1) :: loc
      INTEGER :: ierror
      INTEGER, DIMENSION(mpi_status_size) :: status
      TYPE(buffer_type), DIMENSION(nrbl) :: buff_bl
c-----------------------------------------------------------------------
c     post receives for off-proc communication.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        IF (line(ibl)%ncomm==0) CYCLE
        ALLOCATE(buff_bl(ibl)%link(line(ibl)%nlinks))
        past=line(ibl)%line2bl_parast
        paen=line(ibl)%line2bl_paraen
        DO ip=1,line(ibl)%nlinks
          IF (line(ibl)%node(ip)/=node) THEN
            ppst=line(ibl)%bl_perpst(ip)
            ppen=line(ibl)%bl_perpen(ip)
            ALLOCATE(buff_bl(ibl)%link(ip)%buff(ppst:ppen,nq,past:paen))
            itmp=SIZE(buff_bl(ibl)%link(ip)%buff)
            ir=line(ibl)%req_index(ip)
            CALL mpi_irecv(buff_bl(ibl)%link(ip)%buff(ppst,1,past),itmp,
     $                     mpi_nim_comp,line(ibl)%node(ip),
     $                     line(ibl)%glb_bl(ip),
     $                     mpi_comm_world,line(ibl)%recv_req(ir),
     $                     ierror)
          ENDIF
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     off-processor sends.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        IF (line(ibl)%ncomm==0) CYCLE
        ppst=line(ibl)%perpst
        DO ip=1,line(ibl)%nlinks
          IF (line(ibl)%node(ip)==node) CYCLE
          past=line(ibl)%line2bl_segst(ip)
          ir=line(ibl)%req_index(ip)
          CALL mpi_isend(line_data(ibl)%arr(1,1,past),
     $                   nq*line(ibl)%line2bl_count(ip),
     $                   mpi_nim_comp,line(ibl)%node(ip),
     $                   line(ibl)%glb_bl_recv,mpi_comm_world,
     $                   line(ibl)%send_req(ir),ierror)
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     on-processor copies.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        ppst=line(ibl)%perpst
        ppen=line(ibl)%perpen
        DO ip=1,line(ibl)%nlinks
          IF (line(ibl)%node(ip)==node) THEN
            past=line(ibl)%line2bl_segst(ip)
            paen=line(ibl)%line2bl_segen(ip)
            ibg=line(ibl)%loc_bl(ip)
            IF (line(ibl)%par_dir==1) THEN
              DO ix=past,paen
                DO iq=1,nq
                  bl_data(ibg)%arr(iq,ix+line(ibg)%line2bl_parast-past,
     $                             ppst:ppen)=
     $              line_data(ibl)%arr(:,iq,ix)
                ENDDO
              ENDDO
            ELSE
              DO iq=1,nq
                bl_data(ibg)%arr(iq,ppst:ppen,
     $            line(ibg)%line2bl_parast:line(ibg)%line2bl_paraen)=
     $            line_data(ibl)%arr(:,iq,past:paen)
              ENDDO
            ENDIF
          ENDIF
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     check that data has arrived.  block loop?
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        IF (line(ibl)%ncomm==0) CYCLE
        past=line(ibl)%line2bl_parast
        paen=line(ibl)%line2bl_paraen
        DO itmp=1,line(ibl)%ncomm
          CALL mpi_waitany(line(ibl)%ncomm,line(ibl)%recv_req,ir,
     $                     status,ierror)
          loc=MINLOC(ABS(line(ibl)%req_index-ir))
          ip=loc(1)
          ppst=line(ibl)%bl_perpst(ip)
          ppen=line(ibl)%bl_perpen(ip)
          IF (line(ibl)%par_dir==1) THEN
            DO ix=past,paen
              DO iq=1,nq
                bl_data(ibl)%arr(iq,ix,ppst:ppen)=
     $            buff_bl(ibl)%link(ip)%buff(:,iq,ix)
              ENDDO
            ENDDO
          ELSE
            DO iq=1,nq
              bl_data(ibl)%arr(iq,ppst:ppen,past:paen)=
     $          buff_bl(ibl)%link(ip)%buff(:,iq,:)
            ENDDO
          ENDIF
          DEALLOCATE(buff_bl(ibl)%link(ip)%buff)
        ENDDO
        DEALLOCATE(buff_bl(ibl)%link)
      ENDDO
c-----------------------------------------------------------------------
c     check sends.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        DO itmp=1,line(ibl)%ncomm
          CALL mpi_waitany(line(ibl)%ncomm,line(ibl)%send_req,ir,
     $                     status,ierror)
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE parallel_line_comm_line2bl_comp
c-----------------------------------------------------------------------
c     subprogram 16. qp0_bcast
c     copy and broadcast the symmetric part of a complex array, that
c     has the dimensions for quadrature-point data, to all layers.  The
c     calling routine needs an interface block to pass the pointer.
c-----------------------------------------------------------------------
      SUBROUTINE qp0_bcast(threedata,twodata,nl)
      USE local
      USE pardata
      USE mpi_nim
      IMPLICIT NONE

      COMPLEX(r8), DIMENSION(:,:,:,:) :: threedata
      REAL(r8), DIMENSION(:,:,:) :: twodata
      INTEGER(i4), INTENT(IN) :: nl

      INTEGER :: ndata,ierror
c-----------------------------------------------------------------------
c     first copy the symmetric part of the 3D data then broadcast.
c-----------------------------------------------------------------------
      IF (ilayer==0) twodata=threedata(:,:,:,1)
      IF (nl>1) THEN
        ndata=SIZE(twodata)
        CALL mpi_bcast(twodata,ndata,mpi_nim_real,0,
     $                 comm_mode,ierror)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE qp0_bcast
c-----------------------------------------------------------------------
c     subprogram 17. qpc_bcast
c     copy and broadcast a specified range of a complex array, that
c     has the dimensions for quadrature-point data, to all layers.  The
c     calling routine needs an interface block to pass the pointer.
c
c     the input includes pointers to the complex data, where threedata
c     holds the Fourier components for this layer, and bcastdata is
c     returned with just the istcomp:iencomp components out of the
c     global index, 1:ncomp, where nmodes_total is passed into ncomp.
c-----------------------------------------------------------------------
      SUBROUTINE qpc_bcast(threedata,bcastdata,nl,ncomp,istcomp,iencomp)
      USE local
      USE pardata
      USE mpi_nim
      IMPLICIT NONE

      COMPLEX(r8), DIMENSION(:,:,:,:) :: threedata,bcastdata
      INTEGER(i4), INTENT(IN) :: nl,ncomp,istcomp,iencomp

      INTEGER :: ndata,ierror
      INTEGER(i4) :: il,llo,lhi,nm,rem,ibc,nbc,i0,i1,i,n1,n2,n3
      COMPLEX(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: ctmp
c-----------------------------------------------------------------------
c     the communication is send-ordered, so loop over layers, and
c     assess how much of the specified range is in each layer.
c-----------------------------------------------------------------------
      nm=ncomp/nl
      rem=ncomp-nm*nl
      llo=1
      DO il=0,nl-1
        IF (il<rem) THEN
          lhi=llo+nm
        ELSE
          lhi=llo+nm-1
        ENDIF
        ibc=MAX(istcomp,llo)
        nbc=MIN(iencomp,lhi)-ibc+1
        i0=ibc-istcomp+1
        i1=i0+nbc-1
        IF (nbc>0) THEN
          IF (nl==1) THEN
            bcastdata(:,:,:,i0:i1)=threedata(:,:,:,ibc:ibc+nbc-1)
          ELSE
            n1=SIZE(bcastdata,1)
            n2=SIZE(bcastdata,2)
            n3=SIZE(bcastdata,3)
            ALLOCATE(ctmp(n1,n2,n3,nbc))
            IF (ilayer==il) THEN
              DO i=1,nbc
                ctmp(:,:,:,i)=threedata(:,:,:,ibc+i-llo)
              ENDDO
            ENDIF
            ndata=n1*n2*n3*nbc
            CALL mpi_bcast(ctmp,ndata,mpi_nim_comp,il,
     $                     comm_mode,ierror)
            DO i=1,nbc
              bcastdata(:,:,:,i0+i-1)=ctmp(:,:,:,i)
            ENDDO
            DEALLOCATE(ctmp)
          ENDIF
        ENDIF
        llo=lhi+1
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE qpc_bcast
c-----------------------------------------------------------------------
c     subprogram 18. parallel_block_dealloc
c     deallocate the arrays that are used to hold block decomposition
c     information.
c-----------------------------------------------------------------------

      SUBROUTINE parallel_block_dealloc
      USE pardata
      IMPLICIT NONE

      DEALLOCATE(mode2layer)
      DEALLOCATE(block2proc)
      DEALLOCATE(global2local)
      DEALLOCATE(loc2glob)
      DEALLOCATE(layer2proc)
      DEALLOCATE(block_sizes)

      RETURN
      END SUBROUTINE parallel_block_dealloc

c-----------------------------------------------------------------------
c     subprogram 19. parallel_seam_dealloc
c     deallocate seam communication data stuctures
c-----------------------------------------------------------------------

      SUBROUTINE parallel_seam_dealloc
      USE fields
      USE seam_storage_mod
      USE mpi_nim
      USE pardata
      IMPLICIT NONE

      INTEGER(i4) :: isend,irecv

      DO isend=1,nsend
        DEALLOCATE(send(isend)%block)
        DEALLOCATE(send(isend)%vertex)
        DEALLOCATE(send(isend)%image)
        DEALLOCATE(send(isend)%data)
        DEALLOCATE(send(isend)%cdata)
      ENDDO

      DO irecv=1,nrecv
        DEALLOCATE(recv(irecv)%block)
        DEALLOCATE(recv(irecv)%vertex)
        DEALLOCATE(recv(irecv)%order)
        DEALLOCATE(recv(irecv)%data)
        DEALLOCATE(recv(irecv)%cdata)
      ENDDO

      IF (nrecv > 0) DEALLOCATE(recv_request)
      IF (nsend > 0) DEALLOCATE(send)
      IF (nrecv > 0) DEALLOCATE(recv)
      IF (nself > 0) DEALLOCATE(self)

      RETURN
      END SUBROUTINE parallel_seam_dealloc

c-----------------------------------------------------------------------
c     subprogram 20. parallel_seg_dealloc
c     deallocate segment communication data stuctures
c-----------------------------------------------------------------------

      SUBROUTINE parallel_seg_dealloc
      USE fields
      USE seam_storage_mod
      USE mpi_nim
      USE pardata
      IMPLICIT NONE

      INTEGER(i4) :: isend,irecv

      DO isend=1,nsendseg
        DEALLOCATE(sendseg(isend)%block)
        DEALLOCATE(sendseg(isend)%segment)
        DEALLOCATE(sendseg(isend)%data)
        DEALLOCATE(sendseg(isend)%cdata)
        DEALLOCATE(sendseg(isend)%mat_data)
        DEALLOCATE(sendseg(isend)%cmat_data)
      ENDDO

      DO irecv=1,nrecvseg
        DEALLOCATE(recvseg(irecv)%block)
        DEALLOCATE(recvseg(irecv)%segment)
        DEALLOCATE(recvseg(irecv)%data)
        DEALLOCATE(recvseg(irecv)%cdata)
        DEALLOCATE(recvseg(irecv)%mat_data)
        DEALLOCATE(recvseg(irecv)%cmat_data)
      ENDDO

      IF (nrecvseg > 0) DEALLOCATE(recvseg_request)
      if (nsendseg > 0) DEALLOCATE(sendseg)
      if (nrecvseg > 0) DEALLOCATE(recvseg)

      RETURN
      END SUBROUTINE parallel_seg_dealloc

