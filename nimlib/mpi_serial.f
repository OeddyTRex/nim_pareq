c this file is to be used on serial machines.  it contains two parts.
c the first is a module of parameters and the second is a set of
c dummy routines to replace mpi routines.

c module for interface to F77 MPI include file
c defines all MPI variables via parameter statements
c use this module to define machine-specific MPI datatypes
c   so that NIMROD source will not have to change for a new
c   machine or if system.f SELECTED_KINDS are changed

      MODULE mpi_nim

      USE local

      INTEGER(i4), PARAMETER :: mpi_nim_int = 1_i4
      INTEGER(i4), PARAMETER :: mpi_nim_real = 1_i4
      INTEGER(i4), PARAMETER :: mpi_nim_comp = 1_i4
      INTEGER(i4), PARAMETER :: mpi_nim_logical = 1_i4
      INTEGER(i4), PARAMETER :: mpi_nim_char = 1_i4

      INTEGER(i4), PARAMETER :: mpi_comm_world = 1_i4
      INTEGER(i4), PARAMETER :: mpi_max = 1_i4
      INTEGER(i4), PARAMETER :: mpi_min = 1_i4
      INTEGER(i4), PARAMETER :: mpi_sum = 1_i4
      INTEGER(i4), PARAMETER :: mpi_land = 1_i4
      INTEGER(i4), PARAMETER :: mpi_lor = 1_i4
      INTEGER(i4), PARAMETER :: mpi_source = 1_i4
      INTEGER(i4), PARAMETER :: mpi_any_source = 1_i4
      INTEGER(i4), PARAMETER :: mpi_status_size = 1_i4
      INTEGER(i4), PARAMETER :: mpi_request_null = 1_i4

      END MODULE mpi_nim


c stubs for MPI calls - use on serial machine to replace real MPI library
c  so NIMROD will still compile
c except for mpi_comm_rank, mpi_comm_size and mpi_allreduce
c  these are all no-operation rouines


      subroutine mpi_init(ierror)
      integer ierror

      return
      end


      subroutine mpi_finalize(ierror)
      integer ierror

      return
      end


c return processor id = 0

      subroutine mpi_comm_rank(comm,rank,ierror)
      integer comm,rank,ierror

      rank = 0

      return
      end

c return # of processors = 1

      subroutine mpi_comm_size(comm,size,ierror)
      integer comm,size,ierror

      size = 1

      return
      end


      subroutine mpi_bcast(buffer,count,datatype,root,comm,ierror)
      integer buffer,count,datatype,root,comm,ierror

      return
      end


      subroutine mpi_barrier(comm,ierror)
      integer comm,ierror

      return
      end


      subroutine mpi_allreduce(sendbuf,recvbuf,count,
     $     datatype,op,comm,ierror)
      integer sendbuf,recvbuf,count,datatype,op,comm,ierror
      dimension sendbuf(count),recvbuf(count)

      integer i

      do 10 i=1,count
        recvbuf(i)=sendbuf(i)
   10 continue

      return
      end


      subroutine mpi_reduce_scatter(sendbuf,recvbuf,recvcounts,
     $     datatype,op,comm,ierror)
      integer sendbuf,recvbuf,recvcounts,datatype,op,comm,ierror

      return
      end


      subroutine mpi_send(buf,count,datatype,dest,tag,comm,ierror)
      integer buf,count,datatype,dest,tag,comm,ierror

      return
      end


      subroutine mpi_recv(buf,count,datatype,source,tag,comm,
     $     status,ierror)
      integer buf,count,datatype,source,tag,comm,status,ierror

      return
      end


      subroutine mpi_irecv(buf,count,datatype,source,tag,comm,
     $     request,ierror)
      integer buf,count,datatype,source,tag,comm,request,ierror

      return
      end


      subroutine mpi_isend(buf,count,datatype,source,tag,comm,
     $     request,ierror)
      integer buf,count,datatype,source,tag,comm,request,ierror

      return
      end


      subroutine mpi_test(request,flag,status,ierror)
      logical flag
      integer request,status,ierror

      return
      end


      subroutine mpi_comm_split(datatype1,ilayer,ii,datatype2,ierror)
      integer datatype1,datatype2,ilayer,ierror,ii

      return
      end


      subroutine mpi_wait(recv_request,status,ierror)
      integer nrecv,ierror
      dimension recv_reqest(1),statuses(1)

      return
      end


      subroutine mpi_waitall(nrecv,recv_request,statuses,ierror)
      integer nrecv,ierror
      dimension recv_reqest(nrecv),statuses(1,1)

      return
      end


      subroutine mpi_waitany(nrecv,recv_request,irecv,status,ierror)
      integer nrecv,irecv,ierror
      dimension recv_reqest(nrecv),status(1)

      return
      end



      subroutine mpi_testany(nrecv,recv_request,irecv,flag,status,
     $                       ierror)
      integer nrecv,irecv,ierror
      dimension recv_reqest(nrecv),status(1)
      logical flag

      return
      end



      subroutine mpi_gather(sendbuf,counts,datatypes,
     $        recvbuf,countr,datatyper,irecv,comm,ierror)
      integer sendbuf,recvbuf,counts,countr,datatypes,datatyper,
     $        comm,irecv,ierror
      dimension sendbuf(counts),recvbuf(counts)

      integer i

      do 10 i=1,counts
        recvbuf(i)=sendbuf(i)
   10 continue

      return
      end



      subroutine mpi_allgather(sendbuf,counts,datatypes,
     $        recvbuf,countr,datatyper,comm,ierror)
      integer sendbuf,recvbuf,counts,countr,datatypes,datatyper,
     $        comm,ierror
      dimension sendbuf(counts),recvbuf(counts)

      integer i

      do 10 i=1,counts
        recvbuf(i)=sendbuf(i)
   10 continue

      return
      end



      subroutine mpi_allgatherv(sendbuf,counts,datatypes,
     $        recvbuf,countr,displs,datatyper,comm,ierror)
      integer sendbuf,recvbuf,counts,countr,datatypes,datatyper,
     $        displs,comm,ierror
      dimension sendbuf(counts),recvbuf(counts)

      integer i

      do 10 i=1,counts
        recvbuf(i)=sendbuf(i)
   10 continue

      return
      end



      subroutine mpi_alltoallv(sendbuf,scounts,sdispls,datatypes,
     $        recvbuf,rcounts,rdispls,datatyper,comm,ierror)
      integer sendbuf,scounts,sdispls,recvbuf,rcounts,rdispls,
     $        datatypes,datatyper,comm,ierror
      dimension sendbuf(scounts),recvbuf(rcounts)

      integer i

      do 10 i=1,min(rcounts,scounts)
        recvbuf(i)=sendbuf(i)
   10 continue

      return
      end



      subroutine mpi_abort(comm,code,ierror)
      integer comm,code,ierror

      return
      end



      subroutine comm_create_(nsend,procsend,comm,nrecv,plan)
      USE local
      integer nsend,nrecv,procsend
      real(r8) plan

      nrecv=nsend
      plan=REAL(nsend)

      return
      end



      subroutine comm_destroy_(plan)
      USE local
      real(r8) plan

      return
      end



      subroutine comm_do_(plan,sendbuf,n,recvbuf)
      USE local
      integer n,sendbuf(*),recvbuf(*)
      real(r8) plan

      integer i

      do 10 i=1,nint(plan)
        recvbuf(i)=sendbuf(i)
   10 continue

      return
      end
