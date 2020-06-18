c-----------------------------------------------------------------------
c     file plot_data.f
c     contains dump file and data management for nimplot.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0.  plot_data.
c     1.  acquire_data.
c     2.  fields_at_angle.
c     3.  flux_comps.
c     4.  cart_comps.
c     5.  add_eq_fields.
c     6.  eval_at_angle.
c-----------------------------------------------------------------------
c     subprogram 0. plot_data.
c     module declarations.
c-----------------------------------------------------------------------
      MODULE plot_data
      USE local
      USE fields
      USE edge
      USE global
      IMPLICIT NONE

      REAL(r8) :: phi,dphi
      INTEGER(i4) :: mode_plot_start,mode_plot_end,iphi
      CHARACTER(1) :: flag,add_eq,do_flux='n',do_jfromb='y'
      CHARACTER(16) :: selection_type
      CHARACTER(128) :: last_dump='initial_no_dump'
      LOGICAL :: read_file,inquire,inquire_flux

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. acquire_data.
c     acquires data from a dump file or from an existing data structure.
c-----------------------------------------------------------------------
      SUBROUTINE acquire_data(dump_data,dump_step,dump_time)
      USE input
      USE pardata
      USE rblock
      USE tblock
      USE dump

      CHARACTER(*), INTENT(IN) :: dump_data
      INTEGER(i4), INTENT(OUT) :: dump_step
      REAL(r8), INTENT(OUT) :: dump_time

      INTEGER(i4) :: ibl,imode,max_nqty,offm
      INTEGER(i4), SAVE :: im=0
      LOGICAL, SAVE :: paralloc=.false.,line_init=.false.
      LOGICAL :: new_dump
c-----------------------------------------------------------------------
c     read dump file if needed.  deallocate previous data if it exists.
c     option c corrupts the data.
c-----------------------------------------------------------------------
      dump_file=dump_data
      new_dump=.false.
      IF (read_file.AND.dump_file/=last_dump) THEN
        new_dump=.true.
        IF (last_dump/='initial_no_dump') CALL deallocate_data
        CALL dump_read(nmodes,keff,dump_time,dump_step)
        last_dump=dump_file
      ENDIF
c-----------------------------------------------------------------------
c     compute # of cells in the polar region.
c-----------------------------------------------------------------------
      mx=0
      my=0
      DO ibl=nybl,nrbl,nybl
        mx=mx+rb(ibl)%mx
      ENDDO
      DO ibl=1,nybl
        my=my+rb(ibl)%my
      ENDDO
c-----------------------------------------------------------------------
c     set parallel variables to serial values.
c-----------------------------------------------------------------------
      nbl_total=nbl
      nrbl_total=nrbl
      IF (.NOT.paralloc) THEN
        ALLOCATE(block2proc(nbl))
        ALLOCATE(global2local(nbl))
        ALLOCATE(layer2proc(nbl))
        ALLOCATE(block_sizes(1,nbl))
        ALLOCATE(loc2glob(nbl))
        block2proc=0
        layer2proc=0
        DO ibl=1,nbl
          global2local(ibl)=ibl
          loc2glob(ibl)=ibl
          block_sizes(1,ibl)=seam(ibl)%nvert
        ENDDO
c-----------------------------------------------------------------------
c       determine the total poloidal dimension for pseudo-spectral
c       coding and ffts.
c-----------------------------------------------------------------------
        ALLOCATE(mps_block(nbl))
        DO ibl=1,nbl
          IF (ibl<=nrbl) THEN
            mps_block(ibl)=rb(ibl)%mx*rb(ibl)%my
          ELSE
            mps_block(ibl)=tb(ibl)%mcell
          ENDIF
        ENDDO
        nprocs_layer=nprocs
        node_layer=0
        ilayer=0
        nmodes_total=nmodes
        mode_lo=1
        mode_hi=nmodes
        paralloc=.true.
      ENDIF
c-----------------------------------------------------------------------
c     inquire for configuration space angle or Fourier component
c     selection.
c-----------------------------------------------------------------------
      not_all: IF (flag/='a') THEN
        IF (inquire) THEN
          add_eq='n'
          IF (flag=='c') THEN
            WRITE(nim_wr,20)
     $       'Please enter the desired toroidal angle (in radians/2pi',
     $       'or z/per_length for linear geometry).'
            WRITE(nim_wr,11)
            READ(nim_rd,*) phi
            phi=MODULO(phi,1._r8)
            IF (phi<0) phi=phi+1
            IF (selection_type=='data') THEN
              WRITE(nim_wr,100)
     $         'Add the equilibrium fields to the perturbed fields, ',
     $         'and zero the equilibrium arrays for these plots? ',
     $         '(',"'",'y',"'",' or ',"'",'n',"'",')'
              WRITE(nim_wr,11)
              READ(nim_rd,*) add_eq
            ENDIF
          ELSE IF (flag=='p'.OR.flag=='t') THEN
            WRITE(nim_wr,20)
     $       'Please enter the desired periodic coorindate increment',
     $       '(in radians/2pi or z/per_length for linear geometry).'
            WRITE(nim_wr,11)
            READ(nim_rd,*) dphi
            IF (phi<0) dphi=dphi+1
            IF (selection_type=='data') THEN
              WRITE(nim_wr,100)
     $         'Add the equilibrium fields to the perturbed fields, ',
     $         'and zero the equilibrium arrays for these plots? ',
     $         '(',"'",'y',"'",' or ',"'",'n',"'",')'
              WRITE(nim_wr,11)
              READ(nim_rd,*) add_eq
            ENDIF
          ELSE
            WRITE(nim_wr,20)
     $      'The following is a list of mode indices and wavenumbers',
     $      '(n for toroidal geometry or k=2*pi*n for linear geometry).'
            DO im=1,nmodes
              WRITE(nim_wr,'(a,i3,a,i3,a,es10.4)') '>>>    ',im-1,
     $         ' => k(',im-1,') = ',keff(im)
            ENDDO
            WRITE(nim_wr,10) 'Select the index for the desired mode.'
            WRITE(nim_wr,11)
            READ(nim_rd,*) im
            IF (selection_type=='data'.AND.keff(im+1)==0) THEN
              WRITE(nim_wr,100)
     $         'Add the equilibrium fields to the perturbed fields, ',
     $         'and zero the equilibrium arrays for these plots? ',
     $         '(',"'",'y',"'",' or ',"'",'n',"'",')'
              WRITE(nim_wr,11)
              READ(nim_rd,*) add_eq
            ENDIF
          ENDIF
        ENDIF
      ELSE not_all
        IF (selection_type=='data'.AND.inquire) THEN
          DO imode=1,nmodes
            IF (keff(imode)==0) THEN
              WRITE(nim_wr,100)
     $         'Add the equilibrium fields to the perturbed fields, ',
     $         'and zero the equilibrium arrays for these plots? ',
     $         '(',"'",'y',"'",' or ',"'",'n',"'",')'
              WRITE(nim_wr,11)
              READ(nim_rd,*) add_eq
            ENDIF
          ENDDO
        ENDIF
      ENDIF not_all
c-----------------------------------------------------------------------
c     initialize communication and block data, and
c     compute current density as a vertex quantity.
c-----------------------------------------------------------------------
      IF (new_dump) THEN
        DO ibl=1,nrbl
          CALL rblock_set(ngr,poly_degree,integration_formula,rb(ibl))
          CALL rblock_basis_set(rb(ibl),(/poly_degree/),
     $                          (/-1_i4/),(/-1_i4/),(/-1_i4/),
     $                          met_spl,geom)
        ENDDO
        DO ibl=nrbl+1,nbl
          CALL tblock_set(tb(ibl))
          CALL tri_linear_get_areas(tb(ibl)%tgeom)
          CALL tblock_basis_set(tb(ibl),geom)
        ENDDO
        CALL variable_alloc(nmodes)
        CALL pointer_init
        max_nqty=MAX(9_i4,6_i4*nmodes)
        CALL edge_init(max_nqty,9_i4)
        CALL boundary_init(geom)
        max_nqty=MAX(2_i4,MAX(9_i4,6_i4*nmodes)*(poly_degree-1))
        offm=poly_degree**2+poly_degree
        CALL edge_segment_init(max_nqty,3_i4,offm)
        IF (solver=='gl_diaga'.AND..NOT.line_init) THEN
          CALL parallel_line_init(rb,nrbl,nbl,poly_degree)
          line_init=.true.
        ENDIF
        CALL block_create_tang(poly_degree)
        CALL mass_mat_init
      ENDIF
c-----------------------------------------------------------------------
c     find nodal current density if needed.
c-----------------------------------------------------------------------
      IF (selection_type=='data'.AND.do_jfromb=='y') THEN
        IF (flag=='o') THEN
          mode_plot_start=im+1
          mode_plot_end=im+1
          CALL jfromb(mode_plot_start,mode_plot_end)
        ELSE
          mode_plot_start=1
          mode_plot_end=nmodes
          CALL jfromb(1_i4,nmodes)
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     add the equilibrium field to the n=0 part of the solution if
c     appropriate and requested.
c-----------------------------------------------------------------------
      IF (selection_type=='data'.AND.add_eq=='y') THEN
        CALL add_eq_fields(geom,poly_degree)
        last_dump='xxx_bogus_xxx'
      ENDIF
c-----------------------------------------------------------------------
c     select a single mode, all, or do an FFT.
c-----------------------------------------------------------------------
      IF (flag=='c'.AND.selection_type=='data') THEN
        CALL fields_at_angle
        last_dump='xxx_bogus_xxx'
        mode_plot_start=1
        mode_plot_end=1
      ELSE IF ((flag=='p'.OR.flag=='t').AND.selection_type=='data') THEN
        phi=MIN(dphi*iphi,1._r8)
        phi=MODULO(phi,1._r8)
        CALL fields_at_angle
        IF (flag=='t') CALL cart_comps
        last_dump='xxx_bogus_xxx'
        mode_plot_start=1
        mode_plot_end=1
      ELSE IF (flag=='o') THEN
        mode_plot_start=im+1
        mode_plot_end=im+1
      ELSE
        mode_plot_start=1
        mode_plot_end=nmodes
      ENDIF
c-----------------------------------------------------------------------
c     flux surface component option.  this corrupts the data, so don't
c     allow it to be reused.
c-----------------------------------------------------------------------
      IF (selection_type=='data'.AND.flag/='t') THEN
        IF (inquire_flux) THEN
          WRITE(nim_wr,110)
     $     'Convert poloidal field components to flux surface normal',
     $     'and tangential components? ',
     $     '(',"'",'y',"'",' or ',"'",'n',"'",')',
     $     'This assumes that the grid is aligned with the flux',
     $     'surfaces, and the transformation is done in rblocks only!!!'
          READ(nim_rd,*) do_flux
        ENDIF
        IF (do_flux=='y') THEN
          CALL flux_comps(poly_degree)
          last_dump='xxx_bogus_xxx'
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     formats:
c-----------------------------------------------------------------------
  10  format(/'>>> ',a)
  20  format(2(/'>>> ',a))
  11  format('>>>? ',$)
 100  format(2(/'>>> ',a),/'>>> ',10a)
 110  format(/'>>> ',a,/'>>> ',10a,2(/'>>> ',a),/'>>>? ',$)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE acquire_data
c-----------------------------------------------------------------------
c     subprogram 2. fields_at_angle.
c     set the first mode data to the values at the specified angle in
c     configuration space.
c-----------------------------------------------------------------------
      SUBROUTINE fields_at_angle

c-----------------------------------------------------------------------
c     compute each dependent variable at the appropriate angle.
c-----------------------------------------------------------------------
      CALL eval_at_angle(be,phi)
      CALL eval_at_angle(ja,phi)
      CALL eval_at_angle(ve,phi)
      CALL eval_at_angle(pres,phi)
      CALL eval_at_angle(prese,phi)
      CALL eval_at_angle(nd,phi)
      CALL eval_at_angle(conc,phi)
      CALL eval_at_angle(tele,phi)
      CALL eval_at_angle(tion,phi)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE fields_at_angle
c-----------------------------------------------------------------------
c     subprogram 3. flux_comps.
c     convert b, j and v to flux surface components.  this is based
c     on analyze routines in files get_data.f from Alan Glasser.
c-----------------------------------------------------------------------
      SUBROUTINE flux_comps(poly_degree)
      
      INTEGER(i4), INTENT(IN) :: poly_degree

      INTEGER(i4) :: ib,im,nm,mxb,myb,ix,iy,iq,ibasis,ix0,iy0,ibn
      REAL(r8) :: dx,dy,rzmax,yloc
      REAL(r8), DIMENSION(2) :: norm,tang,temp
      COMPLEX(r8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: jatmp,betmp,
     $             vetmp
c-----------------------------------------------------------------------
c     loop over rblocks only.  the double evaluation of lagrrz allows
c     averaging the potentially discontinous y-derivative.  the usual
c     conformance of block dimensions is assumed.
c-----------------------------------------------------------------------
      DO ib=1,nrbl
        nm=SIZE(rb(ib)%ve%fs,4)
        mxb=rb(ib)%mx
        myb=rb(ib)%my
        rzmax=MAXVAL(ABS(rb(ib)%rz%fs))
        ALLOCATE(jatmp(3,0:mxb,0:myb,nm,poly_degree**2))
        ALLOCATE(betmp(3,0:mxb,0:myb,nm,poly_degree**2))
        ALLOCATE(vetmp(3,0:mxb,0:myb,nm,poly_degree**2))
        DO ibasis=1,poly_degree**2
          ix0=rb(ib)%ve%ix0(ibasis)
          iy0=rb(ib)%ve%iy0(ibasis)
          dx=rb(ib)%ve%dx(ibasis)
          dy=rb(ib)%ve%dy(ibasis)
          DO iy=iy0,myb
            DO ix=ix0,mxb
              yloc=iy-iy0+dy-1.e-7
              IF (yloc<0) THEN 
                ibn=seam(ib)%segment(1)%ptr(1)
                yloc=yloc+rb(ibn)%my
                CALL lagr_quad_eval(rb(ibn)%rz,ix-ix0+dx,yloc,1_i4)
                tang=rb(ibn)%rz%fy
              ELSE
                CALL lagr_quad_eval(rb(ib)%rz,ix-ix0+dx,yloc,1_i4)
                tang=rb(ib)%rz%fy
              ENDIF
              yloc=iy-iy0+dy+1.e-7
              IF (yloc>myb) THEN
                yloc=yloc-myb
                ibn=seam(ib)%segment(mxb+myb+1)%ptr(1)
                CALL lagr_quad_eval(rb(ibn)%rz,ix-ix0+dx,yloc,1_i4)
                tang=0.5*(tang+rb(ibn)%rz%fy)
              ELSE
                CALL lagr_quad_eval(rb(ib)%rz,ix-ix0+dx,yloc,1_i4)
                tang=0.5*(tang+rb(ib)%rz%fy)
              ENDIF
              temp(1)=SQRT(SUM(tang**2))
              IF (temp(1)<rzmax*1.e-13) THEN
                CALL lagr_quad_eval(rb(ib)%rz,ix-ix0+dx+1.e-7_r8,
     $                              iy-iy0+dy,1_i4)
                tang(:)=rb(ib)%rz%fy(:)
                temp(1)= SQRT(SUM(tang(:)**2))
              ENDIF
              tang(1)= tang(1)/temp(1)
              tang(2)= tang(2)/temp(1)
              norm(1)= tang(2)
              norm(2)=-tang(1)
              CALL lagr_quad_eval(rb(ib)%ve,ix-ix0+dx,iy-iy0+dy,0_i4)
              CALL lagr_quad_eval(rb(ib)%be,ix-ix0+dx,iy-iy0+dy,0_i4)
              CALL lagr_quad_eval(rb(ib)%ja,ix-ix0+dx,iy-iy0+dy,0_i4)
              DO im=1,nm
                vetmp(1,ix,iy,im,ibasis)=SUM(norm*rb(ib)%ve%f(1:2,im))
                vetmp(2,ix,iy,im,ibasis)=SUM(tang*rb(ib)%ve%f(1:2,im))
                vetmp(3,ix,iy,im,ibasis)=rb(ib)%ve%f(3,im)
                betmp(1,ix,iy,im,ibasis)=SUM(norm*rb(ib)%be%f(1:2,im))
                betmp(2,ix,iy,im,ibasis)=SUM(tang*rb(ib)%be%f(1:2,im))
                betmp(3,ix,iy,im,ibasis)=rb(ib)%be%f(3,im)
                jatmp(1,ix,iy,im,ibasis)=SUM(norm*rb(ib)%ja%f(1:2,im))
                jatmp(2,ix,iy,im,ibasis)=SUM(tang*rb(ib)%ja%f(1:2,im))
                jatmp(3,ix,iy,im,ibasis)=rb(ib)%ja%f(3,im)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        DO ibasis=1,poly_degree**2
          ix0=rb(ib)%ve%ix0(ibasis)
          iy0=rb(ib)%ve%iy0(ibasis)
          CALL lagr_quad_basis_assign_arr(rb(ib)%ve,
     $      vetmp(:,ix0:mxb,iy0:myb,:,ibasis),ibasis)
          CALL lagr_quad_basis_assign_arr(rb(ib)%be,
     $      betmp(:,ix0:mxb,iy0:myb,:,ibasis),ibasis)
          CALL lagr_quad_basis_assign_arr(rb(ib)%ja,
     $      jatmp(:,ix0:mxb,iy0:myb,:,ibasis),ibasis)
        ENDDO
        DEALLOCATE(jatmp,betmp,vetmp)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE flux_comps
c-----------------------------------------------------------------------
c     subprogram 4. cart_comps.
c     convert the b, j and v at a given angle (what is loaded in the
c     first Fourier component storage) to Cartesian components
c     with consistent z axes for making 3D plots.
c-----------------------------------------------------------------------
      SUBROUTINE cart_comps
      
      INTEGER(i4) :: ib,im,nm,mxb,myb,ix,iy,iq
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: temp
      REAL(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: temp2
c-----------------------------------------------------------------------
c     loop over all blocks.
c-----------------------------------------------------------------------
      DO ib=1,nbl
        mxb=SIZE(ve(ib)%arr,2)
        myb=SIZE(ve(ib)%arr,3)
        ALLOCATE(temp(3,mxb,myb))
        temp=ve(ib)%arr(:,:,:,1)
        ve(ib)%arr(1,:,:,1)= temp(1,:,:)*COS(twopi*phi)
     $                      -temp(3,:,:)*SIN(twopi*phi)
        ve(ib)%arr(2,:,:,1)=-temp(1,:,:)*SIN(twopi*phi)
     $                      -temp(3,:,:)*COS(twopi*phi)
        ve(ib)%arr(3,:,:,1)= temp(2,:,:)
        temp=be(ib)%arr(:,:,:,1)
        be(ib)%arr(1,:,:,1)= temp(1,:,:)*COS(twopi*phi)
     $                      -temp(3,:,:)*SIN(twopi*phi)
        be(ib)%arr(2,:,:,1)=-temp(1,:,:)*SIN(twopi*phi)
     $                      -temp(3,:,:)*COS(twopi*phi)
        be(ib)%arr(3,:,:,1)= temp(2,:,:)
        temp=ja(ib)%arr(:,:,:,1)
        ja(ib)%arr(1,:,:,1)= temp(1,:,:)*COS(twopi*phi)
     $                      -temp(3,:,:)*SIN(twopi*phi)
        ja(ib)%arr(2,:,:,1)=-temp(1,:,:)*SIN(twopi*phi)
     $                      -temp(3,:,:)*COS(twopi*phi)
        ja(ib)%arr(3,:,:,1)= temp(2,:,:)
        DEALLOCATE(temp)
        IF (ASSOCIATED(ve(ib)%arrh)) THEN
          ALLOCATE(temp2(3,SIZE(ve(ib)%arrh,2),mxb-1,myb))
          temp2=ve(ib)%arrh(:,:,:,:,1)
          ve(ib)%arrh(1,:,:,:,1)= temp2(1,:,:,:)*COS(twopi*phi)
     $                           -temp2(3,:,:,:)*SIN(twopi*phi)
          ve(ib)%arrh(2,:,:,:,1)=-temp2(1,:,:,:)*SIN(twopi*phi)
     $                           -temp2(3,:,:,:)*COS(twopi*phi)
          ve(ib)%arrh(3,:,:,:,1)= temp2(2,:,:,:)
          temp2=be(ib)%arrh(:,:,:,:,1)
          be(ib)%arrh(1,:,:,:,1)= temp2(1,:,:,:)*COS(twopi*phi)
     $                           -temp2(3,:,:,:)*SIN(twopi*phi)
          be(ib)%arrh(2,:,:,:,1)=-temp2(1,:,:,:)*SIN(twopi*phi)
     $                           -temp2(3,:,:,:)*COS(twopi*phi)
          be(ib)%arrh(3,:,:,:,1)= temp2(2,:,:,:)
          temp2=ja(ib)%arrh(:,:,:,:,1)
          ja(ib)%arrh(1,:,:,:,1)= temp2(1,:,:,:)*COS(twopi*phi)
     $                           -temp2(3,:,:,:)*SIN(twopi*phi)
          ja(ib)%arrh(2,:,:,:,1)=-temp2(1,:,:,:)*SIN(twopi*phi)
     $                           -temp2(3,:,:,:)*COS(twopi*phi)
          ja(ib)%arrh(3,:,:,:,1)= temp2(2,:,:,:)
          DEALLOCATE(temp2)
        ENDIF
        IF (ASSOCIATED(ve(ib)%arrv)) THEN
          ALLOCATE(temp2(3,SIZE(ve(ib)%arrv,2),mxb,myb-1))
          temp2=ve(ib)%arrv(:,:,:,:,1)
          ve(ib)%arrv(1,:,:,:,1)= temp2(1,:,:,:)*COS(twopi*phi)
     $                           -temp2(3,:,:,:)*SIN(twopi*phi)
          ve(ib)%arrv(2,:,:,:,1)=-temp2(1,:,:,:)*SIN(twopi*phi)
     $                           -temp2(3,:,:,:)*COS(twopi*phi)
          ve(ib)%arrv(3,:,:,:,1)= temp2(2,:,:,:)
          temp2=be(ib)%arrv(:,:,:,:,1)
          be(ib)%arrv(1,:,:,:,1)= temp2(1,:,:,:)*COS(twopi*phi)
     $                           -temp2(3,:,:,:)*SIN(twopi*phi)
          be(ib)%arrv(2,:,:,:,1)=-temp2(1,:,:,:)*SIN(twopi*phi)
     $                           -temp2(3,:,:,:)*COS(twopi*phi)
          be(ib)%arrv(3,:,:,:,1)= temp2(2,:,:,:)
          temp2=ja(ib)%arrv(:,:,:,:,1)
          ja(ib)%arrv(1,:,:,:,1)= temp2(1,:,:,:)*COS(twopi*phi)
     $                           -temp2(3,:,:,:)*SIN(twopi*phi)
          ja(ib)%arrv(2,:,:,:,1)=-temp2(1,:,:,:)*SIN(twopi*phi)
     $                           -temp2(3,:,:,:)*COS(twopi*phi)
          ja(ib)%arrv(3,:,:,:,1)= temp2(2,:,:,:)
          DEALLOCATE(temp2)
        ENDIF
        IF (ASSOCIATED(ve(ib)%arri)) THEN
          ALLOCATE(temp2(3,SIZE(ve(ib)%arri,2),mxb-1,myb-1))
          temp2=ve(ib)%arri(:,:,:,:,1)
          ve(ib)%arri(1,:,:,:,1)= temp2(1,:,:,:)*COS(twopi*phi)
     $                           -temp2(3,:,:,:)*SIN(twopi*phi)
          ve(ib)%arri(2,:,:,:,1)=-temp2(1,:,:,:)*SIN(twopi*phi)
     $                           -temp2(3,:,:,:)*COS(twopi*phi)
          ve(ib)%arri(3,:,:,:,1)= temp2(2,:,:,:)
          temp2=be(ib)%arri(:,:,:,:,1)
          be(ib)%arri(1,:,:,:,1)= temp2(1,:,:,:)*COS(twopi*phi)
     $                           -temp2(3,:,:,:)*SIN(twopi*phi)
          be(ib)%arri(2,:,:,:,1)=-temp2(1,:,:,:)*SIN(twopi*phi)
     $                           -temp2(3,:,:,:)*COS(twopi*phi)
          be(ib)%arri(3,:,:,:,1)= temp2(2,:,:,:)
          temp2=ja(ib)%arri(:,:,:,:,1)
          ja(ib)%arri(1,:,:,:,1)= temp2(1,:,:,:)*COS(twopi*phi)
     $                           -temp2(3,:,:,:)*SIN(twopi*phi)
          ja(ib)%arri(2,:,:,:,1)=-temp2(1,:,:,:)*SIN(twopi*phi)
     $                           -temp2(3,:,:,:)*COS(twopi*phi)
          ja(ib)%arri(3,:,:,:,1)= temp2(2,:,:,:)
          DEALLOCATE(temp2)
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cart_comps
c-----------------------------------------------------------------------
c     subprogram 5. add_eq_fields.
c     add the *_eq fields to the n=0 component of the solution.
c-----------------------------------------------------------------------
      SUBROUTINE add_eq_fields(geom,poly_degree)
      USE fft_mod
      USE global
      USE physdat

      CHARACTER(*), INTENT(IN) :: geom
      INTEGER(i4), INTENT(IN) :: poly_degree

      INTEGER(i4) :: ib,im,ix,iy,vlim,mxb,myb,ix0,iy0,ibasis
      REAL(r8) :: dx,dy,te,ti
      COMPLEX(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: jatmp,betmp,
     $             vetmp,prtmp,petmp,ndtmp,tetmp,titmp,w1tmp
c-----------------------------------------------------------------------
c     set vector component index limit for be and ja.
c-----------------------------------------------------------------------
      IF (geom=='tor') THEN
        vlim=2
      ELSE
        vlim=3
      ENDIF
c-----------------------------------------------------------------------
c     find the n=0 mode.
c-----------------------------------------------------------------------
      mode_loop: DO im=1,nmodes
        IF (keff(im)/=0) CYCLE mode_loop
c-----------------------------------------------------------------------
c       add the equilibrium fields to all basis functions in rblocks,
c       then zero-out equilibrium fields.
c
c       be_eq is added to work1 and saved as cylindrical components
c       of total B in case be is changed to flux components.
c-----------------------------------------------------------------------
        DO ib=1,nrbl
          mxb=rb(ib)%mx
          myb=rb(ib)%my
          ALLOCATE(jatmp(3,0:mxb,0:myb,poly_degree**2))
          ALLOCATE(betmp(3,0:mxb,0:myb,poly_degree**2))
          ALLOCATE(vetmp(3,0:mxb,0:myb,poly_degree**2))
          ALLOCATE(w1tmp(3,0:mxb,0:myb,poly_degree**2))
          ALLOCATE(prtmp(1,0:mxb,0:myb,poly_degree**2))
          ALLOCATE(petmp(1,0:mxb,0:myb,poly_degree**2))
          ALLOCATE(ndtmp(1,0:mxb,0:myb,poly_degree**2))
          ALLOCATE(tetmp(1,0:mxb,0:myb,poly_degree**2))
          ALLOCATE(titmp(1,0:mxb,0:myb,poly_degree**2))
          DO ibasis=1,poly_degree**2
            ix0=rb(ib)%ve%ix0(ibasis)
            iy0=rb(ib)%ve%iy0(ibasis)
            dx=rb(ib)%ve%dx(ibasis)
            dy=rb(ib)%ve%dy(ibasis)
            DO iy=iy0,myb
              DO ix=ix0,mxb
                CALL lagr_quad_eval(rb(ib)%be,ix-ix0+dx,iy-iy0+dy,0_i4)
                CALL lagr_quad_eval(rb(ib)%work1,ix-ix0+dx,iy-iy0+dy,
     $                              0_i4)
                CALL lagr_quad_eval(rb(ib)%ja,ix-ix0+dx,iy-iy0+dy,0_i4)
                CALL lagr_quad_eval(rb(ib)%ve,ix-ix0+dx,iy-iy0+dy,0_i4)
                CALL lagr_quad_eval(rb(ib)%pres,ix-ix0+dx,iy-iy0+dy,
     $                              0_i4)
                CALL lagr_quad_eval(rb(ib)%prese,ix-ix0+dx,iy-iy0+dy,
     $                              0_i4)
                CALL lagr_quad_eval(rb(ib)%nd,ix-ix0+dx,iy-iy0+dy,0_i4)
                CALL lagr_quad_eval(rb(ib)%tele,ix-ix0+dx,iy-iy0+dy,
     $                              0_i4)
                CALL lagr_quad_eval(rb(ib)%tion,ix-ix0+dx,iy-iy0+dy,
     $                              0_i4)
                CALL lagr_quad_eval(rb(ib)%be_eq,ix-ix0+dx,iy-iy0+dy,
     $                              0_i4)
                CALL lagr_quad_eval(rb(ib)%ja_eq,ix-ix0+dx,iy-iy0+dy,
     $                              0_i4)
                CALL lagr_quad_eval(rb(ib)%ve_eq,ix-ix0+dx,iy-iy0+dy,
     $                              0_i4)
                CALL lagr_quad_eval(rb(ib)%pres_eq,ix-ix0+dx,iy-iy0+dy,
     $                              0_i4)
                CALL lagr_quad_eval(rb(ib)%prese_eq,ix-ix0+dx,
     $                              iy-iy0+dy,0_i4)
                CALL lagr_quad_eval(rb(ib)%nd_eq,ix-ix0+dx,iy-iy0+dy,
     $                              0_i4)
                betmp(1:vlim,ix,iy,ibasis)=rb(ib)%be%f(1:vlim,im)
     $                                    +rb(ib)%be_eq%f(1:vlim)
                w1tmp(1:vlim,ix,iy,ibasis)=rb(ib)%work1%f(1:vlim,im)
     $                                    +rb(ib)%be_eq%f(1:vlim)
                jatmp(1:vlim,ix,iy,ibasis)=rb(ib)%ja%f(1:vlim,im)
     $                                    +rb(ib)%ja_eq%f(1:vlim)
                IF (geom=='tor') THEN
                  CALL lagr_quad_eval(rb(ib)%rz,ix-ix0+dx,iy-iy0+dy,
     $                                0_i4)
                  IF (rb(ib)%rz%f(1)>2*TINY(0._r8)) THEN
                    betmp(3,ix,iy,ibasis)=rb(ib)%be%f(3,im)
     $                +rb(ib)%be_eq%f(3)/rb(ib)%rz%f(1)
                    w1tmp(3,ix,iy,ibasis)=rb(ib)%work1%f(3,im)
     $                +rb(ib)%be_eq%f(3)/rb(ib)%rz%f(1)
                  ELSE
                    betmp(3,ix,iy,ibasis)=0
                    w1tmp(3,ix,iy,ibasis)=0
                  ENDIF
                  jatmp(3,ix,iy,ibasis)=rb(ib)%ja%f(3,im)
     $              +rb(ib)%ja_eq%f(3)*rb(ib)%rz%f(1)
                ENDIF
                vetmp(:,ix,iy,ibasis)=rb(ib)%ve%f(:,im)
     $                               +rb(ib)%ve_eq%f(:)
                prtmp(:,ix,iy,ibasis)=rb(ib)%pres%f(:,im)
     $                               +rb(ib)%pres_eq%f(:)
                petmp(:,ix,iy,ibasis)=rb(ib)%prese%f(:,im)
     $                               +rb(ib)%prese_eq%f(:)
                ndtmp(:,ix,iy,ibasis)=rb(ib)%nd%f(:,im)
     $                               +rb(ib)%nd_eq%f(:)
                tetmp(:,ix,iy,ibasis)=rb(ib)%tele%f(:,im)
     $                +rb(ib)%prese_eq%f(:)/(kboltz*rb(ib)%nd_eq%f(1))
                titmp(:,ix,iy,ibasis)=rb(ib)%tion%f(:,im)
     $                +(rb(ib)%pres_eq%f(:)-rb(ib)%prese_eq%f(:))*
     $                 zeff/(kboltz*rb(ib)%nd_eq%f(1))
              ENDDO
            ENDDO
          ENDDO
          DO ibasis=1,poly_degree**2
            ix0=rb(ib)%ve%ix0(ibasis)
            iy0=rb(ib)%ve%iy0(ibasis)
            DO iy=iy0,myb
              DO ix=ix0,mxb
                CALL lagr_quad_basis_assign_loc(rb(ib)%be,
     $            betmp(:,ix,iy,ibasis),ibasis,ix,iy,im)
                CALL lagr_quad_basis_assign_loc(rb(ib)%work1,
     $            w1tmp(:,ix,iy,ibasis),ibasis,ix,iy,im)
                CALL lagr_quad_basis_assign_loc(rb(ib)%ja,
     $            jatmp(:,ix,iy,ibasis),ibasis,ix,iy,im)
                CALL lagr_quad_basis_assign_loc(rb(ib)%ve,
     $            vetmp(:,ix,iy,ibasis),ibasis,ix,iy,im)
                CALL lagr_quad_basis_assign_loc(rb(ib)%pres,
     $            prtmp(:,ix,iy,ibasis),ibasis,ix,iy,im)
                CALL lagr_quad_basis_assign_loc(rb(ib)%prese,
     $            petmp(:,ix,iy,ibasis),ibasis,ix,iy,im)
                CALL lagr_quad_basis_assign_loc(rb(ib)%nd,
     $            ndtmp(:,ix,iy,ibasis),ibasis,ix,iy,im)
                CALL lagr_quad_basis_assign_loc(rb(ib)%tele,
     $            tetmp(:,ix,iy,ibasis),ibasis,ix,iy,im)
                CALL lagr_quad_basis_assign_loc(rb(ib)%tion,
     $            titmp(:,ix,iy,ibasis),ibasis,ix,iy,im)
              ENDDO
            ENDDO
          ENDDO
          DEALLOCATE(jatmp,betmp,vetmp,prtmp,petmp,ndtmp,tetmp,titmp,
     $               w1tmp)
c-----------------------------------------------------------------------
c         zero-out the equilibrium fields, so that other plotting
c         options still work correctly.
c-----------------------------------------------------------------------
          rb(ib)%be_eq=0
          rb(ib)%ve_eq=0
          rb(ib)%ja_eq=0
          rb(ib)%pres_eq=0
          rb(ib)%prese_eq=0
          rb(ib)%nd_eq=0
        ENDDO
c-----------------------------------------------------------------------
c       add the equilibrium fields to all basis functions in tblocks.
c-----------------------------------------------------------------------
        DO ib=nrbl+1,nbl
          be(ib)%arr(1:vlim,:,:,im)=be(ib)%arr(1:vlim,:,:,im)
     $                             +be_eq(ib)%arr(1:vlim,:,:)
          work1(ib)%arr(1:vlim,:,:,im)=work1(ib)%arr(1:vlim,:,:,im)
     $                                +be_eq(ib)%arr(1:vlim,:,:)
          ja(ib)%arr(1:vlim,:,:,im)=ja(ib)%arr(1:vlim,:,:,im)
     $                             +ja_eq(ib)%arr(1:vlim,:,:)
          IF (geom=='tor') THEN
            DO iy=0,tb(ib)%mvert
              IF (tb(ib)%tgeom%xs(iy)>2*TINY(0._r8)) THEN
                be(ib)%arr(3,iy,0,im)=be(ib)%arr(3,iy,0,im)
     $            +be_eq(ib)%arr(3,iy,0)/tb(ib)%tgeom%xs(iy)
                work1(ib)%arr(3,iy,0,im)=work1(ib)%arr(3,iy,0,im)
     $            +be_eq(ib)%arr(3,iy,0)/tb(ib)%tgeom%xs(iy)
              ELSE
                be(ib)%arr(3,iy,0,im)=0
                work1(ib)%arr(3,iy,0,im)=0
              ENDIF
            ENDDO
            ja(ib)%arr(3,:,0,im)=ja(ib)%arr(3,:,0,im)
     $        +ja_eq(ib)%arr(3,:,0)*tb(ib)%tgeom%xs(:)
          ENDIF
          ve(ib)%arr(:,:,:,im)=ve(ib)%arr(:,:,:,im)+ve_eq(ib)%arr
          pres(ib)%arr(:,:,:,im)=pres(ib)%arr(:,:,:,im)+pres_eq(ib)%arr
          prese(ib)%arr(:,:,:,im)=prese(ib)%arr(:,:,:,im)
     $                           +prese_eq(ib)%arr
          nd(ib)%arr(:,:,:,im)=nd(ib)%arr(:,:,:,im)+nd_eq(ib)%arr
          tele(ib)%arr(:,:,:,im)=tele(ib)%arr(:,:,:,im)
     $                +prese_eq(ib)%arr/(kboltz*nd_eq(ib)%arr)
          tion(ib)%arr(:,:,:,im)=tion(ib)%arr(:,:,:,im)
     $                +(pres_eq(ib)%arr-prese_eq(ib)%arr)*
     $                 zeff/(kboltz*nd_eq(ib)%arr)
c-----------------------------------------------------------------------
c         zero-out the equilibrium fields, so that other plotting
c         options still work correctly.
c-----------------------------------------------------------------------
          tb(ib)%be_eq=0
          tb(ib)%ve_eq=0
          tb(ib)%ja_eq=0
          tb(ib)%pres_eq=0
          tb(ib)%prese_eq=0
          tb(ib)%nd_eq=0
        ENDDO
      ENDDO mode_loop
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE add_eq_fields
c-----------------------------------------------------------------------
c     subprogram 6. eval_at_angle.
c     evaluate some quantity at a given angle, and place that value
c     in the first Fourier component.
c
c     note that phi is in radians/2*pi.  [0-1]
c-----------------------------------------------------------------------
      SUBROUTINE eval_at_angle(vec,phi)

      TYPE(cvector_type), DIMENSION(:), POINTER :: vec
      REAL(r8), INTENT(IN) :: phi

      INTEGER(i4) :: ibl,im
      REAL(r8) :: angle
c-----------------------------------------------------------------------
c     all blocks, grid centered bases first.
c
c     there are no separate Fourier layers in nimplot, so we may assume
c     that n=0 is located with index 1.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        IF (ASSOCIATED(vec(ibl)%arr)) THEN
          DO im=2,nmodes
            angle=twopi*(im-1)*phi
            vec(ibl)%arr(:,:,:,1)=vec(ibl)%arr(:,:,:,1)
     $        +2*(COS(angle)*REAL(vec(ibl)%arr(:,:,:,im),r8)
     $           -SIN(angle)*AIMAG(vec(ibl)%arr(:,:,:,im)))
          ENDDO
          vec(ibl)%arr(:,:,:,1)=REAL(vec(ibl)%arr(:,:,:,1),r8)
        ENDIF
c-----------------------------------------------------------------------
c       horizontal side-centered bases.
c-----------------------------------------------------------------------
        IF (ASSOCIATED(vec(ibl)%arrh)) THEN
          DO im=2,nmodes
            angle=twopi*(im-1)*phi
            vec(ibl)%arrh(:,:,:,:,1)=vec(ibl)%arrh(:,:,:,:,1)
     $        +2*(COS(angle)*REAL(vec(ibl)%arrh(:,:,:,:,im),r8)
     $           -SIN(angle)*AIMAG(vec(ibl)%arrh(:,:,:,:,im)))
          ENDDO
          vec(ibl)%arrh(:,:,:,:,1)=REAL(vec(ibl)%arrh(:,:,:,:,1),r8)
        ENDIF
c-----------------------------------------------------------------------
c       vertical side-centered bases.
c-----------------------------------------------------------------------
        IF (ASSOCIATED(vec(ibl)%arrv)) THEN
          DO im=2,nmodes
            angle=twopi*(im-1)*phi
            vec(ibl)%arrv(:,:,:,:,1)=vec(ibl)%arrv(:,:,:,:,1)
     $        +2*(COS(angle)*REAL(vec(ibl)%arrv(:,:,:,:,im),r8)
     $           -SIN(angle)*AIMAG(vec(ibl)%arrv(:,:,:,:,im)))
          ENDDO
          vec(ibl)%arrv(:,:,:,:,1)=REAL(vec(ibl)%arrv(:,:,:,:,1),r8)
        ENDIF
c-----------------------------------------------------------------------
c       interior centered bases.
c-----------------------------------------------------------------------
        IF (ASSOCIATED(vec(ibl)%arri)) THEN
          DO im=2,nmodes
            angle=twopi*(im-1)*phi
            vec(ibl)%arri(:,:,:,:,1)=vec(ibl)%arri(:,:,:,:,1)
     $        +2*(COS(angle)*REAL(vec(ibl)%arri(:,:,:,:,im),r8)
     $           -SIN(angle)*AIMAG(vec(ibl)%arri(:,:,:,:,im)))
          ENDDO
          vec(ibl)%arri(:,:,:,:,1)=REAL(vec(ibl)%arri(:,:,:,:,1),r8)
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE eval_at_angle
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE plot_data
