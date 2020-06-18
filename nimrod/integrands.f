c-----------------------------------------------------------------------
c     file integrands.f
c     this is now just the module that collects the modules for the
c     integrands needed for the physics kernel, which have been
c     separated to avoid unnecessary compiler/optimization time.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     global module referencing mat, rhs, and dot modules.
c-----------------------------------------------------------------------
      MODULE integrands

      USE integrands_rhs
      USE integrands_dot
      USE integrands_mat
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE integrands
