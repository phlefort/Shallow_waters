!-------------------
!  type double precision portable
!-------------------

module mod_precision

  implicit none

  integer, parameter :: pr=selected_real_kind(15,70)

end module mod_precision
