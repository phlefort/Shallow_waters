module mod_flux
  use mod_precision
  use mod_variables
  
  implicit none
!ne change pas de la version du prof avec maillage structure
!reste a determiner le u dans la fonction
!n normeale a larrete
!fhp correspond à la première composante de Fn+
!fhm correspond à la première composante de Fn-
!fup correspond à la deuxième et troisième composante de Fu+ avec u=un et u=ut
!fum correspond à la deuxième et troisième composante de Fu- avec u=un et u=ut
contains

  function fhp(u,a,b) 

    implicit none

    !---------------------- in ----------------------------!
    real(pr) :: u,a,b
    !------------------------------------------------------!

    !---------------------- out ---------------------------!
    real(pr) :: fhp
    !------------------------------------------------------!

    if (u>=b) then
       fhp=2._pr*u*b;
    elseif (u<=-b) then
       fhp=0._pr;
    else
       fhp=(u+b)**2/2._pr;
    end if

    fhp=2._pr*a*b*fhp;

  end function fhp

  function fhm(u,a,b)

    implicit none

    !---------------------- in ----------------------------!
    real(pr) :: u,a,b
    !------------------------------------------------------!

    !---------------------- out ---------------------------!
    real(pr) :: fhm
    !------------------------------------------------------!

    if (u>=b) then
       fhm=0._pr
    elseif (u<=-b) then
       fhm=2._pr*u*b;
    else
       fhm=-(u-b)**2/2._pr;
    end if

    fhm=2._pr*a*b*fhm;

  end function fhm

  function fup(u,a,b) 

    implicit none

    !---------------------- in ----------------------------!
    real(pr) :: u,a,b
    !------------------------------------------------------!

    !---------------------- out ---------------------------!
    real(pr) :: fup
    !------------------------------------------------------!

    if (u>=b) then
       fup=2._pr*b/3._pr*(3._pr*u**2+b**2);
    elseif (u<=-b) then
       fup=0._pr;
    else
       fup=(u+b)**3/3._pr;
    end if

    fup=2._pr*a*b*fup;

  end function fup

  function fum(u,a,b) result(f)

    implicit none

    !---------------------- in ----------------------------!
    real(pr) :: u,a,b
    !------------------------------------------------------!

    !---------------------- out ---------------------------!
    real(pr) :: f
    !------------------------------------------------------!

    if (u>=b) then
       f=0._pr;
    elseif (u<=-b) then
       f=2._pr*b/3._pr*(3._pr*u**2+b**2);
    else
       f=-(u-b)**3/3._pr;
    end if

    f=2._pr*a*b*f;
  end function fum



end module mod_flux
