!----------------------------------------------------------
!ligne 55
!----------------------------------------------------------

program masca

  !----------- sous programmes ------------------------------------ 
  use mod_precision   !- (double ou simple précision)
  use mod_variables
!  use mod_vitesse     !- (champ de vitesse)
!  use mod_ci          !- (donnée initiale)
!  use mod_aplusmoins  !- (calcul de a+ et a-)
  use mod_io    !- résultats de calculs
  use mod_calcul
  use mod_flux
  !----------------------------------------------------------------

  implicit none

  !------------ variables -----------------------------------------
  type(maillage)::mesh
  type(modele)  ::var
  integer :: i,iter,itermax
  integer :: nmax,nfich
  character(len=100) :: fich
  real(pr) :: vmax,tmax,norme_vitesse,stationnaire_erreur
  real(pr), dimension(:), allocatable :: Ux_stat, Uy_stat
  logical::champs_vitesse,stationnaire_debit
  !----------------------------------------------------------------


  !------------- programme ----------------------------------------
  print*,"*************************************************************************************************************"
  print*,"********************************************DEBUT DU CALCUL**************************************************"
  print*,"*************************************************************************************************************"


  !_________________________________________________!
  !             INITIALISATION                      !
  !_________________________________________________! 
  call lecture_donnees(var,mesh)
  call lecture_maillage(mesh,var)

  ! ****Initialisation du champs de vitesse dans le domaine*** 
  allocate(Ux_stat(mesh%nb_trig))
  allocate(Uy_stat(mesh%nb_trig))
  var%iter=0
  itermax=1000000
  stationnaire_erreur=0.000000000000001
  norme_vitesse=1.
  nfich=0
  var%stationnaire=.false.
  call condition_initiale(var,mesh)
  stationnaire_debit=.true.
  !call debit(var,mesh)
  
  
  if (var%lecture_stationnaire=="non") then
    do while (var%iter<var%iter_stationnairemax.and.norme_vitesse>stationnaire_erreur.and.stationnaire_debit)
      var%iter=var%iter+1
      Ux_stat(:)=var%Ux(:)
      Uy_stat(:)=var%Uy(:)
      
      !Actualisation du pas de temps
      var%dt=mesh%deltah/maxval(sqrt(var%ux**2+var%uy**2)+sqrt((3./2.)*var%g*var%h))
      var%t=var%t+var%dt

      call flux(var,mesh)  

      !changement en non conservatif
      call var_non_conservatives(var,mesh)

      ! Norme de la variation du champs du vitesse 
      norme_vitesse=dot_product(var%Ux(1:mesh%nb_trig)-Ux_stat(1:mesh%nb_trig),var%Ux(1:mesh%nb_trig)-Ux_stat(1:mesh%nb_trig))+&
                  & dot_product(var%Uy(1:mesh%nb_trig)-Uy_stat(1:mesh%nb_trig),var%Uy(1:mesh%nb_trig)-Uy_stat(1:mesh%nb_trig))

      if (var%iter==var%pas_affichage*nfich+1) then 
        call sortie_visit(nfich,mesh,var)
        
        print*,"@@ Stationnaire n° iteration=",var%iter,"pour t=",var%t," un dt de",&
               &  var%dt," et une variation de vitesse de ",norme_vitesse
        call debits(var,mesh)
        nfich=nfich+1
      end if 

      do i=1,mesh%nb_trig
        var%Un(i,1)=var%Unp1(i,1)
        var%Un(i,2)=var%Unp1(i,2)
        var%Un(i,3)=var%Unp1(i,3)
      end do 
      stationnaire_debit=.not.debit(var,mesh)
    end do 
    deallocate(Ux_stat,Uy_stat)
    call sauvegarde_stationnaire(mesh,var)
    var%stationnaire=.true.

  else if (var%lecture_stationnaire=="oui") then
    call lecture_stationnaire(mesh,var)
    var%stationnaire=.true.
  end if

  !_________________________________________________!
  !                    Mascaret                     !
  !_________________________________________________! 

  var%iter=0
  nfich=0
  var%t=0.
    do while (var%t<var%tmax)!.and.norme_vitesse>stationnaire_erreur)
      var%iter=var%iter+1
      
      !Actualisation du pas de temps
      var%dt=mesh%deltah/maxval(sqrt(var%ux**2+var%uy**2)+sqrt((3./2.)*var%g*var%h))
      if (var%t+var%dt>var%tmax) var%dt=var%tmax-var%t
      var%t=var%t+var%dt

      call flux(var,mesh)  

      !changement en non conservatif
      call var_non_conservatives(var,mesh)

      if (var%iter==var%pas_affichage*nfich+1) then 
        call sortie_visit(nfich,mesh,var)
        
        print*,"@@ mascaret n° iteration=",var%iter,"pour t=",var%t," un dt de",&
               &  var%dt
        call debits(var,mesh)
        nfich=nfich+1
      end if 

      do i=1,mesh%nb_trig
        var%Un(i,1)=var%Unp1(i,1)
        var%Un(i,2)=var%Unp1(i,2)
        var%Un(i,3)=var%Unp1(i,3)
      end do 
    end do 

  !----------------------------------------------------------------
  print*,"*************************************************************************************************************"
  print*,"**********************************************FIN DU CALCUL**************************************************"
  print*,"*************************************************************************************************************"

end program masca
