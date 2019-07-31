module mod_variables

  use mod_precision

  implicit none
  type maillage
    integer :: nb_noeuds
    integer :: nb_trig
    integer :: nb_cotes
    real(pr):: deltah ! pas du maillage
    integer :: b_nord, b_sud,b_est, b_ouest !refs pour identifier les CL

   !donner par emc2
    real(pr) , dimension(:,:), pointer :: p       => null() ! coordonnées des noeuds 
    integer , dimension(:,:), pointer :: s       => null() ! numero des sommets pour chaque triangle
    integer , dimension(:), pointer :: reft    => null() ! non utiliser ici
    integer , dimension(:), pointer :: refs    => null() ! reference de chaque noeud pour les CL

   !calculer par mod_maillage
    integer , dimension(:,:), pointer :: e       => null() ! numéros des noeuds extrémités de chaque arête
    integer , dimension(:,:), pointer :: ar      => null() ! numéros des arêtes de chaque triangle
    integer , dimension(:,:), pointer :: trig    => null() ! numéro des triangles ayant l'arête i en commun
    integer , dimension(:,:), pointer :: flag    => null() ! flag(n1,n2) contient en fin de boucle le numéro de l'arête liant les sommets de numéros n1 et n2
    real(pr) , dimension(:), pointer :: aire    => null() ! aires des triangles
    real(pr) , dimension(:), pointer :: diam    => null() ! pas du maillage
    real(pr) , dimension(:), pointer :: xc      => null() !--- calcul des centres des triangles--
    real(pr) , dimension(:), pointer :: yc      => null() !--------------------------------------
    real(pr) , dimension(:,:), pointer :: normale => null() ! normales aux arêtes (definies comme pointant du triangle 1 vers le triangle 2)
    
  end type maillage

  type modele
    real(pr) :: dt 
    real(pr) :: cfl
    real(pr) :: tmax
    integer  :: pas_affichage
    character(len=100) :: fich
    real(pr) ::a
    real(pr) ::g
    real(pr) ::t
    real(pr) :: V_courant,h_courant
    logical :: stationnaire
    character(len=4) :: lecture_stationnaire
    character(len=100) ::fichier_stationnaire
    integer ::iter,iter_stationnairemax
    real(pr) , dimension(:) , pointer :: h => null()
    real(pr) , dimension(:) , pointer :: ux,uy => null()
    real(pr) , dimension(:) , pointer :: b=> null()!a devient un reel
    real(pr) , dimension(:) , pointer :: fond => null()!correspond au grand B du cours
    real(pr) , dimension(:,:) , pointer ::Un =>null() !U conservatif
    real(pr) , dimension(:,:) , pointer ::Unp1 =>null()                                                      !U=(h  )
                                                     !  (hux)
                                                     !  (huy)
    real(pr) , dimension(:) , pointer ::vnp,vnm,fh,fun => null()
    real(pr) , dimension(:) , pointer ::fut,fux,fuy => null()
    real(pr) :: h_mascaret,V_mascaret

  end type modele
end module mod_variables
