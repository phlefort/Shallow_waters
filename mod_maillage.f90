module mod_maillage
  
  use mod_precision
  use mod_variables

  implicit none

contains 
  subroutine lecture_maillage(mesh,var)
    implicit none

    !------------ variables -----------------------------------------
    type(maillage)::mesh
    type(modele)  ::var
    integer :: i,j,n,k,arete,tri1,tri2,iA,iB,n1,n2,n3
    character(len=100) :: fich
    real(pr) :: AB,BC,AC,fb0,fb1,cfl,invdt,dt,R,fj,fi,vn,vnplus,vnmoins,temps,l
    real(pr), dimension(1:2) :: A,B,C,v,nij
    !----------------------------------------------------------------

    !--- lecture du maillage (format am_fmt donné par le mailleur EMC2)
    open(1,file=var%fich,form='formatted',status='old')
    read (1,*) mesh%nb_noeuds,mesh%nb_trig    ! nombre de noeuds (sommets) et nombre de triangles
    allocate(mesh%p(mesh%nb_noeuds,2))        ! coordonnées de chaque noeud
    allocate(mesh%s(mesh%nb_trig,3))          ! numéros des sommets de chaque triangle
    allocate(mesh%reft(mesh%nb_trig))         ! référence de chaque triangle (pour les CL) : non utilisé ici
    allocate(mesh%refs(mesh%nb_noeuds))       ! référence de chaque noeud (pour les CL)
    read (1,*) ((mesh%s(j,i),i=1,3),j=1,mesh%nb_trig)
    read (1,*) ((mesh%p(j,i),i=1,2),j=1,mesh%nb_noeuds)
    read (1,*) ( mesh%reft(i),i=1,mesh%nb_trig)
    read (1,*) ( mesh%refs(i),i=1,mesh%nb_noeuds)
    close(1)
    print*,""
    print*,"----------------Données du maillage---------------------------"
    print*,"nb noeuds :",mesh%nb_noeuds,", nb triangles :",mesh%nb_trig

    !--- calcul de la connectivité
  
    !1)-calcul du nombre d'arêtes mesh%nb_cotes
    mesh%nb_cotes=0
    allocate(mesh%flag(1:mesh%nb_noeuds,1:mesh%nb_noeuds))      
            ! mesh%flag(n1,n2) contient en fin de boucle le numéro de l'arête liant les sommets de numéros n1 et n2
    mesh%flag=0
    do i=1,mesh%nb_trig
       n1=mesh%s(i,1); n2=mesh%s(i,2); n3=mesh%s(i,3)
       if (mesh%flag(n1,n2)==0.or.mesh%flag(n2,n1)==0) then ! alors arete non numerotee
          mesh%nb_cotes=mesh%nb_cotes+1
          mesh%flag(n1,n2)=mesh%nb_cotes; mesh%flag(n2,n1)=mesh%nb_cotes
       end if
       if (mesh%flag(n2,n3)==0.or.mesh%flag(n3,n2)==0) then ! alors arete non numerotee
          mesh%nb_cotes=mesh%nb_cotes+1
          mesh%flag(n2,n3)=mesh%nb_cotes; mesh%flag(n3,n2)=mesh%nb_cotes
       end if
       if (mesh%flag(n3,n1)==0.or.mesh%flag(n1,n3)==0) then ! alors arete non numerotee
          mesh%nb_cotes=mesh%nb_cotes+1
          mesh%flag(n3,n1)=mesh%nb_cotes; mesh%flag(n1,n3)=mesh%nb_cotes
       end if
    end do

    !2)-remplissage de mesh%e
    allocate(mesh%e(1:mesh%nb_cotes,2))    ! numéros des noeuds extrémités de chaque arête
    do i=1,mesh%nb_trig
       n1=mesh%s(i,1); n2=mesh%s(i,2); n3=mesh%s(i,3)
       mesh%e(mesh%flag(n1,n2),1)=n1; mesh%e(mesh%flag(n1,n2),2)=n2
       mesh%e(mesh%flag(n2,n3),1)=n2; mesh%e(mesh%flag(n2,n3),2)=n3
       mesh%e(mesh%flag(n3,n1),1)=n3; mesh%e(mesh%flag(n3,n1),2)=n1
    end do

    !3)-remplissage de mesh%ar et mesh%trig
    allocate(mesh%ar(1:mesh%nb_trig,3))    ! numéros des arêtes de chaque triangle
    allocate(mesh%trig(1:mesh%nb_cotes,2)) ! numéro des triangles ayant l'arête i en commun
    mesh%trig=0
    do i=1,mesh%nb_trig
       n1=mesh%s(i,1); n2=mesh%s(i,2); n3=mesh%s(i,3)
       !-arete (n1,n2)
       mesh%ar(i,1)=mesh%flag(n1,n2)
       if (mesh%trig(mesh%ar(i,1),1)/=0) then
          mesh%trig(mesh%ar(i,1),2)=i
       else
          mesh%trig(mesh%ar(i,1),1)=i
       end if
       !-arete (n2,n3)
       mesh%ar(i,2)=mesh%flag(n2,n3)
       if (mesh%trig(mesh%ar(i,2),1)/=0) then
          mesh%trig(mesh%ar(i,2),2)=i
       else
          mesh%trig(mesh%ar(i,2),1)=i
       end if
       !-arete (n3,n1)
       mesh%ar(i,3)=mesh%flag(n3,n1)
       if (mesh%trig(mesh%ar(i,3),1)/=0) then
          mesh%trig(mesh%ar(i,3),2)=i
       else
          mesh%trig(mesh%ar(i,3),1)=i
       end if

    end do

    !--- calcul des aires des triangles et du "pas" du maillage
    allocate(mesh%aire(mesh%nb_trig),mesh%diam(mesh%nb_trig))
     do i=1,mesh%nb_trig
       A=mesh%p(mesh%s(i,1),:)
       B=mesh%p(mesh%s(i,2),:)
       C=mesh%p(mesh%s(i,3),:)
       mesh%aire(i)=0.5*abs((B(1)-A(1))*(C(2)-A(2))-(C(1)-A(1))*(B(2)-A(2)))
       AB=sqrt(sum((B-A)**2))
       BC=sqrt(sum((C-B)**2))
       AC=sqrt(sum((C-A)**2))
       mesh%diam(i)=2._pr*mesh%aire(i)/(AB+BC+AC)
    end do
    mesh%deltah=minval(mesh%diam)
    print*,"pas du maillage :",mesh%deltah

    !--- calcul des centres des triangles
    allocate(mesh%xc(mesh%nb_trig),mesh%yc(mesh%nb_trig))
    do i=1,mesh%nb_trig
       mesh%xc(i)=1._pr/3._pr*(mesh%p(mesh%s(i,1),1)+mesh%p(mesh%s(i,2),1)+mesh%p(mesh%s(i,3),1))
       mesh%yc(i)=1._pr/3._pr*(mesh%p(mesh%s(i,1),2)+mesh%p(mesh%s(i,2),2)+mesh%p(mesh%s(i,3),2))
    end do

    !--- calcul des normales aux arêtes (definies comme pointant du triangle 1 vers le triangle 2)
    !-- on utilise le rangement des noeuds du triangle 1 dans le sens trigo pour orienter
    !-- l'arête dans le sens trigo,
    !-- puis faire une rotation de -pi/2 de l'arête
    allocate(mesh%normale(1:mesh%nb_cotes,1:2))
    do arete=1,mesh%nb_cotes

       iA=mesh%e(arete,1)
       iB=mesh%e(arete,2)

       tri1=mesh%trig(arete,1)

       if (mesh%s(tri1,1)/=iA.and.mesh%s(tri1,1)/=iB) then
          iA=mesh%s(tri1,2)
          iB=mesh%s(tri1,3)
       elseif (mesh%s(tri1,2)/=iA.and.mesh%s(tri1,2)/=iB) then
          iA=mesh%s(tri1,3)
          iB=mesh%s(tri1,1)
       else
          iA=mesh%s(tri1,1)
          iB=mesh%s(tri1,2)
       end if
        mesh%normale(arete,1)=mesh%p(iB,2)-mesh%p(iA,2)
       mesh%normale(arete,2)=-(mesh%p(iB,1)-mesh%p(iA,1))
       
       l=sqrt(mesh%normale(arete,1)**2+mesh%normale(arete,2)**2)
       mesh%normale(arete,1)=(mesh%p(iB,2)-mesh%p(iA,2))/l
       mesh%normale(arete,2)=-(mesh%p(iB,1)-mesh%p(iA,1))/l

    end do
    print*,"--------------------------------------------------------------"
    print*,""
  end subroutine lecture_maillage

end module mod_maillage
