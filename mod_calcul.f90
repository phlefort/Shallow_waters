module mod_calcul
  use mod_precision
  use mod_variables
  use mod_flux
  implicit none

contains
  subroutine condition_initiale(var,mesh)
    implicit none 
    type(maillage)::mesh
    type(modele) ::var
    integer::i

    allocate(var%h(1:mesh%nb_trig))
    allocate(var%Ux(mesh%nb_trig))
    allocate(var%Uy(mesh%nb_trig))
    allocate(var%fond(mesh%nb_trig))!le fond est plat il restera a zero
    allocate(var%b(mesh%nb_cotes))
    allocate(var%Un(mesh%nb_trig,3))!au temps n
    allocate(var%Unp1(mesh%nb_trig,3))!au temps n+1
    allocate(var%vnp(1:mesh%nb_cotes),var%vnm(1:mesh%nb_cotes),var%fh(1:mesh%nb_cotes),var%fun(1:mesh%nb_cotes))
    allocate(var%fut(1:mesh%nb_cotes),var%fux(1:mesh%nb_cotes),var%fuy(1:mesh%nb_cotes))
    var%fond=0
    var%g=10.
    var%ux=0
    var%uy=0
    var%h =1.
    var%a=1./(6*var%g)

    do i=1,mesh%nb_trig
      var%Un(i,1)=var%h(i)
      var%Un(i,2)=var%h(i)*var%Ux(i)
      var%Un(i,3)=var%h(i)*var%Uy(i)  

      var%Unp1(i,1)=var%h(i)
      var%Unp1(i,2)=var%h(i)*var%Ux(i)
      var%Unp1(i,3)=var%h(i)*var%Uy(i)    
   end do 

    open(unit=1,file='SORTIES/gnuplot.100.h')
    do i=1,mesh%nb_trig
      write(1,*) mesh%xc(i) ,mesh%yc(i) ,var%h(i)
    end do
    close(1)
  end subroutine condition_initiale

  subroutine var_non_conservatives(var,mesh)
    implicit none 
    type(maillage)::mesh
    type(modele) :: var
    integer:: i
!
    do i=1,mesh%nb_trig
if (var%h(i)/=0._pr) then
      var%h(i)=var%Unp1(i,1)
      var%Ux(i)=var%Unp1(i,2)/var%Unp1(i,1)
      var%Uy(i)=var%Unp1(i,3)/var%Unp1(i,1)
else
      var%Ux(i)=0
      var%Uy(i)=0
end if
    end do
  end subroutine var_non_conservatives

  function mascaret(var,mesh)
    implicit none
    type(maillage)::mesh
    type(modele) :: var 
    real(pr),dimension(1:3)::mascaret   

  end function mascaret

  subroutine flux(var,mesh)
    implicit none 
    type(maillage)::mesh
    type(modele) :: var
    integer:: arete,i,k,tri1,tri2
    real(pr) :: vn,vt,Rh,Run,Rut,sf,dt,b1,b2,Rux,Ruy,AB
    real(pr), dimension(1:2) :: v,t,A,B
    real(pr), dimension(1:3) :: U1,U2

 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!calcul des flux sur toutes les aretes!!!!!!!!!!!!!!!!!!
    !--- calcul de la vitesse normale aux arêtes

    do arete=1,mesh%nb_cotes
      !-- indices des triangles voisins
      tri1=mesh%trig(arete,1)
      tri2=mesh%trig(arete,2)

      !vecteur t
      t(1)=-mesh%normale(arete,2)
      t(2)=mesh%normale(arete,1) 

      !-- valeurs dans les triangles voisins
      U1(1)=var%h(tri1)
      U1(2)=var%Ux(tri1)*mesh%normale(arete,1)+var%Uy(tri1)*mesh%normale(arete,2)!un
      U1(3)=var%Ux(tri1)*t(1)+var%Uy(tri1)*t(2)!ut

      if (tri2==0) then ! alors arete est sur le bord
        If (mesh%refs(mesh%e(arete,1))==mesh%b_ouest.or.mesh%refs(mesh%e(arete,2))==mesh%b_ouest) then ! CL d'entrée !!!!!!!!!!!! refs sur les noeuds
          U2(1)=var%h_courant
          U2(2)=-var%V_courant
          U2(3)=U1(3)

        else if (mesh%refs(mesh%e(arete,1))==mesh%b_nord.or.mesh%refs(mesh%e(arete,2))==mesh%b_nord) then 
          U2(1)=U1(1)
          U2(2)=-U1(2)
          U2(3)=U1(3)

        else if (mesh%refs(mesh%e(arete,1))==mesh%b_sud.or.mesh%refs(mesh%e(arete,2))==mesh%b_sud) then ! CL de bords
          U2(1)=U1(1)
          U2(2)=-U1(2)
          U2(3)=U1(3)

        else if (mesh%refs(mesh%e(arete,1))==mesh%b_est.or.mesh%refs(mesh%e(arete,2))==mesh%b_est) then ! CL de sortie
          U2(1)=U1(1)
          U2(2)=U1(2)
          U2(3)=U1(3)
          if (var%stationnaire) then
          if(abs(var%t)<1) then
          U2(1)=var%h_mascaret
          else
          U2(1)=1
          end if
          end if
          if (var%stationnaire) U2(2)=U1(2)
          if (var%stationnaire) U2(3)= U1(3)
        end if        
      else
        U2(1)=var%h(tri2)
        U2(2)=var%Ux(tri2)*mesh%normale(arete,1)+var%Uy(tri2)*mesh%normale(arete,2)
        U2(3)=var%Ux(tri2)*t(1)+var%Uy(tri2)*t(2)        
      end if

      b1=sqrt(1.5*var%g*U1(1))
      b2=sqrt(1.5*var%g*U2(1))

      !flux

      var%fh(arete)=fhp(U1(2),var%a,b1)+fhm(U2(2),var%a,b2)
      var%fun(arete)=fup(U1(2),var%a,b1)+fum(U2(2),var%a,b2)
      var%fut(arete)=U1(3)*fhp(U1(2),var%a,b1)+U2(3)*fhm(U2(2),var%a,b2)
      var%fux(arete)=var%fun(arete)*mesh%normale(arete,1)-var%fut(arete)*mesh%normale(arete,2)
      var%fuy(arete)=var%fun(arete)*mesh%normale(arete,2)+var%fut(arete)*mesh%normale(arete,1)
        

    end do


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !-- boucle sur les triangles
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    boucle_triangle:do i=1,mesh%nb_trig

      !-- boucle sur les arêtes du triangle i
      Rh=0
      Rux=0!
      Ruy=0!

      boucle_arete:do k=1,3

        !-- numero global de l'arête
        arete=mesh%ar(i,k)

        !-- numeros globaux des deux triangles contenant arête
        tri1=mesh%trig(arete,1)
        tri2=mesh%trig(arete,2)
        A=mesh%P(mesh%e(arete,1),:)
        B=mesh%P(mesh%e(arete,2),:)
	AB=sqrt(sum((B-A)**2))

        !-- on cherche lequel est i
        if (tri1==i) then  ! (contribution + du flux)
          sf=1
        else               ! (contribution - du flux)
          sf=-1
        end if

        !-- mise a  jour du residu
        Rh=Rh+sf*var%fh(arete)*AB
        Rux=Rux+sf*var%fux(arete)*AB
        Ruy=Ruy+sf*var%fuy(arete)*AB
      
      end do boucle_arete

      var%Unp1(i,1)=var%Un(i,1)-var%dt/mesh%aire(i)*Rh
      var%Unp1(i,2)=var%Un(i,2)-var%dt/mesh%aire(i)*Rux
      var%Unp1(i,3)=var%Un(i,3)-var%dt/mesh%aire(i)*Ruy
      
    end do boucle_triangle
end subroutine flux
function debit(var,mesh)
 implicit none 
    type(maillage)::mesh
    type(modele) :: var
    integer:: arete,n1,n2,tri1,tri2
    real(pr) :: vn,vt,Rh,Run,Rut,sf,dt,b1,b2,Rux,Ruy,AB,h1,h2,V1,V2,CD
    real(pr), dimension(1:2) :: v,t,A,B,debi,C,D
    real(pr), dimension(1:3) :: U1,U2
    logical :: debit
 n1=0
 n2=0
 V1=0
 V2=0
 AB=0
 CD=0
 h1=0
 h2=0
 debi=0
  debit=.false.

 do arete=1,mesh%nb_cotes
      tri1=mesh%trig(arete,1)
      tri2=mesh%trig(arete,2)
        IF (tri2==0) then
        If (mesh%refs(mesh%e(arete,1))==mesh%b_ouest.or.mesh%refs(mesh%e(arete,2))==mesh%b_ouest) then ! CL d'entrée !!!!!!!!!!!! refs sur les noeuds
         n1=n1+1
         A=mesh%P(mesh%e(arete,1),:)
         B=mesh%P(mesh%e(arete,2),:)
	 AB=sqrt(sum((B-A)**2))+AB
         V1=var%Ux(tri1)*mesh%normale(arete,1)+var%Uy(tri1)*mesh%normale(arete,2)+V1
         h1=h1+var%h(tri1)
 debi(1)=h1*AB*V1+debi(1)
!print*,mesh%normale(arete,1),mesh%normale(arete,2)
         

        else if (mesh%refs(mesh%e(arete,1))==mesh%b_est.or.mesh%refs(mesh%e(arete,2))==mesh%b_est) then ! CL de sortie
         n2=n2+1
         V2=var%Ux(tri1)*mesh%normale(arete,1)+var%Uy(tri1)*mesh%normale(arete,2)+V2
         h2=h2+var%h(tri1)
         C=mesh%P(mesh%e(arete,1),:)
         D=mesh%P(mesh%e(arete,2),:)
	 CD=sqrt(sum((D-C)**2))+CD
     debi(2)=h2*CD*V2+debi(2)
!print*,mesh%normale(arete,1),mesh%normale(arete,2)
           
end if   
end if
end do
 

debi(1)=debi(1)/n1
debi(2)=-debi(2)/n2

if (abs(debi(1)-debi(2))<=0.001) debit=.true.



end function debit
subroutine debits(var,mesh)
 implicit none 
    type(maillage)::mesh
    type(modele) :: var
    integer:: arete,n1,n2,tri1,tri2
    real(pr) :: vn,vt,Rh,Run,Rut,sf,dt,b1,b2,Rux,Ruy,AB,h1,h2,V1,V2,CD
    real(pr), dimension(1:2) :: v,t,A,B,debi,C,D
    real(pr), dimension(1:3) :: U1,U2

 n1=0
 n2=0
 V1=0
 V2=0
 AB=0
 CD=0
 h1=0
 h2=0
 debi=0


 do arete=1,mesh%nb_cotes
      tri1=mesh%trig(arete,1)
      tri2=mesh%trig(arete,2)
        IF (tri2==0) then
        If (mesh%refs(mesh%e(arete,1))==mesh%b_ouest.or.mesh%refs(mesh%e(arete,2))==mesh%b_ouest) then ! CL d'entrée !!!!!!!!!!!! refs sur les noeuds
         n1=n1+1
         A=mesh%P(mesh%e(arete,1),:)
         B=mesh%P(mesh%e(arete,2),:)
	 AB=sqrt(sum((B-A)**2))+AB
         V1=var%Ux(tri1)*mesh%normale(arete,1)+var%Uy(tri1)*mesh%normale(arete,2)+V1
         h1=h1+var%h(tri1)
 debi(1)=h1*AB*V1+debi(1)
!print*,mesh%normale(arete,1),mesh%normale(arete,2)
         

        else if (mesh%refs(mesh%e(arete,1))==mesh%b_est.or.mesh%refs(mesh%e(arete,2))==mesh%b_est) then ! CL de sortie
         n2=n2+1
         V2=var%Ux(tri1)*mesh%normale(arete,1)+var%Uy(tri1)*mesh%normale(arete,2)+V2
         h2=h2+var%h(tri1)
         C=mesh%P(mesh%e(arete,1),:)
         D=mesh%P(mesh%e(arete,2),:)
	 CD=sqrt(sum((D-C)**2))+CD
     debi(2)=h2*CD*V2+debi(2)
!print*,mesh%normale(arete,1),mesh%normale(arete,2)
           
end if   
end if
end do
 

debi(1)=debi(1)/n1
debi(2)=-debi(2)/n2



!print*,"@@ debit entree",debi(1),"debit sortie",debi(2)

end subroutine debits
  

end module mod_calcul
