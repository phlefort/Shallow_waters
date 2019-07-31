!-----------------------------------------------
! stockage des résultats (format visit)
!-----------------------------------------------


module mod_io
  use mod_precision
  use mod_maillage
  implicit none

  contains
  
    subroutine sortie_visit(nfich,mesh,var)

      !------------------------------------------------------
      !      sorties au format vtk (ascii)
      !      lisible par visit 
      !------------------------------------------------------

      !-------------------- modules -------------------------!
      use mod_precision
!      use mod_vitesse
      !------------------------------------------------------!

      implicit none

      !---------------------- in ----------------------------!
      type(maillage)::mesh
      type(modele) :: var
      integer, intent(in) :: nfich
      !------------------------------------------------------!

      !--------------- variables locales --------------------!
      integer :: i
      character(len=100) :: fich
      character(len=100) :: numero
      !------------------------------------------------------!

      write(numero,*) nfich
      if (nfich<=9) then
         numero='00'//trim(adjustl(numero))
      elseif (9<nfich.and.nfich<=99) then
         numero='0'//trim(adjustl(numero))
      end if

      if (var%stationnaire) then 
        fich='SORTIES/visit/mascaret.'//trim(adjustl(numero))//'.vtk'
      else 
        fich='SORTIES/visit/stationnaire.'//trim(adjustl(numero))//'.vtk' 
      end if   
      open(unit=1,file=fich)



      write(1,'(1A26)') '# vtk DataFile Version 2.0'
      write(1,*) 'convection 2D'
      write(1,*) 'ASCII'
      write(1,*) 'DATASET UNSTRUCTURED_GRID'

      write(1,*) 'POINTS',mesh%nb_noeuds,' double'
      do i=1,mesh%nb_noeuds
         write(1,*) mesh%p(i,1),mesh%p(i,2),0.0_pr
      end do

      write(1,*) 'CELLS ',mesh%nb_trig,mesh%nb_trig*4
      do i=1,mesh%nb_trig
         write(1,*) 3,mesh%s(i,1)-1, mesh%s(i,2)-1, mesh%s(i,3)-1
      end do

      write(1,*) 'CELL_TYPES ',mesh%nb_trig
      do i=1,mesh%nb_trig
         write(1,*) 5
      end do

      write(1,*) 'CELL_DATA',mesh%nb_trig
      write(1,*) 'SCALARS hauteur double'
      write(1,*) 'LOOKUP_TABLE default'
      do i=1,mesh%nb_trig
         write(1,*) var%h(i)+var%fond(i)
      end do

      write(1,*) 'SCALARS fond double'
      write(1,*) 'LOOKUP_TABLE default'
      do i=1,mesh%nb_trig
        write(1,*) var%fond(i)
      end do

      write(1,*) 'VECTORS u double'
      do i=1,mesh%nb_trig
        write(1,*) var%ux(i),var%uy(i),0.0_pr
      end do
 
      close(1)

    end subroutine sortie_visit

    subroutine sauvegarde_stationnaire(mesh,var)
      implicit none
      !---------------------- in ----------------------------!
      type(maillage)::mesh
      type(modele) ::var
      !------------------------------------------------------!

      !--------------- variables locales --------------------!
      integer :: i
      character(len=300) :: fich
      character(len=100) :: V_courant,h_courant
      !------------------------------------------------------!

      write(V_courant,*) var%V_courant
      write(h_courant,*) var%h_courant

      fich='SORTIES/champs_stationnaire/riviere.V'//trim(adjustl(V_courant))//'_h'//trim(adjustl(h_courant))//'.h'
      open(unit=19,file=fich)
      do i=1,mesh%nb_trig
         write(19,*) var%h(i)
      end do
      close(19)

      fich='SORTIES/champs_stationnaire/riviere.V'//trim(adjustl(V_courant))//'_h'//trim(adjustl(h_courant))//'.ux'
      open(unit=11,file=fich)
      do i=1,mesh%nb_trig
         write(11,*) var%Ux(i)
      end do
      close(11)

      fich='SORTIES/champs_stationnaire/riviere.V'//trim(adjustl(V_courant))//'_h'//trim(adjustl(h_courant))//'.uy'
      open(unit=53,file=fich)
      do i=1,mesh%nb_trig
         write(53,*) var%Uy(i)
      end do
      close(53)

      fich='SORTIES/champs_stationnaire/riviere.V'//trim(adjustl(V_courant))//'_h'//trim(adjustl(h_courant))//'.datas'
      open(unit=67,file=fich)
        write(67,*) var%h_courant
        write(67,*) var%V_courant

      close(67)      

    end subroutine sauvegarde_stationnaire

  subroutine lecture_donnees(var,mesh)
    implicit none 
    type(modele)::var
    type(maillage)::mesh
 

    open(unit=37,file="datas")
      read(37,*) var%tmax
      read(37,*) var%fich
      read(37,*) var%pas_affichage
      read(37,*) mesh%b_ouest , mesh%b_nord ,mesh%b_est,mesh%b_sud
      read(37,*) var%h_courant, var%V_courant
      read(37,*) var%iter_stationnairemax
      read(37,*) var%lecture_stationnaire, var%fichier_stationnaire
      read(37,*) var%h_mascaret, var%V_mascaret
    close(37)
  end subroutine lecture_donnees

    subroutine lecture_stationnaire(mesh,var)
      implicit none
      !---------------------- in ----------------------------!
      type(maillage)::mesh
      type(modele) ::var
      !------------------------------------------------------!

      !--------------- variables locales --------------------!
      integer :: i
      character(len=300) :: fich
      character(len=200) :: prefix
      !------------------------------------------------------!

      write(prefix,*) var%fichier_stationnaire

      fich='SORTIES/champs_stationnaire/'//trim(adjustl(prefix))//'.h'
      open(unit=33,file=fich)
      do i=1,mesh%nb_trig
         read (33,*) var%h(i)
      end do
      close(33)

      fich='SORTIES/champs_stationnaire/'//trim(adjustl(prefix))//'.ux'
      open(unit=11,file=fich)
      do i=1,mesh%nb_trig
         read (11,*) var%Ux(i)
      end do
      close(11)

      fich='SORTIES/champs_stationnaire/'//trim(adjustl(prefix))//'.uy'
      open(unit=53,file=fich)
      do i=1,mesh%nb_trig
         read (53,*) var%Uy(i)
      end do
      close(53)

      fich='SORTIES/champs_stationnaire/'//trim(adjustl(prefix))//'.datas'
      open(unit=67,file=fich)
        read(67,*) var%h_courant
        read(67,*) var%V_courant
      close(67)      
  
      do i=1,mesh%nb_trig
        var%Un(i,1)=var%h(i)
        var%Un(i,2)=var%h(i)*var%Ux(i)
        var%Un(i,3)=var%h(i)*var%Uy(i)
      end do 

    end subroutine lecture_stationnaire

  end module mod_io
