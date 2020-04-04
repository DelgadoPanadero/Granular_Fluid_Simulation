      program granular
c
c     El objetivo de este programa es simular colisiones inelásticas entre esferas duras.
c
      implicit none
                                                             !RECOMENDACIONES
      integer, parameter :: npart=100000
      integer, parameter :: npares=100*npart                !npart*100
      integer, parameter :: nchoques=10000000               !npart*100
      integer, parameter :: nmed=nchoques/100               !muchas medidas
      integer, parameter :: ninter=50                        !no muy alto < 100
      real*8, parameter    :: alfa=0.5d0

      real*8 vx(npart),vy(npart),vz(npart),vmax,h,modulo3
      real*8 suma,modulo2,sumax,sumay,sumaz,temp,maxwell,v2,v4,vv
      real*8 momentox,momentoy,momentoz,fi2,fi4,sumafi2,sumafi4,vsum(3)
      integer i,validos,ni,nj,j
      integer px(0:2*ninter),py(0:2*ninter),pz(0:2*ninter),p(0:2*ninter)
      real*8 vij(3),a,b,pi,mu2,mu4,delta,sonine



      print*,'-----PROGRAMA DE SIMULACION DE CHOQUES INELASTICOS-----'
      print*,'Parametros iniciales:'
      print*,'    numero de particulas:      ' ,npart
      print*,'    numero de choques:         ',nchoques
      print*,'    coeficiente de restitucion:  ',alfa


      pi=3.14159265358979323846d0


      open(14,file='vini.txt')
      open(15,file='vfin.txt')
      open(16,file='histograma.txt')
      open(17,file='maxwell.txt')

!CREAMOS VELOCIDADES ALEATORIAS Y LAS NORMALIZAMOS A UNA ESFERA

      ! CALL init_random_seed()

      do i=1,npart

        modulo2=2.d0

        do while (modulo2.gt.1.d0)

	      call random_number(vx(i))
          call random_number(vy(i))
          call random_number(vz(i))

          vx(i)=vx(i)*2.d0-1.d0
          vy(i)=vy(i)*2.d0-1.d0
          vz(i)=vz(i)*2.d0-1.d0

          modulo2=vx(i)**2+vy(i)**2+vz(i)**2

        end do

        modulo2=dsqrt(modulo2)

        vx(i)=vx(i)/modulo2
        vy(i)=vy(i)/modulo2
        vz(i)=vz(i)/modulo2

        modulo2=vx(i)**2+vy(i)**2+vz(i)**2
        !print*,modulo2

        write(14,*)vx(i), vy(i), vz(i)

      end do
      close(14)

!IMPONEMOS MOMENTO TOTAL NULO

      sumax=0.d0
      sumay=0.d0
      sumaz=0.d0

      do i=1,npart

        sumax=sumax+vx(i)
        sumay=sumay+vy(i)
        sumaz=sumaz+vz(i)

      end do

      momentox=sumax/npart
      momentoy=sumay/npart
      momentoz=sumaz/npart

c       Nuevas velocidades:
      sumax=0.d0
      sumay=0.d0
      sumaz=0.d0

      do i=1,npart

        vx(i)=vx(i)-momentox
        vy(i)=vy(i)-momentoy
        vz(i)=vz(i)-momentoz
        !print*, "vx",i,":",vx(i)
        !print*, "vy",i,":",vy(i)
        !print*, "vz",i,":",vz(i)
        sumax=sumax+vx(i)
        sumay=sumay+vy(i)
        sumaz=sumaz+vz(i)

      end do

      print*,"Momento nuevo:", sumax,sumay,sumaz



c-----------------------------------------------------------------------------
      CALL REESCALAR(vx,vy,vz,npart,temp)
c-----------------------------------------------------------------------------

      suma=0.d0

      do i=1,npart

        suma=vx(i)**2+vy(i)**2+vz(i)**2+suma

      end do


      suma=suma/npart
      temp=2d0*suma/3d0

      print*,'Temperatura tras el reescalado =',temp

!INICIALIZAMOS VARIALBES PARA LOS CHOQUES

      vmax=100.d0
      h=vmax/ninter
      v2=0.d0
      v4=0.d0
      validos=0

      do i=0,2*ninter

        p(i)=0
        px(i)=0
        py(i)=0
        pz(i)=0

      end do

      print*,'INICIANDO CHOQUES DE CALENTAMIENTO:'


c-----------------------------------------------------------------------------
      CALL CHOQUES (vx,vy,vz,npart,nchoques/2,ninter,h,alfa,temp,nmed,
     &validos,p,px,py,pz,v2,v4)
c-----------------------------------------------------------------------------


      print*,'Temperatura tras el calentamiento:',temp


c-----------------------------------------------------------------------------
      CALL REESCALAR(vx,vy,vz,npart,temp)
c-----------------------------------------------------------------------------

      print*,'Reescalamos:',temp


!SACAMOS LAS VARIALBES PARA HACER EL HISTOGRAMA

      vmax=0.d0

      do i=1,npart

        maxwell=4d0*pi*(1d0/(pi*temp**2))**(3d0/2d0)*(i*h)**2d0*
     &exp(-(i*h)**2/temp**2)

        if(modulo2.gt.vmax) then

          vmax=6*modulo2 !chapuza

        end if

      end do

!INICIALIZAMOS VARIALBES PARA LOS CHOQUES

      h=vmax/ninter
      v2=0.d0
      v4=0.d0
      validos=0

      do i=0,2*ninter

        p(i)=0
        px(i)=0
        py(i)=0
        pz(i)=0

      end do



      print*,'BUCLE DE MEDIDA: A CHOCAR!!!!'

      ni=int(nchoques/(npart*100))

      do i=1,ni

c-----------------------------------------------------------------------------
      CALL CHOQUES (vx,vy,vz,npart,npart*100,ninter,h,alfa,temp,nmed,
     &validos,p,px,py,pz,v2,v4)
c-----------------------------------------------------------------------------

      end do

      ni=mod(nchoques,(npart*100))

c-----------------------------------------------------------------------------
      CALL CHOQUES (vx,vy,vz,npart,ni,ninter,h,alfa,temp,nmed,
     &validos,p,px,py,pz,v2,v4)
c-----------------------------------------------------------------------------

      suma=0.d0

      do i=1,npart

        suma=vx(i)**2+vy(i)**2+vz(i)**2+suma

      end do
      suma=suma/npart
      temp=2.d0*suma/3.d0

      print*,'Temperatura tras los choques:',temp


c-----------------------------------------------------------------------------
      CALL REESCALAR(vx,vy,vz,npart,temp)
c-----------------------------------------------------------------------------

      print*,'Reescalamos:',temp
      print*,""


      do i=1,npart

        modulo2=dsqrt(vx(i)**2+vy(i)**2+vz(i)**2)
        write(15,*) vx(i),vy(i),vz(i),modulo2

      end do

      close(15)



      print*,"Tasa de aceptacion: ", float(validos)/nchoques
      print*,"Promedio de v2: ", v2/(nchoques/nmed)
      print*,"Promedio de v4: ", v4/(nchoques/nmed)
      print*,""



      suma=0.d0

      do i=1,npart

        suma=vx(i)**2+vy(i)**2+vz(i)**2+suma

      end do

      suma=suma/npart
      temp=2.d0*suma/3.d0
      temp=dsqrt(temp)

      print*,'Temperatura =',temp
      print*,""

      fi4=0.d0
      fi2=0.d0

      sumafi2=0.d0
      sumafi4=0.d0

!CALCULO DE LOS MOMENTOS DE LA DISTRIBUCIÓN

      do j=1,npares

      ni=1
      nj=ni

        do while (ni.eq.nj)

          call random_number(a)
          call random_number(b)
          ni=floor(npart*a)+1
          nj=floor(npart*b)+1

        end do

        vij(1)=(vx(ni)-vx(nj))
        vij(2)=(vy(ni)-vy(nj))
        vij(3)=(vz(ni)-vz(nj))

        vsum(1)=0.5d0*(vx(ni)+vx(nj))
        vsum(2)=0.5d0*(vy(ni)+vy(nj))
        vsum(3)=0.5d0*(vz(ni)+vz(nj))

        modulo2=vij(1)**2+vij(2)**2+vij(3)**2
        modulo2=dsqrt(modulo2)

        modulo3=vsum(1)**2+vsum(2)**2+vsum(3)**2
        modulo3=dsqrt(modulo3)

        vv=vij(1)*vsum(1)+vij(2)*vsum(2)+vij(3)*vsum(3)

        fi2=modulo2**3*(1d0-alfa**2)*pi/8d0
        sumafi2=fi2+sumafi2

        fi4=modulo2*pi/4d0*(5d0/3d0*(1d0-alfa**2)*modulo2**2*modulo3**2+
     &(1d0-alfa**2)*(2+alfa**2)/12d0*modulo2**4+(3d0-alfa)*(1d0+alfa)*
     &(vv**2-modulo2**2*modulo3**2/3d0))

        sumafi4=fi4+sumafi4

      end do

!---------------------------------------------------------------------------------------
!RESULTADO DE LOS MOMENTOS
!---------------------------------------------------------------------------------------

      print*,"Momento 2: ", sumafi2/npares/temp**3
      print*,"Momento 4: ", sumafi4/npares/temp**5
      mu2=dsqrt(2d0*pi)*(1d0-alfa**2)+
     &dsqrt(2d0*pi)*(1d0-alfa**2)*3d0/16d0*
     &((4d0*v4/(15d0*(nchoques/nmed)))-1d0)
      print*,"Momento 2 teorico: ",mu2
      mu4=mu2*4d0/3d0*v4/(nchoques/nmed)
      print*,"Momento 4 teorico: ",mu4
      print*,""

!---------------------------------------------------------------------------------------
!RESULTADO DEL COEFICIENTE a2
!---------------------------------------------------------------------------------------

      print*,"Valores teoricos de a2=a2(alfa,d):"
      print*,"   a2(p=2): ", 16d0*(1d0-alfa)*(1d0-2d0*alfa**2)/
     &(9d0+24d0*3d0-alfa*(41d0-8d0*3d0)+30d0*(1d0-alfa)*alfa**2)
      print*,"   a2(GT): ", 16d0*(1d0-alfa)*(1d0-2*alfa**2)/
     &(25d0+24d0*3d0-alfa*(57d0-8d0*3d0)-2d0*(1d0-alfa)*alfa**2)
      print*,"   a2(GS): ", 16d0*(1d0-alfa)*(1d0-2*alfa**2)/
     &(401d0-alfa*337d0+190d0*(1d0-alfa)*alfa**2)
      print*,""

!------------------------------RELACIÓN TEÓRICA-----------------------------------------
!La predicción del termostato gaussiano  (GT) y la aproximación para p=2 deben coincidir
!en el intervalo de alfa (0,5-1) y deben tener diferencias significativas en (0-0,3)
!
!El valor de a2(GS) fue una predicción hecha por Goldstein, la cual es erronea. Vamos a
!comprobar que esto es así viendo que da menos de lo medido a2(gs)<0,04
!---------------------------------------------------------------------------------------


      print*,"Valores medidos de a2:"
      print*,"   a2(c^4): ", (4d0*v4/(15d0*(nchoques/nmed)))-1d0

      if (alfa.lt.1d0) then

      print*,"   a2(mu): ", sumafi4/(5d0*sumafi2*temp**2)-1d0

      end if




!ESCRIBIMOS LOS MOMENTOS ALMACENADOS

      do i=0,2*ninter		!No puede haber exponencial para velocidades negativas porque usa modulos

        write(16,*) ((i-ninter)*h),px(i)/((nchoques/nmed)*(npart*h)),
     &py(i)/((nchoques/nmed)*(npart*h)),
     &pz(i)/((nchoques/nmed)*(npart*h))

      end do

c       MAXWELL:
      do i=1,ninter

        maxwell=4d0*pi*(1d0/pi)**(3d0/2d0)*((i+0.5d0)*h)**2d0*exp(-
     &((i+0.5d0)*h)**2)
      
        sonine=float(p(i))/((nchoques/nmed)*(npart*h))   !Aunque hemos usado la variable sonine, esto es  f(c)

      
        delta=((sonine/maxwell)-1d0)/(16d0*(1d0-alfa)*(1d0-2*alfa**2)/
     &(25d0+24d0*3d0-alfa*(57d0-8d0*3d0)-2d0*(1d0-alfa)*alfa**2))

        sonine=0.5d0*((i*h)**2)**2-5.d0/2.d0*(i*h)**2+15d0/8d0


        write(17,*) (i*h),maxwell,p(i)/((nchoques/nmed)*(npart*h)),
     &delta,sonine

      end do

      close(17)
      close(16)

      stop
      end








c-----------------------------------------------------------------------------------------------------------------
c-----------------------------------------------------------------------------------------------------------------

      SUBROUTINE CHOQUES (vx,vy,vz,npart,nchoques,ninter,h,alfa,temp,
     &nmed,validos,p,px,py,pz,v2,v4)

      !ENTRADA:  vx,vy,vz,npart,nchoques,ninter,alfa,h,temp,nmed
      !SALIDA:   vx,vy,vz,validos,p,px,py,pz,v2,v4

      !CALL REESCALAR(vx,vy,vz,npart,temp)
      !CALL HISTOGRAMA(vx,vy,vz,npart,ninter,h,p,px,py,pz)

      integer ni,nj,i,k,validos,npart,nchoques,ninter,nmed
      real*8 vx(npart),vy(npart),vz(npart),vij(3),sigma(3)
      integer px(0:2*ninter),py(0:2*ninter),pz(0:2*ninter),p(0:2*ninter)
      real*8 alfa,temp,wmax,a,b,modulo2,v2,v4,suma2,suma4,h,tirar


      wmax=0.d0

      do k=1,nchoques

      ni=1
      nj=ni

        do while (ni.eq.nj)

          call random_number(a)
          call random_number(b)
          ni=floor(npart*a)+1
          nj=floor(npart*b)+1

        end do

        modulo2=2d0

        do while (modulo2.gt.1d0)

          CALL random_number(sigma(1))
          CALL random_number(sigma(2))
          CALL random_number(sigma(3))

          sigma(1)=sigma(1)*2.d0-1.d0
          sigma(2)=sigma(2)*2.d0-1.d0
          sigma(3)=sigma(3)*2.d0-1.d0





          modulo2= sigma(1)**2+sigma(2)**2+sigma(3)**2

        end do

        modulo2=1/dsqrt(modulo2)

        do i=1,3

          sigma(i)=sigma(i)*modulo2

        end do

        vij(1)=(vx(ni)-vx(nj))
        vij(2)=(vy(ni)-vy(nj))
        vij(3)=(vz(ni)-vz(nj))
        wij=vij(1)*sigma(1)+vij(2)*sigma(2)+vij(3)*sigma(3)

        if(wij.gt.0.d0) then !     Si wij es negativo es que las particulas se están alejando

          if (wmax.lt.abs(wij)) then

            wmax=abs(wij)

          end if

	      call random_number(a)

          if(a.lt.abs(wij)/wmax) then

           !   tirar=vx(ni)*vx(ni)+vy(ni)*vy(ni)+vz(ni)*vz(ni)
            !tirar=tirar+vx(nj)*vx(nj)+vy(nj)*vy(nj)+vz(nj)*vz(nj)
            !print*, tirar

            vx(ni)=vx(ni)-0.5d0*(1.d0+alfa)*wij*sigma(1)
            vy(ni)=vy(ni)-0.5d0*(1.d0+alfa)*wij*sigma(2)
            vz(ni)=vz(ni)-0.5d0*(1.d0+alfa)*wij*sigma(3)
            vx(nj)=vx(nj)+0.5d0*(1.d0+alfa)*wij*sigma(1)
            vy(nj)=vy(nj)+0.5d0*(1.d0+alfa)*wij*sigma(2)
            vz(nj)=vz(nj)+0.5d0*(1.d0+alfa)*wij*sigma(3)

            validos=validos+1
            !tirar=vx(ni)*vx(ni)+vy(ni)*vy(ni)+vz(ni)*vz(ni)
            !tirar=tirar+vx(nj)*vx(nj)+vy(nj)*vy(nj)+vz(nj)*vz(nj)
            !print*, tirar
            !read*

          endif

        endif

!REESCALAMOS LAS VELOCIDADES CON LA TEMPERATURA

        if(mod(k,nmed).eq.0) then

c-------------------------------------------------------------------------------------
          CALL REESCALAR(vx,vy,vz,npart,temp)
c-------------------------------------------------------------------------------------
          suma2=0.d0
          suma4=0.d0

          do i=1,npart

            modulo2=vx(i)**2+vy(i)**2+vz(i)**2
            suma2=modulo2+suma2
            suma4=modulo2**2+suma4

      !      if (sqrt(modulo2).gt.vmax) pause 'Velocidad fuera del
  !   &histograma(|v|>vmax)'

          end do

          suma2=suma2/npart
          suma4=suma4/npart

          v2=v2+suma2
          v4=v4+suma4

          print*,'medida #',k/nmed,'   de',nchoques/nmed
          print*,'v2',suma2,'   v4 ',suma4

c-------------------------------------------------------------------------------------
          CALL HISTOGRAMA(vx,vy,vz,npart,ninter,h,p,px,py,pz)
c-------------------------------------------------------------------------------------

        end if

      end do

      return
      end





c-----------------------------------------------------------------------------------------------------------------
c-----------------------------------------------------------------------------------------------------------------

      SUBROUTINE HISTOGRAMA(vx,vy,vz,npart,ninter,h,p,px,py,pz)
      !ENTRADA: vx,vy,vz,npart,h
      !SALIDA:  p,px,py,pz

      REAL*8 vx(npart),vy(npart),vz(npart),v(npart),h
      INTEGER i,j,npart,ninter
      integer px(0:2*ninter),py(0:2*ninter),pz(0:2*ninter),p(0:2*ninter)


      do i=1,npart

        v(i)=dsqrt(vx(i)**2+vy(i)**2+vz(i)**2)
        j=floor((v(i))/h)

        if ((j.lt.0).OR.(j.gt.2*ninter)) then

          print*,'ERROR, PARTICULA FUERRA DEL HISTOGRAMA p=',j

        else

          p(j)=p(j)+1

        end if

      end do

      do i=1,npart

        j=nint((vx(i))/h)+ninter


        if ((j.lt.0).OR.(j.gt.2*ninter)) then

          print*,'ERROR, PARTICULA FUERRA DEL HISTOGRAMA px=',j

        else

          px(j)=px(j)+1

        end if


      end do

      do i=1,npart

        j=nint((vy(i))/h)+ninter

        if ((j.lt.0).OR.(j.gt.2*ninter)) then

          print*,'ERROR, PARTICULA FUERRA DEL HISTOGRAMA py=',j

        else

          py(j)=py(j)+1

        end if

      end do

      do i=1,npart

        j=nint((vz(i))/h)+ninter

        if ((j.lt.0).OR.(j.gt.2*ninter)) then

          print*,'ERROR, PARTICULA FUERRA DEL HISTOGRAMA pz=',j

        else

          pz(j)=pz(j)+1

        end if


      end do

      return
      end



c-----------------------------------------------------------------------------------------------------------------
c-----------------------------------------------------------------------------------------------------------------

      SUBROUTINE REESCALAR(vx,vy,vz,npart,temp)
      !ENTRADA: vx,vy,vz,npart
      !SALIDA: temp,suma

      REAL*8 suma,vx(npart),vy(npart),vz(npart),temp
      INTEGER i,npart
      REAL*8 sumax,sumay,sumaz,momentox,momentoy,momentoz


!IMPONEMOS MOMENTO TOTAL NULO

      sumax=0.d0
      sumay=0.d0
      sumaz=0.d0

      do i=1,npart

        sumax=sumax+vx(i)
        sumay=sumay+vy(i)
        sumaz=sumaz+vz(i)

      end do

!      print*,"Momento viejo:", sumax,sumay,sumaz

      momentox=sumax/npart
      momentoy=sumay/npart
      momentoz=sumaz/npart


c       Nuevas velocidades:
      sumax=0.d0
      sumay=0.d0
      sumaz=0.d0

      do i=1,npart

        vx(i)=vx(i)-momentox
        vy(i)=vy(i)-momentoy
        vz(i)=vz(i)-momentoz
        !print*, "vx",i,":",vx(i)
        !print*, "vy",i,":",vy(i)
        !print*, "vz",i,":",vz(i)
        sumax=sumax+vx(i)
        sumay=sumay+vy(i)
        sumaz=sumaz+vz(i)

      end do

!      print*,"Momento nuevo:", sumax,sumay,sumaz
       !hasta aqui

      suma=0.d0

      do i=1,npart

        suma=vx(i)**2+vy(i)**2+vz(i)**2+suma

      end do

      suma=suma/npart
      temp=2*suma/3		                        !3 es por la dimensión del espacio
      !print*,'Temperatura =',temp

      temp=1d0/dsqrt(temp)                         !buenas costumbres

      do i=1,npart

        vx(i)=vx(i)*temp
        vy(i)=vy(i)*temp
        vz(i)=vz(i)*temp

      end do

      return
      end
