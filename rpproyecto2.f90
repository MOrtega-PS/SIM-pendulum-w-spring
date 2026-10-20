MODULE constantes 
IMPLICIT NONE
SAVE
INTEGER,PARAMETER::ki=SELECTED_REAL_KIND(30,3000)
END MODULE constantes

Module varglob
  IMPLICIT NONE
  real(SELECTED_REAL_KIND(30,3000))::g,m,r0,kr
end module varglob


program proyectorp
  use constantes
  use varglob
  implicit none
  integer :: i,n,aux2
  real(ki), allocatable:: rs1(:), rs2(:), k(:,:),os1(:),os2(:)
  real(ki),allocatable:: xs(:),ys(:)
  real(ki):: a,b,h,ts,aux,En

  
!!!!diccionario de variables1!!!!!!!!!!1!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !a,b limites del intervalo
  !h tamaño de paso
  !n numero de subintervalos
  !i variable de conteo
  !rs1 Puntos de la solucion para la velocidad radial
  !rs2 Puntos de la solucion para la posicion radial
  !os1 Puntos de la solucion para el angulo
  !os2 Puntos de la solucion para la velocidad angular
  !k Valores de la k para runge kuta
  !aux variable auxiliar en caso de necesitarla
  !xs Posicion en el eje x
  !ys Posicion en el eje y
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!se le pide al usuario los val iniciales!!!!!!!!!!!!!!!
  !print*,"Establezca la a y despues la b"
  !read*,a,b
  !print*,"Establezca el numero de subintervalos"
  !read*,n
  a = 0
  b = 30
  n = 100000
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!! se inicializan variables y arrays!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  h = (b-a)/n!se calcula el paso
  allocate(rs1(n+1)) !vr(t) velocidad radial
  allocate(rs2(n+1)) !r(t) posicion radial
  allocate(os1(n+1)) !w(t) velocidad angular
  allocate(os2(n+1)) !0(t) angulo
  allocate(xs(n+1))  !posicion en x
  allocate(ys(n+1))  !posicion en y
  allocate(k(4,4))
  !!!!!se preguntan los val iniciales de las Ecuaciones
  print*, "Establezca la posicion radial y angular inicial (en ese orden)"
  read*,aux
  rs2(1) = aux
  r0 = aux !posicion angular inicial
  read*,aux
  os2(1) = aux
  print*, "Establezca la velocidad radial y angular inicial (en ese orden)"
  read*,aux
  rs1(1) = aux
  read*,aux
  os1(1) = aux
  !print*, "Establezca la masa y la constante k del resorte (en ese orden)"
  !read*,m,kr
  g = 9.81 !se establece la gravedad
  !!!!!valores automaticos para pruebas rapidas
  !rs1(1) = 0.5
  !os1(1) = 0.5
  !rs2(1) = 2
  !os2(1) = 1
  !r0 = 0.5
  m = 2
  k = 3
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !se abre el archivo donde se guarda la solucion
  open(unit=2,file ="data_polares",status="replace")
  open(unit=3,file ="data_cartesianas",status="replace")
  
  !!!!!!!!!!!!!!!!!!!!metodo runge kuta para y !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do i =1,n+1
     ts=a+i*h !guardas la t que se "usa"
     !!!!!!!!!!!!!!!!!!!!!!!k1!!!!!!!!!!!!!!!!!!!!!!!!!
     k(1,1) = h*r1(ts,rs1(i), rs2(i), os2(i), os1(i))
     k(1,2) = h*r2(ts,rs1(i), rs2(i))
     k(1,3) = h*o1(ts,rs1(i), rs2(i), os2(i), os1(i))
     k(1,4) = h*o2(ts,os1(i), os2(i))
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     !!!!!!!!!!!!!!!!!!!!!!!k2!!!!!!!!!!!!!!!!!!!!!!!!!
     k(2,1) = h*r1(ts + 0.5*h, rs1(i) + 0.5*k(1,1), rs2(i) + 0.5*k(1,2), &
          & + os2(i) + 0.5*k(1,3), os1(i) + 0.5*k(1,4))
     k(2,2) = h*r2(ts + 0.5*h, rs1(i) + 0.5*k(1,1), rs2(i) + 0.5*k(1,2))
     k(2,3) = h*o1(ts + 0.5*h, rs1(i) + 0.5*k(1,1), rs2(i) + 0.5*k(1,2), &
          & + os2(i) + 0.5*k(1,3), os1(i) + 0.5*k(1,4))
     k(2,4) = h*o2(ts + 0.5*h, os1(i) + 0.5*k(1,3), os2(i) + 0.5*k(1,4))
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     !!!!!!!!!!!!!!!!!!!!!!!k3!!!!!!!!!!!!!!!!!!!!!!!!!
     k(3,1) = h*r1(ts + 0.5*h, rs1(i) + 0.5*k(2,1), rs2(i) + 0.5*k(2,2), &
          & + os2(i) + 0.5*k(2,3), os1(i) + 0.5*k(2,4))
     k(3,2) = h*r2(ts + 0.5*h, rs1(i) + 0.5*k(2,1), rs2(i) + 0.5*k(2,2))
     k(3,3) = h*o1(ts + 0.5*h, rs1(i) + 0.5*k(2,1), rs2(i) + 0.5*k(2,2), &
          & + os2(i) + 0.5*k(2,3), os1(i) + 0.5*k(2,4))
     k(3,4) = h*o2(ts + 0.5*h, os1(i) + 0.5*k(2,3), os2(i) + 0.5*k(2,4))
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     !!!!!!!!!!!!!!!!!!!!!!!k4!!!!!!!!!!!!!!!!!!!!!!!!!
     k(4,1) = h*r1(ts + h, rs1(i) + k(3,1), rs2(i) + k(3,2), &
          & + os2(i) + k(3,3), os1(i) + k(3,4))
     k(4,2) = h*r2(ts + h, rs1(i) + k(3,1), rs2(i) + k(3,2))
     k(4,3) = h*o1(ts + h, rs1(i) + k(3,1), rs2(i) + k(3,2), &
          & + os2(i) + k(3,3), os1(i) + k(3,4))
     k(4,4) = h*o2(ts + h, os1(i) + k(3,3), os2(i) + k(3,4))
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     !!!!!!!!!!!!!!!!!!!!!!!!calculo de la y siguiente!!!!!!!!!!!!!!!!!!!!!!!
     rs1(i+1) = rs1(i)+(1.0/6.0)*(k(1,1) + 2*k(2,1) + 2*k(3,1)+k(4,1))
     rs2(i+1) = rs2(i)+(1.0/6.0)*(k(1,2) + 2*k(2,2) + 2*k(3,2)+k(4,2))
     os1(i+1) = os1(i)+(1.0/6.0)*(k(1,3) + 2*k(2,3) + 2*k(3,3)+k(4,3))
     os2(i+1) = os2(i)+(1.0/6.0)*(k(1,4) + 2*k(2,4) + 2*k(3,4)+k(4,4))
     !rs1(i+1) = r1(ts+h,0.0_ki,0.0_ki,0.0_ki,0.0_ki)
     !rs2(i+1) = r2(ts+h,0.0_ki,0.0_ki)
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  enddo
  !!!!!!!!!!!Fin del metodo de runge kuta para y!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!! Conversion de polares a cartesianas!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do i =0,n+1
     xs(i) = rs2(i)*sin(os2(i))
     ys(i) = -rs2(i)*cos(os2(i))
  end do
  
  !!!!!!!!!!! Fin de conversion a cartesianas!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !!!!!!!!!!!!!!escritura de las tabla cartesianas!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ts = a
  do i =0,1000
     ts = ts+h*(n/1000)
     write(3,*) ts, xs((n/1000)*i+1),ys((n/1000)*i+1)
  enddo
  !!!!!!!!!!!!!fin de la escritura de las tablas cartesianas!!!!!!!!!!!!!!!!!
  
!!!!!!!!!!!!!!escritura de las tablas!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ts = a
  do i =0,1000
     ts = ts+h*(n/1000)
     aux2 = (n/1000)*i+1
     En = 0.5*m*(rs2(aux2))**2+0.5*m*((rs1(aux2))**2)*(os2(aux2))**2-m*g*ys(aux2)
     write(2,*) ts, rs1(aux2), rs2(aux2), os1(aux2), os2(aux2), m*rs2(aux2), ((rs1(aux2))**2)*m*(os2(aux2)), en
  enddo
  !!!!!!!!!!!!!fin de la escritura de las tablas!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  
  !!!!!!!!!!!!!!!!!!!!!!!!!graficacion!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!grafica 1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  open(unit= 5, file = "datos.plt", status = "unknown")
  write(5,*)'m = "./data_polares"'!indicaciones del archivo de datos
  write(5,*) 'set title "Trayectorias"'
  write(5,*)'plot m using 1:5 with lines title "pos angular"'
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!grafica 2!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  open(unit= 6, file = "datos2.plt", status = "unknown")
  write(6,*)'m = "./data_polares"'!indicaciones del archivo de datos
  write(6,*) 'set title "Trayectorias"'
  write(6,*)'plot m using 1:3 with lines title "pos radial"'
  call system("gnuplot -p datos.plt")!llamar al sistema para usar gnuplot
  call system("gnuplot -p datos2.plt")
  !!!!!!!!!!!!!!!!!!!!!!!!Fin de graficacion!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
contains
  function r1(t,vr,r,o,w)
    use constantes
    use varglob
    implicit none
    real(ki):: r1,t,vr,r,w,o
    g = 9.81_ki
    !r1 =  r*w**2 + g*cos(o) - kr*(r-r0)/m
    !r1 = 0 !caso trivial, caso vel constante
    !r1= 2*t !funcion creciente vel= 2
    r1 = -0.25*sin(11*t)
  end function r1

  function r2(t,vr,r)
    use constantes
    implicit none
    real(ki)::t,vr,r2,r,aux
    aux = 11
    !r2 = vr
    !r2= 2 !caso de velocidad constante r(t) = .5*t
    !r2 = 0 !caso trivial
    r2 = aux*0.25*cos(aux*t)
  end function r2

  function o1(t,vr,r,o,w)
    use constantes
    use varglob
    implicit none
    real(ki):: r,vr,o,w,o1,t
    o1 =  (-2*vr*w)/r - (g*sin(o))/r
  end function o1

  function o2(t,w,o)
    use constantes
    implicit none
    real(ki)::o2,t,w,o
    o2 = w
  end function o2
  
end program proyectorp
