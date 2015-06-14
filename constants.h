#include<cstdio>
#include<cstdlib>
#include<cmath>
#include<cstring>
#define QCONSTANTS

typedef float real;

#define VERBOSE 0
#define min(x,y) (x<y?x:y)
#define max(x,y) (x>y?x:y)

//;*************************************************************************************************************************
;//Constantes fisicas
//;//*************************************************************************************************************************
int tipo=1;//0:Lvar, 1:Lcon

real Scte = 1.0763e4		;//[MJ m-2 yr-1] incoming solar radiation, aprox. igual que S=341.3 [W m-2] TrenberthEt09, -REVISADO
real Lv = 2.501			;//[MJ kg-1] calor latente de vaporizacion/condensacion (gas/liquido) - REVISADO

//;//*******************************************************************************************************************
;//De la atmosfera
//;//*******************************************************************************************************************
real cpa = 10.			;//[MJ m-2 K-1] calor especifico de la atmosfera
real lapse = 0.0065		;//[K m-1] (Ojo: este valor tambien esta incorporado en la funcion FunLE)
real Ztoai[] = {12.0e3, 12.0e3, 12.0e3, 12.0e3, 10.0e3}	;//[m] [t,g,d,o,cr]

//;//*******************************************************************************************************************
;//De la superficie (land surface, t,g,d)
//;//*******************************************************************************************************************
real cps = 10.						;//[MJ m-2 K-1] calor especifico de la superficie land (equivalente a 2.5 m de agua)
								;//OJO, este valor esta incluido en FunBioticReg
real Wmi[] = {346.2, 86.4, 17.4}		;//[kg m-2] capacidad de campo en [t,g,d]

//;//*******************************************************************************************************************
;//Del oceano
//;//*******************************************************************************************************************
real cpo = 1.64e4	;//[MJ m-2 K-1] calor especifico del oceano (equivalente a 4000 m de agua)
real Moc = 3.7e6		;//[kg m-2] masa por unidad de area, estimado a partir de Rudolf&Rubel05 y TrenberthEt07 (NGK05 usa 3.8e6)

//;//*******************************************************************************************************************
;//De la criosfera
//;//*******************************************************************************************************************
real cpc = 10.		;//[MJ m-2 K-1] calor especifico de la criosfera (equivalente a 2.5 m de agua)
real s1 = 0.435		;//s1 y s2, constantes en la variacion latitudinal de S - ajustada - REVISADA
real s2 = 0.5780335	;//obligatoriamente s2 = (1.-s1)/0.977452 -calculada de s1 - REVISADA
real Mcr = 1.7e6		;//[kg m-2] masa por unidad de area, estimado a partir de Rudolf&Rubel05 (NGK05 usa 0.6e6)
real c1 = 0.509		;//c1 constante en la variacion latitudinal de T - ajustada - REVISADA

//;//*******************************************************************************************************************
;//De la biota
//;//*******************************************************************************************************************
real tada = 0.0	;//[adim] tasa de adaptacion. Valores tipicos: tada=0.0 (sin adaptacion), tada=0.5 (con adaptacion)
real Toptt0 = 295.0	;//[K] temperatura optima inicial. Igual a Tst inicial. 299=273+22, 25 es aprox. la temperatura promedio anual en la Amazonia
real Toptg0 = 290.0	;//[K] 291=273+18, 18 es una temperatura promedio anual razonable en grasslands, Ver por ejemplo CoxEt04, fig2.

//;//*************************************************************************************************************************
;//Constantes numericas
//;//*************************************************************************************************************************
real Tfr = 273.;		//;[K] temperatura de la linea de equilibrio (limite entre ac y ab), NGC05. Tambien esta en FunLE y FunQs
real pinum = 3.141593;	//;pi

real nodata = -9999.;
real amin = 1.e-6		;//area minima en sentido numerico. Por debajo de este valor cualquier area se considera nula (=0.)
				;//equivale a amin*Ae km2 con Ae=5.10e8 km2 (area superficial de la Tierra)
//;//*************************************************************************************************************************
;//Tiempos: integracion, horizonte de modelacion, tamanos de paso
//;//*************************************************************************************************************************
;//escala geologica (millones de anos), en la que cambia el forzamiento solar
real tgini = -3000.		;//[Myr] -3000 tiempo geologico de inicio, 0. es la actualidad
real tgfin = 3000.		;//[Myr] 3000 tiempo geologico de final, 0. es la actualidad (Para llegar hasta 0 poner tgfin = deltatL, o sea un paso adelante
real deltatL = 1000.		;//[Myr] 1000 tamano de paso del tiempo geologico, de la Luminosidad solar (relevante que sea pequeno cuando tada <> 0.0)

;//escala miles de anos, durante la cual el forzamiento solar es aprox. constante, la de buscar estados de equilibrio
real tmod = 5000.		;//[yr] 5000 horizonte de modelacion para cada forzamiento solar constante. Es el tiempo durante el cual se integra para encontrar los estados de equilibrio correspondientes a algun forzamiento solar fijo.
real deltatINT = 1.		;//[yr] tamano de paso del tiempo de modelacion, INTegracion

real tequ = 2000.		;//[yr] 2000 numero de anos que para un L dado describen el estado de equilibrio. Menoro igual que tmod.

;//escala anual
real deltatRK4=1./(365.*2.)	;//[yr] tamano de paso del metodo RK4. 1/2 dia, aprox 0.00137 yr (NGC05 usan 0.001)

;//;//iLparciales = [5,10,15,20,25,30,35]	;//valores de iL donde se escriben los archivos parciales. Ajustar dependiendo de deltaL
real iLparciales[] = {3,6,9,12};


//;//*************************************************************************************************************************
;//Definicion de fracciones de area dinamicas, DAF, en el siguiente orden: [t,g,d,o,cr]
//;//*************************************************************************************************************************
char DAFs[][3] = {"t","g","d","o","cr","ab","ac"}		;//cr = ab+ac

//;//*******************************************************************************************************************
;//Condiciones iniciales. Escritas para cada DAF en el siguiente orden: [t,g,d,o,ac,ab].
;//Estan basadas en promedios globales actuales (ver hoja de excel "ParametrosEarthDAFM")
//;//*******************************************************************************************************************
real CIareai[] = {0.087, 0.087, 0.087, 0.708, 0.031}	;//[adim] DAF iniciales. La suma de todas tiene que ser 1.0. - REVISADAS
real CITsi[] = {Toptt0, Toptg0, 282.0, 289.1, 247.0,}	;//[K] basadas en promedios globales actuales (0 y cr) - REVISADAS

real CITai[5] ;//[K] +10m para inducir flujo inicial de calor sensible
real CIWsi[3] ;//[kg m-2], humedad del suelo en t,g,d. Carece de sentido en o y cr

real CIWai[] = {1037., 1015., 1004., 1025., 1001.} 	;//[kg m-2] agua precipitable total en la atmosfera, Wa=Wvap+Wpre
							;//donde Wvap: remanente en la atmosfera, Wpre: excedente que se precipita
real CIBRi[] = {0.,0.}	;// OJO, ESTO PUEDE CAMBIAR SI SE REINICIA EN PROGRAMA CON UNA CI DE TEMPERATURA DIFERENTE A Ts=Topt
;//;//CIWvapi = [37.0788, 15.3501, 4.34366, 22.7554, 0.628646]		;//[kg m-2], agua precipitable, remanente en la atmosfera
;//;//CIWprei = [999.921, 999.650, 999.656, 1002.24, 1000.37]	;//[kg m-2], agua precipitable, excedente que se vuelve precipitacion
				                                ;//el area se entiende como area-1 (adimensional)
//;//*************************************************************************************************************************
