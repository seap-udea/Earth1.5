#include<constants.h>
#include<earth15-alpha.cpp>

using namespace std;
int main(void)
{
  //INITIALIZATION
  for(int i=0;i<5;i++) CITai[i]=CITsi[i]-lapse*(Ztoai[i]/2.+10.);
  for(int i=0;i<3;i++) CIWsi[i]=Wmi[i]/2;
  printVector("Tai = ",CITai,5,1);
  printVector("Tai = ",CIWsi,5,1);

  real areai_t1[5],at_t1,ag_t1,ao_t1,acr_t1,Tsi_t1[5],Tai_t1[5],Wsi_t1[3],Wai_t1[5],criosaT_t1[6];
  int t,n;
  real acr=CIareai[4];
  real crios[2];
  real Tsac,Taac,Tsab,Taab,inteab,inteac;
  real aac,aab;
  real tetaac;
  if(acr>0){
    real Tscr=CITsi[4];
    real Tacr=CITai[4];
    real tetacr=acos(1.-acr);
    real intecr=FunInte(0.,tetacr);
    printf("teta,integral = %e,%e\n",tetacr,intecr);
    real c2 = (1.-c1)*(1.-cos(tetacr))/intecr;
    FunCrios(Tscr,acr,tetacr,c1,c2,amin,crios);//;tiene incorporada la condicion de amin
    aac = crios[0];
    aab = crios[1];

    if(aac>0){
      tetaac = acos(1.-aac);
      inteac = FunInte(0.,tetaac);
      //;ec (9) modelo criosfera en papel
      Tsac = Tscr*c1 + Tscr*c2*inteac/(1.-cos(tetaac));
      //;ec (9) modelo criosfera en papel
      Taac = Tacr*c1 + Tacr*c2*inteac/(1. - cos(tetaac));	
    }else{
      inteac = 0.;		//;necesaria para zona de ablacion
      Tsac = nodata;	//;no existe zona de acumulacion
      Taac = nodata;
    }

    if(aab>0.){
      //;si aac gt 0. tetaac ya fue calculada	
      if(aac==0){tetaac = 0.;}
      //;da lo mismo que inteab = FunInte(tetaac,tetacr)
      inteab = intecr - inteac;		
      //;ec (9) modelo criosfera en papel
      Tsab = Tscr*c1 + Tscr*c2*inteab/(cos(tetaac) - cos(tetacr));
      //;ec (9) modelo criosfera en papel
      Taab = Tacr*c1 + Tacr*c2*inteab/(cos(tetaac) - cos(tetacr));
    }else{
      Tsab = nodata;	//;no existe zona de ablacion
      Taab = nodata;
    }

  }else{
    aac = nodata;
    aab = nodata;
    Tsac = nodata;
    Taac = nodata;
    Tsab = nodata;
    Taab = nodata;
  }
  if(VERBOSE) printf("tetaac = %e\n",tetaac);
  if(VERBOSE) printf("aab = %lf, aac = %lf\n",aab,aac);
  if(VERBOSE) printf("Tsab = %lf, Taab = %lf, Tsac = %lf, Taac = %lf\n",Tsab,Taab,Tsac,Taac);

  real CIDafmTsi[5];//;ajuste BES
  real CIDafmTai[5];//;ajuste BEA
  real CIDafmWsi[3];//;ajuste BAS (solo contiene t,g,d)
  real CIDafmWai[5];//;ajuste BAA

  real DafmTsi[5];//;ajuste BES
  real* DafmTsiret;
  real DafmTai[5];//;ajuste BEA
  real* DafmTairet;
  real DafmWsi[3];//;ajuste BAS (solo contiene t,g,d)
  real* DafmWsiret;
  real DafmWai[5];//;ajuste BAA
  real* DafmWairet;
  
  real areaiDafm[6];
  real dareaidtDafm[6];
  real TaiDafm[6];
  real TsiDafm[6];
  real cpTsiDafm[6];
  real cpTaiDafm[6];
  real WaiDafm[6];

  real CIcriosaT[6]={aab, aac, Tsab, Tsac, Taab, Taac};
  real criosaT[6];

  //;calor especifico de la superficie  i=t,g,d,o,cr
  real cpsi[5] = {cps, cps, cps, cpo, cpc};
  
  //;numero de iteraciones del forzamiento solar (Luminosidad)
  int niterL = round((tgfin-tgini)/deltatL);
  //;numero de iteraciones de INTegracion
  int niterINT = round(tmod/deltatINT);
  //;numero de iteraciones en los RK4
  int niterRK4 = round(deltatINT/deltatRK4);				

  if(VERBOSE) printf("%f %f\n",deltatINT,deltatRK4);
  if(VERBOSE) printf("%d %d %d\n",niterL,niterINT,niterRK4);

  //;matrices para almacenar resultados con forzamiento solar variables
  real** vecL=fltarr(1,niterL+1);
  
  real** matAreaMED = fltarr(7, niterL+1);
  real** matTsiMED = fltarr(7, niterL+1);
  real** matTaiMED = fltarr(7, niterL+1);
  real** matWsiMED = fltarr(3, niterL+1);
  real** matWaiMED = fltarr(5, niterL+1);
  real** matWviMED = fltarr(7, niterL+1);
  real** matBRiMED = fltarr(2, niterL+1);
  
  real** matAreaMAX = fltarr(7, niterL+1);
  real** matTsiMAX = fltarr(7, niterL+1);
  real** matTaiMAX = fltarr(7, niterL+1);
  real** matWsiMAX = fltarr(3, niterL+1);
  real** matWaiMAX = fltarr(5, niterL+1);
  real** matWviMAX = fltarr(7, niterL+1);
  real** matBRiMAX = fltarr(2, niterL+1);
  
  real** matAreaMIN = fltarr(7, niterL+1);
  real** matTsiMIN = fltarr(7, niterL+1);
  real** matTaiMIN = fltarr(7, niterL+1);
  real** matWsiMIN = fltarr(3, niterL+1);
  real** matWaiMIN = fltarr(5, niterL+1);
  real** matWviMIN = fltarr(7, niterL+1);
  real** matBRiMIN = fltarr(2, niterL+1);

  //;Para almacenar series de variables de estado. Cada columna para una DAF en el siguiente orden: [t,g,d,o,cr]
  real** vecArea = fltarr(7,tmod+1);		//;7 para incluir ac y ab, [t,g,d,o,cr,ac,ab]. +1 porque comienza en t=0
  real** vecTsi = fltarr(7,tmod+1);		//;[t,g,d,o,cr,ac,ab]
  real** vecTai = fltarr(7,tmod+1);		//;[t,g,d,o,cr,ac,ab]
  real** vecWsi = fltarr(3,tmod+1);		//;solo contienen t,g,d
  real** vecWai = fltarr(5,tmod+1);		//;[t,g,d,o,cr]
  real** vecBRi = fltarr(2,tmod+1);		//;[t,g]


  real at,Tst,Wst,Wmt,Betat,tmort;
  real ag,Tsg,Wsg,Wmg,Betag,tmorg;
  
  int iL;
  real tg,L,Toptt,Toptg,Tsetq,Tsgeq;

  real areai[5];		
  real Tsi[5];
  real Tai[5];
  real Wsi[5];
  real Wai[5];
  real Tsteq;

  real tetacr,inte,Sol,Scr;

  real Wvi[5];
  real Preci[5];
  real ETRi[5];
  real Qsi[5];
  real SHsi[5];
  real LHsi[5];
  real LHai[5];
  real Rsi[5];
  real Rai[5];
  real ROCsi[5];
  real Cloudi[5];
  real albci[5];
  real RHi[5];
  real BRi[5];
  real BRalbi[5];
  real* FEai;

  real Wv = nodata;
  real Prec = nodata;
  real ETR = nodata;
  real SHs = nodata;
  real LHs = nodata;
  real LHa = nodata;
  real Rs = nodata;
  real Ra = nodata;
  real ac = nodata;
  real albc = nodata;
  real RH = nodata;

  real Wsat,qv0,RH0,Ta0;
  real Wacr,Wvcr,Ztoacr,fac1,fac2,facab,facac,Wvab,Wvac;
  real ROCs;

  real las,k1at,k1ag,Topt,bs;
  real *BioReg;
  real BRalb,BRetr,deltaETR,BR;

  real* Qsiret;
  real ao;
  real k1ao,k1acr,wmi;

  real* k1Tsi = fltarr1(5);
  real* k1Tai = fltarr1(5);
  real* k1Wai = fltarr1(5);
  real* k1Wsi = fltarr1(3);

  real* k2Tsi = fltarr1(5);
  real* k2Tai = fltarr1(5);
  real* k2Wai = fltarr1(5);
  real* k2Wsi = fltarr1(3);

  real* k3Tsi = fltarr1(5);
  real* k3Tai = fltarr1(5);
  real* k3Wai = fltarr1(5);
  real* k3Wsi = fltarr1(3);

  real* k4Tsi = fltarr1(5);
  real* k4Tai = fltarr1(5);
  real* k4Wai = fltarr1(5);
  real* k4Wsi = fltarr1(3);

  real* dTsidtRK4 = fltarr1(5);
  real* dTaidtRK4 = fltarr1(5);
  real* dWaidtRK4 = fltarr1(5);
  real* dareaidtRK4 = fltarr1(4);
  real* dWsidtRK4 = fltarr1(3);

  //DECLARACION
  real ad,Tscr,Tacr,intecr,c2;
  real k2at,k2ag,k2ao,k2acr;
  real k3at,k3ag,k3ao,k3acr;
  real k4at,k4ag,k4ao,k4acr;

  real dadt1[4],dadt2[2];

  FILE *fl=fopen("results.dat","w");
  //fprintf(fl,"# 0:1 1:area t		      2:area g                3:area d                4:area o                5:area cr               6:Ts t		      8:Ts g                  9:Ts d                  10:Ts o                 11:Ts cr                12:Ta t		      13:Ta g                 14:Ta d                 15:Ta o                 16:Ta cr                17:Wa t		      18:Wa g                 19:Wa d                 20:Wa o                 21:Wa cr                \n");

  for(iL=0;iL<=niterL;iL++){//;Loop del forzamiento solar

    tg = (tgini + deltatL*iL)/1000.;		//;[Gyr] tiempo geologico
    L = 1./(1.-(0.38*tg/4.55));			//;[adim] luminosidad,
                                                // Caldeira and Kasting (1992)

    if(VERBOSE) printf("L = %.17e\n",L);

    //;******************
    if(tada>0.){//then begin
		//;En este bloquecito se introduce la adaptacion biologica al modificar Toptt. Esto se complementa con FunBeta
		//;temperatura optima [K] de reproduccion de las especies biologicas, sujeta a adaptacion
      if(iL==0){//then begin
	Toptt = Toptt0;
	Toptg = Toptg0;
      }else{
	Toptt = Toptt + tada*(Tsteq - Toptt);	//;Tsteq y Tsgeq se calculan al final del periodo de modelacion (loopINT)
	Toptt = Toptt + tada*(Tsgeq - Toptg);
      }
    }else{
      Toptt = Toptt0;
      Toptg = Toptg0;
    }

    // THEY SHOULD BE PUT IN ZERO
    /*
    double vecArea = fltarr(7,tmod+1);		//;7 para incluir ac y ab, [t,g,d,o,cr,ac,ab]. +1 porque comienza en t=0
    double vecTsi = fltarr(7,tmod+1);		//;[t,g,d,o,cr,ac,ab]
    double vecTai = fltarr(7,tmod+1);		//;[t,g,d,o,cr,ac,ab]
    double vecWsi = fltarr(3,tmod+1);		//;solo contienen t,g,d
    double vecWai = fltarr(5,tmod+1);		//;[t,g,d,o,cr]
    double vecBRi = fltarr(2,tmod+1);		//;[t,g]
    */
    
    //;Condiciones iniciales: son vectores con los valores para: [t,g,d,o,cr]
    copyVector(areai,CIareai,5);	//;en el procesamiento no se incluyen ac y ab sino la criosfera como un todo y ac y ab se deducen
    copyVector(Tsi,CITsi,5);		//;en los balances la criosfera se trata como un todo
    copyVector(Tai,CITai,5);	
    copyVector(Wsi,CIWsi,3);
    copyVector(Wai,CIWai,5);
    //;contiene [aab, aac, Tsab, Tsac, Taab, Taac] que se requieren en algunas funciones
    copyVector(criosaT,CIcriosaT,6);		
    aab = criosaT[0];
    aac = criosaT[1];
    Tsab = criosaT[2];
    Tsac = criosaT[3];
    Taab = criosaT[4];
    Taac = criosaT[5];

    //;Ajustes DAFM iniciales (todos iguales a cero porque al inicio no ha habido cambios en las areas)
    //BES: Balance Energia Superficial
    zeroVector(CIDafmTsi,5);//;ajuste BES
    zeroVector(CIDafmTai,5);//;ajuste BEA
    zeroVector(CIDafmWsi,3);//;ajuste BAS (solo contiene t,g,d)
    zeroVector(CIDafmWai,5);//;ajuste BAA

    //;BLOQUE1-inicio: Calculando variables auxiliares iniciales
    //;En este bloque no hacen falta los ajuste DAFM por tratarse de la condicion inicial donde todavia no han cambiado
    //;las areas
    
    //;constante solar regional
    acr = areai[4];
    tetacr = acos(1.-acr);
    inte = FunInte(tetacr,pinum/2.);		//;integral

    printVector("criosaT:",criosaT,6);

    //Scte: Mj / m^2 yr^2 = 341.4 W/m^2
    Sol = Scte*s1 + (Scte*s2 / cos(tetacr) )*inte;	//;constante solar para ocean-land
    if(acr>0.){
      Scr = (Scte-Sol*(1.-acr))/acr;
    }else{
      Scr = 0.;
    }
    
    if(VERBOSE) printf("acr = %.17e\n",acr);

    if(VERBOSE) printf("tetacr = %.17e, Sol = %.17e, Scr = %.17e\n",tetacr,Sol,Scr);

    //;flujos horizontales de energia
    FEai = FunFEa(areai, Tai);		//;usa los vectores completos
    printVector("FEai:",FEai,5,1);

    //;devuelve DeltaEai = vector en orden ['t','g','d','o','cr']
    //;revisa conservacion de la energia, si no se cumple devuelve mensaje de error
    //;tiene incorporada la condicion: si area=0 entonces todo=nodata

    char DAF[10];
    int i;
    real S,area,Ts,Ta,Wa,Ztoa;
    real *atmos,*clouds,*radnet,*laten,*Qsi;//7 atmospheric parameters
    real *FWai;
    real ac,albc,Tc,emc;

    printVector("Areai = ",areai,5);

    for(i=0;i<=4;i++){

      strcpy(DAF,DAFs[i]); //;DAFs= ['t','g','d','o','cr']
      //strcpy(DAF,DAFs[4]); //;DAFs= ['t','g','d','o','cr']
      
      if(strcmp(DAF,"cr")==0){
	S=Scr;
      }else{
	S=Sol;
      }
      area=areai[i];
      if(VERBOSE) printf("MAIN:Area = %f\n",area);

      if(area>0.){
	Ts = Tsi[i];
	Ta = Tai[i];
	Wa = Wai[i];
	Ztoa = Ztoai[i];
	
	//;atmosfera, resultados = [Wv, Wpre, Wsat, RH, qv0, RH0, Ta0]
	if(VERBOSE) 
	  printf("MAIN: FunAtm %s %e %e %e %e\n",DAFs[i],Ta,Wa,lapse,Ztoa);
	
	atmos=FunAtm(Ta, Wa, lapse, Ztoa, DAF);
	printVector("MAIN:Atmos:",atmos,7);
	
	Wv = atmos[0];		//;[kg m-2] agua precipitable que permanece en la atmosfera
	Prec = atmos[1];		//;[kg m-2 yr-1] excedente de Wa que se vuelve precipitacion. Prec=Wpre/1yr
	Wsat = atmos[2];		//;[kg m-2] capacidad maxima de almacenamiento de agua en la atmosfera
	RH = atmos[3];		//;humedad relativa media en la columna
	qv0 = atmos[4];		//;humedad especifica del aire en la superficie
	RH0 = atmos[5];		//;humedad relativa del aire en la superficie
	Ta0 = atmos[6];		//;temperatura del aire en la superficie
	
	//;nubes, resultados = [ac, albc, Tc, emc]
	if(VERBOSE) 
	  printf("MAIN: FunCloud %s %e %e %e %e %e %e\n",DAFs[i],Prec, RH, Ztoa, lapse, qv0, Ta0);
	clouds = FunCloud(Prec, RH, Ztoa, lapse, qv0, Ta0, DAF);

	printVector("Clouds:",clouds,4);

	ac = clouds[0];		//;area total de nubes
	albc = clouds[1];	//;albedo promedio de las nubes
	Tc = clouds[2];		//;temperatura promedio de las nubes
	emc = clouds[3];		//;emisividad promedio de las nubes
	
	//;Net surface radiation
	if(VERBOSE) printf("S=%e,Scr=%e,L=%e,Ts=%e,Ta=%e\n",S,Scr,L,Ts,Ta);
	radnet = FunRs(S,L,Ts,Ta,RH,clouds,criosaT,DAF);		//;Rs en [MJ m-2 yr-1]

	Rs = radnet[0];
	ROCs = radnet[1];
	  
	//;Sensible heat - superficie
	SHs = FunSH(Ts, Ta0, DAF);				//;[MJ m-2 yr-1]
	if(VERBOSE) printf("SHs = %e\n",SHs);

	//;Calor latente, evaporacion, transpiracion, sublimacion - superficie	
	laten = FunLE(Ta0, Ztoa, Rs*1.e6, Prec, qv0, RH0, criosaT, DAF);	//;Rs en [J m-2 yr-1]
	//;//;//;laten = FunLE(Ts,lapse,Rs*1.e6,Prec,qv10,criosaT,DAF);	//;Rs en [J m-2 yr-1]
	LHs = laten[0];							//;[MJ m-2 yr-1]
	ETR = laten[1];							//;[kg m-2 yr-1]
	if(VERBOSE) printf("LHs = %e, ETR = %e\n",LHs,ETR);
	    
	//;radiacion neta en la atmosfera
	Ra = FunRa(S, L, Ts, Ta, RH, clouds, criosaT, DAF);		//;[MJ m-2 yr-1]
	if(VERBOSE) printf("Ra = %e\n",Ra);
	
	//;calor latente en la atmosfera			
	LHa = FunLHa(Prec);						//;[MJ m-2 yr-1]
	if(VERBOSE) printf("LHa = %e\n",LHa);

      }

      if(VERBOSE or 0) printf("%s %e %e %e %e %e %e %e %e %e %e %e\n",
	     DAFs[i],
	     Wv,
	     Prec,
	     ETR,
	     SHs,
	     LHs,
	     LHa,
	     Rs,
	     Ra,
	     ac,
	     albc,
	     RH);

      Wvi[i] = Wv;
      Preci[i] = Prec;
      ETRi[i] = ETR;
      SHsi[i] = SHs;
      LHsi[i] = LHs;
      LHai[i] = LHa;
      Rsi[i] = Rs;
      Rai[i] = Ra;
      Cloudi[i] = ac;
      albci[i] = albc;
      RHi[i] = RH;
      //break;
    } //End for DAFs
    
    //;escorrentia superficial, melting
    //;[km m-2 yr] devuelve un vector
    Qsi = FunQs(Preci, ETRi, Wsi, Wmi, areai, criosaT);
    printVector("Qsi:",Qsi,5,1);
    
    //;flujo horizontal de humedad atmosferica
    FWai = FunFWa(areai,Tai,RHi);	//;[km m-2 yr] devuelve un vector
    printVector("FWai:",FWai,5,1);

    //;*********************************************************************************************************
    //;particion del vapor de agua en la criosfera (acr diferentes de 0. en la condicion inicial)
    if(acr>0.){
      Wacr = Wai[4];
      Wvcr = Wvi[4];					//;Wacr = Wvcr+Wprecr (conservacion)
      Ztoacr = Ztoai[4];
      atmos = FunAtm(Taab, Wacr, lapse, Ztoacr, DAF);	//;atmosfera, resultados = [Wv, Wpre, Wsat, RH, qv0, RH0, Ta0]
      fac1 = atmos[0];
      atmos = FunAtm(Taac, Wacr, lapse, Ztoacr, DAF);	//;atmosfera, resultados = [Wv, Wpre, Wsat, RH, qv0, RH0, Ta0]
      fac2 = atmos[0];
      facab = fac1/(fac1+fac2);
      facac = fac2/(fac1+fac2);
      Wvab = facab*Wvcr;
      Wvac = facac*Wvcr;
    }else{
      fprintf(stdout,"acr=0.0 en la condicion inicial\n");
      Wvab=nodata;
      Wvac=nodata;
    }

    //;Almacenando variables de estado iniciales
    fillRow(vecArea,0,7,(float[7]){areai[0],areai[1],areai[2],areai[3],areai[4],aab,aac});
    fillRow(vecTsi,0,7,(float[7]){Tsi[0],Tsi[1],Tsi[2],Tsi[3],Tsi[4],Tsab,Tsac});
    fillRow(vecTai,0,7,(float[7]){Tai[0],Tai[1],Tai[2],Tai[3],Tai[4],Taab,Taac});
    fillRow(vecWsi,0,3,Wsi);
    fillRow(vecWai,0,3,Wai);
    //;condicion inicial [0.,0.] si se inicia con Ts=Topt. Cambia si se cae una corrida y se reinicia.
    fillRow(vecBRi,0,3,CIBRi);	

    //;BLOQUE1-final: Calculando variables auxiliares iniciales
    //;*********************************************************************************************************
    
    //;Aqui, antes de comenzar el loop del tiempo, estan definidas las condiciones iniciales con valores de
    //;areai, Tsi, Tai, Wsi, Wai (5 vectores), definidos previamente como condiciones iniciales. Estas condiciones
    //;iniciales son las mismas cada que se cambia el forzamiento solar, es decir, cada que se comienza con un nuevo
    //;valor de luminosidad (L). Tambien estan implicitamente definidos los ajustes DAFM iniciales.
    
    //;*********************************************************************************************************
    //;BLOQUE 2 - INTEGRACIONES MEDIANTE RK4
    
    //;vectores (reciclables) para almacenar temporalmente los valores de las variables auxiliares en las 6 DAFs.
    zeroVector(Wvi,5);
    zeroVector(Preci,5);
    zeroVector(ETRi ,5);
    zeroVector(Qsi,5);
    zeroVector(SHsi,5);
    zeroVector(LHsi,5);
    zeroVector(LHai,5);
    zeroVector(Rsi,5);
    zeroVector(ROCsi,5);
    zeroVector(Rai,5);
    zeroVector(Cloudi,5);
    zeroVector(albci,5);
    zeroVector(RHi,5);
    zeroVector(BRi,5);
    zeroVector(BRalbi,5);

    //;mostrando condiciones iniciales
    if(VERBOSE or 0) printf("INITIAL CONDITIONS = \n");
    printVector("Tsi = ",Tsi,5,1);
    printVector("Tai = ",Tai,5,1);
    printVector("Wsi = ",Wsi,3,1);
    printVector("Wai = ",Wai,5,1);
    printVector("Wvi = ",Wvi,5,1);
    printVector("areai = ",areai,5,1);
    
    //;Loop de INTegracion
    niterINT=50;

    for(t=0;t<=niterINT-1;t++){
      //;******************************************************************************
      //;Runge Kutta 4 - Inicio		
      //;OJO: dentro de este loop de RK4 no se guardan valores. Notese que son con muy alta frecuencia.
      //;los valores se guardan es en el loop de integracion, cada deltatINT
      //;Loop RK4
      printf("t = %d\n",t);

      if(VERBOSE or 0) printf("%d\n",niterRK4);

      for(n=0;n<=niterRK4;n++){
	if(VERBOSE or 0) printf("n = %d\n",n);
	//;*************************************************************************************
	//;BLOQUE RESERVA: se guarda copia de las variables de estado en el instante t1
	//;porque los vectores areai, Tsi, ... se modifican para encontrar los valores en el
	//;instante siguiente (t2).
	//;*************************************************************************************			
	copyVector(areai_t1,areai,5);
	at_t1 = areai_t1[0];
	ag_t1 = areai_t1[1];
	ao_t1 = areai_t1[3];
	acr_t1 = areai_t1[4];
	copyVector(Tsi_t1,Tsi,5);
	copyVector(Tai_t1,Tai,5);
	copyVector(Wsi_t1,Wsi,5);
	copyVector(Wai_t1,Wai,5);
	copyVector(criosaT_t1,criosaT,6);
	//;Los valores con subindice t1 son los que se usan en la ec.final de RK4

	//;*************************************************************************************
	//;BLOQUE K1: obtencion de k1y=dydt1=F(y), donde y es cada una de las varibles de estado 
	//;*************************************************************************************
	//;constante solar regional
	acr = areai[4];
	tetacr = acos(1.-acr);
	inte = FunInte(tetacr,pinum/2.);					//;integral
	Sol = Scte*s1 + (Scte*s2 / cos(tetacr) )*inte;	//;constante solar para ocean-land
	if(acr>0.){Scr = (Scte-Sol*(1.-acr))/acr;}else{Scr = 0.;}
	//;flujos horizontales de energia
	FEai = FunFEa(areai, Tai);
	if(VERBOSE or 0) printf("Sol = %e\n",Sol);
	printVector("FEai = ",FEai,5,1);
	
	for(i=0;i<=4;i++){
	  
	  strcpy(DAF,DAFs[i]); //;DAFs= ['t','g','d','o','cr']
	  //strcpy(DAF,DAFs[4]); //;DAFs= ['t','g','d','o','cr']
	  
	  if(strcmp(DAF,"cr")==0){S=Scr;}else{S=Sol;}
	  area=areai[i];
	  if(VERBOSE) printf("RK4:Area = %f\n",area);

	  if(area>0.){
	    Ts = Tsi[i];
	    Ta = Tai[i];
	    Wa = Wai[i];
	    Ztoa = Ztoai[i];
	  
	    //;atmosfera, resultados = [Wv, Wpre, Wsat, RH, qv0, RH0, Ta0]
	    if(VERBOSE) 
	      printf("MAIN: FunAtm %s %e %e %e %e\n",DAFs[i],Ta,Wa,lapse,Ztoa);
	
	    atmos=FunAtm(Ta, Wa, lapse, Ztoa, DAF);
	    printVector("MAIN:Atmos:",atmos,7);
	    
	    //;atmosfera, resultados = [Wv, Wpre, Wsat, RH, qv0, RH0, Ta0]
	    atmos = FunAtm(Ta, Wa, lapse, Ztoa, DAF);
	    
	    Wv = atmos[0];		//;[kg m-2] agua precipitable que permanece en la atmosfera
	    Prec = atmos[1];		//;[kg m-2 yr-1] excedente de Wa que se vuelve precipitacion. Prec=Wpre/1yr
	    Wsat = atmos[2];		//;[kg m-2] capacidad maxima de almacenamiento de agua en la atmosfera
	    RH = atmos[3];		//;humedad relativa media en la columna
	    qv0 = atmos[4];		//;humedad especifica del aire en la superficie
	    RH0 = atmos[5];		//;humedad relativa del aire en la superficie
	    Ta0 = atmos[6];		//;temperatura del aire en la superficie
	    
	    //;nubes, resultados = [ac, albc, Tc, emc]
	    if(VERBOSE) 
	      printf("MAIN: FunCloud %s %e %e %e %e %e %e\n",DAFs[i],Prec, RH, Ztoa, lapse, qv0, Ta0);
	    clouds = FunCloud(Prec, RH, Ztoa, lapse, qv0, Ta0, DAF);
	    
	    printVector("Clouds:",clouds,4);
	    
	    ac = clouds[0];		//;area total de nubes
	    albc = clouds[1];	//;albedo promedio de las nubes
	    Tc = clouds[2];		//;temperatura promedio de las nubes
	    emc = clouds[3];		//;emisividad promedio de las nubes
	    
	    //;Net surface radiation
	    if(VERBOSE) printf("S=%e,Scr=%e,L=%e,Ts=%e,Ta=%e\n",S,Scr,L,Ts,Ta);
	    radnet = FunRs(S,L,Ts,Ta,RH,clouds,criosaT,DAF);		//;Rs en [MJ m-2 yr-1]
	    
	    Rs = radnet[0];
	    ROCs = radnet[1];
	    
	    //;Sensible heat - superficie
	    SHs = FunSH(Ts, Ta0, DAF);				//;[MJ m-2 yr-1]
	    if(VERBOSE) printf("SHs = %e\n",SHs);

	    //;Calor latente, evaporacion, transpiracion, sublimacion - superficie	
	    laten = FunLE(Ta0, Ztoa, Rs*1.e6, Prec, qv0, RH0, criosaT, DAF);	//;Rs en [J m-2 yr-1]
	    //;//;//;laten = FunLE(Ts,lapse,Rs*1.e6,Prec,qv10,criosaT,DAF);	//;Rs en [J m-2 yr-1]
	    LHs = laten[0];							//;[MJ m-2 yr-1]
	    ETR = laten[1];							//;[kg m-2 yr-1]
	    if(VERBOSE) printf("LHs = %e, ETR = %e\n",LHs,ETR);
	    
	    //;radiacion neta en la atmosfera
	    Ra = FunRa(S, L, Ts, Ta, RH, clouds, criosaT, DAF);		//;[MJ m-2 yr-1]
	    if(VERBOSE) printf("Ra = %e\n",Ra);
	    
	    //;calor latente en la atmosfera			
	    LHa = FunLHa(Prec);						//;[MJ m-2 yr-1]
	    if(VERBOSE) printf("LHa = %e\n",LHa);

	  }					
					//;AQUI ESTABA BIOTICREG
	  
	  else{
	    Wv = nodata;
	    Prec = nodata;
	    ETR = nodata;
	    SHs = nodata;
	    LHs = nodata;
	    LHa = nodata;
	    Rs = nodata;
	    ROCs = nodata;
	    Ra = nodata;
	    ac = nodata;
	    albc = nodata;
	    RH = nodata;
	  }

	  if(VERBOSE or 0) printf("%s %e %e %e %e %e %e %e %e %e %e %e\n",
		 DAFs[i],
		 Wv,
		 Prec,
		 ETR,
		 SHs,
		 LHs,
		 LHa,
		 Rs,
		 Ra,
		 ac,
		 albc,
		 RH);
	  
	  Wvi[i] = Wv;
	  Preci[i] = Prec;
	  ETRi[i] = ETR;
	  SHsi[i] = SHs;
	  LHsi[i] = LHs;
	  LHai[i] = LHa;
	  Rsi[i] = Rs;
	  ROCsi[i] = ROCs;
	  Rai[i] = Ra;
	  Cloudi[i] = ac;
	  albci[i] = albc;
	  RHi[i] = RH;

	}
	  
	FWai = FunFWa(areai,Tai,RHi);	//;[km m-2 yr] devuelve un vector
	printVector("FWai:",FWai,5,1);

	//;ajustes DAFM
	if((t==0)&&(n==0)){
	  //;condiciones iniciales. Los vectores son de la forma: [t,g,d,o,cr]. Excepto Wsi que es [t,g,d].
	  copyVector(DafmTsi,CIDafmTsi,5);		//;[MJ m-2 yr-1]
	  copyVector(DafmTai,CIDafmTai,5);		//;[MJ m-2 yr-1]
	  copyVector(DafmWsi,CIDafmWsi,3);		//;[kg m-2 yr-1]	
	  copyVector(DafmWai,CIDafmWai,5);		//;[kg m-2 yr-1]
	}else{
	  //;dareaidtDafm esta calculado al final del loop RK4
	  //;los ajustes Dafm se hacen teniendo en cuenta la separacion de la criosfera.
	  //;Esto se hace internamente en la funcion FunDafm
	  //;Por eso los vectores de entrada son de la forma: [t,g,d,o,ab,ac]. Excepto Wsi que es [t,g,d].
	  //;Sin embargo, en los vectores de salida se integra la criosfera y por eso
	  //;los resultados de FunDafm son AjusDafm = [KHt, KHg, KHd, KHo, KHcr]
	  printVector("areaiDafm = ",areaiDafm,5,1);
	  printVector("dareaidtDafm = ",dareaidtDafm,5,1);
	  printVector("cpsi = ",cpsi,5,1);
	  printVector("TsiDafm = ",TsiDafm,5,1);

	  //;la cantidad sujeta a conservacion es cps*Ts
	  cpTsiDafm[0]=cpsi[0]*TsiDafm[0];
	  cpTsiDafm[1]=cpsi[1]*TsiDafm[1];
	  cpTsiDafm[2]=cpsi[2]*TsiDafm[2];
	  cpTsiDafm[3]=cpsi[3]*TsiDafm[3];
	  cpTsiDafm[4]=cpsi[4]*TsiDafm[4];
	  cpTsiDafm[5]=cpc*TsiDafm[5];

	  printVector("cpTsiDafm = ",cpTsiDafm,6,1);

	  DafmTsiret = FunDafm(areaiDafm, dareaidtDafm,cpTsiDafm);
	  copyVector(DafmTsi,DafmTsiret,5);
	  
	  cpTaiDafm[0]=cpa*TaiDafm[0];
	  cpTaiDafm[1]=cpa*TaiDafm[1];
	  cpTaiDafm[2]=cpa*TaiDafm[2];
	  cpTaiDafm[3]=cpa*TaiDafm[3];
	  cpTaiDafm[4]=cpa*TaiDafm[4];
	  cpTaiDafm[5]=cpa*TaiDafm[5];
	  DafmTairet = FunDafm(areaiDafm, dareaidtDafm, cpTaiDafm);
	  copyVector(DafmTai,DafmTairet,5);
	  
	  DafmWsiret = FunDafmWs(areaiDafm, dareaidtDafm, Wsi);		//;OJO, es diferente a la funcion FunDafm
	  copyVector(DafmWsi,DafmWsiret,3);
	  
	  DafmWairet = FunDafm(areaiDafm, dareaidtDafm, WaiDafm);
	  copyVector(DafmWai,DafmWairet,5);

	  //;Aunque la criosfera como un todo se contraiga (dacrdt<0), puede que el ajuste dafm sea diferente de cero
	  //;porque la zona de acumulacion se expanda.
	  printVector("DafmTsi = ",DafmTsi,5,1);
	  printVector("DafmTai = ",DafmTai,5,1);
	  printVector("DafmWsi = ",DafmWsi,5,1);
	  printVector("DafmWai = ",DafmWai,5,1);
	}

	//;crecimiento/decrecimiento de las regiones (DAFs) ocupadas por especies vegetales (at y ag)
	at = areai[0];
	if(at>0.){
	  Tst = Tsi[0];
	  Wst = Wsi[0];
	  Wmt = Wmi[0];
	  Betat = FunBeta(Tst, Toptt,Toptt0);
	  tmort = FunTmor(Wst, Wmt);
	}else{
	  Betat = 0.;
	  tmort = 0.;
	}
	ag = areai[1];
	if(ag>0.){
	  Tsg = Tsi[1];
	  Wsg = Wsi[1];
	  Wmg = Wmi[1];
	  Betag = FunBeta(Tsg, Toptg,Toptg0);
	  tmorg = FunTmor(Wsg, Wmg);
	}else{
	  Betag = 0.;
	  tmorg = 0.;
	}
	if(VERBOSE or 0) printf("Betat,tmort = %e,%e\n",Betat,tmort);
	if(VERBOSE or 0) printf("Betag,tmorg = %e,%e\n",Betag,tmorg);

	las = areai[2]/(areai[0]+areai[1]+areai[2]);	//;land available space (fraction)
	k1at = at*(las*Betat-tmort);
	k1ag = ag*(las*Betag-tmorg);

	//;*****************************************************************
	//;Biotic Regulation, resultados = [BRalb, BRetr, deltaETR]
	for(i=0;i<=1;i++){
	  strcpy(DAF,DAFs[i]);
	  if(strcmp(DAF,"t")==0){
	    Topt=Toptt;
	    bs = las*Betat-tmort;
	  }else{
	    Topt=Toptg;
	    bs = las*Betag-tmorg;
	  }
	  if(VERBOSE or 0) printf("Input BioReg : %e %e %e %e %e %e %e %e\n",
		 Tsi[i],Topt,ROCsi[i],Rsi[i],LHsi[i],Preci[i],ETRi[i],bs);
	  BioReg = FunBioticReg(Tsi[i],Topt,ROCsi[i],Rsi[i],LHsi[i],Preci[i],ETRi[i],bs);
	  printVector("BioReg = ",BioReg,3,1);

	  BRalb = BioReg[0];
	  BRetr = BioReg[1];
	  deltaETR = BioReg[2];

	  BR = BRalb + BRetr;	//;total
				//;ajuste de flujos
	  ETR = ETR + deltaETR;	//;deltaETR contiene el signo
	  LHs = LHs + BRetr;	//;BRetr contiene el signo
				//;el ajuste BRalb se aplica directamente en la ecuacion de BES
	  BRi[i] = BR;
	  BRalbi[i] = BRalb;
	  if(VERBOSE or 0) printf("BR,ETR,LHs: %e %e %e\n",BR,ETR,LHs);
	}

	
	//;escorrentia superficial, melting
	//;[km m-2 yr] devuelve un vector
	Qsiret = FunQs(Preci, ETRi, Wsi, Wmi, areai, criosaT);
	copyVector(Qsi,Qsiret,5);
	printVector("Qsi = ",Qsi,5,1);
	
	//;crecimiento/decrecimiento del oceano (ao > 0. siempre)
	ao = areai[3];

	//;[yr-1] Qsi(3)=Qso<0 (entrada) contiene caudal (land) y 
	//;derretimiento (cryosphere). Estos calcuos estan en FunQs
	k1ao = ao*(Preci[3]-ETRi[3]-Qsi[3])/Moc;
		
	//;crecimiento/decrecimiento de la criosfera (incluye ambas zonas: ablacion y acumulacion)
	if(acr>0.){k1acr = acr*(Preci[4]-ETRi[4]-Qsi[4])/Mcr;}else{k1acr = 0.;}

	//;OJO, aqui NO se guardan variables como en el BLOQUE1

	if(VERBOSE or 0) printf("ao,k1ao,k1acr = %e, %e, %e\n",ao,k1ao,k1acr);

	//;*********************************************
	//;primeras derivadas en RK4 (los valores de k1)
	real k1areai[] = {k1at, k1ag, k1ao, k1acr};
	//;[yr-1] daidt1=k1ai=F(ai). Vector: [t,g,o,cr]. d se obtiene por conservacion y ab y ac se obtienen dividiendo la criosfera con la funcion FunCrios

	zeroVector(k1Tsi,5);
	zeroVector(k1Tai,5);
	zeroVector(k1Wai,5);
	zeroVector(k1Wsi,3);

	for(i=0;i<=4;i++){
	  area = areai[i];
	  if(area>0.){
	    //;[K yr-1] dTsdt1=k1Ts=F(Ts). Es un vector: [t,g,d,o,cr]. BRetr esta en LHs
	    k1Tsi[i] = (1./cpsi[i])*(Rsi[i]-SHsi[i]-LHsi[i]+BRalbi[i]+DafmTsi[i]);	
	    //;[K yr-1] dTadt1=k1Ta=F(Ta). Es un vector: [t,g,d,o,cr]
	    k1Tai[i] = (1./cpa)*(Rai[i]+SHsi[i]+LHai[i]+FEai[i]+DafmTai[i]);
	    //;[kg m-2 yr-1] dWadt1=k1Wa=F(Wa). Es un vector: [t,g,d,o,cr]
	    k1Wai[i] = ETRi[i]-Preci[i]+FWai[i]+DafmWai[i];
	  }else{
	    k1Tsi[i] = 0.;
	    k1Tai[i] = 0.;
	    k1Wai[i] = 0.;
	  }
	}

	for(i=0;i<=2;i++){
	  area = areai[i];
	  if(area>0.){
	    //;[kg m-2 yr-1] dWsdt1=k1Ws=F(Ws). Vector: [t,g,d]
	    k1Wsi[i] = Preci[i]-ETRi[i]-Qsi[i]+DafmWsi[i];
	  }else{	
	    k1Wsi[i] = 0.;
	  }
	}
	//;*********************************************

	printVector("k1Tsi = ",k1Tsi,5,1);
	printVector("k1Tai = ",k1Tai,5,1);
	printVector("k1Wai = ",k1Wai,5,1);
	printVector("k1Wsi = ",k1Wsi,3,1);
	
	//;********************************************************************************************
	//;BLOQUE K2: obtencion de k2y=dydt2=F(y+k1*h/2), donde y es cada una de las varibles de estado 
	//;********************************************************************************************
	//;Este bloque es una copia del bloque k1 pero evaluando todo en valores modificados de las variables
	//;de estado. Los nuevos valores son de la forma (p.e.) Tsi+k1Tsi*h/2

	//;Variables modificadas (equivalen a y+k1*h/2)
	for(i=0;i<=4;i++){
	  if(Tsi_t1[i]!=nodata){Tsi[i] = Tsi_t1[i] + k1Tsi[i]*deltatRK4/2.;}else{Tsi[i] = nodata;}
	  if(Tai_t1[i]!=nodata){Tai[i] = Tai_t1[i] + k1Tai[i]*deltatRK4/2.;}else{Tai[i] = nodata;}
	  if(i<2){
	    if(Wsi_t1[i]!=nodata){Wsi[i] = max(0., Wsi_t1[i] + k1Wsi[i]*deltatRK4/2.);}else{Wsi[i] = nodata;}
	  }
	  if(Wai_t1[i]!=nodata){Wai[i] = max(0., Wai_t1[i] + k1Wai[i]*deltatRK4/2.);}else{Wai[i] = nodata;}
	}
	at = max(0., at_t1 + k1at*deltatRK4/2.);//;trees
	ag = max(0., ag_t1 + k1ag*deltatRK4/2.);		//;grasses
	ao = ao_t1 + k1ao*deltatRK4/2.;					//;ocean
	acr = max(0., acr_t1 + k1acr*deltatRK4/2.);	//;cryosphere
	ad = 1.-(at+ag+ao+acr);							//;bare land
	areai[0]=at;areai[1]=ag;areai[2]=ad;areai[3]=ao;areai[4]=acr;

	//;Division de la criosfera en las condiciones modificadas (temperatura).
	acr = areai[4];
	if(acr>0.){
	  Tscr = Tsi[4];
	  Tacr = Tai[4];
	  tetacr = acos(1.-acr);
	  intecr = FunInte(0.,tetacr);
	  c2 = (1.-c1)*(1.-cos(tetacr))/intecr;
	  //;tiene incorporada la condicion de amin
	  FunCrios(Tscr,acr,tetacr,c1,c2,amin,crios);	
	  aac = crios[0];
	  aab = crios[1];
	}else{
	  aac = 0.;
	  aab = 0.;
	}
	if(aac>0.){
	  tetaac = acos(1.-aac);
	  inteac = FunInte(0.,tetaac);
	  //;ec (9) modelo criosfera en papel
	  Tsac = Tscr*c1 + Tscr*c2*inteac/(1. - cos(tetaac));	
	  //;ec (9) modelo criosfera en papel
	  Taac = Tacr*c1 + Tacr*c2*inteac/(1. - cos(tetaac));	
	}else{
	  inteac = 0.;		//;necesaria para zona de ablacion
	  Tsac = nodata;	//;no existe zona de acumulacion
	  Taac = nodata;
	}
	if(aab>0.){
	  if(aac==0.){tetaac = 0.;}//;si aac gt 0. tetaac ya fue calculada	
	  //;da lo mismo que inteab = FunInte(tetaac,tetacr)
	  inteab = intecr - inteac;
	  //;ec (9) modelo criosfera en papel
	  Tsab = Tscr*c1 + Tscr*c2*inteab/(cos(tetaac) - cos(tetacr));
	  //;ec (9) modelo criosfera en papel
	  Taab = Tacr*c1 + Tacr*c2*inteab/(cos(tetaac) - cos(tetacr));	
	}else{
	  Tsab = nodata;	//;no existe zona de ablacion
	  Taab = nodata;
	}
	criosaT[0]=aab;criosaT[1]=aac;criosaT[2]=Tsab;criosaT[3]=Tsac;criosaT[4]=Taab;
	criosaT[5]=Taac;//;actualizado
			//;control
	if(minVector(areai,5)<0.){
	  fprintf(stderr,"ERROR en areai, aparece valor negativo\n");
	  exit(1);
	}
	
	printVector("criosaT = ",criosaT,6,1);

	//;A partir de aqui bloque K2 es igual a bloque K1 (Notese que la diferencia esta en las variables modificadas)
	
	//;constante solar regional
	acr = areai[4];
	tetacr = acos(1.-acr);
	inte = FunInte(tetacr,pinum/2.);//;integral
	Sol = Scte*s1 + (Scte*s2 / cos(tetacr) )*inte;	//;constante solar para ocean-land
	if(acr>0.){Scr = (Scte-Sol*(1.-acr))/acr;}else{Scr = 0.;}

	//;flujos horizontales de energia
	FEai = FunFEa(areai, Tai);

	for(i=0;i<=4;i++){//;Loop por DAFs	

	  strcpy(DAF,DAFs[i]); //;DAFs= ['t','g','d','o','cr']
	  if(strcmp(DAF,"cr")==0){S=Scr;}else{S=Sol;}
	  
	  area=areai[i];

	  if(area>0.){
	    Ts = Tsi[i];
	    Ta = Tai[i];
	    Wa = Wai[i];
	    Ztoa = Ztoai[i];

	    atmos=FunAtm(Ta, Wa, lapse, Ztoa, DAF);

	    Wv = atmos[0];		//;[kg m-2] agua precipitable que permanece en la atmosfera
	    Prec = atmos[1];		//;[kg m-2 yr-1] excedente de Wa que se vuelve precipitacion. Prec=Wpre/1yr
	    Wsat = atmos[2];		//;[kg m-2] capacidad maxima de almacenamiento de agua en la atmosfera
	    RH = atmos[3];		//;humedad relativa media en la columna
	    qv0 = atmos[4];		//;humedad especifica del aire en la superficie
	    RH0 = atmos[5];		//;humedad relativa del aire en la superficie
	    Ta0 = atmos[6];		//;temperatura del aire en la superficie
	    
	    clouds = FunCloud(Prec, RH, Ztoa, lapse, qv0, Ta0, DAF);

	    ac = clouds[0];		//;area total de nubes
	    albc = clouds[1];	//;albedo promedio de las nubes
	    Tc = clouds[2];		//;temperatura promedio de las nubes
	    emc = clouds[3];		//;emisividad promedio de las nubes
	    
	    radnet = FunRs(S,L,Ts,Ta,RH,clouds,criosaT,DAF);		//;Rs en [MJ m-2 yr-1]

	    Rs = radnet[0];
	    ROCs = radnet[1];

	    SHs = FunSH(Ts, Ta0, DAF);				//;[MJ m-2 yr-1]
	    laten = FunLE(Ta0, Ztoa, Rs*1.e6, Prec, qv0, RH0, criosaT, DAF);	//;Rs en [J m-2 yr-1
	    LHs = laten[0];							//;[MJ m-2 yr-1]
	    ETR = laten[1];							//;[kg m-2 yr-1]
	    
	    //;radiacion neta en la atmosfera
	    Ra = FunRa(S, L, Ts, Ta, RH, clouds, criosaT, DAF);		//;[MJ m-2 yr-1]
	    LHa = FunLHa(Prec);						//;[MJ m-2 yr-1]

	  }else{		
	    Wv = nodata;
	    Prec = nodata;
	    ETR = nodata;		
	    SHs = nodata;
	    LHs = nodata;
	    LHa = nodata;
	    Rs = nodata;
	    ROCs = nodata;
	    Ra = nodata;
	    ac = nodata;
	    albc = nodata;
	    RH = nodata;
	  }

	  Wvi[i] = Wv;
	  Preci[i] = Prec;
	  ETRi[i] = ETR;
	  SHsi[i] = SHs;
	  LHsi[i] = LHs;
	  LHai[i] = LHa;
	  Rsi[i] = Rs;
	  Rai[i] = Ra;
	  ROCsi[i] = ROCs;
	  Cloudi[i] = ac;
	  albci[i] = albc;
	  RHi[i] = RH;

	}

	//;flujo horizontal de humedad atmosferica
	FWai = FunFWa(areai,Tai,RHi);	//;[km m-2 yr] devuelve un vector

	if((t==0)&&(n==0)){
	  //;condiciones iniciales. Los vectores son de la forma: [t,g,d,o,cr]. Excepto Wsi que es [t,g,d].
	  copyVector(DafmTsi,CIDafmTsi,5);		//;[MJ m-2 yr-1]
	  copyVector(DafmTai,CIDafmTai,5);		//;[MJ m-2 yr-1]
	  copyVector(DafmWsi,CIDafmWsi,3);		//;[kg m-2 yr-1]	
	  copyVector(DafmWai,CIDafmWai,5);		//;[kg m-2 yr-1]
	}else{
	  cpTsiDafm[0]=cpsi[0]*TsiDafm[0];
	  cpTsiDafm[1]=cpsi[1]*TsiDafm[1];
	  cpTsiDafm[2]=cpsi[2]*TsiDafm[2];
	  cpTsiDafm[3]=cpsi[3]*TsiDafm[3];
	  cpTsiDafm[4]=cpsi[4]*TsiDafm[4];
	  cpTsiDafm[5]=cpc*TsiDafm[5];

	  DafmTsiret = FunDafm(areaiDafm, dareaidtDafm,cpTsiDafm);
	  copyVector(DafmTsi,DafmTsiret,5);
	  
	  cpTaiDafm[0]=cpa*TaiDafm[0];
	  cpTaiDafm[1]=cpa*TaiDafm[1];
	  cpTaiDafm[2]=cpa*TaiDafm[2];
	  cpTaiDafm[3]=cpa*TaiDafm[3];
	  cpTaiDafm[4]=cpa*TaiDafm[4];
	  cpTaiDafm[5]=cpa*TaiDafm[5];
	  DafmTairet = FunDafm(areaiDafm, dareaidtDafm, cpTaiDafm);
	  copyVector(DafmTai,DafmTairet,5);
	  
	  DafmWsiret = FunDafmWs(areaiDafm, dareaidtDafm, Wsi);		//;OJO, es diferente a la funcion FunDafm
	  copyVector(DafmWsi,DafmWsiret,3);
	  
	  DafmWairet = FunDafm(areaiDafm, dareaidtDafm, WaiDafm);
	  copyVector(DafmWai,DafmWairet,5);
	}	  

	//;crecimiento/decrecimiento de las regiones (DAFs) ocupadas por especies vegetales (at y ag)
	at = areai[0];
	if(at>0.){
	  Tst = Tsi[0];
	  Wst = Wsi[0];
	  Wmt = Wmi[0];
	  Betat = FunBeta(Tst, Toptt,Toptt0);
	  tmort = FunTmor(Wst, Wmt);
	}else{
	  Betat = 0.;
	  tmort = 0.;
	}
	ag = areai[1];
	if(ag>0.){
	  Tsg = Tsi[1];
	  Wsg = Wsi[1];
	  Wmg = Wmi[1];
	  Betag = FunBeta(Tsg, Toptg,Toptg0);
	  tmorg = FunTmor(Wsg, Wmg);
	}else{
	  Betag = 0.;
	  tmorg = 0.;
	}
	las = areai[2]/(areai[0]+areai[1]+areai[2]);	//;land available space (fraction)
	k2at = at*(las*Betat-tmort);
	k2ag = ag*(las*Betag-tmorg);
	
	//;*****************************************************************
	//;Biotic Regulation, resultados = [BRalb, BRetr, deltaETR]
	for(i=0;i<=1;i++){
	  strcpy(DAF,DAFs[i]);
	  if(strcmp(DAF,"t")==0){
	    Topt=Toptt;
	    bs = las*Betat-tmort;
	  }else{
	    Topt=Toptg;
	    bs = las*Betag-tmorg;
	  }
	  if(VERBOSE or 0) printf("Input BioReg : %e %e %e %e %e %e %e %e\n",
		 Tsi[i],Topt,ROCsi[i],Rsi[i],LHsi[i],Preci[i],ETRi[i],bs);
	  BioReg = FunBioticReg(Tsi[i],Topt,ROCsi[i],Rsi[i],LHsi[i],Preci[i],ETRi[i],bs);
	  printVector("BioReg = ",BioReg,3,1);

	  BRalb = BioReg[0];
	  BRetr = BioReg[1];
	  deltaETR = BioReg[2];

	  BR = BRalb + BRetr;	//;total
				//;ajuste de flujos
	  ETR = ETR + deltaETR;	//;deltaETR contiene el signo
	  LHs = LHs + BRetr;	//;BRetr contiene el signo
				//;el ajuste BRalb se aplica directamente en la ecuacion de BES
	  BRi[i] = BR;
	  BRalbi[i] = BRalb;
	  if(VERBOSE or 0) printf("BR,ETR,LHs: %e %e %e\n",BR,ETR,LHs);
	}

	//;escorrentia superficial, melting
	//;[km m-2 yr] devuelve un vector
	Qsiret = FunQs(Preci, ETRi, Wsi, Wmi, areai, criosaT);
	copyVector(Qsi,Qsiret,5);
	printVector("Qsi = ",Qsi,5,1);
	
	//;crecimiento/decrecimiento del oceano (ao > 0. siempre)
	ao = areai[3];

	//;[yr-1] Qsi(3)=Qso<0 (entrada) contiene caudal (land) y 
	//;derretimiento (cryosphere). Estos calcuos estan en FunQs
	k2ao = ao*(Preci[3]-ETRi[3]-Qsi[3])/Moc;


	//;crecimiento/decrecimiento de la criosfera (incluye ambas zonas: ablacion y acumulacion)
	if(acr>0.){k2acr = acr*(Preci[4]-ETRi[4]-Qsi[4])/Mcr;}else{k2acr = 0.;}



	//;*********************************************
	//;segundas derivadas en RK4 (los valores de k2)

	real k2areai[] = {k2at, k2ag, k2ao, k2acr};
	//;[yr-1] daidt1=k2ai=F(ai). Vector: [t,g,o,cr]. d se obtiene por conservacion y ab y ac se obtienen dividiendo la criosfera con la funcion FunCrios

	zeroVector(k2Tsi,5);
	zeroVector(k2Tai,5);
	zeroVector(k2Wai,5);
	zeroVector(k2Wsi,3);

	for(i=0;i<=4;i++){
	  area = areai[i];
	  if(area>0.){
	    //;[K yr-1] dTsdt1=k2Ts=F(Ts). Es un vector: [t,g,d,o,cr]. BRetr esta en LHs
	    k2Tsi[i] = (1./cpsi[i])*(Rsi[i]-SHsi[i]-LHsi[i]+BRalbi[i]+DafmTsi[i]);	
	    //;[K yr-1] dTadt1=k2Ta=F(Ta). Es un vector: [t,g,d,o,cr]
	    k2Tai[i] = (1./cpa)*(Rai[i]+SHsi[i]+LHai[i]+FEai[i]+DafmTai[i]);
	    //;[kg m-2 yr-1] dWadt1=k2Wa=F(Wa). Es un vector: [t,g,d,o,cr]
	    k2Wai[i] = ETRi[i]-Preci[i]+FWai[i]+DafmWai[i];
	  }else{
	    k2Tsi[i] = 0.;
	    k2Tai[i] = 0.;
	    k2Wai[i] = 0.;
	  }
	}

	for(i=0;i<=2;i++){
	  area = areai[i];
	  if(area>0.){
	    //;[kg m-2 yr-1] dWsdt1=k2Ws=F(Ws). Vector: [t,g,d]
	    k2Wsi[i] = Preci[i]-ETRi[i]-Qsi[i]+DafmWsi[i];
	  }else{	
	    k2Wsi[i] = 0.;
	  }
	}

	printVector("k2Tsi = ",k2Tsi,5,1);
	printVector("k2Tai = ",k2Tai,5,1);
	printVector("k2Wai = ",k2Wai,5,1);
	printVector("k2Wsi = ",k2Wsi,3,1);

	//;*********************************************
	//;********************************************************************************************
	//;BLOQUE K3: obtencion de k3y=dydt3=F(y+k2*h/2), donde y es cada una de las varibles de estado 
	//;********************************************************************************************
	//;Este bloque es una copia del bloque k1 pero evaluando todo en valores modificados de las variables
	//;de estado. Los nuevos valores son de la forma (p.e.) Tsi+k2Tsi*h/2

	//;Variables modificadas (equivalen a y+k2*h/2)
	for(i=0;i<=4;i++){
	  if(Tsi_t1[i]!=nodata){Tsi[i] = Tsi_t1[i] + k2Tsi[i]*deltatRK4/2.;}else{Tsi[i] = nodata;}
	  if(Tai_t1[i]!=nodata){Tai[i] = Tai_t1[i] + k2Tai[i]*deltatRK4/2.;}else{Tai[i] = nodata;}
	  if(i<2){
	    if(Wsi_t1[i]!=nodata){Wsi[i] = max(0., Wsi_t1[i] + k2Wsi[i]*deltatRK4/2.);}else{Wsi[i] = nodata;}
	  }
	  if(Wai_t1[i]!=nodata){Wai[i] = max(0., Wai_t1[i] + k2Wai[i]*deltatRK4/2.);}else{Wai[i] = nodata;}
	}
	at = max(0., at_t1 + k2at*deltatRK4/2.);//;trees
	ag = max(0., ag_t1 + k2ag*deltatRK4/2.);		//;grasses
	ao = ao_t1 + k2ao*deltatRK4/2.;					//;ocean
	acr = max(0., acr_t1 + k2acr*deltatRK4/2.);	//;cryosphere
	ad = 1.-(at+ag+ao+acr);							//;bare land
	areai[0]=at;areai[1]=ag;areai[2]=ad;areai[3]=ao;areai[4]=acr;

	//;Division de la criosfera en las condiciones modificadas (temperatura).
	acr = areai[4];
	if(acr>0.){
	  Tscr = Tsi[4];
	  Tacr = Tai[4];
	  tetacr = acos(1.-acr);
	  intecr = FunInte(0.,tetacr);
	  c2 = (1.-c1)*(1.-cos(tetacr))/intecr;
	  //;tiene incorporada la condicion de amin
	  FunCrios(Tscr,acr,tetacr,c1,c2,amin,crios);	
	  aac = crios[0];
	  aab = crios[1];
	}else{
	  aac = 0.;
	  aab = 0.;
	}
	if(aac>0.){
	  tetaac = acos(1.-aac);
	  inteac = FunInte(0.,tetaac);
	  //;ec (9) modelo criosfera en papel
	  Tsac = Tscr*c1 + Tscr*c2*inteac/(1. - cos(tetaac));	
	  //;ec (9) modelo criosfera en papel
	  Taac = Tacr*c1 + Tacr*c2*inteac/(1. - cos(tetaac));	
	}else{
	  inteac = 0.;		//;necesaria para zona de ablacion
	  Tsac = nodata;	//;no existe zona de acumulacion
	  Taac = nodata;
	}
	if(aab>0.){
	  if(aac==0.){tetaac = 0.;}//;si aac gt 0. tetaac ya fue calculada	
	  //;da lo mismo que inteab = FunInte(tetaac,tetacr)
	  inteab = intecr - inteac;
	  //;ec (9) modelo criosfera en papel
	  Tsab = Tscr*c1 + Tscr*c2*inteab/(cos(tetaac) - cos(tetacr));
	  //;ec (9) modelo criosfera en papel
	  Taab = Tacr*c1 + Tacr*c2*inteab/(cos(tetaac) - cos(tetacr));	
	}else{
	  Tsab = nodata;	//;no existe zona de ablacion
	  Taab = nodata;
	}
	criosaT[0]=aab;criosaT[1]=aac;criosaT[2]=Tsab;criosaT[3]=Tsac;criosaT[4]=Taab;
	criosaT[5]=Taac;//;actualizado
			//;control
	if(minVector(areai,5)<0.){
	  fprintf(stderr,"ERROR en areai, aparece valor negativo\n");
	  exit(1);
	}
	
	printVector("criosaT = ",criosaT,6,1);

	//;A partir de aqui bloque K2 es igual a bloque K1 (Notese que la diferencia esta en las variables modificadas)
	
	//;constante solar regional
	acr = areai[4];
	tetacr = acos(1.-acr);
	inte = FunInte(tetacr,pinum/2.);//;integral
	Sol = Scte*s1 + (Scte*s2 / cos(tetacr) )*inte;	//;constante solar para ocean-land
	if(acr>0.){Scr = (Scte-Sol*(1.-acr))/acr;}else{Scr = 0.;}

	//;flujos horizontales de energia
	FEai = FunFEa(areai, Tai);

	for(i=0;i<=4;i++){//;Loop por DAFs	

	  strcpy(DAF,DAFs[i]); //;DAFs= ['t','g','d','o','cr']
	  if(strcmp(DAF,"cr")==0){S=Scr;}else{S=Sol;}
	  
	  area=areai[i];

	  if(area>0.){
	    Ts = Tsi[i];
	    Ta = Tai[i];
	    Wa = Wai[i];
	    Ztoa = Ztoai[i];

	    atmos=FunAtm(Ta, Wa, lapse, Ztoa, DAF);

	    Wv = atmos[0];		//;[kg m-2] agua precipitable que permanece en la atmosfera
	    Prec = atmos[1];		//;[kg m-2 yr-1] excedente de Wa que se vuelve precipitacion. Prec=Wpre/1yr
	    Wsat = atmos[2];		//;[kg m-2] capacidad maxima de almacenamiento de agua en la atmosfera
	    RH = atmos[3];		//;humedad relativa media en la columna
	    qv0 = atmos[4];		//;humedad especifica del aire en la superficie
	    RH0 = atmos[5];		//;humedad relativa del aire en la superficie
	    Ta0 = atmos[6];		//;temperatura del aire en la superficie
	    
	    clouds = FunCloud(Prec, RH, Ztoa, lapse, qv0, Ta0, DAF);

	    ac = clouds[0];		//;area total de nubes
	    albc = clouds[1];	//;albedo promedio de las nubes
	    Tc = clouds[2];		//;temperatura promedio de las nubes
	    emc = clouds[3];		//;emisividad promedio de las nubes
	    
	    radnet = FunRs(S,L,Ts,Ta,RH,clouds,criosaT,DAF);		//;Rs en [MJ m-2 yr-1]

	    Rs = radnet[0];
	    ROCs = radnet[1];

	    SHs = FunSH(Ts, Ta0, DAF);				//;[MJ m-2 yr-1]
	    laten = FunLE(Ta0, Ztoa, Rs*1.e6, Prec, qv0, RH0, criosaT, DAF);	//;Rs en [J m-2 yr-1
	    LHs = laten[0];							//;[MJ m-2 yr-1]
	    ETR = laten[1];							//;[kg m-2 yr-1]
	    
	    //;radiacion neta en la atmosfera
	    Ra = FunRa(S, L, Ts, Ta, RH, clouds, criosaT, DAF);		//;[MJ m-2 yr-1]
	    LHa = FunLHa(Prec);						//;[MJ m-2 yr-1]

	  }else{		
	    Wv = nodata;
	    Prec = nodata;
	    ETR = nodata;		
	    SHs = nodata;
	    LHs = nodata;
	    LHa = nodata;
	    Rs = nodata;
	    ROCs = nodata;
	    Ra = nodata;
	    ac = nodata;
	    albc = nodata;
	    RH = nodata;
	  }

	  Wvi[i] = Wv;
	  Preci[i] = Prec;
	  ETRi[i] = ETR;
	  SHsi[i] = SHs;
	  LHsi[i] = LHs;
	  LHai[i] = LHa;
	  Rsi[i] = Rs;
	  Rai[i] = Ra;
	  ROCsi[i] = ROCs;
	  Cloudi[i] = ac;
	  albci[i] = albc;
	  RHi[i] = RH;

	}

	//;flujo horizontal de humedad atmosferica
	FWai = FunFWa(areai,Tai,RHi);	//;[km m-2 yr] devuelve un vector

	if((t==0)&&(n==0)){
	  //;condiciones iniciales. Los vectores son de la forma: [t,g,d,o,cr]. Excepto Wsi que es [t,g,d].
	  copyVector(DafmTsi,CIDafmTsi,5);		//;[MJ m-2 yr-1]
	  copyVector(DafmTai,CIDafmTai,5);		//;[MJ m-2 yr-1]
	  copyVector(DafmWsi,CIDafmWsi,3);		//;[kg m-2 yr-1]	
	  copyVector(DafmWai,CIDafmWai,5);		//;[kg m-2 yr-1]
	}else{
	  cpTsiDafm[0]=cpsi[0]*TsiDafm[0];
	  cpTsiDafm[1]=cpsi[1]*TsiDafm[1];
	  cpTsiDafm[2]=cpsi[2]*TsiDafm[2];
	  cpTsiDafm[3]=cpsi[3]*TsiDafm[3];
	  cpTsiDafm[4]=cpsi[4]*TsiDafm[4];
	  cpTsiDafm[5]=cpc*TsiDafm[5];

	  DafmTsiret = FunDafm(areaiDafm, dareaidtDafm,cpTsiDafm);
	  copyVector(DafmTsi,DafmTsiret,5);
	  
	  cpTaiDafm[0]=cpa*TaiDafm[0];
	  cpTaiDafm[1]=cpa*TaiDafm[1];
	  cpTaiDafm[2]=cpa*TaiDafm[2];
	  cpTaiDafm[3]=cpa*TaiDafm[3];
	  cpTaiDafm[4]=cpa*TaiDafm[4];
	  cpTaiDafm[5]=cpa*TaiDafm[5];
	  DafmTairet = FunDafm(areaiDafm, dareaidtDafm, cpTaiDafm);
	  copyVector(DafmTai,DafmTairet,5);
	  
	  DafmWsiret = FunDafmWs(areaiDafm, dareaidtDafm, Wsi);		//;OJO, es diferente a la funcion FunDafm
	  copyVector(DafmWsi,DafmWsiret,3);
	  
	  DafmWairet = FunDafm(areaiDafm, dareaidtDafm, WaiDafm);
	  copyVector(DafmWai,DafmWairet,5);
	}	  

	//;crecimiento/decrecimiento de las regiones (DAFs) ocupadas por especies vegetales (at y ag)
	at = areai[0];
	if(at>0.){
	  Tst = Tsi[0];
	  Wst = Wsi[0];
	  Wmt = Wmi[0];
	  Betat = FunBeta(Tst, Toptt,Toptt0);
	  tmort = FunTmor(Wst, Wmt);
	}else{
	  Betat = 0.;
	  tmort = 0.;
	}
	ag = areai[1];
	if(ag>0.){
	  Tsg = Tsi[1];
	  Wsg = Wsi[1];
	  Wmg = Wmi[1];
	  Betag = FunBeta(Tsg, Toptg,Toptg0);
	  tmorg = FunTmor(Wsg, Wmg);
	}else{
	  Betag = 0.;
	  tmorg = 0.;
	}
	las = areai[2]/(areai[0]+areai[1]+areai[2]);	//;land available space (fraction)
	k3at = at*(las*Betat-tmort);
	k3ag = ag*(las*Betag-tmorg);
	
	//;*****************************************************************
	//;Biotic Regulation, resultados = [BRalb, BRetr, deltaETR]
	for(i=0;i<=1;i++){
	  strcpy(DAF,DAFs[i]);
	  if(strcmp(DAF,"t")==0){
	    Topt=Toptt;
	    bs = las*Betat-tmort;
	  }else{
	    Topt=Toptg;
	    bs = las*Betag-tmorg;
	  }
	  if(VERBOSE or 0) printf("Input BioReg : %e %e %e %e %e %e %e %e\n",
		 Tsi[i],Topt,ROCsi[i],Rsi[i],LHsi[i],Preci[i],ETRi[i],bs);
	  BioReg = FunBioticReg(Tsi[i],Topt,ROCsi[i],Rsi[i],LHsi[i],Preci[i],ETRi[i],bs);
	  printVector("BioReg = ",BioReg,3,1);

	  BRalb = BioReg[0];
	  BRetr = BioReg[1];
	  deltaETR = BioReg[2];

	  BR = BRalb + BRetr;	//;total
				//;ajuste de flujos
	  ETR = ETR + deltaETR;	//;deltaETR contiene el signo
	  LHs = LHs + BRetr;	//;BRetr contiene el signo
				//;el ajuste BRalb se aplica directamente en la ecuacion de BES
	  BRi[i] = BR;
	  BRalbi[i] = BRalb;
	  if(VERBOSE or 0) printf("BR,ETR,LHs: %e %e %e\n",BR,ETR,LHs);
	}

	//;escorrentia superficial, melting
	//;[km m-2 yr] devuelve un vector
	Qsiret = FunQs(Preci, ETRi, Wsi, Wmi, areai, criosaT);
	copyVector(Qsi,Qsiret,5);
	printVector("Qsi = ",Qsi,5,1);
	
	//;crecimiento/decrecimiento del oceano (ao > 0. siempre)
	ao = areai[3];

	//;[yr-1] Qsi(3)=Qso<0 (entrada) contiene caudal (land) y 
	//;derretimiento (cryosphere). Estos calcuos estan en FunQs
	k3ao = ao*(Preci[3]-ETRi[3]-Qsi[3])/Moc;


	//;crecimiento/decrecimiento de la criosfera (incluye ambas zonas: ablacion y acumulacion)
	if(acr>0.){k3acr = acr*(Preci[4]-ETRi[4]-Qsi[4])/Mcr;}else{k3acr = 0.;}



	//;*********************************************
	//;segundas derivadas en RK4 (los valores de k3)

	real k3areai[] = {k3at, k3ag, k3ao, k3acr};
	//;[yr-1] daidt1=k3ai=F(ai). Vector: [t,g,o,cr]. d se obtiene por conservacion y ab y ac se obtienen dividiendo la criosfera con la funcion FunCrios

	zeroVector(k3Tsi,5);
	zeroVector(k3Tai,5);
	zeroVector(k3Wai,5);
	zeroVector(k3Wsi,3);

	for(i=0;i<=4;i++){
	  area = areai[i];
	  if(area>0.){
	    //;[K yr-1] dTsdt1=k3Ts=F(Ts). Es un vector: [t,g,d,o,cr]. BRetr esta en LHs
	    k3Tsi[i] = (1./cpsi[i])*(Rsi[i]-SHsi[i]-LHsi[i]+BRalbi[i]+DafmTsi[i]);	
	    //;[K yr-1] dTadt1=k3Ta=F(Ta). Es un vector: [t,g,d,o,cr]
	    k3Tai[i] = (1./cpa)*(Rai[i]+SHsi[i]+LHai[i]+FEai[i]+DafmTai[i]);
	    //;[kg m-2 yr-1] dWadt1=k3Wa=F(Wa). Es un vector: [t,g,d,o,cr]
	    k3Wai[i] = ETRi[i]-Preci[i]+FWai[i]+DafmWai[i];
	  }else{
	    k3Tsi[i] = 0.;
	    k3Tai[i] = 0.;
	    k3Wai[i] = 0.;
	  }
	}

	for(i=0;i<=2;i++){
	  area = areai[i];
	  if(area>0.){
	    //;[kg m-2 yr-1] dWsdt1=k3Ws=F(Ws). Vector: [t,g,d]
	    k3Wsi[i] = Preci[i]-ETRi[i]-Qsi[i]+DafmWsi[i];
	  }else{	
	    k3Wsi[i] = 0.;
	  }
	}

	printVector("k3Tsi = ",k3Tsi,5,1);
	printVector("k3Tai = ",k3Tai,5,1);
	printVector("k3Wai = ",k3Wai,5,1);
	printVector("k3Wsi = ",k3Wsi,3,1);

	//;********************************************************************************************
	//;BLOQUE K4: obtencion de k4y=dydt4=F(y+k3*h), donde y es cada una de las varibles de estado 
	//;********************************************************************************************
	//;Este bloque es una copia del bloque k1 pero evaluando todo en valores modificados de las variables
	//;de estado. Los nuevos valores son de la forma (p.e.) Tsi+k3Tsi*h
	for(i=0;i<=4;i++){
	  if(Tsi_t1[i]!=nodata){Tsi[i] = Tsi_t1[i] + k3Tsi[i]*deltatRK4;}else{Tsi[i] = nodata;}
	  if(Tai_t1[i]!=nodata){Tai[i] = Tai_t1[i] + k3Tai[i]*deltatRK4;}else{Tai[i] = nodata;}
	  if(i<2){
	    if(Wsi_t1[i]!=nodata){Wsi[i] = max(0., Wsi_t1[i] + k3Wsi[i]*deltatRK4);}else{Wsi[i] = nodata;}
	  }
	  if(Wai_t1[i]!=nodata){Wai[i] = max(0., Wai_t1[i] + k3Wai[i]*deltatRK4);}else{Wai[i] = nodata;}
	}
	at = max(0., at_t1 + k3at*deltatRK4);//;trees
	ag = max(0., ag_t1 + k3ag*deltatRK4);		//;grasses
	ao = ao_t1 + k3ao*deltatRK4;					//;ocean
	acr = max(0., acr_t1 + k3acr*deltatRK4);	//;cryosphere
	ad = 1.-(at+ag+ao+acr);							//;bare land
	areai[0]=at;areai[1]=ag;areai[2]=ad;areai[3]=ao;areai[4]=acr;

	//;Division de la criosfera en las condiciones modificadas (temperatura).
	acr = areai[4];
	if(acr>0.){
	  Tscr = Tsi[4];
	  Tacr = Tai[4];
	  tetacr = acos(1.-acr);
	  intecr = FunInte(0.,tetacr);
	  c2 = (1.-c1)*(1.-cos(tetacr))/intecr;
	  //;tiene incorporada la condicion de amin
	  FunCrios(Tscr,acr,tetacr,c1,c2,amin,crios);	
	  aac = crios[0];
	  aab = crios[1];
	}else{
	  aac = 0.;
	  aab = 0.;
	}
	if(aac>0.){
	  tetaac = acos(1.-aac);
	  inteac = FunInte(0.,tetaac);
	  //;ec (9) modelo criosfera en papel
	  Tsac = Tscr*c1 + Tscr*c2*inteac/(1. - cos(tetaac));	
	  //;ec (9) modelo criosfera en papel
	  Taac = Tacr*c1 + Tacr*c2*inteac/(1. - cos(tetaac));	
	}else{
	  inteac = 0.;		//;necesaria para zona de ablacion
	  Tsac = nodata;	//;no existe zona de acumulacion
	  Taac = nodata;
	}
	if(aab>0.){
	  if(aac==0.){tetaac = 0.;}//;si aac gt 0. tetaac ya fue calculada	
	  //;da lo mismo que inteab = FunInte(tetaac,tetacr)
	  inteab = intecr - inteac;
	  //;ec (9) modelo criosfera en papel
	  Tsab = Tscr*c1 + Tscr*c2*inteab/(cos(tetaac) - cos(tetacr));
	  //;ec (9) modelo criosfera en papel
	  Taab = Tacr*c1 + Tacr*c2*inteab/(cos(tetaac) - cos(tetacr));	
	}else{
	  Tsab = nodata;	//;no existe zona de ablacion
	  Taab = nodata;
	}
	criosaT[0]=aab;criosaT[1]=aac;criosaT[2]=Tsab;criosaT[3]=Tsac;criosaT[4]=Taab;
	criosaT[5]=Taac;//;actualizado
			//;control
	if(minVector(areai,5)<0.){
	  fprintf(stderr,"ERROR en areai, aparece valor negativo\n");
	  exit(1);
	}
	
	printVector("criosaT = ",criosaT,6,1);
	
	//;constante solar regional
	acr = areai[4];
	tetacr = acos(1.-acr);
	inte = FunInte(tetacr,pinum/2.);//;integral
	Sol = Scte*s1 + (Scte*s2 / cos(tetacr) )*inte;	//;constante solar para ocean-land
	if(acr>0.){Scr = (Scte-Sol*(1.-acr))/acr;}else{Scr = 0.;}

	//;flujos horizontales de energia
	FEai = FunFEa(areai, Tai);

	for(i=0;i<=4;i++){//;Loop por DAFs	

	  strcpy(DAF,DAFs[i]); //;DAFs= ['t','g','d','o','cr']
	  if(strcmp(DAF,"cr")==0){S=Scr;}else{S=Sol;}
	  
	  area=areai[i];

	  if(area>0.){
	    Ts = Tsi[i];
	    Ta = Tai[i];
	    Wa = Wai[i];
	    Ztoa = Ztoai[i];

	    atmos=FunAtm(Ta, Wa, lapse, Ztoa, DAF);

	    Wv = atmos[0];		//;[kg m-2] agua precipitable que permanece en la atmosfera
	    Prec = atmos[1];		//;[kg m-2 yr-1] excedente de Wa que se vuelve precipitacion. Prec=Wpre/1yr
	    Wsat = atmos[2];		//;[kg m-2] capacidad maxima de almacenamiento de agua en la atmosfera
	    RH = atmos[3];		//;humedad relativa media en la columna
	    qv0 = atmos[4];		//;humedad especifica del aire en la superficie
	    RH0 = atmos[5];		//;humedad relativa del aire en la superficie
	    Ta0 = atmos[6];		//;temperatura del aire en la superficie
	    
	    clouds = FunCloud(Prec, RH, Ztoa, lapse, qv0, Ta0, DAF);

	    ac = clouds[0];		//;area total de nubes
	    albc = clouds[1];	//;albedo promedio de las nubes
	    Tc = clouds[2];		//;temperatura promedio de las nubes
	    emc = clouds[3];		//;emisividad promedio de las nubes
	    
	    radnet = FunRs(S,L,Ts,Ta,RH,clouds,criosaT,DAF);		//;Rs en [MJ m-2 yr-1]

	    Rs = radnet[0];
	    ROCs = radnet[1];

	    SHs = FunSH(Ts, Ta0, DAF);				//;[MJ m-2 yr-1]
	    laten = FunLE(Ta0, Ztoa, Rs*1.e6, Prec, qv0, RH0, criosaT, DAF);	//;Rs en [J m-2 yr-1
	    LHs = laten[0];							//;[MJ m-2 yr-1]
	    ETR = laten[1];							//;[kg m-2 yr-1]
	    
	    //;radiacion neta en la atmosfera
	    Ra = FunRa(S, L, Ts, Ta, RH, clouds, criosaT, DAF);		//;[MJ m-2 yr-1]
	    LHa = FunLHa(Prec);						//;[MJ m-2 yr-1]

	  }else{		
	    Wv = nodata;
	    Prec = nodata;
	    ETR = nodata;		
	    SHs = nodata;
	    LHs = nodata;
	    LHa = nodata;
	    Rs = nodata;
	    ROCs = nodata;
	    Ra = nodata;
	    ac = nodata;
	    albc = nodata;
	    RH = nodata;
	  }

	  Wvi[i] = Wv;
	  Preci[i] = Prec;
	  ETRi[i] = ETR;
	  SHsi[i] = SHs;
	  LHsi[i] = LHs;
	  LHai[i] = LHa;
	  Rsi[i] = Rs;
	  Rai[i] = Ra;
	  ROCsi[i] = ROCs;
	  Cloudi[i] = ac;
	  albci[i] = albc;
	  RHi[i] = RH;

	}

	//;flujo horizontal de humedad atmosferica
	FWai = FunFWa(areai,Tai,RHi);	//;[km m-2 yr] devuelve un vector

	if((t==0)&&(n==0)){
	  //;condiciones iniciales. Los vectores son de la forma: [t,g,d,o,cr]. Excepto Wsi que es [t,g,d].
	  copyVector(DafmTsi,CIDafmTsi,5);		//;[MJ m-2 yr-1]
	  copyVector(DafmTai,CIDafmTai,5);		//;[MJ m-2 yr-1]
	  copyVector(DafmWsi,CIDafmWsi,3);		//;[kg m-2 yr-1]	
	  copyVector(DafmWai,CIDafmWai,5);		//;[kg m-2 yr-1]
	}else{
	  cpTsiDafm[0]=cpsi[0]*TsiDafm[0];
	  cpTsiDafm[1]=cpsi[1]*TsiDafm[1];
	  cpTsiDafm[2]=cpsi[2]*TsiDafm[2];
	  cpTsiDafm[3]=cpsi[3]*TsiDafm[3];
	  cpTsiDafm[4]=cpsi[4]*TsiDafm[4];
	  cpTsiDafm[5]=cpc*TsiDafm[5];

	  DafmTsiret = FunDafm(areaiDafm, dareaidtDafm,cpTsiDafm);
	  copyVector(DafmTsi,DafmTsiret,5);
	  
	  cpTaiDafm[0]=cpa*TaiDafm[0];
	  cpTaiDafm[1]=cpa*TaiDafm[1];
	  cpTaiDafm[2]=cpa*TaiDafm[2];
	  cpTaiDafm[3]=cpa*TaiDafm[3];
	  cpTaiDafm[4]=cpa*TaiDafm[4];
	  cpTaiDafm[5]=cpa*TaiDafm[5];
	  DafmTairet = FunDafm(areaiDafm, dareaidtDafm, cpTaiDafm);
	  copyVector(DafmTai,DafmTairet,5);
	  
	  DafmWsiret = FunDafmWs(areaiDafm, dareaidtDafm, Wsi);		//;OJO, es diferente a la funcion FunDafm
	  copyVector(DafmWsi,DafmWsiret,3);
	  
	  DafmWairet = FunDafm(areaiDafm, dareaidtDafm, WaiDafm);
	  copyVector(DafmWai,DafmWairet,5);
	}	  

	//;crecimiento/decrecimiento de las regiones (DAFs) ocupadas por especies vegetales (at y ag)
	at = areai[0];
	if(at>0.){
	  Tst = Tsi[0];
	  Wst = Wsi[0];
	  Wmt = Wmi[0];
	  Betat = FunBeta(Tst, Toptt,Toptt0);
	  tmort = FunTmor(Wst, Wmt);
	}else{
	  Betat = 0.;
	  tmort = 0.;
	}
	ag = areai[1];
	if(ag>0.){
	  Tsg = Tsi[1];
	  Wsg = Wsi[1];
	  Wmg = Wmi[1];
	  Betag = FunBeta(Tsg, Toptg,Toptg0);
	  tmorg = FunTmor(Wsg, Wmg);
	}else{
	  Betag = 0.;
	  tmorg = 0.;
	}
	las = areai[2]/(areai[0]+areai[1]+areai[2]);	//;land available space (fraction)
	k4at = at*(las*Betat-tmort);
	k4ag = ag*(las*Betag-tmorg);
	
	//;*****************************************************************
	//;Biotic Regulation, resultados = [BRalb, BRetr, deltaETR]
	for(i=0;i<=1;i++){
	  strcpy(DAF,DAFs[i]);
	  if(strcmp(DAF,"t")==0){
	    Topt=Toptt;
	    bs = las*Betat-tmort;
	  }else{
	    Topt=Toptg;
	    bs = las*Betag-tmorg;
	  }
	  if(VERBOSE or 0) printf("Input BioReg : %e %e %e %e %e %e %e %e\n",
		 Tsi[i],Topt,ROCsi[i],Rsi[i],LHsi[i],Preci[i],ETRi[i],bs);
	  BioReg = FunBioticReg(Tsi[i],Topt,ROCsi[i],Rsi[i],LHsi[i],Preci[i],ETRi[i],bs);
	  printVector("BioReg = ",BioReg,3,1);

	  BRalb = BioReg[0];
	  BRetr = BioReg[1];
	  deltaETR = BioReg[2];

	  BR = BRalb + BRetr;	//;total
				//;ajuste de flujos
	  ETR = ETR + deltaETR;	//;deltaETR contiene el signo
	  LHs = LHs + BRetr;	//;BRetr contiene el signo
				//;el ajuste BRalb se aplica directamente en la ecuacion de BES
	  BRi[i] = BR;
	  BRalbi[i] = BRalb;
	  if(VERBOSE or 0) printf("BR,ETR,LHs: %e %e %e\n",BR,ETR,LHs);
	}

	//;escorrentia superficial, melting
	//;[km m-2 yr] devuelve un vector
	Qsiret = FunQs(Preci, ETRi, Wsi, Wmi, areai, criosaT);
	copyVector(Qsi,Qsiret,5);
	printVector("Qsi = ",Qsi,5,1);
	
	//;crecimiento/decrecimiento del oceano (ao > 0. siempre)
	ao = areai[3];

	//;[yr-1] Qsi(3)=Qso<0 (entrada) contiene caudal (land) y 
	//;derretimiento (cryosphere). Estos calcuos estan en FunQs
	k4ao = ao*(Preci[3]-ETRi[3]-Qsi[3])/Moc;


	//;crecimiento/decrecimiento de la criosfera (incluye ambas zonas: ablacion y acumulacion)
	if(acr>0.){k4acr = acr*(Preci[4]-ETRi[4]-Qsi[4])/Mcr;}else{k4acr = 0.;}



	//;*********************************************
	//;segundas derivadas en RK4 (los valores de k4)

	real k4areai[] = {k4at, k4ag, k4ao, k4acr};
	//;[yr-1] daidt1=k4ai=F(ai). Vector: [t,g,o,cr]. d se obtiene por conservacion y ab y ac se obtienen dividiendo la criosfera con la funcion FunCrios

	zeroVector(k4Tsi,5);
	zeroVector(k4Tai,5);
	zeroVector(k4Wai,5);
	zeroVector(k4Wsi,3);

	for(i=0;i<=4;i++){
	  area = areai[i];
	  if(area>0.){
	    //;[K yr-1] dTsdt1=k4Ts=F(Ts). Es un vector: [t,g,d,o,cr]. BRetr esta en LHs
	    k4Tsi[i] = (1./cpsi[i])*(Rsi[i]-SHsi[i]-LHsi[i]+BRalbi[i]+DafmTsi[i]);	
	    //;[K yr-1] dTadt1=k4Ta=F(Ta). Es un vector: [t,g,d,o,cr]
	    k4Tai[i] = (1./cpa)*(Rai[i]+SHsi[i]+LHai[i]+FEai[i]+DafmTai[i]);
	    //;[kg m-2 yr-1] dWadt1=k4Wa=F(Wa). Es un vector: [t,g,d,o,cr]
	    k4Wai[i] = ETRi[i]-Preci[i]+FWai[i]+DafmWai[i];
	  }else{
	    k4Tsi[i] = 0.;
	    k4Tai[i] = 0.;
	    k4Wai[i] = 0.;
	  }
	}

	for(i=0;i<=2;i++){
	  area = areai[i];
	  if(area>0.){
	    //;[kg m-2 yr-1] dWsdt1=k4Ws=F(Ws). Vector: [t,g,d]
	    k4Wsi[i] = Preci[i]-ETRi[i]-Qsi[i]+DafmWsi[i];
	  }else{	
	    k4Wsi[i] = 0.;
	  }
	}

	printVector("k4Tsi = ",k4Tsi,5,1);
	printVector("k4Tai = ",k4Tai,5,1);
	printVector("k4Wai = ",k4Wai,5,1);
	printVector("k4Wsi = ",k4Wsi,3,1);


	//;********************************************************************************************
	//;FINALIZADO EL CALCULO DE k1,k2,k3 y k4 del metodo RK4 para todas las variables
	//;********************************************************************************************
	
	//;********************************************************************************************
	//;BLOQUE: actualizacion de variables de estado
	//;********************************************************************************************

	//;derivadas promedio segun RK4
	for(i=0;i<5;i++){
	  dTsidtRK4[i] = (k1Tsi[i]  + 2.*k2Tsi[i]  + 2.*k3Tsi[i]  + k4Tsi[i])/6.;
	  dTaidtRK4[i] = (k1Tai[i]  + 2.*k2Tai[i]  + 2.*k3Tai[i]  + k4Tai[i])/6.;
	  dWaidtRK4[i] = (k1Wai[i]  + 2.*k2Wai[i]  + 2.*k3Wai[i]  + k4Wai[i])/6.;
	  dareaidtRK4[i] = (k1areai[i]  + 2.*k2areai[i]  + 2.*k3areai[i]  + k4areai[i])/6.;
	  if(i<3) dWsidtRK4[i] = (k1Wsi[i]  + 2.*k2Wsi[i]  + 2.*k3Wsi[i] + k4Wsi[i])/6.;
	}

	if(VERBOSE or 0) printf("\n");
	printVector("dTsidt : ",dTsidtRK4,5,1);
	printVector("dTaidt : ",dTaidtRK4,5,1);
	printVector("dWaidt : ",dWaidtRK4,5,1);
	printVector("dWsidt : ",dWsidtRK4,3,1);
	printVector("dareaidt : ",dareaidtRK4,4,1);
	
	//;*********************************************
	//;Variables de estado actualizadas, es decir, en el tiempo t+deltatRK4
	for(i=0;i<=4;i++){
	  if(Tsi_t1[i]!=nodata){Tsi[i] = Tsi_t1[i] + dTsidtRK4[i]*deltatRK4;}else{Tsi[i] = nodata;}
	  if(Tai_t1[i]!=nodata){Tai[i] = Tai_t1[i] + dTaidtRK4[i]*deltatRK4;}else{Tai[i] = nodata;}
	  if(i<2){
	    if(Wsi_t1[i]!=nodata){Wsi[i] = max(0., Wsi_t1[i] + dWsidtRK4[i]*deltatRK4);}else{Wsi[i] = nodata;}
	  }
	  if(Wai_t1[i]!=nodata){Wai[i] = max(0., Wai_t1[i] + dWaidtRK4[i]*deltatRK4);}else{Wai[i] = nodata;}
	}
	at = max(0., at_t1 + dareaidtRK4[0]*deltatRK4);//;trees
	ag = max(0., ag_t1 + dareaidtRK4[1]*deltatRK4);		//;grasses
	ao = ao_t1 + dareaidtRK4[2]*deltatRK4;					//;ocean
	acr = max(0., acr_t1 + dareaidtRK4[3]*deltatRK4);	//;cryosphere
	ad = 1.-(at+ag+ao+acr);							//;bare land
	areai[0]=at;areai[1]=ag;areai[2]=ad;areai[3]=ao;areai[4]=acr;

	//;limitando la minima area significativa
	for(i=0;i<=4;i++){
	  if(areai[i]<amin){
	    areai[i] = 0.;
	    Tsi[i] = nodata;
	    Tai[i] = nodata;
	    if(i<=2){Wsi[i] = nodata;}
	    Wai[i] = nodata;
	  }
	}

	if(VERBOSE or 0) printf("\n");
	if(VERBOSE or 0) printf("Paso n: %d\n",n);
	printVector("Tsi : ",Tsi,5,1);
	printVector("Tai : ",Tai,5,1);
	printVector("Wai : ",Wai,5,1);
	printVector("Wsi : ",Wsi,3,1);
	printVector("areai : ",areai,5,1);
	
	//;control
	if(minVector(areai,5)<0.){
	  fprintf(stderr,"ERROR en areai, aparece valor negativo\n");
	  exit(1);
	}

	//;Division de la criosfera con las variables actualizadas (temperatura).
	acr = areai[4];
	if(acr>0.){
	  Tscr = Tsi[4];
	  Tacr = Tai[4];
	  tetacr = acos(1.-acr);
	  intecr = FunInte(0.,tetacr);
	  c2 = (1.-c1)*(1.-cos(tetacr))/intecr;
	  //;tiene incorporada la condicion de amin
	  FunCrios(Tscr,acr,tetacr,c1,c2,amin,crios);	
	  aac = crios[0];
	  aab = crios[1];
	}else{
	  aac = 0.;
	  aab = 0.;
	}
	if(aac>0.){
	  tetaac = acos(1.-aac);
	  inteac = FunInte(0.,tetaac);
	  //;ec (9) modelo criosfera en papel
	  Tsac = Tscr*c1 + Tscr*c2*inteac/(1. - cos(tetaac));	
	  //;ec (9) modelo criosfera en papel
	  Taac = Tacr*c1 + Tacr*c2*inteac/(1. - cos(tetaac));	
	}else{
	  inteac = 0.;		//;necesaria para zona de ablacion
	  Tsac = nodata;	//;no existe zona de acumulacion
	  Taac = nodata;
	}
	if(aab>0.){
	  if(aac==0.){tetaac = 0.;}//;si aac gt 0. tetaac ya fue calculada	
	  //;da lo mismo que inteab = FunInte(tetaac,tetacr)
	  inteab = intecr - inteac;
	  //;ec (9) modelo criosfera en papel
	  Tsab = Tscr*c1 + Tscr*c2*inteab/(cos(tetaac) - cos(tetacr));
	  //;ec (9) modelo criosfera en papel
	  Taab = Tacr*c1 + Tacr*c2*inteab/(cos(tetaac) - cos(tetacr));	
	}else{
	  Tsab = nodata;	//;no existe zona de ablacion
	  Taab = nodata;
	}
	criosaT[0]=aab;criosaT[1]=aac;criosaT[2]=Tsab;criosaT[3]=Tsac;criosaT[4]=Taab;
	criosaT[5]=Taac;//;actualizado

	//;Wvi actualizada (Wv = Wa-Wpre interviene en los ajustes Dafm)
	zeroVector(Wvi,5);

	for(i=0;i<=4;i++){
	  area = areai[i];
	  if(area>0.){
	    strcpy(DAF,DAFs[i]);		//;DAFs= ['t','g','d','o','cr']		
	    Ta = Tai[i];
	    Wa = Wai[i];
	    Ztoa = Ztoai[i];
	    atmos = FunAtm(Ta, Wa, lapse, Ztoa, DAF);	//;atmosfera, resultados = [Wv, Wpre, Wsat, RH, qv0, RH0, Ta0]
	    Wvi[i] = atmos[0];		//;[kg m-2] agua precipitable que permanece en la atmosfera
	  }else{
	    Wvi[i] = nodata;
	  }
	}

	//;particion del vapor de agua en la criosfera
	if(acr>0.){
	  Wvcr = Wvi[4];					//;Wacr = Wvcr+Wprecr (conservacion)
	  Ztoacr = Ztoai[4];
	  if(aab>0.){
	    atmos = FunAtm(Taab, Wacr, lapse, Ztoacr, DAF);	//;atmosfera, resultados = [Wv, Wpre, Wsat, RH, qv0, RH0, Ta0]
	    fac1 = atmos[0];
	  }else{
	    fac1 = 0.;
	  }
	  if(aac>0.){
	    atmos = FunAtm(Taac, Wacr, lapse, Ztoacr, DAF);	//;atmosfera, resultados = [Wv, Wpre, Wsat, RH, qv0, RH0, Ta0]
	    fac2 = atmos[0];
	  }else{
	    fac2 = 0.;
	  }
	  facab = fac1/(fac1+fac2);
	  facac = fac2/(fac1+fac2);
	  Wvab = facab*Wvcr;
	  Wvac = facac*Wvcr;
	}else{
	  Wvab = nodata;
	  Wvac = nodata;
	}

	//;Biotic Regulation actualizada
	Betat = FunBeta(Tsi[0], Toptt,Toptt0);
	tmort = FunTmor(Wsi[0], Wmi[0]);
	Betag = FunBeta(Tsi[1], Toptg,Toptg0);
	tmorg = FunTmor(Wsi[1], Wmi[1]);
	las = areai[2]/(areai[0]+areai[1]+areai[2]);

	for(i=0;i<=1;i++){
	  strcpy(DAF,DAFs[i]);
	  if(strcmp(DAF,"t")==0){
	    Topt=Toptt;
	    bs = las*Betat-tmort;
	  }else{
	    Topt=Toptg;
	    bs = las*Betag-tmorg;
	  }
	  BioReg = FunBioticReg(Tsi[i], Topt, ROCsi[i], Rsi[i], LHsi[i], Preci[i], ETRi[i], bs);
	  BRalb = BioReg[0];
	  BRetr = BioReg[1];
	  deltaETR = BioReg[2];
	  BR = BRalb + BRetr;	//;total
				//;ajuste de flujos
	  ETR = ETR + deltaETR;	//;deltaETR contiene el signo
	  LHs = LHs + BRetr;	//;BRetr contiene el signo
				//;el ajuste BRalb se aplica directamente en la ecuacion de BES
	  BRi[i] = BR;
	  //;//;BRalbi[i] = BRalb
	}

	//;verificacion de compatibilidad de Taab y Wvab (o Taac y Wvac). Pruebas de escritorio
	//;//;atmos = FunAtm(Taab, Wvab, lapse, Ztoacr, DAF)			
	//;//;Wvab2 = atmos(0)	//;da lo mismo que Wvab calculado antes
	//;//;atmos = FunAtm(Taab, Wvac, lapse, Ztoacr, DAF)			
	//;//;Wvac2 = atmos(0)	//;da lo mismo que Wvac calculado antes
	
	//;calculo de daidt que se usa en los ajustes Dafm de la proxima iteracion
	for(i=0;i<4;i++){
	  dadt1[i] = (areai_t1[i] - areai[i]);			//;[t,g,d,o] contraccion o expansion de las DAFs
	}
	dadt2[0] = (criosaT_t1[0] - criosaT[0]);//;[ab, ac]
	dadt2[1] = (criosaT_t1[1] - criosaT[1]);//;[ab, ac]
	
	dareaidtDafm[0] = dadt1[0]; 
	dareaidtDafm[1] = dadt1[1]; 
	dareaidtDafm[2] = dadt1[2]; 
	dareaidtDafm[3] = dadt1[3]; 
	dareaidtDafm[4] = dadt2[0]; 
	dareaidtDafm[5] = dadt2[1]; 

	TsiDafm[0] = Tsi[0];//;[t,g,d,o,ab,ac]
	TsiDafm[1] = Tsi[1];
	TsiDafm[2] = Tsi[2];
	TsiDafm[3] = Tsi[3];
	TsiDafm[4] = Tsab;
	TsiDafm[5] = Tsac;

	TaiDafm[0] = Tai[0];//;[t,g,d,o,ab,ac]
	TaiDafm[1] = Tai[1];
	TaiDafm[2] = Tai[2];
	TaiDafm[3] = Tai[3];
	TaiDafm[4] = Taab;
	TaiDafm[5] = Taac;

	WaiDafm[0] = Wvi[0];//;[t,g,d,o,ab,ac]
	WaiDafm[1] = Wvi[1];
	WaiDafm[2] = Wvi[2];
	WaiDafm[3] = Wvi[3];
	WaiDafm[4] = Wvab;
	WaiDafm[5] = Wvac;

	areaiDafm[0] = areai[0];//;[t,g,d,o,ab,ac]
	areaiDafm[1] = areai[1];
	areaiDafm[2] = areai[2];
	areaiDafm[3] = areai[3];
	areaiDafm[4] = aab;
	areaiDafm[5] = aac;

	if(VERBOSE or 0) printf("\n");
	printVector("dareaidtDafm = ",dareaidtDafm,6,1);
	printVector("areaiDafm = ",areaiDafm,6,1);
	printVector("TsiDafm = ",TsiDafm,6,1);
	printVector("TaiDafm = ",TaiDafm,6,1);
	printVector("WaiDafm = ",WaiDafm,6,1);

	//if(n>100) exit(0);
      }//End for n RK4

      printVector("dareaidtDafm = ",dareaidtDafm,5,2);
      printVector("areaiDafm = ",areaiDafm,5,2);
      printVector("TsiDafm = ",TsiDafm,5,2);
      printVector("TaiDafm = ",TaiDafm,5,2);
      printVector("WaiDafm = ",WaiDafm,5,2);
      
      fprintf(fl,"%5d ",t);
      fprintVector(fl,areaiDafm,5);
      fprintVector(fl,TsiDafm,5);
      fprintVector(fl,TaiDafm,5);
      fprintVector(fl,WaiDafm,5);
      fprintf(fl,"\n");
      fflush(fl);

    }//End for t
    
    break;
    
  }//End for L 

  return 0;
}
