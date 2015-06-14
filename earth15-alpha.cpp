#ifndef QCONSTANTS
#include<constants.h>
#endif
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>

#define expint gsl_sf_expint_E1

int fillRow(real** matrix,int irow,int ncol,real* nrow)
{
  for(int i=0;i<ncol;i++)
    matrix[irow][i]=nrow[i];
  return 0;
}


real aguaPrecipitable(int npart,real Ztoa,
		       real p0,real Ta0,real teta,real qv0,
		       real lapse,real Hv)
{
  real suma = 0.; //;sumador
  real deltaZ = Ztoa/npart;
  real z1,z2,Tz1,Tz2,Pz1,Pz2,qvz1,qvz2,qvmed,deltaP;

  for(int i=0;i<npart;i++){
    z1 = i*deltaZ;
    z2 = (i+1)*deltaZ;
    Tz1 = Ta0 - lapse*z1;
    Tz2 = Ta0 - lapse*z2;
    Pz1 = p0*pow((Tz1/Ta0),teta);
    Pz2 = p0*pow((Tz2/Ta0),teta);
    qvz1 = qv0*exp(-z1/Hv);
    qvz2 = qv0*exp(-z2/Hv);
    qvmed = (qvz1+qvz2)/2.;
    deltaP = Pz2-Pz1;
    suma = suma + qvmed*deltaP;
  }
  if(VERBOSE>1) printf("AGUAPRECIPITABLE: Suma = %.17e\n",suma);
  return suma;
}


real total(real *vec,int n)
{
  real suma=0.;
  for(int i=0;i<n;i++){
    suma+=vec[i];
  }
  return suma;
}

int printVector(const char *msg,real* vec,int n,int verbose=0)
{

  if(verbose==1) verbose=0;
  if(VERBOSE || verbose) printf("%s[",msg);
  for(int i=0;i<n;i++){
    if(VERBOSE || verbose) printf("%.17e ",vec[i]);
  }
  if(VERBOSE || verbose) printf("]\n");
  return 0;
}

int fprintVector(FILE *fl,real* vec,int n)
{
  for(int i=0;i<n;i++) fprintf(fl,"%.17e ",vec[i]);
  return 0;
}

int zeroVector(real* targ,int n)
{
  for(int i=0;i<n;i++) targ[i]=0.0;
  return 0;
}

int copyVector(real* targ,real* source,int n)
{
  for(int i=0;i<n;i++) targ[i]=source[i];
  return 0;
}

real** fltarr(int m,int n)
{
  real** flt=(real **)calloc(n,sizeof(real*));
  for(int i=0;i<n;i++){
    flt[i]=(real *)calloc(m,sizeof(real));
  }
  return flt;
}

real* fltarr1(int n)
{
  real* flt=(real *)calloc(n,sizeof(real));
  return flt;
}


real FunInte(real teta1,real teta2)
{
  real n=100.		;//n tiene que ser un numero PAR (even number)
  real a=teta1;
  real b=teta2;
  real h=(b-a)/n;

  ;//extremos
  real fa = sqrt(a)*sin(a);
  real fb = sqrt(b)*sin(b);

  ;//suma impares
  real sumimp = 0.;
  real x,fx;
  for(int j=1;j<=n-1;j+=2){
    x=a+j*h;
    fx=sqrt(x)*sin(x);
    sumimp = sumimp + fx;
  }

  ;//suma pares
  real sumpar = 0.;
  for(int j=2;j<=n-2;j+=2){
    x=a+j*h;
    fx=sqrt(x)*sin(x);
    sumpar = sumpar + fx;
  }
  real inte = (h/3.)*(fa + 4.*sumimp + 2.*sumpar + fb );
  return inte;
}

real* FunFEa(real areai[5],real Tai[5])
{
/*
;*************************************************************************************************************************
;flujos horizontales de energia en la atmosfera
;*************************************************************************************************************************
function FunFEa, areai, Tai

;JFSalazar, may 2010
;vectores en orden: [t,g,d,o,cr]
;convencion: si FEai > 0 la region gana energia, si FEai < 0 la region pierde energia
*/

  //;parametros
  real Rearth = 6370597.;	//;[m] radio de la Tierra, Re = sqrt(Ae/(4*pinum)) con Ae=5.10e14 m2
  real DTz = 15.8;		//;[MJ m-2 yr-1 K-1] coeficiente de difusion zonal de energia, 
                        // basado en Ledley 91 (ver excel)
			//;Nota: DTm y DTz no se entienden de la misma manera, de ahi la diferencia en el orden de magnitud.
			//;Consistentes con el orden de magnitud de los flujos en Fig. 6 de Ledley 1991.
			//;ver tambien fig 12.5 CW99 p.337
  //;ensayando
  real DTm = 7.5;	 //;el mismo pa no afectar la criosfera


  //;constantes numericas
  real ceronum = 1.e-2;		//;cero numerico. Si no se fija puede haber problemas en la verificacion de la conservacion de energia
  real nodata = -9999.;

  real at = areai[0];
  real ag = areai[1];
  real ad = areai[2];
  real ao = areai[3];
  real acr = areai[4];

  real Tt = Tai[0];
  real Tg = Tai[1];
  real Td = Tai[2];
  real To = Tai[3];
  real Tcr = Tai[4];

  //;medidas de las areas a traves de las cuales se dan los flujos meridionales hacia/desde la criosfera
  //L*: Longitud del borde de la criosfera y el oceano
  real tetacr = acos(1.-acr);
  real Lcr = sin(tetacr);		//;si acr=0. entonces Lcr=0. (no hay intercambios con la criosfera porque esta no existe)
  
  real Locr,Ldcr,Lgcr,Ltcr;

  if(ad>0.){
    Locr = Lcr*ao/(ao+ad);
    Ldcr = Lcr*ad/(ao+ad);
  }else{
    Locr = Lcr*ao/(ao+at+ag);
    Ltcr = Lcr*at/(ao+at+ag);
    Lgcr = Lcr*ag/(ao+at+ag);
  }
  
  //;*************************************************************************************************************************
  
  double aref,Eto,Etg,Etn,FEt;
  //;trees
  if(at>0){
    //;ocean
    aref=min(ao,at);
    Eto=DTz*(To-Tt)*aref;
    //;grass
    aref = min(ag,at);
    //;notese que si aref=0. el flujo es 0. aref es una medida del area a traves de la cual se da la transferencia
    Etg = DTz*(Tg-Tt)*aref;		
    //;desert or cryosphere	
    if(ad>0){
      aref = min(ad,at);
      Etn = DTz*(Td-Tt)*aref;
    }else{
      //;el intercambio es con la criosfera y es meridional, no zonal
      Etn = DTm*(Tcr-Tt)*Ltcr;		//;Ltcr=0. si acr=0.
    }
    //;total
    FEt = Eto + Etg + Etn;	//;flujo. Para obtener la convergencia en t basta dividir por at (esto se hace al final)
  }else{
    FEt = 0.;
  }


  //;*************************************************************************************************************************
  //;grass
  double Ego,Egt,Egn,FEg;
  if(ag>0.){
    //;ocean
    aref = min(ao,ag);
    Ego = DTz*(To-Tg)*aref;
    //;trees
    aref = min(at,ag);
    Egt = DTz*(Tt-Tg)*aref;
    //;desert or cryosphere	
    if(ad>0.){
      aref = min(ad,ag);
      Egn = DTz*(Td-Tg)*aref;
    }else{
      //;el intercambio es con la criosfera y es meridional, no zonal
      Egn = DTm*(Tcr-Tg)*Lgcr;
    }
    //;total
    FEg = Ego + Egt + Egn;	//;flujo. Para obtener la convergencia en g basta dividir por ag (esto se hace al final)
  }else{
    FEg = 0.;
  }
  
  //;*************************************************************************************************************************
  //;desert
  real Edo,Edt,Edg,Edcr,FEd;
  if(ad>0.){
    //;ocean (flujo zonal)
    aref = min(ao,ad);
    Edo = DTz*(To-Td)*aref;
    //;trees  (flujo zonal)
    aref = min(at,ad);
    Edt = DTz*(Tt-Td)*aref;
    //;grass (flujo zonal)
    aref = min(ag,ad);
    Edg = DTz*(Tg-Td)*aref;
    //;cryosphere (flujo meridional)
    Edcr = DTm*(Tcr-Td)*Ldcr;
    //;total
    FEd = Edo + Edt + Edg + Edcr;
  }else{
    FEd = 0.;
  }

  //;*************************************************************************************************************************
  //;ocean (ao > 0. siempre)
  //;trees (zonal)
  real Eot,Eog,Eod,Eocr,FEo;

  aref = min(at,ao);
  Eot = DTz*(Tt-To)*aref;
  //;grass (zonal)
  aref = min(ag,ao);
  Eog = DTz*(Tg-To)*aref;
  //;desert (zonal)
  aref = min(ad,ao);
  Eod = DTz*(Td-To)*aref;
  //;cryosphere (meridional)
  Eocr = DTm*(Tcr-To)*Locr;			
  //;total
  FEo = Eot + Eog + Eod + Eocr;

  //;*************************************************************************************************************************
  //;cryosphere (todos son flujos meridionales)
  real Ecro,Ecrt,Ecrg,Ecrn,FEcr;
  if(acr>0.){
    //;ocean
    Ecro = DTm*(To-Tcr)*Locr;
    if(ad>0.){
      Ecrn = DTm*(Td-Tcr)*Ldcr;
    }else{
      //;trees
      Ecrt = DTm*(Tt-Tcr)*Ltcr;
      //;grass
      Ecrg = DTm*(Tg-Tcr)*Lgcr;
      //;total from/to land
      Ecrn = Ecrt+Ecrg;
    }
    FEcr = Ecro + Ecrn;
  }else{
    FEcr = 0.;
  }
  
  //;*************************************************************************************************************************
  //;resultados
  real* FEai=fltarr1(5);
  FEai[0]=FEt;
  FEai[1]=FEg; 
  FEai[2]=FEd;
  FEai[3]=FEo;
  FEai[4]=FEcr;

  //;conservacion de la energia (aplica sobre los flujos, no sobre la convergencia)
  int control = fabs(total(FEai,5));
  if(control>ceronum){
    fprintf(stderr,"ERROR en FunDeltaEa: se viola la conservacion de la energia\n");
    exit(1);
  }

  //;Convirtiendo los flujos a convergencia sobre cada region.
  //;Convergencia/divergencia total sobre cada region (por unidad de area de cada region)
  FEai[0]=areai[0]>0?FEai[0]/areai[0]:nodata;
  FEai[1]=areai[1]>0?FEai[1]/areai[1]:nodata;
  FEai[2]=areai[2]>0?FEai[2]/areai[2]:nodata;
  FEai[3]=areai[3]>0?FEai[3]/areai[3]:nodata;
  FEai[4]=areai[4]>0?FEai[4]/areai[4]:nodata;

  return FEai;
}


/*

*/

real* FunAtm(real Ta,real Wa,real lapse,real Ztoa,char DAF[]){

  real g = 9.81;			//; [m s-2] aceleracion de la gravedad - REVISADO
  real Rv = 461.0;			//; [J kg-1 K-1]	specific gas constant for vapour (A00, p.36,216, CW99 p.438) - REVISADO
  real R = 287.0;			//; [J kg-1 K-1]	specific gas constant for dry air (A00, p.51,216 CW99 p.437) - REVISADO
  real epsilon = R/Rv;	  //; [adim] aprox. 0.622 cociente (masa molecular de vapor de agua)/(masa molecular promedio de aire humedo o seco) 
                          //  (A00, p.37) - REVISADO
  real etp = 611.;	  //; [Pa] presion de vapor de referencia (Clausius-Clapeyron), triple punto (CW99, p.102) - REVISADO
  real Ttp = 273.16;	  //; [K] temperatura de referencia (Clausius-Clapeyron), triple punto (CW99, p.102) - REVISADO
  real Lv = 2.501e6;	  //; [J kg-1] calor latente de vaporizacion/condensacion (gas/liquido) - REVISADO
  real Ld = 0.333e6;	  //; [J kg-1] calor latente de congelacion/fusion (solido/liquido) (A00p216)  - REVISADO 
  real Ls = 2.834e6;	  //; [J kg-1] calor latente de sublimacion/ (gas/solido) Ls=Lv+Ld (A00p.216) - REVISADO
  real p0 = 1.013e5;	  //; [Pa] presion atmosferica en la superficie (CW99 p.438) - REVISADO
  real Hv = 2000.;	  //; [m] water scale height (Pet00)

  //;*******************************************************************************************************************
  //;Parametros (constantes a lo largo del tiempo)
  //;//;//;RH0i = [0.90, 0.60, 0.30, 1.00, 1.00]	//;humedad relativa en la superficie para ['t','g','d','o','cr'], valor de referencia
  real RH0i[] = {1.00, 0.80, 0.60, 1.00, 1.00};	//;humedad relativa en la superficie para ['t','g','d','o','cr'], valor de referencia
  real RHcrit = 0.85;		//;[adim] humedad reativa critica (de toda la columna) por encima de la cual todo se vuelve 
					//;precipitacion (BjornsonEt97=0.85, Pet00=0.95) - REVISADO
  real* atmos=fltarr1(7);
  

  //;*******************************************************************************************************************
  //;Constantes fisicas
  //;*******************************************************************************************************************
  //;Constantes numericas
  int npart = 100;		//;numero de particiones de la atmosfera para obtener el agua precipitable

  //;*******************************************************************************************************************
  //;BLOQUE 1: Calculos iniciales
  //;*******************************************************************************************************************
  
  //;Temperatura del aire en la superficie (Tao)
  real Ta0 = Ta + lapse*Ztoa/2.;	

  //;humedad de saturacion en la columna (Wsat)
  real teta = g/(lapse*R);
  real e1 = etp*exp(Lv/(Rv*Ttp));	//;Pa
  real x1 = Lv/(Rv*Ta0);
  real x2 = Lv/(Rv*(Ta0-lapse*Ztoa));
  real Wsat =  (e1/(Rv*lapse))*( expint(x1) - expint(x2) );
  //;interesante deduccion, ver en documento tesis
  
  //;humedad relativa en la superficie (RHo)
  real RH0;
  if(strcmp(DAF,"t")==0){RH0=RH0i[0];}
  if(strcmp(DAF,"g")==0){RH0=RH0i[1];}
  if(strcmp(DAF,"d")==0){RH0=RH0i[2];}
  if(strcmp(DAF,"o")==0){RH0=RH0i[3];}
  if(strcmp(DAF,"cr")==0){RH0=RH0i[4];}


  //;humedad especifica en la superficie (qv0 y qvsat0)
  real es0 = etp*exp((Lv/Rv)*(1./Ttp - 1./Ta0));
  real qvsat0 = epsilon*es0/p0;	//;de saturacion
  real qv0 = qvsat0*RH0;

  if(VERBOSE>1) printf("FUNATM: qv0 = %.17e\n",qv0);
  //qv0=0.9;
  //Wa=1045.0;
  

  //;*******************************************************************************************************************
  //;BLOQUE 2: calculo del agua precipitable (integracion numerica)
  //;calculo del agua precipitable en la columna atmosferica a partir de suponer una funcion qv=qv(Z)=qv0*exp(-Z/Hv)
  //;*******************************************************************************************************************
  //;la que habria con el perfil de humedad supuesto
  real Wv = -aguaPrecipitable(npart,Ztoa,p0,Ta0,teta,qv0,lapse,Hv)/g;	

  if(VERBOSE>1) printf("FUNATM: suma,Wa,Wv = %.17e,%.17e,%.17e\n",-aguaPrecipitable(npart,Ztoa,p0,Ta0,teta,qv0,lapse,Hv),Wa,Wv);
  //exit(0);

  //;*******************************************************************************************************************
  //;BLOQUE 3: ajuste por disponibilidad de agua (Wa)
  //;*******************************************************************************************************************

  real ak,bk,ck,fak,fbk,fck,suma,cond1;
  if(Wa<Wv){
    //;en este caso hay que ajustar el valor de qv0, y entonces el de RH0, para que Wv=Wa. El nuevo qv0 se haya como
    //;una raiz de f(qv0) = Wa-Wv que se sabe esta entre 0 y qv0inicial.

    //;**************************************************
    //;Inicio: Metodo de la falsa posicion (Regula Falsi)  1a-vez	

    //;intervalo inicial que contiene la raiz (siempre)
    ak = 0.0;		//;f(0.) = Wa > 0.
    bk = qv0;		//;f(qv0,inicial) = Wa-Wv < 0.
    
    //;inicio de las iteraciones
    real crite = 10.;			//;cualquier valor mas grande que 0.01 (la tolerancia fijada) para iniciar el while
    real numit = 0.;			//;sumador para prevenir loop infinito
    real deltaZ = Ztoa/npart;	        //;tamano de paso vertical

    while(crite>0.01){		//;se fija la maxima diferencia abs(Wa-Wv) permitida en 0.01

      //;paso 1: evaluar funcion en los extremos del intervalo
      qv0 = ak;
      fak=Wa-(-aguaPrecipitable(npart,Ztoa,p0,Ta0,teta,qv0,lapse,Hv)/g);

      qv0 = bk;
      fbk=Wa-(-aguaPrecipitable(npart,Ztoa,p0,Ta0,teta,qv0,lapse,Hv)/g);

      cond1 = fak*fbk;
      if(cond1>0.){
	if(VERBOSE>1) printf("FUNATM: ERROR3 en FunAtm: fak*fbk > 0 resolviendo f(qv0) = Wa-Wv = 0\n");
	exit(0);
      }

      //;paso 2: estimacion de falsa posicion
      ck=(fbk*ak-fak*bk)/(fbk-fak);	//;formula falsa posicion

      qv0 = ck;
      fck=Wa-(-aguaPrecipitable(npart,Ztoa,p0,Ta0,teta,qv0,lapse,Hv)/g);
      if(VERBOSE>1) printf("%f,%f,%f,%f,%f\n",aguaPrecipitable(npart,Ztoa,p0,Ta0,teta,qv0,lapse,Hv),g,fck,fbk*ak-fak*bk,Wa);

      if(VERBOSE>1) printf("fak, fbk, ck, fck = %f,%f,%f,%f\n",fak,fbk,ck,fck);

      //;paso 3: verificar criterio
      crite=fabs(fck);
      //;se fija la maxima diferencia abs(Wa-Wv) permitida en 0.01		
      if(crite<=0.01){				
	//;//;//;qv0 = qv0	//;ya se encontro el valor de qv0, es el ultimo hallado
	Wv = Wa;
	RH0 = qv0/qvsat0;	//;nueva humedad especifica en la superficie
      }else{
	if(fak*fck<0.){
	  bk=ck;
	}else{
	  if(fck*fbk<0.){
	    ak = ck;
	  }else{
	    if(VERBOSE>1) printf("ERROR4 en FunAtm: la raiz (qv0) no esta en el intervalo\n");	
	    exit(1);
	  }
	}
      }
      numit = numit + 1;
      if(numit>1000){
	if(VERBOSE>1) printf("ERROR5 en FunAtm: No encontro raiz (qv0) en 1000 iteraciones\n");
	exit(0);
      }
    }
    if(VERBOSE>1) printf("FUNATM: qv0,RH0=%f,%f\n",qv0,RH0);
  }

  if(VERBOSE>1) printf("FUNATM: Wa,Wv,Wsat,Wv/Wsat=%f,%f,%f,%f,%f\n",Wa,Wv,Wsat,Wv/Wsat,RHcrit);
  real RH,Wpre;

  //;************************************************
  //;BLOQUE 4: calculo de la precipitacion y la humedad relativa de la columna. Ajuste de la atmosfera 
  //;si excede RHcrit.
  //;//;jumpA: //;no hace falta

  if(Wv/Wsat<=RHcrit){
    RH = Wv/Wsat;	//;de toda la columna
    Wpre = Wa-Wv;
  }else{
    //;ajuste de la atmosfera liberando agua (precipitacion) hasta reducir RH hasta RHcrit
    RH = RHcrit;
    Wv = RHcrit*Wsat;		//;este va a ser el nuevo valor de Wv
    Wpre = Wa-Wv;
    //;como Wv se redujo entonces, consecuentemente se debe reducir qv0 y RH0
    //;en este caso hay que ajustar el valor de qv0, y entonces el de RH0, para que Wv=RHcrit*Wsat. El nuevo qv0 se haya como
    //;una raiz de f(qv0) = Wv - Wv* donde Wv* es el valor de Wv que se obtiene con qv0 dado. Se sabe que la raiz esta entre 
    //;0 y qv0inicial.
    
    //;**************************************************
    //;Inicio: Metodo de la falsa posicion (Regula Falsi)  2da-vez	
    
    //;intervalo inicial que contiene la raiz (siempre)
    ak = 0.0;		//;f(0.) = Wv > 0. En este caso Wv se trata como una constante, lo que varia con qv0 es Wv*
    bk = qv0;		//;f(qv0,inicial) = Wv-Wv* < 0 porque se esta tratando el caso Wv<Wv*.

    //;inicio de las iteraciones
    real crite = 10.;			//;cualquier valor mas grande que 0.01 (la tolerancia fijada) para iniciar el while
    int numit = 0.;			//;sumador para prevenir loop infinito
    real deltaZ = Ztoa/npart;	//;tamano de paso vertical

    while(crite>0.01){		//;se fija la maxima diferencia abs(Wa-Wv) permitida en 0.01

      //;paso 1: evaluar funcion en los extremos del intervalo
      qv0 = ak;
      fak=Wv-(-aguaPrecipitable(npart,Ztoa,p0,Ta0,teta,qv0,lapse,Hv)/g);

      if(VERBOSE>1) printf("FUNATM: Suma = %e\n",-aguaPrecipitable(npart,Ztoa,p0,Ta0,teta,qv0,lapse,Hv)/g);

      qv0 = bk;
      fbk=Wv-(-aguaPrecipitable(npart,Ztoa,p0,Ta0,teta,qv0,lapse,Hv)/g);

      cond1 = fak*fbk;
      if(VERBOSE>1) printf("FUNATM: ak = %e, bk = %e, fak, fbk = %e, %e\n",ak,bk,fak,fbk);

      if(cond1>0.){
	if(VERBOSE>1) printf("FUNATM: ERROR3 en FunAtm: fak*fbk > 0 resolviendo f(qv0) = Wa-Wv = 0\n");
	exit(0);
      }

      //;paso 2: estimacion de falsa posicion
      ck=(fbk*ak-fak*bk)/(fbk-fak);	//;formula falsa posicion

      qv0 = ck;
      fck=Wv-(-aguaPrecipitable(npart,Ztoa,p0,Ta0,teta,qv0,lapse,Hv)/g);
      if(VERBOSE>1) printf("%f,%f,%f,%f,%f\n",aguaPrecipitable(npart,Ztoa,p0,Ta0,teta,qv0,lapse,Hv),g,fck,fbk*ak-fak*bk,Wa);

      if(VERBOSE>1) printf("fak, fbk, ck, fck = %f,%f,%f,%f\n",fak,fbk,ck,fck);

      //;paso 3: verificar criterio
      crite=fabs(fck);
      //;se fija la maxima diferencia abs(Wa-Wv) permitida en 0.01		
      if(crite<=0.01){				
	//;//;//;qv0 = qv0	//;ya se encontro el valor de qv0, es el ultimo hallado
	RH0 = qv0/qvsat0;	//;nueva humedad especifica en la superficie
      }else{
	if(fak*fck<0.){
	  bk=ck;
	}else{
	  if(fck*fbk<0.){
	    ak = ck;
	  }else{
	    if(VERBOSE>1) printf("ERROR4 en FunAtm: la raiz (qv0) no esta en el intervalo\n");	
	    exit(1);
	  }
	}
      }
      numit = numit + 1;
      if(numit>1000){
	if(VERBOSE>1) printf("ERROR5 en FunAtm: No encontro raiz (qv0) en 1000 iteraciones\n");
	exit(0);
      }
    }
    if(VERBOSE>1) printf("qv0,RH0=%f,%f\n",qv0,RH0);
  }

  if(VERBOSE>1) printf("FUNATM: Wa,Wv,Wsat=%f,%f,%f,%f\n",Wa,Wv,Wsat,Wv/Wsat);
  
  real resultados[] = {Wv, Wpre, Wsat, RH, qv0, RH0, Ta0};
  copyVector(atmos,resultados,7);

  return atmos;
}


real* FunCloud(real Prec, real RH, real Ztoa, real lapse, real qv0, real Ta0,char  DAF[]){

  real* clouds=fltarr1(4);

  //;Constantes fisicas
  real g = 9.81;			//; [m2 s-1] aceleracion de la gravedad - REVISADO
  real Rv = 461.0;			//; [J kg-1 K-1]	specific gas constant for vapour (A00, p.36,216, CW99 p.438) - REVISADO
  real R = 287.0;			//; [J kg-1 K-1]	specific gas constant for dry air (A00, p.51,216 CW99 p.437) - REVISADO
  real epsilon = R/Rv;		//; [adim] aprox. 0.622 cociente (masa molecular de vapor de agua)/(masa molecular promedio de aire humedo o seco) (A00, p.37) - REVISADO
  real etp = 611.;			//; [Pa] presion de vapor de referencia (Clausius-Clapeyron), triple punto (CW99, p.102) - REVISADO
  real Ttp = 273.16;		//; [K] temperatura de referencia (Clausius-Clapeyron), triple punto (CW99, p.102) - REVISADO
  real Lv = 2.501e6;		//; [J kg-1] calor latente de vaporizacion/condensacion (gas/liquido) - REVISADO
  real Ld = 0.333e6;		//; [J kg-1] calor latente de congelacion/fusion (solido/liquido) (A00p216)  - REVISADO 
  real Ls = 2.834e6;		//; [J kg-1] calor latente de sublimacion/ (gas/solido) Ls=Lv+Ld (A00p.216) - REVISADO
  real p0 = 1.013e5;		//; [Pa] presion atmosferica en la superficie (CW99 p.438) - REVISADO
  real teta = g/(lapse*R);

  //;Constantes numericas
  real nodata = -9999.;

  //;Parametros de las nubes. En orden: [t,g,d,o,cr]
  real eCloi[] = {0.1, 0.1, 0.1, 0.1, 0.1};		        //; [adim] exponente ley potencial entre nubosidad y prec, NGC05 
  real pCloi[] = {0.40, 0.30, 0.20, 0.36, 0.20};			//; [adim] prefactor de la ley potencial entre nubosidad y prec
  real mAlbc = 0.12;				//;[adim] pendiente de la linea albedo=f(altitud)
  real bAlbc = 0.01;				//;[adim] "intercepto" de la linea albedo=f(altitud)
  real Kfstrocean = 0.87;			//;parametros en la particion de las nubes en str y cum
  real Kfstrland = 0.64;			//;notese que son del orden de magnitud de Pet00 (0.80 + algo)
  real Ptopland = 49000.;			//;[Pa] tope de las nubes sobre land (Warren&Hahn03 = WH03, table 2)
  real Ptopocean = 62000.;			//;[Pa] tope de las nubes sobre ocean (Warren&Hahn03 = WH03, table 2)
  real Ptopcryos = 62000.;			//;[Pa] tope de las nubes sobre cryosphere, se asumio igual que over ocean
  //;//;//;kEmis = 20.					//;parametro de la funcion emisividad=f(RH)

  //;Calculos
  real eClo,pClo,Ptop,Ztop;
  if(strcmp(DAF,"t")==0){
    eClo = eCloi[0];
    pClo = pCloi[0];
    Ptop = Ptopland;
  }
  if(strcmp(DAF,"g")==0){
    eClo = eCloi[1];
    pClo = pCloi[1];
    Ptop = Ptopland;
  }
  if(strcmp(DAF,"d")==0){
    eClo = eCloi[2];
    pClo = pCloi[2];
    Ptop = Ptopland;
  }
  if(strcmp(DAF,"o")==0){
    eClo = eCloi[3];
    pClo = pCloi[3];
    Ptop = Ptopocean;
  }
  if(strcmp(DAF,"cr")==0){
    eClo = eCloi[4];
    pClo = pCloi[4];
    Ptop = Ptopcryos;
  }

  //;*******************************************************************************************************************
  //;altura del tope de las nubes (Ztop) 
  //;la suposicion inicial es que:
  //;el tope de las nubes se encuentra siempre a una misma presion asi: land: 490 hPa, ocean: 620 hPa, ice: 490 hPa (WH03)
  Ztop = (Ta0/lapse)*(1.- pow((Ptop/p0),(1./teta)));
  if(Ztop>Ztoa){Ztop = Ztoa;}	//;acotando
  if(VERBOSE>1) printf("Ptop,Ztop = %e,%e\n",Ptop,Ztop);
  
  real Zlcl,ak,bk,es0,fak,z,Tz,esz,pz,qvsat0,qvsatz,fbk,cond1;
  real ck,crite,fck;
  int numit;
  if(strcmp(DAF,"cr")==0){
    Zlcl = nodata;		//;en cr no hace falta calcularlo
  }else{

    //;Zlcl se halla como la raiz de la funcion f(z)=qv0/qvsat(z) - 1.00. O sea que f(Zlcl) = 0.

    //;**************************************************
    //;Inicio: Metodo de la falsa posicion (Regula Falsi)
    
    //;intervalo inicial. Inicialmente se supone el maximo posible que es (0., Ztop)
    ak = 0.0;
    bk = Ztop;
      
    //;vericando si existe la raiz
    //;ak=0.0
    es0 = etp*exp((Lv/Rv)*(1./Ttp - 1./Ta0));
    qvsat0 = epsilon*es0/p0;		//;humedad de saturacion en la superficie
    fak = qv0/qvsat0 - 1.00;

    z = bk;
    Tz = Ta0 - lapse*z;
    esz = etp*exp((Lv/Rv)*(1./Ttp - 1./Tz));
    pz = p0*pow((Tz/Ta0),teta);
    qvsatz = epsilon*esz/pz;		//;humedad de saturacion a la altura z
    fbk = qv0/qvsatz - 1.00;

    cond1 = fak*fbk;
    if(cond1 > 0.){
      //;no existe raiz
      Zlcl = nodata;
    }else{
      //;si existe raiz
      //;inicio de las iteraciones
      crite = 10.;			//;cualquier valor mas grande que 0.001 (la tolerancia fijada) para iniciar el while
      numit = 0.;			//;sumador para prevenir loop infinito

      while(crite>0.01){//;se fija la maxima diferencia abs(qv0/qvsat(z) - 1.00.) permitida en 0.001
	//;tenia 0.001 pero por razones practicas (costo computacional) 0.01 suficiente
	//;ojo, esta tolerancia tambien esta abajo donde se verifica

	//;paso 1: evaluar funcion en los extremos del intervalo
	z = ak;
	Tz = Ta0 - lapse*z;
	esz = etp*exp((Lv/Rv)*(1./Ttp - 1./Tz));
	pz = p0*pow((Tz/Ta0),teta);
	qvsatz = epsilon*esz/pz;		//;humedad de saturacion a la altura z
	fak = qv0/qvsatz - 1.00;

	z = bk;
	Tz = Ta0 - lapse*z;
	esz = etp*exp((Lv/Rv)*(1./Ttp - 1./Tz));
	pz = p0*pow((Tz/Ta0),teta);
	qvsatz = epsilon*esz/pz;		//;humedad de saturacion a la altura z
	fbk = qv0/qvsatz - 1.00;


	//;paso 2: estimacion de falsa posicion
	ck = (fbk*ak - fak*bk)/(fbk - fak);	//;formula falsa posicion
	z = ck;
	Tz = Ta0 - lapse*z;
	esz = etp*exp((Lv/Rv)*(1./Ttp - 1./Tz));
	pz = p0*pow((Tz/Ta0),teta);
	qvsatz = epsilon*esz/pz;		//;humedad de saturacion a la altura z
	fck = qv0/qvsatz - 1.00;

	if(VERBOSE>1) printf("fak,fbk,fck=%e,%e,%e\n",fak,fbk,fck);

	//;paso 3: verificar criterio
	crite = fabs(fck);
	if(crite<=0.01){
	  Zlcl = z;
	    //;paso 4: ajustar el intervalo si no se cumple el criterio	
	}else{
	  if(fak*fck<0.){
	    bk=ck;
	  }else{
	    if(fck*fbk<0.){
	      ak = ck;
	    }else{
	      fprintf(stderr,"ERROR en FunCloud: la raiz (Zlcl) no esta en el intervalo\n");
	      exit(1);
	    }
	  }
	}
	numit = numit + 1;
	if(numit>2000){
	  fprintf(stderr,"ERROR en FunCloud: No encontro raiz (qv0) en 2000 iteraciones\n");
	}
      }
      if(VERBOSE>1) printf("numit = %d, ck = %e, Zlcl = %e\n",numit,ck,Zlcl);
      //;Fin: Metodo de la falsa posicion (Regula Falsi)
      //;*************************************************
    }
  }
    
  //;*******************************************************************************************************************
  //;area total de nubes y particion en stratiform and cumuliform
  //;ac = fstr*ac + fcum*ac, donde fstr+fcum = 1.0
  
  real ac;
  ac = pClo*pow(Prec,eClo);	//;total
  if(ac>1.0){ac=1.0;}	//;acotamiento
  
  real fstr,Kfstr;
  if(strcmp(DAF,"cr")==0){
    fstr=1.0; 
  }else{
    if(Zlcl==nodata){
      //;si no hay Zlcl entonces no hay cumuliform
      fstr=1.0 ;
    }else{
      //Relation between relative humidity and startiform clouds
      if(strcmp(DAF,"o")==0){Kfstr=Kfstrocean;}
      else{Kfstr=Kfstrland;}
      fstr=Kfstr*pow(RH,1.5);	//;forma basada en Pet00
    }
  }

  real fcum,acstr,accum;
  fcum = 1.-fstr;
  acstr = ac*fstr;		//;stratiform
  accum = ac*fcum;		//;cumuliform

  if(VERBOSE>1) printf("str,cum = %e,%e\n",acstr,accum);

  //;altura y albedo de las nubes segun el tipo

  real albstr,Zstr,Zcum,albcum;
  if(acstr>0){
    Zstr = Ztop;
    albstr = bAlbc + mAlbc*Zstr/1000.;		//; /1000. porque Zstr aqui tiene que estar en km
    if(albstr>0.85){albstr = 0.85;}			//;acotamiento
    if(albstr<0.25){albstr = 0.25;}
  }else{
    Zstr = nodata;
    albstr = nodata;
  }
  
  if(accum>0.){
    Zcum = (Zlcl+Ztop)/2.;	//;asi el albedo depende del desarrollo vertical
    //;esto hace que sean mas bajas que Zstr y por ende tengan menor albedo.
    albcum = bAlbc + mAlbc*Zcum/1000.;		//; /1000. porque Zcum aqui tiene que estar en km
    if(albcum>0.85){albcum = 0.85;}	//;acotamiento
    if(albcum<0.25){albcum = 0.25;}
  }else{
    Zcum = nodata;
    albcum = nodata;
  }
  if(VERBOSE>1) printf("Albedos str,cum = %e,%e\n",albstr,albcum);

  //;valores promedio regionales (en toda la DAF considerando las nubes como un todo: ac)
  real albc,Zc,Tc,emc;
  if(ac>0.){
    //;si algun albedo es nodata entonces se anual al multiplicar por cero, no problem
    albc=(albstr*acstr+albcum*accum)/ac;	
    Zc =(Zstr*acstr+Zcum*accum)/ac;
    Tc = Ta0 - lapse*Zc;
    emc = 1.0;	//;valor fijo
  }else{
    albc=nodata;
    //;en los balances de energia queda multiplicado por ac=0 y por ende se anula
    Tc=nodata;
    emc=nodata;
  }

  clouds[0]=ac;
  clouds[1]=albc;
  clouds[2]=Tc;
  clouds[3]=emc;

  if(VERBOSE>1) printf("ac,albc,Tc,emc = %e,%e,%e,%e\n",ac,albc,Tc,emc);
  return clouds;
}

/*
*/

real* FunRs(real S,real L,real Ts,real Ta,real RH,real* clouds,real* criosaT,char* DAF)
{
  //;JFSalazar, may 2010
  //;Rs en [MJ m-2 yr-1] donde [MJ] = [10^6 J]
  
  //;constantes fisicas
  //; [MJ m-2 K-4 yr-1] contante de Stefan-Boltzmann. Igual que 5.67e-8 W m-2 K-4, CW99 p.437 -REVISADO
  real* resultados=fltarr1(2);

  real sigma= 1.79e-6;			
    
  //;De las nubes
  real ac = clouds[0];		//;area total de nubes
  real albc = clouds[1];	//;albedo promedio de las nubes
  real Tc = clouds[2];		//;temperatura promedio de las nubes
  real emc = clouds[3];		//;emisividad promedio de las nubes

  //;De la atmosfera
  real aba = 0.23;			//;[adim] absorcion de ROC, aba=67/342 (TrenberthEt09) -REVISADO
  real emg = 0.25;			//;[adim] efecto invernadero de gases diferentes al vapor de agua (tambien esta en Ra)

  //;De la superficie
  real albcr=0;
  if(strcmp(DAF,"cr")==0)
    albcr = (criosaT[0]*0.40 + criosaT[1]*0.80)/( criosaT[0]+criosaT[1]);

  if(VERBOSE>1) printf("albcr = %e\n",albcr);
  
  real albsi[] = {0.15, 0.25, 0.20, 0.10};	//;[adim] albedo de la superficie, MH, en orden: [t,g,d,o] Ojo: no incluye criosfera
  real albs;
  if(strcmp(DAF,"t")==0) albs = albsi[0];
  if(strcmp(DAF,"g")==0) albs = albsi[1];
  if(strcmp(DAF,"d")==0) albs = albsi[2];
  if(strcmp(DAF,"o")==0) albs = albsi[3];
  if(strcmp(DAF,"cr")==0) albs = albcr;

  real Tae = (Ta + ac*Tc)/(1.+ac);

  if(VERBOSE>1) printf("Tc=%e,Tae=%e\n",Tc,Tae);
  
  real emv,ema;
  if(strcmp(DAF,"cr")==0){
    emv=0.0;
  }else{
    //;debido al bajo contenido de humedad atmosferica en la criosfera. Tambien en FunRs
    //;d no hacerlo asi el efecto invernadero de vapor de agua en cr seria muy alto, irr
    emv = 1.-exp(-RH);
  }
							
  ema = emv*(1.-ac) + emc*ac + emg;
  ema = max(1.0,ema);
  
  real ROCs = S*L*(1.-aba)*(1.-albs)*(1.-ac*albc);
  real Rs = ROCs - sigma*gsl_pow_int(Ts,4) + ema*sigma*gsl_pow_int(Tae,4);

  if(VERBOSE>1) printf("Rs = %e, ROCs = %e\n",Rs,ROCs);
  
  resultados[0]=Rs;
  resultados[1]=ROCs;

  return resultados;
}



/*
//;*************************************************************************************************************************
//;Calcula el flujo de calor sensible desde la superficie
//;*************************************************************************************************************************
function FunSH, Ts, Ta0, DAF
*/

real FunSH(real Ts,real Ta0,char* DAF)
{
  //;Valor revisado (del orden de magnitud de 10 W m-2 K-1 como en BjornssonEt97)
  real Hs = 315.;	//;[MJ m-2 yr-1 K-1]
  real SHs = Hs*(Ts-Ta0);	//;dirigido desde la superficie hacia la atmosfera si Ts>Ta0

  return SHs;
}



/*
//;JFSalazar, may 2010
//;SHs en [MJ m-2 yr-1] donde [MJ] = [10^6 J]

//;Sensible heat transfer coefficients [MJ m-2 yr-1 K-1], donde [MJ] = [10^6 J]
//;//;if (DAF eq 'o') or (DAF eq 'cr') then Hs = 5822030769.e-6					
//;//;if (DAF eq 'd') or (DAF eq 't') or (DAF eq 'g') then Hs = 13099569231.e-6

//;//;ensayando
//;//;if (DAF eq 'd') or (DAF eq 't') or (DAF eq 'g') then Hs = 13099569231.e-7		//;por ahora parece bien
*/

/*
//;*************************************************************************************************************************
//;Calcula el flujo de calor latente desde la superficie, y la evaporacion y transpiracion segun el caso
//;*************************************************************************************************************************
function FunLE, Ta0, Ztoa, Rs, Prec, qv0, RH0, criosaT, DAF
*/
real* FunLE(real Ta0,real Ztoa,real Rs,real Prec,real qv0,real RH0,real* criosaT,char* DAF)
{
  //;JFSalazar, may 2010
  //;retorna Evap en [kg m-2 yr-1] y LH en [MJ m-2 yr-1] 
  real* resultados=fltarr1(2);
 
  //;constantes fisicas
  //;[J kg-1 K-1]	specific gas constant for vapour (A00, p.36,216, CW99 p.438) - REVISADO
  real Rv = 461.0;			
  //; [J kg-1 K-1]	specific gas constant for dry air (A00, p.51,216 CW99 p.437) - REVISADO
  real R = 287.0;
  //; [adim] aprox. 0.622 cociente (masa molecular de vapor de agua)/(masa molecular promedio de aire humedo o seco) (A00, p.37) - REVISADO
  real epsilon = R/Rv;		

  //;[Pa] presion de vapor de referencia (Clausius-Clapeyron), triple punto (CW99, p.102) - REVISADO
  real etp = 611.;
  real Ttp = 273.16;		//;[K] temperatura de referencia (Clausius-Clapeyron), triple punto (CW99, p.102) - REVISADO
  real p0 = 1.013e5;		//;[Pa] presion atmosferica en la superficie (CW99 p.438) - REVISADO
  real Lv = 2.501e6;		//;[J kg-1] calor latente de vaporizacion/condensacion (gas/liquido) (A00p216) - REVISADO 
  real Ld = 0.333e6;		//;[J kg-1] calor latente de congelacion/fusion (solido/liquido) (A00p216)  - REVISADO 
  real Ls = 2.834e6;		//;[J kg-1] calor latente de sublimacion/ (gas/solido) Ls=Lv+Ld (A00p.216) - REVISADO
  real Hv = 2000.;			//;[m] water scale height (Pet00)

  //;parametros
  //;//;//;eCho = 1.8
  //;[adim] exponente en la ecuacion de evaporacion real de Choudhury 1999 - descartada esta funcion
  //;[K m-1] (Ojo: este valor tambien esta al principio del programa)
  lapse = 0.0065;			
  //;[adim] plant available water coefficient (ZhangEt99)
  real wbioi[] = {2.0, 0.5, 0.1};		
  //;[K] //;temperatura anual efectiva de acumulacion de hielo - revisar TAMBIEN ESTA EN FunQS y FunCrios
  real Tfr = 273.;			

  real Evap,LHs,Hl,Taab0,Taac0,qvsat0,es0,qv100,ETR;

  if(strcmp(DAF,"t")==0 || strcmp(DAF,"g")==0 || strcmp(DAF,"d")==0){
    //;evapotranspiracion desde tierra, ecuacion de Zhang et al, 1999.
    if((Rs>0.)&&(Prec>0.)){
      //;[kg m-2 yr-1] radiacion neta expresada como lamina de agua, evapotranspiracion potencial
      real Rslam = Rs/Lv;	
      //;//;//;Evap = Prec / (1.+ (Prec/Rslam)^eCho)^(1./eCho) 	//;Choudhury, 1999
      real wbio;
      if(strcmp(DAF,"t")==0) wbio = wbioi[0];
      if(strcmp(DAF,"g")==0) wbio = wbioi[1];
      if(strcmp(DAF,"d")==0) wbio = wbioi[2];
      real relEP = Rslam/Prec;
      Evap = Prec*(1.+ wbio*relEP)/(1.+wbio*relEP + 1./relEP);		//;ZhangEt99, Notese que ETR<Prec siempre
    }else{
      Evap = 0.;
    }

    LHs = Lv*Evap;		//;[J m-2 yr-1]

    resultados[0] = LHs*1.e-6;
    //;[MJ m-2 yr-1] , [kg m-2 yr-1]
    resultados[1] = Evap;			

  }else if(strcmp(DAF,"o")==0){
    //;evaporacion desde oceano
    Hl=3.0e6;	//;ensayando (3.20e6)
    //;//;Hl=2.25e6	//;[kg m-2 yr-1] coeficiente de transferencia de calor latente
    //;ajustado para que en la CI LHo=3062 MJ m-2 yr-1 (observado)
    //;orden de magnitud: Hl=vel*CDe*dens*86400*365, con vel 0(10), Cde 0(10-2, 10-3), dens O(0)
    real es0 = etp*exp((Lv/Rv)*(1./Ttp - 1./Ta0));	//;se asume saturacion justo en la interfase oceano-atmosfera (CW99 p.252)
    qvsat0 = epsilon*es0/p0;
    //;//;//;qvsat0 = 0.98*qvsat0			//;CW99 p.252, asumiendo salinity=35psu - decidi ignorar este detalle
    //;//;qv0 = qvsat0*RH0	//;es igual que el valor qv0 ingresado
    qv100 = qv0*exp(-100./Hv);		//;qv100<qvsat0 porque RH0=1.0 (over ocean), luego qv0=qvsat0, y qv100<qv0
    LHs = Hl*Lv*(qvsat0-qv100);		//;[J m-2 yr-1] como siempre qv100<qvsat0 entonces siempre LHs>0. (desde oceano)
    Evap = LHs/Lv;					//;[kg m-2 yr-1]
    
    resultados[0] = LHs*1.e-6;
    resultados[1] = Evap;	//;[MJ m-2 yr-1] , [kg m-2 yr-1]
  }else if(strcmp(DAF,"cr")==0){
    real aab = criosaT[0];
    real aac = criosaT[1];
    real Taab = criosaT[4];
    real Taac = criosaT[5];
    Taab0 = Taab + lapse*Ztoa/2.;
    Taac0 = Taac + lapse*Ztoa/2.;

    Hl=3.20e5;	//;[kg m-2 yr-1] Hl=315360. coeficiente de transferencia de calor latente
    //;Hl=vel*CDe*dens*86400*365, con vel=10 m/s (velocidad viento), CDE=10^-3 (CW99p.254), dens=1kg/m3

    //;qv0 y RH0 son sobre toda la criosfera. Se requieren qvab0 y qvac0. 
    //;Se asume RH0 igual en toda la criosfera
    
    //;evaporacion desde zona de ablacion (superficie libre de agua liquida)
    es0 = etp*exp((Lv/Rv)*(1./Ttp - 1./Taab0)); //;se asume saturacion al nivel de la superficie de agua como en oceano
    qvsat0 = epsilon*es0/p0;
    qv0 = qvsat0*RH0;	//;especificamente en la zona de ablacion porque es0 se calculo con Taab0

    //;RH0 viene ajustada desde el modelo de atmosfera para que sea compatible con ella (la atmosfera)
    qv100 = qv0*exp(-100./Hv);		//;qv100<qvsat0 porque RH0=1.0 (over ocean), luego qv0=qvsat0, y qv100<qv0

    //;[J m-2 yr-1] como siempre qv100<qvsat0 entonces siempre LHs>0. (desde ablacion)    
    real LHsab_evap = Hl*Lv*(qvsat0-qv100);		
    real EVAPab = LHsab_evap/Lv;					//;[kg m-2 yr-1]

    //;fusion desde la zona de ablacion (se usa calor latente de fusion, Ld)
    //;heat loss due to icemelt
    real Qab = 25.+50.*(Taab0 - Tfr);	//;[kg m-2 yr-1], NGC05 usa 550 y 1100	TAMBIEN ESTA EN FunQS
    real LHsab_melt = Ld*Qab;	//;[J m-2 yr-1]

    real LHsab = LHsab_evap + LHsab_melt;		//;[J m-2 yr-1] total desde zona de ablacion

    //;sublimacion desde zona de acumulacion (se usa calor latente de sublimacion, Ls)
    real es0 = etp*exp((Lv/Rv)*(1./Ttp - 1./Taac0)); //;se asume saturacion al nivel de la superficie de agua como en oceano
    real qvsat0 = epsilon*es0/p0;
    real qv0 = qvsat0*RH0;
    real qv100 = qv0*exp(-100./Hv);
    real LHsac = Hl*Ls*(qvsat0-qv100);	//;[J m-2 yr-1] Ls en vez de Lv
    real SUBLIac = LHsac/Ls;				//;[kg m-2 yr-1], Ls porque es sublimacion

    //;totales criosfera
    LHs = (aab*LHsab + aac*LHsac)/(aab+aac);			//;[J m-2 yr-1]

    if(VERBOSE) printf("FunLE: aab = %e, LHsab = %e, aac = %e, LHsac = %e , LHs = %e\n",aab,LHsab,aac,LHsac,LHs);

    Evap = (aab*EVAPab + aac*SUBLIac)/(aab+aac);		//;[kg m-2 yr-1]
    
    resultados[0] = LHs*1.e-6;
    resultados[1] = Evap;		//;[MJ m-2 yr-1] , [kg m-2 yr-1]
  }
  
  return resultados;
}
  
/*
//;*************************************************************************************************************************
//;Calcula la radiacion neta en la atmosfera
//;*************************************************************************************************************************
function FunRa, S, L, Ts, Ta, RH, clouds, criosaT, DAF
*/
real FunRa(real S,real L,real Ts,real Ta,real RH,real* clouds,real* criosaT,char* DAF)
{
  //;JFSalazar, may 2010
  //;Ra en [MJ m-2 yr-1] donde [MJ] = [10^6 J]
  
  //;constantes fisicas
  //; [MJ m-2 K-4 yr-1] contante de Stefan-Boltzmann. Igual que 5.67e-8 W m-2 K-4, CW99 p.437 -REVISADO
  real sigma = 1.79e-6;			

  //;De las nubes
  real ac = clouds[0];		//;area total de nubes
  real albc = clouds[1];	//;albedo promedio de las nubes
  real Tc = clouds[2];		//;temperatura promedio de las nubes
  real emc = clouds[3];		//;emisividad promedio de las nubes

  //;De la superficie
  real albcr,albs;
  if(strcmp(DAF,"cr")==0){
    albcr = ( criosaT[0]*0.40 + criosaT[1]*0.80)/( criosaT[0]+criosaT[1] );
  }
  //;//;//;criosaT(0)=aab y criosaT(1)=aac

  real albsi[] = {0.15, 0.25, 0.20, 0.10};	//;[adim] albedo de la superficie, MH, en orden: [t,g,d,o] Ojo: no incluye criosfera
  if(strcmp(DAF,"t")==0) albs = albsi[0];
  if(strcmp(DAF,"g")==0) albs = albsi[1];
  if(strcmp(DAF,"d")==0) albs = albsi[2];
  if(strcmp(DAF,"o")==0) albs = albsi[3];
  if(strcmp(DAF,"cr")==0) albs = albcr;

  //;De la atmosfera
  real aba = 0.23;			//;[adim] absorcion de ROC, aba=67/342 (TrenberthEt09) -REVISADO
  real emg = 0.25;			//;[adim] efecto invernadero de gases diferentes al vapor de agua (tambien esta en FunRs)
  real emv;
  if(strcmp(DAF,"cr")==0){
    emv=0.0;
  }else{
    //;d no hacerlo asi el efecto invernadero de vapor de agua en cr seria muy alto, irreal
    //;debido al bajo contenido de humedad atmosferica en la criosfera. Tambien en FunRs    
    emv = 1. - exp(-RH);
  }
  real Tae = (Ta + ac*Tc)/(1.+ac);
  real ema = emv*(1.-ac) + emc*ac + emg;
  ema = max(1.0, ema);

  real Ra = S*L*aba*(1+albs-albs*aba)*(1.-ac*albc) + ema*sigma*gsl_pow_int(Ts,4) - 2.*ema*sigma*gsl_pow_int(Tae,4);

  return Ra;
}

/*
//;*************************************************************************************************************************
//;Calcula calor latente en la atmosfera
//;*************************************************************************************************************************
function FunLHa, Prec
*/
real FunLHa(real Prec)
{
  //;JFSalazar, may 2010
  //;LHa en [MJ m-2 yr-1]
  
  //;constantes fisicas
  real Lv = 2.501;		//;[MJ kg-1] calor latente de vaporizacion/condensacion (gas/liquido) (A00p216) - REVISADO 

  real LHa = Prec*Lv;	//;[MJ m-2 yr-1]
  return LHa;
}

/*
//;*************************************************************************************************************************
//;Calcula la escorrentia superficial o melting segun el caso
//;*************************************************************************************************************************
function FunQs, Preci, ETRi, wsi, wmi, areai, criosaT
*/
real* FunQs(real* Preci,real* ETRi,real* wsi,real* wmi,real* areai,real* criosaT)
{
  //;devuelve el caudal superficial, sea producto de la lluvia o ice melting from the cryosphere
  //;convencion: si Qsi > 0 la region recibe agua, si Qsi < 0 la region pierde agua
  
  //;parametros del suelo
  //;//;wmi //; //;capacidad de campo, es dato de entrada a la funcion
  
  //;parametros criosfera
  //;[K] //;temperatura anual efectiva de acumulacion de hielo - revisar	
  //TAMBIEN ESTA EN FunLE y FunCrios
  real Tfr = 273.;

  //;constantes numericas
  //;cero numerico (tolerancia en la conservacion de la masa: 0.1 mm = 0.1 kg m-2)
  real ceronum = 0.1;		
    
  real at = areai[0];
  real ag = areai[1];
  real ad = areai[2];
  real ao = areai[3];
  real acr = areai[4];

  real aab = criosaT[0];
  real aac = criosaT[1];
  real Tsab = criosaT[2];
  //;from land (t,g,d)

  real* Qslandi = fltarr1(3);
  real* Qsi = fltarr1(5);
  real area,dwsdt,wsnew;


  for(int i=0;i<=2;i++){
    area = areai[i];
    if(area>0){
      if(wsi[i]<wmi[i]){
	//;[kg m-2 yr-1] = [mm yr-1] mayor que 0. siempre porque over land P-E>0 por definicion 
	dwsdt = Preci[i] - ETRi[i];	
	//;ajuste por capacidad de campo
	wsnew = wsi[i] + dwsdt;	//;dwsdt = dwsdt*1yr
	if(wsnew<=wmi[i]){Qslandi[i]=0.;}
	//;[kg m-2 yr-1] Qslandi(i)=(wsnew-wmi(i))/1yr
	else{Qslandi[i]=wsnew-wmi[i];}
      }else{
	//;en este caso wsi=wmi (nunca wsi>wmi)
	Qslandi[i] = Preci[i] - ETRi[i];			//;implica dwsdt = 0.
      }
    }else{
      Qslandi[i] = nodata;
    }
  }

  real Qt = Qslandi[0];
  real Qg = Qslandi[1];
  real Qd = Qslandi[2];
  real Qab,Qcr;

  //;from cryosphere (melting from ablation zone)
  if(acr>0.){
    if(aab>0.){
      Qab = 25.+ 50.*(Tsab - Tfr);	//;[kg m-2 yr-1], NGC05 TAMBIEN ESTA EN FunLE
      //;Tsab > Tfr siempre porque Tfr es la temperatura minima en la zona de ablacion
      //;revisar ajustar las constantes 550 y 1100 de NGC95. Las dividi por 10
      Qcr = Qab*aab/acr;		//;por unidad de area de la criosfera completa 
    }else{
      Qcr = 0.;
    }
  }else{
    Qcr = nodata;
  }

  //;to ocean (por conservacion de masa) (areao no es cero nunca, si lo fuera aparece error en otra funcion)
  //;[kg m-2 yr-1] se vuelve positivo porque entra al volumen de control
  real Qo = -(Qt*at+Qg*ag+Qd*ad+Qcr*acr)/ao;		
  //;necesarios los factores de ponderacion porque lo que se conservan
  //;son las cantidades extensivas, no las intensivas con respecto al area.
  //;los terminos donde ai=0 se anualan incluso si Qi=nodata.
  Qsi[0]=Qt;
  Qsi[1]=Qg;
  Qsi[2]=Qd;
  Qsi[3]=Qo;
  Qsi[4]=Qcr;
  //{Qt,Qg,Qd,Qo,Qcr};		//;[kg m-2 yr-1]
  
  //;conservacion de la masa
  real cero = fabs(areai[0]*Qsi[0]+
		   areai[1]*Qsi[1]+
		   areai[2]*Qsi[2]+
		   areai[3]*Qsi[3]+
		   areai[4]*Qsi[4]
		   );		//;el vector de caudales absolutos es areai*Qsi (cantidad extensiva no relativa al area)
  //;si algun area=0. entonces Qs=nodata y el producto es igual a cero
  if(cero>ceronum){
    fprintf(stderr,"ERROR en FunQs: se viola la conservacion de la masa");
    exit(1);
  }
  
  return Qsi;
}

/*
//;JFSalazar, may 2010
*/

/*
//;*************************************************************************************************************************
//;Flujo horizontal de humedad en la atmosfera
//;*************************************************************************************************************************
function FunFWa, areai, Tai, RHi
*/
real* FunFWa(real* areai,real* Tai,real* RHi)
{
  //;JFSalazar, may 2010
  //;se supone que el flujo horizontal de humedad en la atmosfera es proporcional a la diferencia de temperaturas (presiones)
  //;el flujo desde oceano siempre esta dirigido hacia afuera (ocean to land-ice). Ver Wang & Mysak 2000 p. 1153
  //;Esta funcion se basa esencialmente en FunFEa (transporte de energia)
  //;inicialmente se calcula el flujo absoluto entre regiones (en kg yr-1 MJ-1) y al final se convierten los flujos
  //;totales en cantidades intensivas por unidad de area de cada region (convergencia/divergencia total) .
  //;vectores en orden: [t,g,d,o,cr], la criosfera se trata como un todo
  //;convencion: si FWai > 0 la region gana agua, si FWai < 0 la region pierde agua

  //;parametros
  //;//;Aearth = 5.10		//;[10^14 m2] area superficial de la Tierra
  //;//;Rearth = 6370597.	//;[m] radio de la Tierra, Re = sqrt(Ae/(4*pinum)) con Ae=5.10e14 m2
  real DWm = 0.7588;		//;[10^14 kg yr-1 K-1] coeficiente de difusion meridional de humedad (calibrado, excel)
					//;DWm = DWm/5.1 = 3.87/5.1 donde 5.1x10^14 m2 = Aearth
  real DWz = 66.0784;		//;[10^14 kg yr-1 K-1] coeficiente de difusion zonal de agua, basado en Rudolf&Rubel 2005 (ver excel)
  //;DWz = DWz/5.1 = 337/5.1 donde 5.1x10^14 m2 = Aearth
  //;Nota: DWm y DWz no se entienden de la misma manera, de ahi la diferencia en el orden de magnitud.
  //;Consistentes con el orden de magnitud de los flujos en Fig. 11.1 de Rudolf&Rubel 2005.
  //;ver tambien fig 12.5 CW99 p.337
  
  //;ensayando
  //;//;DWm = 1.
  //;//;DWz = 100.
  
  //;Biotic pump (arbitrario, t>g>d)
  real Bpt = 2.;
  real Bpg = 1.;
  real Bpd = 0.5;
  
  //;constantes numericas
  real ceronum = 0.1;		//;cero numerico. Si no se fija puede haber problemas en la verificacion de la conservacion de energia
  real nodata = -9999.;

  real at = areai[0];
  real ag = areai[1];
  real ad = areai[2];
  real ao = areai[3];
  real acr = areai[4];

  real Tt = Tai[0];
  real Tg = Tai[1];
  real Td = Tai[2];
  real To = Tai[3];
  real Tcr = Tai[4];

  real RHt = RHi[0];
  real RHg = RHi[1];
  real RHd = RHi[2];
  real RHo = RHi[3];
  real RHcr = RHi[4];

  real aref,Wto,RH;
  real Wtn,Wtg,Wgo,Wgt,Wgn,Wdo,Wdt,Wdcr,Wod,Wog,Wocr,Wdg,Wot,Wcro,Wcrn,Wcrt,Wcrg;
  real FWt,FWg,FWo,FWcr,FWd;
  real control;
  int i;

  //;medidas de las areas a traves de las cuales se dan los flujos meridionales hacia/desde la criosfera
  real tetacr = acos(1.-acr);
  real Lcr = sin(tetacr);		//;si acr=0. entonces Lcr=0. (no hay intercambios con la criosfera porque esta no existe)
  real Locr,Ldcr,Ltcr,Lgcr;
  if(ad>0.){
    Locr = Lcr*ao/(ao+ad);
    Ldcr = Lcr*ad/(ao+ad);
  }else{
    Locr = Lcr*ao/(ao+at+ag);
    Ltcr = Lcr*at/(ao+at+ag);
    Lgcr = Lcr*ag/(ao+at+ag);
  }
  if(VERBOSE>1) printf("Locr,Ldcr=%e,%e\n",Locr,Ldcr);
  
  //;*************************************************************************************************************************
  //;trees
  if(at>0.){
    //;ocean
    aref = min(ao,at);
    Wto = DWz*fabs(To-Tt)*aref*RHo*Bpt;	//;valor absoluto porque siempre va de o hacia t. Por eso tambien RHo
    //;grass
    aref = min(ag,at);
    if(Tt>=Tg){RH=RHt;}else{RH=RHg;}
    //;notese que si aref=0. el flujo es 0. aref es una medida del area a traves de la cual se da la transferencia
    Wtg = DWz*(Tg-Tt)*aref*RH;		
    //;desert or cryosphere	
    if(ad>0.){
      aref = min(ad,at);
      if(Tt>=Td){RH=RHt;}else{RH=RHd;}
      Wtn = DWz*(Td-Tt)*aref*RH;
    }else{
      //;el intercambio es con la criosfera y es meridional, no zonal
      if(Tt>=Tcr){RH=RHt;}else{RH=RHcr;}
      Wtn = DWm*(Tcr-Tt)*Ltcr*RH;
    }
    //;flujo. Para obtener la convergencia en t basta dividir por at (esto se hace al final)
    FWt = Wto + Wtg + Wtn;	
  }else{
    FWt = 0.;
  }

  //;*************************************************************************************************************************
  //;grass
  if(ag>0.){
    //;ocean
    aref = min(ao,ag);
    Wgo = DWz*fabs(To-Tg)*aref*RHo*Bpg;
    //;trees
    aref = min(at,ag);
    if(Tg>=Tt){RH=RHg;}else{RH=RHt;}
    Wgt = DWz*(Tt-Tg)*aref*RH;
    //;desert or cryosphere	
    if(ad>0.){
      aref = min(ad,ag);
      if(Tg>=Td){RH=RHg;}else{RH=RHd;}
      Wgn = DWz*(Td-Tg)*aref*RH;
    }else{
      //;el intercambio es con la criosfera y es meridional, no zonal
      if(Tg>=Tcr){RH=RHg;}else{RH=RHcr;}
      Wgn = DWm*(Tcr-Tg)*Lgcr*RH;
    }
    //;total
    FWg = Wgo + Wgt + Wgn;	//;flujo. Para obtener la convergencia en g basta dividir por ag (esto se hace al final)
  }else{
    FWg = 0.;
  }

  //;*************************************************************************************************************************
  //;desert
  if(ad>0.){
    //;ocean (flujo zonal)
    aref = min(ao,ad);
    Wdo = DWz*fabs(To-Td)*aref*RHo*Bpd;
    //;trees  (flujo zonal)
    aref = min(at,ad);
    if(Td>=Tt){RH=RHd;}else{RH=RHt;}
    Wdt = DWz*(Tt-Td)*aref*RH;
    //;grass (flujo zonal)
    aref = min(ag,ad);
    if(Td>=Tg){RH=RHd;}else{RH=RHg;}
    Wdg = DWz*(Tg-Td)*aref*RH;
    //;cryosphere (flujo meridional)
    if(Td>=Tcr){RH=RHd;}else{RH=RHcr;}
    Wdcr = DWm*(Tcr-Td)*Ldcr*RH;			
    //;total
    FWd = Wdo + Wdt + Wdg + Wdcr;
  }else{
    FWd = 0.;
  }

  //;*************************************************************************************************************************
  //;ocean (ao > 0. siempre)
  //;trees (zonal)
  aref = min(at,ao);
  Wot = -DWz*fabs(Tt-To)*aref*RHo*Bpt;
  //;grass (zonal)
  aref = min(ag,ao);
  Wog = -DWz*fabs(Tg-To)*aref*RHo*Bpg;
  //;desert (zonal)
  aref = min(ad,ao);
  Wod = -DWz*fabs(Td-To)*aref*RHo*Bpd;
  //;cryosphere (meridional)
  Wocr = -DWm*fabs(Tcr-To)*Locr*RHo;
  //;total
  FWo = Wot + Wog + Wod + Wocr;

  //;*************************************************************************************************************************
  //;cryosphere (todos son flujos meridionales)
  if(acr>0.){
    //;ocean
    Wcro = DWm*fabs(To-Tcr)*Locr*RHo;
    if(ad>0.){
      if(Tcr>=Td){RH=RHcr;}else{RH=RHd;}
      Wcrn = DWm*(Td-Tcr)*Ldcr*RH;
    }else{
      //;trees
      if(Tcr>=Tt){RH=RHcr;}else{RH=RHt;}
      Wcrt = DWm*(Tt-Tcr)*Ltcr*RH;
      //;grass
      if(Tcr>=Tg){RH=RHcr;}else{RH=RHg;}
      Wcrg = DWm*(Tg-Tcr)*Lgcr*RH;
      //;total from/to land
      Wcrn = Wcrt+Wcrg;
    }
    FWcr = Wcro + Wcrn;
  }else{
    FWcr = 0.;
  }

  if(VERBOSE>1) printf("FWt,FWg,FWd,FWo,FWcr = %e,%e,%e,%e,%e\n",FWt,FWg,FWd,FWo,FWcr);

  //;*************************************************************************************************************************
  //;resultados
  real* FWai=fltarr1(5);
  FWai[0]=FWt;
  FWai[1]=FWg;
  FWai[2]=FWd;
  FWai[3]=FWo;
  FWai[4]=FWcr;

  //;//;//;Idea: remover esta modificacion de FWai porque hace crecer anormalmente a FWai. Pensar en involucrar umed
  //;//;WvL = total(Wvi*areai)
  //;//;FWai = FWai*WvL/WvL1

  //;conservacion de la energia (aplica sobre los flujos, no sobre la convergencia)
  control = fabs(FWai[0]+FWai[1]+FWai[2]+FWai[3]+FWai[4]);
  if(control>ceronum){
    fprintf(stderr,"ERROR en FunFWa: se viola la conservacion de la masa\n");
    exit(1);
  }

  //;Convirtiendo los flujos a convergencia sobre cada region.
  //;Convergencia/divergencia total sobre cada region (por unidad de area de cada region)
  //ceros = where(areai eq 0.)
  //nocero = where(areai ne 0.)
  for(i=0;i<5;i++){
    if(areai[i]==0){
      FWai[i]=nodata;
    }else{
      FWai[i]=FWai[i]/areai[i];
    }
  }
  return FWai;
}

/*
//;*************************************************************************************************************************
//;Calcula la tasa de mortalidad de las especies biologicas (trees y grasses).
//;*************************************************************************************************************************
function FunTmor, Ws, Wm
*/
real FunTmor(real Ws, real Wm)
{
  //;//;//;ensayando tasa constante
  //;//;//;tmor=0.3
  //;//;//;return, tmor
  
  //;parametros
  real Wsopt = Wm;			//;humedad del suelo optima
  //;"factor de forma", constante arbitraria, tiene que ver con la forma de la funcion tmor=f(ms)
  real ff = 1.19;			
  real tmormin = 0.2;
  real tmor = min(1.0,tmormin+gsl_pow_int(((Ws - Wsopt)/Wsopt)/ff,2));
  return tmor;
}

/*
//;*************************************************************************************************************************
//;Calcula la tasa de crecimiento de las especies biologicas (trees y grasses).
//;*************************************************************************************************************************
function FunBeta,Ts,Topt,Topt0
*/
real FunBeta(real Ts,real Topt,real Topt0)
{
  //;JFSalazar, may 2010
  //;Involucra implicitamente capacidad de adaptacion en escala geologica a traves del valor Topt
  
  
  //;parametros
  //;//;Topt0 = 295.5		//;[K] temperatura optima de referencia, igual al valor inicial que a su vez es igual a Tst inicial
  //;=(298+293)/2
  
  real KBmax = 60.;		//;60 ad hoc (estaba 67)constante que define la amplitud del intervalo de valores de temperatura donde cabe la adaptacion
				//;este intervalo es [Topti0-KBmax, Topti0+KBmax].
  real KB = 17.5;		//;constante que define la amplitud del intervalo de valores de temperatura donde Beta varia
				//;parabolicamente alrededor de Topt(t). Se deja igual que en el Daisyworld.
  
  real Bmax = 1.-gsl_pow_int(((Topt0-Topt)/KBmax),2.);	//;no hace falta imponer aqui la condicion Bmaxt>=0 porque esta mas adelante en Beta
  //;notese que sin adaptacion (tada=0.) Bmaxt=Bmaxg=1.0 siempre porque Toptt y Toptg 
  //;no cambian sino que se mantienen iguales a sus valores iniciales Toptt0 y Toptg0
  
  real Beta = max(0.,Bmax*(1.-gsl_pow_int(((Ts-Topt)/KB),2.)));
  
  return Beta;
}

//*
/*
//;*************************************************************************************************************************
//;Encuentra particion de la criosfera en ablation and accumulation zones, resultado = [Aab, Aac, Tab, Tac]
//;*************************************************************************************************************************
*/
int FunCrios(real Tscr,real acr,real tetacr,real c1,real c2,real amin,real crios[2])
{
  real minTscr,maxTscr,aac,aab,tetaac;//;

  if(acr>0){
    minTscr = c1*Tscr;////;minima temperatura en la criosfera
    maxTscr = (c1+c2*sqrt(tetacr))*Tscr;////;maxima temperatura en la criosfera

    if((Tfr>minTscr) && (Tfr<maxTscr)){
      ////;la criosfera tiene ambas regiones: ac y ab
      tetaac = ((1./c2)*(Tfr/Tscr - c1))*((1./c2)*(Tfr/Tscr - c1));////;colatitud de la linea de equilibrio que limita la zona de acumulacion
      aac = 1.-cos(tetaac);
    }else{
      if(Tfr>=maxTscr){
	////;toda la criosfera es de acumulacion
	aac = acr;
      }else{
	////;toda la criosfera es de ablacion
	aac = 0.;
      }
    }
    
    ////;amin es un cero numerico para evitar numeros demasiado pequenos
    if(aac<amin){aac=0.;}

    aab = acr-aac;
    if(aab<amin){aab=0.;}

  }else{
    aab = 0.;
    aac = 0.;
  }
  crios[0]=aac;
  crios[1]=aab;
  return 0;
}
//*/

/*
//;*************************************************************************************************************************
//;Calcula los ajustes DAFM
//;*************************************************************************************************************************
function FunDafm, areai, dareaidt, Hi
*/
real* FunDafm(real* areai,real* dareaidt,real* Hi)
{
  //;JFSalazar, may 2010
  //;sirve para hallar Ajustes de Tsi, Tai y Wai. NO para Wsi
  //;esta supuesto que ad > 0. siempre como en NGC05 (pendiente verificar efecto de esta suposicion).
  //;la criosfera se trata separadamente en ab y ac. Por eso los vectores de entrada contienen: [t,g,d,o,ab,ac]
  //;Sin embargo en los resultados se promedia sobre toda la criosfera, por eso la salida  es: [t,g,d,o,cr]
  real Ht = Hi[0];
  real Hg = Hi[1];
  real Hd = Hi[2];
  real Ho = Hi[3];
  real Hab = Hi[4];
  real Hac = Hi[5];

  real datdt = dareaidt[0];
  real dagdt = dareaidt[1];
  real daodt = dareaidt[3];
  real daabdt = dareaidt[4];
  real daacdt = dareaidt[5];

  real at = areai[0];
  real ag = areai[1];
  real ad = areai[2];
  real ao = areai[3];
  real aab = areai[4];
  real aac = areai[5];

  real KHd = 0.; //;para ir sumando
  real KHt,KHg,KHo,KHcr,KHac,dacrdt,daabdt1,daabdt2;
  real KHab,KHab1,KHab2,acr;

  //;ocean-land
  if(at>0){
    if(datdt<=0.){
      KHt = 0.;
      //;igual que KHd=KHd-(-datdt)*(Hd-Ht)/ad=KHd-(daddt)*(Hd-Ht)/ad con daddt>0
      KHd = KHd + datdt*(Hd-Ht)/ad;
    }else{
      KHt = -datdt*(Ht-Hd)/at;
    }
  }else{
    KHt = 0.;
  }
  
  if(ag>0.){
    if(dagdt<=0.){
      KHg = 0.;
      KHd = KHd + dagdt*(Hd-Hg)/ad;
    }else{
      KHg = -dagdt*(Hg-Hd)/ag;
    }
  }else{
    KHg = 0.;
  }

  if(daodt<=0.){
    KHo = 0.;
    KHd = KHd + daodt*(Hd-Ho)/ad;
  }else{
    KHo = -daodt*(Ho-Hd)/ao;	//;ao no puede ser nula
  }

  //;cryosphere. Hay que especificar hacia donde se expanden/contraen las zonas de ablacion/acumulacion
  acr = aab+aac;
  if(acr>0.){
    if(aac>0.){
      if(daacdt<=0.){
	KHac = 0.;
	if(aab==0.){KHd = KHd + daacdt*(Hd-Hac)/ad;}	//;como no hay ab la contraccion de acc implica expansion de ad
      }else{
	//;aac se expande hacia aab (si existe), o hacia ad (si no existe aab)
	if(aab>0.){KHac = -daacdt*(Hac-Hab)/aac;}else{KHac = -daacdt*(Hac-Hd)/aac;}
      }
    }else{
      KHac = 0.;
    }

    dacrdt = daacdt + daabdt;

    if(aab>0.){
      if(daabdt<=0.){
	KHab = 0.;
      }else{
	if(daacdt>0.){
	  KHab = -daabdt*(Hab-Hd)/aab;		//;aab se expande hacia ad, crece toda la criosfera
	}else{
	  //;aab se expande y aac se contrae
	  if(dacrdt<=0.){
	    //;la zona de ablacion se expande unicamente hacia la de acumulacion porque la criosfera no crece
	    KHab = -daabdt*(Hab-Hac)/aab;
	  }else{
	    //;la zona de ablacion se expande hacia acumulacion y hacia bare land (d) al mismo tiempo.
	    daabdt1 = -daacdt;			//;todo lo que se reduce acumulacion lo ocupa ablacion
	    daabdt2 = daabdt - daabdt1;	//;el resto se expande hacia bare land (d)
	    KHab1 = -daabdt1*(Hab-Hac)/aab;	//;hacia acumulacion
	    KHab2 = -daabdt2*(Hab-Hd)/aab;	//;hacia bare land (d)
	    KHab = KHab1 + KHab2;
	  }
	}
      }
      //;ajuste KHd	
      if(dacrdt<0.){
	//;la criosfera se contrae y como aab>0 entonces lo hace necesariamente sobre la zona de ablacion
	KHd = KHd + dacrdt*(Hd-Hab)/ad;
      }
    }else{
      KHab = 0.;
    }
    KHcr = (aab*KHab + aac*KHac)/acr; 	//;promediando sobre toda la criosfera
  }else{
    KHcr = 0.;
  }
  
  real* AjusDafm=fltarr1(5);
  AjusDafm[0]=KHt;
  AjusDafm[1]=KHg;
  AjusDafm[2]=KHd;
  AjusDafm[3]=KHo;
  AjusDafm[4]=KHcr;

  return AjusDafm;
}

/*
//;*************************************************************************************************************************
//;Calcula el ajuste Dafm para Ws
//;*************************************************************************************************************************
function FunDafmWs, areai, dareaidt, Hi
*/
real* FunDafmWs(real* areai,real* dareaidt,real* Hi)
{
  //;JFSalazar, may 2010
  //;este ajuste solo tiene sentido en las DAFs land (t,g,d) donde opera el BAS
  
  real Ht = Hi[0];
  real Hg = Hi[1];
  real Hd = Hi[2];

  real datdt = dareaidt[0];
  real dagdt = dareaidt[1];
  
  real at = areai[0];
  real ag = areai[1];
  real ad = areai[2];

  real KHd = 0.; //;para ir sumando
  real KHt,KHg;

  //;ocean-land
  if(at>0.){
    if(datdt<=0.){
      KHt = 0.;
      KHd = KHd + datdt*(Hd-Ht)/ad;		//;igual que KHd=KHd-(-datdt)*(Hd-Ht)/ad=KHd-(daddt)*(Hd-Ht)/ad con daddt>0
    }else{
      KHt = -datdt*(Ht-Hd)/at;
    }
  }else{
    KHt = 0.;
  }

  if(ag>0.){
    if(dagdt<=0.){
      KHg = 0.;
      KHd = KHd + dagdt*(Hd-Hg)/ad;
    }else{
      KHg = -dagdt*(Hg-Hd)/ag;
    }
  }else{
    KHg = 0.;
  }
  real* AjusDafm=fltarr1(3);
  AjusDafm[0]=KHt;
  AjusDafm[0]=KHg;
  AjusDafm[0]=KHd;
  
  return AjusDafm;
}

/*
//;*************************************************************************************************************************
//;Calcula el forzamiento debido a la regulacion biotica de la temperatura local
//;*************************************************************************************************************************
function FunBioticReg, Ts, Topt, ROCs, Rs, LHs, Prec, ETR, bs
*/
real* FunBioticReg(real Ts,real Topt,real ROCs,real Rs,real LHs,real Prec,real ETR,real bs)
{
  //;JFSalazar, jun 2010
  //;basado en la idea de "biotic regulation of the environment" by GGM 2000
  
  //;constantes fisicas
  //;[MJ m-2 K-1] calor especifico de la superficie land (equivalente aprox a 2.5 m de agua)  
  real cps = 10.;
  //;[MJ kg-1] calor latente de vaporizacion/condensacion (gas/liquido) - REVISADO
  real Lv = 2.501;

  //;parametros
  //;biotic regulation capacity. 1.0 by default, (bc = 0.0 control experiment)
  real bc = 5.0;
  real sigE,sigW,limE,limW,limB,absDeltaETR,deltaETR,BRetr,BRalb;
  real* resultados=fltarr1(3);

  if(bc==0.){
    zeroVector(resultados,3);
    return resultados;
  }

  if(bs<0.){
    if(Ts>Topt){
      sigE = -1.;
      sigW = 1.;
    }else{
      sigE = 1.;
      sigW = -1.;
    }
    if(sigW==1.){
      //;Rs en [MJ m-2 yr-1]
      limE = max(0.,Rs/Lv);
      limW = Prec-ETR;
      limB = bc*ETR/10.;
      absDeltaETR = min(limE,min(limW,limB));	//;aumento
    }else{
      absDeltaETR = bc*ETR/10.;	//;reduccion
    }
    deltaETR = sigW*absDeltaETR;
    BRetr = sigE*Lv*absDeltaETR;
    BRalb = sigE*bc*ROCs/10.;
  }else{
    deltaETR = 0.;
    BRetr = 0.;
    BRalb = 0.;
  }
  resultados[0]=BRalb;
  resultados[1]=BRetr;
  resultados[2]=deltaETR;
  return resultados;
}

real minVector(real* vector,int n)
{
  real min=1E100;
  for(int i=0;i<n;i++)
    min=min<vector[i]?min:vector[i];
  return min;
}
