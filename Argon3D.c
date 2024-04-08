#include<stdio.h>          /*  ESTE ES UN PROGRAMA DE DINAMICA MOLECULAR      */
#include<stdlib.h>         /*  PARA UN GRUPO DE PARTICULAS EN TRES DIMEN       */ 
#include<conio.h>          /*  SIONES QUE INTERACTUAN CON UN POTENCIAL DE     */
#include<math.h>           /*  DE LENNAR-JONES                                */ 
#include<malloc.h>                                                                

void INIT_CELDA( ) ;         /*  OTORGA POSICIONES Y VELOCIDADES A LOS ATOMOS*/
void CONDC_FRONT( );         /*  RELIZA EL CONTROL DE FRONTERAS*/
void MOVE_MENT( );           /*  MUEVE LOS ATOMOS              */
void VEL_cmNULL( ) ;         /*  ADJUDICA VELOCIDADES CUYA SUMA VECTORIAL ES NULA     */
void VELOCITIES( ) ;         /*  ADJUDICA VELOCIDADES DE ACUERDO A LA TEMPERATURA DE SIMULACION            */
void D_MAXWELL( );
void CALCULA_FUERZAS( int a );     /* CALCULA LAS PROYECCIONES DE LAS FUERZAS DEPENDIENDO DE a SOBRE CADA ATOMO */
void ESCALA(long int a, double e); /* ESTA FUNCION ESCALA LAS VELOCIDADES*/
double ENER_POTC(  );              /*  CALCULA LA EERGIA POTENCIAL DE LOS ATOMOS            */
double ENER_CINT(  );          /*  CALCULA LA ENERGIA CINETICA DE LOS ATOMOS            */
void ACTUALIZA_TablaVerlet();

double aleat( );             /*  GENERA NUMEROS PSEUDOALEATORIOS  */          
                                                                                         /*   DNM2LJ.c  ACHIVO PRINCIPAL        */
                                                                                         /*   DNM2LJ1.c  INIT_CELDA( )          */
const double sigma = 1.0;          /*  DIAMETRO DEL ATOMO EN UNIDADES REDUCIDAS    */    /*   DNM2LJ2.c  ENER_POTC()            */ 
const double Rcut = 2.5 ;          /*  DISTANCIA DE TRUNCAMIENTO DE LA INTERACCION */    /*   DNM2LJ3.c  ENER_CINT()            */
const int N = 4        ;         /*  NUMERO DE PARTICULAS */                          /*   DNM2LJ4.c  D_MAXWELL()            */  
                                
const double h = 0.032  ;          /*  PASO DE TIEMPO O PASO DE LA SIMULACION*/          /*   DNM2LJ7.c  CONDC_FRONT()          */
const double T = 2.53   ;          /*  TEMPERATURA EN UNIDADES REDUCIDAS     */          /*   DNM2LJ8.c  VELOCITIES( )          */ 
static double Rm = 3.3 ;
double V                ;          /*  VELOCIDAD EN UNIDADES REDUCIDAS*/                 /*   DNM2LJ9.c  CALCULA_FUERZAS(int )     */
double A                ;          /*  POSICION DONDE SE COLOCAN LOS ATOMOS  */          /*   DNM2LJ6.c  MOVE_MENT( long int a )*/ 
double L                ;         /*  LONGITUD DE LA CELDA EN UNIDADES DE SIGMA*/       /*   DNM2LJ5.c  VEL_cmNULL()           */                                               
static int M[3];
long int PASOS;

double X[ 1 ];                 /*  SEMILLA DEL GENERADOR DE NUMEROS PSEUDOALEATORIOS    */
long double *RVFiFf; /****************  ESTE ARREGLO CONTIENE LAS POSICIONES VELOCIDADES ****************/
long int *TablaVerlet;                /*  Y FUERZAS QUE ACTUAN SOBRE LAS PARTICULAS            */
                               /*  EN EL MOMENTO INICIAL Y EN EL MOMENTO FINAL          */ 
         

int main()
{
   long int i;
   int j;
   double Ec,Ep;
   double beta; 
  
   FILE *F;
   FILE *G;
   FILE *H;
   FILE *I;
 
   RVFiFf = malloc( 12*pow(N,4)*sizeof(long double) );                            /* RESERVA MEMORIA PARA LOS DATOS DE ENTRADA*/
   if(RVFiFf!=0)printf("MEMORIA RESERVADA PARA EL ARREGLO: OK\n");   /* CHEQUEA LA RESERVACION DE MEMORIA        */
   else printf("MEMORIA NO RESERVADA PARA EL ARREGLO: FOUL\n");
   
   TablaVerlet = malloc(151*pow(N,4)*sizeof( int));
   if(TablaVerlet!=0)printf("MEMORIA RESERVADA PARA TablaVerlet: OK\n");
   else printf("MEMORIA NO RESERVADA PARA TablaVerlet: FOUL\n");
   
                                
   printf("NUMERO DE PASOS DE LA SIMULACION");
   scanf("%ld",&PASOS);
 
   X[0] = 123457.0    ;
   V = 0.3668;                                              /* VELOCIDAD MEDIA*/
   L = 7.38;
   A = L/N      ;
   Ec = 0;
   Ep = 0;
   M[0] = 0 ; 
   M[1] = -1;
   M[2] =  1;
   
    
   INIT_CELDA();
   
   F = fopen("Argon0032_Ec_450_2000.txt","w");
   G = fopen("Argon0032_Ep_450_2000.txt","w");
   H = fopen("Argon0032_E_450_2000.txt","w");
   I = fopen("Argon0032_T_450_2000.txt","w");
  
   for(i=1;i<PASOS+1;i++)
   {
    Ec = ENER_CINT(  )  ;
    Ep = ENER_POTC(  )  ;
    MOVE_MENT(  );
       
    fprintf(F,"%ld\t%1.2lf\n",i,24*Ec);
    fprintf(G,"%ld\t%1.2lf\n",i,4*Ep);
    fprintf(H,"%ld\t%1.2lf\n",i,4*Ep+24*Ec);
    fprintf(I,"%ld\t%1.2lf\n",i,24*Ec/(1.5*N)); 
    
    if(fmod(i,25)==0)ACTUALIZA_TablaVerlet();
    if(i <= 450)ESCALA(i, Ec);    
   
    printf( "%ld\n",i )  ;
   }
   
   fclose(F);
   fclose(G);
   fclose(H);
   fclose(I);

   free(RVFiFf);
   

   system("PAUSE");
 
   return 0;
}

/******************************************************************************
 ******************************************************************************
 ******************************************************************************/

double aleat()
{
   unsigned k;
   double m;
 
   m = pow(2,31)-1;
   k=16807U;

   X[0]=fmod((k*X[0]),m);

   return(X[0]/m);

}



/******************************************************************************/ 


void CALCULA_FUERZAS( int a )               /*SI a=0 FUERZA DE LA n POSICION, SI a=2 FUERZA DE LA n+1 POSICION*/
{
   int i,j,k,l,m;
   double RijAlCuadrado;
   
   for(i = 0; i < (int)pow(N,4); i++)RVFiFf[12*i+6+a] = RVFiFf[12*i+7+a] = RVFiFf[12*i+8+a] = 0;
 
   for(i = 0; i < (int)pow(N,4); i++)
   {
   for(j = 0; j < TablaVerlet[151*i + 150]; j++) /* TablaVerlet[151*i + 150] contiene la cantidad de particulas  */
   {                                             /* Que interaccionan con la particula i  */
   
    for(k = 0; k < 3; k++)
    {
    for(m = 0; m < 3; m++)
    {
    for(l = 0; l < 3; l++)                       /* TablaVerlet[151*i + j] contiene el numero de cada particula */ 
    {                                            /* Que interacciona con la particula i  */
     RijAlCuadrado =
     pow( RVFiFf[12*i     ] - (RVFiFf[12*TablaVerlet[151*i + j]     ] + M[ k ]*L), 2) +
     pow( RVFiFf[12*i + 1 ] - (RVFiFf[12*TablaVerlet[151*i + j] + 1 ] + M[ m ]*L), 2) +
     pow( RVFiFf[12*i + 2 ] - (RVFiFf[12*TablaVerlet[151*i + j] + 2 ] + M[ l ]*L), 2) ;

     if(RijAlCuadrado <= Rcut*Rcut)/*factor 0.3*pow(10,12)*/
     {
      RVFiFf[12*i+6+a] += ( RVFiFf[12*i    ] - (RVFiFf[12*TablaVerlet[151*i + j]    ] + M[ k ]*L) )*( pow(sigma/RijAlCuadrado, 7) - 
                                                                                                  0.5*pow(sigma/RijAlCuadrado, 4));
      RVFiFf[12*i+7+a] += ( RVFiFf[12*i + 1] - (RVFiFf[12*TablaVerlet[151*i + j] + 1] + M[ m ]*L) )*( pow(sigma/RijAlCuadrado, 7) - 
                                                                                                  0.5*pow(sigma/RijAlCuadrado, 4)); 
      RVFiFf[12*i+8+a] += ( RVFiFf[12*i + 2] - (RVFiFf[12*TablaVerlet[151*i + j] + 2] + M[ l ]*L) )*( pow(sigma/RijAlCuadrado, 7) - 
                                                                                                  0.5*pow(sigma/RijAlCuadrado, 4));                                                                           
     }
    }
    }
    }
       
   }      
   }
}



/******************************************************************************/



 void INIT_CELDA( )
 {
   int i,j,k,l;
   int Posicion_n;
   
   Posicion_n = 0;
   l = 0;
   
     for(i = 0; i < 12*(int)pow(N,4); i++)RVFiFf[i] = 0;                                       /* INICIALIZA A CERO EL ARREGLO DE DATOS   */
    
/**********************DISTRIBUCION DE MOLECULAS EN LA CELDA*******************/  /*   DNM2LJ.c  FUNCION PRINCIPAL       */
                                                                                  /*   DNM2LJ1.c  INIT_CELDA( )          */
     for(i = 0; i < N; i ++)                                     
     {                                                                           
     for(j = 0; j < N; j ++)                                     
     {                                                                         
     for(k = 0; k < N; k ++,l ++)                /*PUNTOS NEGROS 1*/
     {                                                                            
      RVFiFf[ 12*l     ] = 0.75*A + k*A ;                                            
      RVFiFf[ 12*l + 1 ] = 0.25*A + j*A ;                                             
      RVFiFf[ 12*l + 2 ] = 0.25*A + i*A ;                                                                             
     }
     }
     }
     
/******************************************************************************/     
     
     for(i = 0; i < N; i ++)                                     
     {                                                                           
     for(j = 0; j < N; j ++)                                     
     {                                          /*PUNTOS NEGROS 2*/                               
     for(k = 0; k < N; k ++,l ++)
     {                                                                            
      RVFiFf[ 12*l     ] = 0.25*A + k*A ;                                            
      RVFiFf[ 12*l + 1 ] = 0.25*A + j*A ;                                             
      RVFiFf[ 12*l + 2 ] = 0.75*A + i*A ;                                                                            
     }
     }
     }
     
/******************************************************************************/ 
  
     for(i = 0; i < N; i ++)                                     
     {                                                                           
     for(j = 0; j < N; j ++)                                     
     {                                        /*CIRCULOS ROJOS 1*/                                 
     for(k = 0; k < N; k ++,l ++)
     {                                                                            
      RVFiFf[ 12*l     ] = 0.25*A + i*A ;                                            
      RVFiFf[ 12*l + 1 ] = 0.75*A + j*A ;                                             
      RVFiFf[ 12*l + 2 ] = 0.25*A + k*A ;                                                                            
     }
     }
     }
     
/******************************************************************************/      
     
     for(i = 0; i < N; i ++)                                     
     {                                                                           
     for(j = 0; j < N; j ++)                                     
     {                                        /*CIRCULOS ROJOS 2*/                                 
     for(k = 0; k < N; k ++,l ++)
     {                                                                            
      RVFiFf[ 12*l     ] = 0.75*A + i*A ;                                            
      RVFiFf[ 12*l + 1 ] = 0.75*A + j*A ;                                             
      RVFiFf[ 12*l + 2 ] = 0.75*A + k*A ;                                                                            
     }
     }
     }
     
  /*  VELOCITIES( );*/
    
    D_MAXWELL();

    VEL_cmNULL( );    /***ANULAMIENTO DE LA VELOCIDAD DEL CENTRO DE MASA****/    
     
    ACTUALIZA_TablaVerlet();
   
    CALCULA_FUERZAS( Posicion_n  );

    
    
}



/******************************************************************************/




double ENER_POTC( )
{
   int i,j,m,l,k;
   double RijAlCuadrado,E;
                                                                                   /*   DNM2LJ.c  ACHIVO PRINCIPAL        */
                                                                                   /*   DNM2LJ1.c  INIT_CELDA( )          */
   E      =  0;                                                                      /*   DNM2LJ2.c  ENER_POTC()            */
                                                                                   /*   DNM2LJ3.c  ENER_CINT()            */
                                                                                   /*   DNM2LJ8.c  D_MAXWELL()            */
                                                                                   /*   DNM2LJ5.c  VEL_cmNULL()           */
                                                                                   /*   DNM2LJ6.c  MOVE_MENT( long int a )*/  
                                                                                   /*   DNM2LJ7.c  CONDC_FRONT()          */
                                                                                   /*   DNM2LJ8.c  VELOCITIES( )          */
   for(i =     0; i < (int)pow(N,4) - 1; i ++)                                                 /*   DNM2LJ9.c  CALCULA_FUERZAS(int )     */
   {
   for(j = i + 1; j < (int)pow(N,4)    ; j ++)
   {
    for(k = 0; k < 3; k ++)
    {
    for(m = 0; m < 3; m ++)
    {
    for(l = 0; l < 3; l++)
    {      
      RijAlCuadrado   =
      pow( RVFiFf[ 12*i     ] - (RVFiFf[ 12*j     ] + M[ k ]*L), 2) +
      pow( RVFiFf[ 12*i + 1 ] - (RVFiFf[ 12*j + 1 ] + M[ m ]*L), 2) +
      pow( RVFiFf[ 12*i + 2 ] - (RVFiFf[ 12*j + 2 ] + M[ l ]*L), 2)  ;

      if(RijAlCuadrado <= Rcut*Rcut)E += ( pow(sigma/RijAlCuadrado,6) - pow(sigma/RijAlCuadrado,3) );
 
    }
    }
    }   
   
   }
   }
  
 return E; 
}



/******************************************************************************/



double ENER_CINT( )
{
   int i; 

                                                                                  /*   DNM2LJ.c  ACHIVO PRINCIPAL        */
   double E  = 0                ;                                                 /*   DNM2LJ1.c  INIT_CELDA( )          */
                                                                                  /*   DNM2LJ2.c  ENER_POTC()            */
                                                                                  /*   DNM2LJ3.c  ENER_CINT()            */ 
                                                                                  /*   DNM2LJ4.c  D_MAXWELL()            */
   for(i=0; i<(int)pow(N,4); i++)                                                             /*   DNM2LJ5.c  VEL_cmNULL()           */ 
   {                                                                              /*   DNM2LJ6.c  MOVE_MENT( long int a )*/
    E += RVFiFf[12*i + 3]*RVFiFf[12*i + 3]+                                       /*   DNM2LJ7.c  CONDC_FRONT()          */
         RVFiFf[12*i + 4]*RVFiFf[12*i + 4]+                                       /*   DNM2LJ8.c  VELOCITIES( )          */ 
         RVFiFf[12*i + 5]*RVFiFf[12*i + 5];
   }                                                                              /*   DNM2LJ9.c  CALCULA_FUERZAS(int )     */

   return (E) ;
}



/******************************************************************************/



void MOVE_MENT(  )
{
    int    i,m,l,k,Posicion_nmas1;                                                    /*   DNM2LJ.c  ACHIVO PRINCIPAL        */
                                                                                  /*   DNM2LJ1.c  INIT_CELDA( )          */ 
                                                                                  /*   DNM2LJ2.c  ENER_POTC()            */ 
    Posicion_nmas1 = 3; /*CALCULO DE FUERZAS EN LA n+1 POSICION*/                   /*   DNM2LJ3.c  ENER_CINT()            */
                                                                                  /*   DNM2LJ4c  D_MAXWELL()            */ 
                                                                                  /*   DNM2LJ5.c  VEL_cmNULL()           */ 
                                                                                  /*   DNM2LJ6.c  MOVE_MENT( long int a )*/
                                                                                  /*   DNM2LJ7.c  CONDC_FRONT()          */
                                                                                  /*   DNM2LJ8.c  VELOCITIES( )          */
                                                                                  /*   DNM2LJ9.c  CALCULA_FUERZAS(int )  */
 
/*****************************************CALCULO DE LA (n+1) POSICION***************************************************/
 
    for(i = 0; i < (int)pow(N,4); i ++)
    {
     RVFiFf[12*i   ] = RVFiFf[12*i   ] + h*RVFiFf[12*i + 3] + 0.5*h*h*RVFiFf[12*i + 6];           
     RVFiFf[12*i +1] = RVFiFf[12*i +1] + h*RVFiFf[12*i + 4] + 0.5*h*h*RVFiFf[12*i + 7];    /**CALCULO DE LA (n+1) POSICION**/
     RVFiFf[12*i +2] = RVFiFf[12*i +2] + h*RVFiFf[12*i + 5] + 0.5*h*h*RVFiFf[12*i + 8];           
    }
          
/******************************************************************************/     
       
    CONDC_FRONT();                                                                /*CHEQUEA LAS CONDICIONES DE FRONTERA*/
   
/******************************************************************************/

    CALCULA_FUERZAS( Posicion_nmas1  );                                           /*CALCULO DE LA FUERZA EN LA (n+1) POSICION*/

/******************************************************************************/

    for(i = 0; i < (int)pow(N,4); i ++)                                           /*CALCULO DE LA VELOCIDAD EN LA (n+1) POSICION*/ 
    {
     RVFiFf[12*i + 3] = RVFiFf[12*i + 3] + 0.5*h*( RVFiFf[12*i + 6] + RVFiFf[12*i + 9 ] );    
     RVFiFf[12*i + 4] = RVFiFf[12*i + 4] + 0.5*h*( RVFiFf[12*i + 7] + RVFiFf[12*i + 10] );
     RVFiFf[12*i + 5] = RVFiFf[12*i + 5] + 0.5*h*( RVFiFf[12*i + 8] + RVFiFf[12*i + 11] );       
    }                            
 
/******************************************************************************/ 

    for(i = 0; i < (int)pow(N,4); i ++)
    {                                                                             /*COLOCA LAS FUERZAS DE LA (n+1) POSICION*/  
     RVFiFf[ 12*i + 6 ] = RVFiFf[ 12*i + 9  ];                                          
     RVFiFf[ 12*i + 7 ] = RVFiFf[ 12*i + 10 ];                                    /*EN EL LUGAR DE LAS FUERZAS DE LA (n) POSICION*/ 
     RVFiFf[ 12*i + 8 ] = RVFiFf[ 12*i + 11 ];
    }                            
 
}



/******************************************************************************/ 



void ESCALA(long int a, double e)                                      /* ESTA FUNCION ESCALA LAS VELOCIDADES*/
{
   int i;
   long int s;
   long double beta;
  
   if( fmod(a,50) == 0 )
   {
    beta = sqrt( (T*(pow(N,4)-1))/(16*e) );    
 
    for(i = 0; i < (int)pow(N,4); i ++)
    {
     RVFiFf[ 12*i +3 ] *= beta;
     RVFiFf[ 12*i +4 ] *= beta;
     RVFiFf[ 12*i +5 ] *= beta;
    }
     
   } 
}



/******************************************************************************/



void VELOCITIES( )
{
   int i;
   double ALEATORIO1,ALEATORIO11;
   double ALEATORIO2,ALEATORIO22;
   double ALEATORIO3,ALEATORIO33;

   for( i=0; i<(int)pow(N,4); i++ )
   {
     do
      {
        
       ALEATORIO1        = aleat();
       ALEATORIO2        = aleat();
       ALEATORIO3        = aleat();
    
      }while(fabs(ALEATORIO1)+fabs(ALEATORIO2)+fabs(ALEATORIO3)>1);
      
      ALEATORIO11      = ALEATORIO1 - 0.5;
      ALEATORIO22      = ALEATORIO2 - 0.5;
      ALEATORIO33      = ALEATORIO3 - 0.5;
      
      RVFiFf[ 12*i+3 ] = ( ALEATORIO11/fabs(ALEATORIO11) )*fabs(ALEATORIO1)*V;
      RVFiFf[ 12*i+4 ] = ( ALEATORIO22/fabs(ALEATORIO22) )*fabs(ALEATORIO2)*V;
      RVFiFf[ 12*i+5 ] = ( ALEATORIO33/fabs(ALEATORIO33) )*fabs(ALEATORIO3)*V;
   }
  

}




/******************************************************************************/ 



void CONDC_FRONT()
{
   int i;
 
   for( i = 0; i < (int)pow(N,4); i++ )                                                         /*   DNM2LJ.c  ACHIVO PRINCIPAL        */
   {                                                                                /*   DNM2LJ1.c  INIT_CELDA( )          */
                                                                                    /*   DNM2LJ2.c  ENER_POTC()            */
     if( RVFiFf[ 12*i  ] > L  )RVFiFf[ 12*i   ] -= L ;                                 /*   DNM2LJ3.c  ENER_CINT()            */ 
     if( RVFiFf[ 12*i  ] < 0  )RVFiFf[ 12*i   ] += L ;                                 /*   DNM2LJ8.c  D_MAXWELL()            */ 
                                                                                    /*   DNM2LJ5.c  VEL_cmNULL()           */
                                                                                    /*   DNM2LJ6.c  MOVE_MENT( long int a )*/
     if( RVFiFf[ 12*i+1] > L  )RVFiFf[ 12*i+1 ] -= L ;                                 /*   DNM2LJ7.c  CONDC_FRONT()          */
     if( RVFiFf[ 12*i+1] < 0  )RVFiFf[ 12*i+1 ] += L ;                                 /*   DNM2LJ8.c  VELOCITIES( )          */ 
                                                                                    /*   DNM2LJ9.c  CALCULA_FUERZAS(int )     */
     if( RVFiFf[ 12*i+2 ] > L  )RVFiFf[ 12*i+2 ] -= L ;                                                                              
     if( RVFiFf[ 12*i+2 ] < 0  )RVFiFf[ 12*i+2 ] += L ;
  } 

}



/******************************************************************************/ 



void VEL_cmNULL()
{
   int i;                                                                           /*   DNM2LJ.c  ACHIVO PRINCIPAL        */
   long double Vcmx,Vcmy,Vcmz;                                                                /*   DNM2LJ1.c  INIT_CELDA( )          */
                                                                                  /*   DNM2LJ2.c  ENER_POTC()            */
   Vcmx = 0;                                                                        /*   DNM2LJ3.c  ENER_CINT()            */
   Vcmy = 0;                                                                        /*   DNM2LJ4.c  D_MAXWELL()            */
   Vcmz = 0;
                                                                                  /*   DNM2LJ5.c  VEL_cmNULL()           */
   for(i=0; i<(int)pow(N,4); i++)                                                               /*   DNM2LJ6.c  MOVE_MENT( long int a )*/
   {                                                                                /*   DNM2LJ7.c  CONDC_FRONT()          */
    Vcmx += RVFiFf[12*i + 3];                                                        /*   DNM2LJ8.c  VELOCITIES( )          */
    Vcmy += RVFiFf[12*i + 4];                                                        /*   DNM2LJ9.c  CALCULA_FUERZAS(int )     */
    Vcmz += RVFiFf[12*i + 5];
   }
 
   for(i=0; i<(int)pow(N,4); i++)
   {
    RVFiFf[12*i + 3] -= Vcmx/pow(N,4); 
    RVFiFf[12*i + 4] -= Vcmy/pow(N,4);
    RVFiFf[12*i + 5] -= Vcmz/pow(N,4);
   }

}



/******************************************************************************/ 



void ACTUALIZA_TablaVerlet()
{
 int i,j,k,l,m,No_Particulas;
 double RijCuadrado;
 
 No_Particulas = 0;
 
 for(i=0; i<(int)151*pow(N,4); i++)TablaVerlet[i] = 0;
 
 
   for(i = 0; i < (int)pow(N,4); i++)
   {
     No_Particulas = 0;    
         
   for(j = 0; j < (int)pow(N,4); j++)
   {
    if(i == j)continue;
    
    for(k = 0; k < 3; k++)
    {
    for(m = 0; m < 3; m++)
    {
    for(l = 0; l < 3; l++)
    {    
    
     RijCuadrado =
     pow( RVFiFf[12*i     ] - (RVFiFf[12*j     ] + M[ k ]*L), 2) +
     pow( RVFiFf[12*i + 1 ] - (RVFiFf[12*j + 1 ] + M[ m ]*L), 2) +
     pow( RVFiFf[12*i + 2 ] - (RVFiFf[12*j + 2 ] + M[ l ]*L), 2) ;

     if(RijCuadrado <= Rm*Rm)     /*factor 0.3*pow(10,12)*/
     {
       No_Particulas ++;            
       TablaVerlet[151*i + No_Particulas-1] = j;                                                                  
     }
      
   }
   }
   }
             
  }  
       TablaVerlet[151*i + 150] = No_Particulas;  
  }
   
     
}

void D_MAXWELL()
{
 int i;
 long double U1,U2,U3,V1,V2,V3,X1,X2,X3,S,P,ARG;

 for(i=0; i<(int)N*N*N*N; i++)
 {
 do{
    U1 = aleat();
    U2 = aleat();
    U3 = aleat();

    V1 = 2*U1 - 1;
    V2 = 2*U2 - 1;
    V3 = 2*U3 - 1;

    S = V1*V1 + V2*V2 + V3*V3;

    }while(S>=1);

    ARG = -2*log(S)/S;

    X1 = V1*sqrt(ARG);
    X2 = V2*sqrt(ARG);
    X3 = V3*sqrt(ARG);
    
    RVFiFf[12*i + 3] = 0.1*X1;
    RVFiFf[12*i + 4] = 0.1*X2;
    RVFiFf[12*i + 5] = 0.1*X3;


 }


}
