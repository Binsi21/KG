//	Program developed by
//	
//	Informatika Fakultatea
//	Euskal Herriko Unibertsitatea
//	http://www.ehu.eus/if
//
// to compile it: gcc dibujar-triangulos-y-objetos.c -lGL -lGLU -lglut -lm
//
// 
//


#include <GL/glut.h>
#include <stdio.h>
#include <string.h>
#include "cargar-triangulo.h"
#include <math.h>

#define PI 3.1415926535897932384626433

typedef struct mlist
    {
    double m[16];
    struct mlist *hptr;
    } mlist;
    
typedef struct triobj
    {
    hiruki *triptr;
    int num_triangles;
    mlist *mptr;
    struct triobj *hptr;
  //  double mesa [16];
    } triobj;

// testuraren informazioa
// información de textura

extern int load_ppm(char *file, unsigned char **bufferptr, int *dimxptr, int * dimyptr);
unsigned char *bufferra;
int dimx,dimy;

int indexx;
hiruki *triangulosptr;
triobj *foptr;
triobj *sel_ptr;
int denak;
int lineak;
int objektuak;
char aldaketa;
int ald_lokala;
//kamerarentzako
int backculling;
int normalak;
triobj *sel_ptr_k;
triobj *foptr_k;
int kam;
int objikuspuntua;
int proiekzioa; //paraleloa=0 perspektiba=1
char modua; //h=hegaldia(default) a=analisia ·G tekla·
double modelview [16]={1.0, 0.0, 0.0, 0.0,
                        0.0, 1.0, 0.0, 0.0,
                        0.0, 0.0, 1.0, 0.0,
                        0.0, 0.0, 0.0, 0.0};
double mesa [16]={1.0, 0.0, 0.0, 0.0,
                        0.0, 1.0, 0.0, 0.0,
                        0.0, 0.0, 1.0, 0.0,
                        0.0, 0.0, 0.0, 0.0};
double mperspektiba [16];

char fitxiz[100];




void objektuari_aldaketa_sartu_ezk(double m[16])
{
}



void objektuari_aldaketa_sartu_esk(double m[16])
{
}



// TODO
// funtzio honek u eta v koordenatuei dagokien pointerra itzuli behar du.
// debe devolver el pointer correspondiente a las coordenadas u y v
unsigned char * color_textura(float u, float v)
{
int desplazamendua;
char * lag;
int x_x, y_y;

x_x=u*(dimx-1);
x_x=x_x % dimx;

y_y=v*(dimy-1);
y_y=y_y % dimy;
y_y=dimy-y_y-1;

desplazamendua = dimx*y_y + x_x;
lag = (unsigned char *)bufferra;  // pixel on the left and top
return(lag+ (3*desplazamendua));
}


// TODO
// lerroa marrazten du, baina testuraren kodea egokitu behar da
// dibuja una linea pero hay que codificar la textura
void  dibujar_linea_z(int linea,float c1x, float c1z, float c1u,float c1v,float c2x,float c2z,float c2u,float c2v)
{
float xkoord,zkoord;
float u,v;
unsigned char r,g,b;
unsigned char *colorv;
float t, delt_t;

glBegin( GL_POINTS );

    // TODO
    // color_textura funtzioa ondo kodetu
    // programar de forma correcta la función color_textura
    delt_t=1.0/(c2x-c1x);
    for (t=1.0;t>=0.0;t-=delt_t)
    //for (xkoord = c1x,zkoord =c1z, u = c1u, v=c1v; xkoord <= c2x; xkoord ++)
    {      

        xkoord = c1x*t +(1-t)*c2x;
        zkoord = c1z*t +(1-t)*c2z;
        u = c1u*t +(1-t)*c2u;
        v = c1v*t +(1-t)*c2v;

        colorv=  color_textura(u, v); 
        r=colorv[0];
        g=colorv[1];
        b=colorv[2];    
        glColor3ub(r,g,b);
        glVertex3f(xkoord, linea, zkoord );
        

        // TODO 
        // zkoord, u eta v berriak kalkulatu eskuineko puntuarentzat
        // calcular zkoord, u y v del siguiente pixel
        
    }
   

glEnd();
}


void print_matrizea(double *m)
{
int i;


for (i = 0;i<4;i++)
   printf("%lf, %lf, %lf, %lf\n",m[i*4],m[i*4+1],m[i*4+2],                         m[i*4+3]);
}


// TODO
// aurrerago egitekoa
// para más adelante
void mxp(punto *pptr, double m[16], punto p)
{


pptr->x=m[0]*p.x+m[1]*p.y+m[2]*p.z+m[3];
pptr->y=m[4]*p.x+m[5]*p.y+m[6]*p.z+m[7];
pptr->z=m[8]*p.x+m[9]*p.y+m[10]*p.z+m[11];
pptr->u = p.u;
pptr->v = p.v;


}

void dibujar_triangulo(triobj *optr, int i)
{
hiruki *tptr;

punto *pgoiptr, *pbeheptr, *perdiptr;
float x1,h1,z1,u1,v1,x2,h2,z2,u2,v2,x3,h3,z3,u3,v3;
float c1x,c1z,c1u,c1v,c2x,c2z,c2u,c2v;
int linea;
float cambio1,cambio1z,cambio1u,cambio1v,cambio2,cambio2z,cambio2u,cambio2v;
punto p1,p2,p3, pp1, pp2, pp3;
punto p4, pn4, pp4;
float bcx,bcy,bcz,bc,bcu,bcv,bcuv,bccos;
float p1lag, p2lag, p3lag, p4lag, pmlag, pnlag;
float w,b,d;
punto pm, pm2, pm3, pnorm, pnorm2, pnorm3;
double lag[3];
double lag2[3];
double lag3[3];
float zatitzaile;
if (i >= optr->num_triangles) return;
tptr = optr->triptr+i;
//optr->mptr->m



kalkulatu_vnormala((tptr));
pn4.x=(tptr->p1.x + 40*tptr->normala[0]);
pn4.y=(tptr->p1.y + 40*tptr->normala[1]);
pn4.z=(tptr->p1.z + 40*tptr->normala[2]);
mxp(&p4,modelview,pn4);

mxp(&p1,modelview,tptr->p1);
mxp(&p2,modelview,tptr->p2);
mxp(&p3,modelview,tptr->p3);


pm.x=sel_ptr_k->mptr->m[3];  pm.y=sel_ptr_k->mptr->m[7];  pm.z=sel_ptr_k->mptr->m[11]; 



pnorm3.x=tptr->normala[0];   pnorm3.y=tptr->normala[1];   pnorm3.z=tptr->normala[2];
mxp(&pnorm,modelview,pnorm3);

if(proiekzioa==1)
{
    // //no deberia de pintar rayas bugeadas pero lo hace¿?
    // w=modelview[12]*tptr->p1.x+modelview[13]*tptr->p1.y +modelview[14]*tptr->p1.z+modelview[15];
    // b=modelview[12]*tptr->p2.x+modelview[13]*tptr->p2.y +modelview[14]*tptr->p2.z+modelview[15];
    // d=modelview[12]*tptr->p3.x+modelview[13]*tptr->p3.y +modelview[14]*tptr->p3.z+modelview[15];
    // if ((w<0)||(b<0)||(d<0))
    // {
    //     return;
    // }
    p1lag=-p1.z;
    p2lag=-p2.z;
    p3lag=-p3.z;

    mxp(&pp1,mperspektiba,p1);
    mxp(&pp2,mperspektiba,p2);
    mxp(&pp3,mperspektiba,p3);

    p1.x=500*pp1.x/p1lag;
    p1.y=500*pp1.y/p1lag;
    p1.z=-500*pp1.z/p1lag;

    p2.x=500*pp2.x/p2lag;
    p2.y=500*pp2.y/p2lag;
    p2.z=-500*pp2.z/p2lag;

    p3.x=500*pp3.x/p3lag;
    p3.y=500*pp3.y/p3lag;
    p3.z=-500*pp3.z/p3lag;

    p4lag=-p4.z;
    mxp(&pp4,mperspektiba,p4);
    p4.x=500*pp4.x/p4lag;
    p4.y=500*pp4.y/p4lag;
    p4.z=-500*pp4.z/p4lag;

    //back-culling
    pmlag=-pm.z;
    mxp(&pm2,mperspektiba,pm);
    pm.x=500*pm2.x/pmlag;
    pm.y=500*pm2.y/pmlag;
    pm.z=-500*pm2.z/pmlag;
    
    lag[0]=p2.x-p1.x; lag[1]=p2.y-p1.y; lag[2]=p2.z-p1.z;
    lag2[0]=p3.x-p1.x; lag2[1]=p3.y-p1.y; lag2[2]=p3.z-p1.z;

    lag3[0]=lag[1]*lag2[2] - lag[2]*lag2[1];    lag3[1]=-(lag[0]*lag2[2] - lag[2]*lag2[0]);      lag3[2]=lag[0]*lag2[1] - lag[1]*lag2[0];

    zatitzaile=sqrt((lag3[0])*(lag3[0]) + (lag3[1])*(lag3[1]) + (lag3[2])*(lag3[2]));

    pnorm.x=lag3[0]/zatitzaile;
    pnorm.y=lag3[1]/zatitzaile;
    pnorm.z=lag3[2]/zatitzaile;
    pnlag=-pnorm.z;
    mxp(&pnorm2,mperspektiba,pnorm);
    pnorm.x=500*pnorm2.x/pnlag;
    pnorm.y=500*pnorm2.y/pnlag;
    pnorm.z=-500*pnorm2.z/pnlag;
}

if (normalak == 1){
    glBegin( GL_LINES );
    glVertex3f(p1.x, p1.y, p1.z);
    glVertex3f(p4.x, p4.y, p4.z);
    glEnd();
}

//back-culling
if (lineak == 1)
{ 
    if (proiekzioa==1)   {
        bcx=pm.x - p4.x;   bcy=pm.y - p4.y;   bcz=pm.z - p4.z;
    } else  {
        bcx=sel_ptr_k->mptr->m[2];   bcy=sel_ptr_k->mptr->m[6];   bcz=sel_ptr_k->mptr->m[10];
    } 
    
    bcuv=(bcx*pnorm.x + bcy*pnorm.y + bcz*pnorm.z);
    bcu=sqrt(bcx*bcx + bcy*bcy + bcz*bcz);
    bcv=sqrt(pnorm.x*pnorm.x + pnorm.y*pnorm.y + pnorm.z*pnorm.z);
    bccos=bcuv/(bcu*bcv);
    bc=(bcu*bcv*bccos);
    if (bc<=0)
   {
        if(backculling!=1)
        { 
            glColor3ub(255,0,0);
            glBegin(GL_POLYGON);
            glVertex3d(p1.x, p1.y, p1.z);
            glVertex3d(p2.x, p2.y, p2.z);
            glVertex3d(p3.x, p3.y, p3.z);
            glEnd();
            
        } 
        return; 
   } else{
        glColor3ub(255,255,255);
        glBegin(GL_POLYGON);
        glVertex3d(p1.x, p1.y, p1.z);
        glVertex3d(p2.x, p2.y, p2.z);
        glVertex3d(p3.x, p3.y, p3.z);
        glEnd();
        return;
   } 
    
}


punto *goikoaptr, *behekoaptr, *erdikoaptr, e1, e2;

if(p1.y > p2.y){
    goikoaptr = &(p1); 
    behekoaptr = &(p2); 
}else{
    goikoaptr = &(p2); 
    behekoaptr = &(p1); 
}
if(p3.y > goikoaptr->y){
    erdikoaptr = goikoaptr; 
    goikoaptr = &(p3);    
}else{ 
    if(p3.y < behekoaptr->y){
        erdikoaptr = behekoaptr; 
        behekoaptr = &(p3); 
    }else{
        erdikoaptr = &(p3); 
    }
}

for (i=goikoaptr->y;i>erdikoaptr->y;i--){
    ebaketa_kalkulatu(goikoaptr, erdikoaptr, i, &e1);
    ebaketa_kalkulatu(goikoaptr, behekoaptr, i, &e2);
    if(e1.x<e2.x)
        dibujar_linea_z(i, e1.x, e1.z, e1.u, e1.v, e2.x, e2.z, e2.u, e2.v);
    else{
        dibujar_linea_z(i, e2.x, e2.z, e2.u, e2.v, e1.x, e1.z, e1.u, e1.v);
    }
}

for (i=behekoaptr->y;i<=erdikoaptr->y;i++){
    ebaketa_kalkulatu(behekoaptr, erdikoaptr, i, &e1);
    ebaketa_kalkulatu(behekoaptr, goikoaptr, i, &e2);
    if(e1.x<e2.x)
        dibujar_linea_z(i, e1.x, e1.z, e1.u, e1.v, e2.x, e2.z, e2.u, e2.v);
    else{
        dibujar_linea_z(i, e2.x, e2.z, e2.u, e2.v, e1.x, e1.z, e1.u, e1.v);
    }
}



}

void kalkulatu_vnormala(hiruki * hirukia){
    //N = ( (V2 – V1) x (V3 – V1) ) / || (V2 – V1) x (V3 – V1) ||
    punto pn1,pn2,pn3;
    double lag[3];
    double lag2[3];
    double lag3[3];
    double zatitzaile;

    pn1=hirukia->p1; pn2=hirukia->p2; pn3=hirukia->p3;

    lag[0]=pn2.x-pn1.x; lag[1]=pn2.y-pn1.y; lag[2]=pn2.z-pn1.z;
    lag2[0]=pn3.x-pn1.x; lag2[1]=pn3.y-pn1.y; lag2[2]=pn3.z-pn1.z;

    lag3[0]=lag[1]*lag2[2] - lag[2]*lag2[1];    lag3[1]=-(lag[0]*lag2[2] - lag[2]*lag2[0]);      lag3[2]=lag[0]*lag2[1] - lag[1]*lag2[0];

    zatitzaile=sqrt((lag3[0])*(lag3[0]) + (lag3[1])*(lag3[1]) + (lag3[2])*(lag3[2]));

    hirukia->normala[0]=lag3[0]/zatitzaile;
    hirukia->normala[1]=lag3[1]/zatitzaile;
    hirukia->normala[2]=lag3[2]/zatitzaile;
}


void ebaketa_kalkulatu(punto * goikoaptr, punto * behekoaptr, int h, punto * eptr){
    eptr->x=-((goikoaptr->x - behekoaptr->x)*(goikoaptr->y - h)/(goikoaptr->y - behekoaptr->y))+goikoaptr->x;
    eptr->y=h;
    eptr->z=-((goikoaptr->z - behekoaptr->z)*(goikoaptr->y - h)/(goikoaptr->y - behekoaptr->y))+goikoaptr->z;
    eptr->u=-((goikoaptr->u - behekoaptr->u)*(goikoaptr->y - h)/(goikoaptr->y - behekoaptr->y))+goikoaptr->u;
    eptr->v=-((goikoaptr->v - behekoaptr->v)*(goikoaptr->y - h)/(goikoaptr->y - behekoaptr->y))+goikoaptr->v;
}

static void marraztu(void)
{
float u,v;
int i,j;
triobj *auxptr;

/*
unsigned char* colorv;
unsigned char r,g,b;
*/

  // marrazteko objektuak behar dira
  // no se puede dibujar sin objetos
if (foptr ==0) return;

// clear viewport...
if (objektuak == 1) glClear( GL_DEPTH_BUFFER_BIT|GL_COLOR_BUFFER_BIT );
    else 
      {
      if (denak == 0) glClear( GL_DEPTH_BUFFER_BIT|GL_COLOR_BUFFER_BIT );
      }

glMatrixMode(GL_PROJECTION);
glLoadIdentity();
// la penultimaestaba en -500 hay que ponerla en 0??
glOrtho(-500.0, 500.0, -500.0, 500.0, -500.0, 500.0);
//printf("kameraren matrizea: \n");      
//print_matrizea(sel_ptr_k->mptr->m);

                
kalkulatu_mesa(sel_ptr_k->mptr->m,&mesa[0]);
//perspektiba_aldatu(&mperspektiba[0]);
//printf("KAMERAREB matrizea: \n");      
//print_matrizea(sel_ptr_k->mptr->m);

triangulosptr = sel_ptr->triptr;
if (objektuak == 1)
    {
    if (denak == 1)
        {
        for (auxptr =foptr; auxptr != 0; auxptr = auxptr->hptr)
            {
                //printf("mesa matrizea: \n");      
                //print_matrizea(mesa);

                //printf("modelview matrizea: \n");
            biderkatu_matrizeak(&modelview[0], &mesa[0], auxptr->mptr->m);
            //print_matrizea(modelview);

            for (i =0; i < auxptr->num_triangles; i++)
                {
                dibujar_triangulo(auxptr,i);
                }
            }
        
        for (auxptr =foptr_k; auxptr != 0; auxptr = auxptr->hptr)
            {
            biderkatu_matrizeak(&modelview[0], &mesa[0], auxptr->mptr->m);
            for (i =0; i < auxptr->num_triangles; i++)
                {
                dibujar_triangulo(auxptr,i);
                }
            }
        }
    else
        {
        
        for (i =0; i < sel_ptr->num_triangles; i++)
            {
            dibujar_triangulo(sel_ptr,i);
            }
        }
    }
else
    {
     dibujar_triangulo(sel_ptr,indexx);
    }
glFlush();
}



void read_from_file(char *fitx, triobj **selptrptr, triobj **foptrptr)
{
int i,j,retval;
triobj *optr;
unsigned char *kolptr;
// hiruki *tptr;

    //printf("%s fitxategitik datuak hartzera\n",fitx);
    optr = (triobj *)malloc(sizeof(triobj));
    kolptr=0;
    retval = cargar_triangulos_color(fitx, &(optr->num_triangles), &(optr->triptr), &kolptr);
    //retval = cargar_triangulos_color(fitx, &(optr->num_triangles), &(optr->triptr), &kol);
    if ((retval !=9)&&(retval!=15)  ) 
         {
         printf("%s fitxategitik datuak hartzerakoan arazoak izan ditut\n    Problemas al leer\n",fitxiz);
         free(optr);
         }
       else
         {
         triangulosptr = optr->triptr;
         //printf("objektuaren matrizea...\n");
         optr->mptr = (mlist *)malloc(sizeof(mlist));
         for (i=0; i<16; i++) optr->mptr->m[i] =0;
         optr->mptr->m[0] = 1.0;
         optr->mptr->m[5] = 1.0;  
         optr->mptr->m[10] = 1.0;  
         optr->mptr->m[15] = 1.0;   
        
        // for (j=0;j<(optr->num_triangles);j++)
        // { 
        //     tptr = optr->triptr+j;
        //     kalkulatu_vnormala((tptr));
        //     j++;
        // } 

         optr->hptr = *selptrptr;
         *selptrptr = optr;
         *foptrptr = optr;
         }
     printf("datuak irakurrita\nLecura finalizada\n");
}

void biderkatu_matrizeak(double *emaitzaptr, double *mlag, double *eskptr)
{
    int i,j,k;
    double era=0.0;

    for (i=0; i<4; i++)
    {
        for (j=0; j<4; j++)
        {
            era=0.0;
            for(k=0; k<4; k++)
            {
                era+=mlag[4*i+k]*eskptr[4*k+j];
            }
            emaitzaptr[i*4+j]=era; 
        }
        
    }
}

void kalkulatu_mesa(double matriz[16], double *mesa)
{
    mesa[0]=matriz[0];   mesa[1]=matriz[4];     mesa[2]=matriz[8];      mesa[3]=-(matriz[0]*matriz[3] + matriz[4]*matriz[7] + matriz[8]*matriz[11]);
    mesa[4]=matriz[1];   mesa[5]=matriz[5];     mesa[6]=matriz[9];      mesa[7]=-(matriz[1]*matriz[3] + matriz[5]*matriz[7] + matriz[9]*matriz[11]);
    mesa[8]=matriz[2];   mesa[9]=matriz[6];     mesa[10]=matriz[10];    mesa[11]=-(matriz[2]*matriz[3] + matriz[6]*matriz[7] + matriz[10]*matriz[11]);
    mesa[12]=0;          mesa[13]=0;            mesa[14]=0;             mesa[15]=1;
}


void x_aldaketa(int dir)
{
    double biraketa = 0.2618; //15 gradu
    double mlag[16];
    double mlag2[16];
    double mlag3[16];
    double mlag4[16];
    double mlag5[16];
    double mlag6[16];
    double mbirak[16];
    double Xx;
    double Xy;
    double Xz;
    mlist * mleb;      
    
    mlag[0]=1.0;    mlag[1]=0.0;    mlag[2]=0.0;    mlag[3]=0.0;
    mlag[4]=0.0;                                    mlag[7]=0.0;
    mlag[8]=0.0;                                    mlag[11]=0.0;
    mlag[12]=0.0;   mlag[13]=0.0;   mlag[14]=0.0;   mlag[15]=1.0;
 
    mlag2[0]= 1.0;    mlag2[1]= 0.0;    mlag2[2]= 0.0;
    mlag2[4]= 0.0;    mlag2[5]= 1.0;    mlag2[6]= 0.0;    mlag2[7]= 0.0;
    mlag2[8]= 0.0;    mlag2[9]= 0.0;    mlag2[10]= 1.0;   mlag2[11]= 0.0;
    mlag2[12]= 0.0;   mlag2[13]= 0.0;   mlag2[14]= 0.0;   mlag2[15]= 1.0;

    mlag3[0]= 1.0;    mlag3[1]= 0.0;    mlag3[2]= 0.0;    
    mlag3[4]= 0.0;    mlag3[5]= 1.0;    mlag3[6]= 0.0;    mlag3[7]= 0.0;
    mlag3[8]= 0.0;    mlag3[9]= 0.0;    mlag3[10]= 1.0;   mlag3[11]= 0.0;
    mlag3[12]= 0.0;   mlag3[13]= 0.0;   mlag3[14]= 0.0;   mlag3[15]= 1.0;

    mlag4[0]= 1.0;    mlag4[1]= 0.0;    mlag4[2]= 0.0;    
    mlag4[4]= 0.0;    mlag4[5]= 1.0;    mlag4[6]= 0.0;    
    mlag4[8]= 0.0;    mlag4[9]= 0.0;    mlag4[10]= 1.0;   
    mlag4[12]= 0.0;   mlag4[13]= 0.0;   mlag4[14]= 0.0;   mlag4[15]= 1.0;  


    Xx=(sel_ptr_k->mptr->m[0]); Xy=(sel_ptr_k->mptr->m[4]); Xz=(sel_ptr_k->mptr->m[8]);

    mbirak[0]= cos(biraketa)+(1-cos(biraketa))*(Xx*Xx);      mbirak[1]= (1-cos(biraketa))*(Xx*Xy)-Xz*sin(biraketa);    mbirak[2]= (1-cos(biraketa))*(Xx*Xz)+Xy*sin(biraketa);    mbirak[3]= 0.0;  
    mbirak[4]= (1-cos(biraketa))*(Xx*Xy)+Xz*sin(biraketa);   mbirak[5]= cos(biraketa)+(1-cos(biraketa))*(Xy*Xy);       mbirak[6]= (1-cos(biraketa))*(Xy*Xz)-Xx*sin(biraketa);    mbirak[7]= 0.0;  
    mbirak[8]= (1-cos(biraketa))*(Xx*Xz)-Xy*sin(biraketa);   mbirak[9]= (1-cos(biraketa))*(Xy*Xz)+Xx*sin(biraketa);    mbirak[10]= cos(biraketa)+(1-cos(biraketa))*(Xz*Xz);      mbirak[11]= 0.0;  
    mbirak[12]= 0.0;   mbirak[13]= 0.0;   mbirak[14]= 0.0;   mbirak[15]= 1.0;  

    if (dir==1)
    {
        mlag[5]=cos(biraketa);
        mlag[6]=-sin(biraketa);
        mlag[9]=sin(biraketa);
        mlag[10]=cos(biraketa);
        mlag2[3]= 3.3;
        mlag3[3]= 2.0;
    }
    else
    {
        biraketa = 2*PI - biraketa;
        mlag[5]=cos(biraketa);
        mlag[6]=-sin(biraketa);
        mlag[9]=sin(biraketa);
        mlag[10]=cos(biraketa);
        mlag2[3]= -3.3;
        mlag3[3]= -2.0;
    }
    mleb = (mlist *)malloc(sizeof(mlist));
    if(kam==0)
    {
        if(aldaketa=='r')
        {
            biderkatu_matrizeak(mleb->m, sel_ptr->mptr->m, mlag);
        }else{ //translazioa
            if (ald_lokala==1)
            {
                biderkatu_matrizeak(mleb->m, sel_ptr->mptr->m, mlag2);
            }else{
                biderkatu_matrizeak(mleb->m, mlag2, sel_ptr->mptr->m);
            }
        } 
        mleb->hptr = sel_ptr->mptr;
        sel_ptr->mptr = mleb;    
    }else{ //kamara aldatu
        if(modua=='a')
        {
            mlag4[3]=-(sel_ptr->mptr->m[3]);
            mlag4[7]=-(sel_ptr->mptr->m[7]);
            mlag4[11]=-(sel_ptr->mptr->m[11]);
            biderkatu_matrizeak(mlag5, mbirak, mlag4);
            mlag4[3]=(sel_ptr->mptr->m[3]);
            mlag4[7]=(sel_ptr->mptr->m[7]);
            mlag4[11]=(sel_ptr->mptr->m[11]);
            biderkatu_matrizeak(mlag6, mlag4, mlag5);            
            biderkatu_matrizeak(mleb->m, mlag6, sel_ptr_k->mptr->m);
        }
        else if(modua=='h')
        {
            if(aldaketa=='r')
            {
                biderkatu_matrizeak(mleb->m, sel_ptr_k->mptr->m, mlag);
            }else{ //translazioa
                if (ald_lokala==1)
                {
                    biderkatu_matrizeak(mleb->m, sel_ptr_k->mptr->m, mlag3);
                }else{
                    biderkatu_matrizeak(mleb->m, mlag3, sel_ptr_k->mptr->m);
                }
            } 
        }
        mleb->hptr = sel_ptr_k->mptr;
        sel_ptr_k->mptr = mleb;    
    }
}


void y_aldaketa(int dir)
{
    double biraketa = 0.2618; //15 gradu
    double mlag[16];
    double mlag2[16];
    double mlag3[16];
    double mlag4[16];
    double mlag5[16];
    double mlag6[16];
    double mbirak[16];
    double Yx;
    double Yy;
    double Yz;
    mlist * mleb;

                    mlag[1]=0.0;                    mlag[3]=0.0;
    mlag[4]=0.0;    mlag[5]=1.0;    mlag[6]=0.0;    mlag[7]=0.0;
                    mlag[9]=0.0;                    mlag[11]=0.0;
    mlag[12]=0.0;   mlag[13]=0.0;   mlag[14]=0.0;   mlag[15]=1.0;

    mlag2[0]= 1.0;  mlag2[1]= 0.0;  mlag2[2]= 0.0;  mlag2[3]= 0.0;
    mlag2[4]= 0.0;  mlag2[5]= 1.0;  mlag2[6]= 0.0;
    mlag2[8]= 0.0;  mlag2[9]= 0.0;  mlag2[10]= 1.0; mlag2[11]= 0.0;
    mlag2[12]= 0.0; mlag2[13]= 0.0; mlag2[14]= 0.0; mlag2[15]= 1.0;

    mlag3[0]= 1.0;    mlag3[1]= 0.0;    mlag3[2]= 0.0;    mlag3[3]= 0.0;
    mlag3[4]= 0.0;    mlag3[5]= 1.0;    mlag3[6]= 0.0;    
    mlag3[8]= 0.0;    mlag3[9]= 0.0;    mlag3[10]= 1.0;   mlag3[11]= 0.0;
    mlag3[12]= 0.0;   mlag3[13]= 0.0;   mlag3[14]= 0.0;   mlag3[15]= 1.0;

    mlag4[0]= 1.0;    mlag4[1]= 0.0;    mlag4[2]= 0.0;    
    mlag4[4]= 0.0;    mlag4[5]= 1.0;    mlag4[6]= 0.0;    
    mlag4[8]= 0.0;    mlag4[9]= 0.0;    mlag4[10]= 1.0;   
    mlag4[12]= 0.0;   mlag4[13]= 0.0;   mlag4[14]= 0.0;   mlag4[15]= 1.0;  


    Yx=(sel_ptr_k->mptr->m[1]); Yy=(sel_ptr_k->mptr->m[5]); Yz=(sel_ptr_k->mptr->m[9]);

    mbirak[0]= cos(biraketa)+(1-cos(biraketa))*(Yx*Yx);      mbirak[1]= (1-cos(biraketa))*(Yx*Yy)-Yz*sin(biraketa);    mbirak[2]= (1-cos(biraketa))*(Yx*Yz)+Yy*sin(biraketa);    mbirak[3]= 0.0;  
    mbirak[4]= (1-cos(biraketa))*(Yx*Yy)+Yz*sin(biraketa);   mbirak[5]= cos(biraketa)+(1-cos(biraketa))*(Yy*Yy);       mbirak[6]= (1-cos(biraketa))*(Yy*Yz)-Yx*sin(biraketa);    mbirak[7]= 0.0;  
    mbirak[8]= (1-cos(biraketa))*(Yx*Yz)-Yy*sin(biraketa);   mbirak[9]= (1-cos(biraketa))*(Yy*Yz)+Yx*sin(biraketa);    mbirak[10]= cos(biraketa)+(1-cos(biraketa))*(Yz*Yz);      mbirak[11]= 0.0;  
    mbirak[12]= 0.0;   mbirak[13]= 0.0;   mbirak[14]= 0.0;   mbirak[15]= 1.0;  

    if (dir==1)
    {
        mlag[0]=cos(biraketa);
        mlag[2]=sin(biraketa);
        mlag[8]=-sin(biraketa);
        mlag[10]=cos(biraketa);
        mlag2[7]= 3.3;
        mlag3[7]= 2.0;
    }
    else
    {
        biraketa = 2*PI - biraketa;
        mlag[0]=cos(biraketa);
        mlag[2]=sin(biraketa);
        mlag[8]=-sin(biraketa);
        mlag[10]=cos(biraketa);
        mlag2[7]= -3.3;
        mlag3[7]= -2.0;
    }
    mleb = (mlist *)malloc(sizeof(mlist));
    if(kam==0)
    {
        if(aldaketa=='r')
        {
            biderkatu_matrizeak(mleb->m, sel_ptr->mptr->m, mlag);
        }else{
            biderkatu_matrizeak(mleb->m, sel_ptr->mptr->m, mlag2);
        }  
        mleb->hptr = sel_ptr->mptr;
        sel_ptr->mptr = mleb;
   }else{ //kamara aldatu
        if(modua=='a')
        {
            mlag4[3]=-(sel_ptr->mptr->m[3]);
            mlag4[7]=-(sel_ptr->mptr->m[7]);
            mlag4[11]=-(sel_ptr->mptr->m[11]);
            biderkatu_matrizeak(mlag5, mbirak, mlag4);
            mlag4[3]=(sel_ptr->mptr->m[3]);
            mlag4[7]=(sel_ptr->mptr->m[7]);
            mlag4[11]=(sel_ptr->mptr->m[11]);
            biderkatu_matrizeak(mlag6, mlag4, mlag5);
            biderkatu_matrizeak(mleb->m, mlag6, sel_ptr_k->mptr->m);
        }
        else if(modua=='h')
        {
            if(aldaketa=='r')
            {
                biderkatu_matrizeak(mleb->m, sel_ptr_k->mptr->m, mlag);
            }else{ //translazioa
                if (ald_lokala==1)
                {
                    biderkatu_matrizeak(mleb->m, sel_ptr_k->mptr->m, mlag3);
                }else{
                    biderkatu_matrizeak(mleb->m, mlag3, sel_ptr_k->mptr->m);
                } 
            }
        }
        mleb->hptr = sel_ptr_k->mptr;
        sel_ptr_k->mptr = mleb;    
    }    
}



void z_aldaketa(int dir)
{   
    double biraketa = 0.2618; //15 gradu
    double mlag[16];
    double mlag2[16];
    double mlag3[16];
    mlist * mleb;

                                    mlag[2]=0.0;    mlag[3]=0.0;
                                    mlag[6]=0.0;    mlag[7]=0.0;
    mlag[8]=0.0;    mlag[9]=0.0;    mlag[10]=1.0;   mlag[11]=0.0;
    mlag[12]=0.0;   mlag[13]=0.0;   mlag[14]=0.0;   mlag[15]=1.0;

    mlag2[0]= 1.0;    mlag2[1]= 0.0;    mlag2[2]= 0.0;    mlag2[3]= 0.0;
    mlag2[4]= 0.0;    mlag2[5]= 1.0;    mlag2[6]= 0.0;    mlag2[7]= 0.0;
    mlag2[8]= 0.0;    mlag2[9]= 0.0;    mlag2[10]= 1.0;
    mlag2[12]= 0.0;   mlag2[13]= 0.0;   mlag2[14]= 0.0;   mlag2[15]= 1.0;

    mlag3[0]= 1.0;    mlag3[1]= 0.0;    mlag3[2]= 0.0;    mlag2[3]= 0.0;
    mlag3[4]= 0.0;    mlag3[5]= 1.0;    mlag3[6]= 0.0;    mlag3[7]= 0.0;
    mlag3[8]= 0.0;    mlag3[9]= 0.0;    mlag3[10]= 1.0;   
    mlag3[12]= 0.0;   mlag3[13]= 0.0;   mlag3[14]= 0.0;   mlag3[15]= 1.0;

    if (dir==1)
    {
        mlag[0]=cos(biraketa);
        mlag[1]=-sin(biraketa);
        mlag[4]=sin(biraketa);
        mlag[5]=cos(biraketa);

        mlag2[11]= 3.3;
        mlag3[11]= 2.0;
    }
    else
    {
        biraketa = 2*PI - biraketa;
        mlag[0]=cos(biraketa);
        mlag[1]=-sin(biraketa);
        mlag[4]=sin(biraketa);
        mlag[5]=cos(biraketa);

        mlag2[11]= -3.3;
        mlag3[11]= -2.0;
    }
    mleb = (mlist *)malloc(sizeof(mlist));
    
    if(kam==0)
    {
        if(aldaketa=='r')
        {
            biderkatu_matrizeak(mleb->m, sel_ptr->mptr->m, mlag);
        }else{
            biderkatu_matrizeak(mleb->m, sel_ptr->mptr->m, mlag2);
        }
        mleb->hptr = sel_ptr->mptr;
        sel_ptr->mptr = mleb;
   }else{ //kamara aldatu
        if(modua=='a')
        {
            if (ald_lokala==1)
            {
                biderkatu_matrizeak(mleb->m, sel_ptr_k->mptr->m, mlag3);
            }else{
                biderkatu_matrizeak(mleb->m, mlag3, sel_ptr_k->mptr->m);
            }
        } 
        mleb->hptr = sel_ptr_k->mptr;
        sel_ptr_k->mptr = mleb;    
    }    
}

void eskalatu(int dir)        

{
    double a;
    double mlag[16];
    mlist * mleb;
    if(dir==0)
   {
        a=1.5;
   } else{
        a=0.66;
   }
   int i;
   for (i=0; i<16; i++)
  {
    mlag[i] = 0;
  } 
    mlag[0] = a;
    mlag[5] = a;
    mlag[10] = a;
    mlag[15] = 1;
    
    mleb = (mlist *)malloc(sizeof(mlist));
    if(kam==0){
        biderkatu_matrizeak(mleb->m, sel_ptr->mptr->m, mlag);
        sel_ptr->mptr = mleb;
    }else{
        biderkatu_matrizeak(mleb->m, sel_ptr_k->mptr->m, mlag);
        sel_ptr_k->mptr = mleb;
    }
    
}

void undo()
{
    double im[16];
    im[0]=1.0;    im[1]=0.0;    im[2]=0.0;    im[3]=0.0;
    im[4]=0.0;    im[5]=1.0;    im[6]=0.0;    im[7]=0.0;
    im[8]=0.0;    im[9]=0.0;    im[10]=1.0;   im[11]=0.0;
    im[12]=0.0;   im[13]=0.0;   im[14]=0.0;   im[15]=1.0;
    if(kam==0)
    {
        if(sel_ptr->mptr->hptr->m != im)
            sel_ptr->mptr = sel_ptr->mptr->hptr;
    }else{
        if(sel_ptr_k->mptr->hptr->m != im)
            sel_ptr_k->mptr = sel_ptr_k->mptr->hptr;
    }
        
}

void perspektiba_aldatu(double *mperspektiba)
{
    float n,r,t;
    float l,b;
    float f;
    n=5.0; r=5.0; t=5.0;
    l=-5.0; b=-5.0;
    f=500.0;
    mperspektiba[0]=(2*n)/(r-l);     mperspektiba[1]=0.0;    mperspektiba[2]=(r+l)/(r-l);   mperspektiba[3]=0.0;
    mperspektiba[4]=0.0;            mperspektiba[5]=(2*n)/(t-b);  mperspektiba[6]=(t+b)/(t-b);     mperspektiba[7]=0.0;
    mperspektiba[8]=0.0;    mperspektiba[9]=0.0;    mperspektiba[10]=(-(f+n))/(f-n);   mperspektiba[11]=(-2*(f+n))/(f-n); 
    mperspektiba[12]=0.0;   mperspektiba[13]=0.0;   mperspektiba[14]=-1.0;   mperspektiba[15]=0.0;
}

void objektuari_begiratu()
{
    double At[3];
    double E[3];
    double Vup[3];
    double Xc[3];
    double Yc[3];
    double Zc[3];
    double klag[16];
    double lag[3];
    double zatitzaile;
    int i;
   
    At[0]=sel_ptr->mptr->m[3];  At[1]=sel_ptr->mptr->m[7];  At[2]=sel_ptr->mptr->m[11];
    E[0]=sel_ptr_k->mptr->m[3];  E[1]=sel_ptr_k->mptr->m[7];  E[2]=sel_ptr_k->mptr->m[11]; 
    Vup[0]=sel_ptr_k->mptr->m[1];  Vup[1]=sel_ptr_k->mptr->m[5];  Vup[2]=sel_ptr_k->mptr->m[9];  

    lag[0]=E[0]-At[0]; lag[1]=E[1]-At[1]; lag[2]=E[2]-At[2];
    zatitzaile=sqrt((lag[0])*(lag[0]) + (lag[1])*(lag[1]) + (lag[2])*(lag[2]));
    Zc[0]=lag[0]/zatitzaile; Zc[1]=lag[1]/zatitzaile; Zc[2]=lag[2]/zatitzaile; 

    lag[0]=Vup[1]*Zc[2]-Vup[2]*Zc[1];   lag[1]=-(Vup[0]*Zc[2]-Vup[2]*Zc[0]);    lag[2]=Vup[0]*Zc[1]-Vup[1]*Zc[0]; 
    zatitzaile=sqrt((lag[0])*(lag[0]) + (lag[1])*(lag[1]) + (lag[2])*(lag[2]));
    Xc[0]=lag[0]/zatitzaile; Xc[1]=lag[1]/zatitzaile; Xc[2]=lag[2]/zatitzaile; 

    Yc[0]=Zc[1]*Xc[2]-Zc[2]*Xc[1];   Yc[1]=-(Zc[0]*Xc[2]-Zc[2]*Xc[0]);    Yc[2]=Zc[0]*Xc[1]-Zc[1]*Xc[0]; 

    klag[0]=Xc[0];  klag[1]=Yc[0];  klag[2]=Zc[0];  klag[3]=sel_ptr_k->mptr->m[3];
    klag[4]=Xc[1];  klag[5]=Yc[1];  klag[6]=Zc[1];  klag[7]=sel_ptr_k->mptr->m[7];
    klag[8]=Xc[2];  klag[9]=Yc[2];  klag[10]=Zc[2]; klag[11]=sel_ptr_k->mptr->m[11];
    klag[12]=0.0;   klag[13]=0.0;   klag[14]=0.0;   klag[15]=1.0;

    for (i=0; i<16; i++)
    {
        sel_ptr_k->mptr->m[i]=klag[i];
    }
    
}



// This function will be called whenever the user pushes one key
static void teklatua (unsigned char key, int x, int y)
{
int retval;
int i;
FILE *obj_file;

switch(key)
	{
	case 13: 
	        if (foptr != 0)  // objekturik ez badago ezer ez du egin behar
	                         // si no hay objeto que no haga nada
	            {
	            indexx ++;  // azkena bada lehenengoa bihurtu
		                // pero si es el último? hay que controlarlo!
		    if (indexx == sel_ptr->num_triangles) 
		        {
		        indexx = 0;
		        if ((denak == 1) && (objektuak == 0))
		            {
		            glClear( GL_DEPTH_BUFFER_BIT|GL_COLOR_BUFFER_BIT );
		            glFlush();
		            }
		        }
		    }
		break;
	case 'd':
		if (denak == 1) denak = 0;
		    else denak = 1;
		break;
    case 'n':
		if (normalak == 1) normalak = 0;
		    else normalak = 1;
		break;
    case 'c':
        if(objikuspuntua==0){
            if (kam==1) kam=0;
                else kam=1;
        }
        break;
    case 'C':
        if(sel_ptr_k==sel_ptr) 
        {
            sel_ptr_k=foptr_k;//si esto no va hay que crear klag
            objikuspuntua=0;
        }
        else
        {
            sel_ptr_k=sel_ptr;
            objikuspuntua=1;
        }
        break;
    case 'G':
        if (modua=='h') {modua='a'; printf("Orain analisi moduan nago\n"); objektuari_begiratu();}
            else {modua='h'; printf("Orain hegaldi moduan nago\n");}
        break; 
	case 'o':
		if (objektuak == 1) objektuak = 0;
		    else objektuak = 1;
		break;
	case 'l':
		if (lineak == 1) lineak = 0;
		    else lineak = 1;
		break;
    case 'b':
		if (backculling == 1) backculling = 0;
		    else backculling = 1;
		break;
	case 't':
	    aldaketa = 't';
		break;
    case 'p':
        if (proiekzioa == 1) proiekzioa = 0;
		    else proiekzioa = 1;
        break;
	case 'r':
		aldaketa = 'r';
		break;
	case 'g':
		if (ald_lokala == 1)
        {
            ald_lokala = 0;
            printf("Aldaketak orain globalak izango dira\n");
        } 
		else 
        {
            ald_lokala = 1;
            printf("Aldaketak orain lokalak izango dira\n");
        }
		break;
        case 'x':
                x_aldaketa(1);
                break;
        case 'y':
                y_aldaketa(1);
                break;
        case 'z':
                z_aldaketa(1);
                break;
        case 'X':
                x_aldaketa(0);
                break;
        case 'Y':
                y_aldaketa(0);
                break;
        case 'Z':
                z_aldaketa(0);
                break;
        case 'a':
                eskalatu(0);
                break;
        case 'A':
                eskalatu(1);
                break;
        case 'u':
                undo();
                break;
	case 'f':
	        /*Ask for file*/
	        printf("idatzi fitxategi izena\n"); 
	        scanf("%s", &(fitxiz[0]));
            read_from_file(fitxiz, &sel_ptr, &foptr);
	        indexx = 0;
                break;
       /* case 'S':  // save to file
	        printf("idatzi fitxategi izena\n"); 
	        scanf("%s", &(fitxiz[0]));
                if ((obj_file = fopen(fitxiz, "w")) == NULL)
                         {
                         printf("ezin fitxategia ireki\n");
                         }
                     else
                         {
                         for (i =0; i < sel_ptr->num_triangles; i++)
                            {
                            fprintf(obj_file,"t %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",
                                 sel_ptr->triptr[i].p1.x-250, sel_ptr->triptr[i].p1.y-250, sel_ptr->triptr[i].p1.z, 
                                 sel_ptr->triptr[i].p1.u, sel_ptr->triptr[i].p1.v,
                                 sel_ptr->triptr[i].p2.x-250, sel_ptr->triptr[i].p2.y-250, sel_ptr->triptr[i].p2.z, 
                                 sel_ptr->triptr[i].p2.u, sel_ptr->triptr[i].p2.v,
                                 sel_ptr->triptr[i].p3.x-250, sel_ptr->triptr[i].p3.y-250, sel_ptr->triptr[i].p3.z, 
                                 sel_ptr->triptr[i].p3.u, sel_ptr->triptr[i].p3.v );
                            }
                         fclose(obj_file);
                         }
                break; */
        case 9: /* <TAB> */
            if (kam==0)
            {
                if (foptr != 0) // objekturik gabe ez du ezer egin behar
                                // si no hay objeto no hace nada
                    {
                    sel_ptr = sel_ptr->hptr;
                    /*The selection is circular, thus if we move out of the list we go back to the first element*/
                    if (sel_ptr == 0) sel_ptr = foptr;
                    indexx =0; // the selected polygon is the first one
                    }
            }else{
                if (foptr_k != 0) 
                    {
                    sel_ptr_k = sel_ptr_k->hptr;
                    if (sel_ptr_k == 0) sel_ptr_k = foptr_k;
                    indexx =0;
                    }
            if (modua=='a') objektuari_begiratu();  //hay que probar si va esto + nose si va aqui
            }
            break;
	case 27:  // <ESC>
		exit( 0 );
		break;
	default:
		printf("%d %c\n", key, key );
	}

// The screen must be drawn to show the new triangle
glutPostRedisplay();
}

int main(int argc, char** argv)
{
int retval;

	printf(" Triangeluak: barneko puntuak eta testura\n Triángulos con puntos internos y textura \n");
	printf("Press <ESC> to finish\n");
	glutInit(&argc,argv);
	glutInitDisplayMode ( GLUT_RGB|GLUT_DEPTH );
	glutInitWindowSize ( 500, 500 );
	glutInitWindowPosition ( 100, 100 );
	glutCreateWindow( "KBG/GO praktika" );
	
	glutDisplayFunc( marraztu );
	glutKeyboardFunc( teklatua );
	/* we put the information of the texture in the buffer pointed by bufferra. The dimensions of the texture are loaded into dimx and dimy */ 
        retval = load_ppm("testura.ppm", &bufferra, &dimx, &dimy);
        if (retval !=1) 
            {
            printf("Ez dago texturaren fitxategia (testura.ppm)\n");
            exit(-1);
            }
        
	glClearColor( 0.0f, 0.0f, 0.7f, 1.0f );
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glEnable(GL_DEPTH_TEST); // activar el test de profundidad (Z-buffer)
        denak = 1;
        lineak =1;
        objektuak = 1;
        foptr = 0;
        sel_ptr = 0;
        sel_ptr_k=0;
        foptr_k=0;
        kam=0;
        modua = 'h';
        objikuspuntua=0;
        proiekzioa=0;
        aldaketa = 'r';
        ald_lokala = 1;
        normalak=0; 
        backculling=0;

        mperspektiba[0]=1.0; mperspektiba[1]=0.0; mperspektiba[2]=0.0; mperspektiba[3]=0.0; 
        mperspektiba[4]=0.0; mperspektiba[5]=1.0; mperspektiba[6]=0.0; mperspektiba[7]=0.0; 
        mperspektiba[8]=0.0; mperspektiba[9]=0.0; mperspektiba[10]=-1.02; mperspektiba[11]=-10.1; 
        mperspektiba[12]=0.0; mperspektiba[13]=0.0; mperspektiba[14]=-1.0; mperspektiba[15]=0.0; 
                       
        if (argc>1)
                read_from_file(argv[1], &sel_ptr, &foptr);  
            else{ read_from_file("borobila.txt", &sel_ptr, &foptr); sel_ptr->mptr->m[3]=-120;}
        read_from_file("borobila.txt", &sel_ptr, &foptr);
        sel_ptr->mptr->m[3]= 120;
        read_from_file("c.txt", &sel_ptr_k, &foptr_k);
        sel_ptr_k->mptr->m[11]= 200;
	glutMainLoop();

	return 0;   
}
