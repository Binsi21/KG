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
triobj *sel_ptr_k;
triobj *foptr_k;
int kam;

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


void print_matrizea(char *str)
{
int i;

printf("%s\n",str);
for (i = 0;i<4;i++)
   printf("%lf, %lf, %lf, %lf\n",sel_ptr->mptr->m[i*4],sel_ptr->mptr->m[i*4+1],sel_ptr->mptr->m[i*4+2],
                                 sel_ptr->mptr->m[i*4+3]);
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
punto p1,p2,p3;

if (i >= optr->num_triangles) return;
tptr = optr->triptr+i;
mxp(&p1,optr->mptr->m,tptr->p1);
mxp(&p2,optr->mptr->m,tptr->p2);
mxp(&p3,optr->mptr->m,tptr->p3);
if (lineak == 1)
        {
        glBegin(GL_POLYGON);
        glVertex3d(p1.x, p1.y, p1.z);
        glVertex3d(p2.x, p2.y, p2.z);
        glVertex3d(p3.x, p3.y, p3.z);
        glEnd();
        return;
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
glOrtho(-500.0, 500.0, -500.0, 500.0, 0.0, 500.0);


triangulosptr = sel_ptr->triptr;
if (objektuak == 1)
    {
    if (denak == 1)
        {
        for (auxptr =foptr; auxptr != 0; auxptr = auxptr->hptr)
            {
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



void read_from_file(char *fitx)
{
int i,retval;
triobj *optr;

    //unsigned char *kol;
    //kol[0]=0;
    //kol[1]=0;
    //kol[2]=0;

    //printf("%s fitxategitik datuak hartzera\n",fitx);
    optr = (triobj *)malloc(sizeof(triobj));
    retval = cargar_triangulos(fitx, &(optr->num_triangles), &(optr->triptr));
    //retval = cargar_triangulos_color(fitx, &(optr->num_triangles), &(optr->triptr), &kol);
    if (retval !=1) 
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
        
         if(fitx=="k.txt"){
            foptr_k = optr;
            sel_ptr_k = optr;
            kam=1;
         }else{
            foptr = optr;
            sel_ptr = optr;
            kam=0;
         }
         }
     printf("datuak irakurrita\nLecura finalizada\n");
}

void biderkatu_matrizeak(double emaitzaptr[16], double mlag[16], double eskptr[16])
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

void kalkulatu_mesa(double matriz[16], double mesa[16])
{
    mesa[0]=matriz[0];   mesa[1]=matriz[4];     mesa[2]=matriz[8];      mesa[3]=matriz[0]*matriz[3] + matriz[1]*matriz[7] + matriz[2]*matriz[11];
    mesa[4]=matriz[1];   mesa[5]=matriz[5];     mesa[6]=matriz[9];      mesa[7]=matriz[4]*matriz[3] + matriz[5]*matriz[7] + matriz[6]*matriz[11];
    mesa[8]=matriz[2];   mesa[9]=matriz[6];     mesa[10]=matriz[10];    mesa[11]=matriz[8]*matriz[3] + matriz[9]*matriz[7] + matriz[10]*matriz[11];
    mesa[12]=matriz[12];  mesa[13]=matriz[13];  mesa[14]=matriz[14];    mesa[15]=matriz[15];
}


void x_aldaketa(int dir)
{
    double biraketa = 0.2618; //15 gradu
    double mlag[16];
    double mlag2[16];
    mlist * mleb;      
    
    mlag[0]=1.0;    mlag[1]=0.0;    mlag[2]=0.0;    mlag[3]=0.0;
    mlag[4]=0.0;                                    mlag[7]=0.0;
    mlag[8]=0.0;                                    mlag[11]=0.0;
    mlag[12]=0.0;   mlag[13]=0.0;   mlag[14]=0.0;   mlag[15]=1.0;
 
    mlag2[0]= 1.0;    mlag2[1]= 0.0;    mlag2[2]= 0.0;
    mlag2[4]= 0.0;    mlag2[5]= 1.0;    mlag2[6]= 0.0;    mlag2[7]= 0.0;
    mlag2[8]= 0.0;    mlag2[9]= 0.0;    mlag2[10]= 1.0;   mlag2[11]= 0.0;
    mlag2[12]= 0.0;   mlag2[13]= 0.0;   mlag2[14]= 0.0;   mlag2[15]= 1.0;

    if (dir==1)
    {
        mlag[5]=cos(biraketa);
        mlag[6]=-sin(biraketa);
        mlag[9]=sin(biraketa);
        mlag[10]=cos(biraketa);
        mlag2[3]= 3.3;
    }
    else
    {
        biraketa = 2*PI - biraketa;
        mlag[5]=cos(biraketa);
        mlag[6]=-sin(biraketa);
        mlag[9]=sin(biraketa);
        mlag[10]=cos(biraketa);
        mlag2[3]= -3.3;
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
        if(aldaketa=='r')
        {
            biderkatu_matrizeak(mleb->m, sel_ptr_k->mptr->m, mlag);
        }else{ //translazioa
            if (ald_lokala==1)
            {
                biderkatu_matrizeak(mleb->m, sel_ptr_k->mptr->m, mlag2);
            }else{
                biderkatu_matrizeak(mleb->m, mlag2, sel_ptr_k->mptr->m);
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
    mlist * mleb;

                    mlag[1]=0.0;                    mlag[3]=0.0;
    mlag[4]=0.0;    mlag[5]=1.0;    mlag[6]=0.0;    mlag[7]=0.0;
                    mlag[9]=0.0;                    mlag[11]=0.0;
    mlag[12]=0.0;   mlag[13]=0.0;   mlag[14]=0.0;   mlag[15]=1.0;

    mlag2[0]= 1.0;  mlag2[1]= 0.0;  mlag2[2]= 0.0;  mlag2[3]= 0.0;
    mlag2[4]= 0.0;  mlag2[5]= 1.0;  mlag2[6]= 0.0;
    mlag2[8]= 0.0;  mlag2[9]= 0.0;  mlag2[10]= 1.0; mlag2[11]= 0.0;
    mlag2[12]= 0.0; mlag2[13]= 0.0; mlag2[14]= 0.0; mlag2[15]= 1.0;
    if (dir==1)
    {
        mlag[0]=cos(biraketa);
        mlag[2]=sin(biraketa);
        mlag[8]=-sin(biraketa);
        mlag[10]=cos(biraketa);
        mlag2[7]= 3.3;
    }
    else
    {
        biraketa = 2*PI - biraketa;
        mlag[0]=cos(biraketa);
        mlag[2]=sin(biraketa);
        mlag[8]=-sin(biraketa);
        mlag[10]=cos(biraketa);
        mlag2[7]= -3.3;
    }
    mleb = (mlist *)malloc(sizeof(mlist));

    if(aldaketa=='r')
    {
        biderkatu_matrizeak(mleb->m, sel_ptr->mptr->m, mlag);
    }else{
        biderkatu_matrizeak(mleb->m, sel_ptr->mptr->m, mlag2);
    }  
    mleb->hptr = sel_ptr->mptr;
    sel_ptr->mptr = mleb;
    
}



void z_aldaketa(int dir)
{   
    double biraketa = 0.2618; //15 gradu
    double mlag[16];
    double mlag2[16];
    mlist * mleb;

                                    mlag[2]=0.0;    mlag[3]=0.0;
                                    mlag[6]=0.0;    mlag[7]=0.0;
    mlag[8]=0.0;    mlag[9]=0.0;    mlag[10]=1.0;   mlag[11]=0.0;
    mlag[12]=0.0;   mlag[13]=0.0;   mlag[14]=0.0;   mlag[15]=1.0;

    mlag2[0]= 1.0;    mlag2[1]= 0.0;    mlag2[2]= 0.0;    mlag2[3]= 0.0;
    mlag2[4]= 0.0;    mlag2[5]= 1.0;    mlag2[6]= 0.0;    mlag2[7]= 0.0;
    mlag2[8]= 0.0;    mlag2[9]= 0.0;    mlag2[10]= 1.0;
    mlag2[12]= 0.0;   mlag2[13]= 0.0;   mlag2[14]= 0.0;   mlag2[15]= 1.0;
    if (dir==1)
    {
        mlag[0]=cos(biraketa);
        mlag[1]=-sin(biraketa);
        mlag[4]=sin(biraketa);
        mlag[5]=cos(biraketa);

        mlag2[11]= 3.3;
    }
    else
    {
        biraketa = 2*PI - biraketa;
        mlag[0]=cos(biraketa);
        mlag[1]=-sin(biraketa);
        mlag[4]=sin(biraketa);
        mlag[5]=cos(biraketa);

        mlag2[11]= -3.3;
    }
    mleb = (mlist *)malloc(sizeof(mlist));
    
    
    if(aldaketa=='r')
    {
        biderkatu_matrizeak(mleb->m, sel_ptr->mptr->m, mlag);
    }else{
        biderkatu_matrizeak(mleb->m, sel_ptr->mptr->m, mlag2);
    }
    mleb->hptr = sel_ptr->mptr;
    sel_ptr->mptr = mleb;
    
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
    biderkatu_matrizeak(mleb->m, sel_ptr->mptr->m, mlag);
    sel_ptr->mptr = mleb;
}

void undo()
{
    double im[16];
    im[0]=1.0;    im[1]=0.0;    im[2]=0.0;    im[3]=0.0;
    im[4]=0.0;    im[5]=1.0;    im[6]=0.0;    im[7]=0.0;
    im[8]=0.0;    im[9]=0.0;    im[10]=1.0;   im[11]=0.0;
    im[12]=0.0;   im[13]=0.0;   im[14]=0.0;   im[15]=1.0;
    if(sel_ptr->mptr->hptr->m != im)
        sel_ptr->mptr = sel_ptr->mptr->hptr;
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
	case 'o':
		if (objektuak == 1) objektuak = 0;
		    else objektuak = 1;
		break;
	case 'l':
		if (lineak == 1) lineak = 0;
		    else lineak = 1;
		break;
	case 't':
	    aldaketa = 't';
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
	        read_from_file(fitxiz);
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
            if (foptr != 0) // objekturik gabe ez du ezer egin behar
                            // si no hay objeto no hace nada
                {
                sel_ptr = sel_ptr->hptr;
                /*The selection is circular, thus if we move out of the list we go back to the first element*/
                if (sel_ptr == 0) sel_ptr = foptr;
                indexx =0; // the selected polygon is the first one
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
        denak = 0;
        lineak =0;
        objektuak = 0;
        foptr = 0;
        sel_ptr = 0;
        sel_ptr_k=0;
        foptr_k=0;
        kam=0;
        aldaketa = 'r';
        ald_lokala = 1;
        if (argc>1) read_from_file(argv[1]);
            else read_from_file("adibideak.txt");
	glutMainLoop();

	return 0;   
}
