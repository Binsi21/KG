
typedef struct punto
{
float x, y, z, u,v;
} punto;

typedef struct hiruki
{
punto p1,p2,p3;
double normala[3]; //bektore normala ikusteko
} hiruki;

int cargar_triangulos(char *fitxiz, int *hkopptr, hiruki **hptrptr);

void kalkulatu_vnormala(hiruki * hirukia);

void ebaketa_kalkulatu(punto * gptr, punto * bptr, int h, punto *eptr);

void kalkulatu_mesa(double matriz[16], double mesa[16]);

int cargar_triangulos_color(char *fitxiz, int *hkopptr, hiruki **hptrptr, unsigned char **rgbptr);

void biderkatu_matrizeak(double emaitzaptr[16], double mlag[16], double eskptr[16]);

void objektuari_begiratu();

void perspektiba_aldatu(double *mperspektiba);
