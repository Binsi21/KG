
typedef struct punto
{
float x, y, z, u,v;
} punto;

typedef struct hiruki
{
punto p1,p2,p3;
} hiruki;

int cargar_triangulos(char *fitxiz, int *hkopptr, hiruki **hptrptr);

void ebaketa_kalkulatu(punto * gptr, punto * bptr, int h, punto *eptr);

void kalkulatu_mesa(double matriz[16], double mesa[16]);

int cargar_triangulos_color(char *fitxiz, int *hkopptr, hiruki **hptrptr, unsigned char **rgbptr);
