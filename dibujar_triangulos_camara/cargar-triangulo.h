
typedef struct punto
{
float x, y, z, u,v;
} punto;

typedef struct hiruki
{
punto p1,p2,p3;
} hiruki;

int cargar_triangulos_color(char *fitxiz, int *hkopptr, hiruki **hptrptr, unsigned char **rgbptr);
int cargar_triangulos(char *fitxiz, int *hkopptr, hiruki **hptrptr);

