//	Program developed by
//
//	Informatika Fakultatea
//	Euskal Herriko Unibertsitatea
//	http://www.ehu.eus/if
//
// to compile it: gcc dibujar-triangulos-y-objetos.c -lGL -lGLU -lglut

#include <GL/glut.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include "cargar-triangulo.h"

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

typedef struct my_camera
{
    double Mcsr[16];
    mlist *mptr;
} my_camera;

// testuraren informazioa
// información de textura

extern int load_ppm(char *file, unsigned char **bufferptr, int *dimxptr, int *dimyptr);
void print_matrizea(char *str, double to_print[16]);
void establecer_camara(double atx, double aty, double atz);
unsigned char *bufferra;
int dimx, dimy;
int indexx;
// Direcciones de los objetos que NO son cáramas.
triobj *foptr;
triobj *sel_ptr;

// Direcciones de los objetos que SÍ son cámaras.
my_camera *cam_ptr;

// Aquí se guardan los valores para las decisiones que se toman.
int denak;
int lineak;
int analisis;
int camara;
int persp;
int objektuak;
char aldaketa;
int ald_lokala;
unsigned char *colorv;
double Mp[16] = {0.0};
double at[3] = {0.0};
double Mmodelview[16] = {0.0};

char fitxiz[100];

void obtener_rotacion_rodrigues(double x, double y, double z, double angulo, double m[16])
{
    // Matriz de rotación rodrigues

    m[0] = cos(angulo) + (1 - cos(angulo)) * x * x;
    m[1] = (1 - cos(angulo)) * x * y - sin(angulo) * z;
    m[2] = (1 - cos(angulo)) * x * z + sin(angulo) * y;

    // Segunda fila
    m[4] = (1 - cos(angulo)) * x * y + sin(angulo) * z;
    m[5] = cos(angulo) + (1 - cos(angulo)) * y * y;
    m[6] = (1 - cos(angulo)) * y * z - sin(angulo) * x;

    // Tercera fila
    m[8] = (1 - cos(angulo)) * x * z - sin(angulo) * y;
    m[9] = (1 - cos(angulo)) * y * z + sin(angulo) * x;
    m[10] = cos(angulo) + (1 - cos(angulo)) * z * z;

    m[15] = 1.0;
    print_matrizea("rotacion rodrigues", m);
}
void mxm(double resultado[16], double operando_izquierdo[16], double operando_derecho[16])
{
    int i, j, k;
    double res;
    for (i = 0; i < 4; i++)
    {
        for (j = 0; j < 4; j++)
        {
            res = 0.0;
            for (k = 0; k < 4; k++)
            {
                res += operando_izquierdo[i * 4 + k] * operando_derecho[k * 4 + j];
            }
            // printf("res: %f\n", res);
            resultado[i * 4 + j] = res;
        }
    }
}

// TODO
// funtzio honek u eta v koordenatuei dagokien pointerra itzuli behar du.
// debe devolver el pointer correspondiente a las coordenadas u y v
unsigned char *color_textura(float u, float v)
{
    char *lag;

    int desplazamiento_u, desplazamiento_v;

    desplazamiento_u = trunc(u * dimx);
    desplazamiento_v = trunc((1 - v) * dimy);

    lag = (unsigned char *)bufferra; // pixel on the left and top
    return (lag + 3 * (desplazamiento_v * dimx + desplazamiento_u));
}

// TODO
// lerroa marrazten du, baina testuraren kodea egokitu behar da
// dibuja una linea pero hay que codificar la textura
void dibujar_linea_z(int linea, float c1x, float c1z, float c1u, float c1v, float c2x, float c2z, float c2u, float c2v)
{
    float xkoord, zkoord;
    float u, v;
    float difu, difv, difz;
    unsigned char r, g, b;

    glBegin(GL_POINTS);

    // Calculamos la pendiente de la recta que corta el triángulo, si c1x = c2x evitamos dividir por cero.
    if (c1x != c2x)
    {
        difv = (c2v - c1v) / (c2x - c1x);
        difu = (c2u - c1u) / (c2x - c1x);
        difz = (c2z - c1z) / (c2x - c1x);
    }
    else
    {
        difv = 0;
        difu = 0;
        difz = 0;
    }

    for (xkoord = c1x, zkoord = c1z, u = c1u, v = c1v; xkoord <= c2x; xkoord++)
    {
        // TODO
        // color_textura funtzioa ondo kodetu
        // programar de forma correcta la función color_textura

        colorv = color_textura(u, v);
        r = colorv[0];
        g = colorv[1];
        b = colorv[2];
        glColor3ub(r, g, b);
        glVertex3f(xkoord, linea, zkoord);
        // TODO
        // zkoord, u eta v berriak kalkulatu eskuineko puntuarentzat
        // calcular zkoord, u y v del siguiente pixel
        u += difu;
        v += difv;
        zkoord += difz;
    }
    glEnd();
}

void print_matrizea(char *str, double to_print[16])
{
    int i;
    printf("%s\n", str);
    for (i = 0; i < 4; i++)
        printf("%lf, %lf, %lf, %lf\n", to_print[i * 4], to_print[i * 4 + 1], to_print[i * 4 + 2],
               to_print[i * 4 + 3]);
}
void normalizar(double vector[3])
{
    double length = 0.0;
    int i;
    for (i = 0; i < 3; i++)
    {
        length += vector[i] * vector[i];
    }
    length = sqrt(length);
    for (i = 0; i < 3; i++)
    {
        vector[i] = vector[i] / length;
    }
}

// A partir de la matriz de la cánara, se calcula la nueva Mcsr. Hay que llamar a esta función cada vez que se realice una transformación en la cámara.

void calcular_mcsr()
{
    int i, j;

    for (i = 0; i < 16; i++)
    {
        cam_ptr->Mcsr[i] = 0.0;
    }

    cam_ptr->Mcsr[3] = -(cam_ptr->mptr->m[3] * cam_ptr->mptr->m[0] + cam_ptr->mptr->m[7] * cam_ptr->mptr->m[4] + cam_ptr->mptr->m[11] * cam_ptr->mptr->m[8]);
    cam_ptr->Mcsr[7] = -(cam_ptr->mptr->m[3] * cam_ptr->mptr->m[1] + cam_ptr->mptr->m[7] * cam_ptr->mptr->m[5] + cam_ptr->mptr->m[11] * cam_ptr->mptr->m[9]);
    cam_ptr->Mcsr[11] = -(cam_ptr->mptr->m[3] * cam_ptr->mptr->m[2] + cam_ptr->mptr->m[7] * cam_ptr->mptr->m[6] + cam_ptr->mptr->m[11] * cam_ptr->mptr->m[10]);
    // Poner los valores de Mc de fila en columna en Mcsr y la ultima fila de la Mcsr 0,0,0,1
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            cam_ptr->Mcsr[i * 4 + j] = cam_ptr->mptr->m[j * 4 + i];
        }
    }
    cam_ptr->Mcsr[15] = 1.0;
}

// TODO
// aurrerago egitekoa
// para más adelante
void mxp(punto *pptr, double m[16], punto p)
{
    // print_matrizea("");
    pptr->x = p.x * m[0] + p.y * m[1] + p.z * m[2] + m[3];
    pptr->y = p.x * m[4] + p.y * m[5] + p.z * m[6] + m[7];
    pptr->z = p.x * m[8] + p.y * m[9] + p.z * m[10] + m[11];
    pptr->w = p.x * m[12] + p.y * m[13] + p.z * m[14] + m[15];
    pptr->u = p.u;
    pptr->v = p.v;
}

void calcula_punto_corte(punto *punto_superior, punto *punto_inferior, float i, punto *corte)
{
    float m = 0, m2 = 0, m3 = 0, m4 = 0, delta_y;
    delta_y = punto_inferior->y - punto_superior->y;
    if (punto_inferior->x != punto_superior->x)
    {
        m = (delta_y) / (punto_inferior->x - punto_superior->x);
        corte->x = ((i - punto_superior->y) / m) + punto_superior->x;
    }
    else
    {
        corte->x = punto_superior->x;
    }
    m2 = (delta_y) / (punto_inferior->u - punto_superior->u);

    m3 = (delta_y) / (punto_inferior->v - punto_superior->v);
    m4 = (delta_y) / (punto_inferior->z - punto_superior->z);

    corte->u = ((i - punto_superior->y) / m2) + punto_superior->u;
    corte->v = ((i - punto_superior->y) / m3) + punto_superior->v;
    corte->z = ((i - punto_superior->y) / m4) + punto_superior->z;
}

// Función para encontrar los puntos máximo, intermedio y mínimo en tptr. Los guarda, respectivamente, en pgoiptr, perdiptr y pbeheptr.
void encontrar_max_min(punto **pgoiptr, punto **perdiptr, punto **pbeheptr)
{

    punto *aux_max;
    punto *aux_min;
    punto *aux_med;
    if ((*pgoiptr)->y > (*perdiptr)->y)
    {
        aux_max = *pgoiptr;
        aux_min = *perdiptr;
    }
    else
    {
        aux_max = *perdiptr;
        aux_min = *pgoiptr;
    }
    if ((*pbeheptr)->y > aux_max->y)
    {
        aux_med = aux_max;
        aux_max = *pbeheptr;
    }
    else if ((*pbeheptr)->y < aux_min->y)
    {
        aux_med = aux_min;
        aux_min = *pbeheptr;
    }
    else
    {
        aux_med = *pbeheptr;
    }
    *pgoiptr = aux_max;
    *perdiptr = aux_med;
    *pbeheptr = aux_min;
}

void rellenar_triangulo(punto *pgoiptr, punto *perdiptr, punto *pbeheptr)
{
    float i = 0;
    punto corte1;
    punto corte2;
    // Dibujamos la mitad superior del triángulo.

    for (i = pgoiptr->y; i > perdiptr->y; i--)
    {
        calcula_punto_corte(pgoiptr, pbeheptr, i, &corte1);
        calcula_punto_corte(pgoiptr, perdiptr, i, &corte2);

        if (corte1.x <= corte2.x)
            dibujar_linea_z(i, corte1.x, corte1.z, corte1.u, corte1.v, corte2.x, corte2.z, corte2.u, corte2.v);
        else
            dibujar_linea_z(i, corte2.x, corte2.z, corte2.u, corte2.v, corte1.x, corte1.z, corte1.u, corte1.v);
    }
    // Dibujamos la mitad inferior del triángulo.
    for (i = perdiptr->y; i > pbeheptr->y; i--)
    {
        calcula_punto_corte(perdiptr, pbeheptr, i, &corte1);
        calcula_punto_corte(pgoiptr, pbeheptr, i, &corte2);
        if (corte1.x <= corte2.x)
            dibujar_linea_z(i, corte1.x, corte1.z, corte1.u, corte1.v, corte2.x, corte2.z, corte2.u, corte2.v);
        else
            dibujar_linea_z(i, corte2.x, corte2.z, corte2.u, corte2.v, corte1.x, corte1.z, corte1.u, corte1.v);
    }
}

void mpxptr(punto *pt)
{
    mxp(pt, Mp, *pt);
    // printf("punto: %f, %f, %f, %f, %f, %f\n", pt->x, pt->y, pt->z, pt->u, pt->v, pt->w);
    if (pt->w == 0)
        pt->w = -1.0;
    pt->x = (pt->x * 500) / pt->w;
    pt->y = (pt->y * 500.0) / pt->w;
    pt->z = (-pt->z * 500.0) / pt->w;
    pt->w = 1.0;
}

void dibujar_triangulo(triobj *optr, int i)
{
    hiruki *tptr;
    punto *pgoiptr, *pbeheptr, *perdiptr;
    punto *pgoiptr2, *pbeheptr2, *perdiptr2;
    punto corte1, corte2;
    int start1, star2;
    double aux[16];
    float t = 1, s = 1, q = 1;
    float decremento_t = 1, decremento_s = 1, decremento_q = 1;
    float c1x, c1z, c1u, c1v, c2x, c2z, c2u, c2v;
    int linea;
    float cambio1, cambio1z, cambio1u, cambio1v, cambio2, cambio2z, cambio2u, cambio2v;
    punto p1, p2, p3;

    if (i >= optr->num_triangles)
        return;
    tptr = optr->triptr + i;
    // print_matrizea("Mcsr", cam_ptr->Mcsr);
    // print_matrizea("Moptr", optr->mptr->m);
    // printf("Posicion del objeto:\n");
    // printf("x: %f, y: %f, z: %f\n", optr->mptr->m[3], optr->mptr->m[7], optr->mptr->m[11]);
    // mxm(Modelview, sel_cptr->Mcsr, optr->mptr->m);
    // print_matrizea("Modelview", Mmodelview);

    // printf("p1: %f, %f, %f, %f, %f, %f\n", optr->triptr[i].p1.x, optr->triptr[i].p1.y, optr->triptr[i].p1.z, optr->triptr[i].p1.u, optr->triptr[i].p1.v, optr->triptr[i].p1.w);
    // printf("p2: %f, %f, %f, %f, %f, %f\n", optr->triptr[i].p2.x, optr->triptr[i].p2.y, optr->triptr[i].p2.z, optr->triptr[i].p2.u, optr->triptr[i].p2.v, optr->triptr[i].p2.w);
    // printf("p3: %f, %f, %f, %f, %f, %f\n", optr->triptr[i].p3.x, optr->triptr[i].p3.y, optr->triptr[i].p3.z, optr->triptr[i].p3.u, optr->triptr[i].p3.v, optr->triptr[i].p3.w);
    mxp(&p1, Mmodelview, tptr->p1);
    mxp(&p2, Mmodelview, tptr->p2);
    mxp(&p3, Mmodelview, tptr->p3);

    if (persp)
    {
        //     printf("p1: %f, %f, %f, %f, %f, %f\n", p1.x, p1.y, p1.z, p1.u, p1.v, p1.w);
        //     printf("p2: %f, %f, %f, %f, %f, %f\n", p2.x, p2.y, p2.z, p2.u, p2.v, p2.w);
        //     printf("p3: %f, %f, %f, %f, %f, %f\n", p3.x, p3.y, p3.z, p3.u, p3.v, p3.w);
        mpxptr(&p1);
        mpxptr(&p2);
        mpxptr(&p3);
    }
    // printf("p1: %f, %f, %f, %f, %f, %f\n", p1.x, p1.y, p1.z, p1.u, p1.v, p1.w);
    // printf("p2: %f, %f, %f, %f, %f, %f\n", p2.x, p2.y, p2.z, p2.u, p2.v, p2.w);
    // printf("p3: %f, %f, %f, %f, %f, %f\n", p3.x, p3.y, p3.z, p3.u, p3.v, p3.w);

    if (lineak == 1)
    {
        glBegin(GL_POLYGON);
        glVertex3d(p1.x, p1.y, p1.z);
        glVertex3d(p2.x, p2.y, p2.z);
        glVertex3d(p3.x, p3.y, p3.z);
        glEnd();
        return;
    }

    // Encontramos el punto máximo, mínimo, y medio del triángulo.

    pgoiptr = &p1;
    perdiptr = &p2;
    pbeheptr = &p3;

    encontrar_max_min(&pgoiptr, &perdiptr, &pbeheptr);
    // printf("pgoi: %f, %f, %f, %f, %f\n", pgoiptr->x, pgoiptr->y, pgoiptr->z, pgoiptr->u, pgoiptr->v);
    // printf("perdi: %f, %f, %f, %f, %f\n", perdiptr->x, perdiptr->y, perdiptr->z, perdiptr->u, perdiptr->v);
    // printf("pbehe: %f, %f, %f, %f, %f\n", pbeheptr->x, pbeheptr->y, pbeheptr->z, pbeheptr->u, pbeheptr->v);

    // Llamada a función rellenar triángulo.

    rellenar_triangulo(pgoiptr, perdiptr, pbeheptr);
}

static void marraztu(void)
{
    int i, j;
    triobj *auxptr;
    triobj *cauxptr;
    double aux[16];

    /*
    unsigned char* colorv;
    unsigned char r,g,b;
    */

    // marrazteko objektuak behar dira
    // no se puede dibujar sin objetos
    if (foptr == 0)
        return;
    // clear viewport...
    if (objektuak == 1)
    {
        glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
        GLenum error = glGetError();
        if (error != GL_NO_ERROR)
        {
            fprintf(stderr, "OpenGL error after glClear: %s\n", gluErrorString(error));
        }
    }
    else
    {
        if (denak == 0)
            glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
    }
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    if (persp == 1)
    {
        glOrtho(-500.0, 500.0, -500.0, 500.0, -5, 500.0);
    }
    else
    {
        glOrtho(-500.0, 500.0, -500.0, 500.0, -500.0, 500.0);
    }

    for (i = 0; i < sel_ptr->num_triangles; i++)
    {
        // printf("p1: %f, %f, %f, %f, %f, %f\n", sel_ptr->triptr[i].p1.x, sel_ptr->triptr[i].p1.y, sel_ptr->triptr[i].p1.z, sel_ptr->triptr[i].p1.u, sel_ptr->triptr[i].p1.v, sel_ptr->triptr[i].p1.w);
        // printf("p2: %f, %f, %f, %f, %f, %f\n", sel_ptr->triptr[i].p2.x, sel_ptr->triptr[i].p2.y, sel_ptr->triptr[i].p2.z, sel_ptr->triptr[i].p2.u, sel_ptr->triptr[i].p2.v, sel_ptr->triptr[i].p2.w);
        // printf("p3: %f, %f, %f, %f, %f, %f\n", sel_ptr->triptr[i].p3.x, sel_ptr->triptr[i].p3.y, sel_ptr->triptr[i].p3.z, sel_ptr->triptr[i].p3.u, sel_ptr->triptr[i].p3.v, sel_ptr->triptr[i].p3.w);
    }
    if (objektuak == 1)
    {
        if (denak == 1)
        {

            for (auxptr = foptr; auxptr != 0; auxptr = auxptr->hptr)
            {
                mxm(Mmodelview, cam_ptr->Mcsr, auxptr->mptr->m);

                for (i = 0; i < auxptr->num_triangles; i++)
                {
                    dibujar_triangulo(auxptr, i);
                }
            }
        }
        else
        {
            for (i = 0; i < sel_ptr->num_triangles; i++)
            {
                dibujar_triangulo(sel_ptr, i);
            }
        }
    }
    else
    {
        dibujar_triangulo(sel_ptr, indexx);
    }
    glFlush();
}

void read_from_file(char *fitx)
{
    int i, retval;
    triobj *optr;
    my_camera *cam_ptr;

    // printf("%s fitxategitik datuak hartzera\n",fitx);
    optr = (triobj *)malloc(sizeof(triobj));
    cam_ptr = (my_camera *)malloc(sizeof(my_camera));
    retval = cargar_triangulos_color(fitx, &(optr->num_triangles), &(optr->triptr), &(colorv));
    if (retval != 9 & retval != 15)
    {
        printf("%d\n", retval);
        printf("%s fitxategitik datuak hartzerakoan arazoak izan ditut\n    Problemas al leer\n", fitxiz);
        free(optr);
    }
    else
    {
        // for (i = 0; i < optr->num_triangles; i++)
        // {
        //     printf("p1: %f, %f, %f, %f, %f, %f\n", optr->triptr[i].p1.x, optr->triptr[i].p1.y, optr->triptr[i].p1.z, optr->triptr[i].p1.u, optr->triptr[i].p1.v, optr->triptr[i].p1.w);
        //     printf("p2: %f, %f, %f, %f, %f, %f\n", optr->triptr[i].p2.x, optr->triptr[i].p2.y, optr->triptr[i].p2.z, optr->triptr[i].p2.u, optr->triptr[i].p2.v, optr->triptr[i].p2.w);
        //     printf("p3: %f, %f, %f, %f, %f, %f\n", optr->triptr[i].p3.x, optr->triptr[i].p3.y, optr->triptr[i].p3.z, optr->triptr[i].p3.u, optr->triptr[i].p3.v, optr->triptr[i].p3.w);
        // }
        // printf("objektuaren matrizea...\n");
        optr->mptr = (mlist *)malloc(sizeof(mlist));

        for (i = 0; i < 16; i++)
            optr->mptr->m[i] = 0;
        optr->mptr->m[0] = 1.0;
        optr->mptr->m[5] = 1.0;
        optr->mptr->m[10] = 1.0;
        optr->mptr->m[15] = 1.0;
        optr->mptr->hptr = 0;
        // printf("objektu zerrendara doa informazioa...\n");
        optr->hptr = foptr;
        foptr = optr;
        sel_ptr = optr;
    }
    printf("datuak irakurrita\nLectura finalizada\n");
}

void transformacion_principal(double m[16])
{
    // Modo analisis a la derecha, modo vuelo a la izquierda
    double resultado[16];
    double aux[16];
    double Mat[16];
    int i;
    int j;

    print_matrizea("Me ha llegado esta bazofia:", m);
    for (i = 0; i < 16; i++)
    {
        Mat[i] = 0;
        resultado[i] = 0;
        aux[i] = 0;
    }
    mlist *new_m = (mlist *)malloc(sizeof(mlist));
    // printf("aldaera lokala: %d\n", ald_lokala);
    if (ald_lokala == 1)
    {
        if (camara == 1)
        {
            mxm(new_m->m, cam_ptr->mptr->m, m);
            new_m->hptr = cam_ptr->mptr;
            cam_ptr->mptr = new_m;
            print_matrizea("Matriz de la camara:", cam_ptr->mptr->m);
            calcular_mcsr();
            print_matrizea("Mcsr:", cam_ptr->Mcsr);
            return;
        }
        // print_matrizea("Objeto a multiplicar por la izquierda", sel_ptr->mptr->m);
        mxm(new_m->m, sel_ptr->mptr->m, m);
        // print_matrizea("Objeto a multiplicar por la derecha", m);
    }
    else
    {
        if (camara == 1)
        {
            mxm(new_m->m, m, cam_ptr->mptr->m);
            new_m->hptr = cam_ptr->mptr;
            cam_ptr->mptr = new_m;
            calcular_mcsr();
            return;
        }
        mxm(new_m->m, m, sel_ptr->mptr->m);
    }
    new_m->hptr = sel_ptr->mptr;
    sel_ptr->mptr = new_m;
}
void x_aldaketa(int dir)
{
    /**
     * Creamos la matriz con la que vamos a multiplicar para lograr la rotación,
     * y llamamos a la función transformación principal.
     */
    double m[16];
    double Mat[16];
    double aux[16];
    int exponent = 0;
    int i = 0;
    double angulo = 0;
    double x = 0;
    double y = 0;
    double z = 0;

    for (i = 0; i < 16; i++)
    {
        m[i] = 0;
        Mat[i] = 0;
        aux[i] = 0;
    }

    m[0] = 1;
    m[5] = 1;
    m[10] = 1;
    m[15] = 1;

    angulo = 0.075 * dir;

    if (aldaketa == 't')
    {
        m[3] = dir * 5;
        print_matrizea("Traslacion: ", m);
        transformacion_principal(m);
        return;
    }

    if (camara == 1 & ald_lokala == 0)
    {
        x = cam_ptr->mptr->m[0];
        y = cam_ptr->mptr->m[4];
        z = cam_ptr->mptr->m[8];

        Mat[0] = 1;
        Mat[5] = 1;
        Mat[10] = 1;
        Mat[15] = 1;

        Mat[3] = -sel_ptr->mptr->m[3];
        Mat[7] = -sel_ptr->mptr->m[7];
        Mat[11] = -sel_ptr->mptr->m[11];
        obtener_rotacion_rodrigues(x, y, z, angulo, m);
        mxm(aux, m, Mat);
        Mat[3] = -Mat[3];
        Mat[7] = -Mat[7];
        Mat[11] = -Mat[11];
        mxm(m, Mat, aux);
        // transformacion_principal(m);
        // return;
    }
    else
    {
        m[5] = cos(angulo);
        m[6] = -sin(angulo);
        m[9] = sin(angulo);
        m[10] = cos(angulo);
    }

    transformacion_principal(m);
}

void y_aldaketa(int dir)
{
    double m[16];
    double aux[16];
    double Mat[16];
    int exponent = 0;
    int i = 0;
    double x = 0;
    double y = 0;
    double z = 0;
    double angulo = 0;

    for (i = 0; i < 16; i++)
    {
        m[i] = 0;
        aux[i] = 0;
        Mat[i] = 0;
    }

    m[0] = 1;
    m[5] = 1;
    m[10] = 1;
    m[15] = 1;

    if (aldaketa == 't')
    {
        m[7] = dir * 5;
        transformacion_principal(m);
        return;
    }

    angulo = dir * 0.075;

    if (camara == 1 & ald_lokala == 0)
    {

        // getchar();

        x = cam_ptr->mptr->m[1];
        y = cam_ptr->mptr->m[5];
        z = cam_ptr->mptr->m[9];

        Mat[3] = -sel_ptr->mptr->m[3];
        Mat[7] = -sel_ptr->mptr->m[7];
        Mat[11] = -sel_ptr->mptr->m[11];
        obtener_rotacion_rodrigues(x, y, z, angulo, m);
        print_matrizea("rotacion rodrigues", m);
        mxm(aux, m, Mat);
        Mat[3] = -Mat[3];
        Mat[7] = -Mat[7];
        Mat[11] = -Mat[11];
        mxm(m, Mat, aux);
        print_matrizea("rotacion final en modo analisi", m);
    }
    else
    {
        m[0] = cos(angulo);
        m[2] = sin(angulo);
        m[8] = -sin(angulo);
        m[10] = cos(angulo);
    }
    transformacion_principal(m);
}
void z_aldaketa(int dir)
{

    double m[16];
    double Mat[16];
    double aux[16];
    int exponent = 0;
    int i = 0;
    double angulo;
    double x, y, z;
    for (i = 0; i < 16; i++)
    {
        m[i] = 0;
    }

    m[0] = 1;
    m[5] = 1;
    m[10] = 1;
    m[15] = 1;

    if (aldaketa == 't')
    {
        m[11] = dir * 5;
        transformacion_principal(m);
        return;
    }

    angulo = dir * 0.075;
    if (camara == 1 & ald_lokala == 0)
    {
        x = cam_ptr->mptr->m[2];
        y = cam_ptr->mptr->m[6];
        z = cam_ptr->mptr->m[10];

        Mat[0] = 1;
        Mat[5] = 1;
        Mat[10] = 1;
        Mat[15] = 1;

        Mat[3] = -sel_ptr->mptr->m[3];
        Mat[7] = -sel_ptr->mptr->m[7];
        Mat[11] = -sel_ptr->mptr->m[11];
        obtener_rotacion_rodrigues(x, y, z, angulo, m);

        mxm(aux, m, Mat);
        Mat[3] = -Mat[3];
        Mat[7] = -Mat[7];
        Mat[11] = -Mat[11];
        mxm(m, Mat, aux);
        print_matrizea("rotacion final en modo analisi", m);

        obtener_rotacion_rodrigues(x, y, z, angulo, m);
    }
    else
    {
        m[0] = cos(angulo);
        m[1] = -sin(angulo);
        m[4] = sin(angulo);
        m[5] = cos(angulo);
    }
    print_matrizea("Aldaketa", m);
    transformacion_principal(m);
}

void s_aldaketa(int dir)
{
    double m[16];
    int exponent = 0;
    int aldaketa = 0;
    int i = 0;
    for (i = 0; i < 16; i++)
    {
        m[i] = 0;
    }

    if (dir == 0)
    {
        m[0] = 0.950;
        m[5] = 0.950;
        m[10] = 0.950;
        m[15] = 1;
        transformacion_principal(m);
        return;
    }
    m[0] = 1.05;
    m[5] = 1.05;
    m[10] = 1.05;
    m[15] = 1;
    transformacion_principal(m);
}

void undo()
{
    if (camara == 1)
    {
        if (cam_ptr->mptr->hptr != NULL)
        {
            cam_ptr->mptr = cam_ptr->mptr->hptr;
        }
        calcular_mcsr(cam_ptr);
        return;
    }
    if (sel_ptr->mptr->hptr != NULL)
    {
        sel_ptr->mptr = sel_ptr->mptr->hptr;
    }
}

// This function will be called whenever the user pushes one key
static void teklatua(unsigned char key, int x, int y)
{
    int retval;
    int i;
    FILE *obj_file;

    switch (key)
    {
    case 13:
        if (foptr != 0) // objekturik ez badago ezer ez du egin behar
                        // si no hay objeto que no haga nada
        {
            indexx++; // azkena bada lehenengoa bihurtu
                      // pero si es el último? hay que controlarlo!
            if (indexx == sel_ptr->num_triangles)
            {
                indexx = 0;
                if ((denak == 1) && (objektuak == 0))
                {
                    glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
                    glFlush();
                }
            }
        }
        break;
    case 'd':
        if (denak == 1)
            denak = 0;
        else if (denak == 0)
            denak = 1;
        break;
    case 'c':
        if (camara == 1)
            camara = 0;
        else
        {
            if (ald_lokala == 0)
                ald_lokala = 1;
            camara = 1;
        }
        break;
    case 'o':
        if (objektuak == 1)
            objektuak = 0;
        else
            objektuak = 1;
        break;
    case 'l':
        if (lineak == 1)
            lineak = 0;
        else
            lineak = 1;
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
            establecer_camara(sel_ptr->mptr->m[3], sel_ptr->mptr->m[7], sel_ptr->mptr->m[11]);
            print_matrizea("Matriz de la camara mirando hacia el objeto seleccionado:", cam_ptr->mptr->m);
            print_matrizea("Mcsr: ", cam_ptr->Mcsr);
        }
        else
        {
            ald_lokala = 1;
        }
        break;
    case 's':
        s_aldaketa(0);
        break;
    case 'S':
        s_aldaketa(1);
        break;
    case 'x':
        x_aldaketa(1);
        break;
    case 'y':
        y_aldaketa(1);
        break;
    case 'z':
        printf("z\n");
        z_aldaketa(-1);
        break;
    case 'X':
        x_aldaketa(-1);
        break;
    case 'Y':
        y_aldaketa(-1);
        break;
    case 'Z':
        z_aldaketa(1);
        break;
    case 'u':
        undo();
        break;
    case 'p':
        if (persp == 1)
        {
            persp = 0;
        }
        else
        {
            persp = 1;
        }
        break;
    case 'P':
        if (persp == 0)
        {
            persp = 1;
        }
        else
        {
            persp = 0;
        }
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
    case 9:             /* <TAB> */
        if (foptr != 0) // objekturik gabe ez du ezer egin behar
                        // si no hay objeto no hace nada
        {
            sel_ptr = sel_ptr->hptr;
            /*The selection is circular, thus if we move out of the list we go back to the first element*/
            if (sel_ptr == 0)
                sel_ptr = foptr;
            indexx = 0; // the selected polygon is the first one
        }
        break;
    case 27: // <ESC>
        exit(0);
        break;
    default:
        printf("%d %c\n", key, key);
    }

    // The screen must be drawn to show the new triangle
    glutPostRedisplay();
}

void perspectiva()
{

    double l, r, b, t, n, f;

    l = -5.0;
    r = 5.0;
    b = -5.0;
    t = 5.0;
    n = 5.0;
    f = 50000.0;

    Mp[0] = (2 * n) / (r - l);
    Mp[2] = (r + l) / (r - l);
    Mp[5] = (2 * n) / (t - b);
    Mp[6] = (t + b) / (t - b);
    Mp[10] = -(f + n) / (f - n);
    Mp[11] = -(2 * f * n) / (f - n);
    Mp[14] = -1;

    printf("Perspectiva:\n");
}

// Establece la cámara dando un punto de atención y recalcula la matriz de cambio de sistema de referencia.

void establecer_camara(double atx, double aty, double atz)
{
    int i, j;
    double Z_cam[3];
    double Pos_cam[3];
    double X_cam[3];
    double Y_cam[3];
    double up[3] = {0.0, 1.0, 0.0};

    // Establecemos todo a cero, lo preparamos todo para los nuevos cambios.

    for (i = 0, j = 0; i < 16; i++)
    {
        cam_ptr->mptr->m[i] = 0;
        cam_ptr->Mcsr[i] = 0;
        if (j < 3)
        {
            Z_cam[j] = 0;
            Pos_cam[j] = 0;
            X_cam[j] = 0;
            Y_cam[j] = 0;
        }
    }

    cam_ptr->mptr->m[15] = 1;
    cam_ptr->Mcsr[15] = 1;

    // Empezamos por decidir la posición de la cámara.

    cam_ptr->mptr->m[3] = 0;
    cam_ptr->mptr->m[7] = 0;
    cam_ptr->mptr->m[11] = 200; // Para poder verlo todo desde la distancia.

    // Calculamos la columna z de la matriz de la cámara.

    Z_cam[0] = cam_ptr->mptr->m[3] - atx;
    Z_cam[1] = cam_ptr->mptr->m[7] - aty;
    Z_cam[2] = cam_ptr->mptr->m[11] - atz;

    normalizar(Z_cam);

    cam_ptr->mptr->m[2] = Z_cam[0];
    cam_ptr->mptr->m[6] = Z_cam[1];
    cam_ptr->mptr->m[10] = Z_cam[2];

    // Calculamos la columna x de la matriz de la cámara.

    X_cam[0] = up[1] * Z_cam[2] - up[2] * Z_cam[1];
    X_cam[1] = -(up[0] * Z_cam[2] - up[2] * Z_cam[0]);
    X_cam[2] = up[0] * Z_cam[1] - up[1] * Z_cam[0];

    normalizar(X_cam);

    cam_ptr->mptr->m[0] = X_cam[0];
    cam_ptr->mptr->m[4] = X_cam[1];
    cam_ptr->mptr->m[8] = X_cam[2];

    // Calculamos la columna y de la matriz de la cámara.

    Y_cam[0] = Z_cam[1] * X_cam[2] - Z_cam[2] * X_cam[1];
    Y_cam[1] = -(Z_cam[0] * X_cam[2] - Z_cam[2] * X_cam[0]);
    Y_cam[2] = Z_cam[0] * X_cam[1] - Z_cam[1] * X_cam[0];

    cam_ptr->mptr->m[1] = Y_cam[0];
    cam_ptr->mptr->m[5] = Y_cam[1];
    cam_ptr->mptr->m[9] = Y_cam[2];

    // Mco[3] *= -1;
    // Mco[7] *= -1;
    // Mco[11] *= -1;

    calcular_mcsr(cam_ptr);
}

int main(int argc, char **argv)
{
    int retval;
    double at[3] = {0.0, 0.0, 0.0};
    printf(" Triangeluak: barneko puntuak eta testura\n Triángulos con puntos internos y textura \n");
    printf("Press <ESC> to finish\n");
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(1000, 1000);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("Practica GCC");
    glutDisplayFunc(marraztu);
    glutKeyboardFunc(teklatua);
    /* we put the information of the texture in the buffer pointed by bufferra. The dimensions of the texture are loaded into dimx and dimy */
    retval = load_ppm("testura.ppm", &bufferra, &dimx, &dimy);
    if (retval != 1)
    {
        printf("Ez dago texturaren fitxategia (testura.ppm)\n");
        exit(-1);
    }
    glClearColor(0.0, 0.0, 0.0, 0.0);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glEnable(GL_DEPTH_TEST); // activar el test de profundidad (Z-buffer)
    denak = 1;
    lineak = 1;
    persp = 1;
    objektuak = 1;
    camara = 1;
    foptr = 0;
    sel_ptr = 0;
    aldaketa = 'r';
    ald_lokala = 1;
    perspectiva();
    printf("Preparando la camara...\n");
    cam_ptr = (my_camera *)malloc(sizeof(my_camera));
    cam_ptr->mptr = (mlist *)malloc(sizeof(mlist));
    establecer_camara(0, 0, 0);
    calcular_mcsr();
    print_matrizea("Camara estado inicial:", cam_ptr->mptr->m);
    read_from_file("k.txt");
    sel_ptr->mptr->m[3] = -200;
    sel_ptr->mptr->m[11] = 200;
    if (argc > 1)
    {
        read_from_file(argv[1]);
    }
    else
    {
        read_from_file("z.txt");
    }

    sel_ptr->mptr->m[3] = 200;
    sel_ptr->mptr->m[11] = 200;
    glutMainLoop();
    return 0;
}
