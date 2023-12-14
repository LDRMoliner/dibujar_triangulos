//	Program developed by
//
//	Informatika Fakultatea
//	Euskal Herriko Unibertsitatea
//	http://www.ehu.eus/if
//
// to compile it: gcc dibujar-triangulos-y-objetos.c -lGL -lGLU -lglut
//
//
//

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
    struct triobj *coptr;
    double Mcsr[16];
    struct my_camera *hptr;
} my_camera;

// testuraren informazioa
// información de textura

extern int load_ppm(char *file, unsigned char **bufferptr, int *dimxptr, int *dimyptr);
unsigned char *bufferra;
int dimx, dimy;
int indexx;
// Direcciones de los objetos que NO son cáramas.
triobj *foptr;
triobj *sel_ptr;

// Direcciones de los objetos que SÍ son cámaras.
my_camera *focptr;
my_camera *sel_cptr;

// Aquí se guardan los valores para las decisiones que se toman.
int denak;
int lineak;
int cam;
int persp;
int objektuak;
char aldaketa;
int ald_lokala;
unsigned char *colorv;
double Mp[16] = {0.0};
double Mmodelview[16] = {0.0};

char fitxiz[100];

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

void calcular_mcsr(my_camera *cam)
{
    int i, j;
    cam->Mcsr[3] = -(cam->coptr->mptr->m[0] * cam->coptr->mptr->m[3] + cam->coptr->mptr->m[4] * cam->coptr->mptr->m[7] + cam->coptr->mptr->m[8] * cam->coptr->mptr->m[11]);
    cam->Mcsr[7] = -(cam->coptr->mptr->m[1] * cam->coptr->mptr->m[3] + cam->coptr->mptr->m[5] * cam->coptr->mptr->m[7] + cam->coptr->mptr->m[9] * cam->coptr->mptr->m[11]);
    cam->Mcsr[11] = -(cam->coptr->mptr->m[2] * cam->coptr->mptr->m[3] + cam->coptr->mptr->m[6] * cam->coptr->mptr->m[7] + cam->coptr->mptr->m[10] * cam->coptr->mptr->m[11]);
    
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            cam->Mcsr[i * 4 + j] = cam->coptr->mptr->m[i * 4 + j];   
        }
    }
}
// TODO
// aurrerago egitekoa
// para más adelante
void mxp(punto *pptr, double m[16], punto p)
{
    // print_matrizea("");
    pptr->x = p.x * m[0]  + p.y * m[1]  +  p.z * m[2] + m[3];
    pptr->y = p.x * m[4] + p.y * m[5] + p.z * m[6] + m[7];
    pptr->z = p.x * m[8] + p.y * m[9] + p.z * m[10] + m[11];
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

void dibujar_triangulo(triobj *optr, int i)
{
    hiruki *tptr;
    punto *pgoiptr, *pbeheptr, *perdiptr;
    punto *pgoiptr2, *pbeheptr2, *perdiptr2;
    punto corte1, corte2;
    double Modelview[16];
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

    calcular_mcsr(sel_cptr);
    print_matrizea("Mcsr", sel_cptr->Mcsr);
    mxm(Modelview, sel_cptr->Mcsr, optr->mptr->m);
    print_matrizea("Modelview", Modelview);
    printf("p1: %f, %f, %f, %f, %f\n", tptr->p1.x, tptr->p1.y, tptr->p1.z, tptr->p1.u, tptr->p1.v);
    printf("p2: %f, %f, %f, %f, %f\n", tptr->p2.x, tptr->p2.y, tptr->p2.z, tptr->p2.u, tptr->p2.v);
    printf("p3: %f, %f, %f, %f, %f\n", tptr->p3.x, tptr->p3.y, tptr->p3.z, tptr->p3.u, tptr->p3.v);

    // mxm(Mmodelview, my_cam.Mcsr, optr->mptr->m);

    // print_matrizea("Mmodelview", Mmodelview);
    mxp(&p1, Modelview, tptr->p1);
    mxp(&p2, Modelview, tptr->p2);
    mxp(&p3, Modelview, tptr->p3);

    if (persp)
    {
        print_matrizea("Mp", Mp);
        mxp(&p1, Mp, p1);
        mxp(&p2, Mp, p2);
        mxp(&p3, Mp, p3);
        p1.x = p1.x / p1.w * 500.0;
        // p1.x = p1.x/p1.w;
        p1.y = p1.y / p1.w * 500.0;
        // p1.y = p1.y/p1.w;
        p1.z = -p1.z / p1.w;
        p1.w = 1.0;
        p2.x = p2.x / p2.w * 500.0;
        // p2.x = p2.x/p2.w;
        p2.y = p2.y / p2.w * 500.0;
        // p2.y = p2.y/p2.w;
        p2.z = -p2.z / p2.w;
        // p2.z = -p2.z/p2.w;
        p2.w = 1.0;
        p3.x = p3.x / p3.w * 500.0;
        // p3.x = p3.x/p3.w;
        p3.y = p3.y / p3.w * 500.0;
        // p3.y = p3.y/p3.w;
        p3.z = -p3.z / p3.w;
        p3.w = 1.0;
    }
    printf("p1 actual: %f, %f, %f, %f, %f\n", p1.x, p1.y, p1.z, p1.u, p1.v);
    printf("p2 actual: %f, %f, %f, %f, %f\n", p2.x, p2.y, p2.z, p2.u, p2.v);
    printf("p3 actual: %f, %f, %f, %f, %f\n", p3.x, p3.y, p3.z, p3.u, p3.v);
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
    glOrtho(-500.0, 500.0, -500.0, 500.0, 0.0, 500.0);

    if (objektuak == 1)
    {
        if (denak == 1)
        {

            for (auxptr = foptr; auxptr != 0; auxptr = auxptr->hptr)
            {
                for (i = 0; i < auxptr->num_triangles; i++)
                {
                    dibujar_triangulo(auxptr, i);
                }
            }
            printf("%d", focptr->coptr == 0);
            for (auxptr = focptr->coptr; auxptr != 0; auxptr = auxptr->hptr)
            {
                for (i = 0; i < auxptr->num_triangles; i++)
                {
                    printf("A dibujar\n");
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
            printf("Agur\n");
            for (i = 0; i < sel_cptr->coptr->num_triangles; i++)
            {
                dibujar_triangulo(sel_cptr->coptr, i);
            }
        }
    }
    else
    {
        dibujar_triangulo(sel_ptr, indexx);
        dibujar_triangulo(sel_cptr->coptr, indexx);
    }
    glFlush();
}

void read_from_file(char *fitx, int is_camera)
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
        if (is_camera)
        {
            cam_ptr->coptr = optr;
            cam_ptr->hptr = focptr;
            focptr = cam_ptr;
            sel_cptr = cam_ptr;
            for (i = 0; i < 16; i++)
                sel_cptr->Mcsr[i] = 0.0;
            sel_cptr->Mcsr[0] = 1.0;
            sel_cptr->Mcsr[5] = 1.0;
            sel_cptr->Mcsr[10] = 1.0;
            sel_cptr->Mcsr[15] = 200.0;

            print_matrizea("Mcsr", focptr->Mcsr);
        }
        else
        {
            optr->hptr = foptr;
            foptr = optr;
            sel_ptr = optr;
        }
    }
    printf("datuak irakurrita\nLectura finalizada\n");
}

void transformacion_principal(double m[16])
{
    double resultado[16];
    mlist *new_m = (mlist *)malloc(sizeof(mlist));
    int i = 0;

    if (ald_lokala == 1)
    {
        mxm(resultado, sel_ptr->mptr->m, m);
    }
    else
    {
        mxm(resultado, m, sel_ptr->mptr->m);
    }

    for (i = 0; i < 16; i++)
    {
        new_m->m[i] = resultado[i];
    }
    new_m->hptr = sel_ptr->mptr;
    sel_ptr->mptr = new_m;
    print_matrizea("", sel_ptr->mptr->m);
}
void x_aldaketa(int dir)
{
    /**
     * Creamos la matriz con la que vamos a multiplicar para lograr la rotación,
     * y llamamos a la función transformación principal.
     */
    double m[16];
    int exponent = 0;
    int i = 0;
    for (i = 0; i < 16; i++)
    {
        m[i] = 0;
    }

    m[0] = 1;
    m[15] = 1;

    exponent = pow(-1, !dir);

    if (aldaketa == 't')
    {
        m[10] = 1;
        m[5] = 1;
        m[3] = pow(-1, !dir) * 5;
        transformacion_principal(m);
        return;
    }

    m[5] = cos(exponent * 0.075);
    m[6] = -sin(exponent * 0.075);
    m[9] = sin(exponent * 0.075);
    m[10] = cos(exponent * 0.075);
    transformacion_principal(m);
}

void y_aldaketa(int dir)
{
    double m[16];
    int exponent = 0;
    int i = 0;
    for (i = 0; i < 16; i++)
    {
        m[i] = 0;
    }

    m[5] = 1;
    m[15] = 1;

    exponent = pow(-1, !dir);

    if (aldaketa == 't')
    {
        m[10] = 1;
        m[0] = 1;
        m[7] = pow(-1, !dir) * 5;
        transformacion_principal(m);
        return;
    }

    m[0] = cos(exponent * 0.075);
    m[2] = sin(exponent * 0.075);
    m[8] = -sin(exponent * 0.075);
    m[10] = cos(exponent * 0.075);
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

void z_aldaketa(int dir)
{

    double m[16];
    int exponent = 0;
    int i = 0;
    for (i = 0; i < 16; i++)
    {
        m[i] = 0;
    }

    m[0] = 1;
    m[15] = 1;

    exponent = pow(-1, !dir);

    if (aldaketa == 't')
    {
        return;
    }
    m[0] = cos(exponent * 0.075);
    m[1] = -sin(exponent * 0.075);
    m[4] = sin(exponent * 0.075);
    m[5] = cos(exponent * 0.075);
    m[10] = 1;
    print_matrizea("Aldaketa", m);
    transformacion_principal(m);
}

void undo()
{
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
        else
            denak = 1;
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
            ald_lokala = 0;
        else
            ald_lokala = 1;
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
        read_from_file(fitxiz, 0);
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
    n = 1.0;
    f = 1000.0;

    Mp[0] = (2 * n) / (r - l);
    Mp[2] = (r + l) / (r - l);
    Mp[5] = (2 * n) / (t - b);
    Mp[6] = (t + b) / (t - b);
    Mp[10] = -(f + n) / (f - n);
    Mp[11] = -(2 * f * n) / (f - n);
    Mp[14] = -1;

    printf("Perspectiva:\n");
    print_matrizea(" ", Mp);
}

int main(int argc, char **argv)
{
    int retval;

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

    glClearColor(0.0f, 0.0f, 0.7f, 1.0f);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glEnable(GL_DEPTH_TEST); // activar el test de profundidad (Z-buffer)
    denak = 1;
    lineak = 1;
    persp = 0;
    objektuak = 1;
    cam = 0;
    foptr = 0;
    sel_ptr = 0;
    aldaketa = 'r';
    ald_lokala = 1;
    perspectiva();
    // inicializar_camara();
    printf("Leyendo camara...\n");
    read_from_file("cam.txt", 1);
    read_from_file("k.txt", 0);
    sel_ptr->mptr->m[3] = -200;
    if (argc > 1)
    {
        read_from_file(argv[1], 0);
    }
    else
    {
        read_from_file("z.txt", 0);
    }

    sel_ptr->mptr->m[3] = 200;
    glutMainLoop();
    return 0;
}
