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

// testuraren informazioa
// información de textura

extern int load_ppm(char *file, unsigned char **bufferptr, int *dimxptr, int *dimyptr);
unsigned char *bufferra;
int dimx, dimy;
int indexx;
hiruki *triangulosptr;
triobj *foptr;
triobj *sel_ptr;
int denak;
int lineak;
int objektuak;
char aldaketa;
int ald_lokala;

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
unsigned char *color_textura(float u, float v)
{
    char *lag;

    int desplazamiento_u, desplazamiento_v;

    desplazamiento_u = trunc (u * dimx);
    desplazamiento_v = trunc ((1-v)*dimy);

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
    float difu, difv;
    unsigned char r, g, b;
    unsigned char *colorv;

    glBegin(GL_POINTS);

    if (c1x != c2x)
    {
        difv = (c2v - c1v) / (c2x - c1x);
        difu = (c2u - c1u) / (c2x - c1x);
    }
    else
    {
        difv = 0;
        difu = 0;
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
    }
    glEnd();
}

void print_matrizea(char *str)
{
    int i;

    printf("%s\n", str);
    for (i = 0; i < 4; i++)
        printf("%lf, %lf, %lf, %lf\n", sel_ptr->mptr->m[i * 4], sel_ptr->mptr->m[i * 4 + 1], sel_ptr->mptr->m[i * 4 + 2],
               sel_ptr->mptr->m[i * 4 + 3]);
}

// TODO
// aurrerago egitekoa
// para más adelante
void mxp(punto *pptr, double m[16], punto p)
{
    pptr->x = p.x;
    pptr->y = p.y;
    pptr->z = p.z;
    pptr->u = p.u;
    pptr->v = p.v;
}

void calcula_punto_corte(punto *punto_superior, punto *punto_inferior, int i, punto *corte)
{
    float m = 0, m2 = 0, m3 = 0, m4 = 0;

    m = (punto_inferior->y - punto_superior->y) / (punto_inferior->x - punto_superior->x);
    m2 = (punto_inferior->y - punto_superior->y) / (punto_inferior->u - punto_superior->u);
    m3 = (punto_inferior->y - punto_superior->y) / (punto_inferior->v - punto_superior->v);
    m4 = (punto_inferior->y - punto_superior->y) / (punto_inferior->z - punto_superior->z);
    corte->x = ((i - punto_superior->y) / m) + punto_superior->x;
    corte->u = ((i - punto_superior->y) / m2) + punto_superior->u;
    corte->v = ((i - punto_superior->y) / m3) + punto_superior->v;
    corte->z = ((i - punto_superior->y) / m4) + punto_superior->z;
}

void encontrar_max_min(hiruki *tptr, float y1, float y2, float y3, punto **pgoiptr, punto **perdiptr, punto **pbeheptr)
{
    if (tptr->p1.y > tptr->p2.y)
    {
        if (tptr->p1.y > tptr->p3.y)
        {
            *pgoiptr = &(tptr->p1);
            if (tptr->p2.y > tptr->p3.y)
            {
                *perdiptr = &(tptr->p2);
                *pbeheptr = &(tptr->p3);
            }
            else
            {
                *perdiptr = &(tptr->p3);
                *pbeheptr = &(tptr->p2);
            }
        }
        else
        {
            *pgoiptr = &(tptr->p3);
            *perdiptr = &(tptr->p1);
            *pbeheptr = &(tptr->p2);
        }
    }
    else if (tptr->p2.y > tptr->p3.y)
    {
        *pgoiptr = &(tptr->p2);
        *pbeheptr = &(tptr->p1);
        *perdiptr = &(tptr->p3);
    }
    else
    {
        *pgoiptr = &(tptr->p3);
        *pbeheptr = &(tptr->p1);
        *perdiptr = &(tptr->p2);
    }
}

void rellenar_triangulo(punto *pgoiptr, punto *perdiptr, punto *pbeheptr)
{
    int i;
    punto corte1;
    punto corte2;

    for (i = pgoiptr->y; i > perdiptr->y; i--)
    {
        calcula_punto_corte(pgoiptr, pbeheptr, i, &corte1);
        calcula_punto_corte(pgoiptr, perdiptr, i, &corte2);
        if (corte1.x <= corte2.x)
            dibujar_linea_z(i, corte1.x, corte1.z, corte1.u, corte1.v, corte2.x, corte2.z, corte2.u, corte2.v);
        else
            dibujar_linea_z(i, corte2.x, corte2.z, corte2.u, corte2.v, corte1.x, corte1.z, corte1.u, corte1.v);
    }

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
    int start1, star2;
    float t = 1, s = 1, q = 1;
    float decremento_t = 1, decremento_s = 1, decremento_q = 1;
    float c1x, c1z, c1u, c1v, c2x, c2z, c2u, c2v;
    int linea;
    float cambio1, cambio1z, cambio1u, cambio1v, cambio2, cambio2z, cambio2u, cambio2v;
    punto p1, p2, p3;

    if (i >= optr->num_triangles)
        return;
    tptr = optr->triptr + i;
    mxp(&p1, optr->mptr->m, tptr->p1);
    mxp(&p2, optr->mptr->m, tptr->p2);
    mxp(&p3, optr->mptr->m, tptr->p3);
    if (lineak == 1)
    {
        glBegin(GL_POLYGON);
        glVertex3d(p1.x, p1.y, p1.z);
        glVertex3d(p2.x, p2.y, p2.z);
        glVertex3d(p3.x, p3.y, p3.z);
        glEnd();
        return;
    }
   
    //Encontramos el punto máximo, mínimo, y medio del triángulo.

    encontrar_max_min(tptr, tptr->p1.y, tptr->p2.y, tptr->p3.y, &pgoiptr, &perdiptr, &pbeheptr);

    // Ahora vamos a calcular los puntos de corte para dibujar líneas.

    rellenar_triangulo(pgoiptr, perdiptr, pbeheptr);
    
}

static void marraztu(void)
{
    float u, v;
    int i, j;
    triobj *auxptr;
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
        glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
    else
    {
        if (denak == 0)
            glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
    }

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-500.0, 500.0, -500.0, 500.0, -500.0, 500.0);

    triangulosptr = sel_ptr->triptr;
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

    // printf("%s fitxategitik datuak hartzera\n",fitx);
    optr = (triobj *)malloc(sizeof(triobj));
    retval = cargar_triangulos(fitx, &(optr->num_triangles), &(optr->triptr));
    if (retval != 1)
    {
        printf("%s fitxategitik datuak hartzerakoan arazoak izan ditut\n    Problemas al leer\n", fitxiz);
        free(optr);
    }
    else
    {
        triangulosptr = optr->triptr;
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
    printf("datuak irakurrita\nLecura finalizada\n");
}

void x_aldaketa(int dir)
{
}

void y_aldaketa(int dir)
{
}

void z_aldaketa(int dir)
{
}

void undo()
{
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

int main(int argc, char **argv)
{
    int retval;

    printf(" Triangeluak: barneko puntuak eta testura\n Triángulos con puntos internos y textura \n");
    printf("Press <ESC> to finish\n");
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(500, 500);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("KBG/GO praktika");

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
    denak = 0;
    lineak = 0;
    objektuak = 0;
    foptr = 0;
    sel_ptr = 0;
    aldaketa = 'r';
    ald_lokala = 1;
    if (argc > 1)
        read_from_file(argv[1]);
    else
        read_from_file("adibideak.txt");
    glutMainLoop();

    return 0;
}
