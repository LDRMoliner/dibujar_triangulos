#include <stdio.h>
#include <string.h>
#include <math.h>
#include "cargar-triangulo.h"

int print_matrizea(char *str, double to_print[16])
{
    int i;
    printf("%s\n", str);
    for (i = 0; i < 4; i++)
        printf("%lf, %lf, %lf, %lf\n", to_print[i * 4], to_print[i * 4 + 1], to_print[i * 4 + 2],
               to_print[i * 4 + 3]);
    return 0;
}

int obtener_rotacion_rodrigues(double x, double y, double z, double angulo, double m[16])
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
    return 0;
}

int mxm(double resultado[16], double operando_izquierdo[16], double operando_derecho[16])
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
    return 0;
}

int normalizar(double vector[3])
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
    return 0;
}

int mxp(punto *pptr, double m[16], punto p)
{
    // print_matrizea("");
    pptr->x = p.x * m[0] + p.y * m[1] + p.z * m[2] + p.w * m[3];
    pptr->y = p.x * m[4] + p.y * m[5] + p.z * m[6] + p.w * m[7];
    pptr->z = p.x * m[8] + p.y * m[9] + p.z * m[10] + p.w * m[11];
    pptr->w = p.x * m[12] + p.y * m[13] + p.z * m[14] + p.w * m[15];
    pptr->u = p.u;
    pptr->v = p.v;
    return 0;
}

// void equation_plane(float x1, float y1, 
//                     float z1, float x2,
//                     float y2, float z2, 
//                     float x3, float y3, float z3)
// {   
//     double v1[3];
//     float a1 = x2 - x1;
//     float b1 = y2 - y1;
//     float c1 = z2 - z1;
//     float a2 = x3 - x1;
//     float b2 = y3 - y1;
//     float c2 = z3 - z1;
//     float a = b1 * c2 - b2 * c1;
//     float b = a2 * c1 - a1 * c2;
//     float c = a1 * b2 - b1 * a2;
//     float d = (- a * x1 - b * y1 - c * z1);
   
//     v1[0] = a;
//     v1[1] = b;
//     v1[2] = c;

//     normalizar(v1);
//      printf("equation of plane is %.2f x + %.2f"
//         " y + %.2f z + %.2f = 0.",a,b,c,d);
//     printf("normal vector is %.2f %.2f %.2f\n", v1[0], v1[1], v1[2]);
//     return;
// }