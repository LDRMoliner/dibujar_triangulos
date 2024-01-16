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

// Matriz de rotación.
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
    // print_matrizea("rotacion rodrigues", m);
    return 0;
}

// Multiplicar matriz por matriz.
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

// Normalizamos cualquier punto vector pasado por parámetro.
int normalizar_vector(punto *vector)
{
    double length = 0.0;

    length += vector->x * vector->x + vector->y * vector->y + vector->z * vector->z;
    length = sqrt(length);

    vector->x = vector->x / length;
    vector->y = vector->y / length;
    vector->z = vector->z / length;

    return 0;
}

// Normalizamos cualquier vector pasado por parámetro.
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

// Matriz para multiplicar tanto puntos como vectores.
int mxp(punto *pptr, double m[16], punto p)
{
    // print_matrizea("");
    pptr->x = p.x * m[0] + p.y * m[1] + p.z * m[2] + p.w * m[3];
    pptr->y = p.x * m[4] + p.y * m[5] + p.z * m[6] + p.w * m[7];
    pptr->z = p.x * m[8] + p.y * m[9] + p.z * m[10] + p.w * m[11];
    // Distinguimos entre punto y vector de ésta forma.
    if (pptr->w != 0.0)
        pptr->w = p.x * m[12] + p.y * m[13] + p.z * m[14] + p.w * m[15];
    pptr->u = p.u;
    pptr->v = p.v;

    // printf("Punto: %f, %f, %f, %f\n", pptr->x, pptr->y, pptr->z, pptr->w);
    return 0;
}

// Calculamos normal del polígono conformado por los tres puntos pasados por parámetro. Mejor explicado en la documentación.
punto calcular_normal(punto p1, punto p2, punto p3)
{
    punto normal;
    double v1[4];
    double v2[4];
    double N[4];

    v1[0] = p2.x - p1.x;
    v1[1] = p2.y - p1.y;
    v1[2] = p2.z - p1.z;
    v1[3] = 0;

    v2[0] = p3.x - p1.x;
    v2[1] = p3.y - p1.y;
    v2[2] = p3.z - p1.z;
    v2[3] = 0;

    N[0] = v1[1] * v2[2] - v1[2] * v2[1];
    N[1] = -(v1[0] * v2[2] - v2[0] * v1[2]);
    N[2] = v1[0] * v2[1] - v2[0] * v1[1];

    // printf("Perpendicular entre: %f, %f, %f y %f, %f, %f: \n", v1[0], v1[1], v1[2], v2[0], v2[1], v2[2]);
    printf("%f, %f, %f\n", N[0], N[1], N[2]);
    normalizar(N);
    // printf("Normalizado: %f, %f, %f\n", N[0], N[1], N[2]);

    normal.x = N[0];
    normal.y = N[1];
    normal.z = N[2];
    normal.w = 0;

    return normal;
}