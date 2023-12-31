#ifndef OPERACIONES_VECTORIALES_H
#define OPERACIONES_VECTORIALES_H

#include "cargar-triangulo.h"

int mxm(double resultado[16], double operando_izquierdo[16], double operando_derecho[16]);
int normalizar(double vector[3]);
punto calcular_normal(punto p1, punto p2, punto p3);
int obtener_rotacion_rodrigues(double x, double y, double z, double angulo, double m[16]);
int print_matrizea(char *str, double to_print[16]);
int mxp(punto *pptr, double m[16], punto p);
#endif

