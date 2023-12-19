#include <stdio.h>
#include <string.h>
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

int main() {
    int sourceArray[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
    double destinationArray[16] = {0,0};  // Ensure that the destination array is large enough
    double prueba[16] = {0,0};

    mxm (prueba, sourceArray, sourceArray);
    
    // Copy the first four elements from sourceArray to destinationArray
    memcpy(destinationArray + 4, sourceArray, sizeof(int) * 4);

    // Print the result
    printf("Source Array: ");
    for (int i = 0; i < sizeof(sourceArray) / sizeof(sourceArray[0]); ++i) {
        printf("%d ", sourceArray[i]);
    }

    printf("\nDestination Array: ");
    for (int i = 0; i < sizeof(destinationArray) / sizeof(destinationArray[0]); ++i) {
        printf("%d ", destinationArray[i]);
    }

    return 0;
}
