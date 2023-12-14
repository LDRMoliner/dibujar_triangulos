#include <stdio.h>
#include <string.h>

int main() {
    int sourceArray[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
    int destinationArray[16] = {0,0};  // Ensure that the destination array is large enough

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
