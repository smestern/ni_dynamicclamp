#include <stdio.h>
#include "interface_c.c"

int main() {
    init_ni(1e-4, 0.1, 0.5);
    // Test case 1
    double I1[] = {0.1, 0.2, 0.3, 0.4, 0.5};
    double out1[5];
    run_step_loop(I1, out1, 5);
    printf("Test case 1:\n");
    for (int i = 0; i < 5; i++) {
        printf("out[%d] = %lf\n", i, out1[i]);
    }
    
    // Test case 2
    double I2[] = {0.5, 0.4, 0.3, 0.2, 0.1};
    double out2[5];
    run_step_loop(I2, out2, 5);
    printf("Test case 2:\n");
    for (int i = 0; i < 5; i++) {
        printf("out[%d] = %lf\n", i, out2[i]);
    }
    
    // Test case 3
    double I3[] = {0.0, 0.0, 0.0, 0.0, 0.0};
    double out3[5];
    run_step_loop(I3, out3, 5);
    printf("Test case 3:\n");
    for (int i = 0; i < 5; i++) {
        printf("out[%d] = %lf\n", i, out3[i]);
    }
    
    return 0;
}