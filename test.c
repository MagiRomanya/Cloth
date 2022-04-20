#include <stdio.h>

typedef struct {
    float x, y;
}Particle;

int test(Particle particles[]){
    float* y[4];
    for (int i=0; i<2; i+=2) {
        y[i] = &particles[i].x;
        y[i+1] = &particles[i].y;
    }
    ++*y[0];
    return 0;
}

int main(void)
{
    Particle particles[2];
    for (int i=0; i < 2; i++) {
        particles[i].x = 0;
        particles[i].y = 0;
    }
    test(particles);
    test(particles);
    test(particles);
    test(particles);
    test(particles);
    printf("X=%f\n",particles[0].x);
    printf("X=%f\n",particles[1].x);
    return 0;
}
