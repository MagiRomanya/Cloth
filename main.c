#include "raylib.h"
#include <stdio.h>
#include "math.h"

#define screenWidth 1200
#define screenHeight 900

typedef struct {
  float x, y, vx, vy;
  float r;
  int id;
}Particle;


float InvSqrt(float x)
{
        float xhalf = 0.5f * x;
        int i = *(int*)&x;            // store floating-point bits in integer
        i = 0x5f3759df - (i >> 1);    // initial guess for Newton's method
        x = *(float*)&i;              // convert new bits into float
        x = x*(1.5f - xhalf*x*x);     // One round of Newton's method
        return x;
}

int Springs(float x, const Particle particles[], float f[], int size);

Particle genParticle(float x, float y, float vx, float vy, int* Nparticles){
  Particle particle;
  particle.id = *Nparticles;
  (*Nparticles)++;
  particle.x = x;
  particle.y = y;
  particle.vx = vx;
  particle.vy = vy;
  particle.r = 10.f;
  return particle;
}
void drawParticles(Particle particles[], int n){
  Color colors[]= {
    DARKPURPLE,
  };
  int n_colors=sizeof(colors)/sizeof(RED);
  for (int i =0; i<n; i++){
    DrawCircle(particles[i].x, particles[i].y, particles[i].r, colors[i%n_colors]);
    //DrawText(TextFormat("%i", particles[i].id), particles[i].x, particles[i].y, 10, LIGHTGRAY);
  }
}

// BUG: The performance using this routine is different from the previous one. This should not happen. The integration method is the same.
// TODO: Test if a difference exists between the routines.
// TODO: Test this runge kutta in a simple scenario: an harmonic oscilator. Compeare the results with the last runge kutta routine.
void rungekutta4(const Particle particles[], float x, int size, int (*func)(float, const Particle*, float*, int), float h, Particle p2[])
{
	int i;
	float K1[size], K2[size], K3[size], K4[size];

  float* y[size]; // links to p2 coordintes
  for (i=0; i<size; i+=4) {
    // Particles2 == Particles  (there are 4 coordiantes per particle)
    p2[i/4].x = particles[i/4].x;
    p2[i/4].y = particles[i/4].y;
    p2[i/4].vx = particles[i/4].vx;
    p2[i/4].vy = particles[i/4].vy;

    // Linking y[] to p2 coordinates
    y[i] = &p2[i/4].x;
    y[i+1] = &p2[i/4].y;
    y[i+2] = &p2[i/4].vx;
    y[i+3] = &p2[i/4].vy;
  }
	// Sets K1 = f(x,y)
	func(x,particles,K1,size);

	// Sets K2 = f(x + h/2, y + h/2*K1)
	for (i=0; i<size; i++)
	{
		*y[i] += h/2*K1[i];
	}
	func(x + h/2, p2, K2, size);

	// Sets K3 = f(x + h/2, y + h/2*K2)
	for (i=0; i<size; i++)
	{
    *y[i] += h/2*K2[i];
	}
	func(x + h/2, p2, K3, size);

	// Sets K4 = f(x + h, y + h*K3)
	for (i=0; i<size; i++)
	{
		*y[i] += h*K3[i];
	}
	func(x + h, p2, K4, size);

	// Computes the final answer
	for (i=0; i<size; i++)
	{
		*y[i] += h/6*(K1[i] + 2*(K2[i]+K3[i]) + K4[i]);
	}
}


Vector2 spring(Particle p1, Particle p2, float K0, float l0)
{
  Vector2 force;
  Vector2 radialVec = (Vector2){p2.x-p1.x, p2.y-p1.y}; //
  float dist = radialVec.x*radialVec.x + radialVec.y*radialVec.y; // Distance squared

  force.x = K0*(1-l0*InvSqrt(dist))*radialVec.x;
  force.y = K0*(1-l0*InvSqrt(dist))*radialVec.y;

  force.x = K0*radialVec.x;
  force.y = K0*radialVec.y;
  return force;
}

int FirstNeighbors();

void genParticleGrid(Particle particles[] ,int* Nparticles, int nx, int ny, float deltax, float deltay){
  float x0 = 180.0f;
  float y0 = 50.0f;
  for (int i=0; i < nx; i++){
    for (int j=0; j < ny; j++) {
      particles[i*ny + j] = genParticle(x0+deltax*i, y0 +deltay*j, 0.0f, 0.0f, Nparticles);
    }
  }
  FirstNeighbors();
}
#define Nx 40
#define Ny 40
#define DeltaX 20.f
#define DeltaY 20.f

bool lattice[Nx*Ny][Nx*Ny];

int FirstNeighbors(){
  // Loop in the grid: i,j are the coordinates in the grid (x->i, y->j)
  for (int i=0; i<Nx-1; i++) {
    for (int j=0; j<Ny-1; j++) {
      // particle index = j + Nx*i
      lattice[j+Nx*i][j+1+Nx*i] = true;
      lattice[j+Nx*i][j+Nx*(i+1)] = true;
    }
  lattice[Ny-1 +Nx*i][Ny-1 + Nx*(i+1)] = true;
  lattice[i +Nx*(Nx-1)][i+1 + Nx*(Nx-1)] = true;
  }
  return 0;

}
int Springs(float x, const Particle particles[], float f[], int size){
  const float K = 25.0f;
  const float dump = 2.0e0f;
  const float d0 = DeltaX;
  const float g0 = -10; // Gravetat
  Vector2 force;
  int i,j;

  // Set the velocities
  for (i=0; i<size; i+=4) {
    f[i] = particles[i/4].vx;
    f[i+1] = particles[i/4].vy;
    f[i+2] = 0.0f;
    f[i+3] = 0.0f;
  }
  // Set the accelerations
  // One body interactions:
  for (i = 0; i<size; i+=4) {
    f[i+2] += -dump*particles[i/4].vx;
    f[i+3] += -g0 -dump*particles[i/4].vy;
  }
  // Two body interaction
  for (i=0; i<size; i+=4) {
    for (j=i; j<size; j+=4) {
      if (lattice[i/4][j/4]) {
        force = spring(particles[i/4], particles[j/4], K, d0/2);
        f[i+2] += force.x;
        f[i+3] += force.y;
        f[j+2] -= force.x;
        f[j+3] -= force.y;
      }
    }
  }
  f[0+2] = 0;
  f[0+3] = 0;
  f[4*(0+Nx*(Nx-1))+2] = 0;
  f[4*(0+Nx*(Nx-1))+3] = 0;
  return 0;
}

int main(void){
  printf("Hello World!\n");

  InitWindow(screenWidth, screenHeight, "Cloth");

  Particle particles[50*50];
  Particle particles2[50*50];
  float time = 0.0f;
  float h = 1e-2f;
  int Nparticles = 0;
  genParticleGrid(particles, &Nparticles, Nx, Ny, DeltaX, DeltaY);
  printf("Nparticles =%i, size=%i\n",Nparticles, Nparticles*4);

  for (int i=0; i<Nparticles; i++) {
    particles2[i] = particles[i];
  }

  int pause = 0;
  SetTargetFPS(60);
  while (!WindowShouldClose()) {
    if (IsKeyPressed(KEY_SPACE)) {
      if (pause) pause = 0;
      else pause = 1;
    }
    if (!pause) {
      rungekutta4(particles, time, Nparticles*4, Springs, h, particles2);
      for (int i=0; i<Nparticles; i++) {
        particles[i] = particles2[i];
      }
      rungekutta4(particles, time, Nparticles*4, Springs, h, particles2);
      for (int i=0; i<Nparticles; i++) {
        particles[i] = particles2[i];
      }
    }
    BeginDrawing();
        ClearBackground(RAYWHITE);
        drawParticles(particles, Nparticles);
    EndDrawing();
  }
  CloseWindow();
}
