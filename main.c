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
int gravity(float x, const float y[], float f[], int size)
{
	int i, j, k;
	const int SIZE = size;
	const int SIZE2 = size/2; // SIZE/2
	const int N=SIZE/4; // We have N particles, DIM coordinates and DIM velocities => number of dimentions of y[] f[]
	const float G = 1.0e3f; // Gravitational constant multiplied by m
  const int DIM = 2;
	float sum2, dist;


	// Initialize f
	for (i=0; i<SIZE2; i++)
	{
		f[i] = y[SIZE2+i]; // dz/dx = v_z
		f[SIZE2+i] = 0;
	}
	for (i=0; i<SIZE2 -DIM; i+=DIM)
	{
		for (j=i+DIM; j<SIZE2; j+=DIM)
		{
			sum2 = 0;
			// CALCULATE THE DISTANCE
			for (k=0; k<DIM; k++)
			{
				sum2 +=(y[j+k]-y[i+k])*(y[j+k]-y[i+k]); // x**2 + y**2 + z**2
			}
			dist = InvSqrt(sum2);
			// COMPUTES THE ACCELERATION
			for (k=0; k<DIM; k++)
			{
				double f_module = G*(dist*dist*dist)*(y[j+k]-y[i+k]);
				f[SIZE2+i+k] += f_module;
				f[SIZE2+j+k] += -f_module;
			}
		}
	}
	return 0;
}
int Springs(float x, const float y[], float f[], int size);
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
    PINK,
    RED,
    MAROON,
    GREEN,
    LIME,
    BLUE,
    DARKBLUE,
    SKYBLUE,
    PURPLE,
    VIOLET
  };
  int n_colors=sizeof(colors)/sizeof(RED);
  for (int i =0; i<n; i++){
    DrawCircle(particles[i].x, particles[i].y, particles[i].r, colors[i%n_colors]);
    //DrawText(TextFormat("%i", particles[i].id), particles[i].x, particles[i].y, 10, LIGHTGRAY);
  }
}
// Calculates a step of RK4
//
//	EXAMPLE: Harmonic oscilator
//
//	diferential eq  dv_z/dt = -kz     (where v_z = dz/dt)
//		        dz/dt = v_z
//
//	=> y = (v_z, z)	=> dy/dx = f(x,y) = (dv_z/dz, dz/dx) = ( -kz, v_z))
//	   x = t
// y[]       Dependent variables (position & velociteies)
// x         Independent variable
// size      Number of dependent variables
// func      Function which returns the value of the derivative
// h         Step size
// y1[]      Updated dependent variables for the next step
//
void rungekutta4(const float y[], float x, int size, int (*func)(float, const float*, float*, int), float h, float* y1)
{
	int i;
	float K1[size], K2[size], K3[size], K4[size], ycache[size];

	// Sets K1 = f(x,y)
	func(x,y,K1,size);

	// Sets K2 = f(x + h/2, y + h/2*K1)
	for (i=0; i<size; i++)
	{
		ycache[i] = y[i] + h/2*K1[i];
	}
	func(x + h/2, ycache, K2, size);

	// Sets K3 = f(x + h/2, y + h/2*K2)
	for (i=0; i<size; i++)
	{
		ycache[i] = y[i] + h/2*K2[i];
	}
	func(x + h/2, ycache, K3, size);

	// Sets K4 = f(x + h, y + h*K3)
	for (i=0; i<size; i++)
	{
		ycache[i] = y[i] + h*K3[i];
	}
	func(x + h, ycache, K4, size);

	// Computes the final answer
	for (i=0; i<size; i++)
	{
		y1[i] = y[i] + h/6*(K1[i] + 2*(K2[i]+K3[i]) + K4[i]);
	}
}

void calculateNextStep_rk4(Particle particles[], int n,float time, float dt){
  float y[n*4];
  float y1[n*4];
  for (int i =0; i<n*4; i++) y1[i] = 0.0f;
  // Creates the y[] output vector
  for (int i=0; i<n; i++){
    y[i*2] = particles[i].x;
    y[i*2+1] = particles[i].y;
  }
  for (int i=0; i<n; i++){
    y[n*2+i*2] = particles[i].vx;
    y[n*2+i*2+1] = particles[i].vy;
  }
  rungekutta4(y, time, 4*n, Springs, dt, y1);
  // Assigns the y1[] to the future coordinates and velocities
  for (int i=0; i<n; i++){
    particles[i].x = y1[i*2];
    particles[i].y = y1[i*2+1];
  }
  for (int i=0; i<n; i++){
    particles[i].vx = y1[n*2+i*2];
    particles[i].vy = y1[n*2+i*2+1];
  }
}

void genParticleGrid(Particle particles[] ,int* Nparticles, int nx, int ny, float deltax, float deltay){
  float x0 = 180.0f;
  float y0 = 50.0f;
  for (int i=0; i < nx; i++){
    for (int j=0; j < ny; j++) {
      particles[i*ny + j] = genParticle(x0+deltax*i, y0 +deltay*j, 0.0f, 0.0f, Nparticles);
    }
  }
}
#define Nx 40
#define Ny 40
#define DeltaX 20.f
#define DeltaY 20.f
#define index(i,j) ((i)*Ny + j)

int Springs(float x, const float y[], float f[], int size){
  const float K = 25.0f;
  const float dump = 1.0e-1f;
  const float d0 = DeltaX;
  const float g0 = -5;
  for (int i = 0; i<size; i++) f[i] = 0.0f;
  for (int i = 0; i<size/2; i++) { // Velocity is the time derivative of sapce coordinates
    f[i] = y[size/2 + i];
  }
  for (int i = 1; i<Nx-1; i+=1) {
    for (int j = 1; j<Ny-1; j+=1) {
      // X coordinate
      f[size/2 + index(i,j)*2] += K*( (y[index(i+1,j)*2] - y[index(i,j)*2] ) + (y[index(i-1,j)*2] - y[index(i,j)*2] ) ) - (y[index(i,j)*2+size/2])*dump;
      f[size/2 + index(i,j)*2] += K*( (y[index(i,j+1)*2] - y[index(i,j)*2] ) + (y[index(i,j-1)*2] - y[index(i,j)*2] ) ) - (y[index(i,j)*2+size/2])*dump;
      // Y coordinate
      f[size/2 + index(i,j)*2+1] += K*( (y[index(i+1,j)*2+1] - y[index(i,j)*2+1] ) + (y[index(i-1,j)*2+1] - y[index(i,j)*2+1] ) ) - (y[index(i,j)*2+size/2+1])*dump - g0;
      f[size/2 + index(i,j)*2+1] += K*( (y[index(i,j+1)*2+1] - y[index(i,j)*2+1] ) + (y[index(i,j-1)*2+1] - y[index(i,j)*2+1] ) ) - (y[index(i,j)*2+size/2+1])*dump - g0;
    }
  }
  for (int j=1; j<Ny-1; j++){
    int i = 0;
    // LEFT COLUMN
      // X coordinate
      f[size/2 + index(i,j)*2] += K*( (y[index(i+1,j)*2] - y[index(i,j)*2] ) ) - (y[index(i,j)*2+size/2])*dump;
      f[size/2 + index(i,j)*2] += K*( (y[index(i,j+1)*2] - y[index(i,j)*2] ) + (y[index(i,j-1)*2] - y[index(i,j)*2] ) ) - (y[index(i,j)*2+size/2])*dump;
      // Y coordinate
      f[size/2 + index(i,j)*2+1] += K*( (y[index(i+1,j)*2+1] - y[index(i,j)*2+1] ) ) - (y[index(i,j)*2+size/2+1])*dump - g0;
      f[size/2 + index(i,j)*2+1] += K*( (y[index(i,j+1)*2+1] - y[index(i,j)*2+1] ) + (y[index(i,j-1)*2+1] - y[index(i,j)*2+1] ) ) - (y[index(i,j)*2+size/2+1])*dump -g0;
    // RIGHT COLUMN
      i = Nx-1;
      // X coordinate
      f[size/2 + index(i,j)*2] += K*( (y[index(i-1,j)*2] - y[index(i,j)*2] ) ) - (y[index(i,j)*2+size/2])*dump;
      f[size/2 + index(i,j)*2] += K*( (y[index(i,j+1)*2] - y[index(i,j)*2] ) + (y[index(i,j-1)*2] - y[index(i,j)*2] ) ) - (y[index(i,j)*2+size/2])*dump;
      // Y coordinate
      f[size/2 + index(i,j)*2+1] += K*( (y[index(i-1,j)*2+1] - y[index(i,j)*2+1] ) ) - (y[index(i,j)*2+size/2+1])*dump -g0;
      f[size/2 + index(i,j)*2+1] += K*( (y[index(i,j+1)*2+1] - y[index(i,j)*2+1] ) + (y[index(i,j-1)*2+1] - y[index(i,j)*2+1] ) ) - (y[index(i,j)*2+size/2+1])*dump -g0;
  }
  for (int i=1; i<Nx-1;i++){
    int j = Ny-1;
      // X coordinate
      f[size/2 + index(i,j)*2] += K*( (y[index(i+1,j)*2] - y[index(i,j)*2] ) + (y[index(i-1,j)*2] - y[index(i,j)*2] ) ) - (y[index(i,j)*2+size/2])*dump;
      f[size/2 + index(i,j)*2] += K*( (y[index(i,j-1)*2] - y[index(i,j)*2] ) ) - (y[index(i,j)*2+size/2])*dump;
      // Y coordinate
      f[size/2 + index(i,j)*2+1] += K*( (y[index(i+1,j)*2+1] - y[index(i,j)*2+1] ) + (y[index(i-1,j)*2+1] - y[index(i,j)*2+1] ) ) - (y[index(i,j)*2+size/2+1])*dump -g0;
      f[size/2 + index(i,j)*2+1] += K*( (y[index(i,j-1)*2+1] - y[index(i,j)*2+1] ) ) - (y[index(i,j)*2+size/2+1])*dump -g0;

    j = 0;
      // X coordinate
      f[size/2 + index(i,j)*2] += K*( (y[index(i+1,j)*2] - y[index(i,j)*2] ) + (y[index(i-1,j)*2] - y[index(i,j)*2] ) ) - (y[index(i,j)*2+size/2])*dump;
      f[size/2 + index(i,j)*2] += K*( (y[index(i,j+1)*2] - y[index(i,j)*2] ) ) - (y[index(i,j)*2+size/2])*dump;
      // Y coordinate
      f[size/2 + index(i,j)*2+1] += K*( (y[index(i+1,j)*2+1] - y[index(i,j)*2+1] ) + (y[index(i-1,j)*2+1] - y[index(i,j)*2+1] ) ) - (y[index(i,j)*2+size/2+1])*dump - g0;
      f[size/2 + index(i,j)*2+1] += K*( (y[index(i,j+1)*2+1] - y[index(i,j)*2+1] ) ) - (y[index(i,j)*2+size/2+1])*dump - g0;
  }
  // BOOTOM LEFT
  int i = 0;
  int j = Ny-1;
      // X coordinate
      f[size/2 + index(i,j)*2] += K*( (y[index(i+1,j)*2] - y[index(i,j)*2] ) ) - (y[index(i,j)*2+size/2])*dump;
      f[size/2 + index(i,j)*2] += K*( (y[index(i,j-1)*2] - y[index(i,j)*2] ) ) - (y[index(i,j)*2+size/2])*dump;
      // Y coordinate
      f[size/2 + index(i,j)*2+1] += K*( (y[index(i+1,j)*2+1] - y[index(i,j)*2+1] ) ) - (y[index(i,j)*2+size/2+1])*dump - g0;
      f[size/2 + index(i,j)*2+1] += K*( (y[index(i,j-1)*2+1] - y[index(i,j)*2+1] ) ) - (y[index(i,j)*2+size/2+1])*dump - g0;
  // BOTTOM RIGHT
  i = Nx-1;
      // X coordinate
      f[size/2 + index(i,j)*2] += K*( (y[index(i-1,j)*2] - y[index(i,j)*2] ) ) - (y[index(i,j)*2+size/2])*dump;
      f[size/2 + index(i,j)*2] += K*( (y[index(i,j-1)*2] - y[index(i,j)*2] ) ) - (y[index(i,j)*2+size/2])*dump;
      // Y coordinate
      f[size/2 + index(i,j)*2+1] += K*( (y[index(i-1,j)*2+1] - y[index(i,j)*2+1] ) ) - (y[index(i,j)*2+size/2+1])*dump - g0;
      f[size/2 + index(i,j)*2+1] += K*( (y[index(i,j-1)*2+1] - y[index(i,j)*2+1] ) ) - (y[index(i,j)*2+size/2+1])*dump - g0;
  return 0;
}

int main(void){
  printf("Hello World!\n");

  InitWindow(screenWidth, screenHeight, "Cloth");

  Particle particles[50*50];
  int Nparticles = 0;
  genParticleGrid(particles, &Nparticles, Nx, Ny, DeltaX, DeltaY);
  //particles[index(1,1)].x += 20;
  //particles[index(1,1)].y += 20;

  int pause = 1;
  SetTargetFPS(60);
  while (!WindowShouldClose()) {
    if (IsKeyPressed(KEY_SPACE)) {
      if (pause) pause = 0;
      else pause = 1;
    }
    if (!pause) {
      calculateNextStep_rk4(particles, Nparticles, 0.0f, 0.1f);
    }
    BeginDrawing();
        ClearBackground(RAYWHITE);
        drawParticles(particles, Nparticles);
    EndDrawing();
  }
  CloseWindow();
}
