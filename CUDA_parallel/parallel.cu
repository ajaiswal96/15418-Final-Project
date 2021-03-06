
using namespace std;

#include "parallel.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <float.h>
#include <assert.h>
#include <fstream>
#include <iostream>

#include <cuda.h>
#include <cuda_runtime.h>
#include <driver_functions.h>


#define GET_X(x) (2*x)
#define GET_Y(y) (2*y + 1)

#define THREADS_PER_BLOCK 1024
#define DAMPING 0.75
#define OVERLAP_COEFF 1.3

typedef struct Parameters {
  string fileName;
  int numFrames; // Number of frames
  int stepsPerFrame; //Steps per frame
  float size; // Particle size
  float dt; //Time step
  float density_ref; //Reference density
  float k; // Bulk modulus
  float viscocity;
  float g; /* Gravity strength */
} Parameters;

typedef struct CurrState {
  int numParticles;
  float mass;
  float* densities; //array of all the densities of the particles
  float* positions;
  float* velocities_half;
  float* velocities_full;
  float* accelerations;
  int numBlocks;

  int* CUDA_numParticles; //Array of size 1
  int CUDA_tempInt; //Array of size 1
  float CUDA_tempFloat; //Array of size 1
  float* CUDA_velocities_half;
  float* CUDA_velocities_full;
  float* CUDA_accelerations;
  float* CUDA_positions;
  float* CUDA_densities;
} CurrState;

static inline int nextPow2(int n)
{
    n--;
    n |= n >> 1;
    n |= n >> 2;
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16;
    n++;
    return n;
}

CurrState* allocState(int n) {
  CurrState* currState = (CurrState*)calloc(1,sizeof(CurrState));
  currState->densities = (float*)calloc(n, sizeof(float));
  currState->positions = (float*)calloc(2*n, sizeof(float));
  currState->velocities_full = (float*)calloc(2*n, sizeof(float));
  currState->velocities_half = (float*)calloc(2*n, sizeof(float));
  currState->accelerations = (float*)calloc(2*n, sizeof(float));
  currState->numParticles = n;
  currState->numBlocks = 2*currState->numParticles / THREADS_PER_BLOCK + 1;

  cudaMalloc((void**)&(currState->CUDA_positions), (2*n*(sizeof(float))));
  cudaMalloc((void**)&(currState->CUDA_densities), n*sizeof(float));
  cudaMalloc((void**)&(currState->CUDA_numParticles), sizeof(int));
  cudaMalloc((void**)&(currState->CUDA_velocities_full),(2*n*(sizeof(float))));
  cudaMalloc((void**)&(currState->CUDA_velocities_half),(2*n*(sizeof(float))));
  cudaMalloc((void**)&(currState->CUDA_accelerations),(2*n*(sizeof(float))));

  cudaMemset(currState->CUDA_positions, 0, 2*n*sizeof(float));
  cudaMemset(currState->CUDA_velocities_full, 0, 2*n*sizeof(float));
  cudaMemset(currState->CUDA_velocities_half, 0, 2*n*sizeof(float));
  cudaMemset(currState->CUDA_accelerations, 0, 2*n*sizeof(float));
  cudaMemset(currState->CUDA_numParticles, 0, sizeof(int));

  return currState;
}

void freeState(CurrState* currState) {
  free(currState->densities);
  free(currState->positions);
  free(currState->velocities_half);
  free(currState->velocities_full);
  free(currState->accelerations);
  free(currState);

  cudaFree(currState->CUDA_positions);
  //cudaFree(currState->CUDA_tempFloat);
  //cudaFree(currState->CUDA_tempInt);
  cudaFree(currState->CUDA_numParticles);
  cudaFree(currState->CUDA_velocities_full);
  cudaFree(currState->CUDA_velocities_half);
  cudaFree(currState->CUDA_accelerations);
}

__global__ void kernel_velocityStep(float* CUDA_velocities_half, float* CUDA_velocities_full, float* CUDA_accelerations, int CUDA_numParticles, float* CUDA_positions, float CUDA_dt) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i >= 2* CUDA_numParticles) {
    //printf("FUCKME");
    return;
  }
  CUDA_velocities_half[i] += CUDA_accelerations[i] * CUDA_dt;
  CUDA_velocities_full[i] = CUDA_velocities_half[i] + CUDA_accelerations[i] * CUDA_dt / 2;
  CUDA_positions[i] += CUDA_velocities_half[i] * CUDA_dt;
}

__global__ void kernel_calculateDensity(int CUDA_numParticles, float CUDA_innerConstant, float CUDA_outerConstant, float CUDA_size_2, float* CUDA_densities, float* CUDA_positions) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i >= CUDA_numParticles) {
    return;
  }

  CUDA_densities[i] += CUDA_innerConstant;
  for (unsigned int j = 0; j < CUDA_numParticles; j++) {
    if (j != i) {
    float dx = CUDA_positions[GET_X(i)] - CUDA_positions[GET_X(j)];
    float dy = CUDA_positions[GET_Y(i)] - CUDA_positions[GET_Y(j)];
    float r_2 = dx * dx + dy * dy;
    float z = CUDA_size_2 - r_2;
    float z_3 = z * z * z;
    if (z > 0) {
      float densities_ij = CUDA_outerConstant * z_3;
      CUDA_densities[i] += densities_ij;
    }
    }
  }
}

__global__ void kernel_calculateAcceleration(int CUDA_numParticles, float CUDA_g, float CUDA_size, float CUDA_C0, float CUDA_Cp, float CUDA_Cv, float CUDA_density_ref, float* CUDA_positions, float* CUDA_accelerations, float* CUDA_densities, float* CUDA_velocities_full){
  int i = blockIdx.x * blockDim.x+ threadIdx.x;
  if (i >= CUDA_numParticles) {
    return;
  }
  float size_2 = CUDA_size * CUDA_size;

  CUDA_accelerations[GET_X(i)] = 0;
  CUDA_accelerations[GET_Y(i)] = -CUDA_g;

  float currDensity_i = CUDA_densities[i];

  for (unsigned int j = 0; j < CUDA_numParticles; j++) {
    if (j!=i) {
      float dx = CUDA_positions[GET_X(i)] - CUDA_positions[GET_X(j)];
      float dy = CUDA_positions[GET_Y(i)] - CUDA_positions[GET_Y(j)];
      float r_2 = dx * dx + dy * dy;
      if (r_2 < size_2) {
        const float currDensity_j = CUDA_densities[j];
        float q = sqrt(r_2)/CUDA_size;
        float u = 1-q;
        float w0 = CUDA_C0 * u/(currDensity_j * currDensity_i);
        float wp = w0 * CUDA_Cp * (currDensity_i + currDensity_j - 2 * CUDA_density_ref) * u/q;
        float wv = w0 * CUDA_Cv;
        float dvx = CUDA_velocities_full[GET_X(i)] - CUDA_velocities_full[GET_X(j)];
        float dvy = CUDA_velocities_full[GET_Y(i)] - CUDA_velocities_full[GET_Y(j)];
        if (i > j) {
          CUDA_accelerations[GET_X(i)] -= (wp * dx + wv * dvx);
          CUDA_accelerations[GET_Y(i)] -= (wp * dy + wv * dvy);
        }
        else {
          CUDA_accelerations[GET_X(i)] += (wp * dx + wv * dvx);
          CUDA_accelerations[GET_Y(i)] += (wp * dy + wv * dvy);
        }
      }

    }
  }

}

//now we need to compute densities. We are going to compute densities once for
//each ij ji pair since they are the same.
void calculateDensity(Parameters* params, CurrState* currState) {
  int numParticles = currState->numParticles;
  float size = params->size;
  float mass = currState->mass;

  float* positions = currState->positions;
  float* densities = currState->densities;

  memset(densities, 0, numParticles * sizeof(float));

  float size_2 = size * size;
  float size_8 = size_2 * size_2 * size_2 * size_2;
  float outerConstant = (4 * mass) / (M_PI * size_8);
  float innerConstant = (4 * mass) / (M_PI * size_2);

  float* CUDA_positions = currState->CUDA_positions;
  float* CUDA_densities = currState->CUDA_densities;

  //cudaMemcpy(currState->CUDA_densities, currState->densities, sizeof(float) * numParticles, cudaMemcpyHostToDevice);
  //cudaMemcpy(currState->CUDA_positions, currState->positions, sizeof(float) * numParticles * 2, cudaMemcpyHostToDevice);

  int numBlocks = numParticles / THREADS_PER_BLOCK + 1;

  kernel_calculateDensity<<<numBlocks,THREADS_PER_BLOCK>>>(numParticles, innerConstant, outerConstant, size_2, CUDA_densities, CUDA_positions);

  //cudaMemcpy(currState->densities, currState->CUDA_densities, sizeof(float) * numParticles, cudaMemcpyDeviceToHost);
  //cudaMemcpy(currState->positions, currState->CUDA_positions, sizeof(float) * numParticles * 2, cudaMemcpyDeviceToHost);

}

void calculateAcceleration(Parameters* params, CurrState* currState) {
  float size = params->size;
  float g = params->g;
  float k = params->k;

  float viscocity = params->viscocity;
  float density_ref = params->density_ref;
  int numParticles = currState->numParticles;
  float mass = currState->mass;


  float* accelerations = currState->accelerations;
  float* densities = currState->densities;
  float* positions = currState->positions;
  float* velocities = currState->velocities_full;

  float size_2 = size * size;
  float size_4 = size_2 * size_2;
  
  calculateDensity(params, currState);
  
  float C0 = mass / (M_PI * size_4);
  float Cp = 15 * k;
  float Cv = -40 * viscocity;

  float* CUDA_positions = currState->CUDA_positions;
  float* CUDA_accelerations = currState->CUDA_accelerations;
  float* CUDA_densities = currState->CUDA_densities;
  float* CUDA_velocities_full = currState->CUDA_velocities_full;

  //cudaMemcpy(currState->CUDA_positions, currState->positions, sizeof(float) * numParticles * 2, cudaMemcpyHostToDevice);
  //cudaMemcpy(currState->CUDA_accelerations, currState->accelerations, sizeof(float) * numParticles * 2, cudaMemcpyHostToDevice);
  //cudaMemcpy(currState->CUDA_densities, currState->densities, sizeof(float) * numParticles, cudaMemcpyHostToDevice);
  //cudaMemcpy(currState->CUDA_velocities_full, currState->velocities_full, sizeof(float) * numParticles * 2, cudaMemcpyHostToDevice);

  int numBlocks = numParticles / THREADS_PER_BLOCK + 1;

  kernel_calculateAcceleration<<<numBlocks, THREADS_PER_BLOCK>>>(numParticles, g, size, C0, Cp, Cv, density_ref, CUDA_positions, CUDA_accelerations, CUDA_densities, CUDA_velocities_full);


  //cudaMemcpy(currState->accelerations, currState->CUDA_accelerations, sizeof(float) * numParticles * 2, cudaMemcpyDeviceToHost);

}

void reflect(int axis, float barrier, float* positions, float* velocities_full, float* velocities_half) {
  const float damping = DAMPING;

  if (velocities_full[axis] == 0) { //this means the particle has stopped
    return;
  }

  float dt = (positions[axis] - barrier) / velocities_full[axis];

  positions[0] -= velocities_full[0] * (1-damping) * dt;
  positions[1] -= velocities_full[1] * (1-damping) * dt;

  //reflect the positions
  positions[axis] = 2 * barrier - positions[axis];

  //reflect the velocities
  velocities_full[axis] = - velocities_full[axis];
  velocities_half[axis] = - velocities_half[axis];

  //damp the velocities
  velocities_half[0] *= damping;
  velocities_half[1] *= damping;
  velocities_full[0] *= damping;
  velocities_full[1] *= damping;
}

__global__ void boundaryCheckKernel(int CUDA_numParticles, float damping, float* CUDA_positions, float* CUDA_velocities_full, float* CUDA_velocities_half) {
 int i = blockIdx.x * blockDim.x+ threadIdx.x;
 if (i >= CUDA_numParticles) {
    return;
  }
  const float XMIN = 0.0;
  const float YMIN = 0.0;
  const float XMAX = 1.0;
  const float YMAX = 1.0;
  if (CUDA_positions[2*i] < XMIN) {
    if (CUDA_velocities_full[2*i] != 0) { //this means the particle has stopped

    float dt = (CUDA_positions[2*i] - XMIN) / CUDA_velocities_full[2*i];

    CUDA_positions[2*i] -= CUDA_velocities_full[2*i] * (1-damping) * dt;
    CUDA_positions[2*i+1] -= CUDA_velocities_full[2*i+1] * (1-damping) * dt;

    //reflect the positions
    CUDA_positions[2*i] = 2 * XMIN - CUDA_positions[2*i];

    //reflect the velocities
    CUDA_velocities_full[2*i] = - CUDA_velocities_full[2*i];
    CUDA_velocities_half[2*i] = - CUDA_velocities_half[2*i];

    //damp the velocities
    CUDA_velocities_half[2*i] *= damping;
    CUDA_velocities_half[2*i+1] *= damping;
    CUDA_velocities_full[2*i] *= damping;
    CUDA_velocities_full[2*i+1] *= damping;
    }
  } else if (CUDA_positions[2*i] > XMAX) {
    if (CUDA_velocities_full[2*i] != 0) { //this means the particle has stopped

    float dt = (CUDA_positions[2*i] - XMAX) / CUDA_velocities_full[2*i];

    CUDA_positions[2*i] -= CUDA_velocities_full[2*i] * (1-damping) * dt;
    CUDA_positions[2*i+1] -= CUDA_velocities_full[2*i+1] * (1-damping) * dt;

    //reflect the positions
    CUDA_positions[2*i] = 2 * XMAX - CUDA_positions[2*i];

    //reflect the velocities
    CUDA_velocities_full[2*i] = - CUDA_velocities_full[2*i];
    CUDA_velocities_half[2*i] = - CUDA_velocities_half[2*i];

    //damp the velocities
    CUDA_velocities_half[2*i] *= damping;
    CUDA_velocities_half[2*i+1] *= damping;
    CUDA_velocities_full[2*i] *= damping;
    CUDA_velocities_full[2*i+1] *= damping;
    }
  }


  if (CUDA_positions[2*i+1] < YMIN) {
    if (CUDA_velocities_full[2*i+1] != 0) { //this means the particle has stopped

    float dt = (CUDA_positions[2*i+1] - YMIN) / CUDA_velocities_full[2*i+1];

    CUDA_positions[2*i] -= CUDA_velocities_full[2*i] * (1-damping) * dt;
    CUDA_positions[2*i+1] -= CUDA_velocities_full[2*i+1] * (1-damping) * dt;

    //reflect the positions
    CUDA_positions[2*i+1] = 2 * YMIN - CUDA_positions[2*i+1];

    //reflect the velocities
    CUDA_velocities_full[2*i+1] = - CUDA_velocities_full[2*i+1];
    CUDA_velocities_half[2*i+1] = - CUDA_velocities_half[2*i+1];

    //damp the velocities
    CUDA_velocities_half[2*i] *= damping;
    CUDA_velocities_half[2*i+1] *= damping;
    CUDA_velocities_full[2*i] *= damping;
    CUDA_velocities_full[2*i+1] *= damping;  
    }
  } else if (CUDA_positions[2*i+1] > YMAX) {
    if (CUDA_velocities_full[2*i+1] != 0) { //this means the particle has stopped

    float dt = (CUDA_positions[2*i+1] - YMAX) / CUDA_velocities_full[2*i+1];

    CUDA_positions[2*i] -= CUDA_velocities_full[2*i] * (1-damping) * dt;
    CUDA_positions[2*i+1] -= CUDA_velocities_full[2*i+1] * (1-damping) * dt;

    //reflect the positions
    CUDA_positions[2*i+1] = 2 * YMAX - CUDA_positions[2*i+1];

    //reflect the velocities
    CUDA_velocities_full[2*i+1] = - CUDA_velocities_full[2*i+1];
    CUDA_velocities_half[2*i+1] = - CUDA_velocities_half[2*i+1];

    //damp the velocities
    CUDA_velocities_half[2*i] *= damping;
    CUDA_velocities_half[2*i+1] *= damping;
    CUDA_velocities_full[2*i] *= damping;
    CUDA_velocities_full[2*i+1] *= damping;
    }
  }
}

void boundaryCheck(CurrState* currState) {

  int numParticles = currState->numParticles;
  float* velocities_full = currState->velocities_full;
  float* velocities_half = currState->velocities_half;
  float* positions = currState->positions;

  float* CUDA_positions = currState->CUDA_positions;
  float* CUDA_velocities_half = currState->CUDA_velocities_half;
  float* CUDA_velocities_full = currState->CUDA_velocities_full;

  //cudaMemcpy(CUDA_positions, positions, sizeof(float) * numParticles * 2, cudaMemcpyHostToDevice);
  //cudaMemcpy(CUDA_velocities_full, velocities_full, sizeof(float) * numParticles * 2, cudaMemcpyHostToDevice);
  //cudaMemcpy(CUDA_velocities_half, velocities_half, sizeof(float) * numParticles * 2, cudaMemcpyHostToDevice);
 
  int numBlocks = numParticles/THREADS_PER_BLOCK + 1;
  boundaryCheckKernel<<<numBlocks, THREADS_PER_BLOCK>>>(numParticles, DAMPING, CUDA_positions, CUDA_velocities_full, CUDA_velocities_half);

  //cudaMemcpy(positions, CUDA_positions, sizeof(float) * numParticles * 2, cudaMemcpyDeviceToHost);
  //cudaMemcpy(velocities_full, CUDA_velocities_full, sizeof(float) * numParticles * 2, cudaMemcpyDeviceToHost);
  //cudaMemcpy(velocities_half, CUDA_velocities_half, sizeof(float) * numParticles * 2, cudaMemcpyDeviceToHost);

}

void velocityStep (CurrState* currState, double dt) {
  int numParticles = currState->numParticles;
  float* accelerations = currState->accelerations;
  float* velocities_full = currState->velocities_full;
  float* velocities_half = currState->velocities_half;
  float* positions = currState->positions;
  int numBlocks = currState->numBlocks;

  float* CUDA_velocities_half = currState->CUDA_velocities_half;
  float* CUDA_velocities_full = currState->CUDA_velocities_full;
  float* CUDA_accelerations = currState->CUDA_accelerations;
  float* CUDA_positions = currState->CUDA_positions;

  //cudaMemcpy(currState->CUDA_velocities_half, currState->velocities_half, sizeof(float) * 2 * numParticles, cudaMemcpyHostToDevice);
  //cudaMemcpy(currState->CUDA_positions, currState->positions, sizeof(float) * 2 * numParticles, cudaMemcpyHostToDevice);
  //cudaMemcpy(currState->CUDA_velocities_full, currState->velocities_full, sizeof(float) * 2 * numParticles, cudaMemcpyHostToDevice);
  //cudaMemcpy(currState->CUDA_accelerations, currState->accelerations, sizeof(float) * 2 * numParticles, cudaMemcpyHostToDevice);

  kernel_velocityStep<<<numBlocks, THREADS_PER_BLOCK>>>(CUDA_velocities_half, CUDA_velocities_full, CUDA_accelerations, numParticles, CUDA_positions, dt);

  //cudaMemcpy(currState->velocities_half, currState->CUDA_velocities_half,sizeof(float) * 2 * numParticles, cudaMemcpyDeviceToHost);
  //cudaMemcpy(currState->velocities_full, currState->CUDA_velocities_full, sizeof(float) * 2 * numParticles, cudaMemcpyDeviceToHost);
  //cudaMemcpy(currState->positions, currState->CUDA_positions, sizeof(float)*2*numParticles, cudaMemcpyDeviceToHost);

  boundaryCheck(currState);
}

__global__ void kernel_velocityStart(float* CUDA_velocities_half, float* CUDA_velocities_full, float* CUDA_accelerations, int CUDA_numParticles, float* CUDA_positions, float CUDA_dt) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i >= 2* CUDA_numParticles) {
    return;
  }
  CUDA_velocities_half[i] = CUDA_velocities_full[i] + CUDA_accelerations[i] * CUDA_dt/2;
  CUDA_velocities_full[i] += CUDA_accelerations[i] * CUDA_dt;
  CUDA_positions[i] += CUDA_velocities_half[i] * CUDA_dt;
}


void velocityStart (CurrState* currState, double dt) {
  int numParticles = currState->numParticles;
  float* accelerations = currState->accelerations;
  float* velocities_full = currState->velocities_full;
  float* velocities_half = currState->velocities_half;
  float* positions = currState->positions;

  float* CUDA_velocities_half = currState->CUDA_velocities_half;
  float* CUDA_velocities_full = currState->CUDA_velocities_full;
  float* CUDA_accelerations = currState->CUDA_accelerations;
  float* CUDA_positions = currState->CUDA_positions;
  
  int numBlocks = 2* numParticles/THREADS_PER_BLOCK + 1;
  
  kernel_velocityStart<<<numBlocks, THREADS_PER_BLOCK>>>(CUDA_velocities_half, CUDA_velocities_full, CUDA_accelerations, numParticles, CUDA_positions, dt);
  /*for (unsigned int i = 0; i < 2*numParticles; ++i) {
    velocities_half[i] = velocities_full[i] + accelerations[i] * dt / 2;
  }

  for (unsigned int i = 0; i < 2*numParticles; ++i) {
    velocities_full[i] += accelerations[i] * dt;
  }

  for (unsigned int i = 0; i < 2*numParticles; ++i) {
    positions[i] += velocities_half[i] * dt;
  }*/

  boundaryCheck(currState);
}

typedef int (*initialFluidShape_fun)(float, float);

int cornerBoxInit(float x, float y) {
  return ((x < 0.5) && (y < 0.5));
}

int sphereDropInit(float x, float y) {
  float dx = (x-0.5);
  float dy = (y-0.3);
  float r2 = dx*dx + dy*dy;
  return (r2 < 0.25 * 0.25);
}

CurrState* initialParticlePlacement(Parameters* params, initialFluidShape_fun shapeFun) {
  float size = params->size;
  float adjSize = size/OVERLAP_COEFF;
  int newCount = 0;
  int iterations = 1.0f/adjSize + 1;
  int inRegionCount = 0;
  for (unsigned int i = 0; i < iterations; ++i){
    for (unsigned int j = 0; j < iterations; ++j) {
      float x = i * adjSize;
      float y = j * adjSize;
      inRegionCount += shapeFun(x,y);
    }
  }
  CurrState* currState = allocState(inRegionCount);
  int position = 0;
  for (unsigned int i = 0; i < iterations; ++i) {
    for (unsigned int j = 0; j < iterations; ++j) {
      float x = i * adjSize;
      float y = j * adjSize;
      if (shapeFun(x,y)) {
        currState->positions[GET_X(position)] = x;
        currState->positions[GET_Y(position)] = y;

        currState->velocities_full[GET_X(position)] = 0;
        currState->velocities_full[GET_Y(position)] = 0;
        ++position;
      }
    }
  }
  return currState;
}


void normalizeMasses(CurrState* currState, Parameters* params) {
  currState->mass = 1.f;
  calculateDensity(params, currState);
  float density_ref = params->density_ref;
  float cDensity = 0.f;
  float cDensity_2 = 0.f;
  int numParticles = currState->numParticles;
  for (int i = 0; i < numParticles; ++i) {
    float density = currState->densities[i];
    cDensity += density;
    cDensity_2 += (density * density);
  }
  currState->mass *= (density_ref * cDensity / cDensity_2);
}

CurrState* initParticles(Parameters* params) {
  CurrState* currState = initialParticlePlacement(params, sphereDropInit);
  normalizeMasses(currState, params);
  return currState;
}

void initParams(Parameters* params) {
  params->fileName = "output.out";
  params->numFrames = 400;
  params->stepsPerFrame = 100;
  params->dt = 1e-4;
  params->size = 5e-2;
  params->k = 1e3;
  params->density_ref = 1000;
  params->viscocity = 0.1;
  params->g = 9.8;
}

void errorCheck(CurrState* currState) {
  int numParticles = currState->numParticles;
  float currX;
  float currY;

  for (unsigned int i = 0; i < numParticles; ++i) {
    int xIndex = GET_X(i);
    int yIndex = GET_Y(i);

    currX = currState->positions[xIndex];
    currY = currState->positions[yIndex];

    assert(currX >=0 || currX <=1);
    assert(currY >=0 || currY <=1);
  }
}

int run_main() {

  //printf("main called\n");

  Parameters params;
  initParams(&params);
  CurrState* currState = initParticles(&params);
  string fileName = params.fileName;
  int numFrames = params.numFrames;
  int stepsPerFrame = params.stepsPerFrame;
  float timeStep = params.dt;
  int numParticles = currState->numParticles;

  cudaMemcpy(currState->CUDA_positions, currState->positions, sizeof(float) * numParticles * 2, cudaMemcpyHostToDevice);
  cudaMemcpy(currState->CUDA_accelerations, currState->accelerations, sizeof(float) * numParticles * 2, cudaMemcpyHostToDevice);
  cudaMemcpy(currState->CUDA_densities, currState->densities, sizeof(float) * numParticles, cudaMemcpyHostToDevice);
  cudaMemcpy(currState->CUDA_velocities_full, currState->velocities_full, sizeof(float) * numParticles * 2, cudaMemcpyHostToDevice);
  cudaMemcpy(currState->CUDA_velocities_half, currState->velocities_half, sizeof(float) * numParticles * 2, cudaMemcpyHostToDevice);

  calculateAcceleration(&params, currState);

  velocityStart(currState, timeStep);

  cudaMemcpy(currState->positions, currState->CUDA_positions, sizeof(float) * numParticles * 2, cudaMemcpyDeviceToHost);
  /* Write to file */
  ofstream data_file;
  data_file.open("simulation_data.txt", ios::out);
  data_file << params.size << "\n";
  data_file << numFrames * stepsPerFrame << "\n";
  for (int i=0; i < numParticles; i++) {
    data_file << currState->positions[GET_X(i)] << "\n";
    data_file << currState->positions[GET_Y(i)] << "\n";
  }
  data_file << "DONE WITH AN ITERATION\n";
  //data_file.close();
  /* End of write */

  //errorCheck(currState);

  //iterate through all the frames in the image
  for (unsigned int frame = 1; frame < numFrames; ++frame) {
    //iterate through all the steps per frame
    for (unsigned int j = 0; j < stepsPerFrame; ++j) {
      calculateAcceleration(&params, currState);
      //printf("calculaate acceleration returned on loop %d\n",j);
      velocityStep(currState, timeStep);

      cudaMemcpy(currState->positions, currState->CUDA_positions, sizeof(float) * numParticles * 2, cudaMemcpyDeviceToHost);
        /* Write to file */
        //ofstream data_file;
        //data_file.open("simulation_data.txt", ios::out | ios::app);
        for (int i=0; i < numParticles; i++) {
          data_file << currState->positions[GET_X(i)] << "\n";
          data_file << currState->positions[GET_Y(i)] << "\n";
        }
        data_file << "DONE WITH AN ITERATION\n";

        //errorCheck(currState);
      }
      cout << frame << "\n";
    }

    data_file.close();
    freeState(currState);
    return 0;
}

void parallel() {
  run_main();
}
