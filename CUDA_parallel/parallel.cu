//The main SPH algorithm used in this code was taken from Cornell CS5520

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

#define THREADS_PER_BLOCK 128
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
  int* CUDA_tempInt; //Array of size 1
  float* CUDA_tempFloat; //Array of size 1
  float* CUDA_velocities_half;
  float* CUDA_velocities_full;
  float* CUDA_accelerations;
  float* CUDA_positions;
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
  currState->numBlocks = currState->numParticles / THREADS_PER_BLOCK + 1;

  if (cudaMalloc((void**)&(currState->CUDA_positions), (2*n*(sizeof(float)))) == cudaSuccess) std::cout << "1\n";
  if (cudaMalloc((void**)&(currState->CUDA_tempFloat), sizeof(float)) == cudaSuccess) std::cout << "2\n";
  if (cudaMalloc((void**)&(currState->CUDA_tempInt), sizeof(int)) == cudaSuccess) std::cout << "3\n";
  if (cudaMalloc((void**)&(currState->CUDA_numParticles), sizeof(int)) == cudaSuccess) std::cout << "4\n";
  if (cudaMalloc((void**)&(currState->CUDA_velocities_full),(2*n*(sizeof(float)))) == cudaSuccess) std::cout << "5\n";
  if (cudaMalloc((void**)&(currState->CUDA_velocities_half),(2*n*(sizeof(float)))) == cudaSuccess) std::cout << "6\n";
  if (cudaMalloc((void**)&(currState->CUDA_accelerations),(2*n*(sizeof(float)))) == cudaSuccess) std::cout << "7\n";

  cudaMemset(currState->CUDA_positions, 0, 2*n*sizeof(float));
  cudaMemset(currState->CUDA_velocities_full, 0, 2*n*sizeof(float));
  cudaMemset(currState->CUDA_velocities_half, 0, 2*n*sizeof(float));
  cudaMemset(currState->CUDA_accelerations, 0, 2*n*sizeof(float));
  cudaMemset(currState->CUDA_tempFloat, 0, sizeof(float));
  cudaMemset(currState->CUDA_tempInt, 0, sizeof(int));
  cudaMemset(currState->CUDA_numParticles, 0, sizeof(int));

  cout << &(currState->CUDA_velocities_full) << " " << &(currState->velocities_full) << "\n";
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
  cudaFree(currState->CUDA_tempFloat);
  cudaFree(currState->CUDA_tempInt);
  cudaFree(currState->CUDA_numParticles);
  cudaFree(currState->CUDA_velocities_full);
  cudaFree(currState->CUDA_velocities_half);
  cudaFree(currState->CUDA_accelerations);
}

__global__ void kernel_velocityStep(float* CUDA_velocities_half, float* CUDA_velocities_full, float* CUDA_accelerations, int* CUDA_numParticles, float* CUDA_positions, float* CUDA_dt) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i == 4) printf("DT is %d and numParticles is %d, velocities half is %d, velocities full is %d, position is %d, acceleration is %d\n", CUDA_dt[0], CUDA_numParticles[0], CUDA_velocities_half[i], CUDA_velocities_full[i], CUDA_positions[i], CUDA_accelerations[i]);
  if (i >= 2* CUDA_numParticles[0]) return;
  CUDA_velocities_half[i] += CUDA_accelerations[i] * CUDA_dt[0];
  CUDA_velocities_full[i] = CUDA_velocities_half[i] + CUDA_accelerations[i] * CUDA_dt[0] / 2;
  CUDA_positions[i] += CUDA_velocities_half[i] * CUDA_dt[0];  
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

  for (unsigned int i = 0; i < numParticles; ++i) {
    densities[i] += innerConstant;
    for (unsigned int j = i + 1; j < numParticles; ++j) {
      float dx = positions[GET_X(i)] - positions[GET_X(j)];
      float dy = positions[GET_Y(i)] - positions[GET_Y(j)];
      float r_2 = dx * dx + dy * dy;
      float z = size_2 - r_2;
      float z_3 = z * z * z;
      if (z > 0) {
        float densities_ij = outerConstant * z_3;
        densities[i] += densities_ij;
        densities[j] += densities_ij;
      }
    }
  }
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
  //printf("initialization complete in calculateAcceleration\n");
  calculateDensity(params, currState);
  //printf("calculateDensity Returns\n");
  for (unsigned int i = 0; i < numParticles; ++i) {
    accelerations[GET_X(i)] = 0;
    accelerations[GET_Y(i)] = -g;
  }
  float C0 = mass / (M_PI * size_4);
  float Cp = 15 * k;
  float Cv = -40 * viscocity;

  for (unsigned int i = 0; i < numParticles; ++i) {
    const float currDensity_i = densities[i];
    for (unsigned int j = i+1; j < numParticles; ++j) {
      float dx = positions[GET_X(i)] - positions[GET_X(j)];
      float dy = positions[GET_Y(i)] - positions[GET_Y(j)];
      float r_2 = dx * dx + dy * dy;
      if (r_2 < size_2) {
        const float currDensity_j = densities[j];
        float q = sqrt(r_2)/size;
        float u = 1-q;
        float w0 = C0 * u/(currDensity_j * currDensity_i);
        float wp = w0 * Cp * (currDensity_i + currDensity_j - 2 * density_ref) * u/q;
        float wv = w0 * Cv;
        float dvx = velocities[GET_X(i)] - velocities[GET_X(j)];
        float dvy = velocities[GET_Y(i)] - velocities[GET_Y(j)];
        accelerations[GET_X(i)] += (wp * dx + wv * dvx);
        accelerations[GET_Y(i)] += (wp * dy + wv * dvy);
        accelerations[GET_X(j)] -= (wp * dx + wv * dvx);
        accelerations[GET_Y(j)] -= (wp * dy + wv * dvy);
      }
    }
  }
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

void boundaryCheck(CurrState* currState) {
  const float XMIN = 0.0;
  const float YMIN = 0.0;
  const float XMAX = 1.0;
  const float YMAX = 1.0;

  int numParticles = currState->numParticles;
  float* velocities_full = currState->velocities_full;
  float* velocities_half = currState->velocities_half;
  float* positions = currState->positions;

  for (unsigned int i = 0; i < numParticles; ++i, positions+=2,  velocities_full+=2, velocities_half+=2) {

    if (positions[0] < XMIN) {
      reflect (0, XMIN, positions, velocities_full, velocities_half);
    }

    if (positions[0] > XMAX) {
      reflect (0, XMAX, positions, velocities_full, velocities_half);
    }

    if (positions[1] < YMIN) {
      reflect (1, YMIN, positions, velocities_full, velocities_half);
    }

    if (positions[1] > YMAX) {
      reflect (1, YMAX, positions, velocities_full, velocities_half);
    }
  }
}

void velocityStep (CurrState* currState, double dt) {
  int numParticles = currState->numParticles;
  float* accelerations = currState->accelerations;
  float* velocities_full = currState->velocities_full;
  float* velocities_half = currState->velocities_half;
  float* positions = currState->positions;
  int numBlocks = currState->numBlocks;

  float* CUDA_tempFloat = currState->CUDA_tempFloat;
  int* CUDA_numParticles = currState->CUDA_numParticles;
  float* CUDA_velocities_half = currState->CUDA_velocities_half;
  float* CUDA_velocities_full = currState->CUDA_velocities_full;
  float* CUDA_accelerations = currState->CUDA_accelerations;
  float* CUDA_positions = currState->CUDA_positions;

  std::cout << &(currState->CUDA_velocities_full) << " " << &(currState->velocities_full) << "\n";
  std::cout<< currState->velocities_half[4] <<" " << currState->velocities_full[4] << " " << currState->accelerations[4] << "\n";

  if(cudaMemcpy(currState->CUDA_velocities_half, currState->velocities_half, sizeof(float) * 2 * numParticles, cudaMemcpyHostToDevice) == cudaSuccess) std::cout << "YES1\n";
  if(cudaMemcpy(currState->CUDA_positions, currState->positions, sizeof(float) * 2 * numParticles, cudaMemcpyHostToDevice) == cudaSuccess) std::cout << "YES2\n";
  if(cudaMemcpy(currState->CUDA_velocities_full, currState->velocities_full, sizeof(float) * 2 * numParticles, cudaMemcpyHostToDevice) == cudaSuccess) std::cout << "YES3\n";
  if(cudaMemcpy(currState->CUDA_accelerations, currState->accelerations, sizeof(float) * 2 * numParticles, cudaMemcpyHostToDevice) == cudaSuccess) std::cout << "YES4\n";
  if(cudaMemcpy(currState->CUDA_numParticles, &(currState->numParticles), sizeof(int), cudaMemcpyHostToDevice) == cudaSuccess) std::cout << "YES5\n";
  if(cudaMemcpy(currState->CUDA_tempFloat, &dt, sizeof(float), cudaMemcpyHostToDevice) == cudaSuccess) std::cout << "YES6\n";

  kernel_velocityStep<<<numBlocks, THREADS_PER_BLOCK>>>(CUDA_velocities_half, CUDA_velocities_full, CUDA_accelerations, CUDA_numParticles, CUDA_positions, CUDA_tempFloat);

/*
  for (unsigned int i = 0; i < 2*numParticles; ++i) {
    velocities_half[i] += accelerations[i] * dt;
  }

  for (unsigned int i = 0; i < 2*numParticles; ++i) {
    velocities_full[i] = velocities_half[i] + accelerations[i] * dt / 2;
  }

  for (unsigned int i = 0; i < 2*numParticles; ++i) {
    positions[i] += velocities_half[i] * dt;
  }
  */

  float tempdt = -1;
  float tempnum = -1;
  cudaMemcpy(currState->velocities_half, currState->CUDA_velocities_half,sizeof(float) * 2 * numParticles, cudaMemcpyDeviceToHost);
  cudaMemcpy(currState->velocities_full, currState->CUDA_velocities_full, sizeof(float) * 2 * numParticles, cudaMemcpyDeviceToHost);
  //cudaMemcpy(currState->accelerations, currState->CUDA_accelerations, sizeof(float) * 2 * numParticles, cudaMemcpyDeviceToHost);
  cudaMemcpy(currState->positions, currState->CUDA_positions, sizeof(float)*2*numParticles, cudaMemcpyDeviceToHost);
  cudaMemcpy(&tempdt, currState->CUDA_tempFloat, sizeof(float), cudaMemcpyDeviceToHost);
  cudaMemcpy(&tempnum, currState->CUDA_numParticles, sizeof(float), cudaMemcpyDeviceToHost);

  //std::cout<< currState->velocities_half[4] <<" " << currState ->velocities_full[4] << " " << currState->accelerations[4] << " " << tempdt << " " << *(&dt) << " " << numParticles << " " << tempnum << "\n";
  boundaryCheck(currState);
}

void velocityStart (CurrState* currState, double dt) {
  int numParticles = currState->numParticles;
  float* accelerations = currState->accelerations;
  float* velocities_full = currState->velocities_full;
  float* velocities_half = currState->velocities_half;
  float* positions = currState->positions;

  for (unsigned int i = 0; i < 2*numParticles; ++i) {
    velocities_half[i] = velocities_full[i] + accelerations[i] * dt / 2;
  }

  for (unsigned int i = 0; i < 2*numParticles; ++i) {
    velocities_full[i] += accelerations[i] * dt;
  }

  for (unsigned int i = 0; i < 2*numParticles; ++i) {
    positions[i] += velocities_half[i] * dt;
  }

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
  printf("inRegionCount is %d\n", inRegionCount);

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

  calculateAcceleration(&params, currState);
  velocityStart(currState, timeStep);

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

  errorCheck(currState);

  //iterate through all the frames in the image
  for (unsigned int frame = 1; frame < numFrames; ++frame) {
    //iterate through all the steps per frame
    for (unsigned int j = 0; j < stepsPerFrame; ++j) {
      calculateAcceleration(&params, currState);
      //printf("calculaate acceleration returned on loop %d\n",j);
      velocityStep(currState, timeStep);

        /* Write to file */
        //ofstream data_file;
        //data_file.open("simulation_data.txt", ios::out | ios::app);
        for (int i=0; i < numParticles; i++) {
          data_file << currState->positions[GET_X(i)] << "\n";
          data_file << currState->positions[GET_Y(i)] << "\n";
        }
        data_file << "DONE WITH AN ITERATION\n";

        errorCheck(currState);
      }
    }

    data_file.close();
    freeState(currState);
    return 0;
}

void parallel() {
  run_main();
}
