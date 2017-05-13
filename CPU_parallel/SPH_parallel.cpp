//The main SPH algorithm used in this code was taken from Cornell CS5520

using namespace std;

#include "parallel.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <float.h>
#include <string.h>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <omp.h>

#define GET_X(x) (2*x)
#define GET_Y(y) (2*y + 1)

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
} CurrState;


CurrState* allocState(int n) {
  CurrState* currState = (CurrState*)calloc(1,sizeof(CurrState));
  if(currState == NULL) {
      std::cout << "WTF1\n";
  }
  currState->densities = (float*)calloc(n, sizeof(float));
  if(currState->densities == NULL) {
      std::cout << "WTF2\n";
  }
  currState->positions = (float*)calloc(2*n, sizeof(float));
  if(currState->positions == NULL) {
      std::cout << "WTF3\n";
  }
  currState->velocities_full = (float*)calloc(2*n, sizeof(float));
  if(currState->velocities_full == NULL) {
      std::cout << "WTF4\n";
  }
  currState->velocities_half = (float*)calloc(2*n, sizeof(float));
  if(currState->velocities_half == NULL) {
      std::cout << "WTF5\n";
  }
  currState->accelerations = (float*)calloc(2*n, sizeof(float));
  if(currState->accelerations == NULL) {
      std::cout << "WTF6\n";
  }
  currState->numParticles = n;
  /*
  for(int i = 0; i < n; i++){
    float cX = currState->positions[2*i];
    float cY = currState->positions[2*i + 1];
    printf("%f, %f\n",cX, cY);
  }
  printf("DONE\n");
  */
}

void freeState(CurrState* currState) {
  free(currState->densities);
  free(currState->positions);
  free(currState->velocities_half);
  free(currState->velocities_full);
  free(currState->accelerations);
  free(currState);
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

#pragma omp parallel for 
  for (unsigned int i = 0; i < numParticles; ++i) {
    densities[i] += innerConstant;
    for (unsigned int j = 0; j < numParticles; ++j) {
        if (j==i) continue;
        float dx = positions[GET_X(i)] - positions[GET_X(j)];
        float dy = positions[GET_Y(i)] - positions[GET_Y(j)];
        float r_2 = dx * dx + dy * dy;
        float z = size_2 - r_2;
        float z_3 = z * z * z;
        if (z > 0) {
            float densities_ij = outerConstant * z_3;
            densities[i] += densities_ij;
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
#pragma omp parallel for
  for (unsigned int i = 0; i < numParticles; ++i) {
    accelerations[GET_X(i)] = 0;
    accelerations[GET_Y(i)] = -g;
  }
  float C0 = mass / (M_PI * size_4);
  float Cp = 15 * k;
  float Cv = -40 * viscocity;

#pragma omp parallel for 
  for (unsigned int i = 0; i < numParticles; ++i) {
    const float currDensity_i = densities[i];
    for (unsigned int j = 0; j < numParticles; ++j) {
        if (j!=i) {
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
                if (i < j){
                    accelerations[GET_X(i)] += (wp * dx + wv * dvx);
                    accelerations[GET_Y(i)] += (wp * dy + wv * dvy);
                }
                else {
                    accelerations[GET_X(i)] -= (wp * dx + wv * dvx);
                    accelerations[GET_Y(i)] -= (wp * dy + wv * dvy);
                }
            }
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

  #pragma omp parallel for
  for (unsigned int i = 0; i < 2*numParticles; ++i) {
    velocities_half[i] += accelerations[i] * dt;
    velocities_full[i] = velocities_half[i] + accelerations[i] * dt / 2;
    positions[i] += velocities_half[i] * dt;
  }

  boundaryCheck(currState);
}

void velocityStart (CurrState* currState, double dt) {
  int numParticles = currState->numParticles;
  float* accelerations = currState->accelerations;
  float* velocities_full = currState->velocities_full;
  float* velocities_half = currState->velocities_half;
  float* positions = currState->positions;

  #pragma omp parallel for 
  for (unsigned int i = 0; i < 2*numParticles; ++i) {
    velocities_half[i] = velocities_full[i] + accelerations[i] * dt / 2;
    velocities_full[i] += accelerations[i] * dt;
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


  CurrState* currState = allocState(inRegionCount);
  int position = 0;
  for (unsigned int i = 0; i < iterations; ++i) {
    for (unsigned int j = 0; j < iterations; ++j) {
      float x = i * adjSize;
      float y = j * adjSize;
//      std::cout << shapeFun(x,y) << "\n";
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
  params->size = 5e-3;
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

    //check that x and y are in bounds and stop program execution if not
    //printf("currX is %f\n",currX);
    //printf("currY is %f\n",currY);

    assert(currX >=0 || currX <=1);
    assert(currY >=0 || currY <=1);
  }
}

int run_main() {

  Parameters params;
  initParams(&params);
  CurrState* currState = initParticles(&params);
  string fileName = params.fileName;
  int numFrames = params.numFrames;
  int stepsPerFrame = params.stepsPerFrame;
  float timeStep = params.dt;
  int numParticles = currState->numParticles;

  std::cout << "NUM PARTICLES IS " << numParticles << "\n";
  calculateAcceleration(&params, currState);
  velocityStart(currState, timeStep);

  /* Write to file */
  /*ofstream data_file;
  data_file.open("simulation_data.txt", ios::out);
  data_file << params.size << "\n";
  data_file << numFrames * stepsPerFrame << "\n";
  for (int i=0; i < numParticles; i++) {
    data_file << currState->positions[GET_X(i)] << "\n";
    data_file << currState->positions[GET_Y(i)] << "\n";
  }
  data_file << "DONE WITH AN ITERATION\n";
  data_file.close();*/
  /* End of write */

  errorCheck(currState);

  clock_t time_start = clock();

  //iterate through all the frames in the image
  for (unsigned int frame = 1; frame < numFrames; ++frame) {
    //iterate through all the steps per frame
    for (unsigned int j = 0; j < stepsPerFrame; ++j) {
      calculateAcceleration(&params, currState);
      //printf("calculaate acceleration returned on loop %d\n",j);
      velocityStep(currState, timeStep);

        /* Write to file */
        /*ofstream data_file;
        data_file.open("simulation_data.txt", ios::out | ios::app);
        for (int i=0; i < numParticles; i++) {
            if (i == 0) {
                printf("x: %f, %f\n", currState->positions[GET_X(i)], currState->positions[GET_Y(i)]);
                printf("d: %f, %f\n", currState->densities[GET_X(i)], currState->densities[GET_Y(i)]);
                printf("a: %f, %f\n", currState->accelerations[GET_X(i)], currState->accelerations[GET_Y(i)]);
                printf("v: %f, %f\n", currState->velocities_full[GET_X(i)], currState->velocities_full[GET_Y(i)]);
            }
    	  data_file << currState->positions[GET_X(i)] << "\n";
	  data_file << currState->positions[GET_Y(i)] << "\n";
  	}
  	data_file << "DONE WITH AN ITERATION\n";
  	data_file.close();*/
  	/* End of write */

      //printf("379\n");
      errorCheck(currState);
      //printf("iterating through frame %d, step %d\n", frame, j);
    }
    std::cout << frame << "\n";
  }
  clock_t time_end = clock();
  std::cout << "TIME IS " << ((float)time_end-time_start)/CLOCKS_PER_SEC << "\n";
  freeState(currState);
}

void parallel() {
  run_main();
}
