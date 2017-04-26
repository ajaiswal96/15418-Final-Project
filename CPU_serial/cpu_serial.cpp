using namespace std;

#define VISCOCITY 1
#define DIFFRATE 1
#define TIMESTEP 1
#define TOTAL_COLS 1000
#define TOTAL_ROWS 2000

struct fluidParticle {
  float density;
  float vxf;
  float vyf;
  float vx0;
  float vy0;
  float temp;
};
typedef struct fluidParticle fluidParticle;

struct fluidGrid {
  int numCols;
  int numRows;
  float timestep; //the length of a timestep
  float diffRate; //the rate of diffusion
  float viscosity; //thickness of fluid

  fluidParticle **allParticles; //2D-array of all the particles on the canvas
};
typedef struct fluidGrid fluidGrid;




fluidGrid* initCanvas(){
  //allocate space for the canvas
  fluidGrid* canvas = malloc(sizeof(fluidGrid)); //allocate space for the canvas

  //initialize constants in canvas
  canvas->numCols = TOTAL_COLS;
  canvas->numRows = TOTAL_ROWS;
  canvas->timestep = TIMESTEP;
  canvas->diffRate = DIFFRATE;
  canvas->viscosity = VISCOCITY;

  /*
  canvas->allParticles = (fluidParticle***) calloc(SIZE, sizeof(fluidParticle**))
  for (unsigned int i = 0; i < SIZE; i++) {
    canvas->allParticles[i] = (fluidParticle*) calloc(SIZE,sizeof(fluidParticle*));
  }
  */

  //allocate 2d array of all the aprticles
  canvas->allParticles = new fluidParticle* [canvas->numRows];
  for (unsigned int i = 0; i < canvas->numRows; i++) {
    canvas->allParticles[i] = new fluidParticle [canvas->numCols];
  }

  //initialize the values for each particle in the grid
  for (unsigned int row = 0; i < canvas->numRows; row++) {
    for (unsigned int col = 0; j < canvas->numCols; col++) {
      canvas->allParticles[row][col].density = 0;
      canvas->allParticles[row][col].vxf = 0;
      canvas->allParticles[row][col].vyf = 0;
      canvas->allParticles[row][col].vx0 = 0;
      canvas->allParticles[row][col].vy0 = 0;
      canvas->allParticles[row][col].temp = 0;
    }
  }
  return canvas;
}

void freeCanvas(fluidGrid *canvas){
  for (unsigned int row = 0; i < canvas->numRows; row++) {
    delete [] canvas->allParticles[row];
  }
  delete [] canvas->allParticles;
  free(canvas);
}

//This is just a helper function to change the velocity of a particle in the grid
static void addVelocity(fluidGrid *canvas, int row, int col, float drow, float dcol) {
  canvas->allParticles[row][col].vxf += dcol;
  canvas->allParticles[row][col].vyf += drow;
}

//This is just a helper function that adds density to a particle in the grid
static void addDensity(fluidGrid *canvas, int row, int col, float dden) {
  canvas->allParticles[row][col].density += dden
}

//this is an advect function that is going to be called for each element
static void advect(int b, float* d, float* d0, float* vx0, float* vy0, float dt, int rows, int cols, int currRow, int currCol, fluidGrid* canvas) {
  float dtx = dt*(cols-2);
  float dty = dt*(rows-2);
  float cRow = (float)currRow;
  float cCol = (float)currCol;

  float tmp1 = dtx * (*vx0);
  float tmp2 = dty * (*vy0);
  float x = cCol - tmp1;
  float y = cRow - tmp2;

  if (x < 0.5f) {
    x = 0.5f;
  }

  if (x > (float)cols + 0.5f) {
    x = (float)cols + 0.5f;
  }

  float c0 = floorf(x);
  float c1 = c0 + 1.0f;

  if (y < 0.5f) {
    y = 0.5f;
  }

  if (y > (float)rows + 0.5f) {
    y = (float)rows + 0.5f;
  }

  float r0 = floorf(y);
  float r1 = r0 + 1.0f;

  float s1 = x - c0;
  float s0 = 1.0f - s1;

  float t1 = y-r0;
  float t0 = 1.0f - t1;

  int r0_i = (int)r0;
  int r1_i = (int)r1;
  int s0_i = (int)s0;
  int s1_i = (int)s1;
  // s refers to col and t refers to row
  *d =
  c0 * r0 * (canvas -> allParticles[]) +
  c0 * r1 * (shit) +
  c1 * r0 * (shit) +
  c1 * r1 * (shit);

//now we need to update;
}

static void nextStep(fluidGrid *canvas) {
  float dt = canvas->timestep;
  float viscosity = canvas->viscocity;
  int rows = canvas->numRows;
  int cols = canvas->numCols;
  //diffuse for vx (column)
  if (rows <=2 || cols <=2) {
    printf("ERROR: not enough rows or cols\n");
    return;
  }
  float diffuse_vxa = dt * viscocity * (rows-2) * (cols-2);
  int iter_diff_vxa =

}

void navierStokes(int a, int b, int c, fluidGrid* canvas, int iters, int maxRows, int maxCols) {
  for (unsigned int i = 0; i < iters; i++) {
    for (unsigned int row = 0; row < maxRows-1; row++) {
      for (unsigned int col = 0; col < maxCols -1; col++) {
        canvas->allParticles[row][col].
      }
    }
  }
}
