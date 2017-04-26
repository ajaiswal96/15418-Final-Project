### Team

Anubhav Jaiswal  (ajaiswal)

Rohit Pillai  (rrpillai)

## Checkpoint Update

We have updated our schedule to reflect where we currently are. 

### Work Completed Thus Far 
The majority of our time was spent on understanding the Navier Stokes algorithm and trying to figure out the best way to manipulate and engineer it to fit our needs. We did not have much difficulty realizing how to write a sequential algorithm, though we focused on having functions operate on individual elements rather than arrays so we can parallelize it easily. Since each element needs to know information on all the elements around it, it would be trivial to have functions operate on the whole array. It is much more difficult to have per-element functions.

We have mostly completed a sequential implementation of the fluid simulator. The biggest challenge was understanding the Navier Stokes equation and using data structures that made it both simple to display a visual output, and run all the steps of the Navier Stokes equation. There are currently bugs in the sequential code that prevent a correct output from being generated. The bugs have to do with the way we broke down the Navier Stokes equation and we are looking into ways in which we can rewrite the code. 

We are also deciding on a way to display our output image. Initially we thought about just using openGL and making calls within the program to display images at different stages. We realized this probably wouldn't be the best way to display an animation as there are a lot of overheads invovled. We are now deciding between using the SimpleImage library and the Box2D image library. SimpleImage is much lighter and easier to use though may not have all the functionality we need.

### How We are Doing Overall



### Summary

We are going to implement a parallel version of a fluid simulator to run on Nvidia GTX 1080 GPU's on the Gates machines and compare that to a parallel implementation on a CPU, using a sequential CPU implementation as baseline. We will be using the Navier-Stokes equation to run the CFD calculations for the simulation.

### Background

We are implementing a parallel computational fluid dynamics (CFD) simulator to run on GPU's. For the sake of simplicity, we are assuming the fluids to be homogeneous and incompressible. We will be representing the fluid as a simulation grid, which is indexed in by their position, and stores a time and velocity (since velocity and location change based on time). In our simulation we need to account for four things: advection, pressure, diffusion, and external forces. Advection refers to the movement of particles from one region to another. Pressure accounts to the fact that particles closer to a force will move with greater velocity than those further away. Diffusion accounts for viscocity of fluids. External forces refer to the actual forces that act on the fluid. 

We plan to take advantage of the fact that GPU's are optimized to render textures, where a texture is just a 2D vector storing multiple values. We can thus store our simulation grids as these textures. Instead of having each pixel in the texture map correspond to a particle, we are going to have a grid of particles and store each grid component as a pixel texture. The inherent structure of CFD simulations thus makes it best suited to be run through CUDA kernels on the GPU. This way instead of looping through a 2D array and updating values, we can update blocks in parallel. 

<img src="15418-Final-Project/assets/17857867_1329636113793992_1520634717_n.png" alt="Navier-Stokes" width="600" height="200">

The Navier-Stokes equation for incompressible homogeneous fluids forms the basis of a lot of CFD and is used to describe the motion of a fluid. Given the initial conditions of the fluid (which could be parameters in our implementation), we can solve this equation at each time step to find the state of the fluid at the next step. There are 4 terms in the equation that correspond to the 4 components (advection, pressure, diffusion and external forces) as described above. 
 
### Challenges 
One of the challenges we have to deal with is a varied density of particles in a grid cell due to advection. For example if we have a force pushing to the right applied to particles on the left, The grid cells on the left will have fewer particles in them as time moves on, and the grid cells on the right will have more. This will lead to divergent execution due to the workload imbalance, hindering peak performance from a parallel implementation. Another issue involves dealing with dependencies because not all the steps will be able to occur in parallel. One run through the fluid velocity grid is necessary to compute the effects of advection; another for the affects of pressure on the grid that accounts for advection; another for the affects of diffusion on the previous grid; and finally one to compute the affects of the external forces on the previous grid. Finding out how to compute these vector grids in parallel and sequentially imposing them on each other will not be simple due to inherent dependencies in their equations. Making efficient use of shared memory within CUDA kernels is crucial for speedup, thus we will have to find an efficient way to map particles onto CUDA blocks. 

### Resources
- We plan to use the the GTX 1080 GPU’s on the GHC machines that we used for assignment 2, and thus all our code will most likely be in CUDA and C++.
- We may also need libraries to visualize the fluid simulation and display the liquid flow as time progresses. 
- All of our preliminary research has been based off a textbook published by NVIDIA that describes how to do fluid simulation. (http://http.developer.nvidia.com/GPUGems/gpugems_ch38.html) 
- We don’t have an existing code base, and so we will be starting from scratch
- For the CPU implementation, we plan to use one of the Gates machines as well since they have 8 core 3.2 GHz Intel Core i7 processors, which is the fastest processor that we have access to.

### Goals 
- We definitely plan to achieve a parallelized version of the fluid simulator using the Navier Stokes equation, on the GTX 1080 GPU. We also want to make a parallel version of this to run on the CPU so that we can see how much speedup we get from the GPU version over this CPU version.
- We also want to make a sequential version to run on the CPU to get baseline results with no optimizations.
- We definitely also want to have some form of GUI to visualize the outputs of our fluid simulation. 
- Something we hope to achieve is simulating free surface boundaries between 2 different fluids (for example, between air and water). Free surface boundaries are just the points of contact of the 2 fluids. This is different than just a single fluid simulation because when the 2 fluids interact, we will have to take into account their different properties, which will require a lot more computation.
- Another very far fetched goal would be to convert this 2-D fluid simulation to a 3-D fluid simulation, with all the vectors and equations being in 3-D.
- With respect to our demo, we hope to have a visual representation of our 2-D fluid simulation that shows a fluid varying as time progresses.
- We also plan to have speedup graphs that show the speedup that we get from our parallel CUDA version over both CPU versions that we will implement as well. 
- We hope to see significant speedup (>4x) going from the parallel CPU version to the parallel GPU version on larger input sizes as the overhead of parallelism will hide benefits on smaller sized inputs. 
- We also hope to see both the parallel CPU and GPU versions greatly outperform the baseline sequential CPU implementation.


### Platform Choice 
- We will be working in C++ and using libraries for the GUI. We would also like to use math libraries to have access to data structures to efficiently store our velocity, time, and position values in the velocity grid. 
- We choose to use C++ because there are a lot of libraries we can use and it works with CUDA. 
- We are using the Nvidia GPU's and CUDA because the structure of fluid simulation makes it suitable to be run on CUDA thread blocks. 
- We will need to use openMP to parallelize the CPU implementation to write the parallel CPU implementation 

### Schedule
**April 10:** Finish Proposal 

**April 15:** Finish Serial Implementation of Fluid Simulation for the CPU

**April 17:** Finish setting up GUI to see outputs of the simulation 

**April 19:** Finish parallel implementation of Fluid Simulation for the GPU

**April 25:** Have first iteration of parallel GPU implementation working  

**May 1:** Improve the parallel GPU implementation to obtain peak performance

**May 6:** If substantial speedup achieved, work on free surface boundaries and 3d implementation, otherwise try to acheive greater speedup

**May 8:** Analyze the performance outputs and write final report
