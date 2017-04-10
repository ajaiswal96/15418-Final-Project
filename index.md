## Proposal

### Team

Anubhav Jaiswal  (ajaiswal)

Rohit Pillai  (rrpillai)

### Summary

We are going to implement a parallel version of a fluid simulator to run on Nvidia GTX 1080 GPU's on the Gates machines and compare that to a sequential implementation on a CPU.

### Background

We are implementing a parallel computational fluid dynamics simulator to run on GPU's. For the sake of simplicity, we are assuming the fluids to be homogeneous and incompressible. We plan to take advantage of the fact that GPU's are optimized to render textures, and a texture is just a 2D vector storing multiple values. We can thus store our simulation grids as these textures. Instead of having each pixel in the texture map correspond to a particle, we are going to have a grid of particles and store each grid componenet as a pixel texture. The inherent structure of CFD simulations thus makes it best suited to be run through CUDA kernels on the GPU. 

We will be representing the fluid as a simulation grid, which is indexed in by their position, and stores a time and velocity (since velocity and location change based on time). In our simulation we need to account for four things: advection, pressure, diffusion, and external forces. Advection refers to the movement of particles from one region to another. Pressure accounts to the fact that particles closer to a force will move with greater velocity than those further away. Diffusion accounts for viscocity of fluids. External forces refer to the actual forces that act on the fluid. 

<img src="https://github.com/ajaiswal96/15418-Final-Project/blob/master/assets/17857867_1329636113793992_1520634717_n.png" alt="Navier-Stokes Equation" width="600" height="200">

The Navier-Stokes equation for incompressible homogeneous fluids forms the basis of a lot of CFD and is used to describe the motion of a fluid. Given the initial conditions of the fluid (which could be parameters in our implementation), we can solve this equation at each time step to find the state of the fluid at the next step. There are 4 terms in the equation that correspond to the 4 components (advection, pressure, diffusion and external forces) as described above. 
 
### Challenges 
One of the challenges we have to deal with is a varied density of particles in a grid cell due to advection. For example if we have a force pushing to the right applied to particles on the left, The grid cells on the left will have fewer particles in them as time moves on, and the grid cells on the right will have more. This will lead to divergent execution due to the workload imbalance, hindering peak performance from a parallel implementation. Another issue involves dealing with dependencies because not all the steps will be able to occur in parallel. One run through the fluid velocity grid is necessary to compute the effects of advection; another for the affects of pressure on the grid that accounts for advection; another for the affects of diffusion on the previous grid; and finally one to compute the affects of the external forces on the previous grid. Finding out how to compute these vector grids in parallel and sequentially imposing them on each other will not be simple due to inherent dependencies in their equations. 

### Resources
- We plan to use the the GTX 1080 GPU’s on the GHC machines that we used for assignment 2, and thus all our code will most likely be in CUDA and C++.
- We may also need libraries to visualize the fluid simulation and display the liquid flow as time progresses. 
- All of our preliminary research has been based off a textbook published by NVIDIA that describes how to do fluid simulation. (http://http.developer.nvidia.com/GPUGems/gpugems_ch38.html) 
- We don’t have an existing code base, and so we will be starting from scratch
- For the sequential implementation, we plan to use one of the Gates machines as well since they have 8 core 3.2 GHz Intel Core i7 processors, which is the fastest processor that we have access to.

### Goals 
- We definitely plan to achieve a parallelized version of the fluid simulator using the Navier Stokes equation, on the GTX 1080 GPU. We also want to make a sequential version of this so that we can see how much speedup we get from the parallel version over this sequential version.
- Something we hope to achieve is simulating free surface boundaries between 2 different fluids (for example, between air and water). Free surface boundaries are just the points of contact of the 2 fluids. This is different than just a single fluid simulation because when the 2 fluids interact, we will have to take into account their different properties, which will require a lot more computation.
- Another very far fetched goal would be to convert this 2-D fluid simulation to a 3-D fluid simulation, with all the vectors and equations being in 3-D.
- With respect to our demo, we hope to have a visual representation of our 2-D fluid simulation that shows a fluid varying as time progresses.
- We also plan to have speedup graphs that show the speedup that we get from our parallel CUDA version over the sequential version that we will implement as well. 

### Platform Choice 
- We will be working in C++ and using libraries for the GUI. We would also like to use math libraries to have access to data structures to efficiently store our velocity, time, and position values in the velocity grid. 

### Schedule
**April 10:** Finish Proposal 

**April 17:** Finish Serial Implementation of Fluid Simulation for the CPU

**April 19:** Finish setting up GUI to see outputs of the simulation

**April 25:** Have first iteration of parallel GPU implementation working  

**May 1:** Improve the parallel GPU implementation to obtain peak performance

**May 6:** If substantial speedup achieved, work on free surface boundaries and 3d implementation, otherwise try to acheive greater speedup

**May 8:** Analyze the performance outputs and write final report
