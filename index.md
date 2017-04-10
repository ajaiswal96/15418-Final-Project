## Proposal

### Team

Anubhav Jaiswal  (ajaiswal)

Rohit Pillai  (rrpillai)

### Summary

We are going to implement a parallel version of a fluid simulator to run on Nvidia GTX 1080 GPU's on the Gates machines and compare that to a sequential implementation on a CPU.

### Background

We are implementing a parallel computational fluid dynamics simulator to run on GPU's. For the sake of simplicity, we are assuming the fluids to be homogeneous and incompressible. We plan to take advantage of the fact that GPU's are optimized to render textures, and a texture is just a 2D vector storing multiple values. We can thus store our simulation grids as these textures. Instead of having each pixel in the texture map correspond to a particle, we are going to have a grid of particles and store each grid componenet as a pixel texture. The inherent structure of CFD simulations thus makes it best suited to be run through CUDA kernels on the GPU. 

We will be representing the fluid as a simulation grid, which is indexed in by their position, and stores a time and velocity (since velocity and location change based on time). In our simulation we need to account for four things: advection, pressure, diffusion, and external forces. Advection refers to the movement of particles from one region to another. Pressure accounts to the fact that particles closer to a force will move with greater velocity than those further away. Diffusion accounts for viscocity of fluids. External forces refer to the actual forces that act on the fluid. 

<img src="https://github.com/ajaiswal96/15418-Final-Project/blob/master/assets/Screen%20Shot%202017-04-10%20at%202.56.03%20PM.png" alt="Navier-Stokes" width="600" height="200">

The Navier-Stokes equation for incompressible homogeneous fluids forms the basis of a lot of CFD and is used to describe the motion of a fluid. Given the initial conditions of the fluid (which could be parameters in our implementation), we can solve this equation at each time step to find the state of the fluid at the next step. There are 4 terms in the equation that correspond to the 4 components (advection, pressure, diffusion and external forces) as described above. 
 
### Challenges 
One of the challenges we have to deal with is a varied density of particles in a grid cell due to advection. For example if we have a force pushing to the right applied to particles on the left, The grid cells on the left will have fewer particles in them as time moves on, and the grid cells on the right will have more. This will create a workload imbalance, hindering peak performance from a parallel implementation. Another issue involves dealing with dependencies because not all the steps will be able to occur in parallel. One run through the fluid velocity grid is necessary to compute the effects of advection; another for the affects of pressure on the grid that accounts for advection, another for the grid that 

### Resources
- We plan to use the the GTX 1080 GPU’s on the GHC machines that we used for assignment 2, and thus all our code will most likely be in CUDA.
- We may also need libraries to visualize the fluid simulation and display the liquid flow as time progresses. 
- All of our preliminary research has been based off a textbook published by NVIDIA that describes how to do fluid simulation. (http://http.developer.nvidia.com/GPUGems/gpugems_ch38.html) 
- We don’t have an existing code base, and so we will be starting from scratch

```markdown
Syntax highlighted code block

# Summary
## Header 2
### Header 3

- Bulleted
- List

1. Numbered
2. List

**Bold** and _Italic_ and `Code` text

[Link](url) and ![Image](src)
```

For more details see [GitHub Flavored Markdown](https://guides.github.com/features/mastering-markdown/).

### Jekyll Themes

Your Pages site will use the layout and styles from the Jekyll theme you have selected in your [repository settings](https://github.com/ajaiswal96/15418-Final-Project/settings). The name of this theme is saved in the Jekyll `_config.yml` configuration file.

### Support or Contact

Having trouble with Pages? Check out our [documentation](https://help.github.com/categories/github-pages-basics/) or [contact support](https://github.com/contact) and we’ll help you sort it out.
