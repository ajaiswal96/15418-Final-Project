## Proposal

### Team

Anubhav Jaiswal  (ajaiswal)


Rohit Pillai  (rrpillai)

### Summary

We are going to implement a parallel version of a fluid simulator to run on Nvidia GTX 1080 GPU's on the Gates machines and compare that to a sequential implementation on a CPU.

### Background

We are implementing a parallel computational fluid dynamics simulator to run on GPU's. For the sake of simplicity, we are assuming the fluids to be homogeneous and incompressible. We plan to take advantage of the fact that GPU's are optimized to render textures, and a texture is just a 2D vector storing multiple values. We can thus store our simulation grids as these textures. Instead of having each pixel in the texture map correspond to a particle, we are going to have a grid of particles and store each grid componenet as a pixel texture. The inherent structure of CFD simulations thus makes it best suited to be run through CUDA kernels on the GPU. 

We will be representing the fluid as a simulation grid, which is indexed in by their position, and stores a time and velocity (since velocity and location change based on time). 

![Navier-Stokes](https://github.com/ajaiswal96/15418-Final-Project/blob/master/assets/Screen%20Shot%202017-04-10%20at%202.56.03%20PM.png | width=200)


We plan to use the GPU by representing the 
We will be using the Navier-Stokes equation to 
 

Homogeneous Incompressible

Navier-Stokes



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
