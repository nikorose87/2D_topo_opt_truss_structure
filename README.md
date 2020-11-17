# 2D Topology Optimization of truss structure

<div class=caption><b> Caption: </b> Truss structure problem, proposed by Christensen in his book. See the reference for more information. </div>

![](https://github.com/nikorose87/2D_topo_opt_truss_structure/blob/master/Structure.png)


As an academic project, I developed a 2D topology optimization of a truss structure (see above) in python, the algorithm took and performed many examples available in the web. The topology optimization was done under the nested formulation, It means that the process is static. Many examples of this problem have been done, however, this is able to analyze structures of any size, and shows a nice visualization with VTK library and Matplotlib. I hope this work could be useful for you.

Five functions were implemented and/or modified in order to accomplish the process, being inspired mostly by the codes proposed in the book of (Ferreira, 2009):

- _structure2_ Build the element matrix in an array, containing the coordinates in each element.
- _indexation_ Returns the element matrix but with the node indexation.
- _solution_ Function to find solution in terms of global displacements. The code was based on the book of (Ferreira, 2009)
- _formStiffness2Dtruss_ Build the Stiffness matrix of the entire problem.
- _outputDisplacementsReactions_ Calculates the force reactions once the displacements are obtained.
- _stresses2Dtruss_ Calculates the stresses in each element.

Another class was proposed for visualizing the answer, which are:

- _displayTruss:_     Function for calculating stiffness matrix, displacements, stresses of 2D trusses
    with the VTK libraries. Original code proposed by (Sukhbinder, 2015)
- _visualizeAreas:_ A function to plot the results in terms of areas (width of the line) and in terms of Stresses (Blue is compression, red is tension).

Those classes are in the file [FEAtruss.py](https://github.com/nikorose87/2D_topo_opt_truss_structure/blob/master/FEAtrussDP3Opt.py), and the main file is called [main_FEA_viz_opt1.py](https://github.com/nikorose87/2D_topo_opt_truss_structure/blob/master/main_FEA_viz_opt1.py). If you want to change the number of rows and/or columns, change the following in the code:

```
rows = 8
cols = 4
```

You can change the properties of the material in the code as well.

Now, we are ready to execute the code as follows:

<div class=caption><b> Caption: </b> Topology Opt. of a 2D truss Structure at 51 iterations.</div>

![](https://github.com/nikorose87/2D_topo_opt_truss_structure/blob/master/output_20_2.png)

<div class=caption><b> Caption: </b> Topology Opt. of a 2D truss Structure at 206 iterations.</div>

![](https://github.com/nikorose87/2D_topo_opt_truss_structure/blob/master/output_24_3.png)

The results of stresses and deformations should be seen as this:

![](https://github.com/nikorose87/2D_topo_opt_truss_structure/blob/master/deformations.png)

## References:

- Sukhbinder, S. (2015). SukhbinderSingh.com. Retrieved July 15, 2018, from https://sukhbinder.wordpress.com/2015/03/10/displaying-truss/#more-2877.

- Ferreira, A. J. M. (2009). MATLAB Codes for Finite Element Analysis, 157. https://doi.org/10.1007/978-1-4020-9200-8

- Christensen, P. W., & Klabring, A. (2009). An introduction to structural optimization. https://doi.org/10.1007/978-1-4020-8666-3
