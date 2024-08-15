=== Overview over the Code ===

To run the code: Execute run.m
Make sure to have casadi libraries available in the correct location

The main file is "bending_model". Here the NLP is implemented. The
following functions are used here:

assemble_A3.m -> Assembles the sparse matrix for the hinge conditions

discrete_nabla.m -> implementation of the discrete nabla operator

init_configuration.m -> generates the trivial configuration

plot_exo_skeleton.m -> plots a given configuration

stiffness_and_load -> calculates the stiffness matrix for the energy discretization

transform3 -> calculates the necessary affine transformations for the reference plates



To run the test files, they need to be moved to the main directory.