=== Overview ===

To run the code, execute run.m
Make sure to have casadi library available in the correct location 
(see bidirectional_bending_nonlinear.m Line 31)

The main file is bidirectional_bending_nonlinear.m. Here the NLP is implemented.
The following functions are used here:

- assemble_A calculates the matrix for the hinge conditions
- init_configuration returns the DOFs for the trivial configuration
- plot_skeleton_square plots a given configuration

Furthermore the follwing functions were implemented:
- calculate_neighbor explicitly calculates the rotation matrix of a 
neighboring cross such that the hinge condition is satisfied

- calculate_neighbor is used in exact_sol_neighbor where it is used
to plot a 2x2 configuration

- omega is the implementation of the omega function in the thesis

- omega_plot creates plots of the omega functions

- cross_like_structure illustrates one cross using a plot

- thee_by_three_non_existence provides a visual illustration of the
non existence proof in the thesis for the rigid case

- measure_angles numerically measures the hinge angles of a 2x2 configuration
