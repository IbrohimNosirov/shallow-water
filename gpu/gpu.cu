#include <cuda.h>
#include <cuda_runtime.h>

#include <math.h>

#include "../common/common.hpp"
#include "../common/solver.hpp"

/**
 * This is your initialization function! We pass in h0, u0, and v0, which are
 * your initial height, u velocity, and v velocity fields. You should send these
 * grids to the GPU so you can do work on them there, and also these other fields.
 * Here, length and width are the length and width of the domain, and nx and ny are
 * the number of grid points in the x and y directions. H is the height of the water
 * column, g is the acceleration due to gravity, and dt is the time step size.
 * The rank and num_procs variables are unused here, but you will need them
 * when doing the MPI version.
 */
void init(double *h0, double *u0, double *v0, double length_, double width_, int nx_, int ny_, double H_, double g_, double dt_, int rank_, int num_procs_)
{
    // @TODO: your code here
}

/**
 * This is your step function! Here, you will actually numerically solve the shallow
 * water equations. You should update the h, u, and v fields to be the solution after
 * one time step has passed.
 */
void step()
{
    // @TODO: Your code here
}

/**
 * This is your transfer function! You should copy the h field back to the host
 * so that the CPU can check the results of your computation.
 */
void transfer(double *h_host)
{
    // @TODO: Your code here
}

/**
 * This is your finalization function! You should free all of the memory that you
 * allocated on the GPU here.
 */
void free_memory()
{
    // @TODO: Your code here
}