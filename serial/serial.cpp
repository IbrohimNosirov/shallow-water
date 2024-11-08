#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>

#include "../common/common.hpp"
#include "../common/solver.hpp"

int nx, ny;
int nh = 2;
int nx_halo, ny_halo;

#define h(i, j) h[(i) * (ny_halo) + (j)]
#define u(i, j) u[(i) * (ny) + (j)]
#define v(i, j) v[(i) * (ny_halo) + (j)]

double *h, *u, *v, *dh, *du, *dv, *dh1, *du1, *dv1, *dh2, *du2, *dv2;
double H, g, dx, dy, dt;

void init(double *h0, double *u0, double *v0,
          double length_, double width_,
          int nx_, int ny_,
          double H_, double g_, double dt_, int rank_, int num_procs_)
{
    nx_halo = nx_ + nh;
    ny_halo = ny_ + nh;
    nx = nx_;
    ny = ny_;

    // We allocate memory
    h = (double *)calloc((nx_halo)*(ny_halo), sizeof(double));
    u = (double *)calloc((nx_halo + 1)*ny, sizeof(double));
    v = (double *)calloc(nx*(ny_halo+1), sizeof(double));

    dh = (double *)calloc(nx * ny, sizeof(double));
    du = (double *)calloc(nx * ny, sizeof(double));
    dv = (double *)calloc(nx * ny, sizeof(double));

    dh1 = (double *)calloc(nx * ny, sizeof(double));
    du1 = (double *)calloc(nx * ny, sizeof(double));
    dv1 = (double *)calloc(nx * ny, sizeof(double));

    dh2 = (double *)calloc(nx * ny, sizeof(double));
    du2 = (double *)calloc(nx * ny, sizeof(double));
    dv2 = (double *)calloc(nx * ny, sizeof(double));

    // nx + 1 because scenarios.cpp is like this. 
    for (int i = 0; i < nx+1; i++) {
        for (int j = 0; j < ny+1; j++) {
            h(i,j) = h0[i*(ny+1) + j];
            u(i,j) = u0[i*(ny+1) + j];
            v(i,j) = v0[i*(ny + 2) + j];
        }
    }

    //Debug
    //for (int i = 0; i < nx; i++) {
    //    for (int j = 0; j < ny; j++) {
    //        h(i,j) = h0[i*ny + j] + i*ny + j + 1;
    //        u(i,j) = u0[i*ny + j] + i*ny + j + 1;
    //        v(i,j) = v0[i*(ny + 2) + j] + i*(ny+2) + j + 1;
    //    }
    //}

    // Debug.
    //std::cout << "Matrix h after allocation:" << std::endl;
    //for (int i = 0; i < nx_halo; i++) {
    //    for (int j = 0; j < ny_halo; j++) {
    //        std::cout << h(i,j) << " ";
    //    }
    //    std::cout << std::endl;
    //}

    //std::cout << "Matrix u after allocation:" << std::endl;
    //for (int i = 0; i < nx_halo+1; i++) {
    //    for (int j = 0; j < ny; j++) {
    //        std::cout << u(i,j) << " ";
    //    }
    //    std::cout << std::endl;
    //}

    //std::cout << "Matrix v after allocation:" << std::endl;
    //for (int i = 0; i < nx; i++) {
    //    for (int j = 0; j < ny+2; j++) {
    //        std::cout << v(i,j) << " ";
    //    }
    //    std::cout << std::endl;
    //}

    H = H_;
    g = g_;

    dx = length_ / nx_;
    dy = width_ / ny_;

    dt = dt_;
}

void compute_derivatives()
{
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            dh(i, j) = -H * (du_dx(i, j) + dv_dy(i, j));
            du(i, j) = -g * dh_dx(i, j);
            dv(i, j) = -g * dh_dy(i, j);
        }
    }
}

void multistep(double a1, double a2, double a3)
{
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            h(i, j) += (a1*dh(i,j) + a2*dh1(i,j) + a3*dh2(i,j))*dt;
            u(i + 1, j) += (a1*du(i,j) + a2*du1(i,j) + a3*du2(i,j))*dt;
            v(i, j + 1) += (a1*dv(i,j) + a2*dv1(i,j) + a3*dv2(i,j))*dt;
        }
    }
}

// first goes to last
void compute_ghost_horizontal()
{

    for (int j = 0; j < ny; j++)
    {
        h(nx, j) = h(0, j);
    }

    //for (int i = nx; i < nx_halo; i++) {
    //    for (int j = 0; j < ny; j++) {
    //        h(i, j) = h(nx_halo-(i+1), j);
    //        std::cout << "h(nx_halo), j:" << h(nx_halo, j) << std::endl;
    //    }
    //}
    
}

// first goes to last
void compute_ghost_vertical()
{
    for (int i = 0; i < nx; i++)
    {
        h(i, ny) = h(i, 0);
    }
}

// last goes to first
void compute_boundaries_horizontal()
{
    for (int j = 0; j < ny; j++)
    {
        u(0, j) = u(nx, j);
    }
}

// last goes to first
void compute_boundaries_vertical()
{
    for (int i = 0; i < nx; i++)
    {
        v(i, 0) = v(i, ny);
    }
}

void swap_buffers()
{
    double *tmp;

    tmp = dh2;
    dh2 = dh1;
    dh1 = dh;
    dh = tmp;

    tmp = du2;
    du2 = du1;
    du1 = du;
    du = tmp;

    tmp = dv2;
    dv2 = dv1;
    dv1 = dv;
    dv = tmp;
}

int t = 0;

void step()
{
    std::cout << "Matrix h before ghost:" << std::endl;
    for (int i = 0; i < nx_halo; i++) {
        for (int j = 0; j < ny_halo; j++) {
            std::cout << h(i,j) << " ";
        }
        std::cout << std::endl;
    }
    compute_ghost_horizontal();
    compute_ghost_vertical();

    std::cout << "Matrix h after ghost:" << std::endl;
    for (int i = 0; i < nx_halo; i++) {
        for (int j = 0; j < ny_halo; j++) {
            std::cout << h(i,j) << " ";
        }
        std::cout << std::endl;
    }

    compute_derivatives();

    double a1, a2, a3;

    if (t == 0)
    {
        a1 = 1.0;
    }
    else if (t == 1)
    {
        a1 = 3.0 / 2.0;
        a2 = -1.0 / 2.0;
    }
    else
    {
        a1 = 23.0 / 12.0;
        a2 = -16.0 / 12.0;
        a3 = 5.0 / 12.0;
    }

    multistep(a1, a2, a3);

    //std::cout << "Matrix u before boundary:" << std::endl;
    //for (int i = 0; i < nx_halo+1; i++) {
    //    for (int j = 0; j < ny; j++) {
    //        std::cout << u(i,j) << " ";
    //    }
    //    std::cout << std::endl;
    //}
    compute_boundaries_horizontal();
    //std::cout << "Matrix u after boundary:" << std::endl;
    //for (int i = 0; i < nx_halo+1; i++) {
    //    for (int j = 0; j < ny; j++) {
    //        std::cout << u(i,j) << " ";
    //    }
    //    std::cout << std::endl;
    //}

    //std::cout << "Matrix v before boundary:" << std::endl;
    //for (int i = 0; i < nx; i++) {
    //    for (int j = 0; j < ny_halo+1; j++) {
    //        std::cout << v(i,j) << " ";
    //    }
    //    std::cout << std::endl;
    //}
    compute_boundaries_vertical();
    //std::cout << "Matrix v after boundary:" << std::endl;
    //for (int i = 0; i < nx; i++) {
    //    for (int j = 0; j < ny_halo+1; j++) {
    //        std::cout << v(i,j) << " ";
    //    }
    //    std::cout << std::endl;
    //}

    swap_buffers();

    t++;
}

void transfer(double *h)
{
    return;
}

void free_memory()
{
    free(dh);
    free(du);
    free(dv);

    free(dh1);
    free(du1);
    free(dv1);

    free(dh2);
    free(du2);
    free(dv2);
}
