#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "../common/common.hpp"
#include "../common/solver.hpp"

// Here we hold the number of cells we have in the x and y directions
int nx, ny;
int n_halo = 2;

double *h, *u, *v, *dh, *du, *dv, *dh1, *du1, *dv1, *dh2, *du2, *dv2;
double H, g, dx, dy, dt;

void init(double *h0, double *u0, double *v0, double length_, double width_,
int nx_, int ny_, double H_, double g_, double dt_, int rank_, int num_procs_)
{
    nx = nx_ + 2*n_halo;
    ny = ny_ + 2*n_halo;

    // we allocate space for h, u, and v
    h = (double *)calloc(nx * ny, sizeof(double));
    u = (double *)calloc(nx * ny, sizeof(double));
    v = (double *)calloc(nx * ny, sizeof(double));

    // We allocate memory for the derivatives
    dh = (double *)calloc(nx * ny, sizeof(double));
    du = (double *)calloc(nx * ny, sizeof(double));
    dv = (double *)calloc(nx * ny, sizeof(double));

    dh1 = (double *)calloc(nx * ny, sizeof(double));
    du1 = (double *)calloc(nx * ny, sizeof(double));
    dv1 = (double *)calloc(nx * ny, sizeof(double));

    dh2 = (double *)calloc(nx * ny, sizeof(double));
    du2 = (double *)calloc(nx * ny, sizeof(double));
    dv2 = (double *)calloc(nx * ny, sizeof(double));

    // build out halo
    for (int i = 0; i < nx_ ; i++) {
        for (int j = 0; j < ny_; j++) {
            int index = (i + n_halo) * ny + (j + n_halo);
            h[index] = h0[i * ny_ + j];
            u[index] = u0[i * ny_ + j];
            v[index] = v0[i * ny_ + j];
        }
    }

    H = H_;
    g = g_;

    dx = length_ / nx_;
    dy = width_ / ny_;

    dt = dt_;
}

inline int index(int i, int j) {
    return (i + n_halo) * ny + (j + n_halo);
}

void compute_dh()
{
    for (int i = n_halo; i < nx - n_halo; i++)
    {
        for (int j = n_halo; j < ny - n_halo; j++)
        {
            dh[index(i,j)] = -H * (du_dx(i, j) + dv_dy(i, j));
        }
    }
}

void compute_du_dv()
{
    for (int i = n_halo; i < nx - n_halo; i++)
    {
        for (int j = n_halo; j < ny - n_halo; j++)
        {
            du[index(i, j)] = -g * dh_dx(i, j);
            dv[index(i, j)] = -g * dh_dy(i, j);
        }
    }
}

void multistep(double a1, double a2, double a3)
{
    for (int i = n_halo; i < nx - n_halo; i++)
    {
        for (int j = n_halo; j < ny - n_halo; j++)
        {
            h[index(i, j)] += (a1 * dh[index(i, j)] + a2 * dh1[index(i, j)] +
                            a3 * dh2[index(i, j)]) * dt;
            // Remove the +1 from i+1 since index() already handles the offset
            u[index(i+1, j)] += (a1 * du[index(i, j)] + a2 * du1[index(i, j)]
                            + a3 * du2[index(i, j)]) * dt;
            v[index(i, j+1)] += (a1 * dv[index(i, j)] + a2 * dv1[index(i, j)]
                            + a3 * dv2[index(i, j)]) * dt;
        }
    }
}

// take a look at this function more closely
void compute_ghost_horizontal()
{
    for (int j = n_halo; j < ny - n_halo; j++)
    {
        h[index(nx - n_halo, j)] = h[index(n_halo, j)];
        h[index(n_halo - 1, j)] = h[index(nx - n_halo - 1, j)];
    }
}

void compute_ghost_vertical()
{
    for (int i = n_halo; i < nx - n_halo; i++)
    {
        h[index(i, ny - n_halo)] = h[index(i, n_halo)];
        h[index(i, n_halo - 1)] = h[index(i, ny - n_halo - 1)];
    }
}

void compute_boundaries_horizontal()
{
    for (int j = n_halo; j < ny - n_halo; j++)
    {
        u[index(n_halo - 1, j)] = u[index(nx - n_halo - 1, j)];
        u[index(nx - n_halo, j)] = u[index(n_halo, j)];
    }
}

void compute_boundaries_vertical()
{
    for (int i = n_halo; i < nx - n_halo; i++)
    {
        v[index(i, n_halo - 1)] = v[index(i, ny - n_halo - 1)];
        v[index(i, ny - n_halo)] = v[index(i, n_halo)];
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
int N = 3;

void step()
{
    if (t % N == 0) {
    // First, we compute our ghost cells as we need them for our derivatives
        compute_ghost_horizontal();
        compute_ghost_vertical();
    }

    // Next, we compute the derivatives of our fields
    compute_dh();
    compute_du_dv();

    // We set the coefficients for our multistep method
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

    // Finally, we compute the next time step using our multistep method
    multistep(a1, a2, a3);

    // We compute the boundaries for our fields, as they are (1) needed for the
    // next time step, and (2) aren't explicitly set in our multistep method
    if (t % N == 0) {
        compute_boundaries_horizontal();
        compute_boundaries_vertical();
    }

    // We swap the buffers for our derivatives so that we can use the
    // derivatives from the previous time steps in our multistep method, then
    // increment the time step counter
    swap_buffers();

    t++;
}

// Since all of our memory is already on the CPU, and specifically in the
// height field, we don't need to transfer anything
void transfer(double *h)
{
    return;
}

// We free all of the memory that we allocated. We didn't create the initial
// height or velocity fields, so we don't need to free them. They are the
// responsibility of the calling code.
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
