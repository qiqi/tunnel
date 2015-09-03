#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>

// ------------------------------------------------------------------------ //
//                                 ARRY2D                                   //
// ------------------------------------------------------------------------ //

struct Array2D
{
    double * ptr;
    int nVar, nx, ny, nSkip;

    double * operator ()(int ix, int iy) const {
        return ptr + nVar * (iy + ix * nSkip);
    }
};

Array2D subarray(const Array2D & base, int iVar,
        int ix0, int ix1, int iy0, int iy1)
{
    double * ptr = base(ix0, iy0) + iVar;
    return Array2D { ptr, base.nVar, ix1 - ix0, iy1 - iy0, base.nSkip };
}

void swap(Array2D & w0, Array2D & w1) 
{
    assert(w0.nVar == w1.nVar);
    assert(w0.nx == w1.nx);
    assert(w0.ny == w1.ny);
    assert(w0.nSkip == w1.nSkip);

    double * tmp = w0.ptr;
    w0.ptr = w1.ptr;
    w1.ptr = tmp;
}

#define FOR_K4 for(int k = 0; k < 4; ++k)

#define FOR_IXY(w) for(int ix=0; ix<w.nx; ++ix) for(int iy=0; iy<w.ny; ++iy)

// ------------------------------------------------------------------------ //
//                           EULER EQUATION                                 //
// ------------------------------------------------------------------------ //

void rhs(double * resid, const double * wPtr[5], double dx)
{
    const double gamma = 1.4;
    const int S = 1, N = 2, E = 3, W = 4;
    
    double rho[5], u[5], v[5], p[5];
    for (int i = 0; i < 5; ++ i) {
        rho[i] = wPtr[i][0] * wPtr[i][0];
        u[i] = wPtr[i][1] / wPtr[i][0];
        v[i] = wPtr[i][2] / wPtr[i][0];
        p[i] = wPtr[i][3];
    }

    double mass = (rho[E] * u[E] - rho[W] * u[W]) / (2*dx)
               + (rho[N] * v[N] - rho[S] * v[S]) / (2*dx);
    double momentum_x = (rho[E]*u[E]*u[E] - rho[W]*u[W]*u[W]
                      + rho[0]*u[0] * (u[E] - u[W])) / (4*dx)
                     + (rho[N]*u[N]*v[N] - rho[S]*u[S]*v[S]
                      + rho[0]*v[0] * (u[N] - u[S])) / (4*dx)
                     + (p[E] - p[W]) / (2*dx);
    double momentum_y = (rho[E]*u[E]*v[E] - rho[W]*u[W]*v[W]
                      + rho[0]*u[0] * (v[E] - v[W])) / (4*dx)
                     + (rho[N]*v[N]*v[N] - rho[S]*v[S]*v[S]
                      + rho[0]*v[0] * (v[N] - v[S])) / (4*dx)
                     + (p[N] - p[S]) / (2*dx);
    double energy = gamma * (p[E]*u[E] - p[W]*u[W]
                          + p[N]*v[N] - p[S]*v[S]) / (2*dx)
                 + (gamma - 1) * (u[0] * (p[E] - p[W])
                                + v[0] * (p[N] - p[S])) / (2*dx);

    double r = wPtr[0][0];
    resid[0] = 0.5 * mass / r;
    resid[1] = momentum_x / r;
    resid[2] = momentum_y / r;
    resid[3] = energy;
}

// ------------------------------------------------------------------------ //
//                               UNIT TEST                                  //
// ------------------------------------------------------------------------ //

void ddtConserved(Array2D w, double dx) 
{
    const double gamma = 1.4;

    double dmass = 0, dmomentum_x = 0, dmomentum_y = 0, denergy = 0;
    FOR_IXY(w) {
        if (iy != 1) continue;
        // const int S = 1, N = 2, E = 3, W = 4;
        double * w0 = w(ix, iy);
        const double * wPtr[5] = {w0, w(ix, iy+1), w(ix, iy-1),
                                     w(ix+1, iy), w(ix-1, iy)};
        double res[4];
        rhs(res, wPtr, dx);

        double r = w0[0], u = w0[1] / w0[0], v = w0[2] / w0[0];

        double ddt_rho = -2 * res[0] * r;
        double ddt_rhou = -res[1] * r + 0.5 * u * ddt_rho;
        double ddt_rhov = -res[2] * r + 0.5 * v * ddt_rho;
        double ddt_p = -res[3];
        double ddt_rhou2 = 2 * u * ddt_rhou - u * u * ddt_rho;
        double ddt_rhov2 = 2 * v * ddt_rhov - v * v * ddt_rho;

        dmass += ddt_rho;
        dmomentum_x += ddt_rhou;
        dmomentum_y += ddt_rhov;
        denergy += ddt_p / (gamma - 1) + 0.5 * (ddt_rhou2 + ddt_rhov2);
    }
    printf("%e %e %e %e\n", dmass, dmomentum_x, dmomentum_y, denergy);
}

// ------------------------------------------------------------------------ //
//                              TIME STEPPING                               //
// ------------------------------------------------------------------------ //

void initialize(Array2D w) {
    const double gamma = 1.4, R = 287, T0 = 300, p0 = 101325, M0 = .2;

    double rho0 = p0 / (R * T0);
    double r0 = sqrt(rho0);
    double c0 = sqrt(gamma * R * T0);
    double u0 = c0 * M0;

    FOR_IXY(w) {
        w(ix, iy)[0] = r0 * (1 + 0.1 * rand() / RAND_MAX);
        w(ix, iy)[1] = r0 * u0 * (1 + 0.1 * rand() / RAND_MAX);
        w(ix, iy)[2] = r0 * u0 * 0.1 * rand() / RAND_MAX;
        w(ix, iy)[3] = p0 * (1 + 0.1 * rand() / RAND_MAX);
    }
}

void fillGhost(Array2D w)
{
    for (int ix = 0; ix < w.nx; ++ ix) {
        memcpy(w(ix, w.ny), w(ix, 0), sizeof(double) * 4);
        memcpy(w(ix, -1), w(ix, w.ny-1), sizeof(double) * 4);
    }
    for (int iy = 0; iy < w.ny; ++ iy) {
        memcpy(w(w.nx, iy), w(0, iy), sizeof(double) * 4);
        memcpy(w(-1, iy), w(w.nx-1, iy), sizeof(double) * 4);
    }
}

void stepFE(Array2D w1, Array2D w0, double dx, double dt) 
{
    assert(w1.nx == w0.nx);
    assert(w1.ny == w0.ny);

    FOR_IXY(w0) {
        double * w0i = w0(ix, iy), * w1i = w1(ix, iy);
        // const int S = 1, N = 2, E = 3, W = 4;
        const double * wPtr[5] = {w0i, w0(ix, iy+1), w0(ix, iy-1),
                                      w0(ix+1, iy), w0(ix-1, iy)};
        rhs(w1i, wPtr, dx);

        FOR_K4 w1i[k] = w0i[k] + dt * w1i[k];
    }
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    const int nx = 250, ny = 50;
    Array2D base { new double[(nx + 2) * (ny + 2) * 8], 8, nx+2, ny+2, ny+2 };
    Array2D w0 = subarray(base, 0, 1, nx+1, 1, ny+1);
    Array2D w1 = subarray(base, 4, 1, nx+1, 1, ny+1);

    memset(base.ptr, 0, sizeof(double) * (nx + 2) * (ny + 2) * 8);

    initialize(w0);
    fillGhost(w0);
    ddtConserved(w0, 0.01);
    
    for (int iStep = 0; iStep < -10; ++ iStep) {
        fillGhost(w0);
        stepFE(w1, w0, 0.01, 0.001);
        ddtConserved(w1, 0.01);
        swap(w1, w0);
    }

    MPI_Finalize();
    return 0;
}

