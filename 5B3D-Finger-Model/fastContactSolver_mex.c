#include "mex.h"
#include "matrix.h"
#include <math.h>
#include <string.h>

/* --- Helpers ------------------------------------------------------------ */

static double norm2_4(const double v[4]) {
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2] + v[3]*v[3]);
}

/* Solve SPD 4x4: (A) x = b via Cholesky. A is row-major 4x4, overwritten. Returns 1 on success */
static int solve_spd4(double A[16], const double b[4], double x[4]) {
    double L_mat[16];
    for (int i=0; i<16; ++i) L_mat[i] = A[i];
    
    // Factorization
    for (int i=0; i<4; ++i) {
        for (int j=0; j<=i; ++j) {
            double sum = L_mat[i*4 + j];
            for (int k=0; k<j; ++k) sum -= L_mat[i*4 + k]*L_mat[j*4 + k];
            if (i==j) {
                if (sum <= 0.0) return 0;
                L_mat[i*4 + j] = sqrt(sum);
            } else {
                L_mat[i*4 + j] = sum / L_mat[j*4 + j];
            }
        }
        for (int j=i+1; j<4; ++j) L_mat[i*4 + j] = 0.0;
    }
    
    // Forward substitution
    double y[4];
    for (int i=0; i<4; ++i) {
        double sum = b[i];
        for (int k=0; k<i; ++k) sum -= L_mat[i*4 + k]*y[k];
        y[i] = sum / L_mat[i*4 + i];
    }
    
    // Backward substitution
    for (int i=3; i>=0; --i) {
        double sum = y[i];
        for (int k=i+1; k<4; ++k) sum -= L_mat[k*4 + i]*x[k];
        x[i] = sum / L_mat[i*4 + i];
    }
    return 1;
}

/* --- Residual ----------------------------------------------------------- */

static void residual_only(const double z[4], double t1, double t2,
                          double L, double t, double xc, double yc0,
                          double r[4])
{
    const double p1=z[0], p2=z[1], a=z[2], b=z[3];
    const double c1 = cos(t1),  s1 = sin(t1);
    const double c12= cos(t1+t2), s12 = sin(t1+t2);
    
    const double x1 = p1*c1 - t*sin(t1);
    const double y1 = p1*s1 + t*c1;
    const double x2 = L*c1 + p2*c12 - t*sin(t1+t2);
    const double y2 = L*s1 + p2*s12 + t*c12;
    
    const double yc = yc0 + b;
    const double Y1 = y1 - yc;
    const double Y2 = y2 - yc;
    const double a2 = a*a, b2 = b*b;
    
    r[0] = ((x1 - xc)*(x1 - xc))/a2 + (Y1*Y1)/b2 - 1.0;
    r[1] = ((x2 - xc)*(x2 - xc))/a2 + (Y2*Y2)/b2 - 1.0;
    r[2] = ((x1 - xc)/a2)*c1  + (Y1/b2)*s1;
    r[3] = ((x2 - xc)/a2)*c12 + (Y2/b2)*s12;
}

/* --- Newton Solver Core ------------------------------------------------- */

static void newton_contact_one(double L, double t, double xc, double yc0,
                               double t1, double t2, const double z_start[4],
                               double z_out[4], int *converged_out)
{
    const int maxIter = 250;
    const double tolR = 1e-6;
    const double tolStep = 1e-8;
    const double EPS = 1e-6;

    double y[4];
    for (int i=0; i<4; ++i) {
        double v = z_start[i] > EPS ? z_start[i] : EPS;
        y[i] = log(v);
    }

    int converged = 0;
    for (int k=0; k<maxIter; ++k) {
        double z[4] = { exp(y[0]), exp(y[1]), exp(y[2]), exp(y[3]) };
        const double p1=z[0], p2=z[1], a=z[2], b=z[3];
        const double c1 = cos(t1),  s1 = sin(t1);
        const double c12= cos(t1+t2), s12 = sin(t1+t2);

        const double x1 = p1*c1 - t*sin(t1);
        const double y1p = p1*s1 + t*c1;
        const double x2 = L*c1 + p2*c12 - t*sin(t1+t2);
        const double y2p = L*s1 + p2*s12 + t*c12;

        const double yc = yc0 + b;
        const double Y1 = y1p - yc;
        const double Y2 = y2p - yc;
        const double a2=a*a, b2=b*b, a3=a2*a, b3=b2*b;

        double r[4];
        r[0] = ((x1 - xc)*(x1 - xc))/a2 + (Y1*Y1)/b2 - 1.0;
        r[1] = ((x2 - xc)*(x2 - xc))/a2 + (Y2*Y2)/b2 - 1.0;
        r[2] = ((x1 - xc)/a2)*c1  + (Y1/b2)*s1;
        r[3] = ((x2 - xc)/a2)*c12 + (Y2/b2)*s12;

        if (norm2_4(r) < tolR) {
            for (int i=0; i<4; ++i) z_out[i] = z[i];
            converged = 1;
            break;
        }

        double Jz[16] = {0};
        const double dx1_dp1=c1, dy1_dp1=s1;
        const double dx2_dp2=c12, dy2_dp2=s12;

        // Jacobian Jz wrt z=[p1 p2 a b]
        Jz[0*4+0] = 2.0*(x1-xc)/a2 * dx1_dp1 + 2.0*Y1/b2 * dy1_dp1;
        Jz[0*4+2] = -2.0*(x1-xc)*(x1-xc)/a3;
        Jz[0*4+3] = (-2.0*Y1)/b2 - 2.0*(Y1*Y1)/b3;

        Jz[1*4+1] = 2.0*(x2-xc)/a2 * dx2_dp2 + 2.0*Y2/b2 * dy2_dp2;
        Jz[1*4+2] = -2.0*(x2-xc)*(x2-xc)/a3;
        Jz[1*4+3] = (-2.0*Y2)/b2 - 2.0*(Y2*Y2)/b3;

        Jz[2*4+0] = (dx1_dp1/a2)*c1 + (dy1_dp1/b2)*s1;
        Jz[2*4+2] = -2.0*(x1-xc)/a3 * c1;
        Jz[2*4+3] = (-1.0/b2)*s1 + (-2.0*Y1/b3)*s1;

        Jz[3*4+1] = (dx2_dp2/a2)*c12 + (dy2_dp2/b2)*s12;
        Jz[3*4+2] = -2.0*(x2-xc)/a3 * c12;
        Jz[3*4+3] = (-1.0/b2)*s12 + (-2.0*Y2/b3)*s12;

        double Jy[16];
        for (int c=0; c<4; ++c) {
            const double s = z[c];
            for (int r_=0; r_<4; ++r_) Jy[r_*4 + c] = Jz[r_*4 + c]*s;
        }

        double g[4] = {0,0,0,0};
        double H[16] = {0};
        for (int i=0; i<4; ++i) {
            double gi = 0.0;
            for (int k=0; k<4; ++k) gi += Jy[k*4 + i] * r[k];
            g[i] = gi;
            for (int j=i; j<4; ++j) {
                double hij = 0.0;
                for (int k=0; k<4; ++k) hij += Jy[k*4 + i] * Jy[k*4 + j];
                H[i*4 + j] = hij;
                H[j*4 + i] = hij;
            }
        }

        double tr = H[0] + H[5] + H[10] + H[15];
        double mu = 1e-9 * (tr/4.0 + 1.0);
        for (int d=0; d<4; ++d) H[d*4 + d] += mu;

        double rhs[4] = {-g[0], -g[1], -g[2], -g[3]};
        double dy[4];
        if (!solve_spd4(H, rhs, dy)) {
            for (int i=0; i<4; ++i) {
                double diag = (i==0?H[0]:(i==1?H[5]:(i==2?H[10]:H[15])));
                dy[i] = -g[i] / (diag > 0 ? diag : 1.0);
            }
        }

        double step = 1.0;
        const double c_armijo = 1e-4;
        double rnorm = norm2_4(r);
        for (int bt=0; bt<8; ++bt) {
            double y_try[4];
            for (int i=0; i<4; ++i) y_try[i] = y[i] + step*dy[i];
            double zt[4];
            for (int i=0; i<4; ++i) {
                double v = exp(y_try[i]);
                zt[i] = v > EPS ? v : EPS;
            }
            double rt[4];
            residual_only(zt, t1, t2, L, t, xc, yc0, rt);
            if (norm2_4(rt) < (1.0 - c_armijo*step)*rnorm) {
                for (int i=0; i<4; ++i) y[i] = log(zt[i]);
                break;
            } else {
                step *= 0.5;
            }
        }

        double steplen2 = step*step*(dy[0]*dy[0] + dy[1]*dy[1] + dy[2]*dy[2] + dy[3]*dy[3]);
        if (sqrt(steplen2) < tolStep) {
            for (int i=0; i<4; ++i) z_out[i] = exp(y[i]);
            converged = 1;
            break;
        }
    }
    if (!converged) {
        for (int i=0; i<4; ++i) z_out[i] = exp(y[i]);
    }
    *converged_out = converged;
}

/* --- MEX Entry Point ---------------------------------------------------- */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs < 3) {
        mexErrMsgIdAndTxt("fastContactsolver:nrhs", "Need param, t1, t2 (and optional z0).");
    }

    /* Extract Params */
    const mxArray *a_f = mxGetField(prhs[0],0,"a");
    const mxArray *t_f = mxGetField(prhs[0],0,"t");
    const mxArray *O_f = mxGetField(prhs[0],0,"O");
    if (!a_f || !t_f || !O_f) mexErrMsgIdAndTxt("fastContactsolver:param", "Missing struct fields.");

    const double L   = mxGetScalar(a_f);
    const double t_p = mxGetScalar(t_f);
    const double *O  = mxGetPr(O_f);
    const double xc  = O[0];
    const double yc0 = O[1];

    /* Extract Vectors */
    mwSize N = mxGetNumberOfElements(prhs[1]);
    const double *t1p = mxGetPr(prhs[1]);
    const double *t2p = mxGetPr(prhs[2]);

    /* Handle Initial Guess and Warm Start Mode */
    double z_initial_static[4];
    int useWarmStart = 0;

    if (nrhs < 4 || mxIsEmpty(prhs[3])) {
        // DEFAULT MODE: Warm Start Enabled
        z_initial_static[0]=0.02; z_initial_static[1]=0.02; 
        z_initial_static[2]=0.03; z_initial_static[3]=0.03;
        useWarmStart = 1;
    } else {
        // USER MODE: Warm Start Disabled (use static z0 for all array elements)
        const double *z0p = mxGetPr(prhs[3]);
        for (int i=0; i<4; ++i) z_initial_static[i] = z0p[i];
        useWarmStart = 0;
    }

    /* Prepare Outputs */
    plhs[0] = mxCreateDoubleMatrix(N, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(N, 1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(N, 1, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(N, 1, mxREAL);
    double *p1_out = mxGetPr(plhs[0]);
    double *p2_out = mxGetPr(plhs[1]);
    double *a_out  = mxGetPr(plhs[2]);
    double *b_out  = mxGetPr(plhs[3]);

    /* Active Guess Variable */
    double z_active[4];
    for (int k=0; k<4; ++k) z_active[k] = z_initial_static[k];

    /* Main Loop */
    for (mwSize i=0; i<N; ++i) {
        double zi[4];
        int ok1 = 0;

        newton_contact_one(L, t_p, xc, yc0, t1p[i], t2p[i], z_active, zi, &ok1);

        if (!ok1) {
            /* Fallback Logic */
            const double base[4] = {0.01,0.01,0.02,0.02};
            const double defg[4] = {0.02,0.02,0.03,0.03};
            double z_reset[4];
            for (int k=0; k<4; ++k) {
                double cand = 0.5*z_active[k] + 0.5*defg[k];
                z_reset[k] = (cand > base[k]) ? cand : base[k];
            }
            int ok2 = 0;
            newton_contact_one(L, t_p, xc, yc0, t1p[i], t2p[i], z_reset, zi, &ok2);
        }

        p1_out[i] = zi[0];
        p2_out[i] = zi[1];
        a_out[i]  = zi[2];
        b_out[i]  = zi[3];

        /* Apply logic: Only update z_active if warm starting is permitted */
        if (useWarmStart) {
            for (int k=0; k<4; ++k) z_active[k] = zi[k];
        } else {
            /* Keep using the same user-provided z0 */
            for (int k=0; k<4; ++k) z_active[k] = z_initial_static[k];
        }
    }
}