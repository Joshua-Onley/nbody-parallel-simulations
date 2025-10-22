#include <iostream>
#include <cstdlib>
#include <random>
#include <chrono>
// #include <omp.h>
// #include <mpi.h>

// #define G 6.6742e-11 // Gravitational constant in SI units
#define G 1.0 // Gravitational constant

inline void calcForce(
        double pAx, double pAy, double pAz, // position of body A [in]
        double pBx, double pBy, double pBz, // position of body B [in]
        double mA, double mB, // masses of body A and body b [in]
        double* Fx, double* Fy, double* Fz // force [out]
        ) {
    double dx = pBx - pAx;
    double dy = pBy - pAy;
    double dz = pBz - pAz;
    double rSquared = dx*dx + dy*dy + dz*dz;
    double r = sqrt(rSquared);
    double F = (G * mA * mB) / (rSquared * r);
    *Fx = F * dx;
    *Fy = F * dy;
    *Fz = F * dz;
}

int main() {
    srand((unsigned) time(nullptr));
    double t = 0.0; // initial time
    double dt = 0.1; // time-step size
    double T = 2000.0; // final time
    int N = 840; // number of bodies
    auto force = new double[N][3]; // each of the forces

    // random number generator seeded with current clock time
    std::default_random_engine engine(
        std::chrono::system_clock::now().time_since_epoch().count()
    );

    auto mass = new double[N]; // each of the body masses
    double mMean = 1.0, mStdDev = 0.0; // mean and standard deviation of mass values
    std::normal_distribution<double> massDistribution(mMean,mStdDev);
    for(int n = 0; n < N; n++) {
        mass[n] = massDistribution(engine);
    }

    auto p = new double[N][3]; // each of the body positions
    double pMean = 0.0, pStdDev = 1.0; // mean and standard deviation of position values
    std::normal_distribution<double> positionDistribution(pMean,pStdDev);
    for(int n = 0; n < N; n++) {
        p[n][0] = positionDistribution(engine);
        p[n][1] = positionDistribution(engine);
        p[n][2] = positionDistribution(engine);
    }

    auto v = new double[N][3]; // each of the body velocities
    double vMean = 0.0, vStdDev = 1.0; // mean and standard deviation of velocity values
    std::normal_distribution<double> velocityDistribution(vMean,vStdDev);
    for (int n = 0; n < N; n++) {
        v[n][0] = velocityDistribution(engine);
        v[n][1] = velocityDistribution(engine);
        v[n][2] = velocityDistribution(engine);
    }

    // mark start time
    auto start = std::chrono::high_resolution_clock::now();

    while (t < T) {

        for (int i=0; i<N; i++) {
            force[i][0] = 0.0;
            force[i][1] = 0.0;
            force[i][2] = 0.0;
        }
        for (int i=0; i<N; i++) {
            for (int j=i+1; j<N; j++) {
                double Fx, Fy, Fz;
                calcForce(p[i][0], p[i][1], p[i][2], p[j][0], p[j][1], p[j][2],
                          mass[i], mass[j], &Fx, &Fy, &Fz);
                
                force[i][0] += Fx;
                force[i][1] += Fy;
                force[i][2] += Fz;
                force[j][0] += Fx;
                force[j][1] += Fy;
                force[j][2] += Fz;
                
            }
        }
        for (int i=0; i<N; i++) {
            v[i][0] += force[i][0]/mass[i] * dt;
            v[i][1] += force[i][1]/mass[i] * dt;
            v[i][2] += force[i][2]/mass[i] * dt;

            p[i][0] += v[i][0]*dt;
            p[i][1] += v[i][1]*dt;
            p[i][2] += v[i][2]*dt;
        }
        t += dt;
    }

    // measure runtime
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    printf("Time taken: %f seconds", duration.count()/1000000.0);


    delete[](mass);
    delete[](force);
    delete[](p);
    delete[](v);

}




