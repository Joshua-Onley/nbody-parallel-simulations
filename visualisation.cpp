#include <iostream>
#include <cstdlib>
#include <random>
#include <chrono>
#include <thread>
#include <SFML/Graphics.hpp>

//#define G 6.6742e-11 // Gravitational constant in SI units
#define G 1.0 // Gravitational constant

inline void calcForce(
    double pAx, double pAy, double pAz, // position of body A
    double pBx, double pBy, double pBz, // position of body B
    double mA, double mB, // masses of body A and B
    double* Fx, double* Fy, double* Fz // force [out]
) {
    double dx = pBx - pAx;
    double dy = pBy - pAy;
    double dz = pBz - pAz;

    // distance between the two bodies
    double rSquared = dx*dx + dy*dy + dz*dz;
    double r = sqrt(rSquared);

    // Gravitational force formula
    double F = (G * mA * mB) / (rSquared * r);

    *Fx = F * dx;
    *Fy = F * dy;
    *Fz = F * dz;
}

sf::Color randomColour() {
    // random RGB values for the colour
    int r = rand() % 250;
    int g = rand() % 5;
    int b = rand() % 5;
    return sf::Color(r, g, b);
}

int main() {
    srand((unsigned) time(nullptr));
    double t = 0.0; // initial time
    double dt = 0.001; // time-step size
    double T = 2000.0; // final time
    int N = 900; // number of bodies
    auto force = new double[N][3]; // forces acting on each body

    // Random number generator
    std::default_random_engine engine(
        std::chrono::system_clock::now().time_since_epoch().count()
    );

    auto mass = new double[N];
    double mMean = 1.0, mStdDev = 0.0; // All masses will be 1
    std::normal_distribution<double> massDistribution(mMean, mStdDev);
    for (int n = 0; n < N; n++) {
        mass[n] = massDistribution(engine);
    }

    auto p = new double[N][3]; // positions of the bodies
    double pMean = 0.0, pStdDev = 1.0;
    std::normal_distribution<double> positionDistribution(pMean, pStdDev);
    for (int n = 0; n < N; n++) {
        p[n][0] = positionDistribution(engine);
        p[n][1] = positionDistribution(engine);
        p[n][2] = positionDistribution(engine);
    }

    auto v = new double[N][3]; // velocities of the bodies
    double vMean = 0.0, vStdDev = 1.0;
    std::normal_distribution<double> velocityDistribution(vMean, vStdDev);
    for (int n = 0; n < N; n++) {
        v[n][0] = velocityDistribution(engine);
        v[n][1] = velocityDistribution(engine);
        v[n][2] = velocityDistribution(engine);
    }

    // Create SFML window
    sf::RenderWindow window(sf::VideoMode(800, 600), "N-Body Simulation");

    // Create circles for each body
    sf::CircleShape particleShape(5); // Radius of 5 pixels for each particle
    std::vector<sf::CircleShape> particles(N, particleShape);

    // Assigns random colours to each particle
    for (int i = 0; i < N; i++) {
        particles[i].setFillColor(randomColour());
    }

    // Simulation loop
    while (t < T && window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed) {
                window.close();
            }
        }

        // Zero out the force array
        for (int i = 0; i < N; i++) {
            force[i][0] = 0.0;
            force[i][1] = 0.0;
            force[i][2] = 0.0;
        }

        // Calculate gravitational forces
        for (int i = 0; i < N; i++) {
            for (int j = i + 1; j < N; j++) {
                double Fx, Fy, Fz;
                calcForce(p[i][0], p[i][1], p[i][2], p[j][0], p[j][1], p[j][2],
                          mass[i], mass[j], &Fx, &Fy, &Fz);
                force[i][0] += Fx;
                force[i][1] += Fy;
                force[i][2] += Fz;
                force[j][0] -= Fx;
                force[j][1] -= Fy;
                force[j][2] -= Fz;
            }
        }

        // Forward (Symplectic) Euler
        for (int i = 0; i < N; i++) {
            v[i][0] += force[i][0] / mass[i] * dt;
            v[i][1] += force[i][1] / mass[i] * dt;
            v[i][2] += force[i][2] / mass[i] * dt;

            p[i][0] += v[i][0] * dt;
            p[i][1] += v[i][1] * dt;
            p[i][2] += v[i][2] * dt;
        }

        // Clear the screen
        window.clear();

        // Update positions of bodies in the window
        for (int i = 0; i < N; i++) {
            particles[i].setPosition(p[i][0] * 100 + 400, p[i][1] * 100 + 300);
            window.draw(particles[i]);
        }

        window.display();
        t += dt;
        std::this_thread::sleep_for(std::chrono::milliseconds(10)); // as this number decreases, the animation speeds up
    }

    delete[] mass;
    delete[] force;
    delete[] p;
    delete[] v;
}
