#include <cmath>
#include <iostream>
#include <SDL2/SDL.h>
#include <cstdlib>
#include <vector>
#include <map>

namespace consts {
    const float HMol = 1.008e-3;
    const float HeMol = 4.003e-3;
    const float kb = 1.380649e-23;
    const float Na = 6.02214076e23;
    const float r = kb*Na;
    const float QH_He = 4.3e-12;
    
    const float HCv = 10730;

    const float g = -9.8;
    const double EARTH_ORBIT = 1.5e11;
    const double EARTH_MASS = 6e24;
    const double SUN_MASS = 2e30;
    // const double G = 6.67e-11;
    const double G = 6.67e-11; // in 2D
    const double PI = M_PI;
    // const double GRID_WIDTH = 1.4E9;
    // const double GRID_HEIGHT = 1.4E9;
    const int GRID_WIDTH = 800;
    const int GRID_HEIGHT = 800;
    const float dt = 0.01;
    const int N_ROWS = 100;
    const int N_COLS = 100;
    const float delta_x = GRID_WIDTH / N_COLS;
    const float delta_y = GRID_HEIGHT / N_ROWS;

    // const float N_TAIL = 100;
    const float YEAR = 365.24 * 24 * 60 * 60;
    const float SCALE_TIME = 1;

    
}

// cell properties
const uint MASS = 0;
const uint TEMPERATURE = 1;
const uint PRESSURE = 2;
const uint SMOKEMASS = 3;
const uint E = 4;
const uint R = 5;
const uint G = 6;
const uint B = 7;
const uint GRAVITY = 8;
const uint HYDROGEN = 9;
const uint HELIUM = 10;
// mouse modes
const uint SMOKE = 0;
const uint VELOCITY = 1;
const uint SOURCE = 2;
// source types
const uint SMOKEGUN = 0;
const uint POINTSOURCE = 1;
const uint FAN = 2;

const uint MAXSOURCES = 10;

struct position {
    public:
        float x, y;
};

std::vector<std::string> split(std::string s, const std::string delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    size_t last = 0; 
    size_t next = 0; 
    while ((next = s.find(delimiter, last)) != std::string::npos) {   
        token = s.substr(last, next-last);
        tokens.push_back(token);   
        last = next + delimiter.length(); 
    } 
    tokens.push_back(s.substr(last));

    return tokens;
}

// float coolingFunction(float T) {
//     if (T < 1e4) {
//         return 0; // Negligible atomic cooling below 10,000 K
//     } else if (T < 1e5) {
//         return 1e-22 * T; // Line cooling scales linearly
//     } else if (T < 1e6) {
//         return 1e-23 * T * T; // Bremsstrahlung and line cooling
//     } else {
//         return 1e-24 * sqrt(T); // Bremsstrahlung dominant
//     }
// }

float combinedCoolingFunction(float rho, float T) {
    // Constants
    const float C_mol = 1e-27; // Molecular cooling coefficient
    const float C_dust = 1e-25; // Dust cooling coefficient
    const float C_line = 1e-22; // Atomic line cooling coefficient

    float Lambda_mol = 0.0f;
    float Lambda_line = 0.0f;
    float Lambda_dust = 0.0f;
    float Lambda_brem = 0.0f;

    // Molecular cooling (10 K < T < 2000 K)
    if (T >= 10 && T < 2000) {
        Lambda_mol = C_mol * std::pow(T, 4.5);
    }

    // Atomic line cooling (10,000 K < T < 1,000,000 K)
    if (T >= 1e4 && T < 1e6) {
        Lambda_line = (T < 1e5) ? C_line * T : C_line * std::pow(T, 2);
    }

    // Dust cooling (T < 1000 K)
    if (T < 1000) {
        Lambda_dust = C_dust * std::sqrt(T);
    }

    // Bremsstrahlung cooling
    if (T > 1e6) {
        Lambda_brem = 1e-24 * std::sqrt(T);
    }

    // Total cooling
    return rho * rho * (Lambda_mol + Lambda_line + Lambda_dust);
}

float coolingFunction(float rho, float T, float Z) {
    // Constants
    const float C_H = 7.5e-19;     // Hydrogen cooling coefficient
    const float C_He = 9.1e-27;    // Helium cooling coefficient
    const float C_ff = 1.42e-27;   // Free-free cooling coefficient
    const float C_mol = 1e-27;     // Molecular cooling coefficient
    const float C_dust = 1e-25;    // Dust cooling coefficient
    const float Z_solar = 0.02;    // Solar metallicity

    // Hydrogen and Helium cooling
    float Lambda_H_He = 0.0f;
    if (T > 10000) {
        Lambda_H_He = C_H * exp(-118400 / T) + C_He * exp(-473600 / T);
    }

    // Metal cooling
    float Lambda_metal = 0.0f;
    if (T > 1e4 && T <= 1e5) {
        Lambda_metal = Z / Z_solar * 1e-23 * T;
    } else if (T > 1e5 && T <= 1e6) {
        Lambda_metal = Z / Z_solar * 1e-22;
    } else if (T > 1e6) {
        Lambda_metal = Z / Z_solar * 1e-23 / T;
    }

    // Free-free cooling
    float Lambda_ff = 0.0f;
    if (T > 1e7) {
        Lambda_ff = C_ff * sqrt(T);
    }

    // Molecular cooling
    float Lambda_mol = 0.0f;
    if (T < 2000) {
        Lambda_mol = C_mol * pow(T, 4.5);
    }

    // Dust cooling
    float Lambda_dust = 0.0f;
    if (T < 1000) {
        Lambda_dust = C_dust * sqrt(T);
    }

    // Total net cooling rate
    float Lambda_net = Lambda_H_He + Lambda_metal + Lambda_ff + Lambda_mol + Lambda_dust;

    // Return cooling rate per unit volume
    return rho * rho * Lambda_net; // Cooling rate in erg/cm^3/s
}

float fusionRate(float rho, float T, float X) {
    const float epsilon0 = 1.07e-19; // J m^3 kg^-2 s^-1
    float T6 = T / 1e6; // Temperature in millions of Kelvin
    if (T6 < 1) return 0; // No fusion below 1 million K
    return epsilon0 * rho * X * X * pow(T6, 4);
}

class Neighbors;
class VelocityVector;
class VelocityBox;
class VelocityGrid;
class FluidCell;
class Source;
class FluidGrid {
    public:
        FluidGrid(float width, float height, int r, int c, float dt, bool comp);

        FluidCell *getCell(int i, int j);

        void freeGrid();

        std::map<uint, float> sampleCellAtPoint(float x, float y);
        
        std::map<uint, float> sampleCellAtPoint(position xy);

        std::map<uint, float> propsAtij(int i, int j);

        void force(int i, int j, float fx, float fy);

        void diffuse(int iters);

        std::map<char, float> getxyFromij(int i, int j);

        void massContinuity();

        void solveGravPotential(int iters);

        void compressibleMomentumUpdate();

        void energyUpdate();

        void projection(int iters);

        float calculateFlux(float quantity, int i, int j);

        std::map<uint,float> calculateFlux(std::map<uint,float> propDict, int i, int j, float v);

        void advect();

        void FFSLadvect();

        void setVelocities(bool setE, bool noRecalc);

        void update(SDL_Event event);

        std::vector<FluidCell*> getActive();

        float getdt() {
            return dt;
        }
        unsigned int nActive;
        VelocityGrid *vGrid;
        // if 0, mouse left click adds mass; if 1, mouse left click dragging
        // changes velocities; if 2, mouse adds smokegun; if 3, mouse adds point
        // source; if 4, mouse adds fan
        uint mouseVelFlag = 0;
        int buttonHeld = 0;
        Sint32 prevMouseX = 0;
        Sint32 prevMouseY = 0;
        // so that simulator can display the pressure if this is 1
        uint pressureDisplay = 0; 
        uint temperatureDisplay = 0;
        uint densityDisplay = 1;
        uint gravityFlag = 0;
        uint gravPotentialDisplay = 0;
        uint energyDisplay = 0;
        bool compressible;
    private:
        float width, height;
        float cellWidth, cellHeight;
        float dt, maxV, maxT, totalMass;
        float massFactor = 1;
        int rows, cols, cells;
        FluidCell **grid;
        FluidCell **newGrid;
        std::vector<FluidCell*> activeCells;
        float SCALE_H, SCALE_W;
        std::vector<uint> color = {1,1,0};
        VelocityGrid *new_vGrid;
        Source *sourceArray;
        std::vector<Source*> sourceList;
        float viscosity = 0;
};
  
  class Simulator;


