#include <cmath>
#include <iostream>

namespace consts {
    const float HMol = 2.408e-3;
    const float HeMol = 4e-3;
    const float kb = 1.380649e-23;
    const float Na = 6.02214076e23;
    const float R = kb*Na;

    const double EARTH_ORBIT = 1.5e11;
    const double EARTH_MASS = 6e24;
    const double SUN_MASS = 2e30;
    const double G = 6.67e-11;
  // const double G = 6.67e-2;
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

struct position {
    public:
        float x, y;
};

class Neighbors;
class VelocityVector;
class VelocityBox;
class VelocityGrid;
class FluidCell;
class FluidGrid;
class Simulator;


