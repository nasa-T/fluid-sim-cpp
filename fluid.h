#include <cmath>
#include <iostream>
#include <SDL2/SDL.h>
#include <cstdlib>
#include <vector>
#include <map>

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

class Neighbors;
class VelocityVector;
class VelocityBox;
class VelocityGrid;
class FluidCell;
class Source;
class FluidGrid {
    public:
        FluidGrid(float width, float height, int r, int c, float dt);

        FluidCell *getCell(int i, int j);

        void freeGrid();

        std::map<uint, float> sampleCellAtPoint(float x, float y);

        void projection(int iters);

        void advect();

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
    private:
        float width, height;
        float cellWidth, cellHeight;
        float dt, maxV;
        int rows, cols, cells;
        FluidCell **grid;
        FluidCell **newGrid;
        std::vector<FluidCell*> activeCells;
        float SCALE_H, SCALE_W;
        
        VelocityGrid *new_vGrid;
        Source *sourceArray;
        std::vector<Source*> sourceList;
};
  
  class Simulator;


