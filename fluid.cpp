#include <cmath>
#include <iostream>
#include "fluid.h"

#include <map>
#include <string>


#include <algorithm>

#include <GL/glew.h>

#include <GLFW/glfw3.h>
GLFWwindow* window;

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/norm.hpp>
using namespace glm;

void print(char *string) {
    std::cout << string << std::endl;
}

class VelocityVector {
    public:
        VelocityVector(float vx_ = 0, float vy_ = 0, float ax_ = 0, float ay_ = 0): vx(vx_), vy(vy_), ax(ax_), ay(ay_) {}
        VelocityVector operator+(VelocityVector& other) {
            VelocityVector added(this->vx + other.vx, this->vy + other.vy, other.ax, other.ay);
            return added;
        }
        float getVx() {
            return vx;
        }
        float getVy() {
            return vy;
        }
        float getAx() {
            return ax;
        }
        float getAy() {
            return ay;
        }
        void setVx(float v) {
            vx = v;
        }
        void setVy(float v) {
            vy = v;
        }
        void accelerate() {
            vx += ax * consts::dt;
            vy += ay * consts::dt;
        }
        float getMag() {
            return std::sqrt(getVx()*getVx()+getVy()*getVy());
            // return pow(vx*vx+vy*vy,1/2);
        }
        // std::ostream& operator<<(std::ostream &s, VelocityVector &vec) {
        //     return s << "vx: " << vec.getVx() << " vy: " << vec.getVy();
        // }
    private:
        float vx, vy, ax, ay;
};

class VelocityBox {
    public:
        float top, bottom, left, right;
        VelocityBox(float top = 0, float bottom = 0, float left = 0, float right = 0): top(top), bottom(bottom), left(left), right(right) {}
        void setVelocity(VelocityVector v) {
            left = v.getVx();
            right = v.getVx();
            top = v.getVy();
            bottom = v.getVy();
        }
};

class VelocityGrid {
    public:
        VelocityGrid(float width, float height, int rows, int cols): width(width), height(height), rows(rows), cols(cols) {
            int i, j;
            vyArray = (float **) malloc(sizeof(float *) * (rows + 1));
            vxArray = (float **) malloc(sizeof(float *) * rows);
            for (i = 0; i < rows + 1; i++) {
                vyArray[i] = (float *) malloc(sizeof(float) * cols);
                for (j = 0; j < cols; j++) {
                    vyArray[i][j] = 0;
                }
            }
            for (i = 0; i < rows; i++) {
                vxArray[i] = (float *) malloc(sizeof(float) * (cols + 1));
                for (j = 0; j < cols + 1; j++) {
                    vxArray[i][j] = 0;
                }
            }
            cellHeight = height/rows;
            cellWidth = width/cols;
        }

        VelocityBox getVelocityBox(int i, int j) {
            return VelocityBox(vyArray[i][j], vyArray[i+1][j], vxArray[i][j], vxArray[i][j+1]);
        }
        VelocityVector sampleVelocityAtPoint(float x, float y) {
            int vyi = ceil((height - y) / cellHeight);
            int vyj = floor((x-cellWidth/2) / cellWidth);
            // the physical positions relative to nearest box of vys
            float cvyX, cvyY; 
            float cvxX, cvxY; // same for vxs
            cvyX = (fmod(x+cellWidth/2, cellWidth))/cellWidth;
            cvyY = fmod(y, cellHeight)/cellHeight;
            // will throw an error for when x and y are off the grid entirely
            float vy10, vy00, vy01, vy11 = 0;
            if ((vyj > -1)) {
                vy10 = vyArray[vyi][vyj];
            }
            if ((vyj > -1) && (vyi > 0)) {
                vy00 = vyArray[vyi-1][vyj];
            }
            if ((vyj < cols-1) && (vyi > 0)) {
                vy01 = vyArray[vyi-1][vyj+1];
            }
            if (vyj < cols-1) {
                vy11 = vyArray[vyi][vyj+1];
            }

            float vy = (1-cvyX)*(1-cvyY)*vy10 + cvyX*(1-cvyY)*vy11 + (1-cvyX)*cvyY*vy00 + cvyX*cvyY*vy01;

            int vxi = ceil((height - (y+cellHeight/2)) / cellHeight);
            int vxj = floor(x / cellWidth);

            cvxX = fmod(x, cellWidth)/cellWidth;
            cvxY = (fmod(y+cellHeight/2, cellHeight))/cellHeight;

            float vx10, vx00, vx01, vx11 = 0;
            if ((vxi < rows)) {
                vx10 = vxArray[vxi][vxj];
            }
            if ((vxi > 0)) {
                vx00 = vxArray[vxi-1][vxj];
            }
            if ((vxi > 0)) {
                vx01 = vxArray[vxi-1][vxj+1];
            }
            if (vxi < rows) {
                vx11 = vxArray[vxi][vxj+1];
            }
            float vx = (1-cvxX)*(1-cvxY)*vx10 + cvxX*(1-cvxY)*vx11 + (1-cvxX)*cvxY*vx00 + cvxX*cvxY*vx01;
            return VelocityVector(vx, vy);
        }
        VelocityVector getVelocityVector(int i, int j) {
            return VelocityVector(vxArray[i][j], vyArray[i][j]);
        }
        float getVy(int i, int j) {
            return vyArray[i][j];
        }
        float getVx(int i, int j) {
            return vxArray[i][j];
        }
        void setVy(int i, int j, float vy) {
            vyArray[i][j] = vy;
        }
        void setVx(int i, int j, float vx) {
            vxArray[i][j] = vx;
        }
        float getWidth() {
            return width;
        }
        float getHeight() {
            return height;
        }
        float getCellWidth() {
            return cellWidth;
        }
        float getCellHeight() {
            return cellHeight;
        }
        std::map<char, float> getxyFromijVx(int i, int j) {
            std::map<char, float> xy;
            xy['x'] = cellWidth*j;
            xy['y'] = height - (i*cellHeight + cellHeight/2);
            return xy;
        }
        std::map<char, float> getxyFromijVy(int i, int j) {
            std::map<char, float> xy;
            xy['x'] = cellWidth*j + cellWidth/2;
            xy['y'] = height - cellHeight*i;
            return xy;
        }
        void freeGrid() {
            int i, j;
            for (i = 0; i < rows+1; i++) {
                free(vyArray[i]);
                if (i < rows) free(vxArray[i]);
            }
            free(vyArray);
            free(vxArray);
        }
    private:
        float width, height;
        float cellWidth, cellHeight;
        int rows, cols;
        float **vyArray, **vxArray;
};

class Neighbors {
    public:
        FluidCell *top, *bottom, *left, *right;
        Neighbors() {
            top = NULL;
            bottom = NULL;
            left  = NULL;
            right = NULL;
        }
};


class FluidCell {
    public:
        FluidCell(float mass, float width, float height, float temp, float red, float green, float blue, bool comp): mass(mass), width(width), height(height), temperature(temp), compressible(comp) {
            // std::vector<uint> color
            // size is the physical area
            size = width * height;
            density = mass/size;
            // pressure = 0;
            // pressure = density/consts::HMol * consts::r * temperature;
            float cv = consts::r/consts::HMol/(1.4-1);
            e = cv*temperature*density;
            pressure = (1.4-1)*e;
            // r = color[0];
            // g = color[1];
            // b = color[2];
            // r = red;
            // g = green;
            // b = blue;
            // setColor((!!red) * mass, (!!green) * mass, (!!blue) * mass);
            setColor(red, green, blue);
            velocity = VelocityVector();
            // velocity.setVx(1);
            // velocity.setVx((std::rand() % (int)(2*width)) - width/2);
            // velocity.setVy((std::rand() % (int)(2*height)) - height/2);
            // vBounds = VelocityBox(0,0,0,0);
            // vBounds.setVelocity(velocity);
            neighbors = Neighbors();
        }
        float getCv() {
            return consts::r/consts::HMol/(1.4-1);
        }
        float getMass(bool smoke=true) {
            if (smoke) return smokeMass;
            else return mass;
        }
        float getSize() {
            return size;
        }
        float getDensity() {
            return mass/size;
        }
        float getTemp(bool recalculate=true) {
            if (recalculate) {
                float cv = consts::r/consts::HMol/(1.4-1);
                temperature = gete(false)/(getDensity()*cv);
            }
            return temperature;
        }
        void setTemp(float t) {
            temperature = t;
        }
        float getPressure(bool recalculate=true) {
            // float oldPressure = pressure;
            if (compressible && recalculate) {
                // pressure = getDensity()/consts::HMol * consts::r * getTemp();
                pressure = (1.4 - 1)*gete(false);
            }
            return pressure;
        }
        void setPressure(float p) {
            if (!!pressure) {
                float newTemperature = temperature * p/pressure;
                // printf("newtemp: %f\n",newTemperature);
                // temperature = newTemperature;
            }
            pressure = p;
        }
        void setMass(float m, bool smoke=true) {
            if (m >= 0) {
                if (smoke) smokeMass = m;
                else mass = m;
            }
            // setColor(mass*r, mass*g, mass*b);
        }
        void addMass(float m, bool smoke=true) {
            setMass(mass + m, smoke);
        }
        std::vector<float> getColor() {
            return std::vector<float> {r, g, b};
        }

        float gete(bool recalculate=true) {
            if (recalculate) {
                e = getCv()*getTemp()*getDensity();
                // e = E - 1/2*getDensity()*velocity.getMag()*velocity.getMag();
            }
            return e;
        }

        void sete(float energy) {
            e = energy;
        }
        
        float getE(bool recalculate=false) {
            // if (recalculate) {
            //     E = e + 1/2*getDensity()*velocity.getMag()*velocity.getMag();
            // }
            return E;
        }

        void setE(float energy) {
            E = energy;
        }

        void setColor(float red, float green, float blue, bool source=false) {
            float tot = red + green + blue;
            // if (source) {
            //     r = red/tot;
            //     g = green/tot;
            //     b = blue/tot;
            //     return;
            // }
            // if ((tot == 0) || (mass == 0)) {
            //     r = 0;
            //     g = 0;
            //     b = 0;
            // } else {
            //     r = red/tot;
            //     g = green/tot;
            //     b = blue/tot;
            // }
            r = red;
            g = green;
            b = blue;
        }

        FluidCell *getRight() {
            return neighbors.right;
        }
        FluidCell *getLeft() {
            return neighbors.left;
        }
        FluidCell *getTop() {
            return neighbors.top;
        }
        FluidCell *getBottom() {
            return neighbors.bottom;
        }

        void setRight(FluidCell *neighbor) {
            neighbors.right = neighbor;
        }
        void setLeft(FluidCell *neighbor) {
            neighbors.left = neighbor;
        }
        void setTop(FluidCell *neighbor) {
            neighbors.top = neighbor;
        }
        void setBottom(FluidCell *neighbor) {
            neighbors.bottom = neighbor;
        }

        void transferMass(FluidCell *other, float m) {
            if (other != NULL) { // halt at the boundaries; NOTHING GETS OUT
                if (m > 0) {
                    if (this->getMass() < m) {
                        m = this->getMass();
                    }
                    float otherMass = other->getMass();
                    other->setMass(otherMass + m);
                    this->setMass(mass - m);
                } else if (m < 0) {
                    m = -m;
                    if (other->getMass() < m) {
                        m = other->getMass();
                    }
                    float otherMass = other->getMass();
                    other->setMass(otherMass - m);
                    this->setMass(mass + m);
                }
            }
        }
        // void addVelocity(VelocityVector vel2) {
        //     vel = vel + vel2;
        // }
        void setVelocity(float vx, float vy) {
            velocity.setVx(vx);
            velocity.setVy(vy);
            // setE(e + 1/2*getDensity()*velocity.getMag()*velocity.getMag());
            // getE(true);
            // velocity = vel2;
            // vBounds.setVelocity(vel2);
        }
        void setVelocity(VelocityVector v) {
            velocity = v;
            // setE(e + 1/2*getDensity()*velocity.getMag()*velocity.getMag());
            // getE(true);
        }
        VelocityVector getVelocity() {
            return velocity;
        }

        // void solveVelocity() {
        //     // pressure difference first
            
        // }

        char isActive() {
            char b, t, l, r = 0;
            if (neighbors.bottom != NULL) b = !!(neighbors.bottom->getMass());
            if (neighbors.top != NULL) t = !!(neighbors.top->getMass());
            if (neighbors.left != NULL) l = !!(neighbors.left->getMass());
            if (neighbors.right != NULL) r = !!(neighbors.right->getMass());
            // printf("%d ", !!mass || b || t || l || r);
            return !!mass || b || t || l || r;
        }
        
        void setLoc(int r, int c) {
            row = r;
            col = c;
        }

        int getRow() {
            return row;
        }

        int getCol() {
            return col;
        }

    private:
        float mass, size, density, temperature, pressure, width, height, smokeMass, e, E;
        bool compressible;
        int row, col;
        VelocityVector velocity;
        // VelocityBox vBounds;
        Neighbors neighbors;
        float r,g,b;
        // std::vector<uint> color;
};
// class FluidGrid {
//     public:
//         FluidGrid(float width, float height, int r, int c): width(width), height(height), rows(r), cols(c) {}
//         FluidCell *getCell(int i, int j) {
//             return &grid[i][j];
//         }
//         private:
//             float width, height;
//             float cellWidth, cellHeight;
//             float dt, maxV;
//             int rows, cols, cells;
//             FluidCell **grid;
//             FluidCell **newGrid;
//             std::vector<FluidCell*> activeCells;
//             float SCALE_H, SCALE_W;
            
//             VelocityGrid *new_vGrid;
//             std::vector<Source*> sourceList;
// };
class Source {
    public:
        Source(float _x, float _y, FluidGrid *gd, VelocityGrid *vgrid, uint sourceType, std::vector<uint> color):  grid(gd), vGrid(vgrid), type(sourceType) {
            r = color[0];
            g = color[1];
            b = color[2];
            // x(_x), y(_y),
            // int i, j;
            // xs = _x;
            // ys = _y;
            setX(_x);
            setY(_y);
            switch (sourceType) {
                case SMOKEGUN:
                    i = floor((vgrid->getHeight() - ys) / vgrid->getCellHeight());
                    j = floor(xs / vgrid->getCellWidth());
                    grid->getCell(i,j)->addMass(100);
                    grid->getCell(i,j)->setColor(r,g,b);
                    break;
                case FAN:
                    break;
                case POINTSOURCE:
                    i = floor((vgrid->getHeight() - ys) / vgrid->getCellHeight());
                    j = floor(xs / vgrid->getCellWidth());
                    grid->getCell(i,j)->addMass(100);
                    grid->getCell(i,j)->setColor(r,g,b);
                    setVx(vgrid->getWidth()/4);
                    setVy(vgrid->getHeight()/4);
                    // vGrid->setVx(i,j,-vgrid->getWidth()/4);
                    // vGrid->setVx(i,j+1,vgrid->getWidth()/4);
                    // vGrid->setVy(i,j,vgrid->getWidth()/4);
                    // vGrid->setVy(i+1,j,-vgrid->getWidth()/4);
                    break;
            }
            // printf("source created\n");
            // printf("source loc: %f, %f\n", getX(), getY());
            // vgrid->setVx(i,j,vx);
            // vgrid->setVx(i,j+1,vx);
            // vgrid->setVy(i,j,vy);
            // vgrid->setVy(i+1,j,vy);
        }
        void setX(float x) {
            xs = x;
        }
        void setY(float y) {
            ys = y;
        }
        float getX() {
            return xs;
        }
        float getY() {
            return ys;
        }
        float getVx() {
            return vx;
        }
        float getVy() {
            return vy;
        }
        void setVx(float v) {
            vx = v;
            float mass = grid->getCell(i,j)->getMass(false);
            if (type == SMOKEGUN) {
                // if (v > 0) {
                //     grid->force(i,j+1,vx*100,0);
                // } else {
                //     grid->force(i,j,vx*100,0);
                // }
                grid->force(i,j,vx*100*mass,0);
                // vGrid->setVx(i,j,vx);
            } else if (type == POINTSOURCE) {
                vGrid->setVx(i,j,-vx);
            }
            // vGrid->setVx(i,j+1,vx);
        }
        void setVy(float v) {
            vy = v;
            // vGrid->setVy(i,j,vy);
            float mass = grid->getCell(i,j)->getMass(false);
            if (type == SMOKEGUN) {
                // if (v > 0) {
                //     grid->force(i,j,0,vy*100);
                // } else {
                //     grid->force(i+1,j,0,vy*100);
                // }
                grid->force(i,j,0,vy*100*mass);
                // vGrid->setVy(i+1,j,vy);
            } else if (type == POINTSOURCE) {
                vGrid->setVy(i+1,j,-vy);
            }
        }
        void addMass() {
            float v = std::sqrt(vx*vx+vy*vy);
            float cellSize = vGrid->getCellHeight() * vGrid->getCellWidth();
            grid->getCell(i,j)->setColor(r,g,b,true);
            if (getType() == POINTSOURCE) {
                grid->getCell(i,j)->addMass(0.5*v/cellSize);
            } else if (getType() == SMOKEGUN) {
                grid->getCell(i,j)->addMass(20*v*v/cellSize);
            }
            
        }
        uint getType() {
            return type;
        }
    private:
        uint r,g,b;
        uint type;
        FluidGrid *grid;
        VelocityGrid *vGrid;
        float vx, vy = 0;
        float xs, ys;
        int i, j;
};

FluidGrid::FluidGrid(float width, float height, int r, int c, float dt, bool comp): width(width), height(height), rows(r), cols(c), dt(dt), compressible(comp) {
    int i, j;
    cells = rows * cols;
    grid = (FluidCell **) malloc(sizeof(FluidCell*)*r);
    newGrid = (FluidCell **) malloc(sizeof(FluidCell*)*r);
    vGrid = (VelocityGrid *) malloc(sizeof(VelocityGrid));
    vGrid[0] = VelocityGrid(width, height, r, c);
    new_vGrid = (VelocityGrid *) malloc(sizeof(VelocityGrid));
    new_vGrid[0] = VelocityGrid(width, height, r, c);
    sourceArray = (Source *) malloc(MAXSOURCES*sizeof(Source));
    cellWidth = width / c;
    cellHeight = height / r;
    // in meters per pixel
    SCALE_H = height / consts::GRID_HEIGHT;
    SCALE_W = width / consts::GRID_WIDTH;
    maxV = std::max(cellHeight/5, cellWidth/5);
    dt = 0.016;
    for (i = 0; i < rows; i++) {
        grid[i] = (FluidCell *) malloc(sizeof(FluidCell)*c);
        newGrid[i] = (FluidCell *) malloc(sizeof(FluidCell)*c);
        for (j = 0; j < cols; j++) {
            //random mass for now
            // float mass = std::rand() % 256;
            std::vector<uint> color = {0, 0, 0};
            grid[i][j] = FluidCell(1000, cellWidth, cellHeight, 400, 0,0,0,comp);
            // if (i == rows/2) {
            //     grid[i][j] = FluidCell(1000, cellWidth, cellHeight, 300, 0,0,0);
            // } else if (i == rows-1) {
            //     grid[i][j] = FluidCell(1000, cellWidth, cellHeight, 300, 0,0,0);
            // }
            // } else if (i == rows/2 && j == cols/2+1) {
            //     grid[i][j] = FluidCell(0, cellWidth, cellHeight, 100);
            //     // grid[i][j].setVelocity(VelocityVector(10,0));
            // } else {
            //     grid[i][j] = FluidCell(0, cellWidth, cellHeight, 100);
            // }
            if (i != 0) {
                grid[i][j].setTop(&grid[i - 1][j]);
                grid[i - 1][j].setBottom(&grid[i][j]);
            }
            if (j != 0) {
                grid[i][j].setLeft(&grid[i][j - 1]);
                grid[i][j - 1].setRight(&grid[i][j]);
            }
            // if (grid[i][j].getVelocity().getMag() > maxV) {
            //     maxV = grid[i][j].getVelocity().getMag();
            // }
            grid[i][j].setLoc(i, j);
            newGrid[i][j] = grid[i][j];
            // new_vGrid->setVx(i,j,vGrid->getVx(i,j));
            // new_vGrid->setVy(i,j,vGrid->getVy(i,j));
            new_vGrid->setVx(i,j,0);
            new_vGrid->setVy(i,j,0);
        }
    }
    setVelocities(true,false);
    // 0.7/maxV;
    nActive = activeCells.size();
}

FluidCell *FluidGrid::getCell(int i, int j) {
    return &grid[i][j];
}

void FluidGrid::freeGrid() {
    int i, j;
    for (i = 0; i < rows; i++) {
        free(grid[i]);
    }
    vGrid->freeGrid();
    free(vGrid);
    free(grid);
}

std::map<uint, float> FluidGrid::sampleCellAtPoint(float x, float y) {
    // y = round(y);
    // float cHeight = round(cellHeight);
    // printf("%f\n", round(height - y - cellHeight/2));
    int i = ceil((height - y - cellHeight/2) / cellHeight);
    int j = floor(round(x - cellWidth/2) / cellWidth);
    if ((i == rows/2) && (j == cols/2)) {
        // printf("%d: %f %d: %f\n", j,x, i,y);
        // printf("%f %f\n", prevY, physY);
    }
    // switch (prop) {
    //     case MASS:

    // }

    // the physical positions relative to nearest box of cells
    float X, Y; 
    X = fmod(x+cellWidth/2, cellWidth)/cellWidth;
    Y = fmod(y+cellHeight/2, cellHeight)/cellHeight;

    // will throw an error for when x and y are off the grid entirely
    std::map<uint, float> p10, p00, p01, p11;
    if ((j > -1) && (i < rows)) {
        // printf("%d %d\n", j, i);
        // p10 = vyArray[i][j];
        FluidCell *cell = getCell(i,j);
        p10[SMOKEMASS] = cell->getMass();
        p10[MASS] = cell->getMass(false);
        p10[TEMPERATURE] = cell->getTemp();
        p10[PRESSURE] = cell->getPressure();
        p10[E] = cell->getE();
        // p10[R] = cell->getColor()[0]*p10[MASS];
        // p10[G] = cell->getColor()[1]*p10[MASS];
        // p10[B] = cell->getColor()[2]*p10[MASS];
        p10[R] = cell->getColor()[0];
        p10[G] = cell->getColor()[1];
        p10[B] = cell->getColor()[2];
    } else {
        p10[SMOKEMASS] = 0;
        p10[MASS] = 0;
        p10[TEMPERATURE] = 0;
        p10[PRESSURE] = 0;
        p10[E] = 0;
        p10[R] = 0;
        p10[G] = 0;
        p10[B] = 0;
    }
    if ((j > -1) && (i > 0)) {
        // p00 = vyArray[i-1][j];
        FluidCell *cell = getCell(i-1,j);
        p00[SMOKEMASS] = cell->getMass();
        p00[MASS] = cell->getMass(false);
        p00[TEMPERATURE] = cell->getTemp();
        p00[PRESSURE] = cell->getPressure();
        p00[E] = cell->getE();
        // p00[R] = cell->getColor()[0]*p00[MASS];
        // p00[G] = cell->getColor()[1]*p00[MASS];
        // p00[B] = cell->getColor()[2]*p00[MASS];
        p00[R] = cell->getColor()[0];
        p00[G] = cell->getColor()[1];
        p00[B] = cell->getColor()[2];
    } else {
        p00[SMOKEMASS] = 0;
        p00[MASS] = 0;
        p00[TEMPERATURE] = 0;
        p00[PRESSURE] = 0;
        p00[E] = 0;
        p00[R] = 0;
        p00[G] = 0;
        p00[B] = 0;
    }
    if ((j < cols-1) && (i > 0)) {
        // p01 = vyArray[i-1][j+1];
        FluidCell *cell = getCell(i-1,j+1);
        p01[SMOKEMASS] = cell->getMass();
        p01[MASS] = cell->getMass(false);
        p01[TEMPERATURE] = cell->getTemp();
        p01[PRESSURE] = cell->getPressure();
        p01[E] = cell->getE();
        // p01[R] = cell->getColor()[0]*p01[MASS];
        // p01[G] = cell->getColor()[1]*p01[MASS];
        // p01[B] = cell->getColor()[2]*p01[MASS];
        p01[R] = cell->getColor()[0];
        p01[G] = cell->getColor()[1];
        p01[B] = cell->getColor()[2];
    } else {
        p01[SMOKEMASS] = 0;
        p01[MASS] = 0;
        p01[TEMPERATURE] = 0;
        p01[PRESSURE] = 0;
        p01[E] = 0;
        p01[R] = 0;
        p01[G] = 0;
        p01[B] = 0;
    }
    if ((j < cols-1) && (i < rows)) {
        // p11 = vyArray[i][j+1];
        FluidCell *cell = getCell(i,j+1);
        p11[SMOKEMASS] = cell->getMass();
        p11[MASS] = cell->getMass(false);
        p11[TEMPERATURE] = cell->getTemp();
        p11[PRESSURE] = cell->getPressure();
        p11[E] = cell->getE();
        // p11[R] = cell->getColor()[0]*p11[MASS];
        // p11[G] = cell->getColor()[1]*p11[MASS];
        // p11[B] = cell->getColor()[2]*p11[MASS];
        p11[R] = cell->getColor()[0];
        p11[G] = cell->getColor()[1];
        p11[B] = cell->getColor()[2];
    } else {
        p11[SMOKEMASS] = 0;
        p11[MASS] = 0;
        p11[TEMPERATURE] = 0;
        p11[PRESSURE] = 0;
        p11[E] = 0;
        p11[R] = 0;
        p11[G] = 0;
        p11[B] = 0;
    }
    // printf("%f %f %f %f\n", p10[MASS], p00[MASS], p01[MASS], p11[MASS]);
    std::map<uint, float> props;
    props[SMOKEMASS] = (1-X)*(1-Y)*p10[SMOKEMASS] + X*(1-Y)*p11[SMOKEMASS] + (1-X)*Y*p00[SMOKEMASS] + X*Y*p01[SMOKEMASS];
    
    props[MASS] = (1-X)*(1-Y)*p10[MASS] + X*(1-Y)*p11[MASS] + (1-X)*Y*p00[MASS] + X*Y*p01[MASS];
    // printf("%f\n", props[MASS]);
    props[TEMPERATURE] = (1-X)*(1-Y)*p10[TEMPERATURE] + X*(1-Y)*p11[TEMPERATURE] + (1-X)*Y*p00[TEMPERATURE] + X*Y*p01[TEMPERATURE];
    props[PRESSURE] = (1-X)*(1-Y)*p10[PRESSURE] + X*(1-Y)*p11[PRESSURE] + (1-X)*Y*p00[PRESSURE] + X*Y*p01[PRESSURE];
    props[E] = (1-X)*(1-Y)*p10[E] + X*(1-Y)*p11[E] + (1-X)*Y*p00[E] + X*Y*p01[E];
    props[R] = (1-X)*(1-Y)*p10[R] + X*(1-Y)*p11[R] + (1-X)*Y*p00[R] + X*Y*p01[R];
    props[G] = (1-X)*(1-Y)*p10[G] + X*(1-Y)*p11[G] + (1-X)*Y*p00[G] + X*Y*p01[G];
    props[B] = (1-X)*(1-Y)*p10[B] + X*(1-Y)*p11[B] + (1-X)*Y*p00[B] + X*Y*p01[B];
    // props[R] = getCell(i,j)->getColor()[0];
    // props[G] = getCell(i,j)->getColor()[1];
    // props[B] = getCell(i,j)->getColor()[2];
    // if ((i == rows/2) && (j == cols/2))
    // printf("%f %f %f %f\n", p10[R],p00[R],p01[R],p11[R]);

    return props;
}

std::map<char, float> FluidGrid::getxyFromij(int i, int j) {
    std::map<char, float> xy;
    xy['x'] = j*cellWidth + cellWidth/2;
    xy['y'] = height - (i*cellHeight + cellHeight/2);
    return xy;
}

void FluidGrid::diffuse(int iters) {
    if (viscosity > 0) {
        // diffuse velocities, densities, and temperatures using gauss-seidel
        int i, j;
        float x, y;
        VelocityVector c, b, t, l, r; // center bottom top left right
        std::map<char, float> cp, bp, tp, lp, rp;
        float vxC, vxB, vxT, vxR, vxL;
        float vyC, vyB, vyT, vyR, vyL;
        for (int n = 0; n < iters; n++) {
            for (i = 0; i < rows; i++) {
                for (j = 0; j < cols; j++) {
                    // x = j*cellWidth + cellWidth/2;
                    // y = height - (i*cellHeight + cellHeight/2);
                    cp = getxyFromij(i,j);
                    c = vGrid->sampleVelocityAtPoint(cp['x'],cp['y']);
                    if (i - 1 == -1) {
                        t = c;
                    } else {
                        tp = getxyFromij(i-1,j);
                        t = vGrid->sampleVelocityAtPoint(tp['x'],tp['y']);
                    }
                    if (j - 1 == -1) {
                        l = c;
                    } else {
                        lp = getxyFromij(i,j-1);
                        l = vGrid->sampleVelocityAtPoint(lp['x'],lp['y']);
                    }
                    if (i + 1 == rows) {
                        b = c;
                    } else {
                        bp = getxyFromij(i+1,j);
                        b = vGrid->sampleVelocityAtPoint(bp['x'],bp['y']);
                    }
                    if (j + 1 == cols) {
                        r = c;
                    } else {
                        rp = getxyFromij(i,j+1);
                        r = vGrid->sampleVelocityAtPoint(rp['x'],rp['y']);
                    }
                    // solving vxs:
                    vxC = c.getVx();
                    vxB = b.getVx();
                    vxT = t.getVx();
                    vxR = r.getVx();
                    vxL = l.getVx();
                    // float dvx = dt*viscosity*((vxL + vxR - 2*vxC)/(cellWidth*cellWidth)+(vxB + vxT - 2*vxC)/(cellHeight*cellHeight));
                    float vxNew = (vxC+viscosity*((vxL+vxT+vxR+vxB)/4))/(1+viscosity);

                    // solving vys:
                    vyC = c.getVy();
                    vyB = b.getVy();
                    vyT = t.getVy();
                    vyR = r.getVy();
                    vyL = l.getVy();
                    float vyNew = (vyC+viscosity*((vyL+vyT+vyR+vyB)/4))/(1+viscosity);
                    if (i == 0) {
                        vGrid->setVy(i,j,0);
                    } else {
                        vGrid->setVy(i,j,vyNew);
                    }
                    if (i == rows-1) {
                        vGrid->setVy(i+1,j,0);
                    } else {
                        vGrid->setVy(i+1,j,vyNew);
                    }

                    if (j == 0) {
                        vGrid->setVx(i,j,0);
                    } else {
                        vGrid->setVx(i,j,vxNew);
                    }
                    if (j == cols-1) {
                        vGrid->setVx(i,j+1,0);
                    } else {
                        vGrid->setVx(i,j+1,vxNew);
                    }
                    // vGrid->setVx(i,j,vxNew);
                    // vGrid->setVx(i,j+1,vxNew);
                    // vGrid->setVy(i,j,vyNew);
                    // vGrid->setVy(i+1,j,vyNew);
                    // if (round(vyNew) != 0) printf("%d %d, %f\n", i,j,vyNew);
                }
            }
        }
    }
}

void FluidGrid::massContinuity() { // handled by advection, maybe this function can be reused for something else
    int i, j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            std::map<char,float> xyLeft = vGrid->getxyFromijVx(i,j);
            std::map<char,float> xyTop = vGrid->getxyFromijVy(i,j);
            std::map<char,float> xyRight = vGrid->getxyFromijVx(i,j+1);
            std::map<char,float> xyBottom = vGrid->getxyFromijVy(i+1,j);
            float massL = sampleCellAtPoint(xyLeft['x'],xyLeft['y'])[MASS];
            float massT = sampleCellAtPoint(xyTop['x'],xyTop['y'])[MASS];
            float massR = sampleCellAtPoint(xyRight['x'],xyRight['y'])[MASS];
            float massB = sampleCellAtPoint(xyBottom['x'],xyBottom['y'])[MASS];
            float vxL = vGrid->getVx(i,j);
            float vyT = vGrid->getVy(i,j);
            float vxR = vGrid->getVx(i,j+1);
            float vyB = vGrid->getVy(i+1,j);
            float divMass = (-vxL*massL + vyT*massT + vxR*massR - vyB*massB);
            FluidCell *c = getCell(i,j);
            c->setMass(c->getMass(false) - divMass*dt,false);
        }
    }
}

void FluidGrid::compressibleMomentumUpdate() {
    int i, j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            FluidCell *cell = getCell(i,j);
            float mass = cell->getMass(false);
            if (gravityFlag) force(i,j,0,mass*(consts::g));
            // velocities surrounding this cell
            float left = vGrid->getVx(i, j);
            float right = vGrid->getVx(i, j+1);
            float top = vGrid->getVy(i, j);
            float bottom = vGrid->getVy(i+1, j);

            float pressure = cell->getPressure();

            float E = cell->getE();
            // // pressure-work term
            // float div = (right-left)/cellWidth + (top-bottom)/cellHeight;
            // float pdiv = pressure * div;
            // // if (i==1 && j==cols/2) printf("pdiv:%f\n",pdiv);
            // E += -pdiv*getdt();
            
            // gravitational-work term
            // if (gravityFlag) {
            //     float gravWork = cell->getVelocity().getVy()*(consts::g)*cell->getDensity();
            //     E += gravWork*getdt();
            // }
            // cell->setE(E);

            // viscosity term
            std::map<char,float> xyM, xyB, xyL, xyT, xyR;
            float laplacianVx, laplacianVy;
            if ((i == 0) || (i == rows-1) || (j == 0) || (j == cols)) {
                laplacianVx = 0;
            } else {
                laplacianVx = (vGrid->getVx(i,j-1) - 2*vGrid->getVx(i,j) + vGrid->getVx(i,j+1))/(cellWidth*cellWidth) + (vGrid->getVx(i-1,j) - 2*vGrid->getVx(i,j) + vGrid->getVx(i+1,j))/(cellHeight*cellHeight);
            }
            if ((i == 0) || (i == rows) || (j == 0) || (j == cols-1)) {
                laplacianVy = 0;
            } else {
                laplacianVy = (vGrid->getVy(i,j-1) - 2*vGrid->getVy(i,j) + vGrid->getVy(i,j+1))/(cellWidth*cellWidth) + (vGrid->getVy(i-1,j) - 2*vGrid->getVy(i,j) + vGrid->getVy(i+1,j))/(cellHeight*cellHeight);
            }
            left += laplacianVx*viscosity*getdt();
            top += laplacianVy*viscosity*getdt();
            new_vGrid->setVx(i,j,left);
            new_vGrid->setVy(i,j,top);
            if (i < rows-1) {
                float pressureB = getCell(i+1,j)->getPressure();
                float gradY = (pressure-pressureB)/(2*cellHeight);

                // if (i == 1 && j ==cols/2) printf("gradY:%f dens:%f press:%f\n",gradY, cell->getDensity(), cell->getPressure());
                std::map<char,float> xyVy = vGrid->getxyFromijVy(i+1,j);
                float density = sampleCellAtPoint(xyVy['x'],xyVy['y'])[MASS]/(cellHeight*cellWidth);
                // float newBVy = bottom-gradY;
                float newBVy = bottom-(1/density)*gradY*dt;
                // if (i == rows-2 && j ==cols/2) printf("newVy:%f; gradY:%f; 1/dens:%f\n",newBVy,gradY,1/density);
                new_vGrid->setVy(i+1,j,newBVy);
            }
            if (j < cols-1) {
                float pressureR = getCell(i,j+1)->getPressure();
                float gradX = (pressureR-pressure)/(2*cellWidth);
                if (!std::isfinite(gradX)) {
                    std::cerr << "NaN detected in grad X at cell (" << i << ", " << j << ") pressure(i,j): " << pressure << " density: " << cell->getDensity() << "\n";
                    assert(false); // Debug immediately
                }
                // if (i == rows/2 && j ==cols/2) printf("gradX:%f\n",gradX);
                std::map<char,float> xyVx = vGrid->getxyFromijVx(i,j+1);
                float density = sampleCellAtPoint(xyVx['x'],xyVx['y'])[MASS]/(cellHeight*cellWidth);
                float newRVx = right-(1/density)*gradX*dt;
                if (!std::isfinite(newRVx)) {
                    std::cerr << "NaN detected in vxR at cell (" << i << ", " << j << ") old vxR: " << right << " density: " << density << "\n";
                    assert(false); // Debug immediately
                }
                new_vGrid->setVx(i,j+1,newRVx);
                // if (i == rows/2 && j ==cols/2) printf("newVx:%f\n",newRVx);
            }
        }
    }
    vGrid = new_vGrid;
}

void FluidGrid::energyUpdate() {
    int i, j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            FluidCell *cell = getCell(i,j);
            float mass = cell->getMass(false);
            // velocities surrounding this cell
            float left = vGrid->getVx(i, j);
            float right = vGrid->getVx(i, j+1);
            float top = vGrid->getVy(i, j);
            float bottom = vGrid->getVy(i+1, j);

            float pressure = cell->getPressure();
            float E = cell->getE();
            // pressure-work term
            float div = (right-left)/cellWidth + (top-bottom)/cellHeight;
            float pdiv = pressure * div;
            if (E - pdiv * getdt() >= 0) {
                E += -pdiv * getdt();
            } else {
                std::cerr << "Warning: Pressure work term causing negative energy.\n";
                E = 0.0f; // Prevent negative energy
            }
            // if (i==1 && j==cols/2) printf("pdiv:%f\n",pdiv);
            // E += -pdiv*getdt();
            
            // gravitational-work term
            if (gravityFlag) {
                float gravWork = cell->getVelocity().getVy()*(consts::g)*cell->getDensity();
                // E += gravWork*getdt();
                if (E + gravWork * getdt() >= 0) {
                    E += gravWork * getdt();
                } else {
                    std::cerr << "Warning: Gravitational work term causing negative energy.\n";
                    E = 0.0f; // Prevent negative energy
                }
            }
            cell->setE(E);
        }
    }
}

void FluidGrid::projection(int iters) {
    // for incompressible fluids only
    int i, j;
    for (int n = 0; n < iters; n++) {
        for (i = 0; i < rows; i++) {
            for (j = 0; j < cols; j++) {
                FluidCell *cell = getCell(i,j);
                float mass = cell->getMass(false);
                if (n == 0) {
                    if (gravityFlag) force(i,j,0,mass*(consts::g));
                    // float pressure = cell->getDensity()/consts::Hmol * consts::kb * consts::Na * cell->getTemp();
                    cell->setPressure(0);
                }
                int sLeft = !!j;
                int sRight = !!(cols - 1 - j);
                int sTop = !!i;
                int sBottom = !!(rows - 1 - i);                        

                int s = sLeft + sRight + sTop + sBottom;

                float left = vGrid->getVx(i, j);
                float right = vGrid->getVx(i, j+1);
                float top = vGrid->getVy(i, j);
                float bottom = vGrid->getVy(i+1, j);
                
                float d = (-left*sLeft + right*sRight + top*sTop - bottom*sBottom)*1;
                // if (i == rows/2 && j ==cols/2) printf("div:%f\n", d);
                float newPressure = cell->getPressure() - d/s * cell->getDensity()*std::sqrt(cellWidth*cellHeight)/dt;
                cell->setPressure(newPressure);

                
                vGrid->setVx(i, j, left+d*sLeft/s);
                vGrid->setVx(i, j+1, right-d*sRight/s);
                vGrid->setVy(i, j, top-d*sTop/s);
                vGrid->setVy(i+1, j, bottom+d*sBottom/s);

            }
        }
    }
}

void FluidGrid::advect() {
    int i, j;
    for (i = 0; i < rows+1; i++) {
        for (j = 0; j < cols+1; j++) {
            if ((i < rows) && (j < cols)) {
                FluidCell cell = *getCell(i,j);
                if (cell.getDensity() == 0) {
                    printf("%d %d\n", i,j);
                }
                float physX = cellWidth*j+cellWidth/2;
                float physY = height - (cellHeight*i+cellHeight/2); // I want y pointed up
                VelocityVector vCell = vGrid->sampleVelocityAtPoint(physX, physY);

                float prevX = physX - dt*vCell.getVx();
                float prevY = physY - dt*vCell.getVy();

                std::map<uint, float> prevProps;
                // if (round((vCell.getVx()+vCell.getVy())*dt/std::sqrt(cellHeight*cellWidth)) == 0) {
                if (round((vCell.getVx()+vCell.getVy())) == 0) {
                    prevProps[SMOKEMASS] = cell.getMass();
                    prevProps[MASS] = cell.getMass(false);
                    prevProps[TEMPERATURE] = cell.getTemp();
                    prevProps[PRESSURE] = cell.getPressure();
                    prevProps[E] = cell.getE();
                    prevProps[R] = cell.getColor()[0];
                    prevProps[G] = cell.getColor()[1];
                    prevProps[B] = cell.getColor()[2];
                } else {
                    prevX = std::max(0.0f, std::min(prevX, width));
                    prevY = std::max(0.0f, std::min(prevY, height));
                    prevProps = sampleCellAtPoint(prevX, prevY);
                    // if ((0 < prevX) && (prevX < width) && (0 < prevY) && (prevY < height)) {
                        // prevProps = sampleCellAtPoint(prevX, prevY);
                    // } 
                }
                
                if (prevProps[E] < 0) {
                    std::cerr << "Warning: Negative energy detected in advection.\n";
                    prevProps[E] = 0.0f;
                }
                if (prevProps[MASS] == 0.0f) {
                    std::cerr << "Warning: 0 mass detected in advection.\n";
                    printf("prevX: %f prevY: %f\n", prevX,prevY);
                    prevProps[MASS] = 0.0f;
                }
                
                newGrid[i][j].setMass(prevProps[SMOKEMASS]);
                newGrid[i][j].setMass(prevProps[MASS],false);
                newGrid[i][j].setTemp(prevProps[TEMPERATURE]);
                
                newGrid[i][j].setE(prevProps[E]);
                newGrid[i][j].setColor(prevProps[R], prevProps[G], prevProps[B]);
                if (!std::isfinite(prevProps[E])) {
                    std::cerr << "NaN detected in energy at cell (" << i << ", " << j << ") " << "density: " << cell.getDensity() << "\n";
                    assert(false); // Debug immediately
                }
                // if ((prevX < 0) || (prevY < 0) || (prevX > width) || (prevY > height)) {
                //     newGrid[i][j].setMass(0);
                //     newGrid[i][j].setMass(0,false);
                //     newGrid[i][j].setTemp(0);
                //     newGrid[i][j].setE(0);
                //     newGrid[i][j].setColor(0,0,0);
                //     // newGrid[i][j].setMass(cell.getMass());
                //     // newGrid[i][j].setMass(cell.getMass(false),false);
                //     // newGrid[i][j].setTemp(cell.getTemp());
                //     // newGrid[i][j].setE(cell.getE());
                //     // std::vector<float> c = cell.getColor();
                //     // newGrid[i][j].setColor(c[0],c[1],c[2]);
                // } else {
                //     newGrid[i][j].setMass(prevProps[SMOKEMASS]);
                //     newGrid[i][j].setMass(prevProps[MASS],false);
                //     newGrid[i][j].setTemp(prevProps[TEMPERATURE]);
                //     newGrid[i][j].setE(prevProps[E]);
                //     newGrid[i][j].setColor(prevProps[R], prevProps[G], prevProps[B]);
                // }
            }

            float newVx, newVy;
            // velocity self-advection:
            if (i < rows) {
                // vxs of vgrid
                float vxPhysX = cellWidth*j;
                float vxPhysY = height - cellHeight*i - cellHeight/2;
                
                VelocityVector vVx = vGrid->sampleVelocityAtPoint(vxPhysX, vxPhysY);
                float prevVxX = vxPhysX - dt*vVx.getVx();
                float prevVxY = vxPhysY - dt*vVx.getVy();
                prevVxX = std::max(0.0f, std::min(prevVxX, width - 1));
                prevVxY = std::max(0.0f, std::min(prevVxY, height - 1));
                if ((0 < prevVxX) && (prevVxX < width) && (0 < prevVxY) && (prevVxY < height)) {
                    newVx = vGrid->sampleVelocityAtPoint(prevVxX, prevVxY).getVx();
                    if ((j == 0) || (j == cols)) {
                        new_vGrid->setVx(i,j,0);
                    } else {
                        new_vGrid->setVx(i,j,newVx);
                    }

                    // find the max velocity in x direction for dt
                    if (newVx > maxV) {
                        maxV = newVx;
                    }
                }
            }
            
            if (j < cols) {
                // vys of vgrid
                float vyPhysX = cellWidth*j + cellWidth/2;
                float vyPhysY = height - cellHeight*i;
                VelocityVector vVy = vGrid->sampleVelocityAtPoint(vyPhysX, vyPhysY);
                float prevVyX = vyPhysX - dt*vVy.getVx();
                float prevVyY = vyPhysY - dt*vVy.getVy();
                prevVyX = std::max(0.0f, std::min(prevVyX, width - 1));
                prevVyY = std::max(0.0f, std::min(prevVyY, height - 1));
                if ((0 < prevVyX) && (prevVyX < width) && (0 < prevVyY) && (prevVyY < height)) {
                    newVy = vGrid->sampleVelocityAtPoint(prevVyX, prevVyY).getVy();
                    if ((i == 0) || (i == rows)) {
                        new_vGrid->setVy(i,j,0);
                    } else {
                        new_vGrid->setVy(i,j,newVy);
                    }
                    // find the max velocity in y direction for dt
                    if (newVy > maxV) {
                        maxV = newVy;
                    }
                }
            }
        }
    }
    grid = newGrid;
    vGrid = new_vGrid;
}

void FluidGrid::setVelocities(bool setE=false, bool noRecalc=false) {
    // so that each cell knows about its velocity
    float totMass = 0;
    int i,j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            std::map<char, float> xy = getxyFromij(i,j);
            VelocityVector v = vGrid->sampleVelocityAtPoint(xy['x'],xy['y']);
            FluidCell *cell = getCell(i,j);
            cell->setVelocity(v);
            if (setE && !noRecalc) {
                cell->setE(cell->gete()+0.5*cell->getDensity()*v.getMag()*v.getMag());
                // if (i==1 && j ==cols/2) printf("initial V: %f\n",v.getMag());
            } else if (!setE && !noRecalc) {
                float newE = cell->getE();
                float KE = 0.5 * cell->getDensity() * (v.getMag() * v.getMag());
                if (newE >= KE) {
                    cell->sete(newE - KE);
                } else {
                    std::cerr << "Warning: Total energy less than kinetic energy.\n";
                    cell->sete(0.0f); // Prevent negative internal energy
                }
                // cell->sete(cell->getE() - 0.5*cell->getDensity()*(v.getMag()*v.getMag()));
                // if (i==rows-1 && j ==cols/2) printf("E:%f e:%f ∆E:%f KE:%f vx:%f vy:%f\n",cell->getE(),cell->gete(), cell->getE()-cell->gete(), 0.5*cell->getDensity()*v.getMag()*v.getMag(), v.getVx(), v.getVy());
                // recalculates temperature and pressure
                cell->getTemp(true);
                cell->getPressure(true);
            }
            if (v.getMag() > maxV) maxV = v.getMag();
            if (!std::isfinite(cell->getDensity())) {
                std::cerr << "NaN detected in density at cell (" << i << ", " << j << ")\n";
                assert(false); // Debug immediately
            }
            if (!std::isfinite(cell->getE())) {
                std::cerr << "NaN detected in energy at cell (" << i << ", " << j << ") velocity: " << v.getMag() << " density: " << cell->getDensity() << "\ne: " << cell->gete() << "\n";
                printf("vx: %f, vy: %f\n", v.getVx(), v.getVy());
                printf("setE=%d noRe=%d\n",setE,noRecalc);
                assert(false); // Debug immediately
            }
            if (!std::isfinite(cell->gete())) {
                std::cerr << "NaN detected in internal energy at cell (" << i << ", " << j << ") velocity: " << v.getMag() << " density: " << cell->getDensity() << "\ne: " << cell->gete() << "\n";
                printf("vx: %f, vy: %f\n", v.getVx(), v.getVy());
                printf("setE=%d noRe=%d\n",setE,noRecalc);
                assert(false); // Debug immediately
            }
            totMass += cell->getMass(false);
        }
    }
    printf("tot mass: %f\n",totMass);
}

void FluidGrid::update(SDL_Event event) {
    int i, j;
    
    if ((event.type == SDL_MOUSEBUTTONDOWN) || (event.type != SDL_MOUSEBUTTONUP)) {

        
        SDL_MouseButtonEvent buttonEvent = event.button;
        if (buttonEvent.button == SDL_BUTTON_RIGHT) {
            mouseVelFlag = !mouseVelFlag;
        } else if (buttonEvent.button == SDL_BUTTON_LEFT) {
            // printf("button left clicked\n");
            buttonHeld += 1;
            Sint32 x = buttonEvent.x;
            Sint32 y = buttonEvent.y;
            if ((prevMouseX == 0) || (prevMouseY == 0)) {
                prevMouseX = x;
                prevMouseY = y;
            }
            i = y * SCALE_H / cellHeight;
            j = x * SCALE_W / cellWidth;
            if (mouseVelFlag == 1) {
                float vxMouse = (x - prevMouseX)*SCALE_W;
                float vyMouse = -(y - prevMouseY)*SCALE_H;
                vGrid->setVx(i, j+1, vxMouse);
                vGrid->setVx(i, j, vxMouse);
                vGrid->setVy(i, j, vyMouse);
                vGrid->setVy(i+1, j, vyMouse);
            } else if (mouseVelFlag == 0) {
                FluidCell *clickedCell = getCell(i, j);
                clickedCell->addMass(100);
                clickedCell->setColor(color[0],color[1],color[2]);
            } else if ((mouseVelFlag == 2) && (buttonHeld < 2) && (sourceList.size() < MAXSOURCES)) {
                // printf("loc: %f %f\n", x*SCALE_W, height-y*SCALE_H);
                // Source *s;
                // int sourceArrayLen = sizeof(*sourceArray)/sizeof(Source);
                sourceArray[sourceList.size()] = Source(x*SCALE_W, height-y*SCALE_H, this, vGrid, SMOKEGUN, color);
                // printf("%f %f\n", sourceArray[0].getX(), sourceArray[0].getY());
                sourceList.push_back(&sourceArray[sourceList.size()]);
                // printf("loc: %f %f\n", sourceList[0]->getX(), sourceList[0]->getY());
                // printf("pushed\n");
                
            } else if ((mouseVelFlag == 3) && (buttonHeld < 2) && (sourceList.size() < MAXSOURCES)) {
                // Source s = Source(x*SCALE_W, height-y*SCALE_H, this, vGrid, POINTSOURCE);
                sourceArray[sourceList.size()] = Source(x*SCALE_W, height-y*SCALE_H, this, vGrid, POINTSOURCE, color);
                // printf("pushing\n");
                sourceList.push_back(&sourceArray[sourceList.size()]);
                // printf("pushed\n");
            }
        } else if (event.key.type == SDL_KEYDOWN) {
            // printf("key clicked\n");
            switch (event.key.keysym.sym) {
                case SDLK_1:
                    mouseVelFlag = 2;
                    // printf("%d\n", mouseVelFlag);
                    break;
                case SDLK_2:
                    mouseVelFlag = 3;
                    // printf("%d\n", mouseVelFlag);
                    break;
                case SDLK_3:
                    mouseVelFlag = 4;
                    break;
                case SDLK_r:
                    color = std::vector<uint> {1,0,0};
                    break;
                case SDLK_g:
                    color = std::vector<uint> {0,1,0};
                    break;
                case SDLK_b:
                    color = std::vector<uint> {0,0,1};
                    break;
                case SDLK_y:
                    color = std::vector<uint> {1,1,0};
                    break;
                case SDLK_p:
                    densityDisplay = 0;
                    temperatureDisplay = 0;
                    pressureDisplay = !(pressureDisplay);
                    printf("%f\n", vGrid->getVy((rows)/2,(cols)/2));
                    break;
                case SDLK_d:
                    pressureDisplay = 0;
                    temperatureDisplay = 0;
                    densityDisplay = !(densityDisplay);
                    break;
                case SDLK_t:
                    densityDisplay = 0;
                    pressureDisplay = 0;
                    temperatureDisplay = !(temperatureDisplay);
                    break;
                case SDLK_0:
                    gravityFlag = !(gravityFlag);
                    printf("gravity\n");
            }


        } 
    } else if (event.type == SDL_MOUSEBUTTONUP) {
        buttonHeld = 0;
        if (mouseVelFlag == 2) {
            SDL_MouseButtonEvent buttonEvent = event.button;
            Sint32 x = buttonEvent.x;
            Sint32 y = buttonEvent.y;
            float vxMouse = (x - prevMouseX)*SCALE_W;
            float vyMouse = -(y - prevMouseY)*SCALE_H;
            Source *s = sourceList.back();

            s->setVx(vxMouse);
            s->setVy(vyMouse);
        }
        prevMouseX = 0;
        prevMouseY = 0;
    }
    for (int i = 0; i < sourceList.size(); i++) {
        Source *s = sourceList[i];
        if ((s->getVx() != 0) || (s->getVy() != 0)) {
            s->setVx(s->getVx());
            s->setVy(s->getVy());
            if ((s->getType() == SMOKEGUN) || (s->getType() == POINTSOURCE)) {
                s->addMass();
            }
        }
    }
    nActive = activeCells.size();
    // activeCells.clear();
    // int i, j;
    // for (i = 0; i < rows; i++) {
    //     for (j = 0; j < cols; j++) {
    //         if (getCell(i, j)->isActive()) {
    //             activeCells.push_back(getCell(i, j));
    //         }
    //     }
    // }
    // for (i = 0; i < nActive; i++) {
    //     FluidCell *cell = activeCells[i];
    //     float initialMass = cell->getMass();
        
    //     cell->transferMass(cell->getBottom(), 5);
    //     totalMass += cell->getMass();
    //     if (!cell->isActive()) activeCells.erase(activeCells.begin() + i);
    //     if (!initialMass && !!cell->getMass()) {
    //         if ((cell->getBottom() != NULL) && (cell->getBottom()->getMass() == 0)) activeCells.push_back(cell->getBottom());
    //         if ((cell->getTop() != NULL) && (cell->getTop()->getMass() == 0)) activeCells.push_back(cell->getTop());
    //         if ((cell->getLeft() != NULL) && (cell->getLeft()->getMass() == 0)) activeCells.push_back(cell->getLeft());
    //         if ((cell->getRight() != NULL) && (cell->getRight()->getMass() == 0)) activeCells.push_back(cell->getRight());
    //         // if ((cell->getBottom() != NULL) && (cell->getBottom()->isActive())) activeCells.push_back(cell->getBottom());
    //         // if ((cell->getTop() != NULL) && (cell->getTop()->isActive())) activeCells.push_back(cell->getTop());
    //         // if ((cell->getLeft() != NULL) && (cell->getLeft()->isActive())) activeCells.push_back(cell->getLeft());
    //         // if ((cell->getRight() != NULL) && (cell->getRight()->isActive())) activeCells.push_back(cell->getRight());
    //     }
        
    // }
    nActive = activeCells.size();
    
    // diffuse(0);
    FluidCell *c = getCell(rows-1,cols/2);
    // printf("dens:%f press:%f vyBottom:%f\n",c->getDensity(),c->getPressure(), vGrid->getVy(rows, cols/2));
    // FluidCell *c = getCell(44,0);
    if (compressible) {
        compressibleMomentumUpdate();
        setVelocities(true, false);
        energyUpdate();
    } else projection(5);
    // printf("projection — dens:%f press:%f\n  vxL:%f vxR:%f vyB:%f vyT:%f\n", c->getDensity(),c->getPressure(),vGrid->getVx(rows-1,0),vGrid->getVx(rows-1,1),vGrid->getVy(rows,0),vGrid->getVy(rows-1,0));
    // printf("projection — dens:%f press:%f\n  vxL:%f vxR:%f vyB:%f vyT:%f\n", c->getDensity(),c->getPressure(),vGrid->getVx(44,0),vGrid->getVx(44,1),vGrid->getVy(45,0),vGrid->getVy(44,0));
    setVelocities(false, true);
    advect();
    setVelocities(false, false);
    // printf("advection — dens:%f press:%f\n  vxL:%f vxR:%f vyB:%f vyT:%f\n", c->getDensity(),c->getPressure(),vGrid->getVx(rows-1,0),vGrid->getVx(rows-1,1),vGrid->getVy(rows,0),vGrid->getVy(rows-1,0));
    // printf("advection — dens:%f press:%f\n  vxL:%f vxR:%f vyB:%f vyT:%f\n", c->getDensity(),c->getPressure(),vGrid->getVx(44,0),vGrid->getVx(44,1),vGrid->getVy(45,0),vGrid->getVy(44,0));
    
    // printf("%f\n", std::min(cellWidth, cellHeight)/maxV);
    // dt = 0.7 * std::min(cellWidth, cellHeight)/maxV;
    dt = std::min(0.016,0.1 * std::min(cellWidth, cellHeight)/maxV);
    viscosity = 0.01*(cellWidth*cellHeight)/dt;
    // for (i = 0; i < rows; i++) {
    //     for (j = 0; j < cols; j++) {
    //         //random mass for now
    //         FluidCell *cell = getCell(i, j);
    //         // cell->transferMass(cell->getRight(), 5);
    //         cell->transferMass(cell->getBottom(), 5);
    //         totalMass += cell->getMass();
    //         // cell->setMass(cell->getMass()-5);
    //     }
    // }
}

std::vector<FluidCell*> FluidGrid::getActive() {
    return activeCells;
}

void FluidGrid::force(int i, int j, float fx, float fy) {
    FluidCell *cell = getCell(i,j);
    float mass = cell->getMass(false);
    if (mass == 0) return; // nothing to force if mass is zero
    float vxL = vGrid->getVx(i,j);
    float vxR = vGrid->getVx(i,j+1);
    float vyU = vGrid->getVy(i,j);
    float vyB = vGrid->getVy(i+1,j);
    float ax = fx / mass;
    float ay = fy / mass;
    if (fx < 0) {
        vGrid->setVx(i,j,vxL+ax*dt);
    } else {
        vGrid->setVx(i,j+1,vxR+ax*dt);
    }
    // if (j == 0) {
    //     vGrid->setVx(i,j,0); 
    // } else
     if (j == cols-1) {
        vGrid->setVx(i,j+1,0);
    }
    
    // vGrid->setVx(i,j+1,vxR+ax*dt);
    if (fy < 0) {
        vGrid->setVy(i+1,j,vyB+ay*dt);
    } else {
        vGrid->setVy(i,j,vyU+ay*dt);
    }
    // if (i == 0) {
    //     vGrid->setVy(i,j,0); 
    // } else 
    if (i == rows-1) {
        vGrid->setVy(i+1,j,0);
    }
    // vGrid->setVy(i+1,j,vyB+ay*dt);
}

class Simulator {
    public:
        Simulator(float width, float height, int c, int r, bool comp): width(width), height(height), rows(r), cols(c), compressible(comp) {
            dt = 0.016;
            // in meters per pixel
            SCALE_H = height / consts::GRID_HEIGHT;
            SCALE_W = width / consts::GRID_WIDTH;
            cells = rows * cols;
            grid = (FluidGrid *) malloc(sizeof(FluidGrid));
            grid[0] = FluidGrid(width, height, rows, cols, dt, compressible);
            
            SDL_Init(SDL_INIT_VIDEO);       // Initializing SDL as Video
            SDL_CreateWindowAndRenderer(consts::GRID_WIDTH, consts::GRID_HEIGHT, 0, &window, &renderer);
            SDL_SetRenderDrawColor(renderer, 0, 0, 0, 0);      // setting draw color
            SDL_RenderClear(renderer);
            
        }

        void drawGridLines() {
            int i, j;
            SDL_SetRenderDrawColor(renderer, 255, 0, 0, 0);
            for (i = 1; i < rows; i++) {
                SDL_RenderDrawLine(renderer, 0, i*consts::GRID_HEIGHT/rows, consts::GRID_WIDTH, i*consts::GRID_HEIGHT/rows);
            }
            for (j = 1; j < cols; j++) {
                SDL_RenderDrawLine(renderer, j*consts::GRID_WIDTH/cols, 0, j*consts::GRID_WIDTH/cols, consts::GRID_HEIGHT);
            }
            // Clear the newly created window
            SDL_RenderPresent(renderer);
        }

        void drawCells() {
            int i, j;
            // save a max and min pressure for the current iteration so 
                // boundaries don't update before the full render
            float thisMaxPressure = 0;
            float thisMinPressure = 0;
            float thisMaxMass = 255;
            float thisMaxDensity = 0;
            float thisMaxTemperature = 0;
            float thisMinTemperature = 0;
            for (i = 0; i < rows; i++) {
                for (j = 0; j < cols; j++) {
                    
            // for (i = 0; i < grid->nActive; i++) {
                    float mass = grid->getCell(i,j)->getMass();
                    if (mass > maxMass) {
                        maxMass = mass;
                    }
                    float pressure = grid->getCell(i,j)->getPressure();
                    float density = grid->getCell(i,j)->getDensity();
                    float temperature = grid->getCell(i,j)->getTemp();

                    if (density > thisMaxDensity) thisMaxDensity = density;
                    // if (thisMaxDensity > maxDensity) maxDensity = thisMaxDensity;
                    if (pressure > thisMaxPressure) {
                        thisMaxPressure = pressure;
                    } else if (pressure < thisMinPressure) {
                        thisMinPressure = pressure;
                    }
                    if (temperature > thisMaxTemperature) thisMaxTemperature = temperature;
                    // if (pressure > maxPressure) {
                    //     maxPressure = pressure;
                    // } else if (pressure < minPressure) {
                    //     minPressure = pressure;
                    // }
                    SDL_Rect rect{j*consts::GRID_WIDTH/cols,i*consts::GRID_HEIGHT/rows,(j+1)*consts::GRID_WIDTH/cols,(i+1)*consts::GRID_HEIGHT/rows};
                    if (grid->pressureDisplay) {
                        // float scaledP_R = 255 * (pressure-minPressure)/(maxPressure-minPressure);
                        float scaledP_B = 255 * (pressure)/(minPressure); // outwards pressure
                        float scaledP_R = 255 * (pressure)/(maxPressure); // inwards pressure
                        // if (scaledP_R > 255) printf("%d %d\n", i,j);
                        if (scaledP_R < 0) scaledP_R = 0;
                        if (scaledP_B < 0) scaledP_B = 0;
                        SDL_SetRenderDrawColor(renderer, scaledP_R, 0, scaledP_B, 255);
                        // SDL_SetRenderDrawColor(renderer, mass, mass, 0, 255);
                        // SDL_RenderFillRect(renderer, &rect);
                    } else if (grid->densityDisplay) {
                        // FluidCell *cell = grid->getCell(i, j);
                        
                        float scaled_dens = density/maxDensity * 255;
                        if (scaled_dens > 255) scaled_dens = 255;
                        SDL_SetRenderDrawColor(renderer, scaled_dens, scaled_dens, 0, 255);
                        // SDL_RenderFillRect(renderer, &rect);
                    } else if (grid->temperatureDisplay) {
                        float scaled_temp = temperature/maxTemperature * 255;
                        if (grid->getCell(i,j)->getTemp() < 0) {
                            // std::cerr << "Temperature (" << i << ", " << j << ")\n";
                            // assert(false); // Debug immediately
                        }
                        SDL_SetRenderDrawColor(renderer, 0, scaled_temp, 0, 255);
                    } else {
                        // float mass = grid->getActive()[i]->getMass();
                        // FluidCell *cell = grid->getActive()[i];
                        std::vector<float> color = grid->getCell(i,j)->getColor();
                        
                        if ((i == rows/2) && (j == cols/2)) {
                            // std::cout << color[0] << std::endl;
                            // printf("%f %f %f\n", pressure, mass, grid->getCell(i,j)->getMass(false));
                        }
                        // SDL_Rect rect{cell->getCol()*consts::GRID_WIDTH/cols,cell->getRow()*consts::GRID_HEIGHT/rows,(cell->getCol()+1)*consts::GRID_WIDTH/cols,(cell->getRow()+1)*consts::GRID_HEIGHT/rows};
                        // if (mass > 255) {
                        //     mass = 255;
                        // }
                        // if (!cell->isActive()) SDL_SetRenderDrawColor(renderer, 0, 255, 0, 255);
                        uint red = std::min((float)255, 255*color[0]*std::sqrt(mass/maxMass));
                        uint green = std::min((float)255, 255*color[1]*std::sqrt(mass/maxMass));
                        uint blue = std::min((float)255, 255*color[2]*std::sqrt(mass/maxMass));
                        SDL_SetRenderDrawColor(renderer, red, green, blue, 255);
                        // SDL_SetRenderDrawColor(renderer, mass, mass, 0, 255);
                        // SDL_RenderFillRect(renderer, &rect);
                    }
                    SDL_RenderFillRect(renderer, &rect);
                }
            }
            SDL_RenderPresent(renderer);
            maxPressure = thisMaxPressure;
            minPressure = thisMinPressure;
            maxDensity = thisMaxDensity;
            maxTemperature = thisMaxTemperature;
            minTemperature = thisMinTemperature;
            // printf("maxtemp:%f\n",maxTemperature);
            // thisMaxDensity = 1;
            // thisMinPressure = -1;
            // thisMaxPressure = 1;
            
        }
        
        void step(SDL_Event event) {
            grid->update(event);
            drawCells();
            SDL_Delay(grid->getdt()*1000);  // setting some Delay
        }

        void freeSim() {
            grid->freeGrid();
            
            free(grid);
        }
        float SCALE_H, SCALE_W;
        float width, height;
        FluidGrid *grid;
    private:
        bool compressible;
        float dt;
        int rows, cols, cells;
        SDL_Renderer *renderer = NULL;
        SDL_Window *window = NULL;
        SDL_Surface *screenSurface;
        float maxMass = 255;
        float maxPressure = 1;
        float minPressure = 0;
        float maxDensity = 1;
        float maxTemperature = 1;
        float minTemperature = 0;
};


int main(int argv, char **argc) {
    if (argv > 6) std::srand((unsigned) atoi(argc[6]));
    else std::srand((unsigned) std::time(NULL));
    float width, height;
    int c, r;
    bool comp;
    if (argv < 6) {
        width = atoi(argc[1]);
        height = atoi(argc[2]);
        c = atoi(argc[3]);
        r = atoi(argc[4]);
        comp = false;
        // sim.drawCells();
    } else {
        width = atoi(argc[2]);
        height = atoi(argc[3]);
        c = atoi(argc[4]);
        r = atoi(argc[5]);
        if (atoi(argc[1])) comp = true; else comp = false;
        // comp = atoi(argc[1]);
        // Simulator sim(atoi(argc[2]), atoi(argc[3]), atoi(argc[4]), atoi(argc[5]), (bool)atoi(argc[1]));
        
    }

    // Simulator sim(atoi(argc[1]), atoi(argc[2]), atoi(argc[3]), atoi(argc[4]), false);
    Simulator sim(width, height, c, r, comp);
    sim.drawCells();
    SDL_Event event;
    // 
    // sim.step(event);
    // sim.drawGridLines();
    uint pause = 0;
    while(!(event.type == SDL_QUIT)){
        if (!pause) sim.step(event);
        SDL_PollEvent(&event);  // Catching the poll event.
        if (event.key.type == SDL_KEYDOWN) {
            // printf("key clicked\n");
            if (event.key.keysym.sym == SDLK_SPACE) pause = !pause;
        }
        if (event.type == SDL_MOUSEBUTTONDOWN) {
            SDL_MouseButtonEvent buttonEvent = event.button;
            Sint32 x = buttonEvent.x;
            Sint32 y = buttonEvent.y;
            float physY = sim.height - y * sim.SCALE_H;
            float physX = x * sim.SCALE_W;
            VelocityVector v = sim.grid->vGrid->sampleVelocityAtPoint(physX,physY);
            // FluidCell *clickedCell = sim.grid->getCell(i, j);
            // printf("%f %f\n", v.getVx(), v.getVy());
        }
    }
    sim.freeSim();
    return 1;
}

