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
            return pow(vx*vx+vy*vy,1/2);
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
        FluidCell(float mass, float width, float height, float temp, float red, float green, float blue): mass(mass), width(width), height(height), temperature(temp) {
            // std::vector<uint> color
            // size is the physical area
            size = width * height;
            density = mass/size;
            pressure = 0;
            // pressure = density/Hmol * consts::kb * temperature;
            // r = color[0];
            // g = color[1];
            // b = color[2];
            // r = red;
            // g = green;
            // b = blue;
            // setColor((!!red) * mass, (!!green) * mass, (!!blue) * mass);
            setColor(red, green, blue);
            // velocity = VelocityVector();
            // velocity.setVx(1);
            // velocity.setVx((std::rand() % (int)(2*width)) - width/2);
            // velocity.setVy((std::rand() % (int)(2*height)) - height/2);
            // vBounds = VelocityBox(0,0,0,0);
            // vBounds.setVelocity(velocity);
            neighbors = Neighbors();
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
        float getTemp() {
            return temperature;
        }
        void setTemp(float t) {
            temperature = t;
        }
        float getPressure() {
            return pressure;
        }
        void setPressure(float p) {
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
            // velocity = vel2;
            // vBounds.setVelocity(vel2);
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
        float mass, size, density, temperature, pressure, width, height, smokeMass;
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
            if (type == SMOKEGUN) {
                if (v > 0) {
                    grid->force(i,j+1,vx*100,0);
                } else {
                    grid->force(i,j,vx*100,0);
                }
                // vGrid->setVx(i,j,vx);
            } else if (type == POINTSOURCE) {
                vGrid->setVx(i,j,-vx);
            }
            // vGrid->setVx(i,j+1,vx);
        }
        void setVy(float v) {
            vy = v;
            // vGrid->setVy(i,j,vy);
            if (type == SMOKEGUN) {
                if (v > 0) {
                    grid->force(i,j,0,vy*100);
                } else {
                    grid->force(i+1,j,0,vy*100);
                }
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

FluidGrid::FluidGrid(float width, float height, int r, int c, float dt): width(width), height(height), rows(r), cols(c), dt(dt) {
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
    for (i = 0; i < rows; i++) {
        grid[i] = (FluidCell *) malloc(sizeof(FluidCell)*c);
        newGrid[i] = (FluidCell *) malloc(sizeof(FluidCell)*c);
        for (j = 0; j < cols; j++) {
            //random mass for now
            // float mass = std::rand() % 256;
            std::vector<uint> color = {0, 0, 0};
            grid[i][j] = FluidCell(0, cellWidth, cellHeight, 0, 0,0,0);
            if (i < rows/2) {
                grid[i][j] = FluidCell(100, cellWidth, cellHeight, 300, 0,0,0);
                // std::vector<uint> color = {100, 100, 0};
                // grid[i][j] = FluidCell(100, cellWidth, cellHeight, 100, 1,1,0);
                // vGrid->setVy(i,j, 50);
                // vGrid->setVx(i,j, -5);
                // vGrid->setVy(i-1,j, 5);
                // grid[i][j].setVelocity(VelocityVector(10,0));
            }
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
            if (grid[i][j].isActive()) {
                activeCells.push_back(&grid[i][j]);
            }
            // if (grid[i][j].getVelocity().getMag() > maxV) {
            //     maxV = grid[i][j].getVelocity().getMag();
            // }
            grid[i][j].setLoc(i, j);
            newGrid[i][j] = grid[i][j];
            new_vGrid->setVx(i,j,vGrid->getVx(i,j));
            new_vGrid->setVy(i,j,vGrid->getVy(i,j));
        }
    }
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
    free(grid);
}

std::map<uint, float> FluidGrid::sampleCellAtPoint(float x, float y) {
    // y = round(y);
    // float cHeight = round(cellHeight);
    // printf("%f\n", round(height - y - cellHeight/2));
    // need switch statements for sampled property "prop"
    int i = ceil((height - y - cellHeight/2) / cellHeight);
    int j = floor((x - cellWidth/2) / cellWidth);
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

void FluidGrid::projection(int iters) {
    int i, j;
    for (int n = 0; n < iters; n++) {
        for (i = 0; i < rows; i++) {
            for (j = 0; j < cols; j++) {
                FluidCell *cell = getCell(i,j);
                
                // printf("%d %d\n", i, j);
                int sLeft = !!j;
                int sRight = !!(cols - 1 - j);
                int sTop = !!i;
                int sBottom = !!(rows - 1 - i);                        
                // printf("%d %d %d %d\n", sLeft, sRight, sTop, sBottom);
                int s = sLeft + sRight + sTop + sBottom;
                // printf("%d\n", s);
                // if ((i > 0) && (i < rows -1) && (j > 0) && (j < cols - 1)) {
                float left = vGrid->getVx(i, j);
                float right = vGrid->getVx(i, j+1);
                float top = vGrid->getVy(i, j);
                float bottom = vGrid->getVy(i+1, j);
                float mass = cell->getMass(false);
                
                // if (n == 0) {
                    
                // }
                float d = (-left*sLeft + right*sRight + top*sTop - bottom*sBottom)*1.5;
                // if ((i==rows/2) && (j==cols/2) && (n == 0)) printf("%f\n", d);
                if (n == 0) {
                    float pressure = cell->getDensity()/consts::Hmol * consts::kb * consts::Na * cell->getTemp();
                    
                    
                    cell->setPressure(pressure);
                    // cell->setPressure(0);
                    force(i,j,0,mass*(-100));
                    if (i==rows/2 && j == cols/2) {
                        // printf("mid:%f %f, vy: %f\n", cell->getDensity(), pressure, vGrid->getVy(i,j));
                    }
                    else if (i==rows-1 && j == cols/2) {
                        // printf("bottom:%f %f, vy: %f\n", cell->getDensity(), pressure, vGrid->getVy(i,j));
                    }
                } else {
                    float newPressure = cell->getPressure() + d/s * cell->getDensity()*std::sqrt(cellWidth*cellHeight)/dt;
                    
                    cell->setPressure(newPressure);
                }
                
                // if ((i == rows/2+1) && (j == cols/2)) {
                //     printf("%f %f %f %f\n", left, right, top, bottom);
                //     printf("%d %d %d %d\n", sLeft, sRight, sTop, sBottom);
                //     printf("%f\n", d);
                // }
                vGrid->setVx(i, j, left+d*sLeft/s);
                vGrid->setVx(i, j+1, right-d*sRight/s);
                vGrid->setVy(i, j, top-d*sTop/s);
                vGrid->setVy(i+1, j, bottom+d*sBottom/s);
                if (n == iters-1) {
                    // if (cell->getPressure() > 0) {

                    // }
                    // float T1 = cell->getTemp();
                    // cell->setTemp(cell->getPressure()/(consts::kb*consts::Na*cell->getDensity()/consts::Hmol));
                    // float T2 = cell->getTemp();
                    // float m2 = mass*T2/T1;
                    // cell->setMass(m2, false);
                    
                }
                
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
                float physX = cellWidth*j+cellWidth/2;
                float physY = height - (cellHeight*i+cellHeight/2); // I want y pointed up
                VelocityVector vCell = vGrid->sampleVelocityAtPoint(physX, physY);
                float prevX = physX - dt*vCell.getVx();
                float prevY = physY - dt*vCell.getVy();
                std::map<uint, float> prevProps;
                if (round(vCell.getVx()+vCell.getVy()) == 0) {
                    prevProps[SMOKEMASS] = cell.getMass();
                    prevProps[MASS] = cell.getMass(false);
                    prevProps[TEMPERATURE] = cell.getTemp();
                    prevProps[PRESSURE] = cell.getPressure();
                    prevProps[R] = cell.getColor()[0];
                    prevProps[G] = cell.getColor()[1];
                    prevProps[B] = cell.getColor()[2];
                } else {
                    prevProps = sampleCellAtPoint(prevX, prevY);
                }

                if ((prevX < 0) || (prevY < 0) || (prevX > width) || (prevY > height)) {
                    newGrid[i][j].setMass(0);
                    newGrid[i][j].setMass(0,false);
                    newGrid[i][j].setColor(0,0,0);
                } else {
                    newGrid[i][j].setMass(prevProps[SMOKEMASS]);
                    newGrid[i][j].setMass(prevProps[MASS],false);
                    // if ((i == rows/2) && (j == cols/2)) printf("%f %f %f\n", prevProps[R], prevProps[G], prevProps[B]);
                    newGrid[i][j].setColor(prevProps[R], prevProps[G], prevProps[B]);
                }
            }
            if (i < rows) {
                // vxs of vgrid
                float vxPhysX = cellWidth*j;
                float vxPhysY = height - cellHeight*i - cellHeight/2;
                VelocityVector vVx = vGrid->sampleVelocityAtPoint(vxPhysX, vxPhysY);
                float prevVxX = vxPhysX - dt*vVx.getVx();
                float prevVxY = vxPhysY - dt*vVx.getVy();
                float newVx = vGrid->sampleVelocityAtPoint(prevVxX, prevVxY).getVx();
                new_vGrid->setVx(i,j,newVx);

                // find the max velocity in x direction for dt
                if (newVx > maxV) {
                    maxV = newVx;
                }
            }
            if (j < cols) {
                // vys of vgrid
                float vyPhysX = cellWidth*j + cellWidth/2;
                float vyPhysY = height - cellHeight*i;
                VelocityVector vVy = vGrid->sampleVelocityAtPoint(vyPhysX, vyPhysY);
                float prevVyX = vyPhysX - dt*vVy.getVx();
                float prevVyY = vyPhysY - dt*vVy.getVy();
                float newVy = vGrid->sampleVelocityAtPoint(prevVyX, prevVyY).getVy();
                new_vGrid->setVy(i,j,newVy);
                
                // find the max velocity in x direction for dt
                if (newVy > maxV) {
                    maxV = newVy;
                }
            }
        }
    }
    grid = newGrid;
    vGrid = new_vGrid;
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
                    pressureDisplay = !(pressureDisplay);
                    printf("%f\n", vGrid->getVy((rows)/2,(cols)/2));
                    break;
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
    
    projection(2);
    advect();
    // printf("%f\n", std::min(cellWidth, cellHeight)/maxV);
    dt = 0.7 * std::min(cellWidth, cellHeight)/maxV;
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
    // vGrid->setVx(i,j+1,vxR+ax*dt);
    if (fy < 0) {
        vGrid->setVy(i+1,j,vyB+ay*dt);
    } else {
        vGrid->setVy(i,j,vyU+ay*dt);
    }
    // vGrid->setVy(i+1,j,vyB+ay*dt);
}

class Simulator {
    public:
        Simulator(float width, float height, int c, int r): width(width), height(height), rows(r), cols(c) {
            dt = 0.016;
            // in meters per pixel
            SCALE_H = height / consts::GRID_HEIGHT;
            SCALE_W = width / consts::GRID_WIDTH;
            cells = rows * cols;
            grid = (FluidGrid *) malloc(sizeof(FluidGrid));
            grid[0] = FluidGrid(width, height, rows, cols, dt);
            
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
            for (i = 0; i < rows; i++) {
                for (j = 0; j < cols; j++) {
                    
            // for (i = 0; i < grid->nActive; i++) {
                    float mass = grid->getCell(i,j)->getMass();
                    if (mass > maxMass) {
                        maxMass = mass;
                    }
                    float pressure = grid->getCell(i,j)->getPressure();
                    if (pressure > maxPressure) {
                        // maxPressure = pressure;
                    } else if (pressure < minPressure) {
                        minPressure = pressure;
                    }
                    SDL_Rect rect{j*consts::GRID_WIDTH/cols,i*consts::GRID_HEIGHT/rows,(j+1)*consts::GRID_WIDTH/cols,(i+1)*consts::GRID_HEIGHT/rows};
                    if (grid->pressureDisplay) {
                        // float scaledP_R = 255 * (pressure-minPressure)/(maxPressure-minPressure);
                        float scaledP_R = -255 * (pressure)/(maxPressure-minPressure); // inwards pressure
                        float scaledP_B = 255 * (pressure)/(maxPressure-minPressure); // outwards pressure

                        if (scaledP_R < 0) scaledP_R = 0;
                        if (scaledP_B < 0) scaledP_B = 0;
                        SDL_SetRenderDrawColor(renderer, scaledP_R, 0, scaledP_B, 255);
                        // SDL_SetRenderDrawColor(renderer, mass, mass, 0, 255);
                        SDL_RenderFillRect(renderer, &rect);
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
                        SDL_RenderFillRect(renderer, &rect);
                    }
            // }
                }
            }
            SDL_RenderPresent(renderer);
            minPressure = 0;
            maxPressure = 1;
        }
        
        void step(SDL_Event event) {
            grid->update(event);
            drawCells();
            SDL_Delay(grid->getdt()*1);  // setting some Delay
        }

        void freeSim() {
            grid->freeGrid();
            free(grid);
        }
        float SCALE_H, SCALE_W;
        float width, height;
        FluidGrid *grid;
    private:
        float dt;
        int rows, cols, cells;
        SDL_Renderer *renderer = NULL;
        SDL_Window *window = NULL;
        SDL_Surface *screenSurface;
        float maxMass = 255;
        float maxPressure = 1;
        float minPressure = 0;
        
};


int main(int argv, char **argc) {
    if (argv > 5) std::srand((unsigned) atoi(argc[5]));
    else std::srand((unsigned) std::time(NULL));
    Simulator sim(atoi(argc[1]), atoi(argc[2]), atoi(argc[3]), atoi(argc[4]));
    sim.drawCells();
    SDL_Event event;
    // 
    // sim.step(event);
    // sim.drawGridLines();

    while(!(event.type == SDL_QUIT)){
        
        sim.step(event);
        SDL_PollEvent(&event);  // Catching the poll event.
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

