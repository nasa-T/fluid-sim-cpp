#include <cmath>
#include <iostream>
#include <fstream>
#include "fluid.h"

#include <map>
#include <string>

#include <omp.h>
#include <algorithm>

// #include <GL/glew.h>

// #include <GLFW/glfw3.h>
// GLFWwindow* window;

// #include <glm/glm.hpp>
// #include <glm/gtc/matrix_transform.hpp>
// #include <glm/gtx/norm.hpp>
// using namespace glm;

#include <chrono>

void print(char *string) {
    std::cout << string << std::endl;
}

class VelocityVector {
    public:
        VelocityVector(double vx_ = 0, double vy_ = 0, double ax_ = 0, double ay_ = 0): vx(vx_), vy(vy_), ax(ax_), ay(ay_) {}
        VelocityVector operator+(VelocityVector& other) {
            VelocityVector added(this->vx + other.vx, this->vy + other.vy, other.ax, other.ay);
            return added;
        }
        double getVx() {
            return vx;
        }
        double getVy() {
            return vy;
        }
        double getAx() {
            return ax;
        }
        double getAy() {
            return ay;
        }
        void setVx(double v) {
            vx = v;
        }
        void setVy(double v) {
            vy = v;
        }
        void accelerate() {
            vx += ax * consts::dt;
            vy += ay * consts::dt;
        }
        double getMag() {
            return std::sqrt(getVx()*getVx()+getVy()*getVy());
            // return pow(vx*vx+vy*vy,1/2);
        }
        // std::ostream& operator<<(std::ostream &s, VelocityVector &vec) {
        //     return s << "vx: " << vec.getVx() << " vy: " << vec.getVy();
        // }
    private:
        double vx, vy, ax, ay;
};

class VelocityBox {
    public:
        double top, bottom, left, right;
        VelocityBox(double top = 0, double bottom = 0, double left = 0, double right = 0): top(top), bottom(bottom), left(left), right(right) {}
        void setVelocity(VelocityVector v) {
            left = v.getVx();
            right = v.getVx();
            top = v.getVy();
            bottom = v.getVy();
        }
};

class VelocityGrid {
    public:
        VelocityGrid(double width, double height, int rows, int cols): width(width), height(height), rows(rows), cols(cols) {
            int i, j;
            vyArray = (double **) malloc(sizeof(double *) * (rows + 1));
            vxArray = (double **) malloc(sizeof(double *) * rows);
            for (i = 0; i < rows + 1; i++) {
                vyArray[i] = (double *) malloc(sizeof(double) * cols);
                for (j = 0; j < cols; j++) {
                    vyArray[i][j] = 0;
                }
            }
            for (i = 0; i < rows; i++) {
                vxArray[i] = (double *) malloc(sizeof(double) * (cols + 1));
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
        VelocityVector sampleVelocityAtPoint(double x, double y) {
            int vyi = ceil((height - y) / cellHeight);
            int vyj = floor((x-cellWidth/2) / cellWidth);
            // the physical positions relative to nearest box of vys
            double cvyX, cvyY; 
            double cvxX, cvxY; // same for vxs
            cvyX = (fmod(x+cellWidth/2, cellWidth))/cellWidth;
            cvyY = fmod(y, cellHeight)/cellHeight;
            // will throw an error for when x and y are off the grid entirely
            double vy10, vy00, vy01, vy11 = 0;
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

            double vy = (1-cvyX)*(1-cvyY)*vy10 + cvyX*(1-cvyY)*vy11 + (1-cvyX)*cvyY*vy00 + cvyX*cvyY*vy01;

            int vxi = ceil((height - (y+cellHeight/2)) / cellHeight);
            int vxj = floor(x / cellWidth);

            cvxX = fmod(x, cellWidth)/cellWidth;
            cvxY = (fmod(y+cellHeight/2, cellHeight))/cellHeight;

            double vx10, vx00, vx01, vx11 = 0;
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
            double vx = (1-cvxX)*(1-cvxY)*vx10 + cvxX*(1-cvxY)*vx11 + (1-cvxX)*cvxY*vx00 + cvxX*cvxY*vx01;
            return VelocityVector(vx, vy);
        }
        VelocityVector getVelocityVector(int i, int j) {
            return VelocityVector(vxArray[i][j], vyArray[i][j]);
        }
        double getVy(int i, int j) {
            return vyArray[i][j];
        }
        double getVx(int i, int j) {
            return vxArray[i][j];
        }
        void setVy(int i, int j, double vy) {
            vyArray[i][j] = vy;
        }
        void setVx(int i, int j, double vx) {
            vxArray[i][j] = vx;
        }
        double getWidth() {
            return width;
        }
        double getHeight() {
            return height;
        }
        double getCellWidth() {
            return cellWidth;
        }
        double getCellHeight() {
            return cellHeight;
        }
        std::map<char, double> getxyFromijVx(int i, int j) {
            std::map<char, double> xy;
            xy['x'] = cellWidth*j;
            xy['y'] = height - (i*cellHeight + cellHeight/2);
            return xy;
        }
        std::map<char, double> getxyFromijVy(int i, int j) {
            std::map<char, double> xy;
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
        double width, height;
        double cellWidth, cellHeight;
        int rows, cols;
        double **vyArray, **vxArray;
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
        FluidCell(double mass, double width, double height, double temp, double red, double green, double blue, bool comp): mass(mass), width(width), height(height), temperature(temp), compressible(comp) {
            // std::vector<uint> color
            // size is the physical area
            size = width * height;
            density = mass/size;
            // pressure = 0;
            // pressure = density/consts::HMol * consts::r * temperature;
            
            // helium = 0.0;
            // metals = 0.0;
            // hydrogen = 1.0-helium-metals;
            setX(1.0);
            setY(0.0);
            setZ(0.0);
            double nonMet = (hydrogen+helium);
            double cv = consts::r/(consts::HMol*hydrogen/nonMet+consts::HeMol*helium/nonMet)/(1.4-1);
            // double cv = consts::r/(consts::HMol*hydrogen)/(1.4-1);
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
            double gravPotential = 0;
            
        }
        double getCv() {
            double nonMet = (hydrogen+helium);
            return consts::r/(consts::HMol*hydrogen/nonMet+consts::HeMol*helium/nonMet)/(1.4-1);
            // return consts::r/(consts::HMol*hydrogen)/(1.4-1);
        }
        double getMass(bool smoke=true) {
            if (smoke) return smokeMass;
            else return mass;
        }
        double getSize() {
            return size;
        }
        double getDensity() {
            return mass/size;
        }
        double getTemp(bool recalculate=true) {
            if (recalculate && compressible) {
                // double cv = consts::r/consts::HMol/(1.4-1);
                if (round(mass) == 0) {
                    temperature = 0;
                } else {
                    temperature = gete(false)/(getDensity()*getCv());
                }
            }
            return temperature;
        }
        void setTemp(double t) {
            temperature = t;
        }
        double getPressure(bool recalculate=true) {
            // double oldPressure = pressure;
            if (compressible && recalculate) {
                // pressure = getDensity()/consts::HMol * consts::r * getTemp();
                pressure = (1.4 - 1)*gete(false);
            }
            return pressure;
        }
        void setPressure(double p) {
            if (!!pressure) {
                // double newTemperature = temperature * p/pressure;
                // printf("newtemp: %f\n",newTemperature);
                // temperature = newTemperature;
            }
            pressure = p;
        }
        void setMass(double m, bool smoke=true) {
            if (m > 0) {
                if (smoke) smokeMass = m;
                else mass = m;
            } else {
                mass = 0;
            }
            // setColor(mass*r, mass*g, mass*b);
        }
        void addMass(double m, bool smoke=true) {
            setMass(mass + m, smoke);
        }
        std::vector<double> getColor() {
            return std::vector<double> {r, g, b};
        }

        double getX() {
            return hydrogen;
        }

        void setX(double h) {
            if (round(mass) > 0) {
                hydrogen = h;
            } else {
                hydrogen = 1.0;
            }

        }

        double getY() {
            return helium;
        }

        void setY(double he) {
            if (round(mass) > 0) {
                helium = he;
            } else {
                helium = 0;
            }
        }

        double getZ() {
            return metals;
        }

        void setZ(double Z) {
            if (round(mass) > 0) {
                metals = Z;
            } else {
                metals = 0;
            }
        }

        void HtoHe(double dX) {
            if (getX() - dX < 0) {
                setX(0);
                setY(1);
            } else {
                setX(getX() - dX);
                setY(1 - getX());
            }
        }

        double gete(bool recalculate=true) {
            if (recalculate) {
                // if (round(getDensity()*1000) == 0) e = 0;
                // else 
                e = getCv()*getTemp()*getDensity();
                // e = E - 1/2*getDensity()*velocity.getMag()*velocity.getMag();
            }
            return e;
        }

        void sete(double energy) {
            if (mass > 0) {
                e = energy;
            } else {
                e = 0;
            }
        }
        
        double getE(bool recalculate=false) {
            // if (recalculate) {
            //     E = e + 1/2*getDensity()*velocity.getMag()*velocity.getMag();
            // }
            return E;
        }

        void setE(double energy) {
            if (mass > 0) {
                E = energy;
            } else {
                E = 0;
            }
        }

        void setColor(double red, double green, double blue, bool source=false) {
            double tot = red + green + blue;
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

        double getGravPotential() {
            return gravPotential;
        }

        void setGravPotential(double u) {
            gravPotential = u;
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

        void transferMass(FluidCell *other, double m) {
            if (other != NULL) { // halt at the boundaries; NOTHING GETS OUT
                if (m > 0) {
                    if (this->getMass() < m) {
                        m = this->getMass();
                    }
                    double otherMass = other->getMass();
                    other->setMass(otherMass + m);
                    this->setMass(mass - m);
                } else if (m < 0) {
                    m = -m;
                    if (other->getMass() < m) {
                        m = other->getMass();
                    }
                    double otherMass = other->getMass();
                    other->setMass(otherMass - m);
                    this->setMass(mass + m);
                }
            }
        }
        // void addVelocity(VelocityVector vel2) {
        //     vel = vel + vel2;
        // }
        void setVelocity(double vx, double vy) {
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
        double mass, size, density, temperature, pressure, width, height, smokeMass, e, E, gravPotential, hydrogen, helium, metals;
        bool compressible;
        int row, col;
        VelocityVector velocity;
        // VelocityBox vBounds;
        Neighbors neighbors;
        double r,g,b;
        // std::vector<uint> color;
};
// class FluidGrid {
//     public:
//         FluidGrid(double width, double height, int r, int c): width(width), height(height), rows(r), cols(c) {}
//         FluidCell *getCell(int i, int j) {
//             return &grid[i][j];
//         }
//         private:
//             double width, height;
//             double cellWidth, cellHeight;
//             double dt, maxV;
//             int rows, cols, cells;
//             FluidCell **grid;
//             FluidCell **newGrid;
//             std::vector<FluidCell*> activeCells;
//             double SCALE_H, SCALE_W;
            
//             VelocityGrid *new_vGrid;
//             std::vector<Source*> sourceList;
// };
class Source {
    public:
        Source(double _x, double _y, FluidGrid *gd, VelocityGrid *vgrid, uint sourceType, std::vector<uint> color):  grid(gd), vGrid(vgrid), type(sourceType) {
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
        void setX(double x) {
            xs = x;
        }
        void setY(double y) {
            ys = y;
        }
        double getX() {
            return xs;
        }
        double getY() {
            return ys;
        }
        double getVx() {
            return vx;
        }
        double getVy() {
            return vy;
        }
        void setVx(double v) {
            // vx = ((v > 0) - (v < 0)) * 1000;
            // printf("%d %d: %f\n", i,j,vx);
            vx = v;
            double mass = grid->getCell(i,j)->getMass(false);
            if (type == SMOKEGUN) {
                // if (v > 0) {
                //     grid->force(i,j+1,vx*100,0);
                // } else {
                //     grid->force(i,j,vx*100,0);
                // }
                grid->force(i, j + 1, vx * 25 * mass, 0);
                grid->force(i, j - 1, vx * 25 * mass, 0);
                grid->force(i,j,vx*50*mass,0);
                // vGrid->setVx(i,j,vx);
            } else if (type == POINTSOURCE) {
                vGrid->setVx(i,j,-vx);
            }
            // vGrid->setVx(i,j+1,vx);
        }
        void setVy(double v) {
            vy = v;
            // vy = 0;
            // vGrid->setVy(i,j,vy);
            double mass = grid->getCell(i,j)->getMass(false);
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
            double v = std::sqrt(vx*vx+vy*vy);
            double cellSize = vGrid->getCellHeight() * vGrid->getCellWidth();
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
        double vx, vy = 0;
        double xs, ys;
        int i, j;
};

FluidGrid::FluidGrid(double width, double height, int r, int c, double dt, bool comp): width(width), height(height), rows(r), cols(c), dt(dt), compressible(comp) {
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
            // double mass = (std::rand() % 100) * 1e16;
            double mass = 2e17;
            std::vector<uint> color = {0, 0, 0};
            grid[i][j] = FluidCell(mass, cellWidth, cellHeight, 50, 0,0,0,comp);
            if ((i == rows/2 && (j == cols/2 || j == cols/2-1)) || (i == rows/2-1 && (j == cols/2 || j == cols/2-1))) {
                // grid[i][j] = FluidCell(1e21, cellWidth, cellHeight, 50, 0,0,0,comp);
                grid[i][j] = FluidCell(5e17, cellWidth, cellHeight, 50, 0,0,0,comp);
            }
            if (i == 25 && j == 25) {
                // grid[i][j] = FluidCell(1e21, cellWidth, cellHeight, 50, 0,0,0,comp);
                // grid[i][j] = FluidCell(1e18, cellWidth, cellHeight, 50, 0,0,0,comp);
            }
            if (i < rows/5 || i > 4*rows/5 || j < cols/5 || j > 4*cols/5) {
                // grid[i][j] = FluidCell(1e22, cellWidth, cellHeight, 10, 0,0,0,comp);
            }
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
            double randVx = (std::rand() % (int)cellWidth/5)-cellWidth/10;
            double randVy = (std::rand() % (int)cellHeight/5)-cellHeight/10;
            // vGrid->setVx(i,j,randVx);
            // vGrid->setVy(i,j,randVy);
            // new_vGrid->setVx(i,j,randVx);
            // new_vGrid->setVy(i,j,randVy);
            vGrid->setVx(i,j,0);
            vGrid->setVy(i,j,0);
            new_vGrid->setVx(i,j,0);
            new_vGrid->setVy(i,j,0);
        }
    }
    if (compressible) setVelocities(true,false);
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

std::map<uint, double> FluidGrid::sampleCellAtPoint(double x, double y) {
    // y = round(y);
    // double cHeight = round(cellHeight);
    // printf("%f\n", round(height - y - cellHeight/2));
    int i,j;
    if (round(1000*(height - y - cellHeight/2)) == 0) {
        i = 0;
    } else {
        i = ceil((height - y - cellHeight/2) / cellHeight);
    }
    if (round(1000*(x - cellWidth/2)) == 0) {
        j = 0;
    } else {
        j = floor((x - cellWidth/2) / cellWidth);
    }
    
    if ((i == rows/2) && (j == cols/2)) {
        // printf("%d: %f %d: %f\n", j,x, i,y);
        // printf("%f %f\n", prevY, physY);
    }
    // switch (prop) {
    //     case MASS:

    // }

    // the physical positions relative to nearest box of cells
    double X, Y; 
    X = fmod(x+cellWidth/2, cellWidth)/cellWidth;
    Y = fmod(y+cellHeight/2, cellHeight)/cellHeight;

    // will throw an error for when x and y are off the grid entirely
    std::map<uint, double> p10, p00, p01, p11;
    if ((j > -1) && (i < rows)) {
        // printf("%d %d\n", j, i);
        // p10 = vyArray[i][j];
        FluidCell *cell = getCell(i,j);
        p10[SMOKEMASS] = cell->getMass();
        p10[MASS] = cell->getMass(false);
        p10[TEMPERATURE] = cell->getTemp();
        p10[PRESSURE] = cell->getPressure();
        p10[E] = cell->getE(false);
        p10[GRAVITY] = cell->getGravPotential();
        // p10[R] = cell->getColor()[0]*p10[MASS];
        // p10[G] = cell->getColor()[1]*p10[MASS];
        // p10[B] = cell->getColor()[2]*p10[MASS];
        p10[R] = cell->getColor()[0];
        p10[G] = cell->getColor()[1];
        p10[B] = cell->getColor()[2];
        p10[HYDROGEN] = cell->getX();
        p10[HELIUM] = cell->getY();
    } else {
        if (compressible && 0) {
            FluidCell *cell = getCell(i-1,j+1);
            p10[SMOKEMASS] = cell->getMass();
            p10[MASS] = cell->getMass(false);
            p10[TEMPERATURE] = cell->getTemp();
            p10[PRESSURE] = cell->getPressure();
            p10[E] = cell->getE(false);
            p10[GRAVITY] = cell->getGravPotential();
            p10[R] = cell->getColor()[0];
            p10[G] = cell->getColor()[1];
            p10[B] = cell->getColor()[2];
            p10[HYDROGEN] = cell->getX();
            p10[HELIUM] = cell->getY();
        } else {
            p10[SMOKEMASS] = 0;
            p10[MASS] = 0;
            p10[TEMPERATURE] = 0;
            p10[PRESSURE] = 0;
            p10[E] = 0;
            p10[R] = 0;
            p10[G] = 0;
            p10[B] = 0;
            p10[HYDROGEN] = 0;
            p10[HELIUM] = 0;
        }
    }
    if ((j > -1) && (i > 0)) {
        // p00 = vyArray[i-1][j];
        FluidCell *cell = getCell(i-1,j);
        p00[SMOKEMASS] = cell->getMass();
        p00[MASS] = cell->getMass(false);
        p00[TEMPERATURE] = cell->getTemp();
        p00[PRESSURE] = cell->getPressure();
        p00[E] = cell->getE(false);
        p00[GRAVITY] = cell->getGravPotential();
        // p00[R] = cell->getColor()[0]*p00[MASS];
        // p00[G] = cell->getColor()[1]*p00[MASS];
        // p00[B] = cell->getColor()[2]*p00[MASS];
        p00[R] = cell->getColor()[0];
        p00[G] = cell->getColor()[1];
        p00[B] = cell->getColor()[2];
        p00[HYDROGEN] = cell->getX();
        p00[HELIUM] = cell->getY();
    } else {
        if (compressible && 0) {
            FluidCell *cell = getCell(i,j+1);
            p00[SMOKEMASS] = cell->getMass();
            p00[MASS] = cell->getMass(false);
            p00[TEMPERATURE] = cell->getTemp();
            p00[PRESSURE] = cell->getPressure();
            p00[E] = cell->getE(false);
            p00[GRAVITY] = cell->getGravPotential();
            p00[R] = cell->getColor()[0];
            p00[G] = cell->getColor()[1];
            p00[B] = cell->getColor()[2];
            p00[HYDROGEN] = cell->getX();
            p00[HELIUM] = cell->getY();
        } else {
            p00[SMOKEMASS] = 0;
            p00[MASS] = 0;
            p00[TEMPERATURE] = 0;
            p00[PRESSURE] = 0;
            p00[E] = 0;
            p00[R] = 0;
            p00[G] = 0;
            p00[B] = 0;
            p00[HYDROGEN] = 0;
            p00[HELIUM] = 0;
        }
    }
    if ((j < cols-1) && (i > 0)) {
        // p01 = vyArray[i-1][j+1];
        FluidCell *cell = getCell(i-1,j+1);
        p01[SMOKEMASS] = cell->getMass();
        p01[MASS] = cell->getMass(false);
        p01[TEMPERATURE] = cell->getTemp();
        p01[PRESSURE] = cell->getPressure();
        p01[E] = cell->getE(false);
        p01[GRAVITY] = cell->getGravPotential();
        // p01[R] = cell->getColor()[0]*p01[MASS];
        // p01[G] = cell->getColor()[1]*p01[MASS];
        // p01[B] = cell->getColor()[2]*p01[MASS];
        p01[R] = cell->getColor()[0];
        p01[G] = cell->getColor()[1];
        p01[B] = cell->getColor()[2];
        p01[HYDROGEN] = cell->getX();
        p01[HELIUM] = cell->getY();
    } else {
        if (compressible && 0) {
            FluidCell *cell = getCell(i,j);
            p01[SMOKEMASS] = cell->getMass();
            p01[MASS] = cell->getMass(false);
            p01[TEMPERATURE] = cell->getTemp();
            p01[PRESSURE] = cell->getPressure();
            p01[E] = cell->getE(false);
            p01[GRAVITY] = cell->getGravPotential();
            p01[R] = cell->getColor()[0];
            p01[G] = cell->getColor()[1];
            p01[B] = cell->getColor()[2];
            p01[HYDROGEN] = cell->getX();
            p01[HELIUM] = cell->getY();
        } else {
            p01[SMOKEMASS] = 0;
            p01[MASS] = 0;
            p01[TEMPERATURE] = 0;
            p01[PRESSURE] = 0;
            p01[E] = 0;
            p01[R] = 0;
            p01[G] = 0;
            p01[B] = 0;
            p01[HYDROGEN] = 0;
            p01[HELIUM] = 0;
        }
    }
    if ((j < cols-1) && (i < rows)) {
        // p11 = vyArray[i][j+1];
        FluidCell *cell = getCell(i,j+1);
        p11[SMOKEMASS] = cell->getMass();
        p11[MASS] = cell->getMass(false);
        p11[TEMPERATURE] = cell->getTemp();
        p11[PRESSURE] = cell->getPressure();
        p11[E] = cell->getE(false);
        p11[GRAVITY] = cell->getGravPotential();
        // p11[R] = cell->getColor()[0]*p11[MASS];
        // p11[G] = cell->getColor()[1]*p11[MASS];
        // p11[B] = cell->getColor()[2]*p11[MASS];
        p11[R] = cell->getColor()[0];
        p11[G] = cell->getColor()[1];
        p11[B] = cell->getColor()[2];
        p11[HYDROGEN] = cell->getX();
        p11[HELIUM] = cell->getY();
    } else {
        if (compressible && 0) {
            FluidCell *cell = getCell(i-1,j);
            p11[SMOKEMASS] = cell->getMass();
            p11[MASS] = cell->getMass(false);
            p11[TEMPERATURE] = cell->getTemp();
            p11[PRESSURE] = cell->getPressure();
            p11[E] = cell->getE(false);
            p11[GRAVITY] = cell->getGravPotential();
            p11[R] = cell->getColor()[0];
            p11[G] = cell->getColor()[1];
            p11[B] = cell->getColor()[2];
            p11[HYDROGEN] = cell->getX();
            p11[HELIUM] = cell->getY();
        } else {
            p11[SMOKEMASS] = 0;
            p11[MASS] = 0;
            p11[TEMPERATURE] = 0;
            p11[PRESSURE] = 0;
            p11[E] = 0;
            p11[R] = 0;
            p11[G] = 0;
            p11[B] = 0;
            p11[HYDROGEN] = 0;
            p11[HELIUM] = 0;
        }
    }
    // printf("%f %f %f %f\n", p10[MASS], p00[MASS], p01[MASS], p11[MASS]);
    std::map<uint, double> props;
    props[SMOKEMASS] = (1-X)*(1-Y)*p10[SMOKEMASS] + X*(1-Y)*p11[SMOKEMASS] + (1-X)*Y*p00[SMOKEMASS] + X*Y*p01[SMOKEMASS];
    
    props[MASS] = (1-X)*(1-Y)*p10[MASS] + X*(1-Y)*p11[MASS] + (1-X)*Y*p00[MASS] + X*Y*p01[MASS];
    // printf("%f\n", props[MASS]);
    props[TEMPERATURE] = (1-X)*(1-Y)*p10[TEMPERATURE] + X*(1-Y)*p11[TEMPERATURE] + (1-X)*Y*p00[TEMPERATURE] + X*Y*p01[TEMPERATURE];
    props[PRESSURE] = (1-X)*(1-Y)*p10[PRESSURE] + X*(1-Y)*p11[PRESSURE] + (1-X)*Y*p00[PRESSURE] + X*Y*p01[PRESSURE];
    props[E] = (1-X)*(1-Y)*p10[E] + X*(1-Y)*p11[E] + (1-X)*Y*p00[E] + X*Y*p01[E];
    props[R] = (1-X)*(1-Y)*p10[R] + X*(1-Y)*p11[R] + (1-X)*Y*p00[R] + X*Y*p01[R];
    props[G] = (1-X)*(1-Y)*p10[G] + X*(1-Y)*p11[G] + (1-X)*Y*p00[G] + X*Y*p01[G];
    props[B] = (1-X)*(1-Y)*p10[B] + X*(1-Y)*p11[B] + (1-X)*Y*p00[B] + X*Y*p01[B];
    props[HYDROGEN] = (1-X)*(1-Y)*p10[HYDROGEN] + X*(1-Y)*p11[HYDROGEN] + (1-X)*Y*p00[HYDROGEN] + X*Y*p01[HYDROGEN];
    props[HELIUM] = (1-X)*(1-Y)*p10[HELIUM] + X*(1-Y)*p11[HELIUM] + (1-X)*Y*p00[HELIUM] + X*Y*p01[HELIUM];
    // props[R] = getCell(i,j)->getColor()[0];
    // props[G] = getCell(i,j)->getColor()[1];
    // props[B] = getCell(i,j)->getColor()[2];
    // if ((i == rows/2) && (j == cols/2))
    // printf("%f %f %f %f\n", p10[R],p00[R],p01[R],p11[R]);

    return props;
}

std::map<uint, double> FluidGrid::sampleCellAtPoint(position xy) {
    double x = xy.x;
    double y = xy.y;
    int i,j;
    if (round(1000*(height - y - cellHeight/2)) == 0) {
        i = 0;
    } else {
        i = ceil((height - y - cellHeight/2) / cellHeight);
    }
    if (round(1000*(x - cellWidth/2)) == 0) {
        j = 0;
    } else {
        j = floor((x - cellWidth/2) / cellWidth);
    }
    
    if ((i == rows/2) && (j == cols/2)) {
        // printf("%d: %f %d: %f\n", j,x, i,y);
        // printf("%f %f\n", prevY, physY);
    }
    // switch (prop) {
    //     case MASS:

    // }

    // the physical positions relative to nearest box of cells
    double X, Y; 
    X = fmod(x+cellWidth/2, cellWidth)/cellWidth;
    Y = fmod(y+cellHeight/2, cellHeight)/cellHeight;

    // will throw an error for when x and y are off the grid entirely
    std::map<uint, double> p10, p00, p01, p11;
    if ((j > -1) && (i < rows)) {
        // printf("%d %d\n", j, i);
        // p10 = vyArray[i][j];
        FluidCell *cell = getCell(i,j);
        p10[SMOKEMASS] = cell->getMass();
        p10[MASS] = cell->getMass(false);
        p10[TEMPERATURE] = cell->getTemp();
        p10[PRESSURE] = cell->getPressure();
        p10[E] = cell->getE(false);
        p10[GRAVITY] = cell->getGravPotential();
        // p10[R] = cell->getColor()[0]*p10[MASS];
        // p10[G] = cell->getColor()[1]*p10[MASS];
        // p10[B] = cell->getColor()[2]*p10[MASS];
        p10[R] = cell->getColor()[0];
        p10[G] = cell->getColor()[1];
        p10[B] = cell->getColor()[2];
        p10[HYDROGEN] = cell->getX();
        p10[HELIUM] = cell->getY();
    } else {
        if (compressible && 0) {
            FluidCell *cell = getCell(i-1,j+1);
            p10[SMOKEMASS] = cell->getMass();
            p10[MASS] = cell->getMass(false);
            p10[TEMPERATURE] = cell->getTemp();
            p10[PRESSURE] = cell->getPressure();
            p10[E] = cell->getE(false);
            p10[GRAVITY] = cell->getGravPotential();
            p10[R] = cell->getColor()[0];
            p10[G] = cell->getColor()[1];
            p10[B] = cell->getColor()[2];
            p10[HYDROGEN] = cell->getX();
            p10[HELIUM] = cell->getY();
        } else {
            p10[SMOKEMASS] = 0;
            p10[MASS] = 0;
            p10[TEMPERATURE] = 0;
            p10[PRESSURE] = 0;
            p10[E] = 0;
            p10[R] = 0;
            p10[G] = 0;
            p10[B] = 0;
            p10[HYDROGEN] = 0;
            p10[HELIUM] = 0;
        }
    }
    if ((j > -1) && (i > 0)) {
        // p00 = vyArray[i-1][j];
        FluidCell *cell = getCell(i-1,j);
        p00[SMOKEMASS] = cell->getMass();
        p00[MASS] = cell->getMass(false);
        p00[TEMPERATURE] = cell->getTemp();
        p00[PRESSURE] = cell->getPressure();
        p00[E] = cell->getE(false);
        p00[GRAVITY] = cell->getGravPotential();
        // p00[R] = cell->getColor()[0]*p00[MASS];
        // p00[G] = cell->getColor()[1]*p00[MASS];
        // p00[B] = cell->getColor()[2]*p00[MASS];
        p00[R] = cell->getColor()[0];
        p00[G] = cell->getColor()[1];
        p00[B] = cell->getColor()[2];
        p00[HYDROGEN] = cell->getX();
        p00[HELIUM] = cell->getY();
    } else {
        if (compressible && 0) {
            FluidCell *cell = getCell(i,j+1);
            p00[SMOKEMASS] = cell->getMass();
            p00[MASS] = cell->getMass(false);
            p00[TEMPERATURE] = cell->getTemp();
            p00[PRESSURE] = cell->getPressure();
            p00[E] = cell->getE(false);
            p00[GRAVITY] = cell->getGravPotential();
            p00[R] = cell->getColor()[0];
            p00[G] = cell->getColor()[1];
            p00[B] = cell->getColor()[2];
            p00[HYDROGEN] = cell->getX();
            p00[HELIUM] = cell->getY();
        } else {
            p00[SMOKEMASS] = 0;
            p00[MASS] = 0;
            p00[TEMPERATURE] = 0;
            p00[PRESSURE] = 0;
            p00[E] = 0;
            p00[R] = 0;
            p00[G] = 0;
            p00[B] = 0;
            p00[HYDROGEN] = 0;
            p00[HELIUM] = 0;
        }
    }
    if ((j < cols-1) && (i > 0)) {
        // p01 = vyArray[i-1][j+1];
        FluidCell *cell = getCell(i-1,j+1);
        p01[SMOKEMASS] = cell->getMass();
        p01[MASS] = cell->getMass(false);
        p01[TEMPERATURE] = cell->getTemp();
        p01[PRESSURE] = cell->getPressure();
        p01[E] = cell->getE(false);
        p01[GRAVITY] = cell->getGravPotential();
        // p01[R] = cell->getColor()[0]*p01[MASS];
        // p01[G] = cell->getColor()[1]*p01[MASS];
        // p01[B] = cell->getColor()[2]*p01[MASS];
        p01[R] = cell->getColor()[0];
        p01[G] = cell->getColor()[1];
        p01[B] = cell->getColor()[2];
        p01[HYDROGEN] = cell->getX();
        p01[HELIUM] = cell->getY();
    } else {
        if (compressible && 0) {
            FluidCell *cell = getCell(i,j);
            p01[SMOKEMASS] = cell->getMass();
            p01[MASS] = cell->getMass(false);
            p01[TEMPERATURE] = cell->getTemp();
            p01[PRESSURE] = cell->getPressure();
            p01[E] = cell->getE(false);
            p01[GRAVITY] = cell->getGravPotential();
            p01[R] = cell->getColor()[0];
            p01[G] = cell->getColor()[1];
            p01[B] = cell->getColor()[2];
            p01[HYDROGEN] = cell->getX();
            p01[HELIUM] = cell->getY();
        } else {
            p01[SMOKEMASS] = 0;
            p01[MASS] = 0;
            p01[TEMPERATURE] = 0;
            p01[PRESSURE] = 0;
            p01[E] = 0;
            p01[R] = 0;
            p01[G] = 0;
            p01[B] = 0;
            p01[HYDROGEN] = 0;
            p01[HELIUM] = 0;
        }
    }
    if ((j < cols-1) && (i < rows)) {
        // p11 = vyArray[i][j+1];
        FluidCell *cell = getCell(i,j+1);
        p11[SMOKEMASS] = cell->getMass();
        p11[MASS] = cell->getMass(false);
        p11[TEMPERATURE] = cell->getTemp();
        p11[PRESSURE] = cell->getPressure();
        p11[E] = cell->getE(false);
        p11[GRAVITY] = cell->getGravPotential();
        // p11[R] = cell->getColor()[0]*p11[MASS];
        // p11[G] = cell->getColor()[1]*p11[MASS];
        // p11[B] = cell->getColor()[2]*p11[MASS];
        p11[R] = cell->getColor()[0];
        p11[G] = cell->getColor()[1];
        p11[B] = cell->getColor()[2];
        p11[HYDROGEN] = cell->getX();
        p11[HELIUM] = cell->getY();
    } else {
        if (compressible && 0) {
            FluidCell *cell = getCell(i-1,j);
            p11[SMOKEMASS] = cell->getMass();
            p11[MASS] = cell->getMass(false);
            p11[TEMPERATURE] = cell->getTemp();
            p11[PRESSURE] = cell->getPressure();
            p11[E] = cell->getE(false);
            p11[GRAVITY] = cell->getGravPotential();
            p11[R] = cell->getColor()[0];
            p11[G] = cell->getColor()[1];
            p11[B] = cell->getColor()[2];
            p11[HYDROGEN] = cell->getX();
            p11[HELIUM] = cell->getY();
        } else {
            p11[SMOKEMASS] = 0;
            p11[MASS] = 0;
            p11[TEMPERATURE] = 0;
            p11[PRESSURE] = 0;
            p11[E] = 0;
            p11[R] = 0;
            p11[G] = 0;
            p11[B] = 0;
            p11[HYDROGEN] = 0;
            p11[HELIUM] = 0;
        }
    }
    // printf("%f %f %f %f\n", p10[MASS], p00[MASS], p01[MASS], p11[MASS]);
    std::map<uint, double> props;
    props[SMOKEMASS] = (1-X)*(1-Y)*p10[SMOKEMASS] + X*(1-Y)*p11[SMOKEMASS] + (1-X)*Y*p00[SMOKEMASS] + X*Y*p01[SMOKEMASS];
    
    props[MASS] = (1-X)*(1-Y)*p10[MASS] + X*(1-Y)*p11[MASS] + (1-X)*Y*p00[MASS] + X*Y*p01[MASS];
    // printf("%f\n", props[MASS]);
    props[TEMPERATURE] = (1-X)*(1-Y)*p10[TEMPERATURE] + X*(1-Y)*p11[TEMPERATURE] + (1-X)*Y*p00[TEMPERATURE] + X*Y*p01[TEMPERATURE];
    props[PRESSURE] = (1-X)*(1-Y)*p10[PRESSURE] + X*(1-Y)*p11[PRESSURE] + (1-X)*Y*p00[PRESSURE] + X*Y*p01[PRESSURE];
    props[E] = (1-X)*(1-Y)*p10[E] + X*(1-Y)*p11[E] + (1-X)*Y*p00[E] + X*Y*p01[E];
    props[R] = (1-X)*(1-Y)*p10[R] + X*(1-Y)*p11[R] + (1-X)*Y*p00[R] + X*Y*p01[R];
    props[G] = (1-X)*(1-Y)*p10[G] + X*(1-Y)*p11[G] + (1-X)*Y*p00[G] + X*Y*p01[G];
    props[B] = (1-X)*(1-Y)*p10[B] + X*(1-Y)*p11[B] + (1-X)*Y*p00[B] + X*Y*p01[B];
    props[HYDROGEN] = (1-X)*(1-Y)*p10[HYDROGEN] + X*(1-Y)*p11[HYDROGEN] + (1-X)*Y*p00[HYDROGEN] + X*Y*p01[HYDROGEN];
    props[HELIUM] = (1-X)*(1-Y)*p10[HELIUM] + X*(1-Y)*p11[HELIUM] + (1-X)*Y*p00[HELIUM] + X*Y*p01[HELIUM];
    // props[R] = getCell(i,j)->getColor()[0];
    // props[G] = getCell(i,j)->getColor()[1];
    // props[B] = getCell(i,j)->getColor()[2];
    // if ((i == rows/2) && (j == cols/2))
    // printf("%f %f %f %f\n", p10[R],p00[R],p01[R],p11[R]);

    return props;
}

std::map<uint, double> FluidGrid::propsAtij(int i, int j) {
    std::map<uint, double> props;
    uint onGrid = ((i < rows) && (j < cols)) && ((i > -1) && (j > -1));
    if (onGrid) {
        FluidCell *cell = getCell(i,j);
        
        props[SMOKEMASS] = cell->getMass();
        props[MASS] = cell->getMass(false);
        props[TEMPERATURE] = cell->getTemp();
        props[PRESSURE] = cell->getPressure();
        props[E] = cell->getE(false);
        props[GRAVITY] = cell->getGravPotential();
        props[R] = cell->getColor()[0];
        props[G] = cell->getColor()[1];
        props[B] = cell->getColor()[2];
        props[HYDROGEN] = cell->getX();
        props[HELIUM] = cell->getY();
    } else {
        props[SMOKEMASS] = 0;
        props[MASS] = 0;
        props[TEMPERATURE] = 0;
        props[PRESSURE] = 0;
        props[E] = 0;
        props[R] = 0;
        props[G] = 0;
        props[B] = 0;
        props[HYDROGEN] = 0;
        props[HELIUM] = 0;
    }
    return props;
}

std::map<char, double> FluidGrid::getxyFromij(int i, int j) {
    std::map<char, double> xy;
    xy['x'] = j*cellWidth + cellWidth/2;
    xy['y'] = height - (i*cellHeight + cellHeight/2);
    return xy;
}

void FluidGrid::diffuse(int iters) {
    if (viscosity > 0) {
        // diffuse velocities, densities, and temperatures using gauss-seidel
        int i, j;
        double x, y;
        VelocityVector c, b, t, l, r; // center bottom top left right
        std::map<char, double> cp, bp, tp, lp, rp;
        double vxC, vxB, vxT, vxR, vxL;
        double vyC, vyB, vyT, vyR, vyL;
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
                    // double dvx = dt*viscosity*((vxL + vxR - 2*vxC)/(cellWidth*cellWidth)+(vxB + vxT - 2*vxC)/(cellHeight*cellHeight));
                    double vxNew = (vxC+viscosity*((vxL+vxT+vxR+vxB)/4))/(1+viscosity);

                    // solving vys:
                    vyC = c.getVy();
                    vyB = b.getVy();
                    vyT = t.getVy();
                    vyR = r.getVy();
                    vyL = l.getVy();
                    double vyNew = (vyC+viscosity*((vyL+vyT+vyR+vyB)/4))/(1+viscosity);
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
            std::map<char,double> xyLeft = vGrid->getxyFromijVx(i,j);
            std::map<char,double> xyTop = vGrid->getxyFromijVy(i,j);
            std::map<char,double> xyRight = vGrid->getxyFromijVx(i,j+1);
            std::map<char,double> xyBottom = vGrid->getxyFromijVy(i+1,j);
            double massL = sampleCellAtPoint(xyLeft['x'],xyLeft['y'])[MASS];
            double massT = sampleCellAtPoint(xyTop['x'],xyTop['y'])[MASS];
            double massR = sampleCellAtPoint(xyRight['x'],xyRight['y'])[MASS];
            double massB = sampleCellAtPoint(xyBottom['x'],xyBottom['y'])[MASS];
            double vxL = vGrid->getVx(i,j);
            double vyT = vGrid->getVy(i,j);
            double vxR = vGrid->getVx(i,j+1);
            double vyB = vGrid->getVy(i+1,j);
            double divMass = (-vxL*massL + vyT*massT + vxR*massR - vyB*massB);
            FluidCell *c = getCell(i,j);
            c->setMass(c->getMass(false) - divMass*dt,false);
        }
    }
}

void FluidGrid::solveGravPotential(int iters) {
    int i, j, n;
    for (n = 0; n < iters; n++) {
        for (i = 0; i < rows; i++) {
            for (j = 0; j < cols; j++) {
                if (i == 0 || j == 0 || i == rows-1 || j == cols-1) {
                    getCell(i,j)->setGravPotential(0);
                    continue;
                }
                double uL, uR, uT, uB = 0;
                double dens = getCell(i,j)->getDensity();
                double w2 = cellWidth*cellWidth;
                double h2 = cellHeight*cellHeight;
                if (i > 0) uT = getCell(i-1,j)->getGravPotential();
                if (i < rows-1) uB = getCell(i+1,j)->getGravPotential();
                if (j > 0) uL = getCell(i,j-1)->getGravPotential();
                if (j < cols-1) uR = getCell(i,j+1)->getGravPotential();
                double u = (uL + uR)/(2*(1+w2/h2)) + (uT + uB)/(2*(h2/w2+1)) - 2 * consts::PI * consts::G * dens / (1/w2 + 1/h2);
                getCell(i,j)->setGravPotential(u);
            }
        }
    }
    // for (i = 0; i < rows; ++i) {
    //     getCell(i, 0)->setGravPotential(0); // Left boundary
    //     getCell(i, cols - 1)->setGravPotential(0); // Right boundary
    // }
    // for (j = 0; j < cols; ++j) {
    //     getCell(0, j)->setGravPotential(0); // Top boundary
    //     getCell(rows - 1, j)->setGravPotential(0); // Bottom boundary
    // }
}

void FluidGrid::compressibleMomentumUpdate() {
    int i, j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            FluidCell *cell = getCell(i,j);
            double mass = cell->getMass(false);
            if (gravityFlag) force(i,j,0,mass*(consts::g));
            // velocities surrounding this cell
            double left = vGrid->getVx(i, j);
            double right = vGrid->getVx(i, j+1);
            double top = vGrid->getVy(i, j);
            double bottom = vGrid->getVy(i+1, j);

            double pressure = cell->getPressure();

            double E = cell->getE(false);
            
            

            // // pressure-work term
            // double div = (right-left)/cellWidth + (top-bottom)/cellHeight;
            // double pdiv = pressure * div;
            // // if (i==1 && j==cols/2) printf("pdiv:%f\n",pdiv);
            // E += -pdiv*getdt();
            
            // gravitational-work term
            // if (gravityFlag) {
            //     double gravWork = cell->getVelocity().getVy()*(consts::g)*cell->getDensity();
            //     E += gravWork*getdt();
            // }
            // cell->setE(E);
            // auto start = std::chrono::steady_clock::now();
            // viscosity term
            std::map<char,double> xyM, xyB, xyL, xyT, xyR;
            double laplacianVx, laplacianVy;
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

            // auto visc = std::chrono::steady_clock::now();

            double grav = getCell(i,j)->getGravPotential();

            
            // if (i < rows-1) {
            //     double pressureB = getCell(i+1,j)->getPressure();
            //     double gravB = getCell(i+1,j)->getGravPotential();
            //     double gradYP = (pressure-pressureB)/(2*cellHeight);
            //     double gradYG = -(grav-gravB)/(2*cellHeight);
            //     // if (i == rows/2 && j ==cols/2) printf("gradYG:%f gradYP: %f dens:%f GP:%f press:%f\n", gradYG, gradYP, cell->getDensity(), cell->getGravPotential(), cell->getPressure(false));
            //     std::map<char,double> xyVy = vGrid->getxyFromijVy(i+1,j);
            //     double mass = sampleCellAtPoint(xyVy['x'],xyVy['y'])[MASS];
            //     double density = mass/(cellHeight*cellWidth);
            //     // double newBVy = bottom-gradY;
            //     double newBVy = bottom-(1/density)*gradYP*dt+gradYG*dt;
            //     // if (i == rows-2 && j ==cols/2) printf("newVy:%f; gradY:%f; 1/dens:%f\n",newBVy,gradY,1/density);
            //     new_vGrid->setVy(i+1,j,newBVy);
            // }
            double density;
            double pressureB, gravB;
            if (i < rows-1) {
                pressureB = getCell(i+1,j)->getPressure();
                gravB = getCell(i+1,j)->getGravPotential();
            } else {
                gravB = 0;
                pressureB = 0;
            }
            double gradYP = (pressure-pressureB)/(2*cellHeight);
            double gradYG = -(grav-gravB)/(2*cellHeight);
            std::map<char,double> xyVy = vGrid->getxyFromijVy(i+1,j);
            // mass = sampleCellAtPoint(xyVy['x'],xyVy['y'])[MASS];
            if (i < rows-1) {
                mass = (cell->getMass(false) + getCell(i+1,j)->getMass(false))/2;
            } else {
                mass = cell->getMass(false)/2;
            }
            density = mass/(cellHeight*cellWidth);
            double newBVy;
            if (round(density) == 0) {
                newBVy = 0;
            } else {
                newBVy = bottom-(1/density)*gradYP*dt+gradYG*dt;
            }
            new_vGrid->setVy(i+1,j,newBVy);
            // auto GPy = std::chrono::steady_clock::now();

            // if (j < cols-1) {
            //     double pressureR = getCell(i,j+1)->getPressure();
            //     double gravR = getCell(i,j+1)->getGravPotential();
            //     double gradXP = (pressureR-pressure)/(2*cellWidth);
            //     double gradXG = -(gravR-grav)/(2*cellWidth);
            //     if (!std::isfinite(gradXP)) {
            //         std::cerr << "NaN detected in grad X at cell (" << i << ", " << j << ") pressure(i,j): " << pressure << " density: " << cell->getDensity() << "\n";
            //         assert(false); // Debug immediately
            //     }
            //     // if (i == rows/2 && j ==cols/2) printf("gradX:%f\n",gradX);
            //     std::map<char,double> xyVx = vGrid->getxyFromijVx(i,j+1);
            //     double density = sampleCellAtPoint(xyVx['x'],xyVx['y'])[MASS]/(cellHeight*cellWidth);
            //     double newRVx = right-(1/density)*gradXP*dt+gradXG*dt;
            //     if (!std::isfinite(newRVx)) {
            //         std::cerr << "NaN detected in vxR at cell (" << i << ", " << j << ") old vxR: " << right << " density: " << density << "\n";
            //         assert(false); // Debug immediately
            //     }
            //     new_vGrid->setVx(i,j+1,newRVx);
            //     // if (i == rows/2 && j ==cols/2) printf("newVx:%f\n",newRVx);
            // }
            double pressureR, gravR;
            if (j < cols-1) {
                pressureR = getCell(i,j+1)->getPressure(false);
                gravR = getCell(i,j+1)->getGravPotential();
            } else {
                pressureR = 0;
                gravR = 0;
            }
            double gradXP = (pressureR-pressure)/(2*cellWidth);
            double gradXG = -(gravR-grav)/(2*cellWidth);
            // if (i == rows/2 && j ==cols/2) printf("gradX:%f\n",gradX);
            std::map<char,double> xyVx = vGrid->getxyFromijVx(i,j+1);
            // mass = sampleCellAtPoint(xyVx['x'],xyVx['y'])[MASS];
            if (j < cols-1) {
                mass = (cell->getMass(false) + getCell(i,j+1)->getMass(false))/2;
            } else {
                mass = cell->getMass(false)/2;
            }
            density = mass/(cellHeight*cellWidth);
            double newRVx;
            if (round(density) == 0) {
                newRVx = 0;
            } else {
                newRVx = right-(1/density)*gradXP*dt+gradXG*dt;
            }
            new_vGrid->setVx(i,j+1,newRVx);
            // auto GPx = std::chrono::steady_clock::now();

            // std::cout << "visc: ";
            // std::cout << std::chrono::duration_cast<std::chrono::microseconds>(visc - start).count() << "\n";
            // std::cout << "vy: ";
            // std::cout << std::chrono::duration_cast<std::chrono::microseconds>(GPy - visc).count() << "\n";
            // std::cout << "vx: ";
            // std::cout << std::chrono::duration_cast<std::chrono::microseconds>(GPx - GPy).count() << "\n";
            // if (i == rows/2 && j ==cols/2) printf("newVx:%f\n",newRVx);
        
        }
    }
    vGrid = new_vGrid;
}

void FluidGrid::energyUpdate() {
    int i, j;
    // #pragma omp parallel for collapse(2)
    // #pragma omp parallel for schedule(dynamic, 4)
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            FluidCell *cell = getCell(i,j);
            double mass = cell->getMass(false);
            // velocities surrounding this cell
            double left = vGrid->getVx(i, j);
            double right = vGrid->getVx(i, j+1);
            double top = vGrid->getVy(i, j);
            double bottom = vGrid->getVy(i+1, j);

            double pressure = cell->getPressure();
            double E = cell->getE(false);
            // nuclear fusion
            double fusionEnergy = fusionRate(cell->getDensity(), cell->getTemp(false), cell->getX()) * getdt();
            if (round(fusionEnergy) > 0 && i == rows/2 && j==cols/2) printf("fusionEnergy %f\n", fusionEnergy);
            E += fusionEnergy;

            // Update hydrogen and helium abundances
            double dX = fusionEnergy / (consts::QH_He * cell->getDensity());
            cell->HtoHe(dX);
            
            // pressure-work term
            double div = (right-left)/cellWidth + (top-bottom)/cellHeight;
            double pdiv = pressure * div;
            if (E - pdiv * getdt() >= 0) {
                E += -pdiv * getdt();
            } else {
                // std::cerr << "Warning: Pressure work term causing negative energy.\n";
                E = 0.0f; // Prevent negative energy
            }
            // if (i==1 && j==cols/2) printf("pdiv:%f\n",pdiv);
            // E += -pdiv*getdt();
            double density = cell->getDensity();
            double T = cell->getTemp(false);
            double Z = cell->getZ();
            double lambda = coolingFunction(density,T,Z);
            if (E + lambda * getdt() >= 0) {
                E += -lambda*getdt();
            } else {
                E = 0.0f;
            }
            // gravitational-work term
            double grav = getCell(i,j)->getGravPotential();
            double gravB = getCell(i+1,j)->getGravPotential();
            double gy = -(grav-gravB)/(2*cellHeight);
            double gravR = getCell(i,j+1)->getGravPotential();
            double gx = -(gravR-grav)/(2*cellWidth);
            double gravWork = (cell->getVelocity().getVx()*gx+cell->getVelocity().getVy()*gy)*cell->getDensity();
            if (E + gravWork * getdt() >= 0) {
                E += gravWork * getdt();
            } else {
                // std::cerr << "Warning: Gravitational work term causing negative energy.\n";
                E = 0.0f; // Prevent negative energy
            }
            // if (gravityFlag) {
            //     double gravWork = cell->getVelocity().getVy()*(consts::g)*cell->getDensity();
            //     // E += gravWork*getdt();
            //     if (E + gravWork * getdt() >= 0) {
            //         E += gravWork * getdt();
            //     } else {
            //         std::cerr << "Warning: Gravitational work term causing negative energy.\n";
            //         E = 0.0f; // Prevent negative energy
            //     }
            // }
            cell->setE(E);
        }
    }
}

void FluidGrid::projection(int iters) {
    // for incompressible fluids only
    int i, j;
    for (int n = 0; n < iters+1; n++) {
        // for (i = 0; i < rows; i++) {
        //     for (j = 0; j < cols; j++) {
        for (i = 0; i < rows; i++) {
            for (j = cols-1; j >= 0; j--) {
                FluidCell *cell = getCell(i,j);
                double mass = cell->getMass(false);
                if (n == 0) {
                    if (gravityFlag) force(i,j,0,mass*(consts::g));
                    // double pressure = cell->getDensity()/consts::Hmol * consts::kb * consts::Na * cell->getTemp();
                    cell->setPressure(0);
                }
                if (n > 0) {
                    int sLeft = !!j;
                    int sRight = !!(cols - 1 - j);
                    int sTop = !!i;
                    int sBottom = !!(rows - 1 - i);                        

                    int s = sLeft + sRight + sTop + sBottom;

                    double left = vGrid->getVx(i, j);
                    double right = vGrid->getVx(i, j+1);
                    double top = vGrid->getVy(i, j);
                    double bottom = vGrid->getVy(i+1, j);
                    
                    double d = (-left*sLeft + right*sRight + top*sTop - bottom*sBottom)*1.5;
                    // if (n == iters-1 && i == rows/2 && j == cols/2) printf("divergence: %f\n", d);
                    // if (i == rows/2 && j ==cols/2) printf("div:%f\n", d);
                    double newPressure = cell->getPressure(false) - d/s * cell->getDensity()*std::sqrt(cellWidth*cellHeight)/dt;
                    cell->setPressure(newPressure);

                    
                    vGrid->setVx(i, j, left+d*sLeft/s);
                    vGrid->setVx(i, j+1, right-d*sRight/s);
                    vGrid->setVy(i, j, top-d*sTop/s);
                    vGrid->setVy(i+1, j, bottom+d*sBottom/s);
                }
            }
        }
    }
}

double FluidGrid::calculateFlux(double quantity, int i, int j) {
    double vxL = vGrid->getVx(i,j);
    double vxR = vGrid->getVx(i,j+1);
    double vyT = vGrid->getVy(i,j);
    double vyB = vGrid->getVy(i+1,j);

    double flux = quantity*(vxL*cellHeight-vxR*cellHeight-vyT*cellWidth+vyB*cellWidth);
    return flux;
}

// std::map<uint,double> FluidGrid::calculateFlux(std::map<uint,double> propDict, int i, int j, double v) {
//     std::map<uint,double> newPropDict;
//     for (auto const& [property, quantity] : propDict) {
//         newPropDict[property] = quantity*v;
//     }
//     return newPropDict;
// }

void FluidGrid::advect() {
    int i, j;
    for (i = 0; i < rows+1; i++) {
        for (j = 0; j < cols+1; j++) {
            if ((i < rows) && (j < cols)) {
                FluidCell cell = *getCell(i,j);
                double physX = cellWidth*j+cellWidth/2;
                double physY = height - (cellHeight*i+cellHeight/2); // I want y pointed up
                VelocityVector vCell = vGrid->sampleVelocityAtPoint(physX, physY);

                double prevX = physX - dt*vCell.getVx();
                double prevY = physY - dt*vCell.getVy();
                // if (i==rows/2 && j==5) printf("Vx-right: %f, prevX: %f\n",vCell.getVx(), prevX);
                // else if (i==rows/2 && j==cols-6) printf("Vx-left: %f, prevX: %f\n", vCell.getVx(), prevX);
                std::map<uint, double> prevProps;
                // if (round((vCell.getVx()+vCell.getVy())*dt/std::sqrt(cellHeight*cellWidth)) == 0) {
                if (round((vCell.getVx()+vCell.getVy())*100) == 0) {
                    prevProps[SMOKEMASS] = cell.getMass();
                    prevProps[MASS] = cell.getMass(false);
                    prevProps[TEMPERATURE] = cell.getTemp();
                    prevProps[PRESSURE] = cell.getPressure();
                    prevProps[E] = cell.getE(false);
                    prevProps[R] = cell.getColor()[0];
                    prevProps[G] = cell.getColor()[1];
                    prevProps[B] = cell.getColor()[2];
                    prevProps[HYDROGEN] = cell.getX();
                    prevProps[HELIUM] = cell.getY();
                } else {
                    // prevX = std::max(0.0f, std::min(prevX, width));
                    // prevY = std::max(0.0f, std::min(prevY, height));
                    if (round(prevX) < -0.0f) prevX = -prevX;
                    if (round(prevX) > width) prevX = 2 * width - prevX;

                    if (round(prevY) < -0.0f) prevY = -prevY;
                    if (round(prevY) > height) prevY = 2 * height - prevY;
                    prevProps = sampleCellAtPoint(prevX, prevY);
                    // if ((0 < prevX) && (prevX < width) && (0 < prevY) && (prevY < height)) {
                    //     prevProps = sampleCellAtPoint(prevX, prevY);
                    // } 
                }
                
                if (prevProps[E] < 0) {
                    // std::cerr << "Warning: Negative energy detected in advection.\n";
                    prevProps[E] = 0.0f;
                }

                newGrid[i][j].setMass(prevProps[SMOKEMASS]);
                // newGrid[i][j].setMass(prevProps[MASS]*massFactor,false);
                newGrid[i][j].setMass(prevProps[MASS],false);
                // newGrid[i][j].setTemp(prevProps[TEMPERATURE]);
                
                // newGrid[i][j].setE(prevProps[E]*massFactor);
                newGrid[i][j].setE(prevProps[E]);
                newGrid[i][j].setColor(prevProps[R], prevProps[G], prevProps[B]);
                newGrid[i][j].setX(prevProps[HYDROGEN]);
                newGrid[i][j].setY(prevProps[HELIUM]);

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
                //     // std::vector<double> c = cell.getColor();
                //     // newGrid[i][j].setColor(c[0],c[1],c[2]);
                // } else {
                //     newGrid[i][j].setMass(prevProps[SMOKEMASS]);
                //     newGrid[i][j].setMass(prevProps[MASS],false);
                //     newGrid[i][j].setTemp(prevProps[TEMPERATURE]);
                //     newGrid[i][j].setE(prevProps[E]);
                //     newGrid[i][j].setColor(prevProps[R], prevProps[G], prevProps[B]);
                // }
            }

            double newVx, newVy;
            // velocity self-advection:
            if (i < rows) {
                // vxs of vgrid
                double vxPhysX = cellWidth*j;
                double vxPhysY = height - cellHeight*i - cellHeight/2;
                
                VelocityVector vVx = vGrid->sampleVelocityAtPoint(vxPhysX, vxPhysY);
                double prevVxX = vxPhysX - dt*vVx.getVx();
                double prevVxY = vxPhysY - dt*vVx.getVy();
                if (round(prevVxX) < -0.0f) prevVxX = -prevVxX;
                if (round(prevVxX) > width) prevVxX = 2 * width - prevVxX;
                if (round(prevVxY) < -0.0f) prevVxY = -prevVxY;
                if (round(prevVxY) > height) prevVxY = 2 * height - prevVxY;
                // prevVxX = std::max(0.0f, std::min(prevVxX, width));
                // prevVxY = std::max(0.0f, std::min(prevVxY, height));
                if ((0 < prevVxX) && (prevVxX < width) && (0 < prevVxY) && (prevVxY < height)) {
                    newVx = vGrid->sampleVelocityAtPoint(prevVxX, prevVxY).getVx();
                    new_vGrid->setVx(i,j,newVx);
                    // find the max velocity in x direction for dt
                    if (std::abs(newVx) > maxV) {
                        maxV = std::abs(newVx);
                    }
                }
            }
            
            if (j < cols) {
                // vys of vgrid
                double vyPhysX = cellWidth*j + cellWidth/2;
                double vyPhysY = height - cellHeight*i;
                VelocityVector vVy = vGrid->sampleVelocityAtPoint(vyPhysX, vyPhysY);
                double prevVyX = vyPhysX - dt*vVy.getVx();
                double prevVyY = vyPhysY - dt*vVy.getVy();

                if (round(prevVyX) < -0.0f) prevVyX = -prevVyX;
                if (round(prevVyX) > width) prevVyX = 2 * width - prevVyX;
                if (round(prevVyY) < -0.0f) prevVyY = -prevVyY;
                if (round(prevVyY) > height) prevVyY = 2 * height - prevVyY;

                // prevVyX = std::max(0.0f, std::min(prevVyX, width));
                // prevVyY = std::max(0.0f, std::min(prevVyY, height));
                if ((0 < prevVyX) && (prevVyX < width) && (0 < prevVyY) && (prevVyY < height)) {
                    newVy = vGrid->sampleVelocityAtPoint(prevVyX, prevVyY).getVy();
                    // if ((i == 0) || (i == rows)) {
                    //     new_vGrid->setVy(i,j,0);
                    // } else {
                    //     new_vGrid->setVy(i,j,newVy);
                    // }
                    new_vGrid->setVy(i,j,newVy);
                    // find the max velocity in y direction for dt
                    if (std::abs(newVy) > maxV) {
                        maxV = std::abs(newVy);
                    }
                }
            }
            // if (i < rows) {
            //     if (j == 0) {
            //         if (round(100*vGrid->getVx(i,j)) != 0) printf("left bound vx: %f at (%d, %d)\n",vGrid->getVx(i,j),i,j);
            //     }
            //     if (j == cols) {
            //         if (round(100*vGrid->getVx(i,j)) != 0) printf("right bound vx: %f at (%d, %d)\n",vGrid->getVx(i,j),i,j);
            //     }
            // }
        }
    }
    grid = newGrid;
    vGrid = new_vGrid;
    // for (i = 0; i < rows; i++) {
    //     vGrid->setVx(i, 0, vGrid->getVx(i, 1));
    //     vGrid->setVx(i, cols, vGrid->getVx(i, cols-1));
    //     // double avgVx = 0.5 * (vGrid->getVx(i, 0) + vGrid->getVx(i, cols));
    //     // vGrid->setVx(i, 0, avgVx);
    //     // vGrid->setVx(i, cols, avgVx);
    //     // vGrid->setVx(i, 0, 0);
    //     // vGrid->setVx(i, cols, 0);
    // }
    // for (j = 0; j < cols; j++) {
    //     vGrid->setVy(0, j, vGrid->getVy(1, j));
    //     vGrid->setVy(rows, j, vGrid->getVy(rows-1, j));
    //     // vGrid->setVy(0, j, 0);
    //     // vGrid->setVy(rows, j, 0);
    // }

}

void FluidGrid::FFSLadvect() {
    int i, j;
    #pragma omp parallel for collapse(2)
    // #pragma omp parallel for schedule(dynamic, 4)
    for (i = 0; i < rows+1; i++) {
        for (j = 0; j < cols+1; j++) {
            if ((i < rows) && (j < cols)) {
                FluidCell cell = *getCell(i,j);
                std::map<char,double> xyT = vGrid->getxyFromijVy(i,j);
                std::map<char,double> xyL = vGrid->getxyFromijVx(i,j);
                std::map<char,double> xyB = vGrid->getxyFromijVy(i+1,j);
                std::map<char,double> xyR = vGrid->getxyFromijVx(i,j+1);
                // boundaries
                int sLeft = !!j;
                int sRight = !!(cols - 1 - j);
                int sTop = !!i;
                int sBottom = !!(rows - 1 - i);       
                double vyT = vGrid->getVy(i,j)*sTop;
                double vxL = vGrid->getVx(i,j)*sLeft;
                double vyB = vGrid->getVy(i+1,j)*sBottom;
                double vxR = vGrid->getVx(i,j+1)*sRight;


                position prevT, prevL, prevB, prevR;

                prevT.x = xyT['x'];
                prevT.y = xyT['y'] - dt*vyT;

                prevL.x = xyL['x'] - dt*vxL;
                prevL.y = xyL['y'];

                prevB.x = xyB['x'];
                prevB.y = xyB['y'] - dt*vyB;

                prevR.x = xyR['x'] - dt*vxR;
                prevR.y = xyR['y'];

                // if (i==rows/2 && j==5) printf("Vx-right: %f, prevX: %f\n",vCell.getVx(), prevX);
                // else if (i==rows/2 && j==cols-6) printf("Vx-left: %f, prevX: %f\n", vCell.getVx(), prevX);
                std::map<uint, double> prevPropsT, prevPropsL, prevPropsB, prevPropsR, deltaProps;
                
                // prevPropsT = sampleCellAtPoint(prevT);
                // prevPropsL = sampleCellAtPoint(prevL);
                // prevPropsB = sampleCellAtPoint(prevB);
                // prevPropsR = sampleCellAtPoint(prevR);

                prevPropsT = vyT > 0 ? propsAtij(i,j) : propsAtij(i-1,j);
                prevPropsL = vxL > 0 ? propsAtij(i,j-1) : propsAtij(i,j);
                prevPropsB = vyB > 0 ? propsAtij(i+1,j) : propsAtij(i,j);
                prevPropsR = vxR > 0 ? propsAtij(i,j) : propsAtij(i,j+1);
                double cellVolume = cellWidth*cellHeight;

                for(std::map<uint,double>::iterator it = prevPropsT.begin(); it != prevPropsT.end(); ++it) {
                    uint property = it->first;
                    if ((property == HYDROGEN) || (property == HELIUM)) {
                        // I want to directly advect the masses of hydrogen and helium
                        deltaProps[property] = (-prevPropsT[property]*prevPropsT[MASS]*vyT*cellWidth + prevPropsL[property]*prevPropsL[MASS]*vxL*cellHeight + prevPropsB[property]*prevPropsB[MASS]*vyB*cellWidth - prevPropsR[property]*prevPropsR[MASS]*vxR*cellHeight)*getdt()/cellVolume;
                    } else {
                        deltaProps[property] = (-prevPropsT[property]*vyT*cellWidth + prevPropsL[property]*vxL*cellHeight + prevPropsB[property]*vyB*cellWidth - prevPropsR[property]*vxR*cellHeight)*getdt()/cellVolume;
                    }
                }


                
                double nSMass = cell.getMass(true) + deltaProps[SMOKEMASS];
                double nMass = cell.getMass(false) + deltaProps[MASS] > 0 ? cell.getMass(false) + deltaProps[MASS] : 0;
                double nE = cell.getE(false) + deltaProps[E] > 0 ? cell.getE(false) + deltaProps[E] : 0;
                // if (i == 48 && j == 0) {
                //     printf("oldE: %f delE: %f newE: %f\n", cell.getE(false), deltaProps[E], nE);
                // }
                double nHydrogen = cell.getX()*cell.getMass(false) + deltaProps[HYDROGEN] > 0 ? cell.getX()*cell.getMass(false) + deltaProps[HYDROGEN] : 0;
                double nHelium = cell.getY()*cell.getMass(false) + deltaProps[HELIUM] > 0 ? cell.getY()*cell.getMass(false) + deltaProps[HELIUM] : 0;
                nHydrogen = nHydrogen/nMass;
                nHelium = nHelium/nMass;
                if (nHydrogen > 1) nHydrogen = 1;
                if (nHelium > 1) nHelium = 1;
                
                // newGrid[i][j].setMass(prevProps[SMOKEMASS]);
                // newGrid[i][j].setMass(prevProps[MASS],false);
                // // newGrid[i][j].setTemp(prevProps[TEMPERATURE]);
                
                // newGrid[i][j].setE(prevProps[E]);
                // newGrid[i][j].setColor(prevProps[R], prevProps[G], prevProps[B]);
                // newGrid[i][j].setX(prevProps[HYDROGEN]);
                // newGrid[i][j].setY(prevProps[HELIUM]);

                newGrid[i][j].setMass(nSMass);
                newGrid[i][j].setMass(nMass,false);
                // newGrid[i][j].setTemp(prevProps[TEMPERATURE]);
                
                newGrid[i][j].setE(nE);
                newGrid[i][j].setColor(deltaProps[R], deltaProps[G], deltaProps[B]);
                newGrid[i][j].setX(nHydrogen);
                newGrid[i][j].setY(nHelium);

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
                //     // std::vector<double> c = cell.getColor();
                //     // newGrid[i][j].setColor(c[0],c[1],c[2]);
                // } else {
                //     newGrid[i][j].setMass(prevProps[SMOKEMASS]);
                //     newGrid[i][j].setMass(prevProps[MASS],false);
                //     newGrid[i][j].setTemp(prevProps[TEMPERATURE]);
                //     newGrid[i][j].setE(prevProps[E]);
                //     newGrid[i][j].setColor(prevProps[R], prevProps[G], prevProps[B]);
                // }
            }

            double newVx, newVy;
            // velocity self-advection:
            if (i < rows) {
                // vxs of vgrid
                double vxPhysX = cellWidth*j;
                double vxPhysY = height - cellHeight*i - cellHeight/2;
                
                VelocityVector vVx = vGrid->sampleVelocityAtPoint(vxPhysX, vxPhysY);
                double prevVxX = vxPhysX - dt*vVx.getVx();
                double prevVxY = vxPhysY - dt*vVx.getVy();
                if (round(prevVxX) < -0.0f) prevVxX = -prevVxX;
                if (round(prevVxX) > width) prevVxX = 2 * width - prevVxX;
                if (round(prevVxY) < -0.0f) prevVxY = -prevVxY;
                if (round(prevVxY) > height) prevVxY = 2 * height - prevVxY;
                // prevVxX = std::max(0.0f, std::min(prevVxX, width));
                // prevVxY = std::max(0.0f, std::min(prevVxY, height));
                if ((0 < prevVxX) && (prevVxX < width) && (0 < prevVxY) && (prevVxY < height)) {
                    newVx = vGrid->sampleVelocityAtPoint(prevVxX, prevVxY).getVx();
                    // if ((j == 0) || (j == cols)) {
                    //     new_vGrid->setVx(i,j,0);
                    // } else {
                    //     new_vGrid->setVx(i,j,newVx);
                    // }
                    // if (j == 0) {
                    //     new_vGrid->setVx(i, j, -new_vGrid->getVx(i, j + 1));
                    // } else if (j == cols - 1) {
                    //     new_vGrid->setVx(i, j + 1, -new_vGrid->getVx(i, j));
                    // } else {
                    //     new_vGrid->setVx(i,j,newVx);
                    // }
                    new_vGrid->setVx(i,j,newVx);
                    // find the max velocity in x direction for dt
                    if (std::abs(newVx) > maxV) {
                        maxV = std::abs(newVx);
                    }
                }
                if (j == 0 || j == cols) {
                    // if (i == rows/2) printf("vx at j==%d: %f\n", j, new_vGrid->getVx(i,j));
                    
                    // new_vGrid->setVx(i,j,0);
                }
                // if (i == rows/2 && (j == 1 || j==cols-1)) printf("vx at j==%d: %f\n", j, new_vGrid->getVx(i,j));
            }
            
            if (j < cols) {
                // vys of vgrid
                double vyPhysX = cellWidth*j + cellWidth/2;
                double vyPhysY = height - cellHeight*i;
                VelocityVector vVy = vGrid->sampleVelocityAtPoint(vyPhysX, vyPhysY);
                double prevVyX = vyPhysX - dt*vVy.getVx();
                double prevVyY = vyPhysY - dt*vVy.getVy();

                if (round(prevVyX) < -0.0f) prevVyX = -prevVyX;
                if (round(prevVyX) > width) prevVyX = 2 * width - prevVyX;
                if (round(prevVyY) < -0.0f) prevVyY = -prevVyY;
                if (round(prevVyY) > height) prevVyY = 2 * height - prevVyY;

                // prevVyX = std::max(0.0f, std::min(prevVyX, width));
                // prevVyY = std::max(0.0f, std::min(prevVyY, height));
                if ((0 < prevVyX) && (prevVyX < width) && (0 < prevVyY) && (prevVyY < height)) {
                    newVy = vGrid->sampleVelocityAtPoint(prevVyX, prevVyY).getVy();
                    // if ((i == 0) || (i == rows)) {
                    //     new_vGrid->setVy(i,j,0);
                    // } else {
                    //     new_vGrid->setVy(i,j,newVy);
                    // }
                    new_vGrid->setVy(i,j,newVy);
                    // find the max velocity in y direction for dt
                    if (std::abs(newVy) > maxV) {
                        maxV = std::abs(newVy);
                    }
                }
                // if (i == 0 || i == rows) {
                //     new_vGrid->setVy(i,j,0);
                // }
            }
        }
    }
    grid = newGrid;
    vGrid = new_vGrid;
}

void FluidGrid::setVelocities(bool setE=false, bool noRecalc=false) {
    // so that each cell knows about its velocity
    double totMass = 0;
    double newMaxV = 0.1;
    double newMaxT = 0;
    int i,j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            std::map<char, double> xy = getxyFromij(i,j);
            VelocityVector v = vGrid->sampleVelocityAtPoint(xy['x'],xy['y']);
            FluidCell *cell = getCell(i,j);
            cell->setVelocity(v);
            if (setE && !noRecalc) {
                cell->setE(cell->gete(false)+0.5*cell->getDensity()*v.getMag()*v.getMag());
                // if (i==1 && j ==cols/2) printf("initial V: %f\n",v.getMag());
            } else if (!setE && !noRecalc) {
                double newE = cell->getE(false);
                double KE = 0.5 * cell->getDensity() * (v.getMag() * v.getMag());
                // if (i==48 && j==0) {
                //     printf("E: %f KE: %f new e: %f\n", newE, KE, newE-KE);
                // }
                if (newE >= KE) {
                    cell->sete(newE - KE);
                    // if (i==48 && j==0) {
                    //     printf("e: %f\n", cell->gete(false));
                    // }
                } else {
                    // std::cerr << "Warning: Total energy less than kinetic energy.\n";
                    cell->sete(0.0f); // Prevent negative internal energy
                }
                // cell->sete(cell->getE() - 0.5*cell->getDensity()*(v.getMag()*v.getMag()));
                // if (i==rows-1 && j ==cols/2) printf("E:%f e:%f ∆E:%f KE:%f vx:%f vy:%f\n",cell->getE(),cell->gete(), cell->getE()-cell->gete(), 0.5*cell->getDensity()*v.getMag()*v.getMag(), v.getVx(), v.getVy());
                // recalculates temperature and pressure
                cell->getTemp(true);
                cell->getPressure(true);
            }
            if (v.getMag() > newMaxV) newMaxV = v.getMag();
            if (!std::isfinite(cell->getDensity())) {
                std::cerr << "NaN detected in density at cell (" << i << ", " << j << ")\n";
                // assert(false); // Debug immediately
            }
            if (!std::isfinite(cell->getE(false))) {
                std::cerr << "NaN detected in energy at cell (" << i << ", " << j << ") velocity: " << v.getMag() << " density: " << cell->getDensity() << "\ne: " << cell->gete(false) << "\n";
                printf("vx: %f, vy: %f\n", v.getVx(), v.getVy());
                printf("setE=%d noRe=%d\n",setE,noRecalc);
                // assert(false); // Debug immediately
            }
            if (!std::isfinite(cell->gete(false))) {
                std::cerr << "NaN detected in internal energy at cell (" << i << ", " << j << ") velocity: " << v.getMag() << " density: " << cell->getDensity() << "\ne: " << cell->gete(false) << "\n";
                printf("vx: %f, vy: %f\n", v.getVx(), v.getVy());
                printf("setE=%d noRe=%d\n",setE,noRecalc);
                // assert(false); // Debug immediately
            }
            if (cell->getTemp(false) > newMaxT) {
                newMaxT = cell->getTemp(false);
            }
            totMass += cell->getMass(false);
        }
    }
    massFactor = totalMass/totMass;
    totalMass = totMass;
    maxV = newMaxV;
    maxT = newMaxT;
    // printf("tot mass: %f\n",totMass);
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
                double vxMouse = (x - prevMouseX)*SCALE_W;
                double vyMouse = -(y - prevMouseY)*SCALE_H;
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
                    gravPotentialDisplay = 0;
                    energyDisplay = 0;
                    pressureDisplay = !(pressureDisplay);
                    printf("%f\n", vGrid->getVy((rows)/2,(cols)/2));
                    break;
                case SDLK_d:
                    pressureDisplay = 0;
                    temperatureDisplay = 0;
                    gravPotentialDisplay = 0;
                    energyDisplay = 0;
                    densityDisplay = !(densityDisplay);
                    printf("dt: %f, dens: %f, massFactor: %f, totMass: %f\n", 0.1 * std::min(cellWidth, cellHeight)/maxV, getCell(rows/2, cols/2)->getDensity(), massFactor, totalMass);
                    break;
                case SDLK_t:
                    densityDisplay = 0;
                    pressureDisplay = 0;
                    gravPotentialDisplay = 0;
                    energyDisplay = 0;
                    temperatureDisplay = !(temperatureDisplay);
                    printf("maxT: %f\n", maxT);
                    break;
                case SDLK_0:
                    gravityFlag = !(gravityFlag);
                    break;
                case SDLK_u:
                    densityDisplay = 0;
                    pressureDisplay = 0;
                    temperatureDisplay = 0;
                    energyDisplay = 0;
                    gravPotentialDisplay = !(gravPotentialDisplay);
                    printf("gravPot: %f\n",getCell(rows/2,cols/2)->getGravPotential());
                    break;
                case SDLK_e:
                    densityDisplay = 0;
                    pressureDisplay = 0;
                    temperatureDisplay = 0;
                    gravPotentialDisplay = 0;
                    energyDisplay = !(energyDisplay);
                    printf("energy: %f\n",getCell(rows/2,cols/2)->gete(false));
                    break;
            }


        } 
    } else if (event.type == SDL_MOUSEBUTTONUP) {
        buttonHeld = 0;
        if (mouseVelFlag == 2) {
            SDL_MouseButtonEvent buttonEvent = event.button;
            Sint32 x = buttonEvent.x;
            Sint32 y = buttonEvent.y;
            double vxMouse = (x - prevMouseX)*SCALE_W;
            double vyMouse = -(y - prevMouseY)*SCALE_H;
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
    //     double initialMass = cell->getMass();
        
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
    FluidCell *c = getCell(rows/2,cols/2);
    if (round(c->getMass(true)) > 0) {
        // printf("smoke:%f\n", c->getMass(true));
    }
    // printf("dens:%f press:%f vyBottom:%f\n",c->getDensity(),c->getPressure(), vGrid->getVy(rows, cols/2));
    // FluidCell *c = getCell(44,0);

    dt = 0.1 * std::min(cellWidth, cellHeight)/maxV;
    if (compressible) {
        auto start = std::chrono::steady_clock::now();
        solveGravPotential(4);
        auto grav = std::chrono::steady_clock::now();
        compressibleMomentumUpdate();
        auto momentum = std::chrono::steady_clock::now();
        setVelocities(true, false);
        auto vel1 = std::chrono::steady_clock::now();
        dt = 0.1 * std::min(cellWidth, cellHeight)/maxV;
        energyUpdate();
        auto energy = std::chrono::steady_clock::now();
        setVelocities(false, true);
        dt = 0.1 * std::min(cellWidth, cellHeight)/maxV;
        auto vel2 = std::chrono::steady_clock::now();
        FFSLadvect();
        auto advect = std::chrono::steady_clock::now();
        setVelocities(false, false);
        // std::cout << "grav: ";
        // std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(grav - start).count()/1000.0f << "\n";
        // std::cout << "momentum: ";
        // std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(momentum - grav).count()/1000.0f << "\n";
        // std::cout << "energy: ";
        // std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(energy - vel1).count()/1000.0f << "\n";
        // std::cout << "advect: ";
        // std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(advect - vel2).count()/1000.0f << "\n";
        
        
        
        dt = 0.1 * std::min(cellWidth, cellHeight)/maxV;
    } else {
        projection(10);
        advect();
        setVelocities(false, true);
        // std::cout << "Left boundary velocity: " << vGrid->getVx(rows/2, 0) << "\n";
        // std::cout << "Right boundary velocity: " << vGrid->getVx(rows/2, cols) << "\n";
        // setVelocities(false, true);
    }
    // printf("projection — dens:%f press:%f\n  vxL:%f vxR:%f vyB:%f vyT:%f\n", c->getDensity(),c->getPressure(),vGrid->getVx(rows-1,0),vGrid->getVx(rows-1,1),vGrid->getVy(rows,0),vGrid->getVy(rows-1,0));
    // printf("projection — dens:%f press:%f\n  vxL:%f vxR:%f vyB:%f vyT:%f\n", c->getDensity(),c->getPressure(),vGrid->getVx(44,0),vGrid->getVx(44,1),vGrid->getVy(45,0),vGrid->getVy(44,0));
    
    // printf("advection — dens:%f press:%f\n  vxL:%f vxR:%f vyB:%f vyT:%f\n", c->getDensity(),c->getPressure(),vGrid->getVx(rows-1,0),vGrid->getVx(rows-1,1),vGrid->getVy(rows,0),vGrid->getVy(rows-1,0));
    // printf("advection — dens:%f press:%f\n  vxL:%f vxR:%f vyB:%f vyT:%f\n", c->getDensity(),c->getPressure(),vGrid->getVx(44,0),vGrid->getVx(44,1),vGrid->getVy(45,0),vGrid->getVy(44,0));
    
    // printf("%f\n", std::min(cellWidth, cellHeight)/maxV);
    
    // dt = 0.001;
    // if (maxV > 0) dt = std::min(0.016,0.4 * std::min(cellWidth, cellHeight)/maxV);
    // else dt = 0.016;
    // viscosity = 0.01*(cellWidth*cellHeight)/dt;
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

void FluidGrid::force(int i, int j, double fx, double fy) {
    FluidCell *cell = getCell(i,j);
    double mass = cell->getMass(false);
    if (round(mass*100) == 0) return; // nothing to force if mass is zero
    double vxL = vGrid->getVx(i,j);
    double vxR = vGrid->getVx(i,j+1);
    double vyU = vGrid->getVy(i,j);
    double vyB = vGrid->getVy(i+1,j);
    double ax = fx / mass;
    // if (i==10 && (j == 5 || j==cols-6)) printf("ax: %f\n",ax);
    double ay = fy / mass;
    if (fx < 0) {
        vGrid->setVx(i,j,vxL+ax*dt);
    } else {
        vGrid->setVx(i,j+1,vxR+ax*dt);
    }

    // if (i==10 && (j == 5 || j==cols-6)) printf("vxL: %f, vxR: %f\n",vGrid->getVx(i,j),vGrid->getVx(i,j+1));
    // vGrid->setVx(i,j,vxL+ax*dt);
    // vGrid->setVx(i,j+1,vxR+ax*dt);
    if (j == 0) {
        vGrid->setVx(i,j,0); 
    } else if (j == cols-1) {
        vGrid->setVx(i,j+1,0);
    }
    
    // vGrid->setVx(i,j+1,vxR+ax*dt);
    if (fy < 0) {
        vGrid->setVy(i+1,j,vyB+ay*dt);
    } else {
        vGrid->setVy(i,j,vyU+ay*dt);
    }
    // vGrid->setVy(i+1,j,vyB+ay*dt);
    // vGrid->setVy(i,j,vyU+ay*dt);
    if (i == 0) {
        vGrid->setVy(i,j,0); 
    } else if (i == rows-1) {
        vGrid->setVy(i+1,j,0);
    }
    // vGrid->setVy(i+1,j,vyB+ay*dt);
}

class Simulator {
    public:
        Simulator(double width, double height, int c, int r, bool comp, bool display=true, std::string fn=""): width(width), height(height), rows(r), cols(c), compressible(comp), display(display), filename(fn) {
            dt = 0.016;
            // in meters per pixel
            SCALE_H = height / consts::GRID_HEIGHT;
            SCALE_W = width / consts::GRID_WIDTH;
            cells = rows * cols;
            grid = (FluidGrid *) malloc(sizeof(FluidGrid));
            grid[0] = FluidGrid(width, height, rows, cols, dt, compressible);
            output.open(filename);
            if (display) {
                SDL_Init(SDL_INIT_VIDEO);       // Initializing SDL as Video
                SDL_CreateWindowAndRenderer(consts::GRID_WIDTH, consts::GRID_HEIGHT, 0, &window, &renderer);
                SDL_SetRenderDrawColor(renderer, 0, 0, 0, 0);      // setting draw color
                SDL_RenderClear(renderer);
            }
            
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
                // max/min boundaries don't update before the full render
            double thisMaxPressure = 0;
            double thisMinPressure = 0;
            double thisMaxMass = 0;
            double thisMaxDensity = 0;
            double thisMaxTemperature = 0;
            double thisMinTemperature = 0;
            double thisMinGP = 0;
            double thisMaxe = 0;
            int imaxT, jmaxT;
            for (i = 0; i < rows; i++) {
                for (j = 0; j < cols; j++) {
                    
                    double mass = grid->getCell(i,j)->getMass();
                    // if (mass > maxMass) {
                    //     maxMass = mass;
                    // }
                    if (mass > thisMaxMass) thisMaxMass = mass;
                    double pressure = grid->getCell(i,j)->getPressure(false);
                    double density = grid->getCell(i,j)->getDensity();
                    double temperature = grid->getCell(i,j)->getTemp(false);
                    double gravPotential = grid->getCell(i,j)->getGravPotential();
                    double e = grid->getCell(i,j)->gete(false);

                    if (density > thisMaxDensity) thisMaxDensity = density;
                    // if (thisMaxDensity > maxDensity) maxDensity = thisMaxDensity;
                    if (pressure > thisMaxPressure) {
                        thisMaxPressure = pressure;
                    } else if (pressure < thisMinPressure) {
                        thisMinPressure = pressure;
                    }
                    if (temperature > thisMaxTemperature) {
                        thisMaxTemperature = temperature;
                        imaxT = i;
                        jmaxT = j;
                    }
                    if (gravPotential < thisMinGP) thisMinGP = gravPotential;
                    if (e > thisMaxe) thisMaxe = e;
                    // if (pressure > maxPressure) {
                    //     maxPressure = pressure;
                    // } else if (pressure < minPressure) {
                    //     minPressure = pressure;
                    // }
                    SDL_Rect rect{j*consts::GRID_WIDTH/cols,i*consts::GRID_HEIGHT/rows,(j+1)*consts::GRID_WIDTH/cols,(i+1)*consts::GRID_HEIGHT/rows};
                    double maxBit = 255.0;
                    double zero = 0;
                    if (grid->pressureDisplay) {
                        // double scaledP_R = 255 * (pressure-minPressure)/(maxPressure-minPressure);
                        double scaledP_B = std::min(maxBit,255 * (pressure)/(minPressure)); // outwards pressure
                        double scaledP_R = std::min(maxBit,255 * (pressure)/(maxPressure)); // inwards pressure
                        // if (scaledP_R > 255) printf("%d %d\n", i,j);
                        if (scaledP_R < 0) scaledP_R = 0;
                        if (scaledP_B < 0) scaledP_B = 0;

                        SDL_SetRenderDrawColor(renderer, scaledP_R, 0, scaledP_B, 255);
                        // SDL_SetRenderDrawColor(renderer, mass, mass, 0, 255);
                        // SDL_RenderFillRect(renderer, &rect);
                    } else if (grid->densityDisplay) {
                        // FluidCell *cell = grid->getCell(i, j);
                        
                        double scaled_dens = density/maxDensity * 255;
                        if (scaled_dens > 255) scaled_dens = 255;
                        SDL_SetRenderDrawColor(renderer, scaled_dens, scaled_dens*(1-grid->getCell(i,j)->getY()), 0, 255);
                        // SDL_RenderFillRect(renderer, &rect);
                    } else if (grid->temperatureDisplay) {
                        double scaled_temp = std::min(maxBit,temperature/maxTemperature * 255);
                        if (grid->getCell(i,j)->getTemp() < 0) {
                            // std::cerr << "Temperature (" << i << ", " << j << ")\n";
                            // assert(false); // Debug immediately
                        }
                        SDL_SetRenderDrawColor(renderer, 0, scaled_temp, 0, 255);
                    } else if (grid->gravPotentialDisplay) {
                        double scaled_pot = std::max(zero,std::min(maxBit,gravPotential/minGP * 255));
                        // if (i == rows/2 && j == cols/2) printf("scaledGP:%f\n",gravPotential/minGP);
                        SDL_SetRenderDrawColor(renderer, 255, scaled_pot, 0, 255);

                    } else if (grid->energyDisplay) {
                        double scaled_e = std::max(zero,std::min(maxBit,e/maxe * 255));
                        SDL_SetRenderDrawColor(renderer, scaled_e, 0, 0, 255);
                    } else {
                        // double mass = grid->getActive()[i]->getMass();
                        // FluidCell *cell = grid->getActive()[i];
                        std::vector<double> color = grid->getCell(i,j)->getColor();
                        
                        if ((i == rows/2) && (j == cols/2)) {
                            // std::cout << color[0] << std::endl;
                            // printf("%f %f %f\n", pressure, mass, grid->getCell(i,j)->getMass(false));
                        }
                        // SDL_Rect rect{cell->getCol()*consts::GRID_WIDTH/cols,cell->getRow()*consts::GRID_HEIGHT/rows,(cell->getCol()+1)*consts::GRID_WIDTH/cols,(cell->getRow()+1)*consts::GRID_HEIGHT/rows};
                        // if (mass > 255) {
                        //     mass = 255;
                        // }
                        // if (!cell->isActive()) SDL_SetRenderDrawColor(renderer, 0, 255, 0, 255);
                        uint red = std::min((double)255, 255*color[0]*std::sqrt(mass/maxMass));
                        uint green = std::min((double)255, 255*color[1]*std::sqrt(mass/maxMass));
                        uint blue = std::min((double)255, 255*color[2]*std::sqrt(mass/maxMass));
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
            maxMass = thisMaxMass;
            maxTemperature = thisMaxTemperature;
            minTemperature = thisMinTemperature;
            minGP = thisMinGP; // grav potential
            maxe = thisMaxe;
            // printf("maxtemp:%f\n",maxTemperature);
            // thisMaxDensity = 1;
            // thisMinPressure = -1;
            // thisMaxPressure = 1;
            FluidCell *cellMaxT = grid->getCell(imaxT, jmaxT);
            // printf("maxT loc: (%d,%d) temp: %f, pressure: %f\n  density: %f, e: %f\n",imaxT,jmaxT, cellMaxT->getTemp(false), cellMaxT->getPressure(false), cellMaxT->getDensity(), cellMaxT->gete(false));
            
        }
        
        void step(SDL_Event event, bool save=false) {
            grid->update(event);
            dt = grid->getdt();
            t += grid->getdt();
            if (save) {
                saveSim();
            } 
            if (display) {
                drawCells();
            }
            // drawCells();
            // SDL_Delay(grid->getdt()*1000);  // setting some Delay
            if (display) {
                SDL_Delay(1);
            }
        }

        void saveSim() {
            int i,j;
            for (i = 0; i < rows; i++) {
                for (j = 0; j < cols; j++) {   
                    FluidCell cell = *grid->getCell(i,j);
                    double mass = cell.getMass(false);
                    double pressure = cell.getPressure(false);
                    double density = cell.getDensity();
                    double temperature = cell.getTemp(false);
                    double gravPotential = cell.getGravPotential();
                    double e = cell.gete(false);
                    output << "//" << density << "," << cell.getX() << "," << cell.getY() << "," << cell.getZ() << "," << temperature << "," << pressure << "," << e << "," << gravPotential;
                }
                output << "//\n";
            }
            output << "//t: " << t;
            output << "//\n";
        }

        void freeSim() {
            grid->freeGrid();
            free(grid);
            output.close();
        }

        double getT() {
            // time passed
            return t;
        }

        double SCALE_H, SCALE_W;
        double width, height;
        FluidGrid *grid;
    private:
        bool compressible, display;
        double dt;
        int rows, cols, cells;
        SDL_Renderer *renderer = NULL;
        SDL_Window *window = NULL;
        SDL_Surface *screenSurface;
        std::ofstream output = std::ofstream();
        std::string filename = "";
        double t = 0;
        double maxMass = 255;
        double maxPressure = 1;
        double minPressure = 0;
        double maxDensity = 1;
        double maxTemperature = 1;
        double minTemperature = 0;
        double minGP = 0;
        double maxe = 0;
};


int main(int argv, char **argc) {
    if (argv > 8) std::srand((unsigned) atoi(argc[8]));
    else std::srand((unsigned) std::time(NULL));
    double width, height;
    int c, r;
    bool comp, save, display;
    if (argv < 8) {
        width = atol(argc[1]);
        height = atol(argc[2]);
        c = atoi(argc[3]);
        r = atoi(argc[4]);
        comp = false;
        display = (bool)(!!atoi(argc[5]));
        save = (bool)(!!atoi(argc[6]));
        // sim.drawCells();
    } else {
        width = atol(argc[2]);
        height = atol(argc[3]);
        c = atoi(argc[4]);
        r = atoi(argc[5]);
        if (atoi(argc[1])) comp = true; else comp = false;
        display = (bool)(!!atoi(argc[6]));
        save = (bool)(!!atoi(argc[7]));
        // comp = atoi(argc[1]);
        // Simulator sim(atoi(argc[2]), atoi(argc[3]), atoi(argc[4]), atoi(argc[5]), (bool)atoi(argc[1]));
        
    }
    int maxT = 2000000;
    // std::string filename = "/Volumes/Give Me Space/outputGrid_3.txt";
    std::string filename = "outputGrid_2000000.txt";
    // Simulator sim(atoi(argc[1]), atoi(argc[2]), atoi(argc[3]), atoi(argc[4]), false);
    Simulator sim(width, height, c, r, comp, display, filename);
    // sim.drawCells();
    SDL_Event event;
    // 
    // sim.step(event);
    // sim.drawGridLines();
    uint pause = 0;
    uint now = (unsigned) std::time(NULL);
    const std::chrono::time_point<std::chrono::steady_clock> start =
        std::chrono::steady_clock::now();
    
    double progress = 0.0;
    
    while(!(event.type == SDL_QUIT) and sim.getT() < maxT){
        // if (display && !save) maxT = sim.getT()+100;
        if (!pause) {
            sim.step(event, save);
            // printf("%f\n", sim.getT());   
        }
        // sim.drawGridLines();
        SDL_PollEvent(&event);  // Catching the poll event.
        // if (event.key.type == SDL_KEYDOWN) {
        //     // printf("key clicked\n");
        //     if (event.key.keysym.sym == SDLK_SPACE) pause = !pause;
        // }
        int barWidth = 70;

        std::cout << "[";
        int pos = barWidth * progress;
        for (int i = 0; i < barWidth; ++i) {
            if (i < pos) std::cout << "=";
            else if (i == pos) std::cout << ">";
            else std::cout << " ";
        }
        std::cout << "] " << int(progress * 100.0) << " %\r";
        std::cout.flush();
        progress = sim.getT()/maxT; // for demonstration only
    }
    std::cout << "\n";
    const auto end = std::chrono::steady_clock::now();
    std::cout << "real time elapsed: ";
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()/1000.0f << "\nsim time elapsed: " << sim.getT() << "\n";
    // std::cout << "real time elapsed: " << (uint) (end - start) << "\nsim time elapsed: " << sim.getT() << "\n";
    // std::cout
    //     << "Slow calculations took "
    //     << std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    // printf("real time elapsed: %d\nsim time elapsed: %f\n", (end-start, sim.getT());
    sim.freeSim();

    // #pragma omp parallel
    // {
    //     std::cout << "Hello from thread " << omp_get_thread_num() << std::endl;
    // }

    return 1;
}

