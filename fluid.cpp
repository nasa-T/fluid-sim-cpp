#include <cmath>
#include <iostream>
#include "fluid.h"
#include <SDL2/SDL.h>
#include <cstdlib>
#include <map>
#include <string>
#include <vector>
#include <vector>
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
        FluidCell(float mass, float width, float height, float temp): mass(mass), width(width), height(height), temperature(temp) {
            // size is the physical area
            size = width * height;
            density = mass/size;
            pressure = density/mass * consts::kb * temperature;
            // velocity = VelocityVector();
            // velocity.setVx(1);
            // velocity.setVx((std::rand() % (int)(2*width)) - width/2);
            // velocity.setVy((std::rand() % (int)(2*height)) - height/2);
            // vBounds = VelocityBox(0,0,0,0);
            // vBounds.setVelocity(velocity);
            neighbors = Neighbors();
        }

        float getMass() {
            return mass;
        }
        float getSize() {
            return size;
        }
        float getDensity() {
            return density;
        }
        float getTemp() {
            return temperature;
        }
        float getPressure() {
            return pressure;
        }
        void setMass(float m) {
            if (m >= 0) mass = m;
        }
        void addMass(float m) {
            setMass(mass + m);
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
        float mass, size, density, temperature, pressure, width, height;
        int row, col;
        VelocityVector velocity;
        // VelocityBox vBounds;
        Neighbors neighbors;
};

class FluidGrid {
    public:
        FluidGrid(float width, float height, int r, int c): width(width), height(height), rows(r), cols(c) {
            int i, j;
            cells = rows * cols;
            grid = (FluidCell **) malloc(sizeof(FluidCell*)*r);
            newGrid = (FluidCell **) malloc(sizeof(FluidCell*)*r);
            vGrid = (VelocityGrid *) malloc(sizeof(VelocityGrid));
            vGrid[0] = VelocityGrid(width, height, r, c);
            new_vGrid = (VelocityGrid *) malloc(sizeof(VelocityGrid));
            new_vGrid[0] = VelocityGrid(width, height, r, c);
            cellWidth = width / c;
            cellHeight = height / r;
            // in meters per pixel
            SCALE_H = height / consts::GRID_HEIGHT;
            SCALE_W = width / consts::GRID_WIDTH;
            maxV = 0;
            for (i = 0; i < rows; i++) {
                grid[i] = (FluidCell *) malloc(sizeof(FluidCell)*c);
                newGrid[i] = (FluidCell *) malloc(sizeof(FluidCell)*c);
                for (j = 0; j < cols; j++) {
                    //random mass for now
                    // float mass = std::rand() % 256;
                    grid[i][j] = FluidCell(0, cellWidth, cellHeight, 100);
                    if (i == rows/2+1 && j == cols/2) {
                        // grid[i][j] = FluidCell(100, cellWidth, cellHeight, 100);
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
                    if (grid[i][j].getVelocity().getMag() > maxV) {
                        maxV = grid[i][j].getVelocity().getMag();
                    }
                    grid[i][j].setLoc(i, j);
                    newGrid[i][j] = grid[i][j];
                    new_vGrid->setVx(i,j,vGrid->getVx(i,j));
                    new_vGrid->setVy(i,j,vGrid->getVy(i,j));
                }
            }
            dt = 0.25;
            // 0.7/maxV;
            nActive = activeCells.size();
        }

        FluidCell *getCell(int i, int j) {
            return &grid[i][j];
        }

        void freeGrid() {
            int i, j;
            for (i = 0; i < rows; i++) {
                free(grid[i]);
            }
            free(grid);
        }

        std::map<uint, float> sampleCellAtPoint(float x, float y) {
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
                p10[MASS] = cell->getMass();
                p10[TEMPERATURE] = cell->getTemp();
                p10[PRESSURE] = cell->getPressure();
            } else {
                p10[MASS] = 0;
                p10[TEMPERATURE] = 0;
                p10[PRESSURE] = 0;
            }
            if ((j > -1) && (i > 0)) {
                // p00 = vyArray[i-1][j];
                FluidCell *cell = getCell(i-1,j);
                p00[MASS] = cell->getMass();
                p00[TEMPERATURE] = cell->getTemp();
                p00[PRESSURE] = cell->getPressure();
            } else {
                p00[MASS] = 0;
                p00[TEMPERATURE] = 0;
                p00[PRESSURE] = 0;
            }
            if ((j < cols-1) && (i > 0)) {
                // p01 = vyArray[i-1][j+1];
                FluidCell *cell = getCell(i-1,j+1);
                p01[MASS] = cell->getMass();
                p01[TEMPERATURE] = cell->getTemp();
                p01[PRESSURE] = cell->getPressure();
            } else {
                p01[MASS] = 0;
                p01[TEMPERATURE] = 0;
                p01[PRESSURE] = 0;
            }
            if ((j < cols-1) && (i < rows)) {
                // p11 = vyArray[i][j+1];
                FluidCell *cell = getCell(i,j+1);
                p11[MASS] = cell->getMass();
                p11[TEMPERATURE] = cell->getTemp();
                p11[PRESSURE] = cell->getPressure();
            } else {
                p11[MASS] = 0;
                p11[TEMPERATURE] = 0;
                p11[PRESSURE] = 0;
            }
            // printf("%f %f %f %f\n", p10[MASS], p00[MASS], p01[MASS], p11[MASS]);
            std::map<uint, float> props;
            props[MASS] = (1-X)*(1-Y)*p10[MASS] + X*(1-Y)*p11[MASS] + (1-X)*Y*p00[MASS] + X*Y*p01[MASS];
            // printf("%f\n", props[MASS]);
            props[TEMPERATURE] = (1-X)*(1-Y)*p10[TEMPERATURE] + X*(1-Y)*p11[TEMPERATURE] + (1-X)*Y*p00[TEMPERATURE] + X*Y*p01[TEMPERATURE];
            props[PRESSURE] = (1-X)*(1-Y)*p10[PRESSURE] + X*(1-Y)*p11[PRESSURE] + (1-X)*Y*p00[PRESSURE] + X*Y*p01[PRESSURE];

            return props;
        }

        void projection(int iters) {
            int i, j;
            for (int n = 0; n < iters; n++) {
                for (i = 0; i < rows; i++) {
                    for (j = 0; j < cols; j++) {
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
                        float d = -left + right + top - bottom;

                        // if ((i == rows/2+1) && (j == cols/2)) {
                        //     printf("%f %f %f %f\n", left, right, top, bottom);
                        //     printf("%d %d %d %d\n", sLeft, sRight, sTop, sBottom);
                        //     printf("%f\n", d);
                        // }
                        vGrid->setVx(i, j, left+d*sLeft/s);
                        vGrid->setVx(i, j+1, right-d*sRight/s);
                        vGrid->setVy(i, j, top-d*sTop/s);
                        vGrid->setVy(i+1, j, bottom+d*sBottom/s);
                        // } else if ()
                    }
                }
            }
        }

        void advect() {
            int i, j;
            for (i = 0; i < rows+1; i++) {
                for (j = 0; j < cols+1; j++) {
                    if ((i < rows) && (j < cols)) {
                        FluidCell cell = *getCell(i,j);
                        if ((i == rows/2) && (j == cols/2)) {
                            // printf("%f\n", cell.getMass());
                            // printf("%f %f\n", prevY, physY);
                        }
                        // VelocityVector v = cell.getVelocity();
                        float physX = cellWidth*j+cellWidth/2;
                        float physY = height - (cellHeight*i+cellHeight/2); // I want y pointed up
                        // printf("%f %f\n", physX, physY);
                        VelocityVector vCell = vGrid->sampleVelocityAtPoint(physX, physY);
                        if ((i == rows/2) && (j == cols/2)) {
                            // printf("%f\n", cell.getMass());
                            // printf("%f %f\n", vCell.getVx(), vCell.getVy());
                        }
                        // printf("%f %f\n", vCell.getVx(), vCell.getVy());
                        float prevX = physX - dt*vCell.getVx();
                        float prevY = physY - dt*vCell.getVy();
                        std::map<uint, float> prevProps;
                        if (round(vCell.getVx()+vCell.getVy()) == 0) {
                            prevProps[MASS] = cell.getMass();
                            prevProps[TEMPERATURE] = cell.getTemp();
                            prevProps[PRESSURE] = cell.getPressure();
                        } else {
                            // printf("%f %f\n", prevX, prevY);
                            prevProps = sampleCellAtPoint(prevX, prevY);
                        }
                        // printf("%f\n", prevProps[MASS]);
                        if ((i == rows/2) && (j == cols/2)) {
                            // printf("%f\n", prevProps[MASS]);
                            // printf("%f %f\n", prevY, physY);
                        }
                        // printf("%f ", prevX);
                        // std::cout << physX << " ";
                        if ((prevX < 0) || (prevY < 0) || (prevX > width) || (prevY > height)) {
                            newGrid[i][j].setMass(0);

                        } else {
                            // int prev_i = prevY / cellHeight;
                            // int prev_j = prevX / cellWidth;
                            if ((i == rows/2) && (j == cols/2 + 1)) {
                                // printf("%f; %f, %f: %d, %d\n", cell.getMass(), prevX, prevY, prev_j, prev_i);
                                // printf("%f\n", v.getVx());
                            }
                            newGrid[i][j].setMass(prevProps[MASS]);
                            // newGrid[i][j].setMass(grid[prev_i][prev_j].getMass());
                            //  = grid[prev_i][prev_j];
                            // newGrid[i][j].setLoc(i,j);
                        }
                    }
                    if (i < rows) {
                        // vxs of vgrid
                        float vxPhysX = cellWidth*j;
                        float vxPhysY = height - cellHeight*i - cellHeight/2;
                        VelocityVector vVx = vGrid->sampleVelocityAtPoint(vxPhysX, vxPhysY);
                        float prevVxX = vxPhysX - dt*vVx.getVx();
                        float prevVxY = vxPhysY - dt*vVx.getVy();
                        new_vGrid->setVx(i,j,vGrid->sampleVelocityAtPoint(prevVxX, prevVxY).getVx());
                    }
                    if (j < cols) {
                        // vys of vgrid
                        float vyPhysX = cellWidth*j + cellWidth/2;
                        float vyPhysY = height - cellHeight*i;
                        VelocityVector vVy = vGrid->sampleVelocityAtPoint(vyPhysX, vyPhysY);
                        float prevVyX = vyPhysX - dt*vVy.getVx();
                        float prevVyY = vyPhysY - dt*vVy.getVy();
                        new_vGrid->setVy(i,j,vGrid->sampleVelocityAtPoint(prevVyX, prevVyY).getVy());
                        
                    }
                }
            }
            grid = newGrid;
            vGrid = new_vGrid;
            // printf("%f\n", getCell(rows/2,cols/2)->getVelocity().getVx());
        }

        void update(SDL_Event event) {
            int i, j;
            if ((event.type == SDL_MOUSEBUTTONDOWN) || (event.type != SDL_MOUSEBUTTONUP)) {
                SDL_MouseButtonEvent buttonEvent = event.button;
                if (buttonEvent.button == SDL_BUTTON_RIGHT) {
                    mouseVelFlag = !mouseVelFlag;
                } else if (buttonEvent.button == SDL_BUTTON_LEFT) {
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
                        // printf("%f %f\n", vxMouse, vyMouse);
                        vGrid->setVx(i, j+1, vxMouse);
                        vGrid->setVx(i, j, vxMouse);
                        vGrid->setVy(i, j, vyMouse);
                        vGrid->setVy(i+1, j, vyMouse);
                        // if (vxMouse > 0) {
                        //     vGrid->setVx(i, j+1, vxMouse);
                        // } else {
                        //     vGrid->setVx(i, j, vxMouse);
                        // }
                        // if (vyMouse > 0) {
                        //     vGrid->setVy(i, j, vyMouse);
                        // } else {
                        //     vGrid->setVy(i+1, j, vyMouse);
                        // }
                    } else if (mouseVelFlag == 0) {
                        FluidCell *clickedCell = getCell(i, j);
                        clickedCell->addMass(100);
                    }
                }
                // FluidCell *clickedCell = getCell(i, j);
                // printf("%f\n", clickedCell->getMass());
                // if (!clickedCell->isActive()) activeCells.push_back(clickedCell);
                // clickedCell->addMass(100);
            } else if (event.type == SDL_MOUSEBUTTONUP) {
                prevMouseX = 0;
                prevMouseY = 0;
            }
            // float totalMass = 0;
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
            
            projection(4);
            advect();
            // print("advected");
            // printf("%f\n", maxV);
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

        std::vector<FluidCell*> getActive() {
            return activeCells;
        }
        unsigned int nActive;
        VelocityGrid *vGrid;
        // if 0, mouse left click adds mass; if 1, mouse left click dragging changes velocities
        uint mouseVelFlag = 0;
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
};

class Simulator {
    public:
        Simulator(float width, float height, int c, int r): width(width), height(height), rows(r), cols(c) {
            // in meters per pixel
            SCALE_H = height / consts::GRID_HEIGHT;
            SCALE_W = width / consts::GRID_WIDTH;
            cells = rows * cols;
            grid = (FluidGrid *) malloc(sizeof(FluidGrid));
            grid[0] = FluidGrid(width, height, rows, cols);
            
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
                    // float mass = grid->getActive()[i]->getMass();
                    // FluidCell *cell = grid->getActive()[i];
                    
                    SDL_Rect rect{j*consts::GRID_WIDTH/cols,i*consts::GRID_HEIGHT/rows,(j+1)*consts::GRID_WIDTH/cols,(i+1)*consts::GRID_HEIGHT/rows};
                    // SDL_Rect rect{cell->getCol()*consts::GRID_WIDTH/cols,cell->getRow()*consts::GRID_HEIGHT/rows,(cell->getCol()+1)*consts::GRID_WIDTH/cols,(cell->getRow()+1)*consts::GRID_HEIGHT/rows};
                    if (mass > 255) {
                        mass = 255;
                    }
                    // if (!cell->isActive()) SDL_SetRenderDrawColor(renderer, 0, 255, 0, 255);
                    SDL_SetRenderDrawColor(renderer, mass, mass, 0, 255);
                    SDL_RenderFillRect(renderer, &rect);
            // }
                }
            }
            SDL_RenderPresent(renderer);
        }
        
        void step(SDL_Event event) {
            grid->update(event);
            drawCells();

            // int i, j;
            // for (i = 0; i < rows; i++) {
            //     for (j = 0; j < cols; j++) {
            //         FluidCell *cell = grid->getCell(i, j);
            //         cell->
            //     }
            // }
        }

        void freeSim() {
            grid->freeGrid();
            free(grid);
        }
        float SCALE_H, SCALE_W;
        float width, height;
        FluidGrid *grid;
    private:
        
        int rows, cols, cells;
        SDL_Renderer *renderer = NULL;
        SDL_Window *window = NULL;
        SDL_Surface *screenSurface;
        
        
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
        SDL_Delay(20);  // setting some Delay
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

