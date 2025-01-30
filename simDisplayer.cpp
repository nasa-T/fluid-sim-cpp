#include <cmath>
#include <iostream>
#include <fstream>

#include <GL/glew.h>

#include <GLFW/glfw3.h>
GLFWwindow* window;

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/norm.hpp>
using namespace glm;

#include <SDL2/SDL.h>
#include <vector>

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

std::string strip(std::string s, std::string ends) {
    std::string token;

    if (s.find(ends) == 0) {
        
        token = s.substr(ends.length(), s.length()-ends.length());
        
    } else {
        token = s;
    }
    if (token.substr(token.length()-ends.length(), ends.length()) == ends) {

        token = token.substr(0, token.length()-ends.length());
    }

    return token;
}

int main(int argv, char **argc) {
    // std::ifstream inputFile("/Volumes/Give Me Space/outputGrid1000000.txt");

    int width = 800;
    int height = 800;

    SDL_Renderer *renderer = NULL;
    SDL_Window *window = NULL;
    SDL_Surface *screenSurface;

    SDL_Init(SDL_INIT_VIDEO);       // Initializing SDL as Video
    SDL_CreateWindowAndRenderer(width, height, 0, &window, &renderer);
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 0);      // setting draw color
    SDL_RenderClear(renderer);

    std::ifstream inputFile("/Volumes/Give Me Space/outputGrid10000000.txt");
    std::vector<std::string> lines;
    std::string line;
    int rows, cols;
    int i,j = 0;
    if (inputFile.is_open()) {
        // std::string line;
        
        while (std::getline(inputFile, line)) {
            line = strip(line,"//");
            if (line.substr(0,1) == "t") break;
            lines.push_back(line);
            // if (line) lines.push_back(line);
        }
        rows = lines.size();
        cols = split(lines[0], "//").size();
        inputFile.clear();
        inputFile.seekg(0);
        
    
        float maxTemperature = 100;
        float thisMaxTemperature = 0;
        float maxDensity = 1;
        float thisMaxDensity = 0;
        float maxPressure = 1;
        float thisMaxPressure = 0;
        float maxe = 1;
        float thisMaxe = 0;
        float minGP = 1;
        float thisMinGP = 0;

        uint pressureDisplay = 0;
        uint temperatureDisplay = 0;
        uint gravPotentialDisplay = 0;
        uint energyDisplay = 0;

        SDL_Event event;

        i = 0;
        while (!(event.type == SDL_QUIT) and std::getline(inputFile, line)) {
            line = strip(line, "//");
            std::vector<std::string> row;
            if (line.substr(0,1) == "t") {
                SDL_PollEvent(&event);
                if (event.key.type == SDL_KEYDOWN) {
                    switch (event.key.keysym.sym) {
                        case SDLK_p:
                            temperatureDisplay = 0;
                            gravPotentialDisplay = 0;
                            energyDisplay = 0;
                            pressureDisplay = !(pressureDisplay);
                            printf("P: %f\n", maxPressure);
                            break;
                        case SDLK_d:
                            pressureDisplay = 0;
                            temperatureDisplay = 0;
                            gravPotentialDisplay = 0;
                            energyDisplay = 0;
                            break;
                        case SDLK_t:
                            pressureDisplay = 0;
                            gravPotentialDisplay = 0;
                            energyDisplay = 0;
                            temperatureDisplay = !(temperatureDisplay);
                            printf("T: %f\n", maxTemperature);
                            break;
                        case SDLK_u:
                            pressureDisplay = 0;
                            temperatureDisplay = 0;
                            energyDisplay = 0;
                            gravPotentialDisplay = !(gravPotentialDisplay);
                            printf("G: %f\n", minGP);
                            break;
                        case SDLK_e:
                            pressureDisplay = 0;
                            temperatureDisplay = 0;
                            gravPotentialDisplay = 0;
                            energyDisplay = !(energyDisplay);
                            printf("e: %f\n", maxe);
                            break;
                    }
                }
                i = 0;
                maxTemperature = thisMaxTemperature;
                thisMaxTemperature = 0;
                maxDensity = thisMaxDensity;
                thisMaxDensity = 0;
                maxe = thisMaxe;
                thisMaxe = 0;
                maxPressure = thisMaxPressure;
                thisMaxPressure = 0;
                minGP = thisMinGP;
                thisMinGP = 0;
                // SDL_Delay(16);
                SDL_RenderPresent(renderer);
                continue;
            }
            row = split(line, "//");
            // for (j = 0; j < cols; j++) {
            //     SDL_Rect rect{j*width/cols,i*height/rows,(j+1)*width/cols,(i+1)*height/rows};
            //     std::vector<std::string> cell = split(row[j], ",");
            //     float temperature = std::stof(cell[4]);
            //     if (thisMaxTemperature < temperature) {
            //         thisMaxTemperature = temperature;
            //     }
            //     float scaled_temp = std::min(255.0f,temperature/maxTemperature * 255);

            //     SDL_SetRenderDrawColor(renderer, 0, scaled_temp, 0, 255);
            //     SDL_RenderFillRect(renderer, &rect);
            // }
            for (j = 0; j < cols; j++) {
                std::vector<std::string> cell = split(row[j], ",");
                float pressure = (float)std::stod(cell[5]);
                float density = (float)std::stod(cell[0]);
                float X = (float)std::stod(cell[1]);
                float Y = (float)std::stod(cell[2]);
                float Z = (float)std::stod(cell[3]);
                float temperature = (float)std::stod(cell[4]);
                float gravPotential = (float)std::stod(cell[7]);
                float e = (float)std::stod(cell[6]);


                if (density > thisMaxDensity) thisMaxDensity = density;
                // if (thisMaxDensity > maxDensity) maxDensity = thisMaxDensity;
                if (pressure > thisMaxPressure) {
                    thisMaxPressure = pressure;
                }
                if (temperature > thisMaxTemperature) {
                    thisMaxTemperature = temperature;
                }
                if (gravPotential < thisMinGP) thisMinGP = gravPotential;
                if (e > thisMaxe) thisMaxe = e;
                SDL_Rect rect{j*width/cols,i*height/rows,(j+1)*width/cols,(i+1)*height/rows};
                if (pressureDisplay) {
                    float scaledP_R = std::max(0.0f,std::min(255.0f,255 * (pressure)/(maxPressure))); // inwards pressure
                    SDL_SetRenderDrawColor(renderer, scaledP_R, 0, 255, 255);
                } else if (temperatureDisplay) {
                    float scaled_temp = std::max(0.0f,std::min(255.0f,temperature/maxTemperature * 255));
                    SDL_SetRenderDrawColor(renderer, 0, scaled_temp, 0, 255);
                } else if (gravPotentialDisplay) {
                    float scaled_pot = std::max(0.0f,std::min(255.0f,gravPotential/minGP * 255));
                    SDL_SetRenderDrawColor(renderer, 255, scaled_pot, 0, 255);
                } else if (energyDisplay) {
                    float scaled_e = std::max(0.0f,std::min(255.0f,e/maxe * 255));
                    SDL_SetRenderDrawColor(renderer, scaled_e, 0, 0, 255);
                } else {
                    float scaled_dens = std::max(0.0f,std::min(255.0f,std::sqrt(density/maxDensity) * 255));
                    SDL_SetRenderDrawColor(renderer, scaled_dens, scaled_dens*(1-Y), 0, 255);
                }
                SDL_RenderFillRect(renderer, &rect);
            }
            i++;
        }
        inputFile.close();
    } else {
        std::cerr << "Error: Unable to open file." << std::endl;
    }

    // printf("%s\n",lines[0].c_str());
    // printf("%s\n",strip(lines[0], "//").c_str());
    // for (int i = 0; i < lines.size(); i++) {
    //     printf("%s",lines[i].c_str());
    // }
    



    return 1;
}