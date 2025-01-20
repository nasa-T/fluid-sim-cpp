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

std::ifstream inputFile("/Volumes/Give Me Space/outputGrid1000000.txt");

std::vector<std::string> lines;

if (inputFile.is_open()) {
    std::string line;
    while (std::getline(inputFile, line)) {
        lines.push_back(line);
    }
    inputFile.close();
} else {
    std::cerr << "Error: Unable to open file." << std::endl;
}



int main(int argv, char **argc) {
    for (int i = 0; i < lines.size(); i++) {
        printf(lines[i]);
    }
    // SDL_Renderer *renderer = NULL;
    // SDL_Window *window = NULL;
    // SDL_Surface *screenSurface;

    // SDL_Init(SDL_INIT_VIDEO);       // Initializing SDL as Video
    // SDL_CreateWindowAndRenderer(800, 800, 0, &window, &renderer);
    // SDL_SetRenderDrawColor(renderer, 0, 0, 0, 0);      // setting draw color
    // SDL_RenderClear(renderer);



    return 1;
}