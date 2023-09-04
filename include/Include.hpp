#pragma once

#include <unordered_map>
#include <unordered_set>
#include <filesystem>
#include <windows.h>
#include <stdexcept>
#include <iostream>
#include <iostream>
#include <optional>
#include <direct.h>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <variant>
#include <cstdlib>
#include <cstring>
#include <numeric>
#include <vector>
#include <chrono>
#include <thread>
#include <future>
#include <math.h>
#include <array>
#include <tuple>
#include <any>
#include <map>
#include <set>

#define GLFW_INCLUDE_VULKAN
//#include <vulkan/vulkan.h>
#include <GLFW/glfw3.h>
#include "glm/glm.hpp"

using namespace std;
//using namespace glm;

#define PI          3.141592653589793
#define TWO_PI      6.283185307179586
#define INVERTED_PI 0.318309886183791
#define DEG_RAD     0.017453292519943
#define RED         dvec4( 1   ,0  , 0   , 1 )
#define GREEN       dvec4( 0  , 1  , 0   , 1 )
#define BLUE        dvec4( 0  , 0  , 1   , 1 )
#define WHITE       dvec4( 1  , 1  , 1   , 1 )
#define BLACK       dvec4( 0  , 0  , 0   , 1 )
#define GRAY        dvec4( 0.5, 0.5, 0.5 , 1 )