#pragma once

#include "../include/Include.hpp"

std::string get_file_contents(const char* filename);

struct Shader {
	GLuint ID;

	Shader(const char* vertexFile, const char* fragmentFile);

	void Activate();
	void Delete();
	void compileErrors(unsigned int shader, const char* type);
};