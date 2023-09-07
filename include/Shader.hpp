#pragma once

#include "../include/Include.hpp"

string get_file_contents(const char* filename);

struct Shader_Program {
	GLuint ID;
	string Frag_Source, Program_Name;

	Shader_Program(const string& i_name) { Program_Name = i_name; };

	void Init(const char* fragmentFile);
	void ReCompile();

	void Activate();
	void Delete();
	void compileErrors(unsigned int shader, const char* type);
};