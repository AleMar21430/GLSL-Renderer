#pragma once

#include "../include/Include.hpp"

struct EBO {
	GLuint ID;

	EBO(GLuint* indices, GLsizeiptr size);

	void Bind();
	void Unbind();
	void Delete();
};