#pragma once

#include "../include/Include.hpp"

struct VBO {
	GLuint ID;

	VBO(GLfloat* vertices, GLsizeiptr size);

	void Bind();
	void Unbind();
	void Delete();
};