#pragma once

#include "../include/Include.hpp"

struct FBO {
	GLuint ID;

	FBO();

	void Bind();
	void Unbind();
	void Delete();
};