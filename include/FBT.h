#pragma once

#include "../include/Include.hpp"

struct FBT {
	GLuint ID;

	FBT(const int& i_width, const int& i_height);

	void Resize(const int& i_width, const int& i_height);
	void Bind(const GLenum& i_texture_id = GL_TEXTURE0);
	void Unbind();
	void Delete();
};