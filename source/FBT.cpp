#include "../include/FBT.h"

FBT::FBT(const int& i_width, const int& i_height) {
	glGenTextures(1, &ID);
	glBindTexture(GL_TEXTURE_2D, ID);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, i_width, i_height, 0, GL_RGBA, GL_FLOAT, NULL);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, ID, 0);

	glBindTexture(GL_TEXTURE_2D, 0);
}

void FBT::Resize(const int& i_width, const int& i_height) {
	GLuint new_ID;
	glGenTextures(1, &new_ID);
	glBindTexture(GL_TEXTURE_2D, new_ID);

	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, i_width, i_height, 0, GL_RGBA, GL_FLOAT, NULL);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, ID, 0);

	glBindTexture(GL_TEXTURE_2D, 0);
	ID = new_ID;
}

void FBT::Bind(const GLenum& i_texture_id) {
	glActiveTexture(i_texture_id);
	glBindTexture(GL_TEXTURE_2D, ID);
}

void FBT::Unbind() {
	glBindTexture(GL_TEXTURE_2D, 0);
}

void FBT::Delete() {
	glDeleteTextures(1, &ID);
}