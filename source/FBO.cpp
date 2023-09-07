#include "../include/FBO.hpp"

void FBO::Init() {
	glGenFramebuffers(1, &ID);
}

void FBO::Bind() {
	glBindFramebuffer(GL_FRAMEBUFFER, ID);
}

void FBO::Unbind() {
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

void FBO::Delete() {
	glDeleteFramebuffers(1, &ID);
}