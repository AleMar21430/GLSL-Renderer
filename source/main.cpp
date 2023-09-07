#include "../include/Include.hpp"

#include"../include/VAO.h"
#include"../include/VBO.h"
#include"../include/EBO.h"
#include"../include/FBO.h"
#include"../include/FBT.h"
#include"../include/Shader.h"

GLfloat vertices[] = {
	-1.0f, -1.0f, 0.0f, 0.0f,
	 1.0f, -1.0f, 1.0f, 0.0f,
	 1.0f,  1.0f, 1.0f, 1.0f,
	-1.0f,  1.0f, 0.0f, 1.0f,
};
GLuint indices[] = {
	0, 1, 2,
	2, 3, 0
};
uint16_t Width = 860;
uint16_t Height = 540;

void framebuffer_size_callback(GLFWwindow* window, int new_width, int new_height) {
	glViewport(0, 0, new_width, new_height);
	Width = new_width;
	Height = new_height;
}

int main() {
	glfwInit();

	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	GLFWwindow* window = glfwCreateWindow(Width, Height, "Path Tracer", NULL, NULL);

	if (window == NULL) {
		std::cout << "Failed to create GLFW window" << std::endl;
		glfwTerminate();
		return -1;
	}

	glfwMakeContextCurrent(window);
	gladLoadGL();

	glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

	glViewport(0, 0, Width, Height);

	// Generates Shader object using shaders defualt.vert and default.frag
	Shader Buffer_A("./resources/Vert.glsl", "./resources/Frag.glsl");

	// VERTICES //
	VAO VAO_main;
	VAO_main.Bind();
	VBO VBO_main(vertices, sizeof(vertices));
	EBO Faces(indices, sizeof(indices));
	VAO_main.LinkAttrib(VBO_main, 0, 2, GL_FLOAT, 4 * sizeof(float), (void*)0);
	VAO_main.LinkAttrib(VBO_main, 1, 2, GL_FLOAT, 4 * sizeof(float), (void*)(2 * sizeof(float)));
	VAO_main.Unbind();
	VBO_main.Unbind();
	Faces.Unbind();

	// FBOs //
	FBO FBO_main;
	FBO_main.Bind();
	FBT FBT_main(Width, Height);
	FBO_main.Unbind();

	FBT last_frame(Width, Height);

	size_t current_frame = 0;

	glClearColor(0.05f, 0.05f, 0.05f, 1.0f);
	while (!glfwWindowShouldClose(window)) {
		glClear(GL_COLOR_BUFFER_BIT);
		Buffer_A.Activate();
		VAO_main.Bind();

		glUniform1f(glGetUniformLocation(Buffer_A.ID, "iTime"), glfwGetTime());
		glUniform1i(glGetUniformLocation(Buffer_A.ID, "iFrame"), current_frame);
		glUniform2f(glGetUniformLocation(Buffer_A.ID, "iResolution"), Width, Height);
		last_frame.Bind();

		glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

		glCopyTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, 0, 0, Width, Height);
		last_frame.Unbind();

		// Swap the back buffer with the front buffer
		glfwSwapBuffers(window);
		glfwPollEvents();
		current_frame++;
	}

	VAO_main.Delete();
	VBO_main.Delete();
	Faces.Delete();
	Buffer_A.Delete();
	glfwDestroyWindow(window);
	glfwTerminate();
	return 0;
}