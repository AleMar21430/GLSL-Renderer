#include "./include/Include.hpp"

#include"./include/VAO.h"
#include"./include/VBO.h"
#include"./include/EBO.h"
#include"./include/FBO.h"
#include"./include/FBT.h"
#include"./include/Shader.h"

#define STB_IMAGE_IMPLEMENTATION
#include "./include/stb_image.h"

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
size_t current_frame = 0;
double iTime = 0.0;

FBT buffer_tex_a, last_frame_tex;
FBO FBO_main;
Shader_Program Buffer_A("Buffer A");
Shader_Program Main_Image("Main Image");

void framebuffer_size_callback(GLFWwindow* window, int width, int height) {
	glViewport(0, 0, width, height);
	Width = width;
	Height = height;
	current_frame = 0;
	iTime = glfwGetTime();
	FBO_main.Bind();
	buffer_tex_a.Resize(width, height);
	FBO_main.Unbind();
	last_frame_tex.Resize(width, height);
}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
	if (key == GLFW_KEY_R && action == GLFW_PRESS) {
		Buffer_A.ReCompile();
		Main_Image.ReCompile();

		FBO_main.Bind();
		buffer_tex_a.Resize(Width, Height);
		FBO_main.Unbind();
		last_frame_tex.Resize(Width, Height);

		current_frame = 0;
		iTime = glfwGetTime();
	}
}

int main() {
	glfwInit();

	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	GLFWwindow* window = glfwCreateWindow(Width, Height, "GLSL Shader", NULL, NULL);

	GLFWimage icons[1];
	icons[0].pixels = stbi_load("./resources/Icon.png", &icons[0].width, &icons[0].height, 0, 4);
	glfwSetWindowIcon(window, 1, icons);
	stbi_image_free(icons[0].pixels);

	if (window == NULL) {
		std::cout << "Failed to create GLFW window" << std::endl;
		glfwTerminate();
		return -1;
	}

	glfwMakeContextCurrent(window);
	gladLoadGL();

	glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
	glfwSetKeyCallback(window, key_callback);

	glViewport(0, 0, Width, Height);

	// Generates Shader object using shaders defualt.vert and default.frag
	Buffer_A.Init("./resources/Shaders/Buffer_A.glsl");
	Main_Image.Init("./resources/Shaders/Main_Image.glsl");

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
	FBO_main.Init();
	FBO_main.Bind();
	buffer_tex_a.Init(Width, Height);
	FBO_main.Unbind();

	last_frame_tex.Init(Width, Height);

	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	while (!glfwWindowShouldClose(window)) {
		double Time = glfwGetTime() - iTime;
		FBO_main.Bind();
		glClear(GL_COLOR_BUFFER_BIT);
		Buffer_A.Activate();
		VAO_main.Bind();

		glUniform1f(glGetUniformLocation(Buffer_A.ID, "iTime"), float(Time));
		glUniform1ui(glGetUniformLocation(Buffer_A.ID, "iFrame"), current_frame);
		glUniform2f(glGetUniformLocation(Buffer_A.ID, "iResolution"), Width, Height);
		last_frame_tex.Bind(GL_TEXTURE0);

		glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

		glCopyTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, 0, 0, Width, Height);
		last_frame_tex.Unbind();

		FBO_main.Unbind();

		Main_Image.Activate();

		glUniform1f(glGetUniformLocation(Main_Image.ID, "iTime"), float(Time));
		glUniform1ui(glGetUniformLocation(Main_Image.ID, "iFrame"), int(current_frame));
		glUniform2f(glGetUniformLocation(Main_Image.ID, "iResolution"), Width, Height);
		buffer_tex_a.Bind(GL_TEXTURE0);

		glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

		// Swap the back buffer with the front buffer
		glfwSwapBuffers(window);
		glfwPollEvents();
		current_frame++;
	}

	VAO_main.Delete();
	VBO_main.Delete();
	Faces.Delete();
	FBO_main.Delete();
	Buffer_A.Delete();
	Main_Image.Delete();
	buffer_tex_a.Delete();
	last_frame_tex.Delete();
	glfwDestroyWindow(window);
	glfwTerminate();
	return 0;
}