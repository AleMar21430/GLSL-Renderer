#include"../include/Shader.h"

string get_file_contents(const char* filename) {
	ifstream in(filename, ios::binary);
	if (in) {
		string contents;
		in.seekg(0, ios::end);
		contents.resize(in.tellg());
		in.seekg(0, ios::beg);
		in.read(&contents[0], contents.size());
		in.close();
		return(contents);
	}
	throw(errno);
}

string getSubstringAfterDelimiter(const string& input, const string& delimiter) {
	size_t delimiterPos = input.find(delimiter);

	if (delimiterPos != string::npos) {
		return input.substr(delimiterPos + delimiter.length());
	}
	else {
		return input;
	}
}

string insertLinesIfDelimiterFound(const string& input, const string& delimiter, const string& lineToInsert) {
	string result;
	size_t pos = 0;
	size_t prevPos = 0;

	while ((pos = input.find(delimiter, prevPos)) != string::npos) {
		result += input.substr(prevPos, pos - prevPos);
		result += lineToInsert;
		prevPos = pos + delimiter.length();
	}
	result += input.substr(prevPos);
	return result;
}

void Shader_Program::Init(const char* fragmentFile) {
	Frag_Source = fragmentFile;
	ReCompile();
}

void Shader_Program::ReCompile() {
	glDeleteProgram(ID);

	string vertexCode = get_file_contents("./resources/Shaders/Vert.glsl");
	string fragmentCode = get_file_contents(Frag_Source.c_str());

	const char* vertexSource = vertexCode.c_str();
	const char* fragmentSource = fragmentCode.c_str();

	GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
	glShaderSource(vertexShader, 1, &vertexSource, NULL);
	glCompileShader(vertexShader);
	compileErrors(vertexShader, ("VERTEX " + Program_Name).c_str());

	GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
	glShaderSource(fragmentShader, 1, &fragmentSource, NULL);
	glCompileShader(fragmentShader);
	compileErrors(fragmentShader, ("FRAGMENT " + Program_Name).c_str());

	ID = glCreateProgram();
	glAttachShader(ID, vertexShader);
	glAttachShader(ID, fragmentShader);
	glLinkProgram(ID);
	compileErrors(ID, ("PROGRAM " + Program_Name).c_str());

	glDeleteShader(vertexShader);
	glDeleteShader(fragmentShader);
}

void Shader_Program::Activate() {
	glUseProgram(ID);
}

void Shader_Program::Delete() {
	glDeleteProgram(ID);
}

void Shader_Program::compileErrors(unsigned int shader, const char* type) {
	GLint hasCompiled;
	char infoLog[1024];
	if (type != "PROGRAM") {
		glGetShaderiv(shader, GL_COMPILE_STATUS, &hasCompiled);
		if (hasCompiled == GL_FALSE) {
			glGetShaderInfoLog(shader, 1024, NULL, infoLog);
			cout << "SHADER_COMPILATION_ERROR for:" << type << "\n" << infoLog << endl;
		}
	}
	else {
		glGetProgramiv(shader, GL_LINK_STATUS, &hasCompiled);
		if (hasCompiled == GL_FALSE) {
			glGetProgramInfoLog(shader, 1024, NULL, infoLog);
			cout << "SHADER_LINKING_ERROR for:" << type << "\n" << infoLog << endl;
		}
	}
}