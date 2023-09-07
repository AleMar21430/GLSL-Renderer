#version 460 core

in vec2 fragCoord;
in vec3 fragTexCoord;

out vec4 fragColor;

void main() {
	fragColor = vec4(fragCoord.x, 0.0, fragCoord.y, 1.0);
}