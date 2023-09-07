#version 460 core

uniform float iTime;
uniform int iFrame;
uniform vec2 iResolution;

uniform sampler2D iRawFrame;

in vec2 fragCoord;
in vec2 fragTexCoord;

out vec4 fragColor;

void main() {
	fragColor = vec4(tanh(6.0 * texture(iRawFrame, fragTexCoord).xyz), 1);
}