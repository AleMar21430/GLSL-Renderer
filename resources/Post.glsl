#version 460 core

uniform float iTime;
uniform int iFrame;
uniform vec2 iResolution;

uniform sampler2D iRawFrame;

in vec2 fragCoord;
in vec2 fragTexCoord;

out vec4 fragColor;

void main() {
	fragColor = texture(iRawFrame, fragTexCoord) * vec4(1.0, 0.25, 0.25, 1.0);
}