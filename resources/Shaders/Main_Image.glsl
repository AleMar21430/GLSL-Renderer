#version 460 core

uniform float iTime;
uniform int iFrame;
uniform vec2 iResolution;

uniform sampler2D iRawFrame;

in vec2 fragCoord;
in vec2 fragTexCoord;

out vec4 fragColor;

vec3 XYZtosRGB(vec3 XYZ)
{
	mat3 m = mat3(
		3.2404542, -1.5371385, -0.4985314,
		-0.9692660, 1.8760108, 0.0415560,
		0.0556434, -0.2040259, 1.0572252);

	return XYZ * m;
}

vec3 ColorGrade(vec3 vColor) {
	vec3 vHue = vec3(1.0, .7, .2);

	vec3 vGamma = 1.0 + vHue * 0.6;
	vec3 vGain = vec3(.9) + vHue * vHue * 8.0;

	vColor *= 1.5;

	float fMaxLum = 100.0;
	vColor /= fMaxLum;
	vColor = pow(vColor, vGamma);
	vColor *= vGain;
	vColor *= fMaxLum;
	return pow(tanh(vColor), vec3(0.57));
}

vec3 tone(vec3 c)
{
	c = XYZtosRGB(c);
	return ColorGrade(c);
}

void main() {
	vec4 acc = texture(iRawFrame, fragTexCoord);
	fragColor = vec4(tone(0.03 * acc.xyz / acc.w), 1.0);
	//fragColor = texture(iRawFrame, fragTexCoord);
}