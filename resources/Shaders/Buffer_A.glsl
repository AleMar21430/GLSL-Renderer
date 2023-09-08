#version 460 core

uniform float iTime;
uniform uint iFrame;
uniform vec2 iResolution;

uniform sampler2D iLastFrame;

in vec2 fragCoord;
in vec2 fragTexCoord;

out vec4 fragColor;

#define TWO_PI   6.28318530718
#define PI       3.14159265359
#define DEG_RAD  0.01745329252
#define MAX_DIST 500.0

#define SPP 2
#define SAMPLES 256

// GENERIC FUNCTIONS START

uvec4 white_noise_seed;
uvec2 pixel;

uvec4 hash(uvec4 v) {
	v = v * 1664525u + 1013904223u;
	v.x += v.y*v.w; v.y += v.z*v.x; v.z += v.x*v.y; v.w += v.y*v.z;
	v = v ^ (v>>16u);
	v.x += v.y*v.w; v.y += v.z*v.x; v.z += v.x*v.y; v.w += v.y*v.z;
	return v;
}

float rand() { return float(hash(white_noise_seed).x)   / float(0xffffffffu); }
vec2 rand2() { return vec2 (hash(white_noise_seed).xy)  / float(0xffffffffu); }
vec3 rand3() { return vec3 (hash(white_noise_seed).xyz) / float(0xffffffffu); }
vec4 rand4() { return vec4 (hash(white_noise_seed))     / float(0xffffffffu); }

void rng_initialize(vec2 pix, uint frame) {
	pixel = uvec2(pix);
	white_noise_seed = uvec4(pixel, frame, uint(pixel.x) + uint(pixel.y));
}

// GENERIC FUNCTIONS END

void main() {
	rng_initialize(gl_FragCoord.xy, iFrame);

	if (iFrame == 1) {
		fragColor = vec4(rand3(), 1.0);
	}
	else {
		fragColor = vec4(rand3(), 1.0);
	}
}