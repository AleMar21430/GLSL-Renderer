#version 460 core

uniform float iTime;
uniform int iFrame;
uniform vec2 iResolution;

uniform sampler2D iLastFrame;

in vec2 fragCoord;
in vec2 fragTexCoord;

out vec4 fragColor;

#define time (iTime*0.1 + 5.0)

mat2 rot(float a) {
	float ca = cos(a);
	float sa = sin(a);
	return mat2(ca, sa, -sa, ca);
}

// Box SDF
float box(vec3 p, vec3 s) {
	p = abs(p) - s;
	return max(p.x, max(p.y, p.z));
}

// capsule SDF
float caps(vec3 p, vec3 p1, vec3 p2, float s) {
	vec3 pa = p - p1;
	vec3 pb = p2 - p1;
	float prog = dot(pa, pb) / dot(pb, pb);
	prog = clamp(prog, 0., 1.);
	return length(p1 + pb * prog - p) - s;
}

// to switch between the 2 scenes
int scene = 0;

// first SDF function with the refractive geometry
float map(vec3 p) {

	// sphere in cube
	vec3 p2 = p;
	float t = time * 0.1;
	p2.yz *= rot(t);
	p2.yx *= rot(t * 1.3);
	float d4 = max(box(p2, vec3(3)), 1.2 - length(p));

	// KIFS with spheres and cubes
	vec3 pb = p;
	float d2 = 10000.;
	for (float i = 0.; i < 3.; ++i) {
		float t = time * 0.03 + i;
		p.yz *= rot(t + i);
		p.yx *= rot(t * 1.3);
		d2 = min(d2, length(p) - 0.47);
		p = abs(p);
		p -= 0.9;
	}


	float d = box(p, vec3(0.4, 0.4, 0.4));

	d = min(d, d2);
	//d=d2;

	if (scene == 0) d = d4;
	//d=d4;

	return d;
}

// we will have a maximum of 10 bounces of the laser
#define pcount 10
vec3 points[pcount];
int pid = 1;

float atm = 0.;
// second SDF function with both refractive geometry and the laser geometry
float map2(vec3 p) {

	// get refractive geometry
	float d = map(p);

	// loop over laser's collisions, insert it into the SDF and accumulate laser's light
	float d2 = 10000.;
	for (int i = 0; i < pid - 1; ++i) {

		// one capsule per laser part
		float d3 = caps(p, points[i], points[i + 1], 0.01);

		// I use the smoothstep to have more contrast between lit and unlit part of the geometry
		atm += 0.013 / (0.05 + abs(d3)) * smoothstep(4., 0.3, d3);

		d2 = min(d2, d3);
	}


	return min(abs(d), d2);
}

float rnd(vec2 uv) {
	return fract(dot(sin(uv * 452.714 + uv.yx * 547.524), vec2(352.887)));
}

void main()
{
	vec2 uv = vec2(fragCoord.x - 0.5, fragCoord.y - 0.5);
	//uv -= 0.5;
	uv /= vec2(iResolution.y / iResolution.x, 1);

	// change scene every 10s
	scene = int(mod(floor(time / 10.), 2.));

	// first raymarching to get the laser collisions
	vec3 s2 = vec3(10, 0, 0);
	vec3 r2 = normalize(vec3(-1, sin(float(time)) * 0.1, 0));

	// get a random refractive index different per pixel
	float ior = (rnd(uv + fract(iTime * .1)) - 0.5);
	//ior = (fract(gl_FragCoord.y/3.)-.5);
	float id = ior * 2.;
	// compute index of refraction associated color 
	vec3 diff = 1.3 - vec3(1. + id, 0.45 + abs(id), 1. - id);

	vec3 p2 = s2;
	points[0] = p2;
	pid = 1;
	float side = 1.;
	for (int i = 0; i < 60; ++i) {
		// we use only refractive geometry SDF
		float d = abs(map(p2));
		if (d < 0.001) {
			// when laser collide with something, store the collision point into the list
			points[pid] = p2;
			pid += 1;
			if (pid >= pcount - 1) break;

			// and we compute the refracted direction
			vec2 off = vec2(0.01, 0);
			vec3 n2 = side * normalize(d - vec3(map(p2 - off.xyy), map(p2 - off.yxy), map(p2 - off.yyx)));
			//r2=reflect(r2,n2);
			vec3 r3 = refract(r2, n2, 1. - side * (0.3 + 0.1 * ior));
			if (dot(r3, r3) < 0.5) r3 = reflect(r2, n2);
			r2 = r3;
			side = -side;
			d = 0.1;
		}
		if (d > 100.0) break;
		p2 += r2 * d;
	}
	points[pid] = p2 + r2 * 1000.;
	++pid;

	// second raymarching, what we will see on screen
	vec3 s = vec3(0, 0, -10);
	vec3 r = normalize(vec3(uv, 1));

	// dithering noise to eliminate some banding
	float mumu = mix(rnd(-uv + fract(float(iFrame) * .1)), 1., 0.9);
	vec3 p = s;
	float side2 = 1.;
	for (int i = 0; i < 90; ++i) {
		float d = abs(map2(p));
		if (d < 0.001) {
			// collision, compute refracted direction
			// we use same per-pixel color and index of refraction than laser's raymarch
			// it's not 100% correct but hey, it looks good
			vec2 off = vec2(0.01, 0);
			vec3 n = side2 * normalize(d - vec3(map(p - off.xyy), map(p - off.yxy), map(p - off.yyx)));
			vec3 r3 = refract(r, n, 1. - side2 * (0.3 + 0.1 * ior));
			if (dot(r3, r3) < 0.5) r3 = reflect(r, n);
			r = r3;

			side2 = -side2;
			d = 0.1;
			//break;
		}
		if (d > 100.0) break;
		p += r * d * mumu;
	}

	vec3 col = vec3(0);
	col += diff * atm;

	/*
	// original shader was doing gamma here but in shadertoy we can do it in post-effect and it's better
	col=smoothstep(0.01,0.9,col);
	col=pow(col, vec3(0.4545));
	*/

	// use previous frame to have a feedback effect and try to reduce the noise a bit
	vec4 prev = texture(iLastFrame, fragTexCoord);
	prev.w *= 0.98;
	// you can try putting 0.95 instead for lot less noise but also more blurry image
	vec4 col1 = vec4(prev.xyz * prev.w, prev.w) + vec4(col, 1.0);
	fragColor = vec4(col1.xyz / col1.w, col1.w);
}