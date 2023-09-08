#version 460 core

uniform float iTime;
uniform int iFrame;
uniform vec2 iResolution;

uniform sampler2D iLastFrame;

in vec2 fragCoord;
in vec2 fragTexCoord;

out vec4 fragColor;


/*
Stanislav Pidhorskyi 2020

This one is inspired by the book by Peter Shirley "Ray Tracing: The Rest of Your Life".


I implemented a spectral path tracer, where each ray has a random wavelength. Surprisingly, it's not much noisier than conventional RGB path tracing.
Wavelength is sampled uniformly from 400.0nm to 700.0nm, but of course, it would be better to utilize importance sampling here and sample more those wavelengths that contribute to the RGB color more.

Specter can be converted to the color values by computing dot products with color matching functions.
I used CIE color-matching functions that I took from John Walker website: https://www.fourmilab.ch/documents/specrend/
These functions are tabulated and correspond to the CIE tristimulus values X, Y, and Z.
The lookup table would not be too large, but it seems to be a significant source of slowdown, so I decided to approximate them with Gaussians.
But before with dive into the approximation, it worth noting that white color does not have a uniform specter.
In order to get sources of white light, we need to use a specific specter, like Illuminant D65.
I took the values for D65 from here: https://www.waveformlighting.com/tech/calculate-illuminant-d-spd-and-cie-1931-xy-from-color-temperature/
I computed the point-wise product of D65 specter and CIE color-matching functions so that we could pretend that now the white color has a uniform specter.
The result functions I approximated with Gaussians and a mixture of two Gaussians for the Red color:


Red:    	8233.31 * f(l | mu=593.95, sigma=34.00) + 1891.26 * f(l | mu=448.89, sigma=18.785),
Green:      10522.64 * f(l | 555.38, 40.80),
Blue:       111254.78 * f(l | mu=452.98, sigma=21.57)

where f(l) = 1/(sigma sqrt(2 * pi)) * exp(-1/2 * ((l - mu) / sigma)^2)

Code for computing the coefficients can be found here: https://gist.github.com/podgorskiy/99c283773f7cee8e71386aa8ef622fdf

Represented in this way, it would be actually not hard to compute a CDF and perform importance sampling, but I did not do it. I just used uniform sampling for wavelengths.

To compute the XYZ values given a specter, we need to compute the integral of the point-wise product of the specter with the color matching functions.
The integral is computed in the Monte Carlo way by sampling wavelength at random and then accumulating the XYZ values. The averaged XYZ value is then converted to RGB

In order to actually appreciate the spectral path tracing, we need to make the refraction index to be a function of the wavelength. Otherwise, the output of the spectral raytracer won't differ from the conventional one.

I use the Sellmeier equation https://en.wikipedia.org/wiki/Sellmeier_equation. The material of the prism is Flint glass, parameters I found here: https://refractiveindex.info/?shelf=glass&book=SCHOTT-F&page=F2.

To sample from the cubemap, we need to do an inverse operation, go from RGB values to a specter. Sure, that's an ill-posed problem, but we could approach it with some basis of orthogonal functions and Non-negative least squares (NNLS), but that would be too costly. Instead, I used a very naive way and assumed that the CIE color-matching functions are orthogonal (which are not) and used just a dot product.
Because of such a naive assumption, colors are off.

And finally, I wrote the code for a ray-prism intersection because I could not find one.

----------------------------------------

Reuses some code from https://www.shadertoy.com/view/lssBD7 that implements Peter Shirley's "Ray Tracing in One Weekend" pathtracer.

Post-processing from Inigo Quilezh's `Happy Jumping` ttps://www.shadertoy.com/view/3lsSzf


*/

#define PI 			3.1415926535
#define MAXFLOAT	99999.99

#define MAXDEPTH 	8
#define NUMSAMPLES 	64


// Do not unroll.
#define DNU(X)  min(int(iResolution.x) << 12, X)

const mat3 XYZ_2_RGB = (mat3(
	3.2404542, -0.9692660, 0.0556434,
	-1.5371385, 1.8760108, -0.2040259,
	-0.4985314, 0.0415560, 1.0572252
));

const mat3 RGB_2_XYZ = (mat3(
	0.4124564, 0.2126729, 0.0193339,
	0.3575761, 0.7151522, 0.1191920,
	0.1804375, 0.0721750, 0.9503041
));

uint seedA = 0x9c127997U;
uint seedB = 0x140b75b2U;

uint urand()
{
	uint x = seedA;
	uint y = seedB;
	seedA = y;
	x ^= x << 23;
	seedB = x ^ y ^ (x >> 17) ^ (y >> 26);
	uint n = seedB + y;
	return n * (n * n * 15731U + 789221U) + 1376312589U;
}

float rand1()
{
	return uintBitsToFloat((urand() >> 9U) | 0x3f800000U) - 1.0;
}

vec2 rand2()
{
	float a = uintBitsToFloat((urand() >> 9U) | 0x3f800000U) - 1.0;
	float b = uintBitsToFloat((urand() >> 9U) | 0x3f800000U) - 1.0;
	return vec2(a, b);
}


float rand(float a, float b)
{
	return rand1() * (b - a) + a;
}

// Random unit vector, uniformly distributed on a sphere
vec3 random_unit_vector()
{
	float a = rand(0.0, 2. * PI);
	float z = rand(-1.0, 1.0);
	float r = sqrt(1.0 - z * z);
	return vec3(r * cos(a), r * sin(a), z);
}

// Random, uniformely distributed point on unit disk
vec2 random_in_unit_disk()
{
	vec2 uv = rand2();
	float theta = 6.283185 * uv.x;
	return sqrt(uv.y) * vec2(cos(theta), sin(theta));
}

vec3 random_cosine_direction(in vec3 normal)
{
	return normalize(normal + random_unit_vector());
}

float gaussian(float x, float mu, float sigma)
{
	return 1.0 / (sigma * sqrt(2.0 * PI)) * exp(-(x - mu) * (x - mu) / (2. * sigma * sigma));
}


// The CIE color matching functions were taken from  https://www.fourmilab.ch/documents/specrend
// The tabulated functions then were approximated with gaussians (for G and B) and with a mixture of two gaussiuns (R).
vec3 wavelength2XYZ(float l)
{
	return vec3(
		8233.31 * gaussian(l, 593.95, 34.00) + 1891.26 * gaussian(l, 448.89, 18.785),
		10522.64 * gaussian(l, 555.38, 40.80),
		11254.78 * gaussian(l, 452.98, 21.57)
	);
}


// Very very crude convertion from XYZ color to spectrum. It is wrong in many ways, but I don't know a faster way to do it.
// The assumption is that the three color matching functions are orthogonal, but they are not.
float XYZ2WavelengthApprox(float l, vec3 color)
{
	return dot(wavelength2XYZ(l), color) / 100.0;
}

float glass_crown_n(float l)
{
	l /= 1000.0;
	float l2 = l * l;
	vec3 k = vec3(1.03961212, 0.231792344, 1.01046945) * l2 /
		(l2 - vec3(0.0060007, 0.020018, 103.560));
	return 1.0 * (sqrt(1. + dot(k, vec3(1.0))) - 1.52) + 1.52;
}


float glass_flint_n(float l)
{
	l /= 1000.0;
	float l2 = l * l;
	vec3 k = vec3(1.34533359, 0.209073176, 0.937357162) * l2 /
		(l2 - vec3(0.00997743871, 0.0470450767, 111.886764));
	return sqrt(1. + dot(k, vec3(1.0)));
}

struct Ray
{
	vec3 origin;
	vec3 direction;
	float wave_length;
};

struct Material
{
	int   materialType;
	vec3  albedo;
	float fuzz;
	float refractionIndex;
};

struct IntersectInfo
{
	vec2 t;
	vec3  p;
	vec3  normal;
	Material mat;
};

struct Transform
{
	mat3 rotation;
	vec3 position;
	vec3 scale;
};

struct Object
{
	Transform transform;
	Material mat;
};

vec3 InvertseTransformPosition(Transform tr, vec3 p)
{
	return (tr.rotation * (p - tr.position)) / tr.scale;
}

vec3 DirectTransformPosition(Transform tr, vec3 p)
{
	return (p * tr.scale) * tr.rotation + tr.position;
}

vec3 InvertseTransformDirection(Transform tr, vec3 d)
{
	return normalize((tr.rotation * d) / tr.scale);
}

vec3 DirectTransformDirection(Transform tr, vec3 d)
{
	return normalize((d * tr.scale) * tr.rotation);
}

IntersectInfo Box(Object box, Ray ray, vec2 t)
{
	Transform tr = box.transform;
	IntersectInfo rec;
	rec.t.x = MAXFLOAT;
	vec3 o = InvertseTransformPosition(tr, ray.origin);
	vec3 d = InvertseTransformDirection(tr, ray.direction);
	vec3 m = 1.0 / d;
	vec3 k = abs(m);
	vec3 a = -m * o - k;
	vec3 b = a + k * 2.0;
	float near = max(max(a.x, a.y), a.z);
	float far = min(min(b.x, b.y), b.z);
	if (near > far)
	{
		return rec;
	}

	vec3 localPosFar = o + far * d;
	vec3 localPosNear = o + near * d;
	vec3 posFar = DirectTransformPosition(tr, localPosFar);
	vec3 posNear = DirectTransformPosition(tr, localPosNear);

	float dmax = dot(posFar - ray.origin, ray.direction);
	float dmin = dot(posNear - ray.origin, ray.direction);

	if (dmin < t.y && dmin > t.x)
	{
		rec.t = vec2(dmin, dmax);
		rec.p = posNear;
		vec3 f = step(0.999, abs(localPosNear));
		rec.normal = normalize(localPosNear * f * tr.rotation);
		rec.mat = box.mat;
	}
	else
		if (dmax < t.y && dmax > t.x)
		{
			rec.t = vec2(dmax, MAXFLOAT);
			rec.p = posFar;
			vec3 f = step(0.999, abs(localPosFar));
			rec.normal = normalize(localPosFar * f * tr.rotation);
			rec.mat = box.mat;
		}
	return rec;
}


IntersectInfo Prism(Object box, Ray ray, vec2 t)
{
	Transform tr = box.transform;
	IntersectInfo rec;
	rec.t.x = MAXFLOAT;
	vec3 o = InvertseTransformPosition(tr, ray.origin);
	vec3 d = InvertseTransformDirection(tr, ray.direction);
	vec3 m = 1.0 / d;
	vec3 k = abs(m);
	vec3 a = -m * o - k;
	vec3 b = a + k * 2.0;
	float near = max(max(a.x, a.y), a.z);
	float far = min(min(b.x, b.y), b.z);
	if (near > far)
	{
		return rec;
	}

	vec3 pn1 = normalize(vec3(2.0, 1.0, 0.0));
	vec3 pn2 = pn1;
	pn2.x = -pn2.x;
	float pm1 = -1.0 / dot(pn1, d);
	float pm2 = -1.0 / dot(pn2, d);
	float p1 = dot(o, pn1) * pm1 - pn1.y * pm1;
	float p2 = dot(o, pn2) * pm2 - pn2.y * pm2;
	float pb = o.y * m.y + m.y;

	if (pm1 > 0.0)
		near = max(near, p1);
	else
		far = min(far, p1);
	if (pm2 > 0.0)
		near = max(near, p2);
	else
		far = min(far, p2);

	if (near > far)
	{
		return rec;
	}

	vec3 localPosFar = o + far * d;
	vec3 localPosNear = o + near * d;
	vec3 posFar = DirectTransformPosition(tr, localPosFar);
	vec3 posNear = DirectTransformPosition(tr, localPosNear);

	float dmax = dot(posFar - ray.origin, ray.direction);
	float dmin = dot(posNear - ray.origin, ray.direction);

	if (dmin < t.y && dmin > t.x)
	{
		rec.t = vec2(dmin, dmax);
		rec.p = posNear;
		vec3 f = step(0.9999, abs(localPosNear));
		rec.normal = normalize(sign(localPosNear) * f);
		if (all(lessThan(f, vec3(1.0))))
		{
			rec.normal = mix(pn1, pn2, bvec3(localPosNear.x < 0.0));
		}
		rec.normal *= tr.rotation;
		rec.mat = box.mat;
	}
	else
		if (dmax < t.y && dmax > t.x)
		{
			rec.t = vec2(dmax, MAXFLOAT);
			rec.p = posFar;
			vec3 f = step(0.9999, abs(localPosFar));
			rec.normal = normalize(localPosFar * f);
			if (all(lessThan(f, vec3(1.0))))
			{
				rec.normal = mix(pn1, pn2, bvec3(localPosFar.x < 0.0));
			}
			rec.normal *= tr.rotation;
			rec.mat = box.mat;
		}
	return rec;
}



IntersectInfo Sphere(Object sphere, Ray ray, vec2 t)
{
	Transform tr = sphere.transform;
	IntersectInfo rec;
	rec.t.x = MAXFLOAT;
	vec3 o = InvertseTransformPosition(tr, ray.origin);
	vec3 d = InvertseTransformDirection(tr, ray.direction);
	float a = dot(d, d);
	float b = 2.0 * dot(d, o);
	float c = dot(o, o) - 1.0;

	float discriminant = b * b - 4.0 * a * c;

	if (discriminant < 0.0f)
	{
		return rec;
	}

	float D = sqrt(discriminant);
	vec2 temp = vec2(D - b, -D - b) / 2.0 / a;

	vec3 localPosMax = o + temp.x * d;
	vec3 localPosMin = o + temp.y * d;
	vec3 posMax = DirectTransformPosition(tr, localPosMax);
	vec3 posMin = DirectTransformPosition(tr, localPosMin);

	float dmax = dot(posMax - ray.origin, ray.direction);
	float dmin = dot(posMin - ray.origin, ray.direction);

	if (dmin < t.y && dmin > t.x)
	{
		rec.t = vec2(dmin, dmax);
		rec.p = posMin;
		rec.normal = normalize(localPosMin * tr.scale * tr.rotation);
		rec.mat = sphere.mat;
	}
	else
		if (dmax < t.y && dmax > t.x)
		{
			rec.t = vec2(dmax, MAXFLOAT);
			rec.p = posMax;
			rec.normal = normalize(localPosMax * tr.scale * tr.rotation);
			rec.mat = sphere.mat;
		}
	return rec;
}

// Schlick's approximation for approximating the contribution of the Fresnel factor
// in the specular reflection of light from a non-conducting surface between two media
//
// Theta is the angle between the direction from which the incident light is coming and
// the normal of the interface between the two media
float schlick(float cos_theta, float n2)
{
	const float n1 = 1.0f;  // refraction index for air

	float r0s = (n1 - n2) / (n1 + n2);
	float r0 = r0s * r0s;

	return r0 + (1.0f - r0) * pow((1.0f - cos_theta), 5.0f);
}

bool refractVec(vec3 v, vec3 n, float ni_over_nt, out vec3 refracted)
{
	vec3 uv = normalize(v);
	float dt = dot(uv, n);
	float discriminant = 1.0 - ni_over_nt * ni_over_nt * (1.0f - dt * dt);
	refracted = ni_over_nt * (uv - n * dt) - n * sqrt(discriminant);
	return discriminant > 0.0f;
}

vec3 reflectVec(vec3 v, vec3 n)
{
	return v - 2.0f * dot(v, n) * n;
}


struct Camera{
	float apertureRadius;
	vec3 origin;
	mat3 cam;
	mat3 k_inv;
	float focusDist;
};


void Camera_init(out Camera camera, vec3 lookfrom, vec3 lookat, vec3 vup, float vfov, float aspect, float aperture, float focusDist)
{
	camera.apertureRadius = aperture;
	camera.focusDist = focusDist;

	float theta = vfov * PI / 180.0;
	float halfHeight = tan(theta / 2.0);
	float halfWidth = aspect * halfHeight;

	camera.origin = lookfrom;

	vec3 z = normalize(lookfrom - lookat);
	vec3 x = normalize(cross(vup, z));
	vec3 y = normalize(cross(z, x));

	camera.cam = mat3(x, y, z);
	camera.k_inv = mat3(vec3(halfWidth, 0, 0), vec3(0, halfHeight, 0), vec3(0.0, 0.0, -1.));
}


Ray Camera_getRay(Camera camera, vec2 uv)
{
	vec3 rd = camera.apertureRadius * vec3(random_in_unit_disk(), 0.0);

	Ray ray;
	vec2 pixel = 1.0 / iResolution.xy;
	uv += pixel * (rand2() - 0.5);
	uv = 2.0 * uv - 1.0;

	ray.origin = camera.origin;
	ray.direction = camera.cam * (normalize(camera.k_inv * vec3(uv, 1.0)));

	return ray;
}

#define LAMBERT    0
#define METAL      1
#define DIELECTRIC 2
#define LIGHT 3

bool Material_bsdf(IntersectInfo rec, Ray wo, out Ray wi, out float attenuation, out float emission)
{
	int materialType = rec.mat.materialType;

	if (materialType == LAMBERT)
	{
		wi.origin = rec.p;
		wi.direction = random_cosine_direction(rec.normal);
		attenuation = rec.mat.albedo.x;
		return true;
	}
	else if (materialType == LIGHT)
	{
		emission = 1.0;
		return false;
	}
	else
		if (materialType == METAL)
		{
			vec3 reflected = reflect(normalize(wo.direction), rec.normal);

			wi.origin = rec.p;
			wi.direction = reflected + rec.mat.fuzz * random_cosine_direction(rec.normal);

			attenuation = rec.mat.albedo.x;

			return (dot(wi.direction, rec.normal) > 0.0f);
		}
		else
			if (materialType == DIELECTRIC)
			{
				vec3 outward_normal;
				vec3 reflected = reflect(wo.direction, rec.normal);

				float ni_over_nt;

				attenuation = 1.0f;
				vec3 refracted;
				float reflect_prob;
				float cosine;

				float rafractionIndex = glass_flint_n(wo.wave_length);

				if (dot(wo.direction, rec.normal) > 0.0f)
				{
					outward_normal = -rec.normal;
					ni_over_nt = rafractionIndex;

					cosine = dot(wo.direction, rec.normal) / length(wo.direction);
					cosine = sqrt(1.0f - rafractionIndex * rafractionIndex * (1.0f - cosine * cosine));
				}
				else
				{
					outward_normal = rec.normal;
					ni_over_nt = 1.0f / rafractionIndex;
					cosine = -dot(wo.direction, rec.normal) / length(wo.direction);
				}
				if (refractVec(wo.direction, outward_normal, ni_over_nt, refracted))
					reflect_prob = schlick(cosine, rafractionIndex);
				else
					reflect_prob = 1.0f;

				if (rand1() < reflect_prob)
				{
					wi.origin = rec.p;
					wi.direction = reflected;
				}
				else
				{
					wi.origin = rec.p;
					wi.direction = refracted;
				}

				return true;
			}

	return false;
}

Object prism = Object(
	Transform(mat3(1.0), vec3(0.0, 2.0, 0.0), vec3(2.0, 2.0, 2.0)),
	Material(2, vec3(0.193093, 0.510542, 0.613362), 0.200000, 1.500000));

Object table = Object(
	Transform(mat3(1.0), vec3(-4.5, -2.1, 0.0), vec3(5.5, 0.1, 5.0)),
	Material(0, vec3(0.500000, 0.500000, 0.500000), 1.000000, 1.500000));

Object prop1 = Object(
	Transform(mat3(1.0), vec3(0.0, -1.001, -1.9), vec3(1.0, 1.0, 0.1)),
	Material(0, vec3(0.500000, 0.500000, 0.500000), 1.000000, 1.500000));

Object prop2 = Object(
	Transform(mat3(1.0), vec3(0.0, -1.001, 1.9), vec3(1.0, 1.0, 0.1)),
	Material(0, vec3(0.500000, 0.500000, 0.500000), 1.000000, 1.500000));

Object light = Object(
	Transform(mat3(1.0), vec3(-9.0, -1.0, 0.0), vec3(1.0, 1.0, 1.0)),
	Material(3, vec3(0.500000, 0.500000, 0.500000), 1.000000, 1.500000));

Object wall = Object(
	Transform(mat3(1.0), vec3(10.0, 0.0, 0.0), vec3(0.1, 10.0, 10.0)),
	Material(0, vec3(0.500000, 0.500000, 0.500000), 1.000000, 1.500000));

IntersectInfo Add(IntersectInfo a, IntersectInfo b)
{
	if (a.t.x < b.t.x) return a; else return b;
}

bool intersectScene(Ray ray, vec2 t, out IntersectInfo rec)
{
	IntersectInfo hit;
	hit.t.x = MAXFLOAT;

	hit = Add(hit, Box(table, ray, t));
	hit = Add(hit, Box(prop1, ray, t));
	hit = Add(hit, Box(prop2, ray, t));
	hit = Add(hit, Box(light, ray, t));
	hit = Add(hit, Prism(prism, ray, t));

	rec = hit;
	return hit.t.x != MAXFLOAT;
}

float skyColor(Ray ray)
{
	vec3 sky = vec3(0.5);
	sky = RGB_2_XYZ * pow(sky, vec3(2.2));
	return XYZ2WavelengthApprox(ray.wave_length, sky);
}


float radiance(Ray ray)
{
	IntersectInfo rec;
	Material mat;
	float intensity = 1.0;

	for (int i = 0; i < DNU(MAXDEPTH); i++)
	{
		if (intersectScene(ray, vec2(0.001, MAXFLOAT), rec))
		{
			Ray wi;
			float attenuation;
			float emission;

			bool wasScattered = Material_bsdf(rec, ray, wi, attenuation, emission);

			ray.origin = wi.origin;
			ray.direction = wi.direction;

			if (wasScattered && intensity > 0.01)
				intensity *= attenuation;
			else
			{
				return intensity * emission;
			}
		}
		else
		{
			return intensity * skyColor(ray);
		}
	}

	return 0.0;
}

void main()
{
	vec2 mo = vec2(0.);

	uvec2 p = uvec2(fragCoord);
	seedA = p.x + 1920U * p.y + (1920U * 1080U) * uint(iFrame);
	seedB = p.y + 1920U * p.x + (1920U * 1080U) * uint(iFrame);

	//vec3 lookfrom = vec3(13.0, 2.0, 3.0);
	const vec3 lookat = vec3(0.0, 1.5, 0.0);
	float aperture = 1.0;


	vec2 mo1 = vec2(0.59, 0.53);
	vec2 mo2 = vec2(0.5, 0.415);
	vec2 mo3 = vec2(0.5, 0.895);
	vec2 mo4 = vec2(0.4, 0.4);
	vec2 mo5 = vec2(0.1, 0.5);
	mo = mo1;
	float time = iTime - 25.0 * floor(iTime / 25.0);
	mo = mix(mo, mo2, smoothstep(1.0, 4.0, time));
	mo = mix(mo, mo3, smoothstep(5.0, 10.0, time));
	mo = mix(mo, mo4, smoothstep(12.0, 15.0, time));
	mo = mix(mo, mo5, smoothstep(15.0, 20.0, time));
	mo = mix(mo, mo1, smoothstep(20.0, 25.0, time));

	// camera
	float theta = 6.283 * (0.5 * mo.y - 0.5);
	vec3 ro = vec3(cos(6.283 * mo.x) * sin(theta),
		cos(theta),
		sin(6.283 * mo.x) * sin(theta));
	ro *= 35.0;
	float distToFocus = length(lookat - ro);

	Camera camera;
	Camera_init(camera, ro, lookat, vec3(0.0f, 1.0f, 0.0f), 20.0f, float(iResolution.x) / float(iResolution.y), aperture, distToFocus);

	vec3 col = vec3(0.0, 0.0, 0.0);

	vec2 uv = fragCoord;// / vec2(iResolution);

	IntersectInfo rec;
	for (int s = 0; s < DNU(NUMSAMPLES); s++)
	{
		Ray r = Camera_getRay(camera, uv);
		r.wave_length = rand(400.0, 700.0);

		float intensity = radiance(r);

		// vec3 color = cie_colour_match[r.wave_length] * D65[r.wave_length];
		vec3 color = wavelength2XYZ(r.wave_length);

		col += color * intensity;
	}
	col = XYZ_2_RGB * col;
	col /= float(NUMSAMPLES);
	col /= 40.0;
	col = clamp(col, vec3(0.0), vec3(1.0));

	// compress
	col = 1.35 * col / (1.0 + col);

	// gamma
	col = pow(col, vec3(0.4545));

	//// s-surve
	col = clamp(col, 0.0, 1.0);
	col = col * col * (3.0 - 2.0 * col);

	// output
	fragColor = vec4(col, 1.0);
}
