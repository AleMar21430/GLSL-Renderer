#version 460 core

uniform float iTime;
uniform int iFrame;
uniform vec2 iResolution;

uniform sampler2D iLastFrame;

in vec2 fragCoord;
in vec2 fragTexCoord;

out vec4 fragColor;



#define PI 3.1415926535


struct Circle
{
	vec2 ctr;
	float rad;
};

struct Ray
{
	vec2 org;
	vec2 dir;
};

struct Segment
{
	vec2 A;
	vec2 B;
};

struct Bezier3
{
	vec2 A;
	vec2 B;
	vec2 C;
};

struct Bezier4
{
	vec2 A;
	vec2 B;
	vec2 C;
	vec2 D;
};

struct Grid
{
	float inter;
};

bool intersect(Circle c, Ray r, out vec2 t)
{
	float A = dot(r.dir, r.dir);
	float B = dot(r.dir, -c.ctr + r.org);
	float C = dot(c.ctr, c.ctr) + dot(r.org, r.org) - 2. * dot(c.ctr, r.org) - c.rad * c.rad;

	float delta = B * B - A * C;

	if (delta < 0.0)
		return false;

	t = (vec2(-B) + vec2(sqrt(delta)) * vec2(-1., 1.)) / A;
	return true;
}

vec2 getPt(Ray r, float t)
{
	return r.org + t * r.dir;
}

vec2 getPt(Circle c, float t)
{
	return c.ctr + c.rad * vec2(cos(t * PI * 2.), sin(t * PI * 2.));
}

vec2 getPt(Segment c, float t)
{
	return mix(c.A, c.B, t);
}

vec2 getPt(Bezier3 c, float t)
{
	return mix(mix(c.A, c.B, t), mix(c.B, c.C, t), t);
}

vec2 getPt(Bezier4 c, float t)
{
	return mix(mix(mix(c.A, c.B, t), mix(c.B, c.C, t), t), mix(mix(c.B, c.C, t), mix(c.C, c.D, t), t), t);
}

float dist(Ray r, vec2 p)
{
	return abs(dot(r.dir.yx * vec2(-1., 1.), p - r.org) / dot(r.dir, r.dir));
}

float dist(vec2 p0, vec2 p1)
{
	return length(p1 - p0);
}

float dist(Circle c, vec2 p)
{
	return abs(length(c.ctr - p) - c.rad);
}

float dist(Segment s, vec2 p)
{
	vec2 pa = p - s.A, ba = s.B - s.A;
	float h = clamp(dot(pa, ba) / dot(ba, ba), 0.0, 1.0);
	return length(pa - ba * h);

}

float dist(Grid g, vec2 p)
{
	vec2 d = mod(p, g.inter);
	d = min(d, g.inter - d);
	return min(d.x, d.y);
}


// Solve cubic equation for roots
vec3 solveCubic(float a, float b, float c)
{
	float p = b - a * a / 3.0, p3 = p * p * p;
	float q = a * (2.0 * a * a - 9.0 * b) / 27.0 + c;
	float d = q * q + 4.0 * p3 / 27.0;
	float offset = -a / 3.0;
	if (d >= 0.0) {
		float z = sqrt(d);
		vec2 x = (vec2(z, -z) - q) / 2.0;
		vec2 uv = sign(x) * pow(abs(x), vec2(1.0 / 3.0));
		return vec3(offset + uv.x + uv.y);
	}
	float v = acos(-sqrt(-27.0 / p3) * q / 2.0) / 3.0;
	float m = cos(v), n = sin(v) * 1.732050808;
	return vec3(m + m, -n - m, n - m) * sqrt(-p / 3.0) + offset;
}

// Find the signed distance from a point to a bezier curve
float dist(Bezier3 B, vec2 p)
{
	vec2 a = B.B - B.A, b = B.A - B.B * 2.0 + B.C, c = a * 2.0, d = B.A - p;
	vec3 k = vec3(3. * dot(a, b), 2. * dot(a, a) + dot(d, b), dot(d, a)) / dot(b, b);
	vec3 t = clamp(solveCubic(k.x, k.y, k.z), 0.0, 1.0);
	vec2 pos = B.A + (c + b * t.x) * t.x;
	float dis = length(pos - p);
	pos = B.A + (c + b * t.y) * t.y;
	dis = min(dis, length(pos - p));
	pos = B.A + (c + b * t.z) * t.z;
	dis = min(dis, length(pos - p));
	return dis;
}

float dist(Bezier4 BB, vec2 p)
{
	vec2 A = BB.A, B = BB.D;

	float ppt;
	float At = .0;
	float Bt = 1.;

	vec2 pp;

	float dis = dist(pp, p);

	for (int i = 0; i < 5; ++i)
	{
		ppt = (At + Bt) * .5;
		pp = getPt(BB, ppt);

		if (dist(Segment(A, pp), p) < dist(Segment(pp, B), p))
		{
			Bt = ppt;
			B = getPt(BB, Bt);
		}
		else
		{
			At = ppt;
			A = getPt(BB, At);
		}

	}

	return min(dist(Segment(A, pp), p), dist(Segment(B, pp), p));
}


#define DRAW(O,P,C,W, CC)	CC=mix(CC,C,mix(1.,0.,clamp(dist(O,P)*iResolution.y/2. -  W, -1., 1.)*.5+.5))




#define lineWidth 	(0.)

vec3 hsv2rgb_smooth(in vec3 c)
{
	vec3 rgb = clamp(abs(mod(c.x * 6.0 + vec3(0.0, 4.0, 2.0), 6.0) - 3.0) - 1.0, 0.0, 1.0);

	rgb = rgb * rgb * (3.0 - 2.0 * rgb); // cubic smoothing	

	return c.z * mix(vec3(1.0), rgb, c.y);
}

void main()
{

	int N = int(iResolution.y);
	float dy = 1. / (float(N) - 1.);
	//float k = float(iFrame) * dy;
	float k = float(iFrame) / 1000.;
	k += float(iFrame / 1000) / 3000.;

	vec2 p = vec2(fragCoord.x - 0.5, fragCoord.y - 0.5) * 2.5;
	vec3 col = vec3(0.);

	vec2 p0 = vec2(1., .0);
	vec2 p1 = vec2(.5, sin(k * PI * 2. * 5.));
	vec2 p2 = vec2(-.5, sin(k * PI * 2. * 3.) * .2);
	//vec2 p3 = vec2(-1. ,.0);
	vec2 p3 = vec2(-1., cos(k * PI * 2. * 7.) * .7);

	Bezier4 b = Bezier4(p0, p1, p2, p3);


	vec3 cc = hsv2rgb_smooth(vec3(k, 1., 1.));
	//vec3 cc = vec3(1.);
	DRAW(b, p, cc, lineWidth, col);

	col *= 40.;

	/*if (iFrame < 1) {
		fragColor = vec4(col, 1.);
	}
	else {*/

	vec2 pacc = vec2(fragTexCoord);
	fragColor = (texture(iLastFrame, pacc) * float(iFrame) + vec4(0.1, 0.1, 0.1, 1.0) + vec4(col, 1.) * 15.0) / float(iFrame + 1);
	//fragColor = vec4(col, 1.);
	//}
}