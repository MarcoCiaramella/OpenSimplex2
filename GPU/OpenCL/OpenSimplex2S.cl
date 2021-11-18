#define PMASK 2047




typedef struct {
    int xsv, ysv;
	double dx, dy;
} LatticePoint2D;

typedef struct {
    double dxr, dyr, dzr;
	int xrv, yrv, zrv;
    int nextOnFailure;
    int nextOnSuccess;
} LatticePoint3D;

typedef struct {
    int xsv, ysv, zsv, wsv;
	double dx, dy, dz, dw;
} LatticePoint4D;

typedef struct {
    double dx, dy;
} Grad2;

typedef struct {
    double dx, dy, dz;
} Grad3;

typedef struct {
    double dx, dy, dz, dw;
} Grad4;




int get_index(const unsigned int size, const unsigned int num_dimensions){
	int index = get_global_id(0);
	if (index*num_dimensions + num_dimensions < size){
		return index;
	}
	return -1;
}

double get_x(double* buffer, int index, const unsigned int num_dimensions){
	return buffer[index*num_dimensions];
}

double get_y(double* buffer, int index, const unsigned int num_dimensions){
	return buffer[index*num_dimensions + 1];
}

double get_z(double* buffer, int index, const unsigned int num_dimensions){
	return buffer[index*num_dimensions + 2];
}

double get_w(double* buffer, int index, const unsigned int num_dimensions){
	return buffer[index*num_dimensions + 3];
}

int fast_floor(double x){
	int xi = (int)x;
	return x < xi ? xi - 1 : xi;
}

double _noise2_Base(__global short* perm, __global Grad2* permGrad2, __global LatticePoint2D* LOOKUP_2D, double xs, double ys){
	double value = 0;

	// Get base points and offsets
	int xsb = fast_floor(xs);
	int ysb = fast_floor(ys);
	double xsi = xs - xsb;
	double ysi = ys - ysb;

	// Index to point list
	int a = (int)(xsi + ysi);
	int index =
		(a << 2) |
		(int)(xsi - ysi / 2 + 1 - a / 2.0) << 3 |
		(int)(ysi - xsi / 2 + 1 - a / 2.0) << 4;

	double ssi = (xsi + ysi) * -0.211324865405187;
	double xi = xsi + ssi;
	double yi = ysi + ssi;

	// Point contributions
	for (int i = 0; i < 4; i++){
		__global LatticePoint2D *c = &(LOOKUP_2D[index + i]);

		double dx = xi + c->dx;
		double dy = yi + c->dy;
		double attn = 2.0 / 3.0 - dx * dx - dy * dy;
		if (attn <= 0)
			continue;

		int pxm = (xsb + c->xsv) & PMASK;
		int pym = (ysb + c->ysv) & PMASK;
		Grad2 grad = permGrad2[perm[pxm] ^ pym];
		double extrapolation = grad.dx * dx + grad.dy * dy;

		attn *= attn;
		value += attn * attn * extrapolation;
	}

	return value;
}

__kernel void noise2(
	__global short* perm,
    __global Grad2* permGrad2,
    __global LatticePoint2D* LOOKUP_2D,
	__global double* input,
	const unsigned int size,
	__global double* output){

	int index = get_index(size, 2);
    if (index >= 0){

		double x = get_x(input, index, 2);
		double y = get_y(input, index, 2);

		// Get points for A2* lattice
		double s = 0.366025403784439 * (x + y);
		double xs = x + s;
		double ys = y + s;

		output[index] = _noise2_Base(perm, permGrad2, LOOKUP_2D, xs, ys);
	}
}

__kernel void noise2_XBeforeY(
	__global short* perm,
    __global Grad2* permGrad2,
    __global LatticePoint2D* LOOKUP_2D,
	__global double* input,
	const unsigned int size,
	__global double* output){

	int index = get_index(size, 2);
    if (index >= 0){

		double x = get_x(input, index, 2);
		double y = get_y(input, index, 2);

		// Skew transform and rotation baked into one.
		double xx = x * 0.7071067811865476;
		double yy = y * 1.224744871380249;

		output[index] = _noise2_Base(perm, permGrad2, LOOKUP_2D, yy + xx, yy - xx);
	}
}

double _noise3_BCC(__global short* perm, __global Grad3* permGrad3, __global LatticePoint3D* LOOKUP_3D, double xr, double yr, double zr){

	// Get base and offsets inside cube of first lattice.
	int xrb = fast_floor(xr);
	int yrb = fast_floor(yr);
	int zrb = fast_floor(zr);
	double xri = xr - xrb;
	double yri = yr - yrb;
	double zri = zr - zrb;

	// Identify which octant of the cube we're in. This determines which cell
	// in the other cubic lattice we're in, and also narrows down one point on each.
	int xht = (int)(xri + 0.5);
	int yht = (int)(yri + 0.5);
	int zht = (int)(zri + 0.5);
	int index = (xht << 0) | (yht << 1) | (zht << 2);

	// Point contributions
	double value = 0;
	while (index >= 0){
		__global LatticePoint3D *c = &(LOOKUP_3D[index]);
		double dxr = xri + c->dxr;
		double dyr = yri + c->dyr;
		double dzr = zri + c->dzr;
		double attn = 0.75 - dxr * dxr - dyr * dyr - dzr * dzr;
		if (attn < 0){
			index = c->nextOnFailure;
		}
		else{
			int pxm = (xrb + c->xrv) & PMASK;
			int pym = (yrb + c->yrv) & PMASK;
			int pzm = (zrb + c->zrv) & PMASK;
			Grad3 grad = permGrad3[perm[perm[pxm] ^ pym] ^ pzm];
			double extrapolation = grad.dx * dxr + grad.dy * dyr + grad.dz * dzr;

			attn *= attn;
			value += attn * attn * extrapolation;
			index = c->nextOnSuccess;
		}
	}
	return value;
}

__kernel void noise3_Classic(
	__global short* perm,
    __global Grad3* permGrad3,
    __global LatticePoint3D* LOOKUP_3D,
	__global double* input,
	const unsigned int size,
	__global double* output){

	int index = get_index(size, 3);
    if (index >= 0){
		
		double x = get_x(input, index, 3);
		double y = get_y(input, index, 3);
		double z = get_z(input, index, 3);

		// Re-orient the cubic lattices via rotation, to produce the expected look on cardinal planar slices.
		// If texturing objects that don't tend to have cardinal plane faces, you could even remove this.
		// Orthonormal rotation. Not a skew transform.
		double r = (2.0 / 3.0) * (x + y + z);
		double xr = r - x;
		double yr = r - y;
		double zr = r - z;

		// Evaluate both lattices to form a BCC lattice.
		output[index] = _noise3_BCC(perm, permGrad3, LOOKUP_3D, xr, yr, zr);
	}
}

__kernel void noise3_XYBeforeZ(
	__global short* perm,
    __global Grad3* permGrad3,
    __global LatticePoint3D* LOOKUP_3D,
	__global double* input,
	const unsigned int size,
	__global double* output){

	int index = get_index(size, 3);
    if (index >= 0){
		
		double x = get_x(input, index, 3);
		double y = get_y(input, index, 3);
		double z = get_z(input, index, 3);

		// Re-orient the cubic lattices without skewing, to make X and Y triangular like 2D.
		// Orthonormal rotation. Not a skew transform.
		double xy = x + y;
		double s2 = xy * -0.211324865405187;
		double zz = z * 0.577350269189626;
		double xr = x + s2 - zz;
		double yr = y + s2 - zz;
		double zr = xy * 0.577350269189626 + zz;

		// Evaluate both lattices to form a BCC lattice.
		output[index] = _noise3_BCC(perm, permGrad3, LOOKUP_3D, xr, yr, zr);
	}
}

__kernel void noise3_XZBeforeY(
	__global short* perm,
    __global Grad3* permGrad3,
    __global LatticePoint3D* LOOKUP_3D,
	__global double* input,
	const unsigned int size,
	__global double* output){

	int index = get_index(size, 3);
    if (index >= 0){
		
		double x = get_x(input, index, 3);
		double y = get_y(input, index, 3);
		double z = get_z(input, index, 3);

		// Re-orient the cubic lattices without skewing, to make X and Z triangular like 2D.
		// Orthonormal rotation. Not a skew transform.
		double xz = x + z;
		double s2 = xz * -0.211324865405187;
		double yy = y * 0.577350269189626;
		double xr = x + s2 - yy;
		double zr = z + s2 - yy;
		double yr = xz * 0.577350269189626 + yy;

		// Evaluate both lattices to form a BCC lattice.
		output[index] = _noise3_BCC(perm, permGrad3, LOOKUP_3D, xr, yr, zr);
	}
}

double _noise4_Base(__global short* perm, __global Grad4* permGrad4, __global LatticePoint4D* LOOKUP_4D, double xs, double ys, double zs, double ws){
	const unsigned int sizes[256] = {
		20, 15, 16, 17, 15, 16, 12, 15, 16, 12, 10, 14, 17, 15, 14,
		17, 15, 16, 12, 15, 16, 14, 14, 13, 12, 14, 11, 12, 15, 13,
		12, 14, 16, 12, 10, 14, 12, 14, 11, 12, 10, 11, 10, 13, 14,
		12, 13, 15, 17, 15, 14, 17, 15, 13, 12, 14, 14, 12, 13, 15,
		17, 14, 15, 17, 15, 16, 12, 15, 16, 14, 14, 13, 12, 14, 11,
		12, 15, 13, 12, 14, 16, 14, 14, 13, 14, 16, 16, 10, 14, 16,
		19, 11, 13, 10, 11, 10, 12, 14, 11, 12, 14, 16, 19, 11, 11,
		19, 13, 14, 12, 11, 14, 13, 15, 13, 12, 14, 13, 10, 11, 10,
		12, 11, 14, 13, 14, 10, 13, 13, 16, 12, 10, 14, 12, 14, 11,
		12, 10, 11, 10, 13, 14, 12, 13, 15, 12, 14, 11, 12, 14, 16,
		19, 11, 11, 19, 13, 14, 12, 11, 14, 13, 10, 11, 10, 13, 11,
		19, 13, 14, 10, 13, 16, 14, 13, 14, 14, 16, 14, 12, 13, 15,
		12, 11, 14, 13, 13, 14, 14, 16, 15, 13, 16, 15, 17, 15, 14,
		17, 15, 13, 12, 14, 14, 12, 13, 15, 17, 14, 15, 17, 15, 13,
		12, 14, 13, 10, 11, 10, 12, 11, 14, 13, 14, 10, 13, 13, 14,
		12, 13, 15, 12, 11, 14, 13, 13, 14, 14, 16, 15, 13, 16, 15,
		17, 14, 15, 17, 14, 10, 13, 13, 15, 13, 16, 15, 17, 13, 15, 20};

	const unsigned int offsets[256] = {0, 20, 35, 51, 68, 83, 99, 111, 126, 142, 154, 164, 178,
	195, 210, 224, 241, 256, 272, 284, 299, 315, 329, 343, 356, 368,
	382, 393, 405, 420, 433, 445, 459, 475, 487, 497, 511, 523, 537,
	548, 560, 570, 581, 591, 604, 618, 630, 643, 658, 675, 690, 704,
	721, 736, 749, 761, 775, 789, 801, 814, 829, 846, 860, 875, 892,
	907, 923, 935, 950, 966, 980, 994, 1007, 1019, 1033, 1044, 1056,
	1071, 1084, 1096, 1110, 1126, 1140, 1154, 1167, 1181, 1197, 1213,
	1223, 1237, 1253, 1272, 1283, 1296, 1306, 1317, 1327, 1339, 1353,
	1364, 1376, 1390, 1406, 1425, 1436, 1447, 1466, 1479, 1493, 1505,
	1516, 1530, 1543, 1558, 1571, 1583, 1597, 1610, 1620, 1631, 1641,
	1653, 1664, 1678, 1691, 1705, 1715, 1728, 1741, 1757, 1769, 1779,
	1793, 1805, 1819, 1830, 1842, 1852, 1863, 1873, 1886, 1900, 1912,
	1925, 1940, 1952, 1966, 1977, 1989, 2003, 2019, 2038, 2049, 2060,
	2079, 2092, 2106, 2118, 2129, 2143, 2156, 2166, 2177, 2187, 2200,
	2211, 2230, 2243, 2257, 2267, 2280, 2296, 2310, 2323, 2337, 2351,
	2367, 2381, 2393, 2406, 2421, 2433, 2444, 2458, 2471, 2484, 2498,
	2512, 2528, 2543, 2556, 2572, 2587, 2604, 2619, 2633, 2650, 2665,
	2678, 2690, 2704, 2718, 2730, 2743, 2758, 2775, 2789, 2804, 2821,
	2836, 2849, 2861, 2875, 2888, 2898, 2909, 2919, 2931, 2942, 2956,
	2969, 2983, 2993, 3006, 3019, 3033, 3045, 3058, 3073, 3085, 3096,
	3110, 3123, 3136, 3150, 3164, 3180, 3195, 3208, 3224, 3239, 3256,
	3270, 3285, 3302, 3316, 3326, 3339, 3352, 3367, 3380, 3396, 3411,
	3428, 3441, 3456};

	double value = 0;

	// Get base points and offsets
	int xsb = fast_floor(xs);
	int ysb = fast_floor(ys);
	int zsb = fast_floor(zs);
	int wsb = fast_floor(ws);
	double xsi = xs - xsb;
	double ysi = ys - ysb;
	double zsi = zs - zsb;
	double wsi = ws - wsb;

	// Unskewed offsets
	double ssi = (xsi + ysi + zsi + wsi) * -0.138196601125011;
	double xi = xsi + ssi;
	double yi = ysi + ssi;
	double zi = zsi + ssi;
	double wi = wsi + ssi;

	int index = ((fast_floor(xs * 4) & 3) << 0) | ((fast_floor(ys * 4) & 3) << 2) | ((fast_floor(zs * 4) & 3) << 4) | ((fast_floor(ws * 4) & 3) << 6);

	// Point contributions
	__global LatticePoint4D* c = &(LOOKUP_4D[offsets[index]]);
	for (int i = 0; i < sizes[index]; i++){
		double dx = xi + c->dx;
		double dy = yi + c->dy;
		double dz = zi + c->dz;
		double dw = wi + c->dw;
		double attn = 0.8 - dx * dx - dy * dy - dz * dz - dw * dw;
		if (attn > 0){
			attn *= attn;

			int pxm = (xsb + c->xsv) & PMASK;
			int pym = (ysb + c->ysv) & PMASK;
			int pzm = (zsb + c->zsv) & PMASK;
			int pwm = (wsb + c->wsv) & PMASK;
			Grad4 grad = permGrad4[perm[perm[perm[pxm] ^ pym] ^ pzm] ^ pwm];
			double extrapolation = grad.dx * dx + grad.dy * dy + grad.dz * dz + grad.dw * dw;

			value += attn * attn * extrapolation;
		}
		c++;
	}
	return value;
}

__kernel void noise4_Classic(
	__global short* perm,
    __global Grad4* permGrad4,
    __global LatticePoint4D* LOOKUP_4D,
	__global double* input,
	const unsigned int size,
	__global double* output){

	int index = get_index(size, 4);
    if (index >= 0){
		
		double x = get_x(input, index, 4);
		double y = get_y(input, index, 4);
		double z = get_z(input, index, 4);
		double w = get_w(input, index, 4);

		// Get points for A4 lattice
		double s = 0.309016994374947 * (x + y + z + w);
		double xs = x + s;
		double ys = y + s;
		double zs = z + s;
		double ws = w + s;

		output[index] = _noise4_Base(perm, permGrad4, LOOKUP_4D, xs, ys, zs, ws);
	}
}

__kernel void noise4_XYBeforeZW(
	__global short* perm,
    __global Grad4* permGrad4,
    __global LatticePoint4D* LOOKUP_4D,
	__global double* input,
	const unsigned int size,
	__global double* output){

	int index = get_index(size, 4);
    if (index >= 0){
		
		double x = get_x(input, index, 4);
		double y = get_y(input, index, 4);
		double z = get_z(input, index, 4);
		double w = get_w(input, index, 4);

		double s2 = (x + y) * -0.28522513987434876941 + (z + w) * 0.83897065470611435718;
		double t2 = (z + w) * 0.21939749883706435719 + (x + y) * -0.48214856493302476942;
		double xs = x + s2;
		double ys = y + s2;
		double zs = z + t2;
		double ws = w + t2;

		output[index] = _noise4_Base(perm, permGrad4, LOOKUP_4D, xs, ys, zs, ws);
	}
}

__kernel void noise4_XZBeforeYW(
	__global short* perm,
    __global Grad4* permGrad4,
    __global LatticePoint4D* LOOKUP_4D,
	__global double* input,
	const unsigned int size,
	__global double* output){

	int index = get_index(size, 4);
    if (index >= 0){
		
		double x = get_x(input, index, 4);
		double y = get_y(input, index, 4);
		double z = get_z(input, index, 4);
		double w = get_w(input, index, 4);

		double s2 = (x + z) * -0.28522513987434876941 + (y + w) * 0.83897065470611435718;
		double t2 = (y + w) * 0.21939749883706435719 + (x + z) * -0.48214856493302476942;
		double xs = x + s2;
		double ys = y + t2;
		double zs = z + s2;
		double ws = w + t2;

		output[index] = _noise4_Base(perm, permGrad4, LOOKUP_4D, xs, ys, zs, ws);
	}
}

__kernel void noise4_XYZBeforeW(
	__global short* perm,
    __global Grad4* permGrad4,
    __global LatticePoint4D* LOOKUP_4D,
	__global double* input,
	const unsigned int size,
	__global double* output){

	int index = get_index(size, 4);
    if (index >= 0){
		
		double x = get_x(input, index, 4);
		double y = get_y(input, index, 4);
		double z = get_z(input, index, 4);
		double w = get_w(input, index, 4);

		double xyz = x + y + z;
		double ww = w * 1.118033988749894;
		double s2 = xyz * -0.16666666666666666 + ww;
		double xs = x + s2;
		double ys = y + s2;
		double zs = z + s2;
		double ws = -0.5 * xyz + ww;

		output[index] = _noise4_Base(perm, permGrad4, LOOKUP_4D, xs, ys, zs, ws);
	}
}
