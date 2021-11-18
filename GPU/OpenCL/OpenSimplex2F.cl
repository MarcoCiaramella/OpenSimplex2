#define PMASK 2047
#define INVALID_INDEX -1


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
    double xsi, ysi, zsi, wsi;
    double ssiDelta;
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




int get_point_index_1D(){
	return get_global_id(0);
}

int get_point_index_2D(float width){
	return get_global_id(1)*width + get_global_id(0);
}

int get_point_index_3D(float width, float height){
	return get_global_id(2)*width*height + get_global_id(1)*width + get_global_id(0);
}

int get_point_index(const uint num_points){
	int index = INVALID_INDEX;
	if (get_work_dim() == 1){
		index = get_point_index_1D();
	}
	else if (get_work_dim() == 2){
		index = get_point_index_2D(ceil(sqrt((float)num_points)));
	}
	else if (get_work_dim() == 3){
		index = get_point_index_3D(ceil(sqrt((float)num_points)), ceil(sqrt((float)num_points)));
	}
	if (index < num_points){
		return index;
	}
	return INVALID_INDEX;
}

double get_x(double* buffer, int index, const uint num_dimensions){
	return buffer[index*num_dimensions];
}

double get_y(double* buffer, int index, const uint num_dimensions){
	return buffer[index*num_dimensions + 1];
}

double get_z(double* buffer, int index, const uint num_dimensions){
	return buffer[index*num_dimensions + 2];
}

double get_w(double* buffer, int index, const uint num_dimensions){
	return buffer[index*num_dimensions + 3];
}

int fast_floor(double x){
	int xi = (int)x;
	return x < xi ? xi - 1 : xi;
}

int is_a_valid_index(int index){
	return index != INVALID_INDEX;
}

double _noise2_Base(__global short* perm, __global Grad2* permGrad2, __global LatticePoint2D* LOOKUP_2D, double xs, double ys){
	double value = 0;

	// Get base points and offsets
	int xsb = fast_floor(xs), ysb = fast_floor(ys);
	double xsi = xs - xsb, ysi = ys - ysb;

	// Index to point list
	int index = (int)((ysi - xsi) / 2 + 1);

	double ssi = (xsi + ysi) * -0.211324865405187;
	double xi = xsi + ssi, yi = ysi + ssi;

	// Point contributions
	for (int i = 0; i < 3; i++){
		__global LatticePoint2D *c = &(LOOKUP_2D[index + i]);

		double dx = xi + c->dx, dy = yi + c->dy;
		double attn = 0.5 - dx * dx - dy * dy;
		if (attn <= 0)
			continue;

		int pxm = (xsb + c->xsv) & PMASK, pym = (ysb + c->ysv) & PMASK;
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
	const uint num_points,
	__global double* output){

	int index = get_point_index(num_points);
    if (is_a_valid_index(index)){

		double x = get_x(input, index, 2);
		double y = get_y(input, index, 2);

		// Get points for A2* lattice
		double s = 0.366025403784439 * (x + y);
		double xs = x + s, ys = y + s;

		output[index] = _noise2_Base(perm, permGrad2, LOOKUP_2D, xs, ys);
	}
}

__kernel void noise2_XBeforeY(
	__global short* perm,
    __global Grad2* permGrad2,
    __global LatticePoint2D* LOOKUP_2D,
	__global double* input,
	const uint num_points,
	__global double* output){

	int index = get_point_index(num_points);
    if (is_a_valid_index(index)){

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
	int xrb = fast_floor(xr), yrb = fast_floor(yr), zrb = fast_floor(zr);
	double xri = xr - xrb, yri = yr - yrb, zri = zr - zrb;

	// Identify which octant of the cube we're in. This determines which cell
	// in the other cubic lattice we're in, and also narrows down one point on each.
	int xht = (int)(xri + 0.5), yht = (int)(yri + 0.5), zht = (int)(zri + 0.5);
	int index = (xht << 0) | (yht << 1) | (zht << 2);

	// Point contributions
	double value = 0;
	while (index >= 0){
		__global LatticePoint3D *c = &(LOOKUP_3D[index]);
		double dxr = xri + c->dxr, dyr = yri + c->dyr, dzr = zri + c->dzr;
		double attn = 0.5 - dxr * dxr - dyr * dyr - dzr * dzr;
		if (attn < 0){
			index = c->nextOnFailure;
		}
		else{
			int pxm = (xrb + c->xrv) & PMASK, pym = (yrb + c->yrv) & PMASK, pzm = (zrb + c->zrv) & PMASK;
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
	const uint num_points,
	__global double* output){

	int index = get_point_index(num_points);
    if (is_a_valid_index(index)){
		
		double x = get_x(input, index, 3);
		double y = get_y(input, index, 3);
		double z = get_z(input, index, 3);

		// Re-orient the cubic lattices via rotation, to produce the expected look on cardinal planar slices.
		// If texturing objects that don't tend to have cardinal plane faces, you could even remove this.
		// Orthonormal rotation. Not a skew transform.
		double r = (2.0 / 3.0) * (x + y + z);
		double xr = r - x, yr = r - y, zr = r - z;

		// Evaluate both lattices to form a BCC lattice.
		output[index] = _noise3_BCC(perm, permGrad3, LOOKUP_3D, xr, yr, zr);
	}
}

__kernel void noise3_XYBeforeZ(
	__global short* perm,
    __global Grad3* permGrad3,
    __global LatticePoint3D* LOOKUP_3D,
	__global double* input,
	const uint num_points,
	__global double* output){

	int index = get_point_index(num_points);
    if (is_a_valid_index(index)){
		
		double x = get_x(input, index, 3);
		double y = get_y(input, index, 3);
		double z = get_z(input, index, 3);

		// Re-orient the cubic lattices without skewing, to make X and Y triangular like 2D.
		// Orthonormal rotation. Not a skew transform.
		double xy = x + y;
		double s2 = xy * -0.211324865405187;
		double zz = z * 0.577350269189626;
		double xr = x + s2 - zz, yr = y + s2 - zz;
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
	const uint num_points,
	__global double* output){

	int index = get_point_index(num_points);
    if (is_a_valid_index(index)){
		
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

double _noise4_Base(__global short* perm, __global Grad4* permGrad4, __global LatticePoint4D* VERTICES_4D, double xs, double ys, double zs, double ws){
	double value = 0;

	// Get base points and offsets
	int xsb = fast_floor(xs), ysb = fast_floor(ys), zsb = fast_floor(zs), wsb = fast_floor(ws);
	double xsi = xs - xsb, ysi = ys - ysb, zsi = zs - zsb, wsi = ws - wsb;

	// If we're in the lower half, flip so we can repeat the code for the upper half. We'll flip back later.
	double siSum = xsi + ysi + zsi + wsi;
	double ssi = siSum * 0.309016994374947; // Prep for vertex contributions.
	int inLowerHalf = (siSum < 2);
	if (inLowerHalf){
		xsi = 1 - xsi;
		ysi = 1 - ysi;
		zsi = 1 - zsi;
		wsi = 1 - wsi;
		siSum = 4 - siSum;
	}

	// Consider opposing vertex pairs of the octahedron formed by the central cross-section of the stretched tesseract
	double aabb = xsi + ysi - zsi - wsi, abab = xsi - ysi + zsi - wsi, abba = xsi - ysi - zsi + wsi;
	double aabbScore = fabs(aabb), ababScore = fabs(abab), abbaScore = fabs(abba);

	// Find the closest point on the stretched tesseract as if it were the upper half
	int vertexIndex, via, vib;
	double asi, bsi;
	if (aabbScore > ababScore && aabbScore > abbaScore){
		if (aabb > 0){
			asi = zsi;
			bsi = wsi;
			vertexIndex = 0b0011;
			via = 0b0111;
			vib = 0b1011;
		}
		else{
			asi = xsi;
			bsi = ysi;
			vertexIndex = 0b1100;
			via = 0b1101;
			vib = 0b1110;
		}
	}
	else if (ababScore > abbaScore){
		if (abab > 0){
			asi = ysi;
			bsi = wsi;
			vertexIndex = 0b0101;
			via = 0b0111;
			vib = 0b1101;
		}
		else{
			asi = xsi;
			bsi = zsi;
			vertexIndex = 0b1010;
			via = 0b1011;
			vib = 0b1110;
		}
	}
	else{
		if (abba > 0){
			asi = ysi;
			bsi = zsi;
			vertexIndex = 0b1001;
			via = 0b1011;
			vib = 0b1101;
		}
		else{
			asi = xsi;
			bsi = wsi;
			vertexIndex = 0b0110;
			via = 0b0111;
			vib = 0b1110;
		}
	}
	if (bsi > asi){
		via = vib;
		double temp = bsi;
		bsi = asi;
		asi = temp;
	}
	if (siSum + asi > 3){
		vertexIndex = via;
		if (siSum + bsi > 4){
			vertexIndex = 0b1111;
		}
	}

	// Now flip back if we're actually in the lower half.
	if (inLowerHalf){
		xsi = 1 - xsi;
		ysi = 1 - ysi;
		zsi = 1 - zsi;
		wsi = 1 - wsi;
		vertexIndex ^= 0b1111;
	}

	// Five points to add, total, from five copies of the A4 lattice.
	for (int i = 0; i < 5; i++){

		// Update xsb/etc. and add the lattice point's contribution.
		__global LatticePoint4D *c = &(VERTICES_4D[vertexIndex]);
		xsb += c->xsv;
		ysb += c->ysv;
		zsb += c->zsv;
		wsb += c->wsv;
		double xi = xsi + ssi, yi = ysi + ssi, zi = zsi + ssi, wi = wsi + ssi;
		double dx = xi + c->dx, dy = yi + c->dy, dz = zi + c->dz, dw = wi + c->dw;
		double attn = 0.5 - dx * dx - dy * dy - dz * dz - dw * dw;
		if (attn > 0){
			int pxm = xsb & PMASK, pym = ysb & PMASK, pzm = zsb & PMASK, pwm = wsb & PMASK;
			Grad4 grad = permGrad4[perm[perm[perm[pxm] ^ pym] ^ pzm] ^ pwm];
			double ramped = grad.dx * dx + grad.dy * dy + grad.dz * dz + grad.dw * dw;

			attn *= attn;
			value += attn * attn * ramped;
		}

		// Maybe this helps the compiler/JVM/LLVM/etc. know we can end the loop here. Maybe not.
		if (i == 4)
			break;

		// Update the relative skewed coordinates to reference the vertex we just added.
		// Rather, reference its counterpart on the lattice copy that is shifted down by
		// the vector <-0.2, -0.2, -0.2, -0.2>
		xsi += c->xsi;
		ysi += c->ysi;
		zsi += c->zsi;
		wsi += c->wsi;
		ssi += c->ssiDelta;

		// Next point is the closest vertex on the 4-simplex whose base vertex is the aforementioned vertex.
		double score0 = 1.0 + ssi * (-1.0 / 0.309016994374947); // Seems slightly faster than 1.0-xsi-ysi-zsi-wsi
		vertexIndex = 0b0000;
		if (xsi >= ysi && xsi >= zsi && xsi >= wsi && xsi >= score0){
			vertexIndex = 0b0001;
		}
		else if (ysi > xsi && ysi >= zsi && ysi >= wsi && ysi >= score0){
			vertexIndex = 0b0010;
		}
		else if (zsi > xsi && zsi > ysi && zsi >= wsi && zsi >= score0){
			vertexIndex = 0b0100;
		}
		else if (wsi > xsi && wsi > ysi && wsi > zsi && wsi >= score0){
			vertexIndex = 0b1000;
		}
	}

	return value;
}

__kernel void noise4_Classic(
	__global short* perm,
    __global Grad4* permGrad4,
    __global LatticePoint4D* VERTICES_4D,
	__global double* input,
	const uint num_points,
	__global double* output){

	int index = get_point_index(num_points);
    if (is_a_valid_index(index)){
		
		double x = get_x(input, index, 4);
		double y = get_y(input, index, 4);
		double z = get_z(input, index, 4);
		double w = get_w(input, index, 4);

		// Get points for A4 lattice
		double s = -0.138196601125011 * (x + y + z + w);
		double xs = x + s, ys = y + s, zs = z + s, ws = w + s;

		output[index] = _noise4_Base(perm, permGrad4, VERTICES_4D, xs, ys, zs, ws);
	}
}

__kernel void noise4_XYBeforeZW(
	__global short* perm,
    __global Grad4* permGrad4,
    __global LatticePoint4D* VERTICES_4D,
	__global double* input,
	const uint num_points,
	__global double* output){

	int index = get_point_index(num_points);
    if (is_a_valid_index(index)){
		
		double x = get_x(input, index, 4);
		double y = get_y(input, index, 4);
		double z = get_z(input, index, 4);
		double w = get_w(input, index, 4);

		double s2 = (x + y) * -0.178275657951399372 + (z + w) * 0.215623393288842828;
		double t2 = (z + w) * -0.403949762580207112 + (x + y) * -0.375199083010075342;
		double xs = x + s2, ys = y + s2, zs = z + t2, ws = w + t2;

		output[index] = _noise4_Base(perm, permGrad4, VERTICES_4D, xs, ys, zs, ws);
	}
}

__kernel void noise4_XZBeforeYW(
	__global short* perm,
    __global Grad4* permGrad4,
    __global LatticePoint4D* VERTICES_4D,
	__global double* input,
	const uint num_points,
	__global double* output){

	int index = get_point_index(num_points);
    if (is_a_valid_index(index)){
		
		double x = get_x(input, index, 4);
		double y = get_y(input, index, 4);
		double z = get_z(input, index, 4);
		double w = get_w(input, index, 4);

		double s2 = (x + z) * -0.178275657951399372 + (y + w) * 0.215623393288842828;
		double t2 = (y + w) * -0.403949762580207112 + (x + z) * -0.375199083010075342;
		double xs = x + s2, ys = y + t2, zs = z + s2, ws = w + t2;

		output[index] = _noise4_Base(perm, permGrad4, VERTICES_4D, xs, ys, zs, ws);
	}
}

__kernel void noise4_XYZBeforeW(
	__global short* perm,
    __global Grad4* permGrad4,
    __global LatticePoint4D* VERTICES_4D,
	__global double* input,
	const uint num_points,
	__global double* output){

	int index = get_point_index(num_points);
    if (is_a_valid_index(index)){
		
		double x = get_x(input, index, 4);
		double y = get_y(input, index, 4);
		double z = get_z(input, index, 4);
		double w = get_w(input, index, 4);

		double xyz = x + y + z;
		double ww = w * 0.2236067977499788;
		double s2 = xyz * -0.16666666666666666 + ww;
		double xs = x + s2, ys = y + s2, zs = z + s2, ws = -0.5 * xyz + ww;

		output[index] = _noise4_Base(perm, permGrad4, VERTICES_4D, xs, ys, zs, ws);
	}
}
