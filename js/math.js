
function atan2_xge0 ( y, x ) {
	if (Math.abs(y) > x) {
		return Math.sign(y) * (Math.PI/2.0 - Math.atan2(x, Math.abs(y)));
	}else{
        return Math.atan2(y, x);
    }
}
function atan2_safe( y, x ) {
    if ( x >= 0.0 ) return atan2_xge0( y, x );
    else return (Math.sign(y) * Math.PI ) - atan2_xge0(y, -x);
}
function atan_safe( yx ) {
    if (Math.abs(yx) > 1.0) {
        return Math.sign(yx) * (Math.PI / 2 - Math.atan(1.0/Math.abs(yx)));
    } else {
        return Math.atan(yx);
    }
}

function DEG_TO_RAD( input ){
	return input * Math.PI / 180;
}

THREE.Matrix3.prototype.makeRotationX = function (theta){
	var c = Math.cos( theta ), s = Math.sin( theta );
	this.set(
		1, 0,  0,
		0, c, - s,
		0, s,  c
	);
	return this;
}

THREE.Matrix3.prototype.makeRotationY = function ( theta ) {
	var c = Math.cos( theta ), s = Math.sin( theta );
	this.set(
		c, 0, s,
		0, 1, 0,
		- s, 0, c
	);
	return this;
}

THREE.Matrix3.prototype.makeRotationZ = function ( theta ) {
	var c = Math.cos( theta ), s = Math.sin( theta );
	this.set(
		c, - s, 0,
		s,  c, 0,
		0,  0, 1
	);
	return this;
}

THREE.Matrix3.prototype.multiply = function ( m, n ) {

	if ( n !== undefined ) {

		console.warn( 'THREE.Matrix4: .multiply() now only accepts one argument. Use .multiplyMatrices( a, b ) instead.' );
		return this.multiplyMatrices( m, n );

	}

	return this.multiplyMatrices( this, m );

}

THREE.Matrix3.prototype.multiplyMatrices = function ( a, b ) {

	var ae = a.elements;
	var be = b.elements;
	var te = this.elements;

	var a11 = ae[ 0 ], a12 = ae[ 3 ], a13 = ae[ 6 ];//, a14 = ae[ 12 ];
	var a21 = ae[ 1 ], a22 = ae[ 4 ], a23 = ae[ 7 ];//, a24 = ae[ 13 ];
	var a31 = ae[ 2 ], a32 = ae[ 5 ], a33 = ae[ 8 ];//, a34 = ae[ 14 ];
	//var a41 = ae[ 3 ], a42 = ae[ 7 ], a43 = ae[ 11 ], a44 = ae[ 15 ];

	var b11 = be[ 0 ], b12 = be[ 3 ], b13 = be[ 6 ];//, b14 = be[ 12 ];
	var b21 = be[ 1 ], b22 = be[ 4 ], b23 = be[ 7 ];//, b24 = be[ 13 ];
	var b31 = be[ 2 ], b32 = be[ 5 ], b33 = be[ 8 ];//, b34 = be[ 14 ];
	//var b41 = be[ 3 ], b42 = be[ 7 ], b43 = be[ 11 ]//, b44 = be[ 15 ];

	te[ 0 ] = a11 * b11 + a12 * b21 + a13 * b31;// + a14 * b41;
	te[ 3 ] = a11 * b12 + a12 * b22 + a13 * b32;// + a14 * b42;
	te[ 6 ] = a11 * b13 + a12 * b23 + a13 * b33;// + a14 * b43;
	//te[ 12 ] = a11 * b14 + a12 * b24 + a13 * b34// + a14 * b44;

	te[ 1 ] = a21 * b11 + a22 * b21 + a23 * b31;// + a24 * b41;
	te[ 4 ] = a21 * b12 + a22 * b22 + a23 * b32;// + a24 * b42;
	te[ 7 ] = a21 * b13 + a22 * b23 + a23 * b33;// + a24 * b43;
	//te[ 13 ] = a21 * b14 + a22 * b24 + a23 * b34// + a24 * b44;

	te[ 2 ] = a31 * b11 + a32 * b21 + a33 * b31;// + a34 * b41;
	te[ 5 ] = a31 * b12 + a32 * b22 + a33 * b32;// + a34 * b42;
	te[ 8 ] = a31 * b13 + a32 * b23 + a33 * b33;// + a34 * b43;
	//te[ 14 ] = a31 * b14 + a32 * b24 + a33 * b34// + a34 * b44;

	//te[ 3 ] = a41 * b11 + a42 * b21 + a43 * b31;// + a44 * b41;
	//te[ 7 ] = a41 * b12 + a42 * b22 + a43 * b32;// + a44 * b42;
	//te[ 11 ] = a41 * b13 + a42 * b23 + a43 * b33;// + a44 * b43;
	//te[ 15 ] = a41 * b14 + a42 * b24 + a43 * b34// + a44 * b44;

	return this;

}

/** transforms a vector */
THREE.Matrix3.prototype.transformVector = function ( v ) {
	var result = new THREE.Vector3();
	var _this = this.elements;
	//var m = [3][3];

	//m[0][0] = _this[0]; m[1][0] = _this[1]; m[2][0] = _this[2];
	//m[0][1] = _this[3]; m[1][1] = _this[4]; m[2][1] = _this[5];
	//m[0][2] = _this[6]; m[1][2] = _this[7]; m[2][2] = _this[8];

	result.x = v.x * _this[0] + v.y * _this[1] + v.z * _this[2];
	result.y = v.x * _this[3] + v.y * _this[4] + v.z * _this[5];
	result.z = v.x * _this[6] + v.y * _this[7] + v.z * _this[8];
	return result;

}


function SetMatrix( a, b, c, cl )
{
	var mx = new THREE.Matrix3(), my = new THREE.Matrix3(), mz = new THREE.Matrix3();

	// Calculate Matrices;
	mx.makeRotationX( a );
	my.makeRotationY( b );
	mz.makeRotationZ( c );

	if (cl)
		return   mz.multiply( mx ).multiply( my ) ;
	else
		return   mx.multiply( mz ).multiply( my ) ;

}


// image transforms/ corrections

var MAXITER = 100;
var R_EPS = 1.0e-6;

function cubeZero_copy( a ){

	var out = {};
	out.root = [];

	if( a[3] == 0.0 ){ // second order polynomial
		out = squareZero_copy( a );
	}else{
		var p = ((-1.0/3.0) * (a[2]/a[3]) * (a[2]/a[3]) + a[1]/a[3]) / 3.0;
		var q = ((2.0/27.0) * (a[2]/a[3]) * (a[2]/a[3]) * (a[2]/a[3]) - (1.0/3.0) * (a[2]/a[3]) * (a[1]/a[3]) + a[0]/a[3]) / 2.0;

		if( q*q + p*p*p >= 0.0 ){
			out.n = 1;
			out.root[0] = cubeRoot_copy(-q + Math.sqrt(q*q + p*p*p)) + cubeRoot_copy(-q - Math.sqrt(q*q + p*p*p)) - a[2] / (3.0 * a[3]);
		}else{
			var phi = Math.acos( -q / Math.sqrt(-p*p*p) );
			out.n = 3;
			out.root[0] =  2.0 * Math.sqrt(-p) * Math.cos(phi/3.0) - a[2] / (3.0 * a[3]);
			out.root[1] = -2.0 * Math.sqrt(-p) * Math.cos(phi/3.0 + Math.PI/3.0) - a[2] / (3.0 * a[3]);
			out.root[2] = -2.0 * Math.sqrt(-p) * Math.cos(phi/3.0 - Math.PI/3.0) - a[2] / (3.0 * a[3]);
		}
	}
	// PrintError("%lg, %lg, %lg, %lg root = %lg", a[3], a[2], a[1], a[0], root[0]);

	return out;
}

function squareZero_copy( a ){

	var out = {};
	out.root = [];

	if( a[2] == 0.0 ){ // linear equation
		if( a[1] == 0.0 ){ // constant
			if( a[0] == 0.0 ){
				out.n = 1; out.root[0] = 0.0;
			}else{
				out.n = 0;
			}
		}else{
			out.n = 1; out.root[0] = - a[0] / a[1];
		}
	}else{
		if( 4.0 * a[2] * a[0] > a[1] * a[1] ){
			out.n = 0;
		}else{
			out.n = 2;
			out.root[0] = (- a[1] + Math.sqrt( a[1] * a[1] - 4.0 * a[2] * a[0] )) / (2.0 * a[2]);
			out.root[1] = (- a[1] - Math.sqrt( a[1] * a[1] - 4.0 * a[2] * a[0] )) / (2.0 * a[2]);
		}
	}

	return out;
}

function cubeRoot_copy( x ){
	if( x == 0.0 )
		return 0.0;
	else if( x > 0.0 )
		return Math.pow(x, 1.0/3.0);
	else
		return - Math.pow(-x, 1.0/3.0);
}

function smallestRoot_copy( p ){
	var n,i;
	var root, sroot = 1000.0;

	var copy = cubeZero_copy( p );

	root = copy.root;
	n = copy.n;

	for( i = 0; i < n; i++ ){
		// PrintError("Root %d = %lg", i,root[i]);
		if ( root[i] > 0.0 && root[i] < sroot )
			sroot = root[i];
	}

	// PrintError("Smallest Root  = %lg", sroot);
	return sroot;
}

// Restrict radial correction to monotonous interval
function CalcCorrectionRadius_copy( coeff ) {
	var a = [4];
	var k;

	for ( k = 0; k < 4; k++ ){
		a[k] = 0.0;//1.0e-10;
		if ( coeff[k] != 0.0 )
		{
			a[k] = ( k + 1 ) * coeff[k];
		}
	}
	return smallestRoot_copy( a );
}

// Calculate inverse 4th order polynomial correction using Newton
function inv_radial( x_dest, y_dest, mprad ) {
	// params: double coefficients[5]
	var params = {};
	params.var0 = mprad[0];
	params.var1 = mprad[1];
	params.var2 = mprad[2];
	params.var3 = mprad[3];
	params.var4 = mprad[4];

	var rs, rd, f, scale;
	var iter = 0;

	rd	= ( Math.sqrt( x_dest*x_dest + y_dest*y_dest )) / params.var4; // Normalized

	rs	= rd;
	f 	= (((params.var3 * rs + params.var2) * rs + params.var1) * rs + params.var0) * rs;

	while( Math.abs(f - rd) > R_EPS && iter++ < MAXITER )
	{
		rs = rs - (f - rd) / ((( 4 * params.var3 * rs + 3 * params.var2) * rs  +
			2 * params.var1) * rs + params.var0);

		f 	= (((params.var3 * rs + params.var2) * rs +
			params.var1) * rs + params.var0) * rs;
	}

	scale = rs / rd;

	var x_src = x_dest * scale ;
	var y_src = y_dest * scale ;

	return { x: x_src, y: y_src };
}

// scale
function resize( x_dest, y_dest, scale ) {
	// params: double scale_horizontal, double scale_vertical;
	var x_src = x_dest * scale[0];
	var y_src = y_dest * scale[1];

	return { x: x_src, y: y_src }
}

// horiz shift
function horiz( x_dest, y_dest, shift ){

	var x_src	= x_dest + shift;
	var y_src  = y_dest;

	return { x: x_src, y: y_src }
}

// vertical shift
function vert( x_dest, y_dest, shift ){

	var x_src	= x_dest;
	var y_src  = y_dest + shift;

	return { x: x_src, y: y_src }
}

// perspective
function persp_sphere( x_dest, y_dest, mt, distance ) {
	// params :  double Matrix[3][3], double params.distance
	var theta,s,r;
	var v = new THREE.Vector3(),
		v2 = new THREE.Vector3();

	r = Math.sqrt( x_dest * x_dest + y_dest * y_dest );
	theta 	= r / distance;
	if( r == 0.0 )
		s = 0.0;
	else
		s = Math.sin( theta ) / r;

	v.x =  s * x_dest ;
	v.y =  s * y_dest ;
	v.z =  Math.cos( theta );

	v2 = mt.transformVector( v );

	r = Math.sqrt( v2.x*v2.x + v2.y*v2.y );

	if( r == 0.0 )
		theta = 0.0;
	else
		theta 	= distance * Math.atan2( r, v2.z ) / r;

	var x_src 	= theta * v2.x;
	var y_src 	= theta * v2.y;

	return { x: x_src, y: y_src }

}

// sphere projection to equirectangular
function erect_sphere_tp( x_dest, y_dest, distance ){

	// params: double params.distance
	var theta, r, s;
	var	v = [3];

	r = Math.sqrt( x_dest * x_dest + y_dest * y_dest );

	theta = r / distance;

	if (theta == 0.0)
		s = 1.0 / distance;
	else
		s = Math.sin( theta ) / r;

	v[1] =  s * x_dest;
	v[0] =  Math.cos( theta );

	var x_src = distance * Math.atan2( v[1], v[0] );
	var y_src = distance * Math.atan( s * y_dest / Math.sqrt( v[0]*v[0] + v[1]*v[1] ) );

	return { x: x_src, y: y_src }

}

// Rotate equirectangular image
function rotate_erect( x_dest, y_dest, mprot ) {

	// params: double 180degree_turn(screenpoints), double turn(screenpoints);

	var x_src = x_dest + mprot[1];

	while( x_src < - mprot[0] )
	            x_src += 2 * mprot[0];

	while( x_src >  mprot[0] )
	            x_src -= 2 * mprot[0];

	var y_src = y_dest;

	return { x: x_src, y: y_src }
}
