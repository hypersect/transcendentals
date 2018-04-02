/******************************************************************************
  Transcendentals.cpp

  Copyright (c) 2018 Hypersect
  http://www.hypersect.com/

  This software is provided 'as-is', without any express or implied
  warranty. In no event will the authors be held liable for any damages
  arising from the use of this software.

  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:

  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software. If you use this software
     in a product, an acknowledgment in the product documentation would be
     appreciated but is not required.

  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.

  3. This notice may not be removed or altered from any source
     distribution.
******************************************************************************/


#include "Trancendentals.h"
#include <math.h>
#include <assert.h>

//******************************************************************************
//******************************************************************************
#define RJ_ASSERT(condition) assert(condition)

//******************************************************************************
//******************************************************************************
const float c_Pi	  = 3.1415926535897932384626433832795f;
const float c_TwoPi	  = 6.2831853071795864769252867665590f;
const float c_HalfPi  = 1.5707963267948966192313216916398f;
const float c_SqrtTwo = 1.4142135623730950488016887242097f;

//******************************************************************************
// Wrapper functions that depend on math.h
// Note: CRT implementations should produce equivalent results across platforms
//******************************************************************************
static inline float  IsInfinity(float val)  { return ::isinf(val); }
static inline double IsInfinity(double val) { return ::isinf(val); }
static inline float  Sqrt(float val)        { return ::sqrtf(val); }
static inline double Sqrt(double val)       { return ::sqrt(val); }
static inline float  Abs(float val)         { return ::fabsf(val); }
static inline double Abs(double val)        { return ::fabs(val); }
static inline float  Ceil(float val)        { return ::ceilf(val); }
static inline double Ceil(double val)       { return ::ceil(val); }
static inline float  Floor(float val)       { return ::floorf(val); }
static inline double Floor(double val)      { return ::floor(val); }

//******************************************************************************
// Get the log base 2 of a 32-bit unsigned integer.
// http://graphics.stanford.edu/~seander/bithacks.html#IntegerLogLookup
//******************************************************************************
uint32_t Log2( uint32_t val )
{
	static const uint8_t logTable[256] = 
	{
		0, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
		5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
		5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
		6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
		6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
		6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
		6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
		7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
		7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
		7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
		7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
		7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
		7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
		7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
		7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7
	};

	uint32_t temp;

	temp = val >> 24;
	if (temp)
		return 24 + logTable[temp];

	temp = val >> 16;
	if (temp)
		return 16 + logTable[temp];

	temp = val >> 8;
	if (temp)
		return 8 + logTable[temp];

	return logTable[val];
}

//******************************************************************************
// Get the log base 2 of a 64-bit unsigned integer.
//******************************************************************************
inline uint64_t Log2( uint64_t val )
{
	uint64_t temp;

	temp = val >> 32ull;
	if (temp)
		return 32 + Log2((uint32_t)temp);

	return Log2((uint32_t)val);
}
		
//******************************************************************************
// Log2
// return base 2 logarithm of value
//******************************************************************************

// Compute Log2 for inputs greater than 0.5 to 4.14 digits accuracy
// Performs range reduction to get arg into [0.5,1.0] by using following 
// identities where n is an integer that puts f in [0.5,1.0]
// log2(x) = log2(f * 2^n)
//         = log2(f) + log2(2^n)
//         = log2(f) + n
// based on http://www.ganssle.com/approx-2.htm
static inline float LogBase2_4_14(float arg)
{
	// To prevent overflow, the input needs to be less than 2^63 for the integer cast
	// followed by a potential left shift.
	const float twoPowSixtyThree = 9223372036854775808.0f;
	if (arg >= twoPowSixtyThree)
	{
		// Factor out log2(2^32) and recurse.
		return 63.0f + LogBase2_4_14(arg/twoPowSixtyThree);
	}

	uint64_t intArg = (uint64_t)arg;

	// compute integer power of 2 that is greater than or equal to the input argument
	uint64_t logInt    = Log2(intArg); // logInt    = log2(intArg)
	uint64_t expLogInt = 1ull << logInt;   // expLogInt = 2^logInt
	if (intArg > expLogInt)                // increase power if not bounded
	{	
		RJ_ASSERT(expLogInt << 1ull >= intArg);
		logInt += 1ull;
		expLogInt <<= 1ull;
	}

	float n = (float)logInt;

	// remove the 2^n scale to get arg into [0.5, 1]
	arg = arg / (float)expLogInt;
	
	// The format of the polynomial is P(x)/Q(x)
	const float P00 = -1.45326486f;
	const float P01 =  0.951366714f;
	const float P02 =  0.501994886f;
	const float Q00 =  0.352143751f; 
	
	float poly= (P00 + arg * (P01 + arg*P02)) / (Q00 + arg);

	// merge the fractional and integral logarithms back together
	return poly + n;
}

static inline double LogBase2_4_14(double arg)
{
	// To prevent overflow, the input needs to be less than 2^63 for the integer cast
	// followed by a potential left shift.
	const float twoPowSixtyThree = 9223372036854775808.0;
	if (arg >= twoPowSixtyThree)
	{
		// Factor out log2(2^32) and recurse.
		return 63.0 + LogBase2_4_14(arg/twoPowSixtyThree);
	}
	
	uint64_t intArg = (uint64_t)arg;

	// compute integer power of 2 that is greater than or equal to the input argument
	uint64_t logInt    = Log2(intArg); // logInt    = log2(intArg)
	uint64_t expLogInt = 1ull << logInt;   // expLogInt = 2^logInt
	if (intArg < expLogInt)               // increase power if not bounded
	{	
		RJ_ASSERT(expLogInt << 1ull >= intArg);
		logInt += 1ull;
		expLogInt <<= 1ull;
	}

	double n = (double)logInt;

	// remove the 2^n scale to get arg into [0.5, 1]
	arg = arg / (double)expLogInt;
	
	// The format of the polynomial is P(x)/Q(x)
	const double P00 = -1.45326486;
	const double P01 =  0.951366714;
	const double P02 =  0.501994886;
	const double Q00 =  0.352143751; 
	
	double poly= (P00 + arg * (P01 + arg*P02)) / (Q00 + arg);

	// merge the fractional and integral logarithms back together
	return poly + n;
}

float Log2(float value)
{
	// Logarithm can only be taken of positive non-zero values. 
	RJ_ASSERT(value > 0.0f);

	// Log approximation is optimized for values greater than 0.5
	if(value < 0.5f)
	{
		// convert to and from optimal range using identiyy log(x) = -log(1/x)
		return -LogBase2_4_14(1.0f / value);
	}
	else
	{
		return LogBase2_4_14(value);
	}
};

double Log2(double value)
{
	// Logarithm can only be taken of positive non-zero values. 
	RJ_ASSERT(value > 0.0);

	// Log approximation is optimized for values greater than 0.5
	if(value < 0.5)
	{
		// convert to and from optimal range using identiyy log(x) = -log(1/x)
		return -LogBase2_4_14(1.0 / value);
	}
	else
	{
		return LogBase2_4_14(value);
	}
};

//******************************************************************************
// Get base 10 logarithm
//******************************************************************************
float Log10(float value)  
{
	const float invLog2of10 = 0.30102999566f; 
	return Log2(value) * invLog2of10; 
}

double Log10(double value) 
{ 
	const double invLog2of10 = 0.30102999566; 
	return Log2(value) * invLog2of10; 
}

//******************************************************************************
// Get natural logarithm
//******************************************************************************
float Log(float value)
{
	const float invLog2ofE = 0.69314718056f; 
	return Log2(value) * invLog2ofE;
}

double Log(double value) 
{ 
	const double invLog2ofE = 0.69314718056; 
	return Log2(value) * invLog2ofE; 
}

//******************************************************************************
// Exp2
// return	value raised to a power
//
// Use approximation valid for range [0, 0.5]
// based on http://www.ganssle.com/approx-2.htm
// 
// q = 25.0391066503 + x^2
// p = x*8.6778388279
// result = (q + p) / (q - p)
//
// Note: A faster, but less accurate approximation for [0,1] is
//       result = 0.99992520 + 0.69583356*x  + 0.22606716 * x^2  +  0.078024521 * x^3
//       http://jrfonseca.blogspot.com/2008/09/fast-sse2-pow-tables-or-polynomials.html
//******************************************************************************
static inline float ExpBase2_Positive(float exponent)
{
	RJ_ASSERT(exponent >= 0.0f);

	const float q0 = 25.0391066503f;
	const float p1 = 8.6778388279f;

	// Split into integer and fractional part such that we can compute separately
	// and later recombine using 2^x = 2^ipart * 2^fpart
	uint32_t ipart = (uint32_t)(exponent);
	float    fpart = exponent - ipart;

	float expFpart;
	if (fpart <= 0.5)
	{
		float q = q0 + fpart*fpart;
		float p = p1 * fpart;
		expFpart = (q+p) / (q-p);
	}
	else
	{
		// map into [0,0.5] using following relation:
		// 2^fpart = 2^0.5 * 2^(fpart-0.5)
		fpart -= 0.5f;

		float q = q0 + fpart*fpart;
		float p = p1 * fpart;
		expFpart = c_SqrtTwo * (q+p) / (q-p);			 
	}

	// exploit floating point format to compute 2^ipart
	// 1 sign bit:        0
	// 8 exponent bits:  exponent+127 
	// 23 mantissa bits: 0
	if (ipart > 1023)
	{
		return HUGE_VALF;
	}
	else
	{
		union { double f; uint64_t i; } expIpart;
		expIpart.i = ((ipart+1023ull) << 52ull);

		// Recombine integer and fractional parts for result
		return (float)(expIpart.f * expFpart);
	}
}

static inline double ExpBase2_Positive(double exponent)
{
	RJ_ASSERT(exponent >= 0.0f);

	const double q0 = 25.0391066503f;
	const double p1 = 8.6778388279f;

	// Split into integer and fractional part such that we can compute separately
	// and later recombine using 2^x = 2^ipart * 2^fpart
	uint64_t ipart = (uint64_t)(exponent);
	double   fpart = exponent - ipart;

	double expFpart;
	if (fpart <= 0.5)
	{
		double q = q0 + fpart*fpart;
		double p = p1 * fpart;
		expFpart = (q+p) / (q-p);
	}
	else
	{
		// map into [0,0.5] using following relation:
		// 2^fpart = 2^0.5 * 2^(fpart-0.5)
		fpart -= 0.5f;

		double q = q0 + fpart*fpart;
		double p = p1 * fpart;
		expFpart = c_SqrtTwo * (q+p) / (q-p);			 
	}

	// exploit floating point format to compute 2^x
	// 1 sign bit:       0
	// 11 exponent bits: exponent+1023 
	// 52 mantissa bits: 0
	if (ipart > 1023)
	{
		return HUGE_VAL;
	}
	else
	{
		union { double f; uint64_t i; } expIpart;
		expIpart.i = ((ipart+1023ull) << 52ull);

		// Recombine integer and fractional parts for result
		return expIpart.f * expFpart;
	}
}

float Exp2(float exponent)
{
	if (exponent > 0.0f)
	{
		// check for positive infinity
		if (IsInfinity(exponent))
			return exponent;

		return ExpBase2_Positive(exponent);
	}
	else
	{
		// check for negative infinity
		if (IsInfinity(exponent))
			return 0.0f;

		// 2^-x = 1 / 2^x
		return 1.0f / ExpBase2_Positive(-exponent);
	}	
}

double Exp2(double exponent)
{
	if (exponent > 0.0)
	{
		// check for positive infinity
		if (IsInfinity(exponent))
			return exponent;

		return ExpBase2_Positive(exponent);
	}
	else
	{
		// check for negative infinity
		if (IsInfinity(exponent))
			return 0.0;

		// 2^-x = 1 / 2^x
		return 1.0 / ExpBase2_Positive(-exponent);
	}	
}

//******************************************************************************
// Get e raised to a power
//******************************************************************************
float Exp(float exponent)
{ 
	const float log2ofE = 1.44269504089f; 
	return Exp2(log2ofE*exponent); 
}

double Exp(double exponent) 
{
	const double log2ofE = 1.44269504089; 
	return Exp2(log2ofE*exponent); 
}

//******************************************************************************
//******************************************************************************
float Pow(float base, float exponent)   
{
	// Note: There are lots of edge cases for IEEE compliant pow that we are not properly supporting.
	//       For example: negative exponents, zeros and infinities.
	if (base == 0.0f)
		return 0.0f;
	else if (exponent == 1.0f)
		return base;

	return Exp2(Log2(base)*exponent);
}

//******************************************************************************
//******************************************************************************
double Pow(double base, double exponent) 
{
	// Note: There are lots of edge cases for IEEE compliant pow that we are not properly supporting.
	//       For example: negative exponents, zeros and infinities.
	if (base == 0.0)
		return 0.0;
	else if (exponent == 1.0)
		return base;

	return Exp2(Log2(base)*exponent); 
}

//******************************************************************************
// Sine
// The current implementation is based on http://www.coranac.com/2009/07/sines/
//******************************************************************************
float Sin(float val)
{
	const float inv2Pi = 1.0f / c_TwoPi;

	// map [-pi/2,3pi/2] to [0,1]
	val = val*inv2Pi + 0.25f;

	// wrap into [0,1] by subtracting floor(val)
	val -= Floor(val);

	// map from wrapped [0,1] to [-1,3] and reflect into [-1,1]
	if (val > 0.5f)
		val = 3.0f - 4.0f*val;
	else
		val = 4.0f*val - 1.0f;

	// Fifth degree polynomial antisymetric about y=0 and optimized to reduce error in range [0,1] (aka [0,pi/2])
	// Maximum error is about 0.0002 (at pi/8 and 3pi/8).
	float a = 1.569718634205488058453210320940344688827031497770954769944f; // 4*(3/pi - 9/16)
	float b = 0.639437268410976116906420641880689377654062995541909539888f; // 2*4*(3/pi - 9/16) - 5/2
	float c = 0.069718634205488058453210320940344688827031497770954769944f; // 4*(3/pi - 9/16) - 3/2

	float val2 = val*val;
	return val*(val2*(c*val2 - b) + a);
}

double Sin(double val)
{
	const double inv2Pi = 1.0 / c_TwoPi;

	// map [-pi/2,3pi/2] to [0,1]
	val = val*inv2Pi + 0.25f;

	// wrap into [0,1] by subtracting floor(val)
	val -= Floor(val);
	
	// map from wrapped [0,1] to [-1,3] and reflect into [-1,1]
	if (val > 0.5)
		val = 3.0 - 4.0*val;
	else
		val = 4.0*val - 1.0f;

	// Fifth degree polynomial antisymetric about y=0 and optimized to reduce error in range [0,1] (aka [0,pi/2])
	// Maximum error is about 0.0002 (at pi/8 and 3pi/8).
	double a = 1.569718634205488058453210320940344688827031497770954769944; // 4*(3/pi - 9/16)
	double b = 0.639437268410976116906420641880689377654062995541909539888; // 2*4*(3/pi - 9/16) - 5/2
	double c = 0.069718634205488058453210320940344688827031497770954769944; // 4*(3/pi - 9/16) - 3/2

	double val2 = val*val;
	return val*(val2*(c*val2 - b) + a);
}

//******************************************************************************
// Cosine
//******************************************************************************
float Cos(float val) 
{
	return Sin(c_HalfPi - val);
}

double Cos(double val)
{
	return Sin(c_HalfPi - val); 
}

//******************************************************************************
// Tangent
//******************************************************************************

// Compute tangent to 5.6 digits of accuracy where input is [-1,1] representing [-pi/4,pi/4] radians.
// Based on http://www.ganssle.com/approx.htm
inline static float Tan_5_6(float val)
{
	const float c1 = -3.16783027f;
	const float c2 =  0.134516124f;
	const float c3 = -4.033321984f;

	float val2 = val*val;
	return val*(c1 + c2*val2)/(c3+val2);
}

inline static double Tan_5_6(double val)
{
	const double c1 = -3.16783027;
	const double c2 =  0.134516124;
	const double c3 = -4.033321984;

	double val2 = val*val;
	return val*(c1 + c2*val2)/(c3+val2);
}

float Tan(float val)
{
	val = val / c_TwoPi; // map [0,2pi] to [0,1]
	val -= Floor(val); // wrap into [0,1]

	val *= 8.0f; // map [0,1] to [0,8]
	int octant = (int)(val);
	switch (octant)
	{
	case 0:  return         Tan_5_6( val );
	case 1:  return  1.0f / Tan_5_6( 2.0f - val );
	case 2:  return -1.0f / Tan_5_6( val - 2.0f );
	case 3:  return        -Tan_5_6( 4.0f - val );
	case 4:  return         Tan_5_6( val - 4.0f );
	case 5:  return  1.0f / Tan_5_6( 6.0f - val );
	case 6:  return -1.0f / Tan_5_6( val - 6.0f );
	default: return        -Tan_5_6( 8.0f - val );
	}
}

double Tan(double val)
{
	val = val / c_TwoPi; // map [0,2pi] to [0,1]
	val -= Floor(val); // wrap into [0,1]

	val *= 8.0; // map [0,1] to [0,8]
	int octant = (int)(val);
	switch (octant)
	{
	case 0:  return        Tan_5_6( val );
	case 1:  return  1.0 / Tan_5_6( 2.0 - val );
	case 2:  return -1.0 / Tan_5_6( val - 2.0f );
	case 3:  return       -Tan_5_6( 4.0 - val );
	case 4:  return        Tan_5_6( val - 4.0f );
	case 5:  return  1.0 / Tan_5_6( 6.0 - val );
	case 6:  return -1.0 / Tan_5_6( val - 6.0f );
	default: return       -Tan_5_6( 8.0 - val );
	}
}

//******************************************************************************
// Asin
// TODO: Investigate a polynomial approximation or use a fast inverse square
//       root to speed this up.
//******************************************************************************
float Asin(float val)
{
	if (val >= 1.0f)
		return c_HalfPi;
	else if (val <= -1.0f)
		return -c_HalfPi;
	else
		return Atan( val / Sqrt(1.0f - val*val) );
}
double Asin(double val)
{ 
	if (val >= 1.0)
		return c_HalfPi;
	else if (val <= -1.0)
		return -c_HalfPi;
	else
		return Atan( val / Sqrt(1.0 - val*val) );
}

//******************************************************************************
// Acos
// TODO: Investigate a polynomial approximation or use a fast inverse square
//       root to speed this up.
//******************************************************************************
float Acos(float val)
{
	return c_HalfPi - Asin(val);
}
double Acos(double val)
{
	return c_HalfPi - Asin(val);
}

//******************************************************************************
// Atan
//******************************************************************************

// Compute arctangent to 6.6 digits of accuracy over the output range of [-pi/12, pi/12]
// The current implementation is based on http://www.coranac.com/2009/07/sines/
static inline float Atan_6_6(float x)
{
	const float c1=1.6867629106f;
	const float c2=0.4378497304f;
	const float c3=1.6867633134f;

	float x2 = x*x;

	return (x*(c1 + x2*c2)/(c3 + x2));
}

static inline double Atan_6_6(double x)
{
	const double c1=1.6867629106;
	const double c2=0.4378497304;
	const double c3=1.6867633134;

	double x2 = x*x;

	return (x*(c1 + x2*c2)/(c3 + x2));
}

float Atan(float val)
{
	const float piOverSix       = 0.52359877559829887307710723054658f;
	const float tanPiOverSix    = 0.57735026918962576450914878050196f;
	const float tanPiOverTwelve = 0.26794919243112270647255365849413f;

	// keep arg positive
	float signScale = 1.0f;
	if (val < 0.0f)
	{
		// arctan(-x) = -arctan(x)
		val = -val;
		signScale = -1.0f;
	}

	// keep arg between 0 and 1
	float complementOffset = 0.0;
	float complementScale  = 1.0;
	if (val > 1.0f)
	{
		val = 1.0f/val;
		complementOffset = c_HalfPi;
		complementScale  = -1.0f;
	}

	// keep arg less than or equal to tan(pi/12)
	float regionOffset = 0.0f;
	if (val > tanPiOverTwelve)
	{
		// Map value from [tan(pi/12), 1] to [-tan(pi/12), tan(pi/12)] with:
		//  remappedX = ( x - tan(pi/6) ) / (1 + x*tan(pi/6) )
		// We can then use the following identity to unmap the result
		//  atan(a) - atan(b) = atan( (a-b) / (1 + a*b) )
		// We set b=tan(pi/6) and convert like so:
		//  atan(x) = pi/6 + atan(remappedX)
		val = (val - tanPiOverSix) / (1.0f + val*tanPiOverSix); 
		regionOffset = piOverSix;
	}

	// run the approximation
	float y = Atan_6_6(val);
	
	// remap back to correct space
	y = y + regionOffset; 
	y = complementOffset + complementScale*y; 
	y = signScale * y;

	return (y);
}

double Atan(double val)
{
	const double piOverSix       = 0.52359877559829887307710723054658;
	const double tanPiOverSix    = 0.57735026918962576450914878050196;
	const double tanPiOverTwelve = 0.26794919243112270647255365849413;

	// keep arg positive
	double signScale = 1.0;
	if (val < 0.0)
	{
		// arctan(-x) = -arctan(x)
		val = -val;
		signScale = -1.0;
	}

	// keep arg between 0 and 1
	double complementOffset = 0.0;
	double complementScale  = 1.0;
	if (val > 1.0)
	{
		val = 1.0/val;
		complementOffset = c_HalfPi;
		complementScale  = -1.0;
	}

	// keep arg less than or equal to tan(pi/12)
	double regionOffset = 0.0;
	if (val > tanPiOverTwelve)
	{
		// Map value from [tan(pi/12), 1] to [-tan(pi/12), tan(pi/12)] with:
		//  remappedX = ( x - tan(pi/6) ) / (1 + x*tan(pi/6) )
		// We can then use the following identity to unmap the result
		//  atan(a) - atan(b) = atan( (a-b) / (1 + a*b) )
		// We set b=tan(pi/6) and convert like so:
		//  atan(x) = pi/6 + atan(remappedX)
		val = (val - tanPiOverSix) / (1.0 + val*tanPiOverSix); 
		regionOffset = piOverSix;
	}

	// run the approximation
	double y = Atan_6_6(val);
	
	// remap back to correct space
	y = y + regionOffset; 
	y = complementOffset + complementScale*y; 
	y = signScale * y;

	return (y);
}

//******************************************************************************
// Returns in range [-pi,pi]
// Based on https://www.dsprelated.com/showarticle/1052.php
//******************************************************************************
float Atan2(float y, float x)
{
    if (x != 0.0f)
    {
		// if (y/x) is outsize the [-1,1] range
        if (Abs(x) > Abs(y))
        {
            const float z = y / x;
            if (x > 0.0f)
            {
                // atan2(y,x) = atan(y/x) if x > 0
                return Atan(z);
            }
            else if (y >= 0.0f)
            {
                // atan2(y,x) = atan(y/x) + PI if x < 0, y >= 0
                return Atan(z) + c_Pi;
            }
            else
            {
                // atan2(y,x) = atan(y/x) - PI if x < 0, y < 0
                return Atan(z) - c_Pi;
            }
        }
        else // Use property atan(y/x) = PI/2 - atan(x/y) if |y/x| > 1.
        {
            const float z = x / y;
            if (y > 0.0f)
            {
                // atan2(y,x) = PI/2 - atan(x/y) if |y/x| > 1, y > 0
                return -Atan(z) + c_HalfPi;
            }
            else
            {
                // atan2(y,x) = -PI/2 - atan(x/y) if |y/x| > 1, y < 0
                return -Atan(z) - c_HalfPi;
            }
        }
    }
    else
    {
        if (y > 0.0f) // x = 0, y > 0
        {
            return c_HalfPi;
        }
        else if (y < 0.0f) // x = 0, y < 0
        {
            return -c_HalfPi;
        }
		else
		{
		    return 0.0f; // TODO: investigate checking sign bits of zeros
		}
    }
}

double Atan2(double y, double x)
{
    if (x != 0.0)
    {
		// if (y/x) is outsize the [-1,1] range
        if (Abs(x) > Abs(y))
        {
            const double z = y / x;
            if (x > 0.0)
            {
                // atan2(y,x) = atan(y/x) if x > 0
                return Atan(z);
            }
            else if (y >= 0.0)
            {
                // atan2(y,x) = atan(y/x) + PI if x < 0, y >= 0
                return Atan(z) + c_Pi;
            }
            else
            {
                // atan2(y,x) = atan(y/x) - PI if x < 0, y < 0
                return Atan(z) - c_Pi;
            }
        }
        else // Use property atan(y/x) = PI/2 - atan(x/y) if |y/x| > 1.
        {
            const double z = x / y;
            if (y > 0.0)
            {
                // atan2(y,x) = PI/2 - atan(x/y) if |y/x| > 1, y > 0
                return -Atan(z) + c_HalfPi;
            }
            else
            {
                // atan2(y,x) = -PI/2 - atan(x/y) if |y/x| > 1, y < 0
                return -Atan(z) - c_HalfPi;
            }
        }
    }
    else
    {
        if (y > 0.0) // x = 0, y > 0
        {
            return c_TwoPi;
        }
        else if (y < 0.0f) // x = 0, y < 0
        {
            return -c_TwoPi;
        }
		else
		{
		    return 0.0; // TODO: investigate checking sign bits of zeros
		}
    }
}

