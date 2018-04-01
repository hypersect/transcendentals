/******************************************************************************
  Transcendentals.h

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

#include <cstdint>

//******************************************************************************
// Get the base-2 logarithm of value
//******************************************************************************
uint32_t Log2(uint32_t val);
uint64_t Log2(uint64_t val);
float    Log2(float value);
double   Log2(double value);

//******************************************************************************
// Get base-10 logarithm
//******************************************************************************
float  Log10(float value);
double Log10(double value);

//******************************************************************************
// Get base-e logarithm (natural logarithm)
//******************************************************************************
float  Log(float value);
double Log(double value);

//******************************************************************************
// Get base-2 exponential
//******************************************************************************
float  Exp2(float exponent);
double Exp2(double exponent);

//******************************************************************************
// Get base-e exponentual
//******************************************************************************
float  Exp(float exponent);
double Exp(double exponent);
		
//******************************************************************************
// Get value raised to a power
//******************************************************************************
float  Pow(float base, float exponent);
double Pow(double base, double exponent);
		
//******************************************************************************
// Sin
//******************************************************************************
float  Sin(float val);
double Sin(double val);

//******************************************************************************
// Cosine
//******************************************************************************
float  Cos(float val);
double Cos(double val);

//******************************************************************************
// Tangent
//******************************************************************************
float  Tan(float val);
double Tan(double val);

//******************************************************************************
// Arcsine
//******************************************************************************
float  Asin(float val);
double Asin(double val);

//******************************************************************************
// Arccosine
//******************************************************************************
float  Acos(float val);
double Acos(double val);

//******************************************************************************
// Arctangent
//******************************************************************************
float  Atan(float val);
double Atan(double val);

//******************************************************************************
// Arctangent of y/x
//******************************************************************************
float  Atan2(float y, float x);
double Atan2(double y, double x);
