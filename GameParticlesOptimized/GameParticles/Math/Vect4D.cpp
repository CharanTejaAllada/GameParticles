#include <math.h>
#include <assert.h>

#include "Vect4D.h"


Vect4D::Vect4D()
{
	this->_m = _mm_setr_ps(0.0f, 0.0f, 0.0f, 1.0f);
}

Vect4D::Vect4D(const __m128 x)
{
	this->_m = x;
}

Vect4D::Vect4D(float tx, float ty, float tz, float tw)
{
	this->_m = _mm_set_ps(tw, tz, ty, tx);
}

Vect4D::~Vect4D()
{
	// nothing to delete
}

void Vect4D::norm(Vect4D& out)
{
	__m128 nrm = _mm_sqrt_ps(_mm_dp_ps(_m, _m, 0X77));
	out._m = _mm_div_ps(_m, nrm);
	out.w = 1.0f;

}

void Vect4D::operator += (Vect4D &t)
{
	this->_m=_mm_add_ps(this->_m, t._m);
}

Vect4D Vect4D::operator + (Vect4D &t)
{
	return Vect4D(_mm_add_ps(this->_m, t._m));
}

Vect4D Vect4D::operator - (Vect4D &t)
{
	return Vect4D(this->x - t.x, this->y - t.y, this->z - t.z);
}

Vect4D Vect4D::operator *(float a)
{
	__m128 tmp = _mm_set1_ps(a);
	return Vect4D(_mm_mul_ps(tmp, this->_m));
}

float& Vect4D::operator[](VECT_ENUM)
{
		return x;
}

void Vect4D::Cross(Vect4D& one, Vect4D& two)
{
	two._m = _mm_sub_ps(_mm_mul_ps(_mm_shuffle_ps(_m, _m, _MM_SHUFFLE(3, 0, 2, 1)), _mm_shuffle_ps(one._m, one._m, _MM_SHUFFLE(3, 1, 0, 2))),
	 _mm_mul_ps(_mm_shuffle_ps(_m, _m, _MM_SHUFFLE(3, 1, 0, 2)),
	_mm_shuffle_ps(one._m, one._m, _MM_SHUFFLE(3, 0, 2, 1))));
}

void Vect4D::set(float tx, float ty, float tz, float tw)
{
	this->_m = _mm_set_ps(tw, tz, ty, tx);
}

// End of file