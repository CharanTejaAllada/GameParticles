#include <Math.h>
#include <assert.h>
#include "Vect4d.h"
#include "Matrix.h"



Matrix::Matrix()
{	// constructor for the matrix

}

Matrix::Matrix(Matrix& t)
{ // copy constructor
	this->v0 = t.v0;
	this->v1 = t.v1;
	this->v2 = t.v2;
	this->v3 = t.v3;

}

Matrix::~Matrix()
{
	// nothign to delete
}

void Matrix::setIdentMatrix()
{ // initialize to the identity matrix
	this->v0._m = _mm_set_ps(0.0f, 0.0f, 0.0f, 1.0f);
	this->v1._m = _mm_set_ps(0.0f, 0.0f, 1.0f, 0.0f);
	this->v2._m = _mm_set_ps(0.0f, 1.0f, 0.0f, 0.0f);
	this->v3._m = _mm_set_ps(1.0f, 0.0f, 0.0f, 0.0f);
	//simd in reverse
}

void Matrix::setTransMatrix(Vect4D *t)
{ // set the translation matrix (note: we are row major)
	this->v0._m = _mm_set_ps(0.0f, 0.0f, 0.0f, 1.0f);
	this->v1._m = _mm_set_ps(0.0f, 0.0f, 1.0f, 0.0f);
	this->v2._m = _mm_set_ps(0.0f, 1.0f, 0.0f, 0.0f);
	this->v3._m = _mm_set_ps(1.0f, t->z, t->y, t->x);
	//simd in reverse

}



void Matrix::get(MatrixRowEnum, Vect4D *t)
{ 
	t->_m = this->v3._m;
	//switches removed
	//only this is needed....other conditions cleaned
}


Matrix Matrix::operator*(Matrix& p)
{
	// matrix multiplications
	Matrix tmp;
	float *t1, *t2, *t3;


	t1 = (float*)&this->m0;
	t2 = (float*)&p.m0;
	t3 = &tmp.m0;

	__m128 row1 = _mm_load_ps(&t2[0]);
	__m128 row2 = _mm_load_ps(&t2[4]);
	__m128 row3 = _mm_load_ps(&t2[8]);
	__m128 row4 = _mm_load_ps(&t2[12]);


	for (int i = 0; i<4; i++)
	{
		__m128 brod1 = _mm_set1_ps(t1[4 * i + 0]);
		__m128 brod2 = _mm_set1_ps(t1[4 * i + 1]);
		__m128 brod3 = _mm_set1_ps(t1[4 * i + 2]);
		__m128 brod4 = _mm_set1_ps(t1[4 * i + 3]);
		__m128 row = _mm_add_ps(_mm_add_ps(_mm_mul_ps(brod1, row1), 
		 _mm_mul_ps(brod2, row2)), 
  _mm_add_ps(_mm_mul_ps(brod3, row3), 
	 _mm_mul_ps(brod4, row4)));
		_mm_store_ps(&t3[4 * i], row);
	}

	return tmp;
}



void Matrix::Inverse(Matrix &out)
{
	__m128 m01, m_1, m, m1_3, r0, r1, r2, r3, det, tmp1;
	float *src = &this->m0;


	tmp1 = _mm_setr_ps(1.0f, 0.0f, 0.0f, 0.0f);
	det = tmp1;
	r3 = r2 = r1 = tmp1;
	tmp1 = _mm_loadh_pi(_mm_loadl_pi(tmp1, (__m64*)(src)), (__m64*)(src + 4));
	r1 = _mm_loadh_pi(_mm_loadl_pi(r1, (__m64*)(src + 8)), (__m64*)(src + 12));
	r0 = _mm_shuffle_ps(tmp1, r1, 0x88);
	r1 = _mm_shuffle_ps(r1, tmp1, 0xDD);
	tmp1 = _mm_loadh_pi(_mm_loadl_pi(tmp1, (__m64*)(src + 2)), (__m64*)(src + 6));
	r3 = _mm_loadh_pi(_mm_loadl_pi(r3, (__m64*)(src + 10)), (__m64*)(src + 14));
	r2 = _mm_shuffle_ps(tmp1, r3, 0x88);
	r3 = _mm_shuffle_ps(r3, tmp1, 0xDD);

	tmp1 = _mm_mul_ps(r2, r3);
	tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);
	m01 = _mm_mul_ps(r1, tmp1);
	m_1 = _mm_mul_ps(r0, tmp1);
	tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
	m01 = _mm_sub_ps(_mm_mul_ps(r1, tmp1), m01);
	m_1 = _mm_sub_ps(_mm_mul_ps(r0, tmp1), m_1);
	m_1 = _mm_shuffle_ps(m_1, m_1, 0x4E);
	// -----------------------------------------------
	tmp1 = _mm_mul_ps(r1, r2);
	tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);
	m01 = _mm_add_ps(_mm_mul_ps(r3, tmp1), m01);
	m1_3 = _mm_mul_ps(r0, tmp1);
	tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
	m01 = _mm_sub_ps(m01, _mm_mul_ps(r3, tmp1));
	m1_3 = _mm_sub_ps(_mm_mul_ps(r0, tmp1), m1_3);
	m1_3 = _mm_shuffle_ps(m1_3, m1_3, 0x4E);
	// -----------------------------------------------
	tmp1 = _mm_mul_ps(_mm_shuffle_ps(r1, r1, 0x4E), r3);
	tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);
	r2 = _mm_shuffle_ps(r2, r2, 0x4E);
	m01 = _mm_add_ps(_mm_mul_ps(r2, tmp1), m01);
	m = _mm_mul_ps(r0, tmp1);
	tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
	m01 = _mm_sub_ps(m01, _mm_mul_ps(r2, tmp1));
	m = _mm_sub_ps(_mm_mul_ps(r0, tmp1), m);
	m = _mm_shuffle_ps(m, m, 0x4E);
	// -----------------------------------------------
	tmp1 = _mm_mul_ps(r0, r1);

	tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);
	m = _mm_add_ps(_mm_mul_ps(r3, tmp1), m);
	m1_3 = _mm_sub_ps(_mm_mul_ps(r2, tmp1), m1_3);
	tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
	m = _mm_sub_ps(_mm_mul_ps(r3, tmp1), m);
	m1_3 = _mm_sub_ps(m1_3, _mm_mul_ps(r2, tmp1));
	// -----------------------------------------------
	tmp1 = _mm_mul_ps(r0, r3);
	tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);
	m_1 = _mm_sub_ps(m_1, _mm_mul_ps(r2, tmp1));
	m = _mm_add_ps(_mm_mul_ps(r1, tmp1), m);
	tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
	m_1 = _mm_add_ps(_mm_mul_ps(r2, tmp1), m_1);
	m = _mm_sub_ps(m, _mm_mul_ps(r1, tmp1));
	// -----------------------------------------------
	tmp1 = _mm_mul_ps(r0, r2);
	tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);
	m_1 = _mm_add_ps(_mm_mul_ps(r3, tmp1), m_1);
	m1_3 = _mm_sub_ps(m1_3, _mm_mul_ps(r1, tmp1));
	tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);
	m_1 = _mm_sub_ps(m_1, _mm_mul_ps(r3, tmp1));
	m1_3 = _mm_add_ps(_mm_mul_ps(r1, tmp1), m1_3);
	// -----------------------------------------------
	det = _mm_mul_ps(r0, m01);
	det = _mm_add_ps(_mm_shuffle_ps(det, det, 0x4E), det);
	det = _mm_add_ss(_mm_shuffle_ps(det, det, 0xB1), det);
	tmp1 = _mm_rcp_ss(det);
	det = _mm_sub_ss(_mm_add_ss(tmp1, tmp1), _mm_mul_ss(det, _mm_mul_ss(tmp1, tmp1)));
	det = _mm_shuffle_ps(det, det, 0x00);
	m01 = _mm_mul_ps(det, m01);
	_mm_storel_pi((__m64*)(src), m01);
	_mm_storeh_pi((__m64*)(src + 2), m01);
	m_1 = _mm_mul_ps(det, m_1);
	_mm_storel_pi((__m64*)(src + 4), m_1);
	_mm_storeh_pi((__m64*)(src + 6), m_1);
	m = _mm_mul_ps(det, m);
	_mm_storel_pi((__m64*)(src + 8), m);
	_mm_storeh_pi((__m64*)(src + 10), m);
	m1_3 = _mm_mul_ps(det, m1_3);
	_mm_storel_pi((__m64*)(src + 12), m1_3);
	_mm_storeh_pi((__m64*)(src + 14), m1_3);

	out.v0._m = m01;
	out.v1._m = m_1;
	out.v2._m = m;
	out.v3._m = m1_3;

}

void Matrix::setScaleMatrix(Vect4D *scale)
{


	this->v0._m = _mm_set_ps(0.0f, 0.0f, 0.0f, scale->x);
	this->v1._m = _mm_set_ps(0.0f, 0.0f, scale->y, 0.0f);
	this->v2._m = _mm_set_ps(0.0f, scale->z, 0.0f, 0.0f);
	this->v3._m = _mm_set_ps(1.0f, 0.0f, 0.0f, 0.0f);
	//simd in reverse
}

void Matrix::UpdateRotTransMatrix(Vect4D &s, float &az)
{
	float a = sinf(az), b = (float)cosf(az);

	this->v0._m = _mm_set_ps(0.0f, 0.0f, -a*s.x, b*s.x);
	this->v1._m = _mm_set_ps(0.0f, 0.0f, b*s.y, a*s.y);
	this->v2._m = _mm_set_ps(0.0f, s.z, 0.0f, 0.0f);
	this->v3._m = _mm_set_ps(1.0f, 0.0f, 0.0f, 0.0f);
	//simd in reverse
}

void Matrix::UpdateCameraScale(Matrix &p)
{
	Vect4D t,part1,part2,part3,part4;

	part1 = p.v0;
	part2 = p.v1;
	part3 = p.v2;
	part4 = p.v3;

	t._m = _mm_set_ps1(this->m12);
	part1._m = _mm_mul_ps(t._m, part1._m);

	t._m = _mm_set_ps1(this->m13);
	part2._m = _mm_mul_ps(t._m, part2._m);

	t._m = _mm_set_ps1(this->m14);
	part3._m = _mm_mul_ps(t._m, part3._m);

	t._m = _mm_set_ps1(this->m15);
	part4._m = _mm_mul_ps(t._m, part4._m);

	t._m = _mm_add_ps(part1._m, part2._m);
	t._m = _mm_add_ps(t._m, part3._m);
	t._m = _mm_add_ps(t._m, part4._m);

	this->v3._m = t._m;
}

void Matrix::setRotZMatrix(float rot)
{


	float a = sinf(rot), b = (float)cosf(rot);

	this->v0._m = _mm_set_ps(0.0f, 0.0f, -a, b);
	this->v1._m = _mm_set_ps(0.0f, 0.0f, b, a);
	this->v2._m = _mm_set_ps(0.0f, 1.0f, 0.0f, 0.0f);
	this->v3._m = _mm_set_ps(1.0f, 0.0f, 0.0f, 0.0f);
	//simd in reverse
}

// End of file