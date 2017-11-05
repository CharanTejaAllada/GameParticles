#ifndef MATRIX_H
#define MATRIX_H

// includes
#include "Enum.h"
#include"Vect4D.h"
#include<xmmintrin.h>
#include<smmintrin.h>

// forward declare
class Vect4D;

// class
class Matrix
{
public:

	// local enumerations
	enum MatrixRowEnum
	{
		MATRIX_ROW_0,
		MATRIX_ROW_1,
		MATRIX_ROW_2,
		MATRIX_ROW_3
	};

	Matrix();
	Matrix( Matrix&  t);
	~Matrix();

	void get(MatrixRowEnum row, Vect4D *vOut);

	void setIdentMatrix();
	void setTransMatrix( Vect4D * );
	void setScaleMatrix( Vect4D * );
	void setRotZMatrix( float Z_Radians);
	void UpdateRotTransMatrix(Vect4D &, float &);
	void UpdateCameraScale(Matrix &);

	Matrix operator*(  Matrix &t);
	void Inverse(Matrix &out);

	union
	{
		struct
		{
			Vect4D v0;
			Vect4D v1;
			Vect4D v2;
			Vect4D v3;
		};

		struct
		{
			// ROW 0
			float m0;
			float m1;
			float m2;
			float m3;

			// ROW 1
			float m4;
			float m5;
			float m6;
			float m7;

			// ROW 2
			float m8;
			float m9;
			float m10;
			float m11;

			// ROW 3
			float m12;
			float m13;
			float m14;
			float m15;
		};
	};
};

#endif  // Matrix.h
