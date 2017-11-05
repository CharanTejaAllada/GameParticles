#ifndef PARTICLE_H
#define PARTICLE_H


#include "Math\Vect4D.h"
#include "Math\Align.h"

class Particle:public Align
{
public:
	friend class ParticleEmitter;

	Particle();
	~Particle();
	void Update(const float& time_elapsed);

	Particle *next;
	Particle *prev;
	float	rotation;
	float	life;
	float	rotation_velocity;
	Vect4D	position;
	Vect4D	velocity;

	Vect4D	scale;
};


#endif //PARTICLE_H
