#include "Particle.h"
#include "Math\Matrix.h"

Vect4D z_axis(0.0f, -0.25f, 1.0f);
Vect4D v(3.0f, 4.0f, 0.0f);

Particle::Particle()
{
	// construtor
	this->life = 0.0f;
	this->position.set(0.0f, 0.0f, 0.0f);
	this->velocity.set(0.0f, 0.0f, 0.0f);
	this->scale.set(1.0f, 1.0f, 1.0f);
	this->rotation = 0.0f;
	this->rotation_velocity = 0.5f;

}
//------------------------------------------------------
// unwanted methods removed
//---------------------------------------------------
Particle::~Particle()
{
	// nothing to do
}



void Particle::Update(const float& time_elapsed)
{
	// serious math below - magic secret sauce
	life += time_elapsed;
	position += (velocity * time_elapsed);	
	position.Cross(z_axis, v);
	v.norm(v);
	position += v * 0.05f * life;
	rotation  += rotation_velocity * time_elapsed *2.01f;
}


// End of file


