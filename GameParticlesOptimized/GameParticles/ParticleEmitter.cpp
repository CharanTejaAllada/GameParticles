#include "DO_NOT_MODIFY\Timer.h"
#include "DO_NOT_MODIFY\GlobalTimer.h"
#include "DO_NOT_MODIFY\OpenGLInterface.h"
#include"DO_NOT_MODIFY\Trace.h"
#include "ParticleEmitter.h"
#include "Settings.h"

static const unsigned char squareColors[] =
{
	255,  255,  000,  255,
	000,  255,  255,  255,
	000,  000,  000,  000,
	255,  000,  255,  255,
};

static const float squareVertices[] =
{
	-0.015f,  -0.015f, 0.0f,
	 0.015f,  -0.015f, 0.0f,
	-0.015f,   0.015f, 0.0f,
	 0.015f,   0.015f, 0.0f,
};


ParticleEmitter::ParticleEmitter()
	: start_position(0.0f, 0.0f, 0.0f),
	  start_velocity(0.0f, 1.0f, 0.0f),
	max_life(MAX_LIFE),
	max_particles(NUM_PARTICLES),
	spawn_frequency(0.0000001f),
	last_spawn(globalTimer::getTimerInSec()),
	last_loop(globalTimer::getTimerInSec()),
	last_active_particle(-1),	
	vel_variance(1.0f, 4.0f, 0.4f),
	pos_variance(1.0f, 1.0f, 1.0f),
	scale_variance(2.5f),
	headParticle(0)
{
	// nothing to do
	
}



ParticleEmitter::~ParticleEmitter()
{
	// do nothing
}


void ParticleEmitter::update()
{
	// get current time
	float current_time = globalTimer::getTimerInSec();

	// spawn particles
	float time_elapsed = current_time - last_spawn;

	Particle *newParticle;
	
	while (last_active_particle < max_particles - 1)
	{
		// spawn a particle
		newParticle = new Particle();

		// initialize the particle
		newParticle->position = start_position;
		newParticle->velocity = start_velocity;
		

		// apply the variance
		this->Execute(newParticle->position, newParticle->velocity, newParticle->scale);

		// increment count
		last_active_particle++;
		// add to list
		this->addParticleToList(newParticle);
		
	}


	last_spawn = current_time;
	// total elapsed
	time_elapsed = current_time - last_loop;

	Particle *p = this->headParticle;
	// walk the particles

	while (p != 0)
	{
		p->Update(time_elapsed);

		if ((p->life > max_life) && (last_active_particle > 0))
		{
			// particle to remove
			newParticle = p;

			// need to squirrel it away.
			p = p->next;

			// remove last node
			this->removeParticleFromList(newParticle);

			// update the number of particles
			last_active_particle--;
		}
		else
		{
			// increment to next point
			p = p->next;
		}
	}

	last_loop = current_time;
}

void ParticleEmitter::addParticleToList(Particle *p)
{
	if (this->headParticle == 0)
	{
		this->headParticle = p;
		p->next = 0;
		p->prev = 0;
	}
	else
	{ 
		headParticle->prev = p;
		p->next = headParticle;
		p->prev = 0;
		headParticle = p;
	}

}

void ParticleEmitter::removeParticleFromList(Particle *p)
{
    if (p->prev == 0 && p->next != 0)
	{ // first on the list
		p->next->prev = 0;
		this->headParticle = p->next;
	}
	else if (p->prev != 0 && p->next == 0)
	{ // last on the list 
		p->prev->next = 0;
	}
	else
	{ 
		p->prev->next = p->next;
		p->next->prev = p->prev;
	}

	// bye bye
	delete p;
}


void ParticleEmitter::draw()
{

	Matrix transCamera, transParticle, rotParticle, tempo;


	// get the camera matrix from OpenGL
	glGetFloatv(GL_MODELVIEW_MATRIX, reinterpret_cast<float*>(&transCamera));

	//Vect4D camPosVect;
	// OpenGL goo... don't worrry about this
	glVertexPointer(3, GL_FLOAT, 0, squareVertices);
	glEnableClientState(GL_VERTEX_ARRAY);
	glColorPointer(4, GL_UNSIGNED_BYTE, 0, squareColors);
	glEnableClientState(GL_COLOR_ARRAY);

	Particle *it=this->headParticle;
	rotParticle.setRotZMatrix(it->rotation);
	

	while(it != 0)
	{
		// particle position
		transParticle.setTransMatrix(&it->position);


		tempo = transCamera;
		tempo.m12 = tempo.m12 + tempo.m0 * transParticle.m12;
		tempo.m13 = tempo.m13 + tempo.m0 * transParticle.m13;
		tempo.m14 = tempo.m14 + tempo.m0 * transParticle.m14;

		Matrix tmp = tempo *rotParticle;

		// set the transformation matrix
		glLoadMatrixf(reinterpret_cast<float*>(&(tmp)));

		// draw the trangle strip
		glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
		it = it->next;

	}

}



void ParticleEmitter::Execute(Vect4D& pos, Vect4D& vel, Vect4D& sc)
{

	float var = static_cast<float>(rand() % 1000) * 0.001f;
	float sign = static_cast<float>(rand() % 2);
	float *t_pos = reinterpret_cast<float*>(&pos);
	float *t_var = &pos_variance[x];
	if (sign == 0)
	{
		var *= -1.0;
	}
	*t_pos += *t_var * var;
	var = static_cast<float>(rand() % 1000) * 0.001f;
	sign = static_cast<float>(rand() % 2);
	t_pos++;
	t_var++;
	if (sign == 0)
	{
		var *= -1.0;
	}
	*t_pos += *t_var * var;
	var = static_cast<float>(rand() % 1000) * 0.001f;
	sign = static_cast<float>(rand() % 2);
	t_pos++;
	t_var++;
	if (sign == 0)
	{
		var *= -1.0;
	}
	*t_pos += *t_var * var;

	var = static_cast<float>(rand() % 1000) * 0.001f;
	sign = static_cast<float>(rand() % 2);

	t_pos = &vel[x];
	t_var = &vel_variance[x];
	if (sign == 0)
	{
		var *= -1.0;
	}
	*t_pos += *t_var * var;

	var = static_cast<float>(rand() % 1000) * 0.001f;
	sign = static_cast<float>(rand() % 2);
	t_pos++;
	t_var++;
	if (sign == 0)
	{
		var *= -1.0;
	}
	*t_pos += *t_var * var;

	var = static_cast<float>(rand() % 1000) * 0.001f;
	sign = static_cast<float>(rand() % 2);
	t_pos++;
	t_var++;
	if (sign == 0)
	{
		var *= -1.0;
	}
	*t_pos += *t_var * var;

	var = 2.0f * static_cast<float>(rand() % 1000) * 0.001f;
	sign = static_cast<float>(rand() % 2);

	if (sign == 0)
	{
		var *= -1.0;
	}

	sc = sc * var;
}


