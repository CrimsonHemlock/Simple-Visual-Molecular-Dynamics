#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <assert.h>

#include "Graphics.h"  // Contains the definition of "Particule" and "Graphics"

const double PETIT = 1e-10;

enum  col_type {
  bottom,
  right,
  top,
  left,
  animation,
  particle
}; //different types of collision

typedef struct { // a structure describing each collision -- one might want an array of Events
  enum col_type type;
  int ia,ib;
  double time;
} Event;

void initparticles( Particle *p, int np, double Lmin, double Lmax, double diameter){ // Initialise the position and speed of the particles
  int i;
  for( i=0;i<np;i++){ // drand48() -- random double numbers uniform on (0 1)
    p[i].x = Lmin +diameter/2 + (Lmax-Lmin-diameter)*drand48(); //random positions for intial condition using drand48()
    p[i].y = Lmin +diameter/2 + (Lmax-Lmin-diameter)*drand48();
    p[i].vx = 10*(drand48()-0.5) ;// choose random speeds too using drand48();
    p[i].vy = 10*(drand48()-0.5);
  }
}
int HitWall(Event *e,int Np,Particle *p,double Lmin, double Lmax, double diameter){ // Calculate the times until collisions with the walls
	int i, compteur;

	compteur=0;
	for(i=0;i<Np;i++){
	if( p[i].vy > 0 ){
		e[compteur].type=top; // Collision with bottom wall
		e[compteur].ia=i;
		e[compteur].ib=-1; 
		e[compteur].time=(Lmax-p[i].y-diameter/2.)/(p[i].vy); // Add event to the table of events
		compteur=compteur+1; // Update the events counter
}

	if(p[i].vy <0){
		e[compteur].type=bottom; // Collision with top wall
		e[compteur].ia=i;
		e[compteur].ib=-1; 
		e[compteur].time=-(p[i].y-Lmin-diameter/2.)/(p[i].vy);// Add event to the table of events
		compteur=compteur+1;// Update the events counter
}

	if( p[i].vx > 0 ){
		e[compteur].type=right; // Collision with right wall
		e[compteur].ia=i;
		e[compteur].ib=-1; 
		e[compteur].time=(Lmax-p[i].x-diameter/2.)/(p[i].vx);// Add event to the table of events
		compteur=compteur+1;// Update the events counter
	}

	if (p[i].vx <0){
		e[compteur].type=left; // Collision with left wall
		e[compteur].ia=i;
		e[compteur].ib=-1; 
		e[compteur].time=-(p[i].x-Lmin-diameter/2.)/(p[i].vx);// Add event to the table of events
		compteur=compteur+1;// Update the events counter
}

	}
	return compteur;
}
int ParticleCollision(Event *e,int Np,Particle *p, double diameter, int compteur){ // Calculates times until collisions between the different particles
	int i,j;
	int compt;
	double A,B,C,D,t,T;
	compt=compteur;
	
	for(i=0;i<Np;i++){
		for (j=0;j<i;j++){
				A=pow((p[i].vx-p[j].vx),2.0)+pow((p[i].vy-p[j].vy),2.0);
				B=2*((p[i].x-p[j].x)*(p[i].vx-p[j].vx)+(p[i].y-p[j].y)*(p[i].vy-p[j].vy));
				C=pow((p[i].x-p[j].x),2.0)+pow((p[i].y-p[j].y),2.0)-pow(diameter,2.0);
				D=(pow(B,2.0)-4*A*C);
				t=(-B-sqrt(D))/(2*A);
				T=(-B+sqrt(D))/(2*A); // Solving of the second degree aquation evolving the time until collision
				if (t>0){
					e[compt].type=particle;
					e[compt].ia=i;
					e[compt].ib=j;
					e[compt].time=t; // Add event to the table of events
					compt++; // Increment counter of events
			}
		}
	}
	
	return compt;
}
int FindMin(Event *e, int compteur){ // Find the first event in the event matrix by looking at the minimum time until next events
	double tmin;
	int i,indmin;
	tmin=1e+10;
	indmin=-1;
	for(i=0;i<compteur;i++){
		if(e[i].time<tmin && e[i].time>0){
			tmin=e[i].time;
			indmin=i; 
		}
	}
	return indmin;
}

void MoveParticleWall(Event* e, Particle *p,int indmin, int Np){ // Moves the particles and changes the velocity of the particle that collides with the wall
	int i;
	for(i=0;i<Np;i++){
		p[i].x += (p[i].vx)*(e[indmin].time); 
    	p[i].y += (p[i].vy)*(e[indmin].time) ; // Positions updates
	}
	if(e[indmin].type==top || e[indmin].type==bottom){
		p[e[indmin].ia].vy *= -1;
	}
	if(e[indmin].type==right || e[indmin].type==left){
		p[e[indmin].ia].vx *= -1; // Change the velocity of the particle involved in the collision with a wall
	}
}
void MoveParticleCollision(Event* e, Particle *p, int Np, int indmin, double diameter){ // Moves the particles and changes the velocity of the particle that collides with another particle
	int a,b, i;
	double d, rx, ry, ps;
	
	for(i=0;i<Np;i++){
		p[i].x += (p[i].vx)*(e[indmin].time); 
    	p[i].y += (p[i].vy)*(e[indmin].time) ; // Positions updates of all the particles
	}
	a=e[indmin].ia;
	b=e[indmin].ib;
	d = sqrt(pow(p[a].x-p[b].x, 2.0)+pow(p[a].y-p[b].y, 2.0));
	assert(fabs(d-diameter) < PETIT); // Error test
	
	rx = (p[b].x-p[a].x)/d;
	ry = (p[b].y-p[a].y)/d;
	ps = rx*(p[b].vx-p[a].vx)+ ry*(p[b].vy-p[a].vy);
	p[a].vx += ps*rx;
	p[a].vy += ps*ry;
	p[b].vx -= ps*rx;
	p[b].vy -= ps*ry; // Update of the velocities of the particles involved in the collision
}
void Animate(Event* e, Particle *p, int Np ,int indmin){
	int i;
	for(i=0;i<Np;i++){
		p[i].x += (p[i].vx)*e[indmin].time; 
    	p[i].y +=(p[i].vy)*e[indmin].time ; // Positions updates for the event type animation
    	
	}
}

int main(){
	printf("start the main program\n");
	double Tnow=0; //Create a universal time for the simulation
	int Np=10; //Number of particles
  	double Dt=1./32.; //Time interval limit 1./16.    
  	double diameter=1;//particle size
  	int Pix=800; //Number of pixels for window
  	double Lmax=10, Lmin=0; //Physical dimensions of box
  	int sizee; // Current size of the event matrix
  	Graphics gw(Np,Pix, Lmin ,Lmax,diameter);// Open a window to plot particles in
  	srand48(1);//inititalize random numbers -- to find always the same value // you can replace "1" by time(NULL) 
  	Particle *p= (Particle *) malloc( Np *sizeof( Particle)); //an array of particles
  	initparticles(p,Np,Lmin, Lmax,diameter); //place particles in box
  	Event *e = (Event *) malloc( (4*Np+(Np*(Np-1))/2+1)* sizeof(Event) ); // 4 Np possible collisions with walls
  
  	for (int l=0; l<100000;l++){//long loop
    	sizee=HitWall(e,Np,p,Lmin,Lmax,diameter);//find future collision with wall using exact calculations of collision times
    	int indmin,k;
    	e[sizee].type=animation;
    	e[sizee].ia=-1;
    	e[sizee].ib=-1;
    	e[sizee].time=Dt-fmod(Tnow,Dt);
    	sizee++; // add the event type = animation in the event table
    	sizee=ParticleCollision(e,Np,p,diameter,sizee); // Find future collisions with other particules
		indmin=FindMin(e,sizee);// find the very first collision in the future by looking at e[].time
		if(e[indmin].type==animation){ 		
    		Animate(e,p,Np,indmin); // Move particles if next event is a type animation
    		gw.draw(p); // Update the GUI
    		Tnow+=e[indmin].time; // Update the universal time
    	}
    	else if(e[indmin].type==particle){
    		MoveParticleCollision(e,p,Np, indmin, diameter);// Move particles if next collision is with another particles and change speeds of the particles involved
    		Tnow+=e[indmin].time; // Update the universal time
    	}
    	else{
    		MoveParticleWall(e,p,indmin,Np);// move particles if next collision is with a wall and change speeds of the particle involved
    		Tnow+=e[indmin].time; // Update the universal time
    	}
  }

	free(p);
	free(e); // Freeing the memory
  
  return 0;
}


