/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"
#include "helper_functions.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	// random noise generator
    default_random_engine gen;

	num_particles = 200;

 	// Resize weights vector based on num_particles
//  	weights.resize(num_particles);
    
  	// Resize vector of particles
//  	particles.resize(num_particles);

    normal_distribution<double> N_x_init(x, std[0]);
    normal_distribution<double> N_y_init(y, std[1]);
    normal_distribution<double> N_theta_init(theta, std[2]);

	for (int i = 0; i < num_particles; i++) {
		Particle p;
		p.id = i;
		p.x = N_x_init(gen);
		p.y = N_y_init(gen);
		p.theta = N_theta_init(gen);
		p.weight = 1;

		particles.push_back(p);
		weights.push_back(1);
	}

	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	default_random_engine gen;

	normal_distribution<double> N_x(0, std_pos[0]);
	normal_distribution<double> N_y(0, std_pos[1]);
	normal_distribution<double> N_theta(0, std_pos[2]);
		
  for(int i = 0; i<num_particles; i++){
    double theta_new = particles[i].theta + yaw_rate * delta_t;
    
    if(fabs(yaw_rate) > 0.00001){
      particles[i].x += velocity/yaw_rate * (sin(theta_new)-sin(particles[i].theta));
      particles[i].y += velocity/yaw_rate * (cos(particles[i].theta)-cos(theta_new));
      particles[i].theta += yaw_rate * delta_t;
    }else{
      particles[i].x += velocity * sin(particles[i].theta) * delta_t;
      particles[i].y += velocity * cos(particles[i].theta) * delta_t;
      particles[i].theta = particles[i].theta;      
    }

    particles[i].x += N_x(gen);    
    particles[i].y += N_y(gen);    
    particles[i].theta += N_theta(gen);    
  }
  
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {

	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
	weights.clear();
  for(int i = 0; i<num_particles; i++){
    std::vector<LandmarkObs> predicted;
  
    double x_p = particles[i].x;
    double y_p = particles[i].y;
    double theta_p = particles[i].theta;
    
    //Transform car observations to map coordinates supposing that the particle is the car.
    particles[i].associations.clear();
    particles[i].sense_x.clear();
    particles[i].sense_y.clear();
	double weight = 1.0;
    for(int j = 0; j<observations.size(); j++){
      double o_x = observations[j].x;
      double o_y = observations[j].y;
      double o_x_map = o_x * cos(theta_p) - o_y * sin(theta_p) + x_p;
      double o_y_map = o_x * sin(theta_p) + o_y * cos(theta_p) + y_p;

      if(pow(pow(o_x_map-x_p,2)+pow(o_y_map-y_p,2),0.5) > sensor_range) continue;
      	particles[i].sense_x.push_back(o_x_map);
      	particles[i].sense_y.push_back(o_y_map);
      	double min_range = numeric_limits<double>::max();
      	int min_k=-1;
      	for(int k = 0; k<map_landmarks.landmark_list.size(); k++){
        	double l_x = map_landmarks.landmark_list[k].x_f;
        	double l_y = map_landmarks.landmark_list[k].y_f;       
        	double diff_x = l_x - o_x_map;
        	double diff_y = l_y - o_y_map;
        	double range = pow(pow(diff_x,2)+pow(diff_y,2),0.5);
        	if(range < min_range){
          	min_range = range;
          	min_k = k;
        	}
      	}
      	double l_x = map_landmarks.landmark_list[min_k].x_f;
      	double l_y = map_landmarks.landmark_list[min_k].y_f;

      	particles[i].associations.push_back(map_landmarks.landmark_list[min_k].id_i);

		double s_x = std_landmark[0];
        double s_y = std_landmark[1];
      	weight = weight * exp(-0.5 * (pow((l_x - o_x_map) / s_x,2) + pow((l_y - o_y_map) / s_y,2))) / (2*M_PI*s_x*s_y);
    	} 
    particles[i].weight = weight;
    weights.push_back(weight); 
  }		
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    default_random_engine gen;
    discrete_distribution<int> distribution(weights.begin(), weights.end());
    std::vector<Particle> new_particles;
    
    weights.clear();
    
    for(int i=0; i < num_particles; i++){
        int chosen = distribution(gen);
        new_particles.push_back(particles[chosen]);
        weights.push_back(particles[chosen].weight);
    }
    
    particles = new_particles;
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
