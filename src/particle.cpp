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

using namespace std;

std::default_random_engine gen;


//for the update step
double bivariate_normal(double x, double y, double mu_x, double mu_y, double sig_x, double sig_y) {
	return exp(-((x-mu_x)*(x-mu_x)/(2*sig_x*sig_x) + (y-mu_y)*(y-mu_y)/(2*sig_y*sig_y))) / (2.0*3.14159*sig_x*sig_y);
}


void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	// Define normal distribution for noise
	std::normal_distribution<double> dist_x(x, std[0]);
	std::normal_distribution<double> dist_y(y, std[1]);
	std::normal_distribution<double> dist_theta(theta, std[2]);

	num_particles = 20;
	weights.resize(num_particles);

	// init each particle
	for (int i = 0; i < num_particles; i++){
		Particle p;

		p.id = i;
		p.x = dist_x(gen);
		p.y = dist_y(gen);
		p.theta = dist_theta(gen);
		p.weight = 1;

		weights[i] = 1;

		particles.push_back(p);

	}

	is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/


	for (int i = 0; i < num_particles; i++) {

		if (fabs(yaw_rate) > 0.001) {

			// Update the particle position using motion model when yaw rate != 0
			Particle p;
			p.x = particles[i].x + (velocity / yaw_rate) * (sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta));
			p.y = particles[i].y + (velocity / yaw_rate) * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t));
			p.theta = particles[i].theta + yaw_rate*delta_t;

			// Define nosied distribution to add gaussian noise
			normal_distribution<double> dist_x(p.x, std_pos[0]);
			normal_distribution<double> dist_y(p.y, std_pos[1]);
			normal_distribution<double> dist_theta(p.theta, std_pos[2]);

			// Set particle state
			particles[i].x = dist_x(gen);
			particles[i].y = dist_y(gen);
			particles[i].theta = dist_theta(gen);
		}

		else {
			// Update the particle position using motion model when yaw rate = 0
			Particle p;
			p.x = particles[i].x + velocity*cos(particles[i].theta)*delta_t;
			p.y = particles[i].y + velocity*sin(particles[i].theta)*delta_t;
			p.theta = particles[i].theta + yaw_rate*delta_t;

			// Define nosied distribution to add gaussian noise
			normal_distribution<double> dist_x(p.x, std_pos[0]);
			normal_distribution<double> dist_y(p.y, std_pos[1]);
			normal_distribution<double> dist_theta(p.theta, std_pos[2]);

			// Set particle state
			particles[i].x = dist_x(gen);
			particles[i].y = dist_y(gen);
			particles[i].theta = dist_theta(gen);
		}

	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	//Nearest neighbour calculation

	double curr_dist = 0;
	for (int i = 0; i<observations.size(); i++){

		double min_dist = 999999;
		int closest_landmark_id = -1;

		for (int j = 0; j < predicted.size(); j++) {
			curr_dist = dist(observations[i].x,observations[i].y, predicted[j].x, predicted[j].y);

			if (curr_dist < min_dist){
				min_dist = curr_dist;
				closest_landmark_id = predicted[j].id;
			}
		}
		observations[i].id = closest_landmark_id;
	}

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

	// For all particles....

	for (int i=0; i<particles.size();i++){
		/*** Step 1: Transform observations from vehicle's co-ordinate system to map co-ordinate sytem***/
		vector<LandmarkObs> transformed_observations;

		for (int j=0; j<observations.size(); j++){
			LandmarkObs trans_obv;

			trans_obv.x = observations[j].x * cos(particles[i].theta) - observations[j].y * sin(particles[i].theta) + particles[i].x;
			trans_obv.y = observations[j].x * sin(particles[i].theta) + observations[j].y * cos(particles[i].theta) + particles[i].y;
			trans_obv.id = -1;

			transformed_observations.push_back(trans_obv);
		}

		/***Step 2: Keep only those landmarks which are in sensor range ***/
		vector<LandmarkObs> predicted_landmarks;

		for (int j = 0; j <map_landmarks.landmark_list.size(); j++) {
			double landmark_dist;

			landmark_dist = dist(particles[i].x, particles[i].y, map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f );

			if (landmark_dist<sensor_range){
				LandmarkObs pred_landmark;
				pred_landmark.id = map_landmarks.landmark_list[j].id_i;
				pred_landmark.x = map_landmarks.landmark_list[j].x_f;
				pred_landmark.y = map_landmarks.landmark_list[j].y_f;

				predicted_landmarks.push_back(pred_landmark);

			}
		}

		/*** Step 3: Associate observations to predicted landmarks***/
		dataAssociation(predicted_landmarks, transformed_observations);

		/*** Step 4: Calculate the weights using multi variable gaissian probability***/
		double prob = 1;
		double mvgd;

		for (int j = 0; j < predicted_landmarks.size(); j++) {
			int id_min = -1;
			double min_dist = 99999;

			double px = predicted_landmarks[j].x;
			double py = predicted_landmarks[j].y;

			for (int k = 0; k < transformed_observations.size(); k++) {
				double tx = transformed_observations[k].x;
				double ty = transformed_observations[k].y;
				double curr_dist = dist(px, py, tx, ty);

				if (curr_dist< min_dist){
					min_dist = curr_dist;
					id_min = k;
				}
			}

			if (id_min != -1){
				mvgd = bivariate_normal(px, py, transformed_observations[id_min].x, transformed_observations[id_min].y, std_landmark[0], std_landmark[1]);

				prob = prob * mvgd;
			}
		}

		weights.push_back(prob);
		particles[i].weight = prob;

	}


}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	vector<Particle> resampled_particles;

	// Create a generator to be used for generating random particle index and beta value
	default_random_engine gen;
	
	//Generate random particle index
	uniform_int_distribution<int> particle_index(0, num_particles - 1);
	
	int current_index = particle_index(gen);
	
	double beta = 0.0;
	
	double max_weight_2 = 2.0 * *max_element(weights.begin(), weights.end());
	
	for (int i = 0; i < particles.size(); i++) {
		uniform_real_distribution<double> random_weight(0.0, max_weight_2);
		beta += random_weight(gen);

	  while (beta > weights[current_index]) {
	    beta -= weights[current_index];
	    current_index = (current_index + 1) % num_particles;
	  }
	  resampled_particles.push_back(particles[current_index]);
	}
	particles = resampled_particles;

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