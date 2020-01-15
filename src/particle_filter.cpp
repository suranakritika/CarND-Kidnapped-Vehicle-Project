/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"
#include <math.h>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>
// #include "helper_functions.h" error1


using std::string;
using std::vector;
using std::normal_distribution;
using std::default_random_engine;

using namespace std;
static default_random_engine gen; //Change to static
void ParticleFilter::init(double x, double y, double theta, double std[]) { 
 
  num_particles = 100;  // TODO: Set the number of particles
  //  Normal Distribution for x, y and theta
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);
  for(int i=0; i < num_particles; i++)
  {
    Particle particle;
    particle.id = i;
    particle.x = dist_x(gen);
    particle.y = dist_y(gen);
    particle.theta = dist_theta(gen);
    particle.weight = 1.0;

    particles.push_back(particle); 
    weights.push_back(1.0); 
  }
   is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
   for (unsigned int i = 0; i < particles.size(); i++) {

    if(abs(yaw_rate) > 0.00001) {

      double predicted_x = (velocity/yaw_rate) * (sin(particles[i].theta + (yaw_rate * delta_t))- sin(particles[i].theta));
      double predicted_y = (velocity/yaw_rate) * (cos(particles[i].theta) - cos(particles[i].theta + (yaw_rate * delta_t)));
      double predicted_heading = yaw_rate * delta_t;

      std::normal_distribution<double> predict_dist_x(predicted_x, std_pos[0]);
      std::normal_distribution<double> predict_dist_y(predicted_y, std_pos[1]);
      std::normal_distribution<double> predict_dist_theta(predicted_heading, std_pos[2]);

      particles[i].x += predict_dist_x(gen);
      particles[i].y += predict_dist_y(gen);
      particles[i].theta += predict_dist_theta(gen);

    } else {

      double predicted_x = velocity * delta_t * cos(particles[i].theta);
      double predicted_y = velocity * delta_t * sin(particles[i].theta);

      std::normal_distribution<double> predict_dist_x(predicted_x, std_pos[0]);
      std::normal_distribution<double> predict_dist_y(predicted_y, std_pos[1]);
      std::normal_distribution<double> predict_dist_theta(0, std_pos[2]);

      particles[i].x += predict_dist_x(gen);
      particles[i].y += predict_dist_y(gen);
      particles[i].theta += predict_dist_theta(gen);
    }
  }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  
  for (unsigned int i = 0; i<observations.size(); i++)
  {
    double min_dist = numeric_limits<double>::max();
    int closest_landmark_id = -1;
    for (unsigned int j = 0; j < predicted.size(); j++) 
    {
      double curr_dist = dist(observations[i].x,observations[i].y, predicted[j].x, predicted[j].y);
      if (curr_dist < min_dist)
      {
        min_dist = curr_dist;
        closest_landmark_id = j;
      }
    }
    observations[i].id = closest_landmark_id;
  }
}



void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
    // For all particles....

  for (int i=0; i<num_particles; i++){
    /*** Step 1: Transform observations from vehicle's co-ordinate system to map co-ordinate sytem based on Coordinate Transformation under Rotation Link: https://www.miniphysics.com/coordinate-transformation-under-rotation.html ***/
    vector<LandmarkObs> transformed_observations;
    for (int j=0; j<observations.size(); j++){
      LandmarkObs trans_obv;
      trans_obv.x = observations[j].x * cos(particles[i].theta) - observations[j].y * sin(particles[i].theta) + particles[i].x;
      trans_obv.y = observations[j].x * sin(particles[i].theta) + observations[j].y * cos(particles[i].theta) + particles[i].y;
      trans_obv.id = j;
      transformed_observations.push_back(trans_obv);
    }
    /***Step 2: Keep only those landmarks which are in sensor range ***/
    vector<LandmarkObs> predicted_landmarks;
    LandmarkObs pred_landmark;
    for (unsigned int k = 0; k <map_landmarks.landmark_list.size(); k++) {
      double landmark_dist = dist(particles[i].x, particles[i].y, map_landmarks.landmark_list[k].x_f, map_landmarks.landmark_list[k].y_f );
      if (fabs(landmark_dist)<sensor_range){
       
        pred_landmark.id = map_landmarks.landmark_list[k].id_i;
        pred_landmark.x = map_landmarks.landmark_list[k].x_f;
        pred_landmark.y = map_landmarks.landmark_list[k].y_f;
        predicted_landmarks.push_back(pred_landmark);
      }
    }    
    /*** Step 3: Associate observations to predicted landmarks***/
    dataAssociation(predicted_landmarks, transformed_observations);
    
    /*** Step 3:  Update Weights ***/
    double particle_probability=1.0;
    double sig_x = 2 * pow(std_landmark[0], 2);
    double sig_y = 2 * pow(std_landmark[1], 2);
    double normalizer = (1/(2 * M_PI * std_landmark[0] * std_landmark[1]));
    for( unsigned int l=0; l<transformed_observations.size(); l++)
    {
      double x_obv = transformed_observations[l].x;
      double y_obv = transformed_observations[l].y;

      int idx = transformed_observations[l].id;
      double x_pred = predicted_landmarks[idx].x;
      double y_pred = predicted_landmarks[idx].y;

      particle_probability *= normalizer * exp(-(pow(x_obv - x_pred, 2)/sig_x + pow(y_obv - y_pred, 2) / sig_y ));
    }
  particles[i].weight = particle_probability;
  weights[i] = particle_probability;
  }
}

void ParticleFilter::resample() {

  vector<Particle> resampled_particles;
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


void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}