/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
using std::normal_distribution;

using namespace std;

static default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 100;  // TODO: Set the number of particles

  normal_distribution<double> x_dist(0, std[0]);
  normal_distribution<double> y_dist(0, std[1]);
  normal_distribution<double> theta_dist(0, std[2]);

  for (int i =0; i<num_particles; i++) {
    Particle particle;
    particle.id = i;
    particle.x = x + x_dist(gen);
    particle.y = y + y_dist(gen);
    particle.theta = theta + theta_dist(gen);
    particle.weight = 1.0;

    particles.push_back(particle);
  }

  is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */

  normal_distribution<double> norm_x(0, std_pos[0]);
  normal_distribution<double> norm_y(0, std_pos[1]);
  normal_distribution<double> norm_theta(0, std_pos[2]);

  for (int i = 0; i < num_particles; i++) {

    if (fabs(yaw_rate) < 0.00001) {
      particles[i].x += (velocity * delta_t * cos(particles[i].theta)) + norm_x(gen);
      particles[i].y += (velocity * delta_t * sin(particles[i].theta));
      particles[i].theta += norm_theta(gen);
    } 
    else {
      particles[i].x += (velocity / yaw_rate * (sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta))) + norm_x(gen);
      particles[i].y += (velocity / yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t))) + norm_y(gen);
      particles[i].theta += (yaw_rate * delta_t) + norm_theta(gen);
    }
  }

}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */

  for (auto &observation: observations) {

    double minimum_dist = numeric_limits<double>::max();

    observation.id = -1;

    for (auto &pred: predicted ) {
      double distance = dist(pred.x, pred.y, observation.x, observation.y);

      if (distance < minimum_dist) {
        minimum_dist = distance;
        observation.id = pred.id;
      }

    }
  }

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */

  for (int i = 0; i<particles.size(); i++) {

    vector<LandmarkObs> trans_observations(observations.size());

    for (int j = 0; j<observations.size(); j++) {
      trans_observations[j].x = particles[i].y + cos(particles[i].theta) * observations[j].x - sin(particles[i].theta) * observations[j].y;
      trans_observations[j].y = particles[i].y + sin(particles[i].theta) * observations[j].x + cos(particles[i].theta) * observations[j].y;
      trans_observations[j].id = -1
    }

    vector<LandmarkObs> landmarks;

    for (int j = 0; j < map_landmarks.landmark_list.size(); j++) {

      if (dist(particles[i].x, particles[i].y, map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f) <= sensor_range) {
        landmarks.push_back(LandmarkObs{map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f, map_landmarks.landmark_list[j].id_i});
      }
    }

    dataAssociation(landmarks, trans_observations);

    particles[i].weight = 1.0;

    vector<double> observation_probabilities(trans_observations.size());
    for (int j = 0; j < observations.size(); j++) {
      LandmarkObs closest_landmark;

      closest_landmark.id = -1;
      closest_landmark.x = static_cast<double>(map_landmarks.landmark_list[trans_observations[j].id - 1].x_f);
      closest_landmark.y = static_cast<double>(map_landmarks.landmark_list[trans_observations[j].id - 1].y_f);

      observation_probabilities[i] = (1 / (2 * M_PI * std_landmark[0] * std_landmark[1])) *
                                     exp(-(pow(trans_observations[j].x - closest_landmark.x, 2.0) / 
                                     (2 * pow(std_landmark[0], 2.0)) + pow(trans_observations[j].y - closest_landmark.y, 2.0) / 
                                     (2 * pow(std_landmark[1], 2.0))));

      particles[i].weight *= observation_probabilities[j];
    }

    // set weights
    weights[i] = particles[i].weight;
  }

}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

  vector<Particle> updated_particles;

  uniform_int_distribution<int> uniformintdist(0, num_particles-1);
  auto idx = uniformintdist(gen);

  vector<double> weights;
  for (int i = 0; i < num_particles; i++) {
    weights.push_back(particles[i].weight);
  }

  double max_weight = *max_element(weights.begin(), weights.end());

  uniform_real_distribution<double> uniformrealdist(0.0, max_weight);

  double beta = 0.0;
  for (int i = 0; i < num_particles; i++) {
    beta += uniformrealdist(gen) * 2.0;
    while (beta > weights[idx]) {
      beta -= weights[idx];
      idx = (idx + 1) % num_particles;
    }
    updated_particles.push_back(particles[idx]);
  }

  particles = updated_particles;


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