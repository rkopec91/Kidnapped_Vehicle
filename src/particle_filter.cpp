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

static std::default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {

  num_particles = 100;

  std::normal_distribution<double> x_dist(x, std[0]);
  std::normal_distribution<double> y_dist(y, std[1]);
  std::normal_distribution<double> theta_dist(theta, std[2]);

  for (int i =0; i<num_particles; i++) {
    Particle particle;
    particle.id = i;
    particle.x = x_dist(gen);
    particle.y = y_dist(gen);
    particle.theta = theta_dist(gen);
    particle.weight = 1.0;

    particles.push_back(particle);
  }

  is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[],
                                double velocity, double yaw_rate) {

  std::normal_distribution<double> norm_x(0, std_pos[0]);
  std::normal_distribution<double> norm_y(0, std_pos[1]);
  std::normal_distribution<double> norm_theta(0, std_pos[2]);

  for (int i = 0; i < num_particles; i++) {

    if (fabs(yaw_rate) < 0.00001) {
      particles[i].x += (velocity * delta_t * cos(particles[i].theta)) + norm_x(gen);
      particles[i].y += (velocity * delta_t * sin(particles[i].theta)) + norm_y(gen);
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

  for (int i = 0; i < observations.size(); i++) {

    double minimum_distance = std::numeric_limits<double>::max();

    observations[i].id = -1;

    for (int j = 0; j < predicted.size(); j++ ) {
      double distance = dist(predicted[j].x, predicted[j].y, observations[i].x, observations[i].y);

      if (distance < minimum_distance) {
        minimum_distance = distance;
        observations[i].id = predicted[j].id;
      }

    }
  }

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   const vector<LandmarkObs> &observations,
                                   const Map &map_landmarks) {


  for (int i = 0; i < num_particles; i++) {
    // get landmarks within range
    vector<LandmarkObs> close_landmarks;
    for(unsigned int j = 0; j < map_landmarks.landmark_list.size(); j++) {
      double delta_x = particles[i].x - map_landmarks.landmark_list[j].x_f;
      double delta_y = particles[i].y - map_landmarks.landmark_list[j].y_f;
      if ( pow(delta_x,2.0) + pow(delta_y,2.0) <= pow(sensor_range,2.0) ) {
        close_landmarks.push_back(LandmarkObs{ map_landmarks.landmark_list[j].id_i,
                                                map_landmarks.landmark_list[j].x_f,
                                                map_landmarks.landmark_list[j].y_f });
      }
    }

    // Transform observation coordinates.
    vector<LandmarkObs> observed_marks;
    for(unsigned int j = 0; j < observations.size(); j++) {
      double x = cos(particles[i].theta)*observations[j].x - sin(particles[i].theta)*observations[j].y + particles[i].x;
      double y = sin(particles[i].theta)*observations[j].x + cos(particles[i].theta)*observations[j].y + particles[i].y;
      observed_marks.push_back(LandmarkObs{ observations[j].id, x, y });
    }

    // Observation association to landmark.
    dataAssociation(close_landmarks, observed_marks);

    // reseting the weights
    particles[i].weight = 1.0;
    // compute weights
    for(unsigned int j = 0; j < observed_marks.size(); j++) {

      double x_landmark;
      double y_landmark;
      unsigned int n = 0;
      unsigned int nLandmarks = close_landmarks.size();
      bool looking = true;
      while( looking && n < nLandmarks ) {
        if ( close_landmarks[n].id == observed_marks[j].id) {
          looking = false;
          x_landmark = close_landmarks[n].x;
          y_landmark = close_landmarks[n].y;
        }
        n++;
      }

      // find deltas
      double delta_x = observed_marks[j].x - x_landmark;
      double delta_y = observed_marks[j].y - y_landmark;

      // calculate new weight
      double weight = ( 1/(2*M_PI*std_landmark[0]*std_landmark[1])) * exp( -( pow(delta_x,2.0)/(2*pow(std_landmark[0],2.0)) + (pow(delta_y,2.0)/(2*pow(std_landmark[1],2.0))) ) );
      // update weight for particle
      if (weight == 0) {
        particles[i].weight *= 0.00001;
      } else {
        particles[i].weight *= weight;
      }
    }
  }

}

void ParticleFilter::resample() {

  vector<Particle> updated_particles;

  vector<double> weights;
  for (int i = 0; i < num_particles; i++) {
    weights.push_back(particles[i].weight);
  }

  std::uniform_int_distribution<int> uniform_int_dist(0, num_particles-1);
  auto idx = uniform_int_dist(gen);

  double max_weight = *max_element(weights.begin(), weights.end());

  std::uniform_real_distribution<double> uniform_real_dist(0.0, max_weight);

  double beta = 0.0;
  for (int i = 0; i < num_particles; i++) {
    beta += uniform_real_dist(gen) * 2.0;
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