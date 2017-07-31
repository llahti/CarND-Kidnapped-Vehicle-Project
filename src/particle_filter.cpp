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
#include <random>
#include "map.h"

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(const double x, const double y,
                          const double theta, const double std[]) {
  // Initialize random number generator to generator numbers from random distribution
  // http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
  std::random_device rd;
  std::mt19937 gen(rd());
  std::normal_distribution<> g_x(x, std[0]);  // mean = x, std_dev = std[0]
  std::normal_distribution<> g_y(y, std[1]);
  std::normal_distribution<> g_t(theta, std[2]);

  // clear particle and weights vector
  particles.clear();
  weights.clear();

  // Set the number of particles.
  num_particles = 100;

  // Initialize all particles to first position (based on estimates of
  // x, y, theta and their uncertainties from GPS) and all weights to 1.
  // Add random Gaussian noise to each particle.
  const double w = 1;
  for (int i = 0; i < num_particles; i++) {
      Particle p;
      p.id = i;
      p.x = g_x(gen);  // Generate values from gaussian distribution
      p.y = g_y(gen);
      p.theta = g_t(gen);
      p.weight = w;
      particles.push_back(p);
      weights.push_back(w);
    }

  is_initialized = true;
}

void ParticleFilter::prediction(const double delta_t, const double std_pos[],
                                const double velocity, const double yaw_rate) {
  // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
  //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
  //  http://www.cplusplus.com/reference/random/default_random_engine/
  //
  // Initialize random number generator to generator numbers from random distribution
  // http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
  std::random_device rd;
  std::mt19937 gen(rd());
  std::normal_distribution<> g_x(0, std_pos[0]);
  std::normal_distribution<> g_y(0, std_pos[1]);
  std::normal_distribution<> g_t(0, std_pos[2]);

  const double epsilon = 0.00001;
  // variables for storing common calculation results
  double v_per_yr;
  double yr_dt;
  double vdt;
  // Prevent division by zero
  const double abs_yaw_rate = abs(yaw_rate);
  if (abs_yaw_rate < epsilon) {
    vdt = velocity * delta_t;
  }
  else {
    v_per_yr = velocity / yaw_rate;
    yr_dt = yaw_rate * delta_t;
  }
  // Predict particles new position by using velocity and yaw_rate measurements
  // and then add random Gaussian noise.
  for (Particle& p: particles) {
    // 1. Predict new position and angle
    // Prevent division by zero yaw_rate
    if (abs_yaw_rate<epsilon) {
      p.x += vdt * cos(p.theta);
      p.y += vdt * sin(p.theta);
      // Theta is not changed when yaw_rate is close to zero,
      // we assume that change would be very miniscule and therefore ignored
    } else {
      p.x += v_per_yr * (sin(p.theta+yr_dt) - sin(p.theta));
      p.y += v_per_yr * (cos(p.theta)       - cos(p.theta+yr_dt));
      p.theta = p.theta + yr_dt;
    }
    // 2. Add noise to predicted values
    p.x += g_x(gen);
    p.y += g_y(gen);
    p.theta += g_t(gen);
    }
}


void ParticleFilter::dataAssociation(const std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
  // Vector size is needed several times so save it in variable
  const unsigned obs_size = observations.size();
  const unsigned predicted_size = predicted.size();

  // confirm that size of both vectors is more than 0
  if ((0 < predicted_size) && (0 < obs_size)) {
    // Outer loop: Loop through each observations
    for (int i = 0; i < obs_size; i++) {
      double shortest_dist = dist(predicted[0].x, predicted[0].y, observations[i].x, observations[i].y);
      int landmark_id = predicted[0].id;
      // Inned loop to go through each predicted landmark, calculate distance betweeb
      // observation and landmark and select the shortest distance
      for (int j = 1; j < predicted_size; j++) {
        const double distance = dist(predicted[j].x, predicted[j].y, observations[i].x, observations[i].y);
        if (distance < shortest_dist) {
            shortest_dist = distance;
            landmark_id = predicted[j].id;
          }
        }
      observations[i].id = landmark_id;
    }
  }
  else {cout << "NO DATA TO ASSOCIATE! " << endl; }
}

void ParticleFilter::updateWeights(const double sensor_range, const double std_landmark[],
                const std::vector<LandmarkObs> observations, const Map map_landmarks) {
  // Update the weights of each particle using a mult-variate Gaussian distribution. You can read
  //   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
  // NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
  //   according to the MAP'S coordinate system. You will need to transform between the two systems.
  //   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
  //   The following is a good resource for the theory:
  //   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
  //   and the following is a good resource for the actual equation to implement (look at equation
  //   3.33
  //   http://planning.cs.uiuc.edu/node99.html

  const unsigned size = observations.size();
  std::vector<LandmarkObs> transformed_obs{size};

  for (Particle& p:particles) {
    //
    // 1. Transform observation coordinates
    //
    transform(observations, p, transformed_obs);

    // 2. Find landmarks within the sensor range
    std::vector<LandmarkObs> landmarks_near_particle;
    for (const auto lm:map_landmarks.landmark_list) {
      // Go through each landmark and add to vector if within sensor range
      if (dist(lm.x_f, lm.y_f, p.x, p.y) <= sensor_range) {
        LandmarkObs lmark;
        lmark.id = lm.id_i;
        lmark.x = lm.x_f;
        lmark.y = lm.y_f;
        landmarks_near_particle.push_back(lmark);
        }
      }

    //
    // 2. Find associations
    //
    dataAssociation(landmarks_near_particle, transformed_obs);
    // Extract landmark association information into below vectors
    // and save associations to particle
    p.associations.clear();
    p.sense_x.clear();
    p.sense_y.clear();
    for (const auto l_obs:transformed_obs) {
      p.associations.push_back(l_obs.id);
      p.sense_x.push_back(l_obs.x);
      p.sense_y.push_back(l_obs.y);
    }

    //
    // 3. Update weight
    //
    double weight = 1;
    const double std_x = std_landmark[0];
    const double std_y = std_landmark[1];
    const double p1 = 1/(2 * M_PI * std_x * std_y);
    const double pow_stdx = (2 * pow(std_x, 2));
    const double pow_stdy = (2 * pow(std_y, 2));
    for (int i = 0; i < p.associations.size(); i++) {
      const int id = p.associations[i];
      const double mu_x = p.sense_x[i];
      const double mu_y = p.sense_y[i];
      // This loop will be stopped when matchin id is found
      for (const auto lmark:landmarks_near_particle) {
          if (id == lmark.id) {  // found
              const double dx = lmark.x - mu_x;
              const double dy = lmark.y - mu_y;
              const double p2 = pow(dx, 2) / pow_stdx;
              const double p3 = pow(dy, 2) / pow_stdy;
              const double p_exp = exp(-(p2 + p3));
              const double P =  p1 * p_exp;
              // Update total weight
              weight *= P;
              break;  // continue to next association
            }
        }
    }
    // Note. particle filter weight vector need to be updated separately
    p.weight = weight;
  }

  // Update all weights to weights vector (needed in resampling step)
  weights.clear();
  for (const auto p:particles) {
      weights.push_back(p.weight);
  }
}

void ParticleFilter::resample() {
  // TODO: Resample particles with replacement with probability proportional to their weight.
  // NOTE: You may find std::discrete_distribution helpful here.
  //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
  std::random_device rd;
  std::mt19937 gen(rd());
  std::discrete_distribution<int> d(weights.begin(), weights.end());

  std::vector<Particle> resampled_particles;

  for (int i=0; i < particles.size(); i++) {
      Particle new_particle = particles[d(gen)];
      new_particle.id = i;
      resampled_particles.push_back(new_particle);
    }
  // Update particles with resampled set
  particles.clear();
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
