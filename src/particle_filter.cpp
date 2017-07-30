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

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  // Initialize random number generator to generator numbers from random distribution
  // http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
  std::random_device rd;
  std::mt19937 gen(rd());
  std::normal_distribution<> g_x(0, std[0]);
  std::normal_distribution<> g_y(0, std[1]);
  std::normal_distribution<> g_t(0, std[2]);

  std::cout << "Initializing particle filter with estimated coordinates x: "<< x << " y: " << y << " theta: " << theta << std::endl;

  // clear particle and weights vector
  particles.clear();
  weights.clear();

  // Set the number of particles.
  num_particles = 10;


  // Initialize all particles to first position (based on estimates of
  //   x, y, theta and their uncertainties from GPS) and all weights to 1.
  // Add random Gaussian noise to each particle.
  const double w = 1;
  for (int i = 0; i < num_particles; i++) {
      Particle p;
      p.id = i;
      p.x = x + g_x(gen);
      p.y = y + g_y(gen);
      p.theta = theta + g_t(gen);
      p.weight = w;
      particles.push_back(p);
      weights.push_back(w);
    }

  // print out particle coordinates
  if (false) {
      std::cout << "Initialized particles with random gaussian noise" << std::endl;
      for (auto p : particles) {
          std::cout << "[" << p.id << "]";
          std::cout <<" x:" << p.x << " y: " << p.y << " theta: " << p.theta << " weight: " << p.weight << std::endl;
        }
    }

  is_initialized = true;
  std::cout << "Particle filter initialized!" << std::endl;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
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

  double v_per_yr = velocity / yaw_rate;
  double yr_dt = yaw_rate * delta_t;

  // Predict particles new position by using velocity and yaw_rate measurements
  // and then add random Gaussian noise.
  for (auto p: particles) {
    // print out particle coordinates
    if (false) {
      std::cout << "[" << p.id << "]";
      std::cout <<" x:" << p.x << " y: " << p.y << " theta: " << p.theta << " weight: " << p.weight << std::endl;
      }
    // 1. Predict new position and angle
    p.x = p.x + v_per_yr * (sin(p.theta+yr_dt)-sin(p.theta));
    p.y = p.y + v_per_yr * (cos(p.theta) - cos(p.theta+yr_dt));
    p.theta = p.theta + yr_dt;
    // 2. Add noise to predicted values
    p.x += g_x(gen);
    p.y += g_y(gen);
    p.theta += g_t(gen);

    // print out particle coordinates
    if (false) {
      std::cout << "[" << p.id << "]";
      std::cout <<" x:" << p.x << " y: " << p.y << " theta: " << p.theta << " weight: " << p.weight << std::endl;
      }
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
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

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
