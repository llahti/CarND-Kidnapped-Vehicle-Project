/*
 * particle_filter.h
 *
 * 2D particle filter class.
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#ifndef PARTICLE_FILTER_H_
#define PARTICLE_FILTER_H_

#include "helper_functions.h"

struct Particle {

	int id;
	double x;
	double y;
	double theta;
	double weight;
	std::vector<int> associations;
	std::vector<double> sense_x;
	std::vector<double> sense_y;
};



class ParticleFilter {
	
	// Number of particles to draw
	int num_particles; 
	
	
	
	// Flag, if filter is initialized
	bool is_initialized;
	
	// Vector of weights of all particles
	std::vector<double> weights;
	
public:
	
	// Set of current particles
	std::vector<Particle> particles;

	// Constructor
	// @param M Number of particles
	ParticleFilter() : num_particles(0), is_initialized(false) {}

	// Destructor
	~ParticleFilter() {}

	/**
	 * init Initializes particle filter by initializing particles to Gaussian
	 *   distribution around first position and all the weights to 1.
	 * @param x Initial x position [m] (simulated estimate from GPS)
	 * @param y Initial y position [m]
	 * @param theta Initial orientation [rad]
	 * @param std[] Array of dimension 3 [standard deviation of x [m], standard deviation of y [m]
	 *   standard deviation of yaw [rad]]
	 */
	void init(const double x, const double y,
		  const double theta, const double std[]);

	/**
	 * prediction Predicts the state for the next time step
	 *   using the process model.
	 * @param delta_t Time between time step t and t+1 in measurements [s]
	 * @param std_pos[] Array of dimension 3 [standard deviation of x [m], standard deviation of y [m]
	 *   standard deviation of yaw [rad]]
	 * @param velocity Velocity of car from t to t+1 [m/s]
	 * @param yaw_rate Yaw rate of car from t to t+1 [rad/s]
	 */
	void prediction(const double delta_t, const double std_pos[],
			const double velocity, const double yaw_rate);
	
	/**
	 * dataAssociation Finds which observations correspond to which landmarks (likely by using
	 *   a nearest-neighbors data association).
	 * @param predicted Vector of predicted landmark observations
	 * @param observations Vector of measured? landmark observations
	 */
	void dataAssociation(const std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations);
	
	/**
	 * updateWeights Updates the weights for each particle based on the likelihood of the 
	 *   observed measurements. 
	 * @param sensor_range Range [m] of sensor
	 * @param std_landmark[] Array of dimension 2 [standard deviation of range [m],
	 *   standard deviation of bearing [rad]]
	 * @param observations Vector of landmark observations
	 * @param map Map class containing map landmarks
	 */
	void updateWeights(const double sensor_range, const double std_landmark[],
			   const std::vector<LandmarkObs> observations, const Map map_landmarks);
	
	/**
	 * resample Resamples from the updated set of particles to form
	 *   the new set of particles.
	 */
	void resample();

	/*
	 * Set a particles list of associations, along with the associations calculated world x,y coordinates
	 * This can be a very useful debugging tool to make sure transformations are correct and assocations correctly connected
	 */
	Particle SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y);
	
	std::string getAssociations(Particle best);
	std::string getSenseX(Particle best);
	std::string getSenseY(Particle best);

	/**
	 * initialized Returns whether particle filter is initialized yet or not.
	 */
	const bool initialized() const {
		return is_initialized;
	}

	/**
	 * @brief transform transforms given coordites
	 * @param obs landmark observation from car's coordinate system
	 * @param p particle coordinates and angle to which given observation is transformed (note that struct is type of ground_truth)
	 * @param x_g transformed x-coordinate (in global map coordinates)
	 * @param y_g transformed y-coordinate (in global map coordinates)
	 */
	//inline void transform(const LandmarkObs& obs, const ground_truth& p, double& x_t, double& y_t) {
	inline void transform(const std::vector<LandmarkObs>& observations,
			      const Particle& p,
			      std::vector<LandmarkObs>& transformation) {
	  for (int i = 0; i < observations.size(); i++) {
	    const double x = observations[i].x;
	    const double y = observations[i].y;
	    const double xt = p.x;
	    const double yt = p.y;
	    const double theta = p.theta;
	    // Preserve observation id in transformed coordinates
	    transformation[i].id = observations[i].id;
	    transformation[i].x = (x * cos(theta)) - (y * sin(theta)) + xt;
	    transformation[i].y = (x * sin(theta)) + (y * cos(theta)) + yt;
	  }
	  // std::cout << "Transformed x: " << obs.x << " y: " << obs.y << std::endl;
	}
};



#endif /* PARTICLE_FILTER_H_ */
