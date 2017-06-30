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

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  // TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
  //   x, y, theta and their uncertainties from GPS) and all weights to 1. 
  // Add random Gaussian noise to each particle.
  // NOTE: Consult particle_filter.h for more information about this method (and others in this file).
  // default_random_engine gen;
  random_device seed;
  mt19937 gen(seed());
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);
  
  num_particles = 20;

  //increasing the particles vector to num_particles 
  //increasing the weights vector too


  for (int i = 0; i < num_particles; i++)
  {
    Particle p;
    p.id = i;
    p.x = dist_x(gen);
    p.y = dist_y(gen);
    p.theta = dist_theta(gen);
    p.weight = 1.0;
    //save each p into the vector particles
    particles.push_back(p);
    weights.push_back(1.0);
  }

  
  is_initialized = true;
  // cout << "finish init" << endl;

  // for (int i = 0; i < num_particles; i++)
  // {
  //   cout << "id: " << particles[i].id << " , " << "x: " << particles[i].x << " , " << "y: " << particles[i].y << " , " << "theta: " << particles[i].theta << endl;
  // }
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
  // TODO: Add measurements to each particle and add random Gaussian noise.
  // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
  //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
  //  http://www.cplusplus.com/reference/random/default_random_engine/
  // default_random_engine gen;
  random_device seed;
  mt19937 gen(seed());
  double pred_x, pred_y, pred_theta;

  for (int i = 0; i < num_particles; i++)
  {
    if(yaw_rate == 0)
    {
      pred_x = particles[i].x + velocity * delta_t * cos(particles[i].theta);
      pred_y = particles[i].y + velocity * delta_t * sin(particles[i].theta);
      pred_theta = particles[i].theta;
    }
    else
    {
      pred_theta = particles[i].theta + yaw_rate * delta_t;
      pred_x = particles[i].x + velocity / yaw_rate * (sin(pred_theta) - sin(particles[i].theta));
      pred_y = particles[i].y + velocity / yaw_rate * (cos(particles[i].theta) - cos(pred_theta));

    }
    // add noise 
  normal_distribution<double> dist_x(pred_x, std_pos[0]);
  normal_distribution<double> dist_y(pred_y, std_pos[1]);
  normal_distribution<double> dist_theta(pred_theta, std_pos[2]);

  particles[i].x = dist_x(gen);
  particles[i].y = dist_y(gen);
  particles[i].theta = dist_theta(gen);


  }

}
void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
  // TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
  //   observed measurement to this particular landmark.
  // NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
  //   implement this method and use it as a helper during the updateWeights phase.
    
    
  LandmarkObs landmark;
  vector<LandmarkObs> closest_landmark_list;
  // 
  // find closest distance for each observations to predicted, and save that as the landmark coord
  // if dist is too far, skip
  double distance;
  double dist_threshold;
  // int landmark_id = -1;
  for (int i = 0; i < observations.size(); i++)
  {
    dist_threshold = 99999;
    for (int j = 0; j < predicted.size(); j++)
    {

      distance = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);
      if (distance < dist_threshold)
      {
        // if new distnace is below the threshold, this predicted is the new lankmark, and set threshold equal to new distance, so it will only update if it find a better/closer predicted point. 
        dist_threshold = distance;

        landmark.id = predicted[j].id;
        landmark.x = predicted[j].x;
        landmark.y = predicted[j].y;
      }
    }
    closest_landmark_list.push_back(landmark);
    // observations[i].id = landmark_id;
  }
  //update back to observations
  observations = closest_landmark_list;
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

  for(int i = 0; i < num_particles; i++)
  {
  	
    vector<LandmarkObs> tran_obs_list;
    LandmarkObs Ob;
    for (int j =0; j < observations.size(); j++)
    {
      LandmarkObs transformed_obs;
      Ob = observations[j];
      // transformed_obs.id = Ob.id;
      transformed_obs.x = particles[i].x + cos(particles[i].theta) * Ob.x - sin(particles[i].theta) * Ob.y;
      transformed_obs.y = particles[i].y + sin(particles[i].theta) * Ob.x + cos(particles[i].theta) * Ob.y;
      tran_obs_list.push_back(transformed_obs);
    }


    

    // set inititial weights
    particles[i].weight = 1.0;
    weights[i] = 1.0;

    // from the particle, filter out landmarks which are beyond detection range
    vector<LandmarkObs> filtered_landmark_list;

    // remove lankmark that beyond the sensor range
    for (int m = 0; m < map_landmarks.landmark_list.size(); m++)
    {
      Map::single_landmark_s Lm = map_landmarks.landmark_list[m];

      double distance = dist(particles[i].x, particles[i].y, map_landmarks.landmark_list[m].x_f, map_landmarks.landmark_list[m].y_f);
      if (distance <= sensor_range)
      {
        LandmarkObs in_range_landmark;
        in_range_landmark.id = Lm.id_i;
        in_range_landmark.x = Lm.x_f;
        in_range_landmark.y = Lm.y_f;
        filtered_landmark_list.push_back(in_range_landmark);
      }
    }

    // get the nearest point on map
    vector<LandmarkObs> closest_obs_list = tran_obs_list;


    dataAssociation(filtered_landmark_list, closest_obs_list);

    vector<int> P_associations;
    vector<double> P_sense_x;
    vector<double> P_sense_y;

    const double std_x = std_landmark[0];
    const double std_y = std_landmark[1];
    const double denominator = 2. * M_PI * std_x * std_y;
    const double std_x2 = 2. * std_x * std_x;
    const double std_y2 = 2. * std_y * std_y;
    for(int j = 0; j < closest_obs_list.size(); j++)
    {

      // calc multi-variant gaussian distribution prob
      double x = tran_obs_list[j].x;
      double y = tran_obs_list[j].y;
      double mu_x = closest_obs_list[j].x;
      double mu_y = closest_obs_list[j].y;
      double x_diff = (x - mu_x) * (x - mu_x) / std_x2;
      double y_diff = (y - mu_y) * (y - mu_y) / std_y2;

      double probability = exp(-(x_diff + y_diff))/ denominator;

      if(probability > 0)
      {
        particles[i].weight = particles[i].weight * probability;
      }

      P_associations.push_back(closest_obs_list[j].id);
      P_sense_x.push_back(closest_obs_list[j].x);
      P_sense_y.push_back(closest_obs_list[j].y);
      
    }

    weights[i] = particles[i].weight;
    particles[i] = SetAssociations(particles[i], P_associations, P_sense_x, P_sense_y);
  }
  // print_vector(weights);
  // for (int i = 0; i < num_particles; i++)
  // {
  //   cout << "id: " << particles[i].id << " , " << "x: " << particles[i].x << " , " << "y: " << particles[i].y << " , " << "theta: " << particles[i].theta << " weight: " << particles[i].weight <<endl;
  // }

}


void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	// default_random_engine gen;
	vector<Particle> resample_particles;
	random_device seed;
	mt19937 gen(seed());
	discrete_distribution<int> index(weights.begin(), weights.end());
	resample_particles.clear();
	for (int n = 0; n < num_particles; n++)
	{
		resample_particles.push_back(particles[index(gen)]);
	}
	particles = resample_particles;

	cout << "finish resample" << endl;
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
