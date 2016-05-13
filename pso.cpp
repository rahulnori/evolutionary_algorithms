#include <cmath>

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;

#include <iomanip>
using std::fixed;
using std::setw;
using std::setprecision;

#include <limits>
using std::numeric_limits;

#include <vector>
using std::vector;



double objective_function(vector<double> parameters) {
    double sum_of_squares = 0.0;

    for (uint32_t i = 0; i < parameters.size(); i++) {
        sum_of_squares += parameters[i] * parameters[i];
    }

    //return the negative value because we're doing maximization here
    return -sum_of_squares;
}


int main(int argc, char **argv) {
    if (argc != 4) {
        cerr << "ERROR, incorrect arguments!" << endl;
        cerr << "usage: " << endl;
        cerr << "\t" << argv[0] << " <maximum iterations> <number particles> <number parameters>" << endl;
        exit(1);
    }

    //parse the command line parameters
    uint32_t maximum_iterations = atoi(argv[1]);
    uint32_t number_particles = atoi(argv[2]);
    uint32_t number_parameters = atoi(argv[3]);


    //seed the random number generator so it initializes with different particles and velocities 
    //when we run the program
    srand48(time(NULL));

    //these are parameters to the particle swarm optimization,
    //changing these can make it converge faster or slower
    double velocity_scale = 1.0;
    double inertia = 0.95;
    double global_constant = 1.5;
    double local_constant = 1.5;

    //define the minimum and maximum bounds for the search,
    //you'll need to come up with your own for the chemistry
    //code
    vector<double> minimums(number_parameters, -15.0);
    vector<double> maximums(number_parameters,  15.0);

    //vectors to hold the particles, their velocities and their
    //current fitness
    vector< vector<double> > particles;
    vector< vector<double> > velocity;
    vector< double> particle_fitness;

    //vectors to hold the best found positions for each particle
    //and the fitness of that best found position
    vector< vector<double> > local_best;
    vector< double > local_best_fitness;

    //the best parameters found for the entire search along
    //with its fitness.  We can initialize this to the lowest
    //possible double so that anything found after this will
    //update the global best fitness.
    vector<double> global_best(number_parameters, 0.0);
    double global_best_fitness = -numeric_limits<double>::max();

    //initialize all the particles
    //each particle's local best is the initial particle
    for (uint32_t i = 0; i < number_particles; i++) {
        //initialize each particle randomly within the bounds
        vector<double> current_particle(number_parameters, 0.0);
        vector<double> current_velocity(number_parameters, 0.0);

        //cout << "current particle: ";
        for (uint32_t j = 0; j < number_parameters; j++) {
            //set each parameter for each particle to a random number between the maximum
            //and minimum bounds
            current_particle[j] = ((maximums[j] - minimums[j]) * drand48() + minimums[j]);
            //cout << " " << setw(10) << setprecision(5) << current_particle[j];

            //initialize the velocity for each parameter to another number that could be
            //between the negative difference between the max and mim and the difference
            //between the max and min
            double velocity_max = maximums[j] - minimums[j];
            double velocity_min = -velocity_max;

            //scaling the initial velocities down a bit with the velocity scale
            //parameter can improve the convergence rate
            current_velocity[j] = velocity_scale * ( ((velocity_max - velocity_min) * drand48()) + velocity_min );
        }
        //cout << endl;

        //push these vectors back into the vectors holding
        //all the particles, their fitnesses.
        particles.push_back(current_particle);
        velocity.push_back(current_velocity);

        //since this is the first position of the particle
        //it is also the local best
        local_best.push_back(current_particle);

        /*
        cout << "velocity:         ";
        for (uint32_t j = 0; j < number_parameters; j++) {
            cout << " " << setw(10) << setprecision(5) << velocity[i][j];
        }
        cout << endl;
        */

        //calculate the fitness of the particle and set it to
        //the particle's current fitness and the local best
        //fitness
        double fitness = objective_function(current_particle);
        particle_fitness.push_back(fitness);
        local_best_fitness.push_back(fitness);

        if (fitness > global_best_fitness) {
            //update the global best found position if we've found
            //a new global best while initializing the particles

            global_best_fitness = particle_fitness[i];
            cout << "new global best: " << setw(10) << setprecision(5) << fixed << global_best_fitness << " -- ";
            for (uint32_t j = 0; j < number_parameters; j++) {
                global_best[j] = particles[i][j];
                cout << " " << setw(10) << setprecision(5) << fixed << global_best[j];
            }
            cout << endl;
        }
    }

    cout << "Iterating!" << endl;

    //this does the actual particle swarm optimization after the initialization
    //see the algorithm here: https://en.wikipedia.org/wiki/Particle_swarm_optimization
    for (uint32_t iteration = 0; iteration < maximum_iterations; iteration++) {
        for (uint32_t i = 0; i < number_particles; i++) {
            double global_random = drand48();
            double local_random = drand48();

            for (uint32_t j = 0; j < number_parameters; j++) {
                //using the two random numbers, calculate a new velocity and then
                //use that velocity to update the particle's position.
                velocity[i][j] += global_constant * global_random * (global_best[j] - particles[i][j]) + local_constant * local_random * (local_best[i][j] - particles[i][j]);
                //having an inertia < 1 can make things converge quicker
                particles[i][j] += inertia * velocity[i][j];

                //make sure the parameters of the particle stay within bounds
                if (particles[i][j] > maximums[j]) particles[i][j] = maximums[j];
                if (particles[i][j] < minimums[j]) particles[i][j] = minimums[j];
            }

            //calculate the new fitness for the particle
            particle_fitness[i] = objective_function(particles[i]);

            if (particle_fitness[i] > local_best_fitness[i]) {
                //we've found a new local best position for the particle,
                //update the local best fitness and parameters
                for (uint32_t j = 0; j < number_parameters; j++) {
                    local_best[i][j] = particles[i][j];
                }
                local_best_fitness[i] = particle_fitness[i];

                if (particle_fitness[i] > global_best_fitness) {
                    //we've found a new global best position for the particle
                    //swarm.  print it out and update the global best fitness
                    //and parameters
                    global_best_fitness = particle_fitness[i];
                    cout << "new global best: " << setw(10) << setprecision(5) << fixed << global_best_fitness << " -- ";
                    for (uint32_t j = 0; j < number_parameters; j++) {
                        global_best[j] = particles[i][j];
                        cout << " " << setw(10) << setprecision(5) << fixed << global_best[j];
                    }
                    cout << endl;
                }
            }
        }
    }
}
