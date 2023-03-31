/**
 * @file example_serialization1.cpp
 * @author your name (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2023-03-31
 *
 * @copyright Copyright (c) 2023
 *
 */

#include <mpi.h>     // must have a system with an MPI library
#include <stdio.h>   //printf
#include <stdlib.h>  //malloc

#include <algorithm>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#define MASTER 0  // One process will take care of initialization

class Particle {
 public:
    Particle() = default;
    Particle(double x, double y, double z) : m_x(x), m_y(y), m_z(z) {}
    ~Particle() = default;

    double x() const { return m_x; }
    double y() const { return m_y; }
    double z() const { return m_z; }
    const std::vector<double>& time() const { return m_time; }
    const std::vector<std::vector<double>>& positions() const { return m_positions; }


    void set_x(double x) { m_x = x; }
    void set_y(double y) { m_y = y; }
    void set_z(double z) { m_z = z; }

    void add_time(double t) { m_time.push_back(t); }
    void add_current_position() { m_positions.push_back({m_x, m_y, m_z}); }

 private:
 
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
            ar & m_x;
            ar & m_y;
            ar & m_z;
            ar & m_time;
            ar & m_positions;
    }

    double m_x;
    double m_y;
    double m_z;

    std::vector<double> m_time;
    std::vector<std::vector<double>> m_positions;
};

int main(int argc, char** argv) {
    // Initialize the MPI environment.
    MPI_Status status;

    std::cout << "Hello World!" << std::endl;

    int number_processes;
    int process_rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &number_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);

    std::cout << "Process " << process_rank << " of " << number_processes << " is alive." << std::endl;

    if (process_rank == MASTER) {

        std::cout << "Process " << process_rank << " is the master." << std::endl;
        // Create a single particle.
        Particle p(1.0, 2.0, 5.0);
        p.add_current_position();
        std::cout << "Process " << process_rank << " created particle: " << p.x() << " " << p.y() << " " << p.z() << std::endl;
        p.add_time(0.1);
        p.add_time(1.5);
        p.set_x(2.0);
        p.set_y(3.0);
        p.set_z(4.0);
        p.add_current_position();
        p.add_time(2.0);

        std::cout << "Process " << process_rank << " created particle: " << p.x() << " " << p.y() << " " << p.z() << std::endl;

        // Serialize the particle.
        std::stringstream ss;
        boost::archive::text_oarchive oa(ss);
        oa << p;

        // Send the serialized particle to the other processes.
        int size = ss.str().size();
        for (int i = 1; i < number_processes; i++) {
            MPI_Send(&size, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(ss.str().c_str(), size, MPI_CHAR, i, 0, MPI_COMM_WORLD);
        }
    } else {
        // Receive the serialized particle.
        int size;
        MPI_Recv(&size, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD, &status);
        char* buffer = (char*)malloc(size);
        MPI_Recv(buffer, size, MPI_CHAR, MASTER, 0, MPI_COMM_WORLD, &status);

        // Deserialize the particle.
        std::stringstream ss;
        ss.write(buffer, size);
        boost::archive::text_iarchive ia(ss);
        Particle p;
        ia >> p;

        // Print the particle.
        std::cout << "Process " << process_rank << " received particle: " << p.x() << " " << p.y() << " " << p.z() << std::endl;
        if (process_rank == 1) {
            std::cout << "Process " << process_rank << " received particle: " << p.time()[0] << " " << p.time()[1] << " " << p.time()[2] << std::endl;
            std::cout << "Process " << process_rank << " received particle: " << p.positions()[0][0] << " " << p.positions()[0][1] << " " << p.positions()[0][2] << std::endl;
        }
        
    }

    MPI_Finalize();
}