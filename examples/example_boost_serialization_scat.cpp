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
    Particle() : m_x(0.0), m_y(0.0), m_z(0.0) {}
    Particle(double x, double y, double z) : m_x(x), m_y(y), m_z(z) {}
    Particle(const Particle& p) : m_x(p.m_x), m_y(p.m_y), m_z(p.m_z), m_time(p.m_time), m_positions(p.m_positions) {}
    ~Particle() = default;

    double                                  x() const { return m_x; }
    double                                  y() const { return m_y; }
    double                                  z() const { return m_z; }
    const std::vector<double>&              time() const { return m_time; }
    const std::vector<std::vector<double>>& positions() const { return m_positions; }

    void set_x(double x) { m_x = x; }
    void set_y(double y) { m_y = y; }
    void set_z(double z) { m_z = z; }

    void add_time(double t) { m_time.push_back(t); }
    void add_current_position() { m_positions.push_back({m_x, m_y, m_z}); }

 private:
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive& ar, const unsigned int version) {
        ar& m_x;
        ar& m_y;
        ar& m_z;
        ar& m_time;
        ar& m_positions;
    }

    double m_x;
    double m_y;
    double m_z;

    std::vector<double>              m_time;
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

    // Create a single particle.
    Particle p(1.0, 2.0, 5.0);

    int size = 0;

    char* buffer = (char*)malloc(size);
    if (process_rank == MASTER) {
        p.add_current_position();
        std::cout << "Process " << process_rank << " created particle: " << p.x() << " " << p.y() << " " << p.z() << std::endl;
        p.add_time(0.1);
        p.add_time(1.5);
        p.set_x(2.0);
        p.set_y(3.0);
        p.set_z(4.0);
        p.add_current_position();
        p.add_time(2.0);
        // Serialize the particle.
        std::stringstream             ss;
        boost::archive::text_oarchive oa(ss);
        oa << p;
        std::string str = ss.str();
        size            = str.size();
        buffer          = (char*)malloc(size);
        memcpy(buffer, str.c_str(), size);
        std::cout << "Process " << process_rank << " serialized particle: " << str << " of size: " << size << std::endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);



    // Scatter the serialized particle to all processes.

    // First, send the size of the serialized particle.
    MPI_Bcast(&size, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    std::cout << "Process " << process_rank << " received size: " << size << std::endl;
    if (process_rank != MASTER) {
        buffer = (char*)malloc(size);
    }
    // Then, send the serialized particle.
    MPI_Bcast(buffer, size, MPI_CHAR, MASTER, MPI_COMM_WORLD);
    std::cout << "Process " << process_rank << " received buffer: " << buffer << std::endl;

    MPI_Barrier(MPI_COMM_WORLD);

    // // Deserialize the particle.
    std::cout << "Process " << process_rank << " deserializing particle of size " << size << std::endl;

    MPI_Barrier(MPI_COMM_WORLD);

    std::stringstream             ss;
    ss.write(buffer, size);
    boost::archive::text_iarchive ia(ss);
    Particle p2;
    ia >> p2;
    std::cout << "Process " << process_rank << " deserialized particle: " << p2.x() << " " << p2.y() << " " << p2.z() << std::endl;

    MPI_Barrier(MPI_COMM_WORLD);

    // Finalize the MPI environment.
    MPI_Finalize();
    return 0;
}