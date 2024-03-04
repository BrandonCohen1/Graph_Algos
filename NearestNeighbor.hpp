#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <vector>
#include <iomanip>
#include <chrono>

struct Node {
    int ID;
    double x;
    double y;
};

std::vector<Node*> create_Nodes(std::string pathname){
    // Create ifstream
    std::ifstream inputFile(pathname);

    // Create lines and flag t start reading
    std::string line;
    bool found = false;

    // Create vector of nodes
    std::vector<Node*> nodes;

    // Reading in the file and creating Nodes
    while (std::getline(inputFile, line)) {
        if (line.find("NODE_COORD_SECTION") != std::string::npos) {
            found = true;
            break;
        }
    }

    if (found) {
        int i=0;
        while (std::getline(inputFile, line)) {

            std::istringstream iss(line);

            Node* newNode = new Node;

            // Check if the extraction is successful before using the values
            if (!(iss >> newNode->ID >> newNode->x >> newNode->y)) {
                //std::cerr << "Error reading line: " << line << std::endl;
                break;
            }
            nodes.push_back(newNode);
        }
        
    } else {
        std::cout << "NODE_COORD_SECTION not found in the file." << std::endl;
    }

    // Close file
    inputFile.close();

    return nodes;
}

// Function to calculate the distance between two points
double distance(struct Node* p1, struct Node* p2) {
    float diff_x = p2->x - p1->x;
    float diff_y = p2->y - p1->y;

    float diff_x_sqrd = pow(diff_x, 2);
    float diff_y_sqrd = pow(diff_y, 2);

    return sqrt( (diff_x_sqrd + diff_y_sqrd) );
}

// Function to create the adjacency matrix
std::vector<std::vector<double>> createAdjacencyMatrix(const std::vector<Node*>& nodes) {
    int size = nodes.size();
    // Create adjancency matrix by create a vector of double and then n vectors with 0 value
    std::vector<std::vector<double>> adjacencyMatrix(size, std::vector<double>(size, 0.0));

    // Iterate over the whole adjancency matrix
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (i != j) {
                // If not a diagnoal index, calculate distance between nodes
                adjacencyMatrix[i][j] = distance(nodes[i], nodes[j]);
            } else {
                // Diagonal entries (distance to itself) are set to infinity
                adjacencyMatrix[i][j] = std::numeric_limits<double>::infinity();
            }
        }
    }
    return adjacencyMatrix;
}

// Function to find the nearest neighbor tour
std::vector<int> nearestNeighborTour(const std::vector<std::vector<double>>& adjacencyMatrix) {
    int size = adjacencyMatrix.size();      // Size of adjacency matrix
    std::vector<int> tour;                  // Tour of ID's of the nodes
    tour.push_back(0);                      // The first node will be ID = 1

    std::vector<bool> visited(size, false); // Create a vector of bools that will indicate if a node was visited or not
    visited[0] = true;                      // The first node is visited

    // Iterate over the whole adjacency matrix
    for (int i = 1; i < size; i++) {
        int current = tour.back();          // Keep track of the last node that is added to the tour
        double minDistance = std::numeric_limits<double>::infinity();   // Set minimum distance to largest number
        int nearestNeighbor = 0;            // The first nearest neighbor is ID 1

        for (int j = 0; j < size; j++) {
            if (!visited[j] && adjacencyMatrix[current][j] < minDistance) {
                minDistance = adjacencyMatrix[current][j];
                nearestNeighbor = j;
            }
        }

        tour.push_back(nearestNeighbor);
        visited[nearestNeighbor] = true;
    }

    return tour;
}

// Function to print the tour and calculate the total distance
void printTour(const std::vector<int>& tour, const std::vector<Node*>& nodes) {
    double totalDistance = 0.0;
    int last;
    for (int i = 0; i < tour.size(); ++i) {
        int nodeID = nodes[tour[i]]->ID;
        std::cout << nodeID << " ";
        if (i < tour.size() - 1) {
            double distanceToNext = distance(nodes[tour[i]], nodes[tour[i + 1]]);
            totalDistance += distanceToNext;
            last = i+1;
        }
    }
    std::cout << nodes[tour[0]]->ID << std::endl;
    totalDistance += distance(nodes[tour[last]], nodes[tour[0]]);

    std::cout << "Total Distance: " << totalDistance << std::endl;
}

// Nearest Neighbor algorithm
void nearestNeighbor(std::string pathname) {

    // Read the file and create vector of Nodes
    std::vector<Node*> nodes = create_Nodes(pathname);

    // Start Time of the Algorithm
    auto start = std::chrono::high_resolution_clock::now();

    // Create the adjacency matrix
    std::vector<std::vector<double>> adjacencyMatrix = createAdjacencyMatrix(nodes);

    // Find the nearest neighbor tour
    std::vector<int> tour = nearestNeighborTour(adjacencyMatrix);

    // End time of algorithm
    auto end = std::chrono::high_resolution_clock::now();

    auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    int duration = diff.count();

    // Output the tour and total distance
    printTour(tour, nodes);

    // Time
    std::cout << "Time in ms: " << duration;


    // Delete the Node objects
    for (int i = 0; i < nodes.size(); i++) {
        delete nodes[i];
    }
}

