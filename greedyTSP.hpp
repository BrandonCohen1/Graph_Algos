#ifndef GREEDYTSP_HPP
#define GREEDYTSP_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <chrono>
#include <limits>
#include <queue>


struct Node {
    int ID;
    double x;
    double y;
    bool is_in_matrix = false;
    int matrix_index;

    // 0 = gray edge (no edges)
    // 1 = white edge (1 edge)
    // 2 = black edge (2 edges)
    int num_edges = 0;
};

struct Edge {
    Node* from;
    Node* to;
    double weight;

    // bool operator<(const Edge& other) const {
    bool operator<(const Edge& other) const {       // Used for minheap comparisons
        return weight > other.weight;
    }
};


// Read in the file and create the nodes
std::vector<Node*> createNodes(std::ifstream& inputFile) {
    std::string line;
    bool found = false;

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

        return nodes;
    } else {
        std::cout << "NODE_COORD_SECTION not found in the file." << std::endl;
        return {};
    }
}


// Calculate Euclidean distance
double distance(const Node* p1, const Node* p2) {

    double diff_x = p2->x - p1->x;
    double diff_y = p2->y - p1->y;

    double diff_x_sqrd = pow(diff_x, 2);
    double diff_y_sqrd = pow(diff_y, 2);

    return sqrt(diff_x_sqrd + diff_y_sqrd);
}


void greedyTSP(const std::string& filename) {
    std::ifstream inputFile(filename);

    if (!inputFile.is_open()) {
        std::cout << "Unable to open the file." << std::endl;
        return;
    }

    // Vector holding the Node*
    std::vector<Node*> nodes = createNodes(inputFile);
    // Close the file as it is not needed no more
    inputFile.close();


    // Start of greedy algorithm
    auto start = std::chrono::high_resolution_clock::now();

    // Min heap to store edges
    std::priority_queue<Edge> edges;


    // Double for loop to iterate to get all pairs of nodes
    for (int i = 0; i < nodes.size() - 1; i++) {
        for (int j = i + 1; j < nodes.size(); j++) {
            Edge e;
            e.from = nodes[i];
            e.to = nodes[j];
            e.weight = distance(nodes[i], nodes[j]);
            edges.push(e);
        }
    }

    // Will hold vectors of shortest distances
    std::vector<std::vector<Node*>> nodeMatrix;

    // Distance
    double totalDistance = 0.0;
    
    // Add to the nodeMatrix
    while (!edges.empty()) {
        // Make edge and pop
        Edge e = edges.top();
        edges.pop();

        // Create Node ptrs and index
        Node* from = e.from;
        Node* to = e.to;
        int index_from = from->ID - 1;
        int index_to = to->ID - 1;

        // Add to node matrix if the first node is not in the matrix
        if (!from->is_in_matrix) {
            if (!to->is_in_matrix) {
                // Both nodes are not in the matrix
                // Create a new vector and add it to the matrix
                std::vector<Node*> temp_vec = { from, to };
                nodeMatrix.push_back(temp_vec);

                // Add distance
                totalDistance += e.weight;

                // Update node information
                from->is_in_matrix = true;
                from->matrix_index = nodeMatrix.size() - 1;
                to->is_in_matrix = true;
                to->matrix_index = nodeMatrix.size() - 1;

                // Update number of edges
                from->num_edges++;
                to->num_edges++;
            } else {
                // Only e.from is not in the matrix
                // Find which vector has e.to and add e.from at the front or back
                int matrix_index = to->matrix_index;
                if (to == nodeMatrix[matrix_index].front()) {
                    nodeMatrix[matrix_index].insert(nodeMatrix[matrix_index].begin(), from);
                    totalDistance += e.weight;
                } else {
                    nodeMatrix[matrix_index].push_back(from);
                    totalDistance += e.weight;
                }

                // Update node information
                from->is_in_matrix = true;
                from->matrix_index = matrix_index;

                // Update number of edges
                from->num_edges++;
                to->num_edges++;
            }
        // Add to node matrix if second node is in the matrix
        } else if (!to->is_in_matrix) {
            // Only e.to is not in the matrix
            // Find which vector has e.from and add e.to at the front or back
            int matrix_index = from->matrix_index;
            if (from == nodeMatrix[matrix_index].front()) {
                nodeMatrix[matrix_index].insert(nodeMatrix[matrix_index].begin(), to);
                totalDistance += e.weight;
            } else {
                nodeMatrix[matrix_index].push_back(to);
                totalDistance += e.weight;
            }

            // Update node information
            to->is_in_matrix = true;
            to->matrix_index = matrix_index;

            // Update number of edges
            from->num_edges++;
            to->num_edges++;
        }
        // Both e.from and e.to are already in the matrix so don't do anything
    }

    // Vector of sorted edges of the ends of vectors of the node matrix
    std::vector<Edge> sortededges;

    // Iterate over each vector in nodeMatrix
    for (int i = 0; i < nodeMatrix.size() - 1; i++) {
        // Get the last node in the vector in the node matrix
        Node* firstNode = nodeMatrix[i].front();
        Node* lastNode = nodeMatrix[i].back();

        // Iterate over the other vectors to get the first node
        for (int j = i + 1; j < nodeMatrix.size(); j++) {
            // Get the first node in the next vector
            Node* front = nodeMatrix[j].front();
            Node* last = nodeMatrix[j].back();

            // Distance between the last node and the first node of the next vector
            double weight_first_last = distance(firstNode, last);
            double weight_last_front = distance(lastNode, front);

            // Create an edge connecting the last node to the first node of the next vector
            Edge edge_front_last;
            edge_front_last.from = firstNode;
            edge_front_last.to = last;
            edge_front_last.weight = weight_first_last;


            Edge edge_last_front;
            edge_last_front.from = lastNode;
            edge_last_front.to = front;
            edge_last_front.weight = weight_last_front;

            // Add the edge 
            sortededges.push_back(edge_front_last);
            sortededges.push_back(edge_last_front);
        }
    }

    // Sort the edges
    std::sort(sortededges.begin(), sortededges.end());

    // Iterate over the sorted edges
    for (const Edge& e : sortededges) {
        // Find the indices of vectors containing e.from and e.to
        auto it_from = std::find_if(nodeMatrix.begin(), nodeMatrix.end(),
                                    [e](const std::vector<Node*>& nodes) {
                                        return nodes.front() == e.from || nodes.back() == e.from;
                                    });

        auto it_to = std::find_if(nodeMatrix.begin(), nodeMatrix.end(),
                                [e](const std::vector<Node*>& nodes) {
                                    return nodes.front() == e.to || nodes.back() == e.to;
                                });

        if (it_from != nodeMatrix.end() && it_to != nodeMatrix.end() && it_from != it_to) {
            // Merge the vectors and erase the second vector
            if ((*it_from).front() == e.from) {
                // (from, to)
                (*it_from).insert((*it_from).begin(), (*it_to).begin(), (*it_to).end());
                totalDistance += e.weight;
            } else {
                // (to, from)
                (*it_from).insert((*it_from).begin(), (*it_to).rbegin(), (*it_to).rend());
                totalDistance += e.weight;
            }
            // Erase vector because it is not needed no more
            nodeMatrix.erase(it_to);
        }
    }

    // Print tours
    for(int i=0; i<nodeMatrix[0].size()-1; i++){
        std::cout << nodeMatrix[0][i]->ID << " ";
    }
    std::cout << nodeMatrix[0][0]->ID << std::endl;

    // Add last distance between last node and the first nodeS
    totalDistance += distance(nodeMatrix[0][0], nodeMatrix[0][nodeMatrix.size()-1]);

    auto end = std::chrono::high_resolution_clock::now();
    auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    int duration = diff.count();

    // Distance
    std::cout << "Total Distance: " << totalDistance << std::endl;

    // Time
    std::cout << "Time in ms: " << duration << std::endl;


    for (auto node : nodes) {
        delete node;
    }
}

#endif 