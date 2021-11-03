/*******************************************************************************
//  4TSPga.cpp                      Author: Ian Nobile
//  Section: 50                     Due Date: 20 October 2021
//
//  This program solves a Travelling Salesperson Problem dataset with up to 100 
//  cities using a genetic algorithm that generates and mates randomised tours 
//  and adds mutations to 10% of their offspring in order to maintain a healthy 
//  gene pool. Two crossover methods and two mutative methods are available for
//  comparison.
//
*******************************************************************************/

#include <iostream> // print to console
#include <chrono>   // time the speed of the program
#include <fstream>  // read files
#include <string>   // header buffer for seeking inside a file
#include <cfloat>   // allows the use of FLT_MAX
#include <vector>   // easy arraying
#include <cmath>    // distance formula
#include <algorithm>// easy array shuffling
#include <random>   // used for easy array shuffling
#include <C:\Users\admin\Desktop\coastal_islander\me\uofl\ai\project4\4TSPga\matplotlib-cpp-master\matplotlibcpp.h>

using namespace std;
using namespace std::chrono;
namespace plt = matplotlibcpp;

// class declarations:
class Node {
public:
    int num;
    float x;
    float y;
};

class Graph {
public:
    vector<Node> nodes;
};

class Tour {
public:
    array<int,100> nArray;
    float dist;
    array<bool, 101> legitimate;
};


// function prototypes:
Graph buildGraph(char*);
float distCheck(Node, Node);
float calcTourDist(const array<array<float, 101>, 101> &, Tour);
Tour cx2(const array<array<float, 101>, 101> &, Tour, Tour, bool);
Tour scx(const array<array<float, 101>, 101> &, Tour, Tour, bool);


//------------------------------------------------------------------------------
//  Main Function
//------------------------------------------------------------------------------
int main(int argc, char* argv[]) {
    // check if path was passed as arg:
    if (argc == 1) {
        cout << "Please pass the path to the .TSP file as a command line argument" << endl;
        return 1;
    }
    
    auto start = high_resolution_clock::now();  // start timer

    // begin with a friendly greeting
    cout << "Hello and welcome to the (Genetic Algorithm) Travelling Salesperson Problem Solver" << endl;

    // fill graph
    Graph graph = buildGraph(argv[1]);
    
    // graphing:
    //vector<double> xPlot;
    //vector<double> yPlot;
    //plt::xkcd();
    //for (Node node : graph.nodes) {
    //    xPlot.push_back(node.x);
    //    yPlot.push_back(node.y);
    //    plt::plot(xPlot, yPlot,"-o");
    //    xPlot.pop_back();
    //    yPlot.pop_back();
    //}
    //plt::show();
    

    // calculate all possible distances up front (saves cycles))
    array<array<float, 101>, 101> distRef;
    for (int i = 0; i < 101; i++) {
        for (int j = 0; j < 101; j++) {
            if (i == 0 || j == 0) {
                distRef[i][j] = 0.0;    // improves readability of checks to distRef
            } else {
                distRef[i][j] = distCheck(graph.nodes[i - 1], graph.nodes[j - 1]);
            }
        }
    }

    // create an initial population of 100 totally randomised tours
    array<Tour, 100> pop;
    for (int i = 0;i < pop.size();i++) {
        Tour tour = Tour();
        fill(tour.legitimate.begin(), tour.legitimate.end(), true); // bool array initaialisation
        tour.legitimate[0] = false; // no destination 0

        // shuffle nodes of graph object
        unsigned seed = system_clock::now().time_since_epoch().count();
        shuffle(graph.nodes.begin(), graph.nodes.end(), std::default_random_engine(seed));
        
        for (int j = 0;j < 100;j++) {
            tour.nArray[j] = graph.nodes[j].num;
        }

        //set tour distance and store in pop
        tour.dist = calcTourDist(distRef, tour);
        pop[i] = tour;
    }
    sort(graph.nodes.begin(), graph.nodes.end(), [](Node a, Node b) { return a.num < b.num; }); //re-sort graph

    // prep mutation probability
    srand(system_clock::now().time_since_epoch().count());
    default_random_engine generator;
    std::bernoulli_distribution distribution(0.1); //10% of offspring
    bool mutate = false;

    //Main Loop:
    for (int upperBound = 1;upperBound < 1001;upperBound++) {
        // sort by dist
        sort(pop.begin(), pop.end(), [](Tour a, Tour b) { return a.dist < b.dist; });
        
        // calc average
        float avg = 0;
        for (int i = 0;i < 100;i++) {
            avg += pop[i].dist;
        }
        avg /= 100;
        
        // calc standardDev
        float var = 0;
        float sig = 0;
        for (int i = 0;i < 100;i++) {
            var += pow((pop[i].dist - avg),2);
        }
        var /= 100;
        sig = sqrt(var);

        // report stats:
        cout << "For generation " << upperBound << ":\tMin: " << pop[0].dist << "\tMax " << pop[99].dist << "\tAvg: " << avg << "\tSig: " << sig << endl;
        //xPlot.push_back(upperBound);  // improvement graphing
        //yPlot.push_back(pop[0].dist);

        // randomly mate top 50 fittest and store offspring in newPop
        array<Tour, 100> newPop;
        for (int i = 0; i < 100; i++){
            srand(system_clock::now().time_since_epoch().count());
            int p1Index = rand() % 50, p2Index = rand() % 50;
            mutate = false;
            if (distribution(generator)) { mutate = true; }
            //newPop[i] = cx2(distRef, pop[p1Index], pop[p2Index], mutate); //uncomment for cx2()
            //newPop[i] = scx(distRef, pop[p1Index], pop[p2Index], mutate);   //uncomment for scx()
        }
        pop = newPop;
    }   //end Main Loop


    // one final sorting
    sort(pop.begin(), pop.end(), [](Tour a, Tour b) { return a.dist < b.dist; });
    
    // graph winning tour and report statistics:
    cout << endl << "The shortest route visiting all the points is:" << endl << endl << "\t";
    for (int i = 0; i < 100;i++) {
        cout << pop[0].nArray[i];
        //xPlot.push_back(graph.nodes[pop[0].nArray[i] - 1].x);
        //yPlot.push_back(graph.nodes[pop[0].nArray[i] - 1].y);
        //plt::plot(xPlot, yPlot);
        //plt::pause(1);
        cout << " -> ";
    }
    cout << pop[0].nArray[0];
    //xPlot.push_back(graph.nodes[pop[0].nArray[99] - 1].x);
    //yPlot.push_back(graph.nodes[pop[0].nArray[99] - 1].y);
    //plt::plot(xPlot, yPlot);
    //plt::pause(1);
    
    //xPlot.push_back(graph.nodes[pop[0].nArray[0] - 1].x);
    //yPlot.push_back(graph.nodes[pop[0].nArray[0] - 1].y);
    //plt::plot(xPlot, yPlot);
    //plt::pause(1);
    //plt::show();
    cout << endl << "and back again at a total distance of " << pop[0].dist << endl << endl;
    
    auto stop = high_resolution_clock::now();   // stop timer
    auto duration = duration_cast<microseconds>(stop - start);  // calculate elapsed time
    cout << "Program execution took " << duration.count() / 1000000.0 << "s" << endl << endl;

    system("pause");
    return 0;
}



//------------------------------------------------------------------------------
//  Function Definitions:
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//  Read .TSP file, creates nodes from coordinates and combines all in an
//  undirected graph object
//------------------------------------------------------------------------------
Graph buildGraph(char* argv) {
    // open .TSP file in read-only mode:
    ifstream tspfile;
    tspfile.open(argv, ios::in);
    // ensure file exists:
    if (!tspfile.is_open()) {
        Graph graph;
        return graph;
    }
    Graph graph;
    int dimension = 0;
    string heading = "";
    // advance buffer to the dimension section:
    while (heading.compare("DIMENSION:") != 0) {
        tspfile >> heading;
    }
    tspfile >> dimension;
    // advance buffer to the coordinates section:
    while (heading.compare("NODE_COORD_SECTION") != 0) {
        tspfile >> heading;
    }
    // create nodes and push to graph vector
    Node newNode = Node();
    for (int i = 0;i < dimension;i++) {
        tspfile >> newNode.num;
        tspfile >> newNode.x;
        tspfile >> newNode.y;
        graph.nodes.push_back(newNode);
    }
    // The graph is now created, and we are finished with the .TSP file
    tspfile.close();
    return graph;
}

//------------------------------------------------------------------------------
//  Returns the distance between two graph nodes using the pythagoran formula: 
//  dist = sqrt((x2 - x1)^2 + (y2 - y1)^2)
//------------------------------------------------------------------------------
float distCheck(Node first, Node second) {
    float dx = second.x - first.x;
    float dy = second.y - first.y;
    return sqrt(dx * dx + dy * dy);
}

//------------------------------------------------------------------------------
//  Calculates the length of a tour
//------------------------------------------------------------------------------
float calcTourDist(const array<array<float, 101>, 101> &distRef, Tour tour){
    float minDist = 0.0;
    for (int i = 0;i < 99;i++) {
        minDist += distRef[tour.nArray[i]][tour.nArray[i + 1]];
    }
    minDist += distRef[tour.nArray[99]][tour.nArray[0]];
    return minDist;
}

//------------------------------------------------------------------------------
//  An implementation of SCX, by Zakir H. Ahmed
//------------------------------------------------------------------------------
Tour scx(const array<array<float, 101>, 101> &distRef, Tour p1, Tour p2, bool mutate){
    Tour o = Tour();
    fill(o.legitimate.begin(), o.legitimate.end(), true);
    o.legitimate[0] = false;    // no destination 0 in tsp data
    o.nArray[0] = p1.nArray[0]; // start with first dest of first parent
    o.legitimate[o.nArray[0]] = false;
    int index1 = 0, index2 = 0;
    float dist1 = 0.0, dist2 = 0.0;

    for (int i = 0; i < 99; i++){
        // position runners
        index1 = 0, index2 = 0;
        while (p1.nArray[index1] != o.nArray[i]) { index1++; }
        index1++;
        index1 %= 100;
        while (!o.legitimate[p1.nArray[index1]]) { index1++; index1 %= 100; }
        while (p2.nArray[index2] != o.nArray[i]) { index2++; }
        index2++;
        index2 %= 100;
        while (!o.legitimate[p2.nArray[index2]]) { index2++; index2 %= 100;}
        
        // calculate the distances
        dist1 = distRef[o.nArray[i]][p1.nArray[index1]];
        dist2 = distRef[o.nArray[i]][p2.nArray[index2]];
        
        // append the closer node
        o.nArray[i + 1] = dist1 < dist2 ? p1.nArray[index1] : p2.nArray[index2];
        o.legitimate[o.nArray[i + 1]] = false;
    }
    
    // relegitimise o
    fill(o.legitimate.begin(), o.legitimate.end(), true);
    o.legitimate[0] = false;
    
    // mutation
    if (mutate) {
        index1 = rand() % 100;
        index2 = rand() % 100;
        
        // mutation 1 (pair swap)
        //swap(o.nArray[index1], o.nArray[index2]); // uncomment for single-swap mutation

        // mutation 2 (rsm)
        //if (index1 > index2) { swap(index1, index2); }  // uncomment for reverse sequence mutation
        //while (index1 < index2) {
        //    swap(o.nArray[index1++], o.nArray[index2--]);
        //}
    };

    o.dist = calcTourDist(distRef,o);
    return o;
}

//------------------------------------------------------------------------------
//  An implementation of CX2, by Abid Hussain et al.
//------------------------------------------------------------------------------
Tour cx2(const array<array<float, 101>, 101> &distRef, Tour p1, Tour p2, bool mutate){
    Tour o1 = Tour(), o2 = Tour();
    fill(o1.legitimate.begin(), o1.legitimate.end(), true);
    fill(o2.legitimate.begin(), o2.legitimate.end(), true);
    o1.legitimate[0] = false;    // no destination 0 in tsp data
    o2.legitimate[0] = false;
    
    int i = 0, index1 = 0, index2 = 0;

    while (i < 100) {
        if (i == 0) {
            // Step 2. Select 1st bit from second parent as a 1st bit of first offspring.
            o1.nArray[0] = p2.nArray[0];
        } else {
            // Step 4. Find value of selected bit from Step 3 in P1.
            index1 = 0;
            while (p1.nArray[index1] != o2.nArray[i - 1]) {
                index1++;
            }
            if (p2.nArray[index1] == o1.nArray[0]) { break; } // stop cycles happening
            o1.nArray[i] = p2.nArray[index1];
        }
        // Step 3. Find value of selected bit twice in p1 before appending from p2
        index1 = 0;
        while (p1.nArray[index1] != o1.nArray[i]) {
            index1++;
        }
        index2 = 0;
        while (p1.nArray[index2] != p2.nArray[index1]) {
            index2++;
        }
        if (p2.nArray[index2] == o2.nArray[0]) { break; } // stop cycles happening
        o2.nArray[i++] = p2.nArray[index2];
    }

    // Step 6: Handle cycled crossovers:
    if (i != 100) {
        array<int, 100> remainder1 = p2.nArray;
        array<int, 100> remainder2 = p2.nArray;

        // bad n^2 loop removes o1’s bits from second remainder
        for (int j = 0; j < 100; j++) {
            for (int k = 0; k < 100; k++) {
                if (remainder2[j] == o1.nArray[k]) {
                    remainder2[j] = 0;
                    break;
                }
            }
        }
        // bad n^2 loop removes o2’s bits from first remainder
        for (int j = 0; j < 100; j++) {
            for (int k = 0; k < 100; k++) {
                if (remainder1[j] == o2.nArray[k]) {
                    remainder1[j] = 0;
                    break;
                }
            }
        }

        // append remainders to offsprings
        index1 = 0;
        for (int iter = i; iter < 100; iter++) {
            while (index1 != 100 && remainder1[index1] == 0) {
                index1++;
            }
            if (index1 >= 100) { break; }
            o2.nArray[iter] = remainder1[index1++];
        }
        index2 = 0;
        for (int iter = i; iter < 100; iter++) {
            while (index2 != 100 && remainder2[index2] == 0) {
                index2++;
            }
            if (index2 >= 100) { break; }
            o1.nArray[iter] = remainder2[index2++];
        }
    }   // end step6


    // mutation (10%)
    if (mutate) {
        index1 = rand() % 100;
        index2 = rand() % 100;

        // mutation 1 (pair swap)
        //swap(o1.nArray[index1], o1.nArray[index2]);   //uncomment for single-swap mutation
        //swap(o2.nArray[index1], o2.nArray[index2]);

        // mutation 2 (rsm)
        //if (index1 > index2) { swap(index1, index2); }    //uncomment for reverse sequence mutation
        //while (index1 < index2) {
        //    swap(o1.nArray[index1++], o1.nArray[index2--]);
        //    swap(o2.nArray[index1++], o2.nArray[index2--]);
        //}
    };

    o1.dist = calcTourDist(distRef, o1);
    o2.dist = calcTourDist(distRef, o2);

    // return fitter of two siblings
    return o1.dist < o2.dist ? o1 : o2;
}

