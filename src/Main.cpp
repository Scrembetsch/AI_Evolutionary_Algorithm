#include "Main.h"
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <sstream>
#include <cstring>
#include <map>
#include <list>
#include <queue>
#include <ctime>
#include <random>

using namespace std;

std::mt19937 rng(177013);

const int chessBoardSize = 3; // TODO: Read from Console
const int fieldSize = chessBoardSize * chessBoardSize; // TODO: Please Add after chessBoardSize
int initialPopulationCount = 10;

enum Mode 
{
    ASSOCIATIVE,
    COMPOSITE,
    MULTIMAGIC
};

typedef struct
{
    int arrangement[chessBoardSize * chessBoardSize];
    int cost;
} individual;

typedef vector<individual*> population_type;

population_type population;

Mode simMode;

int fitnessValue(int arrangement[])//low value is good
{



    //int fitness = (chessBoardSize * (chessBoardSize - 1)) / 2;          //initialize to a solution
    ////removing pairs that lie on the same row and on the same diagonal
    //for (int i = 0; i < chessBoardSize; i++)
    //    for (int j = i + 1; j < chessBoardSize; j++)
    //        if ((arrangement[i] == arrangement[j]) || (i - arrangement[i] == j - arrangement[j]) || (i + arrangement[i] == j + arrangement[j]))
    //            fitness--;
    //return fitness;
}

individual* createNode()
{
    individual* newNode = new individual;
    return newNode;
}

// Returns a random value between Min and Max (both inclusive)
int random(int min, int max)
{
    // Mersenne-Twister
    double r = (double)max - (double)min + 1;
    return min + (int)(r * rng() / (std::mt19937::max() + 1.0));
}

void generatePopulation()
{
    const int arraySize = chessBoardSize * chessBoardSize;
    int sampleArrangement[arraySize];
    individual* temp;
    for (int i = 0; i < arraySize; i++)
    {
        sampleArrangement[i] = i + 1;
    }

    // adds entries to population list
    for (int i = 0; i < initialPopulationCount; i++)
    {
        // Permute In Place (Random Shuffle)
        for (int j = 0; j < arraySize; j++)
        {
            std::swap(sampleArrangement[j], sampleArrangement[random(j, arraySize-1)]);
        }
        temp = createNode();
        for (int j = 0; j < arraySize; j++)
        {
            temp->arrangement[j] = sampleArrangement[j];
            //cout << sampleArrangement[j] << " ";
        }
        //cout << std::endl;
        temp->cost = fitnessValue(sampleArrangement);
        population.push_back(temp);
    }
}

int* createInversionSequence(individual* individual)
{
    int* inversionSequence = new int[fieldSize];

    for (int i = 0; i < fieldSize; i++)
    {
        int counter = 0;

        for (int value : individual->arrangement)
        {
            if (value == i + 1)
                break;

            if (value > i + 1)
                counter++;
        }
        inversionSequence[i] = counter;
    }

    return inversionSequence;
}

individual* reproduce(individual* parent1, individual* parent2)
{
    individual* child = createNode();
    int childArrangement[chessBoardSize * chessBoardSize];

    int* inversionSequenceP1;
    inversionSequenceP1 = createInversionSequence(parent1);

    for (int i = 0; i < fieldSize; i++)
    {
        std::cout << inversionSequenceP1[i] << " | ";
    }
    std::cout << std::endl;

    int* inversionSequenceP2;
    inversionSequenceP2 = createInversionSequence(parent2);

    for (int i = 0; i < fieldSize; i++)
    {
        std::cout << inversionSequenceP2[i] << " | ";
    }
    std::cout << std::endl;

    int crossoverPoint = 4 + 1; //random(2, 6) + 1;
    int* inversionSequenceChild = new int[fieldSize];

    std::copy(inversionSequenceP1, inversionSequenceP1 + crossoverPoint, inversionSequenceChild);
    std::copy(inversionSequenceP2 + crossoverPoint, inversionSequenceP2 + fieldSize, inversionSequenceChild + crossoverPoint );

    for (int i = 0; i < fieldSize; i++)
    {
        std::cout << inversionSequenceChild[i] << " | ";
    }

    //int n = chessBoardSize;
    //int c = rand() % n;
    // child->arrangement = (x->arrangement).substr(0, c) + (y->arrangement).substr(c, n - c + 1);
    //child->cost = fitnessValue(child->arrangement);
    //return child;

    return nullptr;
}

individual* mutate(individual* child)
{
    int randomQueen = random(0, chessBoardSize);
    int randomPosition = random(0, chessBoardSize);
    child->arrangement[randomQueen] = randomPosition + 48; //48?
    return child;
}

int randomSelection()
{
    int randomPos = random(0, population.size()) % 2;  // TODO: Adjust
    return randomPos;
}

bool isFit(individual* test)
{
    if (fitnessValue(test->arrangement) == ((chessBoardSize * (chessBoardSize - 1)) / 2))
        return true;
    return false;
}

bool comp(individual* a, individual* b)
{
    return(a->cost > b->cost);
}

individual* GA()
{
    int randomNum1 = 0;
    int randomNum2 = 0;
    individual* individualX = nullptr;
    individual* individualY = nullptr;
    individual* child = nullptr;

    bool found = false;
    while (!found)
    {
        population_type new_population;
        for (unsigned int i = 0; i < population.size(); i++)
        {
            sort(population.begin(), population.end(), comp);

            randomNum1 = randomSelection();
            individualX = population[randomNum1];

            randomNum2 = randomSelection();
            individualY = population[randomNum2];

            child = reproduce(individualX, individualY);

            if (random(0, 1))     //random probability
                child = mutate(child);

            if (isFit(child))
            {
                found = 1;
                return child;
            }
            new_population.push_back(child);
        }
        population = new_population;
    }
    return child;
}

void initialize()
{
    srand(time(0));     // no longer neeeded -> use custom random() function instead of rand()
    //chessBoardSize = 3; // TODO: Read from Commandline
}

int main()
{
    ////Reproduction Testing Stuff
    //individual* Testsubject1 = new individual;;
    //individual* Testsubject2 = new individual;;

    //for (int i = 0; i < fieldSize; i++)
    //{
    //    Testsubject1->arrangement[i] = i + 1;
    //    Testsubject2->arrangement[i] = 9 - i;
    //}
    //
    //for (int value : Testsubject1->arrangement)
    //{
    //    std::cout << value << " | ";
    //}
    //std::cout << std::endl;

    //for (int value : Testsubject2->arrangement)
    //{
    //    std::cout << value << " | ";
    //}
    //std::cout << std::endl;

    //reproduce(Testsubject1, Testsubject2);
    ////Reproduction Testing Stuff End

    initialize();
    clock_t start_time, end_time;           //to keep a track of the time spent in computing
    //map<string, int> solutionsFound;
    int maxSolutions = 92, numFound = 0;       //already known that 92 solutions exist for 8 Queen Problem!
    start_time = clock();
    cout << "*Returns the column number corresponding to the row at the index*" << endl << endl;
    while (numFound != maxSolutions)
    {
        generatePopulation();
        cout << "Generated Population successfully." << endl;
        /*individual* solution = */GA();
        // TODO: Print Solutions
        //if (!solutionsFound[solution->arrangement])
        //{
        //    solutionsFound[solution->arrangement] = 1;
        //    cout << "Possible Solution #" << (++numFound) << ":\t" << solution->arrangement << endl;
        //}
    }
    end_time = clock();

    cout << "Time required for execution: \t" << 1000 * ((double)(end_time - start_time) / CLOCKS_PER_SEC) << " milliseconds." << "\n\n";
    return 0;
}
