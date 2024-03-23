#include <string.h>

#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <vector>
#define INT_MAX 2147483647
#define INT_MIN -2147483647
using namespace std;

int getNumber(string);
void print_edge();
void print_matrix();
void dfsCountConnectedGraph(int node, int numNodes, int edge[][9999],
                            vector<bool>& visited, int count);
int countConnectedGraph(int numNodes, int edge[][9999]);
bool isBipartite(unordered_map<int, vector<int> >& graph);
bool dfsBipartite(int node, int currentColor, unordered_map<int, int>& color,
                  unordered_map<int, vector<int> >& graph);
void paintColor(vector<int>, int);

class Component {
   public:
    int x1, y1, x2, y2;
    int color;
    int id;
    int graphId;
    bool graphConflicted;
    Component(int x1, int y1, int x2, int y2, int color, int id) {
        this->x1 = x1;
        this->y1 = y1;
        this->x2 = x2;
        this->y2 = y2;
        this->color = color;
        this->id = id;
    }
    int top() { return y2; }

    int right() { return x2; }

    int bottom() { return y1; }

    int left() { return x1; }

    Component() {
        this->x1 = 0;
        this->y1 = 0;
        this->x2 = 0;
        this->y2 = 0;
        this->color = 0;
        this->id = 0;
        this->graphId = -1;
        this->graphConflicted = false;
    }
};
class DensityWindow {
   public:
    int x1, y1, x2, y2;
    int greenArea;
    int blueArea;
    double greenDensity;
    double blueDensity;
    DensityWindow(int x1, int y1, int x2, int y2) {
        this->x1 = x1;
        this->y1 = y1;
        this->x2 = x2;
        this->y2 = y2;
    }
    DensityWindow() {
        this->x1 = 0;
        this->y1 = 0;
        this->x2 = 0;
        this->y2 = 0;
    }
};
Component component[1000];
int edge[9999][9999];        // adjency matrix
int numberOfComponents = 1;  // number of components
// int color[99] = {0, 0, 0, 0, 0, 0, 2, 1, 1, 2,
//                  2, 1, 1, 2, 1, 2, 2, 1, 2, 1};  // 0 no color, 1 color
//                  A(green), 2
//                                               // color B(blue)

int main() {
    bool flag = true;
    freopen("input.in", "r", stdin);
    freopen("output.txt", "w", stdout);
    string alpha_input;
    cin >> alpha_input;
    int alpha = getNumber(alpha_input);
    // 取得alpha
    string beta_input;
    cin >> beta_input;
    int beta = getNumber(beta_input);
    // 取得beta
    string omega_input;
    cin >> omega_input;
    int omega = getNumber(omega_input);
    // 取得omega

    int a, b, c, d;
    char cc;
    while (cin >> a >> cc >> b >> cc >> c >> cc >> d) {
        component[numberOfComponents].x1 = a;  // x1 x1<x2 bottom left
        component[numberOfComponents].y1 = b;  // y1 y1<y2 bottom left
        component[numberOfComponents].x2 = c;  // x2 top right
        component[numberOfComponents].y2 = d;  // y2 top right
        // component[numberOfComponents].color = color[numberOfComponents];
        numberOfComponents += 1;
    }
    numberOfComponents--;  // fix the number of components

    // create the edges
    memset(edge, 0, sizeof(edge));
    int edge_num = 0;
    for (int i = 1; i <= numberOfComponents; i++) {
        for (int j = 1; j <= numberOfComponents; j++) {
            int delta = 0;
            if (i == j) continue;
            if (!(component[i].bottom() > component[j].top()) &&
                !(component[i].top() < component[j].bottom())) {
                if (component[i].left() > component[j].right())
                    delta = component[i].left() - component[j].right();
                else if (component[i].right() < component[j].left())
                    delta = component[j].left() - component[i].right();
                else
                    continue;

                if (delta <= alpha) {
                    edge_num++;
                    edge[i][j] = 1;
                    edge[j][i] = 1;
                } else
                    continue;
            } else if (!(component[i].left() > component[j].right()) &&
                       !(component[i].right() < component[j].left())) {
                if (component[i].bottom() > component[j].top())
                    delta = component[i].bottom() - component[j].top();
                else if (component[i].top() < component[j].bottom())
                    delta = component[j].bottom() - component[i].top();
                else
                    continue;

                if (delta <= beta) {
                    edge_num++;
                    edge[i][j] = 1;
                    edge[j][i] = 1;
                } else
                    continue;
            }
        }
    }
    print_edge();
    int numConnectedGraph = countConnectedGraph(numberOfComponents, edge) - 1;
    cout << "numConnectedGraph: " << numConnectedGraph << endl;
    int numConnectedGraphConflict = 0;
    // check conflict
    for (int i = 1; i <= numConnectedGraph; i++) {
        unordered_map<int, vector<int> > subGraph;
        for (int j = 1; j <= numberOfComponents; j++) {
            for (int k = 1; k <= numberOfComponents; k++) {
                if (component[j].graphId == i && edge[j][k] == 1) {
                    subGraph[j].push_back(k);
                    subGraph[k].push_back(j);
                }
            }
        }
        bool isConflict = !isBipartite(subGraph);  // not conflict == Bitpartite
        numConnectedGraphConflict += isConflict;
        for (int j = 1; j <= numberOfComponents; j++) {
            if (component[j].graphId == i)
                component[j].graphConflicted = isConflict;
        }
        subGraph.clear();
    }

    // create bounding box
    int BoundingBox_x1 = INT_MAX, BoundingBox_y1 = INT_MAX,
        BoundingBox_x2 = INT_MIN, BoundingBox_y2 = INT_MIN;
    // x1 x1<x2 bottom left
    // y1 y1<y2 bottom left
    // x2 top right
    // y2 top right
    for (int i = 1; i <= numberOfComponents; i++) {
        if (!component[i].graphConflicted) {
            if (component[i].left() < BoundingBox_x1)
                BoundingBox_x1 = component[i].left();
            if (component[i].bottom() < BoundingBox_y1)
                BoundingBox_y1 = component[i].bottom();
            if (component[i].right() > BoundingBox_x2)
                BoundingBox_x2 = component[i].right();
            if (component[i].top() > BoundingBox_y2)
                BoundingBox_y2 = component[i].top();
        }
    }
    vector<DensityWindow> densitywindow;

    // Iterator for bounding box (color density window (window size = omega))
    for (int i = BoundingBox_y1; i < BoundingBox_y2;) {
        for (int j = BoundingBox_x1; j < BoundingBox_x2;) {
            DensityWindow temp;
            temp.x1 = j;
            temp.y1 = i;
            temp.x2 = j + omega;
            temp.y2 = i + omega;
            densitywindow.push_back(temp);
            // cout << j << " " << i << endl;
            if (j == BoundingBox_x2 - omega) break;
            if (j + omega + omega > BoundingBox_x2) {
                j = BoundingBox_x2 - omega;
            } else
                j += omega;
        }
        if (i == BoundingBox_y2 - omega) break;
        if (i + omega + omega > BoundingBox_y2) {
            i = BoundingBox_y2 - omega;
        } else
            i += omega;
    }
    // initialize color
    srand(time(NULL));
    vector<int> ansInit;
    for (int i = 0; i < numConnectedGraph - numConnectedGraphConflict; i++) {
        ansInit.push_back(rand() % 2 + 1);
        // cout << int(ans[i]) << endl;
    }

    // paintColor(ansInit, numberOfComponents);

    // puts("----------------");

    cout << fixed;
    cout << setprecision(2);
    // cout << calculateFitness(densitywindow, omega) << endl;
    double tInitial = 3000;
    double tFinal = 0.1;
    int nMarkov = 100;
    double alfa = 0.98;
    void simulatedAnnealing(double tInitial, double tFinal, int nMarkov,
                            double alfa, vector<int> ansInit,
                            vector<DensityWindow> densitywindow, int omega);
    simulatedAnnealing(tInitial, tFinal, nMarkov, alfa, ansInit, densitywindow,
                       omega);
    for (int i = 1; i <= numberOfComponents; i++) {
        cout << component[i].color << endl;
    }
}

void simulatedAnnealing(double tInitial, double tFinal, int nMarkov,
                        double alfa, vector<int> ansInit,
                        vector<DensityWindow> densitywindow, int omega) {
    double calculateFitness(vector<DensityWindow> densitywindow, int omega);
    int t = tInitial;
    vector<int> bestAns = ansInit;
    vector<int> currentAns = bestAns;
    paintColor(ansInit, numberOfComponents);
    double bestCost = calculateFitness(densitywindow, omega);
    double currentCost = bestCost;
    double answerCost = 0;
    vector<int> answerAns;
    while (t > tFinal) {
        for (int i = 0; i < nMarkov; i++) {
            // generate a new solution
            currentAns[rand() % bestAns.size()] = rand() % 2 + 1;
            // calculate the cost of the new solution
            paintColor(currentAns, numberOfComponents);
            currentCost = calculateFitness(densitywindow, omega);
            // calculate the cost of the current solution
            // if the new solution is better than the current solution, accept
            // the new solution if the new solution is worse than the current
            // solution, accept the new solution with a probability of
            // e^((cost_new-cost_current)/t)
            // Metropolis rule
            if (currentCost > bestCost) {
                bestCost = currentCost;
                bestAns = currentAns;
            } else if (exp((currentCost - bestCost) / t) >
                       (rand() % 100) / 100) {
                bestCost = currentCost;
                bestAns = currentAns;
            }
            if (answerCost < bestCost) {
                answerCost = bestCost;
                answerAns = bestAns;
            }
        }

        t *= alfa;
    }
    cout << "bestCost: " << answerCost << endl;
    paintColor(answerAns, numberOfComponents);
}
void paintColor(vector<int> ans, int numberOfComponents) {
    void dfsPaintColor(int numberOfComponents, vector<bool>& visited, int node,
                       int color);
    vector<bool> visited(numberOfComponents, false);
    for (int i = 1, s = 0; i <= numberOfComponents; i++) {
        if (visited[i] == false && component[i].graphConflicted == false) {
            visited[i] = true;
            component[i].color = ans[s];
            s++;
            dfsPaintColor(numberOfComponents, visited, i, component[i].color);
        }
    }
}
void dfsPaintColor(int numberOfComponents, vector<bool>& visited, int node,
                   int color) {
    visited[node] = true;
    component[node].color = color;
    for (int neighbor = 1; neighbor <= numberOfComponents; ++neighbor) {
        if (edge[node][neighbor] == 1) {
            if (!visited[neighbor]) {
                dfsPaintColor(numberOfComponents, visited, neighbor,
                              (color + 2) % 2 + 1);
            }
        }
    }
}
double calculateFitness(vector<DensityWindow> densitywindow, int omega) {
    // Calculate the color density of each density window
    for (auto& dw : densitywindow) {
        int greenArea = 0, blueArea = 0;
        for (int j = 1; j <= numberOfComponents; j++) {
            if (component[j].graphConflicted) continue;
            // id 1: all area in window
            if (dw.x1 <= component[j].x1 && component[j].x1 <= dw.x2 &&
                dw.y1 <= component[j].y1 && component[j].y1 <= dw.y2 &&
                dw.x1 <= component[j].x2 && component[j].x2 <= dw.x2 &&
                dw.y1 <= component[j].y2 && component[j].y2 <= dw.y2) {
                if (component[j].color == 1)
                    greenArea += abs(component[j].y2 - component[j].y1) *
                                 abs(component[j].x2 - component[j].x1);
                else if (component[j].color == 2)
                    blueArea += abs(component[j].y2 - component[j].y1) *
                                abs(component[j].x2 - component[j].x1);
            }
            // id 2:top right corner
            else if (dw.x1 <= component[j].x1 && component[j].x1 <= dw.x2 &&
                     dw.y1 <= component[j].y1 && component[j].y1 <= dw.y2 &&
                     dw.y2 < component[j].y2 && dw.x2 < component[j].x2) {
                if (component[j].color == 1)
                    greenArea += abs(dw.y2 - component[j].y1) *
                                 abs(dw.x2 - component[j].x1);
                else if (component[j].color == 2)
                    blueArea += abs(dw.y2 - component[j].y1) *
                                abs(dw.x2 - component[j].x1);
            }
            // id 3:right side
            else if (dw.x1 <= component[j].x1 && component[j].x1 <= dw.x2 &&
                     dw.y1 <= component[j].y1 && component[j].y1 <= dw.y2 &&
                     dw.y1 <= component[j].y2 && component[j].y2 <= dw.y2 &&
                     dw.x2 < component[j].x2) {
                if (component[j].color == 1)
                    greenArea += abs(component[j].y2 - component[j].y1) *
                                 abs(dw.x2 - component[j].x1);
                else if (component[j].color == 2)
                    blueArea += abs(component[j].y2 - component[j].y1) *
                                abs(dw.x2 - component[j].x1);
            }
            // id 4: bottom right side
            else if (dw.x1 <= component[j].x1 && component[j].x1 <= dw.x2 &&
                     dw.y1 <= component[j].y2 && component[j].y2 <= dw.y2 &&
                     dw.x2 < component[j].x2 && component[j].y1 < dw.y1) {
                if (component[j].color == 1)
                    greenArea += abs(dw.x2 - component[j].x1) *
                                 abs(component[j].y2 - dw.y1);
                else if (component[j].color == 2)
                    blueArea += abs(dw.x2 - component[j].x1) *
                                abs(component[j].y2 - dw.y1);
            }
            // id 5:bottom side
            else if (dw.x1 <= component[j].x1 && component[j].x1 <= dw.x2 &&
                     dw.x1 <= component[j].x2 && component[j].x2 <= dw.x2 &&
                     dw.y1 <= component[j].y2 && component[j].y2 <= dw.y2 &&
                     component[j].y1 < dw.y1) {
                if (component[j].color == 1)
                    greenArea += abs(component[j].y2 - dw.y1) *
                                 abs(component[j].x2 - component[j].x1);
                else if (component[j].color == 2)
                    blueArea += abs(component[j].y2 - dw.y1) *
                                abs(component[j].x2 - component[j].x1);
            }
            // id 6:bottom left corner
            else if (dw.x1 <= component[j].x2 && component[j].x2 <= dw.x2 &&
                     dw.y1 <= component[j].y2 && component[j].y2 <= dw.y2 &&
                     component[j].y1 < dw.y1 && component[j].x1 < dw.x1) {
                if (component[j].color == 1)
                    greenArea += abs(component[j].y2 - dw.y1) *
                                 abs(component[j].x2 - dw.x1);
                else if (component[j].color == 2)
                    blueArea += abs(component[j].y2 - dw.y1) *
                                abs(component[j].x2 - dw.x1);
            }
            // id 7:left side
            else if (dw.y1 <= component[j].y1 && component[j].y1 <= dw.y2 &&
                     dw.x1 <= component[j].x2 && component[j].x2 <= dw.x2 &&
                     dw.y1 <= component[j].y2 && component[j].y2 <= dw.y2 &&
                     component[j].x1 < dw.x1) {
                if (component[j].color == 1)
                    greenArea += abs(component[j].x2 - dw.x1) *
                                 abs(component[j].y2 - component[j].y1);
                else if (component[j].color == 2)
                    blueArea += abs(component[j].x2 - dw.x1) *
                                abs(component[j].y2 - component[j].y1);
            }
            // id 8:top left corner
            else if (dw.x1 <= component[j].x2 && component[j].x2 <= dw.x2 &&
                     dw.y1 <= component[j].y1 && component[j].y1 <= dw.y2 &&
                     component[j].x1 < dw.x1 && dw.y2 < component[j].y2) {
                if (component[j].color == 1)
                    greenArea += abs(dw.y2 - component[j].y1) *
                                 abs(component[j].x2 - dw.x1);
                else if (component[j].color == 2)
                    blueArea += abs(dw.y2 - component[j].y1) *
                                abs(component[j].x2 - dw.x1);
            }
            // id 9:top side
            else if (dw.x1 <= component[j].x1 && component[j].x1 <= dw.x2 &&
                     dw.x1 <= component[j].x2 && component[j].x2 <= dw.x2 &&
                     dw.y1 <= component[j].y1 && component[j].y1 <= dw.y2 &&
                     dw.y2 < component[j].y2) {
                if (component[j].color == 1)
                    greenArea += abs(dw.y2 - component[j].y1) *
                                 abs(component[j].x2 - component[j].x1);
                else if (component[j].color == 2)
                    blueArea += abs(dw.y2 - component[j].y1) *
                                abs(component[j].x2 - component[j].x1);
            }
        }
        dw.greenArea = greenArea;
        dw.blueArea = blueArea;
        dw.greenDensity = (double)greenArea / (omega * omega) * 100;
        dw.blueDensity = (double)blueArea / (omega * omega) * 100;
    }
    // puts("----");
    double cost = 30.0;
    for (auto& dw : densitywindow) {
        cost += fabs(70.0 / densitywindow.size() -
                     (fabs(dw.greenDensity - dw.blueDensity) / 5));
    }
    // cout << cost << endl;
    return cost;
}
int getNumber(string st) {
    int ret = -1;
    size_t equalPos = st.find("=");
    string numberPart = st.substr(equalPos + 1);
    ret = atoi(numberPart.c_str());
    return ret;
}
void print_edge() {
    for (int i = 1; i <= numberOfComponents; i++) {
        for (int j = 1; j <= numberOfComponents; j++) {
            // cout << edge[i][j];
            if (i <= j && edge[i][j] == 1)
                cout << "e " << i << ' ' << j << '\n';
        }
        // cout << endl;
    }
}
void print_matrix() {  // print edge
    for (int i = 1; i <= numberOfComponents; i++) {
        for (int j = 1; j <= numberOfComponents; j++) {
            cout << edge[i][j] << ' ';
        }
        cout << '\n';
    }
}
void dfsCountConnectedGraph(int node, int numNodes, int edge[][9999],
                            vector<bool>& visited, int count) {
    visited[node] = true;
    component[node].graphId = count;
    for (int neighbor = 1; neighbor <= numNodes; ++neighbor) {
        if (edge[node][neighbor] == 1) {
            if (!visited[neighbor]) {
                dfsCountConnectedGraph(neighbor, numNodes, edge, visited,
                                       count);
            }
            // dfs(neighbor, numNodes, edge, visited);
        }
    }
}
int countConnectedGraph(int numNodes, int edge[][9999]) {
    vector<bool> visited(numNodes, false);
    int count = 1;

    for (int node = 1; node <= numNodes; ++node) {
        if (!visited[node]) {
            dfsCountConnectedGraph(node, numNodes, edge, visited, count);
            count++;
        }
    }
    return count;
}

/**
 * Determines if a given graph is bipartite.
 * A graph is bipartite if its vertices can be divided into two disjoint sets
 * such that every edge connects a vertex from one set to another. This function
 * uses depth-first search (DFS) to check if the graph is bipartite.
 *
 * @param graph The graph represented as an unordered map, where the key is the
 * node and the value is a vector of its neighboring nodes.
 * @return True if the graph is bipartite, false otherwise.
 */
bool isBipartite(unordered_map<int, vector<int> >& graph) {
    unordered_map<int, int> color;

    for (auto& node : graph) {
        if (color.find(node.first) == color.end()) {
            int currentColor = 1;
            if (!dfsBipartite(node.first, currentColor, color, graph)) {
                return false;
            }
        }
    }

    return true;
}

/**
 * Performs a depth-first search to determine if a graph is bipartite.
 *
 * @param node The current node being visited.
 * @param currentColor The color assigned to the current node.
 * @param color A map that stores the color assigned to each node.
 * @param graph A map that represents the graph structure.
 * @return True if the graph is bipartite, false otherwise.
 */
bool dfsBipartite(int node, int currentColor, unordered_map<int, int>& color,
                  unordered_map<int, vector<int> >& graph) {
    color[node] = currentColor;

    for (int neighbor : graph[node]) {
        if (color.find(neighbor) == color.end()) {
            if (!dfsBipartite(neighbor, -currentColor, color, graph)) {
                return false;
            }
        } else if (color[neighbor] == color[node]) {
            return false;
        }
    }

    return true;
}
