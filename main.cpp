#include <string.h>

#include <algorithm>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <vector>
using namespace std;

int getNumber(string);
void print_edge();
void print_matrix();
void dfs(int node, int numNodes, int edge[][9999], vector<bool>& visited,
         int count);
int countConnectedGraph(int numNodes, int edge[][9999]);
bool isBipartite(unordered_map<int, vector<int> >& graph);
bool dfsBipartite(int node, int currentColor, unordered_map<int, int>& color,
                  unordered_map<int, vector<int> >& graph);

class Component {
   public:
    int x1, y1, x2, y2;
    int connectwith;
    int color;
    int id;
    int graphId;
    bool graphConflicted;
    Component(int x1, int y1, int x2, int y2, int color, int id,
              int connectwith) {
        this->x1 = x1;
        this->y1 = y1;
        this->x2 = x2;
        this->y2 = y2;
        this->color = color;
        this->id = id;
        this->connectwith = connectwith;
    }
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
int edge[9999][9999];  // adjency matrix
// int pos[9999][4];      // every component's position(x1,y1,x2,y2)
int numberOfComponents = 1;  // number of components
// bool visited[9999];  // is component visited?
int color[99] = {0, 0, 0, 0, 0, 0, 2, 1, 1, 2, 2,
                 1, 1, 2, 1, 2, 2, 1, 2, 1};  // 0 no color, 1 color A(green), 2
                                              // color B(blue)

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
        component[numberOfComponents].x1 = a;  // x1 x1<x2 buttom left
        component[numberOfComponents].y1 = b;  // y1 y1<y2 buttom left
        component[numberOfComponents].x2 = c;  // x2 top right
        component[numberOfComponents].y2 = d;  // y2 top right
        component[numberOfComponents].color = color[numberOfComponents];
        numberOfComponents += 1;
    }
    numberOfComponents--;  // fix the number of components

    // create the edges
    memset(edge, 0, sizeof(edge));
    int edge_num = 0;
    for (int i = 1; i <= numberOfComponents; i++) {
        for (int j = 1; j <= numberOfComponents; j++) {
            // //i是要比較的shape
            // //j是要跟他比的shape
            if (j == i) continue;

            if (component[i].x1 <= component[j].x1 &&
                component[j].x1 <=
                    component[i]
                        .x2) {  // 看j shape的兩個點的x有沒有在i shape的範圍內
                int y_distance =
                    min(min(abs(component[i].y1 - component[j].y2),
                            abs(component[i].y1 - component[j].y1)),
                        min(abs(component[i].y2 - component[j].y1),
                            abs(component[i].y2 - component[j].y2)));
                if (y_distance <= beta) {
                    edge_num++;
                    edge[i][j] = 1;
                    edge[j][i] = 1;
                }
            } else if (component[i].x1 <= component[j].x2 &&
                       component[j].x2 <=
                           component[i].x2) {  // 看j shape的兩個點的x有沒有在i
                                               // shape的範圍內
                int y_distance =
                    min(min(abs(component[i].y1 - component[j].y2),
                            abs(component[i].y1 - component[j].y1)),
                        min(abs(component[i].y2 - component[j].y1),
                            abs(component[i].y2 - component[j].y2)));
                if (y_distance <= beta) {
                    edge_num++;
                    edge[i][j] = 1;
                    edge[j][i] = 1;
                }
            } else if (component[i].y1 <= component[j].y1 &&
                       component[j].y1 <=
                           component[i].y2) {  // 看j shape的兩個點的y有沒有在i
                                               // shape的範圍內
                int x_distance =
                    min(min(abs(component[i].x1 - component[j].x1),
                            abs(component[i].x1 - component[j].x2)),
                        min(abs(component[i].x2 - component[j].x1),
                            abs(component[i].x2 - component[j].x2)));
                if (x_distance <= alpha) {
                    edge_num++;
                    edge[i][j] = 1;
                    edge[j][i] = 1;
                }
            } else if (component[i].y1 <= component[j].y2 &&
                       component[j].y2 <=
                           component[i].y2) {  // 看j shape的兩個點的y有沒有在i
                                               // shape的範圍內
                int x_distance =
                    min(min(abs(component[i].x1 - component[j].x1),
                            abs(component[i].x1 - component[j].x2)),
                        min(abs(component[i].x2 - component[j].x1),
                            abs(component[i].x2 - component[j].x2)));
                if (x_distance <= alpha) {
                    edge_num++;
                    edge[i][j] = 1;
                    edge[j][i] = 1;
                }
            }
            /*
            x
            pos[i][0] <= pos[j][0] <= pos[i][2];
            pos[i][0] <= pos[j][2] <= pos[i][2];

            y
            pos[i][1] <= pos[j][1] <= pos[i][3];
            pos[i][1] <= pos[j][3] <= pos[i][3];
            */
        }
    }

    int numConnectedGraph = countConnectedGraph(numberOfComponents, edge) - 1;
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
        for (int j = 1; j <= numberOfComponents; j++) {
            if (component[j].graphId == i)
                component[j].graphConflicted = isConflict;
        }
        subGraph.clear();
    }

    // create bounding box
    int BoundingBox_x1 = INT_MAX, BoundingBox_y1 = INT_MAX,
        BoundingBox_x2 = INT_MIN, BoundingBox_y2 = INT_MIN;
    // x1 x1<x2 buttom left
    // y1 y1<y2 buttom left
    // x2 top right
    // y2 top right
    for (int i = 1; i <= numberOfComponents; i++) {
        if (!component[i].graphConflicted) {
            if (component[i].x1 < BoundingBox_x1)
                BoundingBox_x1 = component[i].x1;
            if (component[i].y1 < BoundingBox_y1)
                BoundingBox_y1 = component[i].y1;
            if (component[i].x2 > BoundingBox_x2)
                BoundingBox_x2 = component[i].x2;
            if (component[i].y2 > BoundingBox_y2)
                BoundingBox_y2 = component[i].y2;
        }
    }
    // cout << BoundingBox_x1 << ' ' << BoundingBox_y1 << ' ' << BoundingBox_x2
    //      << ' ' << BoundingBox_y2 << '\n';
    // 540 0
    // 960 0
    // 540 360
    // 960 360
    DensityWindow DensityWindow[1000];
    int numberOfDensityWindows = 0;

    // Iterator for bounding box (color density window (window size = omega))
    for (int i = BoundingBox_y1; i < BoundingBox_y2;) {
        for (int j = BoundingBox_x1; j < BoundingBox_x2;) {
            DensityWindow[numberOfDensityWindows].x1 = j;
            DensityWindow[numberOfDensityWindows].y1 = i;
            DensityWindow[numberOfDensityWindows].x2 = j + omega;
            DensityWindow[numberOfDensityWindows].y2 = i + omega;
            cout << j << " " << i << endl;
            numberOfDensityWindows++;
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

    //  print_edge();

    //---output with edge format---
    // print_matrix();
    //---output with matrix format---
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
void dfs(int node, int numNodes, int edge[][9999], vector<bool>& visited,
         int count) {
    visited[node] = true;
    component[node].graphId = count;
    for (int neighbor = 1; neighbor <= numNodes; ++neighbor) {
        if (edge[node][neighbor] == 1) {
            if (!visited[neighbor]) {
                dfs(neighbor, numNodes, edge, visited, count);
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
            dfs(node, numNodes, edge, visited, count);
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
