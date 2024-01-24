#include <string.h>

#include <algorithm>
#include <iostream>
#include <sstream>
#include <vector>
using namespace std;

int getNumber(string);
void print_edge();
void print_matrix();
void dfs(int node, int numNodes, int edge[][9999], vector<bool>& visited);
int countConnectedComponents(int numNodes, int edge[][9999]);

int edge[9999][9999];  // adjency matrix
int pos[9999][4];      // every component's position(x1,y1,x2,y2)
int cot = 1;           // number of components
bool visited[9999];    // is component visited?
int color[99] = {0, 0, 0, 0, 0, 0, 2, 1, 1, 2, 2,
                 1, 1, 2, 1, 2, 2, 1, 2, 1};  // 0 no color, 1 color A(green), 2
                                              // color B(blue)
bool isconfict[500] = {0};                    // graph isconfict?
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
        pos[cot][0] = a;  // x1 x1<x2 buttom left
        pos[cot][1] = b;  // y1 y1<y2 buttom left
        pos[cot][2] = c;  // x2 top right
        pos[cot][3] = d;  // y2 top right
        cot += 1;
    }
    // input data
    memset(edge, 0, sizeof(edge));
    int edge_num = 0;
    for (int i = 1; i <= cot; i++) {
        for (int j = 1; j <= cot; j++) {
            // //i是要比較的shape
            // //j是要跟他比的shape
            if (j == i) continue;
            // if(edge[j][i]==1) continue;
            // if(edge[i][j]==1) continue;
            if (pos[i][0] <= pos[j][0] &&
                pos[j][0] <=
                    pos[i][2]) {  // 看j shape的兩個點的x有沒有在i shape的範圍內
                int y_distance = min(
                    min(abs(pos[i][1] - pos[j][3]), abs(pos[i][1] - pos[j][1])),
                    min(abs(pos[i][3] - pos[j][1]),
                        abs(pos[i][3] - pos[j][3])));
                if (y_distance <= beta) {
                    edge_num++;
                    edge[i][j] = 1;
                    edge[j][i] = 1;
                }
            } else if (pos[i][0] <= pos[j][2] &&
                       pos[j][2] <=
                           pos[i][2]) {  // 看j shape的兩個點的x有沒有在i
                                         // shape的範圍內
                int y_distance = min(
                    min(abs(pos[i][1] - pos[j][3]), abs(pos[i][1] - pos[j][1])),
                    min(abs(pos[i][3] - pos[j][1]),
                        abs(pos[i][3] - pos[j][3])));
                if (y_distance <= beta) {
                    edge_num++;
                    edge[i][j] = 1;
                    edge[j][i] = 1;
                }
            } else if (pos[i][1] <= pos[j][1] &&
                       pos[j][1] <=
                           pos[i][3]) {  // 看j shape的兩個點的y有沒有在i
                                         // shape的範圍內
                int x_distance = min(
                    min(abs(pos[i][0] - pos[j][0]), abs(pos[i][0] - pos[j][2])),
                    min(abs(pos[i][2] - pos[j][0]),
                        abs(pos[i][2] - pos[j][2])));
                if (x_distance <= alpha) {
                    edge_num++;
                    edge[i][j] = 1;
                    edge[j][i] = 1;
                }
            } else if (pos[i][1] <= pos[j][3] &&
                       pos[j][3] <=
                           pos[i][3]) {  // 看j shape的兩個點的y有沒有在i
                                         // shape的範圍內
                int x_distance = min(
                    min(abs(pos[i][0] - pos[j][0]), abs(pos[i][0] - pos[j][2])),
                    min(abs(pos[i][2] - pos[j][0]),
                        abs(pos[i][2] - pos[j][2])));
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
    print_edge();
    for (int i = 1; i < cot; i++) {
        cout << i << ": " << color[i] << endl;
    }
    //---output with edge format---
    //
    // print_matrix();
    //---output with matrix format---
    int numConnectedComponents = countConnectedComponents(cot - 1, edge);
    cout << numConnectedComponents;
}

int getNumber(string st) {
    int ret = -1;
    size_t equalPos = st.find("=");
    string numberPart = st.substr(equalPos + 1);
    ret = atoi(numberPart.c_str());
    return ret;
}
void print_edge() {
    for (int i = 1; i <= cot; i++) {
        for (int j = 1; j <= cot; j++) {
            // cout << edge[i][j];
            if (i <= j && edge[i][j] == 1)
                cout << "e " << i << ' ' << j << '\n';
        }
        // cout << endl;
    }
}
void print_matrix() {
    for (int i = 1; i <= cot; i++) {
        for (int j = 1; j <= cot; j++) {
            cout << edge[i][j] << ' ';
        }
        cout << '\n';
    }
}
void dfs(int node, int numNodes, int edge[][9999], vector<bool>& visited) {
    visited[node] = true;
    for (int neighbor = 1; neighbor <= numNodes; ++neighbor) {
        if (edge[node][neighbor] == 1 && !visited[neighbor]) {
            dfs(neighbor, numNodes, edge, visited);
        }
    }
}
int countConnectedComponents(int numNodes, int edge[][9999]) {
    vector<bool> visited(numNodes, false);
    int count = 0;

    for (int node = 1; node <= numNodes; ++node) {
        if (!visited[node]) {
            dfs(node, numNodes, edge, visited);
            count++;
        }
    }
    return count;
}
// check conflicts

// get valid connected components and check boundaries(top right and buttom
// left)

// iterate density window to calulate A and B
//  cost function
