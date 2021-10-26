#include <iostream>
#include <unordered_map>
#include <vector>

using namespace std;

int main(void) {
    vector<unordered_map<uint32_t, float>> sim_values;
    sim_values.resize(10);
    cout << sim_values.size();
    for (int i = 0; i < sim_values.size(); i++) {
        for (int j = 0; j < 3; j++) {
            cout << "Hi" << endl;
            if (sim_values[i].find(j) == sim_values[i].end()) {
                cout << "here" << endl;
            }
        }
    }
    return 0;
}