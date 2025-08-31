#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;


int main() {
    int N, a, b;
    vector<int> chords;
    vector<vector<int>> M;

    // Get chords
    cin >> N;
    chords = vector<int>(N, 0);
    while (cin >> a) {
        if (!(cin >> b)){
            // Get 0
            break;
        }
        chords[a] = b;
        chords[b] = a;
    }
    
    M = vector<vector<int>>(N, vector<int>(N, 0));

    // Create table
    for (int step = 1; step < N; step++) {
        for (int s = 0; s + step < N; s++) {
            int e = s + step;

            M[s][e] = M[s+1][e];
            
            int k = chords[s];  // (s, k)
            if (k > s && k <= e) {
                // if k is in the range
                int l = (s+1 <= k-1) ? M[s+1][k-1] : 0;
                int r = (k+1 <= e)   ? M[k+1][e]   : 0;
                M[s][e] = max(M[s][e], 1 + l + r);
            }
        }
    }
    
    // Output data
    cout << M[0][N-1] << '\n';
    return 0;
}
