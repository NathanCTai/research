#include <vector>
#include <random>
#include <fstream>
#include <iostream>
using namespace std;

int T = 100, P = 21, C = 0, dumbRuns = 1000; float L = 0.0f; 
vector<vector<float>> lambda(T), V(T + 1); // λ_t(p) is T x P; 
vector<float> alpha(T), gamma(T), prices(P); // V_t(x) is (T + 1) x C;
ofstream priceFile("price.csv"),
         vtxFile("vtx.csv"),
         dumbFile("dumb.csv"),
         demandFile("demand.csv"),
         splitFile("split.csv"),
         data("data.json");
random_device rd; unsigned int seed = rd(); mt19937 gen(seed); 
uniform_real_distribution<float> alpha_dist(0.5f, 1.0f), gamma_dist(0.1f, 0.4f), coin(0.0f, 1.0f);

void readData() {
    priceFile << "x,t,opt\n"; vtxFile << "x,t,VTX\n"; 
    demandFile << "t,alpha,gamma\n"; dumbFile << "p,rev\n";
    for (int i = 0; i < T; i++) {
        alpha[i] = alpha_dist(gen); gamma[i] = gamma_dist(gen);
        demandFile << i << "," << alpha[i] << "," << gamma[i] << "\n";
        lambda[i] = vector<float>(P);
        for (int j = 0; j < P; j++) {
            if (i == 0) {prices[j] = float(j);} // reading in prices
            lambda[i][j] = exp(-gamma[i] * prices[j]); 
        }
        L += alpha[i] * exp(-1); 
    } 
    C = int(L * 0.8); // can also be 0.6 and now we know C
    V = vector<vector<float>>(T + 1, vector<float>(C, 0.0f)); 
    data << "{\n\t\"seed\": " << seed << ",\n\t\"T\": " << T << ",\n\t\"C\": " << C << ",\n\t\"P\": " << P << "\n}";
}

void smartdp() {
    for (int t = T - 1; t > 0; t--) { 
        for (int x = 1; x < C; x++) {
            float maximum = 0.0; float optimalprice = 0.0; 
            for (int p = 0; p < P; p++) { 
                float lhs = alpha[t] * lambda[t][p] * (prices[p] + V[t + 1][x - 1]); // SALE
                float rhs = (1 - alpha[t] * lambda[t][p]) * V[t + 1][x]; // NO SALE
                if (lhs + rhs > maximum) {
                    maximum = rhs + lhs; 
                    optimalprice = prices[p];
                } 
            } priceFile << x << "," << t << "," << optimalprice << "\n";
            V[t][x] = maximum;
            vtxFile << x << "," << t << "," << V[t][x] << "\n";
        }
    }
}

void lessdumb() {
    int split = 20, b1, b2 = 0; float best_expec = 0.0f;
    for (int p1 : prices) {
        for (int p2 : prices) {
            float Theta1 = 0; float Theta2 = 0;
            for (int t = 0; t < split; t++) {
                Theta1 += (alpha[t] * lambda[t][p1]);
            }
            for (int t = split; t < T; t++) {
                Theta2 += (alpha[t] * lambda[t][p2]);
            }
            float exp = min(float(C), Theta1);
            float rev = p1 * exp + p2 * min(C - exp, Theta2);
            if (rev > best_expec) {
                b1 = p1; b2 = p2; best_expec = rev;
            }
        }
    }
    float rev = 0;
    for (int k = 0; k < dumbRuns; k++) {
        int qty = C;
        for (int t = 0; t < split; t++) {
            if (coin(gen) < alpha[t]) { // guy walks in
                if (coin(gen) < lambda[t][b1] and qty > 0) {
                    qty -= 1; rev += b1;
                }
            }
        }
        for (int t = split; t < T; t++) {
            if (coin(gen) < alpha[t]) {
                if (coin(gen) < lambda[t][b2] and qty > 0) {
                    qty -= 1; rev += b2;
                }
            }
        }
    }
}

void dumbsim() {
    vector<float> rev(P, 0.0f);
    for (int k = 0; k < dumbRuns; k++) {
        vector<int> qty(P, C);
        for (int t = 0; t < T; t++) {
            for (int p = 0; p < P; p++) {
                if (coin(gen) < alpha[t]) { // guy walks in 
                    if (coin(gen) < lambda[t][p] and qty[p] > 0) { // guy buys something
                        qty[p] -= 1; rev[p] += p;
                    }
                }
            }
        }
    }
    for (int p = 0; p < P; p++) {
        dumbFile << p << "," << rev[p] / dumbRuns << "\n";
    }
}

int main() {
    readData();
    smartdp();
    lessdumb();
    dumbsim();
    return 0;
}