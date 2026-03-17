#include <vector>
#include <random>
#include <fstream>
#include <iostream>
using namespace std;

const int T = 100, P = 21, split = 20, dumbRuns = 1000, mainRuns = 10; int C = 0; 
vector<vector<float>> lambda(T), V(T + 1); // λ_t(p) is T x P; read in theta from the start
vector<float> alpha(T, 0.0f), gamma(T, 0.0f), prices(P, 0.0f), theta1(P, 0.0f), theta2(P, 0.0f);
ofstream data("data.csv"); random_device rd; unsigned int seed = rd(); mt19937 gen(seed); 
uniform_real_distribution<float> alpha_dist(0.5f, 1.0f), gamma_dist(0.1f, 0.4f), coin(0.0f, 1.0f);

void readData() {
    for (int p = 0; p < P; p++) {
        prices[p] = float(p);
        for (int t = 0; t < T; t++) {
            if (p == 0) {alpha[t] = alpha_dist(gen); gamma[t] = gamma_dist(gen); lambda[t].resize(P); }
            lambda[t][p] = exp(-gamma[t] * prices[p]); 
            if (t <  split) { theta1[p] += alpha[t] * lambda[t][p]; }
            if (t >= split) { theta2[p] += alpha[t] * lambda[t][p]; } // new way of calculating theta
        }
    }
    C = int(accumulate(alpha.begin(), alpha.end(), 0.0f) * 0.29); data << seed << "," << C << ",";
    V = vector<vector<float>>(T + 1, vector<float>(C, 0.0f)); 
}

void smartdp() {
    for (int t = T - 1; t > -1; t--) { 
        for (int x = 1; x < C; x++) {
            float maximum = 0.0; float optimalprice = 0.0; 
            for (int p = 0; p < P; p++) { 
                float lhs = alpha[t] * lambda[t][p] * (prices[p] + V[t + 1][x - 1]); // SALE
                float rhs = (1 - alpha[t] * lambda[t][p]) * V[t + 1][x]; // NO SALE
                if (lhs + rhs > maximum) {
                    maximum = rhs + lhs; optimalprice = prices[p];
                } 
            } V[t][x] = maximum;
        }
    } data << V[0][C - 1] << ",";
}

void splitsim() {
    vector<vector<float>> sim_rev(P, vector<float>(P, 0.0f));
    vector<vector<float>> calc_rev(P, vector<float>(P, 0.0f));
    for (int p1 = 0; p1 < P; p1++) {
        for (int p2 = 0; p2 < P; p2++) {
            for (int k = 0; k < dumbRuns; k++) {
                int qty = C - 1; 
                for (int t = 0; t < split; t++)
                    if (coin(gen) < alpha[t] && coin(gen) < lambda[t][p1] && qty > 0)
                        { qty--; sim_rev[p1][p2] += prices[p1]; }
                for (int t = split; t < T; t++)
                    if (coin(gen) < alpha[t] && coin(gen) < lambda[t][p2] && qty > 0)
                        { qty--; sim_rev[p1][p2] += prices[p2]; }
            }
        }
    } int c1 = 0, c2 = 0, s1 = 0, s2 = 0; float best_calc_rev = 0.0f, best_sim_rev = 0.0f;
    for (int p1 = 0; p1 < P; p1++) {
        for (int p2 = 0; p2 < P; p2++) {
            float lhs = min(theta1[p1], float(C)) * prices[p1];
            float rhs = min(float(C - min(theta1[p1], float(C))), theta2[p2]) * prices[p2];
            if (sim_rev[p1][p2] > best_sim_rev)  { best_sim_rev = sim_rev[p1][p2]; s1 = p1; s2 = p2; }
            if (lhs + rhs > best_calc_rev) { best_calc_rev = lhs + rhs; c1 = p1; c2 = p2; }
        }
    }
    data << best_calc_rev << "," << prices[c1] << "," << prices[c2] << "," 
    << best_sim_rev / dumbRuns << "," << prices[s1] << "," << prices[s2] << ",";
}

void dumbsim() {
    vector<float> sim_rev(P, 0.0f), calc_rev(P, 0.0f);
    for (int p = 0; p < P; p++) {
        float theta = 0.0f;
        for (int t = 0; t < T; t++) theta += alpha[t] * lambda[t][p];
        calc_rev[p] = prices[p] * min((float)C, theta);
    }
    for (int k = 0; k < dumbRuns; k++) {
        vector<int> qty(P, C - 1);
        for (int t = 0; t < T; t++) {
            for (int p = 0; p < P; p++) {
                if (coin(gen) < alpha[t] && coin(gen) < lambda[t][p] && qty[p] > 0)
                    { qty[p]--; sim_rev[p] += prices[p]; }
            }
        }
    }
    auto sit = max_element(sim_rev.begin(), sim_rev.end());
    auto cit = max_element(calc_rev.begin(), calc_rev.end());
    data << *cit << "," << prices[cit - calc_rev.begin()] << "," 
    << *sit / dumbRuns << "," << prices[sit - sim_rev.begin()] << "\n";
}

int main() {
    data << "seed,capacity,smart_dp,split_calc_rev,split_calc_p1,split_calc_p2,split_sim_rev,"
         << "split_sim_p1,split_sim_p2,dumb_calc_rev,dumb_calc_p,dumb_sim_rev,dum_sim_p\n";
    for (int r = 0; r < mainRuns; r++) {
        C = 0; seed = rd(); gen.seed(seed);
        lambda.assign(T, vector<float>());
        V.assign(T + 1, vector<float>());
        readData(); smartdp(); splitsim(); dumbsim();
    }
    return 0;
}