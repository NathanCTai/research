#include <vector>
#include <random>
#include <fstream>
#include <iostream>
using namespace std;

const int T = 100, P = 41, Split1 = 20, Split2 = 40, DumbRuns = 1000, MainRuns = 10; int C = 0; 
vector<vector<float>> lambda(T); // λ_t(p) is T x P; read in theta from the start
vector<float> alpha(T), gamma(T), prices(P), theta1(P), theta2(P);
ofstream data("data.csv"); random_device rd; unsigned int seed = rd(); mt19937 gen(seed); 
uniform_real_distribution<float> alpha_dist(0.5f, 1.0f), gamma_dist(0.1f, 0.4f), coin(0.0f, 1.0f);

void readData() {
    fill(theta1.begin(), theta1.end(), 0.0f); fill(alpha.begin(), alpha.end(), 0.0f); 
    fill(theta2.begin(), theta2.end(), 0.0f); fill(gamma.begin(), gamma.end(), 0.0f);
    for (int p = 0; p < P; p++) { prices[p] = float(p) * 0.5;
        for (int t = 0; t < T; t++) {
            if (p == 0) {alpha[t] = alpha_dist(gen); gamma[t] = gamma_dist(gen); lambda[t].resize(P); }
            lambda[t][p] = exp(-gamma[t] * prices[p]); 
            if (t <  Split1) { theta1[p] += alpha[t] * lambda[t][p]; }
            if (t >= Split1) { theta2[p] += alpha[t] * lambda[t][p]; } // new way of calculating theta
        }
    } C = int(accumulate(alpha.begin(), alpha.end(), 0.0f) * 0.29); data << seed << "," << C << ",";
}

void smartdp(int start, int end, int cap) {
    vector<vector<float>> V = vector<vector<float>>(T + 1, vector<float>(cap, 0.0f));
    for (int t = start; t >= end; t--) { 
        for (int x = 1; x < cap; x++) {
            float maximum = 0.0; float optimalprice = 0.0; 
            for (int p = 0; p < P; p++) { 
                float lhs = alpha[t] * lambda[t][p] * (prices[p] + V[t + 1][x - 1]); // SALE
                float rhs = (1 - alpha[t] * lambda[t][p]) * V[t + 1][x]; // NO SALE
                if (lhs + rhs > maximum) {
                    maximum = rhs + lhs; optimalprice = prices[p];
                } 
            } V[t][x] = maximum;
        }
    } data << V[0][cap - 1] << ",";
}

void splitsim(int start, int end, int split, int cap) {
    vector<vector<float>> sim_rev(P, vector<float>(P, 0.0f));
    for (int p1 = 0; p1 < P; p1++) {
        for (int p2 = 0; p2 < P; p2++) {
            for (int k = 0; k < DumbRuns; k++) {
                int qty = cap - 1; 
                for (int t = start; t < split; t++)
                    if (coin(gen) < alpha[t] && coin(gen) < lambda[t][p1] && qty > 0)
                        { qty--; sim_rev[p1][p2] += prices[p1]; }
                for (int t = split; t < end; t++)
                    if (coin(gen) < alpha[t] && coin(gen) < lambda[t][p2] && qty > 0)
                        { qty--; sim_rev[p1][p2] += prices[p2]; }
            }
        }
    } int c1 = 0, c2 = 0, s1 = 0, s2 = 0; float best_calc_rev = 0.0f, best_sim_rev = 0.0f;
    for (int p1 = 0; p1 < P; p1++) {
        for (int p2 = 0; p2 < P; p2++) {
            float lhs = min(theta1[p1], float(cap)) * prices[p1];
            float rhs = min(float(cap - min(theta1[p1], float(cap))), theta2[p2]) * prices[p2];
            if (sim_rev[p1][p2] > best_sim_rev)  { best_sim_rev = sim_rev[p1][p2]; s1 = p1; s2 = p2; }
            if (lhs + rhs > best_calc_rev) { best_calc_rev = lhs + rhs; c1 = p1; c2 = p2; }
        }
    }
    data << best_calc_rev << "," << prices[c1] << "," << prices[c2] << "," 
    << best_sim_rev / DumbRuns << "," << prices[s1] << "," << prices[s2] << ",";
}

void dumbsim(int start, int end, int cap) {
    vector<float> sim_rev(P, 0.0f), calc_rev(P, 0.0f);
    for (int p = 0; p < P; p++) {
        float theta = 0.0f;
        for (int t = start; t < end; t++) theta += alpha[t] * lambda[t][p];
        calc_rev[p] = prices[p] * min((float)C, theta);
    }
    for (int k = 0; k < DumbRuns; k++) {
        vector<int> qty(P, cap - 1);
        for (int t = start; t < end; t++) {
            for (int p = 0; p < P; p++) {
                if (coin(gen) < alpha[t] && coin(gen) < lambda[t][p] && qty[p] > 0)
                    { qty[p]--; sim_rev[p] += prices[p]; }
            }
        }
    }
    auto sit = max_element(sim_rev.begin(), sim_rev.end());
    auto cit = max_element(calc_rev.begin(), calc_rev.end());
    data << *cit << "," << prices[cit - calc_rev.begin()] << "," 
    << *sit / DumbRuns << "," << prices[sit - sim_rev.begin()] << "\n";
}

int main() {
    data << "seed,capacity,smart_dp,split_calc_rev,split_calc_p1,split_calc_p2,split_sim_rev,"
         << "split_sim_p1,split_sim_p2,dumb_calc_rev,dumb_calc_p,dumb_sim_rev,dum_sim_p\n";
    for (int r = 0; r < MainRuns; r++) {
        C = 0; seed = rd(); gen.seed(seed);
        lambda.assign(T, vector<float>());
        readData(); smartdp(0, 100, C); splitsim(0, 20, 100, C); dumbsim(0, 100, C);
    }
    return 0;
}