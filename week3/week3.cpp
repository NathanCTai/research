#include <vector>
#include <random>
#include <fstream>
#include <iostream>
using namespace std;

int T = 100, P = 21, C = 0, dumbRuns = 1000, mainRuns = 3; float L = 0.0f; 
vector<vector<float>> lambda(T), V(T + 1); // λ_t(p) is T x P; 
vector<float> alpha(T), gamma(T), prices(P); // V_t(x) is (T + 1) x C;
ofstream data("data.csv"); random_device rd; unsigned int seed = rd(); mt19937 gen(seed); 
uniform_real_distribution<float> alpha_dist(0.5f, 1.0f), gamma_dist(0.25f, 0.4f), coin(0.0f, 1.0f);

void readData() {
    for (int i = 0; i < T; i++) {
        float gap = 0;
        if (i < 50) {  gap = 0.15; } else { gap = 0;}
        alpha[i] = alpha_dist(gen); gamma[i] = gamma_dist(gen) - gap; lambda[i].resize(P);
        for (int j = 0; j < P; j++) {
            if (i == 0) {prices[j] = float(j);} // reading in prices
            lambda[i][j] = exp(-gamma[i] * prices[j]); 
        }
        L += alpha[i] * exp(-1); 
    } 
    C = int(L * 0.8); data << seed << "," << C << ",";
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
    const int split = 20;
    vector<vector<float>> sim_rev(P, vector<float>(P, 0.0f));
    vector<vector<float>> calc_rev(P, vector<float>(P, 0.0f));
    for (int i = 0; i < P; i++) {
        for (int j = 0; j < P; j++) {
            int p1 = prices[i], p2 = prices[j];
            float theta1 = 0.0f, theta2 = 0.0f;
            for (int t = 0; t < split; t++)  theta1 += alpha[t] * lambda[t][p1];
            for (int t = split; t < T; t++)  theta2 += alpha[t] * lambda[t][p2];
            float exp1 = min((float)C, theta1);
            calc_rev[i][j] = p1 * exp1 + p2 * min((float)C - exp1, theta2);
            for (int k = 0; k < dumbRuns; k++) {
                int qty = C - 1;
                for (int t = 0; t < split; t++)
                    if (coin(gen) < alpha[t] && coin(gen) < lambda[t][p1] && qty > 0)
                        { qty--; sim_rev[i][j] += p1; }
                for (int t = split; t < T; t++)
                    if (coin(gen) < alpha[t] && coin(gen) < lambda[t][p2] && qty > 0)
                        { qty--; sim_rev[i][j] += p2; }
            }
        }
    } int c1 = 0, c2 = 0, s1 = 0, s2 = 0; float best_calc_rev = 0.0f, best_sim_rev = 0.0f;
    for (int i = 0; i < P; i++) {
        for (int j = 0; j < P; j++) {
            if (sim_rev[i][j]  > best_sim_rev)  { best_sim_rev  = sim_rev[i][j];  s1 = i; s2 = j; }
            if (calc_rev[i][j] > best_calc_rev) { best_calc_rev = calc_rev[i][j]; c1 = i; c2 = j; }
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
        L = 0.0f; C = 0; seed = rd(); gen.seed(seed);
        fill(alpha.begin(), alpha.end(), 0.0f);
        fill(gamma.begin(), gamma.end(), 0.0f);
        lambda.assign(T, vector<float>());
        V.assign(T + 1, vector<float>());
        readData(); smartdp(); splitsim(); dumbsim();
    }
    return 0;
}