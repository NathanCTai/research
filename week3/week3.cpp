#include <vector>
#include <random>
#include <fstream>
#include <iostream>
using namespace std;

const int T = 100, P = 21, split = 20, dumbRuns = 1000, mainRuns = 1; int C = 0; 
vector<vector<float>> lambda(T), V(T + 1); // λ_t(p) is T x P; read in theta from the start
vector<float> alpha(T, 0.0f), gamma(T, 0.0f), prices(P, 0.0f), theta1(P, 0.0f), theta2(P, 0.0f);
ofstream data("data.csv"), theta_log("log.csv"); random_device rd; unsigned int seed = rd(); mt19937 gen(seed); 
uniform_real_distribution<float> alpha_dist(0.5f, 1.0f), gamma_dist(0.1f, 0.4f), coin(0.0f, 1.0f);

void readData() {
    for (int p = 0; p < P; p++) {
        prices[p] = float(p);
        for (int t = 0; t < T; t++) {
            if (p == 0) {alpha[t] = alpha_dist(gen); gamma[t] = gamma_dist(gen); lambda[t].resize(P); }
            lambda[t][p] = exp(-gamma[t] * prices[p]); 
            if (t <  split) { theta1[p] += alpha[t] * lambda[t][p]; }
            if (t >= split) { theta2[p] += alpha[t] * lambda[t][p]; }
        }
    }
    for (int p = 0; p < P; p++) {
        printf("p : %d, theta1 : %f, theta2 : %f \n", p, theta1[p], theta2[p]);
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
            theta_log << p1 << "," << p2 << "," << theta1 << "," << theta2 << "," << log(theta1) << "," << log(theta2) << "\n";
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
    theta_log << "p1,p2,theta1,theta2,log_t1,log_t2\n";
    for (int r = 0; r < mainRuns; r++) {
        C = 0; seed = rd(); gen.seed(seed);
        lambda.assign(T, vector<float>());
        V.assign(T + 1, vector<float>());
        readData(); smartdp(); splitsim(); dumbsim();
    }
    return 0;
}