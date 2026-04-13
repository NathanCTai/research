#include <cmath>
#include <vector>
#include <string>
#include <random>
#include <fstream>
#include <iostream>
using namespace std;

const int T = 100, P = 41, Runs = 1000, Sims = 20; int C, Seed = 0; float K = 5.0;
vector<vector<int>> TFPreCompP1DX(T), DPsugIDX(T);
vector<vector<float>> Lambda(T);
vector<float> Alpha(T, 0.0f), Gamma(T), Prices(P);
mt19937 Gen(Seed); vector<string> Categories = {"Uniform", "Linear", "Exp_Gentle", "Exp_Steep", "Concave"};
uniform_real_distribution<float> AlphaDist(0.5f, 1.0f), GammaDist(-0.1f, 0.1f), Coin(0.0f, 1.0f);

void resetData() {
    fill(Alpha.begin(), Alpha.end(), 0.0f);
    fill(Gamma.begin(), Gamma.end(), 0.0f);
    Lambda.assign(T, vector<float>(P, 0.0f));
    TFPreCompP1DX.assign(T, vector<int>(P, 0));
}

void readData(int mode, float mult) {
    for (float t = 0.0; t < T; t = t + 1.0) {
        Alpha[t] = AlphaDist(Gen); int g = T - 1;
        switch(mode) {
            case 0: Gamma[t] = 0.25 + GammaDist(Gen); break; // UNIFORM BASELINE
            case 1: Gamma[t] = 0.1 + 0.3 * t / g + GammaDist(Gen); break; // LINEAR
            case 2: Gamma[t] = 0.1 * exp(log(4) * t / g) + GammaDist(Gen); break; // EXP GENTLE
            case 3: Gamma[t] = 0.1 + 0.3 * (exp(K * t / g - 1) / (exp(K) - 1)) + GammaDist(Gen); break; // EXP STEEP
            case 4: Gamma[t] = 0.1 + 0.3 * (1 - exp(-K * t / g)) / (1 - exp(-K)) + GammaDist(Gen); break; // CONCAVE
        }
        for (int p = 0; p < P; p++) {
            if (t == 0) { Prices[p] = float(p) / 2.0; }
            Lambda[t][p] = exp(-Gamma[t] * Prices[p]);
        }
    } C = int(accumulate(Alpha.begin(), Alpha.end(), 0.0f) * mult * exp(-1));
}

float calcDP() {
    vector<vector<float>> DP(T + 1, vector<float>(C + 1, 0.0f));
    for (int t = T - 1; t > -1; t--) { 
        DPsugIDX[t].assign(C + 1, 0); // DP suggests price INDICES to charge at 
        for (int x = 1; x < C + 1; x++) {
            float maximum = 0.0;
            for (int p = 0; p < P; p++) { 
                float lhs = Alpha[t] * Lambda[t][p] * (Prices[p] + DP[t + 1][x - 1]);
                float rhs = (1 - Alpha[t] * Lambda[t][p]) * DP[t + 1][x];
                if (lhs + rhs > maximum) {
                    maximum = rhs + lhs; DPsugIDX[t][x] = p;
                } 
            } DP[t][x] = maximum;
        }
    } return DP[0][C];
}

float simDP() {
    float rev = 0.0;
    for (int run = 0; run < Runs; run++) {
        int stock = C; 
        for (int t = 0; t < T; t++) {
            int priceIDX = DPsugIDX[t][stock];
            if (Coin(Gen) < Alpha[t] && Coin(Gen) < Lambda[t][priceIDX] && stock > 0) {
                stock -= 1; rev += Prices[priceIDX];
            }
        }
    }
    return rev / Runs;
}

float staticOneFareSim() {
    vector<float> rev(P, 0.0f);
    for (int run = 0; run < Runs; run++) {
        vector<int> stock(P, C - 1);
        for (int t = 0; t < T; t++) {
            for (int p = 0; p < P; p++) {
                if (Coin(Gen) < Alpha[t] && Coin(Gen) < Lambda[t][p] && stock[p] > 0) {
                    stock[p] -= 1; rev[p] += Prices[p];
                }
            }
        }
    }
    return *max_element(rev.begin(), rev.end()) / Runs;
}

void fillTFrecalc(int table) { // Stores the INDICES of the optimal P1 at every time/capacity combination
    for (int t = 0; t < T; t++) {
        TFPreCompP1DX[t].assign(C, 0); int split = 0;
        for (int c = 0; c < C; c++) {
            if (table == 0) { split = T * 0.2 + t * 0.8; }
            else { split = T * 0.5 + t * 0.5; }
            float best_rev = 0; int best_p1_idx = 0;
            for (int p1 = 0; p1 < P; p1++) {
                for (int p2 = 0; p2 < P; p2++) {
                    float theta1 = 0.0f, theta2 = 0.0f;
                    for (int u = t; u < split; u++) theta1 += Alpha[u] * Lambda[u][p1];
                    for (int u = split; u < T; u++) theta2 += Alpha[u] * Lambda[u][p2];
                    float exp1 = min(float(c), theta1);
                    float calc_rev = Prices[p1] * exp1 + Prices[p2] * min(float(c) - exp1, theta2);
                    if (calc_rev > best_rev) { best_p1_idx = p1; best_rev = calc_rev; }
                }
            }
            TFPreCompP1DX[t][c] = best_p1_idx; 
        }
    }
}

float recalcTF(int gap) {
    float rev = 0.0f;
    for (int run = 0; run < Runs; run++) {
        int stock = C;
        for (int s = 0; s < T; s = s + gap) {
            for (int t = s; t < s + gap; t++) {
                int p = TFPreCompP1DX[s][stock];
                if (Coin(Gen) < Alpha[t] && Coin(Gen) < Lambda[t][p] && stock > 0) {
                    stock -= 1; rev += Prices[p];
                }
            }
        }
    }
    return rev / Runs;
}

int main() {
    ofstream data("data.csv");
    data << "mode,mult,";
    data << "TF_Recalc(10:20-80),TF_Recalc(5:20-80),TF_Recalc(10:50-50),TF_Recalc(5:50-50),";
    data << "DP_Calc,DP_Sim,One_Fare_Static_Sim,Seed\n";
    for (int sim = 1; sim <= Sims; sim++) {
        Seed = sim; Gen.seed(Seed);
        for (int mode = 0; mode < 5; mode++) {
            for (float mult = 0.6; mult < 0.9; mult = mult + 0.1f) {
                data << Categories[mode] << "," << mult << ",";
                resetData();
                readData(mode, mult);
                fillTFrecalc(0); // 20-80 split
                data << recalcTF(10) << "," << recalcTF(5) << ",";
                fillTFrecalc(1); // 50-50 split
                data << recalcTF(10) << "," << recalcTF(5) << ",";
                data << calcDP() << "," << simDP() << "," << staticOneFareSim() << "," << Seed << "\n";
            }
        }
    }
    return 0;
}