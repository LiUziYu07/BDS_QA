#include <iostream>
#include <vector>
using namespace std;

vector<double> rungeKutta(vector<double>& vec0, int& speed1, int& speed2, int& speed3, double& h, int& epoch);
double calEquation1(vector<double>& vec, int& speed1, int& speed2, int& speed3, double& h);
double calEquation2(vector<double>& vec, int& speed1, int& speed2, double& h);
double calEquation3(vector<double>& vec, int& speed1, int& speed2, int& speed3, double& h);
double calEquation4(vector<double>& vec, int& speed3, double& h);

int main() {
    vector<double> vec0 = {1, 10, 0, 0};
    int speed1 = 100;
    int speed2 = 600;
    int speed3 = 150;
    double h = 0.00001;
    int epoch = 100000;
    vector<double> ans = rungeKutta(vec0, speed1, speed2, speed3, h, epoch);
    cout << "E: " << ans[0] << "\tS: " << ans[1] << "\tES: " << ans[2] << "\tP: " << ans[3] << endl;
    return 0;
}

/**
 * @brief 
 * Implement Fourth Order RungeKutta Algorithm
 * @param vec0 the inital value of four substances
 * @param speed1 rate k1
 * @param speed2 rate k2
 * @param speed3 rete k3
 * @param h step length of RungeKutta
 * @param epoch iteration times of RungeKutta
 * @return vector<double> final result of four substances
 */
vector<double> rungeKutta(vector<double>& vec0, int& speed1, int& speed2, int& speed3, double& h, int& epoch){
    vector<double> ans(vec0);
    for(int i = 0;i < epoch;i++){
        ans[0] = calEquation1(ans, speed1, speed2, speed3, h);
        ans[1] = calEquation2(ans, speed1, speed2, h);
        ans[2] = calEquation3(ans, speed1, speed2, speed3, h);
        ans[3] = calEquation4(ans, speed3, h);
    }
    return ans;
}

/**
 * @brief 
 * Sovle the second equation about substance E
 * @param vec the value of four substances
 * @param speed1 rate k1
 * @param speed2 rate k2
 * @param speed3 rate k3
 * @param h step length of RungeKutta
 * @return double the value of E(n+1)
 */
double calEquation1(vector<double>& vec, int& speed1, int& speed2, int& speed3, double& h){
    double k1 = (speed2 + speed3) * vec[2] - speed1 * vec[0] * vec[1];
    double k2 = (speed2 + speed3) * vec[2] - speed1 * (vec[0] + static_cast<double>(h / 2) * k1) * vec[1];
    double k3 = (speed2 + speed3) * vec[2] - speed1 * (vec[0] + static_cast<double>(h / 2) * k2) * vec[1];
    double k4 = (speed2 + speed3) * vec[2] - speed1 * (vec[0] + h * k3) * vec[1];
    return static_cast<double>(vec[0] + static_cast<double>(h / 6) * (k1 + 2 * k2 + 2 * k3 + k4));
}

/**
 * @brief 
 * Sovle the second equation about substance S
 * @param vec the value of four substances
 * @param speed1 rate k1
 * @param speed2 rate k2
 * @param h step length of RungeKutta
 * @return double the value of S(n+1)
 */
double calEquation2(vector<double>& vec, int& speed1, int& speed2, double& h){
    double k1 = speed2 * vec[2] - speed1 * vec[1] * vec[0];
    double k2 = speed2 * vec[2] - speed1 * (vec[1] + static_cast<double>(h / 2) * k1) * vec[0];
    double k3 = speed2 * vec[2] - speed1 * (vec[1] + static_cast<double>(h / 2) * k2) * vec[0];
    double k4 = speed2 * vec[2] - speed1 * (vec[1] + h * k3) * vec[0];
    return static_cast<double>(vec[1] + static_cast<double>(h / 6) * (k1 + 2 * k2 + 2 * k3 + k4));
}

/**
 * @brief 
 * Solve the third equation about substance ES
 * @param vec the value of four substances
 * @param speed1 rate k1
 * @param speed2 rate k2
 * @param speed3 rate k3
 * @param h step length of RungeKutta
 * @return double the value of ES(n+1)
 */
double calEquation3(vector<double>& vec, int& speed1, int& speed2, int& speed3, double& h){
    double k1 = speed1 * vec[0] * vec[1] - (speed2 + speed3) * vec[2];
    double k2 = speed1 * vec[0] * vec[1] - (speed2 + speed3) * (vec[2] + static_cast<double>(h / 2) * k1);
    double k3 = speed1 * vec[0] * vec[1] - (speed2 + speed3) * (vec[2] + static_cast<double>(h / 2) * k2);
    double k4 = speed1 * vec[0] * vec[1] - (speed2 + speed3) * (vec[2] + h * k3);
    return static_cast<double>(vec[2] + static_cast<double>(h / 6) * (k1 + 2 * k2 + 2 * k3 + k4));
}

/**
 * @brief 
 * Solve the third equation about substance P
 * @param vec the value of four substances
 * @param speed3 rate k3
 * @param h step length of RungeKutta
 * @return double the value of P(n+1)
 */
double calEquation4(vector<double>& vec, int& speed3, double& h){
    double k1 = speed3 * vec[2];
    double k2 = speed3 * vec[2];
    double k3 = speed3 * vec[2];
    double k4 = speed3 * vec[2];

    return static_cast<double>(vec[3] + static_cast<double>(h / 6) * (k1 + 2 * k2 + 2 * k3 + k4));
}
