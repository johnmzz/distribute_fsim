#include "InitializerMPI.h"

InitializerMPI::InitializerMPI(string method, const vector<string>& label_info_1, const vector<string>& label_info_2) {
    if (method == "IdenticalInitializer" || method == "identical") {
        method_ = "identical";
    } else if (method == "JaccardInitializer" || method == "jaccard") {
        method_ = "jaccard";
    } else if (method == "HammingInitializer" || method == "hamming") {
        method_ = "hamming";
    } else if (method == "EditDistInitializer" || method == "edit") {
        method_ = "edit";
    } else if (method == "JarowinklerInitializer" || method == "jarowinkler") {
        method_ = "jarowinkler";
    } else {
        cout << "Initializer Type Error, " << method << endl;
        exit(0);
    }

    label_info_1_ = vector<string>(label_info_1.size(), "");
    label_info_2_ = vector<string>(label_info_2.size(), "");

    for (uint32_t i = 0; i < label_info_1.size(); i++) {
        label_info_1_[i] = label_info_1[i];
    }

    for (uint32_t i = 0; i < label_info_2.size(); i++) {
        label_info_2_[i] = label_info_2[i];
    }

    VectorFloatType temp;
    for (uint32_t i = 0; i < label_info_1.size(); i++) {
        temp = VectorFloatType(label_info_2.size(), 0.0);
        label_sim_matrix_.push_back(temp);
    }

    for (uint32_t i = 0; i < label_sim_matrix_.size(); i++) {
        for (uint32_t j = 0; j < label_sim_matrix_[0].size(); j++) {
            label_sim_matrix_[i][j] = cal_string_sim(label_info_1[i], label_info_2[j]);
        }
    }
}

int InitializerMPI::my_min(int x, int y, int z) {
    return std::min(std::min(x, y), z);
}

float InitializerMPI::cal_string_sim(const string& s1, const string& s2) {
    if (method_ == "identical") {
        float result = 0.0;
        if (s1 == s2) {
            result = 1.0;
        } else {
            result = 0.0;
        }
        return result;
    }
    // jaccard
    else if (method_ == "jaccard") {
        vector<int> freq1 = vector<int>(256, 0);
        vector<int> freq2 = vector<int>(256, 0);

        for (char s : s1) {
            if (s == '\"' || (int)s >= 256 || (int)s < 0) {
                continue;
            }
            (s >= 'A' && s <= 'Z') ? freq1[s + 32]++ : freq1[s]++;
        }
        for (char s : s2) {
            if (s == '\"' || (int)s >= 256 || (int)s < 0) {
                continue;
            }
            (s >= 'A' && s <= 'Z') ? freq2[s + 32]++ : freq2[s]++;
        }

        float common = 0.0, total = 0.0;
        for (int i = 0; i < 256; i++) {
            common += (freq1[i] >= freq2[i] ? freq2[i] : freq1[i]);
            total += (freq1[i] >= freq2[i] ? freq1[i] : freq2[i]);
        }
        if (total == 0) {
            cout << s1 << ", " << s2 << " jaccard similarity error" << endl;
            exit(0);
        }
        return (float)common / total;
    }
    // hamming
    else if (method_ == "hamming") {
        if (s1.length() == 0 || s2.length() == 0) {
            return 0;
        }

        if (s1 == "-" && s2 != "-") {
            return 0;
        }
        int min, max, dist = 0;
        float sim = 0.0;
        if (s1.length() >= s2.length()) {
            max = s1.length();
            min = s2.length();
        } else {
            max = s2.length();
            min = s1.length();
        }

        for (int i = 0; i < min; i++) {
            if (s1[i] != s2[i]) {
                dist += 1;
            }
        }

        dist = dist + max - min;
        sim = 1 - (float)dist / max;
        return sim;
    }
    // edit
    else if (method_ == "edit") {
        if (s1.length() == 0 || s2.length() == 0) {
            return 0;
        }

        if (s1 == "-" && s2 != "-") {
            return 0;
        }

        int m = s1.length();
        int n = s2.length();
        // Create a table to store results of subproblems
        int dp[m + 1][n + 1];

        // Fill d[][] in bottom up manner
        for (int i = 0; i <= m; i++) {
            for (int j = 0; j <= n; j++) {
                // If first string is empty, only option is to insert all characters of second string
                if (i == 0) {
                    dp[i][j] = j;  // Min. operations = j
                }
                // If second string is empty, only option is to remove all characters of second string
                else if (j == 0) {
                    dp[i][j] = i;  // Min. operations = i
                }
                // If last characters are same, ignore last char and recur for remaining string
                else if (s1[i - 1] == s2[j - 1]) {
                    dp[i][j] = dp[i - 1][j - 1];
                }
                // If the last character is different, consider all possibilities and find the minimum
                else {
                    dp[i][j] = 1 + my_min(dp[i][j - 1], dp[i - 1][j], dp[i - 1][j - 1]);
                }
            }
        }

        float sim = 1 - (float)dp[m][n] / std::max(m, n);
        return sim;
        // return (float) dp[m][n];
    }
    // jarowinkler
    else {
        float m = 0;
        int low, high, range;
        int k = 0, numTrans = 0;

        // Exit early if either are empty
        if (s1.length() == 0 || s2.length() == 0) {
            return 0;
        }
        if (s1 == "-" && s2 != "-") {
            return 0;
        }
        // Exit early if they're an exact match.
        if (s1 == s2) {
            return 1;
        }

        range = (std::max(s1.length(), s2.length()) / 2) - 1;
        int s1Matches[s1.length()] = {};
        int s2Matches[s2.length()] = {};

        for (int i = 0; i < s1.length(); i++) {
            // Low Value;
            if (i >= range) {
                low = i - range;
            } else {
                low = 0;
            }
            // High Value;
            if (i + range <= (s2.length() - 1)) {
                high = i + range;
            } else {
                high = s2.length() - 1;
            }

            for (int j = low; j <= high; j++) {
                if (s1Matches[i] != 1 && s2Matches[j] != 1 && s1[i] == s2[j]) {
                    m += 1;
                    s1Matches[i] = 1;
                    s2Matches[j] = 1;
                    break;
                }
            }
        }

        // Exit early if no matches were found
        if (m == 0) {
            return 0;
        }

        // Count the transpositions.
        for (int i = 0; i < s1.length(); i++) {
            if (s1Matches[i] == 1) {
                int j;
                for (j = k; j < s2.length(); j++) {
                    if (s2Matches[j] == 1) {
                        k = j + 1;
                        break;
                    }
                }

                if (s1[i] != s2[j]) {
                    numTrans += 1;
                }
            }
        }

        float weight = (m / s1.length() + m / s2.length() + (m - (numTrans / 2)) / m) / 3;
        float l = 0;
        float p = 0.1;
        if (weight > 0.7) {
            while (s1[l] == s2[l] && l < 4) {
                l += 1;
            }
            weight += l * p * (1 - weight);
        }
        return weight;
    }
}

float InitializerMPI::get_sim_value(uint32_t i, uint32_t j) {
    return label_sim_matrix_[i][j];
}

void InitializerMPI::print_label_info() {
    cout << "Graph 1 labels: ";
    for (uint32_t i = 0; i < label_info_1_.size(); i++) {
        cout << label_info_1_[i] << " ";
    }
    cout << endl;

    cout << "Graph 2 labels: ";
    for (uint32_t i = 0; i < label_info_2_.size(); i++) {
        cout << label_info_2_[i] << " ";
    }
    cout << endl;
}

void InitializerMPI::print_label_sim() {
    for (uint32_t i = 0; i < label_info_1_.size(); i++) {
        for (uint32_t j = 0; j < label_info_2_.size(); j++) {
            cout << setw(7) << setprecision(3) << label_sim_matrix_[i][j] << " ";
        }
        cout << endl;
    }
}