#include <complex>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <time.h>

using namespace std;

class KMP {
    string pattern;
    vector <int> fail;

    void initKMP(const string &p) {
        pattern = p;
        int m = pattern.size();
        fail.assign(m + 1, -1);
        for (int i = 0, j = -1; i < m; i++) {
            while (j >= 0 && pattern[i] != pattern[j]) j = fail[j];
            j++;
            if (pattern[i + 1] == pattern[j]) fail[i + 1] = fail[j];
            else fail[i + 1] = j;
        }
    }

    void initMP(const string &p) {
        pattern = p;
        int m = pattern.size();
        fail.assign(m + 1, -1);
        for (int i = 0, j = -1; i < m; i++) {
            while (j >= 0 && pattern[i] != pattern[j]) j = fail[j];
            fail[i+1] = ++j;
        }
    }

public:
    KMP(const string &p) { initKMP(p); }

    int period(int i) { return i - fail[i]; }

    vector <int> match(const string &s) {
        int n = s.size();
        int m = pattern.size();
        vector <int> res;
        for (int i = 0, k = 0; i < n; i++) {
            while (k >= 0 && s[i] != pattern[k]) k = fail[k];
            k++;
            if (k == m) res.push_back(i - m + 1);
        }
        return res;
    }
};

vector<complex<double>> fft(vector<complex<double>> a, bool inverse = false) {
    int n = a.size();
    int h = 0;
    for (int i = 0; 1 << i < n; i++) h++;
    for (int i = 0; i < n; i++) {
        int j = 0;
        for (int k = 0; k < h; k++) {
            j |= (i >> k & 1) << (h - 1 - k);
        }
        if (i < j) {
            swap(a[i], a[j]);
        }
    }
    for (int b = 1; b < n; b *= 2) {
        for (int j = 0; j < b; j++) {
            complex<double> w = polar(1.0, (2 * M_PI) / (2 * b) * j * (inverse ? 1 : -1));
            for (int k = 0; k < n; k += b * 2) {
                complex<double> s = a[j + k];
                complex<double> t = a[j + k + b] * w;
                a[j + k] = s + t;
                a[j + k + b] = s - t;
            }
        }
    }
    if (inverse) {
        for (int i = 0; i < n; i++) {
            a[i] /= n;
        }
    }
    return a;
}

vector<complex<double>> fft(vector<double> a, bool inverse = false) {
    vector<complex<double>> a_complex(a.size());
    for (int i = 0; i < a.size(); i++) {
        a_complex[i] = complex<double>(a[i], 0);
    }
    return fft(a_complex, inverse);
}

vector<double> convolve(vector<double> a, vector<double> b) {
    int s = a.size() + b.size() - 1;
    int t = 1;
    while (t < s) t *= 2;
    a.resize(t);
    b.resize(t);
    vector<complex<double>> A = fft(a);
    vector<complex<double>> B = fft(b);
    for (int i = 0; i < t; i++) {
        A[i] *= B[i];
    }
    A = fft(A, true);
    a.resize(s);
    for (int i = 0; i < s; i++) a[i] = A[i].real();
    return a;
}

string generate(int len, int num) {
    string text = "";
    for (int i = 0; i < len; i++) {
        int r = rand() % num;
        text += 'a' + r;
    }
    return text;
}

int query(string text, string pattern, int types) {
    clock_t start = clock();

    std::reverse(pattern.begin(), pattern.end());

    int text_size = text.size();
    int pattern_size = pattern.size();

    vector<vector<double>> cnt(types);
    for (int i = 0; i < types; i++) {
        vector<double> t(text_size, 0.0), p(pattern_size, 0.0);
        for (int j = 0; j < text_size; j++) {
            if (text[j] == 'a' + i) {
                t[j] = 1.0;
            }
        }
        for (int j = 0; j < pattern_size; j++) {
            if (pattern[j] == 'a' + i) {
                p[j] = 1.0;
            }
        }

        cnt[i] = convolve(t, p);
    }

    int ans = 0;
    for (int j = 0; j < cnt[0].size(); j++) {
        double sm = 0.0;
        for (int i = 0; i < types; i++) {
            sm += cnt[i][j];
        }
        if (abs(pattern_size - sm) < 0.5) {
            ans++;
            // cout << "index " << j - pattern_size + 1 << " ~ " << j << " mathced!\n";
        }
    }

    clock_t end = clock();
    std::cout << "fft   : " << (double)(end - start) / CLOCKS_PER_SEC << "sec.\n";

    return ans;
}

int naive(string text, string pattern) {
    clock_t start = clock();

    int ans = 0;
    for (int i = 0; i <= text.size() - pattern.size(); i++) {
        for (int j = 0; j < pattern.size(); j++) {
            if (text[i + j] != pattern[j]) break;
            if (j == pattern.size() - 1) ans++;
        }
    }

    clock_t end = clock();
    std::cout << "naive : " << (double)(end - start) / CLOCKS_PER_SEC << "sec.\n";
    return ans;
}

int kmp(string text, string pattern) {
    clock_t start = clock();

    KMP kmp(pattern);
    vector<int> res = kmp.match(text);

    clock_t end = clock();
    std::cout << "kmp   : " << (double)(end - start) / CLOCKS_PER_SEC << "sec.\n";
    return res.size();
}

int main() {
//    for (int t_len = 1000; t_len <= 100000; t_len *= 10) {
//        for (int p_len = 10; p_len < t_len; p_len *= 10) {
//            string text = generate(t_len, 2);
//            string pattern = generate(p_len, 2);
//            cout << "---------------------------------" << endl;
//            cout << "n=" << t_len << ", m=" << p_len << endl;
//            int ans1 = naive(text, pattern);
//            int ans2 = query(text, pattern, 2);
//            int ans3 = kmp(text, pattern);
//            cout << ans1 << " matched!!" << endl;
//        }
//    }

//    for (int t_len = 1000; t_len <= 100000; t_len *= 10) {
//        for (int p_len = 10; p_len < t_len; p_len *= 10) {
//            string tmp = "ab";
//            string text = "";
//            for (int k = 0; k < t_len / 2; k++) text += tmp;
//            string pattern = "";
//            for (int k = 0; k < p_len / 2; k++) pattern += tmp;
//
//            cout << "---------------------------------" << endl;
//            cout << "n=" << text.size() << ", m=" << pattern.size() << endl;
//            int ans1 = naive(text, pattern);
//            int ans2 = query(text, pattern, 2);
//            int ans3 = kmp(text, pattern);
//            cout << ans1 << " matched!!" << endl;
//        }
//    }

    for (int t_len = 1000; t_len <= 100000; t_len *= 10) {
        for (int p_len = 10; p_len < t_len; p_len *= 10) {
            string tmp = "abcde";
            string text = "";
            for (int k = 0; k < t_len / 5; k++) text += tmp;
            string pattern = "";
            for (int k = 0; k < p_len / 5; k++) pattern += tmp;

            cout << "---------------------------------" << endl;
            cout << "n=" << text.size() << ", m=" << pattern.size() << endl;
            int ans1 = naive(text, pattern);
            int ans2 = query(text, pattern, 5);
            int ans3 = kmp(text, pattern);
            cout << ans1 << " matched!!" << endl;
        }
    }

    return 0;
}