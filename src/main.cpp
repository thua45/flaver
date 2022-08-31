//
// Created by Tinghua Huang on 19/12/30.
//

#include <getopt.h>
#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <algorithm>
#include "cdflib.hpp"

using namespace std;

struct stat {
    string tf_name;
    float pp;
    double tau;
    double estimate;
    double pvalue;
    double fdr;
    bool operator > (const stat & r) const
    { return pvalue > r.pvalue; }
    bool operator < (const stat & r) const
    { return pvalue < r.pvalue; }
};

//double weigher(int x, int n, float ratio);
//double weigher(int x);
//double weigher(int x, float pp);
double weigher(int x, float pp, int n);
//long double v_sum(int x_n);
//long double vsq_sum(int x_n);
long double v_sum(int x_n, float pp);
long double vsq_sum(int x_n, float pp);
//long double v_mean(int x_n,  float ratio);
//long double vsq_mean(int x_n,  float ratio);
int sgn(int x);
//void w_ranked_tau(vector<int> vec_x, vector<int> vec_y, float ratio);
void bubble_asort(vector<int> &nums, vector<int> &order);
void bubble_asort(vector<double> &nums, vector<int> &order);
void bubble_lexsort(vector<int> &vec_x, vector<int> &vec_y, vector<int> &order);
void vec_rank(vector<int> &vec, vector<int> &rank);
void vec_rank(vector<double> &vec, vector<int> &rank);
//void read_fpp(string fpp_file);
void trim_string(string &src);
vector<string> str_split(string& src, string delimit);
double tau_pvalue(double estimate);
//double w_ranked_tau(vector<int> &vec_x, vector<int> &vec_y, float pp);
double w_ranked_tau(vector<int> &rank_x, vector<int> &rank_y, vector<int> &r_order, float pp);
void read_grit_dataset(string grit_file);
void read_rank_data(string rank_file);
void overlap_set(map<string, double> &map_x, map<string, double> &map_y, vector<double> &vec_x, vector<double> &vec_y);
void one_tau_analysis(vector<double> &vec_x, vector<double> &vec_y, float &best_pp, double &best_tau, double &estimate, double &best_pvalue);
void performe_wrktf_analysis();
void calculate_fpp(int x_n);
bool stat_less_p(stat a, stat b);
bool stat_less_fdr(stat a, stat b);
void mcf_die(const std::string & message);
void save_result(string output_file);

char const short_options[] = "hm:i:b:n:z:s:t:p:o:d:";
struct option long_options[] =
        {
            {"help", 0, NULL, 'h'},
            {"rank_set", 1, NULL, 's'},
            {"rank_data", 1, NULL, 'i'},
            {"q_score", 0, NULL, 'q'},
            {"output", 1, NULL, 'o'},
            {"dir", 1, NULL, 'd'},
            {0, 0, 0, 0}
        };

map<float, float> fpp_mp;
vector<map<string, double>> grit_set;
vector<string> tf_names;
map<string, double> rank_data;
vector<stat> result_stats;

int main(int argc, char *argv[]) {

    if (argc == 1)
    {
        printf("flaver -s rank_set -i rank_data [-q q_score] -o output\n");
        exit(1);

    }

    string rset_file = "rank_set.txt";
    string rdata_file = "rank_data.txt";
    string output_file = "output.txt";
    double q_score = 0.05;
    string output_dir = "./";

    int c;
    while((c=getopt_long(argc, argv, short_options, long_options, NULL)) != -1)
    {
        switch (c)
        {
            case 'h':
                printf("flaver -s rank_set -i rank_data [-q q_score] -o output\n");
                exit(0);
                break;
            case 's':
                rset_file = optarg;
                break;
            case 'i':
                rdata_file = optarg;
                break;
            case 'q':
                //cout<<optarg<<endl;
                q_score = atof(optarg);
                break;
            case 'o':
                output_file = optarg;
                break;
            default :
                cout << "invalid option: " << optarg << endl;
                exit(1);
        }
    } 

    /*
    int arry_x[5] = {12, 2, 1, 12, 2};
    vector<int> vec_x(arry_x, arry_x + 5);
    int arry_y[5] = {2, 4, 7, 1, 0};
    vector<int> vec_y(arry_y, arry_y + 5);
    */
    /*
    int arry_x[8] = {1,2,3,4,5,6,7,8};
    vector<int> vec_x(arry_x, arry_x + 8);
    int arry_y[8] = {3,4,1,2,5,7,8,6};
    vector<int> vec_y(arry_y, arry_y + 8);

    w_ranked_tau(vec_x, vec_y, 1.0);
    */

    /*
    int arry_x2[8] = {1,2,3,4,5,6,7,8};
    vector<int> vec_x2(arry_x2, arry_x2 + 8);
    int arry_y2[8] = {1,2,3,4,6,8,7,5};
    vector<int> vec_y2(arry_y2, arry_y2 + 8);
    */

    /*
    w_ranked_tau(vec_x2, vec_y2, 1.0);
    w_ranked_tau(vec_x2, vec_y2, 0.7);
    w_ranked_tau(vec_x2, vec_y2, 0.5);
    w_ranked_tau(vec_x2, vec_y2, 0.3);
    */

    /*
    int arry_x2[20] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};
    vector<int> vec_x2(arry_x2, arry_x2 + 20);
    int arry_y2[20] = {1,2,4,3,6,5,20,13,15,19,17,18,14,16,12,9,11,8,10,7};
    vector<int> vec_y2(arry_y2, arry_y2 + 20);

    for (int i = 1; i <= 10; i++) {
        float ratio = float(i) / 10;
        w_ranked_tau(vec_x2, vec_y2, ratio);
        cout << endl;
    }
    */

    /*
    float tau = 1.0;
    long int x_n = 1;
    float estimate = 1.5 * tau * v_sum(x_n) / sqrt(vsq_sum(x_n));
    cout << x_n << "\t" << estimate << endl;

    cout << endl;

    for (int i = 1; i <= 100; i++) {
        long int x_n = i * 100000000;
        float estimate = 1.5 * tau * v_sum(x_n) / sqrt(vsq_sum(x_n));
        cout << x_n << "\t" << estimate << endl;
    }
    */

    /*
    long int x_n = pow(10, 9);
    float estimate = 1.5 * tau * v_sum(x_n) / sqrt(vsq_sum(x_n));
    cout << x_n << "\t" << estimate << endl;
    */

    /*
    cout << endl;
    for (int i = 1073741825; i <= 1573741824; i = i + 10000000) {
        int x_n = i;
        float estimate = 1.5 * tau * v_sum(x_n) / sqrt(vsq_sum(x_n));
        cout << x_n << "\t" << estimate << endl;
    }
    */
    /*
    long int x_n = pow(10, 9);
    for (int i = 10; i <= 50; i++) {
        float pp = i * -1.0 / 10.0;
        //float estimate = 1.5 * tau * v_sum(x_n) / sqrt(vsq_sum(x_n));
        float est_limit = v_sum(x_n, pp) / sqrt(vsq_sum(x_n, pp));
        cout << pp << "\t" << est_limit << endl;
    }
    */
    /*
    long int x_n = pow(10, 8);
    for (float i = 0.0; i <= 1.0; i=i+0.05) {
        float pp = i * -1.0;
        //float estimate = 1.5 * tau * v_sum(x_n) / sqrt(vsq_sum(x_n));
        long double vsum = v_sum(x_n, pp);
        long double vsqsum = vsq_sum(x_n, pp);
        //float est_limit = v_sum(x_n, pp) / sqrt(vsq_sum(x_n, pp));
        float est_limit = vsum / sqrt(vsqsum);
        cout << pp << "\t" << vsum << "\t" << vsqsum << "\t" << est_limit << endl;
    }
    for (float i = 9.0; i >= 1.0; i=i-0.5) {
        float pp = float(-10.0) / i;
        //float estimate = 1.5 * tau * v_sum(x_n) / sqrt(vsq_sum(x_n));
        //float est_limit = v_sum(x_n, pp) / sqrt(vsq_sum(x_n, pp));
        long double vsum = v_sum(x_n, pp);
        long double vsqsum = vsq_sum(x_n, pp);
        float est_limit = vsum / sqrt(vsqsum);
        cout << pp << "\t" << vsum << "\t" << vsqsum << "\t" << est_limit << endl;
        //cout << pp << "\t" << est_limit << endl;
    }
    */
    /*
    long int x_n = pow(10, 9);
    for (float i = -0.5; i >= -5.01; i=i-0.05) {
        float pp = i;
        //float estimate = 1.5 * tau * v_sum(x_n) / sqrt(vsq_sum(x_n));
        long double vsum = v_sum(x_n, pp);
        long double vsqsum = vsq_sum(x_n, pp);
        //float est_limit = v_sum(x_n, pp) / sqrt(vsq_sum(x_n, pp));
        float est_limit = vsum / sqrt(vsqsum);
        cout << pp << "\t" << vsum << "\t" << vsqsum << "\t" << est_limit << endl;
    }
    */
    /*
    long int x_n = pow(10, 8);
    float pp = -5.0;
    //float estimate = 1.5 * tau * v_sum(x_n) / sqrt(vsq_sum(x_n));
    long double vsum = v_sum(x_n, pp);
    long double vsqsum = vsq_sum(x_n, pp);
    //float est_limit = v_sum(x_n, pp) / sqrt(vsq_sum(x_n, pp));
    float est_limit = vsum / sqrt(vsqsum);
    cout << pp << "\t" << vsum << "\t" << vsqsum << "\t" << est_limit << endl;
    */
    
    /*
    int arry_x2[20] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};
    vector<int> vec_x2(arry_x2, arry_x2 + 20);
    int arry_y2[20] = {1,2,4,3,6,5,20,13,15,19,17,18,14,16,12,9,11,8,10,7};
    vector<int> vec_y2(arry_y2, arry_y2 + 20);
    */
    
    /*
    int arry_x2[20] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};
    vector<int> vec_x2(arry_x2, arry_x2 + 20);
    int arry_y2[20] = {8,12,19,4,3,17,5,6,14,2,11,7,15,20,18,13,1,16,9,10};
    vector<int> vec_y2(arry_y2, arry_y2 + 20);
    */

    //string fpp_file = "est_pplimit_lite.txt";
    //read_fpp(fpp_file);

    /*
    map <float, float>::reverse_iterator mp_riter;
    for (mp_riter = fpp_mp.rbegin(); mp_riter != fpp_mp.rend(); mp_riter++) {
        //cout << mp_riter->first << "\t" << mp_riter->second << endl;
        float pp = mp_riter->first;
        float est_limit = mp_riter->second;
        double tau = w_ranked_tau(vec_x2, vec_y2, pp);
        double estimate = 1.5 * tau * est_limit;
        double pvalue = tau_pvalue(estimate);
        //float estimate = 1.5 * tau * v_sum(x_n) / sqrt(vsq_sum(x_n));
        //long double vsum = v_sum(x_n, pp);
        //long double vsqsum = vsq_sum(x_n, pp);
        //float est_limit = v_sum(x_n, pp) / sqrt(vsq_sum(x_n, pp));
        //float est_limit = vsum / sqrt(vsqsum);
        cout << pp << "\t" << est_limit << "\t" << tau << "\t" << estimate << "\t" << pvalue << endl;
    }
    */
    

    //string fpp_file = "est_pplimit_lite.txt";
    //string rank_file = "Neu-PMA-data-norm-corr-cluster1.rnk";
    //string grit_file = "grit-human_2_gsea-e-6.gmt+pvalue_wtau.txt";
    //read_fpp(fpp_file);
    read_grit_dataset(rset_file);
    //cout << "read dataset finished!" << endl;
    read_rank_data(rdata_file);
    //cout << "read rank data finished!" << endl;
    performe_wrktf_analysis();
    save_result(output_file);

}

void calculate_fpp(int x_n) {
    //long int x_n = pow(10, 9);
    fpp_mp.clear();
    //for (float i = -0.1; i >= -2.01; i = i - 0.1) {
    
    /*
    float pps[13] = {0.03, 0.06, 0.1, 0.15, 0.25, 0.5, 1.0, 2.0, 4.0, 10.0, 15.0, 25.0, 50.0};

    for (int i = 0; i <= 12; i++) {
        float pp = pps[i];
        //float estimate = 1.5 * tau * v_sum(x_n) / sqrt(vsq_sum(x_n));
        long double vsum = v_sum(x_n, pp);
        long double vsqsum = vsq_sum(x_n, pp);
        //float est_limit = v_sum(x_n, pp) / sqrt(vsq_sum(x_n, pp));
        float est_limit = vsum / sqrt(vsqsum);
        fpp_mp[pp] = est_limit;
        //cout << pp << "\t" << vsum << "\t" << vsqsum << "\t" << est_limit << endl;
    }
    */
   for (float fi = -2; fi <= 2.1; fi += 0.2) {
        //cout << fi << endl;
        float pp = pow(10.0, fi);
        //float estimate = 1.5 * tau * v_sum(x_n) / sqrt(vsq_sum(x_n));
        long double vsum = v_sum(x_n, pp);
        long double vsqsum = vsq_sum(x_n, pp);
        //float est_limit = v_sum(x_n, pp) / sqrt(vsq_sum(x_n, pp));
        float est_limit = vsum / sqrt(vsqsum);
        fpp_mp[pp] = est_limit;
        //cout << pp << "\t" << vsum << "\t" << vsqsum << "\t" << est_limit << endl;
    }
}

void performe_wrktf_analysis() {
    cout << "tf\tpp\ttau\testimate\tpvalue" << endl;
    //vector<map<string, double>>::iterator gsv_riter;
    map<string, double>::iterator rkm_riter;
    for (int i = 0; i < grit_set.size(); i++) {
        vector<double> vec_x;
        vector<double> vec_y;
        overlap_set(grit_set[i], rank_data, vec_x, vec_y);
        int x_n = vec_x.size();
        calculate_fpp(x_n);
        float pp;
        double tau;
        double estimate;
        double pvalue;
        one_tau_analysis(vec_x, vec_y, pp, tau, estimate, pvalue);
        stat stat_one;
        stat_one.tf_name = tf_names[i];
        stat_one.pp = pp;
        stat_one.tau = tau;
        stat_one.estimate = estimate;
        stat_one.pvalue = pvalue;
        result_stats.push_back(stat_one);
        cout << tf_names[i] << "\t" << pp << "\t" << tau << "\t" << estimate << "\t" << pvalue << endl;
        //cout << tf_names[i] << "\t" << pp << "\t" << tau << "\t" << estimate << "\t" << pvalue << endl;
    }

}

void save_result(string output_file) {
    sort(result_stats.begin(), result_stats.end(), stat_less_p);
    double fdr;
    unsigned int homo_stats_size = result_stats.size();
    for (unsigned int hi = 0; hi < result_stats.size(); ++hi) {
        fdr = result_stats[hi].pvalue * homo_stats_size / (hi + 1.0);
        if (fdr <= 1.0) {
            result_stats[hi].fdr = fdr;
        }
        else {
            result_stats[hi].fdr = 1.0;
        }
    }

    ofstream file_writer;
    //string stat_file = dir + "/" + session + ".stat.txt";
    string stat_file = output_file;
    file_writer.open(output_file.c_str(), ios::out);
    if (!file_writer) mcf_die("Sorry, couldn't open " + stat_file);
    file_writer << "tf\tpp\ttau\testimate\tpvalue\tfdr" << endl;
    //cout << "tf\tpp\ttau\testimate\tpvalue\tfdr" << endl;
    sort(result_stats.begin(), result_stats.end(), stat_less_fdr);
    for (unsigned int ri = 0; ri < result_stats.size(); ++ri) {
        file_writer << result_stats[ri].tf_name << "\t" << result_stats[ri].pp << "\t" << result_stats[ri].tau << "\t";
        file_writer << result_stats[ri].estimate << "\t" << result_stats[ri].pvalue << "\t" << result_stats[ri].fdr << endl;
    }
    file_writer.close();
}

void one_tau_analysis(vector<double> &vec_x, vector<double> &vec_y, float &pp, double &tau, double &estimate, double &pvalue) {
    unsigned n = 0;
    unsigned int len_x = vec_x.size();
    unsigned int len_y = vec_y.size();
    if (len_x < 1 || len_y < 1 || len_x != len_y) {
        cout << "length error: vec_x:" << len_x << " vec_y:" << len_y << endl;
        exit(1);
    }
    else {
        n = len_x;
    }
    
    vector<int> r_order;
    /*
    bubble_asort(vec_x, r_order);
    for (int i = 0; i < r_order.size(); i++) {
        cout << i << ": " << r_order[i] << endl;
    }*/
    vector<int> rank_x;
    vector<int> rank_y;
    vec_rank(vec_x, rank_x);
    vec_rank(vec_y, rank_y);

    bubble_lexsort(rank_x, rank_y, r_order);

    float best_pp = 0.0;
    //float best_est_limit = 0.0;
    double best_tau = 0.0;
    double best_estimate = 0.0;
    double best_pvalue = 1.0;

    map <float, float>::iterator mp_riter;
    for (mp_riter = fpp_mp.begin(); mp_riter != fpp_mp.end(); mp_riter++) {
        //cout << mp_riter->first << "\t" << mp_riter->second << endl;
        float tmp_pp = mp_riter->first;
        float tmp_est_limit = mp_riter->second;
        double tmp_tau = w_ranked_tau(rank_x, rank_y, r_order, tmp_pp);
        double tmp_estimate = 1.5 * tmp_tau * tmp_est_limit;
        double tmp_pvalue = tau_pvalue(tmp_estimate);
        //float estimate = 1.5 * tau * v_sum(x_n) / sqrt(vsq_sum(x_n));
        //long double vsum = v_sum(x_n, pp);
        //long double vsqsum = vsq_sum(x_n, pp);
        //float est_limit = v_sum(x_n, pp) / sqrt(vsq_sum(x_n, pp));
        //float est_limit = vsum / sqrt(vsqsum);
        if (tmp_pvalue < best_pvalue) {
            best_pp = tmp_pp;
            //best_est_limit = est_limit;
            best_tau = tmp_tau;
            best_estimate = tmp_estimate;
            best_pvalue = tmp_pvalue;
        }
    }
    pp = best_pp;
    tau = best_tau;
    estimate = best_estimate;
    pvalue = best_pvalue;
    //cout << best_pp << "\t" << best_est_limit << "\t" << best_tau << "\t" << best_estimate << "\t" << best_pvalue << endl;

}

void overlap_set(map<string, double> &map_x, map<string, double> &map_y, vector<double> &vec_x, vector<double> &vec_y) {
    vec_x.clear();
    vec_y.clear();
    map<string, double>::iterator mp_riter_x;
    map<string, double>::iterator mp_riter_y;
    for (mp_riter_x = map_x.begin(); mp_riter_x != map_x.end(); mp_riter_x++) {
        for (mp_riter_y = map_y.begin(); mp_riter_y != map_y.end(); mp_riter_y++) {
            if (mp_riter_x->first == mp_riter_y->first) {
                vec_x.push_back(mp_riter_x->second);
                vec_y.push_back(mp_riter_y->second);
            }
        }
    }
}

void read_grit_dataset(string grit_file) {
    tf_names.clear();
    grit_set.clear();
    ifstream reader(grit_file.c_str());
    if (!reader) {
        cout << "Sorry, couldn't open file " + grit_file << endl;
        exit(1);
    };

    string tmp_str;
	vector<string> tmp_vec;
	while (!reader.eof()) {
		getline(reader, tmp_str, '\n');
		trim_string(tmp_str);
        if (tmp_str.length() == 0) {
            continue;
        }
        else if (tmp_str[0] == '#') {
            continue;
        }
		//cout<<tmp_str<<endl;
        tmp_vec = str_split(tmp_str, "\t");
        if (tmp_vec.size() < 2 || tmp_vec[0] == "") {
		    cout << "wrong line: " + tmp_str << endl;
            exit(1);
	    }
        else {
            tf_names.push_back(tmp_vec[0]);
            map<string, double> one_set;
            for (int i = 1; i < tmp_vec.size(); i++) {
                vector<string> one_vec = str_split(tmp_vec[i], "|");
                if (one_vec.size() != 2 || one_vec[0] == "" || one_vec[1] == "") {
                    cout << "wrong format: " << tmp_vec[i] << endl;
                    exit(1);
                }
                one_set[one_vec[0]] = atof(one_vec[1].c_str());
            }
            grit_set.push_back(one_set);
        }
		
	}
	reader.close();
}

void read_rank_data(string rank_file) {
    rank_data.clear();
    ifstream reader(rank_file.c_str());
    if (!reader) {
        cout << "Sorry, couldn't open file " + rank_file << endl;
        exit(1);
    };

    string tmp_str;
	vector<string> tmp_vec;
	while (!reader.eof()) {
		getline(reader, tmp_str, '\n');
        //cout << tmp_str << endl;
		trim_string(tmp_str);
        if (tmp_str.length() == 0) {
            continue;
        }
        else if (tmp_str[0] == '#') {
            continue;
        }
		//cout<<tmp_str<<endl;
        tmp_vec = str_split(tmp_str, "\t");
        //cout<<tmp_str<<", 2"<<endl;
        if (tmp_vec.size() < 2 || tmp_vec[0] == "" || tmp_vec[1] == "") {
		    cout << "wrong line: " + tmp_str << endl;
            exit(1);
	    }
        else {
            //cout << tmp_vec[0] << endl;
            rank_data[tmp_vec[0].c_str()] = atof(tmp_vec[1].c_str());
        }
		
	}
	reader.close();

}

/*
void read_fpp(string fpp_file) {
    fpp_mp.clear();
    ifstream reader(fpp_file.c_str());
    if (!reader) {
        cout << "Sorry, couldn't open file " + fpp_file << endl;
        exit(1);
    };

    string tmp_str;
	vector<string> tmp_vec;
	while (!reader.eof()) {
		getline(reader, tmp_str, '\n');
		trim_string(tmp_str);
        if (tmp_str.length() == 0) {
            continue;
        }
        else if (tmp_str[0] == '#') {
            continue;
        }
		//cout<<tmp_str<<endl;
        tmp_vec = str_split(tmp_str, "\t");
        if (tmp_vec.size() != 2) {
		    cout << "wrong line: " + tmp_str << endl;
            exit(1);
	    }
        else {
            fpp_mp[atof(tmp_vec[0].c_str())] = atof(tmp_vec[1].c_str());
        }
		
	}
	reader.close();

    
    //map<float, float>::reverse_iterator mp_riter;
    //for (mp_riter = fpp_mp.rbegin(); mp_riter != fpp_mp.rend(); mp_riter++) {
    //    cout << mp_riter->first << "\t" << mp_riter->second << endl;
    //}
    
}
*/

/*
double weigher(int x, int n, float ratio) {
    //return float(1.0) / (1.0 + x);
    //float weight = 1.0 / floor(n * ratio); 
    if (float(x) / n <= ratio) {
        return 1.0;
    }
    else {
        return 0.0;
    }
    
}
*/

/*
double weigher(int x) {
    return float(1.0) / (1.0 + x);  
}
*/

/*
double weigher(int x, float pp) {
    if (x < 1.0 || pp > 0.0) {
        cout << "wrong x or pp" << endl;
        exit(1);
    }
    //return pow((1.0 + x), pp);
    return pow(x, pp);  
}
*/

double weigher(int x, float pp, int n) {
    //if (x < 1.0 || pp < 1.0) {
    if (x < 1.0 || pp == 0.0) {
        cout << "wrong x or pp" << endl;
        exit(1);
    }
    //return pow((1.0 + x), pp);
    return pow(1.0 * (n - x + 1) / n, pp);  
}

int sgn(int x) {
    if (x < 0) {
        return -1;
    }
    else if (x == 0) {
        return 0;
    }
    else {
        return 1;
    }
}

/*
long double v_mean(int x_n,  float ratio) {
    if (x_n < 1) {
        cout << "x_n < 1" << endl;
        exit(1);
        //return 0.0;
    }
    long double sum_vi = 0.0;
    for (int x = 1; x <= x_n; x++) {
        sum_vi += weigher(x, x_n, ratio);
    }
    return sum_vi / x_n;
}
*/

/*
long double vsq_mean(int x_n,  float ratio) {
    if (x_n < 1) {
        cout << "x_n < 1" << endl;
        exit(1);
        //return 0.0;
    }
    long double sum_vi = 0.0;
    for (int x = 1; x <= x_n; x++) {
        sum_vi += pow(weigher(x, x_n, ratio), 2);
    }
    return sum_vi / x_n;
}
*/

/*
long double v_sum(int x_n) {
    if (x_n < 1) {
        cout << "x_n < 1" << endl;
        exit(1);
        //return 0.0;
    }
    long double sum_vi = 0.0;
    for (int x = 1; x <= x_n; x++) {
        sum_vi += weigher(x);
    }
    return sum_vi;
}
*/

/*
long double vsq_sum(int x_n) {
    if (x_n < 1) {
        cout << "x_n < 1" << endl;
        exit(1);
        //return 0.0;
    }
    long double sum_vi = 0.0;
    for (int x = 1; x <= x_n; x++) {
        sum_vi += pow(weigher(x), 2);
    }
    return sum_vi;
}
*/

long double v_sum(int x_n, float pp) {
    if (x_n < 1) {
        cout << "x_n < 1" << endl;
        exit(1);
        //return 0.0;
    }
    long double sum_vi = 0.0;
    for (int x = 1; x <= x_n; x++) {
        sum_vi += weigher(x, pp, x_n);
    }
    return sum_vi;
}

long double vsq_sum(int x_n, float pp) {
    if (x_n < 1) {
        cout << "x_n < 1" << endl;
        exit(1);
        //return 0.0;
    }
    long double sum_vi = 0.0;
    for (int x = 1; x <= x_n; x++) {
        sum_vi += pow(weigher(x, pp, x_n), 2);
    }
    return sum_vi;
}

/*
double w_ranked_tau(vector<int> &vec_x, vector<int> &vec_y, float pp) {   
    unsigned n = 0;
    unsigned int len_x = vec_x.size();
    unsigned int len_y = vec_y.size();
    if (len_x < 1 || len_y < 1 || len_x != len_y) {
        cout << "length error: vec_x:" << len_x << " vec_y:" << len_y << endl;
        exit(1);
    }
    else {
        n = len_x;
    }
    
    vector<int> r_order;
    
    //bubble_asort(vec_x, r_order);
    //for (int i = 0; i < r_order.size(); i++) {
    //    cout << i << ": " << r_order[i] << endl;
    //}
    vector<int> rank_x;
    vector<int> rank_y;
    vec_rank(vec_x, rank_x);
    vec_rank(vec_y, rank_y);

    bubble_lexsort(rank_x, rank_y, r_order);
    
    //for (int i = 0; i < r_order.size(); i++) {
    //    cout << i << ": " << r_order[i] << "\t" << rank_x[r_order[i]] << "\t" << rank_y[r_order[i]] << endl;
    //}
    
    long double sum_vi = 0.0;
    long double sum_svi = 0.0;
    float vi = weigher(rank_x[r_order[0]], pp, n);
    sum_vi += vi;
    sum_svi += pow(vi, 2);

    for (int i = 1; i < r_order.size(); i++) {
        //if (rank_x[r_order[i]] != rank_x[r_order[i - 1]] || rank_y[r_order[i]] != rank_y[r_order[i - 1]]) {
        //    float vi = weigher(rank_x[r_order[i]]);
        //    sum_vi += vi;
        //    sum_svi += pow(vi, 2);
        //}
        float vi = weigher(rank_x[r_order[i]], pp, n);
        sum_vi += vi;
        sum_svi += pow(vi, 2);
    }
    float tot = pow(sum_vi, 2) - sum_svi; 

    long double tsgn =  0.0;
    for (int i = 1; i < r_order.size(); i++) {
        for (int j = 0; j < i; j++) {
            float vi = weigher(rank_x[r_order[i]], pp, n);
            float vj = weigher(rank_x[r_order[j]], pp, n);
            tsgn += vi * vj * sgn(rank_y[r_order[i]] - rank_y[r_order[j]]);
        }
    }

    float tau = 0.0;
    if (tot != 0.0) {
        tau = float(2.0) / tot * tsgn;
    }

    //float vm = v_mean(1000000, ratio);
    //cout << "vm: " << vm << endl;

    //float vsqm = vsq_mean(1000000, ratio);
    //cout << "vsqm: " << vsqm << endl;

    //float estimate = sqrt(n) * tau * 3 * vm / sqrt(4 * vsqm);
    //float estimate = sqrt(n) * tau * 3 * vm / sqrt(4 * vsqm);
    //float estimate = 1.5 * tau * sqrt(n * ratio);

    return tau;

}
*/

double w_ranked_tau(vector<int> &rank_x, vector<int> &rank_y, vector<int> &r_order, float pp) {   
    unsigned n = 0;
    unsigned int len_x = rank_x.size();
    unsigned int len_y = rank_y.size();
    if (len_x < 1 || len_y < 1 || len_x != len_y) {
        cout << "length error: vec_x:" << len_x << " vec_y:" << len_y << endl;
        exit(1);
    }
    else {
        n = len_x;
    }
    
    bubble_lexsort(rank_x, rank_y, r_order);
    /*
    for (int i = 0; i < r_order.size(); i++) {
        cout << i << ": " << r_order[i] << "\t" << rank_x[r_order[i]] << "\t" << rank_y[r_order[i]] << endl;
    }
    */

    long double sum_vi = 0.0;
    long double sum_svi = 0.0;
    float vi = weigher(rank_x[r_order[0]], pp, n);
    sum_vi += vi;
    sum_svi += pow(vi, 2);

    for (int i = 1; i < r_order.size(); i++) {
        /*
        if (rank_x[r_order[i]] != rank_x[r_order[i - 1]] || rank_y[r_order[i]] != rank_y[r_order[i - 1]]) {
            float vi = weigher(rank_x[r_order[i]]);
            sum_vi += vi;
            sum_svi += pow(vi, 2);
        }
        */
        float vi = weigher(rank_x[r_order[i]], pp, n);
        sum_vi += vi;
        sum_svi += pow(vi, 2);
    }
    float tot = pow(sum_vi, 2) - sum_svi; 

    long double tsgn =  0.0;
    for (int i = 1; i < r_order.size(); i++) {
        for (int j = 0; j < i; j++) {
            float vi = weigher(rank_x[r_order[i]], pp, n);
            float vj = weigher(rank_x[r_order[j]], pp, n);
            tsgn += vi * vj * sgn(rank_y[r_order[i]] - rank_y[r_order[j]]);
        }
    }

    float tau = 0.0;
    if (tot != 0.0) {
        tau = float(2.0) / tot * tsgn;
    }

    //float vm = v_mean(1000000, ratio);
    //cout << "vm: " << vm << endl;

    //float vsqm = vsq_mean(1000000, ratio);
    //cout << "vsqm: " << vsqm << endl;

    //float estimate = sqrt(n) * tau * 3 * vm / sqrt(4 * vsqm);
    //float estimate = sqrt(n) * tau * 3 * vm / sqrt(4 * vsqm);
    //float estimate = 1.5 * tau * sqrt(n * ratio);

    return tau;

}

double tau_pvalue(double estimate) {
    double pvalue;
    double bound;
    double x = estimate;
    double mean = 0;
    double sd = 1;
    double p;
    double q;
    int status;
    int which = 1;
    cdfnor(&which, &p, &q, &x, &mean, &sd, &status, &bound);
    if ( status != 0 ) {
        pvalue = 1.0;
    }
    else {
        pvalue = 1.0 - p;
    }
    return pvalue;
    //cout << "ratio: " << ratio << endl;
    //cout << "tau: " << tau << endl;
    //cout << "estimate: " << estimate << endl;
    //cout << "pvalue: " << pvalue << endl;
}

/*
void w_ranked_tau(vector<int> vec_x, vector<int> vec_y, float ratio) {   
    unsigned n = 0;
    unsigned int len_x = vec_x.size();
    unsigned int len_y = vec_y.size();
    if (len_x < 1 || len_y < 1 || len_x != len_y) {
        cout << "length error: vec_x:" << len_x << " vec_y:" << len_y << endl;
        exit(1);
    }
    else {
        n = len_x;
    }
    
    vector<int> r_order;
    vector<int> rank_x;
    vector<int> rank_y;
    vec_rank(vec_x, rank_x);
    vec_rank(vec_y, rank_y);

    bubble_lexsort(rank_x, rank_y, r_order);

    long double sum_vi = 0.0;
    long double sum_svi = 0.0;
    float vi = weigher(rank_x[r_order[0]], n, ratio);
    sum_vi += vi;
    sum_svi += pow(vi, 2);

    for (int i = 1; i < r_order.size(); i++) {
        float vi = weigher(rank_x[r_order[i]], n, ratio);
        sum_vi += vi;
        sum_svi += pow(vi, 2);
    }
    float tot = pow(sum_vi, 2) - sum_svi; 

    long double tsgn =  0.0;
    for (int i = 1; i < r_order.size(); i++) {
        for (int j = 0; j < i; j++) {
            float vi = weigher(rank_x[r_order[i]], n, ratio);
            float vj = weigher(rank_x[r_order[j]], n, ratio);
            tsgn += vi * vj * sgn(rank_y[r_order[i]] - rank_y[r_order[j]]);
        }
    }

    float tau = 0.0;
    if (tot != 0.0) {
        tau = float(2.0) / tot * tsgn;
    }

    //float vm = v_mean(1000000, ratio);
    //cout << "vm: " << vm << endl;

    //float vsqm = vsq_mean(1000000, ratio);
    //cout << "vsqm: " << vsqm << endl;

    //float estimate = sqrt(n) * tau * 3 * vm / sqrt(4 * vsqm);
    //float estimate = sqrt(n) * tau * 3 * vm / sqrt(4 * vsqm);
    float estimate = 1.5 * tau * sqrt(n * ratio);

    double pvalue;
    double bound;
    double x = estimate;
    double mean = 0;
    double sd = 1;
    double p;
    double q;
    int status;
    int which = 1;
    cdfnor(&which, &p, &q, &x, &mean, &sd, &status, &bound);
    if ( status != 0 ) {
        pvalue = 1.0;
    }
    else {
        pvalue = 1.0 - p;
    }

    cout << "ratio: " << ratio << endl;
    cout << "tau: " << tau << endl;
    cout << "estimate: " << estimate << endl;
    cout << "pvalue: " << pvalue << endl;

}
*/

void vec_rank(vector<int> &vec, vector<int> &rank) {
    rank.clear();
    if (vec.size() < 1) {
        return;
    }
    vector<int> order;
    bubble_asort(vec, order);
    for (int i = 0; i < vec.size(); i++) {
        rank.push_back(0);
    }
    rank[order[0]] = 1;
	for (int i = 1; i < order.size(); i++) {
        if (vec[order[i]] != vec[order[i - 1]]) {
            rank[order[i]] = i + 1;
        }
        else {
            rank[order[i]] = rank[order[i - 1]];
        }
	}
}

void vec_rank(vector<double> &vec, vector<int> &rank) {
    rank.clear();
    if (vec.size() < 1) {
        return;
    }
    vector<int> order;
    bubble_asort(vec, order);
    for (int i = 0; i < vec.size(); i++) {
        rank.push_back(0);
    }
    rank[order[0]] = 1;
	for (int i = 1; i < order.size(); i++) {
        if (vec[order[i]] != vec[order[i - 1]]) {
            rank[order[i]] = i + 1;
        }
        else {
            rank[order[i]] = rank[order[i - 1]];
        }
	}
}

void bubble_asort(vector<int> &nums, vector<int> &order) {
	int temp = 0;
    order.clear();
    for (int i = 0; i < nums.size(); i++) {
        order.push_back(i);
    }
	for (int i = 0; i < order.size() - 1; i++) {
		for (int j = 0; j < order.size() - 1 - i; j++) {
			if (nums[order[j]] > nums[order[j + 1]]) {
				temp = order[j];
				order[j] = order[j + 1];
				order[j + 1] = temp;
			}
		}
	}
}

void bubble_asort(vector<double> &nums, vector<int> &order) {
	int temp = 0;
    order.clear();
    for (int i = 0; i < nums.size(); i++) {
        order.push_back(i);
    }
	for (int i = 0; i < order.size() - 1; i++) {
		for (int j = 0; j < order.size() - 1 - i; j++) {
			if (nums[order[j]] > nums[order[j + 1]]) {
				temp = order[j];
				order[j] = order[j + 1];
				order[j + 1] = temp;
			}
		}
	}
}

void bubble_lexsort(vector<int> &vec_x, vector<int> &vec_y, vector<int> &order) {
    unsigned int len_x = vec_x.size();
    unsigned int len_y = vec_y.size();
    if (len_x < 1 || len_y < 1 || len_x != len_y) {
        cout << "length error: vec_x:" << len_x << " vec_y:" << len_y << endl;
        exit(1);
    }
	int temp = 0;
    order.clear();
    for (int i = 0; i < len_x; i++) {
        order.push_back(i);
    }
	for (int i = 0; i < order.size() - 1; i++) {
		for (int j = 0; j < order.size() - 1 - i; j++) {
			if (vec_x[order[j]] > vec_x[order[j + 1]]) {
				temp = order[j];
				order[j] = order[j + 1];
				order[j + 1] = temp;
			}
            else if (vec_x[order[j]] == vec_x[order[j + 1]]) {
                if (vec_y[order[j]] > vec_y[order[j + 1]]) {
                    temp = order[j];
                    order[j] = order[j + 1];
                    order[j + 1] = temp;
                }
            }
		}
	}
}

vector<string> str_split(string& src, string delimit) {
	string null_subst = "";
	vector<string> strvec;
	if (src.empty() || delimit.empty()) {
		//cout<<"Split: Empty String\n";
		return strvec;
	}
	int deli_len = delimit.length();
	size_t index = string::npos, last_search_position = 0;
	while ((index = src.find(delimit, last_search_position)) != string::npos) {
        //cout << "del-index: " << index << endl;
        strvec.push_back(src.substr(last_search_position, index - last_search_position));
        last_search_position = index + deli_len;
        /*
		if (index == last_search_position) {
			// a blank block
			strvec.push_back(null_subst);
            last_search_position = index + deli_len;
		} else {
			strvec.push_back(src.substr(last_search_position, index - last_search_position));
			last_search_position = index + deli_len;
		}
        */
	}
	string last_one = src.substr(last_search_position);
	if (!last_one.empty()) {
		strvec.push_back(last_one);
	}
	return strvec;
}

void trim_string(string &src) {
	size_t pr = string::npos;
	pr = src.find('\r');
	if(pr != string::npos) {
		src.erase(pr, src.size() - pr);
		//src = src.substr(0, pr);
	}

}

bool stat_less_p(stat a, stat b)
{
  return a.pvalue < b.pvalue; //升序排列，如果改为return a>b，则为降序
}

bool stat_less_fdr(stat a, stat b)
{
  return a.fdr < b.fdr; //升序排列，如果改为return a>b，则为降序
}

void mcf_die(const std::string & message)
{
  std::cerr << message << std::endl;
  exit(EXIT_FAILURE);
}
