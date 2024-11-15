//
// Created by Tinghua Huang on 02/12/23.
//

#include <getopt.h>
#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>
#include <math.h>
#include <algorithm>
#include <random>
#include <thread>
#include <mutex>

using namespace std;

struct fdata {
    vector<string> gnames;
    vector<double> vec_x, vec_y;
    vector<int> rank_x, rank_y, rank_diff;
    vector<double> weight_x, weight_y, weights;
    vector<double> rank_m;
    vector<int> ranks;
};

struct fstat {
    string tf_name;
    int count_n;
    double logpp;
    char set_dir;
    char list_dir;
    char cor_dir;
    double tau;
    double estimate;
    double pvalue;
    double fdr;
    fdata data;
    bool operator > (const fstat & r) const
    { return pvalue > r.pvalue; }
    bool operator < (const fstat & r) const
    { return pvalue < r.pvalue; }
};

int distribution_bn = 50;

int sgn(int x);
void bubble_asort(vector<int> &nums, vector<int> &order);
void bubble_asort(vector<double> &nums, vector<int> &order);
void bubble_asort_rev(vector<int> &nums, vector<int> &order);
void bubble_asort_rev(vector<double> &nums, vector<int> &order);
void bubble_lexsort(vector<int> &vec_x, vector<int> &vec_y, vector<int> &order);
int vec_min(vector<int> &vec_data);
int vec_max(vector<int> &vec_data);
double vec_min(vector<double> &vec_data);
double vec_max(vector<double> &vec_data);
double vec_sum(vector<double> &vec_data);
double vec_stdev(vector<double> &vec_data);
void vec_rank(vector<double> &vec, vector<int> &rank);
void vec_rank_rev(vector<double> &vec, vector<int> &rank);
void vec_rank(vector<int> &vec, vector<int> &rank);
void vec_rank_rev(vector<int> &vec, vector<int> &rank);
void trim_string(string &src);
vector<string> str_split(string& src, string delimit);
double tau_pvalue(double estimate);
double tau0(vector<int>& rank_x, vector<int>& rank_y);
double wtau2(vector<int> &rank_x, vector<int> &rank_y, vector<double> &weights);
void read_grit_dataset_bed(string grit_file_bed);
void grit_set_abs();
void grit_set_rev();
void read_rank_data(string rank_file);
void grit_list_abs();
void grit_list_rev();
void overlap_set(map<string, double> &map_x, map<string, double> &map_y, vector<string> &gnames, vector<double> &vec_x, vector<double> &vec_y);
void one_tau_analysis(vector<double> &vec_x, vector<double> &vec_y, double &logp, double &tau, double &estimate, double &pvalue, char &cor_dir, fdata &data_tmp);
void performe_threading();
void performe_wrktf_analysis(vector<string> *index_sect);
void performe_wrktf_analysis_auto_p(vector<string> *index_sect);
bool stat_less_p(struct fstat a, struct fstat b);
bool stat_less_fdr(struct fstat a, struct fstat b);
void mcf_die(const std::string &message);
void save_result(string output_file);
void save_data(string output_dfile);
void print_result();
float round_x(float value, unsigned int dec);
void cal_weight_ranks(vector<double> &vec_data, double &pp, vector<double> &weights, vector<int> &data_rank);
void cal_weights(vector<double> &vec_data, double &pp, vector<double> &weights);
void cal_dweights(vector<double> &vec_data, int bn, vector<double> &weights);
//void cal_dweights_kde(vector<double> &vec_data, vector<double> &weights);
void cal_pdf(vector<double>& data_vec, vector<double>& weights);
double box_pdf(double x, double m, double s);
double default_bandwidth(vector<double>& data_vec);
void cal_weights_ori(vector<double> &vec_data, vector<double> &weights);
double weigher(double x1, double pp);
double calculate_limit(vector<double> &weights);
void abs_vec(vector<double> &vec_data, vector<double> &abs_data);
void merge_weights(vector<double> &weight_x, vector<double> &weight_y, vector<double> &weights);
void merge_ranks(vector<int> &rank_x, vector<int> &rank_y, vector<double> &ranks);
void print_datas(vector<double> &vec_x, vector<double> &vec_y, vector<int> &rank_x, vector<int> &rank_y, vector<double> &weight_x, vector<double> &weight_y, vector<double> &weights);
double normalCDF(double x);

char const short_options[] = "hs:i:g:l:r:w:pt:o:d:b:";
struct option long_options[] =
        {
            {"help", 0, NULL, 'h'},
            {"rank_set", 1, NULL, 's'},
            {"rank_data", 1, NULL, 'i'},
            {"tscore_set", 0, NULL, 'g'},
            {"tscore_list", 0, NULL, 'l'},
            {"adj_method", 0, NULL, 'r'},
            {"wt_method", 0, NULL, 'w'},
            {"auto_p", 0, NULL, 'p'},
            {"thread", 0, NULL, 't'},
            {"output", 0, NULL, 'o'},
            {"outdata", 0, NULL, 'd'},
            {"bin", 0, NULL, 'b'},
            {0, 0, 0, 0}
        };

map<float, float> fpp_mp;
map<string, map<string, double>> grit_set;
vector<string> tf_names;
map<string, double> rank_data;
vector<struct fstat> result_stats;
vector<fdata> result_datas;
double tscore_set = 0.0;
double tscore_list = 0.0;
string wt_method = "7";
string adj_method = "11";
unsigned int thread_n = 6;
mutex thread_lock;
bool mute = true;
bool dump = false;
bool auto_p = false;

int main(int argc, char *argv[]) {
    //cout << normalCDF(1.64) << endl;
    //cout << normalCDF(-1.64) << endl;
    //permutation_test();
    if (argc == 1)
    {
        printf("flaver -s gene_set -i gene_list [-t thread] [-o output]\n");
        exit(1);
    }
    string rset_file = "gene_set.txt";
    string rdata_file = "gene_list.txt";
    string output_file = "output.txt";
    string output_dfile = "outdata.txt";
    //double q_score = 0.05;
    //string output_dir = "./";
    int permutation_n = 100;
    string run_type = "precise";
    int c;
    while((c=getopt_long(argc, argv, short_options, long_options, NULL)) != -1)
    {
        switch (c)
        {
            case 'h':
                printf("flaver -s gene_set -i gene_list [-t thread] [-o output]\n");
                exit(0);
                break;
            case 's':
                rset_file = optarg;
                break;
            case 'i':
                rdata_file = optarg;
                break;
            case 'g':
                tscore_set = atof(optarg);
                break;
            case 'l':
                tscore_list = atof(optarg);
                break;
            case 'r': {
                //cout<<optarg<<endl;
                string optarg_str = optarg;
                if (optarg_str.size() == 2) {
                    adj_method = optarg_str;
                } else {
                    cout << "wrong adj_method: " << optarg_str << endl;
                    exit(1);
                }
                break;
            }
            case 'w': {
                wt_method = optarg;
                break;
            }
            //case 'p':
                //cout<<optarg<<endl;
                //permutation_n = atoi(optarg);
                //run_type = "permutation";
                //break;
            case 't':
                thread_n = atoi(optarg);
                break;
            case 'b':
                distribution_bn = atoi(optarg);
                break;
            case 'o':
                output_file = optarg;
                mute = false;
                break;
            case 'd':
                output_dfile = optarg;
                dump = true;
                break;
            case 'p':
                auto_p = true;
                break;
            default :
                cout << "invalid option: " << optarg << endl;
                exit(1);
        }
    } 

    read_grit_dataset_bed(rset_file);
    //cout << "read dataset finished!" << endl;
    read_rank_data(rdata_file);
    //cout << "read rank data finished!" << endl;
    performe_threading();
    //vector<string> gset_indexs;
    //gset_indexs.push_back("MA0471.2 E2F6");
    //performe_wrktf_analysis(&gset_indexs);
    if (mute == false) {
        save_result(output_file);
    } else {
        print_result();
    }
    if (dump == true) {
        save_data(output_dfile);
    }
}

void cal_weight_ranks(vector<double> &vec_data, double &pp, vector<double> &weights, vector<int> &data_rank) {
    //double pp = log(0.5) / log(logp);
    weights.clear();
    unsigned int x_n = vec_data.size();
    if (x_n < 1) {
        cout << "x_n < 1" << endl;
        exit(1);
        //return 0.0;
    }
    vector<double> abs_data;
    abs_vec(vec_data, abs_data);
    vector<int> abs_rank;
    vec_rank_rev(vec_data, data_rank);
    vec_rank_rev(abs_data, abs_rank);
    for (int x = 0; x < x_n; x++) {
        double x1 = (abs_rank[x] - 0.5) / x_n;
        double tmp_wt = weigher(x1, pp);
        weights.push_back(tmp_wt);
    }
}

void cal_weights(vector<double> &vec_data, double &pp, vector<double> &weights) {
    //double pp = log(0.5) / log(logp);
    weights.clear();
    unsigned int x_n = vec_data.size();
    if (x_n < 1) {
        cout << "x_n < 1" << endl;
        exit(1);
        //return 0.0;
    }
    vector<double> abs_data;
    abs_vec(vec_data, abs_data);
    vector<int> abs_rank;
    vec_rank_rev(abs_data, abs_rank);
    for (int x = 0; x < x_n; x++) {
        double x1 = (abs_rank[x] - 0.5) / x_n;
        double tmp_wt = weigher(x1, pp);
        weights.push_back(tmp_wt);
    }
}

double vec_stdev(vector<double> &vec_data) {
    if (vec_data.size() < 2) {
        return 0.0;
    }
    double avg = vec_sum(vec_data) / vec_data.size();
    double s_total = 0.0;
    for (int i=0; i<vec_data.size(); i++) {
        s_total += pow(vec_data[i] - avg, 2.0);
    }
    double stdev = sqrt(s_total / vec_data.size());
    return stdev;
}

double default_bandwidth(vector<double>& data_vec){
  int x_n = data_vec.size();
  if (x_n < 2) {
    return 1.0;
  }
  double sigma = vec_stdev(data_vec);
  double b = sigma * (pow((3.0 * x_n / 4.0), (-1.0 / 5.0)));
  return b;
}

double box_pdf(double x, double m, double s) {
	if(x < m-0.5*s || x > m+0.5*s){
		return (0.0);
	}else{
		return (1.0/s);
	}
}	

void cal_pdf(vector<double>& data_vec, vector<double>& weights) {
  int x_n = data_vec.size();
  double bandwidth = default_bandwidth(data_vec);
  double data_max = vec_max(data_vec);
  double data_min = vec_min(data_vec);
  for(int i = 0; i < x_n; i++) {
    double data_rg_min = data_min > data_vec[i] - 0.5*bandwidth ? data_min : data_vec[i] - 0.5*bandwidth;
    double data_rg_max = data_max < data_vec[i] + 0.5*bandwidth ? data_max : data_vec[i] + 0.5*bandwidth;
    double da = 0.0;
    for(int j = 0; j < x_n; j++) {
      da += box_pdf(data_vec[j], data_vec[i], bandwidth);
    }
    if (data_rg_max - data_rg_min == 0) {
        weights.push_back(0.0);
    } else {
        weights.push_back(da / double(x_n) * bandwidth / (data_rg_max - data_rg_min));
    }
  }
  double bin_xmax = vec_max(weights);
  if (bin_xmax <= 0.0) {
    //cout << "bin_xmax <= 0.0" << endl;
    //exit(1);
    for (int x = 0; x < x_n; x++) {
        weights[x] = 1.0;
    }
  } else {
    for (int x = 0; x < x_n; x++) {
        weights[x] = 1.0 - (weights[x] / bin_xmax);
    }
  }
}

/*
void cal_dweights_kde(vector<double> &vec_data, vector<double> &weights) {
    weights.clear();
    unsigned int x_n = vec_data.size();
    if (x_n < 1) {
        cout << "x_n < 1" << endl;
        exit(1);
    }
    else if (x_n == 1) {
        weights.push_back(1.0);
        return;
    }
    KDE* kde = new KDE();
    kde->set_kernel_type(2);
    for (int x = 0; x < x_n; x++) {
        vector<double> tmp_data;
        tmp_data.push_back(vec_data[x]);
        kde->add_data(tmp_data);
    }
    //cout << "# bandwidth var 1: " << kde->get_bandwidth(0) << endl;
    for (int x = 0; x < x_n; x++) {
        //cout << vec_data[x] << ":" << kde->pdf(vec_data[x])<< endl;
        weights.push_back(kde->pdf(vec_data[x]));
    }
    double bin_xmax = vec_max(weights);
    if (bin_xmax <= 0.0) {
        cout << "bin_xmax <= 0.0" << endl;
        exit(1);
    }
    for (int x = 0; x < x_n; x++) {
        weights[x] = 1.0 - (weights[x]  / bin_xmax);
    }
}
*/

void cal_dweights(vector<double> &vec_data, int bn, vector<double> &weights) {
    //double pp = log(0.5) / log(logp);
    weights.clear();
    unsigned int x_n = vec_data.size();
    if (x_n < 1) {
        cout << "x_n < 1" << endl;
        exit(1);
        //return 0.0;
    }
    vector<int> order;
    bubble_asort_rev(vec_data, order);
    double max_data = vec_max(vec_data);
    double min_data = vec_min(vec_data);
    double max_data_rg = vec_data[order[round((x_n - 1) * 0.05)]];
    double min_data_rg = vec_data[order[round((x_n - 1) * 0.95)]];
    //double bin = (max_data_rg - min_data_rg) / bn;
    double sigma = vec_stdev(vec_data);
    double bin =  sigma * pow(4.0 / (3.0 * x_n), (1.0 / 5.0));

    double bin_start, bin_end;
    int bin_xn, bin_xn_adj, bin_xmax = 0, bin_xmin = 0.0;
    //bin_xtot = 0;
    double pseudo_min = 1e-6;
    for (int x = 0; x < x_n; x++) {
        bin_start = vec_data[x] - (bin / 2.0);
        bin_end = vec_data[x] + (bin / 2.0);
        if (bin_start < min_data) {
            bin_start = min_data;
            bin_end = bin_start + bin;
        }
        if (bin_end > max_data) {
            bin_end = max_data;
            bin_start = bin_end - bin;
        }
        if (bin_end <= bin_start) {
            weights.push_back(pseudo_min);
            continue;
        }
        bin_xn = 0;
        for (int di = 0; di < x_n; di++) {
            if (vec_data[di] >= bin_start && vec_data[di] <= bin_end) {
                bin_xn++;
            }
        }
        weights.push_back(bin_xn + pseudo_min);
    }
    //double weight_avg = vec_sum(weights) / x_n * bn;
    bin_xmax = vec_max(weights);
    //bin_xmin = vec_min(weights);
    for (int x = 0; x < x_n; x++) {
        weights[x] = 1.0 - (weights[x]  / bin_xmax);
        //weights[x] = 1.0 - ((weights[x] - bin_xmin) / (bin_xmax - bin_xmin));
        //weights[x] = fabs(vec_data[x]) / (weights[x] + 0.25);
        //weights[x] = 1.0 / (weights[x] / weight_avg * (double)x_n + 1.0);
        //weights[x] = 1.0 / weights[x];
        //cout << vec_data[x] << "\t" << weights[x] << endl;
    }
    /*
    bin_xmax = vec_max(weights);
    bin_xmin = vec_min(weights);
    for (int x = 0; x < x_n; x++) {
        weights[x] /= bin_xmax;
        cout << x << ": " << weights[x] << endl;
    }
    */
    /*
    double weight_sum = vec_sum(weights);
    for (int x = 0; x < x_n; x++) {
        weights[x] = weights[x] / weight_sum;
    }
    */
}

/*
void cal_dweights(vector<double> &vec_data, int bn, vector<double> &weights) {
    //double pp = log(0.5) / log(logp);
    weights.clear();
    unsigned int x_n = vec_data.size();
    if (x_n < 1) {
        cout << "x_n < 1" << endl;
        exit(1);
        //return 0.0;
    }
    vector<int> order;
    bubble_asort_rev(vec_data, order);
    double max_data = vec_max(vec_data);
    double min_data = vec_min(vec_data);
    double max_data_rg = vec_data[order[round((x_n - 1) * 0.005)]];
    double min_data_rg = vec_data[order[round((x_n - 1) * 0.995)]];
    double bin = (max_data_rg - min_data_rg) / bn;
    double bin_start, bin_end;
    int bin_xn, bin_xn_adj, bin_xmax = 0, bin_xmin = 0.0;
    //bin_xtot = 0;
    double pseudo_min = 1e-4;
    for (int x = 0; x < x_n; x++) {
        bin_start = vec_data[x] - (bin / 2.0);
        bin_end = vec_data[x] + (bin / 2.0);
        if (bin_start < min_data) {
            bin_start = min_data;
        }
        if (bin_end > max_data) {
            bin_end = max_data;
        }
        if (bin_end <= bin_start) {
            weights.push_back(pseudo_min);
            continue;
        }
        bin_xn = 0;
        for (int di = 0; di < x_n; di++) {
            if (vec_data[di] >= bin_start && vec_data[di] <= bin_end) {
                bin_xn++;
            }
        }
        bin_xn_adj = bin_xn / (bin_end - bin_start) * bin;
        weights.push_back(bin_xn_adj + pseudo_min);
    }
    bin_xmax = vec_max(weights);
    bin_xmin = vec_min(weights);
    for (int x = 0; x < x_n; x++) {
        //weights[x] = 1 - (weights[x] / bin_xmax);
        //weights[x] = fabs(vec_data[x]) / (weights[x] + 0.25);
        weights[x] = 1.0 / (weights[x] - bin_xmin + 1.0);
        //cout << vec_data[x] << "\t" << weights[x] << endl;
    }
}
*/

void cal_weights_ori(vector<double> &vec_data, vector<double> &weights) {
    //double pp = log(0.5) / log(logp);
    weights.clear();
    unsigned int x_n = vec_data.size();
    if (x_n < 1) {
        cout << "x_n < 1" << endl;
        exit(1);
        //return 0.0;
    }
    vector<double> vec_data_abs;
    abs_vec(vec_data, vec_data_abs);
    vector<int> order;
    bubble_asort_rev(vec_data_abs, order);
    //double max_data1 = vec_max(vec_data_abs);
    double max_data = vec_data_abs[order[round((x_n - 1) * 0.01)]];
    //cout << max_data1 << ", " << max_data << endl;
    double min_data = vec_min(vec_data_abs);
    if (max_data <= 0.0) {
        cout << "max_data < 0.0" << endl;
        exit(1);
        //return 0.0;
    }
    for (int x = 0; x < x_n; x++) {
        if (vec_data_abs[x] > max_data) {
            weights.push_back(1.0);
        } else {
            double tmp_wt = (vec_data_abs[x] - min_data) / (max_data - min_data);
            weights.push_back(tmp_wt);
        }
    }
}

void abs_vec(vector<double> &vec_data, vector<double> &abs_data) {
    abs_data.clear();
    unsigned int x_n = vec_data.size();
    if (x_n < 1) {
        cout << "x_n < 1" << endl;
        exit(1);
        //return 0.0;
    }
    for (int x = 0; x < x_n; x++) {
        abs_data.push_back(fabs(vec_data[x]));
    }
}

int vec_min(vector<int> &vec_data) {
    unsigned int x_n = vec_data.size();
    if (x_n < 1) {
        cout << "x_n < 1" << endl;
        exit(1);
        //return 0.0;
    }
    int res_min = vec_data[0];
    for (int x = 1; x < x_n; x++) {
        if (vec_data[x] < res_min) {
            res_min = vec_data[x];
        }
    }
    return res_min;
}

int vec_max(vector<int> &vec_data) {
    unsigned int x_n = vec_data.size();
    if (x_n < 1) {
        cout << "x_n < 1" << endl;
        exit(1);
        //return 0.0;
    }
    int res_min = vec_data[0];
    for (int x = 1; x < x_n; x++) {
        if (vec_data[x] > res_min) {
            res_min = vec_data[x];
        }
    }
    return res_min;
}

double vec_min(vector<double> &vec_data) {
    unsigned int x_n = vec_data.size();
    if (x_n < 1) {
        cout << "x_n < 1" << endl;
        exit(1);
        //return 0.0;
    }
    double res_min = vec_data[0];
    for (int x = 1; x < x_n; x++) {
        if (vec_data[x] < res_min) {
            res_min = vec_data[x];
        }
    }
    return res_min;
}

double vec_max(vector<double> &vec_data) {
    unsigned int x_n = vec_data.size();
    if (x_n < 1) {
        cout << "x_n < 1" << endl;
        exit(1);
        //return 0.0;
    }
    bool valid = false;
    double res_min = 0.0;
    for (int x = 0; x < x_n; x++) {
        if (!isfinite(vec_data[x])) {
            continue;
        }
        else if (!valid) {
            res_min = vec_data[x];
            valid = true;
        }
        else if (vec_data[x] > res_min) {
            res_min = vec_data[x];
        }
    }
    if (!valid) {
        cout << "No valid number in vec_max()!" << endl;
        exit(1);
    }
    return res_min;
}

double vec_sum(vector<double> &vec_data) {
    unsigned int x_n = vec_data.size();
    if (x_n < 1) {
        cout << "x_n < 1" << endl;
        exit(1);
        //return 0.0;
    }
    double res_sum = 0.0;
    for (int x = 0; x < x_n; x++) {
        if (!isfinite(vec_data[x])) {
            continue;
        }
        res_sum += vec_data[x];
    }
    return res_sum;
}

double weigher(double x1, double pp) {
    double wrs = 0.0;
    if (x1 > 1 || x1 < 0) {
        cout << "wrong x1: " << x1 << endl;
        exit(1);
    }

    double y1 = pow(1.0 - pow(x1, pp), 1.0 / pp);
    wrs = y1;
    
    return wrs;
    /*
    //if (x < 1.0 || pp < 1.0) {
    if (x < 1.0 || pp == 0.0) {
        cout << "wrong x or pp: " << x << ", " << pp << endl;
        exit(1);
    }
    //return pow((1.0 + x), pp);
    return pow(1.0 * (n - x + 0.5) / n, pp);
    */ 
}

double calculate_limit(vector<double> &weights) {
    unsigned int x_n = weights.size();
    if (x_n < 1) {
        cout << "x_n < 1" << endl;
        exit(1);
        //return 0.0;
    }
    long double sum_vi = 0.0;
    long double sum_viq = 0.0;
    for (int x = 0; x < x_n; x++) {
        sum_vi += weights[x];
        sum_viq += pow(weights[x], 2.0);
    }
    double est_limit = 0.0;
    if (sum_viq > 0.0) {
        est_limit = (double)(sum_vi / sqrt(sum_viq) / sqrt(x_n));
    }
    return est_limit;
}

void one_tau_analysis(vector<double> &vec_x, vector<double> &vec_y, double &logp, double &tau, double &estimate, double &pvalue, char &cor_dir, fdata &data_tmp) {
    int n = 0;
    int len_x = vec_x.size();
    int len_y = vec_y.size();
    if (len_x < 1 || len_y < 1 || len_x != len_y) {
        cout << "length error3: vec_x:" << len_x << " vec_y:" << len_y << endl;
        exit(1);
    }
    else {
        n = len_x;
    }
    double pp = log(0.5) / log(logp);
    //vector<int> *rank_x = &data_tmp.rank_x;
    //vector<int> *rank_y = &data_tmp.rank_y;
    //vector<double> *weight_x = &data_tmp.weight_x;
    //vector<double> *weight_y = &data_tmp.weight_y;
    //cal_weights(rdm_x, pp, weight_x);
    //cal_weights(rdm_y, pp, weight_y);
    if (wt_method == "3") {
        cal_weight_ranks(vec_x, pp, data_tmp.weights, data_tmp.rank_x);
        //cal_weight_ranks(vec_y, pp, data_tmp.weight_y, data_tmp.rank_y);
        vec_rank_rev(vec_y, data_tmp.rank_y);
        //vector<double> weights;
        //data_tmp.weights = data_tmp.weight_x;
        //merge_weights(data_tmp.weight_x, data_tmp.weight_y, data_tmp.weights);
        //if (logp == 0.8) {
        //    print_datas(vec_x, vec_y, rank_x, rank_y, weight_x, weight_y, weights);
        //}
    }
    else if (wt_method == "2") {
        vec_rank_rev(vec_x, data_tmp.rank_x);
        //cal_weight_ranks(vec_x, pp, data_tmp.weight_x, data_tmp.rank_x);
        cal_weight_ranks(vec_y, pp, data_tmp.weights, data_tmp.rank_y);
        //vector<double> weights;
        //data_tmp.weights = data_tmp.weight_y;
        //merge_weights(data_tmp.weight_x, data_tmp.weight_y, data_tmp.weights);
        //if (logp == 0.8) {
        //    print_datas(vec_x, vec_y, rank_x, rank_y, weight_x, weight_y, weights);
        //}
    /*
    } else if (wt_method == "1") {
        vector<double> vecx_abs, vecy_abs;
        abs_vec(vec_x, vecx_abs);
        abs_vec(vec_y, vecy_abs);
        vector<int> data_rank_x, data_rank_y;
        vec_rank(vecx_abs, data_rank_x);
        vec_rank(vecy_abs, data_rank_y);
        vector<double> mix_ranks;
        merge_ranks(data_rank_x, data_rank_y, mix_ranks);
        cal_weight_ranks(mix_ranks, pp, data_tmp.weights, data_tmp.ranks);
        vec_rank_rev(vec_x, data_tmp.rank_x);
        vec_rank_rev(vec_y, data_tmp.rank_y);
    }*/
    } 
    else if (wt_method == "0") {
        vec_rank_rev(vec_x, data_tmp.rank_x);
        vec_rank_rev(vec_y, data_tmp.rank_y);
        //vector<double> data_weight_x, data_weight_y;
        //cal_weights(vec_x, pp, data_weight_x);
        //cal_weights(vec_y, pp, data_weight_y);
        //merge_weights(data_weight_x, data_weight_y, data_tmp.weights);
        //vec_rank_rev(vec_y, data_tmp.ranks);
    }
    else if (wt_method == "1") {
        vec_rank_rev(vec_x, data_tmp.rank_x);
        vec_rank_rev(vec_y, data_tmp.rank_y);
        vector<double> data_weight_x, data_weight_y;
        cal_weights(vec_x, pp, data_weight_x);
        cal_weights(vec_y, pp, data_weight_y);
        merge_weights(data_weight_x, data_weight_y, data_tmp.weights);
        //vec_rank_rev(vec_y, data_tmp.ranks);
    }
    else if (wt_method == "4") {
        vec_rank_rev(vec_x, data_tmp.rank_x);
        vec_rank_rev(vec_y, data_tmp.rank_y);
        vector<double> data_weight_x, data_weight_y;
        cal_weights_ori(vec_x, data_weight_x);
        cal_weights_ori(vec_y, data_weight_y);
        merge_weights(data_weight_x, data_weight_y, data_tmp.weights);
        //vec_rank_rev(vec_y, data_tmp.ranks);
    }
    else if (wt_method == "5") {
        vec_rank_rev(vec_x, data_tmp.rank_x);
        vec_rank_rev(vec_y, data_tmp.rank_y);
        //cal_weight_ranks(vec_x, pp, data_tmp.weight_x, data_tmp.rank_x);
        //cal_weight_ranks(vec_y, pp, data_tmp.weights, data_tmp.rank_y);
        cal_weights_ori(vec_y, data_tmp.weights);
        //vector<double> weights;
        //data_tmp.weights = data_tmp.weight_y;
        //merge_weights(data_tmp.weight_x, data_tmp.weight_y, data_tmp.weights);
        //if (logp == 0.8) {
        //    print_datas(vec_x, vec_y, rank_x, rank_y, weight_x, weight_y, weights);
        //}
    }
    else if (wt_method == "6") {
        //cal_weight_ranks(vec_x, pp, data_tmp.weights, data_tmp.rank_x);
        //cal_weight_ranks(vec_y, pp, data_tmp.weight_y, data_tmp.rank_y);
        vec_rank_rev(vec_x, data_tmp.rank_x);
        vec_rank_rev(vec_y, data_tmp.rank_y);
        cal_weights_ori(vec_x, data_tmp.weights);
        //vector<double> weights;
        //data_tmp.weights = data_tmp.weight_x;
        //merge_weights(data_tmp.weight_x, data_tmp.weight_y, data_tmp.weights);
        //if (logp == 0.8) {
        //    print_datas(vec_x, vec_y, rank_x, rank_y, weight_x, weight_y, weights);
        //}
    }
    else if (wt_method == "7") {
        //cal_weight_ranks(vec_x, pp, data_tmp.weights, data_tmp.rank_x);
        //cal_weight_ranks(vec_y, pp, data_tmp.weight_y, data_tmp.rank_y);
        vec_rank_rev(vec_x, data_tmp.rank_x);
        vec_rank_rev(vec_y, data_tmp.rank_y);
        vector<double> data_weight_x, data_weight_y;
        cal_pdf(vec_x, data_weight_x);
        cal_pdf(vec_y, data_weight_y);
        merge_weights(data_weight_x, data_weight_y, data_tmp.weights);
        //vector<double> weights;
        //data_tmp.weights = data_tmp.weight_x;
        //merge_weights(data_tmp.weight_x, data_tmp.weight_y, data_tmp.weights);
        //if (logp == 0.8) {
        //    print_datas(vec_x, vec_y, rank_x, rank_y, weight_x, weight_y, weights);
        //}
    }
    else if (wt_method == "8") {
        //cal_weight_ranks(vec_x, pp, data_tmp.weights, data_tmp.rank_x);
        //cal_weight_ranks(vec_y, pp, data_tmp.weight_y, data_tmp.rank_y);
        vec_rank_rev(vec_x, data_tmp.rank_x);
        vec_rank_rev(vec_y, data_tmp.rank_y);
        cal_pdf(vec_y, data_tmp.weights);
        //vector<double> data_weight_x, data_weight_y;
        //cal_dweights(vec_x, distribution_bn, data_weight_x);
        //cal_dweights(vec_y, distribution_bn, data_weight_y);
        //merge_weights(data_weight_x, data_weight_y, data_tmp.weights);
        //vector<double> weights;
        //data_tmp.weights = data_tmp.weight_x;
        //merge_weights(data_tmp.weight_x, data_tmp.weight_y, data_tmp.weights);
        //if (logp == 0.8) {
        //    print_datas(vec_x, vec_y, rank_x, rank_y, weight_x, weight_y, weights);
        //}
    }
    else if (wt_method == "9") {
        //cal_weight_ranks(vec_x, pp, data_tmp.weights, data_tmp.rank_x);
        //cal_weight_ranks(vec_y, pp, data_tmp.weight_y, data_tmp.rank_y);
        vec_rank_rev(vec_x, data_tmp.rank_x);
        vec_rank_rev(vec_y, data_tmp.rank_y);
        cal_pdf(vec_x, data_tmp.weights);
        //vector<double> data_weight_x, data_weight_y;
        //cal_dweights(vec_x, distribution_bn, data_weight_x);
        //cal_dweights(vec_y, distribution_bn, data_weight_y);
        //merge_weights(data_weight_x, data_weight_y, data_tmp.weights);
        //vector<double> weights;
        //data_tmp.weights = data_tmp.weight_x;
        //merge_weights(data_tmp.weight_x, data_tmp.weight_y, data_tmp.weights);
        //if (logp == 0.8) {
        //    print_datas(vec_x, vec_y, rank_x, rank_y, weight_x, weight_y, weights);
        //}
    }
    else {
        cout << "unknown weighting method: " << wt_method << endl;
    }
    double limit;
    //limit = calculate_limit(data_tmp.weights);
    //cout << limit << ", " << n << endl;
    if (wt_method == "0") {
        limit = sqrt(n * (n - 1) / (2 * (2 * n + 5)));
    }
    else if (pp == 1.0 && (wt_method == "2" || wt_method == "3")) {
        limit = 0.866025;
        //cout << "pp 1.0" << endl;
    } else {
        limit = calculate_limit(data_tmp.weights);
    }
    if (wt_method == "0") {
        tau = tau0(data_tmp.rank_x, data_tmp.rank_y);
        estimate = 3.0 * tau * limit;
    }
    else {
        tau = wtau2(data_tmp.rank_x, data_tmp.rank_y, data_tmp.weights);
        estimate = 1.5 * tau * limit * sqrt(n);
    }
    if (estimate >= 0) {
        cor_dir = '+';
    } else {
        cor_dir = '-';
    }
    pvalue = tau_pvalue(estimate);

    //cout << best_pp << "\t" << best_est_limit << "\t" << best_tau << "\t" << best_estimate << "\t" << best_pvalue << endl;
}

void merge_weights(vector<double> &weight_x, vector<double> &weight_y, vector<double> &weights) {
    weights.clear();
    unsigned int x_n = weight_x.size();
    unsigned int y_n = weight_y.size();
    if (x_n != y_n) {
        cout << "x_n != y_n!" << endl;
        exit(1);
    }
    for (unsigned int x = 0; x < x_n; x++) {
        double wi = pow(weight_x[x] * weight_y[x], 0.5);
        //double wi = (weight_x[x] + weight_y[x])/ 2.0;
        weights.push_back(wi);
    }
}

void merge_ranks(vector<int> &rank_x, vector<int> &rank_y, vector<double> &ranks) {
    ranks.clear();
    unsigned int x_n = rank_x.size();
    unsigned int y_n = rank_y.size();
    if (x_n != y_n) {
        cout << "x_n != y_n!" << endl;
        exit(1);
    }
    for (unsigned int x = 0; x < x_n; x++) {
        double wi = pow(rank_x[x] * rank_y[x], 0.5);
        ranks.push_back(wi);
    }
}

void print_datas(vector<double> &vec_x, vector<double> &vec_y, vector<int> &rank_x, vector<int> &rank_y, vector<double> &weight_x, vector<double> &weight_y, vector<double> &weights) {
    //weights.clear();
    unsigned int x_n = weight_x.size();
    unsigned int y_n = weight_y.size();
    unsigned int w_n = weights.size();
    if (x_n != y_n || x_n != w_n) {
        cout << "x_n != y_n != w_n!" << endl;
        exit(1);
    }
    cout << "vec_x\tvec_y\trank_x\trank_y\tgset_wt\tglist_wt\tweight\n";
    for (unsigned int x = 0; x < x_n; x++) {
        cout << vec_x[x] << "\t" << vec_y[x] << "\t" << rank_x[x] << "\t" << rank_y[x] << "\t";
        cout << weight_x[x] << "\t" << weight_y[x] << "\t" << weights[x] << endl;
    }
}

void performe_threading() {
    unsigned int gset_n = grit_set.size();
    vector<string> gset_index_pool;
    map<string, map<string, double>>::iterator iter;
    iter = grit_set.begin();
    while (iter != grit_set.end()) {
        string key = iter->first;
        map<string, double> val = iter->second;
        gset_index_pool.push_back(key);
        iter++;
    }
    unsigned int bunch_size = ceil((double)(gset_n) / thread_n);
    vector<vector<string>> gset_indexs;
    for(unsigned int i = 0; i < gset_n; i += bunch_size) {
        unsigned last = min(gset_n, i + bunch_size);
        vector<string> index_sect;
        index_sect.assign(gset_index_pool.begin() + i, gset_index_pool.begin() + last);
        gset_indexs.push_back(index_sect);
    }
    //thread t1(performe_wrktf_analysis, &(gset_indexs[0]));
    //thread t2(performe_wrktf_analysis, &(gset_indexs[1]));
    //t1.join();
    //t2.join();
    if (mute == false) {
        //cout << "gene-set\tpp\tcor_dir\tset_wt\tlist_wt\ttau\testimate\tpvalue" << endl;
        cout << "Gene-set\tCount_N\tCor_dir\tTau\tEstimate\tP-value" << endl;
    }
    vector<thread> threads;
    for (unsigned int thread_i = 0; thread_i < gset_indexs.size(); thread_i++) {
        //vector<int> index_sect;
        if (auto_p) {
            thread t1(performe_wrktf_analysis_auto_p, &(gset_indexs[thread_i]));
            threads.push_back(move(t1));
        }
        else {
            thread t1(performe_wrktf_analysis, &(gset_indexs[thread_i]));
            threads.push_back(move(t1));
        }
        //t1.join();
        //threads.push_back(move(t1));
        //threads[0]->join();
    }
    //cout << "test point!" << endl;
    //threads[0]->join();
    for (unsigned int thread_i = 0; thread_i < threads.size(); thread_i++) {
        threads[thread_i].join();
    }
}

void performe_wrktf_analysis(vector<string> *index_sect) {
    //cout << "gene-set\tpp\tset_dir\tlist_dir\ttau\testimate\tpvalue" << endl;
    //vector<map<string, double>>::iterator gsv_riter;
    //map<string, double>::iterator rkm_riter;
    for (int ii = 0; ii < index_sect->size(); ii++) {
        string gset_name = index_sect->at(ii);
        //vector<double> vec_x;
        //vector<double> vec_y;
        fdata data_tmp;
        overlap_set(grit_set[gset_name], rank_data, data_tmp.gnames, data_tmp.vec_x, data_tmp.vec_y);
        int x_n = data_tmp.vec_x.size();
        if (x_n < 1) {
            continue;
        }
        //cout << "x_n:" << x_n << endl;
        double logp;
        double tau;
        double estimate;
        double pvalue;
        char cor_dir;
        struct fstat stat_best;
        //fdata data_best;
		stat_best.pvalue = 1.1;

        //char set_wt = wt_method[0];
        //char list_wt = wt_method[1];

        //for (float fi = 0.8; fi >= 0.15; fi -= 0.1) {
        //for (float fi = 0.8; fi >= 0.15; fi -= 0.15) {
        //for (int fi = 0; fi < 5; fi++) {
            //logp = round_x(fi, 2);
            //logp = 0.8 - (0.15 * fi);
            logp = 0.5;
            //char set_dir = (*it)[0];
            //char list_dir = (*it)[1];
            //one_tau_analysis(vec_x, vec_y, logpp, tau, estimate, pvalue, set_dir, list_dir);
            one_tau_analysis(data_tmp.vec_x, data_tmp.vec_y, logp, tau, estimate, pvalue, cor_dir, data_tmp);
            //if (pvalue < stat_best.pvalue) {
                //stat_best.tf_name = tf_names[i];
                stat_best.tf_name = gset_name;
                stat_best.count_n = x_n;
                stat_best.logpp = logp;
                //stat_best.set_dir = set_dir;
                //stat_best.list_dir = list_dir;
                stat_best.cor_dir = cor_dir;
                stat_best.tau = tau;
                stat_best.estimate = estimate;
                stat_best.pvalue = pvalue;
                stat_best.data = data_tmp;
            //}
        //}
        thread_lock.lock();
        result_stats.push_back(stat_best);
        if (mute == false) {
            cout << gset_name << "\t" << stat_best.count_n << "\t" << stat_best.cor_dir;
            //cout << "\t" << stat_best.set_wt << "\t" << stat_best.list_wt;
            cout << "\t" << stat_best.tau << "\t" << stat_best.estimate << "\t" << stat_best.pvalue << endl;
        }
        //cout << tf_names[i] << "\t" << pp << "\t" << tau << "\t" << estimate << "\t" << pvalue << endl;
        thread_lock.unlock();
    }
}

void performe_wrktf_analysis_auto_p(vector<string> *index_sect) {
    //cout << "gene-set\tpp\tset_dir\tlist_dir\ttau\testimate\tpvalue" << endl;
    //vector<map<string, double>>::iterator gsv_riter;
    //map<string, double>::iterator rkm_riter;
    for (int ii = 0; ii < index_sect->size(); ii++) {
        string gset_name = index_sect->at(ii);
        //vector<double> vec_x;
        //vector<double> vec_y;
        fdata data_tmp;
        overlap_set(grit_set[gset_name], rank_data, data_tmp.gnames, data_tmp.vec_x, data_tmp.vec_y);
        int x_n = data_tmp.vec_x.size();
        if (x_n < 1) {
            continue;
        }
        //cout << "x_n:" << x_n << endl;
        double logp;
        double tau;
        double estimate;
        double pvalue;
        char cor_dir;
        struct fstat stat_best;
        //fdata data_best;
        stat_best.pvalue = 1.1;

        //char set_wt = wt_method[0];
        //char list_wt = wt_method[1];

        //for (float fi = 0.8; fi >= 0.15; fi -= 0.1) {
        //for (float fi = 0.8; fi >= 0.15; fi -= 0.15) {
        for (int fi = 0; fi < 5; fi++) {
            //logp = round_x(fi, 2);
            logp = 0.8 - (0.15 * fi);
            //logp = 0.5;
            //char set_dir = (*it)[0];
            //char list_dir = (*it)[1];
            //one_tau_analysis(vec_x, vec_y, logpp, tau, estimate, pvalue, set_dir, list_dir);
            one_tau_analysis(data_tmp.vec_x, data_tmp.vec_y, logp, tau, estimate, pvalue, cor_dir, data_tmp);
            if (pvalue < stat_best.pvalue) {
                //stat_best.tf_name = tf_names[i];
                stat_best.tf_name = gset_name;
                stat_best.count_n = x_n;
                stat_best.logpp = logp;
                //stat_best.set_dir = set_dir;
                //stat_best.list_dir = list_dir;
                stat_best.cor_dir = cor_dir;
                stat_best.tau = tau;
                stat_best.estimate = estimate;
                stat_best.pvalue = pvalue;
                stat_best.data = data_tmp;
            }
        }
        thread_lock.lock();
        result_stats.push_back(stat_best);
        if (mute == false) {
            cout << gset_name << "\t" << stat_best.count_n << "\t" << stat_best.cor_dir;
            //cout << "\t" << stat_best.set_wt << "\t" << stat_best.list_wt;
            cout << "\t" << stat_best.tau << "\t" << stat_best.estimate << "\t" << stat_best.pvalue << endl;
        }
        //cout << tf_names[i] << "\t" << pp << "\t" << tau << "\t" << estimate << "\t" << pvalue << endl;
        thread_lock.unlock();
    }
}

void save_result(string output_file) {
    sort(result_stats.begin(), result_stats.end(), stat_less_p);
    double fdr;
    int homo_stats_size = result_stats.size();
    for (int hi = 0; hi < result_stats.size(); ++hi) {
        fdr = result_stats[hi].pvalue * homo_stats_size / (hi + 1.0);
        if (fdr <= 1.0) {
            result_stats[hi].fdr = fdr;
        }
        else {
            result_stats[hi].fdr = 1.0;
        }
    }
    sort(result_stats.begin(), result_stats.end(), stat_less_fdr);
    ofstream file_writer;
    //string stat_file = dir + "/" + session + ".stat.txt";
    string stat_file = output_file;
    file_writer.open(output_file.c_str(), ios::out);
    if (!file_writer) mcf_die("Sorry, couldn't open " + stat_file);
    file_writer << "Gene-set\tCount_N\tCor_dir\tTau\tEstimate\tP-value\tFDR" << endl;
    //cout << "tf\tpp\ttau\testimate\tpvalue\tfdr" << endl;
    //sort(result_stats.begin(), result_stats.end(), stat_less_fdr);
    for (int ri = 0; ri < result_stats.size(); ++ri) {
        file_writer << result_stats[ri].tf_name << "\t" << result_stats[ri].count_n;
        file_writer << "\t" << result_stats[ri].cor_dir;
        //file_writer << "\t" << result_stats[ri].set_wt << "\t" << result_stats[ri].list_wt;
        file_writer << "\t" << result_stats[ri].tau << "\t";
        file_writer << result_stats[ri].estimate << "\t" << result_stats[ri].pvalue << "\t" << result_stats[ri].fdr << endl;
    }
    file_writer.close();
}

void save_data(string output_dfile) {
    ofstream file_writer;
    file_writer.open(output_dfile.c_str(), ios::out);
    vector<int> rank_diff_tmp;
    vector<int> fin_rank_tmp;
    for (int ri = 0; ri < result_stats.size(); ++ri) {
        if (result_stats[ri].pvalue >= 0.05) {
            continue;
        }
        /*
        if (result_stats[ri].cor_dir == '+') {
            rank_diff_tmp.clear();
            for (int xi = 0; xi < result_stats[ri].data.vec_x.size(); xi++) {
                result_stats[ri].data.rank_m.push_back(pow(result_stats[ri].data.rank_x[xi] * result_stats[ri].data.rank_y[xi], 0.5));
                rank_diff_tmp.push_back(abs(result_stats[ri].data.rank_x[xi] - result_stats[ri].data.rank_y[xi]));
            }
            vec_rank(rank_diff_tmp, result_stats[ri].data.rank_diff);
        } else if (result_stats[ri].cor_dir == '-') {
            rank_diff_tmp.clear();
            int max_rank_y = vec_max(result_stats[ri].data.rank_y);
            int rank_y_tmp;
            for (int xi = 0; xi < result_stats[ri].data.vec_x.size(); xi++) {
                rank_y_tmp = max_rank_y - result_stats[ri].data.rank_y[xi] + 1;
                result_stats[ri].data.rank_m.push_back(pow(result_stats[ri].data.rank_x[xi] * rank_y_tmp, 0.5));
                rank_diff_tmp.push_back(abs(result_stats[ri].data.rank_x[xi] - rank_y_tmp));
            }
            vec_rank(rank_diff_tmp, result_stats[ri].data.rank_diff);
        } else {
            cout << "Unknown direction!" << endl;
            exit(1);
        }
        fin_rank_tmp.clear();
        for (int xi = 0; xi < result_stats[ri].data.vec_x.size(); xi++) {
            fin_rank_tmp.push_back(pow(result_stats[ri].data.rank_m[xi] * result_stats[ri].data.rank_diff[xi], 0.5) / result_stats[ri].data.weights[xi]);
            //result_stats[ri].data.ranks.push_back(result_stats[ri].data.rank_m[xi] / result_stats[ri].data.weights[xi]);
        }
        vec_rank(fin_rank_tmp, result_stats[ri].data.ranks);
        vector<int> fin_rank;
        bubble_asort(result_stats[ri].data.ranks, fin_rank);
        file_writer << "# Gene-set: " << result_stats[ri].tf_name << " (" << result_stats[ri].cor_dir << ")" << endl;
        file_writer << "Gene\tVec_x\tVec_y\tRank_x\tRank_y\tWeight_x\tWeight_y\tWeights\tRank_m\tRank_diff\tRanks" << endl;
        for (int fi = 0; fi < fin_rank.size(); fi++) {
            int rri = fin_rank[fi];
            file_writer << result_stats[ri].data.gnames[rri];
            file_writer << "\t" << result_stats[ri].data.vec_x[rri];
            file_writer << "\t" << result_stats[ri].data.vec_y[rri];
            file_writer << "\t" << result_stats[ri].data.rank_x[rri];
            file_writer << "\t" << result_stats[ri].data.rank_y[rri];
            file_writer << "\t" << result_stats[ri].data.weight_x[rri];
            file_writer << "\t" << result_stats[ri].data.weight_y[rri];
            file_writer << "\t" << result_stats[ri].data.weights[rri];
            file_writer << "\t" << result_stats[ri].data.rank_m[rri];
            file_writer << "\t" << result_stats[ri].data.rank_diff[rri];
            file_writer << "\t" << result_stats[ri].data.ranks[rri];
            file_writer << endl;
        }
        file_writer << "//" << endl << endl;
        */
        vector<int> fin_rank;
        bubble_asort_rev(result_stats[ri].data.weights, fin_rank);
        file_writer << "# Gene-set: " << result_stats[ri].tf_name << " (" << result_stats[ri].cor_dir << ")" << endl;
        file_writer << "Gene\tVec_x\tVec_y\tRank_x\tRank_y\tWeights" << endl;
        for (int fi = 0; fi < fin_rank.size(); fi++) {
            int rri = fin_rank[fi];
            file_writer << result_stats[ri].data.gnames[rri];
            file_writer << "\t" << result_stats[ri].data.vec_x[rri];
            file_writer << "\t" << result_stats[ri].data.vec_y[rri];
            file_writer << "\t" << result_stats[ri].data.rank_x[rri];
            file_writer << "\t" << result_stats[ri].data.rank_y[rri];
            file_writer << "\t" << result_stats[ri].data.weights[rri];
            file_writer << endl;
        }
        file_writer << "//" << endl << endl;
    }
    file_writer.close();
}

void print_result() {
    sort(result_stats.begin(), result_stats.end(), stat_less_p);
    double fdr;
    int homo_stats_size = result_stats.size();
    for (int hi = 0; hi < result_stats.size(); ++hi) {
        fdr = result_stats[hi].pvalue * homo_stats_size / (hi + 1.0);
        if (fdr <= 1.0) {
            result_stats[hi].fdr = fdr;
        }
        else {
            result_stats[hi].fdr = 1.0;
        }
    }
    sort(result_stats.begin(), result_stats.end(), stat_less_fdr);
    cout << "Gene-set\tCount_N\tCor_dir\tTau\tEstimate\tP-value\tFDR" << endl;
    //cout << "gene-set\tpp\tcor_dir\tset_wt\tlist_wt\ttau\testimate\tpvalue\tfdr" << endl;
    //cout << "tf\tpp\ttau\testimate\tpvalue\tfdr" << endl;
    sort(result_stats.begin(), result_stats.end(), stat_less_fdr);
    for (int ri = 0; ri < result_stats.size(); ++ri) {
        cout << result_stats[ri].tf_name << "\t" << result_stats[ri].count_n;
        cout << "\t" << result_stats[ri].cor_dir;
        //cout << "\t" << result_stats[ri].set_wt << "\t" << result_stats[ri].list_wt;
        cout << "\t" << result_stats[ri].tau << "\t";
        cout << result_stats[ri].estimate << "\t" << result_stats[ri].pvalue << "\t" << result_stats[ri].fdr << endl;
    }
}

void overlap_set(map<string, double> &map_x, map<string, double> &map_y, vector<string> &gnames, vector<double> &vec_x, vector<double> &vec_y) {
    gnames.clear();
    vec_x.clear();
    vec_y.clear();
    map<string, double>::iterator mp_riter_x;
    map<string, double>::iterator mp_riter_y;
    for (mp_riter_x = map_x.begin(); mp_riter_x != map_x.end(); mp_riter_x++) {
        for (mp_riter_y = map_y.begin(); mp_riter_y != map_y.end(); mp_riter_y++) {
            if (mp_riter_x->first == mp_riter_y->first) {
                gnames.push_back(mp_riter_x->first);
                vec_x.push_back(mp_riter_x->second);
                vec_y.push_back(mp_riter_y->second);
            }
        }
    }
}

void read_grit_dataset_bed(string grit_file_bed) {
    tf_names.clear();
    grit_set.clear();
    ifstream reader(grit_file_bed.c_str());
    if (!reader) {
        cout << "Sorry, couldn't open file " + grit_file_bed << endl;
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
        if (tmp_vec.size() < 6 || tmp_vec[3] == "") {
            cout << "wrong line: " + tmp_str << endl;
            exit(1);
        } else if (fabs(atof(tmp_vec[4].c_str())) < tscore_set) {
            continue;
        }
        else {
            vector<string> one_vec = str_split(tmp_vec[3], "; ");
            if (one_vec.size() < 2) {
                cout << "wrong line: " + tmp_str << endl;
                exit(1);
            } else {
                vector<string> one_vec0 = str_split(one_vec[0], " ");
                vector<string> one_vec1 = str_split(one_vec[1], "|");
                //if (one_vec0.size() < 2 || one_vec1.size() < 2) {
                if (one_vec1.size() < 2) {
                    cout << "wrong line: " + tmp_str << endl;
                    exit(1);
                }
                grit_set[one_vec[0]][one_vec1[1]] = atof(tmp_vec[4].c_str());
                //save2_grit_set(one_vec[0], one_vec1[1], atof(tmp_vec[4].c_str()));
            }
        }
    }
    reader.close();
    if (adj_method[0] == '2') {
        grit_set_abs();
    } else if (adj_method[0] == '3') {
        grit_set_rev();
    }
}

void grit_set_abs() {
    map<string, map<string, double>>::iterator mmp_riter_x;
    for (mmp_riter_x = grit_set.begin(); mmp_riter_x != grit_set.end(); mmp_riter_x++) {
        map<string, double>::iterator mp_riter_x;
        for (mp_riter_x = mmp_riter_x->second.begin(); mp_riter_x != mmp_riter_x->second.end(); mp_riter_x++) {
            mp_riter_x->second = fabs(mp_riter_x->second);
        }
    }
}

void grit_set_rev() {
    map<string, map<string, double>>::iterator mmp_riter_x;
    for (mmp_riter_x = grit_set.begin(); mmp_riter_x != grit_set.end(); mmp_riter_x++) {
        map<string, double>::iterator mp_riter_x;
        for (mp_riter_x = mmp_riter_x->second.begin(); mp_riter_x != mmp_riter_x->second.end(); mp_riter_x++) {
            mp_riter_x->second = mp_riter_x->second * -1.0;
        }
    }
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
            double score = atof(tmp_vec[1].c_str());
            if (fabs(score) < tscore_list) {
                continue;
            }
            rank_data[tmp_vec[0].c_str()] = atof(tmp_vec[1].c_str());
        }
	}
	reader.close();
    if (adj_method[1] == '2') {
        grit_list_abs();
    } else if (adj_method[1] == '3') {
        grit_list_rev();
    }
}

void grit_list_abs() {
    map<string, double>::iterator mp_riter_x;
    for (mp_riter_x = rank_data.begin(); mp_riter_x != rank_data.end(); mp_riter_x++) {
        mp_riter_x->second = fabs(mp_riter_x->second);
    }
}

void grit_list_rev() {
    map<string, double>::iterator mp_riter_x;
    for (mp_riter_x = rank_data.begin(); mp_riter_x != rank_data.end(); mp_riter_x++) {
        mp_riter_x->second = mp_riter_x->second * -1.0;
    }
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

double tau0(vector<int>& rank_x, vector<int>& rank_y) {
    int n = 0;
    int len_x = rank_x.size();
    int len_y = rank_y.size();
    if (len_x < 1 || len_y < 1 || len_x != len_y) {
        cout << "length error1: rank_x:" << len_x << " rank_y:" << len_y << endl;
        exit(1);
    }
    else {
        n = len_x;
    }
    vector<int> r_order;
    bubble_lexsort(rank_x, rank_y, r_order);
    long double tot = 0.0;
    long double tsgn = 0.0;
    for (int i = 1; i < r_order.size(); i++) {
        for (int j = 0; j < i; j++) {
            //float vi = weigher(rank_x[r_order[i]], logpp, n);
            //float vj = weigher(rank_y[r_order[j]], logpp, n);
            //float vi = weights[r_order[i]];
            //float vj = weights[r_order[j]];
            //tot += vi * vj;
            //tsgn += vi * vj * sgn(rank_x[r_order[i]] - rank_x[r_order[j]]) * sgn(rank_y[r_order[i]] - rank_y[r_order[j]]);
            tsgn += sgn(rank_x[r_order[i]] - rank_x[r_order[j]]) * sgn(rank_y[r_order[i]] - rank_y[r_order[j]]);
        }
    }
    double tau = 0.0;
    if (n > 1) {
        tau = 2.0 * tsgn / (double)(n * (n - 1));
    }
    else {
        tau = 0.0;
    }
    return tau;
}

double wtau2(vector<int> &rank_x, vector<int> &rank_y, vector<double> &weights) {   
    int n = 0;
    int len_x = rank_x.size();
    int len_y = rank_y.size();
    if (len_x < 1 || len_y < 1 || len_x != len_y) {
        cout << "length error1: rank_x:" << len_x << " rank_y:" << len_y << endl;
        exit(1);
    }
    else {
        n = len_x;
    }
    vector<int> r_order;
    bubble_lexsort(rank_x, rank_y, r_order);
    long double tot = 0.0;
    long double tsgn =  0.0;
    for (int i = 1; i < r_order.size(); i++) {
        for (int j = 0; j < i; j++) {
            //float vi = weigher(rank_x[r_order[i]], logpp, n);
            //float vj = weigher(rank_y[r_order[j]], logpp, n);
            float vi = weights[r_order[i]];
            float vj = weights[r_order[j]];
            tot += vi * vj;
            tsgn += vi * vj * sgn(rank_x[r_order[i]] - rank_x[r_order[j]]) * sgn(rank_y[r_order[i]] - rank_y[r_order[j]]);
        }
    }
    double tau = 0.0;
    if (tot != 0.0) {
        tau = (double)(tsgn / tot);
    }
    return tau;
}

double tau_pvalue(double estimate) {
    double raw_p_r = 1.0 - (erfc((fabs(estimate) * -1.0) / sqrt(2.0)) / 2.0);
    double raw_p_l = erfc(fabs(estimate) / sqrt(2.0)) / 2.0;
    double pvalue = raw_p_l + raw_p_r;
    if (pvalue > 1.0) {
        pvalue = 1.0;
    }
    return pvalue;
}

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

void vec_rank_rev(vector<int> &vec, vector<int> &rank) {
    rank.clear();
    if (vec.size() < 1) {
        return;
    }
    vector<int> order;
    bubble_asort_rev(vec, order);
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

void vec_rank_rev(vector<double> &vec, vector<int> &rank) {
    rank.clear();
    if (vec.size() < 1) {
        return;
    }
    vector<int> order;
    bubble_asort_rev(vec, order);
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

void bubble_asort_rev(vector<int> &nums, vector<int> &order) {
	int temp = 0;
    order.clear();
    for (int i = 0; i < nums.size(); i++) {
        order.push_back(i);
    }
	for (int i = 0; i < order.size() - 1; i++) {
		for (int j = 0; j < order.size() - 1 - i; j++) {
			if (nums[order[j]] < nums[order[j + 1]]) {
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

void bubble_asort_rev(vector<double> &nums, vector<int> &order) {
	int temp = 0;
    order.clear();
    for (int i = 0; i < nums.size(); i++) {
        order.push_back(i);
    }
	for (int i = 0; i < order.size() - 1; i++) {
		for (int j = 0; j < order.size() - 1 - i; j++) {
			if (nums[order[j]] < nums[order[j + 1]]) {
				temp = order[j];
				order[j] = order[j + 1];
				order[j + 1] = temp;
			}
		}
	}
}

void bubble_lexsort(vector<int> &vec_x, vector<int> &vec_y, vector<int> &order) {
    int len_x = vec_x.size();
    int len_y = vec_y.size();
    if (len_x < 1 || len_y < 1 || len_x != len_y) {
        cout << "length error4: vec_x:" << len_x << " vec_y:" << len_y << endl;
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

bool stat_less_p(struct fstat a, struct fstat b)
{
  return a.pvalue < b.pvalue;
}

bool stat_less_fdr(struct fstat a, struct fstat b)
{
  return a.fdr < b.fdr;
}

void mcf_die(const std::string & message)
{
  std::cerr << message << std::endl;
  exit(EXIT_FAILURE);
}

float round_x(float value, unsigned int dec) {
    double px = pow(10, dec);
    return (round(value * px) / px);
}

double normalCDF(double x) {
    double cdf = erfc(-x / sqrt(2.0)) / 2.0;
    return cdf;
}
