#define FLOAT double

/* Weighted Side EMD Weight Function Type */
enum weightTpye
{
	exp_kernel,
	sin_kernel, 
	spline_kernel
};

int sift(const std::vector<FLOAT> *src, std::vector<FLOAT> *imf);
int stopcondition(const std::vector<FLOAT> *src);
void emd(int, const FLOAT*, FLOAT*, int&, const int, const int, const int);
int wsemd(int, int, int, FLOAT *const, FLOAT *const, FLOAT*, int&, const int, const int, const int);
int extrema(const FLOAT *src, const int n, std::vector<FLOAT> &x_pick, std::vector<FLOAT> &y_pick, std::vector<FLOAT> &x_bottom, std::vector<FLOAT> &y_bottom);
void weightfunction(FLOAT *w, int n, int type, FLOAT alpha);
