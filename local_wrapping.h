#pragma once
#include"global.h"

DisMap local_wrap(Mat img, Mat mask, Mat& wrappedImg, int mode);
DisMap get_local_dis(Mat img, Mat mask, int mode);
void init_border_index(int& tmp_max, int& tmp_start, int& tmp_end, bool& countingFlag);
void longgest_border(const Mat img, const Mat mask, int& begin, int& end, int& direction, int& maxLength);
Mat calc_local_grad(const Mat img);
Mat calc_local_energy(const Mat img, const Mat mask);
vector<int> get_seam(Mat img, Mat mask, const int seamdir, const int begin, const int end);
Mat insert_seam(Mat img, Mat& mask, vector<int> seam, const int seamdir, const int shiftToRight, const int begin, const int end);