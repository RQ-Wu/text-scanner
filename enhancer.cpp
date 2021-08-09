#include"enhancer.h"

Mat sharpen(Mat input, double add_weight)
{
	Mat dst;
	Mat kernel = (Mat_<char>(3, 3) << 0, -1, 0, -1, 5.3, -1, 0, -1, 0);
	filter2D(input, dst, input.depth(), kernel);

	return dst;
}

Mat remove_bgd(Mat input)
{
	Mat blur, usm, dst, res;
	GaussianBlur(input, blur, Size(31, 31), 0);
	divide(input, blur, dst);
	Mat dst2 = sharpen(dst, 0.5);
	dst2.convertTo(res, CV_8UC1, 255);

	return res;
}
Mat remove_shadow(Mat src, int n) {
	vector<Mat> channels;
	split(src, channels);


	for (int i = 0; i < 3; i++)
	{
		//先对原灰度图做最大滤波
		Mat element = getStructuringElement(cv::MorphShapes::MORPH_RECT, cv::Size(n, n));
		int iteration = 1;
		Mat maxFilterMat_A;
		cv::morphologyEx(channels[i], maxFilterMat_A, MORPH_DILATE, element, cv::Point(-1, -1), iteration, cv::BORDER_CONSTANT, cv::morphologyDefaultBorderValue());

		//再对maxFilterMat_A做最小化滤波
		Mat minFilterMat_B;
		cv::morphologyEx(maxFilterMat_A, minFilterMat_B, MORPH_ERODE, element, cv::Point(-1, -1), iteration, cv::BORDER_CONSTANT, cv::morphologyDefaultBorderValue());

		//(先最大再最小滤波) - 原灰度图
		Mat diffMat = minFilterMat_B - channels[i];
		diffMat = 255 - diffMat;
		Mat normalizeMat;
		normalize(diffMat, normalizeMat, 0, 255, NORM_MINMAX);
		channels[i] = normalizeMat.clone();
	}
	merge(channels, src);

	return src;
}