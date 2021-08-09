#include<cmath>
#include<algorithm>
#include"global.h"
#include"local_wrapping.h"
#include"wrap_back.h"
#include"global_wrap.h"
#include"post-processing.h"
#include"display.h"
#include"enhancer.h"

const enum
{
	NOT_SET = 0,
	PROCESS = 1,
	SET = 2
};

const Scalar RED = Scalar(0, 0, 255);
const Scalar BLUE = Scalar(255, 0, 0);
const Scalar GREEN = Scalar(0, 255, 0);
const int USER_KEY = EVENT_FLAG_CTRLKEY;
const int radius = 2;

Mesh mesh, outMesh;
Mat scale_img;
int mode = NORMAL;

bool init = false;
Mat grab_mask, bgdmodel, fgdmodel;
int rectState = NOT_SET, bgdPointState = NOT_SET, fgdPointState = NOT_SET;
Rect rect;
Point pStart, pEnd;
vector<Point> fgdPixels, bgdPixels;

void showImage()
{
	Mat res, binMask;
	binMask = grab_mask & 1;
	if (init)
	{
		scale_img.copyTo(res, binMask);
	}
	else
	{
		res = scale_img.clone();
	}
	// 显示矩形框
	if (rectState == PROCESS || rectState == SET)
	{
		rectangle(res, rect, RED, 2);
	}

	for (int i = 0; i < bgdPixels.size(); i++)
		circle(res, bgdPixels[i], radius, BLUE, -1);
	for (int i = 0; i < fgdPixels.size(); i++)
		circle(res, fgdPixels[i], radius, GREEN, -1);

	imshow("grab", res);
}

void onMouse(int events, int x, int y, int flag, void*)
{
	if (x < 0 || y < 0 || x > scale_img.cols || y > scale_img.rows)	//无效区域
		return;
	// 使用switch稍微加个速
	switch (events)
	{
	case EVENT_LBUTTONDOWN: // 开始设置矩形框或者背景像素
	{
		bool bgd = (USER_KEY & flag) != 0;
		if (bgd && rectState == SET)
		{
			bgdPointState = PROCESS;
		}
		else if (rectState == NOT_SET && !bgd)
		{
			rectState = PROCESS;
			pStart.x = x;
			pStart.y = y;
		}
	}
	break;
	case EVENT_RBUTTONDOWN:	// 开始设置前景像素
	{
		bool fgd = (USER_KEY & flag) != 0;
		if (fgd)
		{
			fgdPointState = PROCESS;
		}
	}
	break;
	case EVENT_MOUSEMOVE:
		if (rectState == PROCESS)
		{
			pEnd.x = x;
			pEnd.y = y;
			rect = Rect(pStart, pEnd);
			showImage();
		}
		else if (bgdPointState == PROCESS)
		{
			Point p = Point(x, y);
			bgdPixels.push_back(p);
			showImage();
			circle(grab_mask, p, radius, BGD, -1);
		}
		else if (fgdPointState == PROCESS)
		{
			Point p = Point(x, y);
			fgdPixels.push_back(p);
			showImage();
			circle(grab_mask, p, radius, FGD, -1);
		}
		break;
	case EVENT_LBUTTONUP: // 结束设置矩形或者背景像素
		if (rectState == PROCESS)
		{
			rectState = SET;
			pEnd.x = x;
			pEnd.y = y;
			rect = Rect(pStart, pEnd);
			grab_mask.setTo(Scalar::all(BGD));//背景
			grab_mask(rect).setTo(Scalar(FGD_IN_RECT));//前景
		}
		else if (bgdPointState == PROCESS)
		{
			bgdPointState = SET;
			Point p = Point(x, y);
			bgdPixels.push_back(p);
			showImage();
			circle(grab_mask, p, radius, BGD, -1);
		}
		break;
	case EVENT_RBUTTONUP: // 结束设置前景像素
		if (fgdPointState == PROCESS)
		{
			fgdPointState = SET;
			Point p = Point(x, y);
			fgdPixels.push_back(p);
			showImage();
			circle(grab_mask, p, radius, FGD, -1);
		}
		break;
	default:
		break;
	}
}

void reset()
{
	if (!grab_mask.empty())
		grab_mask.setTo(Scalar::all(BGD));
	bgdPixels.clear(); fgdPixels.clear();

	init = false;
	rectState = NOT_SET;
	fgdPointState = NOT_SET;
	bgdPointState = NOT_SET;
	//rect = Rect();
}

void runGrabCut()
{
	if (init)
	{
		grabCut(scale_img, grab_mask, rect, bgdmodel, fgdmodel, 2, GC_INIT_WITH_RECT);
	}
	else
	{
		grabCut(scale_img, grab_mask, rect, bgdmodel, fgdmodel, 2, GC_INIT_WITH_RECT);
		init = true;
	}
}

Mat calcgrayhist(const Mat& image)
{
	Mat histogram = Mat::zeros(Size(256, 1), CV_32SC1);
	//图像宽高
	int rows = image.rows;
	int cols = image.cols;
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < rows; j++)
		{
			int index = int(image.at<uchar>(i, j));
			histogram.at<int>(0, index) += 1;
		}
	}
	return histogram;
};

Rect get_max_rect(Mat dst)
{
	const double ratio = 0.03;

	Mat element = getStructuringElement(cv::MORPH_ELLIPSE, cv::Size(5, 5));
	erode(dst, dst, element);
	dilate(dst, dst, element);
	int min_row = dst.rows, max_row = 0, min_col = dst.cols, max_col = 0;
	bool find_min = false;
	for (int row = 0; row < dst.rows; row++)
	{
		for (int col = 0; col < dst.cols; col++)
		{
			if (dst.at<uchar>(row, col) == 255)
			{
				min_row = min(min_row, row);
				min_col = min(min_col, col);
				max_row = max(max_row, row);
				max_col = max(max_col, col);
			}
		}
	}
	int col_dis = max_col - min_col;
	int row_dis = max_row - min_row;

	int start_row = max(min_row - int(ratio * dst.rows), 0);
	int start_col = max(min_col - int(ratio * dst.cols), 0);
	int end_row = min(max_row + int(ratio * dst.rows), dst.rows);
	int end_col = min(max_col + int(ratio * dst.cols), dst.cols);

	Rect rect(start_col, start_row, end_col - start_col, end_row - start_row);

	return rect;
}

Rect get_max_mask(Mat dst)
{
	Mat element = getStructuringElement(cv::MORPH_ELLIPSE, cv::Size(3, 3));
	erode(dst, dst, element);
	dilate(dst, dst, element);
	int min_row = dst.rows, max_row = 0, min_col = dst.cols, max_col = 0;
	bool find_min = false;
	for (int row = 0; row < dst.rows; row++)
	{
		for (int col = 0; col < dst.cols; col++)
		{
			if (dst.at<uchar>(row, col) == 0)
			{
				min_row = min(min_row, row);
				min_col = min(min_col, col);
				max_row = max(max_row, row);
				max_col = max(max_col, col);
			}
		}
	}
	int col_dis = max_col - min_col;
	int row_dis = max_row - min_row;

	int start_row = max(min_row - 10, 0);
	int start_col = max(min_col - 10, 0);
	int end_row = min(max_row + 10, dst.rows);
	int end_col = min(max_col + 10, dst.cols);

	Rect rect(min_col, min_row, max_col - min_col, max_row - min_row);
	return rect;
}

//OTSU
int OTSU(const Mat& image, Mat& OTSU_image)
{
	int rows = image.rows;
	int cols = image.cols;
	Mat histogram = calcgrayhist(image);
	//归一化直方图
	Mat normhist;
	histogram.convertTo(normhist, CV_32FC1, 1.0 / (rows * cols), 0.0);
	//计算累加直方图和一阶累积矩
	Mat zeroaccumulate = Mat::zeros(Size(256, 1), CV_32FC1);
	Mat oneaccumulate = Mat::zeros(Size(256, 1), CV_32FC1);
	for (int i = 0; i < 256; i++)
	{
		if (i == 0)
		{
			zeroaccumulate.at<float>(0, i) = normhist.at<float>(0, i);
			oneaccumulate.at<float>(0, i) = i * normhist.at<float>(0, i);
		}
		else
		{
			zeroaccumulate.at<float>(0, i) = zeroaccumulate.at<float>(0, i - 1)
				+ normhist.at<float>(0, i);
			oneaccumulate.at<float>(0, i) = oneaccumulate.at<float>(0, i - 1)
				+ i * normhist.at<float>(0, i);
		}
	}
	//计算间类方差
	Mat variance = Mat::zeros(Size(256, 1), CV_32FC1);
	float mean = oneaccumulate.at<float>(0, 255);
	for (int i = 0; i < 255; i++)
	{
		if (zeroaccumulate.at<float>(0, i) == 0 || zeroaccumulate.at<float>(0, i) == 1)
			variance.at<float>(0, i) = 0;
		else
		{
			float cofficient = zeroaccumulate.at<float>(0, i) * (1.0 -
				zeroaccumulate.at<float>(0, i));
			variance.at<float>(0, i) = pow(mean * zeroaccumulate.at<float>(0, i)
				- oneaccumulate.at<float>(0, i), 2.0) / cofficient;
		}
	}
	//找到阈值；
	Point maxloc;
	//计算矩阵中最大值
	minMaxLoc(variance, NULL, NULL, NULL, &maxloc);
	int otsuthreshold = maxloc.x;
	threshold(image, OTSU_image, otsuthreshold, 255, THRESH_BINARY);
	return otsuthreshold;
};

int main(int argc, char* argv[]) {
	Mat gray_img, dst, res;
	Mat src_img = imread("C://Users//12//Desktop//test.jpg");
	int size = 150;
	
	int s = src_img.rows * src_img.cols;
	double scale = sqrt(fixPixelNum / double(s));
	resize(src_img, scale_img, Size(src_img.cols * scale, src_img.rows * scale), 0, 0, INTER_NEAREST);

	cvtColor(scale_img, gray_img, COLOR_BGR2GRAY);
	GaussianBlur(gray_img, gray_img, Size(5, 5), 0);
	OTSU(gray_img, dst);
	rect = get_max_rect(dst);
	
	Mat rect_img = scale_img.clone();
	rectangle(rect_img, rect, RED, 2);
	imshow("grab", rect_img);

	grab_mask.create(scale_img.size(), CV_8U);
	setMouseCallback("grab", onMouse, 0);

	while (true)
	{
		char c = (char)waitKey(0);
		if (c == ' ')
		{
			runGrabCut();
			showImage();
		}
		else if (c == 'r')
		{
			reset();
			showImage();
		}
		else if (c == 'g')
		{
			runGrabCut();
			break;
		}
	}
	double t = (double)getTickCount();
	scale_img.copyTo(res, grab_mask & 1);

	Mat wrapped_img = Mat::zeros(scale_img.size(), CV_8UC3);
	Mat mask(Size(scale_img.cols, scale_img.rows), CV_8UC1);

	cout << scale_img.rows << " " << scale_img.cols << endl;

	// local wrap
	//mask = Mask_contour(res);
	mask = grab_mask & 1;
	mask = mask * 255;
	mask = ~mask;
	Rect mask_rect = get_max_mask(mask);
	Mat c_img, c_mask;
	c_img = scale_img(mask_rect);
	c_mask = mask(mask_rect);
	Mat transres;
	transpose(res, transres);
	flip(transres, transres, 0);
	imshow("grab", transres);
	imwrite("trasn.png", transres);
	//waitKey(1);
	DisMap displacementMap = local_wrap(c_img, c_mask, wrapped_img, mode);

	// wrap back
	mesh = init_mesh(c_img.rows, c_img.cols);
	wrap_mesh(mesh, displacementMap, c_mask);
	

	// global wrap
	outMesh = global_wrap(c_img, c_mask, mesh, mode);

	// Stretching Reduction
	double scale_x = 0;
	double scale_y = 0;
	stretching_reduction(mesh, outMesh, scale_x, scale_y);

	// get result
	Mat outputimg = Mat::zeros(c_img.size(), CV_32FC3);
	Mat ouputcnt = Mat::zeros(c_img.size(), CV_32FC3);

#pragma omp parallel for
	for (int row = 0; row < mesh.mesh_rows - 1; row++) {
		for (int col = 0; col < mesh.mesh_cols - 1; col++) {
			meshPos topleft = outMesh.meshPoints[row][col];
			meshPos topright = outMesh.meshPoints[row][col + 1];
			meshPos bottomleft = outMesh.meshPoints[row + 1][col];
			meshPos bottomright = outMesh.meshPoints[row + 1][col + 1];

			VectorXd Vq = VectorXd::Zero(8);
			Vq << topleft.col, topleft.row, topright.col, topright.row, bottomleft.col, bottomleft.row, bottomright.col, bottomright.row;

			meshPos topleft1 = mesh.meshPoints[row][col];
			meshPos topright1 = mesh.meshPoints[row][col + 1];
			meshPos bottomleft1 = mesh.meshPoints[row + 1][col];
			meshPos bottomright1 = mesh.meshPoints[row + 1][col + 1];

			VectorXd Vo = VectorXd::Zero(8);
			Vo << topleft1.col, topleft1.row, topright1.col, topright1.row, bottomleft1.col, bottomleft1.row, bottomright1.col, bottomright1.row;

			double col_len = max(Vq(0), max(Vq(2), max(Vq(4), Vq(6)))) - min(Vq(0), min(Vq(2), min(Vq(4), Vq(6))));
			double row_len = max(Vq(1), max(Vq(3), max(Vq(5), Vq(7)))) - min(Vq(1), min(Vq(3), min(Vq(5), Vq(7))));
			double col_step = 1 / (4 * col_len);
			double row_step = 1 / (4 * row_len);
			
			for (double i = 0; i < 1; i += row_step) {
				for (double j = 0; j < 1; j += col_step) {
					double v1w = 1 - i - j + i * j;
					double v2w = j - i * j;
					double v3w = i - i * j;
					double v4w = i * j;
					MatrixXd matt(2, 8);
					matt << v1w, 0, v2w, 0, v3w, 0, v4w, 0,
						0, v1w, 0, v2w, 0, v3w, 0, v4w;
					VectorXd pout = matt * Vq;
					VectorXd pref = matt * Vo;
					if (int(pout(1)) >= 0 && int(pout(0)) >= 0 && int(pout(1)) < c_img.rows && int(pout(0)) < c_img.cols) {
						Vec3b pixel = c_img.at<Vec3b>(int(pref(1)), int(pref(0)));
						cv::Vec3f pixelf = cv::Vec3f(float(pixel[0]), float(pixel[1]), float(pixel[2]));
						outputimg.at<Vec3f>(int(pout(1)), int(pout(0))) = outputimg.at<Vec3f>(int(pout(1)), int(pout(0))) + pixelf;
						ouputcnt.at<Vec3f>(int(pout(1)), int(pout(0))) += cv::Vec3f(1, 1, 1);
					}
				}
			}
		}
	}
	Mat finaloutput = outputimg / (255 * ouputcnt);
	Mat result = finaloutput(Range(0, c_img.rows * scale_y), Range(0, c_img.cols * scale_x));
	//Mat enhance = sharpen(result, 1);
	Mat transRes;
	transpose(result, transRes);
	flip(transRes, transRes, 0);
	//transpose(transRes, transRes);
	//transpose(transRes, transRes);
	imshow("final", transRes);
	imwrite("final.png", transRes*255);
	transRes = imread("final.png");
	transRes = remove_shadow(transRes, size);
	imshow("final", transRes);
	imwrite("final.png", transRes);
	t = (double)getTickCount() - t;
	cout << t / (getTickFrequency()) << endl;
	waitKey(0);

	return 0;
}