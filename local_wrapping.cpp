#include"local_wrapping.h"
#include"display.h"

DisMap local_wrap(Mat img, Mat mask, Mat& wrappedImg, int mode)
{
	Mat scaled_img, scaled_mask;
	cv::resize(img, scaled_img, cv::Size(0, 0), 0.5, 0.5);
	cv::resize(mask, scaled_mask, cv::Size(0, 0), 0.5, 0.5);
	DisMap disImg = get_local_dis(scaled_img, scaled_mask, mode);
	for (int i = 0; i < scaled_img.rows; i ++)
	{
		for (int j = 0; j < scaled_img.cols; j++)
		{
			position pos = disImg[i][j];
			wrappedImg.at<Vec3b>(i, j) = scaled_img.at<Vec3b>(i + pos.row, j + pos.col);
		}
	}
	return disImg;
}

DisMap get_local_dis(Mat img, Mat mask, int mode)
{
	int rows = img.rows, cols = img.cols;
	Mat grayImg;
	cvtColor(img, grayImg, COLOR_BGR2GRAY);

	// 定义并初始化disMap
	DisMap displacementMap, finalDisplacementMap;
	for (int row = 0; row < rows; row++) {
		vector<position> displacement_row;
		for (int col = 0; col < cols; col++) {
			position c;
			displacement_row.push_back(c);
		}
		displacementMap.push_back(displacement_row);
		finalDisplacementMap.push_back(displacement_row);
	}

	int i = -1;
	double border_time = 0, seam_time = 0, insert_time = 0, dis_time = 0;
	double t = (double)getTickCount();
	while (true)
	{
		i++;
		int begin = 0, end = 0, maxLength = 0, direction;

		// 寻找到最长边的开始点、结束点、方向
		longgest_border(img, mask, begin, end, direction, maxLength);

		// 边缘没有loss pixel，则结束
		if (end - begin == 0)
		{
			return displacementMap;
		}
		else
		{
			int seamdir = direction % 2; // 上下左右转换为垂直还是水平
			int shiftToRight = direction / 2; // 是否向右shift（bottom的方向旋转后还是朝向右）

			// 获取seam
			vector<int> seam = get_seam(img, mask, seamdir, begin, end);

			// 插入seam
			img = insert_seam(img, mask, seam, seamdir, shiftToRight, begin, end);

			if (mode && i % 5 == 0)
			{
				display_local_warp(img);
				imshow("maskk", mask);
			}
			//更新置换矩阵
			for (int row = 0; row < rows; row++) {
				for (int col = 0; col < cols; col++) {
					position tmpdisplacement;
					if (seamdir == VERTICAL && row >= begin && row <= end) {
						int local_row = row - begin;
						if (col > seam[local_row] && shiftToRight) {
							tmpdisplacement.col = -1;
						}
						else {
							if (col < seam[local_row] && !shiftToRight) {
								tmpdisplacement.col = 1;
							}
						}
					}
					else {
						if (seamdir == HORIZONTAL && col >= begin && col <= end) {
							int local_col = col - begin;
							if (row > seam[local_col] && shiftToRight) {
								tmpdisplacement.row = -1;
							}
							else {
								if (row < seam[local_col] && !shiftToRight) {
									tmpdisplacement.row = 1;
								}
							}
						}
					}
					int tmpdisplace_row = row + tmpdisplacement.row;
					int tmpdisplace_col = col + tmpdisplacement.col;
					position displacementOftarget = displacementMap[tmpdisplace_row][tmpdisplace_col];
					finalDisplacementMap[row][col].row = tmpdisplacement.row + displacementOftarget.row;
					finalDisplacementMap[row][col].col = tmpdisplacement.col + displacementOftarget.col;
				}
			}
			for (int row = 0; row < rows; row++) {
				for (int col = 0; col < cols; col++) {
					position& displacement = displacementMap[row][col];
					position finalDisplacement = finalDisplacementMap[row][col];
					displacement.row = finalDisplacement.row;
					displacement.col = finalDisplacement.col;
				}
			}
			//dis_time += ((double)getTickCount() - t) / getTickFrequency();
			//t = (double)getTickCount();
		}
	}
	return displacementMap;
}

void init_border_index(int& tmp_max, int& tmp_start, int& tmp_end, bool& countingFlag)
{
	tmp_max = 0, tmp_start = 0, tmp_end = 0;
	countingFlag = false;
}

void longgest_border(const Mat img, const Mat mask, int& begin, int& end, int& direction, int& maxLength)
{
	int rows = img.rows, cols = img.cols;

	int tmp_max, tmp_start, tmp_end;
	bool countingFlag;
	// left
	init_border_index(tmp_max, tmp_start, tmp_end, countingFlag);
	for (int row = 0; row < rows; row++)
	{
		// 该点不是透明的或者到达最后一个点
		if (mask.at<uchar>(row, 0) == 0 || row == rows - 1)
		{
			if (countingFlag) // 上一个点还在计数，整理结果
			{
				if (mask.at<uchar>(row, 0) != 0)
				{
					tmp_end++;
					tmp_max++;
				}
				if (tmp_max > maxLength)
				{
					maxLength = tmp_max;
					begin = tmp_start;
					end = tmp_end;
					direction = LEFT;
				}
				countingFlag = false;
			}
			tmp_start = tmp_end = row;
			tmp_max = 0;
		}
		else if (mask.at<uchar>(row, 0) == 255)// 遇到的点是透明的
		{
			tmp_end++;
			tmp_max++;
			countingFlag = true;
		}
	}

	// top
	init_border_index(tmp_max, tmp_start, tmp_end, countingFlag);
	for (int col = 0; col < cols; col++)
	{
		// 该点不是透明的或者到达最后一个点
		if (mask.at<uchar>(0, col) == 0 || col == cols - 1)
		{
			if (countingFlag) // 上一个点还在计数，整理结果
			{
				if (mask.at<uchar>(0, col) != 0)
				{
					tmp_end++;
					tmp_max++;
				}
				if (tmp_max > maxLength)
				{
					maxLength = tmp_max;
					begin = tmp_start;
					end = tmp_end;
					direction = TOP;
				}
				countingFlag = false;
			}
			tmp_start = tmp_end = col;
			tmp_max = 0;
		}
		else // 遇到的点是透明的
		{
			tmp_end++;
			tmp_max++;
			countingFlag = true;
		}
	}

	// right
	init_border_index(tmp_max, tmp_start, tmp_end, countingFlag);
	for (int row = 0; row < rows; row++)
	{
		// 该点不是透明的或者到达最后一个点
		if (mask.at<uchar>(row, cols - 1) == 0 || row == rows - 1)
		{
			if (countingFlag) // 上一个点还在计数，整理结果
			{
				if (mask.at<uchar>(row, cols - 1) != 0)
				{
					tmp_end++;
					tmp_max++;
				}
				if (tmp_max > maxLength)
				{
					maxLength = tmp_max;
					begin = tmp_start;
					end = tmp_end;
					direction = RIGHT;
				}
				countingFlag = false;
			}
			tmp_start = tmp_end = row;
			tmp_max = 0;
		}
		else // 遇到的点是透明的
		{
			tmp_end++;
			tmp_max++;
			countingFlag = true;
		}
	}

	// bottom
	init_border_index(tmp_max, tmp_start, tmp_end, countingFlag);
	for (int col = 0; col < cols; col++)
	{
		// 该点不是透明的或者到达最后一个点
		if (mask.at<uchar>(rows - 1, col) == 0 || col == cols - 1)
		{
			if (countingFlag) // 上一个点还在计数，整理结果
			{
				if (mask.at<uchar>(rows - 1, col) != 0)
				{
					tmp_end++;
					tmp_max++;
				}
				if (tmp_max > maxLength)
				{
					maxLength = tmp_max;
					begin = tmp_start;
					end = tmp_end;
					direction = BOTTOM;
				}
				countingFlag = false;
			}
			tmp_start = tmp_end = col;
			tmp_max = 0;
		}
		else // 遇到的点是透明的
		{
			tmp_end++;
			tmp_max++;
			countingFlag = true;
		}
	}
	if (end == rows) end = end - 1;
	if (end == cols) end = end - 1;
}

Mat calc_local_grad(const Mat img)
{
	Mat grayImg;
	cvtColor(img, grayImg, COLOR_BGR2GRAY);

	
	Mat grad_x(grayImg.rows, grayImg.cols, CV_8U, Scalar(0));
	Mat grad_y(grayImg.rows, grayImg.cols, CV_8U, Scalar(0));

	Mat kernel_x = (Mat_<float>(3, 3) << -1, 0, 1, -2, 0, 2, -3, 0, 3); //求水平梯度所使用的卷积核（赋初始值）
	Mat kernel_y = (Mat_<float>(3, 3) << 1, 2, 1, 0, 0, 0, -1, -2, -1); //求垂直梯度所使用的卷积核（赋初始值）
	
	filter2D(grayImg, grad_x, grad_x.depth(), kernel_x);
	filter2D(grayImg, grad_y, grad_y.depth(), kernel_y);

	Mat gradRes(grayImg.rows, grayImg.cols, CV_8U, Scalar(0));
	add(abs(grad_x), abs(grad_y), gradRes);
	
	/*
	Mat grad_U(grayImg.rows, grayImg.cols, CV_8U, Scalar(0));
	Mat grad_L(grayImg.rows, grayImg.cols, CV_8U, Scalar(0));
	Mat grad_R(grayImg.rows, grayImg.cols, CV_8U, Scalar(0));

	Mat kernel_U = (Mat_<float>(3, 3) << 0, 0, 0, -1, 0, 1, 0, 0, 0);
	Mat kernel_L = (Mat_<float>(3, 3) << 0, 1, 0, -1, 0, 0, 0, 0, 0);
	Mat kernel_R = (Mat_<float>(3, 3) << 0, 1, 0, 0, 0, -1, 0, 0, 0);

	filter2D(grayImg, grad_U, grad_U.depth(), kernel_U);
	filter2D(grayImg, grad_L, grad_L.depth(), kernel_L);
	filter2D(grayImg, grad_R, grad_R.depth(), kernel_R);

	Mat gradtemp(grayImg.rows, grayImg.cols, CV_8U, Scalar(0));
	Mat gradRes(grayImg.rows, grayImg.cols, CV_8U, Scalar(0));

	add(abs(grad_U), abs(grad_L), gradtemp);
	add(abs(gradtemp), abs(grad_R), gradRes);
	*/

	return gradRes;
}

Mat calc_local_energy(const Mat img, const Mat mask)
{
	Mat localGrad = calc_local_grad(img);
	int rows = img.rows, cols = img.cols;
	localGrad.convertTo(localGrad, CV_32F);
	// 把loss pixel的地方设置为INF，防止seam通过

	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			if ((int)mask.at<uchar>(i, j) == 255) localGrad.at<float>(i, j) += INF;
		}
	}

	// dp
	// dp[i][j]: 到达i，j点的最小能量值
	Mat dp;
	localGrad.copyTo(dp);
	for (int i = 1; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			if (j == 0) dp.at<float>(i, j) += min(dp.at<float>(i - 1, j), dp.at<float>(i - 1, j + 1));
			else if (j == cols - 1) dp.at<float>(i, j) += min(dp.at<float>(i - 1, j), dp.at<float>(i - 1, j - 1));
			else
				dp.at<float>(i, j) += min(dp.at<float>(i - 1, j - 1), min(dp.at<float>(i - 1, j), dp.at<float>(i - 1, j + 1)));
		}
	}

	return dp;
}

vector<int> get_seam(Mat img, Mat mask, const int seamdir, const int begin, const int end)
{
	// 如果是寻找水平的seam，将图片旋转即可，该函数只处理垂直方向
	if (seamdir == HORIZONTAL)
	{
		transpose(img, img);
		transpose(mask, mask);
	}

	int cols = img.cols, rows = img.rows, range = end - begin + 1;

	//cout << img.at<Vec3b>(p2) << endl;
	//cout << begin << " " << end << endl;
	Mat sub_img = img(Range(begin, end + 1), Range(0, cols));
	Mat sub_mask = mask(Range(begin, end + 1), Range(0, cols));
	Mat dp = calc_local_energy(sub_img, sub_mask);

	vector<int> seam;
	float minIdx = 0, minEn = INF;
	for (int i = 0; i < cols; i++)
	{
		if (dp.at<float>(range - 1, i) < minEn)
		{
			minIdx = i;
			minEn = dp.at<float>(range - 1, i);
		}
	}
	seam.push_back(minIdx);

	for (int i = 1; i < range; i++)
	{
		int cur_row = range - i - 1;
		if (seam[i - 1] == 0)
		{
			if (dp.at<float>(cur_row, seam[i - 1]) < dp.at<float>(cur_row, seam[i - 1] + 1))
			{
				seam.push_back(seam[i - 1]);
			}
			else
			{
				seam.push_back(seam[i - 1] + 1);
			}
		}
		else if (seam[i - 1] == cols - 1)
		{
			if (dp.at<float>(cur_row, seam[i - 1]) < dp.at<float>(cur_row, seam[i - 1] - 1))
			{
				seam.push_back(seam[i - 1]);
			}
			else seam.push_back(seam[i - 1] - 1);
		}
		else
		{
			int tempMinIdx = 0, tempMin = 0;
			if (dp.at<float>(cur_row, seam[i - 1]) < dp.at<float>(cur_row, seam[i - 1] - 1))
			{
				tempMinIdx = 1;
				tempMin = dp.at<float>(cur_row, seam[i - 1]);
			}
			else
			{
				tempMin = dp.at<float>(cur_row, seam[i - 1] - 1);
			}
			if (dp.at<float>(cur_row, seam[i - 1] + 1) < tempMin)
			{
				tempMinIdx = 2;
			}
			seam.push_back(seam[i - 1] + tempMinIdx - 1);
		}
	}
	//reverse(seam.begin(), seam.end());

	return seam;
}

Mat insert_seam(Mat img, Mat &mask, vector<int> seam, const int seamdir, const int shiftToRight, const int begin, const int end)
{
	if (seamdir == HORIZONTAL)
	{
		transpose(img, img);
		transpose(mask, mask);
	}
	int rows = img.rows, cols = img.cols;
	Mat res, maskOrigin;
	img.copyTo(res);
	mask.copyTo(maskOrigin);

	if (!shiftToRight)
	{
		
//#pragma omp parallel for
		for (int i = begin; i <= end; i++)
		{
			int cur_row = end - i;

			// 更新mask
			mask.at<uchar>(i, seam[cur_row]) = 0;

			// 更新seam处的像素
			if (seam[cur_row] == 0) res.at<Vec3b>(i, 0) = img.at<Vec3b>(i, 1);
			else if (seam[cur_row] == cols - 1) res.at<Vec3b>(i, cols - 1) = img.at<Vec3b>(i, cols - 2);
			else
			{
				res.at<Vec3b>(i, seam[cur_row]) =
					img.at<Vec3b>(i, seam[cur_row] - 1) / 2 + img.at<Vec3b>(i, seam[cur_row] + 1) / 2;
			}

			// shift
			//if (seam[cur_row])
			//{
			//	img(cv::Range(i, i + 1), cv::Range(1, seam[cur_row] + 1)).copyTo(res(cv::Range(i, i + 1), cv::Range(0, seam[cur_row])));
			//	mask(cv::Range(i, i + 1), cv::Range(1, seam[cur_row] + 1)).copyTo(mask(cv::Range(i, i + 1), cv::Range(0, seam[cur_row])));
			//}
			

			/*int offset = seam[cur_row];
			uchar* pImage = img.ptr<uchar>(i);
			uchar* pRes = res.ptr<uchar>(i);
			uchar* pOriMask = maskOrigin.ptr<uchar>(i);
			uchar* pMask = mask.ptr<uchar>(i);

			memcpy(pRes, pImage + 1, sizeof(uchar) * offset * 3);
			memcpy(pMask, pOriMask + 1, sizeof(uchar) * offset);*/

			for (int j = 0; j < seam[cur_row]; j ++)
			{
				res.at<Vec3b>(i, j) = img.at<Vec3b>(i, j + 1);
				mask.at<uchar>(i, j) = mask.at<uchar>(i, j + 1);
			}
			
		}
	}
	else
	{
//#pragma omp parallel for
		for (int i = begin; i <= end; i++)
		{
			int cur_row = end - i;

			mask.at<uchar>(i, seam[cur_row]) = 0;
			if (seam[cur_row] == 0) res.at<Vec3b>(i, 0) = img.at<Vec3b>(i, 1);
			else if (seam[cur_row] == cols - 1) res.at<Vec3b>(i, cols - 1) = img.at<Vec3b>(i, cols - 2);
			else
			{
				res.at<Vec3b>(i, seam[cur_row]) =
					img.at<Vec3b>(i, seam[cur_row] - 1) / 2 + img.at<Vec3b>(i, seam[cur_row] + 1) / 2;
			}
			
			
			/*int offset = cols - seam[cur_row];
			uchar* pImage = img.ptr<uchar>(i);
			uchar* pRes = res.ptr<uchar>(i);
			uchar* pOriMask = maskOrigin.ptr<uchar>(i);
			uchar* pMask = mask.ptr<uchar>(i);

			memcpy(pRes + seam[cur_row], pImage + seam[cur_row] - 1, sizeof(uchar) * offset * 3);
			memcpy(pMask + seam[cur_row], pOriMask + seam[cur_row] - 1, sizeof(uchar) * offset);*/
			
			
			for (int j = cols - 1; j > seam[cur_row]; j--)
			{
				// shift
				res.at<Vec3b>(i, j) = img.at<Vec3b>(i, j - 1);
				mask.at<uchar>(i, j) = mask.at<uchar>(i, j - 1);
			}
			
			
		}
	}
	if (seamdir == HORIZONTAL)
	{
		transpose(res, res);
		transpose(mask, mask);
	}

	return res;
}







