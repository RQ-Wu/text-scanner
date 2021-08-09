#include"global_wrap.h"
#include"display.h"

Mesh global_wrap(const Mat src, Mat mask, Mesh mesh, int mode)
{
	// 准备不会在优化过程中发生改变的量如Aq，e^hat
	SpareseMatrixD_Row shapePrepared = getShapePrepared(mesh);
	SpareseMatrixD_Row map_qua2point = init_q(mesh);
	pair<SpareseMatrixD_Row, VectorXd> pair_dvec_B = get_boundary_mat(src, mesh);
	vector<pair<int, double>> idx_theta;
	vector<double> rotate;
	LineInQua lineInQua = init_lines(mask, src, mesh, idx_theta, rotate, mode);
	Mesh meshUpdated;

	for (int iteration = 0; iteration < 10; iteration++)
	{
		//cout << iteration << endl;
		// =======================fix theta and update V======================= //
		vector<pair<MatrixXd, MatrixXd>> BilinearVec;
		vector<bool> badpoints;
		SpareseMatrixD_Row linePrepared = get_line_mat(src, mask, mesh, rotate, lineInQua, BilinearVec, badpoints);

		// 规定系数
		double lambdaLine = 100;
		double lambdaBoundary = 1e8;
		double Nq = (double)(mesh.mesh_cols - 1) * (mesh.mesh_rows - 1);

		// 得到各个能量函数对应的矩阵
		SpareseMatrixD_Row shapeEnergy = sqrt(1 / Nq) * (shapePrepared * map_qua2point);
		SpareseMatrixD_Row lineEnergy = sqrt(lambdaLine / (linePrepared.rows() / 2)) * (linePrepared * map_qua2point);
		SpareseMatrixD_Row boundaryEnergy = sqrt(lambdaBoundary) * pair_dvec_B.first;

		// 能量函数的最优化（向量范数平方的最小化）转为求 AX = b 的最小二乘解
		// 得到矩阵A和b
		SpareseMatrixD_Row A(shapeEnergy.rows() + lineEnergy.rows() + boundaryEnergy.rows(), shapeEnergy.cols());
		A.topRows(shapeEnergy.rows()) = shapeEnergy;
		A.middleRows(shapeEnergy.rows(), lineEnergy.rows()) = lineEnergy;
		A.bottomRows(boundaryEnergy.rows()) = boundaryEnergy;

		VectorXd b = VectorXd::Zero(A.rows());
		b.tail(pair_dvec_B.second.size()) = sqrt(lambdaBoundary) * pair_dvec_B.second;

		// 为求解最小二乘解做矩阵变换
		// Ax = b ==> x^hat = (A^TA)^-1 * A^T * b
		SparseMatrixD A_trans = A.transpose();
		SparseMatrixD  K = A_trans * A;
		VectorXd bA = A_trans * b;

		CSolve* p_A = new CSolve(K);
		VectorXd x = p_A->solve(bA);

		//VectorXd x = K.inverse() * A_trans * b;

		meshUpdated = updateMesh(x, mesh);
		// =======================fix theta and update V======================= //


		// =======================fix V and update theta======================= //
		int qua_rows = mesh.mesh_rows - 1, qua_cols = mesh.mesh_cols - 1;
		int tempNum = 0;
		double deltathetas[50];
		int thetaCounts[50];
		for (int i = 0; i < 50; i++)
		{
			deltathetas[i] = 0, thetaCounts[i] = 0;
		}
		for (int i = 0; i < qua_rows; i++)
		{
			for (int j = 0; j < qua_cols; j++)
			{
				meshPos topleft = meshUpdated.meshPoints[i][j];
				meshPos topright = meshUpdated.meshPoints[i][j + 1];
				meshPos bottomleft = meshUpdated.meshPoints[i + 1][j];
				meshPos bottomright = meshUpdated.meshPoints[i + 1][j + 1];

				VectorXd Vq = VectorXd::Zero(8);
				Vq << topleft.col, topleft.row, topright.col, topright.row, bottomleft.col, bottomleft.row, bottomright.col, bottomright.row;

				vector<lineSeg> lineInThisQua = lineInQua[i][j];
				for (int k = 0; k < lineInThisQua.size(); k++)
				{
					if (isnan(idx_theta[tempNum].second) || badpoints[tempNum] == true)
					{
						//cout << "badpoints " << badpoints[tempNum] << endl;
						tempNum++;
						continue;
					}

					MatrixXd beginWeight = BilinearVec[tempNum].first;
					MatrixXd endWeight = BilinearVec[tempNum].second;

					Vector2d newBegin = beginWeight * Vq;
					Vector2d newEnd = endWeight * Vq;
					double newtheta = atan((newBegin(1) - newEnd(1)) / (newBegin(0) - newEnd(0)));
					if (isnan(newtheta))
					{
						tempNum++;
						continue;
					}
					
					double deltatheta = newtheta - idx_theta[tempNum].second;
					if (deltatheta < (-PI / 2)) deltatheta += PI;
					if (deltatheta > (PI / 2)) deltatheta -= PI;

					deltathetas[idx_theta[tempNum].first] += deltatheta;
					thetaCounts[idx_theta[tempNum].first]++;

					tempNum++;
				}
			}
		}

		for (int i = 0; i < 50; i++) if (thetaCounts[i]) deltathetas[i] /= thetaCounts[i];
		for (int i = 0; i < rotate.size(); i++)
		{
			rotate[i] = deltathetas[idx_theta[i].first];
			//cout << rotate[i] << endl;
		}
	}
	return meshUpdated;
}

SpareseMatrixD_Row getShapePrepared(Mesh mesh)
{
	int qua_rows = mesh.mesh_rows - 1;
	int qua_cols = mesh.mesh_cols - 1;
	MeshPoints meshps = mesh.meshPoints;

	SpareseMatrixD_Row shapePrepared(8 * qua_rows * qua_cols, 8 * qua_rows * qua_cols);
	for (int i = 0; i < qua_rows; i++)
	{
		for (int j = 0; j < qua_cols; j++)
		{
			meshPos p0 = meshps[i][j];
			meshPos p1 = meshps[i][j + 1];
			meshPos p2 = meshps[i + 1][j];
			meshPos p3 = meshps[i + 1][j + 1];

			MatrixXd Aq(8, 4);
			Aq << p0.col, -p0.row, 1, 0,
				p0.row, p0.col, 0, 1,
				p1.col, -p1.row, 1, 0,
				p1.row, p1.col, 0, 1,
				p2.col, -p2.row, 1, 0,
				p2.row, p2.col, 0, 1,
				p3.col, -p3.row, 1, 0,
				p3.row, p3.col, 0, 1;

			MatrixXd AqT = Aq.transpose();
			MatrixXd I = MatrixXd::Identity(8, 8);
			MatrixXd Aq_trans_mul_Aq_reverse = (AqT * Aq).inverse();
			MatrixXd res = (Aq * (Aq_trans_mul_Aq_reverse) * AqT - I);

			int startPos = (i * qua_cols + j) * 8;
			for (int m = 0; m < 8; m++)
			{
				for (int n = 0; n < 8; n++)
				{
					shapePrepared.insert(startPos + m, startPos + n) = res(m, n);
				}
			}
		}
	}
	shapePrepared.makeCompressed();
	return shapePrepared;
}

SpareseMatrixD_Row init_q(Mesh mesh)
{
	int qua_rows = mesh.mesh_rows - 1;
	int qua_cols = mesh.mesh_cols - 1;
	int mesh_rows = mesh.mesh_rows;
	int mesh_cols = mesh.mesh_cols;

	SpareseMatrixD_Row Q(8 * qua_rows * qua_cols, 2 * mesh_rows * mesh_cols);
	for (int row = 0; row < qua_rows; row++) {
		for (int col = 0; col < qua_cols; col++) {
			int quadid = 8 * (row * qua_cols + col);
			int topleftvertexId = 2 * (row * mesh_cols + col);
			Q.insert(quadid, topleftvertexId) = 1;
			Q.insert(quadid + 1, topleftvertexId + 1) = 1;
			Q.insert(quadid + 2, topleftvertexId + 2) = 1;
			Q.insert(quadid + 3, topleftvertexId + 3) = 1;
			Q.insert(quadid + 4, topleftvertexId + 2 * mesh_cols) = 1;
			Q.insert(quadid + 5, topleftvertexId + 2 * mesh_cols + 1) = 1;
			Q.insert(quadid + 6, topleftvertexId + 2 * mesh_cols + 2) = 1;
			Q.insert(quadid + 7, topleftvertexId + 2 * mesh_cols + 3) = 1;
		}
	}
	Q.makeCompressed();

	return Q;
}

pair<SpareseMatrixD_Row, VectorXd> get_boundary_mat(const Mat src, Mesh mesh) {

	int rows = src.rows;
	int cols = src.cols;
	int numMeshRow = mesh.mesh_rows;
	int numMeshCol = mesh.mesh_cols;
	int vertexnum = numMeshRow * numMeshCol;

	//顺序：[x0,y0,x1,y1...] row为y轴，col为x轴
	VectorXd dvec = VectorXd::Zero(vertexnum * 2);//dvec表标记
	VectorXd B = VectorXd::Zero(vertexnum * 2);//B表位置
	for (int i = 0; i < vertexnum * 2; i += numMeshCol * 2) {//left
		dvec(i) = 1;
		B(i) = 0;
	}//x
	for (int i = numMeshCol * 2 - 2; i < vertexnum * 2; i += numMeshCol * 2)
	{//right
		dvec(i) = 1;
		B(i) = cols - 1;
	}//y

	for (int i = 1; i < 2 * numMeshCol; i += 2)
	{//top
		dvec(i) = 1;
		B(i) = 0;
	}

	for (int i = 2 * vertexnum - 2 * numMeshCol + 1; i < vertexnum * 2; i += 2)
	{//bottom
		dvec(i) = 1;
		B(i) = rows - 1;
	}

	//把dvec换成稀疏矩阵对角来存储和后续计算
	SpareseMatrixD_Row diag(dvec.size(), dvec.size());
	for (int i = 0; i < dvec.size(); i++) {
		diag.insert(i, i) = dvec(i);
	}
	diag.makeCompressed();
	return make_pair(diag, B);
};

//===============================================得到关于line的初始化======================================================//
LineInQua init_lines(Mat &mask, Mat img, Mesh mesh, vector<pair<int, double>>& idx_theta, vector<double>& rotate, int mode)
{
	Mat draw_mask;
	mask.copyTo(draw_mask);
	// 消除mask上的line（将mask腐蚀操作）
	clean_line_on_border_mask(mask);

	// 得到在图片内的所有线段
	vector<lineSeg> lines = get_lines(img, mask);

	// 分配每个线段到对应的quad中
	LineInQua lineInQua = allocate_line_to_qua(lines, mesh);
	if (mode)
	{
		draw_line_in_qua(img, lineInQua, mesh);
	}

	// 构造直方图
	double thetaPerBin = PI / 50;
	for (int i = 0; i < mesh.mesh_rows - 1; i++)
	{
		for (int j = 0; j < mesh.mesh_cols - 1; j++)
		{
			vector<lineSeg> lineInThisQua = lineInQua[i][j];
			for (int k = 0; k < lineInThisQua.size(); k++)
			{
				lineSeg curLine = lineInThisQua[k];
				double theta;
				double slop = (curLine.begin.row - curLine.end.row) / (curLine.begin.col - curLine.end.col);
				if (isnan(slop)) theta = PI - 0.0001;
				else theta = atan((curLine.begin.row - curLine.end.row) / (curLine.begin.col - curLine.end.col));
				int whichOne = (int)((theta + PI / 2) / thetaPerBin);

				idx_theta.push_back(make_pair(whichOne, theta));
				rotate.push_back(0);
			}
		}
	}

	return lineInQua;
}

// lsd算法得到所有的line
vector<lineSeg> get_lines(Mat img, Mat mask)
{
	vector<lineSeg> lines;

	// 先调用接口得到out数组
	int rows = img.rows;
	int cols = img.cols;
	Mat gray_img;
	cvtColor(img, gray_img, COLOR_BGR2GRAY);
	double* image = new double[rows * cols];

	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			image[i * cols + j] = gray_img.at<uchar>(i, j);
		}
	}
	double* out;
	int n_out;

	out = lsd(&n_out, image, cols, rows);
	//draw_lsd(out, n_out, img);

	// 遍历并判断是否在mask的边缘上
	for (int i = 0; i < n_out; i++)
	{
		double start_x = out[i * 7];
		double start_y = out[i * 7 + 1];
		double end_x = out[i * 7 + 2];
		double end_y = out[i * 7 + 3];

		if (inside_mask(mask, start_x, start_y, end_x, end_y))
		{
			lineSeg line;
			linePos start, end;

			start.col = start_x;
			start.row = start_y;
			end.col = end_x;
			end.row = end_y;

			line.begin = start;
			line.end = end;

			lines.push_back(line);
		}
	}
	return lines;
}

// 通过对mask进行腐蚀，排除在mask边缘的line
void clean_line_on_border_mask(Mat& mask)
{
	int rows = mask.rows;
	int cols = mask.cols;
	for (int row = 0; row < rows; row++) {
		mask.at<uchar>(row, 0) = 255;
		mask.at<uchar>(row, cols - 1) = 255;
	}
	for (int col = 0; col < cols; col++) {
		mask.at<uchar>(0, col) = 255;
		mask.at<uchar>(rows - 1, col) = 255;
	}
	Mat element = cv::getStructuringElement(cv::MORPH_ELLIPSE, cv::Size(8, 8));
	cv::dilate(mask, mask, element);
	cv::dilate(mask, mask, element);
}

// 判断是否在一个点mask边缘
bool inside_mask(Mat mask, double start_x, double start_y, double end_x, double end_y)
{
	if (mask.at<uchar>(int(start_y), int(start_x)) == 0 || mask.at<uchar>(int(end_y), int(end_x)) == 0) return true;
	else return false;
}

// 将line分配到对应的网格中
LineInQua allocate_line_to_qua(vector<lineSeg> lines, Mesh mesh)
{
	LineInQua lineInQua;
	int qua_cols = mesh.mesh_cols - 1;
	int qua_rows = mesh.mesh_rows - 1;
	MeshPoints meshPoints = mesh.meshPoints;
	int lineState[20000] = { 0 };

	for (int i = 0; i < qua_rows;	i++)
	{
		vector<vector<lineSeg>> row;
		for (int j = 0; j < qua_cols; j++)
		{
			vector<lineSeg> lineInThisQua;
			meshPos leftTop = meshPoints[i][j];
			meshPos rightTop = meshPoints[i][j + 1];
			meshPos leftBottom = meshPoints[i + 1][j];
			meshPos rightBottom = meshPoints[i + 1][j + 1];

			for (int k = 0; k < lines.size(); k++)
			{
				lineSeg curLine = lines[k];
				if (!lineState[k])
				{
					bool startInQua = isPointInQua(leftTop, rightTop, leftBottom, rightBottom, curLine.begin);
					bool endInQua = isPointInQua(leftTop, rightTop, leftBottom, rightBottom, curLine.end);
					int pointNum = 0;
					vector<linePos> IntersectionPoints;
					// 都在
					if (startInQua && endInQua)
					{
						lineInThisQua.push_back(curLine);
						lineState[k] = 1;
					}
					// begin点在
					else if (startInQua)
					{
						//cout << "(" << curLine.begin.row << "," << curLine.begin.col << ")(" << curLine.end.row << "," << curLine.end.col << ")" << endl;
						//cout << "(" << leftTop.row << "," << leftTop.col << ")(" << rightTop.row << "," << rightTop.col << ")(" << leftBottom.row << "," << leftBottom.col << ")(" << rightBottom.row << "," << rightBottom.col << ")" << endl;

						IntersectionPoints = getIntersectionPoints(leftTop, rightTop, leftBottom, rightBottom, curLine);
						if (IntersectionPoints.size())
						{
							lineSeg temp;
							temp.begin = curLine.begin;
							temp.end = IntersectionPoints[0];
							lineInThisQua.push_back(temp);
						}
					}
					// end点在
					else if (endInQua)
					{
						IntersectionPoints = getIntersectionPoints(leftTop, rightTop, leftBottom, rightBottom, curLine);
						if (IntersectionPoints.size())
						{
							lineSeg temp;
							temp.begin = IntersectionPoints[0];
							temp.end = curLine.end;
							lineInThisQua.push_back(temp);
						}
					}
					// 都不在
					else
					{
						IntersectionPoints = getIntersectionPoints(leftTop, rightTop, leftBottom, rightBottom, curLine);
						if (IntersectionPoints.size() == 2)
						{
							lineSeg temp;
							temp.begin = IntersectionPoints[0];
							temp.end = IntersectionPoints[1];
							lineInThisQua.push_back(temp);
						}
					}
				}
			}
			row.push_back(lineInThisQua);
		}
		lineInQua.push_back(row);
	}

	return lineInQua;
}

// 判断点是否在网格中，因为网格是不规则的四边形，需要用到线性规划的知识
bool isPointInQua(meshPos leftTop, meshPos rightTop, meshPos leftBottom, meshPos rightBottom, linePos point)
{
	// 坐标轴有些反常（为了使所有的点都在第一象限）
	double x = point.row;
	double y = point.col;
	double k, b;

	// 点斜式 y = kx + b， 然后通过线性规划实现
	// 是否在最上面的线的下面
	if (rightTop.row == leftTop.row)
	{
		if (x < rightTop.row) return false;
	}
	else
	{
		k = (rightTop.col - leftTop.col) / (rightTop.row - leftTop.row);
		b = leftTop.col - k * leftTop.row;
		if (x < (y - b) / k) return false;
	}

	// 是否在最下面线的上面
	if (rightBottom.row == leftBottom.row)
	{
		if (x > rightBottom.row) return false;
	}
	else
	{
		k = (rightBottom.col - leftBottom.col) / (rightBottom.row - leftBottom.row);
		b = leftBottom.col - k * leftBottom.row;
		if (x > (y - b) / k) return false;
	}

	// 是否在最左边线的右边
	if (leftTop.col == leftBottom.col)
	{
		if (y < leftTop.col) return false;
	}
	else
	{
		k = (leftTop.col - leftBottom.col) / (leftTop.row - leftBottom.row);
		b = leftBottom.col - k * leftBottom.row;
		if (y < k * x + b) return false;
	}

	// 是否在最右边线的左边
	if (rightTop.col == rightBottom.col)
	{
		if (y > rightTop.col) return false;
	}
	else
	{
		k = (rightTop.col - rightBottom.col) / (rightTop.row - rightBottom.row);
		b = rightBottom.col - k * rightBottom.row;
		if (y > k * x + b) return false;
	}

	return true;
}

// 得到line和该网格的截线
vector<linePos> getIntersectionPoints(meshPos leftTop, meshPos rightTop, meshPos leftBottom, meshPos rightBottom, lineSeg line)
{
	vector<linePos> IntersectionPoints;
	linePos p;

	// 与上方边
	if (getIntersection(leftTop, rightTop, line, p))
	{
		if (isBetween(leftTop, rightTop, p) && isBetween(line.begin, line.end, p))
			IntersectionPoints.push_back(p);
	}

	// 与下方边
	if (getIntersection(leftBottom, rightBottom, line, p))
	{
		if (isBetween(leftBottom, rightBottom, p) && isBetween(line.begin, line.end, p))
			IntersectionPoints.push_back(p);
	}

	// 与左方边
	if (getIntersection(leftTop, leftBottom, line, p))
	{
		if (isBetween(leftTop, leftBottom, p) && isBetween(line.begin, line.end, p))
			IntersectionPoints.push_back(p);
	}

	// 与右方边
	if (getIntersection(rightTop, rightBottom, line, p))
	{
		if (isBetween(rightTop, rightBottom, p) && isBetween(line.begin, line.end, p))
			IntersectionPoints.push_back(p);
	}

	return IntersectionPoints;
}

// 得到一条line和对应的边的交点
bool getIntersection(meshPos a, meshPos b, lineSeg line, linePos& p)
{
	// 通过a，b求出 y = k1x + b1 的表达式
	double k1 = (a.col - b.col) / (a.row - b.row);
	double b1 = a.col - k1 * a.row;

	// 通过line的两个端点求出 y = k2x + b2的表达式
	linePos c = line.begin;
	linePos d = line.end;
	double k2 = (c.col - d.col) / (c.row - d.row);
	double b2 = c.col - k2 * c.row;

	if (k1 == k2)
	{
		return false;
	}
	else if (a.row == b.row && c.row == d.row)
	{
		return false;
	}
	else if (a.row == b.row)
	{
		double x = a.row;
		double y = k2 * x + b2;
		p.row = x;
		p.col = y;
		return true;
	}
	else if (c.row == d.row)
	{
		double x = c.row;
		double y = k1 * x + b1;
		p.row = x;
		p.col = y;
		return true;
	}
	else
	{
		// 列等式 k1x + b1 = k2x + b2求出交点的x坐标
		double x = (b2 - b1) / (k1 - k2);
		double y = k1 * x + b1;
		p.row = x;
		p.col = y;
		return true;
	}
}

// 判断交点是否在两点中间
bool isBetween(meshPos a, meshPos b, meshPos target)
{
	if (target.col >= a.col && target.col <= b.col && target.row <= a.row && target.row >= b.row) return true;
	else if (target.col <= a.col && target.col >= b.col && target.row <= a.row && target.row >= b.row) return true;
	else if (target.col <= a.col && target.col >= b.col && target.row >= a.row && target.row <= b.row) return true;
	else if (target.col >= a.col && target.col <= b.col && target.row >= a.row && target.row <= b.row) return true;
	else return false;
}
//========================================================================================================================//

SpareseMatrixD_Row get_line_mat(Mat img, Mat mask, Mesh mesh, vector<double>rotate_theta, LineInQua lineInQua, 
	vector<pair<MatrixXd, MatrixXd>>& BilinearVec, vector<bool>& badpoints)
{
	int qua_rows = mesh.mesh_rows - 1;
	int qua_cols = mesh.mesh_cols - 1;
	int linenum = 0;

	SpareseMatrixD_Row energy_line(2 * rotate_theta.size(), 8 * qua_rows * qua_cols);
	for (int i = 0; i < qua_rows; i++)
	{
		for (int j = 0; j < qua_cols; j++)
		{
			int QuadID = i * qua_cols + j;
			vector<lineSeg> lineInThisQua = lineInQua[i][j];
			MatrixXd C_row_stack(0, 8);
			for (int k = 0; k < lineInThisQua.size(); k++)
			{
				lineSeg temp = lineInThisQua[k];
				linePos begin = temp.begin;
				linePos end = temp.end;

				MatrixXd begin_weight = bilinear_interpolation(begin, i, j, mesh);
				MatrixXd end_weight = bilinear_interpolation(end, i, j, mesh);

				meshPos topleft = mesh.meshPoints[i][j];
				meshPos topright = mesh.meshPoints[i][j + 1];
				meshPos bottomleft = mesh.meshPoints[i + 1][j];
				meshPos bottomright = mesh.meshPoints[i + 1][j + 1];

				VectorXd Vq = VectorXd::Zero(8);
				Vq << topleft.col, topleft.row, topright.col, topright.row, bottomleft.col, bottomleft.row, bottomright.col, bottomright.row;
				Vector2d ans = begin_weight * Vq - Vector2d(begin.col, begin.row);
				Vector2d ans2 = end_weight * Vq - Vector2d(end.col, end.row);

				if (ans.norm() > 1e-10 || ans2.norm() > 1e-10 || isnan(ans(0)) || isnan(ans2(0)))
				{
					//cout << begin.col << " " << begin.row << endl;
					//cout << end.col << " " << end.row << endl;
					//cout << ans << endl;
					//cout << ans2 << endl;
					badpoints.push_back(true);
					BilinearVec.push_back(make_pair(MatrixXd::Zero(2, 8), MatrixXd::Zero(2, 8)));
					linenum++;
					continue;
				}

				badpoints.push_back(false);
				BilinearVec.push_back(make_pair(begin_weight, end_weight));

				double theta = rotate_theta[linenum];
				MatrixXd R(2, 2);
				MatrixXd ehat(2, 1);
				MatrixXd I = MatrixXd::Identity(2, 2);

				R << cos(theta), -sin(theta), 
					sin(theta), cos(theta);
				ehat << begin.col - end.col, begin.row - end.row;

				MatrixXd e_trans = ehat.transpose();
				MatrixXd C = R * ehat * (e_trans * ehat).inverse() * e_trans * R.transpose() - I;
				MatrixXd CW = C * (begin_weight - end_weight);	

				for (int row = 0; row < 2; row++)
				{
					for (int col = 0; col < 8; col++)
					{
						if (isnan(CW(row, col)))
						{
							cout << begin_weight << endl;
						}
						energy_line.insert(linenum * 2 + row, 8 * QuadID + col) = CW(row, col);
					}
				}
				linenum++;
			}
		}
	}
	
	energy_line.makeCompressed();
	return energy_line;
}
// 不规则四边形的双线性插值：
// If we represent the two end points of a line segment as a bilinear interpolation of its quad vertexes Vq, 
// we can write e as a linear function of Vq
// 救命的博客：https://blog.csdn.net/lk040384/article/details/104939742
MatrixXd bilinear_interpolation(linePos linePoint, int row, int col, Mesh mesh)
{
	double u, v;
	MatrixXd mat(2, 8);

	MeshPoints meshPoints = mesh.meshPoints;;
	meshPos A = meshPoints[row][col];
	meshPos B = meshPoints[row][col + 1];
	meshPos C = meshPoints[row + 1][col + 1];
	meshPos D = meshPoints[row + 1][col];

	doublePos E, F, G, H;
	E.row = B.row - A.row, E.col = B.col - A.col;
	F.row = D.row - A.row, F.col = D.col - A.col;
	G.row = A.row - B.row + C.row - D.row, G.col = A.col - B.col + C.col - D.col;
	H.row = linePoint.row - A.row, H.col = linePoint.col - A.col;

	double k0 = H.col * E.row - H.row * E.col;
	double k1 = E.col * F.row - E.row * F.col + H.col * G.row - H.row * G.col;
	double k2 = G.col * F.row - G.row * F.col;

	if (abs(k2) < 0.01)
	{
		v = -k0 / k1;
		u = (H.col - F.col * v) / (E.col + G.col * v);
	}
	else
	{
		double daita = sqrt(k1 * k1 - 4 * k0 * k2);
		double v1 = (-k1 + daita) / (2 * k2);
		double v2 = (-k1 - daita) / (2 * k2);

		if (v1 >= 0 && v1 <= 1) v = v1;
		else v = v2;
		u = (H.col - F.col * v) / (E.col + G.col * v);
	}

	double v1w = 1 - u - v + u * v;
	double v2w = u - u * v;
	double v3w = v - u * v;
	double v4w = u * v;

	mat << v1w, 0, v2w, 0, v3w, 0, v4w, 0,
		0, v1w, 0, v2w, 0, v3w, 0, v4w;

	return mat;
}

Mesh updateMesh(VectorXd x, Mesh& mesh)
{
	Mesh tempMesh;
	tempMesh.mesh_rows = mesh.mesh_rows;
	tempMesh.mesh_cols = mesh.mesh_cols;

	int mesh_rows = mesh.mesh_rows;
	int mesh_cols = mesh.mesh_cols;

	for (int i = 0; i < mesh_rows; i++)
	{
		vector<meshPos> meshRow;
		for (int j = 0; j < mesh_cols; j++)
		{
			meshPos p;
			int id = (i * mesh_cols + j) * 2;
			p.col = x(id);
			p.row = x(id + 1);

			meshRow.push_back(p);
		}
		tempMesh.meshPoints.push_back(meshRow);
	}

	return tempMesh;
}

Mat fill_missing_pixel(Mat& img, const Mat mask)
{
	//row为y轴，col为x轴
	assert(img.rows == mask.rows);
	assert(img.cols == mask.cols);
	Mat mask_2;
	int size_erode = 3;
	Mat element = cv::getStructuringElement(cv::MORPH_RECT, cv::Size(size_erode, size_erode));
	cv::erode(mask, mask_2, element);//255是非图的部分，，腐蚀掉一些这个
	/*cv::imshow("temp", mask_2);
	cv::imshow("temp2", mask);
	cv::waitKey(0);	*/
	for (int row = 0; row < mask.rows; row++)
	{
		for (int col = 0; col < mask.cols; col++)
		{
			if (mask.at<uchar>(row, col) == 255 && mask_2.at<uchar>(row, col) == 0)
			{
				for (int i = 0; i < size_erode; i++)
				{
					int temp_y = row - 2 + i / size_erode;
					int temp_x = col - 2 + i % size_erode;
					if (temp_y >= 0 && temp_y <= mask.rows && temp_x >= 0 && temp_x <= mask.cols)
					{
						if (mask.at<uchar>(temp_y, temp_x) == 0)
						{
							img.at<cv::Vec3b>(row, col) = img.at<cv::Vec3b>(temp_y, temp_x);
							break;
						}
					}
				}
			}
		}
	}
	return img;
}