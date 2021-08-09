#pragma once
#include"global.h"

Mesh global_wrap(const Mat src, Mat mask, Mesh mesh, int mode);
SpareseMatrixD_Row getShapePrepared(Mesh mesh);
SpareseMatrixD_Row init_q(Mesh mesh);
pair<SpareseMatrixD_Row, VectorXd> get_boundary_mat(const Mat src, Mesh mesh);
LineInQua init_lines(Mat& mask, Mat img, Mesh mesh, vector<pair<int, double>>& idx_theta, vector<double>& rotate, int mode);
void clean_line_on_border_mask(Mat& mask);
vector<lineSeg> get_lines(Mat img, Mat mask);
bool inside_mask(Mat mask, double start_x, double start_y, double end_x, double end_y);
LineInQua allocate_line_to_qua(vector<lineSeg> lines, Mesh mesh);
bool isPointInQua(meshPos leftTop, meshPos rightTop, meshPos leftBottom, meshPos rightBottom, linePos point);
vector<linePos> getIntersectionPoints(meshPos leftTop, meshPos rightTop, meshPos leftBottom, meshPos rightBottom, lineSeg line);
bool getIntersection(meshPos a, meshPos b, lineSeg line, linePos& p);
bool isBetween(meshPos a, meshPos b, meshPos target);
SpareseMatrixD_Row get_line_mat(Mat img, Mat mask, Mesh mesh, vector<double>rotate_theta, LineInQua lineInQua,
	vector<pair<MatrixXd, MatrixXd>>& BilinearVec, vector<bool>& badpoints);
MatrixXd bilinear_interpolation(linePos linePoint, int row, int col, Mesh mesh);
Mesh updateMesh(VectorXd x, Mesh& mesh);
SpareseMatrixD_Row row_stack(SpareseMatrixD_Row origin, SpareseMatrixD_Row diag);
MatrixXd row_stack(MatrixXd mat1, MatrixXd mat2);
SpareseMatrixD_Row block_diag(SpareseMatrixD_Row origin, MatrixXd addin, int QuadID, Mesh mesh);
Mat fill_missing_pixel(Mat& img, const Mat mask);