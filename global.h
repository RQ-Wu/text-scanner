#pragma once
#include<opencv2/opencv.hpp>
#include<iostream>
#include<algorithm>
#include<vector>
#include<Eigen/Dense>
#include<Eigen/Sparse>
#include<GL/glut.h>
#include"lsd.h"

using namespace cv;
using namespace std;

#define BGD 0  // 一定是背景
#define FGD 1  // 一定是前景
#define BGD_IN_RECT 2  // 矩形内用户未指定的背景
#define FGD_IN_RECT 3  // 矩形内用户未指定的前景

#define DisMap vector<vector<position>>
#define MeshPoints vector<vector<meshPos>>
#define LineInQua vector<vector<vector<lineSeg>>>

#define LEFT 0
#define TOP 1
#define RIGHT 2
#define BOTTOM 3

#define VERTICAL 0
#define HORIZONTAL 1

#define INF 100000000
#define PI 3.14159265358979323846

#define NORMAL 0
#define DISPLAY 1

typedef Eigen::SparseMatrix<double> SparseMatrixD;
typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SpareseMatrixD_Row;
typedef Eigen::VectorXd VectorXd;
typedef Eigen::MatrixXd MatrixXd;
typedef Eigen::Vector2d Vector2d;
typedef Eigen::Vector2i Vector2i;
typedef Eigen::MatrixXi MatrixXi;
typedef Eigen::Matrix2d Matrix2d;
typedef Eigen::Triplet<double> T;
typedef Eigen::SimplicialCholesky<SparseMatrixD> CSolve;

const int fixPixelNum = 500000;
const int fixVNum = 400;

typedef struct position
{
	int row = 0, col = 0;
}position;

typedef struct meshPos
{
	double row = 0, col = 0;
}meshPos, linePos, doublePos;

typedef struct Mesh
{
	MeshPoints meshPoints;
	int mesh_rows;
	int mesh_cols;
	double rowPerMesh;
	double colPerMesh;
}Mesh;

typedef struct lineSeg
{
	linePos begin;
	linePos end;
}lineSeg;

