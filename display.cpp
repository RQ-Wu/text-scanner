#include"display.h"

void draw_lsd(double out[], int n_out, Mat img)
{
	for (int i = 0; i < n_out; i++)
	{
		int start_x = out[i * 7];
		int start_y = out[i * 7 + 1];
		int end_x = out[i * 7 + 2];
		int end_y = out[i * 7 + 3];

		Point start(start_x, start_y);
		Point end(end_x, end_y);

		line(img, start, end, CV_RGB(0, 255, 0));
	}
	imshow("lsd", img);
	waitKey(0);
}

void draw_line_in_qua(Mat img, LineInQua lineInQua, Mesh mesh)
{
	int qua_cols = mesh.mesh_cols - 1;
	int qua_rows = mesh.mesh_rows - 1;
	MeshPoints meshPoints = mesh.meshPoints;

	Mat meshImg, lineImg;
	img.copyTo(meshImg);
	img.copyTo(lineImg);

	for (int i = 0; i < qua_rows; i++)
	{
		for (int j = 0; j < qua_cols; j++)
		{
			meshPos leftTop = meshPoints[i][j];
			meshPos rightTop = meshPoints[i][j + 1];
			meshPos leftBottom = meshPoints[i + 1][j];
			meshPos rightBottom = meshPoints[i + 1][j + 1];
			//cout << i << " " << rightTop.col << endl;

			// draw qua
			line(meshImg, Point(leftTop.col, leftTop.row), Point(rightTop.col, rightTop.row), CV_RGB(0, 255, 0), 2);
			line(meshImg, Point(leftTop.col, leftTop.row), Point(leftBottom.col, leftBottom.row), CV_RGB(0, 255, 0), 2);
			line(meshImg, Point(rightBottom.col, rightBottom.row), Point(rightTop.col, rightTop.row), CV_RGB(0, 255, 0), 2);
			line(meshImg, Point(leftBottom.col, leftBottom.row), Point(rightBottom.col, rightBottom.row), CV_RGB(0, 255, 0), 2);

			// draw lines
			vector<lineSeg> lineInThisQua = lineInQua[i][j];
			for (int k = 0; k < lineInThisQua.size(); k++)
			{
				lineSeg curLine = lineInThisQua[k];
				line(meshImg, Point(curLine.begin.col, curLine.begin.row), Point(curLine.end.col, curLine.end.row), CV_RGB(255, 0, 0), 2);
			}

			imshow("line-mesh", meshImg);
			//imshow("line", lineImg);
			if (i != qua_rows - 1 || j != qua_cols - 1)
			{
				waitKey(10);
			}
			else
			{
				waitKey(10);
			}
		}
	}
}

void draw_mesh(Mat img, Mesh mesh, string winName)
{
	int qua_cols = mesh.mesh_cols - 1;
	int qua_rows = mesh.mesh_rows - 1;
	MeshPoints meshPoints = mesh.meshPoints;

	Mat meshImg, lineImg;
	img.copyTo(meshImg);
	img.copyTo(lineImg);

	for (int i = 0; i < qua_rows; i++)
	{
		for (int j = 0; j < qua_cols; j++)
		{
			meshPos leftTop = meshPoints[i][j];
			meshPos rightTop = meshPoints[i][j + 1];
			meshPos leftBottom = meshPoints[i + 1][j];
			meshPos rightBottom = meshPoints[i + 1][j + 1];
			// cout << i << " " << rightTop.col << endl;

			// draw qua
			line(meshImg, Point(leftTop.col, leftTop.row), Point(rightTop.col, rightTop.row), CV_RGB(0, 255, 0), 1);
			line(meshImg, Point(leftTop.col, leftTop.row), Point(leftBottom.col, leftBottom.row), CV_RGB(0, 255, 0), 1);
			line(meshImg, Point(rightBottom.col, rightBottom.row), Point(rightTop.col, rightTop.row), CV_RGB(0, 255, 0), 1);
			line(meshImg, Point(leftBottom.col, leftBottom.row), Point(rightBottom.col, rightBottom.row), CV_RGB(0, 255, 0), 1);

			imshow(winName, meshImg);
			//imshow("line", lineImg);
			if (i != qua_rows - 1 || j != qua_cols - 1)
			{
				waitKey(10);
			}
			else
			{
				waitKey(10);
			}
		}
	}
}

void display_local_warp(Mat img)
{
	int rows = img.rows;
	int cols = img.cols;

	Mat disImg;
	resize(img, disImg, Size(cols * 2, rows * 2), 0, 0, INTER_NEAREST);

	imshow("img", disImg);
	
	waitKey(5);
}