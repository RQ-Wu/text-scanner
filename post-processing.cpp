#include"post-processing.h"

void get_scale(Mesh mesh, Mesh meshUpdated, double& x_scale, double& y_scale)
{
	int qua_rows = mesh.mesh_rows - 1;
	int qua_cols = mesh.mesh_cols - 1;

	MeshPoints meshPoints = mesh.meshPoints;
	MeshPoints meshUpdatedPoints = meshUpdated.meshPoints;

	int N = qua_rows * qua_cols;
	double sx = 0, sy = 0;
	for (int i = 0; i < qua_rows; i++)
	{
		for (int j = 0; j < qua_cols; j++)
		{
			int x_min1, x_max1, y_min1, y_max1;
			int x_min2, x_max2, y_min2, y_max2;

			get_min_max(meshPoints, i, j, x_min1, x_max1, y_min1, y_max1);
			get_min_max(meshUpdatedPoints, i, j, x_min2, x_max2, y_min2, y_max2);

			double y_temp_scale = ((double)y_max2 - (double)y_min2) / ((double)y_max1 - (double)y_min1);
			double x_temp_scale = ((double)x_max2 - (double)x_min2) / ((double)x_max1 - (double)x_min1);

			if (y_max1 == y_min1) y_temp_scale = 0;
			if (x_max1 == x_min1) x_temp_scale = 0;

			sx += x_temp_scale;
			sy += y_temp_scale;
		}
	}

	x_scale = 1 / (sx / (double)N);
	y_scale = 1 / (sy / (double)N);
}

void get_min_max(MeshPoints meshPoints, int row, int col, int& x_min, int& x_max, int& y_min, int& y_max)
{
	meshPos leftTop = meshPoints[row][col];
	meshPos rightTop = meshPoints[row][col + 1];
	meshPos leftBottom = meshPoints[row + 1][col];
	meshPos rightBottom = meshPoints[row + 1][col + 1];

	int x[4] = { leftTop.col, rightTop.col, leftBottom.col, rightBottom.col };
	int y[4] = { leftTop.row, rightTop.row, leftBottom.row, rightBottom.row };
	sort(x, x + 4);
	sort(y, y + 4);

	x_min = x[0], x_max = x[3], y_min = y[0], y_max = y[3];
}

void stretching_reduction(Mesh mesh, Mesh& meshUpdated, double& scale_x, double& scale_y)
{
	get_scale(mesh, meshUpdated, scale_x, scale_y);

	cout << "scale is  " << scale_x << " " << scale_y << endl;
	int rows = mesh.mesh_rows;
	int cols = mesh.mesh_cols;

	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			meshUpdated.meshPoints[i][j].row *= scale_y;
			meshUpdated.meshPoints[i][j].col *= scale_x;
		}
	}
}