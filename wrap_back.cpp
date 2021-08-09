#include"wrap_back.h"

Mesh init_mesh(int rows, int cols)
{
	MeshPoints meshPoints;
	Mesh mesh;

	rows /= 2;
	cols /= 2;
	//(M_w - 1) / (M_h - 1) = w / h, M_w * M_h = fixVNum，解方程组
	double nw, nh;
	double daita = pow(rows - cols, 2) + 4 * fixVNum * rows * cols;
	nw = ((rows - cols) + sqrt(daita)) / (2 * rows);
	nh = ((cols - rows) + sqrt(daita)) / (2 * cols);

	//cout << nw << " " << nh << endl;
	int mesh_rows = int(nh);
	int mesh_cols = int(nw);

	double rowPerMesh = double(rows - 1) / (mesh_rows - 1);
	double colPerMesh = double(cols - 1) / (mesh_cols - 1);

	for (int row_mesh = 0; row_mesh < mesh_rows; row_mesh++) {
		vector<meshPos> meshrow;
		for (int col_mesh = 0; col_mesh < mesh_cols; col_mesh++) {
			meshPos coord;
			coord.row = row_mesh * rowPerMesh;
			coord.col = col_mesh * colPerMesh;
			meshrow.push_back(coord);
		}
		meshPoints.push_back(meshrow);
	}

	mesh.meshPoints = meshPoints;
	mesh.mesh_rows = mesh_rows;
	mesh.mesh_cols = mesh_cols;
	mesh.rowPerMesh = rowPerMesh;
	mesh.colPerMesh = colPerMesh;

	return mesh;
}

void wrap_mesh(Mesh& mesh, DisMap displacementMap, Mat mask)
{
	int meshnum_row = mesh.mesh_rows;
	int meshnum_col = mesh.mesh_cols;

	for (int row_mesh = 0; row_mesh < meshnum_row; row_mesh++) {
		for (int col_mesh = 0; col_mesh < meshnum_col; col_mesh++) {
			if (row_mesh == meshnum_row - 1 && col_mesh == meshnum_col - 1) {
				meshPos& meshVertexCoord = mesh.meshPoints[row_mesh][col_mesh];
				position vertexDisplacement = displacementMap[floor(meshVertexCoord.row) - 1][floor(meshVertexCoord.col) - 1];
				meshVertexCoord.row += vertexDisplacement.row;
				meshVertexCoord.col += vertexDisplacement.col;
				meshVertexCoord.row *= 2;
				meshVertexCoord.col *= 2;
				fixMesh(mask, meshVertexCoord);
			}
			else
			{
				meshPos& meshVertexCoord = mesh.meshPoints[row_mesh][col_mesh];
				position vertexDisplacement = displacementMap[floor(meshVertexCoord.row)][floor(meshVertexCoord.col)];
				meshVertexCoord.row += vertexDisplacement.row;
				meshVertexCoord.col += vertexDisplacement.col;
				meshVertexCoord.row *= 2;
				meshVertexCoord.col *= 2;
				fixMesh(mask, meshVertexCoord);
			}
		}
	}
}

void fixMesh(Mat mask, meshPos& meshVertexCoord)
{
	int row1 = meshVertexCoord.row;
	int col1 = meshVertexCoord.col;
	int row2 = meshVertexCoord.row;
	int col2 = meshVertexCoord.col;

	int rows = mask.rows;
	int cols = mask.cols;

	while (meshVertexCoord.row < rows / 2 && row1 < rows && mask.at<uchar>(row1, col1) == 255) row1++;
	while (meshVertexCoord.row > rows / 2 && row1 >= 0 && mask.at<uchar>(row1, col1) == 255) row1--;
	while (meshVertexCoord.col < cols / 2 && col2 < cols && mask.at<uchar>(row2, col2) == 255) col2++;
	while (meshVertexCoord.col > cols / 2 && col2 >= 0 && mask.at<uchar>(row2, col2) == 255) col2--;

	if (abs(meshVertexCoord.row - row1) < abs(meshVertexCoord.col - col2)) meshVertexCoord.row = row1;
	else meshVertexCoord.col = col2;
}