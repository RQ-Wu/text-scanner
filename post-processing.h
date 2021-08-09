#pragma once
#include"global.h"

void get_scale(Mesh mesh, Mesh meshUpdated, double& x_scale, double& y_scale);
void get_min_max(MeshPoints meshPoints, int row, int col, int& x_min, int& x_max, int& y_min, int& y_max);
void stretching_reduction(Mesh mesh, Mesh& meshUpdated, double& scale_x, double& scale_y);