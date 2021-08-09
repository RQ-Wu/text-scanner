#pragma once
#include"global.h"

void draw_lsd(double out[], int n_out, Mat img);
void draw_line_in_qua(Mat img, LineInQua lineInQua, Mesh mesh);
void draw_mesh(Mat img, Mesh mesh, string winName);
void display_local_warp(Mat img);