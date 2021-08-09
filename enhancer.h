#pragma once
#include "global.h"

Mat sharpen(Mat input, double add_weight);
Mat remove_bgd(Mat input);
Mat remove_shadow(Mat src, int n);