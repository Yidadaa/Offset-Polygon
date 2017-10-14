#pragma once
#include <string>
#include <cmath>
#include <iostream>
#include <vector>
#include <minmax.h>
//#include "House.h"

#define MIN_ERR 0.00001 // 定义最小误差，用于相等计算
using std::string;
using std::vector;
using namespace std;

Point::Point()
{
    isNULL = true;
}

Point::Point(float x_val, float y_val, float bulge_val, string id_val) {
    x = x_val;
    y = y_val;
    z = 0;
    bulge = bulge_val;
    id = id_val;
    isNULL = false;
}

Point::Point(float x_val, float y_val) {
    x = x_val;
    y = y_val;
    z = 0;
    bulge = 0;
    id = "No ID";
    isNULL = false;
}

/* 判断该点与另外一个点是否近似相等 */
bool Point::isEqualTo(Point p) {
    return sqrt(pow(x - p.x, 2) + pow(y - p.y, 2)) < MIN_ERR;
}

/* 判断某点是否在某区域内 */
bool Point::isInRegion(Region r) {
    Point zeroPoint(0, 0);
    //Segment line = Segment(zeroPoint, *this);
    // TODO: 完成Segment类的getCurWithRegion函数
    return true;
}

