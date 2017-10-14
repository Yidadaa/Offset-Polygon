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

Segment::Segment(Point sp, Point ep, string id_val) {
    startPoint = sp;
    endPoint = ep;
    id = id_val;
    a = sp.y - ep.y;
    b = sp.x - ep.x;
    c = sp.x * ep.y - sp.y * ep.x;
    distance = sqrt(pow(sp.x - ep.x, 2) + pow(sp.y - ep.y, 2)); // 计算长度
    xRange.min = min(sp.x, ep.x);
    xRange.max = max(sp.x, ep.x);
    yRange.min = min(sp.y, ep.y);
    yRange.max = max(sp.y, ep.y);
}

Segment::Segment(Point sp, Point ep) {
    startPoint = sp;
    endPoint = ep;
    id = "No ID";
    a = sp.y - ep.y;
    b = sp.x - ep.x;
    c = sp.x * ep.y - sp.y * ep.x;
    distance = sqrt(pow(sp.x - ep.x, 2) + pow(sp.y - ep.y, 2)); // 计算长度
    xRange.min = min(sp.x, ep.x);
    xRange.max = max(sp.x, ep.x);
    yRange.min = min(sp.y, ep.y);
    yRange.max = max(sp.y, ep.y);
}


/* 判断是否与另一条线段平行 */
bool Segment::isParalWith(Segment s) {
    return abs(a * s.b - b * s.a) < MIN_ERR;
}

/* 计算与另一条线段的交点 */
Point Segment::getCorWith(Segment s) {
    auto isInRange = [](float n, Range range) { // 用于判断某个值是否在范围内
        return n >= range.min && n <= range.max;
    };
    if (this->isParalWith(s)) {
        return Point(); // 如果是平行的，就不存在交点
    }
    float der = a * s.b - b * s.a;
    float x = (b * s.c - c * s.b) / der;
    float y = (a * s.c - c * s.a) / der; // 计算出交点坐标值
    if (isInRange(x, xRange)
        && isInRange(y, yRange)
        && isInRange(x, s.xRange)
        && isInRange(y, s.yRange)) {
        Point p(x, y);
        return p;
    }
    else {
        return Point(); // 如果交点不在线段范围内，也不作数
    }
}

/* 计算线段与区域的交点 */
vector<Point> Segment::getCorWithRegion(Region r) {
    //return vector<Point>();
}