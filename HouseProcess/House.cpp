#pragma once
#include "House.h"
//#include "stdafx.h"

#define MIN_ERR 0.00001 // 定义最小误差，用于相等计算
/* 计算三角形面积 */
double computeTriArea(Point a, Point b, Point c) {
    const double l1 = Segment(a, b).distance;
    const double l2 = Segment(a, c).distance;
    const double l3 = Segment(c, b).distance;
    const double p = (l1 + l2 + l3) / 2; // 半周长
    return sqrt(p * (p - l1) * (p - l2) * (p - l3));
};

House::House() {
    isNULL = true;
}
House::House(vector<Segment> lines) {
    isNULL = false;
    this->lines = lines;
    regions = this->findRegions();
}

/* 寻找所有的闭合区域 */
vector<Region> House::findRegions() {
    vector<Region> regions; // 用于存放所有的区域
    vector<Segment> tmpLines = this->lines; // 用于存放所有的线

                                            // 将线段的首尾颠倒
    auto reverseSeg = [](Segment s) {
        Point tmp = s.startPoint;
        s.startPoint = s.endPoint;
        s.endPoint = tmp;
        s.startPoint.bulge = -s.startPoint.bulge;
        s.endPoint.bulge = -s.endPoint.bulge;
        return s;
    };

    /* 删除集合中某个墙壁 */
    auto delWall = [](vector<Segment> lines, Segment seg) {
        for (auto i = lines.begin(); i != lines.end(); i++) {
            if ((i->startPoint.isEqualTo(seg.startPoint)
                && i->endPoint.isEqualTo(seg.endPoint)) || i->id == seg.id) {
                lines.erase(i); // 删除与seg相同id的元素
                break;
            }
            //i->startPoint.isEqualTo(seg.startPoint)
            //    && i->endPoint.isEqualTo(seg.endPoint)
        }
        return lines; // 返回删除后的集合
    };

    /* 在集合中寻找与某线段连接的线段 */
    auto findNextWall = [](vector<Segment> lines, Segment seg) {
        for each (Segment s in lines) {
            if (s.startPoint.isEqualTo(seg.endPoint)
                || s.endPoint.isEqualTo(seg.endPoint)
                || s.startPoint.isEqualTo(seg.startPoint)
                || s.endPoint.isEqualTo(seg.startPoint)) {
                return s;
            }
        }
        return Segment();
    };

    while (tmpLines.size() > 0) {
        Segment curWall = tmpLines.at(0);
        vector<Segment> borders; // 用来存放墙壁
        borders.push_back(curWall);
        tmpLines = delWall(tmpLines, curWall);
        while (tmpLines.size() > 0) { // 不断地从所有墙壁中找出首尾相连的墙
            Segment nextWall = findNextWall(tmpLines, curWall);
            if (!nextWall.isNULL) {
                // 找到的线段可能与现在的线段有四种不同的连接情况
                // 要将他们调整成首尾相连的状态
                if (curWall.endPoint.isEqualTo(nextWall.endPoint)) {
                    // curWall:  ---->
                    // nextWall:     <----
                    nextWall = reverseSeg(nextWall);
                } else if (curWall.startPoint.isEqualTo(nextWall.startPoint)) {
                    // curWall:  ---->
                    // nextWall: ---->
                    curWall = reverseSeg(curWall);
                    borders.pop_back();
                    borders.push_back(curWall);
                    // 把当前墙壁反转后重新入栈，一般来说，只会在迭代到第二次时执行到这条
                } else if (curWall.startPoint.isEqualTo(nextWall.endPoint)) {
                    // curWall:  ---->
                    // nextWall: <----
                    curWall = reverseSeg(curWall);
                    borders.pop_back();
                    borders.push_back(curWall);
                    nextWall = reverseSeg(nextWall);
                }
            } else {
                // 找不到下一个相邻墙壁了，也结束
                // TIPs: 可优化，将多余墙壁删除，但有风险，如果区域不是闭环的话，
                //       区域将无法被找到。
                break;
            }
            borders.push_back(nextWall); // 将下一个墙壁推入栈
            int lastSize = tmpLines.size();
            tmpLines = delWall(tmpLines, nextWall); // 将已推入栈的墙壁删除
            if (tmpLines.size() == lastSize) break; // 删除失败，跳出
            curWall = nextWall;
            if (nextWall.endPoint.isEqualTo(borders.at(0).startPoint)) {
                // 形成闭环，结束
                // break;
                // Tips: 这里不跳出，反而会有更好的效果，特别是当数据有问题时。
            }
        }
        if (borders.size() == 0) {
            break; // 没有收集到区域，停止循环
        } else {
            regions.push_back(Region(borders)); // 将收集到的区域入栈
            borders.clear();
        }
    }
    return regions;
}

Region::Region() {
    isNULL = true;
}

Region::Region(vector<Segment> s) {
    borders = s;
    center = this->findCenter(); // 需要更改，不能将私有属性暴露出来
    area = this->computeArea();
    perimeter = this->computePerimeter();
    isNULL = false;
}

/* 查找视觉中心位 */
/* 查找策略，寻找最长切分线，切分线不能与边界有交点，而且中点在区域内
* 取切分线中点作为视觉中心位
*/
Point Region::findCenter() {
    vector<Segment> inLines;
    int borderNum = this->borders.size();
    for (int i = 0; i < int(borderNum / 2); i++) { // 开始遍历所有切分线
        for (int j = i + 1; j < borderNum; j++) {
            Segment s = Segment(
                this->borders.at(i).startPoint,
                this->borders.at(j).startPoint
                );
            // 判断这条切分线是否在边线上
            // TIPs: 这里使算法复杂度上升到了o(n! * n)
            bool isInBorder = false;
            for each (Segment seg in this->borders) {
                if (seg.isParalWith(s)) {
                    isInBorder = true;
                    break;
                };
            }
            if (isInBorder) continue;

            // 判断这条线是否与边线相交
            vector<Point> corPoints = s.getCorWithRegion(*this);
            if (corPoints.size() > 0) continue;

            inLines.push_back(s); // 保存该线
        }
    }
    if (inLines.size() == 0) {
        // TODO: 实在找不到，就用所有线段的加权中心点
        double min_cx = 10000;
        double max_cx = 0;
        double min_cy = 10000;
        double max_cy = 0;
        for each (auto s in this->borders) {
            double x = s.startPoint.x;
            double y = s.startPoint.y;
            min_cx = x < min_cx ? x : min_cx;
            max_cx = x > max_cx ? x : max_cx;
            min_cy = y < min_cy ? y : min_cy;
            max_cy = y > max_cy ? y : max_cy;
        }
        return Point((min_cx + max_cx) / 2, (min_cy + max_cy) / 2); // 没有符合条件的线
    }

    // 开始寻找最佳切分点
    Point bestPoint;
    double maxRatio = 0;
    for each (Segment seg in inLines) {
        // 计算线段横跨矩形的面积
        double l = abs(seg.xRange.max - seg.xRange.min); // 长
        double w = abs(seg.yRange.max - seg.yRange.min); // 宽
        double ratio = l * w;
        Point center = seg.center; // 选取切分点的中点作为最佳视觉中心点
        bool isInRegion = center.isInRegion(*(this));
        if (ratio > maxRatio && isInRegion) {
            maxRatio = ratio;
            bestPoint = center;
        }
    }
    return bestPoint;
}

/* 计算区域面积 */
double Region::computeArea() {
    vector<Point> points; // 区域的所有角点
    double area = 0;

    /* 从点集中删除点 */
    auto delPointFromPoints = [](Point p, vector<Point> points) {
        for (auto i = points.begin(); i != points.end(); i++) {
            if (p.isEqualTo(*i)) {
                points.erase(i);
                break;
            }
        }
        return points;
    };

    double arcArea = 0; // 先计算带有弧边的面积

    for each (Segment s in this->borders) {
        points.push_back(s.startPoint);
        double b = abs(s.startPoint.bulge);
        double p = s.startPoint.bulge > 0 ? 1 : -1; // 区分凸出来还是凹进去
        if (b > MIN_ERR) {
            double alpha = 2 * atan(b); // 二分之一角度
            double a = s.distance / 2;
            double R = a / sin(alpha);
            double b = a / tan(alpha);
            double s = 0.5 * a * b; // 计算三角形面积
            double arc = 0.5 * alpha * pow(R, 2); // 计算扇形面积
            arcArea += p * 2 * (arc - s); // 计算弧面切边面积
        }
    } // 获取所有角点

    while (points.size() > 0) { // 不断地从多边形中选取点，切分成三角形进行消解
        int lastsize = points.size();
        for (int i = 0; i < lastsize; i++) {
            Point sp = points.at(i);
            Point cp = points.at((i + 1) % lastsize);
            Point ep = points.at((i + 2) % lastsize);
            Segment triLine = Segment(sp, ep); // 斜边
            vector<Point> corPoints = triLine.getCorWithRegion(*this);
            bool centerIsInRegion = triLine.center.isInRegion(*this);
            if (corPoints.size() == 0 && centerIsInRegion) {
                area += computeTriArea(sp, cp, ep); // 计算三角形的面积
                points = delPointFromPoints(cp, points); // 将中点删去
                break;
            }
        }
        if (points.size() == lastsize) break; // 如果没有可以选取的点了，那么跳出
    }
    return area + arcArea;
}

/* 计算周长 */
double Region::computePerimeter() {
    double perimeter = 0;
    for each (auto l in this->borders) {
        double b = abs(l.startPoint.bulge);
        if (b > MIN_ERR) {
            // 计算弧度周长
            double alpha = 2 * atan(b); // 得到二分之一角度
            double a = l.distance / 2;
            double R = a / sin(alpha);
            double arc = R * alpha; // 计算半弧长
            perimeter += 2 * arc;
        } else {
            perimeter += l.distance;
        }
    }
    return perimeter;
};

Point::Point() {
    isNULL = true;
}

Point::Point(double x_val, double y_val, double bulge_val, string id_val) {
    x = x_val;
    y = y_val;
    z = 0;
    bulge = bulge_val;
    id = id_val;
    isNULL = false;
}

Point::Point(double x_val, double y_val) {
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
    Point zeroPoint(-1, -1);
    Segment line = Segment(zeroPoint, *this); // 画一条射向区域外的射线
    vector<Point> corPoints = line.getCorWithRegion(r); // 取得射线与区域的交点
    return corPoints.size() % 2 == 1; // 如果交点个数为奇数个，则判定该点在区域内
}


Segment::Segment(Point sp, Point ep, string id_val) {
    startPoint = sp;
    endPoint = ep;
    id = id_val;
    a = sp.y - ep.y;
    b = sp.x - ep.x;
    c = sp.x * ep.y - sp.y * ep.x;
    isNULL = false;
    distance = sqrt(pow(sp.x - ep.x, 2) + pow(sp.y - ep.y, 2)); // 计算长度
    xRange.min = min(sp.x, ep.x);
    xRange.max = max(sp.x, ep.x);
    yRange.min = min(sp.y, ep.y);
    yRange.max = max(sp.y, ep.y);
    center = Point(
        (this->startPoint.x + this->endPoint.x) / 2,
        (this->startPoint.y + this->endPoint.y) / 2
        );
}

Segment::Segment(Point sp, Point ep) {
    startPoint = sp;
    endPoint = ep;
    id = "No ID";
    a = sp.y - ep.y;
    b = sp.x - ep.x;
    c = sp.x * ep.y - sp.y * ep.x;
    isNULL = false;
    distance = sqrt(pow(sp.x - ep.x, 2) + pow(sp.y - ep.y, 2)); // 计算长度
    xRange.min = min(sp.x, ep.x);
    xRange.max = max(sp.x, ep.x);
    yRange.min = min(sp.y, ep.y);
    yRange.max = max(sp.y, ep.y);
    center = Point(
        (this->startPoint.x + this->endPoint.x) / 2,
        (this->startPoint.y + this->endPoint.y) / 2
        );
}

Segment::Segment() {
isNULL: true;
}


/* 判断是否与另一条线段平行 */
bool Segment::isParalWith(Segment s) {
    return abs(a * s.b - b * s.a) < MIN_ERR;
}

/* 计算与另一条线段的交点 */
Point Segment::getCorWith(Segment s) {
    auto isInRange = [](double n, Range range) { // 用于判断某个值是否在范围内
        return n >= range.min - MIN_ERR && n <= range.max + MIN_ERR;
    };
    if (this->isParalWith(s)) {
        return Point(); // 如果是平行的，就不存在交点
    }
    double der = a * s.b - b * s.a;
    double x = (b * s.c - c * s.b) / der;
    double y = (a * s.c - c * s.a) / der; // 计算出交点坐标值
    if (isInRange(x, xRange)
        && isInRange(y, yRange)
        && isInRange(x, s.xRange)
        && isInRange(y, s.yRange)) {
        Point p = Point(x, y);
        return p;
    } else {
        return Point(); // 如果交点不在线段范围内，也不作数
    }
}

/* 计算线段与区域的交点 */
vector<Point> Segment::getCorWithRegion(Region r) {
    vector<Segment> borders = r.borders;
    vector<Point> corPoints; // 交点集合
    auto hasInSet = [](Point p, vector<Point> pset) {
        bool flag = false;
        for each (Point pi in pset) {
            flag = flag || p.isEqualTo(pi);
            if (flag) break;
        }
        return flag;
    };
    for each (Segment s in borders) {
        Point corPoint = this->getCorWith(s);
        Point(1, 1);
        if (!corPoint.isNULL // 非空
            && !corPoint.isEqualTo(this->startPoint) // 不算线段的端点
            && !corPoint.isEqualTo(this->endPoint)
            && !hasInSet(corPoint, corPoints)) { // 避免重复添加
            corPoints.push_back(corPoint); // 将交点保存入集合
        }
    }
    return corPoints;
}

