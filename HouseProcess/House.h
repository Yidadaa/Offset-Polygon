#include <string>
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>


using namespace std;

#ifndef HOUSE_H
#define HOUSE_H

namespace HouseProcess {
    class YFSegment;
    class YFPoint;
    class YFRegion;
    class YFHouse;

    /* 点类 */
    class YFPoint {
    public:
        double x;
        double y;
        double z;
        string id;
        double bulge;
        bool isNULL;
        YFPoint();
        YFPoint(double x, double y, double bulge, string id);
        YFPoint(double x_val, double y_val);
        bool isEqualTo(YFPoint p);
        bool isInRegion(YFRegion r); // 如果点在区域边界上，那么将其算作在区域内
        bool isInRegionWithoutBorder(YFRegion r); // 如果点在区域边界上，那么不将其算作在区域内
    };

    /* 定义一条线段 */
    class YFSegment {
    public:
        // a, b, c 为直线的一般形式的三个参数
        struct Range {
            double min;
            double max;
        };
        double a, b, c;
        string id;
        Range xRange;
        Range yRange;
        YFPoint startPoint;
        YFPoint endPoint;
        YFPoint center;
        bool isNULL;
        double distance; // 直线距离，而非曲线距离
        YFSegment();
        YFSegment(YFPoint sp, YFPoint ep, string id);
        YFSegment(YFPoint sp, YFPoint ep);
        YFPoint getCorWith(YFSegment s);
        vector<YFPoint> getCorWithRegion(YFRegion r);
        bool isParalWith(YFSegment s);
        bool hasPoint(YFPoint p); // 判断某点是否在线段上
    };

    class YFRegion {
    public:
        YFRegion();
        YFRegion(vector<YFSegment> s);
        vector<YFSegment> borders; // 边界集合
        bool isNULL;
        double area;
        double perimeter;
        YFPoint center;
        YFPoint findCenter();
        double computeArea();
        double computePerimeter();
        vector<YFPoint> getCorWithRegion(YFRegion r); // 计算两个区域的交点
    };

    class YFHouse {
    public:
        vector<YFRegion> regions;
        YFHouse();
        YFHouse(vector<YFSegment> line, double outWallThickness);
        double outWallThickness;
        vector<YFRegion> findRegions(vector<YFSegment> lines);
        vector<YFSegment> findOutLines();
        vector<YFSegment> findInnerLiners();
        vector<YFSegment> lines;
        vector<YFSegment> outLines; // 外延轮廓线
        vector<YFSegment> innerLines; // 中墙线
        bool isNULL;
    };

    double computeTriArea(YFPoint a, YFPoint b, YFPoint c);
}
#endif