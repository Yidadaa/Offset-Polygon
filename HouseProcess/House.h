#include <string>
#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

#ifndef HOUSE_H
#define HOUSE_H

namespace HouseProcess {
    class Segment;
    class Point;
    class Region;
    class House;

    double computeTriArea(Point a, Point b, Point c);

    /* 点类 */
    class Point {
    public:
        double x;
        double y;
        double z;
        string id;
        double bulge;
        bool isNULL;
        Point();
        Point(double x, double y, double bulge, string id);
        Point(double x_val, double y_val);
        bool isEqualTo(Point p);
        bool isInRegion(Region r);
    };

    /* 定义一条线段 */
    class Segment {
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
        Point startPoint;
        Point endPoint;
        Point center;
        bool isNULL;
        double distance; // 直线距离，而非曲线距离
        Segment();
        Segment(Point sp, Point ep, string id);
        Segment(Point sp, Point ep);
        Point getCorWith(Segment s);
        vector<Point> getCorWithRegion(Region r);
        bool isParalWith(Segment s);
    };

    class Region {
    public:
        Region();
        Region(vector<Segment> s);
        vector<Segment> borders; // 边界集合
        bool isNULL;
        double area;
        double perimeter;
        Point center;
        Point findCenter();
        double computeArea();
        double computePerimeter();
    };

    class House {
    public:
        vector<Region> regions;
        House();
        House(vector<Segment> lines);
        vector<Region> findRegions();
        vector<Segment> findOutLines();
        vector<Segment> lines;
        vector<Segment> outLines; // 外延轮廓线
        bool isNULL;
    };
}
#endif