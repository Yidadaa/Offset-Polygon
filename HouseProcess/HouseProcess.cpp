// HouseProcess.cpp : 定义控制台应用程序的入口点。
//

#include "House.h"
#include "House.cpp"

using namespace HouseProcess;

void test();
void isPassed(string name, bool t);

int main()
{
    test();
    system("pause");
    return 0;
}
/* 该函数用于测试 */
void test() {
    /*测试Point类*/
    YFPoint p = YFPoint(1, 1, 1, "hello");
    YFPoint p2 = YFPoint(1.0000001, 1.0000001, 1, "fuck");
    YFPoint emptyP = YFPoint();
    bool a = p.isEqualTo(p2) == true;
    isPassed("Point::isEqualTo()", a);
    isPassed("Point::Point()", emptyP.isNULL);
    /*测试Segment类*/
    YFSegment s1 = YFSegment(YFPoint(1, 1), YFPoint(2, 2));
    YFSegment s2 = YFSegment(YFPoint(0, 0), YFPoint(1, 1));
    YFSegment s3 = YFSegment(YFPoint(0, 1), YFPoint(1, 0));
    isPassed("Segment::isParalWith()", s1.isParalWith(s2));
    isPassed("Segment::getCorWith()", s1.getCorWith(s2).isNULL);
    isPassed("Segment::getCorWith()", s2.getCorWith(s3).isEqualTo(YFPoint(0.5, 0.5)));
    vector<YFSegment> segs = {
        YFSegment(YFPoint(0.4, 0.4), YFPoint(0, 1), "line_1"),
        YFSegment(YFPoint(0, 1), YFPoint(1, 0), "line_2"),
        YFSegment(YFPoint(1, 0), YFPoint(0.4, 0.4), "line_3")
    };
    YFRegion r = YFRegion(segs);
    isPassed("Point::isInReigon", YFPoint(0.5, 0.5).isInRegion(r));

    /* 开始测试House类的findRegions方法
     * 分为两种情况，一种是闭合，一种是不闭合
     * 每种情况又分线段的四种排列方式
     */
    vector<YFSegment> segs_2 = {
        YFSegment(YFPoint(0, 1), YFPoint(1, 0), "line_2"),
        YFSegment(YFPoint(0, 0), YFPoint(0, 1), "line_1"),
        YFSegment(YFPoint(0.5, 0), YFPoint(0, 0), "line_3"),
        YFSegment(YFPoint(1, 0), YFPoint(0.5, 0), "line_4")
    };
    YFHouse h = YFHouse(segs_2, 0.35);
    bool hasPassed = false;
    vector<string> conds = { "line_2", "line_1", "line_3", "line_4" };
    hasPassed = hasPassed || h.regions.size() == 1; // 是否只找到了一个区域
    for (int i = 0; i < h.regions.at(0).borders.size(); i++) {
        hasPassed = hasPassed && h.regions.at(0).borders.at(i).id == conds.at(i);
    }
    YFPoint center = h.regions.at(0).center;
    double area = h.regions.at(0).computeArea();
    isPassed("Regin::computeArea()", abs(area - 0.5) < MIN_ERR);
    isPassed("House::findRegions()", hasPassed);
    isPassed("Region::findCenter()", center.isEqualTo(YFPoint(0.25, 0.5)));

    double triArea = computeTriArea(
        YFPoint(0, 0), YFPoint(0, 1), YFPoint(1, 0)
    );

    double nu = 0.1;
    isPassed("computeTriArea()", true);
}

/* 用于输出测试信息 */
void isPassed(string name, bool t) {
    cout << (t ? "passed" : "wrong") << '\t' << name << '\n';
}