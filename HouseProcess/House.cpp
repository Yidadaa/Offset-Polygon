#pragma once
#include "House.h"
//#include "stdafx.h"

#define MIN_ERR 0.000001 // 定义最小误差，用于相等计算
#define PI 3.14159265358 // 预定义π值

namespace HouseProcess {

    /* 计算三角形面积 */
    double computeTriArea(YFPoint a, YFPoint b, YFPoint c) {
        const double l1 = YFSegment(a, b).distance;
        const double l2 = YFSegment(a, c).distance;
        const double l3 = YFSegment(c, b).distance;
        const double p = (l1 + l2 + l3) / 2; // 半周长
        return sqrt(p * (p - l1) * (p - l2) * (p - l3));
    };

    YFHouse::YFHouse() {
        this->isNULL = true;
    }
    YFHouse::YFHouse(vector<YFSegment> lines, double thickness) {
        this->isNULL = false;
        this->lines = lines;
        this->outWallThickness = thickness;
        this->regions = this->findRegions(lines);
        this->outLines = this->findOutLines();
        this->innerLines = this->findInnerLiners();
    }


    /* 寻找所有的闭合区域 */
    vector<YFRegion> YFHouse::findRegions(vector<YFSegment> lines) {
        vector<YFRegion> regions; // 用于存放所有的区域
        vector<YFSegment> tmpLines; // 用于存放所有的线

                                    // 将线段的首尾颠倒
        auto reverseSeg = [](YFSegment s) {
            auto negNum = [](double n) { // 取相反数
                if (abs(n) < MIN_ERR) return n;
                else {
                    return -n;
                }
            };
            YFPoint tmp = s.startPoint;
            s.startPoint = s.endPoint;
            s.endPoint = tmp;
            double tmpBulge = s.startPoint.bulge;
            s.startPoint.bulge = negNum(s.endPoint.bulge); // 首尾反转，弧度交换并且取反
            s.endPoint.bulge = negNum(tmpBulge);
            return s;
        };

        /* 删除集合中某个墙壁 */
        auto delWall = [](vector<YFSegment> lines, YFSegment seg) {
            for (auto i = lines.begin(); i != lines.end(); i++) {
                if ((i->startPoint.isEqualTo(seg.startPoint) && i->endPoint.isEqualTo(seg.endPoint))
                    || (i->endPoint.isEqualTo(seg.startPoint) && i->startPoint.isEqualTo(seg.endPoint))) {
                    lines.erase(i);
                    break;
                }
            }
            return lines; // 返回删除后的集合
        };

        /* 在集合中寻找与某线段连接的线段 */
        auto findNextWall = [](vector<YFSegment> lines, YFSegment seg) {
            for (YFSegment s : lines) {
                if (s.startPoint.isEqualTo(seg.endPoint)
                    || s.endPoint.isEqualTo(seg.endPoint)
                    || s.startPoint.isEqualTo(seg.startPoint)
                    || s.endPoint.isEqualTo(seg.startPoint)) {
                    return s;
                }
            }
            return YFSegment();
        };
        //初始化tmpLines，将所有孤立墙壁剔除
        for (YFSegment seg : lines) {
            int sFlag = 0;
            int eFlag = 0;
            for (YFSegment s : lines) {
                if (s.startPoint.isEqualTo(seg.startPoint)
                    || s.endPoint.isEqualTo(seg.startPoint)) {
                    sFlag++;
                } else if (s.startPoint.isEqualTo(seg.endPoint)
                    || s.endPoint.isEqualTo(seg.endPoint)) {
                    eFlag++;
                }
            }
            // 孤立墙壁必有一个点与其他点都不重叠
            if (sFlag > 0 && eFlag > 0) tmpLines.push_back(seg);
        }

        while (tmpLines.size() > 0) {
            YFSegment curWall = tmpLines.at(0);
            vector<YFSegment> borders; // 用来存放墙壁
            bool isClosure = false; // 是否是闭合区域
            borders.push_back(curWall);
            tmpLines = delWall(tmpLines, curWall);
            while (tmpLines.size() > 0) { // 不断地从所有墙壁中找出首尾相连的墙
                YFSegment nextWall = findNextWall(tmpLines, curWall);
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
                    break;
                }
                borders.push_back(nextWall); // 将下一个墙壁推入栈
                int lastSize = tmpLines.size();
                tmpLines = delWall(tmpLines, nextWall); // 将已推入栈的墙壁删除
                if (tmpLines.size() == lastSize) break; // 删除失败，跳出
                curWall = nextWall;
                if (nextWall.endPoint.isEqualTo(borders.at(0).startPoint)) {
                    // 形成闭环，结束
                    // Tips: 这里不跳出，反而会有更好的效果，特别是当数据有问题时。
                    isClosure = true;
                    break;
                }
            }
            vector<YFSegment> tmpBorders;
            for (auto l : borders) if (l.distance > MIN_ERR) tmpBorders.push_back(l); // 过滤零线段
            borders = tmpBorders;
            YFRegion region(borders);
            if (borders.size() == 0) {
                break; // 没有收集到区域，停止循环
            } else if (isClosure && region.area > MIN_ERR) { // 仅收集闭合区域
                regions.push_back(region); // 将收集到的区域入栈
                borders.clear();
            }
        }
        return regions;
    }

    /* 获取外延轮廓线 */
    vector<YFSegment> YFHouse::findOutLines() {
        vector<YFSegment> outLines;
        vector<YFSegment> lines = this->lines;
        vector<YFRegion> regions = this->regions;
        const double distance = this->outWallThickness;
        const double D = 0.3; // 最小距离

        vector<YFRegion> outRegions;
        vector<YFPoint> allPoints; // 所有用来计算外线的点都放在这里

                                   /* 计算某点的极坐标表示，接收直角坐标，返回极坐标中的theta角 */
        auto computeTheta = [](double x, double y) {
            double theta = 0;
            if (abs(x) < MIN_ERR) {
                // x = 0
                if (abs(y) < MIN_ERR) {
                    // x = 0, y = 0
                    return 0.0;
                } else if (y > 0) {
                    // x = 0, y > 0
                    return PI / 2;
                } else {
                    // x = 0. y < 0
                    return -PI / 2;
                }
            } else if (x > 0) {
                // x > 0
                return atan(y / x);
            } else {
                // x < 0
                if (y < 0) {
                    // y < 0
                    return atan(y / x) - PI;
                } else {
                    // y >= 0
                    return atan(y / x) + PI;
                }
            }
        };

        /* 将多边形向内或者向外扩展一定距离，当d大于0时，向多边形外部扩张；否则向内部收缩 */
        auto zoomRegion = [computeTheta](YFRegion r, double distance) {
            double d = abs(distance);
            bool isZoomOut = distance > 0; // 是否外扩
            vector<YFPoint> outPoints; // 将所有角平分线上的点保存下来
            int borderCount = r.borders.size();
            for (int i = 0; i < borderCount; i++) {
                YFSegment curSeg = r.borders.at(i);
                YFSegment nextSeg = r.borders.at((i + 1) % borderCount);
                YFPoint sp, mp, ep; // 开始点，中间点，结束点
                if (curSeg.endPoint.isEqualTo(nextSeg.startPoint)) {
                    sp = curSeg.startPoint;
                    mp = nextSeg.startPoint;
                    ep = nextSeg.endPoint;
                } else if (curSeg.endPoint.isEqualTo(nextSeg.endPoint)) {
                    sp = curSeg.startPoint;
                    mp = curSeg.endPoint;
                    ep = nextSeg.startPoint;
                } else if (curSeg.startPoint.isEqualTo(nextSeg.startPoint)) {
                    sp = curSeg.endPoint;
                    mp = curSeg.startPoint;
                    ep = nextSeg.endPoint;
                } else if (curSeg.startPoint.isEqualTo(nextSeg.endPoint)) {
                    sp = curSeg.endPoint;
                    mp = curSeg.startPoint;
                    ep = nextSeg.startPoint;
                }
                // 开始计算角平分线上的点，此时以mp为坐标原点进行计算
                double oa[2] = { sp.x - mp.x, sp.y - mp.y };
                double ob[2] = { ep.x - mp.x, ep.y - mp.y };
                double thetaOA = computeTheta(oa[0], oa[1]); // 计算点a的极坐标角度
                double thetaOB = computeTheta(ob[0], ob[1]); // 计算点b的极坐标角度
                auto lenOfVector = [](double vector[]) { return sqrt(pow(vector[0], 2) + pow(vector[1], 2)); }; // 计算向量模长
                double alpha = acos((oa[0] * ob[0] + oa[1] * ob[1]) / (lenOfVector(oa) * lenOfVector(ob))); // oa与ob的夹角
                                                                                                            // 根据a和b的极坐标，可以计算出两者角平分线上的任意一点坐标
                double rho = d / sin(alpha / 2); // 计算点p的极坐标的r
                double testRho = 0.01 / sin(alpha / 2); // testRho用于判断点是否在区域内，如果用d，有时会无法判断
                double theta1 = (thetaOA + thetaOB) / 2; // 角平分线上的点
                double theta2 = theta1 + PI; // 角平分线延长线上的点
                YFPoint p1(rho * cos(theta1) + mp.x, rho * sin(theta1) + mp.y); // 需要加上mp的值进行复原操作
                YFPoint p2(rho * cos(theta2) + mp.x, rho * sin(theta2) + mp.y);
                YFPoint testp1(testRho * cos(theta1) + mp.x, testRho * sin(theta1) + mp.y); // 用于测试的两个点,testp1和testp2
                YFPoint testp2(testRho * cos(theta2) + mp.x, testRho * sin(theta2) + mp.y); // 这两个点代表了p1和p2与区域的位置关系
                                                                                            // 判断p1和p2谁在区域内，如果是外扩，那么就选区域外的点；如果是内缩，那么就选区域内的点
                YFPoint p = testp1.isInRegion(r) == isZoomOut ? p2 : p1; // 使用异或操作进行判断，是选p1还是p2
                p.bulge = mp.bulge; // 保持弧线一致
                outPoints.push_back(p); // 将点p收集起来
            }
            int pointCount = outPoints.size();
            vector<YFSegment> regionOutLines;
            double fullLength = 0.0; // 计算所有线段总长
            for (int i = 0; i < pointCount; i++) {
                auto ep = outPoints.at((i + 1) % pointCount);
                auto s = YFSegment(outPoints.at(i), YFPoint(ep.x, ep.y)); // 去除endPoint的bulge信息
                fullLength += s.distance;
                regionOutLines.push_back(s); // 将所有点连起来，作为外边线
            }
            // 检测是否存在自相交，如果有自相交，删除自相交的部分
            for (int i = 0; i < regionOutLines.size(); i++) {
                auto curLine = regionOutLines.at(i);
                for (int j = 0; j < regionOutLines.size(); j++) {
                    if ((i == 0 && j == regionOutLines.size() - 1) ||
                        (j == 0 && i == regionOutLines.size() - 1) ||
                        (abs(i - j) <= 1)) continue; // 跳过相邻边
                    auto theJline = regionOutLines.at(j);
                    auto corPoint = theJline.getCorWith(curLine); // 计算交点
                    if (corPoint.isNULL) continue; // 交点为空，及时跳出
                                                   // 如果存在交点，就将两交点之间的线段删除
                                                   // 通过计算周长来判断删除哪部分线段
                    int minIndex = i > j ? j : i;
                    int maxIndex = i > j ? i : j;
                    double min2maxLength = 0.0; // 计算min到max之间的线段长度
                    for (int k = 0; k < regionOutLines.size(); k++) {
                        if (k >= minIndex && k <= maxIndex) min2maxLength += regionOutLines.at(k).distance;
                    }
                    int beforeDel = regionOutLines.size();
                    if (min2maxLength < fullLength / 2.0) {
                        // 需要删除(min, max)之间的线段
                        corPoint.bulge = regionOutLines.at(maxIndex).startPoint.bulge; // 保持弧线信息一致
                        regionOutLines.at(minIndex) = YFSegment(regionOutLines.at(minIndex).startPoint, corPoint);
                        regionOutLines.at(maxIndex) = YFSegment(corPoint, regionOutLines.at(maxIndex).endPoint);
                        regionOutLines.erase(regionOutLines.begin() + minIndex + 1, regionOutLines.begin() + maxIndex); // 删除(minIndex, maxIndex)
                    } else {
                        // 需要删除[0, min), (max, N)之间的线段
                        corPoint.bulge = regionOutLines.at(minIndex).startPoint.bulge; // 保持弧线信息一致
                        regionOutLines.at(minIndex) = YFSegment(corPoint, regionOutLines.at(minIndex).endPoint);
                        regionOutLines.at(maxIndex) = YFSegment(regionOutLines.at(maxIndex).startPoint, corPoint);
                        regionOutLines.erase(regionOutLines.begin() + maxIndex + 1, regionOutLines.end()); // 删除(maxIndex, end]
                        regionOutLines.erase(regionOutLines.begin(), regionOutLines.begin() + (minIndex - 1 <= 0 ? 1 : minIndex - 1)); // 删除[0, minIndex)
                    }
                    int afterDel = regionOutLines.size();
                    // 删除线段后，重新开始循环
                    if (beforeDel != afterDel) {
                        // 这里的判断是防止死循环
                        i = 0;
                        break;
                    }
                }
            }
            vector<YFSegment> tmpLines;
            for (auto l : regionOutLines) {
                // 删除零线段
                if (l.distance > MIN_ERR) tmpLines.push_back(l);
            }
            regionOutLines = tmpLines;
            return YFRegion(regionOutLines);
        };

        for (auto r : regions) outRegions.push_back(zoomRegion(r, distance >= D ? distance : D));
        // 将所有区域外扩distance个单位，为了保证区域之间有重叠，设置最小值为D，小于D的另行处理

        /* 将n点集按坐标排序 */
        auto compare = [](YFPoint a, YFPoint b) {
            if (abs(a.x - b.x) < MIN_ERR) return a.y < b.y;
            else return a.x < b.x;
        };

        if (outRegions.size() > 1) {
            // 开始计算所有区域之间的交点
            vector<YFPoint> corPoints;
            for (int i = 0; i < outRegions.size(); i++) {
                auto curRegion = outRegions.at(i);
                for (auto l : curRegion.borders) { // 对于当前区域的每一个边界，都计算其与其他区域的交点
                    vector<YFPoint> corPts; // 用于存放当前边界与其他的交点
                    vector<YFRegion> corRegions; // 用于存放与当前边界有交点的区域
                    for (int j = 0; j < outRegions.size(); j++) {
                        if (j == i) continue;
                        auto nextRegion = outRegions.at(j);
                        auto regionCorPts = l.getCorWithRegion(nextRegion);
                        if (regionCorPts.size() > 0) {
                            corPts.insert(corPts.end(), regionCorPts.begin(), regionCorPts.end()); // 将所有交点收集起来
                            corRegions.push_back(nextRegion); // 将有交点的区域也收集起来
                        }
                    }
                    if (corPts.size() > 0) { // 存在交点时，进行合并运算
                        for (int k = 0; k < corPts.size(); k++) {
                            corPts[k] = YFPoint(corPts.at(k).x, corPts.at(k).y, l.startPoint.bulge, ""); // 保持弧线信息
                        }
                        corPts.push_back(l.startPoint);
                        corPts.push_back(l.endPoint);
                        corPoints.insert(corPoints.end(), corPts.begin(), corPts.end()); // 保存交点
                                                                                         // 对所有点进行排序，排序之后取子线段，如果线段中点在区域外，则保存线段
                        sort(corPts.begin(), corPts.end(), compare);
                        for (int k = 0; k < corPts.size() - 1; k++) { // 将排过序的点顺次连起来
                            YFSegment tmpSeg(corPts.at(k), corPts.at(k + 1));
                            // 判断连起来的线段是否在相交区域外
                            bool isInRegion = false;
                            for (int index = 0; index < outRegions.size(); index++) {
                                // 把在区域内的线去除
                                if (index == i) continue;
                                auto r = outRegions.at(index);
                                isInRegion = isInRegion || tmpSeg.center.isInRegionWithoutBorder(r);
                            }
                            for (auto r : this->regions) {
                                //continue;
                                // 把在内区域的线也去除掉
                                isInRegion = isInRegion || tmpSeg.center.isInRegion(r) || tmpSeg.startPoint.isInRegion(r) || tmpSeg.endPoint.isInRegion(r);
                            }
                            if (!isInRegion) { // 若在区域外，则判定其为边界
                                outLines.push_back(tmpSeg);
                            }
                        }
                    } else {
                        // 如果没有交点，证明这条边界不需要合并
                        bool isInRegion = false;
                        for (int index = 0; index < outRegions.size(); index++) {
                            if (index == i) continue;
                            auto r = outRegions.at(index);
                            isInRegion = isInRegion || l.center.isInRegionWithoutBorder(r);
                        }
                        for (int index = 0; index < this->regions.size(); index++) {
                            //continue;
                            auto r = this->regions.at(index);
                            isInRegion = isInRegion || l.center.isInRegion(r) || l.startPoint.isInRegion(r) || l.endPoint.isInRegion(r);
                        }
                        if (!isInRegion) outLines.push_back(l);
                    }
                }
            }
        } else if (outRegions.size() == 1) {
            // 如果只有一个区域，则不需要计算交点，直接将外部区域的边线返回即可
            outLines = outRegions.begin()->borders;
        }

        if (outLines.size() > 0) {
            // 对outLines中的线段进行去重操作
            vector<YFSegment> tmpOutLines;
            for (auto l : outLines) {
                bool alreadyExist = false;
                for (auto s : tmpOutLines) {
                    // 检测是否有相同线段已经存在了，这里复杂度为o(n^2)，可以用哈希来优化
                    alreadyExist = abs(s.center.x - l.center.x) < MIN_ERR && abs(s.center.y - l.center.y) < MIN_ERR;
                    if (alreadyExist) break;
                }
                if (!alreadyExist) {
                    tmpOutLines.push_back(l);
                }
            }
            outLines = tmpOutLines;
        }
        auto rs = YFHouse().findRegions(outLines); // 利用findRegions函数对外线进行重排序
        auto finalOutRegion = YFRegion();
        int maxSize = 0;
        if (rs.size() > 0) {
            // 如果找到多个外区域，返回最大的那个
            for (int i = 0; i < rs.size(); i++) {
                if (rs[i].borders.size() > maxSize) {
                    finalOutRegion = rs[i];
                    maxSize = rs[i].borders.size();
                }
            }
        }
        if (!finalOutRegion.isNULL && distance < D) {
            finalOutRegion = zoomRegion(finalOutRegion, distance - D);
        }
        if (!finalOutRegion.isNULL) outLines = finalOutRegion.borders;
        return outLines;
    }

    /* 计算中墙线 */
    vector<YFSegment> YFHouse::findInnerLiners() {
        struct LINE {
            YFSegment s;
            bool hasChanged;
        };
        map <string, LINE> inlinesMap; // 用来存放中线
        map <string, string> singleMap; // 用来建立对应关系，<当前线段，lineMap的key>
                                        // 建立一个哈希表，用于快速查找两条线之间的中线
        vector<YFRegion> allLines;
        // 存放所有线段，依次是 各区域边线 - 外线
        allLines = this->regions;
        allLines.push_back(YFRegion(this->outLines));
        // 每一条线在allLines中，都可以由[rIndex][lIndex]索引到
        // 这个索引将用于查找哈希表

        /* 将n点集按坐标排序 */
        auto compare = [](YFPoint a, YFPoint b) {
            if (abs(a.x - b.x) < MIN_ERR) return a.y < b.y;
            else return a.x < b.x;
        };

        /* 从小到大将点按坐标排序 */
        auto sortPoints = [compare](YFPoint a, YFPoint b) {
            vector<YFPoint> arr({ a, b });
            sort(arr.begin(), arr.end(), compare);
            return arr;
        };

        /* 计算探针线段 */
        auto computeSeekers = [sortPoints](YFSegment s, YFRegion r, bool isOutLine) {
            const double delta = 0.01; // 偏移量
            const double d = 0.5; // 探针伸出长度
            auto sorted = sortPoints(s.startPoint, s.endPoint);
            auto sp = sorted[0];
            auto ep = sorted[1];
            auto cp = s.center;
            vector<YFSegment> seekers;
            double cosAlpha = abs(ep.y - sp.y) / s.distance;
            double sinAlpha = abs(ep.x - sp.x) / s.distance;
            YFPoint testPoint(cp.x + delta * cosAlpha, cp.y + delta * sinAlpha); // 用于判定探针方向
            int direction = (testPoint.isInRegion(r) == isOutLine) ? 1 : -1;
            seekers.push_back(YFSegment(sp, YFPoint(sp.x + direction * d * cosAlpha, sp.y + direction * d * sinAlpha)));
            seekers.push_back(YFSegment(ep, YFPoint(ep.x + direction * d * cosAlpha, ep.y + direction * d * sinAlpha)));
            return seekers;
        };

        /* 生成哈希表的键值，rIndex - 所在区域的索引，lIndex - 线段在区域中的索引 */
        auto computeKeys = [](int rIndex1, int lIndex1, int rIndex2, int lIndex2) {
            vector<int> k1({ rIndex1, lIndex1 });
            vector<int> k2({ rIndex2, lIndex2 });
            vector<int> sk;
            string key;
            if (k1[0] > k2[0] || (k1[0] == k2[0] && k1[1] > k2[1])) {
                sk = k1;
                k1 = k2;
                k2 = sk;
            }
            key = to_string(k1[0]) + string("-") + to_string(k1[1]) + string("-")
                + to_string(k2[0]) + string("-") + to_string(k2[1]);
            return key;
        };

        /* 生成singleMap的key */
        auto compute2Key = [](int i, int j) {
            return to_string(i) + string("-") + to_string(j);
        };

        /* 计算两条直线的交点 */
        auto computeCorOfLines = [](YFSegment a, YFSegment b) {
            double der = a.a * b.b - a.b * b.a;
            if (a.isParalWith(b)) return YFPoint(); // a, b 平行不共线，不进行计算
            double x = (a.b * b.c - a.c * b.b) / der;
            double y = (a.a * b.c - a.c * b.a) / der; // 计算出交点坐标值
            return YFPoint(x, y);
        };

        /* 计算两条线段的中线 */
        auto computeMidLine = [sortPoints](YFSegment s, YFSegment l) {
            auto a = sortPoints(s.startPoint, s.endPoint);
            auto b = sortPoints(l.startPoint, l.endPoint);
            auto msp = YFSegment(a[0], b[0]).center;
            auto mep = YFSegment(a[1], b[1]).center;

            return YFSegment(msp, mep);
        };

        /* 将线段排序 */
        auto compareLine = [](YFSegment a, YFSegment b) {
            if (abs(a.center.x - b.center.x) < MIN_ERR) return a.center.y < b.center.y;
            else return a.center.x < b.center.x;
        };

        // 以下算法是o(n^2)的
        int count = allLines.size();
        for (int i = 0; i < count; i++) {
            bool isOutLine = i >= this->regions.size(); // 是否是外边线段
            auto r = allLines[i]; // 当前区域
            for (int j = 0; j < r.borders.size(); j++) {
                auto s = r.borders[j]; // 当前线段
                auto seekers = computeSeekers(s, r, isOutLine); // 一般是两个探针

                struct TMP_LINE {
                    YFSegment nearestLine;
                    vector<int> keys;
                    double minDistance;
                };

                vector<TMP_LINE> opLines(seekers.size(), { YFSegment(), vector<int>(), 100000 });

                for (int k = 0; k < count; k++) { // 遍历所有线段，找出对应的最近的线
                    auto rr = allLines[k];
                    for (int n = 0; n < rr.borders.size(); n++) {
                        for (int ptr = 0; ptr < seekers.size(); ptr++) {
                            // 对于每一个seeker，计算其与其他线段的交点
                            auto seeker = seekers[ptr];
                            auto curSeg = rr.borders[n];
                            double cosBeta = (s.a * curSeg.a + s.b * curSeg.b)
                                / sqrt(s.a * s.a + s.b * s.b)
                                / sqrt(curSeg.a * curSeg.a + curSeg.b * curSeg.b);
                            if (abs(cosBeta) < 0.9) continue; // 两直线夹角不能大于10度
                            auto p = seeker.getCorWith(curSeg);
                            if (!p.isNULL) {
                                YFSegment tmp(seeker.startPoint, p);
                                if (tmp.distance < opLines[ptr].minDistance && tmp.distance > 0.1) {
                                    opLines[ptr].minDistance = tmp.distance;
                                    opLines[ptr].nearestLine = curSeg;
                                    opLines[ptr].keys = { k, n };
                                }
                            }
                        }
                    }
                }
                // 找到最近的线以及对应的索引之后，更新inLinesMap
                struct TMP_MID_LINE {
                    YFSegment midLine;
                    vector<int> keys;
                };
                vector<TMP_MID_LINE> midLines;
                if (opLines[0].nearestLine.isNULL && opLines[1].nearestLine.isNULL) {
                    continue; // 如果没找到对应中线，就不进行操作
                } else if (!opLines[0].nearestLine.isNULL && !opLines[1].nearestLine.isNULL) {
                    // 两个seeker都找到了对应中线，需要对中线进行修剪
                    vector<YFPoint> tmpPts({
                        opLines[0].nearestLine.startPoint,
                        opLines[0].nearestLine.endPoint,
                        opLines[1].nearestLine.startPoint,
                        opLines[1].nearestLine.endPoint,
                    });
                    sort(tmpPts.begin(), tmpPts.end(), compare);

                    double midX = (tmpPts[1].x + tmpPts[2].x) / 2;
                    double midY = (tmpPts[1].y + tmpPts[2].y) / 2;

                    auto m1 = computeMidLine(s, opLines[0].nearestLine);
                    auto m2 = computeMidLine(s, opLines[1].nearestLine);
                    vector<YFSegment> m({ m1, m2 });
                    sort(m.begin(), m.end(), compareLine);
                    bool hasReverse = !m[0].center.isEqualTo(m1.center); // 是否进行重排序了

                    auto computeAnotherPoint = [midX, midY](YFSegment s) {
                        double x;
                        double y;
                        if (abs(s.a) < MIN_ERR) {
                            x = midX;
                            y = (s.a * midX + s.c) / s.b;
                        } else {
                            x = (s.b * midY - s.c) / s.a;
                            y = midY;
                        }
                        return YFPoint(x, y);
                    };

                    auto m1sp = m[0].startPoint;
                    auto m1ep = computeAnotherPoint(m[0]);
                    auto m2sp = computeAnotherPoint(m[1]);
                    auto m2ep = m[1].endPoint;
                    midLines.push_back({ m1, opLines[hasReverse ? 1 : 0].keys });
                    midLines.push_back({ m2, opLines[hasReverse ? 0 : 1].keys });
                } else {
                    // 只有一个seeker找到了对应中线，直接返回对应中线即可
                    auto theLINE = opLines[0].nearestLine.isNULL ?
                        opLines[1] : opLines[0];
                    auto nearestLine = theLINE.nearestLine;
                    auto midLine = computeMidLine(s, nearestLine);
                    midLine.thickness = theLINE.minDistance;
                    midLine.startPoint.bulge = s.startPoint.bulge; // 这里的弧度信息需要判断一下

                    midLines.push_back({ midLine, theLINE.keys });
                }

                // 将处理好的中线存储到对应位置
                for (auto midLineWithKeys : midLines) {
                    auto keys = midLineWithKeys.keys;
                    auto midLine = midLineWithKeys.midLine;
                    if (keys.size() < 2) continue; // 防止越界
                    auto theKey = computeKeys(i, j, keys[0], keys[1]);

                    bool isInRegion = false; // 判断这条线是否在区域内
                    for (auto r : this->regions) {
                        isInRegion = isInRegion ||
                            midLine.startPoint.isInRegionWithoutBorder(r) ||
                            midLine.endPoint.isInRegionWithoutBorder(r) ||
                            midLine.center.isInRegionWithoutBorder(r);
                    }
                    if (isInRegion) continue;

                    LINE wmLine = { midLine, false };
                    auto curLineKey = compute2Key(i, j);
                    auto nearestLineKey = compute2Key(keys[0], keys[1]);

                    if (singleMap.find(curLineKey) == singleMap.end()) singleMap[curLineKey] = theKey; // 储存当前边线
                    if (singleMap.find(nearestLineKey) == singleMap.end()) singleMap[nearestLineKey] = theKey; // 储存对应边线

                    if (inlinesMap.find(theKey) == inlinesMap.end()) {
                        inlinesMap[theKey] = wmLine; // 两条线段之间只需要一条中线
                    }
                }
            }
        }

        for (int i = 0; i < count; i++) {
            auto r = allLines[i];
            for (int j = 0; j < r.borders.size(); j++) {
                int nextIndex = (j + 1) % r.borders.size();
                auto curLineKey = compute2Key(i, j);
                auto nextLineKey = compute2Key(i, nextIndex);
                auto curLine = r.borders.at(j);
                auto nextLine = r.borders.at(nextIndex);

                YFPoint linePoint; // 两条边缘的公共点
                if (curLine.endPoint.isEqualTo(nextLine.startPoint) || curLine.endPoint.isEqualTo(nextLine.endPoint)) {
                    linePoint = curLine.endPoint;
                } else linePoint = curLine.startPoint;

                if (singleMap.find(curLineKey) != singleMap.end() &&
                    singleMap.find(nextLineKey) != singleMap.end()) {
                    // 当两条线段的中线都存在时
                    auto curMidLineKey = singleMap[curLineKey];
                    auto nextMidLineKey = singleMap[nextLineKey];
                    auto curMidLine = inlinesMap[curMidLineKey];
                    auto nextMidLine = inlinesMap[nextMidLineKey];

                    // 开始削减或者延长对应的边
                    if (!curMidLine.s.isParalWith(nextMidLine.s)) {
                        // 只计算不平行的时候的情况
                        YFPoint corPoint = computeCorOfLines(curMidLine.s, nextMidLine.s); // 计算两条中线的交点
                        if (corPoint.isNULL) continue; // 没有交点，直接跳过
                        vector<LINE> tmpArr = { curMidLine, nextMidLine };
                        string key[2] = { curMidLineKey, nextMidLineKey };
                        for (int k = 0; k < tmpArr.size(); k++) {
                            auto innerLine = tmpArr.at(k).s;
                            bool hasChanged = tmpArr.at(k).hasChanged; // 是否已经对其削减过
                                                                       // 如果线段已经被削减过，那么只延长，不削减
                            if (((corPoint.x < innerLine.xRange.max && corPoint.x > innerLine.xRange.min) ||
                                (corPoint.y < innerLine.yRange.max && corPoint.y > innerLine.yRange.min)) && !hasChanged) {
                                // 如果点在直线的范围内，对其进行削减
                                double v1[2] = { innerLine.endPoint.x - corPoint.x, innerLine.endPoint.y - corPoint.y };
                                double v2[2] = { innerLine.startPoint.x - corPoint.x, innerLine.startPoint.y - corPoint.y };
                                YFPoint otherPoint;
                                if (pow(v1[0], 2) + pow(v1[1], 2) < pow(v2[0], 2) + pow(v2[1], 2)) otherPoint = innerLine.startPoint;
                                else otherPoint = innerLine.endPoint;
                                inlinesMap.at(key[k]) = LINE({ YFSegment(corPoint, otherPoint), true });
                            } else {
                                // 如果点在直线范围外，对其进行延长
                                vector<YFPoint> edgePoints = {
                                    innerLine.startPoint, innerLine.endPoint,
                                    corPoint
                                };
                                sort(edgePoints.begin(), edgePoints.end(), compare); // 将两条线的端点进行排序
                                YFSegment combineLine = YFSegment(edgePoints.at(0), edgePoints.at(2));
                                inlinesMap.at(key[k]) = LINE({ combineLine, true }); // 更新对应的中线
                            }
                        }
                    }
                }
            }
        }
        vector<YFSegment> inlines;
        for (auto it = inlinesMap.begin(); it != inlinesMap.end(); it++) {
            inlines.push_back(it->second.s);
        }
        return inlines;
    }

    YFRegion::YFRegion() {
        this->isNULL = true;
    }

    YFRegion::YFRegion(vector<YFSegment> s) {
        this->borders = s;
        this->center = this->findCenter(); // 需要更改，不能将私有属性暴露出来
        this->area = this->computeArea();
        this->perimeter = this->computePerimeter();
        this->isNULL = false;
    }

    /* 查找视觉中心位 */
    /* 查找策略，寻找最长切分线，切分线不能与边界有交点，而且中点在区域内
    * 取切分线中点作为视觉中心位
    */
    YFPoint YFRegion::findCenter() {
        vector<YFSegment> inLines;
        int borderNum = this->borders.size();
        for (int i = 0; i < int(borderNum / 2); i++) { // 开始遍历所有切分线
            for (int j = i + 1; j < borderNum; j++) {
                YFSegment s = YFSegment(
                    this->borders.at(i).startPoint,
                    this->borders.at(j).startPoint
                    );
                // 判断这条切分线是否在边线上
                // TIPs: 这里使算法复杂度上升到了o(n! * n)
                bool isInBorder = false;
                for (YFSegment seg : this->borders) {
                    if (seg.isParalWith(s)) {
                        isInBorder = true;
                        break;
                    };
                }
                if (isInBorder) continue;

                // 判断这条线是否与边线相交
                vector<YFPoint> corPoints = s.getCorWithRegion(*this);
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
            for (auto s : this->borders) {
                double x = s.startPoint.x;
                double y = s.startPoint.y;
                min_cx = x < min_cx ? x : min_cx;
                max_cx = x > max_cx ? x : max_cx;
                min_cy = y < min_cy ? y : min_cy;
                max_cy = y > max_cy ? y : max_cy;
            }
            return YFPoint((min_cx + max_cx) / 2, (min_cy + max_cy) / 2); // 没有符合条件的线
        }

        // 开始寻找最佳切分点
        YFPoint bestPoint;
        double maxRatio = 0;
        for (YFSegment seg : inLines) {
            // 计算线段横跨矩形的面积
            double l = abs(seg.xRange.max - seg.xRange.min); // 长
            double w = abs(seg.yRange.max - seg.yRange.min); // 宽
            double ratio = l * w;
            YFPoint center = seg.center; // 选取切分点的中点作为最佳视觉中心点
            bool isInRegion = center.isInRegion(*(this));
            if (ratio > maxRatio && isInRegion) {
                maxRatio = ratio;
                bestPoint = center;
            }
        }
        return bestPoint;
    }

    /* 计算区域面积 */
    double YFRegion::computeArea() {
        vector<YFPoint> points; // 区域的所有角点
        double area = 0;

        /* 从点集中删除点 */
        auto delPointFromPoints = [](YFPoint p, vector<YFPoint> points) {
            for (auto i = points.begin(); i != points.end(); i++) {
                if (p.isEqualTo(*i)) {
                    points.erase(i);
                    break;
                }
            }
            return points;
        };

        /* 判断某弓形是凸出去还是凹进来 */
        auto arcDirection = [](YFSegment curSeg, YFSegment nextSeg) {
            double vector_1[2] = { curSeg.endPoint.x - curSeg.startPoint.x, curSeg.endPoint.y - curSeg.startPoint.y }; // 向量化线段
            double vector_2[2] = { nextSeg.endPoint.x - nextSeg.startPoint.x, nextSeg.endPoint.y - nextSeg.startPoint.y };
            return vector_1[0] * vector_2[1] > vector_1[1] * vector_2[0]; // 如果为true，那么就是逆时针；如果为false，那么就是逆时针
        };

        /* 判断整个区域的顺逆方向 */
        int cw = 0; // 顺时针方向的角 clockwise
        int anticw = 0; // 逆时针方向的角 anticlockwise
        for (int i = 0; i < this->borders.size() - 1; i++) {
            if (arcDirection(this->borders.at(i), this->borders.at(i + 1))) {
                anticw++;
            } else cw++;
        }

        int arcDirect = anticw > cw ? 1 : -1; // 如果为1，那么总体为逆时针；如果为负，那么总体为-1

        double arcArea = 0; // 先计算带有弧边的面积

        for (int i = 0; i < this->borders.size(); i++) {
            YFSegment s = this->borders.at(i);
            points.push_back(s.startPoint);
            double sb = abs(s.startPoint.bulge);
            double eb = abs(s.endPoint.bulge);
            double b = sb > eb ? sb : eb; // 弧线反转之后，必有一个点的凸度不为0
            if (b > MIN_ERR) {
                int p = 1;
                if (sb > eb) {
                    p = arcDirect * (s.startPoint.bulge > 0 ? 1 : -1); // 如果弧的方向与总体方向一致，那么就应该加上弓形面积
                } else {
                    p = arcDirect * (s.endPoint.bulge > 0 ? -1 : 1); // 否则就减去弓形面积
                }
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
                YFPoint sp = points.at(i);
                YFPoint cp = points.at((i + 1) % lastsize);
                YFPoint ep = points.at((i + 2) % lastsize);
                YFSegment triLine = YFSegment(sp, ep); // 斜边
                vector<YFPoint> corPoints = triLine.getCorWithRegion(*this);
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
    double YFRegion::computePerimeter() {
        double perimeter = 0;
        for (auto l : this->borders) {
            double sb = abs(l.startPoint.bulge);
            double eb = abs(l.endPoint.bulge);
            double b = sb > eb ? sb : eb; // 弧线反转之后，必有一个点的凸度不为0
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

    /* 计算当前区域与指定区域的交点 */
    vector<YFPoint> YFRegion::getCorWithRegion(YFRegion r) {
        vector<YFPoint> regionCorPoints;
        for (YFSegment l : this->borders) {
            vector<YFPoint> lineCorPoints = l.getCorWithRegion(r);
            if (lineCorPoints.size() == 0) continue;
            regionCorPoints.insert(regionCorPoints.end(), lineCorPoints.begin(), lineCorPoints.end()); // 将交点保存起来
        }
        return regionCorPoints;
    }

    YFPoint::YFPoint() {
        this->isNULL = true;
    }

    YFPoint::YFPoint(double x_val, double y_val, double bulge_val, string id_val) {
        this->x = x_val;
        this->y = y_val;
        this->z = 0;
        this->bulge = bulge_val;
        this->id = id_val;
        this->isNULL = false;
    }

    YFPoint::YFPoint(double x_val, double y_val) {
        this->x = x_val;
        this->y = y_val;
        this->z = 0;
        this->bulge = 0;
        this->id = "No ID";
        this->isNULL = false;
    }

    /* 判断该点与另外一个点是否近似相等 */
    bool YFPoint::isEqualTo(YFPoint p) {
        return sqrt(pow(x - p.x, 2) + pow(y - p.y, 2)) < MIN_ERR;
    }

    /* 判断某点是否在某区域内 */
    bool YFPoint::isInRegion(YFRegion r) {
        YFPoint zeroPoint(-1, -1);
        YFSegment line = YFSegment(zeroPoint, *this); // 画一条射向区域外的射线
        vector<YFPoint> corPoints = line.getCorWithRegion(r); // 取得射线与区域的交点
        return corPoints.size() % 2 == 1; // 如果交点个数为奇数个，则判定该点在区域内
    }

    /* 判断点是否在不包含边界的区域内 */
    bool YFPoint::isInRegionWithoutBorder(YFRegion r) {
        bool isOnBorder = false;
        for (auto l : r.borders) {
            bool hasPoint = l.hasPoint(*(this));
            isOnBorder = isOnBorder || hasPoint;
        }
        if (isOnBorder) return false; // 如果在边界上，那么就直接判定为不在区域内
        YFPoint zeroPoint(-1, -1);
        YFSegment line = YFSegment(zeroPoint, *this); // 画一条射向区域外的射线
        vector<YFPoint> corPoints = line.getCorWithRegion(r); // 取得射线与区域的交点
        return corPoints.size() % 2 == 1; // 如果交点个数为奇数个，则判定该点在区域内
    }

    YFSegment::YFSegment(YFPoint sp, YFPoint ep, string id_val) {
        this->startPoint = sp;
        this->endPoint = ep;
        this->id = id_val;
        this->a = sp.y - ep.y;
        this->b = sp.x - ep.x;
        this->c = sp.x * ep.y - sp.y * ep.x;
        this->isNULL = false;
        this->distance = sqrt(pow(sp.x - ep.x, 2) + pow(sp.y - ep.y, 2)); // 计算长度
        this->xRange.min = sp.x < ep.x ? sp.x : ep.x;
        this->xRange.max = sp.x > ep.x ? sp.x : ep.x;
        this->yRange.min = sp.y < ep.y ? sp.y : ep.y;
        this->yRange.max = sp.y > ep.y ? sp.y : ep.y;
        this->center = YFPoint(
            (this->startPoint.x + this->endPoint.x) / 2,
            (this->startPoint.y + this->endPoint.y) / 2
            );
    }

    YFSegment::YFSegment(YFPoint sp, YFPoint ep) {
        this->startPoint = sp;
        this->endPoint = ep;
        this->id = "NO ID";
        this->a = sp.y - ep.y;
        this->b = sp.x - ep.x; // 直线的一般式竟然一直都写错了！
        this->c = sp.x * ep.y - sp.y * ep.x;
        this->isNULL = false;
        this->distance = sqrt(pow(sp.x - ep.x, 2) + pow(sp.y - ep.y, 2)); // 计算长度
        this->xRange.min = sp.x < ep.x ? sp.x : ep.x;
        this->xRange.max = sp.x > ep.x ? sp.x : ep.x;
        this->yRange.min = sp.y < ep.y ? sp.y : ep.y;
        this->yRange.max = sp.y > ep.y ? sp.y : ep.y;
        this->center = YFPoint(
            (this->startPoint.x + this->endPoint.x) / 2,
            (this->startPoint.y + this->endPoint.y) / 2
            );
    }

    YFSegment::YFSegment() {
        this->isNULL = true;
    }


    /* 判断是否与另一条线段平行 */
    bool YFSegment::isParalWith(YFSegment s) {
        return abs(a * s.b - b * s.a) < MIN_ERR;
    }

    /* 计算与另一条线段的交点 */
    YFPoint YFSegment::getCorWith(YFSegment s) {
        auto isInRange = [](double n, Range range) { // 用于判断某个值是否在范围内
            return n >= range.min - MIN_ERR && n <= range.max + MIN_ERR;
        };
        if (this->isParalWith(s)) {
            return YFPoint(); // 如果是平行的，就不存在交点
        }
        double der = a * s.b - b * s.a;
        double x = (b * s.c - c * s.b) / der;
        double y = (a * s.c - c * s.a) / der; // 计算出交点坐标值
        if (isInRange(x, xRange)
            && isInRange(y, yRange)
            && isInRange(x, s.xRange)
            && isInRange(y, s.yRange)) {
            YFPoint p = YFPoint(x, y);
            return p;
        } else {
            return YFPoint(); // 如果交点不在线段范围内，也不作数
        }
    }

    /* 计算线段与区域的交点 */
    vector<YFPoint> YFSegment::getCorWithRegion(YFRegion r) {
        vector<YFSegment> borders = r.borders;
        vector<YFPoint> corPoints; // 交点集合
        auto hasInSet = [](YFPoint p, vector<YFPoint> pset) {
            bool flag = false;
            for (YFPoint pi : pset) {
                flag = flag || p.isEqualTo(pi);
                if (flag) break;
            }
            return flag;
        };
        for (YFSegment s : borders) {
            YFPoint corPoint = this->getCorWith(s);
            YFPoint(1, 1);
            if (!corPoint.isNULL // 非空
                && !corPoint.isEqualTo(this->startPoint) // 不算线段的端点
                && !corPoint.isEqualTo(this->endPoint)
                && !hasInSet(corPoint, corPoints)) { // 避免重复添加
                corPoints.push_back(corPoint); // 将交点保存入集合
            }
        }
        return corPoints;
    }

    /* 判断某线段是否包含某点 */
    bool YFSegment::hasPoint(YFPoint p) {
        auto l = *(this);
        bool isOnLine = abs(l.a * p.x - l.b * p.y + l.c) < MIN_ERR;
        bool isInRange;
        if (abs(l.a) < MIN_ERR) {
            // 如果直线平行于x轴，就判定x的坐标范围
            isInRange = p.x <= l.xRange.max && p.x >= l.xRange.min;
        } else {
            // 如果直线平行于y轴，就判定y的坐标范围
            isInRange = p.y <= l.yRange.max && p.y >= l.yRange.min;
        }
        return isInRange && isOnLine;
    }
}

