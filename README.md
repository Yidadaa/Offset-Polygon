# 使用文档

主要包含`House`，`Segment`，`Region`和`Point`类，实际过程中，要先将墙壁转换为`Segment`类，并将其中的`startPoint`和`endPoint`属性包装成`Point`类，下面是一段示例代码，演示了从 json 文件中读取数据并执行操作的过程。

```C++
#pragma comment(lib, "lib_json.lib"); // jsoncpp的库，仅用于读取json文件
#include "House.h" // 引入Hosue.h头文件
#include <fstream>

int main() {
    ifstream ifs;
    ifs.open("data/data_5.json"); // 打开json文件

    CharReaderBuilder builder; // 调用jsoncpp的reader类，解析json内容。
    Value data;

    if (!ifs.is_open()) return;
    if (!parseFromStream(builder, ifs, &data, false)) return;

    /* 第一步：包装墙壁集合，用来初始化House类 */
    vector<Segment> wallLines; // 存放包装后的墙壁

    auto walls = data["walls"]; // 读取json中的墙壁数据

    for each (auto w in walls) {
        auto sp = w["startPoint"];
        auto ep = w["endPoint"];
        Point spp = Point(sp["x"].asDouble(), sp["y"].asDouble(), sp["bulge"].asDouble(), sp["Id"].asString()); // 将startPoint包装为Point类
        Point epp = Point(ep["x"].asDouble(), ep["y"].asDouble(), ep["bulge"].asDouble(), ep["Id"].asString()); // 将endPoint包装为Point类
        Segment s = Segment(spp, epp, w["Id"].asString()); // 用包装后的Point与Id初始化一条线段
        wallLines.push_back(s); // 将线段存入vector
    }

    /* 第二步：初始化House类 */
    House house(wallLines); // 用转换后的Segment集合构造House对象

    /* 第三步：开始使用*/
    vector<Region> regions = house.regions; // 找到的区域存放在House的regions属性中
    Region r = regions.at(0);

    Point center = r.center; // 区域视觉中心点
    double area = r.area; // 区域面积
    double perimeter = r.perimeter; // 区域周长

    return 0;
}

```